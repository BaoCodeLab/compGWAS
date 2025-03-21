import argparse
import pandas as pd
import numpy as np
import re
import pickle
import sys
import glob
import logging
import os
import copy


class StrToBytes:
    def __init__(self,fileobj):
        self.fileobj = fileobj
    def read(self,size):
        return self.fileobj.read(size).encode()
    def readline(self,size=-1):
        return self.fileobj.readline(size).encode()


def nonCDSgwas(args):
    if args.allDIP != None:
        DIPfiles = glob.glob(args.allDIP+r"/*")
    else:
        raise Exception("Please provide DIP files path!")

    if args.CDSdic != None:
        CDSdic = args.CDSdic
    else:
        raise Exception("Please provide CDS annotation dictionaries of the species to be analyzed!")

    if args.size != None:
        size = int(args.size)
    else:
        raise Exception("Please provide the size of the window!")

    if args.pheno1 != None:
        pheno1 = args.pheno1[0]
        phe1 = args.pheno1[1]
    else:
        raise Exception("Please provide phenotype1's GCA number file!")

    if args.pheno0 != None:
        pheno0 = args.pheno0[0]
        phe0 = args.pheno0[1]
    else:
        raise Exception("Please provide phenotype2's GCA number file!")

    if args.genome != None:
        genome = args.genome
    else:
        raise Exception("Please provide the FASTA format genome file of the species to be analyzed!")

    if args.cols != None:
        ctigindx = args.cols[0]
        locindx = args.cols[1]
        refindx = args.cols[2]
        muttpindx = args.cols[3]
    else:
        raise Exception("Please provide the column index of Loc, Ref, and Muttype of the mutation annotation file!")

    if args.thread != None:
        thread = args.thread

    if args.threshold != None:
        threshold = args.threshold

    if args.outdir != None:
        outdir = args.outdir

    if args.prefix != None:
        prefix = args.prefix

    rootout = outdir + "/" + prefix
    out = rootout + ".window.all"
    out1 = rootout + ".window.filtered"
    out2 = rootout + ".window.filtered.information"
    out3 = rootout + ".LogistRegre.all.CHI2"
    out4 = rootout + ".LogistRegre.all"
    out5 = rootout + ".LogistRegre.filtered"

    if args.Rscript != None:
        Rscr = args.Rscript
    else:
        raise Exception("Please provide the path of Rscript!")

    if args.rscript != None:
        script1 = args.rscript + "allGWAS/GWASlib/CHI2distribution.R"
        script2 = args.rscript + "allGWAS/GWASlib/GWAS.R"
    else:
        raise Exception("Please provide the path of GWAS package!")

    logging.basicConfig(
            format='%(asctime)s \tFile \"%(filename)s" %(levelname)s: \n %(message)s \n',
            datefmt='%a, %d %b %Y %H:%M:%S',
            filename=rootout+'.log', filemode="w", level=logging.INFO)


    anno_dic = {}
    for m in range(0,len(CDSdic)):
        anno_m = open(CDSdic[m],"r")
        anno_m_dic = pickle.load(StrToBytes(anno_m))
        anno_dic.update(anno_m_dic)
        anno_m.close()
    contigs = []
    for i in anno_dic.keys():
        contigs.append(i[0])
    contigs = list(set(contigs))

    genome_file = open(genome,"r")
    genome_tmp = genome_file.read()
    CtigLendic = {}
    for i in contigs:
        tag = i
        re_genome = re.compile(">"+tag+" "+ '.+\n' + '((\w+\n)+)')
        genome_new = list(re_genome.findall(genome_tmp)[0][0])
        while "\n" in genome_new:
            genome_new.remove("\n")
        CtigLendic[tag] = len(genome_new)

    ctig_nonCDS_dic = {}
    for jj in contigs:
        newdic = {}
        for i in anno_dic.keys():
            if i[0] == jj:
                start_i = int(anno_dic[i][0][0])
                newdic[start_i] = i
        keys = sorted(list(newdic.keys()))

        # filter the start site of gene whose both start site and end site are in an another gene
        # creat a dic (anno_dic_2) which contains all CDS regions
        anno_dic_2,filter_start = {},[]

        for mm in keys:
            i = newdic[mm]
            start_i = int(anno_dic[i][0][0])
            end_i = int(anno_dic[i][0][1])
            for n in keys:
                m = newdic[n]
                start_m = int(anno_dic[m][0][0])
                end_m = int(anno_dic[m][0][1])
                if (start_i < start_m and end_i > end_m):
                    filter_start.append(int(start_m))
                if (start_i < start_m and end_i < end_m and end_i >= start_m):
                    filter_start.append(int(start_m))
                    end_i = end_m
                if (start_i > start_m and start_i <= end_m and end_i >= end_m):
                    filter_start.append(int(start_i))
                    start_i = start_m
                if (start_i > start_m and end_i < end_m):
                    filter_start.append(int(start_i))
                    start_i = start_m
                    end_i = end_m

                if (start_i == start_m and end_i == end_m and anno_dic[m] != anno_dic[i]):
                    logging.info('There is an another gene {1}, which possesses the same start and end site with the gene {0}.'.format(i[1],m[1]))
                if (start_i == start_m and end_i < end_m and end_i >= start_m):
                    logging.info('There is an another gene {1}, whose start site is same but end site is bigger than the gene {0}.'.format(i[1],m[1]))
                if (start_i > start_m and start_i <= end_m and end_i == end_m):
                    logging.info('There is an another gene {1}, whose end site is same but start site is smaller than the gene {0}.'.format(i[1],m[1]))
            if (start_i not in filter_start):
                anno_dic_2[start_i] = [start_i,end_i]

        dic2_keys = sorted(list(anno_dic_2.keys()))

        # creat a dic (nonCDS_anno_dic) which contains all nonCDS regions
        nonCDS_anno_dic = {}
        nonCDS_anno_dic[1] = [1,dic2_keys[0]-1]
        for i in range(0,len(dic2_keys)-1):
            left = anno_dic_2[dic2_keys[i]][1]+1
            right = anno_dic_2[dic2_keys[i+1]][0]-1
            nonCDS_anno_dic[left] = [left,right]
        if (anno_dic_2[dic2_keys[-1]][1] < CtigLendic[jj]):
            nonCDS_anno_dic[anno_dic_2[dic2_keys[-1]][1]+1] = [anno_dic_2[dic2_keys[-1]][1]+1,CtigLendic[jj]]

        del anno_dic_2
        del anno_dic

        dataF=pd.DataFrame(nonCDS_anno_dic)
        nonCDSdata = pd.DataFrame(dataF.values.T, index=dataF.columns, columns=dataF.index)
        del nonCDS_anno_dic
        del dataF
        ctig_nonCDS_dic[jj] = nonCDSdata

    # merge all DIP sites of all DIP files
    data,count = {},1
    for i in range(0,len(DIPfiles)):
        file_nm = DIPfiles[i]
        i_file = pd.read_csv(file_nm, "\t")
        indxi = list(map(lambda x,y:str(x)+"_!T!E!S!T!_"+str(y), list(i_file.iloc[:,ctigindx]), list(i_file.iloc[:,locindx])))
        data[count] = pd.Series(list(i_file.iloc[:,refindx]), index=indxi)
        data[count] = data[count].astype('category')
        data[count] = data[count].cat.add_categories(['',"O","P"])
        count=count+1
    dataF=pd.DataFrame(data)
    dataF.columns = ["ALT"] * dataF.columns.size
    dataF.insert(0,'Loc',dataF.index)
    del data
    if args.ref != None:
        dataF["REF"]= np.nan

    # filter CDS indel sites
    data,colnames = {},[]
    for i in range(0,len(DIPfiles)):
        file_i = DIPfiles[i]
        file_nm = file_i.split("/")[-1].split("GCA-")[-1].split(".")[0]
        colnames.append(str(file_nm))
        i_file = pd.read_csv(file_i, "\t")
        indxi = list(map(lambda x,y:str(x)+"_!T!E!S!T!_"+str(y), list(i_file.iloc[:,ctigindx]), list(i_file.iloc[:,locindx])))
        data[file_nm] = pd.Series(list(i_file.iloc[:,muttpindx]), index=indxi)
        data[file_nm] = data[file_nm].astype('category')
        data[file_nm] = data[file_nm].cat.add_categories([''])
    MutType_data=pd.DataFrame(data)
    MutType_data.columns = colnames
    MutType_data.insert(0,'Loc',MutType_data.index)
    MutType_data = MutType_data.fillna('')

    nonCDS_mut_types = ['inter-gene-3-5','inter-gene','pseudo-gene|ncRNA','inter-gene-5-3','inter-gene-5-5','inter-gene-3-3','inter-gene--5','inter-gene--3','inter-gene-5-','inter-gene-3-']
    cds_filter_index = []
    for i in range(0,len(MutType_data)):
        mut_types_i = list(set(MutType_data.iloc[i,:].tolist()))
        inter_nonCDS_i = list(set(mut_types_i).intersection(set(nonCDS_mut_types)))
        if(len(inter_nonCDS_i)==0):
            cds_filter_index.append(MutType_data.iloc[i,0])
    dataF = dataF.drop(cds_filter_index)
    del MutType_data

    # calculate the MAF of each window
    dataF_colnames = dataF.columns.tolist()
    dataF.columns = list(range(0,dataF.shape[1]))
    dataF.index = list(range(0,len(dataF)))
    dataF=dataF.fillna("O")

    for m in dataF.index.tolist():
        for n in range(1,dataF.shape[1]):
            if(dataF.loc[m,n] != "O" ):
                dataF.loc[m,n] = "P"


    dataF1 = dataF.iloc[:,1:]
    allCount=dataF1.apply(pd.value_counts,axis=1)
    allCount = pd.DataFrame(allCount)
    del dataF1
    minor_allele = allCount.idxmin(axis=1)
    dataF.insert(dataF.shape[1],'Minor',minor_allele)
    for i in range(1,dataF.shape[1]-1):
        dataF[i] = np.where((dataF["Minor"] == dataF[i]), np.nan, dataF[i]) 
    dataF = dataF.iloc[:,0:dataF.shape[1]-1]
    dataF.index = dataF.iloc[:,0]
    dataF.columns = dataF_colnames
    INDEX = dataF.iloc[:,0].tolist()
    Ctig = list(map(lambda x:x.split("_!T!E!S!T!_")[0], INDEX))
    LOC = list(map(lambda x:int(x.split("_!T!E!S!T!_")[1]), INDEX))
    dataFF = dataF.drop(columns=["Loc"])
    dataFF.insert(0, "Loc", LOC)
    dataFF.insert(0, "RefName", Ctig)

    colnames = ["RefName","Loc","Sub"]
    for i in range(0,len(DIPfiles)):
        file_i = DIPfiles[i]
        file_nm = file_i.split("/")[-1].split("_GCA-")[-1].split(".")[0]
        colnames.append(str(file_nm))
    if args.ref != None:
        colnames.append(str(args.ref[0]))

    OUT,indx,allele_filter = [],0,[]
    for jj in contigs: 
        min_size = 1
        max_size = int(size)
        step = int(max_size/10)
        DATA,locs,subs = {},[],[]
    
        dataF2 = dataFF[dataFF.iloc[:,0] == jj]
        nonCDSdata = ctig_nonCDS_dic[jj]
        for i in range(1,dataF2.iloc[:,1].max()):
            nonCDSdata_i = nonCDSdata[((min_size <= nonCDSdata.iloc[:,0]) & (nonCDSdata.iloc[:,0] <= max_size)) | ((min_size <= nonCDSdata.iloc[:,1]) & (nonCDSdata.iloc[:,1] <= max_size)) | ((min_size >= nonCDSdata.iloc[:,0]) & (nonCDSdata.iloc[:,1] >= max_size)) ]
            if len(nonCDSdata_i) == 0 or len(nonCDSdata_i) == 1:
                dataF_i = dataF2[(min_size <= dataF2.iloc[:,1]) & (dataF2.iloc[:,1] <= max_size)]
                dataF_i = dataF_i.iloc[:,1:]
                sites_num = len(dataF_i)
                results_i = []
                if(sites_num==0):
                    min_size = min_size + step
                    max_size = max_size + step
                    continue
                else:
                    na_num_list = dataF_i.isna().sum().tolist()
                    major_num = 0               
                    for m in range(1,len(na_num_list)):     
                        if(na_num_list[m]==0):
                            major_num=major_num+1
                    minor_num = int(dataF_i.shape[1]-1)-major_num    
                    if(min(major_num,minor_num)<2):
                        allele_filter.append(indx)
                    for n in range(1,len(na_num_list)):     
                        if((na_num_list[n]!=0 and major_num>=minor_num) or (na_num_list[n]==0 and major_num<minor_num)):
                            results_i.append(2)
                        else:
                            results_i.append(1)
                DATA[indx] = pd.Series(results_i, index=range(0,len(results_i)))
                DATA[indx] = DATA[indx].astype('category')
                locs.append(min_size)
                subs.append(0)
                min_size = min_size + step
                max_size = max_size + step
                indx = indx + 1
                if(min_size > dataF2.iloc[:,1].max()):
                    break
            elif len(nonCDSdata_i) > 1:
                for ii in range(0,len(nonCDSdata_i)):
                    if (nonCDSdata_i.iloc[ii,0] >= min_size and nonCDSdata_i.iloc[ii,1] <= max_size):
                        dataF_i = dataF2[(nonCDSdata_i.iloc[ii,0] <= dataF2.iloc[:,1]) & (dataF2.iloc[:,1] <= nonCDSdata_i.iloc[ii,1])]
                    elif (nonCDSdata_i.iloc[ii,0] < min_size and nonCDSdata_i.iloc[ii,1] <= max_size):
                        dataF_i = dataF2[(min_size <= dataF2.iloc[:,1]) & (dataF2.iloc[:,1] <= nonCDSdata_i.iloc[ii,1])]
                    elif (nonCDSdata_i.iloc[ii,0] >= min_size and nonCDSdata_i.iloc[ii,1] > max_size):
                        dataF_i = dataF2[(nonCDSdata_i.iloc[ii,0] <= dataF2.iloc[:,1]) & (dataF2.iloc[:,1] <= max_size)]
                    elif (nonCDSdata_i.iloc[ii,0] < min_size and nonCDSdata_i.iloc[ii,1] > max_size):
                        dataF_i = dataF2[(min_size <= dataF2.iloc[:,1]) & (dataF2.iloc[:,1] <= max_size)]
                    dataF_i = dataF_i.iloc[:,1:]
            
                    sites_num = len(dataF_i)
                    results_i = []
                    if(sites_num==0):
                        continue
                    else:
                        na_num_list = dataF_i.isna().sum().tolist()
                        major_num = 0                 
                        for m in range(1,len(na_num_list)):     
                            if(na_num_list[m]==0):
                                major_num=major_num+1
                        minor_num = int(dataF_i.shape[1]-1)-major_num    
                        if(min(major_num,minor_num)<2):
                            allele_filter.append(indx)
                        for n in range(1,len(na_num_list)):      
                            if((na_num_list[n]!=0 and major_num>=minor_num) or (na_num_list[n]==0 and major_num<minor_num)):
                                results_i.append(2)
                            else:
                                results_i.append(1)
                    DATA[indx] = pd.Series(results_i, index=range(0,len(results_i)))
                    DATA[indx] = DATA[indx].astype('category')
                    locs.append(min_size)
                    subs.append(ii+1)
                    indx = indx + 1
                min_size = min_size + step
                max_size = max_size + step
                if(min_size > dataF2.iloc[:,1].max()):
                    break

        data=pd.DataFrame(DATA)
        data2 = pd.DataFrame(data.values.T, index=data.columns, columns=data.index)
        data2.insert(0,"Sub",subs)
        data2.insert(0,"Loc",locs)
        data2.insert(0,"RefName",[jj]*len(locs))
        data2.columns = colnames
        OUT.append(data2)    
    
    data2 = pd.concat(OUT)
    header_str = '\t'.join(data2.columns.values)
    np.savetxt(out, data2.values ,delimiter="\t" ,fmt='%s',header= header_str)
    data2 = data2.drop(allele_filter)
    header_str = '\t'.join(data2.columns.values)
    np.savetxt(out1, data2.values ,delimiter="\t" ,fmt='%s',header= header_str)
    geno = copy.copy(data2)

    # get the information table for LogisticRegre
    input1 =  open(pheno1,"r")
    input2 =  open(pheno0,"r")
    GCAs_pheno1 = input1.readlines()
    GCAs_pheno0 = input2.readlines()
    GCAs,PH1_GCAs,PH0_GCAs = {},[],[]
    for i in GCAs_pheno1:
        GCAs[i.split("\n")[0]] = 1
        PH1_GCAs.append(i.split("\n")[0])
    for i in GCAs_pheno0:
        GCAs[i.split("\n")[0]] = 0
        PH0_GCAs.append(i.split("\n")[0])
    if args.ref != None:
        if args.ref[1] == "p":
            GCAs[str(args.ref[0])] = 1
            PH1_GCAs.append(str(args.ref[0]))
        elif args.ref[1] == "P":
            GCAs[str(args.ref[0])] = 0
            PH0_GCAs.append(str(args.ref[0]))
        else:
            raise Exception("Please provide right phenotype of reference strain, p(means phenotype1) or P(means phenotype0)!")

    nrow = len(data2)
    re_GCA = re.compile(r'GCA_(\d+)')
    tab = {}
    if(args.info != None):  
        infotab = args.info[0]
        covar = args.info[1]
        types = args.info[2]
        covar = covar.split(",")
        covar.append("Phenotype")
        input4 = open(infotab, "r")
        infor_tab = input4.readlines()
        cols = infor_tab[0].split("\n")[0].split("\t")
        for i in range(1,len(infor_tab)):
            cols_i = infor_tab[i].split("\n")[0].split("\t")
            GCA_i = re_GCA.findall(cols_i[0])[0]
            tab_i = {}
            for m in range(1, len(cols)):
                tab_i[m] = pd.Series([cols_i[m]]*nrow, index=data2.index)
            if(GCA_i in GCAs.keys()):
                tab_i["pheno"] = pd.Series([GCAs[GCA_i]] * nrow, index=data2.index)
                data_i = pd.DataFrame(tab_i)
                data_i.columns = cols[1:]+["Phenotype"]
                data_i = data_i[covar]
                tab[GCA_i] = data_i
        input4.close()
        Rcols = ",".join(['Allele']+covar)
    else:
        types = "N"
        logging.info('There is no covariate for logistic regression.')
        for m in GCAs.keys():
            tab_i = {}
            tab_i["pheno"] = pd.Series([GCAs[m]] * nrow, index=data2.index)
            data_i = pd.DataFrame(tab_i)
            data_i.columns = ["Phenotype"]
            tab[m] = data_i
        Rcols = 'Allele,Phenotype'
    final = []
    for i in tab.keys():
         tab[i].insert(0, 'Allele', data2[i].values.tolist())
         final.append(tab[i])
    finaltab = pd.concat(final, axis=1)
    RefNmae = data2.iloc[:,0].values.tolist()
    Loc = data2.iloc[:,1].values.tolist()
    Sub = data2.iloc[:,2].values.tolist()
    del data2
    finaltab.insert(0,"Sub",Sub)
    finaltab.insert(0,"Loc",Loc)
    finaltab.insert(0,"RefNmae",RefNmae)
    header_str = '\t'.join(finaltab.columns.values)
    np.savetxt(out2, finaltab.values ,delimiter="\t" ,fmt='%s',header= header_str)
    input1.close()
    input2.close()

    
    oldname = geno.columns.tolist()
    newname = ["RefName","Loc","Sub"]
    for i in range(3,len(oldname)):
        GCA_i =  oldname[i]
        if(GCA_i in PH1_GCAs):
            newname.append(phe1)
        elif(GCA_i in PH0_GCAs):
            newname.append(phe0)
    geno.columns = newname
    geno_phe1,geno_phe0 = geno[phe1],geno[phe0]
    phe1allCount = geno_phe1.apply(pd.value_counts,axis=1)
    phe1allCount = pd.DataFrame(phe1allCount)
    phe1_cols=[]
    for i in range(0,phe1allCount.shape[1]):
        phe1_cols.append(phe1+"."+str(phe1allCount.columns.tolist()[i]))
    phe1allCount.columns = phe1_cols
    phe1allCount = phe1allCount.fillna("0")
    phe1allCount = pd.DataFrame(phe1allCount,dtype="int")
    
    phe0allCount=geno_phe0.apply(pd.value_counts,axis=1)
    phe0allCount = pd.DataFrame(phe0allCount)
    phe0_cols=[]
    for i in range(0,phe0allCount.shape[1]):
        phe0_cols.append(phe0+"."+str(phe0allCount.columns.tolist()[i]))
    phe0allCount.columns = phe0_cols
    phe0allCount = phe0allCount.fillna("0")
    phe0allCount = pd.DataFrame(phe0allCount,dtype="int")
    
    AllCount = pd.concat([phe1allCount,phe0allCount],axis=1)
    AllCount = pd.DataFrame(AllCount)
    AllCount.insert(0,"Sub",geno.iloc[:,2].tolist())
    AllCount.insert(0,"Loc",geno.iloc[:,1].tolist())
    AllCount.insert(0,"RefName",geno.iloc[:,0].tolist())
    colnames = AllCount.columns.tolist()
    header_str = '\t'.join(AllCount.columns.values)
    np.savetxt(out4, AllCount.values ,delimiter="\t" ,fmt='%s',header= header_str)

    os.system('%s %s %s %s %s %s %s nonCDS %d %s' % (Rscr,script1,out4,colnames[3],colnames[4],colnames[5],colnames[6],thread,out3))

    os.system('%s %s %s %s %s nonCDS %d %f %s %s' % (Rscr,script2,out2,Rcols,types,thread,threshold,out4,out5))





