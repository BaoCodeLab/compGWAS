import glob
import pandas as pd
import numpy as np
import re
import os
import logging
import copy


def SNPgwas(args):

    if args.allSNP != None:
        SNPfiles = glob.glob(args.allSNP+r"/*")
    else:
        raise Exception("Please provide SNP files path!")

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
   
    if args.cols != None:
        ctigindx = args.cols[0]
        locindx = args.cols[1]
        aleindx = args.cols[2]
    else:
        raise Exception("Please provide the column index of Loc and Allele of the SNPs annotation file!")

    if args.thread != None:
        thread = args.thread

    if args.threshold != None:
        threshold = args.threshold

    if args.outdir != None:
        outdir = args.outdir
    
    if args.prefix != None:
        prefix = args.prefix

    rootout = outdir + "/" + prefix
    out = rootout + ".information"
    out2 = rootout + ".LogistRegre.all.CHI2"
    out3 = rootout + ".LogistRegre.all"
    out4 = rootout + ".LogistRegre.filtered"

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


    data,count,names = {},1,[]
    for i in range(0,len(SNPfiles)):
        file_i = SNPfiles[i]
        file_nm = file_i.split("/")[-1].split("GCA-")[-1].split(".")[0]
        names.append(file_nm)
        i_file = pd.read_csv(file_i, "\t")
        indxi = list(map(lambda x,y:str(x)+"_!T!E!S!T!_"+str(y), list(i_file.iloc[:,ctigindx]), list(i_file.iloc[:,locindx])))
        data[count] = pd.Series(list(i_file.iloc[:,aleindx]), index=indxi)
        data[count] = data[count].astype('category')
        data[count] = data[count].cat.add_categories([''])
        count=count+1

    logging.info('The account of samples is:\t %d',count-1)
    dataF = pd.DataFrame(data)
    dataF2 = pd.DataFrame(data)
    dataF.columns = ['ALT']*int(count-1)
    dataF.insert(0, 'Loc', dataF.index)
    dataF = dataF.fillna('')
    dataF2 = dataF2.fillna('')
    
    data2 = dataF["ALT"]
    data2["REF"] = pd.Series([""]*len(data2), index=dataF.index).astype('category')
    dataF2["REF"] = pd.Series([""]*len(dataF2), index=dataF2.index).astype('category')
    
    allCunt = data2.apply(pd.value_counts, axis=1)
    NaNCunt = allCunt.isnull().sum(axis=1)
    allCuntMin = allCunt.min(axis=1)
    alleFilt1 = NaNCunt[NaNCunt<3].index.tolist()
    alleFilt2 = allCuntMin[allCuntMin<2].index.tolist()
    alleFilt = list(set(alleFilt1 + alleFilt2))

    dataF = dataF.drop(alleFilt)
    dataF2 = dataF2.drop(alleFilt)
    allCunt = pd.DataFrame(allCunt.drop(alleFilt))
    majAlle = allCunt.idxmax(axis=1)

    dataF2.insert(dataF2.columns.size, 'major', majAlle)
    data3 = pd.DataFrame()
    for i in range(0, len(SNPfiles)):
        file_i = SNPfiles[i]
        file_nm = file_i.split("/")[-1].split("GCA-")[-1].split(".")[0]
        data3.insert(i, file_nm, np.where((dataF2["major"] == dataF2[i+1]), 1, dataF2[i+1]))
    if args.ref != None:
        data3.insert(data3.columns.size,str(args.ref[0]),np.where((dataF2["major"] == dataF2["REF"]), 1, dataF2["REF"]))

    data3 = pd.DataFrame(data3)
    data3 = data3.replace('',2)
    data3 = data3.replace('A',2)
    data3 = data3.replace('T',2)
    data3 = data3.replace('C',2)
    data3 = data3.replace('G',2)
    data3 = data3.astype("int8")
    data3.insert(0, "Loc", dataF["Loc"].values.tolist())
    data3.index = dataF2.index.tolist()
    geno = copy.copy(data3) 
    

    input1 = open(pheno1,"r")
    input2 = open(pheno0,"r")
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
    nrow = len(data3)
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
                tab_i[m] = pd.Series([cols_i[m]]*nrow, index=data3.index)
            if(GCA_i in GCAs.keys()):
                tab_i["pheno"] = pd.Series([GCAs[GCA_i]] * nrow, index=data3.index)
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
            tab_i["pheno"] = pd.Series([GCAs[m]] * nrow, index=data3.index)    
            data_i = pd.DataFrame(tab_i)
            data_i.columns = ["Phenotype"]
            tab[m] = data_i
        Rcols = 'Allele,Phenotype'
    final = []
    for i in tab.keys():
        tab[i].insert(0, 'Allele', data3[i].values.tolist())
        final.append(tab[i])
    finaltab = pd.concat(final, axis=1)
    finalindx = data3.index.tolist()
    finalctig = list(map(lambda x:x.split("_!T!E!S!T!_")[0], finalindx))
    finalloc = list(map(lambda x:int(x.split("_!T!E!S!T!_")[1]), finalindx))
    finaltab.insert(0, "Loc", finalloc)
    finaltab.insert(0, "RefName", finalctig)
    finaltab = finaltab.sort_values(by=["RefName","Loc"],ascending=(True,True))
    header_str = '\t'.join(finaltab.columns.values)
    np.savetxt(out, finaltab.values, delimiter="\t", fmt='%s', header=header_str)
    input1.close()
    input2.close()
    
    
    if args.ref != None:
        names.append(str(args.ref[0]))
    newname = ["Loc"]
    for i in range(0,len(names)):
        GCA_i = names[i]
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
    Locs = AllCount.index.tolist()
    Ctig = list(map(lambda x:x.split("_!T!E!S!T!_")[0], Locs))
    LOC = list(map(lambda x:int(x.split("_!T!E!S!T!_")[1]), Locs))
    AllCount.insert(0,"Loc",LOC)
    AllCount.insert(0,"RefName",Ctig)
    colnames = AllCount.columns.tolist()
    AllCount = AllCount.sort_values(by=["RefName","Loc"],ascending=(True,True))
    header_str = '\t'.join(AllCount.columns.values)
    np.savetxt(out3, AllCount.values ,delimiter="\t" ,fmt='%s',header= header_str)

    os.system('%s %s %s %s %s %s %s SNP %d %s' % (Rscr,script1,out3,colnames[2],colnames[3],colnames[4],colnames[5],thread,out2))

    os.system('%s %s %s %s %s SNP %d %f %s %s' % (Rscr,script2,out,Rcols,types,thread,threshold,out3,out4)) 



