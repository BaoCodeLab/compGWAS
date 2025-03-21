import glob
import pandas as pd
import numpy as np
import re
import os
import pickle
import copy
import logging


def merge_DIPorSNP_ALT(allfiles,ALT_col,Loc_col,RefName_col):
    data,colnames = {},[]
    for i in range(0,len(allfiles)):
        file_i = allfiles[i]
        file_nm = file_i.split("/")[-1].split("GCA-")[-1].split(".")[0]   
        colnames.append(str(file_nm))
        i_file = pd.read_csv(file_i, "\t")
        indxi = list(map(lambda x,y:str(x)+"_!T!E!S!T!_"+str(y), list(i_file.iloc[:,RefName_col]), list(i_file.iloc[:,Loc_col])))
        data[file_nm] = pd.Series(list(i_file.iloc[:,ALT_col]), index=indxi)
        data[file_nm] = data[file_nm].astype('category')
        data[file_nm] = data[file_nm].cat.add_categories([''])
    dataF=pd.DataFrame(data)
    del data
    dataF.columns = colnames
    indxi = dataF.index.tolist()
    Ctig = list(map(lambda x:x.split("_!T!E!S!T!_")[0], indxi))
    LOC = list(map(lambda x:int(x.split("_!T!E!S!T!_")[1]), indxi))
    dataF.insert(0,"Loc",LOC)
    dataF.insert(0,"RefName",Ctig)
    return dataF


class StrToBytes:
    def __init__(self,fileobj):
        self.fileobj = fileobj
    def read(self,size):
        return self.fileobj.read(size).encode()
    def readline(self,size=-1):
        return self.fileobj.readline(size).encode()


def load_anno_dic(anno_dics_list):
    anno_dic = {}
    for m in range(0,len(anno_dics_list)):
        anno_m = open(anno_dics_list[m],"r")
        anno_m_dic = pickle.load(StrToBytes(anno_m))
        anno_dic.update(anno_m_dic)
        anno_m.close()
    return anno_dic


def get_DIPorSNP_geneID(merged_DIP_dataframe,merged_SNP_dataframe,*CDSandPseudo_dic):
    merged_DIP_dataframe.insert(0,"geneID",[""]*len(merged_DIP_dataframe))
    merged_SNP_dataframe.insert(0, "geneID", [""] * len(merged_SNP_dataframe))
    fisrt_filter_geneID = []
    mutated_cds_count = 0
    highmutated_cds_count = 0
    highmutated_cds_dataframe_list = []
    nonmutated_cds_count = 0
    nonmutated_geneID = []
    for i in CDSandPseudo_dic:
        for m in i.keys():
            CDS_contig = m[0]
            CDS_start = int(i[m][0][0])
            CDS_end = int(i[m][0][1])
            CDS_length = CDS_end - CDS_start + 1
            merged_DIP_dataframe_m = merged_DIP_dataframe[(merged_DIP_dataframe.iloc[:, 1] <= CDS_contig) & (CDS_start <= merged_DIP_dataframe.iloc[:, 2]) & (merged_DIP_dataframe.iloc[:, 2] <= CDS_end)]
            merged_SNP_dataframe_m = merged_SNP_dataframe[(merged_SNP_dataframe.iloc[:, 1] <= CDS_contig) & (CDS_start <= merged_SNP_dataframe.iloc[:, 2]) & (merged_SNP_dataframe.iloc[:, 2] <= CDS_end)]
            if ( len(merged_DIP_dataframe_m) == 0 and len(merged_SNP_dataframe_m)==0 ):
                nonmutated_cds_count = nonmutated_cds_count + 1
                nonmutated_geneID.append(m[1])
            else:
                mutated_cds_count = mutated_cds_count + 1
                merged_DIP_dataframe_m_m = merged_DIP_dataframe_m.iloc[:, 3:]
                merged_SNP_dataframe_m_m = merged_SNP_dataframe_m.iloc[:, 3:]
                merged_DIPSNP_dataframe_m = pd.concat([merged_DIP_dataframe_m_m, merged_SNP_dataframe_m_m], axis=1)
                NaN_count = merged_DIPSNP_dataframe_m.isnull().sum(axis=1)
                mut_ratio = (  ( len(NaN_count) * merged_DIPSNP_dataframe_m.shape[1] - sum(NaN_count.values.tolist()) ) / (merged_DIPSNP_dataframe_m.shape[1]/2)  ) /  CDS_length
                if ( mut_ratio > 0.15 ):
                    highmutated_cds_count = highmutated_cds_count + 1
                    merged_DIPSNP_dataframe_m_indx = merged_DIPSNP_dataframe_m.index.tolist()
                    merged_DIPSNP_dataframe_m_ctig = list(map(lambda x:x.split("_!T!E!S!T!_")[0], merged_DIPSNP_dataframe_m_indx))
                    merged_DIPSNP_dataframe_m_loc = list(map(lambda x:int(x.split("_!T!E!S!T!_")[1]), merged_DIPSNP_dataframe_m_indx))
                    merged_DIPSNP_dataframe_m.insert(0, "Loc", merged_DIPSNP_dataframe_m_loc)
                    merged_DIPSNP_dataframe_m.insert(0, "RefName", merged_DIPSNP_dataframe_m_ctig)
                    merged_DIPSNP_dataframe_m.insert(0, "geneID", [""] * len(merged_DIPSNP_dataframe_m))
                    for n in merged_DIPSNP_dataframe_m.index.tolist():
                        merged_DIPSNP_dataframe_m.loc[n, "geneID"] = m[1]
                    highmutated_cds_dataframe_list.append(merged_DIPSNP_dataframe_m)
                else:
                    fisrt_filter_geneID.append(m[1])
                    for n in merged_DIP_dataframe_m.index.tolist():
                        if (merged_DIP_dataframe.loc[n, "geneID"] == ''):
                            merged_DIP_dataframe.loc[n, "geneID"] = m[1]
                        else:
                            merged_DIP_dataframe.loc[n, "geneID"] = merged_DIP_dataframe.loc[n, "geneID"] + ";"+ m[1]
                    for n in merged_SNP_dataframe_m.index.tolist():
                        if (merged_SNP_dataframe.loc[n, "geneID"] == ''):
                            merged_SNP_dataframe.loc[n, "geneID"] = m[1]
                        else:
                            merged_SNP_dataframe.loc[n, "geneID"] = merged_SNP_dataframe.loc[n, "geneID"] + ";"+ m[1]
    if(len(highmutated_cds_dataframe_list)!=0):
        highmutated_cds_dataframe = pd.concat(highmutated_cds_dataframe_list)
        highmutated_cds_dataframe = highmutated_cds_dataframe.fillna('')
    else:
        highmutated_cds_dataframe = []
    merged_DIP_dataframe = merged_DIP_dataframe[merged_DIP_dataframe.iloc[:, 0] != '']
    merged_DIP_dataframe = merged_DIP_dataframe.fillna('')
    merged_SNP_dataframe = merged_SNP_dataframe[merged_SNP_dataframe.iloc[:, 0] != '']
    merged_SNP_dataframe = merged_SNP_dataframe.fillna('')
    return [merged_DIP_dataframe,merged_SNP_dataframe,highmutated_cds_dataframe,nonmutated_cds_count,mutated_cds_count,highmutated_cds_count,fisrt_filter_geneID,nonmutated_geneID]


def CDSseq_rc(CDSseq_list,orientation):
    if(orientation == "+"):
        return CDSseq_list
    elif(orientation == "-"):
        new_CDSseq_list = CDSseq_list[::-1]
        for i in range(0,len(new_CDSseq_list)):
            if (new_CDSseq_list[i] == 'A'):
                new_CDSseq_list[i] = 'T'
            elif (new_CDSseq_list[i] == 'T'):
                new_CDSseq_list[i] = 'A'
            elif (new_CDSseq_list[i] == 'C'):
                new_CDSseq_list[i] = 'G'
            elif (new_CDSseq_list[i] == 'G'):
                new_CDSseq_list[i] = 'C'
        return new_CDSseq_list


def check_CDSinternal_stopCodon(nocodon_CDSseq_list):
    for i in range(0,len(nocodon_CDSseq_list),3):
        seq = "".join( nocodon_CDSseq_list[i:i+3] )
        if (seq in ["TAA","TGA","TAG"]):
            return True
            break
    else:
        return False


def get_codon_index(SEQ_list,targrt_codon_list):
    codon_index = []
    for i in range(len(SEQ_list)-2):
        seq = "".join( SEQ_list[i:i+3] )
        if(seq in targrt_codon_list):
            codon_index.append(i)
    return codon_index


def check_CDS(CDSseq_list,nonCDSseq_list_up,nonCDSseq_list_down,extendlst):
    extendI5,extendI3 = extendlst[1],extendlst[2]
    start_codon_list = list("".join(CDSseq_list))[0:3]
    start_codon = "".join(start_codon_list)
    nocodon_list = list("".join(CDSseq_list))[3:-3]
    stop_codon_list = list("".join(CDSseq_list))[-3:]
    stop_codon = "".join(stop_codon_list)
    count = (len(CDSseq_list) % 3)
    
    check_bool = False
    if ( (count==0) and (start_codon in ["ATG","GTG","TTG","CTG"]) and (stop_codon in ["TAA","TGA","TAG"]) and ( check_CDSinternal_stopCodon(nocodon_list) == False ) ):
        check_bool = True
        MutType = "normal"
    elif ( (count==0) and (start_codon in ["ATG","GTG","TTG","CTG"]) and (stop_codon in ["TAA","TGA","TAG"]) and ( check_CDSinternal_stopCodon(nocodon_list) == True ) ):
        check_bool = True
        MutType = "abnormal"
    else:
        new_CDS_list = CDSseq_list[extendI5:-extendI3] 
        start_codon_TargrtSeq = nonCDSseq_list_up + CDSseq_list[0:extendI5]
        stop_codon_TargrtSeq = CDSseq_list[-extendI3:] + nonCDSseq_list_down
        start_index = get_codon_index(start_codon_TargrtSeq,["ATG","GTG","TTG","CTG"])
        stop_index = get_codon_index(stop_codon_TargrtSeq,["TAA","TGA","TAG"])
        for m in start_index:
            for n in stop_index:
                if ( (m+3)<len(start_codon_TargrtSeq)  ):
                    new_CDSseq_list_final = start_codon_TargrtSeq[m+3:] + new_CDS_list + stop_codon_TargrtSeq[:-(len(stop_codon_TargrtSeq)-n)]
                else:
                    new_CDSseq_list_final = new_CDS_list + stop_codon_TargrtSeq[:-(len(stop_codon_TargrtSeq)-n)]
                if ( (len(new_CDSseq_list_final) % 3 == 0) and ( check_CDSinternal_stopCodon(new_CDSseq_list_final) == False ) ):
                    check_bool = True
                    MutType = "normal"
                    break
            if (check_bool == True):
                break

    if ( check_bool != True ):
        MutType = "abnormal"
    return MutType


def add_CDS_SNPandDIP(geneID,Target_seq_list, CDS_dataframe_cols_list, cds_start_site, cds_end_site, orientation, TargetRegion_ALLDIP_dataframe, TargetRegion_ALLSNP_dataframe, extendlst):
    result,extendO = [],extendlst[0]
    original_Target_seq_list = copy.deepcopy(Target_seq_list)
    targetseq_list = Target_seq_list
    for i in range(3,len(CDS_dataframe_cols_list)):
        strain = CDS_dataframe_cols_list[i]
        if (len(TargetRegion_ALLSNP_dataframe) != 0):
            for m in TargetRegion_ALLSNP_dataframe.index.tolist():
                loc_m = int(m.split("_!T!E!S!T!_")[1])
                if (TargetRegion_ALLSNP_dataframe.loc[m, strain] != ''):
                    targetseq_list[loc_m-1-(cds_start_site-extendO-1)] = TargetRegion_ALLSNP_dataframe.loc[m, strain]
        if (len(TargetRegion_ALLDIP_dataframe) != 0):
            for n in TargetRegion_ALLDIP_dataframe.index.tolist():
                loc_n = int(n.split("_!T!E!S!T!_")[1])
                if (TargetRegion_ALLDIP_dataframe.loc[n, strain] != ''):
                    DIP_n = list(TargetRegion_ALLDIP_dataframe.loc[n, strain])
                    if ("-" in DIP_n):
                        count = len(DIP_n)
                        if ( (loc_n-1-(cds_start_site-extendO-1) + count) <= len(targetseq_list) ):
                            targetseq_list[loc_n-1-(cds_start_site-extendO-1):(loc_n-1-(cds_start_site-extendO-1) + count)] = [''] * count
                        else:
                            targetseq_list[loc_n-1-(cds_start_site-extendO-1):len(targetseq_list)] = [''] * (count - (loc_n-1-(cds_start_site-extendO-1) + count - len(targetseq_list)))
                    else:
                        targetseq_list[loc_n-1-(cds_start_site-extendO-1)] =  TargetRegion_ALLDIP_dataframe.loc[n, strain] + targetseq_list[loc_n-1-(cds_start_site-extendO-1)]
        
        tmp_up = list("".join(targetseq_list[0:extendO]))
        tmp_down = list("".join(targetseq_list[-extendO:]))
        tmp_cds = list("".join(targetseq_list[extendO:-extendO]))
        if ( orientation == "-" ):
            nonCDS_seq_list_up = CDSseq_rc(tmp_down,orientation)
            nonCDS_seq_list_down = CDSseq_rc(tmp_up,orientation)
        else:
            nonCDS_seq_list_up = CDSseq_rc(tmp_up,orientation)
            nonCDS_seq_list_down = CDSseq_rc(tmp_down,orientation)
        CDS_seq_list = CDSseq_rc(tmp_cds,orientation)
        result.append(check_CDS(CDS_seq_list,nonCDS_seq_list_up,nonCDS_seq_list_down,extendlst))
        targetseq_list = copy.deepcopy(original_Target_seq_list)
    return result


def get_MAF01_dataframe(geneID,CDSkey_dic,CDS_dic,MAF):
    data = []
    key = CDSkey_dic[geneID]
    data.append(key[0])
    data.append(CDS_dic[key][0][0])
    data.append(CDS_dic[key][0][1])
    data.append(CDS_dic[key][0][2])
    data.append(CDS_dic[key][0][3])
    data.append(MAF)
    return data


def CDSgwas(args):
    if args.allDIP != None:
        DIPfiles = glob.glob(args.allDIP+r"/*")
    else:
        raise Exception("Please provide DIP files path!")

    if args.allSNP != None:
        SNPfiles = glob.glob(args.allSNP+r"/*")
    else:
        raise Exception("Please provide SNP files path!")

    if args.CDSdic != None:
        CDSdic = args.CDSdic
    else:
        raise Exception("Please provide CDS annotation dictionaries of the species to be analyzed!")
    
    if args.genome != None:
        genome = args.genome
    else:
        raise Exception("Please provide the FASTA format genome file of the species to be analyzed!")

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
    
    if args.length != None:
        etdO = args.length[0]
        etdI5 = args.length[1]
        etdI3 = args.length[2]
        extdlst = [etdO,etdI5,etdI3]
    else:
        raise Exception("Please provide the length extending outward and inward!")

    if args.cols != None:
        ctigindx_S = args.cols[0]
        locindx_S = args.cols[1]
        aleindx_S = args.cols[2]
        ctigindx_I = args.cols[3]
        locindx_I = args.cols[4]
        aleindx_I = args.cols[5]
    else:
        raise Exception("Please provide the column index of Loc and Allele of the mutation annotation file!")

    if args.thread != None:
        thread = args.thread
    
    if args.threshold != None:
        threshold = args.threshold

    if args.outdir != None:
        outdir = args.outdir
    
    if args.prefix != None:
        prefix = args.prefix

    rootout = outdir + "/" + prefix
    out = rootout + ".highlymutated"
    out0 = rootout + ".nonmutated"
    out1 = rootout + ".allDIP"
    out2 = rootout + ".allSNP"
    out3 = rootout + ".wrongDIPs"
    out4 = rootout + ".MutTypes"
    out5 = rootout + ".MAF0and1"
    out6 = rootout + ".information"
    out7 = rootout + ".LogistRegre.all.CHI2"
    out8 = rootout + ".LogistRegre.all"
    out9 = rootout + ".LogistRegre.filtered"

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


    #0. load the annotation dictionaries creat a anno_CDSkey_dic 
    anno_dic_CDS = load_anno_dic(CDSdic)
    anno_CDSkey_dic,contigs = {},[]
    for i in anno_dic_CDS.keys():
        anno_CDSkey_dic[i[1]] = i
        contigs.append(i[0])
    contigs = list(set(contigs))

    #1. get genome list
    genome_file = open(genome,"r")
    genome_tmp = genome_file.read()
    Contigdic = {}
    for i in contigs:
        tag = i
        re_genome = re.compile(">"+tag+" "+ '.+\n' + '((\w+\n)+)')
        genome_new = list(re_genome.findall(genome_tmp)[0][0])
        while "\n" in genome_new:
            genome_new.remove("\n")
        Contigdic[tag] = genome_new
        logging.info('The length of the contig {0} of the species being analyzed is:\t{1}'.format(tag,len(genome_new)))

    #2. get merged allDIP and allSNP ALT dataframe
    ALLDIP_data_nonfilter = merge_DIPorSNP_ALT(DIPfiles,aleindx_I,locindx_I,ctigindx_I)
    logging.info('All InDels have been merged.')
    ALLSNP_data_nonfilter = merge_DIPorSNP_ALT(SNPfiles,aleindx_S,locindx_S,ctigindx_S)
    logging.info('All SNPs have been merged.')

    #3. get each CDS's SNPs and DIPs of all strains and filter the CDS which has "wrong" DIP
    get_DIPorSNP_geneID_result = get_DIPorSNP_geneID(ALLDIP_data_nonfilter,ALLSNP_data_nonfilter,anno_dic_CDS)
    highmutated_dataframe = get_DIPorSNP_geneID_result[2]
    if(len(highmutated_dataframe)!=0):
        header_str = '\t'.join(highmutated_dataframe.columns.values)
        np.savetxt(out, highmutated_dataframe.values ,delimiter="\t" ,fmt='%s',header= header_str)
    else:
        OUT = open(out,"w")
        OUT.write("There is no highly mutated CDS.")
        OUT.close()

    nonmutated_geneID = get_DIPorSNP_geneID_result[7]
    nonmutate_DATA = {}
    for i in nonmutated_geneID:
        result = get_MAF01_dataframe(i, anno_CDSkey_dic,anno_dic_CDS, "non")
        nonmutate_DATA[i] = pd.Series(result,index=["CDS_RefName","CDS_start","CDS_end","CDS_orientation","Gene","Mutation"])
        nonmutate_data=pd.DataFrame(nonmutate_DATA)
    del nonmutate_DATA
    nonmutate_dataF = pd.DataFrame(nonmutate_data.values.T, index=nonmutate_data.columns, columns=nonmutate_data.index)
    nonmutate_dataF.insert(0,"geneID",nonmutate_dataF.index.tolist())
    header_str = '\t'.join(nonmutate_dataF.columns.values)
    np.savetxt(out0, nonmutate_dataF.values ,delimiter="\t" ,fmt='%s',header= header_str)
    del nonmutate_dataF

    allDIP_data = get_DIPorSNP_geneID_result[0]
    allSNP_data = get_DIPorSNP_geneID_result[1]
    ALLDIP_data_nonfilter = ALLDIP_data_nonfilter.fillna('')
    ALLSNP_data_nonfilter = ALLSNP_data_nonfilter.fillna('')
    header_str = '\t'.join(allDIP_data.columns.values)
    np.savetxt(out1, allDIP_data.values ,delimiter="\t" ,fmt='%s',header= header_str)
    header_str = '\t'.join(allSNP_data.columns.values)
    np.savetxt(out2, allSNP_data.values ,delimiter="\t" ,fmt='%s',header= header_str)

    wrong_index = []
    colnames = allDIP_data.columns
    for m in ALLDIP_data_nonfilter.index.tolist():
        for n in range(3,len(colnames)):
            list_mn = list(set(list(ALLDIP_data_nonfilter.loc[m,colnames[n]])))
            for i in list_mn:
                if(i not in ["A","T","C","G","","-","N"]):
                    wrong_index.append(m)
    right_index = list(set(ALLDIP_data_nonfilter.index.tolist()).difference(set(wrong_index)))
    wrong_allDIP_nonfilter_data = ALLDIP_data_nonfilter.drop(right_index)
    header_str = '\t'.join(wrong_allDIP_nonfilter_data.columns.values)
    np.savetxt(out3, wrong_allDIP_nonfilter_data.values ,delimiter="\t" ,fmt='%s',header= header_str)
    ALLDIP_data_nonfilter = ALLDIP_data_nonfilter.drop(wrong_index)

    wrong_geneID = list(set(wrong_allDIP_nonfilter_data["geneID"].values.tolist()))
    wrong_geneID_list = []
    filter_index = []
    for i in wrong_geneID:
        if(i == ''):
            continue
        filter_data = ALLDIP_data_nonfilter[ALLDIP_data_nonfilter.iloc[:,0] == i]
        if(len(filter_data)!=0):
            for m in filter_data.index.tolist():
                filter_index.append(m)
        split_i = i.split(";")
        if(len(split_i)==1):
            wrong_geneID_list.append(split_i[0])
        elif(len(split_i)>1):
            for n in split_i:
                wrong_geneID_list.append(n)
                filter_data = ALLDIP_data_nonfilter[ALLDIP_data_nonfilter.iloc[:,0] == n]
                if(len(filter_data)!=0):
                    for m in filter_data.index.tolist():
                        filter_index.append(m)
    ALLDIP_data_nonfilter = ALLDIP_data_nonfilter.drop(filter_index)

    allDIPandSNP_geneID = get_DIPorSNP_geneID_result[6]
    wrong_geneID_list = list(set(wrong_geneID_list))
    for i in wrong_geneID_list:
        allDIPandSNP_geneID.remove(i)

    #4. add DIPs and SNPs to each CDS of all strains and judge their MutTypes
    data2 = {}
    for i in allDIPandSNP_geneID: 
            CDS_ctig = anno_CDSkey_dic[i][0]
            CDS_start = anno_dic_CDS[anno_CDSkey_dic[i]][0][0]
            CDS_end = anno_dic_CDS[anno_CDSkey_dic[i]][0][1]
            CDS_ori = anno_dic_CDS[anno_CDSkey_dic[i]][0][2]
            ALLDIP_dataframe_i = ALLDIP_data_nonfilter[ (ALLDIP_data_nonfilter.iloc[:, 1] == CDS_ctig) & ( ALLDIP_data_nonfilter.iloc[:, 2] >= CDS_start-etdO) & (ALLDIP_data_nonfilter.iloc[:, 2] <= CDS_end+etdO) ]
            ALLSNP_dataframe_i = ALLSNP_data_nonfilter[ (ALLSNP_data_nonfilter.iloc[:, 1] == CDS_ctig) & ( ALLSNP_data_nonfilter.iloc[:, 2] >= CDS_start-etdO) & (ALLSNP_data_nonfilter.iloc[:, 2] <= CDS_end+etdO) ]
            data_cols = ALLDIP_data_nonfilter.columns
            genome_new = Contigdic[CDS_ctig]
            if( (CDS_start-etdO-1)<0 ):
                up_part = ["N"]*(etdO+1-CDS_start) + genome_new[0:CDS_start-1]
            else:
                up_part = genome_new[CDS_start-etdO-1:CDS_start-1]
            cds_part = genome_new[CDS_start-1:CDS_end]
            if ( (CDS_end+etdO)>len(genome_new) ):
                down_part = genome_new[CDS_end:] + ["N"]*(CDS_end+etdO-len(genome_new))
            else:
                down_part = genome_new[CDS_end:CDS_end+etdO]
            Regionseq_list = up_part+cds_part+down_part
            Regionseq_cdslist = genome_new[CDS_start-1:CDS_end]
            data_i = add_CDS_SNPandDIP(i,Regionseq_list, data_cols, CDS_start, CDS_end, CDS_ori, ALLDIP_dataframe_i, ALLSNP_dataframe_i, extdlst)
            data2[i] = pd.Series(data_i,index=data_cols[3:])

    dataF = pd.DataFrame(data2)
    dataF2 = pd.DataFrame(dataF.values.T, index=dataF.columns, columns=dataF.index)
    dataF2.insert(0,'geneID',dataF2.index)
    del data2
    del dataF
    del allSNP_data
    del allDIP_data
    if args.ref != None:
        dataF2.insert(dataF2.shape[1],str(args.ref[0]),["normal"]*len(dataF2))
        re_pseudo = re.compile(r'pseudo')
        for m in dataF2.index.tolist():
            if (len( re_pseudo.findall(m)) != 0):
                dataF2.loc[m, str(args.ref[0])] = "abnormal"
    
    header_str = '\t'.join(dataF2.columns.values)
    np.savetxt(out4, dataF2.values ,delimiter="\t" ,fmt='%s',header= header_str)
    Muttps = copy.copy(dataF2)

    #5. calculate the MAF of each CDS
    dataF = dataF2.iloc[:,1:]
    strains_count = dataF.shape[1]
    allCount=dataF.apply(pd.value_counts,axis=1)
    allCount = pd.DataFrame(allCount)
    del dataF
    NaN_count = allCount.isnull().sum(axis=1)
    allele_filter1 = NaN_count[NaN_count!=0].index.tolist()   
    allCount = allCount.drop(allele_filter1)
    dataF2 = dataF2.drop(allele_filter1)
    allCount_min = allCount.min(axis=1)
    minor_allele = pd.DataFrame(allCount.idxmin(axis=1))
    major_allele = pd.DataFrame(allCount.idxmax(axis=1))
    allele_filter2 = allCount_min[allCount_min < 2].index.tolist()  
    dataF2 = dataF2.drop(allele_filter2)
    minor_allele = minor_allele.drop(allele_filter2)
    major_allele = major_allele.drop(allele_filter2)

    cols = dataF2.columns
    dataF2.columns = range(dataF2.shape[1])
    dataF2.insert(dataF2.shape[1],'Minor',minor_allele.iloc[:,0].tolist())
    dataF2.insert(dataF2.shape[1],'Major',major_allele.iloc[:,0].tolist())
    for i in range(1,dataF2.shape[1]-2):
        dataF2[i] = np.where((dataF2["Minor"] == dataF2[i]), 2, dataF2[i])  
        dataF2[i] = np.where((dataF2["Major"] == dataF2[i]), 1, dataF2[i])  
    dataF2 = dataF2.iloc[:,0:dataF2.shape[1]-2]
    dataF2.columns = cols
    dataF2 = dataF2.replace('normal',1)

    DATA = {}
    for i in allele_filter1:
        result = get_MAF01_dataframe(i, anno_CDSkey_dic,anno_dic_CDS,0)
        DATA[i] = pd.Series(result,index=["CDS_RefName","CDS_start","CDS_end","CDS_orientation","gene","MAF"])
    for i in allele_filter2:
        result = get_MAF01_dataframe(i, anno_CDSkey_dic,anno_dic_CDS,1)
        DATA[i] = pd.Series(result,index=["CDS_RefName","CDS_start","CDS_end","CDS_orientation","gene","MAF"])
    data=pd.DataFrame(DATA)
    del DATA
    dataF = pd.DataFrame(data.values.T, index=data.columns, columns=data.index)
    dataF.insert(0,"geneID",dataF.index.tolist())
    header_str = '\t'.join(dataF.columns.values)
    np.savetxt(out5, dataF.values ,delimiter="\t" ,fmt='%s',header= header_str)
    del dataF

    #6. get the information table for LogisticRegre
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
    nrow = len(dataF2)
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
                tab_i[m] = pd.Series([cols_i[m]]*nrow, index=dataF2.index)
            if(GCA_i in GCAs.keys()):
                tab_i["pheno"] = pd.Series([GCAs[GCA_i]] * nrow, index=dataF2.index)
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
            tab_i["pheno"] = pd.Series([GCAs[m]] * nrow, index=dataF2.index)
            data_i = pd.DataFrame(tab_i)
            data_i.columns = ["Phenotype"]
            tab[m] = data_i
        Rcols = 'Allele,Phenotype'
    final = []
    for i in tab.keys():
        tab[i].insert(0, 'Allele', dataF2[i].values.tolist())
        final.append(tab[i])
    finaltab = pd.concat(final, axis=1)
    finaltab.insert(0, "geneID", dataF2["geneID"].values.tolist())
    header_str = '\t'.join(finaltab.columns.values)
    np.savetxt(out6, finaltab.values, delimiter="\t", fmt='%s', header=header_str)


    #7. get the chi2 distribution
    oldname = Muttps.columns.tolist() 
    newname = ["geneID"]
    for i in range(1,len(oldname)):
        GCA_i = oldname[i]
        if(GCA_i in PH1_GCAs):
            newname.append(phe1)
        elif(GCA_i in PH0_GCAs):
            newname.append(phe0)
    Muttps.index = Muttps.iloc[:,0].tolist()
    Muttps.columns = newname    
    Muttps_phe1 = Muttps[phe1]
    Muttps_phe0 = Muttps[phe0]

    phe1allCount = Muttps_phe1.apply(pd.value_counts,axis=1)
    phe1allCount = pd.DataFrame(phe1allCount)
    phe1_cols=[]
    for i in range(0,phe1allCount.shape[1]):
        phe1_cols.append(phe1+"."+str(phe1allCount.columns.tolist()[i]))
    phe1allCount.columns = phe1_cols
    phe1allCount = phe1allCount.fillna("0")
    phe1allCount = pd.DataFrame(phe1allCount,dtype="int")

    phe0allCount = Muttps_phe0.apply(pd.value_counts,axis=1)
    phe0allCount = pd.DataFrame(phe0allCount)
    phe0_cols=[]
    for i in range(0,phe0allCount.shape[1]):
        phe0_cols.append(phe0+"."+str(phe0allCount.columns.tolist()[i]))
    phe0allCount.columns = phe0_cols
    phe0allCount = phe0allCount.fillna("0")
    phe0allCount = pd.DataFrame(phe0allCount,dtype="int")

    AllCount = pd.concat([phe1allCount,phe0allCount],axis=1)
    AllCount = pd.DataFrame(AllCount)
    geneIDs = AllCount.index.tolist()
    AllCount.insert(0,"geneID",geneIDs)
    colnames = AllCount.columns.tolist()
    fil = list(set(geneIDs).difference(set(finaltab["geneID"].values.tolist())))
    AllCount = AllCount.drop(fil)
    header_str = '\t'.join(AllCount.columns.values)
    np.savetxt(out8, AllCount.values ,delimiter="\t" ,fmt='%s',header= header_str)
    input1.close()
    input2.close()


    os.system('%s %s %s %s %s %s %s CDS %d %s' % (Rscr,script1,out8,colnames[1],colnames[2],colnames[3],colnames[4],thread,out7))

    os.system('%s %s %s %s %s CDS %d %f %s %s' % (Rscr,script2,out6,Rcols,types,thread,threshold,out8,out9))


