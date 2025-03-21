import pandas as pd
import numpy as np
import glob
import re
import pickle
import os


def LDprun(args):
    if args.allSNP != None:
        SNPfiles = glob.glob(args.allSNP+r"/*")
    else:
        raise Exception("Please provide SNP files path!")

    if args.pheno1 != None:
        pheno1 = args.pheno1
    else:
        raise Exception("Please provide phenotype1's GCA number file!")

    if args.pheno0 != None:
        pheno0 = args.pheno0
    else:
        raise Exception("Please provide phenotype0's GCA number file!")

    if args.genome != None:
        genome = args.genome
    else:
        raise Exception("Please provide the FASTA format genome file of the species to be analyzed!")
    
    if args.logisre != None:
        logisre = args.logisre
    else:
        raise Exception("Please provide the results file of all SNPs logistic regression!")

    if args.threshold != None:
        threshold = args.threshold

    if args.outdir != None:
        outdir = args.outdir

    if args.prefix != None:
        prefix = args.prefix
    
    rootout = outdir + "/" + prefix

    if args.java != None:
        java = args.java
    else:
        raise Exception("Please provide the path of java!")

    if args.Haploview != None:
        Haploview = args.Haploview
    else:
        raise Exception("Please provide the path of Haploview.jar!")

    if args.Argument != None:
        Argument = args.Argument
    else:
        raise Exception("Please provide the parameters of Haploview.jar, excluding input and output parameters!")

    
    genome_file = open(genome,"r")
    genome_tmp = genome_file.read()
    re_ctig = re.compile(">" + '(\S+)' + " ")
    contigs = re_ctig.findall(genome_tmp)

    input1 = open(pheno1,"r")
    input2 = open(pheno0,"r")
    GCAs_pheno1 = input1.readlines()
    GCAs_pheno0 = input2.readlines()
    GCAs = {}
    for i in GCAs_pheno1:
        GCAs[i.split("\n")[0]] = 1
        for i in GCAs_pheno0:
            GCAs[i.split("\n")[0]] = 2
            if args.ref != None:
                if args.ref[1] == "p":
                    GCAs[str(args.ref[0])] = 1
                elif args.ref[1] == "P":
                    GCAs[str(args.ref[0])] = 2
                else:
                    raise Exception("Please provide right phenotype of reference strain, p(means phenotype1) or P(means phenotype0)!")
    input1.close()
    input2.close()

    data_alt,data_ref = {},{}
    for i in range(0,len(SNPfiles)):
        i_file = SNPfiles[i]
        file_nm = i_file.split("/")[-1].split("GCA-")[-1].split(".")[0]
        file_i = pd.read_csv(i_file, "\t")
        INDEX = list(map(lambda x,y:str(x)+"_!T!E!S!T!_"+str(y), list(file_i.iloc[:,0]), list(file_i.iloc[:,1])))
        data_alt[file_nm] = pd.Series(list(file_i.iloc[:,6]), index=INDEX)
        data_alt[file_nm] = data_alt[file_nm].astype('category')
        data_alt[file_nm] = data_alt[file_nm].cat.add_categories([''])
        data_ref[file_nm] = pd.Series(list(file_i.iloc[:,4]), index=INDEX)
        data_ref[file_nm] = data_ref[file_nm].astype('category')
        data_ref[file_nm] = data_ref[file_nm].cat.add_categories([''])
    ALTData = pd.DataFrame(data_alt)
    REFData = pd.DataFrame(data_ref)
    INDEX = ALTData.index.tolist()
    ctig = list(map(lambda x:x.split("_!T!E!S!T!_")[0], INDEX)) 
    LOC = list(map(lambda x:int(x.split("_!T!E!S!T!_")[1]), INDEX))
    ALTData.insert(0, "Loc", LOC)
    ALTData.insert(0, "RefName", ctig)
    REFData.insert(0, "Loc", LOC)
    REFData.insert(0, "RefName", ctig)

    data_i = pd.read_csv(logisre, "\t")
    account = 0
    for jj in contigs:
        data_i = data_i[(data_i.iloc[:,6]<=float(threshold)) & (data_i.iloc[:,0]==jj)]
        if len(data_i) == 0:
            account = account + 1
            continue
        else:
            out = rootout + "." + jj + ".Haploview"
            out1 = rootout + "." + jj + ".Haploview.info"
            out2 = rootout + "." + jj + ".Haploview.data"
            goalLoc = data_i.iloc[:,1].values.tolist()
            
            AltData = ALTData[ALTData.iloc[:,0]==jj]
            AltData = AltData.sort_values(by=["Loc"],ascending=True)
            AltData.index = AltData.iloc[:,1].tolist()
            AltData = AltData.iloc[:,2:]
            RefData = REFData[REFData.iloc[:,0]==jj]
            RefData = RefData.sort_values(by=["Loc"],ascending=True)
            RefData.index = RefData.iloc[:,1].tolist()
            RefData = RefData.iloc[:,2:]

            filtloc = list(set(AltData.index.tolist()).difference(set(goalLoc)))
            AltData = AltData.drop(filtloc)
            RefData = RefData.drop(filtloc) 

            info_data = {}
            info_data["loc"] = pd.Series(AltData.index.tolist(), index=range(1,len(AltData)+1))
            info_out = pd.DataFrame(info_data)
            info_out.insert(0,'index',info_out.index)
            info_out.to_csv(out1,sep="\t",header=False,index=False)

            if args.ref != None:
                AltData[str(args.ref[0])]= np.nan

            RefData = RefData.fillna("")
            for m in range(0,len(AltData)):
                for n in range(0,RefData.shape[1]):
                    if RefData.iloc[m,n] != "":
                        ref_m = RefData.iloc[m,n]
                        break
                AltData.iloc[m,:] = AltData.iloc[m,:].fillna(ref_m)

            BaseData = pd.DataFrame(AltData.values.T, index=AltData.columns, columns=AltData.index)
            for i in range(0,BaseData.shape[1]):
                BaseData.iloc[:,i] = BaseData.iloc[:,i].str.cat(BaseData.iloc[:,i],sep=" ")
            phenos = []
            for i in BaseData.index.tolist():
                phenos.append(GCAs[i])
            BaseData.insert(0,'Pheno',phenos)
            BaseData.insert(0,'Sex',[1]*len(BaseData))
            BaseData.insert(0,'MoID',[0]*len(BaseData))
            BaseData.insert(0,'FaID',[0]*len(BaseData))
            BaseData.insert(0,'InID',BaseData.index.tolist())
            BaseData.insert(0,'Name',range(1,len(BaseData)+1))
            BaseData.to_csv(out2,sep="\t",header=False,index=False)

            os.system('%s -jar %s -out %s -pedfile %s -info %s %s' % (java,Haploview,out,out2,out1,Argument))

    if account == len(contigs):
        raise Exception("There is no snp (P <= 0.05) in the logistic regression result file!")

