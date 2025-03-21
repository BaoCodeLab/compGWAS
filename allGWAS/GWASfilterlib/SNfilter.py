import pickle
import pandas as pd


class StrToBytes:
    def __init__(self,fileobj):
        self.fileobj = fileobj
    def read(self,size):
        return self.fileobj.read(size).encode()
    def readline(self,size=-1):
        return self.fileobj.readline(size).encode()


def check_dist(mat,gendic,threshold5,threshold3):
    out = []
    mat5 = mat[mat["MutType"]=="5"]
    mat3 = mat[mat["MutType"]=="3"]
    mat5.index = list(range(len(mat5)))
    mat3.index = list(range(len(mat3)))
    if len(mat5) == 0:
        out.append(pd.DataFrame())
    else:
        filt5 = []
        for ii in range(len(mat5)):
            dist = abs(mat5.iloc[ii,0] - gendic[mat5.iloc[ii,5]])
            if dist > threshold5:
                filt5.append(ii)
        mat5 = mat5.drop(filt5)
        if len(mat5) == 0:
            out.append(pd.DataFrame())
        else:
            out.append(mat5)
    if len(mat3) == 0:
        out.append(pd.DataFrame())
    else:       
        filt3 = []
        for ii in range(len(mat3)):
            dist = abs(mat3.iloc[ii,0] - gendic[mat3.iloc[ii,5]])
            if dist > threshold3:
                filt3.append(ii)
        mat3 = mat3.drop(filt3)
        if len(mat3) == 0:
            out.append(pd.DataFrame())
        else:
            out.append(mat3)
    return out


def getres(mat):
    Loco,Estimateo,Pvalueo,MutTypeo,Geneo,GeneIDo = [],[],[],[],[],[]
    IDs = list(set(mat["GeneID"].tolist()))
    for ii in IDs:
        mat2 = mat[mat["GeneID"]==ii]
        mat3 = mat2[mat2["Pvalue"]==mat2["Pvalue"].min()]
        Loco.append(mat3.iloc[0,0])
        Estimateo.append(mat3.iloc[0,1])
        Pvalueo.append(mat3.iloc[0,2])
        MutTypeo.append(mat3.iloc[0,3])
        Geneo.append(mat3.iloc[0,4])
        GeneIDo.append(mat3.iloc[0,5])
    return [Loco,Estimateo,Pvalueo,MutTypeo,Geneo,GeneIDo]


def SNfilter(args):
    if args.CDSdic != None:
        CDSdic = args.CDSdic
    else:
        raise Exception("Please provide CDS annotation dictionaries of the species to be analyzed!")

    if args.LogisAnno != None:
        LogisAnno = args.LogisAnno
        Out = LogisAnno + ".SN.filterout"
    else:
        raise Exception("Please provide the annotation file of SNP logistic regression results!")

    if args.anno != None:
        anno = args.anno
    else:
        raise Exception("Please provide the annotation file of all genes of the species to be analyzed!")

    if args.Chi2 != None:
        Chi2 = args.Chi2
    else:
        raise Exception("Please provide the chi2 distribution file of the 2 phenotypes!")

    if args.threshold != None:
        threshold = args.threshold
                
    if args.distance != None:
        distance5 = args.distance[0]
        distance3 = args.distance[1]
    else:
        raise Exception("Please provide the threshold of 5' distance(bp) and 3' distance(bp)!")


    dat = pd.read_csv(anno,"\t")
    dat = dat.fillna('')
    geneid = dat.iloc[:,10].tolist()
    sub1 = dat.iloc[:,16].tolist()
    subanno = {}
    for i in range(len(geneid)):
        if(sub1[i] == ''):
            subanno[geneid[i]] = 'Unknown'
        else:
            subanno[geneid[i]] = sub1[i].split(" //")[0].rstrip(' ')

    anno_dic = {}
    for m in range(0,len(CDSdic)):
        anno_m = open(CDSdic[m],"r")
        anno_m_dic = pickle.load(StrToBytes(anno_m))
        anno_dic.update(anno_m_dic)
        anno_m.close()

    CHI2d = pd.read_csv(Chi2,sep="\t")
    DATA = pd.read_csv(LogisAnno,"\t")
    contigs = list(set(DATA.iloc[:,0].tolist()))
    finalOUT = []
    for jj in contigs:
        chi2flt = []
        chi2d = CHI2d[CHI2d.iloc[:,0]==jj]
        Locs = chi2d.iloc[:,1]
        freq1 = chi2d.iloc[:,2]
        freq2 = chi2d.iloc[:,3]
        freq3 = chi2d.iloc[:,4]
        freq4 = chi2d.iloc[:,5]
        for i in range(len(Locs)):
            freqs = sorted([freq1[i],freq2[i],freq3[i],freq4[i]])
            if freqs[1] < 3:
                chi2flt.append(Locs[i])

        newdic = {}
        for i in anno_dic.keys():
            if i[0] == jj:
                start_i = int(anno_dic[i][0][0])
                newdic[start_i] = i
        keys = sorted(list(newdic.keys()))

        # filter the start site of gene whose both start site and end site are in an another gene
        # creat a dic (anno_dic_2) which contains all CDS regions
        genleft,genright,filter_start = {},{},[]
        for mm in keys:
            i = newdic[mm]
            start_i = int(anno_dic[i][0][0])
            end_i = int(anno_dic[i][0][1])
            geneID_left_i = i[1]
            geneID_right_i = i[1]
            for n in keys:
                m = newdic[n]
                start_m = int(anno_dic[m][0][0])
                end_m = int(anno_dic[m][0][1])
                geneID_left_m = m[1]
                geneID_right_m = m[1]
                if (start_i < start_m and end_i > end_m):
                    filter_start.append(int(start_m))
                if (start_i < start_m and end_i < end_m and end_i >= start_m):
                    filter_start.append(int(start_m))
                    end_i = end_m
                    geneID_right_i = geneID_right_m
                if (start_i > start_m and start_i <= end_m and end_i >= end_m):
                    filter_start.append(int(start_i))
                    start_i = start_m
                    geneID_left_i = geneID_left_m
                if (start_i > start_m and end_i < end_m):
                    filter_start.append(int(start_i))
                    start_i = start_m
                    geneID_left_i = geneID_left_m
                    end_i = end_m
                    geneID_right_i = geneID_right_m
            if (start_i not in filter_start):
                genleft[geneID_left_i] = int(start_i)
                genright[geneID_right_i] = int(end_i)

        filtindx = []
        data = DATA[DATA.iloc[:,0]==jj]
        data = data[(data["MutType"].str.contains("inter-gene")) & (data["Allele2_Pvalue"] <= threshold)]
        data.index = list(range(len(data)))
        for i in range(len(data)):
            if data.iloc[i,1] in chi2flt:
                filtindx.append(i)
                continue

        data = data.drop(filtindx)
        data = data.iloc[:,1:]
        genIDs = list(set(data["GeneID"].tolist()))

        Loc_o,Estimate_o,Pvalue_o,MutType_o,Gene_o,GeneID_o = [],[],[],[],[],[]
        for i in genIDs:
            datmp = data[data["GeneID"]==i]
            tmp_l,Loc_l,Estimate_l,Pvalue_l,MutType_l,Gene_l,GeneID_l = {},[],[],[],[],[],[]
            tmp_r,Loc_r,Estimate_r,Pvalue_r,MutType_r,Gene_r,GeneID_r = {},[],[],[],[],[],[]
            for m in range(len(datmp)):
                Loc_l.append(datmp.iloc[m,0])
                Loc_r.append(datmp.iloc[m,0])
                Estimate_l.append(datmp.iloc[m,1])
                Estimate_r.append(datmp.iloc[m,1])
                Pvalue_l.append(datmp.iloc[m,2])
                Pvalue_r.append(datmp.iloc[m,2])
                MutType_l.append(datmp.iloc[m,3].split("inter-gene-")[-1].split("-")[0])
                MutType_r.append(datmp.iloc[m,3].split("inter-gene-")[-1].split("-")[1])
                Gene_l.append(datmp.iloc[m,4].split("||")[0])
                Gene_r.append(datmp.iloc[m,4].split("||")[1])
                GeneID_l.append(datmp.iloc[m,5].split("||")[0])
                GeneID_r.append(datmp.iloc[m,5].split("||")[1])
            tmp_l["Loc"] = pd.Series(Loc_l, index=list(range(len(Loc_l))))
            tmp_r["Loc"] = pd.Series(Loc_r, index=list(range(len(Loc_r))))
            tmp_l["Estimate"] = pd.Series(Estimate_l, index=list(range(len(Estimate_l))))
            tmp_r["Estimate"] = pd.Series(Estimate_r, index=list(range(len(Estimate_r))))
            tmp_l["Pvalue"] = pd.Series(Pvalue_l, index=list(range(len(Pvalue_l))))
            tmp_r["Pvalue"] = pd.Series(Pvalue_r, index=list(range(len(Pvalue_r))))
            tmp_l["MutType"] = pd.Series(MutType_l, index=list(range(len(MutType_l))))
            tmp_r["MutType"] = pd.Series(MutType_r, index=list(range(len(MutType_r))))
            tmp_l["Gene"] = pd.Series(Gene_l, index=list(range(len(Gene_l))))
            tmp_r["Gene"] = pd.Series(Gene_r, index=list(range(len(Gene_r))))
            tmp_l["GeneID"] = pd.Series(GeneID_l, index=list(range(len(GeneID_l))))
            tmp_r["GeneID"] = pd.Series(GeneID_r, index=list(range(len(GeneID_r))))
            tmpF_l = pd.DataFrame(tmp_l)
            tmpF_r = pd.DataFrame(tmp_r)
            tmpF5_l = check_dist(tmpF_l,genright,distance5,distance3)[0]
            tmpF3_l = check_dist(tmpF_l,genright,distance5,distance3)[1]
            tmpF5_r = check_dist(tmpF_r,genleft,distance5,distance3)[0]
            tmpF3_r = check_dist(tmpF_r,genleft,distance5,distance3)[1]
            tmpF5 = pd.concat( [tmpF5_l,tmpF5_r], axis=0)
            tmpF3 = pd.concat( [tmpF3_l,tmpF3_r], axis=0)
            if len(tmpF5) != 0:
                tmpF5_o = getres(tmpF5)
                for n in range(len(tmpF5_o[0])):
                    Loc_o.append(tmpF5_o[0][n])
                    Estimate_o.append(tmpF5_o[1][n])
                    Pvalue_o.append(tmpF5_o[2][n])
                    MutType_o.append(tmpF5_o[3][n])
                    Gene_o.append(tmpF5_o[4][n])
                    GeneID_o.append(tmpF5_o[5][n])
            if len(tmpF3) != 0:
                tmpF3_o = getres(tmpF3)
                for n in range(len(tmpF3_o[0])):
                    Loc_o.append(tmpF3_o[0][n])
                    Estimate_o.append(tmpF3_o[1][n])
                    Pvalue_o.append(tmpF3_o[2][n])
                    MutType_o.append(tmpF3_o[3][n])
                    Gene_o.append(tmpF3_o[4][n])
                    GeneID_o.append(tmpF3_o[5][n])

        Subsystem1_o = []
        for i in range(len(GeneID_o)):
            Subsystem1_o.append(subanno[GeneID_o[i]])

        OUT = {}
        OUT["RefName"] = pd.Series([jj]*len(Loc_o), index=list(range(len(Loc_o))))
        OUT["Site"] = pd.Series(Loc_o, index=list(range(len(Loc_o))))
        OUT["Allele2_Coefficient"] = pd.Series(Estimate_o, index=list(range(len(Estimate_o))))
        OUT["Allele2_Pvalue"] = pd.Series(Pvalue_o, index=list(range(len(Pvalue_o))))
        OUT["Type"] = pd.Series(MutType_o, index=list(range(len(MutType_o))))
        OUT["Gene"] = pd.Series(Gene_o, index=list(range(len(Gene_o))))
        OUT["GeneID"] = pd.Series(GeneID_o, index=list(range(len(GeneID_o))))
        OUT["Subsystem1"] = pd.Series(Subsystem1_o, index=list(range(len(Subsystem1_o))))
        OUTF = pd.DataFrame(OUT)
        finalOUT.append(OUTF)
    OUTD = pd.concat(finalOUT)
    OUTD.to_csv(Out, sep="\t", header=True, index=False)



