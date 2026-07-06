### version3
import pickle
import logging
import vcf
import re
import os
import glob

contg_dic,gene_dic,CDS_dic,mol_dic,CDS_seq_dic={},{},{},{},{}


class StrToBytes:
    def __init__(self,fileobj):
        self.fileobj = fileobj
    def read(self,size):
        return self.fileobj.read(size).encode()
    def readline(self,size=-1):
        return self.fileobj.readline(size).encode()



def create_contg_dic(contg_filefile):
    global contg_dic
    contg_file=open(contg_filefile,'r')
    contg_dic=pickle.load(StrToBytes(contg_file))

    contg_file.close()
#    return contg_dic


def create_gene_dic(gene_filefile):
    global gene_dic
    gene_file=open(gene_filefile,'r')
    gene_dic=pickle.load(StrToBytes(gene_file))
    gene_file.close()
#    return gene_dic


def create_CDS_dic(CDS_filefile):
    global CDS_dic
    CDS_file=open(CDS_filefile,'r')
    CDS_dic=pickle.load(StrToBytes(CDS_file))
    CDS_file.close()
#    return CDS_dic


def create_mol_dic(mol_filefile):
    global mol_dic
    mol_file=open(mol_filefile,'r')
    mol_dic=pickle.load(StrToBytes(mol_file))
    mol_file.close()
#    return mol_dic


## input CDS sequences for determination of CDS location
def create_CDSseq_dic(CDS_seq_filefile):
    global CDS_seq_dic
    CDS_seq_file=open(CDS_seq_filefile,'r')
    CDS_seq_dic=pickle.load(StrToBytes(CDS_seq_file))
    CDS_seq_file.close()
#    return CDS_seq_dic


cod=    {'CDS_nonSynon':-2,'CDS_synon':-1,'CDS_unanno':0,'UTR':1,'intron-splice-5':2,'intron-splice-3':3,'intron':4,'inter-gene-3-5':5,'inter-gene-5':6,'inter-gene-3':7,'inter-gene':8,'pseudo-gene|ncRNA':9,'inter-gene-5-3':10,'inter-gene-5-5':11,'inter-gene-3-3':12,'inter-gene--5':13,'inter-gene--3':14,'inter-gene-5-':15,'inter-gene-3-':16}

antiCod={-2:'CDS_nonSynon',-1:'CDS_synon',0:'CDS_unanno',1:'UTR',2:'intron-splice-5',3:'intron-splice-3',4:'intron',5:'inter-gene-3-5',6:'inter-gene-5',7:'inter-gene-3',8:'inter-gene',9:'pseudo-gene|ncRNA',10:'inter-gene-5-3',11:'inter-gene-5-5',12:'inter-gene-3-3',13:'inter-gene--5',14:'inter-gene--3',15:'inter-gene-5-',16:'inter-gene-3-'}

nt=['TTT','TTC','TTA','TTG','TCT','TCC','TCA','TCG','TAT','TAC','TAA','TAG','TGT','TGC','TGA','TGG',\
    'CTT','CTC','CTA','CTG','CCT','CCC','CCA','CCG','CAT','CAC','CAA','CAG','CGT','CGC','CGA','CGG',\
    'ATT','ATC','ATA','ATG','ACT','ACC','ACA','ACG','AAT','AAC','AAA','AAG','AGT','AGC','AGA','AGG',\
    'GTT','GTC','GTA','GTG','GCT','GCC','GCA','GCG','GAT','GAC','GAA','GAG','GGT','GGC','GGA','GGG']

aa=['Phe','Phe','Leu','Leu','Ser','Ser','Ser','Ser','Tyr','Tyr','Stp','Stp','Cys','Cys','Stp','Trp',\
    'Leu','Leu','Leu','Leu','Pro','Pro','Pro','Pro','His','His','Gln','Gln','Arg','Arg','Arg','Arg',\
    'Ile','Ile','Ile','Met','Thr','Thr','Thr','Thr','Asn','Asn','Lys','Lys','Ser','Ser','Arg','Arg',\
    'Val','Val','Val','Val','Ala','Ala','Ala','Ala','Asp','Asp','Glu','Glu','Gly','Gly','Gly','Gly']

aa_single=['F','F','L','L','S','S','S','S','Y','Y','*','*','C','C','*','W',\
    'L','L','L','L','P','P','P','P','H','H','Q','Q','R','R','R','R',\
    'I','I','I','M','T','T','T','T','N','N','K','K','S','S','R','R',\
    'V','V','V','V','A','A','A','A','D','D','E','E','G','G','G','G']


coden={}
def coden_dic():
    global coden
    for ii in range(0,len(aa)):
        coden[nt[ii]]=aa[ii]

def char_cmp(str0):
    if str0=='A':
        return 'T'
    elif str0=='T':
        return 'A'
    elif str0=='G':
        return 'C'
    elif str0=='C':
        return 'G'
    else:
        return 'N'

def find_cod(str1,str2,str3,snp_pos,snp_allel): ## str1=contig; str2=geneID; str3=geneAcc
    snpCDSpos=0
    CDS_ky=(str1,str2,str3)
    CDS_seq=CDS_seq_dic_dic[CDS_ky]
    CDS_lg=len(CDS_seq)
    if CDS_lg%3!=0:
        snp_anno,snp_prot,aa_pos='CDS_unanno','',''

        return snp_anno,snp_prot,aa_pos

    CDS_dir=CDS_dic[CDS_ky][0][2]
    for CDS_val in sorted(CDS_dic[CDS_ky]):
        if snp_pos < CDS_val[0]:  ## if snp_pos is located before an exon 
            break
        snpCDSpos+=min(snp_pos,CDS_val[1])-CDS_val[0]+1 # the variation position in the CDS
    trip_st,trip_pos=(snpCDSpos-1)//3*3,snpCDSpos%3 # the coden starting position and the variation position in the coden
    aa_pos=(snpCDSpos-1)//3+1

    trip_ref=CDS_seq[trip_st:trip_st+3]
#   print str2,str3,CDS_lg,snpCDSpos,trip_st,trip_ref
    trip_ref_lst=list(trip_ref)
    trip_allel_lst=list(trip_ref)
    trip_allel_lst[trip_pos-1]=snp_allel

    if CDS_dir=='-':
        aa_pos=CDS_lg/3-(snpCDSpos-1)//3

        trip_ref_lst.reverse()
        trip_allel_lst.reverse()

        trip_ref_lst=map(char_cmp,trip_ref_lst)
        trip_allel_lst=map(char_cmp,trip_allel_lst)

    ref_cod=coden[''.join(trip_ref_lst)]
    allel_cod=coden[''.join(trip_allel_lst)]
#   print CDS_dir,''.join(trip_ref_lst),''.join(trip_allel_lst)

    if ref_cod == allel_cod:
        snp_anno,snp_prot='CDS_synon',''
    else:
        snp_anno,snp_prot='CDS_nonSynon','%s%d%s'  % (ref_cod,aa_pos,allel_cod)

    return snp_anno,snp_prot,int(aa_pos)



def find_cod_multi(str1,str2,str3,snp_pos_lst,snp_allel_lst):

    CDS_ky=(str1,str2,str3)
    CDS_seq=CDS_seq_dic[CDS_ky]
    CDS_lg=len(CDS_seq)
    if CDS_lg%3!=0:
        snp_anno,snp_prot,aa_pos='CDS_unanno','',''
        return snp_anno,snp_prot,aa_pos

    snpCDSpos_lst=[]
    CDS_dir=CDS_dic[CDS_ky][0][2]
		
    for ech_snp_pos in snp_pos_lst:
        snpCDSpos=0
        for CDS_val in sorted(CDS_dic[CDS_ky]):
            if ech_snp_pos < CDS_val[0]:
                break
            snpCDSpos+=min(ech_snp_pos,CDS_val[1])-CDS_val[0]+1
        snpCDSpos_lst+=[snpCDSpos]  ## the position of snp in the CDS

#   trip_st,trip_pos=(snpCDSpos-1)//3*3,snpCDSpos%3 # the coden starting position and the snp frame in the coden

    trip_st=(snpCDSpos_lst[0]-1)//3*3  ## the coden starting position in the CDS, i.e., 6,7,8
    aa_pos=(snpCDSpos_lst[0]-1)//3+1   ## the amino acid position in the protein
#    print(trip_st)

    trip_ref=CDS_seq[trip_st:trip_st+3]
#   print str2,str3,CDS_lg,snpCDSpos,trip_st,trip_ref
    trip_ref_lst=list(trip_ref)
    trip_allel_lst=list(trip_ref)
	
## calculate the snp frame in the codon, and replace the corresponding allele in the reference with the snp
    snpLg=len(snp_pos_lst)
    for tt in range(0,snpLg):
        ech_trip_pos=snpCDSpos_lst[tt]%3  ## the snp frame in the codon, i.e.,1,2,3
        trip_allel_lst[ech_trip_pos-1]=snp_allel_lst[tt] ## replace the allele in the reference with the snp

    if CDS_dir=='-':
        aa_pos=CDS_lg/3-(snpCDSpos_lst[0]-1)//3

        trip_ref_lst.reverse()
        trip_allel_lst.reverse()

        trip_ref_lst=map(char_cmp,trip_ref_lst)
        trip_allel_lst=map(char_cmp,trip_allel_lst)

#    print(trip_allel_lst)
    ref_cod=coden[''.join(trip_ref_lst)]
    allel_cod=coden[''.join(trip_allel_lst)]
    
    if ref_cod == allel_cod:
        snp_anno,snp_prot='CDS_synon',''
    else:
        snp_anno,snp_prot='CDS_nonSynon','%s%d%s'  % (ref_cod,aa_pos,allel_cod)

    return snp_anno,snp_prot,int(aa_pos)



def format_var_file(snp_input,fileType):

    mat=[]
    snp_file=open(snp_input,'r')

    if fileType == "table":
        lines = snp_file.readlines()
        for lin in lines:
            if lin[0] == '#':   # skip the header line
                continue

            lst=lin.split('\t')
            contig,pos=lst[0].split('.')[0],lst[1]
            allel_ref,allel=lst[4],lst[6]
            mutType,refLen,allelNum=lst[2],lst[3],lst[5]
            mat+=[[contig,pos,mutType,refLen,allel_ref,allelNum,allel]]
#            mat+=[[contig,pos,pos,mutType,refLen,allel_ref,allelNum,allel]]


    elif fileType == "vcf":
    
        vcf_reader=vcf.Reader(snp_file)
        for lin in vcf_reader:
#            print(lin.POS,lin.start,lin.end,lin.affected_start,lin.affected_end,lin.is_indel,lin.alleles)            
            contig,pos=str(lin.CHROM),str(lin.POS)
            allel_ref,allel=str(lin.REF),list(map(str,lin.ALT))
            refLen,allelNum=str(len(allel_ref)),str(len(allel))
            mutType = "SNP" if lin.is_indel==False else "DIP" 
            mat+=[[contig,pos,mutType,refLen,allel_ref,allelNum,'|'.join(allel)]]
#            mat+=[[contig,pos,pos,mutType,refLen,allel_ref,allelNum,'|'.join(allel)]]

    snp_file.close()

    return mat



def anno_snp(big_lst,prefix):

    out_fname=prefix+'.anno'
    if os.path.exists(out_fname):
        i=1
        new_out="{0}.{1}".format(out_fname, i)
        while os.path.exists(new_out):
            i+=1
            new_out="{0}.{1}".format(out_fname, i)
        logging.warning("Output file %s exists, writing to %s instead", out_fname, new_out)
        out_file=open(new_out,'w')
    else:
        out_file=open(out_fname,'w')

#    out_file=open(prefix+'.anno','w')
    anno_dic={}

################################################################
#   cds_dx,gene_dx=10,11
#    while True:
#        lin=snp_file.readline().rstrip('\n')
#        if lin=='':
#            break

#       if lin[0] == '#':  ## skip the comment
#            continue

#        lst=lin.split('\t')
#        contig,pos=lst[0].split('.')[0],int(lst[1])
#        allel_ref,allel=lst[3],lst[4].split(',')
#       print("The contig, pos, allele_ref, allele are: ",contig,pos,allel_ref,allel)

###############################################################

    for lst in big_lst:
#        print(lst)
        contig,pos=lst[0],int(lst[1])
#        allel_ref,allel=lst[3],lst[4].split('|')
        allel_ref,allel=lst[4],re.split('[|,:/]{1}',lst[6])
#        print(allel)

        try:
            allel.remove(allel_ref)
            allel_var=allel[0]
        except ValueError:
            allel_var=allel[0]


        anno,gene,geneIDlst,geneAcc,snp_dir,aa_chg,aa_chgPos=[],[],[],[],[],[],[]
        ed0,geneSymb0,geneID0,ori0=0,'','',''

        if contig not in contg_dic:
            anno,gene,geneIDlst,geneAcc,snp_dir,aa_chg,aa_chgPos=['inter-gene'],[''],[''],[''],['0'],[''],[''] # if the dic has no contig, the SNP is inter-genic
            funx=anno[0]
            funx_gene=gene[0]
            funx_geneID=geneIDlst[0]
            funx_geneAcc=geneAcc[0]
            map_dir=snp_dir[0]
            funx_aa=aa_chg[0]
            funx_aaPos=aa_chgPos[0]

#           out_file.write('\t'.join(lst)+'\t'+funx+'\t'+funx_gene+'\t'+funx_geneID+'\t'+map_dir+'\t'+funx_aa+'\n')
			
            new_ky=(contig,pos)
            if new_ky not in anno_dic:
                anno_dic[new_ky]=['\t'.join(lst),funx,funx_gene,funx_geneID,map_dir,funx_aa,funx_geneAcc,allel_var,funx_aaPos]
            elif new_ky in anno_dic:
                logging.warning("the snp (%s,%d) has existed!" % contig,pos)

            continue

#       sort_contgDic=sorted(contg_dic[contig])+[]]
        for contg_val in sorted(contg_dic[contig]):

            flag,flag0=0,0
            st,ed,ori=contg_val[0],contg_val[1],contg_val[2]
            geneID,geneSymb=contg_val[3],contg_val[4]

            if ed0 < pos < st: # intergenic
                the_anno,the_gene,the_geneID,the_geneAcc='inter-gene','','',''

#                if pos-ed0 <= 500 and ori == '+':
#                    the_anno+='-3'
#                    the_gene+=geneSymb0+'||'
#                    the_geneID+=geneID0+'||'

#                if pos-ed0 <= 1000 and ori == '-':
#                    the_anno+='-5'
#                    the_gene+=geneSymb0+'||'
#                    the_geneID+=geneID0+'||'
            
                if ori0 =='+':
                    the_anno+='-3'
                    the_gene+=geneSymb0+'||'
                    the_geneID+=geneID0+'||'
                if ori0 =='-':
                    the_anno+='-5'
                    the_gene+=geneSymb0+'||'
                    the_geneID+=geneID0+'||'
                if ori0 =='':
                    the_anno+='-'
                    the_gene+=geneSymb0+'||'
                    the_geneID+=geneID0+'||'



#                if st-pos <= 1000 and ori == '+':
#                    the_anno+='-5'
#                    the_gene+=geneSymb+'||'
#                    the_geneID+=geneID+'||'

#                if st-pos <= 500 and ori == '-':
#                    the_anno+='-3'
#                    the_gene+=geneSymb+'||'
#                    the_geneID+=geneID+'||'

                if ori =='+':
                    the_anno+='-5'
                    the_gene+=geneSymb+'||'
                    the_geneID+=geneID+'||'
                if ori =='-':
                    the_anno+='-3'
                    the_gene+=geneSymb+'||'
                    the_geneID+=geneID+'||'
                if ori =='':
                    the_anno+='-'
                    the_gene+=geneSymb+'||'
                    the_geneID+=geneID+'||'



                anno+=[cod[the_anno]]
                #anno+=[the_anno]
                gene+=[the_gene]
                geneIDlst+=[the_geneID]
                geneAcc+=[the_geneAcc]
                snp_dir+=['0']
                aa_chg+=['']
                aa_chgPos+=['']

                break # if intergenic,stop searching

            if st <= pos <= ed and geneID not in gene_dic: # the SNP is in some psudo-gene or non-coding RNA

                anno+=[cod['pseudo-gene|ncRNA']]
                #anno+=['pseudo-gene']
                gene+=[geneSymb]
                geneIDlst+=[geneID]
                geneAcc+=['']
                snp_dir+=[ori]
                aa_chg+=['']
                aa_chgPos+=['']
                flag,flag0=1,1
		
            elif st <= pos <= ed and geneID in gene_dic: # in some coding gene

                gene_ed0=st-1
                sort_geneDic=sorted(gene_dic[geneID])+[[ed,ed,ori,'UTR','-']] # include the last UTR-fragment

                for geneID_val in sort_geneDic:
                    gene_st,gene_ed,gene_dir=geneID_val[0],geneID_val[1],geneID_val[2]
                    gene_reg,gene_Acc=geneID_val[3],geneID_val[4]

                    if gene_ed0 < pos < gene_st: # intronic

                        the_anno,the_gene,the_geneID,the_geneAcc='intron',geneSymb,geneID,gene_Acc
                        if gene_st-pos  <= 2:
                            the_anno+='-splice-3'
                        if pos-gene_ed0 <= 2:
                            the_anno+='-splice-5'

                        anno+=[cod[the_anno]]
                        #anno+=[the_anno]
                        gene+=[the_gene]
                        geneIDlst+=[the_geneID]
                        geneAcc+=[the_geneAcc]
                        snp_dir+=[gene_dir]
                        aa_chg+=['']
                        aa_chgPos+=['']
                        #break # if intronic, stop searching
			

                    if gene_st <= pos <= gene_ed: # exonic

                        the_anno,the_gene,the_geneID,the_geneAcc=gene_reg,geneSymb,geneID,gene_Acc

                        if the_anno != 'CDS':
                            anno+=[cod[the_anno]]
                            #anno+=[the_anno]
                            gene+=[the_gene]
                            geneIDlst+=[the_geneID]
                            geneAcc+=[the_geneAcc]
                            snp_dir+=[gene_dir]
                            aa_chg+=['']
                            aa_chgPos+=['']

                        else:
                            CDS_anno,Prot_anno,AApos_anno=find_cod_multi(contig,geneID,gene_Acc,[pos],[allel_var])
                            anno+=[cod[CDS_anno]]
                            gene+=[the_gene]
                            geneIDlst+=[the_geneID]
                            geneAcc+=[the_geneAcc]
                            snp_dir+=[gene_dir]
                            aa_chg+=[Prot_anno]
                            aa_chgPos+=[AApos_anno]


                    gene_ed0=gene_ed

                flag,flag0=1,1


            elif pos > ed: # if the location is in the last inter-gene region of that contig,stop
                flag=0


            if flag==0 and flag0==1:
                break
					

            ed0=ed
            geneSymb0=geneSymb
            geneID0=geneID
            ori0=ori

#       if contig=='NT_167214' and pos==116855:
#           print anno,gene,geneIDlst	

        zipped=list(zip(anno,gene,geneIDlst,snp_dir,aa_chg,geneAcc,aa_chgPos))
        zipped.sort()

        if zipped!=[]:

            zipped0=zipped[0]
            funx=antiCod[zipped0[0]]
            funx_gene=zipped0[1]
            funx_geneID=zipped0[2]
            map_dir=zipped0[3]
            funx_aa=zipped0[4]
            funx_geneAcc=zipped0[5]
            funx_aaPos=zipped0[6]

        else:

            funx='inter-gene'
            funx_gene=''
            funx_geneID=''
            map_dir='0'
            funx_aa=''
            funx_geneAcc=''
            funx_aaPos=''

#       out_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' % ('\t'.join(lst),funx,funx_gene,funx_geneID,map_dir,funx_aa) )

		
        new_ky=(contig,pos)
        if new_ky not in anno_dic:
            anno_dic[new_ky]=['\t'.join(lst),funx,funx_gene,funx_geneID,map_dir,funx_aa,funx_geneAcc,allel_var,funx_aaPos]
        elif new_ky in anno_dic:
            logging.warning("the snp (%s,%d) has existed!"  % (contig,pos))

#       if lst[cds_dx].find('CDS')==-1:
#           lst[cds_dx]=funx
#           lst[gene_dx]=funx_gene
#           lst+=[map_dir]
#       elif lst[cds_dx].find('CDS')!=-1:
#           lst+=[map_dir]
#       out_file.write('\t'.join(lst)+'\n')		

#    snp_file.close()


    anno_kys=list(anno_dic.keys())
    anno_kys.sort()
    kysLg=len(anno_kys)

## initial values for each key and value	
    ech_ky0=anno_kys[0]
    ech_val0=anno_dic[ech_ky0]

    ech_contig0,ech_pos0=ech_ky0[0],ech_ky0[1]
    ech_geneID0,ech_geneAcc0,ech_aaPos0=ech_val0[3],ech_val0[6],ech_val0[8]
    print(ech_contig0,ech_pos0)

## checking neighboring SNPs within the same codon:
    infram_lst,ech_infram=[],[]
    for kk in range(1,kysLg):

        ech_ky=anno_kys[kk]
        ech_contig,ech_pos=ech_ky[0],ech_ky[1]
        ech_val=anno_dic[ech_ky]
        ech_geneID,ech_geneAcc,ech_aaPos=ech_val[3],ech_val[6],ech_val[8]
        
        if ech_aaPos != '' and ech_contig == ech_contig0 and ech_geneID == ech_geneID0 and ech_geneAcc == ech_geneAcc0 and ech_aaPos == ech_aaPos0:
            ech_infram+=[ech_ky0]
            ech_infram+=[ech_ky]
			
        else:
            if ech_infram != []:
                ech_infram=list(set(ech_infram))
                ech_infram.sort()
                infram_lst+=[ech_infram]
                ech_infram=[]

        ech_ky0=ech_ky
        ech_contig0,ech_pos0=ech_ky0[0],ech_ky0[1]
        ech_val0=ech_val
        ech_geneID0,ech_geneAcc0,ech_aaPos0=ech_val0[3],ech_val0[6],ech_val0[8]


## extract the positions and alleles of SNPs within the same codon, and call find_cod_multi for functional identification
    for ech_infram in infram_lst:

        pos_infram,allel_infram=[],[]
        inframLg=len(ech_infram)
        comm_contig=ech_infram[0][0]
        val0=anno_dic[ech_infram[0]]
        comm_geneID,comm_geneAcc=val0[3],val0[6]
        
        for pp in range(0,inframLg):
            ech_ky=ech_infram[pp]
            ech_val=anno_dic[ech_ky]
            ech_pos=ech_ky[1]
            ech_allel=ech_val[7]
            pos_infram+=[ech_pos]
            allel_infram+=[ech_allel]

        CDS_anno,Prot_anno,AApos_anno=find_cod_multi(comm_contig,comm_geneID,comm_geneAcc,pos_infram,allel_infram)
        for pp in range(0,inframLg):
            ech_ky=ech_infram[pp]
            anno_dic[ech_ky][1]=CDS_anno
            anno_dic[ech_ky][5]=Prot_anno

	
    for ech in anno_kys:
        anno_val=anno_dic[ech]
        out_file.write('\t'.join(anno_val[0:6])+'\t'+str(anno_val[-1])+'\n')

    out_file.close()



def snps(args):

    if args.contg:
        contg_filefile = args.contg
        create_contg_dic(contg_filefile)
#        print(contg_dic.keys())
    else:
        raise Exception("Please provide the contig dictionary file !")


    if args.gene:
        gene_filefile = args.gene
        create_gene_dic(gene_filefile)
    else:
        raise Exception("Please provide the gene dictionary file !")


    if args.CDS:
        CDS_filefile = args.CDS
        create_CDS_dic(CDS_filefile)
    else:
        raise Exception("Please provide the CDS dictionary file !")


    if args.CDS_seq:
        CDS_seq_filefile = args.CDS_seq
        create_CDSseq_dic(CDS_seq_filefile)
    else:
        raise Exception("Please provide the CDS-Seq dictionary file !")


    if args.mol:
        mol_filefile = args.mol
        create_mol_dic(mol_filefile)
    else:
        pass



    coden_dic()

#    if args.i:
#        file_list = glob.glob(args.i.rstrip("/")+"/*")
#    else:
#        raise Exception("Please provide the directory of mutation files")


# Handle args.i as either a directory or file(s)
    if args.input:
        file_list = []
        input_paths = args.input if isinstance(args.input, list) else [args.input]
        print(input_paths) 
        for input_path in input_paths:
            input_path = input_path.rstrip("/")
            if os.path.isdir(input_path):
                # If it's a directory, get all files in it
                file_list.extend(glob.glob(os.path.join(input_path, "*")))
            elif os.path.isfile(input_path):
                file_list.append(input_path)
            else:
                logging.warning("Path does not exist: %s" % input_path)
    else:
        raise Exception("Please provide the directory or file(s) of mutation files")


#    if args.table is not None and args.vcf is None:
    if args.format == "table":
#        tables = args.table
        for ech_file in file_list:
            snpMat = format_var_file(ech_file,"table")
            anno_snp(snpMat,ech_file)
#    elif args.vcf is not None and args.table is None:
    elif args.format == "vcf":
#        tables = args.vcf
        for ech_file in file_list:
            snpMat = format_var_file(ech_file,"vcf")
            anno_snp(snpMat,ech_file)

    else:
        raise Exception("Please provide either vcf or tab-delimited SNP files!")

    
    



