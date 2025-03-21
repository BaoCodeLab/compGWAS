import sys
import re
import logging



def str_rc(str1):
    new_str=list(str1[::-1])
    for i in range(len(str1)):
        if new_str[i]=='A':
            new_str[i]='T'
        elif new_str[i]=='T':
            new_str[i]='A'
        elif new_str[i]=='C':
            new_str[i]='G'
        elif new_str[i]=='G':
            new_str[i]='C'

    return ''.join(new_str)



def parseSeq(the_str):
    seq=''
    orig=the_str.find('ORIGIN')
    st=the_str.find('1',orig)
    ed=the_str.find('//',orig)
    lst=the_str[st:ed].split('\n')
    for itm in lst:
        seq+=''.join(itm.split()[1:])
    return seq



def parseGeneLoc(the_str):
    cflag,cdir=0,'+'
    locLst=[]

    keywd=re.compile(r"\n\s+\/")
    endPos=keywd.search(the_str)
    locStr=the_str[:endPos.start()]+'\n'
    region=locStr.split()[0]
    if region=='gene':
        region='GENE'
    locStr=the_str[:endPos.start()]+'\n'

    locStr=locStr.replace('<','')
    locStr=locStr.replace('>','')

    if locStr.find('complement(')!=-1:
        cdir='-'
        cflag=1
    elif locStr.find('join(')!=-1 or locStr.find('order(')!=-1:
        cflag=1
    if cflag==0:
        g1,g2=' ','\n'
    if cflag==1:
        g1,g2='(',')'

    loc_st=locStr.rfind(g1)
    loc_ed=locStr.find(g2,1)
    allLoc=locStr[loc_st+1:loc_ed].split(',')

    for echLoc in allLoc:
        locLst+=[list(map(int,echLoc.split('..')))]
    return cdir,locLst,region



def parseMarkLoc(the_str,keywd0,cflag=0):
    keywd='     '+keywd0
    MarkFlag=the_str.find(keywd)



def gbkFormat(gbk_file,tab_file,seq_file,theOrg):

    big_sz=100000000

    file_in=open(gbk_file,'r')
    file_out=open(tab_file,'w')
    file_seq=open(seq_file,'w')
    

    file_out.write("# taxomony\tr1\tr2\tr3\tr4\tstart\tend\tlength\tstrand\tsymbol\tgeneID\tmol_type\tfigID\tprotID\tfunction\tfigFam\tsubsystem1\tsubsystem2\subsystem3\n")
    gbk0=''
    cunt=0
    org=theOrg
    pat=re.compile(r"\n     (gene|mRNA|CDS|misc_RNA|rRNA|tRNA|ncRNA)     ")


    chrID,chr_st,chr_ed,chr_dir='-','-','-','-'
    while True:
        gbk1=file_in.read(big_sz)
        cunt+=1
        if gbk1=='':
            gbk=gbk0+'LOCUS'
        else:
            gbk=gbk0+gbk1
        st=gbk.find('LOCUS',0)
        while True:
            ed=gbk.find('LOCUS',st+1)
            if ed==-1:
                gbk0=gbk[st:]
                break


#get the scaffold sequence
            seq_st=gbk.find('\nORIGIN',st+1)
            seq_ed=gbk.find('//\n',st+1)
            seq_tmp=gbk[seq_st:seq_ed]+'//'
            the_seq=parseSeq(seq_tmp)


#get the contig accession
            contgAcc_ed=gbk.find('\n',st+1)
            contg=gbk[st:contgAcc_ed].split()[1]


            feat=gbk.find('\nFEATURES',st+1)
            the_blck=gbk[feat:seq_st]


#find each feature
            feat_st=[]
		
            for ech_pos in pat.finditer(the_blck):
                feat_st+=[ech_pos.start()]
            if feat_st==[]:
                logging.warning('the contig %s has no genes' % contg)
            feat_st+=[len(the_blck)]


#initialize the region and location
            prevReg,prevPortID,prevLoc='','',[]
            for tt in range(1,len(feat_st)):
                region=''
                thePrt=the_blck[feat_st[tt-1]:feat_st[tt]]


                geneSymb_st0=thePrt.find('/gene="',0)
                if geneSymb_st0==-1:
                    geneSymb='-'
                else:
                    geneSymb_st=thePrt.find('"',geneSymb_st0)
                    geneSymb_ed=thePrt.find('"',geneSymb_st+1)
                    geneSymb=thePrt[geneSymb_st+1:geneSymb_ed]


                geneID_st0=thePrt.find('/locus_tag=',0)	
                if geneID_st0==-1:
                    geneID='-'
                else:
                    geneID_st=thePrt.find('"',geneID_st0)
                    geneID_ed=thePrt.find('"',geneID_st+1)
                    geneID=thePrt[geneID_st+1:geneID_ed]


                if thePrt.find('/transcript_id="',0)==-1:
                    protID_st0=thePrt.find('/protein_id="',0)
                else:
                    protID_st0=thePrt.find('/transcript_id="',0)
                if protID_st0==-1:
                    protID='-'
                else:
                    protID_st=thePrt.find('"',protID_st0)
                    protID_ed=thePrt.find('"',protID_st+1)
                    protID=thePrt[protID_st+1:protID_ed]


                prod_st0=thePrt.find('/product=',0)
                if prod_st0==-1:
                    prod='-'
                else:
                    prod_st=thePrt.find('"',prod_st0)
                    prod_ed=thePrt.find('"',prod_st+1)
                    prod=thePrt[prod_st+1:prod_ed].replace('\n'+21*' ',' ')


#find the codon start

                contg_dir,theLoc,region=parseGeneLoc(thePrt)
                for regLoc in theLoc:
                    if len(regLoc)==1: # it is 1-base exon
                        regLoc_dx=theLoc.index(regLoc)
                        theLoc[regLoc_dx]=[regLoc[0],regLoc[0]]

#5' UTR
                if region=='CDS' and prevReg=='mRNA':
                    for prevRegLoc in prevLoc:
                        if  prevRegLoc[0] < theLoc[0][0] <= prevRegLoc[1]:
                            file_out.write('%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\t\t\t\n' % (org,chrID,chr_st,chr_ed,chr_dir,contg,prevRegLoc[0],theLoc[0][0]-1,contg_dir,geneSymb,geneID,'UTR','Primary_Assembly',prevProtID,prod) )
                            break
                        elif theLoc[0][0] > prevRegLoc[1]:
                            file_out.write('%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\t\t\t\n' % (org,chrID,chr_st,chr_ed,chr_dir,contg,prevRegLoc[0],prevRegLoc[1],contg_dir,geneSymb,geneID,'UTR','Primary_Assembly',prevProtID,prod) )


#3' UTR			
                if region=='CDS' and prevReg=='mRNA':
                    prevLoc.reverse()
                    for prevRegLoc in prevLoc:
                        if prevRegLoc[0] <= theLoc[-1][-1] < prevRegLoc[1]:
                            file_out.write('%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\t\t\t\n' % (org,chrID,chr_st,chr_ed,chr_dir,contg,theLoc[-1][-1]+1,prevRegLoc[1],contg_dir,geneSymb,geneID,'UTR','Primary_Assembly',prevProtID,prod) )
                            break
                        elif theLoc[-1][-1] < prevRegLoc[0]:
                            file_out.write('%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\t\t\t\n' % (org,chrID,chr_st,chr_ed,chr_dir,contg,prevRegLoc[0],prevRegLoc[1],contg_dir,geneSymb,geneID,'UTR','Primary_Assembly',prevProtID,prod) )


#mRNA/CDS 			
                prt_seq=''
                for regLoc in theLoc:
                    contg_st,contg_ed=regLoc[0],regLoc[1]
                    file_out.write('%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\t\t\t\n' % (org,chrID,chr_st,chr_ed,chr_dir,contg,contg_st,contg_ed,contg_dir,geneSymb,geneID,region,'Primary_Assembly',protID,prod) )
                    prt_seq+=the_seq[contg_st-1:contg_ed]
			
                who_seq=prt_seq.upper()
                if contg_dir=='+' and region=='GENE':
                    file_seq.write('>%s|%s\n%s\n' % (geneID,geneSymb,who_seq) )
                elif contg_dir=='-' and region=='GENE':
                    file_seq.write('>%s|%s\n%s\n' % (geneID,geneSymb,str_rc(who_seq)) )


                prevReg,prevProtID=region,protID
                prevLoc=theLoc

            st=ed

        if gbk1=='':
            break

    file_in.close()
    file_out.close()
    file_seq.close()
		

		
def gbk2seqGene(args):

    if args.gbk == None:
        raise Exception("Please provide the gbk file.")
    else:
        gbk_file=args.gbk

    if args.out == None:
        logging.basicConfig(
                format='%(asctime)s \tFile \"%(filename)s" %(levelname)s: \n %(message)s \n',
                datefmt='%a, %d %b %Y %H:%M:%S',
                filename=args.out+'.log', filemode="w", level=logging.INFO)
        raise Exception("Please provide the file for tab-delimted output.")
    else:
        out_file=args.out

    if args.seq == None:
        raise Exception("Please provide the file for protein sequences output.")
    else:
        seq_file=args.seq

    if args.org == None:
        logging.warning("No organism is provided, default NULL value will be used.")
        Org='-'
    else:
        Org = args.org

    gbkFormat(gbk_file,out_file,seq_file,Org)






