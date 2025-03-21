import os
import glob
import logging


def read_big_file(file0,lines):
    bigSize = 1000000000
    file_lst=[]
    siz=os.path.getsize(file0)
    if siz < bigSize:
        file_lst=[file0]
    else:
        splt_cmd='split --lines=%d %s %s' % (lines,file0,file0+'.sub')
        os.system(splt_cmd)
        file_lst=glob.glob(file0+'.sub*')
    return file_lst

def creat_dic(som_file,flag):
    som_dic={}
    file_lst=read_big_file(som_file,1000000)
    for sub in file_lst:
        sub_file=open(sub,'r')
        sub_file.readline()

        lst=sub_file.read().split('\n')[:-1]
        if flag=='all':
            for ss in lst:
                tt=ss.split('\t')
                try:
                    ky=(tt[0],int(tt[1]),tt[4],tt[6],tt[10])
                except IndexError:
                    logging.error('The index does not match, please re-define your key.')
					
                val=tt
                val[7],val[8],val[9]='-','-','-'
                if ky not in som_dic: # and abs(val[0])>=2.0:
                    som_dic[ky]=val
                elif ky in som_dic: # and abs(val[0])>=2.0:
                    logging.warning('The key ' + ky + ' has existed!')

        else:
            tab=int(flag)
            for ss in lst:
                tt=ss.split('\t')
                try:
                    ky=(tt[0],int(tt[1]))
                except IndexError:
                    logging.error('The index does not match, please re-define your key.')

                val=tt[tab]
                if ky not in som_dic:
                    som_dic[ky]=val
                elif ky in som_dic:
                    logging.warning('the key' + ky + 'has existed!')
					
        sub_file.close()
    return som_dic



##creat the total dictionary
def creat_totDic(input_file):

    tot_dic={}
    totNum = len(input_file)

    if totNum == 2:  ## two parameters
        file_list = glob.glob(input_file[0].rstrip("/")+"/*")
        file_flag=input_file[1]
        
        for file_nm in file_list:
            cunt_rept=0
            i_dic = creat_dic(file_nm,file_flag)
            for i_ky in i_dic.keys():
                if i_ky not in tot_dic:
                    tot_dic[i_ky] = i_dic[i_ky]
                elif i_ky in tot_dic:
                    cunt_rept += 1

        logging.info("The total number of redundant elements is %d" % cunt_rept)

    elif totNum > 2 and totNum % 2 == 0:  ## even number of parameters
        for ii in range(0,totNum,2):
            cunt_rept=0
            file_nm = input_file[ii]
            file_flag = input_file[ii+1]
            i_dic = creat_dic(file_nm,file_flag)

            for i_ky in i_dic.keys():
                if i_ky not in tot_dic:
                    tot_dic[i_ky] = i_dic[i_ky]
                elif i_ky in tot_dic:
                    cunt_rept += 1

        logging.info("The total number of redundant elements is %d" % cunt_rept)

    return tot_dic



##creat the dictionary list for the data to be extracted
def ins_multiDic(insfile):

    dic_lst,val_lg_lst=[],[]
    insNum = len(insfile)

    if insNum == 2:  ## Two parameters
        insfile_list = glob.glob(insfile[0].rstrip("/")+"/*")
        file_flag=insfile[1]

        for file_nm in insfile_list:
            i_dic = creat_dic(file_nm,file_flag)
            dic_lst.append(i_dic)
            try:
                i_dic_ky0 = list(i_dic.keys())[0]
                val_lg_lst.append(len(i_dic[i_dic_ky0]))
            except IndexError:
                val_lg_lst.append(3)
                pass

    elif insNum > 2 and insNum % 2 == 0:

        for jj in range(0,insNum,2):
            file_nm = insfile[jj]
            file_flag = insfile[jj+1]
            i_dic = creat_dic(file_nm,file_flag)
            dic_lst.append(i_dic)
            try:
                i_dic_ky0 = list(i_dic.keys())[0]
                val_lg_lst.append(len(i_dic[i_dic_ky0]))
            except IndexError:
                val_lg_lst.append(3)  # if some dictionary is empty, just randomly assign a column-number for write-out
                pass 

    return dic_lst,val_lg_lst


## write the results to file
def writeCombined(insfile_list,totDic,dicList,val_lg_lst,outFileName):

    #out_str = '# RefName\tLoc\tType\tRefLen\tRef\tAlleleNum\tAlt\tAltRatio\tAltDepth\tDepth\tMutType\tGene\tGeneID\tOri\tAAchange\tAAsite\n'
    out_str = '#RefName\tLoc\tType\tRefLen\tRef\tAlleleNum\tAllele\tAlleleRatio\tAlleleDepth\tDepth\tMutType\tGene\tGeneID\tOri\tAAchange\tAAsite\n'

    tot_keys = list(totDic.keys())
    tot_keys.sort()
    for kyky in tot_keys:
        out_str+='\t'.join(totDic[kyky])+'\n'

    file_out=open(outFileName,'w')
    file_out.write(out_str)
    file_out.close()



def SNPmerge(args):

    if args.comb == None:
        raise Exception("Please provide the files to be combined.")
    elif len(args.comb) % 2 > 0:
        raise Exception("Please provide correct number of parameters for -comb")
    elif len(args.comb) == 2 and os.path.exists(args.comb[0]) == False:
        raise Exception("The directory %s does not exist!" % args.comb[0])
    else:
        tot_dic=creat_totDic(args.comb)


    if args.ins == None:
        raise Exception("Please provide the files to be inserted.")
    elif len(args.ins) % 2 > 0:
        raise Exception("Please provide correct number of parameters for -ins")
    elif len(args.ins) == 2 and os.path.exists(args.ins[0]) == False:
        raise Exception("The directory %s does not exist!" % args.ins[0])
    else:
        dic_lst,val_lg_lst=ins_multiDic(args.ins)


    if args.out == None:
        logging.basicConfig(
                format='%(asctime)s \tFile \"%(filename)s" %(levelname)s: \n %(message)s \n',
                datefmt='%a, %d %b %Y %H:%M:%S',
                filename=args.out+'.log', filemode="w", level=logging.INFO)

        raise Exception("Please provide the output file.")
    else:
        writeCombined(args.ins,tot_dic,dic_lst,val_lg_lst,args.out)




