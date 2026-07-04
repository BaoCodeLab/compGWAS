import argparse
import glob
from allGWAS.preGWASlib import gbk2seqGene
from allGWAS.preGWASlib import parseGpff
from allGWAS.GWASlib.SNPgwas import SNPgwas
from allGWAS.GWASlib.CDSgwas import CDSgwas
from allGWAS.GWASlib.nonCDSgwas import nonCDSgwas
from allGWAS.GWASannolib import SNPmerge
from allGWAS.GWASannolib import SNPCDSanno
from allGWAS.GWASannolib.nonCDSanno import nonCDSanno
from allGWAS.GWASLDlib import LDprun
from allGWAS.GWASLDlib import Block
from allGWAS.GWASLDlib import Screen
from allGWAS.GWASfilterlib.SCfilter import SCfilter
from allGWAS.GWASfilterlib.SNfilter import SNfilter
from allGWAS.GWASfilterlib.Cfilter import Cfilter
from allGWAS.GWASfilterlib.Nfilter import Nfilter



if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(title="GWAS",description="process of GWAS",help="The whole process of GWAS")
    
   

    preGWAS_parser=subparsers.add_parser("preGWAS",help="Get the tab-delimted format gene annotation file , the CDS sequence file, and the annotation dictionaries of the species to be analyzed")
    preGWAS_parser.add_argument("-F","--FLAG", nargs='?', const=1, default="preGWAS", type=str, help="ignore this argument or input string 'preGWAS'")
    preGWAS_parser.add_argument("-g","--gbk",required=True,help="The gbk file for converting to tab-delimted format annotation")
    preGWAS_parser.add_argument("-o","--out",required=True,help="The output file for tab-delimted format annotation file of all genes")
    preGWAS_parser.add_argument("-s","--seq",required=True,help="The protein sequences file from gbk file")
    preGWAS_parser.add_argument("-t","--org",help="The organism you are working on")
    preGWAS_parser.add_argument("-r","--ref",required=True,help="The reference genome file")
    preGWAS_parser.add_argument("-p","--prefix",required=True,help="The prefix of output file")


    
    SNPgwas_parser=subparsers.add_parser("SNPgwas",help="GWAS of SNPs") 
    SNPgwas_parser.add_argument("-S","--allSNP",required=True,help="the path of all strains' SNP files, not include 'SNP file name'")
    SNPgwas_parser.add_argument("-p","--pheno1",nargs=2, required=True,help="input phenotype1 strains' GCAs (only the number between GCA_ and .), and its name of the phenotype p")
    SNPgwas_parser.add_argument("-P","--pheno0",nargs=2, required=True,help="input phenotype0 strains' GCAs (only the number between GCA_ and .), and its name of the phenotype P")
    SNPgwas_parser.add_argument("-i","--info", nargs=3, help="input strain's covariates information table, the covariates for Logistic(seperated by ,), and the types of each covariate(f or n, f means factor and n means numeric, seperated by ,)")
    SNPgwas_parser.add_argument("-f","--ref",nargs=2,type=str,help="input reference strain's GCAs (only the number between GCA_ and .) and its phenotype(p or P)")
    SNPgwas_parser.add_argument("-c","--cols",type=int,nargs=3,required=True,help="input the column index of RefName (Contig), Loc and Allele of SNPs annotation file (start with 0, e.g first column = 0)")
    SNPgwas_parser.add_argument("-t","--thread", nargs='?', const=1, default=1, type=int, help="input the thread used in 'CHI2distribution.R' and 'GWAS.R' (default: 1)")
    SNPgwas_parser.add_argument("-T","--threshold", nargs='?', const=1, default=0.05, type=float, help="input the threshold Pvalue of Allele used to filter logisticRegre results (default: 0.05)")
    SNPgwas_parser.add_argument("-o","--outdir", nargs='?', const=1, default=".", help="output folder (default: .)")
    SNPgwas_parser.add_argument("-O","--prefix", nargs='?', const=1, default="SNP.output", help="filename output prefix of SNP results (default: SNP.output)")
    SNPgwas_parser.add_argument("-R","--Rscript",required=True,help="the path of Rscript, include 'Rscript'")
    SNPgwas_parser.add_argument("-r","--rscript",required=True,help="the path of the GWAS package")
    SNPgwas_parser.set_defaults(func=SNPgwas)



    CDSgwas_parser=subparsers.add_parser("CDSgwas",help="GWAS of CDS")
    CDSgwas_parser.add_argument("-D","--allDIP",required=True,help="the path of all strains' InDel files, not include 'InDel file name'")
    CDSgwas_parser.add_argument("-S","--allSNP",required=True,help="the path of all strains' SNP files, not include 'SNP file name'")
    CDSgwas_parser.add_argument("-C","--CDSdic",nargs='+',required=True,help="input the CDS annotation dictionaries of the species to be analyzed")
    CDSgwas_parser.add_argument("-g","--genome",required=True,help="input the FASTA format genome file of the species to be analyzed")
    CDSgwas_parser.add_argument("-p","--pheno1",nargs=2, required=True,help="input phenotype1 strains' GCAs (only the number between GCA_ and .), and its name of the phenotype p")
    CDSgwas_parser.add_argument("-P","--pheno0",nargs=2, required=True,help="input phenotype0 strains' GCAs (only the number between GCA_ and .), and its name of the phenotype P")
    CDSgwas_parser.add_argument("-i","--info", nargs=3, help="input strain's covariates information table, the covariates for Logistic(seperated by ,), and the types of each covariate(f or n, f means factor and n means numeric, seperated by ,)")
    CDSgwas_parser.add_argument("-f","--ref",nargs=2,type=str,help="input reference strain's GCAs (only the number between GCA_ and .) and its phenotype(p or P)")
    CDSgwas_parser.add_argument("-l","--length",nargs=3,type=int,required=True,help="input the length extending outward, the length extending inward at the 5'end of CDS, and the length extending inward at the 3'end of CDS")
    CDSgwas_parser.add_argument("-c","--cols",type=int,nargs=6,required=True,help="input the column index of RefName (Contig), Loc and Allele of SNPs annotation file and InDels annotation file (start with 0, e.g first column = 0)")
    CDSgwas_parser.add_argument("-t","--thread", nargs='?', const=1, default=1, type=int, help="input the thread used in 'CHI2distribution.R' and 'GWAS.R' (default: 1)")
    CDSgwas_parser.add_argument("-T","--threshold", nargs='?', const=1, default=0.05, type=float, help="input the threshold Pvalue of Allele used to filter logisticRegre results (default: 0.05)")
    CDSgwas_parser.add_argument("-o","--outdir", nargs='?', const=1, default=".", help="output folder (default: .)")
    CDSgwas_parser.add_argument("-O","--prefix", nargs='?', const=1, default="CDS.output", help="filename output prefix of CDS results (default: CDS.output)")
    CDSgwas_parser.add_argument("-R","--Rscript",required=True,help="the path of Rscript, include 'Rscript'")
    CDSgwas_parser.add_argument("-r","--rscript",required=True,help="the path of the GWAS package")
    CDSgwas_parser.set_defaults(func=CDSgwas)



    nonCDSgwas_parser=subparsers.add_parser("nonCDSgwas",help="GWAS of nonCDS")
    nonCDSgwas_parser.add_argument("-D","--allDIP",required=True,help="the path of all strains' InDel files, not include 'InDel file name'")
    nonCDSgwas_parser.add_argument("-C","--CDSdic",nargs='+',required=True,help="input the CDS annotation dictionaries of the species to be analyzed")
    nonCDSgwas_parser.add_argument("-Z","--size",type=int,required=True,help="the size of the window")
    nonCDSgwas_parser.add_argument("-p","--pheno1",nargs=2,required=True,help="input phenotype1 strains' GCAs (only the number between GCA_ and .), and its name of the phenotype p")
    nonCDSgwas_parser.add_argument("-P","--pheno0",nargs=2,required=True,help="input phenotype0 strains' GCAs (only the number between GCA_ and .), and its name of the phenotype P")
    nonCDSgwas_parser.add_argument("-i","--info", nargs=3, help="input strain's covariates information table, the covariates for Logistic(seperated by ,), and the types of each covariate(f or n, f means factor and n means numeric, seperated by ,)")
    nonCDSgwas_parser.add_argument("-f","--ref",nargs=2,type=str,help="input reference strain's GCAs (only the number between GCA_ and .) and its phenotype(p or P)")
    nonCDSgwas_parser.add_argument("-g","--genome",required=True,help="input the FASTA format genome file of the species to be analyzed")
    nonCDSgwas_parser.add_argument("-c","--cols",type=int,nargs=4,required=True,help="input the column index of RefName, Loc, Ref, and Muttype of InDels annotation file (start with 0, e.g first column = 0)")
    nonCDSgwas_parser.add_argument("-t","--thread", nargs='?', const=1, default=1, type=int, help="input the thread used in 'CHI2distribution.R' and 'GWAS.R' (default: 1)")
    nonCDSgwas_parser.add_argument("-T","--threshold", nargs='?', const=1, default=0.05, type=float, help="input the threshold Pvalue of Allele used to filter logisticRegre results (default: 0.05)")
    nonCDSgwas_parser.add_argument("-o","--outdir", nargs='?', const=1, default=".", help="output folder (default: .)")
    nonCDSgwas_parser.add_argument("-O","--prefix", nargs='?', const=1, default="nonCDS.output", help="filename output prefix of nonCDS results (default: nonCDS.output)")
    nonCDSgwas_parser.add_argument("-R","--Rscript",required=True,help="the path of Rscript, include 'Rscript'")
    nonCDSgwas_parser.add_argument("-r","--rscript",required=True,help="the path of the GWAS package")
    nonCDSgwas_parser.set_defaults(func=nonCDSgwas)



    SNPCDSanno_parser=subparsers.add_parser("SNPCDSanno",help="The annotations of SNP or CDS GWAS results") 
    SNPCDSanno_parser.add_argument("-F","--FLAG", nargs='?', const=1, default="SNPCDSanno", type=str, help="ignore this argument or input string 'SNPCDSanno'")
    SNPCDSanno_parser.add_argument("-c","--comb",nargs="+",help="When annotation type is SNP, input the path of the files to be merged and the column to be combined, for example: 'allSNP/ all'")
    SNPCDSanno_parser.add_argument("-i","--ins",nargs="+",help="When annotation type is SNP, input the path of the files to be merged and the column to be combined, for example: 'allSNP/ all'")
    SNPCDSanno_parser.add_argument("-o","--out", nargs='?', const=1, default="./SNP.merged.allanno", help="When annotation type is SNP, input the output annotation file of megered all SNPs (default: ./SNP.merged.allanno)")
    SNPCDSanno_parser.add_argument("-l","--logisre",required=True,help="the results file of SNP or CDS logistic regression")
    SNPCDSanno_parser.add_argument("-a","--anno",help="When annotation type is CDS, input the annotation file of all genes of the species to be analyzed")
    SNPCDSanno_parser.add_argument("-t","--type",required=True,help="the type of annotation, SNP or CDS")



    nonCDSanno_parser=subparsers.add_parser("nonCDSanno",help="The annotations of nonCDS GWAS results")
    nonCDSanno_parser.add_argument("-C","--CDSdic",nargs='+',required=True,help="input the CDS annotation dictionaries of the species to be analyzed")
    nonCDSanno_parser.add_argument("-g","--genome",required=True,help="input the FASTA format genome file of the species to be analyzed")
    nonCDSanno_parser.add_argument("-Z","--size",type=int,required=True,help="the size of the window")
    nonCDSanno_parser.add_argument("-l","--logisre",required=True,help="the results file of nonCDS logistic regression")
    nonCDSanno_parser.set_defaults(func=nonCDSanno)


    
    LDprun_parser=subparsers.add_parser("LDprun",help="Linkage disequilibrium analysis of SNPs") 
    LDprun_parser.add_argument("-F","--FLAG", nargs='?', const=1, default="LDprun", type=str, help="ignore this argument or input string 'LDprun'")
    LDprun_parser.add_argument("-S","--allSNP",required=True,help="the path of all strains' SNP files, not include 'SNP file name'")
    LDprun_parser.add_argument("-p","--pheno1",required=True,help="input phenotype1 strains' GCAs (only the number between GCA_ and .)")
    LDprun_parser.add_argument("-P","--pheno0",required=True,help="input phenotype0 strains' GCAs (only the number between GCA_ and .)")
    LDprun_parser.add_argument("-f","--ref",nargs=2,type=str,help="input reference strain's GCAs (only the number between GCA_ and .) and its phenotype(p or P)")
    LDprun_parser.add_argument("-g","--genome",required=True,help="input the FASTA format genome file of the species to be analyzed")
    LDprun_parser.add_argument("-l","--logisre",required=True,help="the results file of all SNPs logistic regression")
    LDprun_parser.add_argument("-T","--threshold", nargs='?', const=1, default=0.05, type=float, help="input the threshold Pvalue of Allele used to filter SNP logisticRegre results (default: 0.05)")
    LDprun_parser.add_argument("-o","--outdir", nargs='?', const=1, default=".", help="output folder (default: .)")
    LDprun_parser.add_argument("-O","--prefix", nargs='?', const=1, default="SNP.LD.output", help="filename output prefix of SNPs linkage disequilibrium analysis (default: SNP.LD.output)")    
    LDprun_parser.add_argument("-J","--java",required=True,help="the path of java, include 'java'")
    LDprun_parser.add_argument("-H","--Haploview",required=True,help="the path of Haploview.jar, include 'Haploview.jar'")
    LDprun_parser.add_argument("-A","--Argument",nargs='?',const=1,default="-n -dprime -blockoutput GAB -check -chromosome X",help="the parameters of Haploview.jar, excluding input and output parameters (default: '-n -dprime -blockoutput GAB -check -chromosome X')")
    
    LDprun_parser.add_argument("-T2","--threshold2", nargs='?', const=1, default=0.8, type=float, help="the threshold of R^2 to separate the LD snps into blocks (default: 0.8)")
    LDprun_parser.add_argument("-d","--maxdist", nargs='?', const=1, default=10, type=float, help="the max distance of each block, the unit is kb (default: 10)")
    
    LDprun_parser.add_argument("--Chi2",required=True,help="input the chi2 distribution file of the 2 phenotypes")
    LDprun_parser.add_argument("-LA","--LogisAnno",required=True,help="input the annotation file of SNP logistic regression results")  

    

    SCfilter_parser=subparsers.add_parser("SCfilter",help="Filter coding region SNPs GWAS results") 
    SCfilter_parser.add_argument("-s","--screenout",required=True,help="the path of SNPblockID.screenout file(s) of SNP linkage disequilibrium analysis")
    SCfilter_parser.add_argument("-O","--prefix", nargs='?', const=1, default="SNP.LD.output", help="filename output prefix of SNPs linkage disequilibrium analysis (default: SNP.LD.output)")
    SCfilter_parser.add_argument("-a","--anno",required=True,help="the annotation file of all genes of the species to be analyzed")
    SCfilter_parser.add_argument("-T","--threshold", nargs='?', const=1, default=0.05, type=float, help="input the significance threshold Pvalue of Allele used to filter logisticRegre results (default: 0.05)")
    SCfilter_parser.set_defaults(func=SCfilter)

    

    SNfilter_parser=subparsers.add_parser("SNfilter",help="Filter non-coding region SNP GWAS results")
    SNfilter_parser.add_argument("-C","--CDSdic",nargs='+',required=True,help="input the CDS annotation dictionaries of the species to be analyzed")
    SNfilter_parser.add_argument("-LA","--LogisAnno",required=True,help="input the annotation file of SNP logistic regression results")
    SNfilter_parser.add_argument("-a","--anno",required=True,help="the annotation file of all genes of the species to be analyzed")
    SNfilter_parser.add_argument("--Chi2",required=True,help="input the chi2 distribution file of the 2 phenotypes")
    SNfilter_parser.add_argument("-T","--threshold", nargs='?', const=1, default=0.05, type=float, help="input the significance threshold Pvalue of Allele used to filter logisticRegre results (default: 0.05)")
    SNfilter_parser.add_argument("-d","--distance",nargs=2,type=int,help="input the threshold of 5' distance(bp) and 3' distance(bp)")
    SNfilter_parser.set_defaults(func=SNfilter)


    
    Cfilter_parser=subparsers.add_parser("Cfilter",help="Filter CDS GWAS results")
    Cfilter_parser.add_argument("-l","--logisre",required=True,help="the results file of CDS logistic regression")
    Cfilter_parser.add_argument("-a","--anno",required=True,help="the annotation file of all genes of the species to be analyzed")
    Cfilter_parser.add_argument("-c2","--Chi2",required=True,help="input the chi2 distribution file of the 2 phenotypes")
    Cfilter_parser.add_argument("-T","--threshold", nargs='?', const=1, default=0.05, type=float, help="input the significance threshold Pvalue of Allele used to filter logisticRegre results (default: 0.05)")
    Cfilter_parser.set_defaults(func=Cfilter)



    Nfilter_parser=subparsers.add_parser("Nfilter",help="Filter nonCDS GWAS results")
    Nfilter_parser.add_argument("-C","--CDSdic",nargs='+',required=True,help="input the CDS annotation dictionaries of the species to be analyzed")
    Nfilter_parser.add_argument("-LA","--LogisAnno",required=True,help="input the annotation file of SNP logistic regression results")
    Nfilter_parser.add_argument("-a","--anno",required=True,help="the annotation file of all genes of the species to be analyzed")
    Nfilter_parser.add_argument("-c2","--Chi2",required=True,help="input the chi2 distribution file of the 2 phenotypes")
    Nfilter_parser.add_argument("-T","--threshold", nargs='?', const=1, default=0.05, type=float, help="input the significance threshold Pvalue of Allele used to filter logisticRegre results (default: 0.05)")
    Nfilter_parser.add_argument("-d","--distance",nargs=2,type=int,help="input the threshold of 5' distance(bp) and 3' distance(bp)")
    Nfilter_parser.set_defaults(func=Nfilter)



    args = parser.parse_args()

    if "FLAG" in dir(args) and args.FLAG == "preGWAS":
        args.gpff = args.out
        gbk2seqGene.gbk2seqGene(args)
        parseGpff.parseGpff(args)
    elif "FLAG" in dir(args) and args.FLAG == "SNPCDSanno":
        if args.type == "CDS":
            SNPCDSanno.SNPCDSanno(args)
        elif args.type == "SNP":
            args.anno = args.out
            SNPmerge.SNPmerge(args)
            SNPCDSanno.SNPCDSanno(args)
    elif "FLAG" in dir(args) and args.FLAG == "LDprun":
        args.info = glob.glob(args.outdir + "/" + args.prefix + r".*.Haploview.info")
        LDprun.LDprun(args)
        Block.Block(args)
        Screen.Screen(args)
    else:
        args.func(args)        
    print(args)



