import argparse
import glob



if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(title="GWAS",description="process of GWAS",help="The whole process of GWAS")
    

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



