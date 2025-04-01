# compGWAS
An integrated toolkit for comprehensive GWAS analysis of SNPs and Indels to identify the association between phenotypes and genotypes.
* [preGWAS](#preGWAS)

        Prepare the reference genome and annotation files
        
* [GWAS](#GWAS)

        GWAS analysis of SNPs in the coding regions
        GWAS analysis of the synergistic effect of SNPs/Indels on CDS
        GWAS analysis of SNPs/Indels windows in the non-coding regions

* [GWASanno](#GWASanno)

        The annotation of coding genes affected by SNPs/Indels in the GWAS results
        The annotation of SNPs/Indels windows in the non-coding regions in the GWAS results

* [LDprun](#LDprun)

        SNP pruning based on linkage disequilibrium to remove the redundancy among SNPs

* [Filtering](#Filtering)

        Filtering of SNPs in the coding regions based on GWAS and LD results
        Filtering of SNPs in the non-coding regions based on GWAS and LD results
        Filtering of coding genes affected by SNPs/Indels based on GWAS results
        Filtering of Indel windows in the non-coding regions based on GWAS results


# Install
### Requirements
* Python >= 3.7
* R >= 3.5
* Java
* packages for python:

&emsp;&emsp;&emsp;&emsp;numpy 

&emsp;&emsp;&emsp;&emsp;pandas
    
* packages for R:

&emsp;&emsp;&emsp;&emsp;foreach

&emsp;&emsp;&emsp;&emsp;doParallel

&emsp;&emsp;&emsp;&emsp;BaylorEdPsych

* software for linkage disequilibrium analysis:
    
&emsp;&emsp;&emsp;&emsp;Haploview



### Usage 
#### python compGWAS.py [-h] {preGWAS,SNPgwas,CDSgwas,nonCDSgwas,SNPCDSanno,nonCDSanno,LDprun,SCfilter,SNfilter,Cfilter,Nfilter}
#### Command line usage
    -h, --help          Show the help message.
    preGWAS             Get the tab-delimited format gene annotation file, and the annotation dictionaries.
    SNPgwas             GWA analysis of SNPs in the coding regions.
    CDSgwas             GWA analysis of the synergistic effect of SNPs/Indels in the coding regions.
    nonCDSgwas          GWA analysis of SNPs/Indels windows in the non-coding regions.
    SNPCDSanno          The annotations of coding genes affected by SNPs/Indels in the GWAS results.
    nonCDSanno          The annotations of SNPs/Indels windows in the non-coding regions in the GWAS results.
    LDprun              Linkage disequilibrium analysis of SNPs.
    SCfilter            Filtering of SNPs in the coding regions based on GWAS and LD results.
    SNfilter            Filtering of SNPs in the non-coding regions based on GWAS and LD results.
    Cfilter             Filtering of coding genes affected by SNPs/Indels based on GWAS results
    Nfilter             Filtering of Indel windows in the non-coding regions based on GWAS results.


# preGWAS
### Get the tab-delimited format annotation file, and the annotation dictionaries
#### Usage:&emsp;python compGWAS.py preGWAS [-h] [-F [FLAG]] -g GBK -o OUT -s SEQ [-t ORG] -r REF -p PREFIX
#### {arguments}
    -h, --help                          Show this help message and exit.
    -F [FLAG], --FLAG [FLAG]            Ignore this argument or input string 'preGWAS'.
    -g GBK, --gbk GBK                   The gbk file for converting to tab-delimted format annotation.
    -t ORG, --org ORG                   The organism you are working on (default: -).
    -r REF, --ref REF                   The reference genome file.
    -o OUT, --out OUT                   The output file for tab-delimted format annotation file of all genes.
    -s SEQ, --seq SEQ                   The output protein sequences file from gbk file.
    -p PREFIX, --prefix PREFIX          The prefix of output file.
* The tab-delimited format annotation file contains 19 columns. The columns (**r1, r2, r3, r4**) are reserved for users to add other information they are interested in. The last 4 columns (**category0, category1, category2, category3**) are reserved for users to add higher-level functional annotations. In this package, the columns (**category1, category2, category3**) will be used to annotate functional categories of coding genes and the column **category1** will be used to annotate genes after `Filtering`. Here is an example of the annotation file with each line for a gene:

| # taxonomy | r1 | r2 | r3 | r4 | contig | start | end  | strand | symbol | geneID | mol_type | assemblyType               | productID | function                                       | category0                                                                                            | category1                                            | catetgory2 | category3 |
| ---------- | -- | -- | -- | -- | ------ | ----- | ---- | ------ | ------ | ------ | -------- | ------------------ | --------- | ---------------------------------------------- | --------------------------------------------------------------------------------------------------- | ---------------------------------------------------- | --------- | ------------ |
| 1314   | -  | -  | -  | -  | AP53   | 232   | 1587 | +      | dnaA   | AP53_1 | GENE     | complete | AP53_1         | Chromosomal replication initiator protein DnaA |    | DNA Metabolism // DNA replication // DNA-replication |           |                 | 
| 1314   | -  | -  | -  | -  | AP53   | 1742  | 2878 | +      | dnaN   | AP53_2 | GENE     | complete | AP53_2         | DNA polymerase III beta subunit (EC 2.7.7.7)   |    | DNA Metabolism // DNA replication // DNA-replication |           |                 |
| 1314   | -  | -  | -  | -  | AP53   | 2953  | 3150 | +      | -      | AP53_3 | GENE     | complete | AP53_3         | FIG01114317: hypothetical protein              |    | Hypothetical protein                                 |           |                 |


# GWAS
### GWA analysis of SNPs in the coding regions
#### Usage:&emsp;python compGWAS.py SNPgwas [-h] -S ALLSNP -c COLS COLS COLS -p PHENO1 PHENO1 -P PHENO0 PHENO0 [-i INFO INFO INFO] [-f REF REF] [-t [THREAD]] [-T [THRESHOLD]] [-o [OUTDIR]] [-O [PREFIX]] -R RSCRIPT -r RSCRIPT
#### {arguments}
    -h, --help                                  Show this help message and exit.
    -S ALLSNP, --allSNP ALLSNP                  The path of all strains' SNP annotation files, not include the files.
    -c COLS COLS COLS, --cols COLS COLS COLS    Input the column index of RefName (Contig), Loc and Allele of the SNPs annotation files (start with 0, e.g first column = 0).                         
    -p PHENO1 PHENO1, --pheno1 PHENO1 PHENO1    Input the file of phenotype1 strains' GCAs (only the number between GCA_ and .), and the name of the phenotype p.
    -P PHENO0 PHENO0, --pheno0 PHENO0 PHENO0    Input the file of phenotype0 strains' GCAs (only the number between GCA_ and .), and the name of the phenotype P. 
    -i INFO INFO INFO, --info INFO INFO INFO    Input the file of strains' covariates information, the covariates for logistic regression (seperated by ,), and the type of each covariate (f or n, f means factor and n means numeric, seperated by ,).
    -f REF REF, --ref REF REF                   Input the reference strain's GCA (only the number between GCA_ and .) and its phenotype (p or P).
    -t [THREAD], --thread [THREAD]              Input the thread used in 'CHI2distribution.R' and 'GWAS.R' (default: 1).
    -T [THRESHOLD], --threshold [THRESHOLD]     Input the Pvalue threshold of Allele used to filter logistic regression results (default: 0.05).
    -o [OUTDIR], --outdir [OUTDIR]              Output folder (default: .).
    -O [PREFIX], --prefix [PREFIX]              Output filename prefix of SNP results (default: SNP.output).
    -R RSCRIPT, --Rscript RSCRIPT               The path of Rscript, include 'Rscript'.
    -r RSCRIPT, --rscript RSCRIPT               The path of the compGWAS package, include 'compGWAS'.

### GWA analysis of the synergistic effect of SNPs/Indels in the coding regions
#### Usage:&emsp;python compGWAS.py CDSgwas [-h] -D ALLDIP -S ALLSNP -c COLS COLS COLS COLS COLS COLS -C CDSDIC [CDSDIC ...] -g GENOME -p PHENO1 PHENO1 -P PHENO0 PHENO0 [-i INFO INFO INFO] [-f REF REF] -l LENGTH LENGTH LENGTH [-t [THREAD]] [-T [THRESHOLD]] [-o [OUTDIR]] [-O [PREFIX]] -R RSCRIPT -r RSCRIPT
#### {arguments}
    -h, --help                                                                  Show this help message and exit.
    -D ALLDIP, --allDIP ALLDIP                                                  The path of all strains' InDel annotation files, not include the files.
    -S ALLSNP, --allSNP ALLSNP                                                  The path of all strains' SNP annotation files, not include the files.
    -c COLS COLS COLS COLS COLS COLS, --cols COLS COLS COLS COLS COLS COLS      Input the column index of RefName (Contig), Loc and Allele of the SNP annotation files and the InDel annotation files (start with 0, e.g first column = 0).            
    -C CDSDIC [CDSDIC ...], --CDSdic CDSDIC [CDSDIC ...]                        Input the CDS annotation dictionaries of the species to be analyzed.
    -g GENOME, --genome GENOME                                                  Input the FASTA format reference genome file.
    -f REF REF, --ref REF REF                                                   Input reference strain's GCA (only the number between GCA_ and .) and its phenotype (p or P).
    -p PHENO1 PHENO1, --pheno1 PHENO1 PHENO1                                    Input the file of phenotype1 strains' GCAs (only the number between GCA_ and .), and the name of the phenotype p.
    -P PHENO0 PHENO0, --pheno0 PHENO0 PHENO0                                    Input the file of phenotype0 strains' GCAs (only the number between GCA_ and .), and the name of the phenotype P.
    -i INFO INFO INFO, --info INFO INFO INFO                                    Input the file of strains' covariates information, the covariates for logistic regression (seperated by ,), and the type of each covariate (f or n, f means factor and n means numeric, seperated by ,).
    -l LENGTH LENGTH LENGTH, --length LENGTH LENGTH LENGTH                      Input the length extending outward, the length extending inward at the 5'end of CDS, and the length extending inward at the 3'end of CDS. The unit is bp.
    -t [THREAD], --thread [THREAD]                                              Input the thread used in 'CHI2distribution.R' and 'GWAS.R' (default: 1).
    -T [THRESHOLD], --threshold [THRESHOLD]                                     Input the Pvalue threshold of Allele used to filter logistic regression results (default: 0.05).
    -o [OUTDIR], --outdir [OUTDIR]                                              Output folder (default: .).
    -O [PREFIX], --prefix [PREFIX]                                              Output filename prefix of CDS results (default: CDS.output).
    -R RSCRIPT, --Rscript RSCRIPT                                               The path of Rscript, include 'Rscript'.
    -r RSCRIPT, --rscript RSCRIPT                                               The path of the compGWAS package, include 'compGWAS'.

### GWA analysis of SNPs/Indels windows in the non-coding regions
#### Usage:&emsp;python compGWAS.py nonCDSgwas [-h] -D ALLDIP -c COLS COLS COLS COLS -C CDSDIC [CDSDIC ...] -Z SIZE -p PHENO1 PHENO1 -P PHENO0 PHENO0 [-i INFO INFO INFO] [-f REF REF] -g GENOME [-t [THREAD]] [-T [THRESHOLD]] [-o [OUTDIR]] [-O [PREFIX]] -R RSCRIPT -r RSCRIPT
#### {arguments}
    -h, --help                                                  Show this help message and exit.
    -D ALLDIP, --allDIP ALLDIP                                  The path of all strains' InDel annotation files, not include the files.
    -c COLS COLS COLS COLS, --cols COLS COLS COLS COLS          Input the column index of RefName, Loc, Ref, and Muttype of InDels annotation file (start with 0, e.g first column = 0). 
    -C CDSDIC [CDSDIC ...], --CDSdic CDSDIC [CDSDIC ...]        Input the CDS annotation dictionaries of the species to be analyzed.
    -Z SIZE, --size SIZE                                        The size of the sliding window, the unit is bp.
    -p PHENO1 PHENO1, --pheno1 PHENO1 PHENO1                    Input the file of phenotype1 strains' GCAs (only the number between GCA_ and .), and the name of the phenotype p.
    -P PHENO0 PHENO0, --pheno0 PHENO0 PHENO0                    Input the file of phenotype0 strains' GCAs (only the number between GCA_ and .), and the name of the phenotype P.
    -i INFO INFO INFO, --info INFO INFO INFO                    Input the file of strains' covariates information, the covariates for logistic regression (seperated by ,), and the type of each covariate (f or n, f means factor and n means numeric, seperated by ,).
    -g GENOME, --genome GENOME                                  Input the FASTA format reference genome file.
    -f REF REF, --ref REF REF                                   Input reference strain's GCA (only the number between GCA_ and .) and its phenotype (p or P).    
    -t [THREAD], --thread [THREAD]                              Input the thread used in 'CHI2distribution.R' and 'GWAS.R' (default: 1).
    -T [THRESHOLD], --threshold [THRESHOLD]                     Input the Pvalue threshold of Allele used to filter logistic regression results (default: 0.05).
    -o [OUTDIR], --outdir [OUTDIR]                              Output folder (default: .).
    -O [PREFIX], --prefix [PREFIX]                              Output filename prefix of nonCDS results (default: nonCDS.output).
    -R RSCRIPT, --Rscript RSCRIPT                               The path of Rscript, include 'Rscript'.
    -r RSCRIPT, --rscript RSCRIPT                               The path of the compGWAS package, include 'compGWAS'.
* The directory of all strains' SNP annotation files and that of all strains' InDel annotation files need to be provided to the parameters `-S (--allSNP)` and `-D (--allDIP)` respectively. The naming format of these annotation files should be uniform. **The annotation file names should include the complete GCA number of the strain in the format 'GCA-' and 'GCA-' can only appear once in the file name. The naming format is `xxxGCA-12312412.1xxx`, of which xxx is the variable part. We recommend using `_` to connect different attributes of the strain, such as strain name, serotype, and GCA number. To avoid confusion, if there is a character `_` in the attribute, change it to `-`.**A typical example of the file name looks like: `bwa_AP53_EMM58.0_PHE-39584-W1_GCA-901571015.1_varSNP.anno`.
* The following table is the format of the SNP annotation file. The column `RefName` is the contig IDs in the reference genome. The columns `RefName`, `Loc` and `Allele` will be used in SNPgwas and CDSgwas. The `MutType` of SNPs includes 'inter-gene', 'inter-gene--5', 'inter-gene--3', 'inter-gene-5-', 'inter-gene-3-', 'inter-gene-3-5', 'inter-gene-5-3', 'inter-gene-5-5', 'inter-gene-3-3', 'pseudo-gene|ncRNA', 'CDS_nonSynon', and 'CDS_synon'.  

| #RefName | Loc    | Type | RefLen | Ref | AlleleNum | Allele | AlleleRatio | AlleleDepth | Depth | MutType         | Gene          | GeneID                 | Ori | AAchange  | AAsite |
| -------- | ------ | ---- | ------ | --- | --------- | ------ | ----------- | ----------- | ----- | --------------- | ------------- | ---------------------- | --- | --------- | ------ |
| AP53     | 40     | SNP  | 1      | A   | 1         | G      | 100.00      | 1           | 1     | inter-gene--5   | \|\|dnaA\|\|  | \|\|AP53_1\|\|         | 0   |           |        |
| AP53     | 462    | SNP  | 1      | T   | 1         | C      | 100.00      | 1           | 1     | CDS_synon       | dnaA          | AP53_1                 | +   |           | 77     |
| AP53     | 2354   | SNP  | 1      | G   | 1         | A      | 100.00      | 1           | 1     | CDS_nonSynon    | dnaN          | AP53_2                 | +   | Gly205Ser | 205    |
| AP53     | 14785  | SNP  | 1      | G   | 1         | T      | 100.00      | 1           | 1     | inter-gene-3-5  | ftsH\|\|-\|\| | AP53_14\|\|AP53_15\|\| | 0   |           |        |

* The following table is the format of the InDel annotation file. The column `RefName` is the contig IDs in the reference genome. The columns `RefName`, `Loc` and `Allele` will be used in CDSgwas. The columns `RefName`, `Loc`, `Ref`, and `Muttype` will be used in nonCDSgwas. The `MutType` of InDels includes 'inter-gene', 'inter-gene--5', 'inter-gene--3', 'inter-gene-5-', 'inter-gene-3-', 'inter-gene-3-5', 'inter-gene-5-3', 'inter-gene-5-5', 'inter-gene-3-3', 'pseudo-gene|ncRNA', 'frame-shift', 'in-frame', and 'in-frame-Stp'.

| #RefName | Loc    | Type | RefLen | Ref | AlleleNum | Allele | AlleleRatio | AlleleDepth | Depth | MutType         | Gene          | GeneID                    | Ori | AAchange  | AAsite |
| -------- | ------ | ---- | ------ | --- | --------- | ------ | ----------- | ----------- | ----- | --------------- | ------------- | ------------------------- | --- | --------- | ------ |
| AP53     | 3364   | DIP  | 1      | G   | 1         | -      | 100.00      | 1           | 1     | frame-shift     | -             | AP53_4                    | -   |           | 18     |
| AP53     | 29967  | DIP  | 3      | TTC | 1         | ---    | 100.00      | 1           | 1     | inter-gene-3-5  | -\|\|-\|\|    | AP53_rna34\|\|AP53_16\|\| | 0   |           |        |
| AP53     | 708569 | DIP  | 3      | --- | 1         | TTA    | 100.00      | 1           | 1     | in-frame-Stp    | mac;ideS      | AP53_713                  | -   |           | 166    |
| AP53     | 708680 | DIP  | 3      | --- | 1         | TTG    | 100.00      | 1           | 1     | in-frame        | mac;ideS      | AP53_713                  | -   |           | 129    |

* The file of strains' GCAs consists of one column (without column name). The number of lines is equal to the number of strains, and each line contains the GCA number (only the number between `GCA_` and `.`) of one strain. Here is an example of the file: 

        900996065
        901570175
        901570175
        901565655
        900993575

* When no covariate is considered, the parameter `-i (--info)` is not needed. In this case, logistic regression analysis is only performed on Allele (representing genotype). Otherwise, the users should provide three parameters to `-i (--info)`，which are the file of covariate information, the names of covariates (seperated by ,), and the types of covariates (f or n, f means factor and n means numeric, seperated by ,). **The first column of the covariate information file should be the GCA number (only the part before `.`) of each strain, and here `_` should be changed into `-`. The names of covariates are the headers in the covariate information file. The types of the covariates determine whether these covariates are used as character or numeric variables in logistic regression.** Here is an example of the covariate information file:

| Assembly_accession_number | Group | Country   | Year_group |
| ------------------------- | ----- | --------- | ---------- |
| GCA_900992365             | g1    | Australia | time2      |
| GCA_900984885             | g3    | Brazil    | time5      |
| GCA_901564955             | g8    | USA       | time5      |
| GCA_900988925             | g11   | Kenya     | time3      |

* The users should provide the parameter `-f (--ref)` only when the reference strain is included in the strains for GWA analysis. 
* Each GWAS will generate chi-square distribution results of Allele (representing genotype) and logistic regression results. We provide the parallel implimentation in these two analysis with the parameter `-t (--thread)`, and the default number of threads is 1. In order to avoid memory blow-up, we recommend that users choose an appropriate number of threads according to their computing resources.
* The users should provide the CDS annotation dictionaries generated by the program **preGWAS** to the parameter `-C (--CDSdic)`. We recommend using `CDS.dic` and `pseudo.dic` in CDSgwas, while `CDS.dic`, `pseudo.dic` and `RNA.dic` in nonCDSgwas.
* In CDSgwas, base substitutions are performed on each gene of the reference genome to determined the coding status of each gene in each strain based on the variant call results. The users need to provide three lengths (bp) to `-l (--length)`. The first length is used to determine the length of extension in the upstream and downstream of each coding gene in the reference genome. The second length is used to determine the length of extension inward at the 5'end of the coding gene. The third length is used to determine the length of extension inward at the 3'end of the coding gene. 
* In nonCDSgwas, the users need to provide the size (bp) of the sliding window. The length of sliding step is default to one tenth of the window size.


# GWASanno
### The annotation of coding genes affected by SNPs/Indels in the GWAS results.
#### Usage:&emsp;python compGWAS.py SNPCDSanno [-h] [-F [FLAG]] -t TYPE -l LOGISRE [-a ANNO] [-c COMB [COMB ...]] [-i INS [INS ...]] [-o [OUT]]
#### {arguments}
    -h, --help                                          Show this help message and exit.
    -F [FLAG], --FLAG [FLAG]                            Ignore this argument or input string 'SNPCDSanno'.
    -t TYPE, --type TYPE                                The type of annotation, SNP or CDS.
    -l LOGISRE, --logisre LOGISRE                       The results file of SNP or CDS logistic regression.
    -a ANNO, --anno ANNO                                When the annotation type is CDS, input the annotation file of all genes of the species to be analyzed.
    -c COMB [COMB ...], --comb COMB [COMB ...]          When annotation type is SNP, input the path of the files to be merged and the column to be combined, for example: 'allSNP/ all'.
    -i INS [INS ...], --ins INS [INS ...]               When annotation type is SNP, input the path of the files to be merged and the column to be combined, for example: 'allSNP/ all'.
    -o [OUT], --out [OUT]                               When annotation type is SNP, input the output annotation file of megered all SNPs (default: ./SNP.merged.allanno).
* `-t (--type)` and `-l (--logisre)` are compulsory parameters that the users should provide when performing the annotation of SNP or CDS GWAS results. 
* When the annotation type is **CDS**, the users need to provide to `-a (--anno)` the annotation file generated from `preGWAS`. When the annotation type is **SNP**, the users need to provide the correspoinding files to `-c (--comb)`, `-i (--ins)`, and `-o (--out)`. 
* The annotation output file of the SNPgwas results will be required in the program `LDprun` and `Filtering`.

### The annotation of SNPs/Indels windows in the non-coding regions in the GWAS results
#### Usage:&emsp;python compGWAS.py nonCDSanno [-h] -C CDSDIC [CDSDIC ...] -g GENOME -Z SIZE -l LOGISRE
#### {arguments}
    -h, --help                                                  Show this help message and exit.
    -C CDSDIC [CDSDIC ...], --CDSdic CDSDIC [CDSDIC ...]        Input the CDS annotation dictionaries of the species to be analyzed.
    -g GENOME, --genome GENOME                                  Input the FASTA format reference genome file.
    -Z SIZE, --size SIZE                                        The size of the sliding window.
    -l LOGISRE, --logisre LOGISRE                               The results file of nonCDS logistic regression.
* The CDS annotation dictionaries provided by the users for the parameter `-C (--CDSdic)` should be consistent with the dictionaries provided in the nonCDSgwas. 
* The annotation output file of the nonCDSgwas results will be required in the program `Filtering`.


# LDprun
### Linkage disequilibrium analysis of SNPs
#### Usage:&emsp;python compGWAS.py LDprun [-h] [-F [FLAG]] -S ALLSNP -p PHENO1 -P PHENO0 -g GENOME [-f REF REF] -l LOGISRE -LA LOGISANNO [-T [THRESHOLD]] -J JAVA -H HAPLOVIEW [-A [ARGUMENT]] [-T2 [THRESHOLD2]] [-d [MAXDIST]] -c2 CHI2 [-o [OUTDIR]] [-O [PREFIX]] 
#### {arguments}
    -h, --help                                          Show this help message and exit.
    -F [FLAG], --FLAG [FLAG]                            Ignore this argument or input string 'LDprun'.
    -S ALLSNP, --allSNP ALLSNP                          The path of all strains' SNP annotation files, not include the files.
    -p PHENO1, --pheno1 PHENO1                          Input the file of phenotype1 strains' GCAs (only the number between GCA_ and .).
    -P PHENO0, --pheno0 PHENO0                          Input the file of phenotype0 strains' GCAs (only the number between GCA_ and .).
    -g GENOME, --genome GENOME                          Input the FASTA format reference genome file.
    -f REF REF, --ref REF REF                           Input the reference strain's GCA (only the number between GCA_ and .) and its phenotype (p or P).
    -l LOGISRE, --logisre LOGISRE                       The results file of all SNPs logistic regression.
    -LA LOGISANNO, --LogisAnno LOGISANNO                Input the annotation file of SNP logistic regression results.
    -T [THRESHOLD], --threshold [THRESHOLD]             Input the Pvalue threshold of Allele used to filter SNP logistic regression results (default: 0.05).
    -J JAVA, --java JAVA                                The path of java, include 'java'.
    -H HAPLOVIEW, --Haploview HAPLOVIEW                 The path of Haploview.jar, include 'Haploview.jar'.
    -A [ARGUMENT], --Argument [ARGUMENT]                The parameters of Haploview.jar, excluding input and output parameters (default: '-n -dprime -blockoutput GAB -check -chromosome X').
    -T2 [THRESHOLD2], --threshold2 [THRESHOLD2]         The threshold of R^2 to separate the LD snps into blocks (default: 0.8).
    -d [MAXDIST], --maxdist [MAXDIST]                   The max distance of each block, the unit is kb (default: 10).
    -c2 CHI2, --Chi2 CHI2                               Input the SNP chi2 distribution file of the 2 phenotypes.
    -o [OUTDIR], --outdir [OUTDIR]                      Output folder (default: .).
    -O [PREFIX], --prefix [PREFIX]                      Output filename prefix of linkage disequilibrium analysis of SNPs (default: SNP.LD.output).
* The users need to provide the parameter `-f (--ref)` only when the reference strain is included in the strains for GWA analysis.
* For each contig, the program `LDprun` performs linkage disequilibrium analysis on the SNPs whose p-values meet the p-value threshold of `-T (--threshold)`. The output file contains the non-synonymous SNPs that do not belong to any block and the non-synonymous SNPs with the lowest p-value within each block. Those SNPs' the second lowest allele frequency should be greater than 2. 


# Filtering
### Filtering of SNPs in the coding regions based on GWAS and LD results 
#### Usage:&emsp;python compGWAS.py SCfilter [-h] -s SCREENOUT [-O [PREFIX]] -a ANNO [-T [THRESHOLD]]
#### {arguments}
    -h, --help                                  Show this help message and exit.
    -s SCREENOUT, --screenout SCREENOUT         The path of SNPblockID.screenout file(s) of linkage disequilibrium analysis of SNPs.
    -O [PREFIX], --prefix [PREFIX]              Output filename prefix of linkage disequilibrium analysis of SNPs (default: SNP.LD.output).
    -a ANNO, --anno ANNO                        The annotation file of all genes of the species to be analyzed.
    -T [THRESHOLD], --threshold [THRESHOLD]     Input the Pvalue significance threshold of Allele used to filter logistic regression results (default: 0.05).
* The parameter `-O (--prefix)` should be consistent with that in `LDprun`.
* For SNPs in the coding regions, the command `SCfilter` merges all files generated by `LDprun`, removes all SNPs that do not meet the p-value threshold, and only keeps the SNPs with the lowest p-value within each gene.

### Filtering of SNPs in the non-coding regions based on GWAS and LD results
#### Usage:&emsp;python compGWAS.py SNfilter [-h] -C CDSDIC [CDSDIC ...] -LA LOGISANNO -a ANNO -c2 CHI2 [-T [THRESHOLD]] [-d DISTANCE DISTANCE]
#### {arguments}
    -h, --help                                                  Show this help message and exit.
    -C CDSDIC [CDSDIC ...], --CDSdic CDSDIC [CDSDIC ...]        Input the CDS annotation dictionaries of the species to be analyzed.
    -LA LOGISANNO, --LogisAnno LOGISANNO                        Input the annotation file of SNP logistic regression results.
    -a ANNO, --anno ANNO                                        The annotation file of all genes of the species to be analyzed.
    -c2 CHI2, --Chi2 CHI2                                       Input the SNP chi2 distribution file of the 2 phenotypes.
    -T [THRESHOLD], --threshold [THRESHOLD]                     Input the Pvalue significance threshold of Allele used to filter logistic regression results (default: 0.05).
    -d DISTANCE DISTANCE, --distance DISTANCE DISTANCE          Input the threshold of 5' distance (bp) and 3' distance (bp).
* The CDS annotation dictionaries provided by the users to `-C (--CDSdic)` should be consistent with the dictionaries provided in the nonCDSgwas. 
* For SNPs in the non-coding regions, the command `SNfilter` filters the SNPs that do not meet the p-value threshold or whose second lowest allele frequency <= 2. It keeps the SNP with the lowest p-value in each non-coding region, and annotates its upstream and downstream neighboring genes that meet the distance threshold.  

### Filtering of coding genes affected by SNPs/Indels based on GWAS results
#### Usage:&emsp;python compGWAS.py Cfilter [-h] -l LOGISRE -a ANNO -c2 CHI2 [-T [THRESHOLD]]
#### {arguments}
    -h, --help                                  Show this help message and exit.
    -l LOGISRE, --logisre LOGISRE               The results file of CDS logistic regression.
    -a ANNO, --anno ANNO                        The annotation file of all genes of the species to be analyzed.
    -c2 CHI2, --Chi2 CHI2                       Input the CDS chi2 distribution file of the 2 phenotypes.
    -T [THRESHOLD], --threshold [THRESHOLD]     Input the Pvalue significance threshold of Allele used to filter logisticRegre results (default: 0.05).
* The command `Cfilter` performs filtering on CDS logistic regression results and it removes genes that do not meet the p-value threshold or whose second lowest allele frequency <= 2. 

### Filtering of windows in the non-coding regions based on GWAS results
#### Usage:&emsp;python compGWAS.py Nfilter [-h] -C CDSDIC [CDSDIC ...] -LA LOGISANNO -a ANNO -c2 CHI2 [-T [THRESHOLD]] [-d DISTANCE DISTANCE]
#### {arguments}
    -h, --help                                                  Show this help message and exit
    -C CDSDIC [CDSDIC ...], --CDSdic CDSDIC [CDSDIC ...]        Input the CDS annotation dictionaries of the species to be analyzed.
    -LA LOGISANNO, --LogisAnno LOGISANNO                        Input the annotation file of SNP logistic regression results.
    -a ANNO, --anno ANNO                                        The annotation file of all genes of the species to be analyzed.
    -c2 CHI2, --Chi2 CHI2                                       Input the nonCDS chi2 distribution file of the 2 phenotypes.
    -T [THRESHOLD], --threshold [THRESHOLD]                     Input the Pvalue significance threshold of Allele used to filter logisticRegre results (default: 0.05).
    -d DISTANCE DISTANCE, --distance DISTANCE DISTANCE          Input the threshold of 5' distance(bp) and 3' distance(bp).
* The CDS annotation dictionaries provided by the users to `-C (--CDSdic)` should be consistent with the dictionaries provided in the nonCDSgwas.
* The command `Nfilter` filters the windows that do not meet the p-value threshold or whose second lowest allele frequency <= 2. It keeps the windows with the lowest p-value in each non-coding region, and annotates its upstream and downstream neighboring genes that meet the distance threshold.

# Citation
Wang PY, Liang Z, Chen ZS, Bao YJ*, Castellino F*. compGWAS: a new GWAS tool allows revelation of the genetic architecture and risk stratification for the versatile pathogen Streptococcus pyogenes. Under review 2025. https://www.researchsquare.com/article/rs-5049272/v1.
