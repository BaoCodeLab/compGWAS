## The data in this folder is an exmple output of reference processing using preGWAS.

### Prepare the gbk file of your reference as input:
```
AP53_NCBI.gbk
AP53_genome.fasta
```

### Run the preGWAS to process the gbk file:
```
python compGWAS.py preGWAS -g AP53_NCBI.gbk -o AP53_NCBI.gbk2gpff.tab -s AP53.NCBI.gbk.prot.faa [-t AP53] -r AP53_genome.fasta -p AP53.anno2
```

### A series of output files will be generated and will be used for mutation annotation later:
```
AP53_NCBI.gbk2gpff.tab
AP53.anno2.CDS.dic
AP53.anno2.CDSseq.dic
AP53.anno2.contg.dic
AP53.anno2.contg.dic
AP53.anno2.gene.dic
AP53.anno2.mol.dic
AP53.anno2.ncRNA.dic
AP53.anno2.ncrnaSeq.dic
AP53.anno2.pseudo.dic
AP53.anno2.pseudoSeq.dic
AP53.NCBI.gbk.prot.faa
```
