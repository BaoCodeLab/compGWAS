## This folder contains demo files used for running the packages.
* [SNPcall_INV-SSTI](#SNPcall_INV-SSTI)
  
        This folder contains example SNP call results for starting from preAnno. Users should prepare their own mutation call files.
  
* [InDelcall_INV-SSTI](#InDelcall_INV-SSTI)
  
        This folder contains example InDel call results for starting from preAnno. Users should prepare their own mutation call files.
        
* [INV_GCAs](#INV_GCAs)
  
        The file contains the sample IDs for the phenotype INV.
  
* [SSTI_GCAs](#SSTI_GCAs)
  
        The file contains the sample IDs for the phenotype SSTI.
  
* The two ID files are required for running pair-wise GWAS analysis. The IDs should be consistent with the IDs embedded in the file names in the folder [SNPcall_INV-SSTI](#SNPcall_INV-SSTI) or [InDelcall_INV-SSTI](#InDelcall_INV-SSTI).

* [sample_covariates_information](#sample_covariates_information)
  
        This file contains the covariate information with each column representing one covariate. If the users hope to consider confounding effects of covariates, the related information should be prepared by users and provide it when running GWAS analysis.
  
