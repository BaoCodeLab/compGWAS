# compGWAS Installation Guide

A quick and easy installation guide for **compGWAS** - A Python toolkit for comprehensive GWAS analysis of SNPs/Indels.

## Prerequisites

- **Python** > 3.7+
- **numpy** > 1.16.0
- **pandas** >=0.24.0
- **R** > 3.5+
- **R packages**: "foreach", "doParallel", "BaylorEdPsych"
- **Java** = 1.8
- **Haploview** (for Haploview/LD analysis)

## Quick Start (5 Steps)

### Step 1: Clone the Repository

```bash
git clone https://github.com/BaoCodeLab/compGWAS.git
cd compGWAS
```

### Step 2: Install System Dependencies

Before installing compGWAS, ensure you have the required system packages installed:

**Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install python3 python3-pip r-base default-jdk
```

**CentOS/RHEL:**
```bash
sudo yum install python3 python3-pip R java-1.8.0-openjdk
```

**macOS (with Homebrew):**
```bash
brew install python3 r openjdk
```

### Step 3: Install compGWAS

**Standard installation:**
```bash
pip3 install -e .
```

**Or from source directory:**
```bash
python3 setup.py install
```

**Verify installation:**
```bash
compGWAS -h
```

### Step 4. Install R Packages

Inside R console:
```r
install.packages("foreach")
install.packages("doParallel")
install.packages("BaylorEdPsych")
q()
```

### Step 5. Download Haploview (Optional but Recommended)

Haploview is needed for linkage disequilibrium analysis:

```bash
mkdir -p ~/compgwas_tools
cd ~/compgwas_tools
wget https://www.broadinstitute.org/ftp/pub/mpg/haploview/Haploview.jar
cd ..
```

## Test compGWAS

```
# List all available commands
compGWAS -h

# Test preGWAS module
compGWAS preGWAS -h

# Test SNPgwas module
compGWAS SNPgwas -h
```

## Usage

After successful installation, run compGWAS with:

```bash
compGWAS <command> [options]
```

### Available Commands

- **preGWAS** - Prepare reference genome and annotation files
- **SNPgwas** - GWAS analysis of SNPs in coding regions
- **CDSgwas** - GWAS analysis of SNPs/Indels synergistic effects
- **nonCDSgwas** - GWAS analysis in non-coding regions
- **SNPCDSanno** - Annotate coding genes affected by variants
- **nonCDSanno** - Annotate non-coding region variants
- **LDprun** - Linkage disequilibrium pruning
- **SCfilter** - Filter SNPs in coding regions
- **SNfilter** - Filter SNPs in non-coding regions
- **Cfilter** - Filter genes truncated by variants
- **Nfilter** - Filter windows in non-coding regions

### Example Usage

```bash
# Prepare reference genome and annotation files
compGWAS preGWAS -g input.gbk -o output_tab-delimited.txt -s reference_proteins.faa -r genome.fasta -p PREFIX

# Run SNP GWAS analysis
compGWAS SNPgwas -S /path/to/snps/ -c 0 1 6 \
  -p pheno1_IDs.txt pheno1 -P pheno0_IDs.txt pheno0 \
  -f reference_ID pheno1 -t 4 -T 0.001 -o output_dir -O PREFIX
  -R /usr/bin/Rscript -r /path/to/compGWAS/
```

See [README.md](README.md) for detailed usage documentation.

---

## Uninstallation

To remove compGWAS:

```bash
pip3 uninstall compGWAS
```

To remove the cloned repository:

```bash
rm -rf ~/compGWAS
```

## License

compGWAS is released under the GNU General Public License v3.0. See [LICENSE](LICENSE) for details.

---

**Last Updated:** June 2026  
**compGWAS Version:** 1.0.0+
