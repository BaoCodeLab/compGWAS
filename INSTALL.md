# compGWAS Installation Guide

A quick and easy installation guide for **compGWAS** - A Python toolkit for comprehensive GWAS analysis of SNPs/Indels.

## Prerequisites

- **Python** > 3.7+
- **R** > 3.5+
- **Java** = 1.8 (for Haploview/LD analysis)
- **numpy** > 1.16.0 (Python package manager)
- **pandas** >=0.24.0
- **R packages**: "foreach", "doParallel", "BaylorEdPsych"

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
wget https://www.broadinstitute.org/files/shared/diabetes/haploview/Haploview.jar
cd ..
```

## Test compGWAS

```
### Test Individual Modules

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
- **Cfilter** - Filter genes affected by variants
- **Nfilter** - Filter windows in non-coding regions

### Example Usage

```bash
# Prepare genomic data
compGWAS preGWAS -g input.gbk -o output.txt -s sequences.fasta -r genome.fasta -p PREFIX

# Run SNP GWAS analysis
compGWAS SNPgwas -S /path/to/snps/ -c 0 1 2 \
  -p pheno1.txt p1 -P pheno0.txt p0 \
  -R /usr/bin/Rscript -r /path/to/compGWAS/ -t 4
```

See [README.md](README.md) for detailed usage documentation.

---

## Troubleshooting

### Issue: `command not found: compGWAS`

**Solution:** Ensure the installation completed successfully:
```bash
pip3 install -e . --force-reinstall
```

Or use the full path:
```bash
python3 /path/to/compGWAS/compGWAS.py -h
```

### Issue: `ModuleNotFoundError: No module named 'pandas'`

**Solution:**
```bash
pip3 install --upgrade pandas numpy
```

### Issue: `R package not found`

**Solution:** Install R packages with proper repository:
```bash
Rscript -e "install.packages(c('foreach', 'doParallel', 'BaylorEdPsych'), repos='https://cloud.r-project.org')"
```

### Issue: `Rscript: command not found`

**Solution:** Find the location of Rscript:
```bash
which Rscript
```

Then use the full path when running compGWAS:
```bash
compGWAS SNPgwas -R /usr/bin/Rscript [other arguments]
```

### Issue: `java: command not found`

**Solution:** Install Java and verify:
```bash
# Ubuntu/Debian
sudo apt-get install default-jdk

# Check installation
java -version
```

### Issue: Permission denied

**Solution:**
```bash
chmod +x compGWAS.py
```

Or run with explicit Python:
```bash
python3 -m compGWAS [command] [arguments]
```

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

---

## Getting Help

- **View help:** `compGWAS -h`
- **Module-specific help:** `compGWAS <command> -h`
- **GitHub Issues:** https://github.com/BaoCodeLab/compGWAS/issues
- **README:** See [README.md](README.md) for detailed documentation

---


## License

compGWAS is released under the GNU General Public License v3.0. See [LICENSE](LICENSE) for details.

---

**Last Updated:** June 2024  
**compGWAS Version:** 1.0.0+
