# compGWAS Installation Guide

A quick and easy installation guide for **compGWAS** - A Python toolkit for comprehensive GWAS analysis of SNPs/Indels.

## Quick Start (3 Steps)

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

---

## Detailed Installation Steps

### Prerequisites

compGWAS requires:

- **Python** 3.7+
- **R** 3.5+
- **Java** (for Haploview/LD analysis)
- **pip** (Python package manager)

### 1. Install System Dependencies

#### Ubuntu/Debian
```bash
sudo apt-get update
sudo apt-get install -y \
    python3 \
    python3-pip \
    python3-dev \
    r-base \
    r-base-dev \
    default-jdk \
    git
```

#### CentOS/RHEL
```bash
sudo yum install -y \
    python3 \
    python3-pip \
    python3-devel \
    R \
    R-devel \
    java-1.8.0-openjdk \
    git
```

#### macOS
```bash
# Using Homebrew (install from https://brew.sh if needed)
brew install python3 r openjdk git
```

#### Verify installations:
```bash
python3 --version
R --version
java -version
```

### 2. Clone compGWAS Repository

```bash
git clone https://github.com/BaoCodeLab/compGWAS.git
cd compGWAS
```

### 3. Install Python Package Dependencies

The setup.py automatically handles Python dependencies (numpy, pandas):

```bash
pip3 install -e .
```

Or with upgrade flag:
```bash
pip3 install -e . --upgrade
```

### 4. Install R Packages

**Interactive installation:**
```bash
R
```

Inside R console:
```r
install.packages("foreach")
install.packages("doParallel")
install.packages("BaylorEdPsych")
q()
```

**Or non-interactive:**
```bash
Rscript -e "install.packages(c('foreach', 'doParallel', 'BaylorEdPsych'), repos='http://cran.r-project.org')"
```

### 5. Download Haploview (Optional but Recommended)

Haploview is needed for linkage disequilibrium analysis:

```bash
mkdir -p ~/compgwas_tools
cd ~/compgwas_tools
wget https://www.broadinstitute.org/files/shared/diabetes/haploview/Haploview.jar
cd -
```

---

## Installation Verification

### Test All Components

```bash
# Test Python packages
python3 -c "import numpy, pandas; print('✓ Python packages OK')"

# Test R packages
Rscript -e "library(foreach); library(doParallel); library(BaylorEdPsych); print('R packages OK')"

# Test Java
java -version

# Test compGWAS
compGWAS -h
```

### Test Individual Modules

```bash
# List all available commands
compGWAS -h

# Test preGWAS module
compGWAS preGWAS -h

# Test SNPgwas module
compGWAS SNPgwas -h
```

---

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

## Advanced Installation Options

### Using a Virtual Environment (Recommended)

**Create and activate a virtual environment:**

```bash
# Python venv
python3 -m venv compgwas_env
source compgwas_env/bin/activate

# Or with conda
conda create -n compgwas python=3.9
conda activate compgwas
```

**Then install compGWAS:**
```bash
pip3 install -e .
```

### Using Conda (Recommended for HPC)

```bash
# Create conda environment with all dependencies
conda create -n compgwas python=3.9 r-base numpy pandas

# Activate environment
conda activate compgwas

# Install R packages
R --slave -e "install.packages(c('foreach', 'doParallel', 'BaylorEdPsych'))"

# Install compGWAS
pip3 install -e .
```

### Docker Installation

**Build Docker image:**
```bash
docker build -t compgwas:latest .
```

**Run compGWAS in Docker:**
```bash
docker run --rm -v /path/to/data:/data compgwas:latest compGWAS -h
```

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

## Installation on HPC Clusters

### Using SLURM Module System

**Create a batch script `run_compgwas.sh`:**

```bash
#!/bin/bash
#SBATCH --job-name=compgwas
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=16G

# Load modules
module load python/3.9
module load r/4.0
module load java/11

# Activate virtual environment (if using one)
source /path/to/compgwas_env/bin/activate

# Run compGWAS
compGWAS SNPgwas -S /data/snps/ -t 4 [other arguments]
```

**Submit job:**
```bash
sbatch run_compgwas.sh
```

### Using Conda on HPC

```bash
# Install miniconda if not available
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# Create environment
conda create -n compgwas python=3.9 r-base numpy pandas -c conda-forge

# Activate and install
conda activate compgwas
pip3 install -e /path/to/compGWAS
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

## Citation

If you use compGWAS in your research, please cite:

> Wang PY, Liang Z, Chen ZS, Bao YJ*, Castellino F*. compGWAS: a new GWAS tool allows revelation of the genetic architecture and risk stratification for the versatile pathogen Streptococcus pyogenes.

---

## License

compGWAS is released under the GNU General Public License v3.0. See [LICENSE](LICENSE) for details.

---

**Last Updated:** June 2024  
**compGWAS Version:** 1.0.0+
