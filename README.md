# 🧬 EM-based STR Structure Prediction

This project introduces a novel statistical model for Short Tandem Repeat (STR) structure prediction and genotyping based on the Expectation-Maximization (EM) algorithm.

## Directory Structure

* **Root Directory**
    * `src/`
        * `main.py`
        * `em.py`
        * `fun.py`
    * `data/`
    * `output/`
    * `environment.txt`
    * `README.md`

## ⚙️ Setup Instructions

### 1. Install Micromamba

Follow the [official Micromamba installation instructions](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) for your operating system.

### 2. Download the Reference Genome  
Before running the analysis, you must download the **HG002_NA24385_son/
NIST_HiSeq_HG002** reference genome. You can obtain it from trusted sources such as:
- **UCSC Genome Browser** ([link](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/))
- **Ensembl** ([link](https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/))
- **Other genomic databases**

Once downloaded, place the file in the `data/` directory.

### 3. Create the Environment  
Once Micromamba is installed, navigate to the project's root directory in your terminal and run the following commands to create and activate the environment:

```bash
# Create the environment (you can choose a different name instead of 'str-em-env')
micromamba create --name str-em-env --file environment.txt

# Activate the environment
micromamba activate str-em-env
```

### 4. Run the Project

With the environment activated, navigate to the src directory and run the main script:

```bash
cd src
python main.py
```
This will process the input data (as configured in main.py) and generate output in the output/ directory.

## 📊 Output Interpretation

Upon successful execution, the `output/` directory will contain two subfolders:

- **`output/frequencies/`**: Contains files representing the probability distribution across all considered alleles for each STR locus analyzed.
- **`output/alleles/`**: Contains files listing the two most probable alleles identified by the EM model for each STR locus.
