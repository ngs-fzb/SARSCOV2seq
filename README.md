# SARSCOV2seq
Pipeline for the analysis of SARS-CoV-2 tiled amplicon sequencing data

# Installation
```
git clone https://github.com/ngs-fzb/SARSCOV2seq
```
# Requirements

## Conda
Install [Conda](https://conda.io/docs/) or [Miniconda](https://conda.io/miniconda.html)
```
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```
### Pagolin
```
git clone https://github.com/cov-lineages/pangolin
cd pangolin
conda env create -f environment.yml
conda activate pangolin
pip install . 
conda deactivate
```

### MTBseq
```
conda create --name mtbseq
conda activate mtbseq
conda install -c bioconda mtbseq
conda deactivate
```

For more information on MTBseq click here: https://github.com/ngs-fzb/MTBseq_source

Copy reference fasta and annotation to the mtbseq /ref/ directory!

### iVar
```
conda create --name ivar
conda activate ivar
conda install -c bioconda ivar
conda deactivate
```

### Open the sarscov2seq file in a text editor and set paths to match the actual location of the "Ref" and "Scripts" directory on your computer

# Usage
Fastq files need to be in the format required by MTBseq please read the [MANUAL.md](https://github.com/ngs-fzb/MTBseq_source/blob/master/MANUAL.md).

Change to the directory containing your fastq files, copy sarscov2seq to the directory and change parameters as needed afterwards execute:
```
sarscov2seq
```

