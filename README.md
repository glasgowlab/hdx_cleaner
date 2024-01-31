# HDX cleaner and parser

## Main purpose

Clean/parse the raw output from **biotools** and **HDExaminer** to do peptide level or high resolution (single residue to mini-peptides) HDX/MS data anlysis .


## Typical workflow

### 1. PIGEON

To exclude the ambiguous peptides

```
python biotools-cleaner.py --t *MSMMS.csv  --o peptide_pool.csv --ov drop --r rangelist_drop.csv 
```

where:

* --t: msmms.csv from biotools analysis
* --o: peptide pool from HDExaminer
* --r: rangelist of the good peptide,


### 2. HDExaminer maunual QC

To check the MS spectra, exclude the bad peptides


### 3. POST_PIGEON 

To parse the output from HDExaminer using interactive object-oriented notebook 

Main function:
* exclude the bad peptides
* add back exchange rate
* peptide subtraction
* make uptake plot/pymol plot/heatmap (peptide level analysis)
* prepare the input bayesian hdx iso (high resolution protection fators)
* analysis of the PFs estaiomation resultes


## Installation

```
# clone the PIGEON repo
git clone git@github.com:glasgowlab/PIGEON.git

# clone the bayesian hdx iso repo
git clone git@github.com:lucl13/bayesian_hdx.git
```
### Dependencies

```
# install the dependencies
conda install numpy pandas numba
conda install conda-forge::seaborn
conda install conda-forge::mdanalysis
pip install --index-url https://pypi.cs.uni-tuebingen.de/simple/ pyopenms

```

## Docker
```
docker build -f docker/Dockerfile -t pigeon:latest .
docker run -p 8888:8888 -v $(pwd):/home/jovyan/work --rm pigeon
```

click the link in the terminal to open the jupyter notebook
cd to example/docker_example
click and run


> **NOTE:** bayesian_hdx at the above link is a fork of the original repo at https://github.com/salilab/bayesian_hdx with some modifications to the code to make it compatible with the PIGEON workflow and support for the isotopic envelope fitting.


