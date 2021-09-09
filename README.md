# Mos99 NA Deep mutation scanning

## Dependencies ##
* python=3.9
* [snakemake](https://snakemake.readthedocs.io/en/stable/)
* [flash](https://github.com/dstreett/FLASH2)
* [seqtk](https://github.com/lh3/seqtk)
* [cutadapt](https://cutadapt.readthedocs.io/en/stable/)
* [pandas](https://pandas.pydata.org/)
* [biopython](https://github.com/biopython/biopython)
* [DSSP](https://ssbio.readthedocs.io/en/latest/instructions/dssp.html)
## Installation ##
Install dependencies by conda:

```conda create -n NA -c bioconda -c anaconda python=3.9 seqtk flash biopython cutadapt snakemake```

Install DSSP by brew (Mac OS X):

```
brew install dssp
```

Before running the analysis, do:

```
source activate NA(for Mac)
```
## Steps ##
