## Mos99 NA Deep mutation scanning

### Dependencies ###
* [python](https://www.python.org/) (version 3.9)
* [snakemake](https://snakemake.readthedocs.io/en/stable/)
* [flash](https://github.com/dstreett/FLASH2)
* [seqtk](https://github.com/lh3/seqtk)
* [cutadapt](https://cutadapt.readthedocs.io/en/stable/)
* [pandas](https://pandas.pydata.org/)
* [biopython](https://github.com/biopython/biopython)
* [DSSP](https://ssbio.readthedocs.io/en/latest/instructions/dssp.html)

### Input files ###
* All raw reads in fastq format should be placed in fastq/
* 

### Installation ###
Install dependencies by conda:   
```
conda create -n NA -c bioconda -c anaconda -c conda-forge -c salilab \
  python=3.9 \
  seqtk \
  flash \
  biopython \
  cutadapt \
  snakemake \
  dssp
```   

### Data analysis ###
1. Activate conda environment:   
``source activate NA``

2. Using UMI to correct sequencing errors:   
``python script/Dedup_UMI.py fastq NNNNNNN 0.8 2``
