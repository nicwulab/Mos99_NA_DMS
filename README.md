## Mos99 NA Deep mutation scanning

### Dependencies ###
* [python](https://www.python.org/) (version 3.9)
* [snakemake](https://snakemake.readthedocs.io/en/stable/)
* [flash](https://github.com/dstreett/FLASH2)
* [seqtk](https://github.com/lh3/seqtk)
* [cutadapt](https://cutadapt.readthedocs.io/en/stable/)
* [pandas](https://pandas.pydata.org/)
* [biopython](https://github.com/biopython/biopython)

### Input files ###
* All raw reads in fastq format should be placed in fastq/
* 

### Installation ###
1. Install dependencies by conda:   
```
conda create -n NA -c bioconda -c anaconda -c conda-forge \
  python=3.9 \
  seqtk \
  flash \
  biopython \
  cutadapt \
  snakemake
```   

2. Activate conda environment:   
``source activate NA``

### Calling mutations from sequencing data ###
1. Using UMI to correct sequencing errors:   
``python script/Dedup_UMI.py fastq NNNNNNN 0.8 2``

2. Counting mutations:   
``snakemake -s Mos99_pipeline.smk -j 10``

3. Convert counts to fitness:   
``python3 script/count2fitness.py``
    - Output file: [./result/Mos99_fit.csv](./result/Mos99_fit.csv)
