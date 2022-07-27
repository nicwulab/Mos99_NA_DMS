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

### Dependencies installation ###
1. Install dependencies by conda:   
```
conda create -n NA -c bioconda -c anaconda -c conda-forge \
  python=3.9 \
  seqtk \
  flash \
  biopython \
  cutadapt \
  snakemake \
  prody
```   

2. Activate conda environment:   
``source activate NA``

### Calculating mutational fitness from sequencing data ###
1. Using UMI to correct sequencing errors:   
``python3 script/Dedup_UMI.py fastq NNNNNNN 0.8 2``

2. Counting mutations:   
``snakemake -s Mos99_pipeline.smk -j 10``

3. Convert counts to fitness:   
``python3 script/count2fitness.py``
    - Output file: [./result/Mos99_fit.csv](./result/Mos99_fit.csv)

### Data analysis ###
1. Compute mutational tolerance for each residue
``python3 script/Mean_mut_fit_per_resi.py``
    - Input file: [./result/Mos99_fit.csv](./result/Mos99_fit.csv)
    - Output file: [./result/Mos99_mean_mut_fit.tsv](./result/Mos99_mean_mut_fit.tsv)

### Plotting ###
1. Plots for checking data quality
``Rscript script/plot_QC.R``
    - Input file: [./result/Mos99_fit.csv](./result/Mos99_fit.csv)
    - Output file:
      - [./graph/Mos99_sil_non_mis.png](./graph/Mos99_sil_non_mis.png)
      - [./graph/Mos99_rep_cor.png](./graph/Mos99_rep_cor.png)

2. Comparing the data in this study with our previous study ([Wang et al. 2021](https://elifesciences.org/articles/72516))
``Rscript script/plot_cross_valid.R``
    - Input files:
      - [./result/Mos99_fit.csv](./result/Mos99_fit.csv)
      - [./data/Mos99_single_mut_Wang_et_al.tsv](./data/Mos99_single_mut_Wang_et_al.tsv)
    - Output file: [./graph/DMS_cross_validate.png](./graph/DMS_cross_validate.png)

3. Heatmap of mutational fitness
``Rscript script/plot_fitness_heatmap.R``
    - Input file: [./result/Mos99_fit.csv](./result/Mos99_fit.csv)
    - Output file: [./graph/Mos99_fit_heatmap.png](./graph/Mos99_fit_heatmap.png)
