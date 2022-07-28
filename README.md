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
* [./Fasta/Mos99_NA.pep](./Fasta/Mos99_NA.pep): Mos99 NA protein sequence.
* [./Fasta/Human_H3N2_NA_2020.aln](./Fasta/Human_H3N2_NA_2020.aln): Full-length NA protein sequences from human H3N2 downloaded from [GISAID](https://www.gisaid.org/).
* [./data/ASA.table](./data/ASA.table): Amino acid solvent accessibility (ASA) from [Tien et al. 2013](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0080635).
* [./data/sites_info.tsv](./data/sites_info.tsv): Antigenic regions and acive site residues are defined by [Colman et al. 1983](https://www.nature.com/articles/303041a0) and [McAuley et al. 2019](https://www.frontiersin.org/articles/10.3389/fmicb.2019.00039/full), respectively.
* [./data/foldx_msa_transformer.csv](./data/foldx_msa_transformer.csv): Stability effect was predicted using [FoldX](https://academic.oup.com/bioinformatics/article/35/20/4168/5381539) and natural fitness was inferred using [MSA Transformer](https://www.biorxiv.org/content/10.1101/2021.02.12.430858v1).
* [./PDB/Mos99_WT_NA_monomer.pdb](./PDB/Mos99_WT_NA_monomer.pdb)
* [./PDB/Mos99_WT_NA_tetramer.pdb](./PDB/Mos99_WT_NA_tetramer.pdb)
* [./PDB/Mos99_WT_sialic_acid_final.pdb](./PDB/Mos99_WT_sialic_acid_final.pdb)

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

2. Assign residue type and calculate RSA   
``python3 script/pos_type_analysis.py``
    - Input files:
      - [./data/ASA.table](./data/ASA.table)
      - [./result/Mos99_mean_mut_fit.tsv](./result/Mos99_mean_mut_fit.tsv)
      - [./data/sites_info.tsv](./data/sites_info.tsv)
      - [./PDB/Mos99_WT_NA_monomer.pdb](./PDB/Mos99_WT_NA_monomer.pdb)
      - [./PDB/Mos99_WT_NA_tetramer.pdb](./PDB/Mos99_WT_NA_tetramer.pdb)
    - Output file: [./result/position_type_vs_fit.tsv](./result/position_type_vs_fit.tsv)

3. Calculate distance to active site for each residue   
``python3 script/Dist_analysis.py``
    - Input file: [./PDB/Mos99_WT_sialic_acid_final.pdb](./PDB/Mos99_WT_sialic_acid_final.pdb)
    - Output file: [./result/Dist_to_active_site.tsv](./result/Dist_to_active_site.tsv)

4. Calculate natural mutation frequency   
``python3 script/natural_mut_analysis.py``
    - Input files:
      - [./Fasta/Mos99_NA.pep](./Fasta/Mos99_NA.pep)
      - [./Fasta/Human_H3N2_NA_2020.aln](./Fasta/Human_H3N2_NA_2020.aln)
      - [./result/Mos99_fit.csv](./result/Mos99_fit.csv)
    - Output file:
      - [./result/N2_mutation_freq.tsv](./result/N2_mutation_freq.tsv)

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

4. Compare RSA and fit across residue types   
``Rscript script/plot_pos_type_analysis.R``
    - Input file: [./result/position_type_vs_fit.tsv](./result/position_type_vs_fit.tsv)
    - Output files:
      - [./graph/fit_vs_RSA.png](./graph/fit_vs_RSA.png)
      - [./graph/position_type_vs_RSA.png](./graph/position_type_vs_RSA.png)
      - [./graph/position_type_vs_fit.png](./graph/position_type_vs_fit.png)

5. Plot correlation between fitness and distance to active site   
``Rscript script/plot_dist_to_active_site.R``
    - Input files:
      - [./result/position_type_vs_fit.tsv](./result/position_type_vs_fit.tsv)
      - [./result/Dist_to_active_site.tsv](./result/Dist_to_active_site.tsv)
    - Output file: [./graph/fit_vs_dist.png](./graph/fit_vs_dist.png)

6. Plot correlation between fitness and natural mutation frequency   
``Rscript script/plot_natural_mut_fit.R``
    - Input file: [./result/N2_mutation_freq.tsv](./result/N2_mutation_freq.tsv)
    - Ouput file: [./graph/natural_mut_freq_fit.png](./graph/natural_mut_freq_fit.png)

7. Plot DMS fitness vs predicted stability effect using [FoldX](https://academic.oup.com/bioinformatics/article/35/20/4168/5381539) and predicted fitness using [MSA Transformer](https://www.biorxiv.org/content/10.1101/2021.02.12.430858v1)   
``Rscript script/plot_fit_vs_predict.R``
    - Input files:
      - [./result/Mos99_fit.csv](./result/Mos99_fit.csv)
      - [./data/foldx_msa_transformer.csv](./data/foldx_msa_transformer.csv)
    - Output files:
      - [./graph/fit_vs_MSA_transformer.png](./graph/fit_vs_MSA_transformer.png)
      - [./graph/fit_vs_foldX_ddG.png](./graph/fit_vs_foldX_ddG.png)
