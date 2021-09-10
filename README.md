# Introduction

This repository contains all publicly available data, the scripts to run simulations and analysis, and the scripts to generate plots and tables for the paper "Scalable DNA-protein binding changer test for insertion and deletion of bases in the genome".

This code was written and tested primarily on R 4.0.3 on a laptop (Windows 10) equipped with an i7 processor. It was tested on Macbook (macOS Catalina) equipped with an i5 processor as well.

The Binding Changer test runs with the package [atIndel](https://github.com/sunyoungshin/atIndel).

# R package installation

```{r}
library("devtools")
devtools::install_github("sunyoungshin/atIndel")   
```
The leukemia data analysis and simulations require the additional packages: 
`BSgenome.Hsapiens.UCSC.hg19`, `data.table`, `markovchain`, `parallel`, `ggplot2`, `PRROC`, `plotROC`, `SMM`, `seqinr`, `plyr`, `seqRFLP`.

# Data

## For simulation:
The motif data are originally from [JASPAR 2018 vertebrata motif library](http://jaspar2018.genereg.net/downloads/) The three motifs for the simulation studies are stored in `simulation.Rdata`. The sequence data are generated for each of the four simulation studies. The sequence generating R codes are in the four
.R files under the directory `Simulation`.

## For leukemia InDel analysis:
The motif data are obtained from [NB4 cell line ENCODE ChIP-seq data](https://www.encodeproject.org/files/ENCFF001VQK/) with [MEME-ChIP](https://meme-suite.org/meme/tools/meme-chip). The leukemia InDel data were provided by Jian Xu lab from UTSW. The data will be provided upon request. The point of contact is [Dr. Jian Xu](mailto:Jian.Xu@UTSouthwestern.edu). 

# Reproducibility workflow

## For simulation:
Motif MSC, Ddit3::Cebpa and Hes1 are from [JASPAR 2018 vertebrata motif library](http://jaspar2018.genereg.net/downloads/) and they are saved in `motifs for `Simulation.Rdata`.

The directory `Simulation` has five subdirectories for the simulation studies.

### 0. Model parameter estimation
Model parameters were obtained with R codes in `0.Model parameter estimation.R`.
The random sequences used for model parameter estimation are saved in `100000 seqs for parameter estimation.RData`.
The model parameters are saved in `FO model parameters.Rdata`.
`0.Model parameter estimation.R` produces Table D1 & D2.


.R files contain codes for input sequence pairs generation.
.Rdata files provide the input random sequence pairs and output tables and figures.

### 1. Simulation under the null
The random sequence pairs are generated with R codes in `1.Null model-First order input.R` and they are saved in `1.Null-input.Rdata`.

`1.Null model-first order test & summary.R` produces Table 1 & D3, Figure 2, 3 & D2.

### 2. Simulation under the alternative
The random sequence pairs are generated with R codes in `2.Alternative model input.R` and they are saved in `2.AL input.Rdata`.

`2.Alternative model test & summary.R` produces Table 2 & D4, Figure 4 & D3.

### 3. Sensitivity and specificity analysis 
The sequence pairs are generated with R codes in `3.Sensitivity and Specificity input.R` and they are saved
in `3.Sensitivity and Specificity input.Rdata`.

### Subdirectory FIMO SS
The input of FIMO are generated with R codes in `generate fimo fasta.R` and they are saved in the directory. 
There are three subdirectories in this directory. Each is for the corresponding motif. The inputs .fasta and .txt for the motif are 
saved in the subdirectory.
FIMO results are obtained using [FIMO](https://meme-suite.org/meme/tools/fimo) with p-value cutoff 1 and they are saved in the subdirectories.

`3.Sensitivity and Specificity test & summary.R` produces Figure D4.

### 4. More background models
The sequence pairs are generated with R codes in `4.More models input.R` and they are saved in `4.More models input.Rdata`.

`4.More models test & summary.R` produces Table D5 & D6.

## For leukemia InDel analysis:
The directory `Leukemia` is for leukemia InDel analysis.

### 1. MYC
**1.1** The input of MEME-ChIP is a fasta file generated from `NB4_Myc_narrowPeak.bed` with R code in `MYC.R`.

`NB4 MYCcombined.meme.txt` shows the motifs from MEME-ChIP.

`MYC.R` produces Figure 5. 

**1.2** We select the top 3 motifs for our analysis and saved the motifs in `MYC motif.Rdata`. 

### 2. Binding changer analysis
**2.1** The model parameters for leukemia InDel analysis were obtained with R codes in `leukemia model parameter.R` and they are saved in `leukemia parameters.Rdata`.
The input R object _lek_seq_info_ of the BC test are generated with R codes in `leukemia indel seq.R`.

**2.2** `leukemia test & plots.R` ran the BC test and produced Table E9 & E10 and Figure 6. 

### 3. FIMO based analysis
`leukemia FIMO analysis.R` contains the codes for FIMO-based analysis.
The input and output of FIMO are saved in directory `FIMO based analysis`. 

`Running time.R` produces Table D7 & D8.
`v.motif.Rdata` contains the motifs for testing the running time.

