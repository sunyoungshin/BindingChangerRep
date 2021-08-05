Introduction：
This repository contains all the PUBLIC? data, scripts to run simulations and analysis, and scripts to generate plots and tables
for the paper "Scalable DNA-protein binding changer test for insertion and deletion of bases in the genome".

This code was written and tested primarily on R 4.0.3 on a laptop (Windows 10) equipped with an i7 processor.、
It was tested on Macbook (macOS Catalina) equipped with an i5 processor as well.

The Binding Changer test runs with the package "Binding Changer" available at: ?????


R package installation:
library("devtools")
devtools::install_github("??", subdir = "????")

Additionally, the leukemia data analysis and simulations require the additional packages: 
BSgenome.Hsapiens.UCSC.hg19, data.table, markovchain, parallel, ggplot2, PRROC, plotROC, SMM


Data:
For simulation:
The motif data are originally from JASPAR vertebrata motif library (http://jaspar.genereg.net/search?q=&collection=CORE&tax_group=vertebrates)
The three motifs for the simulation studies are stored in simulation.Rdata.
The sequence data are generated for each of the four simulation studies. The sequence generating R codes are in the four
.R files under directory "Simulation codes".


For leukemia InDel analysis:
The motif data are generated from NB4 cell line ENCODE ChIP-seq data (https://www.encodeproject.org/files/ENCFF001VQK/)  
with MEME-ChIP (https://meme-suite.org/meme/tools/meme-chip).
The leukemia InDel data were provided by Jian Xu lab from UTSW. The Lab will provide the data upon request. 
The point of contact is Dr. Jian Xu at Jian.Xu@UTSouthwestern.edu. 



Reproducibility workflow：
For simulation:
0.a Motif MSC, Ddit3::Cebpa and Hes1 from JASPAR vertebrata motif library.
0.b Model parameters were obtained with R codes in 0.Model parameter estimation.R.
The model parameters and motifs are saved in simulation.Rdata.

.R files contain codes for input sequence pairs generation and .Rdata files provide the input random sequence pairs.
1. Simulation under the null
The random sequence pairs could be generated with R codes in 1.Null model-First order.R and they are saved in simulation FO.Rdata.

2. Simulation under the alternative
The random sequence pairs could be generated with R codes in 2.Alternative model.R and they are saved in simulation AL.Rdata.

3. Sensitivity and specificity analysis 
The sequence pairs could be generated with R codes in 3.Sensitivity and Specificity.R and they are saved
in simulation SS.Rdata.
FIMO results are obtained using FIMO (https://meme-suite.org/meme/tools/fimo) with p-value cutoff 1 and they are saved in "FIMO SS".

4. More background models
The sequence pairs could be generated with R codes in 4.More models.R and they are saved in more models.Rdata.


For leukemia InDel analysis:
The input of MEME-ChIP is a fasta file generated from NB4_Myc_narrowPeak.bed with R code in MYC.R.
NB4 MYCcombined.meme.txt shows the motifs from MEME-ChIP.
We select the top 3 motifs for our analysis.
leukemia application.R ran the BC test. 



.R files contain codes for output tables and figures.
0.Model parameter estimation.R produces Table 3-4, Figure 2.
1.Null model-First order.R produces Table 5-8, Figure 3-4.
2.Alternative model.R produces Table 9-14, Figure 5-10.
3.Sensitivity and Specificity.R produces Figure 11-13.
4.More models.R produces Table 15-16.
leukemia application.R produces Table 17-18, Figure 14-15.