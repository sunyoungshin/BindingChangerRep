#install.packages("devtools")
library(devtools)
#install_github("sunyoungshin/atIndel")
library(atIndel)

load("Simulation/0.Model parameter estimation/FO model parameters.Rdata")
load("v.motif.Rdata")
load("Leukemia/Binding changer analysis/lek indel seq.Rdata") 

motif_library<-v.motif[1:10]
seq_info<-lek_seq_info[1:100]
  
ts<-atIndel::indel_motif_scores(motif_library,seq_info) 
t1<-Sys.time()
tq<-atIndel::indel_p_values(motif_library,seq_info,ts$list,prior,trans_mat,2000)
t2<-Sys.time()
t2-t1

t3<-Sys.time()
tq<-atIndel::indel_p_values(motif_library,seq_info,ts$list,prior,trans_mat,2000,num_cores=4)
t4<-Sys.time()
t4-t3

motif_library<-v.motif[1:500]
seq_info<-lek_seq_info[1:20]

ts<-atIndel::indel_motif_scores(motif_library,seq_info) 
t5<-Sys.time()
tq<-atIndel::indel_p_values(motif_library,seq_info,ts$list,prior,trans_mat,2000,num_cores=12)
t6<-Sys.time()
t6-t5

t7<-Sys.time()
tq<-atIndel::indel_p_values(motif_library,seq_info,ts$list,prior,trans_mat,10000,num_cores=12)
t8<-Sys.time()
t8-t7