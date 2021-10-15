#install.packages("data.table")
library(data.table)
#install.packages("devtools")
library(devtools)
#install_github("sunyoungshin/atIndel")
library(atIndel)
#install.packages("ggplot2")
library(ggplot2)

load("MYC motif.Rdata")
load("lek indel seq.Rdata") 
load("leukemia parameters.Rdata")

load("Simulation/1.Null model-First order/1.Null-input.RData")
load("Leukemia/lek indel seq.Rdata")

Lek_score<-atIndel::indel_motif_scores(myMYC,lek_seq_info)
Lek_result_10000<-atIndel::indel_p_values(myMYC,lek_seq_info,Lek_score$list,prior1000,trans1000,10000)

Lek_result_10000<-data.table(Lek_result_10000)
Lek_result2<-Lek_result_10000[,pval_bh:=p.adjust(p_value_change.rank,method="BH"),by=motif]
Lek_result3<-Lek_result2[pval_bh<0.10&(ref.pval<0.05|mutation.pval<0.05),]  #7 indels
Lek_result3$direction<-ifelse(Lek_result3$ref.pval<Lek_result3$mutation.pval,"-","+")

#add Nearest genes of the InDels
Neargene<- read.table("InDel nearest genes.txt", header = FALSE,stringsAsFactors = F, sep = "\t")
colnames(Neargene)<-c("id","Nearest.genes")
f1<-merge(Lek_result3,Neargene,by="id") 
IDinfo<-I3[I3$id %in% Lek_result3$id,c(1:3,12)]
f2<-merge(f1,IDinfo,by="id")  


table1213<-f2[,c(15:17,1:14)]

mmm<-Lek_result_10000[Lek_result_10000$motif=="MYC1",]
m<-abs(nchar(mmm$ref)-nchar(mmm$alt))
table(m)
m2<-ifelse(m<=10&m>=6,"6-10",ifelse(m<=20&m>=11,"11-20",ifelse(m>=21,"21+",m)))
table5<-table(m2)

#######plot the results (Figure 7) 

which(names(lek_seq_info)=="chr15:102265302:2:1") #1676
atIndel::plot_indel_binding(lek_seq_info[1676],Lek_score$list,myMYC[1])

which(names(lek_seq_info)=="chr19:23973237:2:1") #2663
atIndel::plot_indel_binding(lek_seq_info[2663],Lek_score$list,myMYC[1])

which(names(lek_seq_info)=="chr19:48758981:2:1") #2825
atIndel::plot_indel_binding(lek_seq_info[2825],Lek_score$list,myMYC[2])

which(names(lek_seq_info)=="chr21:36194111:2:1") #3514
atIndel::plot_indel_binding(lek_seq_info[3514],Lek_score$list,myMYC[1])

which(names(lek_seq_info)=="chr21:38593248:2:1") #3526
atIndel::plot_indel_binding(lek_seq_info[3526],Lek_score$list,myMYC[1])

which(names(lek_seq_info)=="chr4:83189692:2:1") #4161
atIndel::plot_indel_binding(lek_seq_info[4161],Lek_score$list,myMYC[1])

which(names(lek_seq_info)=="chr7:151860748:1:2") #4802
atIndel::plot_indel_binding(lek_seq_info[4802],Lek_score$list,myMYC[1])

save(table1213,table5,file="leukemia output.Rdata")
