#install.packages("seqinr")
library(seqinr)
#install.packages("plyr")
library(plyr)
#install.packages("seqRFLP")
library(seqRFLP)


load("Simulation/3.Sensitivity and Specificity/3.Sensitivity and Specificity input.Rdata")

seq_info2fasta<-function(seq_info,p1,p2){
  n<-length(seq_info)
  short_seq<-list()
  long_seq<-list()
  for (i in seq_len(n)) {
    t<-n2s(seq_info[i][[1]]$inserted_sequence-1,levels = c("A", "C", "G", "T"), base4 = TRUE)
    
    m<-seq_info[i][[1]]$insertion_len
    seq_len<-length(t)
    width<-(seq_len-m)/2
    
    t1<-t[1:width]
    t2<-t[(width+1):(width+m)]
    t3<-t[(width+m+1):seq_len]
    
    short_seq[[i]]<-paste((c(t1,t3)), collapse="")
    long_seq[[i]]<-paste(t, collapse="")
    
  }
  Y<-c(1:n)
  x1<-data.frame(Y,ldply (short_seq, data.frame))
  x2<-data.frame(Y,ldply (long_seq, data.frame))
  short.fasta = dataframe2fas(x1, file=p1)
  long.fasta = dataframe2fas(x2, file=p2)
}

seq_info2fasta(MSC_50_10000_indel_info,"MSC_50_10000_short.fasta","MSC_50_10000_long.fasta")
seq_info2fasta(MSC_100_10000_indel_info,"MSC_100_10000_short.fasta","MSC_100_10000_long.fasta")
seq_info2fasta(MSC_500_10000_indel_info,"MSC_500_10000_short.fasta","MSC_500_10000_long.fasta")

seq_info2fasta(Ddit3_50_10000_indel_info,"Ddit3_50_10000_short.fasta","Ddit3_50_10000_long.fasta")
seq_info2fasta(Ddit3_100_10000_indel_info,"Ddit3_100_10000_short.fasta","Ddit3_100_10000_long.fasta")
seq_info2fasta(Ddit3_500_10000_indel_info,"Ddit3_500_10000_short.fasta","Ddit3_500_10000_long.fasta")

seq_info2fasta(Hes1_50_10000_indel_info,"Hes1_50_10000_short.fasta","Hes1_50_10000_long.fasta")
seq_info2fasta(Hes1_100_10000_indel_info,"Hes1_100_10000_short.fasta","Hes1_100_10000_long.fasta")
seq_info2fasta(Hes1_500_10000_indel_info,"Hes1_500_10000_short.fasta","Hes1_500_10000_long.fasta")
