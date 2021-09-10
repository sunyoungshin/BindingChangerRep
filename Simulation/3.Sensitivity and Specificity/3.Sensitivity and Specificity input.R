#install.packages("data.table")
library(data.table)
#install.packages("gridExtra")
library(gridExtra)
#install.packages("markovchain")
library(markovchain)
#install.packages("seqinr")
library(seqinr)


load("FO model parameters.Rdata")
load("motifs for Simulation.Rdata")
n_indels<-10000

#attributes(trans_mat)$class<-"matrix"
mcDNA= new("markovchain", states = as.character(seq(4)),transitionMatrix = trans_mat, name = "DNA")

GenerateMotif<-function(motif,size){
  n<-nrow(motif)
  m<- matrix(, nrow =size, ncol = n)
  for(i in 1:n){
    m[,i]<-t(sample(seq(4), size, replace=TRUE, prob = motif[i,]))
  }
  m
}

Indel_prop_indel<-function(p,n_indels,motif_library){
  all_indel_info<-list()
  
  width<-nrow(motif_library[[1]])
  insertion_len <-width
  ma<- vector()
  n1<-round(p*n_indels)
  motif<-GenerateMotif(motif_library[[1]],n1)
  
  for (i in seq_len(n1)) {
    t<-as.integer(markovchainSequence(n=width*2-3+insertion_len,
                                      mcDNA,t0 = sample(seq(4), 1,prob=prior),include.t0 = TRUE, useRCpp = FALSE))
    t1<-t[1:(width-1)]
    t2<-t[(width+insertion_len):length(t)]
    t3<-motif[i,]
    
    a1<-n2s(t1-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
    a2<-n2s(t2-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
    a3<-n2s(t3-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
    a4<-n2s(t-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
    aref<-n2s(t[(width-1)]-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
    
    ref<-paste(aref, collapse="")
    alt<-paste((c(aref,a3)), collapse="")
    
    
    all_indel_info[[i]] <- list(
      inserted_sequence = as.integer(c(t1,t3,t2)),
      insertion_len = insertion_len,
      insertion=1,
      ref=ref,
      alt=alt
    )
  }
  
  for (i in (n1+1):n_indels) {
    t<-as.integer(markovchainSequence(n=width * 2 - 3 + insertion_len,
                                      mcDNA,t0 = sample(seq(4), 1,prob=prior),include.t0 = TRUE, useRCpp = FALSE))
    t1<-t[1:(width-1)]
    t2<-t[width:(width+insertion_len-1)]
    t3<-t[(width+insertion_len):length(t)]
    
    a1<-n2s(t1-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
    a2<-n2s(t2-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
    a3<-n2s(t3-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
    a4<-n2s(t-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
    aref<-n2s(t[(width-1)]-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
    
    
    ref<-paste(aref, collapse="")
    alt<-paste((c(aref,a2)), collapse="")
    
    all_indel_info[[i]] <- list(
      inserted_sequence = t,
      insertion_len = insertion_len,
      insertion=1,
      ref=ref,
      alt=alt
    )
  }
  
  names(all_indel_info)<-c(1:n_indels)
  all_indel_info
}

#MSC 0.005,10000
MSC_50_10000_indel_info<-Indel_prop_indel(0.005,10000,motif_library[1])

#MSC 0.01,10000
MSC_100_10000_indel_info<-Indel_prop_indel(0.01,10000,motif_library[1])

#MSC 0.05,10000
MSC_500_10000_indel_info<-Indel_prop_indel(0.05,10000,motif_library[1])

#Ddit3 0.005,10000
Ddit3_50_10000_indel_info<-Indel_prop_indel(0.005,10000,motif_library[2])

#Ddit3 0.01,10000
Ddit3_100_10000_indel_info<-Indel_prop_indel(0.01,10000,motif_library[2])

#Ddit3 0.05,10000
Ddit3_500_10000_indel_info<-Indel_prop_indel(0.05,10000,motif_library[2])

#Hes1 0.005,10000
Hes1_50_10000_indel_info<-Indel_prop_indel(0.005,10000,motif_library[3])

#Hes1 0.01,10000
Hes1_100_10000_indel_info<-Indel_prop_indel(0.01,10000,motif_library[3])

#Hes1 0.05,10000
Hes1_500_10000_indel_info<-Indel_prop_indel(0.05,10000,motif_library[3])

save(MSC_50_10000_indel_info,MSC_100_10000_indel_info,MSC_500_10000_indel_info,
     Ddit3_50_10000_indel_info,Ddit3_100_10000_indel_info,Ddit3_500_10000_indel_info,
     Hes1_50_10000_indel_info,Hes1_100_10000_indel_info,Hes1_500_10000_indel_info,
     file="3.Sensitivity and Specificity input.Rdata")






