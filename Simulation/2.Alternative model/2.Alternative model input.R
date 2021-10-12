#install.packages("data.table")
library(data.table)
#install.packages("markovchain")
#install.packages('markovchain', repos = 'https://cloud.r-project.org')
library(markovchain)
#install.packages("seqinr")
library(seqinr)

load("FO model parameters.Rdata")
load("motifs for Simulation.Rdata")


#generate sequence pairs
mcDNA= new("markovchain", states = as.character(seq(4)),transitionMatrix = trans_mat, name = "DNA")

GenerateMotif<-function(motif,size){
  n<-nrow(motif)
  m<- matrix(nrow =size, ncol = n)
  for(i in 1:n){
    m[,i]<-t(sample(seq(4), size, replace=TRUE, prob = motif[i,]))
  }
  m
}

all_indel_info2<- list()
all_indel_info3<- list()
GenerateSeq2<-function(n_indels,m,motif){
  motifpool<-GenerateMotif(motif,n_indels)
  width<-nrow(motif)
  all_indel_info<- list()
  for (i in seq_len(n_indels)) {
    t<-as.integer(markovchainSequence(n=width*3-3,
                                      mcDNA,t0 = sample(seq(4), 1,prob=prior),include.t0 = TRUE, useRCpp = FALSE))
    
    t1<-t[1:(width-1)]
    t2<-t[(2*width):length(t)]
    t3<-motifpool[i,1:m]
    t4<-t[(width+m):(2*width-1)] #rest of InDel  
    
    a3<-n2s(t3-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
    a4<-n2s(t4-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
    aref<-n2s(t[(width-1)]-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
    
    ref<-paste(aref, collapse="")
    alt<-paste((c(aref,a3,a4)), collapse="")
    
    all_indel_info[[i]] <- list(
      inserted_sequence = as.integer(c(t1,t3,t4,t2)),
      insertion_len = width,
      insertion=1,
      ref=ref,
      alt=alt
    )
  }
  names(all_indel_info)<-c(1:n_indels)
  all_indel_info2<-all_indel_info
  names(all_indel_info2)<-c((n_indels*(m-1)+1):(n_indels*(m-1)+n_indels))
  all_indel_info3<-c(all_indel_info3,all_indel_info2)
  all_indel_info3
}


#MSC
k1<-mapply(function(x) GenerateSeq2(2000,x,motif_library[1][[1]]),1:10,SIMPLIFY = FALSE)
MSC_AL_indel_info<-do.call(c,k1)

#Ddit3::Cebpa
k1<-mapply(function(x) GenerateSeq2(2000,x,motif_library[2][[1]]),1:12,SIMPLIFY = FALSE)
Ddit3_AL_indel_info<-do.call(c,k1)

#Hes1
k1<-mapply(function(x) GenerateSeq2(2000,x,motif_library[3][[1]]),1:10,SIMPLIFY = FALSE)
Hes1_AL_indel_info<-do.call(c,k1)

save(MSC_AL_indel_info,Ddit3_AL_indel_info,Hes1_AL_indel_info,file="2.AL input.Rdata")