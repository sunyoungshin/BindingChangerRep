#install.packages("data.table")
library(data.table)
#install.packages("markovchain")
library(markovchain)
#install.packages("seqinr")
library(seqinr)

load("FO model parameters.Rdata")
#First order
fo_indel_info<- list()
short_seq<-list()
long_seq<-list()
n_indels<-10000

mcDNA= new("markovchain", states = as.character(seq(4)),transitionMatrix = data.matrix(trans_mat), name = "DNA")
max_width<-max(sapply(motif_library, nrow))
insertion_len <-6

#random seq generation
for (i in seq_len(n_indels)) {
  
  t<-as.integer(markovchainSequence(n=max_width * 2 - 3 + insertion_len,
                                    mcDNA,t0 = sample(seq(4), 1,prob=prior),include.t0 = TRUE, useRCpp = FALSE))
  t1<-t[1:(max_width-1)]
  t2<-t[max_width:(max_width+insertion_len-1)]
  t3<-t[(max_width+insertion_len):length(t)]
  
  a1<-n2s(t1-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
  a2<-n2s(t2-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
  a3<-n2s(t3-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
  a4<-n2s(t-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
  aref<-n2s(t[(max_width-1)]-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
  
  ref<-paste(aref, collapse="")
  alt<-paste((c(aref,a2)), collapse="")
  
  fo_indel_info[[i]] <- list(
    inserted_sequence = t,
    insertion_len = insertion_len,
    insertion=1,
    ref=ref,
    alt=alt
  )
}

names(fo_indel_info)<-c(1:n_indels)
save(fo_indel_info, file="1.Null-input.RData")


