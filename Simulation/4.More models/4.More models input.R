#install.packages("seqinr")
library(seqinr)


load("100000 seqs for parameter estimation.RData")

#Fifth-order
#install.packages("SMM")
library(SMM)
#install.packages("data.table")
library(data.table)

   
Ref_seq<-data.frame(Ref)  #100000 seqs 
Ref_seq<-as.matrix(Ref_seq)

E <- c("A","C","G","T")
S = length(E)


## estimation of seqs from hg19
seq<-list()
for (i in 1:100000) {
  seq[[i]] <- strsplit(Ref_seq[i], "")[[1]]
  
}
fit5order<- estimMk(seq = seq, E = E, k = 5)
Ptr5<-data.frame(fit5order$Ptrans)

#make row sum 1
Ptr5$X1<-round(Ptr5$X1,digits=4)
Ptr5$X2<-round(Ptr5$X2,digits=4)
Ptr5$X3<-round(Ptr5$X3,digits=4)
Ptr5$X4<-round(1-Ptr5$X1-Ptr5$X2-Ptr5$X3,digits=4)
which(rowSums(Ptr5)!=1)  #558
Ptr5[558,]<-c(0.3343,0.2375,0.0418,0.3864) #(0.3344,0.2375,0.0418,0.3863)
#To avoid error with simulMk(), we slightly modified row 558 to have row sum 1.  



#random seqs of 5th order
fifth_indel_info<- list()
short_seq<-list() 
long_seq<-list()
n_indels<-10000

max_width<-max(sapply(motif_library, nrow))
insertion_len <-6

#simulate seqs following 5th order model
simu_seq5<-simulMk(E = E, nbSeq =n_indels, lengthSeq =rep(max_width * 2 - 2 + insertion_len,n_indels), Ptrans = as.matrix(Ptr5), init = fit5order$init, k = 5)
  
for (i in 1:length(simu_seq5)) {
  
  t<-s2n(simu_seq5[[i]], levels = c("A", "C", "G", "T"), base4 = FALSE, forceToLower = FALSE)
  
  t1<-t[1:(max_width-1)]
  t2<-t[max_width:(max_width+insertion_len-1)]
  t3<-t[(max_width+insertion_len):length(t)]
  
  a1<-n2s(t1-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
  a2<-n2s(t2-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
  a3<-n2s(t3-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
  a4<-n2s(t-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
  aref<-n2s(t[(max_width-1)]-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
  
  short_seq[[i]]<-paste((c(a1,a3)), collapse="")
  long_seq[[i]]<-paste(a4, collapse="")
  
  ref<-paste(aref, collapse="")
  alt<-paste((c(aref,a2)), collapse="")
  
  fifth_indel_info[[i]] <- list(
    inserted_sequence = as.integer(t),
    insertion_len = insertion_len,
    insertion=1,
    ref=ref,
    alt=alt
  )
}

names(fifth_indel_info)<-c(1:n_indels)





########
#Zero model
GenerateZO<-function(nofrow,nofcol,prob){
  m<- matrix(, nrow =nofrow, ncol = nofcol)
  for(i in 1:nofrow){
    for(j in 1:nofcol){
      m[i,j]<-sample(seq(4),1, replace=TRUE, prob = prob)
    }
  }
  m
}

ZeroNull<-GenerateZO(10000,max_width*2-2+insertion_len,prior)

Zero_indel_info<-list()

for (i in seq_len(n_indels)) {
  t<-ZeroNull[i,]
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
  
  Zero_indel_info[[i]] <- list(
    inserted_sequence = t,
    insertion_len = insertion_len,
    insertion=1,
    ref=ref,
    alt=alt
  )
}
names(Zero_indel_info)<-c(1:n_indels)


save(fifth_indel_info,Zero_indel_info,file="4.More models input.Rdata")

