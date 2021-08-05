#Fifth-order
library(SMM)
library(data.table)

library(BSgenome.Hsapiens.UCSC.hg19)
genome<-BSgenome.Hsapiens.UCSC.hg19

getRef<- function(ch,sta,lth) {
  gr<- GRanges(seqnames=ch,ranges=IRanges(start=sta, width=lth))
  getSeq(genome, gr)
}

chr_length<-function(x){
  ifelse(x=="chr1",249250621,ifelse(x=="chr2",243199373,ifelse(x=="chr3",198022430,ifelse(x=="chr4",191154276,
                                                                                          ifelse(x=="chr5",180915260,ifelse(x=="chr6",171115067,ifelse(x=="chr7",159138663,ifelse(x=="chrX",155270560,
                                                                                                                                                                                  ifelse(x=="chr8",146364022,ifelse(x=="chr9",141213431,ifelse(x=="chr10",135534747,ifelse(x=="chr11",135006516,
                                                                                                                                                                                                                                                                           ifelse(x=="chr12",133851895,ifelse(x=="chr13",115169878,ifelse(x=="chr14",107349540,
                                                                                                                                                                                                                                                                                                                                          ifelse(x=="chr15",102531392,ifelse(x=="chr16",90354753,ifelse(x=="chr17",81195210,
                                                                                                                                                                                                                                                                                                                                                                                                        ifelse(x=="chr18",78077248,ifelse(x=="chr20",63025520,ifelse(x=="chrY",59373566,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                     ifelse(x=="chr19",59128983,ifelse(x=="chr21",48129895,51304566)))))))))))))))))))))))
}

get_startpos<-function(x){
  l<-chr_length(x)
  e<-l-10000
  sample(10000:e,1)  
}

ch1<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chrX")
ch2<-c("chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15")
ch3<-c("chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrY")
ch<-c(ch1,ch2,ch3)
Chromosome<-sample(ch,100000,replace =TRUE)
Start.sites<-unlist(lapply(Chromosome, get_startpos))

G<-data.frame(Chromosome,Start.sites)

Ref<-getRef(G$Chromosome,G$Start.sites,100)    #100000 seqs 
Ref_seq<-data.frame(Ref)
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


Ptr5[558,]<-c(0.3343,0.2375,0.0418,0.3864)
Ptr5[282,]<-c(0.2338,0.5781,0.0579,0.1302)
Ptr5[359,]<-c(0.114,0.5084,0.2724,0.1052)
Ptr5[382,]<-c(0.277,0.3702,0.0336,0.3192)

fifth_indel_info<- list()
short_seq<-list() 
long_seq<-list()
n_indels<-10000

max_width<-max(sapply(motif_library, nrow))
insertion_len <-6

#simulate seqs following k-th order model
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

Fifth_score<-atIndel::indel_motif_scores(motif_library,fifth_indel_info) 
Fifth_pvalue<-atIndel::indel_p_values(motif_library,fifth_indel_info,Fifth_score$list,atIndel_prior,atIndel_tran,2000)

#summary table
at_result_table<-function(x){
  p<-c(0.01,0.05,0.1)
  atIndel_p<-mcmapply(function(p) emp(p,x$p_value_change.rank),p)
  data.table(p,atIndel_p)
}

emp<-function(p,x){
  length(which(x<p))/length(x)
}

FO5_1<-at_result_table(Fifth_pvalue[Fifth_pvalue$motif=="MSC",])
FO5_1

FO5_2<-at_result_table(Fifth_pvalue[Fifth_pvalue$motif=="Ddit3::Cebpa",])
FO5_2

FO5_3<-at_result_table(Fifth_pvalue[Fifth_pvalue$motif=="Hes1",])
FO5_3



########
#Zero model
max_width<-max(sapply(motif_library, nrow))
insertion_len <-6

GenerateZO<-function(nofrow,nofcol,prob){
  m<- matrix(, nrow =nofrow, ncol = nofcol)
  for(i in 1:nofrow){
    for(j in 1:nofcol){
      m[i,j]<-sample(seq(4),1, replace=TRUE, prob = prob)
    }
  }
  m
}

FimoNull<-GenerateZO(10000,max_width*2-3+insertion_len,atIndel_prior)

Zero_indel_info<-list()

for (i in seq_len(n_indels)) {
  t<-FimoNull[i,]
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

ZO_score<-atIndel::indel_motif_scores(motif_library,Zero_indel_info) 
ZO_pvalue<-atIndel::indel_p_values(motif_library,Zero_indel_info,ZO_score$list,atIndel_prior,atIndel_tran,2000)


#table
FO0_1<-at_result_table(ZO_pvalue[ZO_pvalue$motif=="MSC",])
FO0_1

FO0_2<-at_result_table(ZO_pvalue[ZO_pvalue$motif=="Ddit3::Cebpa",])
FO0_2

FO0_3<-at_result_table(ZO_pvalue[ZO_pvalue$motif=="Hes1",])
FO0_3

save(fifth_indel_info,Fifth_pvalue,Zero_indel_info,ZO_pvalue,file="moremodels.Rdata")

