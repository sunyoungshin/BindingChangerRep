BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)
genome<-BSgenome.Hsapiens.UCSC.hg19

getRef<- function(ch,sta,lth) {
  gr<- GRanges(seqnames=ch,ranges=IRanges(start=sta, width=lth))
  getSeq(genome, gr)
}


ch1<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chrX")
ch2<-c("chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15")
ch3<-c("chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrY")
ch<-c(ch1,ch2,ch3)
Chromosome<-sample(ch,100000,replace =TRUE)

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

Start.sites<-unlist(lapply(Chromosome, get_startpos))

G<-data.frame(Chromosome,Start.sites)


Ref<-getRef(G$Chromosome,G$Start.sites,100)    #100000 seqs 
Ref_seq<-data.frame(Ref)
Ref_seq<-as.matrix(Ref_seq)

codes <- seq(4)
names(codes) <- c("A", "C", "G", "T")
refmat <- sapply(Ref_seq, function(x) codes[strsplit(x, "")[[1]]])

colnames(refmat) <- colnames(refmat) <- NULL
transition<-trans.matrix(refmat,prob=F) 
prior<- apply(transition, 1, sum)
prior <- prior/sum(prior)
trans_mat<- transition/apply(transition, 1, sum)
names(prior) <- colnames(trans_mat) <- rownames(trans_mat) <- c("A", "C", "G", "T")

atIndel_tran<-trans_mat
atIndel_prior<-prior


trans.matrix<-function(X, prob=T)
{
  tt <- table( c(X[,-ncol(X)]), c(X[,-1]) )
  if(prob) tt <- tt / rowSums(tt)
  tt
}








