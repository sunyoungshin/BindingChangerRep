#install.packages("data.table")
library(data.table)
#install.packages("seqinr")
library(seqinr)
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)
genome<-BSgenome.Hsapiens.UCSC.hg19


load("")
#Model parameter estimation 
getRef<-function(ch,sta,lth) {
  gr<- GRanges(seqnames=ch,ranges=IRanges(start=sta, width=lth))
  getSeq(genome, gr)
}

trans.matrix<-function(X, prob=T)
{
  tt <- table( c(X[,-ncol(X)]), c(X[,-1]) )
  if(prob) tt <- tt / rowSums(tt)
  tt
}

Lek_Ref<-getRef(I2$CHROM,I2$START,1000)    
Lek_Ref_seq<-data.frame(Lek_Ref)
Lek_Ref_seq<-as.matrix(Lek_Ref_seq)

codes <- seq(4)
names(codes) <- c("A", "C", "G", "T")
refmat <- sapply(Lek_Ref_seq, function(x) codes[strsplit(x, "")[[1]]])

colnames(refmat) <- colnames(refmat) <- NULL
transition<-trans.matrix(refmat,prob=F) 
prior1000<- apply(transition, 1, sum)
prior1000<- prior1000/sum(prior1000)
trans1000<- transition/apply(transition, 1, sum)
names(prior1000) <- colnames(trans1000) <- rownames(trans1000) <- c("A", "C", "G", "T")

save(prior1000,trans1000,file="leukemia parameters.Rdata")
