#install.packages("motifStack")
library(motifStack)
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)
genome<-BSgenome.Hsapiens.UCSC.hg19

MYCseq<-function(ch,sta,end){
  gr<- GRanges(seqnames=ch,ranges=IRanges(start=sta, end=end))
  a<-getSeq(genome, gr)
  a<-paste(a, sep="", collapse=NULL)
}

#input for MEME-ChIP
MYC_bed<-read.table("NB4_Myc_narrowPeak.bed", header = FALSE,stringsAsFactors = F, sep = "\t")
s<-mapply(MYCseq,MYC_bed$V1,MYC_bed$V2,MYC_bed$V3)
x<-data.frame(c(1:length(s)),s)
MYC.fasta = dataframe2fas(x, file="MYC.fasta")  


#Top 3 motifs from MEME-ChIP results
MYC1<-matrix(1:24, nrow =6, ncol = 4)
MYC1[1,]<-c(0.520285,0.124954,0.3547609999,0.0000000001)
MYC1[2,]<-c(0.0000000001,0.0000000001,0.9999999997,0.0000000001)
MYC1[3,]<-c(0.0000000001,0.0000000001,0.9999999997,0.0000000001)
MYC1[4,]<-c(0.9999999997,0.0000000001,0.0000000001,0.0000000001)
MYC1[5,]<-c(0.9999999997,0.0000000001,0.0000000001,0.0000000001)
MYC1[6,]<-c( 0.377446999,0.0000000001,0.6225529999,0.0000000001)

MYC2<-matrix(1:28, nrow =7, ncol = 4)
MYC2[1,]<-c(0.0000000001,0.9999999997,0.0000000001,0.0000000001)
MYC2[2,]<-c(0.9999999997,0.0000000001,0.0000000001,0.0000000001)
MYC2[3,]<-c(0.0000000001,0.9999999997,0.0000000001,0.0000000001)
MYC2[4,]<-c(0.6977549999,0.0000000001,0.3022449999,0.0000000001)
MYC2[5,]<-c(0.0000000001,0.0000000001,0.0000000001,0.9999999997)
MYC2[6,]<-c(0.0000000001,0.0000000001,0.9999999997,0.0000000001)
MYC2[7,]<-c(0.0000000001,0.2889209999,0.4651300000,0.2459490000)

MYC3<-matrix(1:32, nrow =8, ncol = 4)
MYC3[1,]<-c(0.5249419999,0.0000000001,0.4750579999,0.0000000001)
MYC3[2,]<-c(0.0000000001,0.0000000001,0.0000000001,0.9999999997)
MYC3[3,]<-c(0.0000000001,0.0000000001,0.9999999997,0.0000000001)
MYC3[4,]<-c(0.9999999997,0.0000000001,0.0000000001,0.0000000001)
MYC3[5,]<-c(0.0000000001,0.5277389999,0.4722609999,0.0000000001)
MYC3[6,]<-c(0.0000000001,0.0000000001,0.0000000001,0.9999999997)
MYC3[7,]<-c(0.0000000001,0.9999999997,0.0000000001,0.0000000001)
MYC3[8,]<-c(0.9999999997,0.0000000001,0.0000000001,0.0000000001)

myMYC<-list()
myMYC$MYC1<-MYC1
myMYC$MYC2<-MYC2
myMYC$MYC3<-MYC3

p1<-t(myMYC$MYC1)
rownames(p1)<-c("A","C","G","T")
fig51<-plotMotifLogo(p1,xlab="")


p2<-t(myMYC$MYC2)
rownames(p2)<-c("A","C","G","T")
fig52<-plotMotifLogo(p2,xlab="")

p3<-t(myMYC$MYC3)
rownames(p3)<-c("A","C","G","T")
fig53<-plotMotifLogo(p3,xlab="")

save(myMYC,file="MYC motif.Rdata")
