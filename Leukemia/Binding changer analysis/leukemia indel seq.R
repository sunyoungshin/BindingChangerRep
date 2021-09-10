#install.packages("data.table")
library(data.table)
#install.packages("gridExtra")
library(gridExtra)
#install.packages("devtools")
library(devtools)
install_github("sunyoungshin/BindingChanger")
#install.packages("markovchain")
library(markovchain)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("parallel")
library(parallel)
#install.packages("seqinr")
library(seqinr)


I1<- read.table("input.txt", header = TRUE,stringsAsFactors = F, sep = "\t") 
#input.txt can be obtained by contacting Dr.Jian Xu at UTSW
#Clean_somatic() is a duplicate remover.
Clean_somatic<-function(b){
  b_urows<-row.names(unique(b[c("CHROM","START","END","REF","ALT")]))
  b<- b[b_urows,c("CHROM","START","END","REF","ALT")]
  b
}

#I2 includes 6,760 InDels
I2<-Clean_somatic(I1)  
I2$nref<-nchar(I2$REF)
I2$nalt<-nchar(I2$ALT)
I2$m<-abs(I2$nref-I2$nalt)
I2$Insertion<-ifelse(I2$nref<I2$nalt,1,0)
I2$id<-paste(I2$CHROM,I2$START,I2$nref,I2$nalt,sep=":", collapse=NULL)

#I3 includes 5,737 InDels
I3<-I2[I2$END==I2$START+I2$nref-1,]
I3$id[1082]<-"chr12:1100118:1:5(2)"

genome<-BSgenome.Hsapiens.UCSC.hg19

FixLong<- function(ch,sta,end,alt,nref,nalt,m,wid,Insertion) {
  lth<-2*wid+m
  if(Insertion==0){
    gr<- GRanges(seqnames=ch,ranges=IRanges(start=sta+nalt-wid, width=lth))
    a<-getSeq(genome, gr)
    a<-paste(a, sep="", collapse=NULL)
  }else{
    grbef<- GRanges(seqnames=ch,ranges=IRanges(start=sta+nref-wid, width=wid-nref))
    graft<- GRanges(seqnames=ch,ranges=IRanges(start=end+1, width=wid))
    bef<-getSeq(genome, grbef)
    aft<-getSeq(genome, graft)
    a<-paste(bef,alt,aft, sep="", collapse=NULL)
  }
  a
}

L2<-mapply(FixLong,I3$CHROM,I3$START,I3$END,I3$ALT,I3$nref,I3$nalt,I3$m,100,I3$Insertion)

lek_seq_info<-list()
for (i in 1:length(L2)){
  b<-strsplit(L2[i], "")[[1]]
  insertion_len<-I3[i,9]
  insertion<-I3[i,10]
  ref<-I3[i,4]
  alt<-I3[i,5]
  t<-s2n(b, levels = c("A", "C", "G", "T"), base4 = FALSE, forceToLower = FALSE)
  a<-I2$id[i]
  lek_seq_info[[i]] <- list(
    inserted_sequence = as.integer(t),
    insertion_len = insertion_len,
    insertion=insertion,
    ref=ref,
    alt=alt)
}
names(lek_seq_info)<-I3$id
save(I2,lek_seq_info,file="lek indel seq.Rdata") 
