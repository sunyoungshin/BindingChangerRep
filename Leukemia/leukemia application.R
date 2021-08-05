library(atIndel)
library(seqinr)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)
library(data.table)

I1<- read.table("input.txt", header = TRUE,stringsAsFactors = F, sep = "\t")
Clean_somatic<-function(b){
  b_urows<-row.names(unique(b[c("CHROM","START","END","REF","ALT")]))
  b<- b[b_urows,c("CHROM","START","END","REF","ALT")]
  b
}

I2<-Clean_somatic(I1)
I2$nref<-nchar(I2$REF)
I2$nalt<-nchar(I2$ALT)
I2$m<-abs(I2$nref-I2$nalt)
I2$Insertion<-ifelse(I2$nref<I2$nalt,1,0)
I2$id<-paste(I2$CHROM,I2$START,I2$nref,I2$nalt,sep=":", collapse=NULL)

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

L<-mapply(FixLong,I2$CHROM,I2$START,I2$END,I2$ALT,I2$nref,I2$nalt,I2$m,100,I2$Insertion)
L2<-mapply(FixLong,I3$CHROM,I3$START,I3$END,I3$ALT,I3$nref,I3$nalt,I3$m,100,I3$Insertion)

lek_seq_info<-list()
for (i in 1:length(L2)){
  b<-strsplit(L2[i], "")[[1]]
  insertion_len<-I3[i,8]
  insertion<-I3[i,9]
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


#Model parameter estimation 
getRef<-function(ch,sta,lth) {
  gr<- GRanges(seqnames=ch,ranges=IRanges(start=sta, width=lth))
  getSeq(genome, gr)
}

Ref<-getRef(I2$CHROM,I2$START,1000)    
Ref_seq<-data.frame(Ref)
Ref_seq<-as.matrix(Ref_seq)

codes <- seq(4)
names(codes) <- c("A", "C", "G", "T")
refmat <- sapply(Ref_seq, function(x) codes[strsplit(x, "")[[1]]])

colnames(refmat) <- colnames(refmat) <- NULL
transition<-trans.matrix(refmat,prob=F) 
prior<- apply(transition, 1, sum)
prior<- prior/sum(prior)
trans_mat<- transition/apply(transition, 1, sum)
names(prior) <- colnames(trans_mat) <- rownames(trans_mat) <- c("A", "C", "G", "T")

prior1000<-prior
trans1000<-trans_mat


#MYC motif from MEME-Chip  
#The motifs are the top 3 motifs from NB4 MYCcombined.meme.txt in directory MYC motif.
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


Lek_score<-atIndel::indel_motif_scores(myMYC,lek_seq_info)
Lek_result_10000<-atIndel::indel_p_values(myMYC,lek_seq_info,Lek_score$list,prior1000,trans1000,10000)

Lek_result_10000<-data.table(Lek_result_10000)
Lek_result2<-Lek_result_10000[,pval_bh:=p.adjust(p_value_change.rank,method="BH"),by=motif]
Lek_result3<-Lek_result2[pval_bh<0.10&(ref.pval<0.05|mutation.pval<0.05),]  #7 indels
Lek_result3$direction<-ifelse(Lek_result3$ref.pval<Lek_result3$mutation.pval,"-","+")

#add Nearest genes of the InDels
Neargene<- read.table("0509indel.txt", header = FALSE,stringsAsFactors = F, sep = "\t")
colnames(Neargene)<-c("id","Nearest.genes")
f1<-merge(Lek_result3,Neargene,by="id") 
IDinfo<-I3[I3$id %in% Lek_result3$id,c(1:3,12)]
f2<-merge(f1,IDinfo,by="id")  


final_table<-f2[,c(15:17,1:14)]


#######plot the results   

which(names(lek_seq_info)=="chr15:102265302:2:1") #1676
atIndel::plot_indel_binding(lek_seq_info[1676],Lek_score$list,myMYC[1])

which(names(lek_seq_info)=="chr19:23973237:2:1") #2663
atIndel::plot_indel_binding(lek_seq_info[2663],Lek_score$list,myMYC[1])

which(names(lek_seq_info)=="chr19:48758981:2:1") #2825
atIndel::plot_indel_binding(lek_seq_info[2825],Lek_score$list,myMYC[2])

which(names(lek_seq_info)=="chr21:36194111:2:1") #3514
atIndel::plot_indel_binding(lek_seq_info[3514],Lek_score$list,myMYC[1])

which(names(lek_seq_info)=="chr21:38593248:2:1") #3526
atIndel::plot_indel_binding(lek_seq_info[3526],Lek_score$list,myMYC[1])

which(names(lek_seq_info)=="chr4:83189692:2:1") #4161
atIndel::plot_indel_binding(lek_seq_info[4161],Lek_score$list,myMYC[1])

which(names(lek_seq_info)=="chr7:151860748:1:2") #4802
atIndel::plot_indel_binding(lek_seq_info[4802],Lek_score$list,myMYC[1])





