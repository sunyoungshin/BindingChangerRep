#install.packages("seqRFLP")
library(seqRFLP)
#install.packages("seqinr")
library(seqinr)
#install.packages("plyr")
library(plyr)

load("lek indel seq.Rdata")

ref_seq<-list()
alt_seq<-list()
for (i in 1:length(lek_seq_info)) {
  
  t<-lek_seq_info[i][[1]]$inserted_sequence
  m<-lek_seq_info[i][[1]]$insertion_len
  t1<-t[1:100]
  t2<-t[101:(100+m)]
  t3<-t[(100+m+1):length(t)]

  a1<-n2s(t1-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
  
  a2<-n2s(t3-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
  
  a<-n2s(t-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
  
  if(lek_seq_info[i][[1]]$insertion==0){
  ref_seq[[i]]<-paste((c(a1,a2)), collapse="")
  alt_seq[[i]]<-paste(a, collapse="")
  }else{
    ref_seq[[i]]<-paste(a, collapse="")
    alt_seq[[i]]<-paste((c(a1,a2)), collapse="") 
  }
}

Y<-names(lek_seq_info)
x1<-data.frame(Y,ldply (ref_seq, data.frame))
x2<-data.frame(Y,ldply (alt_seq, data.frame))

#generate fasta files---input of FIMO
ref.fasta = dataframe2fas(x1, file="ref.fasta")
alt.fasta = dataframe2fas(x2, file="alt.fasta")





#Direct p-value from FIMO
FIMO_directp<-function(p1,p2){
  pofshort<-data.table(read.table(p1, sep = '\t', header = TRUE))  
  poflong<-data.table(read.table(p2, sep = '\t', header = TRUE)) 
  fimos<-fimo_table_dir(pofshort)
  colnames(fimos)[2]<-"ref_pval"
  fimol<-fimo_table_dir(poflong)
  colnames(fimol)[2]<-"alt_pval"
  fimoresult<-merge(fimos,fimol,by="id")
  fimoresult
}

fimo_table_dir<-function(x){
  id<-unique(x$sequence_name)
  pvalue<-mcmapply(function(id) findfimominp(x,id),id)
  data.frame(id,pvalue)
}

findfimominp<-function(x,i){
  min(x[x$sequence_name==i,]$p.value)
}

Fimo_p_MYC1<-FIMO_directp("Leukemia/FIMO based analysis/ref_MYC1.tsv","Leukemia/FIMO based analysis/alt_MYC1.tsv")
Fimo_p_MYC1$motif<-"MYC1"

Fimo_p_MYC2<-FIMO_directp("Leukemia/FIMO based analysis/ref_MYC2.tsv","Leukemia/FIMO based analysis/alt_MYC2.tsv")
Fimo_p_MYC2$motif<-"MYC2"

Fimo_p_MYC3<-FIMO_directp("Leukemia/FIMO based analysis/ref_MYC3.tsv","Leukemia/FIMO based analysis/alt_MYC3.tsv")
Fimo_p_MYC3$motif<-"MYC3"

FIMO_p<-rbind(Fimo_p_MYC1,Fimo_p_MYC2,Fimo_p_MYC3)
summary(FIMO_p$ref_pval)
summary(FIMO_p$alt_pval)
#all the p-values are smaller than 0.05.


#Direct q-value from FIMO
FIMO_directq<-function(p1,p2){
  pofshort<-data.table(read.table(p1, sep = '\t', header = TRUE))  
  poflong<-data.table(read.table(p2, sep = '\t', header = TRUE)) 
  fimos<-fimo_table_dirq(pofshort)
  colnames(fimos)[2]<-"ref_pval"
  fimol<-fimo_table_dirq(poflong)
  colnames(fimol)[2]<-"alt_pval"
  fimoresult<-merge(fimos,fimol,by="id")
  fimoresult
}

fimo_table_dirq<-function(x){
  id<-unique(x$sequence_name)
  pvalue<-mcmapply(function(id) findminq(x,id),id)
  data.frame(id,pvalue)
}

findminq<-function(x,i){
  min(x[x$sequence_name==i,]$q.value)
}


Fimo_q_MYC1<-FIMO_directq("Leukemia/FIMO based analysis/ref_MYC1.tsv","Leukemia/FIMO based analysis/alt_MYC1.tsv")
Fimo_q_MYC1$motif<-"MYC1"

Fimo_q_MYC2<-FIMO_directq("Leukemia/FIMO based analysis/ref_MYC2.tsv","Leukemia/FIMO based analysis/alt_MYC2.tsv")
Fimo_q_MYC2$motif<-"MYC2"

Fimo_q_MYC3<-FIMO_directq("Leukemia/FIMO based analysis/ref_MYC3.tsv","Leukemia/FIMO based analysis/alt_MYC3.tsv")
Fimo_q_MYC3$motif<-"MYC3"

FIMO_q<-rbind(Fimo_q_MYC1,Fimo_q_MYC2,Fimo_q_MYC3)
summary(FIMO_q$ref_pval)
summary(FIMO_q$alt_pval)
#all the q-values are larger than 0.05.


#adjusted p-value by sequence
FIMO_pval2<-function(p1,p2){
  pofshort<-data.table(read.table(p1, sep = '\t', header = TRUE))  
  pofshort<-pofshort[,pval_bh:=p.adjust(p.value,method="BH"),by=sequence_name]
  poflong<-data.table(read.table(p2, sep = '\t', header = TRUE)) 
  poflong<-poflong[,pval_bh:=p.adjust(p.value,method="BH"),by=sequence_name]
  fimos<-fimo_table_bh(pofshort)
  colnames(fimos)[2]<-"fimo_ref"
  fimol<-fimo_table_bh(poflong)
  colnames(fimol)[2]<-"fimo_alt"
  fimoresult<-merge(fimos,fimol,by="id")
  fimoresult
}


fimo_table_bh<-function(x){
  id<-unique(x$sequence_name)
  pvalue<-mcmapply(function(id) findminp_bh(x,id),id)
  data.frame(id,pvalue)
}

findminp_bh<-function(x,i){
  min(x[x$sequence_name==i,]$pval_bh)
}

Fimo_MYC1<-FIMO_pval2("Leukemia/FIMO based analysis/ref_MYC1.tsv","Leukemia/FIMO based analysis/alt_MYC1.tsv")
Fimo_MYC1$motif<-"MYC1"
Fimo_MYC2<-FIMO_pval2("Leukemia/FIMO based analysis/ref_MYC2.tsv","Leukemia/FIMO based analysis/alt_MYC2.tsv")
Fimo_MYC2$motif<-"MYC2"
Fimo_MYC3<-FIMO_pval2("Leukemia/FIMO based analysis/ref_MYC3.tsv","Leukemia/FIMO based analysis/alt_MYC3.tsv")
Fimo_MYC3$motif<-"MYC3"

FIMO_p<-rbind(Fimo_MYC1,Fimo_MYC2,Fimo_MYC3)
summary(FIMO_p$fimo_ref-FIMO_p$fimo_alt)

###FIMO cannot detect binding change/ no difference 










