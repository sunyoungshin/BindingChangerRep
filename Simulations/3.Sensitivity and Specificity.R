library(PRROC)
library(plotROC)
library(gridExtra)
library(ggplot2)
library(egg)

n_indels<-10000

mcDNA= new("markovchain", states = as.character(seq(4)),transitionMatrix = atIndel_tran, name = "DNA")

GenerateMotif<-function(motif,size){
  n<-nrow(motif)
  m<- matrix(, nrow =size, ncol = n)
  for(i in 1:n){
    m[,i]<-t(sample(seq(4), size, replace=TRUE, prob = motif[i,]))
  }
  m
}

Indel_prop<-function(p,n_indels,motif_library){
  short_seq<-list()
  long_seq<-list()
  all_indel_info<-list()
  
  width<-nrow(motif_library[[1]])
  insertion_len <-width
  ma<- vector()
  n1<-round(p*n_indels)
  motif<-GenerateMotif(motif_library[[1]],n1)
  
  for (i in seq_len(n1)) {
    t<-as.integer(markovchainSequence(n=width*2-3+insertion_len,
                                      mcDNA,t0 = sample(seq(4), 1,prob=atIndel_prior),include.t0 = TRUE, useRCpp = FALSE))
    t1<-t[1:(width-1)]
    t2<-t[(width+insertion_len):length(t)]
    t3<-motif[i,]
    
    a1<-n2s(t1-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
    a2<-n2s(t2-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
    a3<-n2s(t3-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
    a4<-n2s(t-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
    aref<-n2s(t[(width-1)]-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
    
    ref<-paste(aref, collapse="")
    alt<-paste((c(aref,a3)), collapse="")
    
    short_seq[[i]]<-paste((c(a1,a2)), collapse="")
    long_seq[[i]]<-paste((c(a1,a3,a2)), collapse="")
    
    all_indel_info[[i]] <- list(
      inserted_sequence = as.integer(c(t1,t3,t2)),
      insertion_len = insertion_len,
      insertion=1,
      ref=ref,
      alt=alt
    )
  }
  
  for (i in (n1+1):n_indels) {
    t<-as.integer(markovchainSequence(n=width * 2 - 3 + insertion_len,
                                      mcDNA,t0 = sample(seq(4), 1,prob=atIndel_prior),include.t0 = TRUE, useRCpp = FALSE))
    t1<-t[1:(width-1)]
    t2<-t[width:(width+insertion_len-1)]
    t3<-t[(width+insertion_len):length(t)]
    
    a1<-n2s(t1-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
    a2<-n2s(t2-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
    a3<-n2s(t3-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
    a4<-n2s(t-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
    aref<-n2s(t[(width-1)]-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
    
    short_seq[[i]]<-paste((c(a1,a3)), collapse="")
    long_seq[[i]]<-paste(a4, collapse="")
    
    ref<-paste(aref, collapse="")
    alt<-paste((c(aref,a2)), collapse="")
    
    all_indel_info[[i]] <- list(
      inserted_sequence = t,
      insertion_len = insertion_len,
      insertion=1,
      ref=ref,
      alt=alt
    )
  }
  
  # Y<-c(1:n_indels)
  # x1<-data.frame(Y,ldply (short_seq, data.frame))
  # x2<-data.frame(Y,ldply (long_seq, data.frame))
  # short.fasta = dataframe2fas(x1, file=p1)
  # long.fasta = dataframe2fas(x2, file=p2)
  # 
  names(all_indel_info)<-c(1:n_indels)
  
  s1<-atIndel::indel_motif_scores(motif_library,all_indel_info)
  r<-atIndel::indel_p_values(motif_library,all_indel_info,s1$list,atIndel_prior,atIndel_tran,2000)
  r$id<-as.numeric(r$id)
  rr<-r[order(r$id),] 
  rr
}

#MSC 0.005,10000
MSC_50_10000<-Indel_prop(0.005,10000,motif_library[1])
MSC_50_10000$rank_ratio<-abs(log(MSC_50_10000$p_value_affinity1)-log(MSC_50_10000$p_value_affinity2))
MSC_50_10000$label<-0
MSC_50_10000$label[1:50]<-1

FIMO_pval<-function(p1,p2){
  pofshort<-data.table(read.table(p1, sep = '\t', header = TRUE))  
  pofshort<-pofshort[,pval_bh:=p.adjust(p.value,method="BH"),by=sequence_name]
  poflong<-data.table(read.table(p2, sep = '\t', header = TRUE)) 
  poflong<-poflong[,pval_bh:=p.adjust(p.value,method="BH"),by=sequence_name]
  fimos<-fimo_table_bh(pofshort)
  colnames(fimos)[2]<-"s.pval_bh"
  fimol<-fimo_table_bh(poflong)
  colnames(fimol)[2]<-"l.pval_bh"
  fimoresult<-merge(fimos,fimol,by="id")
  fimoresult$pval_bh<-apply(fimoresult,1,function(x) min(x[2],x[3]))
  fimoresult
}

d1<-FIMO_pval("MSC/0.005 10000/short1.tsv","MSC/0.005 10000/long1.tsv")
d2<-FIMO_pval("MSC/0.005 10000/short2.tsv","MSC/0.005 10000/long2.tsv")
d3<-FIMO_pval("MSC/0.005 10000/short3.tsv","MSC/0.005 10000/long3.tsv")
d4<-FIMO_pval("MSC/0.005 10000/short4.tsv","MSC/0.005 10000/long4.tsv")
d5<-FIMO_pval("MSC/0.005 10000/short5.tsv","MSC/0.005 10000/long5.tsv")
fimo_MSC_50_10000<-rbind(d1,d2,d3,d4,d5)
fimo_MSC_50_10000$rank_ratio<-abs(log(fimo_MSC_50_10000$s.pval_bh)-log(fimo_MSC_50_10000$l.pval_bh))
fimo_MSC_50_10000$label<-0
fimo_MSC_50_10000$label[1:50]<-1

x1<-data.frame(MSC_50_10000$rank_ratio,MSC_50_10000$label,method="BC test")
colnames(x1)<-c("rank_ratio","label","method")
x2<-data.frame(fimo_MSC_50_10000$rank_ratio,fimo_MSC_50_10000$label,method="FIMO-based rank")
colnames(x2)<-c("rank_ratio","label","method")

x3<-rbind(x1,x2)

k1<-ggplot(x3, aes(m =rank_ratio, d = label,color=method))+geom_roc(n.cuts=0,labels=FALSE)+ggtitle("MSC (50/10,000)")+
  theme_bw()+ylab("Sensitivity")+xlab("1-Specificity")+theme(title = element_text(size=12),
                                                             axis.text=element_text(size=12),axis.title=element_text(size=12),
                                                             legend.text = element_text(size=12),legend.title = element_text(size=12),
                                                             panel.border = element_blank(), panel.grid.major = element_blank(),
                                                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text", x = 0.5, y = 0.5, label = "AUC: 0.9994",color="#F8766D",size=5)+
  annotate("text", x = 0.5, y = 0.4, label = "AUC: 0.9995",color="#00BFC4",size=5)
calc_auc(k1)
k1

#MSC 0.01,10000
MSC_100_10000<-Indel_prop(0.01,10000,motif_library[1])
MSC_100_10000$rank_ratio<-abs(log(MSC_100_10000$p_value_affinity1)-log(MSC_100_10000$p_value_affinity2))
MSC_100_10000$label<-0
MSC_100_10000$label[1:100]<-1

d1<-FIMO_pval("MSC/0.01 10000/short1.tsv","MSC/0.01 10000/long1.tsv")
d2<-FIMO_pval("MSC/0.01 10000/short2.tsv","MSC/0.01 10000/long2.tsv")
d3<-FIMO_pval("MSC/0.01 10000/short3.tsv","MSC/0.01 10000/long3.tsv")
d4<-FIMO_pval("MSC/0.01 10000/short4.tsv","MSC/0.01 10000/long4.tsv")
d5<-FIMO_pval("MSC/0.01 10000/short5.tsv","MSC/0.01 10000/long5.tsv")
fimo_MSC_100_10000<-rbind(d1,d2,d3,d4,d5)
fimo_MSC_100_10000$rank_ratio<-abs(log(fimo_MSC_100_10000$s.pval_bh)-log(fimo_MSC_100_10000$l.pval_bh))
fimo_MSC_100_10000$label<-0
fimo_MSC_100_10000$label[1:100]<-1

x1<-data.frame(MSC_100_10000$rank_ratio,MSC_100_10000$label,method="BC test")
colnames(x1)<-c("rank_ratio","label","method")
x2<-data.frame(fimo_MSC_100_10000$rank_ratio,fimo_MSC_100_10000$label,method="FIMO-based rank")
colnames(x2)<-c("rank_ratio","label","method")

k2<-ggplot(x3, aes(m =rank_ratio, d = label,color=method))+geom_roc(n.cuts=0,labels=FALSE)+ggtitle("MSC (100/10,000)")+
  theme_bw()+ylab("Sensitivity")+xlab("1-Specificity")+theme(title = element_text(size=12),
                                                             axis.text=element_text(size=12),axis.title=element_text(size=12),
                                                             legend.text = element_text(size=12),legend.title = element_text(size=12),
                                                             panel.border = element_blank(), panel.grid.major = element_blank(),
                                                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text", x = 0.5, y = 0.5, label = "AUC: 0.9985",color="#F8766D",size=5)+
  annotate("text", x = 0.5, y = 0.4, label = "AUC: 0.9982",color="#00BFC4",size=5)
calc_auc(k2)
k2

#MSC 0.05,10000
MSC_500_10000<-Indel_prop(0.05,10000,motif_library[1])
MSC_500_10000$rank_ratio<-abs(log(MSC_500_10000$p_value_affinity1)-log(MSC_500_10000$p_value_affinity2))
MSC_500_10000$label<-0
MSC_500_10000$label[1:500]<-1

d1<-FIMO_pval("MSC/0.05 10000/short1.tsv","MSC/0.05 10000/long1.tsv")
d2<-FIMO_pval("MSC/0.05 10000/short2.tsv","MSC/0.05 10000/long2.tsv")
d3<-FIMO_pval("MSC/0.05 10000/short3.tsv","MSC/0.05 10000/long3.tsv")
d4<-FIMO_pval("MSC/0.05 10000/short4.tsv","MSC/0.05 10000/long4.tsv")
d5<-FIMO_pval("MSC/0.05 10000/short5.tsv","MSC/0.05 10000/long5.tsv")
fimo_MSC_500_10000<-rbind(d1,d2,d3,d4,d5)
fimo_MSC_500_10000$rank_ratio<-abs(log(fimo_MSC_500_10000$s.pval_bh)-log(fimo_MSC_500_10000$l.pval_bh))
fimo_MSC_500_10000$label<-0
fimo_MSC_500_10000$label[0:500]<-1

x1<-data.frame(MSC_500_10000$rank_ratio,MSC_500_10000$label,method="BC test")
colnames(x1)<-c("rank_ratio","label","method")
x2<-data.frame(fimo_MSC_500_10000$rank_ratio,fimo_MSC_500_10000$label,method="FIMO-based rank")
colnames(x2)<-c("rank_ratio","label","method")

x3<-rbind(x1,x2)

k3<-ggplot(x3, aes(m =rank_ratio, d = label,color=method))+geom_roc(n.cuts=0,labels=FALSE)+ggtitle("MSC (500/10,000)")+
  theme_bw()+ylab("Sensitivity")+xlab("1-Specificity")+theme(title = element_text(size=12),
                                                             axis.text=element_text(size=12),axis.title=element_text(size=12),
                                                             legend.text = element_text(size=12),legend.title = element_text(size=12),
                                                             panel.border = element_blank(), panel.grid.major = element_blank(),
                                                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text", x = 0.5, y = 0.5, label = "AUC: 0.9936",color="#F8766D",size=5)+
  annotate("text", x = 0.5, y = 0.4, label = "AUC: 0.9977",color="#00BFC4",size=5)
calc_auc(k3)
k3



#Ddit3 0.005,10000
Ddit3_50_10000<-Indel_prop(0.005,10000,motif_library[2])
Ddit3_50_10000$rank_ratio<-abs(log(Ddit3_50_10000$p_value_affinity1)-log(Ddit3_50_10000$p_value_affinity2))
Ddit3_50_10000$label<-0
Ddit3_50_10000$label[1:50]<-1

d1<-FIMO_pval("Ddit3/0.005 10000/short1.tsv","Ddit3/0.005 10000/long1.tsv")
d2<-FIMO_pval("Ddit3/0.005 10000/short2.tsv","Ddit3/0.005 10000/long2.tsv")
d3<-FIMO_pval("Ddit3/0.005 10000/short3.tsv","Ddit3/0.005 10000/long3.tsv")
d4<-FIMO_pval("Ddit3/0.005 10000/short4.tsv","Ddit3/0.005 10000/long4.tsv")
d5<-FIMO_pval("Ddit3/0.005 10000/short5.tsv","Ddit3/0.005 10000/long5.tsv")
fimo_Ddit3_50_10000<-rbind(d1,d2,d3,d4,d5)
fimo_Ddit3_50_10000$rank_ratio<-abs(log(fimo_Ddit3_50_10000$s.pval_bh)-log(fimo_Ddit3_50_10000$l.pval_bh))
fimo_Ddit3_50_10000$label<-0
fimo_Ddit3_50_10000$label[1:50]<-1

x1<-data.frame(Ddit3_50_10000$rank_ratio,Ddit3_50_10000$label,method="BC test")
colnames(x1)<-c("rank_ratio","label","method")
x2<-data.frame(fimo_Ddit3_50_10000$rank_ratio,fimo_Ddit3_50_10000$label,method="FIMO-based rank")
colnames(x2)<-c("rank_ratio","label","method")

x3<-rbind(x1,x2)

j1<-ggplot(x3, aes(m =rank_ratio, d = label,color=method))+geom_roc(n.cuts=0,labels=FALSE)+ggtitle("Ddit3::Cebpa (50/10,000)")+
  theme_bw()+ylab("Sensitivity")+xlab("1-Specificity")+theme(title = element_text(size=12),
                                                             axis.text=element_text(size=12),axis.title=element_text(size=12),
                                                             legend.text = element_text(size=12),legend.title = element_text(size=12),
                                                             panel.border = element_blank(), panel.grid.major = element_blank(),
                                                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text", x = 0.5, y = 0.5, label = "AUC: 0.9630",color="#F8766D",size=5)+
  annotate("text", x = 0.5, y = 0.4, label = "AUC: 0.9598",color="#00BFC4",size=5)
calc_auc(j1)
j1

#Ddit3 0.01,10000
Ddit3_100_10000<-Indel_prop(0.01,10000,motif_library[2])
Ddit3_100_10000$rank_ratio<-abs(log(Ddit3_100_10000$p_value_affinity1)-log(Ddit3_100_10000$p_value_affinity2))
Ddit3_100_10000$label<-0
Ddit3_100_10000$label[1:100]<-1

d1<-FIMO_pval("Ddit3/0.01 10000/short1.tsv","Ddit3/0.01 10000/long1.tsv")
d2<-FIMO_pval("Ddit3/0.01 10000/short2.tsv","Ddit3/0.01 10000/long2.tsv")
d3<-FIMO_pval("Ddit3/0.01 10000/short3.tsv","Ddit3/0.01 10000/long3.tsv")
d4<-FIMO_pval("Ddit3/0.01 10000/short4.tsv","Ddit3/0.01 10000/long4.tsv")
d5<-FIMO_pval("Ddit3/0.01 10000/short5.tsv","Ddit3/0.01 10000/long5.tsv")
fimo_Ddit3_100_10000<-rbind(d1,d2,d3,d4,d5)
fimo_Ddit3_100_10000$rank_ratio<-abs(log(fimo_Ddit3_100_10000$s.pval_bh)-log(fimo_Ddit3_100_10000$l.pval_bh))
fimo_Ddit3_100_10000$label<-0
fimo_Ddit3_100_10000$label[1:100]<-1

x1<-data.frame(Ddit3_100_10000$rank_ratio,Ddit3_100_10000$label,method="BC test")
colnames(x1)<-c("rank_ratio","label","method")
x2<-data.frame(fimo_Ddit3_100_10000$rank_ratio,fimo_Ddit3_100_10000$label,method="FIMO-based rank")
colnames(x2)<-c("rank_ratio","label","method")

x3<-rbind(x1,x2)

j2<-ggplot(x3, aes(m =rank_ratio, d = label,color=method))+geom_roc(n.cuts=0,labels=FALSE)+ggtitle("Ddit3::Cebpa (100/10,000)")+
  theme_bw()+ylab("Sensitivity")+xlab("1-Specificity")+theme(title = element_text(size=12),
                                                             axis.text=element_text(size=12),axis.title=element_text(size=12),
                                                             legend.text = element_text(size=12),legend.title = element_text(size=12),
                                                             panel.border = element_blank(), panel.grid.major = element_blank(),
                                                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text", x = 0.5, y = 0.5, label = "AUC: 0.9551",color="#F8766D",size=5)+
  annotate("text", x = 0.5, y = 0.4, label = "AUC: 0.9524",color="#00BFC4",size=5)
calc_auc(j2)
j2

#Ddit3 0.05,10000
Ddit3_500_10000<-Indel_prop(0.05,10000,motif_library[2])
Ddit3_500_10000$rank_ratio<-abs(log(Ddit3_500_10000$p_value_affinity1)-log(Ddit3_500_10000$p_value_affinity2))
Ddit3_500_10000$label<-0
Ddit3_500_10000$label[1:500]<-1

d1<-FIMO_pval("Ddit3/0.05 10000/short1.tsv","Ddit3/0.05 10000/long1.tsv")
d2<-FIMO_pval("Ddit3/0.05 10000/short2.tsv","Ddit3/0.05 10000/long2.tsv")
d3<-FIMO_pval("Ddit3/0.05 10000/short3.tsv","Ddit3/0.05 10000/long3.tsv")
d4<-FIMO_pval("Ddit3/0.05 10000/short4.tsv","Ddit3/0.05 10000/long4.tsv")
d5<-FIMO_pval("Ddit3/0.05 10000/short5.tsv","Ddit3/0.05 10000/long5.tsv")
fimo_Ddit3_500_10000<-rbind(d1,d2,d3,d4,d5)
fimo_Ddit3_500_10000$rank_ratio<-abs(log(fimo_Ddit3_500_10000$s.pval_bh)-log(fimo_Ddit3_500_10000$l.pval_bh))
fimo_Ddit3_500_10000$label<-0
fimo_Ddit3_500_10000$label[1:500]<-1

x1<-data.frame(Ddit3_500_10000$rank_ratio,Ddit3_500_10000$label,method="BC test")
colnames(x1)<-c("rank_ratio","label","method")
x2<-data.frame(fimo_Ddit3_500_10000$rank_ratio,fimo_Ddit3_500_10000$label,method="FIMO-based rank")
colnames(x2)<-c("rank_ratio","label","method")

x3<-rbind(x1,x2)

j3<-ggplot(x3, aes(m =rank_ratio, d = label,color=method))+geom_roc(n.cuts=0,labels=FALSE)+ggtitle("Ddit3::Cebpa (500/10,000)")+
  theme_bw()+ylab("Sensitivity")+xlab("1-Specificity")+theme(title = element_text(size=12),
                                                             axis.text=element_text(size=12),axis.title=element_text(size=12),
                                                             legend.text = element_text(size=12),legend.title = element_text(size=12),
                                                             panel.border = element_blank(), panel.grid.major = element_blank(),
                                                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text", x = 0.5, y = 0.5, label = "AUC: 0.9676",color="#F8766D",size=5)+
  annotate("text", x = 0.5, y = 0.4, label = "AUC: 0.9742",color="#00BFC4",size=5)
calc_auc(j3)
j3

#Hes1 0.005,10000
Hes1_50_10000<-Indel_prop(0.005,10000,motif_library[3])
Hes1_50_10000$rank_ratio<-abs(log(Hes1_50_10000$p_value_affinity1)-log(Hes1_50_10000$p_value_affinity2))
Hes1_50_10000$label<-0
Hes1_50_10000$label[1:50]<-1

d1<-FIMO_pval("Hes1/0.005 10000/short1.tsv","Hes1/0.005 10000/long1.tsv")
d2<-FIMO_pval("Hes1/0.005 10000/short2.tsv","Hes1/0.005 10000/long2.tsv")
d3<-FIMO_pval("Hes1/0.005 10000/short3.tsv","Hes1/0.005 10000/long3.tsv")
d4<-FIMO_pval("Hes1/0.005 10000/short4.tsv","Hes1/0.005 10000/long4.tsv")
d5<-FIMO_pval("Hes1/0.005 10000/short5.tsv","Hes1/0.005 10000/long5.tsv")
fimo_Hes1_50_10000<-rbind(d1,d2,d3,d4,d5)
fimo_Hes1_50_10000$rank_ratio<-abs(log(fimo_Hes1_50_10000$s.pval_bh)-log(fimo_Hes1_50_10000$l.pval_bh))
fimo_Hes1_50_10000$label<-0
fimo_Hes1_50_10000$label[1:50]<-1

x1<-data.frame(Hes1_50_10000$rank_ratio,Hes1_50_10000$label,method="BC test")
colnames(x1)<-c("rank_ratio","label","method")
x2<-data.frame(fimo_Hes1_50_10000$rank_ratio,fimo_Hes1_50_10000$label,method="FIMO-based rank")
colnames(x2)<-c("rank_ratio","label","method")

x3<-rbind(x1,x2)

h1<-ggplot(x3, aes(m =rank_ratio, d = label,color=method))+geom_roc(n.cuts=0,labels=FALSE)+ggtitle("Hes1 (50/10,000)")+
  theme_bw()+ylab("Sensitivity")+xlab("1-Specificity")+theme(title = element_text(size=12),
                                                             axis.text=element_text(size=12),axis.title=element_text(size=12),
                                                             legend.text = element_text(size=12),legend.title = element_text(size=12),
                                                             panel.border = element_blank(), panel.grid.major = element_blank(),
                                                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text", x = 0.5, y = 0.5, label = "AUC: 0.8771",color="#F8766D",size=5)+
  annotate("text", x = 0.5, y = 0.4, label = "AUC: 0.8854",color="#00BFC4",size=5)
calc_auc(h1)
h1

#Hes1 0.01,10000
Hes1_100_10000<-Indel_prop(0.01,10000,motif_library[3])
Hes1_100_10000$rank_ratio<-abs(log(Hes1_100_10000$p_value_affinity1)-log(Hes1_100_10000$p_value_affinity2))
Hes1_100_10000$label<-0
Hes1_100_10000$label[1:100]<-1

d1<-FIMO_pval("Hes1/0.01 10000/short1.tsv","Hes1/0.01 10000/long1.tsv")
d2<-FIMO_pval("Hes1/0.01 10000/short2.tsv","Hes1/0.01 10000/long2.tsv")
d3<-FIMO_pval("Hes1/0.01 10000/short3.tsv","Hes1/0.01 10000/long3.tsv")
d4<-FIMO_pval("Hes1/0.01 10000/short4.tsv","Hes1/0.01 10000/long4.tsv")
d5<-FIMO_pval("Hes1/0.01 10000/short5.tsv","Hes1/0.01 10000/long5.tsv")
fimo_Hes1_100_10000<-rbind(d1,d2,d3,d4,d5)
fimo_Hes1_100_10000$rank_ratio<-abs(log(fimo_Hes1_100_10000$s.pval_bh)-log(fimo_Hes1_100_10000$l.pval_bh))
fimo_Hes1_100_10000$label<-0
fimo_Hes1_100_10000$label[1:100]<-1

x1<-data.frame(Hes1_100_10000$rank_ratio,Hes1_100_10000$label,method="BC test")
colnames(x1)<-c("rank_ratio","label","method")
x2<-data.frame(fimo_Hes1_100_10000$rank_ratio,fimo_Hes1_100_10000$label,method="FIMO-based rank")
colnames(x2)<-c("rank_ratio","label","method")

x3<-rbind(x1,x2)

h2<-ggplot(x3, aes(m =rank_ratio, d = label,color=method))+geom_roc(n.cuts=0,labels=FALSE)+ggtitle("Hes1 (100/10,000)")+
  theme_bw()+ylab("Sensitivity")+xlab("1-Specificity")+theme(title = element_text(size=12),
                                                             axis.text=element_text(size=12),axis.title=element_text(size=12),
                                                             legend.text = element_text(size=12),legend.title = element_text(size=12),
                                                             panel.border = element_blank(), panel.grid.major = element_blank(),
                                                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text", x = 0.5, y = 0.5, label = "AUC: 0.8298",color="#F8766D",size=5)+
  annotate("text", x = 0.5, y = 0.4, label = "AUC: 0.8487",color="#00BFC4",size=5)
calc_auc(h2)
h2

#Hes1 0.05,10000
Hes1_500_10000<-Indel_prop(0.05,10000,motif_library[3])
Hes1_500_10000$rank_ratio<-abs(log(Hes1_500_10000$p_value_affinity1)-log(Hes1_500_10000$p_value_affinity2))
Hes1_500_10000$label<-0
Hes1_500_10000$label[1:500]<-1

d1<-FIMO_pval("Hes1/0.05 10000/short1.tsv","Hes1/0.05 10000/long1.tsv")
d2<-FIMO_pval("Hes1/0.05 10000/short2.tsv","Hes1/0.05 10000/long2.tsv")
d3<-FIMO_pval("Hes1/0.05 10000/short3.tsv","Hes1/0.05 10000/long3.tsv")
d4<-FIMO_pval("Hes1/0.05 10000/short4.tsv","Hes1/0.05 10000/long4.tsv")
d5<-FIMO_pval("Hes1/0.05 10000/short5.tsv","Hes1/0.05 10000/long5.tsv")
fimo_Hes1_500_10000<-rbind(d1,d2,d3,d4,d5)
fimo_Hes1_500_10000$rank_ratio<-abs(log(fimo_Hes1_500_10000$s.pval_bh)-log(fimo_Hes1_500_10000$l.pval_bh))
fimo_Hes1_500_10000$label<-0
fimo_Hes1_500_10000$label[1:500]<-1

x1<-data.frame(Hes1_500_10000$rank_ratio,Hes1_a500_10000$label,method="BC test")
colnames(x1)<-c("rank_ratio","label","method")
x2<-data.frame(fimo_Hes1_500_10000$rank_ratio,fimo_Hes1_500_10000$label,method="FIMO-based rank")
colnames(x2)<-c("rank_ratio","label","method")

x3<-rbind(x1,x2)

h3<-ggplot(x3, aes(m =rank_ratio, d = label,color=method))+geom_roc(n.cuts=0,labels=FALSE)+ggtitle("Hes1 (500/10,000)")+
  theme_bw()+ylab("Sensitivity")+xlab("1-Specificity")+theme(title = element_text(size=12),
                                                             axis.text=element_text(size=12),axis.title=element_text(size=12),
                                                             legend.text = element_text(size=12),legend.title = element_text(size=12),
                                                             panel.border = element_blank(), panel.grid.major = element_blank(),
                                                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text", x = 0.5, y = 0.5, label = "AUC: 0.8394",color="#F8766D",size=5)+
  annotate("text", x = 0.5, y = 0.4, label = "AUC: 0.8586",color="#00BFC4",size=5)
calc_auc(h3)
h3

