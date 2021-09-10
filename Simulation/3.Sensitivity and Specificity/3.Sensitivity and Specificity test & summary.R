#install.packages("data.table")
library(data.table)
#install.packages("gridExtra")
library(gridExtra)
#install.packages("devtools")
library(devtools)
#install_github("sunyoungshin/atIndel")
library(atIndel)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("plotROC")
library(plotROC)
#install.packages("egg")
library(egg)

load("3.Sensitivity and Specificity input.Rdata")
load("FO model parameters.Rdata")
load("motifs for Simulation.Rdata")

BCtest<-function(motif_library,all_indel_info){
s<-atIndel::indel_motif_scores(motif_library,all_indel_info)
r<-atIndel::indel_p_values(motif_library,all_indel_info,s$list,prior,trans_mat,2000)
r$id<-as.numeric(r$id)
rr<-r[order(r$id),] 
rr
}

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

fimo_table_bh<-function(x){
  id<-unique(x$sequence_name)
  pvalue<-mcmapply(function(id) findminp_bh(x,id),id)
  data.frame(id,pvalue)
}

findminp_bh<-function(x,i){
  min(x[x$sequence_name==i,]$pval_bh)
}

#MSC 0.005,10000
MSC_50_10000<-BCtest(motif_library[1],MSC_50_10000_indel_info)
MSC_50_10000$rank_ratio<-abs(log(MSC_50_10000$ref.pval)-log(MSC_50_10000$mutation.pval))
MSC_50_10000$label<-0
MSC_50_10000$label[1:50]<-1


d1<-FIMO_pval("MSC/0.005 10000/short1.tsv","MSC/0.005 10000/long1.tsv")
d2<-FIMO_pval("MSC/0.005 10000/short2.tsv","MSC/0.005 10000/long2.tsv")
d3<-FIMO_pval("MSC/0.005 10000/short3.tsv","MSC/0.005 10000/long3.tsv")
d4<-FIMO_pval("MSC/0.005 10000/short4.tsv","MSC/0.005 10000/long4.tsv")
fimo_MSC_50_10000<-rbind(d1,d2,d3,d4)
fimo_MSC_50_10000$rank_ratio<-abs(log(fimo_MSC_50_10000$s.pval_bh)-log(fimo_MSC_50_10000$l.pval_bh))
fimo_MSC_50_10000$label<-0
fimo_MSC_50_10000$label[1:50]<-1

x1<-data.frame(MSC_50_10000$rank_ratio,MSC_50_10000$label,method="BC test")
colnames(x1)<-c("rank_ratio","label","method")
x2<-data.frame(fimo_MSC_50_10000$rank_ratio,fimo_MSC_50_10000$label,method="FIMO-based rank")
colnames(x2)<-c("rank_ratio","label","method")

MSC50_10000<-rbind(x1,x2)

fig411<-ggplot(MSC50_10000, aes(m =rank_ratio, d = label,color=method))+geom_roc(n.cuts=0,labels=FALSE)+ggtitle("MSC (50/10,000)")+
  theme_bw()+ylab("Sensitivity")+xlab("1-Specificity")+theme(title = element_text(size=12),
                                                             axis.text=element_text(size=12),axis.title=element_text(size=12),
                                                             legend.text = element_text(size=12),legend.title = element_text(size=12),
                                                             panel.border = element_blank(), panel.grid.major = element_blank(),
                                                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text", x = 0.5, y = 0.5, label = "AUC: 0.9994",color="#F8766D",size=5)+
  annotate("text", x = 0.5, y = 0.4, label = "AUC: 0.9995",color="#00BFC4",size=5)
calc_auc(fig411)
fig411

#MSC 0.01,10000
MSC_100_10000<-BCtest(motif_library[1],MSC_100_10000_indel_info)
MSC_100_10000$rank_ratio<-abs(log(MSC_100_10000$ref.pval)-log(MSC_100_10000$mutation.pval))
MSC_100_10000$label<-0
MSC_100_10000$label[1:100]<-1

d1<-FIMO_pval("MSC/0.01 10000/short1.tsv","MSC/0.01 10000/long1.tsv")
d2<-FIMO_pval("MSC/0.01 10000/short2.tsv","MSC/0.01 10000/long2.tsv")
d3<-FIMO_pval("MSC/0.01 10000/short3.tsv","MSC/0.01 10000/long3.tsv")
d4<-FIMO_pval("MSC/0.01 10000/short4.tsv","MSC/0.01 10000/long4.tsv")
fimo_MSC_100_10000<-rbind(d1,d2,d3,d4)
fimo_MSC_100_10000$rank_ratio<-abs(log(fimo_MSC_100_10000$s.pval_bh)-log(fimo_MSC_100_10000$l.pval_bh))
fimo_MSC_100_10000$label<-0
fimo_MSC_100_10000$label[1:100]<-1

x1<-data.frame(MSC_100_10000$rank_ratio,MSC_100_10000$label,method="BC test")
colnames(x1)<-c("rank_ratio","label","method")
x2<-data.frame(fimo_MSC_100_10000$rank_ratio,fimo_MSC_100_10000$label,method="FIMO-based rank")
colnames(x2)<-c("rank_ratio","label","method")

MSC100_10000<-rbind(x1,x2)

fig421<-ggplot(MSC100_10000, aes(m =rank_ratio, d = label,color=method))+geom_roc(n.cuts=0,labels=FALSE)+ggtitle("MSC (100/10,000)")+
  theme_bw()+ylab("Sensitivity")+xlab("1-Specificity")+theme(title = element_text(size=12),
                                                             axis.text=element_text(size=12),axis.title=element_text(size=12),
                                                             legend.text = element_text(size=12),legend.title = element_text(size=12),
                                                             panel.border = element_blank(), panel.grid.major = element_blank(),
                                                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text", x = 0.5, y = 0.5, label = "AUC: 0.9980",color="#F8766D",size=5)+
  annotate("text", x = 0.5, y = 0.4, label = "AUC: 0.9982",color="#00BFC4",size=5)
calc_auc(fig421)
fig421


#MSC 0.05,10000
MSC_500_10000<-BCtest(motif_library[1],MSC_500_10000_indel_info)
MSC_500_10000$rank_ratio<-abs(log(MSC_500_10000$ref.pval)-log(MSC_500_10000$mutation.pval))
MSC_500_10000$label<-0
MSC_500_10000$label[1:500]<-1

d1<-FIMO_pval("MSC/0.05 10000/short1.tsv","MSC/0.05 10000/long1.tsv")
d2<-FIMO_pval("MSC/0.05 10000/short2.tsv","MSC/0.05 10000/long2.tsv")
d3<-FIMO_pval("MSC/0.05 10000/short3.tsv","MSC/0.05 10000/long3.tsv")
d4<-FIMO_pval("MSC/0.05 10000/short4.tsv","MSC/0.05 10000/long4.tsv")

fimo_MSC_500_10000<-rbind(d1,d2,d3,d4)
fimo_MSC_500_10000$rank_ratio<-abs(log(fimo_MSC_500_10000$s.pval_bh)-log(fimo_MSC_500_10000$l.pval_bh))
fimo_MSC_500_10000$label<-0
fimo_MSC_500_10000$label[1:500]<-1

x1<-data.frame(MSC_500_10000$rank_ratio,MSC_500_10000$label,method="BC test")
colnames(x1)<-c("rank_ratio","label","method")
x2<-data.frame(fimo_MSC_500_10000$rank_ratio,fimo_MSC_500_10000$label,method="FIMO-based rank")
colnames(x2)<-c("rank_ratio","label","method")

MSC500_10000<-rbind(x1,x2)

fig431<-ggplot(MSC500_10000, aes(m =rank_ratio, d = label,color=method))+geom_roc(n.cuts=0,labels=FALSE)+ggtitle("MSC (500/10,000)")+
  theme_bw()+ylab("Sensitivity")+xlab("1-Specificity")+theme(title = element_text(size=12),
                                                             axis.text=element_text(size=12),axis.title=element_text(size=12),
                                                             legend.text = element_text(size=12),legend.title = element_text(size=12),
                                                             panel.border = element_blank(), panel.grid.major = element_blank(),
                                                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text", x = 0.5, y = 0.5, label = "AUC: 0.9951",color="#F8766D",size=5)+
  annotate("text", x = 0.5, y = 0.4, label = "AUC: 0.9977",color="#00BFC4",size=5)
calc_auc(fig431)
fig431

#Ddit3 0.005,10000
Ddit3_50_10000<-BCtest(motif_library[2],Ddit3_50_10000_indel_info)
Ddit3_50_10000$rank_ratio<-abs(log(Ddit3_50_10000$ref.pval)-log(Ddit3_50_10000$mutation.pval))
Ddit3_50_10000$label<-0
Ddit3_50_10000$label[1:50]<-1

d1<-FIMO_pval("Ddit3/0.005 10000/short1.tsv","Ddit3/0.005 10000/long1.tsv")
d2<-FIMO_pval("Ddit3/0.005 10000/short2.tsv","Ddit3/0.005 10000/long2.tsv")
d3<-FIMO_pval("Ddit3/0.005 10000/short3.tsv","Ddit3/0.005 10000/long3.tsv")
d4<-FIMO_pval("Ddit3/0.005 10000/short4.tsv","Ddit3/0.005 10000/long4.tsv")
fimo_Ddit3_50_10000<-rbind(d1,d2,d3,d4)
fimo_Ddit3_50_10000$rank_ratio<-abs(log(fimo_Ddit3_50_10000$s.pval_bh)-log(fimo_Ddit3_50_10000$l.pval_bh))
fimo_Ddit3_50_10000$label<-0
fimo_Ddit3_50_10000$label[1:50]<-1

x1<-data.frame(Ddit3_50_10000$rank_ratio,Ddit3_50_10000$label,method="BC test")
colnames(x1)<-c("rank_ratio","label","method")
x2<-data.frame(fimo_Ddit3_50_10000$rank_ratio,fimo_Ddit3_50_10000$label,method="FIMO-based rank")
colnames(x2)<-c("rank_ratio","label","method")

Ddit350_10000<-rbind(x1,x2)

fig412<-ggplot(Ddit350_10000, aes(m =rank_ratio, d = label,color=method))+geom_roc(n.cuts=0,labels=FALSE)+ggtitle("Ddit3::Cebpa (50/10,000)")+
  theme_bw()+ylab("Sensitivity")+xlab("1-Specificity")+theme(title = element_text(size=12),
                                                             axis.text=element_text(size=12),axis.title=element_text(size=12),
                                                             legend.text = element_text(size=12),legend.title = element_text(size=12),
                                                             panel.border = element_blank(), panel.grid.major = element_blank(),
                                                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text", x = 0.5, y = 0.5, label = "AUC: 0.9560",color="#F8766D",size=5)+
  annotate("text", x = 0.5, y = 0.4, label = "AUC: 0.9764",color="#00BFC4",size=5)
calc_auc(fig412)
fig412

#Ddit3 0.01,10000
Ddit3_100_10000<-BCtest(motif_library[2],Ddit3_100_10000_indel_info)
Ddit3_100_10000$rank_ratio<-abs(log(Ddit3_100_10000$ref.pval)-log(Ddit3_100_10000$mutation.pval))
Ddit3_100_10000$label<-0
Ddit3_100_10000$label[1:100]<-1

d1<-FIMO_pval("Ddit3/0.01 10000/short1.tsv","Ddit3/0.01 10000/long1.tsv")
d2<-FIMO_pval("Ddit3/0.01 10000/short2.tsv","Ddit3/0.01 10000/long2.tsv")
d3<-FIMO_pval("Ddit3/0.01 10000/short3.tsv","Ddit3/0.01 10000/long3.tsv")
d4<-FIMO_pval("Ddit3/0.01 10000/short4.tsv","Ddit3/0.01 10000/long4.tsv")
fimo_Ddit3_100_10000<-rbind(d1,d2,d3,d4)
fimo_Ddit3_100_10000$rank_ratio<-abs(log(fimo_Ddit3_100_10000$s.pval_bh)-log(fimo_Ddit3_100_10000$l.pval_bh))
fimo_Ddit3_100_10000$label<-0
fimo_Ddit3_100_10000$label[1:100]<-1

x1<-data.frame(Ddit3_100_10000$rank_ratio,Ddit3_100_10000$label,method="BC test")
colnames(x1)<-c("rank_ratio","label","method")
x2<-data.frame(fimo_Ddit3_100_10000$rank_ratio,fimo_Ddit3_100_10000$label,method="FIMO-based rank")
colnames(x2)<-c("rank_ratio","label","method")

Ddit3100_10000<-rbind(x1,x2)

fig422<-ggplot(Ddit3100_10000<-rbind(x1,x2), aes(m =rank_ratio, d = label,color=method))+geom_roc(n.cuts=0,labels=FALSE)+ggtitle("Ddit3::Cebpa (100/10,000)")+
  theme_bw()+ylab("Sensitivity")+xlab("1-Specificity")+theme(title = element_text(size=12),
                                                             axis.text=element_text(size=12),axis.title=element_text(size=12),
                                                             legend.text = element_text(size=12),legend.title = element_text(size=12),
                                                             panel.border = element_blank(), panel.grid.major = element_blank(),
                                                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text", x = 0.5, y = 0.5, label = "AUC: 0.9553",color="#F8766D",size=5)+
  annotate("text", x = 0.5, y = 0.4, label = "AUC: 0.9651",color="#00BFC4",size=5)
calc_auc(fig422)
fig422


#Ddit3 0.05,10000
Ddit3_500_10000<-BCtest(motif_library[2],Ddit3_500_10000_indel_info)
Ddit3_500_10000$rank_ratio<-abs(log(Ddit3_500_10000$ref.pval)-log(Ddit3_500_10000$mutation.pval))
Ddit3_500_10000$label<-0
Ddit3_500_10000$label[1:500]<-1

d1<-FIMO_pval("Ddit3/0.05 10000/short1.tsv","Ddit3/0.05 10000/long1.tsv")
d2<-FIMO_pval("Ddit3/0.05 10000/short2.tsv","Ddit3/0.05 10000/long2.tsv")
d3<-FIMO_pval("Ddit3/0.05 10000/short3.tsv","Ddit3/0.05 10000/long3.tsv")
d4<-FIMO_pval("Ddit3/0.05 10000/short4.tsv","Ddit3/0.05 10000/long4.tsv")
fimo_Ddit3_500_10000<-rbind(d1,d2,d3,d4)
fimo_Ddit3_500_10000$rank_ratio<-abs(log(fimo_Ddit3_500_10000$s.pval_bh)-log(fimo_Ddit3_500_10000$l.pval_bh))
fimo_Ddit3_500_10000$label<-0
fimo_Ddit3_500_10000$label[1:500]<-1

x1<-data.frame(Ddit3_500_10000$rank_ratio,Ddit3_500_10000$label,method="BC test")
colnames(x1)<-c("rank_ratio","label","method")
x2<-data.frame(fimo_Ddit3_500_10000$rank_ratio,fimo_Ddit3_500_10000$label,method="FIMO-based rank")
colnames(x2)<-c("rank_ratio","label","method")

Ddit3500_10000<-rbind(x1,x2)

fig432<-ggplot(Ddit3500_10000, aes(m =rank_ratio, d = label,color=method))+geom_roc(n.cuts=0,labels=FALSE)+ggtitle("Ddit3::Cebpa (500/10,000)")+
  theme_bw()+ylab("Sensitivity")+xlab("1-Specificity")+theme(title = element_text(size=12),
                                                             axis.text=element_text(size=12),axis.title=element_text(size=12),
                                                             legend.text = element_text(size=12),legend.title = element_text(size=12),
                                                             panel.border = element_blank(), panel.grid.major = element_blank(),
                                                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text", x = 0.5, y = 0.5, label = "AUC: 0.9647",color="#F8766D",size=5)+
  annotate("text", x = 0.5, y = 0.4, label = "AUC: 0.9816",color="#00BFC4",size=5)
calc_auc(fig432)
fig432


#Hes1 0.005,10000
Hes1_50_10000<-BCtest(motif_library[3],Hes1_50_10000_indel_info)
Hes1_50_10000$rank_ratio<-abs(log(Hes1_50_10000$ref.pval)-log(Hes1_50_10000$mutation.pval))
Hes1_50_10000$label<-0
Hes1_50_10000$label[1:50]<-1

d1<-FIMO_pval("Hes1/0.005 10000/short1.tsv","Hes1/0.005 10000/long1.tsv")
d2<-FIMO_pval("Hes1/0.005 10000/short2.tsv","Hes1/0.005 10000/long2.tsv")
d3<-FIMO_pval("Hes1/0.005 10000/short3.tsv","Hes1/0.005 10000/long3.tsv")
d4<-FIMO_pval("Hes1/0.005 10000/short4.tsv","Hes1/0.005 10000/long4.tsv")

fimo_Hes1_50_10000<-rbind(d1,d2,d3,d4)
fimo_Hes1_50_10000$rank_ratio<-abs(log(fimo_Hes1_50_10000$s.pval_bh)-log(fimo_Hes1_50_10000$l.pval_bh))
fimo_Hes1_50_10000$label<-0
fimo_Hes1_50_10000$label[1:50]<-1

x1<-data.frame(Hes1_50_10000$rank_ratio,Hes1_50_10000$label,method="BC test")
colnames(x1)<-c("rank_ratio","label","method")
x2<-data.frame(fimo_Hes1_50_10000$rank_ratio,fimo_Hes1_50_10000$label,method="FIMO-based rank")
colnames(x2)<-c("rank_ratio","label","method")

Hes150_10000<-rbind(x1,x2)

fig413<-ggplot(Hes150_10000, aes(m =rank_ratio, d = label,color=method))+geom_roc(n.cuts=0,labels=FALSE)+ggtitle("Hes1 (50/10,000)")+
  theme_bw()+ylab("Sensitivity")+xlab("1-Specificity")+theme(title = element_text(size=12),
                                                             axis.text=element_text(size=12),axis.title=element_text(size=12),
                                                             legend.text = element_text(size=12),legend.title = element_text(size=12),
                                                             panel.border = element_blank(), panel.grid.major = element_blank(),
                                                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text", x = 0.5, y = 0.5, label = "AUC: 0.8804",color="#F8766D",size=5)+
  annotate("text", x = 0.5, y = 0.4, label = "AUC: 0.8854",color="#00BFC4",size=5)
calc_auc(fig413)
fig413

#Hes1 0.01,10000
Hes1_100_10000<-BCtest(motif_library[3],Hes1_100_10000_indel_info)
Hes1_100_10000$rank_ratio<-abs(log(Hes1_100_10000$ref.pval)-log(Hes1_100_10000$mutation.pval))
Hes1_100_10000$label<-0
Hes1_100_10000$label[1:100]<-1

d1<-FIMO_pval("Hes1/0.01 10000/short1.tsv","Hes1/0.01 10000/long1.tsv")
d2<-FIMO_pval("Hes1/0.01 10000/short2.tsv","Hes1/0.01 10000/long2.tsv")
d3<-FIMO_pval("Hes1/0.01 10000/short3.tsv","Hes1/0.01 10000/long3.tsv")
d4<-FIMO_pval("Hes1/0.01 10000/short4.tsv","Hes1/0.01 10000/long4.tsv")

fimo_Hes1_100_10000<-rbind(d1,d2,d3,d4)
fimo_Hes1_100_10000$rank_ratio<-abs(log(fimo_Hes1_100_10000$s.pval_bh)-log(fimo_Hes1_100_10000$l.pval_bh))
fimo_Hes1_100_10000$label<-0
fimo_Hes1_100_10000$label[1:100]<-1

x1<-data.frame(Hes1_100_10000$rank_ratio,Hes1_100_10000$label,method="BC test")
colnames(x1)<-c("rank_ratio","label","method")
x2<-data.frame(fimo_Hes1_100_10000$rank_ratio,fimo_Hes1_100_10000$label,method="FIMO-based rank")
colnames(x2)<-c("rank_ratio","label","method")

Hes1100_10000<-rbind(x1,x2)

fig423<-ggplot(Hes1100_10000, aes(m =rank_ratio, d = label,color=method))+geom_roc(n.cuts=0,labels=FALSE)+ggtitle("Hes1 (100/10,000)")+
  theme_bw()+ylab("Sensitivity")+xlab("1-Specificity")+theme(title = element_text(size=12),
                                                             axis.text=element_text(size=12),axis.title=element_text(size=12),
                                                             legend.text = element_text(size=12),legend.title = element_text(size=12),
                                                             panel.border = element_blank(), panel.grid.major = element_blank(),
                                                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text", x = 0.5, y = 0.5, label = "AUC: 0.8300",color="#F8766D",size=5)+
  annotate("text", x = 0.5, y = 0.4, label = "AUC: 0.8487",color="#00BFC4",size=5)
calc_auc(fig423)
fig423


#Hes1 0.05,10000
Hes1_500_10000<-BCtest(motif_library[3],Hes1_500_10000_indel_info)
Hes1_500_10000$rank_ratio<-abs(log(Hes1_500_10000$ref.pval)-log(Hes1_500_10000$mutation.pval))
Hes1_500_10000$label<-0
Hes1_500_10000$label[1:500]<-1

d1<-FIMO_pval("Hes1/0.05 10000/short1.tsv","Hes1/0.05 10000/long1.tsv")
d2<-FIMO_pval("Hes1/0.05 10000/short2.tsv","Hes1/0.05 10000/long2.tsv")
d3<-FIMO_pval("Hes1/0.05 10000/short3.tsv","Hes1/0.05 10000/long3.tsv")
d4<-FIMO_pval("Hes1/0.05 10000/short4.tsv","Hes1/0.05 10000/long4.tsv")
fimo_Hes1_500_10000<-rbind(d1,d2,d3,d4)
fimo_Hes1_500_10000$rank_ratio<-abs(log(fimo_Hes1_500_10000$s.pval_bh)-log(fimo_Hes1_500_10000$l.pval_bh))
fimo_Hes1_500_10000$label<-0
fimo_Hes1_500_10000$label[1:500]<-1

x1<-data.frame(Hes1_500_10000$rank_ratio,Hes1_500_10000$label,method="BC test")
colnames(x1)<-c("rank_ratio","label","method")
x2<-data.frame(fimo_Hes1_500_10000$rank_ratio,fimo_Hes1_500_10000$label,method="FIMO-based rank")
colnames(x2)<-c("rank_ratio","label","method")

Hes1500_10000<-rbind(x1,x2)

fig433<-ggplot(Hes1500_10000, aes(m =rank_ratio, d = label,color=method))+geom_roc(n.cuts=0,labels=FALSE)+ggtitle("Hes1 (500/10,000)")+
  theme_bw()+ylab("Sensitivity")+xlab("1-Specificity")+theme(title = element_text(size=12),
                                                             axis.text=element_text(size=12),axis.title=element_text(size=12),
                                                             legend.text = element_text(size=12),legend.title = element_text(size=12),
                                                             panel.border = element_blank(), panel.grid.major = element_blank(),
                                                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text", x = 0.5, y = 0.5, label = "AUC: 0.8404",color="#F8766D",size=5)+
  annotate("text", x = 0.5, y = 0.4, label = "AUC: 0.8586",color="#00BFC4",size=5)
calc_auc(fig433)
fig433


save(fig411,fig412,fig413,fig421,fig422,fig423,fig431,fig432,fig433,file="3.Sensitivity and Specificity output.Rdata")

