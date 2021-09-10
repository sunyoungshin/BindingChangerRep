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
#install.packages("motifStack")
library(motifStack)


load("FO model parameters.Rdata")
load("motifs for Simulation.Rdata")
load("1.Null-input.RData")

FO_score<-atIndel::indel_motif_scores(motif_library,fo_indel_info) 
FO_pvalue<-atIndel::indel_p_values(motif_library,fo_indel_info,FO_score$list,prior,trans_mat,2000)


#summary table
at_result_table<-function(x){
  p<-c(0.01,0.05,0.1)
  change.rank<-mcmapply(function(p) emp(p,x$p_value_change.rank),p)
  data.table(p,change.rank)
}

at_result_table_s<-function(x){
  p<-c(0.01,0.05,0.1)
  change.score<-mcmapply(function(p) emp(p,x$p_value_change.score),p)
  data.table(p,change.score)
}

emp<-function(p,x){
  length(which(x<p))/length(x)
}

tab21<-at_result_table(FO_pvalue[FO_pvalue$motif=="MSC",])

tab22<-at_result_table(FO_pvalue[FO_pvalue$motif=="Ddit3::Cebpa",])

tab23<-at_result_table(FO_pvalue[FO_pvalue$motif=="Hes1",])

tab61<-at_result_table_s(FO_pvalue[FO_pvalue$motif=="MSC",])

tab62<-at_result_table_s(FO_pvalue[FO_pvalue$motif=="Ddit3::Cebpa",])

tab63<-at_result_table_s(FO_pvalue[FO_pvalue$motif=="Hes1",])




#score difference
p81<-ggplot(FO_pvalue[FO_pvalue$motif=="MSC",],aes(x=ref.score-mutation.score))+ ggtitle("MSC")+ylab("Frequency") +
  ylim(0,2500)+xlab("Score difference")+ geom_histogram(color="black", fill="white",bins=50)+theme_bw()+
  theme(title = element_text(size=14),axis.text=element_text(size=14),axis.title=element_text(size=14),
        legend.text = element_text(size=14),legend.title = element_text(size=14),panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

p82<-ggplot(FO_pvalue[FO_pvalue$motif=="Ddit3::Cebpa",],aes(x=ref.score-mutation.score))+ ggtitle("Ddit3::Cebpa")+ylab("Frequency") +
  ylim(0,2500)+xlab("Score difference")+ geom_histogram(color="black", fill="white",bins=50)+theme_bw()+
  theme(title = element_text(size=14),axis.text=element_text(size=14),axis.title=element_text(size=14),
        legend.text = element_text(size=14),legend.title = element_text(size=14),panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

p83<-ggplot(FO_pvalue[FO_pvalue$motif=="Hes1",],aes(x=ref.score-mutation.score))+ ggtitle("Hes1")+ylab("Frequency") +
  ylim(0,2500)+xlab("Score difference")+ geom_histogram(color="black", fill="white",bins=50)+theme_bw()+
  theme(title = element_text(size=14),axis.text=element_text(size=14),axis.title=element_text(size=14),
        legend.text = element_text(size=14),legend.title = element_text(size=14),panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
fig8<-grid.arrange(p81, p82, p83,nrow = 1)

#log rank ratio
p21<-ggplot(FO_pvalue[FO_pvalue$motif=="MSC",],aes(x=log(ref.pval)-log(mutation.pval)))+ ggtitle("MSC")+ylab("Frequency") +
  ylim(0,4000)+xlab("Binding changer test statistic")+ geom_histogram(color="black", fill="white",bins=60)+theme_bw()+
  theme(title = element_text(size=14),axis.text=element_text(size=14),axis.title=element_text(size=14),
        legend.text = element_text(size=14),legend.title = element_text(size=14),panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

p22<-ggplot(FO_pvalue[FO_pvalue$motif=="Ddit3::Cebpa",],aes(x=log(ref.pval)-log(mutation.pval)))+ ggtitle("Ddit3::Cebpa")+ylab("Frequency") +
  ylim(0,4000)+xlim(-20,20)+xlab("Binding changer test statistic")+ geom_histogram(color="black", fill="white",binwidth =0.6)+theme_bw()+
  theme(title = element_text(size=14),axis.text=element_text(size=14),axis.title=element_text(size=14),
        legend.text = element_text(size=14),legend.title = element_text(size=14),panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

p23<-ggplot(FO_pvalue[FO_pvalue$motif=="Hes1",],aes(x=log(ref.pval)-log(mutation.pval)))+ ggtitle("Hes1")+ylab("Frequency") +
  ylim(0,4000)+xlim(-10,10)+xlab("Binding changer test statistic")+ geom_histogram(color="black", fill="white",bins=60)+theme_bw()+
  theme(title = element_text(size=14),axis.text=element_text(size=14),axis.title=element_text(size=14),
        legend.text = element_text(size=14),legend.title = element_text(size=14),panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
fig2<-grid.arrange(p21, p22, p23,nrow = 1)


p1<-t(motif_library$MSC)
rownames(p1)<-c("A","C","G","T")
fig21<-plotMotifLogo(p1,xlab="")


p2<-t(motif_library$`Ddit3::Cebpa`)
rownames(p2)<-c("A","C","G","T")
fig22<-plotMotifLogo(p2,xlab="")

p3<-t(motif_library$Hes1)
rownames(p3)<-c("A","C","G","T")
fig23<-plotMotifLogo(p3,xlab="")

save(tab21,tab22,tab23,tab61,tab62,tab63,p21,p22,p23,fig2,p81,p82,p83,fig8,file="1.Null-output.RData")
