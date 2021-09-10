#install.packages("data.table")
library(data.table)
#install.packages("gridExtra")
library(gridExtra)
#install.packages("devtools")
library(devtools)
#install_github("sunyoungshin/atIndel")
library(atIndel)
#install.packages("markovchain")
library(markovchain)
#install.packages("ggplot2")
library(ggplot2)

load("FO model parameters.Rdata")
load("motifs for Simulation.Rdata")
load("2.AL input.Rdata")

#empirical rejection probability
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

getP<-function(x,m,n_indels){
  start<-n_indels*(m-1)+1
  end<-n_indels*m
  r<-x[start:end,]
  k<-at_result_table(r)
  k$m<-m
  k
}  

getP_s<-function(x,m,n_indels){
  start<-n_indels*(m-1)+1
  end<-n_indels*m
  r<-x[start:end,]
  k<-at_result_table_s(r)
  k$m<-m
  k
}  

#Tests conducted
#MSC
MSC_AL_score<-atIndel::indel_motif_scores(motif_library[1],MSC_AL_indel_info)
MSC_AL_pval<-atIndel::indel_p_values(motif_library[1],MSC_AL_indel_info,MSC_AL_score$list,prior,trans_mat,2000)
MSC_AL_pval$id<-as.numeric(MSC_AL_pval$id)
MSC_AL_pval2<-MSC_AL_pval[order(MSC_AL_pval$id),]


#Ddit3::Cebpa
Ddit3_AL_score<-atIndel::indel_motif_scores(motif_library[2],Ddit3_AL_indel_info)
Ddit3_AL_pval<-atIndel::indel_p_values(motif_library[2],Ddit3_AL_indel_info,Ddit3_AL_score$list,prior,trans_mat,2000)
Ddit3_AL_pval$id<-as.numeric(Ddit3_AL_pval$id)
Ddit3_AL_pval2<-Ddit3_AL_pval[order(Ddit3_AL_pval$id),]


#Hes1
Hes1_AL_score<-atIndel::indel_motif_scores(motif_library[3],Hes1_AL_indel_info)
Hes1_AL_pval<-atIndel::indel_p_values(motif_library[3],Hes1_AL_indel_info,Hes1_AL_score$list,prior,trans_mat,2000)
Hes1_AL_pval$id<-as.numeric(Hes1_AL_pval$id)
Hes1_AL_pval2<-Hes1_AL_pval[order(Hes1_AL_pval$id),]


#Summary-BC score
#MSC
MSC_list <-mapply(function(x) getP(MSC_AL_pval2,x,2000),1:10,SIMPLIFY = FALSE)
MSCpa_tbl<-do.call(rbind.data.frame,MSC_list)
MSCpa_tbl$p<-factor(MSCpa_tbl$p)

fig31<-ggplot(MSCpa_tbl, aes(x=m, y=change.rank, colour=p),size=0.1)+ggtitle("MSC")+ylab("Empirical Power") +ylim(0,1)+
  geom_line(size=0.1) +scale_x_discrete(name ="Inserted base length", limits=c(2,4,6,8,10))+
  geom_point()+theme_bw()+theme(title = element_text(size=14),
                                axis.text=element_text(size=14),axis.title=element_text(size=14),
                                legend.text = element_text(size=14),legend.title = element_text(size=14),
                                panel.border = element_blank(), panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#Ddit3::Cebpa
Ddit3_list <-mapply(function(x) getP(Ddit3_AL_pval2,x,2000),1:12,SIMPLIFY = FALSE)
Ddit3pa_tbl<-do.call(rbind.data.frame,Ddit3_list)
Ddit3pa_tbl$p<-factor(Ddit3pa_tbl$p)

fig32<-ggplot(Ddit3pa_tbl, aes(x=m, y=change.rank, colour=p),size=0.1)+ggtitle("Ddit3::Cebpa")+ylab("Empirical Power") +ylim(0,1)+
  geom_line(size=0.1) +scale_x_discrete(name ="Inserted base length", limits=c(2,4,6,8,10))+
  geom_point()+theme_bw()+theme(title = element_text(size=14),
                                axis.text=element_text(size=14),axis.title=element_text(size=14),
                                legend.text = element_text(size=14),legend.title = element_text(size=14),
                                panel.border = element_blank(), panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#Hes1
Hes1_list <-mapply(function(x) getP(Hes1_AL_pval2,x,2000),1:10,SIMPLIFY = FALSE)
Hes1pa_tbl<-do.call(rbind.data.frame,Hes1_list)
Hes1pa_tbl$p<-factor(Hes1pa_tbl$p)

fig33<-ggplot(Hes1pa_tbl, aes(x=m, y=change.rank, colour=p),size=0.1)+ggtitle("Hes1")+ylab("Empirical Power") +ylim(0,1)+
  geom_line(size=0.1) +scale_x_discrete(name ="Inserted base length", limits=c(2,4,6,8,10))+
  geom_point()+theme_bw()+theme(title = element_text(size=14),
                                axis.text=element_text(size=14),axis.title=element_text(size=14),
                                legend.text = element_text(size=14),legend.title = element_text(size=14),
                                panel.border = element_blank(), panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#summary table
tab31<-MSCpa_tbl[m==10,]
tab32<-Ddit3pa_tbl[m==12,]
tab33<-Hes1pa_tbl[m==10,]


#Summary-pval_score
#MSC
MSC_list_s <-mapply(function(x) getP_s(MSC_AL_pval2,x,2000),1:10,SIMPLIFY = FALSE)
MSCpa_tbl_s<-do.call(rbind.data.frame,MSC_list_s)
MSCpa_tbl_s$p<-factor(MSCpa_tbl_s$p)
fig91<-ggplot(MSCpa_tbl_s, aes(x=m, y=change.score, colour=p),size=0.1)+ggtitle("MSC")+ylab("Empirical Power") +ylim(0,1)+
  geom_line(size=0.1) +scale_x_discrete(name ="Inserted base length", limits=c(2,4,6,8,10))+
  geom_point()+theme_bw()+theme(title = element_text(size=14),
                                axis.text=element_text(size=14),axis.title=element_text(size=14),
                                legend.text = element_text(size=14),legend.title = element_text(size=14),
                                panel.border = element_blank(), panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#Ddit3::Cebpa
Ddit3_list_s <-mapply(function(x) getP_s(Ddit3_AL_pval2,x,2000),1:12,SIMPLIFY = FALSE)
Ddit3pa_tbl_s<-do.call(rbind.data.frame,Ddit3_list_s)
Ddit3pa_tbl_s$p<-factor(Ddit3pa_tbl_s$p)
fig92<-ggplot(Ddit3pa_tbl_s, aes(x=m, y=change.score, colour=p),size=0.1)+ggtitle("Ddit3::Cebpa")+ylab("Empirical Power") +ylim(0,1)+
  geom_line(size=0.1) +scale_x_discrete(name ="Inserted base length", limits=c(2,4,6,8,10))+
  geom_point()+theme_bw()+theme(title = element_text(size=14),
                                axis.text=element_text(size=14),axis.title=element_text(size=14),
                                legend.text = element_text(size=14),legend.title = element_text(size=14),
                                panel.border = element_blank(), panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#Hes1
Hes1_list_s <-mapply(function(x) getP_s(Hes1_AL_pval2,x,2000),1:10,SIMPLIFY = FALSE)
Hes1pa_tbl_s<-do.call(rbind.data.frame,Hes1_list_s)
Hes1pa_tbl_s$p<-factor(Hes1pa_tbl_s$p)
fig93<-ggplot(Hes1pa_tbl_s, aes(x=m, y=change.score, colour=p),size=0.1)+ggtitle("Hes1")+ylab("Empirical Power") +ylim(0,1)+
  geom_line(size=0.1) +scale_x_discrete(name ="Inserted base length", limits=c(2,4,6,8,10))+
  geom_point()+theme_bw()+theme(title = element_text(size=14),
                                axis.text=element_text(size=14),axis.title=element_text(size=14),
                                legend.text = element_text(size=14),legend.title = element_text(size=14),
                                panel.border = element_blank(), panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#summary table
tab71<-MSCpa_tbl_s[m==10,]
tab72<-Ddit3pa_tbl_s[m==12,]
tab73<-Hes1pa_tbl_s[m==10,]

save(fig31,fig32,fig33,fig91,fig92,fig93,tab31,tab32,tab33,tab71,tab72,tab73,file="2.AL output.Rdata")
