library(data.table)
library(gridExtra)
library(atIndel)
library(markovchain)
library(ggplot2)
library(parallel)

#motif library
motif_library<-list()
motif_library$MSC<-v.motif$MSC
MSC<-v.motif$MSC
motif_library$`Ddit3::Cebpa`<-v.motif$`Ddit3::Cebpa`
motif_library$Hes1<-v.motif$Hes1


#generate sequence pairs 
mcDNA= new("markovchain", states = as.character(seq(4)),transitionMatrix = atIndel_tran, name = "DNA")

GenerateMotif<-function(motif,size){
  n<-nrow(motif)
  m<- matrix(nrow =size, ncol = n)
  for(i in 1:n){
    m[,i]<-t(sample(seq(4), size, replace=TRUE, prob = motif[i,]))
  }
  m
}

all_indel_info2<- list()
all_indel_info3<- list()
GenerateSeq<-function(n_indels,m,motif){
  motif<-GenerateMotif(motif,n_indels)
  width<-nrow(motif)
  all_indel_info<- list()
  for (i in seq_len(n_indels)) {
    insertion_len<-m
    t<-as.integer(markovchainSequence(n=width*2-3+insertion_len,
                                      mcDNA,t0 = sample(seq(4), 1,prob=atIndel_prior),include.t0 = TRUE, useRCpp = FALSE))
    
    t1<-t[1:(width-1)]
    t2<-t[(width+insertion_len):length(t)]
    t3<-motif[i,1:m]
    
    a1<-n2s(t1-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
    a2<-n2s(t2-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
    a3<-n2s(t3-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
    a4<-n2s(t-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
    aref<-n2s(t[(width-1)]-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
    
    ref<-paste(aref, collapse="")
    alt<-paste((c(aref,a3)), collapse="")
    
    all_indel_info[[i]] <- list(
      inserted_sequence = as.integer(c(t1,t3,t2)),
      insertion_len = insertion_len,
      insertion=1,
      ref=ref,
      alt=alt
    )
  }
  names(all_indel_info)<-c(1:n_indels)
  all_indel_info2<-all_indel_info
  names(all_indel_info2)<-c((n_indels*(m-1)+1):(n_indels*(m-1)+n_indels))
  all_indel_info3<-c(all_indel_info3,all_indel_info2)
  all_indel_info3
}

#empirical rejection probability
at_result_table<-function(x){
  p<-c(0.01,0.05,0.1)
  atIndel_p<-mcmapply(function(p) emp(p,x$p_value_change.rank),p)
  data.table(p,atIndel_p)
}

at_result_table_s<-function(x){
  p<-c(0.01,0.05,0.1)
  atIndel_p<-mcmapply(function(p) emp(p,x$p_value_change.score),p)
  data.table(p,atIndel_p)
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


#MSC
k1<-mapply(function(x) GenerateSeq(2000,x,motif_library[1][[1]]),1:10,SIMPLIFY = FALSE)
MSC_AL_indel_info<-do.call(c,k1)
MSC_AL_score<-atIndel::indel_motif_scores(motif_library[1],MSC_AL_indel_info)
MSC_AL_pval<-atIndel::indel_p_values(motif_library[1],MSC_AL_indel_info,MSC_AL_score$list,atIndel_prior,atIndel_tran,2000)
MSC_AL_pval$id<-as.numeric(MSC_AL_pval$id)
MSC_AL_pval2<-MSC_AL_pval[order(MSC_AL_pval$id),]

MSC_list <-mapply(function(x) getP(MSC_AL_pval2,x,2000),1:10,SIMPLIFY = FALSE)
MSCpa_tbl<-do.call(rbind.data.frame,MSC_list)
MSCpa_tbl$p<-factor(MSCpa_tbl$p)

ggplot(MSCpa_tbl, aes(x=m, y=atIndel_p, colour=p),size=0.1)+ggtitle("MSC")+ylab("Empirical Power") +ylim(0,1)+
  geom_line(size=0.1) +scale_x_discrete(name ="Inserted base length", limits=c(2,4,6,8,10))+
  geom_point()+theme_bw()+theme(title = element_text(size=14),
                                axis.text=element_text(size=14),axis.title=element_text(size=14),
                                legend.text = element_text(size=14),legend.title = element_text(size=14),
                                panel.border = element_blank(), panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#Ddit3::Cebpa
k1<-mapply(function(x) GenerateSeq(2000,x,motif_library[2]),1:12,SIMPLIFY = FALSE)
Ddit3_AL_indel_info<-do.call(c,k1)

Ddit3_AL_score<-atIndel::indel_motif_scores(motif_library[2],Ddit3_AL_indel_info)
Ddit3_AL_pval<-atIndel::indel_p_values(motif_library[2],Ddit3_AL_indel_info,Ddit3_AL_score$list,atIndel_prior,atIndel_tran,2000)
Ddit3_AL_pval$id<-as.numeric(Ddit3_AL_pval$id)
Ddit3_AL_pval2<-Ddit3_AL_pval[order(Ddit3_AL_pval$id),]

Ddit3_list <-mapply(function(x) getP(motif_library[2],Ddit3_AL_indel_info,x,5),1:12,SIMPLIFY = FALSE)
Ddit3pa_tbl<-do.call(rbind.data.frame,Ddit3_list)
Ddit3pa_tbl$p<-factor(Ddit3pa_tbl$p)

ggplot(Ddit3pa_tbl, aes(x=m, y=atIndel_p, colour=p),size=0.1)+ggtitle("Ddit3::Cebpa")+ylab("Empirical Power") +ylim(0,1)+
  geom_line(size=0.1) +scale_x_discrete(name ="Inserted base length", limits=c(2,4,6,8,10))+
  geom_point()+theme_bw()+theme(title = element_text(size=14),
                                axis.text=element_text(size=14),axis.title=element_text(size=14),
                                legend.text = element_text(size=14),legend.title = element_text(size=14),
                                panel.border = element_blank(), panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#Hes1
k1<-mapply(function(x) GenerateSeq(2000,x,motif_library[3]),1:10,SIMPLIFY = FALSE)
Hes1_AL_indel_info<-do.call(c,k1)

Hes1_AL_score<-atIndel::indel_motif_scores(motif_library[3],Hes1_AL_indel_info)
Hes1_AL_pval<-atIndel::indel_p_values(motif_library[3],Hes1_AL_indel_info,Hes1_AL_score$list,atIndel_prior,atIndel_tran,2000)
Hes1_AL_pval$id<-as.numeric(Hes1_AL_pval$id)
Hes1_AL_pval2<-Ddit3_AL_pval[order(Hes1_AL_pval$id),]

Hes1_list <-mapply(function(x) getP(2000,x,motif_library[3]),1:10,SIMPLIFY = FALSE)
Hes1pa_tbl<-do.call(rbind.data.frame,Hes1_list)
Hes1pa_tbl$p<-factor(Hes1pa_tbl$p)

ggplot(Hes1pa_tbl, aes(x=m, y=atIndel_p, colour=p),size=0.1)+ggtitle("Hes1")+ylab("Empirical Power") +ylim(0,1)+
  geom_line(size=0.1) +scale_x_discrete(name ="Inserted base length", limits=c(2,4,6,8,10))+
  geom_point()+theme_bw()+theme(title = element_text(size=14),
                                axis.text=element_text(size=14),axis.title=element_text(size=14),
                                legend.text = element_text(size=14),legend.title = element_text(size=14),
                                panel.border = element_blank(), panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#summary table
MSCpa_tbl[m==10,]
Ddit3pa_tbl[m==12,]
Hes1pa_tbl[m==10,]

#pval_score
MSC_list_s <-mapply(function(x) getP_s(MSC_AL_pval2,x,2000),1:10,SIMPLIFY = FALSE)
MSCpa_tbl_s<-do.call(rbind.data.frame,MSC_list_s)
MSCpa_tbl_s$p<-factor(MSCpa_tbl_s$p)

Ddit3_list_s <-mapply(function(x) getP_s(Ddit3_AL_pval2,x,2000),1:12,SIMPLIFY = FALSE)
Ddit3pa_tbl_s<-do.call(rbind.data.frame,Ddit3_list_s)
Ddit3pa_tbl_s$p<-factor(Ddit3pa_tbl_s$p)

Hes1_list_s <-mapply(function(x) getP_s(Hes1_AL_pval2,x,2000),1:10,SIMPLIFY = FALSE)
Hes1pa_tbl_s<-do.call(rbind.data.frame,Hes1_list_s)
Hes1pa_tbl_s$p<-factor(Hes1pa_tbl_s$p)


#summary table
MSCpa_tbl_s[m==10,]
Ddit3pa_tbl_s[m==12,]
Hes1pa_tbl_s[m==10,]
ggplot(MSCpa_tbl_s, aes(x=m, y=atIndel_p, colour=p),size=0.1)+ggtitle("MSC")+ylab("Empirical Power") +ylim(0,1)+
  geom_line(size=0.1) +scale_x_discrete(name ="Inserted base length", limits=c(2,4,6,8,10))+
  geom_point()+theme_bw()+theme(title = element_text(size=14),
                                axis.text=element_text(size=14),axis.title=element_text(size=14),
                                legend.text = element_text(size=14),legend.title = element_text(size=14),
                                panel.border = element_blank(), panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggplot(Ddit3pa_tbl_s, aes(x=m, y=atIndel_p, colour=p),size=0.1)+ggtitle("Ddit3::Cebpa")+ylab("Empirical Power") +ylim(0,1)+
  geom_line(size=0.1) +scale_x_discrete(name ="Inserted base length", limits=c(2,4,6,8,10))+
  geom_point()+theme_bw()+theme(title = element_text(size=14),
                                axis.text=element_text(size=14),axis.title=element_text(size=14),
                                legend.text = element_text(size=14),legend.title = element_text(size=14),
                                panel.border = element_blank(), panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggplot(Hes1pa_tbl_s, aes(x=m, y=atIndel_p, colour=p),size=0.1)+ggtitle("Hes1")+ylab("Empirical Power") +ylim(0,1)+
  geom_line(size=0.1) +scale_x_discrete(name ="Inserted base length", limits=c(2,4,6,8,10))+
  geom_point()+theme_bw()+theme(title = element_text(size=14),
                                axis.text=element_text(size=14),axis.title=element_text(size=14),
                                legend.text = element_text(size=14),legend.title = element_text(size=14),
                                panel.border = element_blank(), panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
