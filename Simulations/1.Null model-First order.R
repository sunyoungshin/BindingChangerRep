library(data.table)
library(gridExtra)
library(atIndel)
library(markovchain)
library(ggplot2)
library(parallel)

load("simulation.Rdata")
#First order
fo_indel_info<- list()
short_seq<-list()
long_seq<-list()
n_indels<-10000

mcDNA= new("markovchain", states = as.character(seq(4)),transitionMatrix = atIndel_tran, name = "DNA")
max_width<-max(sapply(motif_library, nrow))
insertion_len <-6

#random seq generation
for (i in seq_len(n_indels)) {
  
  t<-as.integer(markovchainSequence(n=max_width * 2 - 3 + insertion_len,
                                    mcDNA,t0 = sample(seq(4), 1,prob=atIndel_prior),include.t0 = TRUE, useRCpp = FALSE))
  t1<-t[1:(max_width-1)]
  t2<-t[max_width:(max_width+insertion_len-1)]
  t3<-t[(max_width+insertion_len):length(t)]
  
  a1<-n2s(t1-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
  a2<-n2s(t2-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
  a3<-n2s(t3-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
  a4<-n2s(t-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
  aref<-n2s(t[(max_width-1)]-1, levels = c("A", "C", "G", "T"), base4 = TRUE)
  
  ref<-paste(aref, collapse="")
  alt<-paste((c(aref,a2)), collapse="")
  
  fo_indel_info[[i]] <- list(
    inserted_sequence = t,
    insertion_len = insertion_len,
    insertion=1,
    ref=ref,
    alt=alt
  )
}

names(fo_indel_info)<-c(1:n_indels)


FO_score<-atIndel::indel_motif_scores(motif_library,fo_indel_info) 
FO_pvalue<-atIndel::indel_p_values(motif_library,fo_indel_info,FO_score$list,atIndel_prior,atIndel_tran,2000)



#summary table
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

FO_1<-at_result_table(FO_pvalue[FO_pvalue$motif=="MSC",])
FO_1

FO_2<-at_result_table(FO_pvalue[FO_pvalue$motif=="Ddit3::Cebpa",])
FO_2

FO_3<-at_result_table(FO_pvalue[FO_pvalue$motif=="Hes1",])
FO_3


#score difference
p1<-ggplot(FO_pvalue[FO_pvalue$motif=="MSC",],aes(x=ref.score-mutation.score))+ ggtitle("MSC")+ylab("Frequency") +
  ylim(0,2500)+xlab("Score difference")+ geom_histogram(color="black", fill="white",bins=50)+theme_bw()+
  theme(title = element_text(size=14),axis.text=element_text(size=14),axis.title=element_text(size=14),
        legend.text = element_text(size=14),legend.title = element_text(size=14),panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

p2<-ggplot(FO_pvalue[FO_pvalue$motif=="Ddit3::Cebpa",],aes(x=ref.score-mutation.score))+ ggtitle("Ddit3::Cebpa")+ylab("Frequency") +
  ylim(0,2500)+xlab("Score difference")+ geom_histogram(color="black", fill="white",bins=50)+theme_bw()+
  theme(title = element_text(size=14),axis.text=element_text(size=14),axis.title=element_text(size=14),
        legend.text = element_text(size=14),legend.title = element_text(size=14),panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

p3<-ggplot(FO_pvalue[FO_pvalue$motif=="Hes1",],aes(x=ref.score-mutation.score))+ ggtitle("Hes1")+ylab("Frequency") +
  ylim(0,2500)+xlab("Score difference")+ geom_histogram(color="black", fill="white",bins=50)+theme_bw()+
  theme(title = element_text(size=14),axis.text=element_text(size=14),axis.title=element_text(size=14),
        legend.text = element_text(size=14),legend.title = element_text(size=14),panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
grid.arrange(p1, p2, p3,nrow = 1)

#log rank ratio
p1<-ggplot(FO_pvalue[FO_pvalue$motif=="MSC",],aes(x=log(ref.pval)-log(mutation.pval)))+ ggtitle("MSC")+ylab("Frequency") +
  ylim(0,4000)+xlab("Binding changer test statistic")+ geom_histogram(color="black", fill="white",bins=60)+theme_bw()+
  theme(title = element_text(size=14),axis.text=element_text(size=14),axis.title=element_text(size=14),
        legend.text = element_text(size=14),legend.title = element_text(size=14),panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

p2<-ggplot(FO_pvalue[FO_pvalue$motif=="Ddit3::Cebpa",],aes(x=log(ref.pval)-log(mutation.pval)))+ ggtitle("Ddit3::Cebpa")+ylab("Frequency") +
  ylim(0,4000)+xlim(-20,20)+xlab("Binding changer test statistic")+ geom_histogram(color="black", fill="white",binwidth =0.6)+theme_bw()+
  theme(title = element_text(size=14),axis.text=element_text(size=14),axis.title=element_text(size=14),
        legend.text = element_text(size=14),legend.title = element_text(size=14),panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

p3<-ggplot(FO_pvalue[FO_pvalue$motif=="Hes1",],aes(x=log(ref.pval)-log(mutation.pval)))+ ggtitle("Hes1")+ylab("Frequency") +
  ylim(0,4000)+xlim(-10,10)+xlab("Binding changer test statistic")+ geom_histogram(color="black", fill="white",bins=60)+theme_bw()+
  theme(title = element_text(size=14),axis.text=element_text(size=14),axis.title=element_text(size=14),
        legend.text = element_text(size=14),legend.title = element_text(size=14),panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
grid.arrange(p1, p2, p3,nrow = 1)

