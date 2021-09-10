#install.packages("data.table")
library(data.table)
#install.packages("devtools")
library(devtools)
#install_github("sunyoungshin/BindingChanger")
library(BindingChanger)

load("FO model parameters.Rdata")
load("motifs for Simulation.Rdata")
load("4.More models input.Rdata")

#Fifth order
Fifth_score<-atIndel::indel_motif_scores(motif_library,fifth_indel_info) 
Fifth_pvalue<-atIndel::indel_p_values(motif_library,fifth_indel_info,Fifth_score$list,prior,trans_mat,2000)

#summary table
at_result_table<-function(x){
  p<-c(0.01,0.05,0.1)
  change.rank<-mcmapply(function(p) emp(p,x$p_value_change.rank),p)
  data.table(p,change.rank)
}

emp<-function(p,x){
  length(which(x<p))/length(x)
}

tab91<-at_result_table(Fifth_pvalue[Fifth_pvalue$motif=="MSC",])
tab91

tab92<-at_result_table(Fifth_pvalue[Fifth_pvalue$motif=="Ddit3::Cebpa",])
tab92

tab93<-at_result_table(Fifth_pvalue[Fifth_pvalue$motif=="Hes1",])
tab93


#Zero order
ZO_score<-atIndel::indel_motif_scores(motif_library,Zero_indel_info) 
ZO_pvalue<-atIndel::indel_p_values(motif_library,Zero_indel_info,ZO_score$list,prior,trans_mat,2000)


#table
tab81<-at_result_table(ZO_pvalue[ZO_pvalue$motif=="MSC",])
tab81

tab82<-at_result_table(ZO_pvalue[ZO_pvalue$motif=="Ddit3::Cebpa",])
tab82

tab83<-at_result_table(ZO_pvalue[ZO_pvalue$motif=="Hes1",])
tab83

save(tab81,tab82,tab83,tab91,tab92,tab93,file="4.More models output.Rdata")
