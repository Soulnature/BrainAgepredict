library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(dplyr)
library(igraph)
library(GOSemSim)
library(data.table)
library(ggpubr)
load("E:/yang/bag_gwas/results/新建文件夹/WGCNA_data_brain_adult.RData")
bag_gene<-read.table('F:/复旦/brainAge/GWAS/10.17/bag_all_gene.txt')
#bag_gene<-bag_gene$V1
## PPI network ##
TRN_network_data<-read.csv('TRN.csv')
co_expression_network<-fread('gene_cor_expression_no_filter.txt')##
#
write.csv(gene_cor_expression,'gene_cor_expression.csv',row.names = FALSE)
##
#
PPi_intersection<-read.table('PPI_intersection.txt',header = TRUE)
PPi_intersection_score_data<-PPi_intersection[which(PPi_intersection$combined_score>=700),]
###
PPi_intersection_score_data<-PPi_intersection_score_data[,c(1,2,10)]

PPI_info<-fread('./9606.protein.info.v11.5.txt',header = TRUE,sep = '\t')
PPI_info<-PPI_info[which(!duplicated(PPI_info$preferred_name)),c(1:2)]
PPI_info<-as.data.frame(PPI_info)
rownames(PPI_info)<-PPI_info$string_protein_id
###
PPi_intersection_score_data$protein1<-PPI_info[PPi_intersection_score_data$protein1,2]
PPi_intersection_score_data$protein2<-PPI_info[PPi_intersection_score_data$protein2,2]
###
all_protein<-union(PPi_intersection_score_data$protein1,PPi_intersection_score_data$protein2)
###

get_gene_brain<-intersect(all_protein,rownames(WGCNA_adult_data))
####
adult_gene_exp_mean<-apply(WGCNA_adult_data, 1, mean)%>%as.data.frame()
adult_gene_exp_mean<-data.frame(symbol=rownames(adult_gene_exp_mean),exp=adult_gene_exp_mean$.)
rownames(adult_gene_exp_mean)<-adult_gene_exp_mean$symbol
protein_exp<-adult_gene_exp_mean[get_gene_brain,]
protein_exp<-protein_exp[which(protein_exp$exp>0),]  ## 8769 protein #
###
protein_brain_exp<-protein_exp$symbol%>%unique()

##
####  Data  ##
row_seq<-c()
for (i in c(1:nrow(PPi_intersection_score_data))) {
  
  inter_gene<-intersect(union(PPi_intersection_score_data[i,1],PPi_intersection_score_data[i,2]),protein_brain_exp)
  if(length(inter_gene)>1){
    row_seq<-c(row_seq,i)
  }
}
PPi_intersection_score_data<-PPi_intersection_score_data[row_seq,]##summary the STRING
## brain specific TRN PPI and corexpression ##
write.csv(PPi_intersection_score_data,"PPI_network.csv",row.names = FALSE,quote = FALSE,sep = '\t')
## summary the networt information #


#########    ##############
get_split_gene<-function(network,gene_list){
  PPi_intersection_save_label<-apply(network, 1, function(x){
    tem_se<-union(x[1],x[2])
    if(length(intersect(tem_se,gene_list))>0){ 
      return(1)
    }
    else{
      return(0)
    }
  })
  PPi_intersection_label_data<-cbind(network,PPi_intersection_save_label)
  PPi_intersection_label_chosing<-PPi_intersection_label_data[which(PPi_intersection_label_data$PPi_intersection_save_label==1),]
  return(PPi_intersection_label_chosing)
}
bag_PPi<-get_split_gene(PPi_intersection_score_data,bag_gene)
bag_TRN<-get_split_gene(TRN_network_data,bag_gene)
bag_coexp<-get_split_gene(co_expression_network,bag_gene)

#disorder PPI
##
disorder_table<-read.csv('F:/复旦/brainAge/GWAS/disease_gene.csv',fill = TRUE,header = TRUE)
##
get_disorder_gene<-function(disorder_table,disorder_names){
  disorder_gene<-disorder_table[,disorder_names]
  disorder_gene<-disorder_gene[which(disorder_gene!="")]
  return(disorder_gene)
}

hyper_test<-function(disorder_net,bag_net,background){
  #disorder_net<-ADHD
  disorder_net<-disorder_net%>%as.data.frame()
  bag_net<-bag_net%>%as.data.frame()
  ##
  #background<-background_gene
  disorder_node<-union(disorder_net[,1],disorder_net[,2])%>%unique()
  bag_node<-union(bag_net[,1],bag_net[,2])%>%unique()
  intersec_bag_d<-intersect(disorder_node,bag_node)
  #########
  q<-intersec_bag_d
  m<-setdiff(disorder_node,bag_node)%>%unique()
  n<-setdiff(background,m)%>%unique()
  k<-setdiff(bag_node,disorder_node)%>%unique()
  Hyper_test<-phyper(length(q)-1,length(m),length(n),length(k),lower.tail = F)
  print(length(q))
  ## compute the jacard distance #
  jacard_dis<-(length(intersect(disorder_node,bag_node)))/length(union(disorder_node,bag_node))
  print(jacard_dis)
  save_hyper_jaccard<-data.frame(Hyper_test,jacard_dis)
  others_fisher<-setdiff(background,union(disorder_node,bag_node))%>%length()
  fisher_matrix<-rbind(c(length(q),length(m)),c(length(k),others_fisher))
  fisher_test<-fisher.test(fisher_matrix)
  return(data.frame(jacard=save_hyper_jaccard,fisher=fisher_test$p.value))
  ### fisher test ##
}
##### disease gene #
test_func_d<-function(Net_work_art,bag_net,back_ground_gene){
  colnames(disorder_table)
  ADHD<-get_split_gene(Net_work_art,get_disorder_gene(disorder_table,colnames(disorder_table)[1]))
  AD<-get_split_gene(Net_work_art,get_disorder_gene(disorder_table,colnames(disorder_table)[2]))
  MDD<-get_split_gene(Net_work_art,get_disorder_gene(disorder_table,colnames(disorder_table)[3]))
  SCZ<-get_split_gene(Net_work_art,get_disorder_gene(disorder_table,colnames(disorder_table)[4]))
  PD<-get_split_gene(Net_work_art,get_disorder_gene(disorder_table,colnames(disorder_table)[5]))
  ASD<-get_split_gene(Net_work_art,get_disorder_gene(disorder_table,colnames(disorder_table)[6]))
  BIP<-get_split_gene(Net_work_art,get_disorder_gene(disorder_table,colnames(disorder_table)[7]))
  ## Hypergeometric test ## PPi_intersection
  hyper_ADHD<-hyper_test(ADHD,bag_net,back_ground_gene)
  hyper_AD<-hyper_test(AD,bag_net,back_ground_gene)
  hyper_MDD<-hyper_test(MDD,bag_net,back_ground_gene)
  hyper_SCZ<-hyper_test(SCZ,bag_net,back_ground_gene)
  hyper_PD<-hyper_test(PD,bag_net,back_ground_gene)
  hyper_ASD<-hyper_test(ASD,bag_net,back_ground_gene)
  hyper_BIP<-hyper_test(BIP,bag_net,back_ground_gene)
  sta_values<-rbind(hyper_ADHD,hyper_AD,hyper_MDD,hyper_SCZ,hyper_PD,hyper_ASD,hyper_BIP)
  ## combined the enriched p values in the data #
  re_df_co<-data.frame(disorder_name=c('ADHD','AD','MDD','SCZ','PD','ASD','BIP'),statistic_value=sta_values)
  return(re_df_co)
}
background_gene<-rownames(ying_new_data) 
###
PPI_res<-test_func_d(PPi_intersection_score_data,bag_PPi,background_gene)
##
TRN_res<-test_func_d(TRN_network_data,bag_TRN,background_gene)
co_exp_res<-test_func_d(gene_cor_expression,bag_coexp,background_gene)
##plot the enrichment result ##
PPI_res$statistic_value.fisher<-p.adjust(PPI_res$statistic_value.fisher,method = 'bonferroni')
TRN_res$statistic_value.fisher<-p.adjust(TRN_res$statistic_value.fisher,method = 'bonferroni')
co_exp_res$statistic_value.fisher<-p.adjust(co_exp_res$statistic_value.fisher,method = 'bonferroni')
network_label<-c(rep('PPI',nrow(PPI_res)),rep('TRN',nrow(TRN_res)),rep('coexpression',nrow(co_exp_res)))
enrichment_plot<-data.frame(rbind(PPI_res,TRN_res,co_exp_res),group=network_label)
#
enrichment_plot$group<-factor(enrichment_plot$group,levels = c('TRN','PPI','coexpression'))
enrichment_plot$disorder_name<-factor(enrichment_plot,levels = c('',))
#enrichment_plot$statistic_value.jacard.jacard_dis<-factor(enrichment_plot$statistic_value.jacard.jacard_dis,levels = enrichment_plot$statistic_value.jacard.jacard_dis)
ggbarplot(enrichment_plot,x='group',y='statistic_value.jacard.jacard_dis',fill ='disorder_name',
          position =position_dodge())+
  geom_hline(aes(yintercept=mean(PPI_res$statistic_value.jacard.jacard_dis)), colour="#990000", linetype="dashed")+
  geom_hline(aes(yintercept=mean(TRN_res$statistic_value.jacard.jacard_dis)), colour="#990000", linetype="dashed")+
  geom_hline(aes(yintercept=mean(co_exp_res$statistic_value.jacard.jacard_dis)), colour="#990000", linetype="dashed")
####
mean(PPI_res$statistic_value.jacard.jacard_dis)
mean(PPI_res$statistic_value.jacard.jacard_dis)
mean(PPI_res$statistic_value.jacard.jacard_dis)
graph2ppt(file='network_enrichment',width=8,height=6,append=TRUE)



###
## all test was significant  test th#
PPI_res<-PPI_res[order(PPI_res$statistic_value.jacard_dis),]
TRN_res<-TRN_res[order(TRN_res$statistic_value.jacard_dis),]
##

PPI_res[which(PPI_res$statistic_value.Hyper_test==0),2]<-2.2e-16
TRN_res$enrichment_data<-fiser_test
ggbarplot(PPI_res,x='disorder_name',y='statistic_value.jacard_dis',fill='#7F7F7F')
ggbarplot(TRN_res,x='disorder_name',y='statistic_value.jacard_dis',fill='#7F7F7F')
##
fiser_test<-c(0.922,4.78e-10,2.22e-4,2.2e-16,1.2e-06,2.2e-16,5.545e-06)
disorder_d<-c('ADHD','AD','MDD','SCZ','PD','ASD','BIP')
disorder_df<-data.frame(Name=disorder_d,enrichment=fiser_test)
disorder_df$enrichment<-p.adjust(disorder_df$enrichment,method='fdr')
disorder_df$enrichment<--log10(disorder_df$enrichment)
disorder_df<-disorder_df[order(disorder_df$enrichment),]
ggbarplot(disorder_df,x='Name',y='enrichment',fill='#7F7F7F')

###
###
PPi_intersection_label_chosing<-PPi_intersection_label_chosing[which(PPi_intersection_label_chosing$combined_score>700),]
##
PPi_intersection_label_chosing<-PPi_intersection_label_chosing[!duplicated(PPi_intersection_label_chosing),] ##
## ENS id to gene symbol #
rownames(PPI_info)<-PPI_info$string_protein_id
protein1_ens<-PPI_info[PPi_intersection_label_chosing$protein1,2]
protein2_ens<-PPI_info[PPi_intersection_label_chosing$protein2,2]
PPi_intersection_label_chosing$protein1<-protein1_ens
PPi_intersection_label_chosing$protein2<- protein2_ens
##


colnames(disorder_table)
### physical interaction score > 700, mapping the disorder gene in the PPI intersection ##
#

###
mappping_disorder_gene_PPI<-function(network_PPI,disorder_name){
  #network_PPI<-gene_cor_expression
 # disorder_name<-disorder_name[1]
  disorder_gene<-disorder_table[,disorder_name]
  rownames(network_PPI)<-c(1:nrow(network_PPI))
  disorder_gene<-disorder_gene[which(disorder_gene!="")]
  ##map disorder gene in
  inter_gene<-intersect(bag_gene$V1,disorder_gene)
  ## mapping the disorder gene in BAG PPI interaction #
 # intersect(disorder_gene,union(PPi_intersection_label_chosing$protein1,PPi_intersection_label_chosing$protein2))
  multi_gene<-union(disorder_gene,bag_gene$V1)
  ##
 # te_g<-union(network_PPI$start,network_PPI$end_node)%>%unique()
  #intersect(te_g,multi_gene)
  disorder_label<-apply(network_PPI, 1, function(x){
  #  x<-network_PPI[2,]
   # all_gene_data<-union(x[,1],x[,2])
    all_gene_data<-union(x[1],x[2])
    #print(all_gene_data)
    if(length(intersect(all_gene_data,multi_gene))>1){
      x<-1
    }else{
      x<-0
    }
  })
  disorder_subnetwork<-cbind(network_PPI,disorder_label)
  disorder_subnetwork<-disorder_subnetwork[which(disorder_subnetwork$disorder_label==1),]
  ######write the node annotation and the network files ##
  write.table(disorder_subnetwork[,c(1,2,3)],paste(disorder_name,'net.txt',sep = ""),sep = '\t',row.names = FALSE,quote = FALSE)
  ##annotate the node ## the 1 is bag gene, the 2 is the disorder gene #
  all_gene_name<-union(disorder_subnetwork$node1,disorder_subnetwork$node2)%>%unique()
  #####
  bag_gene_ld<-intersect(all_gene_name,bag_gene$V1)
  disorder_ge_l<-intersect(all_gene_name,disorder_gene)
  ###
  bag_gene_ld<-setdiff(bag_gene_ld,inter_gene)
  disorder_ge_l<-setdiff(disorder_ge_l,inter_gene)
  ###
  label_data_gene<-c(rep(1,length(bag_gene_ld)),rep(2,length(disorder_ge_l)),rep(3,length(inter_gene)))
  gene_comb<-c(bag_gene_ld,disorder_ge_l,inter_gene)
  node_ano<-data.frame(gene_name=gene_comb,label_data=label_data_gene)
  colnames(node_ano)<-c('gene',paste(disorder_name,'_anno',sep = ""))
  #####
  write.table(node_ano,paste(disorder_name,'node.txt',sep = ""),sep = '\t',row.names = FALSE,quote = FALSE)

  ##
  return(disorder_subnetwork)
} 
######### PPI  network ##
disorder_name<-colnames(disorder_table)
ADHD<-mappping_disorder_gene_PPI(PPi_intersection_label_chosing,disorder_name[1])
AD<-mappping_disorder_gene_PPI(PPi_intersection_label_chosing,disorder_name[2])
MDD<-mappping_disorder_gene_PPI(PPi_intersection_label_chosing,disorder_name[3])
SCZ<-mappping_disorder_gene_PPI(PPi_intersection_label_chosing,disorder_name[4])
PD<-mappping_disorder_gene_PPI(PPi_intersection_label_chosing,disorder_name[5])
ASD<-mappping_disorder_gene_PPI(PPi_intersection_label_chosing,disorder_name[6])
BIP<-mappping_disorder_gene_PPI(PPi_intersection_label_chosing,disorder_name[7])
#### Co-expression network ###
ADHD<-mappping_disorder_gene_PPI(gene_cor_expression,disorder_name[1])
AD<-mappping_disorder_gene_PPI(gene_cor_expression,disorder_name[2])
MDD<-mappping_disorder_gene_PPI(gene_cor_expression,disorder_name[3])
SCZ<-mappping_disorder_gene_PPI(gene_cor_expression,disorder_name[4])
PD<-mappping_disorder_gene_PPI(gene_cor_expression,disorder_name[5])
ASD<-mappping_disorder_gene_PPI(gene_cor_expression,disorder_name[6])
BIP<-mappping_disorder_gene_PPI(gene_cor_expression,disorder_name[7])







###################### bag_mag_network ###
disorder_name<-colnames(disorder_table)
ADHD_TRN<-mappping_disorder_gene_PPI(bag_mag_network,disorder_name[1])
AD_TRN<-mappping_disorder_gene_PPI(bag_mag_network,disorder_name[2])
MDD_TRN<-mappping_disorder_gene_PPI(bag_mag_network,disorder_name[3])
SCZ_TRN<-mappping_disorder_gene_PPI(bag_mag_network,disorder_name[4])
PD_TRN<-mappping_disorder_gene_PPI(bag_mag_network,disorder_name[5])
ASD_TRN<-mappping_disorder_gene_PPI(bag_mag_network,disorder_name[6])
BIP_TRN<-mappping_disorder_gene_PPI(bag_mag_network,disorder_name[7])








#####
####   ##
PPI_TRN_common<-function(net1,net2){
 # net1<-SCZ
 # net2<-SCZ_TRN
  label_data<-net1[1,]
  for (j in c(1:nrow(net1))) {
    #j=1
    union_data<-union(net1[j,1],net1[j,2])
    for (i in c(1:nrow(net2))) {
      #i=1
      trn_gene<-union(net2[i,1],net2[i,2])
      if(length(intersect(union_data,trn_gene))>1){
        #label_data<-c(label_data,1)
        label_data<-rbind(label_data,net1[j,])
      }
    }
  } 
  label_data<-label_data[-1,]
  print(label_data)
  return(label_data)
}
PPI_TRN_common(ADHD,ADHD_TRN)
PPI_TRN_common(AD,AD_TRN)
PPI_TRN_common(MDD,MDD_TRN)
PPI_TRN_common(SCZ,SCZ_TRN)
PPI_TRN_common(PD,PD_TRN)
PPI_TRN_common(BIP,BIP_TRN)
###
###  ##
PPi_intersection_label_chosing[which(PPi_intersection_label_chosing$protein1=='9606.ENSP00000358404'),]%>%nrow()
PPi_intersection_label_chosing[which(PPi_intersection_label_chosing$protein1=='9606.ENSP00000360493'),]%>%nrow()
###
node<-union(PPi_intersection_label_chosing$protein1,PPi_intersection_label_chosing$protein2)
all_node_data<-PPI_info[node,2]
####
data_bag_data<-intersect(bag_gene,all_node_data)
data_inter<-intersect(all_node_data,mag_gene_sym)
set_other<-setdiff(all_node_data,union(data_bag_data,data_inter))
node_annote_data<-c(rep(1,length(data_bag_data)),rep(2,length(data_inter)),rep(3,length(set_other)))
data_node<-cbind(c(data_bag_data,data_inter,set_other),node_annote_data)
write.table(data_node,'subnet_all.txt',sep = '\t',row.names = FALSE,quote = FALSE)
##### bag gene ##
rownames(PPI_info)<-PPI_info$string_protein_id
node1_name<-PPI_info[PPi_intersection_label_chosing$protein1,2]
node2_name<-PPI_info[PPi_intersection_label_chosing$protein2,2]
####
PPi_intersection_label_chosing$protein1<-node1_name
PPi_intersection_label_chosing$protein2<-node2_name
write.table(PPi_intersection_label_chosing,'PPI_mag_net.txt',sep = '\t',row.names = FALSE,quote = FALSE)
## PPI pathway analysis #
all_pathway_data<-fread('9606.protein.enrichment.terms.v11.5.txt',sep = '\t',header = TRUE)
all_pathway_data_bp<-all_pathway_data[which(all_pathway_data$category=='Biological Process (Gene Ontology)'),]
TCF7L2_pathway<-all_pathway_data_bp[which(all_pathway_data_bp$`#string_protein_id`=='9606.ENSP00000358404',)]
## three gene in one pathway #
TCF7L2_pathway<-as.data.frame(TCF7L2_pathway)
TCF_bp_id<-TCF7L2_pathway$term
id_all_three<-c()
for (i in TCF_bp_id) {
  #i=TCF_bp_id[1]
  bp_protein_data<-all_pathway_data_bp[which(all_pathway_data_bp$term==i),]
  if(length(intersect(bp_protein_data$`#string_protein_id`,three_protein_data))>2){
    id_all_three<-c(id_all_three,i)
  }
}
rownames(TCF7L2_pathway)<-TCF7L2_pathway$term

target_pathway<-TCF7L2_pathway[id_all_three,]
bag_sig_pathway<-read.table('bag_pathway.txt',sep = '\t',header = TRUE,fill = TRUE)
mag_pathway_data<-read.table('mag_pathway.txt',sep = '\t',header = TRUE,fill = TRUE)
bag_sig_bp<-bag_sig_pathway[which(bag_sig_pathway$Category=='GO: Biological Process'),]
mag_pathway_data<-mag_pathway_data[which(mag_pathway_data$Category=='GO: Biological Process'),]
### plot the bag pathway ##
mag_pathway_data<-mag_pathway_data[c(1:20),]
bag_sig_bp<-bag_sig_bp[c(1:20),]
###
intersect(mag_pathway_data$ID,bag_sig_bp$ID)

mag_pathway_data$q.value.FDR.B.H<- -log10(mag_pathway_data$q.value.FDR.B.H)
ggdotchart(mag_pathway_data, x="Name", y="q.value.FDR.B.H",
           palette = c("#00AFBB", "#E7B800", "#FC4E07"), sorting = "descending", 
           add = "segments", rotate = TRUE, dot.size = 10, 
           font.label = list(color="white", size=8, vjust=0.5), 
           ggtheme = theme_pubr())
intersect(id_all_three,)


#####
##brain TRN load, explore the role of associated gene in brain TRN#
RUNX2_data<-TRN_network_data[which(TRN_network_data$node1=='RUNX2'),]
TCF7L2_data<-TRN_network_data[which(TRN_network_data$node1=='TCF7L2'),]
##
bag_tf_data<-intersect(TRN_node$Var1,bag_gene)
TRN_node<-table(TRN_network_data$node1)%>%as.data.frame()
rownames(TRN_node)<-TRN_node$Var1
###
bag_tf_degree<-TRN_node[bag_tf_data,]
other_degree<-TRN_node[setdiff(TRN_node$Var1,bag_gene),]
###



####

