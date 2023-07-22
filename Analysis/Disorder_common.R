###
library(dplyr)
library(UpSetR)         #Upset图（upset 包，适用样本数 2-7）
library(VennDiagram) 
library(export)
library(pheatmap)
# gene node ##
AD_node<-read.table('Alzheimer.s.Diseasenode.txt',header = TRUE)
PD_node<-read.table('Parkinson.s.Diseasenode.txt',header = TRUE)
ASD_node<-read.table('Autism.Spectrum.Disordernode.txt',header = TRUE)
BIP_node<-read.table('Bipolar.Disordernode.txt',header = TRUE)
scz_node<-read.table('Schizophrenianode.txt',header = TRUE)
dep_node<-read.table('Depressionnode.txt',header = TRUE)
ADHD_node<-read.table('Attention.deficit.hyperactivity.disordernode.txt',header = TRUE)
AD_gene<-AD_node[which(AD_node$Alzheimer.s.Disease_anno==2|AD_node$Alzheimer.s.Disease_anno==3),]$gene
pd_gene<-PD_node[which(PD_node$Parkinson.s.Disease_anno==2),]$gene
ASD_gene<-ASD_node[which(ASD_node$Autism.Spectrum.Disorder_anno==2|ASD_node$Autism.Spectrum.Disorder_anno==3),]$gene
BIP_gene<-BIP_node[which(BIP_node$Bipolar.Disorder_anno==2|BIP_node$Bipolar.Disorder_anno==3),]$gene
scz_gene<-scz_node[which(scz_node$Schizophrenia_anno==2|scz_node$Schizophrenia_anno==3),]$gene
dep_gene<-dep_node[which(dep_node$Depression_anno==2|dep_node$Depression_anno==3),]$gene
ADHD_gene<-ADHD_node[which(ADHD_node$Attention.deficit.hyperactivity.disorder_anno==2|ADHD_node$Attention.deficit.hyperactivity.disorder_anno==3),]$gene
## gene edge ##
AD_edge<-read.table('Alzheimer.s.Diseasenet.txt',header = TRUE)
PD_edge<-read.table('Parkinson.s.Diseasenet.txt',header = TRUE)
ASD_edge<-read.table('Autism.Spectrum.Disordernet.txt',header = TRUE)
BIP_edge<-read.table('Bipolar.Disordernet.txt',header = TRUE)
scz_edge<-read.table('Schizophrenianet.txt',header = TRUE)
dep_edge<-read.table('Depressionnet.txt',header = TRUE)
ADHD_edge<-read.table('Attention.deficit.hyperactivity.disordernet.txt',header = TRUE)
###  BAG TFs regulated gene ##
gene_data_tf_f<-function(disorder_net,bag_gene){
  #disorder_net<-AD_edge
  inter_tfs<-intersect(disorder_net$node1,bag_gene$V1)
  save_edge_d<-disorder_net[1,]
  for (i in inter_tfs) {
    te_gene<-disorder_net[which(disorder_net$node1==i),]
    save_edge_d<-rbind(save_edge_d,te_gene)
  }
  save_edge_d<-save_edge_d[-1,]
  ###
  return(unique(save_edge_d$node2))
}
bag_ad_re<-gene_data_tf_f(AD_edge,bag_gene)
bag_pd_re<-gene_data_tf_f(PD_edge,bag_gene)
bag_asd_re<-gene_data_tf_f(ASD_edge,bag_gene)
bag_bip_re<-gene_data_tf_f(BIP_edge,bag_gene)
bag_scz_re<-gene_data_tf_f(scz_edge,bag_gene)
bag_dep_re<-gene_data_tf_f(dep_edge,bag_gene)
bag_adhd_re<-gene_data_tf_f(ADHD_edge,bag_gene)


####
#whole_net<-rbind(AD_edge,PD_edge,ASD_edge,BIP_edge,scz_edge,dep_edge,ADHD_edge)
bag_gene_net<-intersect(union(whole_net$node1,whole_net$node2),bag_gene$V1)
##
#upset_list<-list(bag_ad_re,bag_pd_re,bag_asd_re,bag_bip_re,bag_scz_re,bag_dep_re,bag_adhd_re)
upset_list<-list(AD_gene,pd_gene,ASD_gene,BIP_gene,scz_gene,dep_gene,ADHD_gene)
names(upset_list) <-c('AD','PD','ASD','BIP','SCZ','DEP','ADHD')
###
upset(fromList(upset_list),  # fromList一个函数，用于将列表转换为与UpSetR兼容的数据形式。
      nsets = 100,     # 绘制的最大集合个数
      nintersects = 40, #绘制的最大交集个数，NA则全部绘制
      order.by = "freq", # 矩阵中的交点是如何排列的。 "freq"根据交集个数排序，"degree"根据
      keep.order = F, # 保持设置与使用sets参数输入的顺序一致。默认值是FALSE，它根据集合的大小排序。
      mb.ratio = c(0.6,0.4),   # 左侧和上方条形图的比例关系
      text.scale = 2 # 文字标签的大小
)
## 
graph2ppt(file='TRN_upset.pptx',append=TRUE,width=10,height=8)
##
inter <- get.venn.partitions(upset_list)
inter_c<-inter[which(inter$..count..>0),]
##
cont_symbol<-inter_c[,1:7]
cont_symbol_num<-apply(cont_symbol, 1, function(x){
  x<-as.numeric(x)
})
cont_symbol_num<-t(cont_symbol_num)
#
inter_c[,1:7]<-cont_symbol_num
###
annota_node<-inter_c[1,1:7]
gene<-c()
for (i in c(1:nrow(inter_c))) {
  #i<-10
  gene_number<-inter_c[i,]$..count..
  for (j in c(1:gene_number)) {
    annota_node<-rbind(annota_node,inter_c[i,1:7])
  }
  gene<-c(gene,inter_c[i,]$..values..%>%unlist())
}
###
annota_node<-annota_node[-1,]
annota_node<-cbind(gene,annota_node)
#####
write.table(annota_node,'node_all_data.txt',sep = '\t',row.names = FALSE,quote = FALSE)
write.table(whole_net,'wholenet.txt',sep = '\t',row.names = FALSE,quote = FALSE)
####
node_annote<-apply(annota_node, 1, function(x){
 # x<-annota_node[1,]
  if(sum(as.numeric(x[2:8]))<2){
    x<-0
  }else{
    x<-1
  }
})
####
gene_shared_data<-annota_node[which(node_annote==1),]
###
shared_gene<-gene_shared_data$gene
share_node_index<-apply(whole_net,1, function(x){
  if(length(intersect(union(x[1],x[2]),shared_gene))>1){
    x<-1
  }else{
    x<-0
  }
  
})
shared_net<-whole_net[which(share_node_index==1),]

write.table(shared_net,'share_net.txt',row.names = FALSE,quote = FALSE,sep = '\t')
write.table(gene_shared_data,'gene_node_shared.txt',row.names = FALSE,quote = FALSE,sep ='\t')
###
gene_share_node<-read.table('gene_node_shared.txt',header = TRUE)
#
gene_share_node<-gene_share_node[,-ncol(gene_share_node)]
########
write.table(gene_share_node ,'gene_share_node_noBAG.txt',sep = '\t',row.names = FALSE,quote = FALSE)
#####
TF_mode<-gene_share_node$gene%>%unique()
###
TFs<-intersect(bag_gene$V1,TF_mode)
unTFs<-setdiff(TF_mode,bag_gene$V1)
###  Model  ###
TF_mode<-data.frame(gene=c(TFs,unTFs),mode=c(rep(1,length(TFs)),rep(2,length(unTFs))))
###
write.table(TF_mode,'TF_mode.txt',row.names = FALSE,quote = FALSE,sep = '\t')
#### Pathway  analysis ##
TCF7L2_path<-fread('./TCF7L2.txt')
KLF16_path<-fread('./KLF16.txt')
HES7_path<-fread('./HES7.txt')
p53<-fread('./TP53.txt')
TCF7L2_path<-TCF7L2_path[which(TCF7L2_path$Category=='GO: Biological Process'),]
KLF16_path<-KLF16_path[which(KLF16_path$Category=='GO: Biological Process'),]
HES7_path<-HES7_path[which(HES7_path$Category=='GO: Biological Process'),]
p53<-p53[which(p53$Category=='GO: Biological Process'),]
###
TCF7L2_path<-TCF7L2_path[1:10,]
KLF16_path<-KLF16_path[1:10,]
HES7_path<-HES7_path[1:10,]
###
row_matrix_name<-unique(c(TCF7L2_path$Name,KLF16_path$Name,HES7_path$Name))
####### ##########
matrix_data<-matrix(0,nrow =length(row_matrix_name),ncol = 3)
rownames(matrix_data)<-row_matrix_name
colnames(matrix_data)<-c('TCF7L2','KLF16','HES7')
##,=#
matrix_data[TCF7L2_path$Name,1]<--log10(TCF7L2_path$`q-value FDR B&H`)
matrix_data[KLF16_path$Name,2]<--log10(KLF16_path$`q-value FDR B&H`)
matrix_data[HES7_path$Name,3]<--log10(HES7_path$`q-value FDR B&H`)
###

pheatmap(matrix_data, color = colorRampPalette(colors = c("white","#C00000"))(20),treeheight_col = 0, treeheight_row = 0)





####
