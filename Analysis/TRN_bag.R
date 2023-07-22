library(dplyr)
library(stringi)
library(stringr)
TRN_info<-read.csv('./TRN.csv')
##



##
TRN_number<-intersect(unique(TRN_info$node1),bag_gene$V1)
#####
degree_data<-table(TRN_info$node1)%>%as.data.frame()
##
rownames(degree_data)<-degree_data$Var1
bag_degree<-degree_data[TRN_number,]

mag_gene<-read.table('mag_gene_list.txt')
bag_gene<-read.table('./10.17/bag_all_gene.txt')
### BAG and MAG TRN reult ##
mag_bag_gene<-union(mag_gene$V1,bag_gene$V1)
####
TF_data<-intersect(mag_bag_gene,TRN_info$node1)

bag_mag_network<-TRN_info[1,]
####
for (i in c(1:nrow(TRN_info))) {
  ###
  if(length(intersect(TRN_info[i,1],mag_bag_gene))>0 || length(intersect(TRN_info[i,2],mag_bag_gene))>0 ){
    bag_mag_network<-rbind(bag_mag_network,TRN_info[i,])
  }
}
bag_mag_network<-bag_mag_network[-1,]
#Compute the coefficient of bag and mag  {connected mag/all bag} / (connected other node/ all other gene number)
##
all_node_info<-union(bag_mag_network$node1,bag_mag_network$node2)%>%unique()
##
other_network<-TRN_info[1,]
bag_mag_connected<-TRN_info[1,]
for (j in 1:nrow(bag_mag_network)) {
  #j=1
  row_gene_number<-union(bag_mag_network[j,1],bag_mag_network[j,2])
  #j=1
    if(length(intersect(row_gene_number,mag_bag_gene))>1){
      bag_mag_connected<-rbind(bag_mag_connected,bag_mag_network[j,])
    }else{
      other_network<-rbind(other_network,bag_mag_network[j,])
    }
}
bag_mag_connected<-bag_mag_connected[-1,]
other_network<-other_network[-1,]
## Compute the coefficient ##
innner_node<-union(bag_mag_connected$node1,bag_mag_connected$node2)%>%unique()
out_node<-union(other_network$node1,other_network$node2)%>%unique()
###
innner_edge<-nrow(bag_mag_connected)
other_edge<-nrow(other_network)
### 1.579262 
(length(innner_node)/length(out_node))/(innner_edge/other_edge)
###
KLF3_case_data<-bag_mag_connected[which(bag_mag_connected$node1=='KLF3'),]
KLF3_gene<-c(KLF3_case_data$node1,KLF3_case_data$node2)%>%unique()
write.table(KLF3_gene,'KLF3.txt',row.names = FALSE,quote = FALSE)
##
####
backgound_gene_all<-setdiff(union(TRN_info$node1,TRN_info$node2),mag_bag_gene)%>%length()
#####
mag_gene_connect<-setdiff(union(bag_mag_connected$node1,bag_mag_connected$node2),bag_gene$V1)%>%length()
other_gene_connect<-setdiff(union(other_network$node1,other_network$node2),bag_gene$V1)%>%length()
###### compute the effient ##
coffient<-(mag_gene_connect/nrow(mag_gene))/(other_gene_connect/backgound_gene_all) ## 6.725177
###
connect_gene<-intersect(mag_bag_gene,union(bag_mag_connected$node1,bag_mag_connected$node2)%>%unique())


##
#### save the BAG MAG network #
write.table(bag_mag_connected,'bag_mag_connected',row.names = FALSE,quote = FALSE,sep = '\t')
####### label the node #
all_node_data<-union(bag_mag_network$node1,bag_mag_network$node2)%>%unique()
### 
other_gene<-setdiff(all_node_data,mag_bag_gene)
other_dif<-rep(3,length(other_gene))
bag_label<-rep(1,length(bag_gene$V1))
mag_label<-rep(2,length(mag_gene$V1)-2)
mag_gene_ot<-setdiff(mag_gene$V1,bag_gene$V1)
###
label_all<-c(bag_label,mag_label)
all_gene<-c(bag_gene$V1,mag_gene_ot)
############
anno_node<-data.frame(all_gene=all_gene,label_gene=label_all)
write.table(anno_node,'node_anno_inter.txt',row.names = FALSE,quote = FALSE,sep = '\t')
####### label the node #
other_Tf<-degree_data[setdiff(unique(TRN_info$node1),bag_gene$V1),]
### the BAG and MAG in TRN ##



