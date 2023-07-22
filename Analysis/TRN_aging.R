## Aging pathway interaction analysis ##
brain_TRN_analysis_dta<-fread('./ADNI/PPI_intersection/TRN.csv')
brain_TRN_net<-brain_TRN_analysis_dta[,1:4]
######
intersect(aging_gene_set$Symbol,bag_gene$V1)

##
cell_aging_clock<-fread('cell_aging.csv')%>%as.data.frame()
###
colnames(cell_aging_clock)<-c('Symbol','Gene_Set','state','Reference')
##



##
aging_gene_set<-fread('./age_gene_set.csv',fill = TRUE)
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
#
intersect(mag_gene$V1,aging_gene_set$Symbol)
##
get_interaction_pathway_data<-function(bag_gene,mag_gene,age_set){
 # age_set<-aging_gene_set
  bag_sub_net<-get_split_gene(brain_TRN_net,bag_gene$V1)
  mag_sub_net<-get_split_gene(brain_TRN_net,mag_gene$V1)
  bag_mag_inter<-get_split_gene(mag_sub_net[,1.:3],bag_gene$V1)
  ###########
  pathway_info<-matrix(NA,nrow = length(unique(age_set$Gene_Set)),3)
  rownames(pathway_info)<-unique(age_set$Gene_Set)
  all_class_age<-unique(age_set$Gene_Set)
  for (i in c(1:length(all_class_age))) {
   #i<-1
    age_pathway_da<-age_set[which(age_set$Gene_Set==all_class_age[i]),]$Symbol
    print(intersect(age_pathway_da,bag_gene$V1))
    print(all_class_age[i])
    age_pathway_sub_net<-get_split_gene(brain_TRN_net,age_pathway_da)
    ####
    age_bag_sub_net<-get_split_gene(age_pathway_sub_net[,1:3],bag_gene$V1)
    age_mag_sub_net<-get_split_gene(age_pathway_sub_net[,1:3],mag_gene$V1)
    ####
    pathway_info[i,1]<-length(age_pathway_da)
    pathway_info[i,2]<-nrow(age_bag_sub_net)
    pathway_info[i,3]<-nrow(age_mag_sub_net)
    ##
  }
  ##
  bag_df<-data.frame(node1=c(rep('BAG',12)),node2=rownames(pathway_info),edge=pathway_info[,2])
  mag_df<-data.frame(node1=c(rep('MAG',12)),node2=rownames(pathway_info),edge=pathway_info[,3])
  bag_df<-bag_df[-1,]
  mag_df<-mag_df[-1,]
  ##
  bag_df$edge<-scale(bag_df$edge,center = FALSE)
  mag_df$edge<-scale(mag_df$edge,center = FALSE)
  bag_mag<-data.frame(node1=c('MAG'),node2=c('BAG'),edge=3)
  age_df_data<-rbind(bag_df,mag_df,bag_mag)
  age_df_data<-age_df_data[which(age_df_data$edge>0),]
  return(age_df_data)
}
###
get_inter_gene<-function(bag_gene,mag_gene,age_set,network){
  bag_sub_net<-get_split_gene(brain_TRN_net,bag_gene$V1)
  mag_sub_net<-get_split_gene(brain_TRN_net,mag_gene$V1)
  bag_mag_inter<-get_split_gene(mag_sub_net[,1.:3],bag_gene$V1)
  ####
  symbol_d<-c()
  # created vector with 5 characters
  columns= c("node1","node2","sign","PPi_intersection_save_label")

  save_net_w = data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(save_net_w) = columns
  all_class_age<-unique(age_set$Gene_Set)
  for (i in c(1:length(all_class_age))) {
    #i<-1
    age_pathway_da<-age_set[which(age_set$Gene_Set==all_class_age[i]),]$Symbol
    print(intersect(age_pathway_da,bag_gene$V1))
    print(all_class_age[i])
    age_pathway_sub_net<-get_split_gene(brain_TRN_net,age_pathway_da)
    age_bag_sub_net<-get_split_gene(age_pathway_sub_net[,1:3],bag_gene$V1)
    save_net_w<-rbind(save_net_w,age_bag_sub_net)
    symbol_d<-c(symbol_d,rep(all_class_age[i],nrow(age_bag_sub_net)))
  }
  save_net_w<-cbind(save_net_w,symbol_d)
  ####################################
  ##
  save_net_w<-save_net_w[which(save_net_w$node1=='KLF3'),]
  ##
  kcl_example<-save_net_w[which(save_net_w$symbol_d!='others'),]
  ##
  
  write.table(kcl_example[,1:3],'kcl_example.txt',row.names = FALSE,quote = FALSE,sep = '\t')
  ###
  bag_mag_example<-bag_mag_inter[which(bag_mag_inter$node1=='KLF3'),]
  #####

  
  
  #
  
  
  
  
  
  ###
  
  
}




