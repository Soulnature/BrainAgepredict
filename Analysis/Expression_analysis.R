#########  the gene expression in brain  ##
# Load required libraries
library(tidyverse)
library(data.table)
library(preprocessCore)
library(ggpubr)
library(RColorBrewer)

# Define a function to read and preprocess expression data
preprocess_expression <- function(file_path, gene_files, sample_info_path) {
  # Read expression data
  expression <- fread(file_path, header = FALSE)
  # Extract gene IDs and clean up
  gene_ids <- expression$V1 %>% str_split('[|]') %>% map_chr(2)
  # Join expression data with gene IDs
  expression <- expression %>%
    select(-V1) %>%
    mutate(symbol = gene_ids) %>%
    distinct(symbol, .keep_all = TRUE)
  rownames(expression) <- expression$symbol
  expression <- select(expression, -symbol)
  
  # Read sample information
  sample_info <- read.csv(sample_info_path, header = TRUE)
  colnames(expression) <- sample_info$Sample
  
  # Read gene lists
  bag_gene <- fread(gene_files$bag, sep = '\t', header = FALSE)
  mag_gene <- fread(gene_files$mag, sep = '\t', header = FALSE)
  
  list(expression = expression, bag_gene = bag_gene, mag_gene = mag_gene, sample_info = sample_info)
}

# Define a function to normalize expression data
normalize_expression <- function(expression_data, sample_list, species) {
  sub_data <- expression_data[, sample_list$Sample, with = FALSE]
  sub_data <- as.matrix(sub_data)
  if (species == 'Macaque') {
    sub_data <- sub_data[, -1]
  }
  normalized_data <- normalize.quantiles(sub_data)
  colnames(normalized_data) <- colnames(sub_data)
  rownames(normalized_data) <- rownames(sub_data)
  normalized_data
}

# Define a function to plot expression data
plot_expression <- function(expression_data, gene_list, sample_periods, title) {
  expression_mean <- expression_data[intersect(rownames(expression_data), gene_list$V1),]
  expression_mean <- apply(expression_mean, 1, median)
  boxplot(expression_mean, main = title)
}

# Preprocess the data
data_files <- list(
  expression = './brain_development_RPKM.txt',
  bag = './10.17/bag_all_gene.txt',
  mag = 'mag_gene_list.txt'
)
sample_info_file <- './promoter_expression/RNA-sample.csv'
data <- preprocess_expression(data_files$expression, data_files, sample_info_file)

# Normalize expression data for human and macaque samples
human_expression <- normalize_expression(data$expression, data$sample_info[data$sample_info$Species == 'Human',], 'Human')
macaque_expression <- normalize_expression(data$expression, data$sample_info[data$sample_info$Species == 'Macaque',], 'Macaque')

# Plot expression data
plot_expression(human_expression, data$bag_gene, data$sample_info[data$sample_info$Species == 'Human',], 'Human BAG Gene Expression')
plot_expression(macaque_expression, data$bag_gene, data$sample_info[data$sample_info$Species == 'Macaque',], 'Macaque BAG Gene Expression')

boxplot(Macaque_expression_mean,Macaque_expression_mean)
RUNX2_exp_fea<-human_expression['KLF3',fetal_period$Sample]
RUNX2_exp_new<-human_expression['KLF3',newborn_period$Sample]
RUNX2_exp_adult<-human_expression['KLF3',adult_period$Sample]
##

##

plot_df<-data.frame(exp=c(RUNX2_exp_new,RUNX2_exp_adult),group_d=c(rep('young',length(RUNX2_exp_new)),rep('older',length(RUNX2_exp_adult))))
#
##
plot_df[which(plot_df$group_d=='young'),1]%>%mean()
plot_df[which(plot_df$group_d=='older'),1]%>%mean()

ggboxplot(plot_df,x='group_d',y='exp', add = "jitter",size=1)+ stat_compare_means()

ggboxplot(ToothGrowth, x="supp",
          y = "len", color = "supp",
          palette = "jco", add = "jitter")

boxplot(RUNX2_exp_new,RUNX2_exp_adult)
###
cerebra_f<-human_expression_mean[,newborn_period[which(newborn_period$NCXRegion=='NCX'),]$Sample]
sub_cor_f<-human_expression_mean[,newborn_period[which(newborn_period$NCXRegion!='NCX'),]$Sample]
###
cerebra_a<-apply(human_expression_mean[,adult_period[which(adult_period$NCXRegion=='NCX'),]$Sample], 1, median)
sub_cor_a<-apply(human_expression_mean[,adult_period[which(adult_period$NCXRegion!='NCX'),]$Sample],1,median)
            
plot(as.numeric(human_sample$Period),RUNX2_exp_adult)
boxplot(RUNX2_exp_new,RUNX2_exp_adult)
mean(RUNX2_exp_adult)/mean(RUNX2_exp_new)

###
age_group_ying<-human_sample$Period
###
pfc_exp<-human_sample[which(human_sample$Species=="Human"),1]
###
pfc_example<-human_sample[which(human_sample$Species=="Human"),]
table(pfc_example$Region)
cortex_sample<-pfc_example[which(pfc_example$NCXRegion=='NCX'),]
no_cortex_sample<-pfc_example[which(pfc_example$NCXRegion!='NCX'),]
#################The bag expression in  different brain regions ###############
brain_expression<-human_expression[intersect(rownames(human_expression),bag_gene$V1),pfc_example$Sample]
############# examine the region expression in different periods   ####
fetal_period<-pfc_example[which(pfc_example$Period<=7),]
newborn_period<-pfc_example[which(pfc_example$Period>7&pfc_example$Period<=11),]
adult_period<-pfc_example[which(pfc_example$Period>11&pfc_example$Period<14),]
chosing_region<-Reduce(intersect,list(fetal_period$Region,newborn_period$Region,adult_period$Region))
###
get_region_expression<-function(value_expression,period_sample,sample_name){
  period_sub<-value_expression[,period_sample$Sample]
  period_sub_mean<-apply(period_sub, 1, function(x){
    x<-tapply(x, period_sample$Region, median)
  })
  #####  
  inter_region_exp<-period_sub_mean[sample_name,]
  return(inter_region_exp)  
}
fetal_exp<-get_region_expression(brain_expression,fetal_period,chosing_region)%>%as.vector()
newborn_exp<-get_region_expression(brain_expression,newborn_period,chosing_region)%>%as.vector()
adult_exp<-get_region_expression(brain_expression,adult_period,chosing_region)%>%as.vector()
##
plot_module<-rbind(fetal_exp,newborn_exp,adult_exp)

pheatmap(brain_expression,color = colorRampPalette(c("blue", "white", "red"))(20),cluster_cols = FALSE,
         cluster_rows = FALSE,scale = 'row')
graph2ppt(file='coexpression_data',width=8,height=6,append=TRUE)
##
region_label<-c()
#period_label<-c(rep('fetal',512),rep('infant',512),rep('adult',512))
period_label<-c(rep('Early postnatal',512),rep('Adult',512))
for (i in chosing_region) {
  region_label<-c(region_label,rep(i,32))
}
region_label<-rep(region_label,2)
###
region_exp_df<-data.frame(exp_data=c(newborn_exp,adult_exp),period=period_label,region_name=region_label)
region_exp_df<-data.frame(exp_data=c(newborn_exp,adult_exp),period=period_label,region_name=region_label)
region_exp_df$period<-factor(region_exp_df$period,levels = c('Adult','Early postnatal'))
###
NCX<-pfc_example[which(pfc_example$NCXRegion=='NCX'),]$Region%>%unique()
NCX<-intersect(region_exp_df$region_name,NCX)
other<-setdiff(region_exp_df$region_name,NCX)
###
region_exp_df$region_name<-factor(region_exp_df$region_name,levels = region_exp_cortical$region)
###
ggboxplot(region_exp_df,x='region_name',y='exp_data',fill ='period',palette =c('#FCC8C5','#9FDFE7'),outlier.shape = NA)+
  stat_compare_means(aes(group = period),label = "p.signif")
###
region_exp_cortical<-tapply(region_exp_df$exp_data, region_exp_df$region_name, median)%>%as.data.frame()
region_exp_cortical<-cbind(rownames(region_exp_cortical),region_exp_cortical)
colnames(region_exp_cortical)<-c('region','expression')
####
NCX_order_mean<-region_exp_cortical[NCX,]
sub_mean_order<-region_exp_cortical[other,]
NCX_order_mean<-NCX_order_mean[order(NCX_order_mean$expression,decreasing = TRUE),]
sub_mean_order<-sub_mean_order[order(sub_mean_order$expression,decreasing = TRUE),]
region_exp_cortical<-rbind(NCX_order_mean,sub_mean_order)
region_exp_cortical$region<-factor(region_exp_cortical$region,levels =region_exp_cortical$region)
ggbarplot(region_exp_cortical,x='region',y='expression')

### trajectory ##
period_all<-intersect(cortex_sample$Period,no_cortex_sample$Period)%>%unique()
period_all<-period_all[period_all>7]
#### Compute the correaltion between BAG and MAG gene in different periods##
exp_all_data<-function(expressin_1,expression_2){
  # expressin_1<-bag_tem_expre
   #expression_2<-mag_tem_expre
  expressin_1<-expressin_1[!as.logical(rowSums(expressin_1==0)), ]
  expression_2<- expression_2[!as.logical(rowSums(expression_2==0)), ]
  co_expre_vec<-c()
  # expressin_1_mean<-apply(expressin_1, 2, mean)
 # expression_2_mean<-apply(expression_2, 2, mean)
 # co_expre_vec<-c(co_expre_vec,cor(expressin_1_mean,expression_2_mean))
  for (i in c(1:nrow(expressin_1))) {
   # i=1
    single_gene<-c()
    for (j in c(1:nrow(expression_2))) {
     # j=1
      temp_cor<-cor.test(expressin_1[i,],expression_2[j,],method = 'pearson')
      single_gene<-c(single_gene,temp_cor$p.value)
    }
    co_expre_vec<-c(co_expre_vec,single_gene)
  }
  return(co_expre_vec)
}
brain_expression_all<-human_expression[intersect(rownames(human_expression),union(bag_gene$V1,mag_gene$V1)),pfc_example$Sample]
period_cor_cal<-function(value_expression,period_sample,sample_name){
 # value_expression<-brain_expression_all
  #sample_name<-chosing_region
 # period_sample<-pfc_example
    period_sub_bag<-value_expression[intersect(rownames(value_expression),bag_gene$V1),]
    period_sub_mag<-value_expression[intersect(rownames(value_expression),mag_gene$V1),]
    other_sub_expression<-value_expression[setdiff(rownames(value_expression),mag_gene$V1),]
    save_data<-c()
    other_cor<-c()
    region_index<-c()
    for (i in sample_name) {
     # i=sample_name[9]
      sample_index<-period_sample[which(period_sample$Region==i),1]
      if(length(sample_index)>3){
        period_sub_index_bag<-period_sub_bag[,sample_index]
        period_sub_index_mag<-period_sub_mag[,sample_index]
        ###
        period_sub_index_bag<-period_sub_index_bag[!as.logical(rowSums(period_sub_index_bag==0)), ]
        period_sub_index_mag<-period_sub_index_mag[!as.logical(rowSums(period_sub_index_mag==0)), ]
        ######
        pca_bag<-prcomp(period_sub_index_bag)
        mag_pca<-prcomp(period_sub_index_mag)
        
        p_lest<-cor.test(pca_bag$rotation[,1],mag_pca$rotation[,1])
        
        print(p_lest$p.value)
        print(i)
       # period_sub_bag_mag<-exp_all_data(period_sub_index_bag,period_sub_index_mag)
        region_index<-c(region_index,rep(i,length(period_sub_bag_mag)))
        save_data<-c(save_data,p_lest$p.value)
      }
      }
    re_cor_df<-data.frame(cor_index=save_data,region_name=region_index)
   # re_cor_df$cor_index<--log10(re_cor_df$cor_index)
  
    ggboxplot(re_cor_df,x='region_name',y='cor_index')
    return(re_cor_df)
  }
fetal_exp_cor<-period_cor_cal(brain_expression_all,fetal_period,chosing_region)
newborn_exp_cor<-period_cor_cal(brain_expression_all,newborn_period,chosing_region)
adult_exp_cor<-period_cor_cal(brain_expression_all,adult_period,chosing_region)
#####
period_label_cor<-c(rep('Early postnatal',nrow(newborn_exp_cor)),rep('adult',nrow(adult_exp_cor)))
cor_df<-rbind(fetal_exp_cor,adult_exp_cor)
cor_df<-cbind(cor_df,period_label_cor)
cor_df<-cor_df[which(cor_df$cor_index>0),]
cor_df$cor_index<-abs(cor_df$cor_index)
ggboxplot(cor_df,x='region_name',y='cor_index',fill ='period_label_cor',outlier.shape = NA,palette =c('#FCC8C5','#868686'))
+ stat_compare_means(aes(group = period_label_cor),label = "p.signif")+ylab('Coexpression coefficient')

boxplot(fetal_exp_cor$cor_index,newborn_exp_cor$cor_index,adult_exp_cor$cor_index)

### Compare the expression in cortex and subcortical regions ##
get_region_difference<-function(bag_cortex,structure_name){
  # bag_cortex<-ying_new_data_brain_bag_nocortex
  structure_name<-chosing_region
  bag_expression<-human_expression[intersect(bag_gene$V1,rownames(human_expression)),pfc_example$Sample]
  other_expression<-human_expression[setdiff(rownames(human_expression),bag_gene$V1),pfc_example$Sample]
  bag_expression<-apply(bag_expression, 1, function(x){
    x<-tapply(x, pfc_example$Region, mean)
  })
  other_expression<-apply(other_expression, 1, function(x){
    x<-tapply(x, pfc_example$Region, mean)
  })
  ##
  bag_expression<-bag_expression[structure_name,]
  other_expression<-other_expression[structure_name,]
  plot_df_all<-data.frame(Expression=0,label='region',region='ss')
  for (i in rownames(bag_expression)) {
   # i= rownames(bag_expression)[1]
    re_bag<-bag_expression[i,]
    ot_bag<-other_expression[i,]
    label_name<-c(rep('BAG',length(re_bag)),rep('Others',length(ot_bag)))
    contaion_df<-c(re_bag,ot_bag)
    contaion_df<-data.frame(Expression=contaion_df,label=label_name,region=rep(i,length(contaion_df)))
    plot_df_all<-rbind(plot_df_all,contaion_df)
    ### permutation  test ##
    permutation_group<-c()
    for (i in c(1:10000)) {
      sample_temp<-sample(ot_bag,length(bag_expression))%>%median()
      permutation_group<-c(permutation_group,sample_temp)
    }
    permutation_group_or<-permutation_group[order(permutation_group,decreasing = TRUE)]
    ################
    permutation_group_big<-permutation_group_or[which(permutation_group_or>=median(re_bag))]%>%length()
    ###
    p<-permutation_group_big/length(permutation_group_or)
    print(p)
  }
####
  

  # contaion_df<-contaion_df[which(contaion_df$Expression!=0),]
  plot_df_all<-plot_df_all[-1,]
  plt<-ggboxplot(plot_df_all,x='label',y='Expression',fill = 'region',)+stat_compare_means()
  
  print(plt)
  #cortex_other<-gather(other_expression,sample,exp,1:ncol(other_expression))
  #######
  return(contaion_df)
}
cortex_df<-get_region_difference(ying_new_data_brain_bag_cortex,'cortex')
subcor_df<-get_region_difference(ying_new_data_brain_bag_nocortex,'subcortex')
##group plot ##
label_name<-c(rep('cortex',nrow(cortex_df)),rep('sub_cortical',nrow(subcor_df)))
diff_cor_df<-cbind(rbind(cortex_df,subcor_df),label_name)
ggviolin(diff_cor_df,x='label_name',y='Expression',fill = 'label',palette ='jco' )+stat_compare_means(aes(group = label),label = "p.format")


#### BAG and MAG co-expression ########
get_region_difference<-function(bag_cortex,structure_name){
  #qaqbag_cortex<-ying_new_data_brain_bag_nocortex
  bag_cortex<-human_expression
  structure_name<-chosing_region
  for (i in structure_name) {
    #i=structure_name[1]
    sample_info<-pfc_example[which(pfc_example$Region==i),1]
    bag_expression<-bag_cortex[intersect(bag_gene$V1,rownames(bag_cortex)),sample_info]
    mag_gene_expression<-bag_cortex[intersect(mag_gene$V1,rownames(bag_cortex)),sample_info]
    other_expression<-bag_cortex[setdiff(rownames(bag_cortex),union(bag_gene$V1,mag_gene$V1)),sample_info]
    bag_pc<-prcomp(bag_expression)
    mac_pc<-prcomp(mag_gene_expression)
    ### data ##
    other_cor<-c()
    for (i in c(1:10000)) {
      other_expression_sample<-other_expression[sample(nrow(other_expression),nrow(mag_gene_expression)),]
      other_pc<-prcomp(other_expression_sample)
      ote<-cor(bag_pc$rotation[,1],other_pc$rotation[,1])
      other_cor<-c(other_cor,ote)
    }
   ###
    other_cor_abs<-mean(other_cor)
    tar<-cor(bag_pc$rotation[,1],mac_pc$rotation[,1])

   print(c(tar,ote,i)) 
  }
  ###
  ## Compute the coeffiecent of the BAG and MAG ##
  bag_mag_coexp<-apply(bag_expression, 1, function(x){
    x<-bag_expression[1,]
    co_exp<-c()
    for (i in c(1:nrow(mag_gene_expression))) {
      co_exp<-c(co_exp,cor(as.numeric(x),as.numeric(mag_gene_expression[i,]),method ='pearson'))%>%abs()
    }
    ###
    x<-mean(co_exp)
  })
  bag_other_coexp<-apply(bag_expression, 1, function(x){
    # x<-bag_expression[1,]
    co_exp<-c()
    for (i in c(1:nrow(mag_gene_expression))) {
      co_exp<-c(co_exp,cor(as.numeric(x),as.numeric(other_expression[i,]),method ='pearson'))%>%abs()
    }
    ###
    x<-mean(co_exp)
  })
  #########################
  # subcortical permutation test #
  if(structure_name=='cortical'){
    permutation_test_cortex<-read.table('./cortex_d_permutation_res.txt')
    ####
    
  }else{
    permutation_test_cortex<-read.table('./no_cortex_d_permutation_res.txt')
    
  }
  permutation_test_cortex<-as.numeric(permutation_test_cortex$V1)
  permutation_test_cortex<-permutation_test_cortex[which(!is.na(permutation_test_cortex))]
  ### compute the p values #
  median_bag_mag<-median(bag_mag_coexp)
  Big_values<-permutation_test_cortex[which(permutation_test_cortex>median_bag_mag)]
  print(Big_values/10000)
#
#  subclortical_regions<-read.table('./no_cortex_d_permutation_res.txt')
# boxplot(subclortical_regions,bag_mag_coexp)
  ### plot the data
  coexp_plot_df<-data.frame(group_data=c(rep('Internal',length(bag_mag_coexp)),rep('External',length(bag_other_coexp))),
                            co_exp=c(bag_mag_coexp,bag_other_coexp),structure_name=rep(structure_name,2*length(bag_mag_coexp)))
  ##
 # ggviolin(coexp_plot_df,x='group_data',y='co_exp',fill ='group_data',palette = 'jco' )+ggtitle('bag_mag_exp & other_1000_no_cortex')
  return(coexp_plot_df)
}
##
cortex_coexp_subcortical<-get_region_difference(ying_new_data_brain_bag_nocortex,'subcortical')
cortex_coexp<-get_region_difference(ying_new_data_brain_bag_cortex,'cortical')
####
cortex_expression_df<-rbind(cortex_coexp,cortex_coexp_subcortical)
###
ggviolin(cortex_expression_df,x='structure_name',y='co_exp',fill = 'group_data')+stat_compare_means(aes(group = group_data))

##
# plot the coexpression data expression ##
bag_other_exp<-read.table('permutation_res_cortex.txt',header = TRUE)
bag_other_exp<-as.numeric(bag_other_exp$x)
###
coexp_plot_df<-data.frame(group_data=c(rep('Internal',length(bag_mag_coexp)),rep('External',length(bag_other_exp))),
                          co_exp=c(bag_mag_coexp,bag_other_exp))
##
ggboxplot(coexp_plot_df,x='group_data',y='co_exp',fill = )+ggtitle('bag_mag_exp & other_1000_no_cortex')
graph2ppt(file='cor_expression.pptx',width=8,height=6,append=TRUE)

###
wilcox.test(as.numeric(bag_other_exp$x),bag_mag_coexp)

cortex_bag<-gather(ying_new_data_brain_bag_cortex,sample,exp,1:ncol(ying_new_data_brain_bag_cortex))
subcortex_bag<-gather(ying_new_data_brain_bag_nocortex,sample,exp,1:ncol(ying_new_data_brain_bag_nocortex))
#######
label_name<-c(rep('cortex',nrow(cortex_bag)),rep('sub_cortical',nrow(subcortex_bag)))
contaion_df<-rbind(cortex_bag,subcortex_bag)
contaion_df<-cbind(contaion_df,label_name)
##
contaion_df$sample<-apply(as.matrix(contaion_df$sample),1,function(x){
  x<-str_split(x,'[.]')%>%unlist()
  x<-x[2]
})
colnames(contaion_df)<-c('Brain_region','Expression','Cortical_label')
##
contaion_order<-tapply(contaion_df$Expression, contaion_df$Brain_region, median)
contaion_order<-contaion_order[order(contaion_order,decreasing = TRUE)]
contaion_df$Brain_region<-factor(contaion_df$Brain_region,levels = names(contaion_order))
ggboxplot(contaion_df,x='Brain_region',y='Expression',fill = 'grey')+stat_compare_means()



### dynamic expression -- cortex vs no cortex #
trajectory<-function(expression_data,exp_sample){
  
 # expression_data<-brain_expression
  # exp_sample<-cortex_sample
  exp_tra<-c()
  period_df<-c()
  for (i in period_all) {
    # i=period_all[1]
    tem_per_sam<-exp_sample[which(exp_sample$Period==i),]
    period_data<-tem_per_sam$Period
    period_df<-c(period_df,period_data)
    tem_per_sam_exp<-expression_data[,tem_per_sam$Sample]
    tem_per_sam_exp<-apply(tem_per_sam_exp, 2, mean)
    # tem_per_sam_exp<-apply(tem_per_sam_exp, 1, mean)
    exp_tra<-c(exp_tra,tem_per_sam_exp)
  }
  return(data.frame(period_df,exp_tra))
}
cortex<-trajectory(brain_expression,cortex_sample)
sub_cortex<-trajectory(brain_expression,no_cortex_sample)
###
plot_df_data<-data.frame(Expression=rbind(cortex,sub_cortex),group=c(rep('cortex',nrow(cortex)),rep('sub_cortex',nrow(sub_cortex))))
colnames(plot_df_data)<-c('Periods','Expression','group')
plot_df_data$period <- factor(plot_df_data$Periods,levels = unique(plot_df_data$Periods))
##

ggplot(data=plot_df_data, 
       aes(x=Periods, y=Expression, group=group, color=factor(group),
           fill=factor(group)))+ geom_smooth(method='loess')+theme_bw()+
  theme(panel.grid =element_blank())+
  # theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  theme(panel.border = element_blank()) +   ## 删去外层边框
  theme(axis.line = element_line(size=1, colour = "black"))+
  xlab('periods') +ylab("Normalizated expression")

+
  geom_smooth(se=FALSE, linetype="dashed", size=0.5) +
  #   geom_xspline(spline_shape=0.4, size=0.5)
  +theme_bw()+theme(panel.grid=element_blank())


#
####
p = ggplot(data = deviates, aes(x = month)) +
  geom_point(aes(x = month, y = seasonal_deviate), color = col, size = 1) +
  geom_errorbar(aes(x = month, ymin = seasonal_deviate - sem, ymax = seasonal_deviate + sem), width = 0.5, color = col) + 
  stat_function(fun = cos_func, args = list(a = amplitude, phase = phase, omega = omega, intercept = 0), size = 0.7, color = col) +
  geom_ribbon(data = ci, aes(x = month, ymin = lower_ci, ymax = upper_ci), fill = col, alpha = 0.3) +
  scale_x_continuous(breaks=c(1, 3, 5, 7, 9, 11)) +
  ggtitle(title) +
  xlab("Month") +
  theme_classic() +
  theme(legend.position="none",
        plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_blank()
  )




######### single cell expression ##
meta_data<-read.csv('E:/yang/bag_gwas/results/singlecell/single_cell_meta.csv')
#####
all_bag_cell_type<-as.data.frame(as.matrix(all_bag_cell_type))
cell_types<-meta_data$cluster%>%unique()
#####
cell_exp_name<-c()
cell_exp_value<-c()
type_class_all<-c()
period_info<-c()
cell_type_class<-c("Olig","Micro", "Astro", "OPC", "InN", "ExN", "Purk")
## Plot the cell type  expression ##
for ( cellindex in c(1:length(cell_type_class))) {
  #cell_type<-cell_types[1]
  # cellindex=1
  if(cellindex<=4){
    type_class<-'Glia'
  }else{
    type_class<-'Neuron'
  }
  cell_type<-cell_type_class[cellindex]
  sample_cell_sample_fetal<-meta_data[which(meta_data$cluster==cell_type&meta_data$Period=='P6'),1]
  sample_cell_sample_adult<-meta_data[which(meta_data$cluster==cell_type&meta_data$Period=='P14'),1]
  # sample_cell_sample<-cbind(sample_cell_sample_fetal,sample_cell_sample_adult)
  # cell_type_expression<-all_bag_cell_type[,sample_cell_sample]%>%as.matrix()
  #
  sample_cell_exp_fetal<-all_bag_cell_type[,sample_cell_sample_fetal]%>%as.matrix()
  sample_cell_exp_adult<-all_bag_cell_type[,sample_cell_sample_adult]%>%as.matrix()
  cell_type_expression_fetal<-apply(sample_cell_exp_fetal, 1, function(x){
    x<-x[which(x!=0)]
    x<-mean(x)
  })
  cell_type_expression_adult<-apply(sample_cell_exp_adult, 1, function(x){
    x<-x[which(x!=0)]
    x<-mean(x)
  })
  ##
  cell_exp_name<-c(cell_exp_name,rep(cell_type,length(c(cell_type_expression_fetal,cell_type_expression_adult))))
  cell_exp_value<-c(cell_exp_value,c(cell_type_expression_fetal,cell_type_expression_adult))
  type_class_all<-c(type_class_all,rep(type_class,length(c(cell_type_expression_fetal,cell_type_expression_adult))))
  period_info<-c(period_info,c(rep('P6',length(cell_type_expression_fetal)),rep('P14',length(cell_type_expression_adult))))
}
##
single_df<-data.frame(cell_name=cell_exp_name,cell_exp_val=cell_exp_value,type_class=type_class_all,period=period_info)%>%na.omit()
##
adult_cell<-single_df[which(single_df$period=='P14'),]
early_cell<-single_df[which(single_df$period=='P6'),]
wilcox.test(adult_cell$cell_exp_val,early_cell$cell_exp_val)

##
ggboxplot(single_df,x='cell_name',y='cell_exp_val', color='type_class',facet.by ='period',add='jitter',palette = 'npj',outlier.shape= NA)
graph2ppt(file='cell_type_expression',height=8,width=10,append=TRUE)
###### plot the expresion trajectory of the single cell ####
glia_cell<-c('Olig','Micro','Astro','OPC')
neuro_cell<-c('InN','ExN')
get_glia_cell_type<-meta_data[1,]
for (i in c(1:length(glia_cell))) {
  tem_sample<-meta_data[which(meta_data$cluster==glia_cell[i]),]
  get_glia_cell_type<-rbind(get_glia_cell_type,tem_sample)
}
get_glia_cell_type<-get_glia_cell_type[-1,]
###
get_neruo_cell_type<-meta_data[1,]
for (i in c(1:length(neuro_cell))) {
  tem_sample<-meta_data[which(meta_data$cluster==neuro_cell[i]),]
  get_neruo_cell_type<-rbind(get_neruo_cell_type,tem_sample)
}
get_neruo_cell_type<-get_neruo_cell_type[-1,]
#########
cell_types_mean_exp<-apply(all_bag_cell_type, 2, function(x){
  x<-x[which(x!=0)]
  x<-mean(x)
})%>%as.data.frame()
#######
glia_exp_data<-cell_types_mean_exp[get_glia_cell_type$X,]
neuro_exp_data<-cell_types_mean_exp[get_neruo_cell_type$X,]
glia_exp_data<-data.frame(exp=glia_exp_data,cell_type=get_glia_cell_type$cluster,period=get_glia_cell_type$Period)%>%na.omit()
neuro_exp_data<-data.frame(exp=neuro_exp_data,cell_type=get_neruo_cell_type$cluster,period=get_neruo_cell_type$Period)%>%na.omit()

get_cell_type_expression<-function(glia_exp_data){
  period_tem<-glia_exp_data$period%>%unique()
  cell_type_info<-c()
  period_data<-c()
  for (i in period_tem) {
    # i=period_tem[1]
    period_d<-glia_exp_data[which(glia_exp_data$period==i),]
    period_d<-tapply(period_d$exp, period_d$cell_type, mean)%>%as.vector()
    cell_type_info<-c(cell_type_info,period_d)
    period_data<-c(period_data,rep(i,length(period_d)))
  }
  return(data.frame(exp=cell_type_info,period=period_data))
  
  
}
glia_exp_df<-get_cell_type_expression(glia_exp_data)
neuro_exp_df<-get_cell_type_expression(neuro_exp_data)
######
group=c(rep('glia',nrow(glia_exp_df)),rep('neuro',nrow(neuro_exp_df)))
trajetory_cell_tyep<-cbind(rbind(glia_exp_df,neuro_exp_df),group)
trajetory_cell_tyep$period<-apply(as.matrix(trajetory_cell_tyep$period), 1,function(x){
  x<-str_split(x,'P')%>%unlist()
  x<-x[2]
})
trajetory_cell_tyep<-na.omit(trajetory_cell_tyep)
trajetory_cell_tyep$period<-as.numeric(trajetory_cell_tyep$period)
trajetory_cell_tyep<-trajetory_cell_tyep[which(trajetory_cell_tyep$period>11),]
ggplot(data=trajetory_cell_tyep, 
       aes(x=period, y=exp, group=group, color=factor(group),
           fill=factor(group)))+ geom_smooth(method='loess')+theme_bw()+
  theme(panel.grid =element_blank())+
  # theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  theme(panel.border = element_blank()) +   ## 
  theme(axis.line = element_line(size=1, colour = "black"))+
  xlab('periods') +ylab("Normalizated expression")
+
  geom_smooth(se=FALSE, linetype="dashed", size=0.5) +
  #   geom_xspline(spline_shape=0.4, size=0.5)
  +theme_bw()+theme(panel.grid=element_blank())
####





