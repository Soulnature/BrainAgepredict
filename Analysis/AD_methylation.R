## ad 
library(glmnet)
library(dplyr)
library(stringi)
library(stringr)
library(ggpubr)
data(IlluminaHumanMethylation450kmanifest)
##
annotation_data<-IlluminaHumanMethylation450kmanifest@data
SNP_data<-annotation_data$%>%as.data.frame()
SNP_data$nCpG
load('age_methylation.Rdata')
load('./elas_methy_data.Rdata')
s353_tar<-read.csv('353_cgp.csv')
rownames(s353_tar)<-s353_tar$CpGmarker
# Elat coefficient data #
elas_coef<-read.table('en.coef',header  =TRUE)
inter_sec_CpG<-intersect(elas_coef$probe,s353_tar$CpGmarker)
##
rownames(s353_tar)<-s353_tar$CpGmarker
AD_sample_data<-read.csv('ADNI_bag_methylation_data.csv')
bag_image<-AD_sample_data$AD_pred-AD_sample_data$Age
##
AD_sample_data<-cbind(AD_sample_data,bag_image)
s353_tar_chosing_data<-s353_tar[rownames(methy_data_all_sample),]
#########################################
normal_sample<-AD_sample_data[which(AD_sample_data$Group=='CN'),]
############### data ###############
bag_data<-normal_sample$bag_image
#################
get_id<-function(normal_sample){
  image_name<-normal_sample$Loc%>%as.matrix()
  name_img<-apply(image_name, 1, function(x){
   # x<-image_name[1]
    x_sub<-str_split(x,'[/]')%>%unlist()
    x_sub<-x_sub[6]
    x_sub<-str_split(x_sub,'[_]')%>%unlist()
    x_sub<-x_sub[1]
  })
  return(name_img)
}
normal_image_id<-get_id(normal_sample)
names(bag_data)<-normal_image_id
## analysis the normal data ##
methy_data_all_sample<-elas_methylation_data
methylation_normal<-methy_data_all_sample[,normal_image_id] ### 151 sample ##
## test the association between the bag and methylation 
correlation_test<-function(bag,methylation_data){
  # immediate correlation test 
  association_value<-matrix(0,311,1)
  ##
  for (i in c(1:nrow(methylation_data))) {
   p_test_cor<- cor.test(as.numeric(bag),as.numeric(methylation_data[i,]))
   association_value[i]<-p_test_cor$p.value
  }
  association_value<-p.adjust(association_value,method = 'fdr')
  association_value<-as.matrix(association_value)
  rownames(association_value)<-rownames(methylation_data)
  #association_value<-association_value[which(association_value<=0.05),]
  association_value<-cbind(probe_annotation,association_value)
  return(association_value)                                    
}
perarson_test<-correlation_test(bag_data,methylation_normal)
###### linear regresion #
regression_data<-cbind(t(methylation_normal),bag_data)%>%as.data.frame()%>%na.omit()
sst<-glm(bag_data~.,family = gaussian,regression_data)
#####
summary_data_glm<-summary(sst)
coefficient_sst<-sst$coefficients
coefficient_sst<-coefficient_sst[which(!is.na(coefficient_sst))]%>%names()
coefficients_data<-summary_data_glm$coefficients
coefficients_data<-coefficients_data[coefficient_sst,]
#coefficients_data[,4]<-p.adjust(coefficients_data[,4],method = 'BH')
coefficients_data<-coefficients_data[which(coefficients_data[,4]<=0.05),]
###
coefficients_data<-cbind(s353_tar_chosing_data[rownames(coefficients_data),c(1,10,12,15,16,17,18)],coefficients_data)
##
write.csv(coefficients_data,'methylation_cpg_bag.txt',row.names = FALSE)
write.csv(perarson_test,'pearson_cor_cpg_bag.txt',row.names = FALSE)
gene_id_20<-coefficients_data$Gene_ID
write.table(gene_id_20,'gene_id20.txt',row.names = FALSE,quote = FALSE,col.names = FALSE) 
##
negtative_data<-coefficients_data[which(coefficients_data$`t value`<0),]
positive_data<-coefficients_data[which(coefficients_data$`t value`>0),]
########### compare the expression data for the 12 cpG island #
all_variable_data<-read.csv('allVariable.csv')

methy_data_all_sample<-all_variable_data[,!duplicated(colnames(methy_data_all_sample))]
all_variable_data<-all_variable_data[which(!duplicated(all_variable_data$Sample_Name)),]
normal_data<-all_variable_data[which(all_variable_data$Sample_Group=='CN'),]
MCI_data<-all_variable_data[which(all_variable_data$Sample_Group=='MCI'),]
AD_data<-all_variable_data[which(all_variable_data$Sample_Group=='AD'),]
#########################
normal_methy<-methy_data_all_sample[,as.character(normal_data$Sample_Name)]
MCI_methy<-methy_data_all_sample[,as.character(MCI_data$Sample_Name)]
AD_methy<-methy_data_all_sample[,as.character(AD_data$Sample_Name)]
#########################
normal_methy<-apply(normal_methy, 1, mean)
MCI_methy<-apply(MCI_methy, 1, mean)
AD_methy<-apply(AD_methy, 1, mean)
## ggplot the data #
library(ggplot2)
data_com<-c(normal_methy,MCI_methy,AD_methy)
label_group<-c(rep('CN',length(normal_methy)),rep('MCI',length(MCI_methy)),rep('AD',length(AD_methy)))
##
label_group<-factor(label_group,levels = c('CN','MCI','AD'))
data_fram_lan<-data.frame(x=label_group,methylation=data_com)
ggplot(data_fram_lan,aes(x=x,y=methylation,fill=x))+geom_boxplot()
#######  the brain age and methlation ##
methylation_data_elas<-read.table('methylation_pred_age_elastic.txt')
methylation_data_blue<-read.table('methylation_pred_age_blue.txt')
#########
all_variable_data<-read.csv('allVariable.csv')
######################################
all_variable_data_elas<-cbind(all_variable_data,t(methylation_data_elas))
all_variable_data_blue<-cbind(all_variable_data,t(methylation_data_blue))
all_variable_data_elas<-all_variable_data_elas[which(!duplicated(all_variable_data_elas$Sample_Name)),]
all_variable_data_blue<-all_variable_data_blue[which(!duplicated(all_variable_data_blue$Sample_Name)),]
#############################         ###### 
elas_bag<-all_variable_data_elas$`1`-all_variable_data_elas$age%>%as.matrix()
blue_bag<-all_variable_data_blue$age-all_variable_data_blue$`1`%>%as.matrix()
##
age<-cbind(bag_elas_group$X,bag_elas_group$age,bag_elas_group$`1`,bag_blue_group$`1`)
colnames(age)<-c("ID","age.raw","Elastic_Net","BLUP")
age<-as.data.frame(age)
cor(age$age.raw,age$BLUP)
mean(abs((age$age.raw-age$Elastic_Net)))
library(ggplot2)
library(reshape2)
library(ggpubr)
ggplot(age,aes(x=age.raw,y=Elastic_Net,fill='red')) +geom_smooth(method = 'lm')+geom_point(color='black')+xlab('Chronological')+ylab('Predicted age')+
  ggtitle('Elastic_Net')+theme(title=element_text(,size=20))+theme_bw()+theme(panel.grid=element_blank())
#######
mean(abs(elas_bag))
mean(abs(blue_bag))
bag_elas_group<-cbind(all_variable_data_elas,elas_bag)
bag_blue_group<-cbind(all_variable_data_blue,blue_bag)
######################################################
plot_group_difference<-function(data){
  #data<-bag_elas_group
  normal_data<-data[which(data$Sample_Group=='CN'),13]
  MCI_data<-data[which(data$Sample_Group=='EMCI'),13]
  AD<-data[which(data$Sample_Group=='AD'),13]%>%abs()
  data_com<-c(normal_data,MCI_data,AD)
  label_group<-c(rep('CN',length(normal_data)),rep('MCI',length(MCI_data)),rep('AD',length(AD)))
  ##
  label_group<-factor(label_group,levels = c('CN','MCI','AD'))
  data_fram_lan<-data.frame(x=label_group,methylation=data_com)
  p<- ggplot(data_fram_lan,aes(x=x,y=methylation,fill=x))+geom_boxplot()+theme(title=element_text(,size=20))+xlab('Group')+ylab('Merhylation_bag')
  print(p)
}
### compute the association between mag and bag
AD_sample_data_bag<-AD_sample_data$AD_pred-AD_sample_data$Age
loc_image<-AD_sample_data$Loc%>%as.matrix()
loc_image_id<-apply(loc_image, 1, function(x){
  #x<-loc_image[1]
  x<-str_split(x,'[/]')%>%unlist()
  x<-str_split(x[6],'[_]')%>%unlist()
  x<-x[1]
  
})

###
bag_elas_group<-bag_elas_group[which(!duplicated(bag_elas_group$Sample_Name)),]
bag_blue_group<-bag_blue_group[which(!duplicated(bag_blue_group$Sample_Name)),]
rownames(bag_elas_group)<-bag_elas_group$Sample_Name
rownames(bag_blue_group)<-bag_blue_group$Sample_Name
bag_elas_group<-bag_elas_group[loc_image_id,]
bag_blues_group<-bag_blues_group[loc_image_id,]
AD_sample_data_mag<-bag_elas_group$elas_bag
AD_sample_data_mag_blue<-bag_blue_group$blue_bag
## age difference in AD patients ##
boxplot(AD_sample_data_bag,AD_sample_data_mag_blue)
cor.test(AD_sample_data_bag,AD_sample_data_mag)
##
ad_real_pred<-read.table('ADNI_methy_real.txt')
AD_sample_data_bag<-ad_real_pred$V1-AD_sample_data$Age
########## load the group methylation data ######
get_data<-function(label_name){
  #label_name='AD'
  data_group<-read.table(paste('./methyation_age/methylation_pred_age_ela_',as.character(label_name),'.txt',sep = ""))
  if(label_name=='normal'){
    label_name<-'CN'
  }else if(label_name=='eMCI'){
    label_name<-'EMCI'
  }
  data_group_info<-all_variable_data[which(all_variable_data$Sample_Group==label_name),]
  if(ncol(data_group)>3){
    data_group<-t(data_group)
  }
  bag_data<-as.vector(data_group)-data_group_info$age
  bag_data<-cbind(data_group_info[,c(3,10,11)],bag_data)
  colnames(bag_data)<-c('Sample_name','Age','Sex','MAG')
  return(bag_data)
}
####
normal_mag<-get_data('normal')
LMCI_mag<-get_data('LMCI')
EMCI_mag<-get_data('eMCI')
MCI_mag<-get_data('MCI')
AD_mag<-get_data('AD')
############ mag data plot ##
#AD_sample_data


all_methlation_age_data<-rbind(normal_mag,LMCI_mag,EMCI_mag,MCI_mag,AD_mag)
all_methlation_age_data<-all_methlation_age_data[!duplicated(all_methlation_age_data$Sample_name),]
rownames(all_methlation_age_data)<-all_methlation_age_data$Sample_name
inter_name<-intersect(loc_image_id,all_methlation_age_data$Sample_name)
all_methlation_age_data<-all_methlation_age_data[as.character(inter_name),]
##
AD_sample_data<-cbind(loc_image_id,AD_sample_data)
AD_sample_data<-AD_sample_data[which(!duplicated(AD_sample_data$loc_image_id)),]
rownames(AD_sample_data)<-AD_sample_data$loc_image_id
AD_sample_data_iner<-AD_sample_data[as.character(inter_name),]
##


##
#### 回归图 ##
group_data<-c(rep('MAG',nrow(all_methlation_age_data)),rep('BAG',nrow(AD_sample_data_iner)))
####  adjust the brain age #
adjust_fun_age<-function(brain_age){
 # AD_age<-AD_sample_data_iner$Age
  brain_age<-AD_sample_data_iner$AD_pred
  lm_data_fram<-data.frame(brain_age,AD_age)
  ##
  age_real<-lm(brain_age~AD_age,lm_data_fram)
  rsq(age_real)
  ##
  correct_age<-brain_age+(AD_age-(AD_age*age_real$coefficients[2]+age_real$coefficients[1]))
  ###
  lm_data_fram$cor_age<-correct_age

  ## compute the r square ##
  return(pred_real_gap)
}
##
value_data<-c(all_methlation_age_data$MAG,AD_sample_data_iner$AD_sample_data_bag)
data_bag_cor<-data.frame(MAG=all_methlation_age_data$MAG,BAG=AD_sample_data_iner$AD_sample_data_bag)
ggplot(data_bag_cor,aes(x=MAG,y=BAG,fill='red')) +geom_smooth(method='lm')+geom_point()+xlab('MAG')+ylab('BAG')+
  ggtitle('Correlation_bag_MAG_multi')+theme(title=element_text(,size=20))+theme_bw()+theme(panel.grid=element_blank())
####bag ##
normal_data<-AD_sample_data[which(AD_sample_data$Group=='CN'),8]
MCI_data<-AD_sample_data[which(AD_sample_data$Group=='MCI'),8]
LMCI_data<-AD_sample_data[which(AD_sample_data$Group=='EMCI'),8]
EMCI_data<-AD_sample_data[which(AD_sample_data$Group=='EMCI'),8]
AD_data<-AD_sample_data[which(AD_sample_data$Group=='AD'),8]
###
normal_data<-normal_mag$MAG
MCI_data<-MCI_mag$MAG
LMCI_data<-LMCI_mag$MAG
EMCI_data<-EMCI_mag$MAG
AD_data<-AD_mag$MAG

data_com<-c(normal_data,LMCI_data,EMCI_data,MCI_data,AD_data)
label_group<-c(rep('CN',length(normal_data)),rep('EMCI',length(EMCI_data)),rep('LMCI',length(LMCI_data)),
               rep('MCI',length(MCI_data)),rep('AD',length(AD_data)))
##
label_group<-factor(label_group,levels = c('CN','EMCI','LMCI','MCI','AD'))
data_fram_lan<-data.frame(x=label_group,methylation=data_com)
my_comparisons <- list( c("CN", "AD"), c("CN", "MCI"), c("MCI", "AD") )
data_fram_lan$methylation<--data_fram_lan$methylation
ggboxplot(data_fram_lan, x = "x", y = "methylation",
          color = "x", palette = "jco")+ 
  stat_compare_means(comparisons = my_comparisons,method ='t.test' )+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50)+     # Add global p-value
  theme(title=element_text(,size=20))+xlab('Group')+ylab('BAG_real')
##
graph2ppt(file='ADNI.pptx',width=8,height=6,append=TRUE)
p<- ggplot(data_fram_lan,aes(x=x,y=methylation,fill=x))+geom_boxplot()
######## using the 


