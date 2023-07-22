## bag data ##
library(dplyr)
bag_pheno<-read.csv('UKB_pheno_type/UKB_match.csv',header = TRUE)
prs_cov<-read.table('./PRS/bag_prs_cov.txt')
rownames(bag_pheno)<-bag_pheno$UKB_image_id
bag_pheno_c<-bag_pheno[as.character(prs_cov$V1),]
#####
bag_pheno_c_data<-cbind(prs_cov,bag_pheno_c$Age,bag_pheno_c$Gende)%>%na.omit()
###
write.table(bag_pheno_c_data,'./PRS/bag_prs_cov.txt',row.names = FALSE,col.names = FALSE,quote = FALSE)
#### Bag prs phenotype set  ##
bag_pheno<-read.table('./BAG_phenotype.txt',header = TRUE,sep = '\t')
write.table(bag_pheno,'BAG_phenotype.txt',row.names = FALSE,quote = FALSE,sep = '\t')
