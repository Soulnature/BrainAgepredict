# PRS data preprocess ##
library(stringi)
library(stringr)
library(ggpubr)
library(export)
PRS_dir<-list.files('./PRS_summary/')
data_spec<-data.frame()
for (i in PRS_dir) {
  # i<-PRS_dir[1]
  tem_data<-read.table(paste('./PRS_summary/',i,sep = ""),header =TRUE)
  data_spec<-rbind(data_spec,tem_data)
}
data_spe<-data_spec

data_spe<-cbind(PRS_dir,data_spe)
## deal the data ##
data_spe$PRS_dir<-apply(as.matrix(data_spe$PRS_dir),1,function(x){
  x<-str_split(x,'[.]')%>%unlist()
  x<-x[1]
})
data_spe$P_adjust<-p.adjust(data_spe$P,method = 'BH')
data_spe$P_adjust<--log10(data_spe$P_adjust)
write.csv(data_spe,'PRS.csv',row.names = FALSE)
###
data_spe$P_adjust[10]<-1.01
## add the disorder classification ##
disorder_class<-c('Cognitive-behavioral', 'Neurological disorder','Cognitive-behavioral','Psychotic disorder','aging disorder',
                  'Others','Neurological disorder','Cognitive-behavioral','Psychotic disorder','Psychotic disorder',
                  'Psychotic disorder','Neurological disorder','Cognitive-behavioral',
                  'aging disorder','Cognitive-behavioral','Others','Neurological disorder' )



##
data_spe_order<-data_spe[order(data_spe$P_adjust,decreasing = F),]
data_spe_order$disorder_type<-rev(disorder_class)
mycolor <- c("#A5562D","#E4191C",'#964EA2','#F681BF','#000000')

p<-ggbarplot(data_spe_order, x = "PRS_dir", y = "P_adjust",
             fill  = "disorder_type", 
             sorting ="descending",                  
             add = "segments",                             
             xlab="Traits",
             rotate = TRUE,
             ylab='-log10(FDR)')+
  geom_hline(aes(yintercept=1.3), colour="#990000", linetype="dashed")
p1 <- p + scale_fill_manual(values=alpha(mycolor,0.7))

###
graph2ppt(file="./bag_gc.pptx",append=TRUE,height=8,width=10)
###
gc_files<-list.files('./genetic_correlation/')
test_data_d<-data.frame()
for (i in gc_files) {
  tem_gc<-read.table(paste('./genetic_correlation/',i,sep = ''),fill = TRUE)
  tem_gc<-tem_gc[c(54,55),]
  colnames(tem_gc)<-tem_gc[1,]
  tem_gc<-tem_gc[-1,]
  test_data_d<-rbind(test_data_d,tem_gc)
}
#####
test_data_d<-cbind(test_data_d,gc_files)
test_data_d$gc_files<-apply(as.matrix(test_data_d$gc_files),1,function(x){
  # x<-test_data_d$gc_files[1]
  x<-str_split(x,'[.]')%>%unlist()
  x<-str_split(x[1],'[__]')%>%unlist()
  x<-x[3]
})
##
#####
test_data_d$FDR<--log10(p.adjust(test_data_d$p,method = 'BH'))
test_data_d$se<-as.numeric(test_data_d$se)
test_data_d$rg<-as.numeric(test_data_d$rg)

test_data_d$label_sig_p<-label_sig_p%>%as.character()
#

###
test_data_d$gc_files<-factor(test_data_d$gc_files,levels =data_spe_order$PRS_dir)
test_data_d$disorder_type<-rev(disorder_class)
mycolor <- c("#A5562D","#E4191C",'#964EA2','#F681BF','#000000')
pd <- position_dodge(0.1) # move them .05 to the left and right
#####
ggplot(test_data_d, aes(x=gc_files, y=rg,fill=disorder_type)) + 
  geom_errorbar(aes(ymin=rg-se, ymax=rg+se), colour="black", width=.1, position=pd) +
  geom_point(position=pd, size=6, shape=21) + # 21 is filled circle
  geom_line()+
  xlab("brain_region") +
  ylab("Enrichment") +
  ggtitle("Partitioning of the SNP heritability of the BAG GWA by brain region annotations") +
  theme_bw() +theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(legend.justification=c(1,0),# 
        legend.position=c(1,0)) +              # Position legend in bottom right
  coord_flip()+
  scale_fill_manual(values=alpha(mycolor,0.7))+
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed")
#####

###




