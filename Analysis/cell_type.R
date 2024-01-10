library(data.table)
library(stringi)
library(stringr)
library(ggpubr)
library(export)
library(Matrix)


resd_files_of_first_row<-function(rootdir,metadata_fan){
  # return the tissue df and cell df  
  # rootdir<-'./FANTOM/'
  # metadata_fan<-'./FANTOM.csv'
  all_cell_type<-list.files(rootdir)
  tem_file<-data.frame()
  for (i in all_cell_type) {
    #i<-all_cell_type[1]
    tem_file_data<-read.table(paste(rootdir,i,sep = ""),header = TRUE)
    tem_file<-rbind(tem_file,tem_file_data[1,])
  }
  cell_name<-c() 
  for (j in all_cell_type) {
    # j<-all_cell_type[1]
    cell_tem<-str_split(j,'_')%>%unlist()
    cell_name<-c(cell_name, cell_tem[3])
  }
  tem_file$Category<-cell_name
  all_cell<-tem_file[1:6,]
  tissue<-tem_file[7:nrow(tem_file),]
  ###
  tissue$FDR<- -log(p.adjust(tissue$Enrichment_p,method = 'fdr'))
  all_cell$FDR<- -log10(p.adjust(all_cell$Enrichment_p,method = 'fdr'))
  ########### #############
  meta_data<-read.csv(metadata_fan)
  cell_annotation<-meta_data[which(meta_data$SampleType=='primary cell'),]
  tissue_annot<-meta_data[which(meta_data$SampleType=='tissue'),]
  tisse_name<-apply(tissue_annot, 1, function(x){
    tissue_names<-str_split(x[2],', ')%>%unlist()
    x<-tissue_names[1]
  })
  periods<-apply(tissue_annot, 1, function(x){
    tissue_names<-str_split(x[2],', ')%>%unlist()
    x<-tissue_names[2]
  })
  #####
  tissue_df<-data.frame(state_id=tissue_annot$StateID,tissue_name=tisse_name,period=periods)
  tissue_time<-table(tissue_df$tissue_name)%>%as.data.frame()
  tissue_time<-tissue_time[which(tissue_time$Freq==2),]
  ###
  tissue_df_two_per<-tissue_df[which(tissue_df$tissue_name%in%tissue_time$Var1),]
  tissue$Category<-apply(as.matrix(tissue$Category),1,function(x){
    x<-str_split(x,'[.]')%>%unlist()
    x<-paste('State',x[1],sep = '_')
  })
  ####
  rownames(tissue)<-tissue$Category
  tissue<-tissue[tissue_df_two_per$state_id,]
  
  
  #### combin the information with the original data ##
  tissue<-cbind(tissue,tissue_df_two_per[,2:3])
  tissue<-tissue[order(tissue$Enrichment_p,decreasing =TRUE),] 
  ##
  tissue<-tissue[which(tissue$period=='newborn'|tissue$period=='adult'),]
  #
  tissue<-tissue[which(tissue$tissue_name!='brain'),]
  tissue<-tissue[which(tissue$tissue_name!='occipital lobe'),]
  tissue<-tissue[which(tissue$tissue_name!='temporal lobe'),]
  p1 <- ggplot(tissue,aes(x=tissue_name,
                          y=ifelse(period=="adult",FDR,-FDR)))+coord_flip()
  p2 <- p1+geom_bar(aes(fill=period),colour="white",
                    size=0,width = 0.9,stat = "identity")

  p2 <- p2+labs(fill="",x="",y="FDR",title="")+guides(color = FALSE)
  p3<-p2+scale_y_continuous(expand=expansion(add = c(0.1, 0.1)),
                            limits = c(-5, 5),
                            breaks = c(-4,-3,-2,-1,0,1,2,3,4),
                            label = c('4',"3","2", "1","0","1","2","3",'4'))+ylab('-log10(FDR)')
  p3
  mycolor <- c("#FF9900","#9933FF")
  p4 <- p3 + scale_fill_manual(values=alpha(mycolor,0.7))
  top.mar=0.2
  right.mar=0.1
  bottom.mar=0.2
  left.mar=0.1
  
  mytheme<-theme_classic()+
    theme(text=element_text(family = "sans",colour ="gray30",size = 12),
          axis.text.y = element_text(size = 12),
          axis.line.y=element_blank(),
          axis.ticks.y = element_blank(),
          axis.line = element_line(size = 0.6,colour = "gray30"),
          axis.ticks = element_line(size = 0.6,colour = "gray30"),
          axis.ticks.length = unit(1.5,units = "mm"),
          plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                           units="inches"))
  p5<-p4+mytheme+ geom_hline(yintercept = 1.3, linetype = "dashed", color = "red", size = 0.8)+
    geom_hline(yintercept = -1.3, linetype = "dashed", color = "red", size = 0.8)
  ### other function regions of the bag #
  others_res<-read.csv('./conder_part_res.csv',header = TRUE)
  ###
  ggbarplot(others_res,x='Category',y='FDR',sort.val = 'asc',fill = 'grey')+
    coord_flip()+  geom_hline(yintercept = 1.3, linetype = "dashed", color = "red", size = 0.8)+
    ylab('-log10(FDR)')
  
  others_res$sig<-ifelse( others_res$FDR>1.3,'Sig','No sig')
  pd <- position_dodge(0.1) # move them .05 to the left and right
  ggplot(others_res, aes(x=Category, y=Enrichment,fill=sig)) + 
    geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error), colour="black", width=.1, position=pd) +
    geom_point(position=pd, size=5, shape=21) + # 21 is filled circle
    geom_line()+
    xlab("brain_region") +
    ylab("Enrichment") +
    ggtitle("Partitioning of the SNP heritability of the BAG GWA by brain region annotations") +
    theme_bw() +theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    theme(legend.justification=c(1,0),# 
          legend.position=c(1,0)) +              # Position legend in bottom right
    coord_flip()+
    scale_fill_brewer(palette = "Set1")+
    geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed")
  
  
  return(list(Tissue = tissue, Cell = all_cell))
}

cell_tiisue<-resd_files_of_first_row('./FANTOM/','./FANTOM.csv')
graph2ppt(file='par_h2_FANTOM.pptx',height=8,width=10,append=TRUE)

##Single cell par ##
singlecell_root<-list.files('./single_cell_adult_brain_enhan/')
single_cell_en<-data.frame()
for(i in singlecell_root){
  # x<-singlecell_root[1]
  name_data<-str_split(i,'1_')%>%unlist()
  name_data<-str_split(name_data[2],'[.]')%>%unlist()
  name_data<-name_data[1]
  x_data<-read.table(paste('./single_cell_adult_brain_enhan/',i,sep = ""),header = TRUE)
  x_data<-x_data[1,]
  x_data$Category<-name_data
  single_cell_en<-rbind(single_cell_en,x_data)
  
}

single_cell_en$FDR<--log10(p.adjust(single_cell_en$Enrichment_p,method = 'BH'))
####
rownames(single_cell_en)<-single_cell_en$Category
annot_single_cell<-read.csv('./single_cell_class_annotation.csv',fill = TRUE)
annot_single_cell<-annot_single_cell[which(!duplicated(annot_single_cell$Cell.subclass)),]
annot_single_cell$Cell.subclass[order(annot_single_cell$Cell.subclass)]
rownames(annot_single_cell)<-annot_single_cell$Cell.subclass
using_single_cell<-intersect(annot_single_cell$Cell.subclass,single_cell_en$Category)
single_cell_en<-single_cell_en[using_single_cell,]
annot_single_cell<-annot_single_cell[using_single_cell,]
###
single_cell_en$color<-annot_single_cell$Cell.class.name
single_cell_en$clf<-annot_single_cell$Cell.class.color
##
c<-single_cell_en[which(single_cell_en$Enrichment>0),]

p<-ggbarplot(single_cell_en,x='Category',y='FDR',sort.val = 'asc',fill = 'color')+
  coord_flip()+  geom_hline(yintercept = 1.3, linetype = "dashed", color = "red", size = 0.8)+
  ylab('-log10(FDR)')+scale_fill_manual(values=unique(single_cell_en$clf))
###
p+geom_line(aes(y = mean_g, group = color), 
) # 
mean_enrich<-tapply(single_cell_en$FDR, single_cell_en$color, mean)%>%as.data.frame()
single_cell_en$mean_g<-mean_enrich[single_cell_en$color,]

