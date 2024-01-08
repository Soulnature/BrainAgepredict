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
    #设置图表的图例、坐标轴、标题；
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
      theme(legend.justification=c(1,0),# 这一项很关键,如果没有这个参数,图例会偏移,读者可以试一试
            legend.position=c(1,0)) +              # Position legend in bottom right
      coord_flip()+
      scale_fill_brewer(palette = "Set1")+
      geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed")
   #others_res<-others_res[-c(1,2),]
   #others_res$FDR<--log10(p.adjust(others_res$Enrichment_p,method = 'BH'))
    #############
   #write.csv(others_res,'conder_part_res.csv',row.names = FALSE,quote = FALSE)
    
    
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
single_cell_en$FDR<--log10(single_cell_en$Enrichment_p)
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






### deal the tissue label #
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
#trajetory_cell_tyep[which(trajetory_cell_tyep$period==14),2]<-13
#trajetory_cell_tyep[which(trajetory_cell_tyep$period==15),2]<-14
#trajetory_cell_tyep$period<-factor(trajetory_cell_tyep$period)
ggplot(data=trajetory_cell_tyep, 
       aes(x=period, y=exp, group=group, color=factor(group),
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
####



####

#

