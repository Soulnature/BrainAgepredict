## MR ##
library(data.table)
library(dplyr)

library(forestploter)
MR_exp_files<-fread('TwoSampleMR.BAG_exp.5E-8.all.txt')
MR_out_files<-fread('TwoSampleMR.BAG.out.ld_0.1.5E-8_ld0.1.all.txt')
###
list_name_exp<- unique(MR_exp_files$outcome)
list_name_out<- unique(MR_out_files$exposure)
list_name_exp<-list_name_exp[c(2,3,5,9,12,14,15,18,19,20,22,24,25,27,28,29,30)]
list_name_out<-list_name_out[c(1,2,4,8,10,12,13,15,16,17,19,22,21,24,25,26,27)]
###
sample_mr_exp<-MR_exp_files[1,]
sample_mr_out<-MR_out_files[1,]
for (i in list_name_exp) {
 # temp_d<-MR_files[which(MR_files$outcome==i),]
  temp_d<-MR_exp_files[which(MR_exp_files$outcome==i),]
  sample_mr_exp<-rbind(sample_mr_exp,temp_d)

}
for (i in list_name_out) {
  # temp_d<-MR_files[which(MR_files$outcome==i),]
  temp_d<-MR_out_files[which(MR_out_files$exposure==i),]
  sample_mr_out<-rbind(sample_mr_out,temp_d)
  
}
sample_mr_exp<-sample_mr_exp[-1,]
sample_mr_out<-sample_mr_out[-1,]
########
FDR_exp<-p.adjust(sample_mr_exp$pval,method = 'BH')
FDR_out<-p.adjust(sample_mr_out$pval,method = 'BH')
sample_mr_exp<-cbind(sample_mr_exp,FDR_exp)
sample_mr_out<-cbind(sample_mr_out,FDR_out)
###sort the value
sample_mr_exp<-sample_mr_exp[order(sample_mr_exp$FDR_exp),]
sample_mr_out<-sample_mr_out[order(sample_mr_out$FDR_out),]

###
write.csv(sample_mr_exp,'BAGtoDisorders.csv',row.names = FALSE)
write.csv(sample_mr_out,'DisordertoBAG.csv',row.names = FALSE)
### plot the sig result ##
sample_mr_out<-sample_mr_out[which(sample_mr_out$FDR_out<=0.05),]
sample_mr_chose<-sample_mr_out[,c(2,3,4,5,6,7,10)]
low_ci<-sample_mr_chose$b-sample_mr_chose$se
up_ci<-sample_mr_chose$b+sample_mr_chose$se
###
sample_mr_chose<-cbind(sample_mr_chose,low_ci,up_ci)
sample_mr_chose$` ` <- paste(rep(" ", 15), collapse = " ")
##
all_disease<-sample_mr_chose$exposure%>%unique()
all_tem<-sample_mr_chose[1,]
for (i in all_disease) {
  tem_d<-sample_mr_chose[which(sample_mr_chose$exposure==i),]
  all_tem<-rbind(all_tem,tem_d)
}
all_tem<-all_tem[-1,]
write.csv(all_tem,'all_res1.csv',row.names = FALSE)
all_tem<-read.csv('all_res1.csv')
all_tem$` ` <- paste(rep(" ", 8), collapse = " ")
tm <- forest_theme(base_size = 10,
                   # Confidence interval point shape, line type/color/width
                   ci_pch = 16,
                   ci_col = "#762a83",
                   ci_lty = 1,
                   ci_lwd = 1.5,
                   ci_Theight = 0.2, # Set an T end at the end of CI 
                   # Reference line width/type/color
                   refline_lwd = 1,
                   refline_lty = "dashed",
                   refline_col = "grey20",
                   # Vertical line width/type/color
                   vertline_lwd = 1,
                   vertline_lty = "dashed",
                   vertline_col = "grey20",
                   # Change summary color for filling and borders
                   summary_fill = "#4575b4",
                   summary_col = "#4575b4",
                   # Footnote font size/face/color
                   footnote_cex = 0.6,
                   footnote_fontface = "italic",
                   footnote_col = "blue")

###
#####
 forest(all_tem[, c(1, 2, 3,4,11)],
               est = all_tem$b,
               lower = all_tem$low_ci,
               upper = all_tem$up_ci,
               ci_column = 5,
               ticks_at = c(-1,-0.5,0,0.5,1),
               ticks_digits=2,
               theme = tm)


