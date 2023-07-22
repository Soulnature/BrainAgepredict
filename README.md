## Induction
The code is used for the work of “Deciphering genetic architecture of brain age gap with functional genomic analysis”.
This code include two parts，the first part is the SACN model, and second part is the analysis process for the BAG GWAS.

### SACN mode 
The code for the SACN model, train process and evalution process.
- model.py, Resnet.py. the network architecture of SACN and resnet
- Evalation.py, the evalution process for SACN
- train_c.py, the train and test process
### Analysis
- PRS.R. the PGS analysis of different phenotype, and association with BAG
- MR_plot.R, the plot for MR result
- Expression_analysis.R. The expression analysis of BAG-related gene
- AD_methylation.R, the process of methylation data on ADNI
- TRN_bag.R. the associated analysis of BAG-related gene and MAG-related gene in brain TRN.
- TRN_aging.R. the association between BAG and aging-relate pathway
- Disorder_common.R.  the pleiotropy of BAG-related gene in diverse brain disorders
- PPI.R. the  PPI analysis of BAG-related gene

## Contact information
if you have any question about this code, please contact the author (18210850006@fudan.edu.cn)
