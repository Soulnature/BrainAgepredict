## Induction
The code is used for the work of “Deciphering genetic architecture of brain age gap”.
This code include two parts，the first part is the SACN model, and second part is the analysis process for the BAG GWAS.

### SACN model 
The code for the SACN model, train process and evalution process.
- SACN.py, Resnet.py. the network architecture of SACN and resnet
- Evalation.py, the evalution process for SACN
- fune_tune.py, transfer the pretrained  network to the UKB biobank
- train_c.py, the train and test process
- main.py, the main process for paramater modify
- requirements.txt the environment required to run this model


You can run  following command to start the train process
```python
python main.py --root_dir /public/home/zhaoxz/reCOns/full_age1.csv --network_type SACN --model_name SACN --base SACN --epochs 200 --batch_size 12  --lr 1e-04
# network_type one of [base,Gender_AD,Domain_AD,SACN]
```
### Analysis
- PRS.R. the PGS analysis of different phenotype, and association with BAG
- MR_plot.R, the plot for MR result
- Expression_analysis.R. The expression analysis of BAG-related gene
- AD_methylation.R, the process of methylation data on ADNI
- TRN_bag.R. the associated analysis of BAG-related gene and MAG-related gene in brain TRN.
- TRN_aging.R. the association between BAG and aging-relate pathway
- Disorder_common.R.  the pleiotropy of BAG-related gene in diverse brain disorders
- PPI.R. the PPI analysis of BAG-related gene
- cell_type.R the code for cell type analysis

## Contact information
if you have any question about this code, please contact the author (18210850006@fudan.edu.cn)
