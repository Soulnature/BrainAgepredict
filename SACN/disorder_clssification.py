# -*- coding: utf-8 -*-
"""
Created on Sat Feb 26 16:25:19 2022

@author: admin
"""
from sklearn.svm import SVC
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import  auc,plot_roc_curve
from sklearn.model_selection import cross_val_score,StratifiedKFold
class SVM_model:
    def __init__(self,train_x,train_y,fold):
        super(SVM_model, self).__init__()
        self.x=train_x
        self.y=train_y
        self.fold=fold
    def svm_cross_validation(self):
        model = SVC()
        param_grid = {'C': [1e-3, 1e-2, 1e-1, 1, 10, 100, 1000], 'gamma': [1,0.1, 0.01, 0.001, 0.0001]}
        grid_search = GridSearchCV(model, param_grid, n_jobs=8, verbose=1)
        grid_search.fit(self.x, self.y)
        best_parameters = grid_search.best_estimator_.get_params()
        for para, val in list(best_parameters.items()):
            print(para, val)
        model = SVC(C=best_parameters['C'], gamma=best_parameters['gamma'], probability=False)
        tprs = []
        aucs = []
        mean_fpr = np.linspace(0, 1, 100)
        fig, ax = plt.subplots()
        cv=StratifiedKFold(n_splits=self.fold)
        for i,(train, test) in enumerate(cv.split(self.x, self.y)):
            model.fit(self.x[train], self.y[train])
            viz = plot_roc_curve(model, self.x[test], self.y[test],
                                 name='ROC fold {}'.format(i),
                                 alpha=0.3, lw=1, ax=ax)
            interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
            interp_tpr[0] = 0.0
            tprs.append(interp_tpr)
            aucs.append(viz.roc_auc)
        ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
                label='Chance', alpha=.8)

        mean_tpr = np.mean(tprs, axis=0)
        mean_tpr[-1] = 1.0
        mean_auc = auc(mean_fpr, mean_tpr)
        std_auc = np.std(aucs)
        ax.plot(mean_fpr, mean_tpr, color='b',
                label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
                lw=2, alpha=.8)
        std_tpr = np.std(tprs, axis=0)
        tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
        tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
        ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                        label=r'$\pm$ 1 std. dev.')

        ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05],
               title="PD_classification")
        ax.legend(loc="lower right")
        plt.savefig('PD_classification.pdf')
        plt.show()
        predicted_socre=str(np.mean(mean_auc))+'±'+str(np.std(mean_auc))
        return predicted_socre
    def check_data(self):
        isnan_x = np.isnan(self.x)  # 判断每个元素是不是nan,返回[False,False,False,False,True]
        isnan_y = np.isnan(self.y)
        if True in isnan_x or True in isnan_y:
            self.x[np.isnan(self.x)] = 0
            self.y[np.isnan(self.y)] = 0
            print('success in check data')
def get_data(label_files,disorder_fea,normal_fea,taskname):

    if taskname=='PPMI':
        all_label_sample = pd.read_csv(label_files)
        all_lable = all_label_sample['Group'].values
        all_feathre=np.load(disorder_fea)
        x_train = all_feathre[:, :, 0]
        y_train = all_lable
    else:
        all_label_sample = pd.read_csv(label_files)
        ad_label = all_label_sample.loc[all_label_sample['Group'] == 2]
        normal_label = all_label_sample.loc[all_label_sample['Group'] == 0]
        disorder_label = ad_label['Group'].values
        normal_label = normal_label['Group'].values
        ad_fea = np.load(disorder_fea)
        ad_fea = ad_fea[:, :, 0]
        normal_fea = np.load(normal_fea)
        normal_fea = normal_fea[:, :, 0]
        print(normal_fea.shape)
        x_train = np.vstack((ad_fea, normal_fea))
        y_train = np.hstack((disorder_label, normal_label))
    return x_train,y_train
if __name__=='__main__':
    label_files='/AD_data_input_unsweed.csv'
    disorder_fea='/AD_input_fea.npy'
    normal_fea='/AD_nor_fea.npy'
    train_x,train_y=get_data(label_files,disorder_fea,normal_fea,'AD')
    model_data=SVM_model(train_x,train_y,10)
    score_metrics=model_data.svm_cross_validation()
    print(score_metrics)
