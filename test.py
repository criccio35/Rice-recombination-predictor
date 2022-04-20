#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 12:02:50 2022

@author: camila
"""

import rice_recombination_predictor as rrp
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

size_w = 100_000
model_chr_nbr = '01'

ref_path = 'input_data/IR64/Osat_IR64_chr{0}.fasta'.format(model_chr_nbr) 
qry_path = 'input_data/Azucena/Osat_Azucena_chr{0}.fasta'.format(model_chr_nbr) 
coords_path = 'input_data/coords/IR64_Azucena_chr{0}.coords'.format(model_chr_nbr)
variants_path = 'input_data/snps/IR64_Azucena_chr{0}.snps'.format(model_chr_nbr)
rec_path = 'input_data/recombination/experimental_recombination_chr{0}.csv'.format(model_chr_nbr)
CentO_path = 'input_data/CentO_AA.fasta'
results_path = 'output_data/'

chr_mod = rrp.rice_recombination_predictor(size_w, ref_path, qry_path, 
                                          coords_path, variants_path, 
                                          CentO_path, rec_path)

chr_mod.data_preprocessing()

initial_parameters = [0.5,0.97,1,0.9,1,0,0.002] #(p1,p2,p3,p4,p5,p5,p7)
chr_mod.optimize_model_parameters(initial_parameters)
chr_mod.prediction

chr_mod.plot_landscape()
chr_mod.plot_correlation()

chromosomes = ['01','02','03','04','05','06','07','08','09','10','11','12']

df_pred = pd.DataFrame()
df_eval = pd.DataFrame(columns=['chromosome','R2_idt','R2_pred','pearson_idt','pearson_pred'])

for chr_nbr in chromosomes:
    ref_path = 'input_data/IR64/Osat_IR64_chr{0}.fasta'.format(chr_nbr) 
    qry_path = 'input_data/Azucena/Osat_Azucena_chr{0}.fasta'.format(chr_nbr) 
    coords_path = 'input_data/coords/IR64_Azucena_chr{0}.coords'.format(chr_nbr)
    variants_path = 'input_data/snps/IR64_Azucena_chr{0}.snps'.format(chr_nbr)
    rec_path = 'input_data/recombination/experimental_recombination_chr{0}.csv'.format(chr_nbr)
    
    chr_tmp = rrp.rice_recombination_predictor(size_w, ref_path, qry_path, 
                                                coords_path, variants_path, 
                                                CentO_path, rec_path,
                                                params=chr_mod.params)
    
    print('---Chromosome {0}---'.format(chr_nbr))
    chr_tmp.data_preprocessing()
    chr_tmp.predict_recombination()
    
    df_pred1 = pd.DataFrame({'chr_'+chr_nbr : chr_tmp.prediction})
    df_pred = pd.concat([df_pred,df_pred1], axis=1)

    # Plot landscape and correlation
    fig, ax = plt.subplots(1,2,figsize=(20,5),gridspec_kw={'width_ratios': [3, 1]})
    chr_tmp.plot_landscape(ax=ax[0])
    chr_tmp.plot_correlation(ax=ax[1])
    fig.suptitle('Prediction of chromosome {0} recombination, calibrated with chromosome {1}'.format(chr_nbr,model_chr_nbr),
                 fontsize=15)
    
    # model evaluation
    corr_mod, r2_mod = chr_tmp.prediction_evaluation('model')
    corr_idt, r2_idt = chr_tmp.prediction_evaluation('identity')
    df_eval = df_eval.append({'chromosome':chr_nbr,
                                'R2_idt':r2_idt,
                                'R2_pred':r2_mod,
                                'pearson_idt':corr_idt,
                                'pearson_pred':corr_mod},
                             ignore_index=True)

df_pred.to_csv(results_path+'predictions.csv',header=True,index=True)

# Plot model evaluation (R2)
plt.figure(figsize=(10,5))
ax0 = sns.barplot(data=df_eval,x='chromosome',y='R2_pred',color='tab:blue',label='model')
ax0.bar_label(ax0.containers[0],label_type='edge',padding=10,fmt='%.3f',fontsize=10)
ax1 = sns.barplot(data=df_eval,x='chromosome',y='R2_idt',color='lightgray',label='identity')
ax1.bar_label(ax1.containers[1],label_type='edge',padding=-20,fmt='%.3f',fontsize=10)
ax1.set_xlabel('Tested chromosome')
ax1.set_ylabel('$R^2$')
ax1.set_title('$R^2$ of the model calibrated on chromosome {0}'.format(model_chr_nbr),
              fontsize=15)
ax1.legend(loc='lower right',framealpha=1)
plt.ylim(None,2.5)
    
# Plot model evaluation (pearson)
plt.figure(figsize=(10,5))
ax0 = sns.barplot(data=df_eval,x='chromosome',y='pearson_pred',color='tab:blue',label='model')
ax0.bar_label(ax0.containers[0],label_type='edge',padding=10,fmt='%.3f',fontsize=10)
ax1 = sns.barplot(data=df_eval,x='chromosome',y='pearson_idt',color='lightgray',label='identity')
ax1.bar_label(ax1.containers[1],label_type='edge',padding=-20,fmt='%.3f',fontsize=10)
ax1.set_xlabel('Tested chromosome')
ax1.set_ylabel('correlation (r)')
ax1.set_title('Correlations of the model calibrated on chromosome {0}'.format(model_chr_nbr),
              fontsize=15)
ax1.legend(loc='lower right',framealpha=1)
plt.ylim(None,1)

