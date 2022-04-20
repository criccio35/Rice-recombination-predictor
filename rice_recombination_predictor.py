#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 27 09:42:21 2021

@author: Camila Riccio and Mauricio Peñuela
"""

# Libraries
import pandas as pd
import numpy as np
import CentOFinder as cf
import os
from scipy.stats import pearsonr
from sklearn.metrics import r2_score
from scipy.optimize import minimize

import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
# # Plotting, own library
#from .plots import *

class rice_recombination_predictor():
    '''
    Class for predicting chromosomal recombination in rice (Oryza sativa) 
    from the alignment of the two parental sequences.
    
    :param size_w: Observation window size (base pairs)
    :type size_w: int
    :param ref_path: fasta file path with a chromosome sequence of the 
        reference genome.
    :type ref_path: str
    :param qry_path: fasta file path with a chromosome sequence of the 
        qry genome. The sequence must be of the same chromosome number 
        as ref_path.
    :type qry_path: str
    :param coords_path: file path with extension .coords that is obtained 
        from MUMmer alignment, containing the contig coordinates.
    :type coords_path: str
    :param variants_path: file path with extension .snps that is obtained 
        from MUMmer alignment, containing the contig coordinates.
    :type variants_path: str
    :param CentO_path: fasta file containing the CentO sequence.
    :type CentO_path: str
    :param rec_path: file path with experimental recombination values.
    :type rec_path: str
    :param coords: contig coordinates of the alignment between the reference 
        genome and the query genome.
    :type coords: DataFrame.
    :param variants: variants detected in the alignment between the reference 
        genome and the query genome.
    :type variants: DataFrame.
    :param experimental_rec: Experimental recombination values per window.
    :type experimental_rec: list
    :param wbp_r: reference chromosome bins/windows in base pairs.
    :type wbp_r: list of ints
    :param wbp_q: query chromosome bins/windows in base pairs.
    :type wbp_q: list of ints
    :param features: relative frequency of the following features by window:
        absent bases, bases in inversions, variant bases, identical bases.
    :type features: DataFrame
    :param params: recombination prediction model parameters [p1,p2,p3,p4,p5,p6,p7].
    :type params: list of floats
    :param size_wc: number of windows used to the left and right of the 
        reference and query centromeres to make the linear transition
        from zero to one of the weight function that corrects for 
        recombination in the centromeric region.
    :type size_wc: int
    :param alpha: smoothing factor, 0 < alpha <=1. The lower the smoother.
    :type alpha: float
    :param cases: listing of the case number (0,1,2 or 3) that was applied
        in each window for the prediction of recombination, being 0 the
        non-application of cases.
    :type cases: list of ints
    :param prediction: predicted recombination values per window.
    :type prediction: list of floats
    '''
    
    # Default constructor
    def __init__(self, size_w, ref_path, qry_path,
                 coords_path, variants_path, CentO_path, rec_path=None,
                 coords=None, variants=None, experimental_rec=None,
                 wbp_r = None, wbp_q = None,
                 features = None, params = None,
                 size_wc=50, alpha = 0.1, cases = None, prediction=None):
        
        self.size_w = size_w
        self.ref_path = ref_path
        self.qry_path = qry_path
        self.coords_path = coords_path
        self.variants_path = variants_path
        self.CentO_path = CentO_path
        self.rec_path = rec_path
        self.coords = coords
        self.variants = variants
        self.experimental_rec = experimental_rec
        self.wbp_r = wbp_r
        self.wbp_q = wbp_q
        self.features = features
        self.params = params
        self.size_wc = size_wc
        self.alpha = alpha
        self.cases = cases
        self.prediction = prediction

        
    def _process_coordinates_file(self):
        '''
        This function reads the coordinate file, with extension .coords 
        obtained with the MUMmer program, extracts only the relevant 
        information for recombination prediction, and saves it in a csv format.
        '''
        final = []
        with open(self.coords_path) as file:
            a = file.readlines()
            head = a[3]
            # remove unnecessary characters
            new_head = head.replace(' ','')         
            new_head = new_head.replace('|','')     
            new_head = new_head.replace('[','')
            new_head = new_head.split(']')
            new_head = '\t'.join(new_head)
        
            a2 = a[5:]
            for line in a2:
                c = line.replace('|', '')
                c = c.split()
                c = '\t'.join(c)
                final.append(c)
           
        # export pre-processed file  
        path_head, path_tail = os.path.split(self.coords_path)
        new_coords_path = path_head + '/processed_coords_' + path_tail.split('.')[0] + '.csv' 
        out = open(new_coords_path, "w")
        out.write(new_head)                     
        for line in final:    
            # write line to output file
            out.write(line)
            out.write("\n")
        out.close()
        
        # read new_coords as dataframe
        coords = pd.read_csv(new_coords_path,sep='\t',usecols=['S1', 'E1', 'S2', 'E2', 'LEN1', 'LEN2','%IDY'])
        coords.rename(columns={'S1': 'Sr', 'E1': 'Er', 'S2': 'Sq', 'E2': 'Eq', 'LEN1':'len_r', 'LEN2':'len_q'}, inplace=True)
        coords['inversion_q'] = coords.apply(lambda df: (df.Sq >= df.Eq)*1, axis=1)
        coords = coords.sort_values(by=['Sr','Er'])
        coords = coords.reset_index()
        coords = coords.drop('index',axis=1)
        
        coords.to_csv(new_coords_path,header=True,index=True)
        self.coords = coords
        
    def _process_variants_file(self):
        '''
        This function reads the variants file, with extension .snps 
        obtained with the MUMmer program, extracts only the relevant 
        information for recombination prediction, and saves it in a csv format.
        '''
        final = []
        with open(self.variants_path) as file:
            a = file.readlines()
            head = a[3]
            # remove unnecessary characters
            new_head = head.replace(' ','')
            new_head = new_head.replace('[','')
            new_head = new_head.split(']')
            new_head = ''.join(new_head)
            new_head = new_head.replace('\n','\tTAGS\tTAGS\n')
            final = a[4:]
               
        # export pre-processed file  
        path_head, path_tail = os.path.split(self.variants_path)
        new_variants_path = path_head + '/processed_variants' + path_tail.split('.')[0] + '.csv' 
        out = open(new_variants_path, 'w')
        out.write(new_head)
        for line in final:
            # write line to output file
            out.write(line)
            out.write("\n")
        out.close()
        
        variants = pd.read_csv(new_variants_path,sep='\t',
                        usecols=['P1', 'SUB', 'SUB.1', 'P2'])
        variants.rename(columns={'P1':'rpos', 'SUB':'rbase', 'SUB.1':'qbase', 'P2':'qpos'}, inplace=True)
        variants = variants.drop_duplicates()
        variants = variants.reset_index(drop=True)
        
        variants.to_csv(new_variants_path,header=True,index=True)
        self.variants = variants
            
    def _read_sequence(self,path):
        '''
        Reads a fasta file with multiple rows and  
        returns the whole sequence on a single string.
        
        :param path: path where the fasta file with the genomic sequence is located.
        :type path: str
        
        :return: genomic sequence
        :rtype: str
        '''
        with open(path,'r') as f:
            f.readline()  
            seq = ''
            for line in f.readlines():
                seq += line.strip('\n').upper()
        return seq
    
    
    def _process_sequences(self):
        '''
        Reads the reference and query genomes and divedes them into windows of size size_w.
        '''
        self.ref_seq = self._read_sequence(self.ref_path)
        self.qry_seq = self._read_sequence(self.qry_path)    
        self.wbp_r = [i for i in range(0,self.size_w*(len(self.ref_seq)//self.size_w)+1,self.size_w)]
        self.wbp_q = [i for i in range(0,self.size_w*(len(self.qry_seq)//self.size_w)+1,self.size_w)]
    
    
    def _process_experimental_recombination(self):
        '''
        Reads a file with a .csv extension containing experimental windowed 
        recombination values. Fill the windows with missing data with nan 
        and smooth all values.
        '''
        df = pd.read_csv(self.rec_path)
        exp_rec = df.recombination.tolist()
        
        n = len(self.wbp_r)-1
        
        #expand real recombination if incomplete
        if len(exp_rec) < n:
            exp_rec += [np.nan]*(n-len(exp_rec)) 
        elif len(exp_rec) > n:
            exp_rec = exp_rec[:-(len(exp_rec)-n)]
        
        self.experimental_rec = self._exponential_smoothing(exp_rec)       
    
    
    def data_preprocessing(self):
        '''
        Perform pre-processing of the following files:
        coordinates with extension .coords,
        variants with extension .snps,
        sequences (reference and query) with extension .fasta, and 
        experimental recombination with extension .csv.
        '''
        self._process_coordinates_file()
        self._process_variants_file()
        self._process_sequences()
        if not(self.rec_path is None):
            self._process_experimental_recombination()
        
    
    def _exponential_smoothing(self,x):
        '''
        Usual exponential smoothing function.
        
        :param x: Array of values to be smoothed.
        :type x: np.ndarray
        
        :return: smoothed values.
        :rtype: np.ndarray
        '''
        # if alpha == 1 s=x, i.e. no smoothing
        assert 0 < self.alpha <= 1
        s = [x[0]]
        for t in range(1,len(x)):
            s.append(self.alpha*x[t] + (1-self.alpha)*s[t-1])
        return s
    
    def _my_exponential_smoothing(self,x):
        '''
        Modified exponential smoothing function, starting from zero.
        
        :param x: Array of values to be smoothed.
        :type x: np.ndarray
        
        :return: smoothed values.
        :rtype: np.ndarray
        '''
        # if alpha == 1 s=x, i.e. no smoothing
        assert 0 < self.alpha <= 1
        s = [0]      # modification to exponential smoothing
        for t in range(1,len(x)):
            s.append(self.alpha*x[t] + (1-self.alpha)*s[t-1])
        return s
    
    def _identify_bases(self,start,end,duplicates=False):
        '''
        Creates a list with the positions of the bases belonging to the 
        intervals determined by start and end lists.
        
        :param start: list where each element corresponds to a base pair 
            number that starts an interval.
        :type start: list of ints.
        :param end: list where each element corresponds to a base pair 
            number that ends an interval.
        :type end: list of ints.
        :param duplicates: If False, bases belonging to overlapping intervals 
            are counted only once. Defaults to False.
        :type duplicates: bool.
        
        :return: list of basepairs belonging to the input intervals.
        :rtype: list of ints
        '''
        assert len(start) == len(end)
        ans = []
        for i in range(len(start)):
            ans += [j for j in range(min(start[i],end[i]),max(start[i],end[i]))]
        if duplicates==False:
            ans = list(set(ans))
        return ans
    
    def _window_features(self):
        '''
        This function calculates the relative frequency of the following 
        features for each window in which the reference chromosome has been divided:
        absent bases (bases that do not belong to any contig of the coordinates file),
        identical bases (bases that belong to the contigs and are not inversions or variants),
        variant bases (bases presenting a variant of type snps and indel in the variants file),
        inversions (bases belonging to inverted contigs of the coordinates file).
        '''        
      
        # mapped bases------------------------------------------------------
        mb = self._identify_bases(self.coords.Sr,self.coords.Er)
        mb_af, _ = np.histogram(mb,bins=self.wbp_r)     #absolute frequency
        mb_rf = mb_af/np.array(self.size_w)             #relative frequency
        
        # absent (no-mapped) bases------------------------------------------
        ab = list(set(range(len(self.ref_seq))).difference(mb))
        ab_af, _ = np.histogram(ab,bins=self.wbp_r)     #absolute frequency
        ab_rf = ab_af/np.array(self.size_w)             #relative frequency
        
        # inversions--------------------------------------------------------
        inv = self._identify_bases(self.coords[self.coords.inversion_q==1].Sr.tolist(),
                                  self.coords[self.coords.inversion_q==1].Er.tolist())
        inv_af, _ = np.histogram(inv,bins=self.wbp_r)
        inv_rf = inv_af/np.array(self.size_w)
        
        # variant bases-----------------------------------------------------
        # consider only snps and indels. 
        # drop duplicates due to contigs overlap
        vb_af, _ = np.histogram(list(set(self.variants[self.variants.rbase!='.'].rpos).difference(inv)),
                            bins=self.wbp_r)
        vb_rf = vb_af/np.array(self.size_w)
           
        #--------identity -------------------------------
        idt = np.subtract(mb_af,vb_af)
        idt = np.subtract(idt,inv_af)
        idt_rf = idt/np.array(self.size_w)
        
        features = pd.DataFrame(list(zip(ab_rf, idt_rf, vb_rf, inv_rf)),
                          columns=['ab_rf','idt_rf', 'vb_rf', 'inv_rf'])
        
        self.features = features

        
    def _model_cases(self,p=None):
        '''
        Model step 1: Modifications to the identity values ​​through 3 cases 
        that consider the identity values, variants and absent bases.
        
        :param p: parameters of step 1 of the model [p1,p2,p3,p4,p5,p6,p7]. 
            A value of p1 is subtracted from windows with identity ​​less than p2 and variants greater than p7.
            A value of p3 is added to windows with identity ​​less than p4 and variants less than p7.
            A value of p5 is subtracted from windows with absent bases ​​less than p6 and variants less than p7.
        :type p: list
        
        :return: Modified identity due to the 3 cases, and list of the applied
            case in each window.
        :rtype: list, list
        '''
        if p is None:
            p = self.params
            
        ans = []
        cases = []
        for i in range(len(self.features.idt_rf)):
            if (self.features.idt_rf[i] < p[1]) and (self.features.vb_rf[i] > p[6]):
                ans.append(max(0,self.features.idt_rf[i] - p[0]))
                cases.append(1)
            elif (self.features.idt_rf[i] < p[3]) and (self.features.vb_rf[i] < p[6]):
                ans.append(max(0,self.features.idt_rf[i] + p[2]))
                cases.append(2)
            elif (self.features.ab_rf[i] > p[5]) and (self.features.vb_rf[i] < p[6]):
                ans.append(max(0,self.features.idt_rf[i]- p[4]))
                cases.append(3)
            else:
                ans.append(max(0,self.features.idt_rf[i]))
                cases.append(0) 
                
        return ans, cases
    
    
    def _centromer_correction(self):
        '''
        Calculation of a weight function that approximates the location of 
        the chromosome centromere. The function is zero between the reference 
        and query windows with the largest number of CentO alignments, forward
        and backward of this region the function performs a linear transition
        to the value 1 in a range of size_wc windows.
        
        :return: weight function to correct for centromere recombination.
        :rtype: list
        '''
        ref_cent = cf.CentOFinder(self.CentO_path, self.ref_path, self.size_w)
        ref_cent.detect_centromere(verbose=False)
        
        qry_cent = cf.CentOFinder(self.CentO_path, self.qry_path, self.size_w)
        qry_cent.detect_centromere(verbose=False)
        
        c0 = min(ref_cent.c_window_number, qry_cent.c_window_number)
        c1 = max(ref_cent.c_window_number, qry_cent.c_window_number)
        
        n = len(self.wbp_r)-1
        
        if c1 >= n//4:
            f = []
            for w in range(n):
                if (c0-self.size_wc) <= w < c0:
                    f.append((-1/self.size_wc)*(w-c0))
                elif c0 <= w < c1:
                    f.append(0)
                elif c1 <= w < c1+self.size_wc:
                    f.append((1/self.size_wc)*(w-c1))
                else: #(0 <= w < c0-size_wc or (c1+size_wc)<=w<=n)
                    f.append(1)
    
        else:
            f = [0 if 0<=i<c1 else 1 for i in range(n)]  # g(w) in the paper
        
        return f
        
        
    def predict_recombination(self,p=None):
        '''
        Perfom the prediction of the crossover recombination through a 4-step model,
        using information from the alignment between the reference and query genomes.
        Calculates the windowed features and then uses them to apply the model
        and make the prediction.
        
        :param p: parameters of step 1 of the model [p1,p2,p3,p4,p5,p6,p7]. 
            A value of p1 is subtracted from windows with identity ​​less than p2 and variants greater than p7.
            A value of p3 is added to windows with identity ​​less than p4 and variants less than p7.
            A value of p5 is subtracted from windows with absent bases ​​less than p6 and variants less than p7.
        :type p: list 
        '''
        
        if p is None:
            p = self.params
        
        if self.features is None:
            print('Computing window features...')
            self._window_features()
            
            if self.params is None:
                print('Optimizing parameters...')
            else:
                print('Computing prediction...')
            
        prediction, self.cases = self._model_cases(p)
        prediction = [i*c for i,c in zip(prediction,self._centromer_correction())]
        self.prediction = self._my_exponential_smoothing(prediction)
        
        
    def _tuning_model(self,test_params,metric='r2'):
        '''
        Evaluates recombination prediction with Pearson's correlation coefficient
        or R2 for the parameter optimization process.
        
        :return: 1 minus the evaluation metric. Value to be minimized in 
            parameter optimization.
        :rtype: float
        '''
        self.predict_recombination(p=test_params)
        
        mask = ~np.isnan(self.experimental_rec)
        if metric == 'pearson':
            ans, _ = pearsonr(np.array(self.prediction)[mask],
                              np.array(self.experimental_rec)[mask])
        elif metric == 'r2':
            ans = r2_score(np.array(self.prediction)[mask],
                           np.array(self.experimental_rec)[mask])
        
        else:
            print('Ivalid value for metric parameter')
        
        return 1-ans   
        
    
    def optimize_model_parameters(self,x0,metric='r2'):
        '''
        Optimization of the 7 parameters of step 1 of the model (p1,p2,p3,p4,p5,p6,p7).
        
        :param x0: Initial values ​​of the parameters to start the search for the optimal ones.
        :type x0: list of floats
        :param metric: Evaluation metric of the model prediction for the 
            optimization process. Choosing between 'pearson' and 'r2', Defoult to 'r2'.
        :type metric: str
        '''
        
        fun = lambda x: self._tuning_model(x,metric)
        bnds = ((0, 1), (0, 1), (0, 1), (0,1), (0,1), (0,1), (0,1))
        cons = ({'type': 'ineq', 'fun': lambda x:  x[1] - x[6]},
                {'type': 'ineq', 'fun': lambda x:  x[3] - x[6]},)
        res = minimize(fun, x0, method='SLSQP', bounds=bnds, constraints=cons)
        self.params = list(res.x)
        
        print('Optimum {0}: {1}'.format(metric,1-res.fun))
        print('Optimum parameters (p1,p2,p3,p4,p5,p6,p7): {0}'.format(self.params))
        
        
    def prediction_evaluation(self,approach='model'):
        '''
        Evaluates recombination prediction with Pearson's correlation coefficient
        and R2. It can compare the experimental recombination with the 
        identity or with the final prediction of the model.
        
        :param approach: Recombination prediction approach to compare against
            experimental recombination. Choosing between 'model' and 'identity'.
            Default set to 'model'.
            
        :return: pearson's correlation coefficient, and R2. 
        :rtype: float, float
        '''
        
        if self.experimental_rec is None:
            print('No experimental recombination data, unable to assess prediction.')
        
        else:
            if approach == 'model' or approach == 'identity':
                if approach == 'model':
                    value = self.prediction
                elif approach == 'identity':
                    idt = np.array(self.features.idt_rf.tolist())
                    value = list(self._my_exponential_smoothing(idt))
            
                
                mask = ~np.isnan(self.experimental_rec)
                corr, _ = pearsonr(np.array(value)[mask],
                                   np.array(self.experimental_rec)[mask])
                r2 = r2_score(np.array(value)[mask],
                                   np.array(self.experimental_rec)[mask])
                return corr,r2
            
            else:
                print('Invalid value for approach variable. Choose between "model" or "identity"')
                
    
    
    def plot_landscape(self,ax=None,fontsize_legend=10,fontsize_axtitle=10):
        '''
        Landscape of identity, model prediction and experimental 
        recombination (if available). The colored bars at the bottom of the
        landscapes indicate which case from the first step of the model is 
        applied in each window.
        '''
        if ax is None:
            fig,ax = plt.subplots(figsize=(15,4))
        
        window = np.array(range(len(self.prediction)))
        idt = np.array(self.features.idt_rf.tolist())
        idt_s = self._my_exponential_smoothing(idt)
        pred = np.array(self.prediction)      # already smoothed
        cases = np.array(self.cases)
        if self.experimental_rec is None:
            rec = np.array([np.nan]*len(pred))
        else:
            rec = np.array(self.experimental_rec) # already smoothed
        
        df = pd.DataFrame(list(zip(window,idt_s,pred,rec,cases)),
                          columns=['window','idt','pred_value','exp_rec','cases'])
    
        # ----plot landscapes-------------------------------------------------------------
        if not(np.isnan(rec).all()):
            sns.lineplot(data=df,x='window',y='exp_rec',color='tab:orange',label='Experimental',ax=ax)
        sns.lineplot(data=df,x='window',y='pred_value',color='tab:blue',label='Prediction',ci='sd',ax=ax)
        sns.lineplot(data=df,x='window',y='idt',color='darkgrey',label='Identity',ax=ax)
        
        
        # ----plot cases-------------------------------------------------------------
        chr_cases=df.cases.tolist()
        for i in  range(0, df.window.nunique()):
            if chr_cases[i] == 1:
                ax.axvspan(i, i+1, facecolor='g', alpha=0.5, ymax=0.1)
            elif chr_cases[i] == 2:
                ax.axvspan(i, i+1, facecolor='b', alpha=0.5, ymax=0.1)
            elif chr_cases[i] == 3:
                ax.axvspan(i, i+1, facecolor='r', alpha=0.5, ymax=0.1)
            else:
                ax.axvspan(i, i+1, facecolor='darkgrey', alpha=0.5, ymax=0.1)
                
        #---plot customization----------------------------------------------------------
        legend1 = ax.legend(loc ="upper left",prop={'size': fontsize_legend})
        ax.add_artist(legend1)
        ax.set_xlabel('Chromosome length (X {:,} bp)'.format(self.size_w), fontsize=fontsize_axtitle)  
        ax.set_ylabel('cM/100Kbp', fontsize=fontsize_axtitle)
    
        patch0 = mpatches.Patch(color='darkgrey', alpha=0.5, label='No case')
        patch1 = mpatches.Patch(color='g', alpha=0.5, label='case 1')
        patch2 = mpatches.Patch(color='b', alpha=0.5, label='case 2')
        patch3 = mpatches.Patch(color='r', alpha=0.5, label='case 3')
        legend2 = ax.legend(handles=[patch0,patch1,patch2,patch3],loc='lower right',
                            prop={'size': fontsize_legend})
        ax.add_artist(legend2)
    
    
    def plot_correlation(self,ax=None,fontsize_legend=10,fontsize_axtitle=10):
        '''
        Plots linear relationship between the experimental recombination and
        the prediction of the model. The marker color in the scatter plot
        indicate which case from the first step of the model is 
        applied in each window.
        '''
        if self.experimental_rec is None:
            print('No experimental recombination data, unable to generate correlation plot')
        
        else:
            ax = ax or plt.gca()
            
            colormap = np.array(['darkgrey','g','b','r'])
            scatter = ax.scatter(self.experimental_rec, self.prediction,
                                 s=50, c=colormap[self.cases],alpha=0.5)
        
            patch0 = mpatches.Patch(color='darkgrey', alpha=0.5, label='No case')
            patch1 = mpatches.Patch(color='g', alpha=0.5, label='case 1')
            patch2 = mpatches.Patch(color='b', alpha=0.5, label='case 2')
            patch3 = mpatches.Patch(color='r', alpha=0.5, label='case 3')
            ax.legend(handles=[patch0,patch1,patch2,patch3],loc='lower right',
                      prop={'size': fontsize_legend})
        
            ax.set_xlabel('recombination', fontsize=fontsize_axtitle)  
            ax.set_ylabel('model prediction', fontsize=fontsize_axtitle)
        
    