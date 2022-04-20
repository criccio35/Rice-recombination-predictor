#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  18 12:54:34 2022

@author: camila Riccio and Mauricio Peñuela
"""

# Libraries
import numpy as np
import matplotlib.pyplot as plt

from Bio.Blast.Applications import NcbiblastnCommandline
from io import StringIO
from Bio.Blast import NCBIXML

class CentOFinder():
    '''
    A class used to detect the location of the centromere on rice chromosomes,
    based on the frequency of CentO sequences.
    
    :param CentO_path: CentO fasta file path. 
    :type CentO_path: str
    :param chromosome_path: chromosome fasta file path.
    :type chromosome_path: str
    :param size_w: Observation window size (base pairs)
    :type size_w: int
    :param wbp: chromosome bins/windows in base pairs.
    :type wbp: list of ints
    :param total_windows: Total windows number
    :type total_windows: int
    :param chromosome_length: total number of base pairs on the chromosome
    :type chromosome_length: int
    :param CentO_freq: CentO base pair frequency per chromosome.
    :type CentO_freq: list of ints.
    :param c_window_number: chromosome window number with the highest frequency of CentO alignments.
    :type c_wondow_number: int
    :param c_window_interval: chromosome window interval, in base pairs, with the highest frequency of CentO alignments.
    :type c_window_interval: tuple of ints.
    '''
    
    # Default constructor
    def __init__(self, CentO_path, chromosome_path, size_w,
                 wbp=None,total_windows=None,chromosome_length=None,
                 CentO_freq=None, c_window_number=None,c_window_interval=None):

        self.CentO_path = CentO_path         
        self.chromosome_path = chromosome_path
        self.size_w = size_w
        self.wbp = wbp
        self.total_windows = total_windows
        self.chromosome_length = chromosome_length 
        self.CentO_freq = CentO_freq
        self.c_window_number = c_window_number
        self.c_window_interval = c_window_interval
        
    
    def _read_chr_sequence(self):
        '''
        Returns the whole chromosome sequence of a fasta file in only one line.
        '''
        with open(self.chromosome_path,'r') as f:
            f.readline()  
            seq = ''
            for line in f.readlines():
                seq += line.strip('\n').upper()
        return seq
    
    
    def _CentO_locs(self,start,end,duplicates=False):
        '''
        List the base pair locations that belong to each CentO alignment interval.
        
        :param start: list where each element corresponds to the base number where an alignment between the CentO and the chromosome sequence begins.
        :type start: list of ints.
        :param end: list where each element corresponds to the base number where an alignment between the CentO and the chromosome sequence ends.
        :type end: list of ints.
        :param duplicates: If False, bases belonging to overlapping regions of multiple CentO alignments are counted only once.Defaults to False
        :type duplicates: bool.
        '''
        assert len(start) == len(end)
        ans = []
        for i in range(len(start)):
            ans += [j for j in range(min(start[i],end[i]),max(start[i],end[i]))]
        if duplicates==False:
            ans = list(set(ans))
        return ans
    
    
    def _CentO_frequency(self):
        '''
        Returns the number of base pairs belonging to a CentO alignment
        per window in the chromosome.
        '''
        output = NcbiblastnCommandline(query= self.CentO_path,
                                       subject= self.chromosome_path,
                                       outfmt=5,
                                       gapextend=1,gapopen=1,penalty=-3)()[0]
        
        blast_result_record = NCBIXML.read(StringIO(output))

        start, end = [],[]
        for alignment in blast_result_record.alignments:
            for hsp in alignment.hsps:
                start.append(hsp.sbjct_start-1)
                end.append(hsp.sbjct_end)
        
        chr_sequence = self._read_chr_sequence()
        self.chromosome_length = len(chr_sequence)
        wbp = [i for i in range(0,self.size_w*(self.chromosome_length//self.size_w)+1,self.size_w)]
        CentO_freq_x_window, _ = np.histogram(self._CentO_locs(start,end),bins=wbp)
        
        self.wbp = wbp
        self.total_windows = len(wbp)-1
        self.CentO_freq = CentO_freq_x_window
        
        
    def detect_centromere(self, verbose=True):
        '''
        Computes the frequency of CentO alignments per window, 
        the window number with the highest frequency, 
        and the corresponding 2Mbp and 3Mbp centromeric region prediction.
        
        :param verbose: If True print information about the results, 
            i.e., chromosome length, window size, Total chromosome windows,
            window with highest CentO frequency, window midpoint and 
            centromeric region prediction of 2Mbp and 3Mbp.
            Default to True.
        :type verbose: bool
        '''
        self._CentO_frequency()
        self.c_window_number = np.where(self.CentO_freq == max(self.CentO_freq))[0][0]
        w_start = self.c_window_number*self.size_w
        w_end = (self.c_window_number+1)*self.size_w
        self.c_midpoint = (w_start + w_end)//2
        
        if verbose == True:
        
            print('Chromosome length: {0} base pairs'.format(self.chromosome_length))
            print('Window size: {0}'.format(self.size_w))
            print('Total chromosome whindows: {0}'.format(self.total_windows))
            print('----------------------------------------')
            print('Approximate location of centromere:')
            print('Window number: {0}'.format(self.c_window_number))
            print('Window midpoint (base pair): {0}'.format(self.c_midpoint))
            print('2Mbp centromeric region prediction: [{0},{1}]'.format(self.c_midpoint-1_000_000,
                                                                         self.c_midpoint+1_000_000))
            print('3Mbp centromeric region prediction: [{0},{1}]'.format(self.c_midpoint-1_500_000,
                                                                         self.c_midpoint+1_500_000))
        
    
    def plot_CentO_frequency(self,color='dodgerblue',label='CentO frequency'):
        '''
        Plot the frequency of base pairs belonging to a CentO alignment for
        each window of the chromosome. Also displays the predicted centromeric region
        of 2Mbp and 3Mbp.
        
        :param color: color of the line that represents the frequency of alignments
        :type color: string
        :param label: label for the frequency line of the alignments
        :type label: string
        '''
        fig, ax = plt.subplots(figsize=(15,4))
        ax.axvline(x=np.digitize(self.c_midpoint-1_000_000,self.wbp,right=False)-1,
                   linestyle='--',color='red',label='2Mbp centromeric region')
        ax.axvline(x=np.digitize(self.c_midpoint+1_000_000,self.wbp,right=False)-1,
                   linestyle='--',color='red')
        ax.axvline(x=np.digitize(self.c_midpoint-1_500_000,self.wbp,right=False)-1,
                   linestyle='--',color='green',label='3Mbp centromeric region')
        ax.axvline(x=np.digitize(self.c_midpoint+1_500_000,self.wbp,right=False)-1,
                   linestyle='--',color='green')
        
        ax.plot(self.CentO_freq,label=label,lw=2,color=color)
        
        #---plot customization-----------------------------------------------------------
        ax.set_xlabel('Chromosome length (x{0}bp)'.format(self.size_w), fontsize=20)
        ax.set_ylabel('Count', fontsize=20)
        ax.legend(fontsize=16)
        plt.show()
        