* Written by Camila Riccio, Mauricio Peñuela, Camilo Rocha and Jorge Finke
* Last update: 04/04/22 

Rice recombination predictor
============================

This is a Python3 implementation to predict local chromosomal recombination in rice
using sequence identity, combined with other features derived from
genome alignment (including the number of variants, inversions, absent bases, and
CentO sequences).

Preliminaries
-------------

In order for the user to predict recombination between two parental rice varieties,
arbitrarily select one of them as the reference genome and the other as the query genome.
As an example we will take the IR64 variety as reference genome and Azucena as query, and
we will use the data from chromosome 01 of both varieties to calibrate the prediction model.

We suggest organizing the data as follows:

* **input_data** folder containing:

  * **IR64** folder containig:

    * 12 fasta files with the amino acid sequences of each chromosome of the IR64 variety. 
      Download from `NCBI Genome database <https://www.ncbi.nlm.nih.gov/genome>`_
      searching the accession number RWKJ00000000.

  * **Azucena** folder containig:

    * 12 fasta files with the amino acid sequences of each chromosome of the Azucena variety.
      Download from `NCBI Genome database <https://www.ncbi.nlm.nih.gov/genome>`_
      searching the accession number PKQC000000000.

  * **coords** folder (Empty, the corresponding files will be generated later). 

  * **snps** folder (Empty, the corresponding files will be generated later).

  * **recombination** folder containig:

    * 12 csv files with windowed experimental recombination values for each chromosome.

  * **CentO** fasta file.

* outpud_data folder (Empty, the corresponding files will be generated later).

The alignment process between the two parental chromosomes must be done independently
using the `MUMmer <http://mummer.sourceforge.net/manual/>`_
software, with the commands shown below executed from the **input_data** folder:

Align reference and query fasta files::
  
  nucmer --prefix=IR64_Azucena_chr01 IR64/Osat_IR64_chr01.fasta Azucena/Osat_Azucena_chr01.fasta

Filter the aligned data::
  
  delta-filter -r -q IR64_Azucena_chr01.delta > IR64_Azucena_chr01.filter

Extract contig coordinates from the filtered file::
  
  show-coords -r IR64_Azucena_chr01.filter > coords/IR64_Azucena_chr01.coords

Extract variants from the filtered file::
  
  show-snps -lr -x 1  -T IR64_Azucena_chr01.filter  > snps/IR64_Azucena_chr01.snps

We are now ready to make use of the Rice recombination predictor software.

Setup
------
Clone the repository::

  git clone git@github.com/criccio35/Rice-recombination-predictor


Requirements
------------
Install the requirements by entering the following commands in the terminal:

Install biopython module::

  pip install biopython

Install Basic Local Alignment Search Tool (BLAST)::

  sudo apt update
  sudo apt install ncbi-blast+


How to use
----------

For optimal model performance, you also need experimental recombination data from at least one chromosome.
Depending on the availability of experimental data, 3 cases can be presented:

- No experimental data at all.
- Experimental data for one chromosome only.
- Experimental data for all chromosomes.

We will show the example of the case in which the experimental data is available for all 
chromosomes and and along the way we will explain what modifications would be necessary 
for the other two cases.

Import module::
  
  import rice_recombination_predictor as rrp

Specify the paths of the input files of the chromosome with which the model is to be calibrated::
  
  ref_path = 'input_data/IR64/Osat_IR64_chr01.fasta'
  qry_path = 'input_data/Azucena/Osat_Azucena_chr01.fasta'
  coords_path = 'input_data/coords/IR64_Azucena_chr01.coords'
  variants_path = 'input_data/snps/IR64_Azucena_chr01.snps'
  rec_path = 'input_data/recombination/experimental_recombination_chr01.csv'
  CentO_path = 'input_data/CentO_AA.fasta'
  results_path = 'output_data/'

Input window size::

  size_w = 100_000

Instantiate the rice_recombination_predictor class::

  chr_mod = rrp.rice_recombination_predictor(size_w, ref_path, qry_path, 
                                          coords_path, variants_path, 
                                          CentO_path, rec_path)

Call the method to preprocess the input files::
  
  chr_mod.data_preprocessing()

Define the initial parameters p1,p2,p3,p4,p5,p5 and p7 (see **Software description** file
for model details) and call the method to optimize them::
  
  initial_parameters = [0.5,0.97,1,0.9,1,0,0.002]
  chr_mod.optimize_model_parameters(initial_parameters)

You can access the predicted values with the following attribute of the class::
  
  chr_mod.prediction

You can visualize the results of the prediction through the method::
  
  chr_mod.plot_landscape()

.. image:: output_data/prediction_chr01.png
  :width: 800

You can visualize the performance of the prediction with respect to the experimental data through the method::
  
  chr_mod.plot_correlation()

.. image:: output_data/performance_chr01.png
  :width: 400

You can now use the optimized parameters to predict recombination on the other chromosomes.
If **experimental recombination data is NOT available for the remaining chromosomes** 
you can proceed as follows::

  import matplotlib.pyplot as plt
  import seaborn as sns
  import pandas as pd

  chromosomes = ['01','02','03','04','05','06','07','08','09','10','11','12']

  df_pred = pd.DataFrame()
  df_eval = pd.DataFrame(columns=['chromosome','R2_idt','R2_pred','pearson_idt','pearson_pred'])

  for chr_nbr in chromosomes:
      ref_path = 'input_data/IR64/Osat_IR64_chr{0}.fasta'.format(chr_nbr) 
      qry_path = 'input_data/Azucena/Osat_Azucena_chr{0}.fasta'.format(chr_nbr) 
      coords_path = 'input_data/coords/IR64_Azucena_chr{0}.coords'.format(chr_nbr)
      variants_path = 'input_data/snps/IR64_Azucena_chr{0}.snps'.format(chr_nbr)
      
      chr_tmp = rrp.rice_recombination_predictor(size_w, ref_path, qry_path, 
                                                  coords_path, variants_path, 
                                                  CentO_path, params=chr_mod.params)
      
      print('---Chromosome {0}---'.format(chr_nbr))
      chr_tmp.data_preprocessing()
      chr_tmp.predict_recombination()
      
      df_pred1 = pd.DataFrame({'chr_'+chr_nbr : chr_tmp.prediction})
      df_pred = pd.concat([df_pred,df_pred1], axis=1)
      
      chr_tmp.plot_landscape()
      
  df_pred.to_csv(results_path+'predictions.csv',header=True,index=True)


The above will generate a .csv file where each column corresponds to the prediction of recombination on a chromosome.
At the same time, for each chromosome, a plot will be generated as follows:

.. image:: output_data/prediction_no_exprec_chr12.png
  :width: 800

If **experimental recombination data are available for all chromosomes**
you can proceed as follows::

  import matplotlib.pyplot as plt
  import seaborn as sns
  import pandas as pd

  model_chr_nbr = '01' # chromosome used in parameter optimization
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

The above will generate a .csv file where each column corresponds to the prediction of recombination on a chromosome.
At the same time, for each chromosome, a plot will be generated as follows:

.. image:: output_data/prediction_and_performance_chr02.png
  :width: 800

With the performance data stored in the *df_eval* dataframe, the plots shown below can be generated.
These plots compare the performance of the two approaches to predict recombination: the identity and 
the model (modified identity).

Plot for coefficient of determination R2::
  
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

.. image:: output_data/R2_performance.png
  :width: 600

Plot for pearson correlation coefficient r::
  
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

.. image:: output_data/r_performance.png
  :width: 600

If there is **no experimental data for any chromosome** you can use your own combination of parameters
to make the prediction::
  
  parameters = [0.5,0.97,1,0.9,1,0,0.002]
  chr_test = rrp.rice_recombination_predictor(size_w, ref_path, qry_path, 
                                            coords_path, variants_path, 
                                            CentO_path, params=parameters)
  chr_test.data_preprocessing()
  chr_test.prediction
  chr_test.plot_landscape()


The complete example, for the case where experimental recombination data is available 
for all chromosomes, is in the file **test.py**