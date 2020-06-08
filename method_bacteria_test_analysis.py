import os
import math
import pandas as pd
import regex as re
import numpy as np
import matplotlib.pyplot as plt
import collections
from Bio import SeqIO
from Bio.SeqUtils import GC
from scipy import stats

DIR = r''  # main files directory

genome_dict = dict()  # dict of 30 samples
becterium = ['B.anthracis Genome.fasta',
             'B.fragilis Genome.fasta',
             'B.subtilis Genome.fasta',
             'C.difficile Genome.fasta',
             'E.coli Genome.fasta',
             'M.pneumoniae Genome.fasta',
             'M.smegmatis Genome.fasta',
             'S.aureus Genome.fasta',
             'S.avermitilis Genome.fasta',
             'S.enterica Genome.fasta',
             'V.cholerae Genome.fasta',
             'A.baumannii Genome.fasta',
             'C.botulinum Genome.fasta',
             'C.tetani Genome.fasta',
             'S.cellulosum Genome.fasta',
             'S.coelicolor Genome.fasta',
             'P.fluorescens Genome.fasta',
             'H.pylori Genome.fasta',
             'A.vinelandii Genome.fasta',
             'Y.pestis Genome.fasta',
             'B.cereus Genome.fasta',
             'P.aureginosa Genome.fasta',
             'S.flexneri Genome.fasta',
             'C.pneumoniae Genome.fasta',
             'S.pneumoniae Genome.fasta',
             'S.agalactiae Genome.fasta',
             'E.faecalis Genome.fasta',
             'B.henselae Genome.fasta',
             'B.abortus Genome.fasta',
             'C.trachomatis Genome.fasta']

for file_name in becterium:                     # extract genomes from files
    for seq_record in SeqIO.parse(os.path.join(DIR,file_name), 'fasta'):
        genome_seq = str(seq_record.seq) + str(seq_record.seq.reverse_complement())
        name = file_name.split(' ')[0]
        genome_dict.update({name:genome_seq})
        break
      
def find_matches(genomes, regex):
    """
    funds matches of regex within a genome.
    
    genomes - dictionary
    regex - string
    """
    ta_dict = dict()
    for name, genome in genome_dict.items():
        ta_hits = re.findall(regex, genome.upper(), overlapped = True)
        ta_dict.update({name:ta_hits})
    return ta_dict # All the hits with TA within the bacterial genomes (number of TA)

def count_matches(ta_hits):
    """
    Counts the number of a match appearances in the genome.
    
    ta_hits - dictionary
    """
    df_list = list()   # a frame containing all the counted 14 lenght sequences
    for name, hits in ta_hits.items():
        ta_counted_14 = collections.Counter(hits)
        seq = []
        counts = []
        for k,v in ta_counted_14.items():
            seq.append(k)
            counts.append(v)
        df = {name + ' seq': seq, name + ' counted hits' : counts}
        df_list.append(pd.DataFrame(df))
    return df_list #  df of data frames containing all the counted 14 lenght sequences

# Previous methods of a single side reading
regex_14 = 'TA.{14}'
regex_14_matches = find_matches(genome_dict, regex_14)
counted_14_hits = count_matches(regex_14_matches)  # for each bacteria, counts the identical accuring seq
 
# Our proposed method
regex_26 = '.{12}TA.{12}'
regex_26_matches = find_matches(genome_dict, regex_26)
counted_26_hits = count_matches(regex_26_matches)

def get_ambiguous_hits(data):   # removes percise hits from data
    """
    Extracts the ambiguous and seperates from the 
    percise (single) hits of the data.
    
    data - DataFrame
    """
    sum_amb_hits = dict() # per bacterium 
    amb_hits = dict() 
    for df in data:
        name = str(df.columns).split(' ')[0].split("'")[1]
        hits = df.iloc[:,1][df.iloc[:,1] > 1] # ambiguous hits (not those which appears once)
        sum_amb_hits.update({name:sum(hits)}) 
        amb_hits.update({name:hits})
    total_amb_hits = pd.DataFrame.from_dict(sum_amb_hits, orient='index')
    amb_hits_by_group = pd.DataFrame.from_dict(amb_hits)
    return total_amb_hits, amb_hits_by_group # 2 df, total hits and grouped

total_14, by_group_14 = get_ambiguous_hits(counted_14_hits)
total_26, by_group_26 = get_ambiguous_hits(counted_26_hits)

total_14 = total_14.rename(columns={0: 'Current System\n  (TA+14N)'})
total_26 = total_26.rename(columns={0: 'Proposed System\n(12N+TA+12N)'})
both_total_compared = pd.merge(total_14, total_26, left_index=True,
                                                   right_index=True)
both_total_compared['Ratio (N/M)'] = (both_total_compared['Current System\n  (TA+14N)'] /
                               both_total_compared['Proposed System\n(12N+TA+12N)'])

################## statistical analysis ##################

def get_statistics(data1, data2, lst):
    """
    Gets statistics from data, comparison between both methods.
    
    data1 - DataFrame, method 1
    data2 - DataFrame, method 2
    lst - list, the bacterium info
    """
    variances_dict = dict()
    statistics_dict = dict()
    for bacteria in lst:
        bacteria_name = bacteria.split(' ')[0]
        x = data1[bacteria_name].dropna()
        y = data2[bacteria_name].dropna()
        variances_dict.update({bacteria_name : stats.levene(x,y)})
        statistics_dict.update({bacteria_name : stats.ttest_ind(x,y,
                                                equal_var = False)})
    return variances_dict, statistics_dict

v_dict,s_dict = get_statistics(by_group_14, by_group_26, becterium) # variance, statistics

################## further analysis, table of  bacterial dataframe ##################

ta_over_genome = {k:math.log10(len(v)/2) for k, v in regex_14_matches.items()}  # number of TA over bacterial genome
both_total_compared['TA over genome (log10)'] = pd.DataFrame(
        list(ta_over_genome.values()), index = ta_over_genome.keys()) 

organism, genome_len, gc = [],[],[]
for k, v in genome_dict.items():
    organism.append(k)
    genome_len.append(round((len(v)/2)/(10**6),2))
    gc.append(round(GC(v),2))
data = {'Genome lenght (Mbp)': genome_len, 'GC content (%)':gc}
bacteria_data = pd.DataFrame(data, index = organism, )
both_total_compared = pd.concat([both_total_compared, bacteria_data], axis=1)

# Graphical representation of methods, comparison of total ambiguous insertation
boxplot = both_total_compared.boxplot(column=['Current System\n  (TA+14N)',
                                              'Proposed System\n(12N+TA+12N)'],
                                        meanline=True, showmeans=True,
                                        fontsize='large',
                                        figsize=(10,10))                                              
boxplot.set_title('Comparison of Methods',fontsize=18)
boxplot.set_ylabel('Equivocal Insertions',fontsize=14)
fig = boxplot.get_figure()
fig.savefig('results.pdf')                                             

# Linear regression for TA counts and hits ratio
X = np.array(both_total_compared['TA over genome (log10)'])
Y = np.array(both_total_compared['Ratio (N/M)'])
slope, intercept, r_value, p_value, std_err = stats.linregress(X, Y) # reference
cm = both_total_compared.corr() # correlation matrix to the data

plt.style.use('ggplot')
X_mean = np.mean(X)
Y_mean = np.mean(Y)
num = 0
den = 0
for i in range(len(X)):
    num += (X[i] - X_mean)*(Y[i] - Y_mean)
    den += (X[i] - X_mean)**2
m = num / den
c = Y_mean - m*X_mean
Y_pred = m*X + c #  linear equation

plt.scatter(X, Y)
plt.plot([min(X), max(X)], [min(Y_pred), max(Y_pred)], color='black')
plt.ylabel('Improvement Ratio', fontsize = 'small') 
plt.xlabel('Number of “TA” Within The Genome (log10)', fontsize = 'small')
plt.savefig('graph.pdf')
