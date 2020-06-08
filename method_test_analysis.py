import os
import random
import math
import pandas as pd
import regex as re
import numpy as np
import matplotlib.pyplot as plt
import collections

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.SeqUtils import GC
from scipy import stats


DIR = r'C:\Users\mark\Desktop\Bacterial Genomes'

genome_dict = dict()
nucleotides = ['A', 'T', 'G', 'C']
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
             'T.pallidum Genome.fasta',
             'C.botulinum Genome.fasta',
             'C.tetani Genome.fasta',
             'S.cellulosum Genome.fasta',
             'S.coelicolor Genome.fasta',
             'P.fluorescens Genome.fasta',
             'P.syringae Genome.fasta',
             'A.vinelandii Genome.fasta',
             'Y.pestis Genome.fasta',
             'B.cereus Genome.fasta',
             'P.aureginosa Genome.fasta',
             'S.flexneri Genome.fasta',
             'S.boydii Genome.fasta',
             'S.pneumoniae Genome.fasta',
             'S.agalactiae Genome.fasta',
             'E.faecalis Genome.fasta',
             'B.henselae Genome.fasta',
             'B.abortus Genome.fasta',
             'C.trachomatis Genome.fasta']

for file_name in becterium:
    for seq_record in SeqIO.parse(os.path.join(DIR,file_name), 'fasta'):
        genome_seq = str(seq_record.seq)
        name = file_name.split(' ')[0]
        genome_dict.update({name:genome_seq})
        break
    
for k, v in genome_dict.items():
    print(k, 'GC%:', round(GC(v), 2), 'Genome lenght:', len(v))

def find_matches(genomes, regex):  # genomes = dictionary, regex = string
    ta_dict = dict()
    for name, genome in genome_dict.items():
        ta_hits = re.findall(regex, genome.upper(), overlapped = True)
        ta_dict.update({name:ta_hits})
    
    return ta_dict # All the hits with TA within the bacterial genomes (number of TA)

def count_matches(ta_hits):  # ta_hits = dictionary
    df_list = list()   # a data frame containing all the counted 14 lenght sequences
    for name, hits in ta_hits.items():
        ta_counted_14 = collections.Counter(hits)
        seq = []
        counts = []
        for k,v in ta_counted_14.items():
            seq.append(k)
            counts.append(v)
        df = {name + ' seq': seq, name + ' counted hits' : counts}
        df_list.append(pd.DataFrame(df))
    return df_list # a list of data frames containing all the counted 14 lenght sequences

# Old method of single side reading
regex_14 = 'TA..............'
regex_14_matches = find_matches(genome_dict, regex_14)
counted_14_hits = count_matches(regex_14_matches)  # for each bacteria, counts the identical accuring seq
 
# New method
regex_26 = '............TA............'
regex_26_matches = find_matches(genome_dict, regex_26)
counted_26_hits = count_matches(regex_26_matches)

def get_ambiguous_hits(data):   # removes percise hits from data
    sum_amb_hits = dict() # per bacterium 
    amb_hits = dict() 
    for df in data:
        name = str(df.columns).split(' ')[0].split("'")[1]
        hits = df.iloc[:,1][df.iloc[:,1] > 1] # ambiguous hits (not those who appear once)
        sum_amb_hits.update({name:sum(hits)}) 
        amb_hits.update({name:hits})
    total_amb_hits = pd.DataFrame.from_dict(sum_amb_hits, orient='index')
    amb_hits_by_group = pd.DataFrame.from_dict(amb_hits)
    return total_amb_hits, amb_hits_by_group

total_14, by_group_14 = get_ambiguous_hits(counted_14_hits)
total_26, by_group_26 = get_ambiguous_hits(counted_26_hits)

total_14 = total_14.rename(columns={0: 'TA+14 Hits (N)'})
total_26 = total_26.rename(columns={0: '12+TA+12 Hits (M)'})
both_total_compared = pd.merge(total_14, total_26, left_index=True, right_index=True)
both_total_compared['Ratio (N/M)'] = (both_total_compared['TA+14 Hits (N)'] /
                               both_total_compared['12+TA+12 Hits (M)'])

# statistical analysis
def get_statistics(data1, data2, lst):
    variances_dict = dict()
    statistics_dict = dict()
    for bacteria in lst:
        bacteria_name = bacteria.split(' ')[0]
        x = data1[bacteria_name].dropna()
        y = data2[bacteria_name].dropna()
        variances_dict.update({bacteria_name : stats.levene(x,y)})
        statistics_dict.update({bacteria_name : stats.ttest_ind(x,y, equal_var = False)})
    return variances_dict, statistics_dict

a,b = get_statistics(by_group_14, by_group_26, becterium)

ta_over_genome = {k:math.log10(len(v)) for k, v in regex_14_matches.items()}
both_total_compared['TA over genome (log10)'] = pd.DataFrame(list(ta_over_genome.values()), index = ta_over_genome.keys()) 

# Linear  regression for TA and hits 
X = np.array(both_total_compared['TA over genome (log10)'])
Y = np.array(both_total_compared['Ratio (N/M)'])

slope, intercept, r_value, p_value, std_err = stats.linregress(X, Y)
m, b = np.polyfit(X, Y, 1)
plt.plot(X, Y, 'o')
plt.plot(X, m*X + b)


organism, genome_len, gc = [],[],[]
for k, v in genome_dict.items():
    organism.append(k)
    genome_len.append(round(len(v)/(10**6),2))
    gc.append(round(GC(v),2))
data = {'Genome lenght (Mbp)': genome_len, 'GC content (%)':gc}
bacteria_data = pd.DataFrame(data, index = organism, )
both_total_compared = pd.concat([both_total_compared, bacteria_data], axis=1)




def mean_over_each(array):
    amb_insertions = collections.Counter(array.dropna())    
    total_amb_insertions = 0
    total_amb_seq = 0
    for repeat, hit in amb_insertions.items():
        total_amb_insertions += repeat*hit
        total_amb_seq += hit
    mean_number_of_options = total_amb_seq / total_amb_insertions  # there is need reveresed
    return mean_number_of_options

def get_mean_over_each(dataframe):
    mean_data = list()
    for bac in list(dataframe):
        mean = mean_over_each(dataframe[bac])
        mean_data.append(mean)
    return mean_data

both_total_compared['Average ambiguous insertions (n)'] = get_mean_over_each(by_group_14)
both_total_compared['Average ambiguous insertions (m)'] = get_mean_over_each(by_group_26)        
both_total_compared['Ratio (n/m)'] = (both_total_compared['Average ambiguous insertions (n)'] /
                               both_total_compared['Average ambiguous insertions (m)'])



X = np.array(both_total_compared['TA over genome (log10)'])
Y = np.array(both_total_compared['Ratio (n/m)'])

slope, intercept, r_value, p_value, std_err = stats.linregress(X, Y)
m, b = np.polyfit(X, Y, 1)
plt.plot(X, Y, 'o')
plt.plot(X, m*X + b)

# אם יש יותר היטס אז כמובן שהממוצע הוא קטן יותר!!!!!!!!!!!!!!!!!!!
# צריך לשנות אתזה כי עכשיו זה ניראה כיאלו בשיטה הישנה יש פחות AMB רצפים

# המספר של הרשיו מראה שבממוצע על כל 100 רצפים של השיטה החדשה יהיה לך 179 רצפים בשיטה הישנה. זה באופן יחסי ולא מוחלט כמו התוצאות האחרות! וזה הכי נכון



'0'