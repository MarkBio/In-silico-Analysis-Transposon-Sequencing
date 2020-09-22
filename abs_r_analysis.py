import os
import pandas as pd
import regex as re
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Simple functions that are used in this analysis

# First step of analysis, extracing files from the resulted fasta/fastq files.
def get_sequences(file_name, _format):
    """
    Parsing a file to retrieve its sequences.
    
    file_name - file directory
    _format - fasta/fastq
    """
    all_seq = list()
    for seq_record in SeqIO.parse(os.path.join(DIR,file_name), _format):
        seq = seq_record.seq
        all_seq.append(str(seq))
    return list(set(all_seq)) # stored as a list of unique seq, removing duplicates.

# After extracting the sequences, reverse compliment is required for true hit detection.
def get_sequences_reverse(file_name, _format):
    """
    Parsing a file to retrieve the reverse complement of sequences.
    
    file_name - file directory
    _format - fasta/fastq
    """
    all_seq = list()
    for seq_record in SeqIO.parse(os.path.join(DIR,file_name), _format):
        seq = seq_record.seq
        r_seq = seq.reverse_complement()  # reverse complement of the sequence
        all_seq.append(str(r_seq))
    return list(set(all_seq))    # stored as a list of unique seq, removing reverse duplicates

# Basicly searching for the required pattern, using regex over the data
# In our case: IR + TA + 24N + TA + IR , asa whole
def find_unique(sequences, pattern):
    """
    Finding matches of a regex pattern within each sequence.
    
    sequences - a list of sequences
    pattern - regular expression pattern (as acceptable)
    """
    possible_hits = list()
    for read in sequences:
        matches = re.findall(pattern, read, overlapped = False)
        if len(matches) == 0:
            pass
        else:
            possible_hits.append(matches[0])
    return list(set(possible_hits)) # list of unique matches of the regex

# For the string (seq) manipulation. Cutting of the desired sequence for rearrangement 
def extract_28(unique):
    """
    Based on matches, snipsa a 28 bp sequence, removes invariables.
    
    unique - a list of sequences
    """
    seq_28 = list()
    for hit in unique:
        target_seq = hit[27:55]
        seq_28.append(target_seq)
    return seq_28 # list of short sequences

# The process of rearrangement of the sequence to its original order
def rearrangement(sequences_28):
    """
    Based on short sequences, basic rearrangement of sequences to 
    thier original sequence as was in the genome.
    
    sequences_28 - a list of short sequences
    """
    to_seek = list()
    for i in range(len(sequences_28)):
        before = sequences_28[i]
        half1 = before[2:14]
        half2 = before[14:26]
        after = Seq(half2 + 'TA' + half1)
        record = SeqRecord(after,str(i))
        to_seek.append(record)
    return to_seek  # list of queries 12N+TA+12N, 'TA' is the insertion point, 'N' are genomic 

# Additional file parser
def get_reads(file1, file2, _format):
    """
    Parsing a file for its sequences.
    
    file1 - file directory
    file2 - file directory
    _format - fasta/fastq
    """
    R1 = list()
    for s1 in SeqIO.parse(os.path.join(DIR,file1), _format):
        seq1 = str(s1.seq)
        R1.append(seq1)
    R2 = list()
    for s2 in SeqIO.parse(os.path.join(DIR, file2), _format):
        seq2 = str(s2.seq)
        R2.append(seq2)
    return R1, R2 # as lists of sequences extracted from files (fasta/fastq)

# Final step of mapping to a reference genome, finding hits 
def init_mapper(reads, strand): 
    """
    Basic mapping function, based on regex, finds a match in a genome 
    and retrieves the location of insertion within the genome.
    
    reads - a list of queries 
    strand - the genome sequence as a string (posstive of negative)
    """
    found = list()
    for read in reads:
        for match in re.finditer(read, strand):
            info = '%02d-%02d: %s' % (match.start(), match.end(), match.group(0))
            found.append((read,info))
    return found # the hits found by mapping

################## sequences handling step ##################
    
DIR = ''            #  main files directory 
file_name_1 = ''    
file_name_2 = ''    #  fasta/fastq directory
FORMAT = r'AGACCGGGGACTTATCATCCAACCTGT.{28}ACAGGTTGGATGATAAGTCCCCGGTCT'  # in our case IR + 28N (from genomic seq) + IR
# Or for a more flexible search (based on read quility), format could be: r'TA.{24}TA

# "r" in variable name means reverse
# First analysis through parsing, based on function above
all_sequeces = get_sequences(file_name_1, 'fastq') # Start with parsing and sequence retrieving R1 and R2
all_sequeces_r = get_sequences_reverse(file_name_2, 'fastq')    # all file seq
unique_sequences = find_unique(all_sequeces, FORMAT)
unique_sequences_r = find_unique(all_sequeces_r, FORMAT)        # unique seq
seq_28 = extract_28(unique_sequences)
seq_28_r = extract_28(unique_sequences_r)                       # short not-ordered seq
reads_28 = rearrangement(seq_28)
reads_28_r = rearrangement(seq_28_r)                            # short ordered seq

with open("abs_R1.fasta", "w") as output_handle: # reading newly constructed files
    SeqIO.write(reads_28, output_handle, "fasta")
with open("abs_R2.fasta", "w") as output_handle:
    SeqIO.write(reads_28_r, output_handle, "fasta")
R1_reads, R2_reads = get_reads('abs_R1.fasta', 'abs_R2.fasta', 'fasta') # different DS, as lists from FASTA

# The genome of M.abs (Mycobrowser)
for record in SeqIO.parse(os.path.join(DIR,'Mycobacterium_abscessus_ATCC_19977_genome_v3.fasta'), 'fasta'):
    genome_pos = str(record.seq)
    genome_neg = str(record.seq.reverse_complement())
    break # the raw genome as string

# Mapping to the genome of M.abs, tring to map on both strands, (R1 and R2)
r1_to_pos = init_mapper(R1_reads, genome_pos)   # run on both strands of genome
r1_to_neg = init_mapper(R1_reads, genome_neg)
r2_to_pos = init_mapper(R2_reads, genome_pos)
r2_to_neg = init_mapper(R2_reads, genome_neg)

gene_dict = dict()  # Retrieving The annotated genes of M.abs (Mycobrowser)
for gene_record in SeqIO.parse(os.path.join(DIR,'Mycobacterium_abscessus_ATCC_19977_genes_v3.fasta'), "fasta"):
    key = str(gene_record.id)
    value = str(gene_record.seq).upper()
    gene_dict.update({key:value})

quary_reads = r1_to_neg + r1_to_pos + r2_to_neg + r2_to_pos # to pool all sequences, matches
match_list = list()   
for quary in quary_reads:
    quary_insert = quary[0]
    quary_insert_r = str(Seq(quary[0]).reverse_complement())
    for gene_id, gene in gene_dict.items():
        if quary_insert in gene:
            match_list.append((quary_insert,gene_id))
        if quary_insert_r in gene:
            match_list.append((quary_insert_r,gene_id))
gene_hits = list(set(match_list))  # initial results of mapped hits, based on genes of M.abs (seen above)

################## Attempt-of-annotation step ##################

insert_seq = []         # the insertet sequences
mabs_x = []             # name of locus
gene_coor = []          # gene coordinates as numbers
strand_type = []        # orientation (+/-)
gene_seq = []           # actual gene sequence
pinpoint = []           # insertion coordinates
g_length = []           # actual gene lenght in bp

# Secondary results of parsing and attempts to annotate the genes
for duo in gene_hits: # parsing to create a table of genes and more information
    insert_seq.append(duo[0])
    temp = duo[1].split('|')
    mabs_x.append(temp[0])
    gene_coor.append(temp[3])
    strand_type.append(temp[4])
for locus in mabs_x:  # find the full sequence of the harmed gene (insertion)
    for k,v in gene_dict.items():
        if locus in k:
            gene_seq.append(v)
for i in range(len(insert_seq)):
    pin = [(m.start(0), m.end(0)) for m in re.finditer(insert_seq[i], gene_seq[i])]
    pinpoint.append(pin[0][0]+12)
    g_length.append(len(gene_seq[i]))
# The Dataframe that contain all the gene information, library result. 
mabs_dict = {'Hit Sequence': insert_seq, 'Gene Sequence':gene_seq,
             'Gene length': g_length, 'Insertion Point': pinpoint,
             'Locus': mabs_x, 'Gene Coordinates': gene_coor,
             'Strand': strand_type}
mabs_df = pd.DataFrame(mabs_dict)    # DataFrame of the results

################## non coding sequences ##################

# Same as above, but with no reference to gene, just the position in the genome
not_genes = list()
matches_as_list = list(set([elem[0] for elem in match_list]))
quary_as_list = list(set([q[0] for q in quary_reads]))

for q in quary_as_list:
    if q not in matches_as_list:
        not_genes.append(q)
test = list()
for ng in not_genes:
    ng_r = str(Seq(ng).reverse_complement())
    for gene_id2, gene2 in gene_dict.items():
        if ng in gene2 or ng_r in gene2:
            test.append((ng))
test2 = list()           
for sq in not_genes:
    if sq not in test:
        test2.append(sq)
        test2.append(str(Seq(sq).reverse_complement()))           
ncr = init_mapper(test2, genome_pos)
ncr_r = init_mapper(test2, genome_neg)  #  results of non-coding seq that were hit by transposition
