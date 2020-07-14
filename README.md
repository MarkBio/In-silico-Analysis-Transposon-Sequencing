# Bacterial-Genome-Analysis
method_bacteria_test_analysis.py - in-silico analysis of 30 bacterial genomes. Comparisons of our proposed 
method to current Tn-seq analysis methods.

- Parsing of 30 selected bacterial genomes (names included).
- Functions contain simple analysis of the bacterial genomes; referred as string, and string manipulation. TA dinucleotide 
  were counted and marked. Theoretical hits assessments were performed by regular expression (regex) including overlapping 
  sequences. Matches comparison; 14 nucleotides (current common methods) and 26 nucleotides (proposed method).
- Statistical analysis includes; Levene's test to assess the equality of variances, T-test for the means and
  Mann–Whitney U test (not shown).
- Visual representation; DataFrame of bacterial Data. Boxplot depicting the data, Linear regression depiction of correlation

# M.abs-library-preparation-Analysis
abs_r_analysis.py - in-silico analysis of experimental results. Transposon insertion in M. abs (ATCC 19977). 
Smooth to Rough transition associated locus. 
- Work layout follows a simple process: Parse and extract sequences, trimming, nucleotides rearrangement and then mapping.
- Documentation is presented in the py.files, explanaton & commentary are located by the layout and the function used.
  The stages of each analysis is reported with I/O and the desired results. The analysis is elementary, accessible, very mutable
  and can change based on the desired modifications and outcomes.
- For this analysis GFF, Genome sequences and individual gene sequences files were downloaded
  from https://mycobrowser.epfl.ch/releases
- Defined function were used for a simple analysis of the experimental results; genes associated with smooth to rough transition.
  Results are FASTQ files of our M.abs (ATCC 19977) Rough morphology Library (about 200 rough m.abs colonies pooled together)
  constructed with our Tn system offered as the propsed method. 
- Hits analysis was perfomed on both stands of the gDNA, locating matches with regular expression portraying the precise
  sought out insertion sequence (IR+28N+IR). Rearrangement of the resulted hits to thier original order is a simple string
  manipulation of moving characters. 
- init_mapper fucnction is specified to find the insertion location (in bases) and locus based on the provided genome/gene file. 
- A dataFrame (mabs_df) contains the resulted information.
- Similar analysis for the non-coding sequences (inter-genetic rigions) was perfomed.
- The raw data files (reads) of the sequenced amplicons (FASTQ) are available at:
  https://www.dropbox.com/sh/u2jvlp4oimurmkp/AABMSIqWjYwcvIx3BD-1Gpe8a?dl=0
