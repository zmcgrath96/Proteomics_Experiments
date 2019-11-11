# Database Experiments
The experiments in this folder involve different ways of constructing databases and queries in order to identify hybrid peptides. The basic overview for these experiments:

1. Create a known hybrid from known parents
2. Break each parent into peptides of a varying length
3. Track how each of these queries perform

The idea is that there should be some sort of identifiable trend in score based on the peptide lengths

The two default protiens used are:

1. Human Insulin: https://www.uniprot.org/uniprot/P01308
2. Human Brain-derived neurotrophic factor: https://www.uniprot.org/uniprot/P23560

The default hybrid is:

   ALYLVCGE-LYTSRV  
   INS[37:44]-BDNF[88:93]

The parents and hybrid can be changed in the sequences.json file

## General Methods
### Generating Spectra 
### Scoring

# Flipped Database Query

### Hypothesis
Flipping the database and the query strings should allow us to see how well our Amino Acid (AA) chunks score against a hybrid peptide. We should be able to see any trends in scores for these chunks. These trends should be:
1. For smaller cuts, we should score very high with lots of false positives 
2. For cuts about the size of the contribution to the hybrid, we should see a rise in score above the noise level, peak, then taper back down
3. For cuts larger than the contribution, we should see a rise as we overlap the hybrid, but not peak as well as number 2

The chunk sizes used in the original runs were 3, 4, 6, 8, 10. 

### Methods
1. Create a hybrid sequence from two known parent proteins
2. Put this sequence in FASTA format and load it as the database into a database scoring algorithm
3. For both of the parent proteins do the following:
   1. Use a sliding window to cut protein into smaller sequences
   2. Use each of these sequences as the query against our database of 1
   3. Record how well each of these sequences score against it
   4. Repeat this for varying sizes of windows (¼, ½, 1, >1 size the contribution of each parent protein to the hybrid)

### Results

# Fractionated Database Query
### Hypothesis
Creating varying types of databases made from the parent proteins of the sequences, we should be able to see trends in the score of sequences in the database. These trends should be:
1.  For sequences smaller than the contribution to the hybrid, there should be a lot of noise with no real contender for the source of the hybrid 
2. For sequences about the length of the contribution from the parent protein, we should see a rise in score as we cover the full contribution, then a taper back down
3. For sequences longer than the contribution from the parent, we should see the same as we saw for number 2

The chunk sizes used in the original runs were 3, 4, 6, 8, 10. 
### Methods
1. Create a hybrid sequence from two known parent proteins
2. Create a new database (from the list below)
3. Query the hybrid against each of these databases 
4. Track the scores the hybrid gets across all of the sequences in the database
##### Database List
- Small cuts (about a quarter of the size of the contribution to the hybrid) from sliding window of parent protein 1 mixed with small cuts (about a quarter of the size of the contribution to the hybrid) from sliding window of parent protein 2
- Small cuts (about half of the size of the contribution to the hybrid) from sliding window of parent protein 1 mixed with small cuts (about half of the size of the contribution to the hybrid) from sliding window of parent protein 2
- Medium cuts (about the size of the contribution to the hybrid) from sliding window of parent protein 1 mixed with medium cuts (about the size of the contribution to the hybrid) from sliding window of parent protein 2
- Large cuts (larger than the size of the contribution to the hybrid) from sliding window of parent protein 1 mixed with large cuts (larger than the size of the contribution to the hybrid) from sliding window of parent protein 2


### Results



Mass source: http://education.expasy.org/student_projects/isotopident/htdocs/aa-list.html