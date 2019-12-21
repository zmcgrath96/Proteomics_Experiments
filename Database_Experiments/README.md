# Flipped Database Query

### Hypothesis
The "flipped" aspect to this experiment is only possible because this is all theoretical. We are essentially
using scoring algorithms in a manner in the correct way because we 

Using k-mers from protein sequences should allow us to do two things: 
1. Accurately identify peptides from one parent protein
2. Identify if a peptide could sourced from two parent protiens

Ideally, as our k-mer lengths envelope our spectra from MS/MS, we should see the score of this comparison peak.
Some aggregation function applied to all the k-mers at the same starting subsequence position should result
in an even more identifiable peak, making it easier to identify parent proteins.

For shorter peptides, we would expect their to be a lot of noise, from which more sophisticated searches should 
spawn in order to accurately identify a parent protein. 

This method should be effective for identifying potential hybrid peptides for the following reasons:
* since we deal with k-mers, as k-mers overlap the spectra, we should get good peaks from 2 parents
* knowing the k-mer size allows us to find the junction site as we have information from sets of k-mers

### Methods
1. Given a list of proteins, 
   1. either randomly generate or load peptides from memory
   2. generate hybrid proteins 
   3. generate hybrid peptides from the hybrid proteins
2. Save all peptides to .FASTA files as the future databases
3. From the proteins, generate k-mers with +1 offset (move the window by 1 position)
4. Calculate the theoretical masses of these k-mers from their residual amino acid masses
4. Save these k-mer spectra to .mzML files
5. Use a scoring algorithm<sup>1</sup> to compare these spectra to the peptides in the .FASTA files
6. Load the scores into memory and format them in a human-readable json format
7. Use this to run analyses on
   1. Plot the scores of k-mers for each peptide and plot aggregate of protein k-mers for many proteins against each peptide
   2. Rank how well the correct starting position for k-mers scores and run statistics on this


[1] crux search tool was used. it is an implementation of the SEQUEST scoring algorithm. Link to and citation in the top level README file


### Download and Usage
To download 
```bash
~>git clone https://github.com/zmcgrath96/Proteomics_Experiments.git
```

To run 
```bash
~>cd Proteomics_Experiments/Database_Experiments/src
Proteomics_Experiments/Database_Experiments/src> main.py <options>
```

For help
```bash
~>cd Proteomics_Experiments/Database_Experiments/src
Proteomics_Experiments/Database_Experiments/src> main.py --help
```

Usage
```bash
usage: main.py [-h] [--hybrid-peptide-file HYB_PEP]
               [--hybrid-protein-file HYB_PROT] [--num-peptides NUM_PEPTIDES]
               [--num-hybrids NUM_HYBRIDS] [--aggregate-function AGG_FUNC]
               [--show-all-graphs SHOW_ALL] [--output-dir SAVE_DIR]
               [--min-length MIN_LENGTH] [--max-length MAX_LENGTH]
               [--top-n TOP_N] [--n N] [--measure-func M_FUNC]
               [--peptide-file D_FILE] [--mix-prots MIX]
               P

Entry file for the database experiments

positional arguments:
  P                     Path to a file (csv or fasta) with protein name and
                        sequences

optional arguments:
  -h, --help            show this help message and exit
  --hybrid-peptide-file HYB_PEP
                        Path to a hybrid peptide file. If none is given, new
                        hybrid peptides generated. Default=
  --hybrid-protein-file HYB_PROT
                        Path to a hybrid protein file. If none is given, new
                        hybrid proteins generated. Default=
  --num-peptides NUM_PEPTIDES
                        Number of peptides to generate as the fake sample.
                        Default=50
  --num-hybrids NUM_HYBRIDS
                        Number of hybrid proteins and peptides to generate.
                        Default=10
  --aggregate-function AGG_FUNC
                        Which aggregation function to use for combining k-mer
                        scores. Pick either sum or product. Default=sum
  --show-all-graphs SHOW_ALL
                        Show all the graphs generated. Will save to directory
                        either way. Default=False.
  --output-dir SAVE_DIR
                        Directory to save all figures. Default=./
  --min-length MIN_LENGTH
                        Minimum length peptide to create. Default=3
  --max-length MAX_LENGTH
                        Maximum length peptide to create. Cuts from N terminus
                        (left) side. Default=20
  --top-n TOP_N         When recording how well a peptide scores against a
                        protein, only use the top n proteins. Default=False
  --n N                 n to use if using --top-n. Default=5
  --measure-func M_FUNC
                        Measuring function for determining the top n proteins.
                        Options are: sum, average, max. Default=average
  --peptide-file D_FILE
                        Peptides from a past experiment. Default=None
  --mix-prots MIX       Whether or not to also use huybrid proteins when
                        calculating scores. Default=False
```


Mass source: http://education.expasy.org/student_projects/isotopident/htdocs/aa-list.html