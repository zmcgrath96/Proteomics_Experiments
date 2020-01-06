import os
import argparse
import json
from time import time
import database
from spectra import gen_spectra_files
from scoring import score_peptides
from sequence_generation import peptides
from analysis.plotting import plot_experiment
from analysis.analyze_experiment import analyze
from utils import __file_exists, __make_valid_dir_string, __make_dir
from protein_utils import read_proteins
from sequence_generation import generate_hybrids

''' old hybrid "sequence": "ALYLVCGELYTSRV", 
    second hybid: GFFYTPKEANIR
    IMPORT DEFAULTS
'''
cwd = '/'.join(str(os.path.dirname(os.path.realpath(__file__))).split('/')[:-1])
default_json_file = cwd + '/defaults.json'
defaults = None
with open(default_json_file, 'r') as o:
    defaults = json.load(o)
'''
    END IMPORT DEFAULTS
'''

def main(args):
    '''
    STEPS
    1. Generate hybrid proteins and sequences (or load from file)
    2. Load proteins and generate (or load) peptides from digestion of proteins
    3. Save all peptides as .FASTA files
    4. For each protein, create k-mers and save to mzml files
    5. Score all the k-mers against all the peptides
    6. Perform analysis on these scores
    '''
    prot_file = args.protein_sequences
    hyb_pep = args.hyb_pep
    hyb_prot = args.hyb_prot
    num_peptides = args.num_peptides
    num_hybs = args.num_hybrids
    agg_func = args.agg_func
    show_all = args.show_all
    save_dir = __make_valid_dir_string(args.save_dir)
    __make_dir(save_dir)
    min_length = args.min_length
    max_length = args.max_length
    top_n = args.top_n
    n = args.n
    m_func = args.m_func
    old_digest = args.d_file 
    mix = args.mix
    compress = args.compress

    start_time = time()

    # load in a list of proteins from a source file
    prots = None 
    if '.csv' in prot_file or '.fasta' in prot_file:
        prots, dups = read_proteins.from_csv(prot_file) if '.csv' in prot_file else read_proteins.from_fasta(prot_file) 
        print('Number of duplicates found: {}'.format(len(dups)))
    else:
        raise Exception('Protein file should be csv or fasta. File passed in: {}'.format(prot_file))

    # create peptides
    non_hybrid_peps = peptides.get_peptides(prots=prots, number_peptides=num_peptides, min_length=min_length, max_length=max_length, save_dir=save_dir, peptide_file=old_digest)

    # create hybrids
    hyb_peps, hyb_prots = generate_hybrids.generate_hybrids(hyb_pep_file=hyb_pep, hyb_prot_file=hyb_prot, prots=prots, num_gen=num_hybs, min_length=min_length, max_length=max_length, save_dir=save_dir)

    # combine them for later use
    print('Combining hybrid and non hybrid lists...')
    all_proteins_raw = prots + hyb_prots
    all_proteins_cleaned = prots + [{'name': x['name'] , 'sequence': x['protein']} for x in hyb_prots]
    all_peptides_raw = non_hybrid_peps + hyb_peps
    all_peptides_cleaned = [{'name': x['peptide_name'], 'sequence': x['peptide_sequence']} for x in non_hybrid_peps] + [{'name': x['peptide_name'], 'sequence': x['peptide_sequence']} for x in hyb_peps]
    print('Done')

    # create database files
    print('Generating fasta databases...')
    fasta_databases = database.generate(all_peptides_cleaned, save_dir=save_dir)
    print('Done')

    # create spectrum files
    print('Generating spectra files...')
    spectra_files = gen_spectra_files.generate(all_proteins_cleaned, defaults['window_sizes'], save_dir=save_dir, compress=compress)
    print('Done.')

    # run scoring algorithm on database and k-mers
    print('Scoring...')
    score_output_files = score_peptides.score_peptides(spectra_files, fasta_databases, defaults['crux_cmd'], save_dir, compress=compress)
    print('Done.')

    # save scores to json
    protein_names = []
    print('Analyzing Experiment...')
    exp_json_path = analyze(all_proteins_raw, all_peptides_raw, score_output_files, {**vars(args), **defaults}, predicting_agg_func=agg_func, saving_dir=save_dir, mix_in_hybrids=mix, show_all=show_all)
    print('Done.')

    print('Finished experiment. Time to complete: {} seconds'.format(time() - start_time))

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Entry file for the database experiments')
    parser.add_argument('protein_sequences', type=str, metavar='P', help='Path to a file (csv or fasta) with protein name and sequences')
    parser.add_argument('--hybrid-peptide-file', dest='hyb_pep', type=str, default='', help='Path to a hybrid peptide file. If none is given, new hybrid peptides generated. Default=')
    parser.add_argument('--hybrid-protein-file', dest='hyb_prot', type=str, default='', help='Path to a hybrid protein file. If none is given, new hybrid proteins generated. Default=')
    parser.add_argument('--num-peptides', dest='num_peptides', type=int, default=50, help='Number of peptides to generate as the fake sample. Default=50')
    parser.add_argument('--num-hybrids', dest='num_hybrids', type=int, default=10, help='Number of hybrid proteins and peptides to generate. Default=10')
    parser.add_argument('--aggregate-function', dest='agg_func', type=str, default='sum', help='Which aggregation function to use for combining k-mer scores. Pick either sum or product. Default=sum')
    parser.add_argument('--show-all-graphs', dest='show_all', type=bool, default=False, help='Show all the graphs generated. Will save to directory either way. Default=False.')
    parser.add_argument('--output-dir', dest='save_dir', type=str, default='./', help='Directory to save all figures. Default=./')
    parser.add_argument('--min-length', dest='min_length', type=int, default=3, help='Minimum length peptide to create. Default=3')
    parser.add_argument('--max-length', dest='max_length', type=int, default=20, help='Maximum length peptide to create. Cuts from N terminus (left) side. Default=20')
    parser.add_argument('--top-n', dest='top_n', type=bool, default=False, help='When recording how well a peptide scores against a protein, only use the top n proteins. Default=False')
    parser.add_argument('--n', dest='n', type=int, default=5, help='n to use if using --top-n. Default=5')
    parser.add_argument('--measure-func', dest='m_func', type=str, default='average', help='Measuring function for determining the top n proteins. Options are: sum, average, max. Default=average')
    parser.add_argument('--peptide-file', dest='d_file', type=str, default='', help='Peptides from a past experiment. Default=None')
    parser.add_argument('--mix-prots', dest='mix', type=bool, default=False, help='Whether or not to also use huybrid proteins when calculating scores. Default=False')
    parser.add_argument('--compress', dest='compress', type=bool, default=True, help='Compress spectra files while generating them. Default=True')
    args = parser.parse_args()
    main(args)
    