from utils.utils import gzip_file, make_valid_dir_string, make_dir
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib as mpl
import math
from statistics import mean, median
import sys 

###############################################
#               CONSTANTS
###############################################

start_pos = 'starting_position'
rank = 'ranking'
parent_prot = 'parent_protein'
prot_name = 'name'
seq = 'sequence'
seq_len = 'sequence_length'
plot_width = 14
plot_height = 4

###############################################
#           END CONSTANTS
###############################################

###############################################
#               PRIVATE FUNCTIONS
###############################################
'''__sort_ranks

DESC:
    sort peptide ranks by the protein they're associated with
PRAMS:
    ranks:  list of peptides with positional ranking information in the form of 
            {
                'starting_position': int,
                'ranking': int,
                'parent_protein': str
            }
Outputs:
    ranks: dictionary with rank entries as lists as values of protein name
'''
def __sort_ranks(ranks):
    sorted_ranks = {}
    for rank in ranks:
        pp = rank[parent_prot]
        if pp not in sorted_ranks:
            sorted_ranks[pp] = []
        sorted_ranks[pp].append(rank)
    return sorted_ranks

'''__protein_summary

DESC:  
    give a breif summary of the protein
Inputs:
    ranks: list of dictionaries of the form 
    {
                'starting_position': int, 
                'ranking': int,
                'sequence_length': int
    }
Outputs:
    string of brief summary
'''
def __protein_summary(ranks):
    summary = 'average peptide length: {}    average rank: {}    median rank: {}    number scores > 5: {}'
    filtered_ranks = [r for r in ranks if r[rank] is not None]
    avg_len = mean([r[seq_len] for r in filtered_ranks])
    avg_rank = mean([r[rank] for r in filtered_ranks])
    med_rank = median(r[rank] for r in filtered_ranks)
    num_bad = 0
    for r in ranks:
        if r[rank] > 5:
            num_bad += 1
    return summary.format(avg_len, avg_rank, med_rank, num_bad)


###############################################
#           END PRIVATE FUNCTIONS
###############################################

'''protein_pos_ranks

DESC:
    plot ranks of peptides against protein position
Inputs:
    prot: dictionary of the form {'name': str, 'sequence': str}
    ranks: list of dictionaries of the form
            {
                'starting_position': int, 
                'ranking': int,
                'sequence_length': int
            }
            NOTE: these should all be associated with the parent protein passed
                in. No check is done for this
kwargs:
    save_dir: str name of directory to save plot to. Default='./'
    show: bool whether or not to show the plot generated. Default=False
    compress: bool compress the plot after its generated. Default=True
Outputs:
    None
'''
def protein_pos_ranks(prot, ranks, save_dir='./', show=False, compress=True):
    rank_dist = [[] for _ in range(len(prot[seq]))]
    # TODO: this function is a temporary fix to an issue where ranks are none
    func = lambda r: r if r is not None else 0
    for r in ranks:
        rank_dist[r[start_pos]].append(func(r[rank]))

    for l in rank_dist:
        if len(l) == 0:
            l.append(0)

    # plot median rank as positions
    plt.figure(figsize=(plot_width, plot_height))
    median_ranks = zip([median(x) for x in rank_dist], [x for x in range(len(prot[seq]))])
    to_scatter = [x for x in median_ranks if x[0] != 0]
    x = [s[1] for s in to_scatter]
    y = [s[0] for s in to_scatter]

    plt.scatter(x, y)

    plt.xlabel('protein sequence position')
    plt.xlim(0, len(prot[seq]) + 1)
    plt.xticks([i for i in range(0, len(prot[seq]), 50)], [i for i in range(0, len(prot[seq]), 50)])
    plt.ylabel('median rank')
    plt.ylim(0, max(y) + 1)

    plt.annotate(__protein_summary(ranks), xy=(0.05, 0.95), xycoords='axes fraction')
    plt.title(prot[prot_name])
    
    True and plt.show()
    file_name = save_dir + prot[prot_name] + '.png'
    plt.savefig(file_name)
    plt.close()
    compress and gzip_file(file_name)
    

'''prots_pep_pos_rankings

DESC:
    Organize ranks and proteins to plot rank information against protein positin
Inputs:
    prots: list of proteins of the form {'name': str, 'sequence': str}
    ranks: list of peptides with positional ranking information in the form of 
            {
                'starting_position': int,
                'ranking': int,
                'parent_protein': str
            }
kwargs:
    save_dir: str name of directory to save plots at. Default='./'
    show_all: bool show all the plots generated. Default=False
    compress: bool compress each subplot to save room. Default=True
Outputs:
    None
'''
def prots_pep_pos_rankings(prots, ranks, save_dir='./', show_all=False, compress=True): 
    save_dir = make_valid_dir_string(save_dir + 'protein_position_rankings') 
    make_dir(save_dir)

    ranks = __sort_ranks(ranks)
    num_prots = len(prots)
    for i, prot in enumerate(prots):
        print('Plotting protein {}/{} [{}%]\r'.format(i, num_prots, int( (float(i) / float(num_prots)) * 100)), end='')
        name = prot[prot_name]
        if name not in ranks:
            print('No peptide scores with parent protein {} found. Skipping.'.format(name))
        protein_pos_ranks(prot, ranks[name], save_dir=save_dir, show=show_all, compress=compress)
        