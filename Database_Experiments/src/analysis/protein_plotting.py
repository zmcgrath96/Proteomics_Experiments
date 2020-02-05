import utils 
import matplotlib.pyplot as plt

###############################################
#               CONSTANTS
###############################################

start_pos = 'starting_position'
rank = 'ranking'
parent_prot = 'parent_protein'
prot_name = 'name'
seq = 'sequence'

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
RETURNS:
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

###############################################
#           END PRIVATE FUNCTIONS
###############################################

'''protein_pos_ranks

DESC:
    plot ranks of peptides against protein position
PARAMS:
    prot: dictionary of the form {'name': str, 'sequence': str}
    ranks: list of dictionaries of the form
            {
                'starting_position': int, 
                'ranking': int
            }
            NOTE: these should all be associated with the parent protein passed
                in. No check is done for this
OPTIONAL:
    save_dir: str name of directory to save plot to. Default='./'
    show: bool whether or not to show the plot generated. Default=False
RETURNS:
    None
'''
def protein_pos_ranks(prot, ranks, save_dir='./', show=False):
    rank_dist = [[] for _ in range(len(prot[seq]))]
    for r in ranks:
        rank_dist[r[start_pos]].append(r[rank])

    for l in rank_dist:
        if len(l) == 0:
            l.append(0)

    title = 'position scores vs {}'.format(prot[prot_name])
    plt.figure(figsize=(10, 7))
    plt.violinplot(rank_dist, [i for i in range(len(prot[seq]))])
    plt.xticks([i for i in range(len(prot[seq]))], prot[seq])
    plt.title(title)
    show and plt.show()
    plt.savefig(save_dir + title)
    

'''prots_pep_pos_rankings

DESC:
    Organize ranks and proteins to plot rank information against protein positin
PARAMS:
    prots: list of proteins of the form {'name': str, 'sequence': str}
    ranks: list of peptides with positional ranking information in the form of 
            {
                'starting_position': int,
                'ranking': int,
                'parent_protein': str
            }
OPTIONAL:
    save_dir: str name of directory to save plots at. Default='./'
    show_all: bool show all the plots generated. Default=False
RETURNS:
    None
'''
def prots_pep_pos_rankings(prots, ranks, save_dir='./', show_all=False): 
    save_dir = utils.__make_valid_dir_string(save_dir + 'protein_position_rankings') 
    utils.__make_dir(save_dir)

    ranks = __sort_ranks(ranks)
    for prot in prots:
        name = prot[prot_name]
        if name not in ranks:
            print('No peptide scores with parent protein {} found. Skipping.'.format(name))
        protein_pos_ranks(prot, ranks[name], save_dir=save_dir, show=show_all)
        