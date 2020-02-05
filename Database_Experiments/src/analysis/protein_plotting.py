import utils 
import matplotlib.pyplot as plt
import math

###############################################
#               CONSTANTS
###############################################

start_pos = 'starting_position'
rank = 'ranking'
parent_prot = 'parent_protein'
prot_name = 'name'
seq = 'sequence'
max_x_len = 75
plot_width = 14
subplot_height = 2
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
    compress: bool compress the plot after its generated. Default=True
RETURNS:
    None
'''
def protein_pos_ranks(prot, ranks, save_dir='./', show=False, compress=True):
    rank_dist = [[] for _ in range(len(prot[seq]))]
    for r in ranks:
        rank_dist[r[start_pos]].append(r[rank])

    for l in rank_dist:
        if len(l) == 0:
            l.append(0)

    title = 'position scores vs {}'.format(prot[prot_name])
    # for longer proteins, the length is an issue. Split into subplots
    num_subplots = math.ceil(len(prot[seq]) / max_x_len)
    fig, ax = plt.subplots(nrows=num_subplots, figsize=(plot_width, subplot_height * num_subplots))
    # divide the figure into subplots where the width is limited to constant of max_x_len
    for i in range(num_subplots):
        end = min([((i+1) * max_x_len), len(prot[seq])])
        divided_seq = prot[seq][i*max_x_len:end]
        divided_ranks = rank_dist[i*max_x_len:end]
        ax[i].violinplot(divided_ranks, [j for j in range(len(divided_ranks))], showmeans=True, showmedians=True)
        ax[i].set_xticks([j for j in range(len(divided_ranks))])
        ax[i].set_xticklabels(divided_seq)
        y_label = 'sequence:  {} - {}'.format(i * max_x_len, end - 1)
        ax[i].set_ylabel(y_label, fontsize=8)

    # give a common label for the y axis as the protein name
    fig.text(0.01, 0.5, prot[prot_name], ha='center', va='center', rotation='vertical')
    fig.suptitle('peptide rank distributions vs protein position')
    show and plt.show()
    file_name = save_dir + title + '.png'
    plt.savefig(file_name)
    plt.close()
    utils.__gzip(file_name)
    

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
    compress: bool compress each subplot to save room. Default=True
RETURNS:
    None
'''
def prots_pep_pos_rankings(prots, ranks, save_dir='./', show_all=False, compress=True): 
    save_dir = utils.__make_valid_dir_string(save_dir + 'protein_position_rankings') 
    utils.__make_dir(save_dir)

    ranks = __sort_ranks(ranks)
    num_prots = len(prots)
    for i, prot in enumerate(prots):
        print('Plotting protein {}/{} [{}%]\r'.format(i, num_prots, int( (float(i) / float(num_prots)) * 100)), end='')
        name = prot[prot_name]
        if name not in ranks:
            print('No peptide scores with parent protein {} found. Skipping.'.format(name))
        protein_pos_ranks(prot, ranks[name], save_dir=save_dir, show=show_all, compress=compress)
        