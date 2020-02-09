import utils 
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
max_x_len = 75
plot_width = 14
subplot_height = 2

colors = [
            '#FFFFFF', '#EF7F37', '#EFA737', '#EFD037', '#EFD937', '#E4EF37', 
            '#BDEF37', '#88EF37', '#37EF58', '#37EFA1', '#37EFCE',
            '#37D9EF', '#37B4EF', '#3769EF', '#4D37EF', '#7737EF',
            '#A137EF', '#CE37EF', '#EF37E7', '#EF37B2', '#EF3772'
        ]
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

'''__avg_len_by_pos

DESC:
    calculate the average length of peptides found at a position
PARAMS:
    ranks: list of dictionaries
    prot_length: int length of the protein sequence
RETURNS:
    list of ints of average length of peptide by position
'''
def __avg_len_by_pos(ranks, prot_length):
    list_lens = [[] for _ in range(prot_length)]
    avg_lens = [0 for _ in range(prot_length)]
    for r in ranks:
        list_lens[r[start_pos]].append(r[seq_len])
    for i in range(prot_length):
        if list_lens[i] is None or list_lens[i] == []:
            continue
        avg_lens[i] = mean(list_lens[i])
    return avg_lens

'''__protein_summary

DESC:  
    give a breif summary of the protein
PARAMS:
    ranks: list of dictionaries of the form 
    {
                'starting_position': int, 
                'ranking': int,
                'sequence_length': int
    }
RETURNS:
    string of brief summary
'''
def __protein_summary(ranks):
    summary = 'average peptide length: {}    average rank: {}    median rank: {}'
    filtered_ranks = [r for r in ranks if r[rank] is not None]
    avg_len = mean([r[seq_len] for r in filtered_ranks])
    avg_rank = mean([r[rank] for r in filtered_ranks])
    med_rank = median(r[rank] for r in filtered_ranks)
    return summary.format(avg_len, avg_rank, med_rank)


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
                'ranking': int,
                'sequence_length': int
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
    # NOTE: this function is a temporary fix to an issue where ranks are none
    fun = lambda r: r if r is not None else 0
    for i, r in enumerate(ranks):
        rank_dist[r[start_pos]].append(fun(r[rank]))

    for l in rank_dist:
        if len(l) == 0:
            l.append(0)

    #calculate the averate subsequence length at each position
    avg_lens = __avg_len_by_pos(ranks, len(prot[seq]))

    # for longer proteins, the length is an issue. Split into subplots
    num_subplots = math.ceil(len(prot[seq]) / max_x_len)
    gs = gridspec.GridSpec(num_subplots, 2, width_ratios=[plot_width, .3]) 
    fig = plt.figure(figsize=(plot_width, subplot_height*num_subplots))
    
    # divide the figure into subplots where the width is limited to constant of max_x_len
    for i in range(num_subplots):  
        this_axis = plt.subplot(gs[i, 0])
        # get the data needed for this row
        end = min([((i+1) * max_x_len), len(prot[seq])])
        divided_seq = prot[seq][i*max_x_len:end]
        divided_rank_dist = rank_dist[i*max_x_len:end]
        divided_avg_len = avg_lens[i*max_x_len:end]

        # plot labeling, plotting, coloring
        bplot = this_axis.boxplot(divided_rank_dist)
        for attr in ['boxes', 'whiskers', 'caps', 'fliers', 'medians']:
            for n, item in enumerate(bplot[attr]):
                t = n if attr != 'whiskers' and attr != 'caps' else math.floor(n/2)
                color_to_use = colors[divided_avg_len[t]]
                item.set(color=color_to_use)

        this_axis.set_xticks([j for j in range(len(divided_rank_dist))])
        this_axis.set_xticklabels(divided_seq)
        y_label = 'sequence:  {} - {}'.format(i * max_x_len, end - 1)
        this_axis.set_ylabel(y_label, fontsize=8)
        
    # add color bar
    cmap = mpl.colors.ListedColormap(colors)
    bounds = [i for i in range(21)]
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    color_axis = plt.subplot(gs[:, 1])
    cb = mpl.colorbar.ColorbarBase(color_axis, cmap=cmap, norm=norm)
    cb.set_label('mean peptide length')

    # give a common label for the y axis as the protein name
    fig.text(0.06, 0.5, prot[prot_name], ha='center', va='center', rotation='vertical')
    fig.suptitle('peptide rank distributions vs protein position')
    fig.text(.25, 0.005, __protein_summary(ranks))
    show and plt.show()
    file_name = save_dir + prot[prot_name] + '.png'
    plt.savefig(file_name)
    plt.close()
    compress and utils.__gzip(file_name)
    

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
        