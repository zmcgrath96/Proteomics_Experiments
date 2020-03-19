import matplotlib.pyplot as plt
import itertools, os
from utils.utils import make_valid_dir_string, make_dir, gzip_dir
from file_io import JSON

####################################################
#                   CONSTANTS
####################################################
plot_colors = [
    'r', 'b', 'g', 'y', 'k', 'c', 'm'
]
plot_markers = [
    ' ', '--', '-.', ':'
]
all_line_types = [''.join(x) for x in list(itertools.product(plot_colors, plot_markers))]
all_line_types.sort(key=lambda x: x[1])

hybrid_prefix = 'HYBRID'

json_header = 'experiment_info'
json_header_prots = 'proteins'
json_header_peps = 'peptides'
json_exp = 'experiment'
json_analysis =  'analysis'
json_rankings = 'ranks'
json_rankings_rankings = 'ranks'
json_seq_len = 'sequence_length'
###################################################
#               END CONSTANTS
###################################################
###################################################
#              PRIVATE FUNCTIONS
###################################################

'''__collect_k_rankings

DESC:
    Collect the score rankings for each peptide for each k
Inputs:
    exp: dictionary experiment summary dictionary
Outputs:
    dictionary of the form
    {
        k: {
            l_0: [], 
            l_1: [],
            ...
        }
    }
'''
def __collect_k_rankings(exp):
    rankings = {}
    for _, pep in exp[json_exp].items():
        l = pep[json_analysis][json_rankings][json_seq_len]
        for k, rank in pep[json_analysis][json_rankings][json_rankings_rankings].items():
            if not k in rankings:
                rankings[k] = {}
            if not l in rankings[k]:
                rankings[k][l] = []
            rankings[k][l] += [rank]

    return rankings

###################################################
#           END PRIVATE FUNCTIONS
###################################################

'''plot_subsequence_vs_protein

DESC:
    generate and save a subsequence vs a protein graph
Inputs:
    k_mers: dictionary of lists of floats to plot. entry names used for legend
kwargs:
    title: string to label this plot. Also saving name for the plot. Default=''
    save_dir: string name of directory to save all output under. Default=./
    show_graph: bool whether or not to show the graph. Default=False
    save_raw_json: bool whether or not to save to raw json data. Default=False
'''
def plot_subsequence_vs_protein(k_mers, title='', save_dir='./', show_graph=False, save_raw_json=False):
    save_dir = make_valid_dir_string(save_dir) + title + '/'
    make_dir(save_dir)

    plt.figure(figsize=(10, 7))
    i  = 0
    for mer in k_mers:
        score = k_mers[mer]
        plt.plot([j for j in range(len(score))], score, all_line_types[i].strip(), label=str(mer))
        i += 1

    plt.xlabel('subsequence start position')
    plt.ylabel('k-mer scores')
    plt.title(title)
    plt.legend()
    plt.savefig(save_dir + title)
    show_graph and plt.show()
    save_raw_json and JSON.save_dict(save_dir + title, k_mers)
    plt.close()


'''plot_subsequence

DESC:
    generate and save a subsequence against the aggregation scores
Inputs:
    aggs: dictionary of list of aggregation scores. Entry names used for legend
kwargs:
    title: string to label this plot. Also saving name for the plot. Default=''
    save_dir: string name of directory to save all output under. Default =./
    show_graph: bool whether or not to show the graph. Default=False
    peaks: list of dictionaries of the form
            [{
                'protein_name': str,
                'position': int,
                'score': float
            }]
            these are the points of interest to mark
    save_raw_json: bool save the raw json data for each graph. Default=False
    compress: bool to compress the directory. Default=True
'''
def plot_subsequence(aggs, title='', save_dir='./', show_graph=False, agg_func='sum', sequence_info=None, save_raw_json=False, compress=True):
    save_dir = make_valid_dir_string(save_dir)
    make_dir(save_dir)

    plt.figure(figsize=(10, 7))
    i = 0
    for agg in aggs:
        score = aggs[agg]
        if len(score) == 0:
            continue
        plt.plot([j for j in range(len(score))], score, all_line_types[i].strip(), label=str(agg))
        i += 1

    plt.xlabel('subsequence start position')
    plt.ylabel('{} of k-mer scores'.format(agg_func))
    plt.legend()
    plt.title(title)
    plt.savefig(save_dir + title)
    show_graph and plt.show()
    save_raw_json and JSON.save_dict(save_dir + title, aggs)
    plt.close()
    compress and gzip_dir(save_dir)


'''plot_score_rankings

DESC:
    Generates violin plots of ranking history
Inputs:
    exp: dictionary experiment summary dictionary
kwargs:
    save_dir: str path to directory to save images too. Default=./
    show_all: bool to show or not all plots. Default=False
Outputs:
'''
def plot_score_rankings(exp, save_dir='./', show_all=False):
    # go through each peptide and collect stats
    # x axis: length of the sequnce
    # y axis: violin plot of distribution
    # do it for each k or aggregation
    save_dir = make_valid_dir_string(save_dir)
    make_dir(save_dir)

    rs = __collect_k_rankings(exp)
    for k in rs:
        pos = []
        data =[]
        for t in rs[k]:
            pos.append(int(t)) 
            # clean up. Remove None values from lists
            d = [x for x in rs[k][t] if x is not None]
            data.append(d)
        
        if not isinstance(data, list) or not all(isinstance(x, list) for x in data) or not all(isinstance(x, int) for y in data for x in y): 
            print('ERROR: data is not in correct format ( [[int]] ). \ndata: {}'.format(data))
            continue
        if not isinstance(pos, list) or not all(isinstance(x, int) for x in pos): 
            print('ERROR: pos is not in correnct format ( [int] ). \npos: {}'.format(pos))
            continue
   
        plt.violinplot(data, pos)
        plt.title(k)
        plt.xlabel('subsequence length')
        plt.ylabel('score distribution')
        show_all and plt.show()
        plt.savefig(save_dir + str(k))
        plt.close()