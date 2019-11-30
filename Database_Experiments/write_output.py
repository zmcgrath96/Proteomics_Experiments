from utils import __make_valid_text_name

SUMMARY_HEADER = '''
###############################################
                   SUMMARY 
###############################################
'''

'''__save_plot_data

DESC: 
    save data from a plot to a text file
PARAMS:

'''
def __save_plot_data(file, data, title):
    file = __make_valid_text_name(file)
    if type(data) is not dict: 
        raise Exception('data should be a dictionary')
    with open(file, 'w') as o:
        o.write(str(title) + '\n')
        o.write(SUMMARY_HEADER)
        protein_line = 'proteins: ' + ' '.join(key for key in data) 
        o.write(protein_line)