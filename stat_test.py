import pandas as pd
import scipy.stats as stats 
from os import listdir
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt


p_values = np.zeros((4,4))
num_sig = [''*4 for i in range(4)]
intron_files = listdir('processed_data/proteins_of_interest_excel_files/proportion_intron_content_files')
exon_files = listdir('processed_data/proteins_of_interest_excel_files/proportion_exon_content_files')
intron_heat = listdir('processed_data/analysis/intron_length_heat_map')
exon_heat = listdir('processed_data/analysis/exon_length_heat_map')
total_files = listdir('processed_data/proteins_of_interest_excel_files/total_files')
exon_total = total_files[0]
print(len(intron_heat))
print(len(intron_files))
print(len(exon_heat))
print((len(exon_files)))

exon_total_data = pd.read_excel('processed_data/proteins_of_interest_excel_files/total_files/total_proportion_exon_proteins_of_interest.xlsx')
exon_total_data['origin'].replace({1:'bacterial', 2:'archaeal', 3:'eukaryotic', 4:'universal'})

intron_total_data = pd.read_excel('processed_data/proteins_of_interest_excel_files/total_files/total_intron_content_proteins_of_interest.xlsx')
intron_total_data['origin'].replace({1:'bacterial', 2:'archaeal', 3:'eukaryotic', 4:'universal'})

#print(stats.kruskal(exon_total_data['exon_lengths'],intron_total_data['intron_lengths']))
# p_values[0][0] = None
# p_values[1][0] = stats.kruskal(exon_total_data['proportion_exon'][exon_total_data['origin'] == 'bacterial'], exon_total_data['proportion_exon'][exon_total_data['origin'] == 'archaeal'])[1]
# p_values[2][0] = stats.kruskal(exon_total_data['proportion_exon'][exon_total_data['origin'] == 'bacterial'], exon_total_data['proportion_exon'][exon_total_data['origin'] == 'eukaryotic'])[1]
# p_values[3][0] = stats.kruskal(exon_total_data['proportion_exon'][exon_total_data['origin'] == 'bacterial'], exon_total_data['proportion_exon'][exon_total_data['origin'] == 'universal'])[1]
# #print('archaeal')
# p_values[0][1] = None
# p_values[1][1] = None
# p_values[2][1] = stats.kruskal(exon_total_data['proportion_exon'][exon_total_data['origin'] == 'archaeal'], exon_total_data['proportion_exon'][exon_total_data['origin'] == 'eukaryotic'])[1]
# p_values[3][1] = stats.kruskal(exon_total_data['proportion_exon'][exon_total_data['origin'] == 'archaeal'], exon_total_data['proportion_exon'][exon_total_data['origin'] == 'universal'])[1]
# #print('eukaryotic')
# p_values[0][2] = None
# p_values[1][2] = None
# p_values[2][2] = None
# p_values[3][2] = stats.kruskal(exon_total_data['proportion_exon'][exon_total_data['origin'] == 'eukaryotic'], exon_total_data['proportion_exon'][exon_total_data['origin'] == 'universal'])[1]

# p_values[0][3] = None
# p_values[1][3] = None
# p_values[2][3] = None
# p_values[3][3] = None

# #print(p_values)
# fig = plt.subplots(figsize=(7.25, 5))
# plt.title('Total Proportion of Exons P values by Origin', fontsize=11)
# ax = sns.heatmap(p_values, annot=True, xticklabels = ['Bacterial', 'Archaeal', 'Eukaryotic', 'Universal'], linewidth=0.2, linecolor='#dcdfe3', clip_on=False)
# sns.dark_palette('#69d', reverse=True, as_cmap=True)
# ax.set_yticklabels(['Bacterial', 'Archaeal', 'Eukaryotic', 'Universal'], rotation=0, va='center')
# plt.savefig('total_proportion_exon_heat_map.png')
# plt.close('all')

for f in range(len(exon_files)):
    print(exon_files[f])
    exon_data = pd.read_excel('processed_data/proteins_of_interest_excel_files/proportion_exon_content_files/' + '_'.join(exon_files[f].split('exon')[:-1]) + 'exon_content_proteins_of_interest.xlsx')
    exon_data['origin'].replace({1:'bacterial', 2:'archaeal', 3:'eukaryotic', 4:'universal'})
    if all(exon_data['proportion_exon'][exon_data['origin'] == 'bacterial'] == 1):
        pass
    else:
        p_values[1][0] = stats.kruskal(exon_data['proportion_exon'][exon_data['origin'] == 'bacterial'], exon_data['proportion_exon'][exon_data['origin'] == 'archaeal'])[1]
        p_values[2][0] = stats.kruskal(exon_data['proportion_exon'][exon_data['origin'] == 'bacterial'], exon_data['proportion_exon'][exon_data['origin'] == 'eukaryotic'])[1]
        p_values[3][0] = stats.kruskal(exon_data['proportion_exon'][exon_data['origin'] == 'bacterial'], exon_data['proportion_exon'][exon_data['origin'] == 'universal'])[1]
    if all(exon_data['proportion_exon'][exon_data['origin']=='archaeal'] == 1):
        pass
    else:
        p_values[2][1] = stats.kruskal(exon_data['proportion_exon'][exon_data['origin'] == 'archaeal'], exon_data['proportion_exon'][exon_data['origin'] == 'eukaryotic'])[1]
        p_values[3][1] = stats.kruskal(exon_data['proportion_exon'][exon_data['origin'] == 'archaeal'], exon_data['proportion_exon'][exon_data['origin'] == 'universal'])[1]
    if all(exon_data['proportion_exon'][exon_data['origin'] == 'eukaryotic'] == 1):
        pass
    else:
        p_values[3][2] = stats.kruskal(exon_data['proportion_exon'][exon_data['origin'] == 'eukaryotic'], exon_data['proportion_exon'][exon_data['origin'] == 'universal'])[1]

    p_values[0][2] = None
    p_values[1][2] = None
    p_values[2][2] = None
    p_values[0][1] = None
    p_values[1][1] = None
    p_values[0][0] = None
    p_values[0][3] = None
    p_values[1][3] = None
    p_values[2][3] = None
    p_values[3][3] = None

    for i in range(len(num_sig[0])-1):
        for j in range(len(num_sig[0])-1):
            if p_values[i][j] <= 0.05:
                num_sig[i][j].append('_'.join(exon_files[f].split('exon')[:-1]))
    fig = plt.subplots(figsize=(7.25, 5))
    plt.title('P values of Exon Content by Origin for ' + '_'.join(exon_files[f].split('exon')[:-1]), fontsize=11)
    ax = sns.heatmap(p_values, annot=True, xticklabels = ['Bacterial', 'Archaeal', 'Eukaryotic', 'Universal'], linewidth=0.2, linecolor='#dcdfe3', clip_on=False)
    sns.dark_palette('#69d', reverse=True, as_cmap=True)
    ax.set_yticklabels(['Bacterial', 'Archaeal', 'Eukaryotic', 'Universal'], rotation=0, va='center')
    #plt.savefig('_'.join(exon_files[f].split('num')[:-1]) + 'num_exons_heat_map.png')
    plt.close('all')
print(num_sig)