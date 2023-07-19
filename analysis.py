from tree_processing import bacterial_proteins, archaeal_proteins, universal_proteins, eukaryotic_proteins
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from os import listdir, path

total_bacterial = []
total_archaeal = []
total_eukaryotic = []
total_universal = []
all_data = []

def make_hist_total(full_data, variable):
    plt.title(variable)
    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax1.set(title=('Total Number of Exons'), ylabel='Number of Exons')
    ax1.set_xticklabels(['n: '+ str(len(full_data))])
    bp = plt.boxplot(full_data, patch_artist=True)
    for tick, label in zip(range(1), ax1.get_xticklabels()):
        ax1.text((np.arange(1)+1)[tick], .97, ('m: '+[str(round(s, 3)) for s in bp['medians'][tick].get_ydata()][0]), transform=ax1.get_xaxis_transform(), horizontalalignment='center', size='small')
    plt.setp(bp['medians'], color='black')
    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    plt.savefig(('Total_' + variable + '.png'))
    plt.close('all')

def make_hist_total_by_origin(bact, arch, euk, uni, variable):
    plt.suptitle(variable)
    data = [bact,arch,euk,uni]
    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax1.set(title=('Total Number of Exons by Origin'), xlabel='Origin', ylabel='Number of Exons')
    bp = plt.boxplot(data, patch_artist=True)
    ax1.set_xticklabels([('Bacterial\nn: ' + str(len(data[0]))), ('Archaeal\nn: ' + str(len(data[1]))), ('Eukaryotic\nn: ' + str(len(data[2]))), ('Universal\nn: ' + str(len(data[3])))])
    for tick, label in zip(range(len(data)), ax1.get_xticklabels()):
        ax1.text((np.arange(len(data))+1)[tick], .97, ('m: '+[str(round(s, 3)) for s in bp['medians'][tick].get_ydata()][0]), transform=ax1.get_xaxis_transform(), horizontalalignment='center', size='small')
    plt.setp(bp['medians'], color='black')
    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    plt.savefig(('Total_' + variable + '_by_origin_boxplot.png'))
    plt.close('all')


def make_hist_by_origin(species, variable, data):
    bacterial = []
    archaeal = []
    universal = []
    eukaryotic = []
    bacterial_ids = []
    archaeal_ids = []
    universal_ids = []
    eukaryotic_ids = []
    for i in range(len(data[variable])):
        if variable == 'gc_content':
            if data['protein_id'].values[i] in bacterial_proteins:
                bacterial.append(data[variable].values[i]/100)
                bacterial_ids.append(data['protein_id'].values[i])
            elif data['protein_id'].values[i] in archaeal_proteins:
                archaeal.append(data[variable].values[i]/100)
                archaeal_ids.append(data['protein_id'].values[i])
            elif data['protein_id'].values[i] in universal_proteins:
                universal.append(data[variable].values[i]/100)
                universal_ids.append(data['protein_id'].values[i])
            elif data['protein_id'].values[i] in eukaryotic_proteins:
                eukaryotic.append(data[variable].values[i]/100)
                eukaryotic_ids.append(data['protein_id'].values[i])
        else:
            if data['protein_id'].values[i] in bacterial_proteins:
                bacterial.append(data[variable].values[i])
                archaeal_ids.append(data['protein_id'].values[i])
            elif data['protein_id'].values[i] in archaeal_proteins:
                archaeal.append(data[variable].values[i])
                archaeal_ids.append(data['protein_id'].values[i])
            elif data['protein_id'].values[i] in universal_proteins:
                universal.append(data[variable].values[i])
                universal_ids.append(data['protein_id'].values[i])
            elif data['protein_id'].values[i] in eukaryotic_proteins:
                eukaryotic.append(data[variable].values[i])
                eukaryotic_ids.append(data['protein_id'].values[i])
    plt.subplots_adjust(hspace=0.3)
    data = [bacterial,archaeal,eukaryotic,universal]
    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax1.set(title=('Number of Exons in ' + species), xlabel='Origin', ylabel='Number of Exons')
    bp = plt.boxplot(data, patch_artist=True)
    ax1.set_xticklabels([('Bacterial\nn: ' + str(len(data[0]))), ('Archaeal\nn: ' + str(len(data[1]))), ('Eukaryotic\nn: ' + str(len(data[2]))), ('Universal\nn: ' + str(len(data[3])))])
    for tick, label in zip(range(len(data)), ax1.get_xticklabels()):
        ax1.text((np.arange(len(data))+1)[tick], .97, ('m: '+[str(round(s, 3)) for s in bp['medians'][tick].get_ydata()][0]), transform=ax1.get_xaxis_transform(), horizontalalignment='center', size='small')
    plt.setp(bp['medians'], color='black')
    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    if species == 'GCF_000001405.37_GRCh38.p11':
        plt.savefig((species + '_' + variable + '_box_plot.png'))
    plt.close('all')
    return bacterial, archaeal, eukaryotic, universal, bacterial_ids + archaeal_ids + eukaryotic_ids + universal_ids

def ordinal_lengths(x):
    intron_num = []
    for i in range(len(x)):
        for j in range(x[i]):
            intron_num.append(j)
    return intron_num

def make_comparison(species, x, y):
    # x is the number of introns
    # y is the exon length
    # max is max number of introns (x range)
    # find how many introns in current protein i
    # plot all the intron lengths
    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax1.set(title=('GC Content vs Exon Length of ' + species + ' by Origin'), xlabel='Exon Length', ylabel='GC Content')
    for l in range(len(x)):
    #     intron_nums = ordinal_lengths(x[l])
        plt.scatter(x[l], y[l], s=10)

    ax1.legend(('Bacterial','Archaeal','Eukaryotic','Universal'),title='Origins')
    if species == 'GCF_000001405.37_GRCh38.p11':
        plt.savefig((species + 'gc_content_vs_exon_length.png'))
    plt.close('all')

def make_total_comparison(x, y):
    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax1.set(title=('GC Content vs Exon Length by Origin'), xlabel='Exon Length', ylabel='GC Content')
    for l in range(len(x)):
        plt.scatter(x[l], y[l], s=10)
    ax1.legend(('Bacterial','Archaeal','Eukaryotic','Universal'),title='Origins')
    plt.savefig('Total_gc_content_vs_exon_length.png')
    plt.close('all')

def save_proteins_of_interest_by_species(species, protein_ids, vals, labels):
    if len(vals[0]) > 0:
        data = {'protein_id':protein_ids[0], 'gc_content': vals[0], 'origin': labels[0]}
        df = pd.DataFrame(data, columns = ['protein_id', 'gc_content', 'origin'])
        if species == 'GCF_000001405.37_GRCh38.p11':
            df.to_excel((species + '_num_exons_proteins_of_interest.xlsx'), index=False)

def save_proteins_of_interest_total(variable, protein_ids, vals, labels):
    print(variable)
    print(len(protein_ids))
    print(len(vals[0]))
    print(len(labels[0]))
    data = {'protein_id':protein_ids, 'gc_content': vals[0], 'origin': labels[0]}
    df = pd.DataFrame(data, columns = ['protein_id', 'gc_content', 'origin'])
    df.to_excel(('total_' +variable + '_proteins_of_interest.xlsx'), index=False)

intron_and_exon_files = listdir('processed_data/intron_and_exon_excel_files')
intron_files = listdir('processed_data/intron_excel_files')
exon_files = listdir('processed_data/exon_excel_files')

for f in range(len(exon_files)):
    print('_'.join(exon_files[f].split('_')[:-1]))
    
    if exon_files[f].endswith('.csv'):
        exons = pd.read_csv('processed_data/exon_excel_files/' + exon_files[f])
    else:
        exons = pd.read_excel('processed_data/exon_excel_files/' + exon_files[f])
    if f == 0:
        total_bacterial, total_archaeal, total_eukaryotic, total_universal, total_ids = make_hist_by_origin('_'.join(exon_files[f].split('_')[:-1]),'gc_content', exons)
        length_bacterial, length_archaeal, length_eukaryotic, length_universal, lengths_ids = make_hist_by_origin('_'.join(exon_files[f].split('_')[:-1]),'exon_lengths', exons)
        make_comparison('_'.join(exon_files[f].split('_')[:-1]), [length_bacterial, length_archaeal, length_eukaryotic, length_universal], [total_bacterial, total_archaeal, total_eukaryotic, total_universal])
        num_exon_labels = [(['bacterial'] * len(length_bacterial)) + (['archaeal'] * len(length_archaeal)) + (['eukaryotic'] * len(length_eukaryotic)) + (['universal'] * len(length_universal))]
        total_labels = [ (['bacterial'] * len(total_bacterial)) + (['archaeal'] * len(total_archaeal)) + (['eukaryotic'] * len(total_eukaryotic)) + (['universal'] * len(total_universal))]
        save_proteins_of_interest_by_species('_'.join(exon_files[f].split('_')[:-1]), total_ids, [total_bacterial + total_archaeal + total_eukaryotic + total_universal], total_labels)
        #save_proteins_of_interest_by_species('_'.join(intron_files[f].split('_')[:-1]), total_ids, [total_bacterial + total_archaeal + total_eukaryotic + total_universal], total_labels)
    else:
        temp_b, temp_a, temp_e, temp_u, temp_ids = make_hist_by_origin('_'.join(exon_files[f].split('_')[:-1]),'gc_content', exons)
        length_b, length_a, length_e, length_u, length_ids = make_hist_by_origin('_'.join(exon_files[f].split('_')[:-1]),'exon_lengths', exons)
        #num_introns_bacterial, num_introns_archaeal, num_introns_eukaryotic, num_introns_universal = make_hist_by_origin('_'.join(intron_files[f].split('_')[:-1]),'num_introns', introns_and_exons)
        make_comparison('_'.join(exon_files[f].split('_')[:-1]), [length_b, length_a, length_e, length_u], [temp_b, temp_a, temp_e, temp_u])
        labels = [ (['bacterial'] * len(temp_b)) + (['archaeal'] * len(temp_a)) + (['eukaryotic'] * len(temp_e)) + (['universal'] * len(temp_u))]
        num_exon_labels = [ (['bacterial'] * len(length_b)) + (['archaeal'] * len(length_a)) + (['eukaryotic'] * len(length_e)) + (['universal'] * len(length_u))]
        save_proteins_of_interest_by_species('_'.join(exon_files[f].split('_')[:-1]), total_ids, [temp_b + temp_a + temp_e + temp_u], labels)
        #save_proteins_of_interest_by_species('_'.join(exon_files[f].split('_')[:-1]), gc_temp_ids, [gc_temp_b + gc_temp_a + gc_temp_e + gc_temp_u], gc_labels)
        length_bacterial = length_bacterial + length_b
        length_archaeal = length_archaeal + length_a
        length_eukaryotic = length_eukaryotic + length_e
        length_universal = length_universal + length_u
        total_bacterial = total_bacterial + temp_b
        total_archaeal = total_archaeal + temp_a
        total_eukaryotic = total_eukaryotic + temp_e
        total_universal = total_universal + temp_u
        total_ids = total_ids + temp_ids
        lengths_ids += length_ids

all_data = total_eukaryotic+total_archaeal+total_bacterial+total_universal
total_labels = [ (['bacterial'] * len(total_bacterial)) + (['archaeal'] * len(total_archaeal)) + (['eukaryotic'] * len(total_eukaryotic)) + (['universal'] * len(total_universal))]
save_proteins_of_interest_total('gc_content', total_ids, [total_bacterial + total_archaeal + total_eukaryotic + total_universal], total_labels)
make_total_comparison([length_bacterial, length_archaeal, length_eukaryotic, length_universal], [total_bacterial, total_archaeal, total_eukaryotic, total_universal])
make_hist_total_by_origin(total_archaeal, total_bacterial, total_eukaryotic, total_universal, 'gc_content')
make_hist_total(all_data, 'gc_content')

    