import pandas as pd
from os import listdir

tree_names = listdir('data/trees_with_dup')
supp = []
trees_origin = []
bacterial_proteins = []
archaeal_proteins = []
universal_proteins = []
eukaryotic_proteins = []
num_replications = []

with open('files/SuppD1.txt') as s:
    while True:
        line = s.readline().split('\t')[:-1]
        if len(line) == 0:
             break
        if line[2] != '0' and line[2] != 'No of LECA dupl.':
            num_replications.append(line[2])
            supp.append(line)
supp_transpose = pd.DataFrame(supp).T.values.tolist()

for t in range(len(tree_names)):
    #print(tree_names[t])
    species_and_proteins = []
    species_ids = []
    protein_ids = []
    leca_dups = []
    
    trees_origin.append([tree_names[t],supp_transpose[1][supp_transpose[0].index(tree_names[t])]])
    with open('data/trees_with_dup/' + tree_names[t]) as f:
        data = f.read().replace('(','').replace(')','').split(',')
        for i in range(len(data)):
            species_and_proteins.append(data[i].split(':')[0])
        #print(species_and_proteins)
    tree_ids = pd.read_excel('trees.xlsx').T.values.tolist()
    #print(tree_ids)

    for j in range(len(species_and_proteins)):
        split = species_and_proteins[j].split('_',1)
        if split[0] in tree_ids[0] and split[0] not in species_ids: 
            species_ids.append(tree_ids[1][tree_ids[0].index(split[0])])
            protein_ids.append(split[1])
            if trees_origin[-1][1] == 'Bacterial':
                bacterial_proteins.append(split[1])
            elif trees_origin[-1][1] == 'Eukaryotic':
                eukaryotic_proteins.append(split[1])
            elif trees_origin[-1][1] == 'Archaeal':
                archaeal_proteins.append(split[1])
            else:
                universal_proteins.append(split[1])
# print(len(trees_origin))
# print(len(num_replications))
# print(len(bacterial_proteins))
# print(len(eukaryotic_proteins))
# print(len(archaeal_proteins))
# print(len(universal_proteins))
# print(len(species_and_proteins))
