import gzip
import pandas as pd
from os import listdir
from Bio.Seq import Seq
from Bio.SeqUtils import GC

# requires path of the feature table
# returns start index, end index, genomic accession (used to id correct unfiltered genome), protein id (used to id correct filtered protein sequence)
def get_features(path): 
    starts = []
    ends = []
    protein_ids = []
    genomic_accession = []
    strand = []
    with open(path) as f:
        while True:
            line = (f.readline())
            if len(line) == 0:
                break
            if line.startswith('CDS'):
                vals = " ".join(line.split())
                split_vals = vals.split()
                if all(x in split_vals for x in ['-','+']):
                    if split_vals.index('+') < split_vals.index('-'):
                        strand_idx = split_vals.index('+')
                    else:
                        strand_idx = split_vals.index('-')
                elif '-' in split_vals:
                    strand_idx = split_vals.index('-')
                else:
                    strand_idx = split_vals.index('+')
                
                starts.append(int(split_vals[strand_idx-2]))
                ends.append(int(split_vals[strand_idx-1]))
                genomic_accession.append(split_vals[strand_idx-3])
                strand.append(split_vals[strand_idx])
                protein_ids.append(split_vals[strand_idx+1])

    return starts, ends, protein_ids, genomic_accession, strand

# returns file number for filtered protein sequence
def get_file_num(access_number):
    with open('/Data/Databases_derived/Eukaryotes/January_2018/genomes/eukaryote_list_v1.txt') as f:
        for num,line in enumerate(f):
            if access_number in line:
                temp = line.split('file_')[1]
                return temp.split('\t')[0]

# depending on whether input needs the protein sequence or the genome, either can be returned given the right input
def get_filtered_seq(file_number, id_number): 
    seq = []
    mark = False
    if 'Exon' in file_number:
        f = open(file_number + "_translated_cds.faa")
    elif 'GCF' in file_number:
        g = gzip.open('/Data/Databases_raw/Eukaryotes/January_2018/genomes/' + file_number + '_genomic_refseq.fna.gz', mode='r') #file_number is GCA or GCF identifier
        f = str(g.read()).split('\\n')
    else:
        f = open('/Data/Databases_derived/Eukaryotes/January_2018/genomes/sequences/file_' + file_number) #file_number is an integer between 1 and 151
    for num,line in enumerate(f):
        if id_number in line and 'location=' in line:
            mark = True
        if mark == True:
            #if len(seq)>0:
            if line[0] == '>' and line[1:4] != 'lcl':
                break
            else:
                seq.append(line.strip())
                if 'gbkey' in line:
                    break
        if id_number in str(line):
            mark = True
    return seq

def get_start_end_idxs(frames, exons):
    starts = []
    ends = []
    for i in range(len(exons)):
        if exons[i] in frames[0]:
            starts.append(frames[0].find(exons[i]))
            ends.append(starts[i] + len(exons[i]))
        elif exons[i] in frames[1]:
            starts.append(frames[1].find(exons[i]))
            ends.append(starts[i] + len(exons[i]))
        else:
            starts.append(frames[2].find(exons[i]))
            ends.append(starts[i] + len(exons[i]))
        starts = sorted(starts)
        ends = sorted(ends)
    return starts, ends

def save_results_by_species(table_name, protein_ids, exon_nums, prop_exon, intron_nums, prop_intron):
    data = {'protein_id':protein_ids, 'num_exons': exon_nums, 'proportion_exon':prop_exon, 'num_introns': intron_nums, 'proportion_intron':prop_intron}
    df = pd.DataFrame(data, columns = ['protein_id', 'num_exons', 'proportion_exon', 'num_introns', 'proportion_intron'])
    df.to_excel((table_name + '_introns_and_exons.xlsx'), index=False)

def inspect_introns(gen_access_ids, table_name, protein_ids, intron_pos, intron_lens, gc_content):
    data = {'genome_id':gen_access_ids, 'protein_id':protein_ids, 'intron_position': intron_pos, 'intron_lengths': intron_lens, 'gc_content': gc_content}
    df = pd.DataFrame(data, columns = ['genome_id', 'protein_id', 'intron_position', 'intron_lengths', 'gc_content'])
    if len(df.index) > 1048576:
        df.to_csv((table_name + '_introns.csv'), index=False)
    else:
        df.to_excel((table_name + '_introns.xlsx'), index=False)

def inspect_exons(gen_access_ids, table_name, protein_ids, exon_pos, exon_lens, gc_content):
    data = {'genome_id':gen_access_ids, 'protein_id':protein_ids, 'exon_position': exon_pos, 'exon_lengths': exon_lens, 'gc_content': gc_content}
    df = pd.DataFrame(data, columns = ['genome_id', 'protein_id', 'exon_position', 'exon_lengths', 'gc_content'])
    if len(df.index) > 1048576:
        df.to_csv((table_name + '_exons.csv'), index=False)
    else:
        df.to_excel((table_name + '_exons.xlsx'), index=False)

def get_introns_and_exons(path, protein_id):
    temp = get_filtered_seq(path, protein_id)
    if len(temp) > 0:
        if not temp[0].startswith('>lcl'):
            splits = temp[-1].split('[')
        else:
            splits = temp[0].split('[')
        exon_locations = splits[-2].split('(')[-1][:-2] #list of exon positions, if there's still a parentheses at the end it means it is on the complement strand, if starts with location it is only 1 exon no introns
        return exon_locations
    else:
        print(temp)
        print(protein_id)
        return -1

GCA_table_names = listdir('8_of_9_GCAs_feature_tables')[1:-1]
names = listdir('94_of_121_GCFs_feature_tables')
GCF_table_names = []
for i in range(len(names)):
    if 'GCF' in names[i]:
        GCF_table_names.append(names[i])
table_names = GCA_table_names + GCF_table_names
species_names = []  
strands = []
exon_table_name = []

for i in range(54,len(GCF_table_names)):
    exon_nums = []
    intron_nums = []
    exon_positions = [] 
    intron_positions = [] 
    protein_id_nums_exon_rep = [] 
    protein_id_nums_intron_rep = []
    gen_access_ids_exon_rep = []
    gen_access_ids_intron_rep = []
    exon_lengths = [] 
    intron_lengths = [] 
    exon_gcs = [] 
    intron_gcs = []
    exon_proportions = []
    intron_proportions = []
    start_vals, end_vals, protein_id_nums, gen_access_ids, strands = get_features('94_of_121_GCFs_feature_tables/' + GCF_table_names[i])
    exon_table_name.append("_".join(GCF_table_names[i].split('_')[:-2]))
    species_names.append('_'.join(GCF_table_names[i].split('_')[2:-2]))
    print(GCF_table_names[i])
    for j in range(len(protein_id_nums)):
        #go through each of the files
        protein_id_nums_edit = protein_id_nums
        exon_idxs = get_introns_and_exons('Exon_Positions/' + exon_table_name[i-54], protein_id_nums[j])
        if exon_idxs == -1:
            print(len(protein_id_nums))
            print(start_vals[j])
            print(end_vals[j])
            print(gen_access_ids[j])
            protein_id_nums = protein_id_nums[:j] + protein_id_nums[j+1:]
            start_vals = start_vals[:j] + start_vals[j+1:]
            end_vals = end_vals[:j] + end_vals[j+1:]
            gen_access_ids = gen_access_ids[:j] + gen_access_ids[j+1:]
            strands = strands[:j] + strands[j+1:]
            print(len(protein_id_nums))
            exon_idxs = get_introns_and_exons('Exon_Positions/' + exon_table_name[i-54], protein_id_nums[j])
        #getting the genome
        if (j != 0) & (gen_access_ids[j] != gen_access_ids[j-1]):
            gg = ''.join(get_filtered_seq("_".join(GCF_table_names[i].split('_')[:-2]), gen_access_ids[j])) 
        elif j == 0:
            gg = ''.join(get_filtered_seq("_".join(GCF_table_names[i].split('_')[:-2]), gen_access_ids[j]))
        if j % 10000 == 0:
            print(gen_access_ids[j])
        exon_list = []
        intron_list = []
        exoncount1= len(exon_lengths)
        introncount1 = len(intron_lengths)
        elist = exon_idxs.replace('<','').replace(')','').replace('>','').replace('location=','').split(',')
        for k in range(len(elist)):
            exon_list.append(elist[k].split('..'))
        if len(elist) == 1:
            exon_list[0][0] = int(exon_list[0][0])
            exon_list[0][1] = int(exon_list[0][1])
        for k in range(0,len(exon_list)-1):
            if len(exon_list[k]) == 1:
                exon_list[k].append(exon_list[k][0])
            if len(exon_list[-1]) == 1:
                exon_list[-1].append(exon_list[-1][0])
            if strands[j] == '-':
                intron_list.append([int(exon_list[k][1])+1,int(exon_list[k+1][0])])
                intron_gcs.append(GC(gg[int(exon_list[k][1])+1:int(exon_list[k+1][0])]))
                exon_list[k][0] = int(exon_list[k][0])
                exon_list[k][1] = int(exon_list[k][1])+1
                exon_list[-1][0] = int(exon_list[-1][0])
                exon_list[-1][1] = int(exon_list[-1][1])+1
            else:
                intron_list.append([int(exon_list[k][1])+1,int(exon_list[k+1][0])+1])
                intron_gcs.append(GC(gg[int(exon_list[k][1])+1:int(exon_list[k+1][0])+1]))
                exon_list[k][0] = int(exon_list[k][0])-1
                exon_list[k][1] = int(exon_list[k][1])
                exon_list[-1][0] = int(exon_list[-1][0])-1
                exon_list[-1][1] = int(exon_list[-1][1])
            intron_positions.append([int(exon_list[k][1])+1,int(exon_list[k+1][0])-1])
            intron_lengths.append((int(exon_list[k+1][0])-1)-(int(exon_list[k][1])+1))
            protein_id_nums_intron_rep.append(protein_id_nums[j])
            gen_access_ids_intron_rep.append(gen_access_ids[j])
        exon_nucleotides = []
        exon_proteins = []
        for l in range(len(exon_list)):
            exon_positions.append(exon_list[l])
            exon_lengths.append(int(exon_list[l][1])-int(exon_list[l][0]))
            protein_id_nums_exon_rep.append(protein_id_nums[j])
            if strands[j] == '-':
                seq = Seq(gg[exon_list[l][0]:exon_list[l][1]])
                exon_nucleotides.append(seq.reverse_complement())
                exon_proteins.append(exon_nucleotides[l].translate())
            else:
                exon_nucleotides.append(Seq(gg[exon_list[l][0]:exon_list[l][1]]))
                exon_proteins.append(exon_nucleotides[l].translate())
            exon_gcs.append(GC(gg[exon_list[l][0]:exon_list[l][1]]))
            gen_access_ids_exon_rep.append(gen_access_ids[j])
            
            introncount2 = len(intron_lengths)
            exoncount2 = len(exon_lengths)

        exon_proportions.append(sum(exon_lengths[exoncount1:exoncount2])/(sum(exon_lengths[exoncount1:exoncount2])+sum(intron_lengths[introncount1:introncount2])))
        intron_proportions.append(sum(intron_lengths[introncount1:introncount2])/(sum(exon_lengths[exoncount1:exoncount2])+sum(intron_lengths[introncount1:introncount2])))

        exon_nums.append(len(exon_proteins))
        intron_nums.append(len(intron_list))

    species_names.append('_'.join(table_names[2].split('_')[2:-2]))

    print(len(protein_id_nums))
    print(len(exon_nums))
    print(len(exon_proportions))
    print(len(intron_nums))
    print(len(intron_proportions))
    #save_results_by_species(exon_table_name[i-54], protein_id_nums, exon_nums, exon_proportions, intron_nums, intron_proportions)

    inspect_exons(gen_access_ids_exon_rep, exon_table_name[i-54], protein_id_nums_exon_rep, exon_positions, exon_lengths, exon_gcs)
    inspect_introns(gen_access_ids_intron_rep, exon_table_name[i-54], protein_id_nums_intron_rep, intron_positions, intron_lengths, intron_gcs)
