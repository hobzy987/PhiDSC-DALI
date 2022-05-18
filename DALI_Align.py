import Filter_PDB
from Filter_PDB import List_Wild_PDB1, input_protein_list2
import urllib3
import os
import sys
import numpy as np
if sys.version_info[0] == 2:
    from urllib import urlretrieve
else:
    from urllib.request import urlretrieve
urllib3.disable_warnings()

# pdb file download (extentin .ent.gz)
###############################################################################################

PDB_Protein_Input = Filter_PDB.input_protein_list2[0]
unq_Ston_PDB = Filter_PDB.Wild_PDB


if sys.version_info[0] == 2:
    from urllib import urlretrieve
else:
    from urllib.request import urlretrieve

download_folder = 'PDBIDs'
compressed = False

try:
    os.makedirs(download_folder)
except OSError as e:
    pass
unq_Ston_PDB1 = []
for i in range(len(unq_Ston_PDB)):
    if type(unq_Ston_PDB[i]) == list:
        d = unq_Ston_PDB[i][0]
        unq_Ston_PDB1.append(d)
    else:
        d = unq_Ston_PDB[i]
        unq_Ston_PDB1.append(d)
for pdb_code in unq_Ston_PDB1:
    # for pdb_code in unq_Ston_PDB:
    filename = '%s.ent.gz' % pdb_code
    if compressed:
        filename = '%.gz' % filename
    url = 'https://files.rcsb.org/pub/pdb/data/structures/all/pdb/pdb%s.ent.gz' % pdb_code.lower()
    destination_file = os.path.join(download_folder, filename)

    urlretrieve(url, destination_file)

url = 'https://files.rcsb.org/pub/pdb/data/structures/all/pdb/pdb%s.ent.gz' % PDB_Protein_Input.lower()
filename = '%s.ent.gz' % PDB_Protein_Input
destination_file = os.path.join(download_folder, filename)
urlretrieve(url, destination_file)

# print(hotspot_chang)
##############################################################################################

# Dali import file .DAT
################################################################################################
os.chdir('../PhiDSC-DALI')
DAT_file_Dali = unq_Ston_PDB1
DAT_file_Dali.append(PDB_Protein_Input)
import subprocess, sys
for pdb_dat in DAT_file_Dali:
    imdali = "./DaliLite.v5/bin/import.pl --pdbfile ./PDBIDs/%s.ent.gz --pdbid %s --dat ./DAT" % (pdb_dat, pdb_dat)
    ## run it ##
    p = subprocess.Popen(imdali, shell=True, stderr=subprocess.PIPE)
    p.communicate()

########################################################################################

#Align with Dali lite v5
##########################################################################################
target_list = []
for i in range(len(List_Wild_PDB1)):
    # for j in List_Wild_PDB1[i]:
    chain_comb = List_Wild_PDB1[i][0] + List_Wild_PDB1[i][1]
    target_list.append(chain_comb)
input_protein_list3 = input_protein_list2[0] + input_protein_list2[1]
for i in range(len(target_list)):
    var1 = './DaliLite.v5/bin/dali.pl --cd1 %s --cd2 %s --dat1 ./DAT --dat2 ./DAT --title "output options" --outfmt ' \
           '"summary,alignments,equivalences,transrot" --clean 2> err' % ( target_list[i], input_protein_list3)
    p = subprocess.Popen(var1, shell=True, stderr=subprocess.PIPE)
    p.communicate()

###############################################################################################
def reemovNestings(l):
    output = []
    for i in range(len(l)):
        if len(l[i]) == 1:
            output.append(l[i][0])
        else:
            output.append(l[i])
    return output


Input_Protein = str(input_protein_list2[0].lower() + input_protein_list2[1])
List_Wild_PDB2 = reemovNestings(List_Wild_PDB1)
Dali_Input = []
for i in range(len(List_Wild_PDB1)):
    Dali_Input.append(str(List_Wild_PDB2[i][0] + List_Wild_PDB2[i][1]))

alignment = [[[] for _ in range(2)] for _ in range(len(List_Wild_PDB1))]
end_alignment = [[[] for _ in range(2)] for _ in range(len(List_Wild_PDB1))]
for k in range(len(Dali_Input)):
    with open('%s.txt' % Dali_Input[k]) as f:
        s = f.read()
        w = s.splitlines()
        for j in range(len(w)):
            if len(w[j]) != 0:
                if w[j][0] == 'Q':
                    alignment[k][1].append(w[j].split()[1])
                    end_alignment[k][1].append(w[j].split()[2])
                if w[j][0] == 'S':
                    alignment[k][0].append(w[j].split()[1])
                    end_alignment[k][0].append(w[j].split()[2])
end_alignment1 = [[[] for _ in range(2)] for _ in range(len(end_alignment))]
for k in range(len(end_alignment)):
    if len(end_alignment[k][0]) != 0:
        end_alignment1[k][0].append(end_alignment[k][0][-1])
        end_alignment1[k][1].append(end_alignment[k][1][-1])
    if len(end_alignment[k][0]) == 0:
        del Dali_Input[k]
        continue

alignment1 = [[[] for _ in range(2)] for _ in range(len(alignment))]
for i in range(len(alignment)):
    if len(alignment[i][0]) != 0:
        alignment1[i][0].append("".join(alignment[i][0]))
        alignment1[i][1].append("".join(alignment[i][1]))
    else:
        continue

# print(alignment1)
for i in range(len(alignment1)):
    # for j in range(len(alignment1)):
    alignment1[i] = [x for x in alignment1[i] if x]
    end_alignment1[i] = [x for x in end_alignment1[i] if x]
alignment1 = [x for x in alignment1 if x]
end_alignment1 = [x for x in end_alignment1 if x]


def Reverse_counter(string, end_number):
    counter_list = []
    string = str(string)
    end_number = int(end_number)
    for c in reversed(string):
        if c.isalpha():
            counter_list.append(end_number)
            end_number = end_number - 1
        else:
            counter_list.append(0)
    return counter_list[::-1]


align_matrix_gapped = [[[] for _ in range(2)] for _ in range(len(Dali_Input))]
for i in range(len(alignment1)):
    for j in range(len(alignment1[i])):
        s = Reverse_counter(alignment1[i][j][0], end_alignment1[i][j][0])
        align_matrix_gapped[i][j].append(s)


# print(align_matrix_gapped)


# get_indexes = lambda x, xs: [i for (y, i) in zip(xs, range(len(xs))) if x == y]
def find(q, ch):
    return [i for i, ltr in enumerate(q) if ltr == ch]


def removal_index(string, list_of_index):
    list_of_index = list(sorted(list_of_index))
    sep_str = []

    if max(list_of_index) < len(string):
        sep_str.append(string[0:list_of_index[0]])
        for z in range(len(list_of_index)):
            z = z - 1
            sep_str.append(string[list_of_index[z] + 1:list_of_index[z + 1]])
        sep_str.append(string[list_of_index[-1] + 1:])
        list_of_index.pop(-1)
        result = ''.join(sep_str)
        return result


index_list = []
for i in range(len(alignment1)):
    for j in range(len(alignment1[i][0])):
        w = find(alignment1[i][0][0], '-')
        index_list.append(w)
# print(index_list)

removed_gaps_residue = []
for i in range(len(alignment1)):
    if len(index_list[i]) != 0:
        x = removal_index(alignment1[i][1][0], index_list[i])
        removed_gaps_residue.append(x)
    else:
        removed_gaps_residue.append(alignment1[i][1][0])
removed_gaps_residue1 = []
for i in range(len(alignment1)):
    if len(index_list[i]) != 0:
        x = removal_index(alignment1[i][0][0], index_list[i])
        removed_gaps_residue1.append(x)
    else:
        removed_gaps_residue1.append(alignment1[i][0][0])


# print(removed_gaps_residue, removed_gaps_residue1)
def split(word):
    return [char for char in word]


for i in range(len(removed_gaps_residue)):
    removed_gaps_residue[i] = split(removed_gaps_residue[i])
Matrix_of_residue_align = np.asanyarray(removed_gaps_residue)
Matrix_of_residue_align = np.transpose(Matrix_of_residue_align)
# print(Matrix_of_residue_align)
removed_gaps_position = []
for i in range(len(align_matrix_gapped)):
    removed_gaps_position.append(align_matrix_gapped[i][1][0])
for k in range(len(index_list)):
    if len(index_list[k]) != 0:
        for ele in sorted(index_list[k], reverse=True):
            del removed_gaps_position[k][ele]
    else:
        continue
Matrix_of_position_align = np.asanyarray(removed_gaps_position)
Matrix_of_position_align = np.transpose(Matrix_of_position_align)

len_input_protein = len(removed_gaps_position[0])
print(len_input_protein)
