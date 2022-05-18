import os
# import re
# import numpy as np
import numpy as np
from grap_alignment import filename


get_indexes = lambda x, xs: [i for (y, i) in zip(xs, range(len(xs))) if x == y]


def removal_index(string, list_of_index):
    list_of_index = list(sorted(list_of_index))
    sep_str = []

    if max(list_of_index) < len(string):
        sep_str.append(string[0:list_of_index[0]])
        for z in range(len(list_of_index)):
            z = z-1
            sep_str.append(string[list_of_index[z]+1:list_of_index[z+1]])
        sep_str.append(string[list_of_index[-1]+1:])
        list_of_index.pop(-1)
        result = ''.join(sep_str)
        return result



os.chdir('/Users/mohamadhoballa/PycharmProjects/tmaligners/PDBIDs')
path = '/Users/mohamadhoballa/PycharmProjects/tmaligners/PDBIDs'
# filename = "Alignment-4Q21.txt"
file_name = os.path.join(path, filename)
my_file = open(filename, 'r')
similar: str = my_file.read()
with open(filename, 'r') as file:
    x = file.readlines()
align_matrix = []
for i in range(len(x)):
    x[i] = x[i].strip()
    align_matrix.append(x[i])
align_matrix = list(filter(None, align_matrix))
for x in align_matrix:
    if ':' in x:
        align_matrix.pop(align_matrix.index(x))
# align_matrix.pop()
# print(align_matrix)
input_align_list = []
fam_align_list = []
for i in range(len(align_matrix)):
    if i % 2 == 0:
        input_align_list.append(align_matrix[i])
    else:
        fam_align_list.append(align_matrix[i])


def String_character_counter(string, start):
    string = str(string)
    start = int(start)
    # end = int(end)
    counter_list = []
    # if start != end:
    for char in string:
        if char.isalpha() is True:
            counter_list.append(start)
            start = start + 1
        if not char.isalpha():
            counter_list.append(0)

    return counter_list

Residue_Align_Matrix = []
for i in range(len(align_matrix)):
    w = String_character_counter(align_matrix[i], 1)
    Residue_Align_Matrix.append(w)
residue_input_align = []
residue_fam_align = []
for i in range(len(align_matrix)):
    if i % 2 == 0:
        residue_input_align.append(Residue_Align_Matrix[i])
    else:
        residue_fam_align.append(Residue_Align_Matrix[i])
# print(Residue_Align_Matrix)

index_list = []
for i in range(len(Residue_Align_Matrix)):
    if i % 2 == 0:
        s = get_indexes(0, Residue_Align_Matrix[i])
        index_list.append(s)
fam_align_list1 = []
fam_align_list1.append(input_align_list[0])
for i in range(len(index_list)):
    if len(index_list[i]) != 0:
        d = removal_index(fam_align_list[i], index_list[i])
        fam_align_list1.append(d)
    else:
        fam_align_list1.append(fam_align_list[i])
# print(fam_align_list1)
residue_fam_align1 = []
for i in range(len(index_list)):
    n = [y for x, y in enumerate(residue_fam_align[i]) if x not in index_list[i]]
    residue_fam_align1.append(n)
# print(residue_fam_align1)
end_protein = []
for i in range(len(residue_fam_align1)):
    end_protein.append(max(residue_fam_align1[i]))
# print(end_protein)
def split(word):
    return [char for char in word]


for i in range(len(fam_align_list)):
    fam_align_list[i] = split(fam_align_list[i])
Matrix42 = np.asanyarray(fam_align_list)
Matrix42 = np.transpose(Matrix42)
# print(Matrix42)
len_input_protein = removal_index(input_align_list[0], index_list[0])
