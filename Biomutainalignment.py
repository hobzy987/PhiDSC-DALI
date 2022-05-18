import numpy as np
import operator
from align_matrix import residue_fam_align1
from Biomuta_data_retrieve import temp_biomuta_list3
from align_matrix import input_align_list as Input_Protein_seq
from align_matrix import residue_input_align



def String_character_counter(string, start, end):
    string = str(string)
    start = int(start)
    end = int(end)
    couter_list = []
    if start != end:
        for char in string:
            if char.isalpha() is True or char == '_':
                couter_list.append(start)
                start = start + 1
            if char == '.':
                couter_list.append(0)

    return couter_list


for i in range(len(residue_fam_align1)):
    for j in range(len(residue_fam_align1[i])):
        if residue_fam_align1[i][j] not in temp_biomuta_list3[i]:
            residue_fam_align1[i][j] = 0
        else:
            residue_fam_align1[i][j] = 1
# print(residue_fam_align1)

# Number_Matrix2 = []
len_Input_Protein = len(residue_input_align[0])
Input_Seq_Num = String_character_counter(Input_Protein_seq[0], 1, len_Input_Protein)
# for i in range(len(Number_Matrix1)):
#     a = [0] * (add_preindexes[i] - 1)
#     b = [0] * ((len_Input_Protein[0]) - add_postindexs[i])
#     Number_Matrix2.append(a + Number_Matrix1[i] + b)
# Number_Matrix3 = Number_Matrix2
Number_Matrix = [Input_Seq_Num] + residue_fam_align1
# # Number_Matrix3 = []
# # for i in range(len(Number_Matrix)):
Number_Matrix4 = np.array(Number_Matrix)
Number_Matrix4 = np.transpose(Number_Matrix)
# print(Number_Matrix4)
scores = np.count_nonzero(Number_Matrix4, axis=1)
scores[:] = [x - 1 for x in scores]
#
#
sorted_scores = sorted(enumerate(scores), key=operator.itemgetter(1), reverse=True)
# print(sorted_scores)