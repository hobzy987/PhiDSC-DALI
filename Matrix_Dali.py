from Biomuta_data_retrieve import temp_biomuta_list2
from DALI_Align import Matrix_of_residue_align
from DALI_Align import removed_gaps_position
from Biomuta_data_retrieve import hotspot_chang
import numpy as np


score_matrix = [[] for _ in range(len(removed_gaps_position))]
for i in range(len(removed_gaps_position)):
    for j in range(len(removed_gaps_position[i])):
        if removed_gaps_position[i][j] in temp_biomuta_list2[i]:
            score_matrix[i].append(1)
        if removed_gaps_position[i][j] not in temp_biomuta_list2[i]:
            score_matrix[i].append(0)
        if removed_gaps_position[i][j] in hotspot_chang[i]:
            score_matrix[i][j] = 2
# print(score_matrix)
scorearray = np.array(score_matrix)
score_array = np.transpose(scorearray)
# print(Number_Matrix4)
scores = np.count_nonzero(score_array, axis=1)
Hotspot_Score = np.sum(score_array, axis=1)
# scores[:] = [x - 1 for x in scores]
# print(scores)
score_of_residue = np.concatenate((Matrix_of_residue_align, scores[:, None], Hotspot_Score[:, None]), axis=1)
# np.savetxt('score_of_residue.csv', score_of_residue, delimiter=',', fmt='%s')

print(score_of_residue)
