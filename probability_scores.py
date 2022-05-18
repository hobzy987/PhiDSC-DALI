from DALI_Align import removed_gaps_residue
from Biomuta_data_retrieve import temp_biomuta_list2, len_UniProt_IDs, hotspot_chang
from Matrix_Dali import score_matrix, Hotspot_Score
import numpy as np
import operator

# print(Matrix41)
# print(temp_biomuta_list2)
# print(temp_biomuta_list3)
# print(len(Number_Matrix3))
Mutation_Probability = []
for i in range(len(temp_biomuta_list2)):
    a = len(temp_biomuta_list2[i])
    b = len_UniProt_IDs[i][0]
    x = a / b
    Mutation_Probability.append(x)
Hotspot_Probability = []
for i in range(len(hotspot_chang)):
    if len(hotspot_chang[i]) != 0:
        b = len_UniProt_IDs[i][0]
        c = len(hotspot_chang[i])
        x = c / b
        Hotspot_Probability.append(x)
    else:
        Hotspot_Probability.append(1)
# print(Mutation_Probability)
Non_Mutation_Probability = []
for i in Mutation_Probability:
    x = 1 - i
    Non_Mutation_Probability.append(x)
# print(Non_Mutation_Probability)
Matrix43 = []
for i in range(len(removed_gaps_residue)):
    Matrix43.append(''.join(removed_gaps_residue[i]))
# print(Matrix43)
Probability_Matrix = [[] for _ in range(len(Mutation_Probability))]
for i in range(len(Mutation_Probability)):
    for j in range(len(score_matrix[i])):
        if Matrix43[i][j].isalpha() is True:
            if score_matrix[i][j] == 1:
                Probability_Matrix[i].append(Mutation_Probability[i])
            if score_matrix[i][j] == 2:
                Probability_Matrix[i].append(Hotspot_Probability[i])
            if score_matrix[i][j] == 0:
                Probability_Matrix[i].append(Non_Mutation_Probability[i])
        if Matrix43[i][j].isalpha() is False:
            Probability_Matrix[i].append(1.0)
Probability_Matrix = np.asanyarray(Probability_Matrix)
Probability_Matrix = np.transpose(Probability_Matrix)
# print(Probability_Matrix)
Probability_Score = []
Probability_Score1 = []
for i in Probability_Matrix:
    x = np.prod(i)
    Probability_Score.append(x)
    Probability_Score1.append("{:.2e}".format(x))
filter_scores = []
for i in range(len(Probability_Score1)):
    s = [i+1, Hotspot_Score[i], Probability_Score[i]]
    filter_scores.append(s)

# print(len(Probability_Score))
sorted_Probability_Sores = sorted(enumerate(Probability_Score), key=operator.itemgetter(1), reverse=False)
print(sorted_Probability_Sores)
