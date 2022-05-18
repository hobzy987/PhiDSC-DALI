from probability_scores import Probability_Score1
from probability_scores import sorted_Probability_Sores, filter_scores
from Matrix_Dali import scores
from Matrix_Dali import Matrix_of_residue_align as Matrix42, Hotspot_Score
from DALI_Align import Dali_Input as unq_Ston_PDB1
from DALI_Align import len_input_protein
from Filter_PDB import input_protein_list2
from Biomuta_data_retrieve import HGNC_Fam_List, mutation_position_input
import numpy as np
import pandas as pd
import os


def Sort(sub_li):
    l = len(sub_li)
    for i in range(0, l):
        for j in range(0, l-i-1):
            if (sub_li[j][1] > sub_li[j + 1][1]):
                tempo = sub_li[j]
                sub_li[j]= sub_li[j + 1]
                sub_li[j + 1]= tempo
    return sub_li



d = len(unq_Ston_PDB1)

Threshold1 = "{:.2e}".format(0.05/d)
Threshold = 0.01/d

t =[]
Sorted_Probability_Scores = []
for i in range(len(sorted_Probability_Sores)):
    t.append(list(sorted_Probability_Sores[i]))
for i in range(len(t)):
    t[i][0] = int(t[i][0]) + 1
    t[i][1] = int(t[i][1])
    Sorted_Probability_Scores.append(tuple(t[i]))
# filter_scores = dict(filter_scores)
sorted_Probability_Sores1 = []
for i in range(len(filter_scores)):
    if filter_scores[i][2] <= Threshold and filter_scores[i][1] > (d//2):
        w = [filter_scores[i][0], filter_scores[i][2]]
        # if w[0] not in mutation_position_input:
        sorted_Probability_Sores1.append(w)
sorted_Probability_Sores1 = Sort(sorted_Probability_Sores1)
Final_Result_Matrix = np.append(Matrix42, scores[:, None], 1)
Final_Result_Matrix = np.append(Final_Result_Matrix, Hotspot_Score[:, None], 1)
Final_Result_Matrix = np.append(Final_Result_Matrix, np.array(Probability_Score1)[:, None], 1)
# print(Final_Result_Matrix)

# np.savetxt('results.csv', Final_Result_Matrix, delimiter=',', fmt='%s')

Index = []
Numbers = range(1, len_input_protein+1)
for count in Numbers:
    Index.append(count)


columns = HGNC_Fam_List + ['Scores', 'Hotspot Score', 'Probability scores']


df = pd.DataFrame(Final_Result_Matrix, index=Index, columns=columns)

html = df.to_html()




def wrapStringInHTMLMac(program, url, body, body1, body2):
    import datetime
    from webbrowser import open_new_tab

    now = datetime.datetime.today().strftime("%Y%m%d-%H%M%S")
    filename = program + '.html'
    f = open(filename, 'w')

    wrapper = """<html>
    <head>
    <title>%s output - %s</title>
    </head>
    <body><p>InputPDB: <a href=\"%s\">%s</a></p><p>%s</p></body>
    </head>
    <body1><p>Sorted_Probability(Pos,Pro);Threshold : <a href=\"%s\">%s</a></p><p>%s</p></body1>
     </head>
    <body2><p>Biomuta_mutations;Threshold : <a href=\"%s\">%s</a></p><p>%s</p></body2>
    </html>"""

    whole = wrapper % (program, now, program, program, body, url, url, body1, url, url, body2)
    f.write(whole)
    f.close()

    filename = os. getcwd() + '/' + filename

    open_new_tab(filename)


wrapStringInHTMLMac(input_protein_list2[0], url="{:.2e}".format(Threshold), body=html, body1=sorted_Probability_Sores1, body2=mutation_position_input)
# print(html)
