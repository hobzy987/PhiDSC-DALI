import os
import string

protein_name1 = input('Type The Name of Your Protein: ')
protein_name = ' ' + protein_name1.upper() + '_HUMAN'

os.chdir('../PhiDSC-DALI')
path = "../PhiDSC-DALI"
file_name = os.path.join(path, "similar.txt")
my_file = open("similar.txt", 'r')
similar: str = my_file.read()
with open("similar.txt", 'r') as file:
    x = file.readlines()
for i in range(len(x)):
    if ',' not in x[i]:
        temp_family = x[i]
    if protein_name in x[i]:
        related_family = temp_family
        r = x.index(related_family)
        global protein_list
        for j in range(r + 1, len(x)):
            if ',' not in x[j]:
                break
            f = str(x[r + 1:j + 1])
            protein_list = [word.strip(string.punctuation) for word in f.split(',')]
            # print(protein_list)
human_protein_list = []
for protein in protein_list:
    if '_HUMAN' in protein:
        human_protein_list.append(protein)
Index_Input_Protein = [human_protein_list.index(i) for i in human_protein_list if protein_name in i]
Index_Input_Protein = [int(i) for i in Index_Input_Protein]
uniprot_list = []
for i in range(len(human_protein_list)):
    uniprot_list.append(human_protein_list[i][human_protein_list[i].find('(') + 1:human_protein_list[i].find(')')])
Input_Protein = human_protein_list[Index_Input_Protein[0]]
Input_Protein = Input_Protein[Input_Protein.find('(') + 1:Input_Protein.find(')')]
# print(Input_Protein)
for ID in uniprot_list:
    if ID == Input_Protein:
        uniprot_list.remove(ID)
# print(uniprot_list)
