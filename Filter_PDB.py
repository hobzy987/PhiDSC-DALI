import os
# from List_PDB_File import PDB_List, UniProt_IDs, result
# from List_PDB_File import result
from urllib.request import urlopen
from Readfamily import Input_Protein, uniprot_list
# from List_PDB_File import wild_pdb_id


def protein_length(uniprot_id):
    data = urlopen("https://rest.uniprot.org/uniprotkb/" + uniprot_id + ".txt").read()
    data1 = data.decode()
    data = data1.split('\n')
    len_Protein = [int(s) for s in data[0].split() if s.isdigit()]
    return len_Protein


len_Input_Protein = int(protein_length(Input_Protein)[0])

def remover(l):
    # new list
    newlist = []
    # loop over elements
    for i in l:
        # pdb.set_trace()
        # is element a non-empty list? then call self on it
        if isinstance(i, list) and len(i) != 0:
            newlist.append(remover(i))
        # if not a list
        if not isinstance(i, list):
            newlist.append(i)

    # return newlist
    return newlist


def reemovNestings(l):
    output = []
    for i in range(len(l)):
        if len(l[i]) == 1:
            output.append(l[i][0])
        else:
            output.append(l[i])
    return output


path = "../tmaligners"
file_name = os.path.join(path, "pdb_chain_uniprot.csv")
my_file = open("pdb_chain_uniprot.csv", 'r')
similar: str = my_file.read()
with open("pdb_chain_uniprot.csv", 'r') as file:
    x = file.readlines()
z = []
for i in range(len(x)):
    z.append(x[i].split(','))
# print(z[2][2])
pdb_chain_id = [[] for _ in range(len(uniprot_list))]
input_protein_list = []
for i in range(len(uniprot_list)):
    # for j in range(len(PDB_List[i])):
    for k in range(2, len(z)):
        if uniprot_list[i] == z[k][2]:
            f = [z[k][0].upper(), z[k][1], int(z[k][4]), z[k][2]]
            pdb_chain_id[i].extend([f])
for k in range(2, len(z)):
    if Input_Protein == z[k][2]:
        s = [z[k][0].upper(), z[k][1], int(z[k][4])]
        input_protein_list.append(s)
pdb_chain_id = remover(pdb_chain_id)
# print(input_protein_list)


UniProt_IDs = []
for i in range(len(pdb_chain_id)):
    UniProt_IDs.append(pdb_chain_id[i][0][3])
# print(UniProt_IDs)

uniprot_len_family = []
for i in range(len(UniProt_IDs)):
    s = protein_length(UniProt_IDs[i])
    uniprot_len_family.append(s[0])

os.chdir('../PhiDSC-DALI')
path = "../DALI"
file_name = os.path.join(path, "wild_pdb_id.txt")
my_file = open("wild_pdb_id.txt", 'r')
y = my_file.read()
with open("wild_pdb_id.txt", 'r') as file:
    x = file.readlines()
wild_pdb_id = []
for i in range(len(x)):
    wild_pdb_id = x[i].split(", ")

remover(wild_pdb_id)
# remove_empty_lists(wild_pdb_id)
# print(wild_pdb_id)


temp_pdb = []
List_Wild_PDB = [[] for _ in range(len(pdb_chain_id))]
for i in range(len(pdb_chain_id)):
    if len(pdb_chain_id[i]) > 1:
        for j in range(len(pdb_chain_id[i])):
            if pdb_chain_id[i][j][0] in wild_pdb_id:
                List_Wild_PDB[i].append(pdb_chain_id[i][j])
            # else:
            #     List_Wild_PDB[i].append(0)
    else:
        for j in range(len(pdb_chain_id[i])):
            List_Wild_PDB[i].append(pdb_chain_id[i][j])
for i in range(len(List_Wild_PDB)):
    if len(List_Wild_PDB[i]) == 0:
        for j in range(len(pdb_chain_id[i])):
            List_Wild_PDB[i].append(pdb_chain_id[i][j])
# List_Wild_PDB = reemovNestings(List_Wild_PDB)
# List_Wild_PDB = reemovNestings(List_Wild_PDB)

# List_Wild_PDB = remover(List_Wild_PDB)
# print(len(List_Wild_PDB))
input_protein_list1 = []
for i in range(len(input_protein_list)):
    if input_protein_list[i][0] in wild_pdb_id and input_protein_list[i][2] == len_Input_Protein:
        input_protein_list1.append(input_protein_list[i])
if len(input_protein_list1) == 0:
    input_protein_list1 = input_protein_list
# print(input_protein_list1)
List_Wild_PDB1 = []
temp1 = []
for i in range(len(List_Wild_PDB)):
    if len(List_Wild_PDB[i]) > 1:
        for j in range(len(List_Wild_PDB[i])):
            temp_pdb.append(int(List_Wild_PDB[i][j][2]))
        if int(uniprot_len_family[i]) in temp_pdb:
            w = temp_pdb.index(uniprot_len_family[i])
            List_Wild_PDB1.append(List_Wild_PDB[i][w])
            del temp_pdb[:]
        else:
            for k in range(len(temp_pdb)):
                temp1.append(len_Input_Protein-temp_pdb[k])
            s = min(temp1)
            n = temp1.index(s)
            List_Wild_PDB1.append(List_Wild_PDB[i][n])
            del temp1[:]
            del temp_pdb[:]
    else:
        List_Wild_PDB1.extend(List_Wild_PDB[i])
# print(List_Wild_PDB1)
input_protein_list2 = []
for i in range(len(input_protein_list1)):
    temp_pdb.append(int(input_protein_list1[i][2]))
    s = max(temp_pdb)
    n = temp_pdb.index(s)
    input_protein_list2 = input_protein_list1[n]
# print(input_protein_list2)
# test = []
# for i in range(len(List_Wild_PDB1)):
#     test.append(List_Wild_PDB1[i][0])
for k in range(len(List_Wild_PDB1)):
    if List_Wild_PDB1[k] == input_protein_list2:
        List_Wild_PDB1.remove(List_Wild_PDB1[k])
        break
Wild_PDB = []
for i in range(len(List_Wild_PDB1)):
    if type(List_Wild_PDB1[i][0]) is list:
        Wild_PDB.append(List_Wild_PDB1[i][0][0])
    else:
        Wild_PDB.append(List_Wild_PDB1[i][0])
print(Wild_PDB)


