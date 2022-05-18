import Readfamily
import os


# meta function
def remove_empty_lists(l):
    keep_going = True
    prev_l = l
    while keep_going:
        # call remover on the list
        new_l = remover(prev_l)
        # are they identical objects?
        if new_l == prev_l:
            keep_going = False
        # set prev to new
        prev_l = new_l
    # return the result
    return new_l


# function
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


uniprot_list = Readfamily.uniprot_list
Input_Protein = Readfamily.Input_Protein
path = "/Users/mohamadhoballa/PycharmProjects/tmaligners"
file_name = os.path.join(path, "UniProt_to_PDB.txt")
my_file = open("UniProt_to_PDB.txt", 'r')
y = my_file.read()
with open("UniProt_to_PDB.txt", 'r') as file:
    x = file.readlines()
Uni_Split_PDB = []
for i in range(len(x)):
    h = x[i].split("\t")
    Uni_Split_PDB.append(h)
# print(Uni_Split_PDB)
PDB_List = []
f = []
k = []
UniProt_IDs = []
for ID in uniprot_list:
    for j in range(len(Uni_Split_PDB)):
        if ID in Uni_Split_PDB[j]:
            f.append(Uni_Split_PDB[j][3])
for i in range(len(f)):
    if len(f[i]) > 3:
        UniProt_IDs.append(uniprot_list[i])
    z = f[i].split(";")
    k.append(z)
    if len(k[i]) > 1:
        PDB_List.append(k[i][:-1])
# for j in range(len(PDB_List)):
#     PDB_List[j].pop(-1)
# print(PDB_List)
Input_PDB = []
q = []
for i in range(len(Uni_Split_PDB)):
    if Input_Protein in Uni_Split_PDB[i]:
        q.append(Uni_Split_PDB[i][3])
        z = q[0].split(";")
        Input_PDB = z
        Input_PDB.pop(-1)
# print("The List of PDBs of Your Input is: " + str(Input_PDB))
file_name = os.path.join(path, "wild_pdb_id.txt")
my_file = open("wild_pdb_id.txt", 'r')
y = my_file.read()
with open("wild_pdb_id.txt", 'r') as file:
    x = file.readlines()
wild_pdb_id = []
for i in range(len(x)):
    wild_pdb_id = x[i].split(", ")

remover(wild_pdb_id)
remove_empty_lists(wild_pdb_id)


# print(wild_pdb_id)
List_PDB = [[] for _ in range(len(PDB_List))]
for j in range(len(PDB_List)):
    for k in range(len(PDB_List[j])):
        if PDB_List[j][k] in wild_pdb_id:
            List_PDB[j].append(PDB_List[j][k])
    if not List_PDB[j]:
        List_PDB[j].append(PDB_List[j])
        List_PDB = reemovNestings(List_PDB)

# print(List_PDB)

result = []
for element in Input_PDB:
    if element in wild_pdb_id:
        result.append(element)
# print("The List of Wild PDBs of Your Input is: " + str(result))
