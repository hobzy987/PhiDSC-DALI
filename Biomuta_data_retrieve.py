import os
from urllib.request import urlopen
import requests
from xml.etree.ElementTree import fromstring
# from DALI_Align import Dali_Input as unq_Ston_PDB
from Filter_PDB import UniProt_IDs, Input_Protein
from Readfamily import Input_Protein as input_uniprot
# from align_matrix import end_protein
# from Bio import SeqIO
# import urllib


def protein_length(uniprot_id):
    data = urlopen("http://www.uniprot.org/uniprot/" + uniprot_id + ".txt").read()
    data1 = data.decode()
    data = data1.split('\n')
    len_Protein = [int(s) for s in data[0].split() if s.isdigit()]
    return len_Protein


pdb_mapping_url = 'http://www.rcsb.org/pdb/rest/das/pdb_uniprot_mapping/alignment'
uniprot_url = 'http://www.uniprot.org/uniprot/{}.xml'


def get_uniprot_accession_id(response_xml):
    root = fromstring(response_xml)
    return next(
        el for el in root.getchildren()[0].getchildren()
        if el.attrib['dbSource'] == 'UniProt'
    ).attrib['dbAccessionId']


def get_uniprot_protein_name(uniport_id):
    uinprot_response = requests.get(
        uniprot_url.format(uniport_id)
    ).text
    return fromstring(uinprot_response).find(
        './/{http://uniprot.org/uniprot}recommendedName/{http://uniprot.org/uniprot}fullName'
    ).text


def map_pdb_to_uniprot(pdb_id):
    pdb_mapping_response = requests.get(
        pdb_mapping_url, params={'query': pdb_id}
    ).text
    uniprot_id = get_uniprot_accession_id(pdb_mapping_response)
    uniprot_name = get_uniprot_protein_name(uniprot_id)
    return {
        'pdb_id': pdb_id,
        'uniprot_id': uniprot_id,
        'uniprot_name': uniprot_name
    }


def remove_dup_list(x):
    return list(dict.fromkeys(x))


def Union(lst1, lst2):
    final_list = lst1 + lst2
    remove_dup_list(final_list)
    return final_list


# UniProt_IDs = []
# for i in unq_Ston_PDB:
#     x = map_pdb_to_uniprot(i)
#     UniProt_IDs.append(x['uniprot_id'])






len_UniProt_IDs = []
for i in UniProt_IDs:
    len_UniProt_IDs.append(protein_length(i))

path = "../PhiDSC-DALI"
os.chdir('../PhiDSC-DALI')
file_name = os.path.join(path, "biomuta-master.csv")
my_file = open("biomuta-master.csv", 'r')
similar: str = my_file.read()
with open("biomuta-master.csv", 'r') as file:
    x = file.readlines()
temp_biomuta_list =[[] for _ in range(len(UniProt_IDs))]
for i in range(len(UniProt_IDs)):
    for j in range(len(x)):
        if UniProt_IDs[i] in x[j]:
            temp_biomuta_list[i].append(x[j].split('","'))
# input_uniprot = map_pdb_to_uniprot(Input_Protein)['uniprot_id']
mutation_position_input = []
freq_mutation = []
for k in range(len(x)):
    if input_uniprot in x[k]:
        z = x[k].split('","')
        key = int(z[10])
        value = int(z[13])
        g = {'key': key,  'value': value}
        freq_mutation.append(g)
        mutation_position_input.append(key)
mutation_position_input = sorted(remove_dup_list(mutation_position_input))
print(mutation_position_input)
from collections import defaultdict
result = defaultdict(int)

for d in freq_mutation:
    result[d['key']] += int(d['value'])
freq_mutation = [{'key': key, 'value': value} for key, value in result.items()]
temp_biomuta_list1 = [[] for _ in range(len(UniProt_IDs))]
occerance = [{} for _ in range(len(UniProt_IDs))]
for i in range(len(UniProt_IDs)):
    for j in range(len(temp_biomuta_list[i])):
        temp_biomuta_list1[i].append(int(temp_biomuta_list[i][j][10]))
        if len(occerance[i]) == 0:
            occerance[i].update({temp_biomuta_list[i][j][10]: int(temp_biomuta_list[i][j][13])})
        elif temp_biomuta_list[i][j][10] in occerance[i]:
            occerance[i][temp_biomuta_list[i][j][10]] = int(occerance[i][temp_biomuta_list[i][j][10]]) + int(temp_biomuta_list[i][j][13])
        else:
            occerance[i].update({temp_biomuta_list[i][j][10]: int(temp_biomuta_list[i][j][13])})


temp_biomuta_list2 = []
for i in range(len(temp_biomuta_list1)):
    x = remove_dup_list(temp_biomuta_list1[i])
    temp_biomuta_list2.append(sorted(x))
# print(temp_biomuta_list2)
my_file1 = open("hotspots_v2.csv", 'r')
similar: str = my_file1.read()
with open("hotspots_v2.csv", 'r') as file1:
    y = file1.readlines()
temp_hotspot1 = []
for i in range(len(y)):
    temp_hotspot1.append(y[i].split(','))
my_file2 = open("hugo_to_uniprot.txt", 'r')
similar: str = my_file2.read()
with open("hugo_to_uniprot.txt", 'r') as file2:
    z = file2.readlines()
HGNC_To_Uni = []
for i in range(len(z)):
    HGNC_To_Uni.append(z[i].split('\t'))
HGNC_Fam_List = []
HGNC_Input_List = []
for k in range(len(UniProt_IDs)):
    for i in range(len(HGNC_To_Uni)):
        if UniProt_IDs[k] in HGNC_To_Uni[i][9]:
            HGNC_Fam_List.append(HGNC_To_Uni[i][1])
        if Input_Protein in HGNC_To_Uni[i][9]:
            HGNC_Input_List.append(HGNC_To_Uni[i][1])
# print(HGNC_Fam_List)
hotspot_chang = [[] for _ in range(len(HGNC_Fam_List))]
hotspot_chang_input = []
for i in range(len(HGNC_Fam_List)):
    for j in range(len(temp_hotspot1)):
        if HGNC_Fam_List[i] == temp_hotspot1[j][0]:
            hotspot_chang[i].append(int(temp_hotspot1[j][1]))
        if HGNC_Input_List[0] == temp_hotspot1[j][0]:
            hotspot_chang_input.append(int(temp_hotspot1[j][1]))
hotspot_chang_input = sorted(remove_dup_list(hotspot_chang_input))
for i in range(len(hotspot_chang)):
    hotspot_chang[i] = sorted(remove_dup_list(hotspot_chang[i]))
    # if len(hotspot_chang[i]) == 0:
    #     hotspot_chang[i].append(0)
for i in range(len(temp_biomuta_list2)):
    temp_biomuta_list2[i] = sorted(remove_dup_list(Union(temp_biomuta_list2[i], hotspot_chang[i])))
mutation_position_input = sorted(remove_dup_list(Union(mutation_position_input, hotspot_chang_input)))
# print(hotspot_chang)
