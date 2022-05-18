import Filter_PDB

PDB_Protein_Input = Filter_PDB.input_protein_list2[0]
unq_Ston_PDB = Filter_PDB.Wild_PDB

import os
import sys


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
    url = 'https://files.rcsb.org/pub/pdb/data/structures/all/pdb/pdb%s.ent.gz' %pdb_code.lower()
    destination_file = os.path.join(download_folder, filename)
    urlretrieve(url, destination_file)

url = 'https://files.rcsb.org/pub/pdb/data/structures/all/pdb/pdb%s.ent.gz' % PDB_Protein_Input.lower()
filename = '%s.ent.gz' % PDB_Protein_Input
destination_file = os.path.join(download_folder, filename)
urlretrieve(url, destination_file)

