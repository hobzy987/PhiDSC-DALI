import os
import re
import urllib

from robobrowser.forms import form
from Download_PDBs import PDB_Protein_Input
from Download_PDBs import unq_Ston_PDB1
from robobrowser import RoboBrowser

br = RoboBrowser()
br.open('https://zhanglab.ccmb.med.umich.edu/TM-align/')
upload_form = br.get_form()
os.chdir('/Users/mohamadhoballa/PycharmProjects/tmaligners/PDBIDs')
filename = "Alignment-%s.txt" % PDB_Protein_Input
f = open(filename, "w+")
for pdb in unq_Ston_PDB1:
    upload_form['pdb1'].value = open('/Users/mohamadhoballa/PycharmProjects/tmaligners/PDBIDs/%s.pdb' % PDB_Protein_Input, 'r')
    upload_form['pdb2'].value = open('/Users/mohamadhoballa/PycharmProjects/tmaligners/PDBIDs/%s.pdb' % pdb, 'r')
    # upload_form["email"] = "mhdhoballa@gmail.com"
    br.submit_form(upload_form)
    src = str(br.parsed())
    start = 'TM-align/tmp/'
    end = '.html'
    alignment_code = ((src.split(start))[1].split(end)[0])
    urlpage = 'https://zhanglab.ccmb.med.umich.edu/TM-align/tmp/%s.html' % alignment_code
    page = urllib.request.urlopen(urlpage)
    mybytes = page.read()
    mystr = mybytes.decode("utf8")
    start1 = '(":" denotes aligned residue pairs of d < 5.0 A, "." denotes other aligned residues)'
    end1 = '</pre'
    result = ((mystr.split(start1))[1].split(end1)[0])
    page.close()
    result = os.linesep.join([s for s in result.splitlines() if s])
    with open(filename, "a") as text_file:
        text_file.write(result)
    text_file.close()

    # print(result)
    # print(s)

