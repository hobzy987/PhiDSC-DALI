import os
# from Filter_PDB import List_Wild_PDB1, input_protein_list2
#
#
# # target_list = []
# # for i in range(len(List_Wild_PDB1)):
# #     # for j in List_Wild_PDB1[i]:
# #     chain_comb = List_Wild_PDB1[i][0] + List_Wild_PDB1[i][1]
# #     target_list.append(chain_comb)
# # input_protein_list3 = input_protein_list2[0] + input_protein_list2[1]
# # for i in range(len(target_list)):
# #     var1 = './DaliLite.v5/bin/dali.pl --cd1 %s --cd2 %s --dat1 ./DAT --dat2 ./DAT --title "output options" --outfmt ' \
# #            '"summary,alignments,equivalences,transrot" --clean 2> err' % ( target_list[i], input_protein_list3)
# #     p = subprocess.Popen(var1, shell=True, stderr=subprocess.PIPE)
# #     p.communicate()
# def reemovNestings(l):
#     output = []
#     for i in range(len(l)):
#         if len(l[i]) == 1:
#             output.append(l[i][0])
#         else:
#             output.append(l[i])
#     return output
#
#
# Input_Protein = str(input_protein_list2[0].lower() + input_protein_list2[1])
# List_Wild_PDB2 = reemovNestings(List_Wild_PDB1)
# Dali_Input = []
# for i in range(len(List_Wild_PDB1)):
#     Dali_Input.append(str(List_Wild_PDB2[i][0] + List_Wild_PDB2[i][1]))
#
# alignment = [[[] for _ in range(2)] for _ in range(len(List_Wild_PDB1))]
# end_alignment = [[[] for _ in range(2)] for _ in range(len(List_Wild_PDB1))]
# for k in range(len(Dali_Input)):
#     with open('%s.txt' % Dali_Input[k]) as f:
#         s = f.read()
#         w = s.splitlines()
#         for j in range(len(w)):
#             if len(w[j]) != 0:
#                 if w[j][0] == 'Q':
#                     alignment[k][1].append(w[j].split()[1])
#                     end_alignment[k][1].append(w[j].split()[2])
#                 if w[j][0] == 'S':
#                     alignment[k][0].append(w[j].split()[1])
#                     end_alignment[k][0].append(w[j].split()[2])
print (os. getcwd())
