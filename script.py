import math
import numpy as np
from os import listdir
from os.path import isfile, join


import numpy as np
import matplotlib.pyplot as plt


#N = 50
#x = np.random.rand(N)
#y = np.random.rand(N)
#colors = np.random.rand(N)
#area = np.pi * (15 * np.random.rand(N))**2  # 0 to 15 point radiuses




def class_by_atomname(atomname):
   if (atomname == 'N'):
      return 0
   if (atomname == 'CA'):
      return 1
   if (atomname == 'C' or atomname == 'O'):
      return 2
   return 3

def get_atom_type(atom,aa):
	if( atom == 'N'):
		return 'N'
	if( atom == 'CA'):
		if(aa == 'GLY'):
			return 'GCA'
		else:
			return 'CA'
	if( atom == 'C'):
		return 'C'
	if( atom == 'O'):
		return 'O'
	if( atom == 'CB'):
		if (aa== 'SER'):
			return 'SO'
		else:
			return 'CB'
	if (atom == 'CG'):
		if (aa == 'PRO'):
			return 'CB'
		if (aa == 'ASP'):
			return 'DO'
		if (aa == 'ASN'):
			return 'NN'
		if (aa == 'HIS'):
			return 'HN'
		if (aa == 'ARG' or aa == 'GLN' or aa == 'GLU' or aa == 'LEU' or aa == 'LYS' or aa == 'MET' or aa == 'PHE' or aa == 'TRP' or aa =='TYR'):
			return 'FC'
	if (atom == 'CG1'):
		if ( aa == 'ILE'):
			return 'FC'
		if ( aa == 'VAL'):
			return 'LC'
	if (atom == 'CG2'):
		if (aa=='THR'):
			return 'FC'
		if (aa == 'ILE' or aa == 'VAL'):
			return 'LC'
	if (atom == 'CD'):
		if (aa == 'LYS'):
			return 'KC'
		if (aa == 'GLU'):
			return 'DO'
		if (aa =='GLN'):
			return 'NN'
		if (aa == 'ARG'):
			return 'RNE'
		if (aa == 'ILE'):
			return 'LC'
		if (aa == 'PRO'):
			return 'CB'
	if (atom == 'CD1' or atom == 'CD2'):
		if(aa == 'PHE'):
			return 'FC'
		if (aa == 'TRP' or aa == 'TYR'):
			return 'FC'
		if (aa == 'LEU' or 'ILE'):
			return 'LC'
		if (aa == 'HIS'):
			return 'HN'
	if (atom == 'CE'):
		if ( aa == 'LYS'):
			return 'KN'
		if ( aa == 'MET'):
			return 'LC'
	if (atom == 'CE1'):
		if (aa == 'HIS'):
			return 'HN'
		if (aa == 'TYR'):
			return 'YC'
		if (aa == 'PHE'):
			return 'FC'
	if (atom == 'CE2'):
		if(aa == 'TYR'):
			return 'YC'
		if (aa == 'PHE'):
			return 'FC'
		if (aa == 'TRP'):
			return 'FC'
	if (atom == 'CE3'):
		if (aa == 'TRP'):
			return 'FC'
	if (atom == 'CZ'):
		if (aa == 'ARG'):
			return 'RNH'
		if (aa == 'TYR'):
			return 'YC'
		if (aa == 'PHE'):
			return 'FC'
	if (atom == 'CZ2' or atom == 'CZ3'):
		if(aa == 'TRP'):
			return 'FC'
	if (atom == 'CH2'):
		if(aa == 'TRP'):
			return 'FC'
	if (atom == 'NZ'):
		if (aa == 'LYS'):
			return 'KN'
	if (atom == 'NH1' or atom == 'NH2'):
		if(aa == 'ARG'):
			return 'RNH'
	if (atom == 'ND2'):
		if (aa == 'ASN'):
			return 'NN'
	if (atom == 'ND1'):
		if(aa == 'HIS'):
			return 'HN'
	if (atom == 'NE'):
		if(aa=='ARG'):
			return 'RNE'
	if (atom =='NE1'):
		if(aa=='TRP'):
			return 'HN'
	if (atom =='NE2'):
		if(aa=='HIS'):
			return 'HN'
		if(aa=='GLN'):
			return 'NN'
	if (atom == 'OD1'):
		if(aa == 'ASP'):
			return 'DO'
		if(aa == 'ASN'):
			return 'NN'
	if (atom == 'OD2'):
		if(aa== 'ASP'):
			return 'DO'
	if (atom == 'OE1'):
		if (aa == 'GLU'):
			return 'DO'
		if (aa== 'GLN'):
			return 'NN'
		if (aa == 'HIS'):
			return 'HN'
	if (atom == 'OE2'):
		if(aa =='GLU'):
			return 'DO'
	if (atom == 'OG'):
		if (aa == 'SER'):
			return 'SO'
	if (atom == 'OG1'):
		if (aa== 'THR'):
			return 'SO'
	if (atom == 'OH'):
		if(aa=='TYR'):
			return 'SO'
	if (atom == 'SD'):
		if(aa == 'MET'):
			return 'FC'
	if (atom == 'SG'):
		if (aa == 'CYS'):
			return 'CS'
	return 'rest', atom , aa

w, h = 18, 18
Matrix = [[0 for x in range(w)] for y in range(h)]
Matrix[0][0] = -0.724
Matrix[0][1] = -0.903
Matrix[0][2] = -0.722
Matrix[0][3] = -0.322
Matrix[0][4] = -0.331
Matrix[0][5] = -0.603
Matrix[0][6] = 1.147
Matrix[0][7] = 0.937
Matrix[0][8] = 0.752
Matrix[0][9] = 0.502
Matrix[0][10] = 0.405
Matrix[0][11] = 0.495
Matrix[0][12] = -0.116
Matrix[0][13] = -0.092
Matrix[0][14] = -0.412
Matrix[0][15] = -0.499
Matrix[0][16] = -1.005
Matrix[0][17] = -2.060
Matrix[1][0] = -0.119
Matrix[1][1] = -0.842
Matrix[1][2] = -0.850
Matrix[1][3] = -0.332
Matrix[1][4] = -0.365
Matrix[1][5] = -0.531
Matrix[1][6] = 1.139
Matrix[1][7] = 0.918
Matrix[1][8] = 0.852
Matrix[1][9] = 0.452
Matrix[1][10] = 0.445
Matrix[1][11] = 0.426
Matrix[1][12] = -0.044
Matrix[1][13] = -0.165
Matrix[1][14] = -0.539
Matrix[1][15] = -0.590
Matrix[1][16] = -1.123
Matrix[1][17] = -2.028
Matrix[2][0] = -0.008
Matrix[2][1] = -0.077
Matrix[2][2] = -0.704
Matrix[2][3] = -0.371
Matrix[2][4] = -0.308
Matrix[2][5] = -0.499
Matrix[2][6] = 1.140
Matrix[2][7] = 0.846
Matrix[2][8] = 0.788
Matrix[2][9] = 0.490
Matrix[2][10] = 0.408
Matrix[2][11] = 0.418
Matrix[2][12] = -0.045
Matrix[2][13] = -0.089
Matrix[2][14] = -0.415
Matrix[2][15] = -0.478
Matrix[2][16] = -0.963
Matrix[2][17] = -2.003
Matrix[3][0] = 0.048
Matrix[3][1] = 0.097
Matrix[3][2] = -0.011
Matrix[3][3] = -0.016
Matrix[3][4] = 0.205
Matrix[3][5] = -0.059
Matrix[3][6] = 1.414
Matrix[3][7] = 1.075
Matrix[3][8] = 1.152
Matrix[3][9] = 0.831
Matrix[3][10] = 0.758
Matrix[3][11] = 0.649
Matrix[3][12] = 0.383
Matrix[3][13] = 0.241
Matrix[3][14] = -0.073
Matrix[3][15] = -0.166
Matrix[3][16] = -0.650
Matrix[3][17] = -1.650
Matrix[4][0] = -0.059
Matrix[4][1] = -0.034
Matrix[4][2] = -0.047
Matrix[4][3] = 0.122
Matrix[4][4] = 0.182
Matrix[4][5] = -0.009
Matrix[4][6] = 1.502
Matrix[4][7] = 1.385
Matrix[4][8] = 1.192
Matrix[4][9] = 0.912
Matrix[4][10] = 0.872
Matrix[4][11] = 0.943
Matrix[4][12] = 0.434
Matrix[4][13] = 0.392
Matrix[4][14] = -0.062
Matrix[4][15] = -0.021
Matrix[4][16] = -0.342
Matrix[4][17] = -1.212
Matrix[5][0] = -0.007
Matrix[5][1] = 0.124
Matrix[5][2] = 0.088
Matrix[5][3] = 0.184
Matrix[5][4] = 0.135
Matrix[5][5] = -0.469
Matrix[5][6] = 1.310
Matrix[5][7] = 1.010
Matrix[5][8] = 0.859
Matrix[5][9] = 0.671
Matrix[5][10] = 0.591
Matrix[5][11] = 0.565
Matrix[5][12] = 0.122
Matrix[5][13] = 0.000
Matrix[5][14] = -0.438
Matrix[5][15] = -0.533
Matrix[5][16] = -1.034
Matrix[5][17] = -1.800
Matrix[6][0] = 0.000
Matrix[6][1] = 0.051
Matrix[6][2] = -0.017
Matrix[6][3] = -0.088
Matrix[6][4] = -0.098
Matrix[6][5] = 0.036
Matrix[6][6] = 3.018
Matrix[6][7] = 2.911
Matrix[6][8] = 1.157
Matrix[6][9] = 2.635
Matrix[6][10] = 1.848
Matrix[6][11] = 2.699
Matrix[6][12] = 1.587
Matrix[6][13] = 1.557
Matrix[6][14] = 1.004
Matrix[6][15] = 1.340
Matrix[6][16] = 1.252
Matrix[6][17] = 0.700
Matrix[7][0] = -0.106
Matrix[7][1] = -0.066
Matrix[7][2] = -0.207
Matrix[7][3] = -0.323
Matrix[7][4] = -0.111
Matrix[7][5] = -0.161
Matrix[7][6] = -0.003
Matrix[7][7] = 2.811
Matrix[7][8] = 0.989
Matrix[7][9] = 2.439
Matrix[7][10] = 1.622
Matrix[7][11] = 2.525
Matrix[7][12] = 1.361
Matrix[7][13] = 1.277
Matrix[7][14] = 0.685
Matrix[7][15] = 0.938
Matrix[7][16] = 0.814
Matrix[7][17] = 0.411
Matrix[8][0] = 0.126
Matrix[8][1] = 0.284
Matrix[8][2] = 0.151
Matrix[8][3] = 0.171
Matrix[8][4] = 0.112
Matrix[8][5] = 0.105
Matrix[8][6] = -1.341
Matrix[8][7] = -1.405
Matrix[8][8] = -1.978
Matrix[8][9] = 0.695
Matrix[8][10] = 1.386
Matrix[8][11] = 0.641
Matrix[8][12] = 1.002
Matrix[8][13] = 0.642
Matrix[8][14] = 0.618
Matrix[8][15] = 0.894
Matrix[8][16] = 0.792
Matrix[8][17] = -0.029
Matrix[9][0] = 0.069
Matrix[9][1] = 0.079
Matrix[9][2] = 0.047
Matrix[9][3] = 0.045
Matrix[9][4] = 0.027
Matrix[9][5] = 0.111
Matrix[9][6] = 0.331
Matrix[9][7] = 0.240
Matrix[9][8] = -1.088
Matrix[9][9] = 1.589
Matrix[9][10] = 1.395
Matrix[9][11] = 1.515
Matrix[9][12] = 1.022
Matrix[9][13] = 0.948
Matrix[9][14] = 0.457
Matrix[9][15] = 0.707
Matrix[9][16] = 0.586
Matrix[9][17] = 0.021
Matrix[9][0] = 0.117
Matrix[10][1] = 0.216
Matrix[10][2] = 0.110
Matrix[10][3] = 0.115
Matrix[10][4] = 0.131
Matrix[10][5] = 0.175
Matrix[10][6] = -0.311
Matrix[10][7] = -0.434
Matrix[10][8] = -0.253
Matrix[10][9] = -0.050
Matrix[10][10] = 1.300
Matrix[10][11] = 1.283
Matrix[10][12] = 0.931
Matrix[10][13] = 0.829
Matrix[10][14] = 0.484
Matrix[10][15] = 0.633
Matrix[10][16] = 0.389
Matrix[10][17] = -0.046
Matrix[11][0] = 0.203
Matrix[11][1] = 0.203
Matrix[11][2] = 0.116
Matrix[11][3] = 0.002
Matrix[11][4] = 0.197
Matrix[11][5] = 0.145
Matrix[11][6] = 0.535
Matrix[11][7] = 0.465
Matrix[11][8] = -1.002
Matrix[11][9] = 0.066
Matrix[11][10] = -0.022
Matrix[11][11] = 1.309
Matrix[11][12] = 0.921
Matrix[11][13] = 0.777
Matrix[11][14] = 0.265
Matrix[11][15] = 0.484
Matrix[11][16] = 0.274
Matrix[11][17] = -0.253
Matrix[12][0] = -0.033
Matrix[12][1] = 0.097
Matrix[12][2] = 0.028
Matrix[12][3] = 0.111
Matrix[12][4] = 0.064
Matrix[12][5] = 0.077
Matrix[12][6] = -0.201
Matrix[12][7] = -0.324
Matrix[12][8] = -0.266
Matrix[12][9] = -0.052
Matrix[12][10] = 0.002
Matrix[12][11] = -0.013
Matrix[12][12] = 0.559
Matrix[12][13] = 0.319
Matrix[12][14] = 0.301
Matrix[12][15] = 0.212
Matrix[12][16] = -0.155
Matrix[12][17] = -0.949
Matrix[13][0] = 0.421
Matrix[13][1] = 0.406
Matrix[13][2] = 0.413
Matrix[13][3] = 0.399
Matrix[13][4] = 0.451
Matrix[13][5] = 0.385
Matrix[13][6] = 0.199
Matrix[13][7] = 0.022
Matrix[13][8] = -0.197
Matrix[13][9] = 0.304
Matrix[13][10] = 0.330
Matrix[13][11] = 0.273
Matrix[13][12] = 0.190
Matrix[13][13] = -0.301
Matrix[13][14] = -0.159
Matrix[13][15] = -0.132
Matrix[13][16] = -0.338
Matrix[13][17] = -1.324
Matrix[14][0] = 0.107
Matrix[14][1] = 0.039
Matrix[14][2] = 0.094
Matrix[14][3] = 0.092
Matrix[14][4] = 0.004
Matrix[14][5] = -0.046
Matrix[14][6] = -0.348
Matrix[14][7] = -0.563
Matrix[14][8] = -0.214
Matrix[14][9] = -0.181
Matrix[14][10] = -0.008
Matrix[14][11] = -0.233
Matrix[14][12] = 0.179
Matrix[14][13] = 0.149
Matrix[14][14] = -0.314
Matrix[14][15] = -0.487
Matrix[14][16] = -0.964
Matrix[14][17] = -1.410
Matrix[15][0] = 0.207
Matrix[15][1] = 0.175
Matrix[15][2] = 0.218
Matrix[15][3] = 0.185
Matrix[15][4] = 0.232
Matrix[15][5] = 0.045
Matrix[15][6] = 0.175
Matrix[15][7] = -0.123
Matrix[15][8] = 0.249
Matrix[15][9] = 0.256
Matrix[15][10] = 0.326
Matrix[15][11] = 0.173
Matrix[15][12] = 0.277
Matrix[15][13] = 0.362
Matrix[15][14] = 0.023
Matrix[15][15] = -0.687
Matrix[15][16] = -1.240
Matrix[15][17] = -1.784
Matrix[16][0] = 0.293
Matrix[16][1] = 0.235
Matrix[16][2] = 0.325
Matrix[16][3] = 0.294
Matrix[16][4] = 0.504
Matrix[16][5] = 0.137
Matrix[16][6] = 0.680
Matrix[16][7] = 0.345
Matrix[16][8] = 0.740
Matrix[16][9] = 0.728
Matrix[16][10] = 0.675
Matrix[16][11] = 0.556
Matrix[16][12] = 0.503
Matrix[16][13] = 0.749
Matrix[16][14] = 0.130
Matrix[16][15] = 0.040
Matrix[16][16] = -1.873
Matrix[16][17] = -2.402
Matrix[17][0] = 0.173
Matrix[17][1] = 0.264
Matrix[17][2] = 0.190
Matrix[17][3] = 0.288
Matrix[17][4] = 0.569
Matrix[17][5] = 0.306
Matrix[17][6] = 1.063
Matrix[17][7] = 0.876
Matrix[17][8] = 0.853
Matrix[17][9] = 1.097
Matrix[17][10] = 1.175
Matrix[17][11] = 0.964
Matrix[17][12] = 0.642
Matrix[17][13] = 0.697
Matrix[17][14] = 0.618
Matrix[17][15] = 0.430
Matrix[17][16] = 0.405
Matrix[17][17] = -3.742
print Matrix


Con = [[0 for x in range(4)] for y in range(4)]
#print "con matrix" ,Con

Con[0][0] = 3
Con[0][1] = 3
Con[0][2] = 2
Con[0][3] = 2
Con[1][0] = 3
Con[1][1] = 3
Con[1][2] = 3
Con[1][3] = 3
Con[2][0] = 4
Con[2][1] = 3
Con[2][2] = 3
Con[2][3] = 3
Con[3][0] = 3
Con[3][1] = 3
Con[3][2] = 2
Con[3][3] = 2

print Con[1][2]
print "Con matrix" , Con


amino_acids = ['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE','LEU', 'LYS','MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

types = ['N', 'CA', 'C', 'O' , 'GCA', 'CB' , 'KN', 'KC' , 'DO' , 'RNH' , 'NN', 'RNE', 'SO' , 'HN', 'YC' , 'FC' , 'LC', 'CS']

mutant_temp = {'155l' : 54.49, '156l':51.13, '157l':54.84, '158l':55.02, '159l':51.14, '160l':51.07, '161l':52.14, '162l':50.94, '163l':51.01, '164l':51.12, '165l':55.29, '166l':51.22}

x = []
y = []

print mutant_temp

from Bio.PDB.PDBParser import PDBParser

parser = PDBParser()


def get_temp(file):
	return mutant_temp[file.split('.')[0]]


blaber_files = [f for f in listdir('blaber_pdbs/') if isfile(join('blaber_pdbs/', f))]

if (blaber_files[0] == '.DS_Store'):
	blaber_files.pop(0)
print blaber_files

for blaber_pdb_filename in blaber_files:
	structure = parser.get_structure(blaber_pdb_filename.split('.')[0], 'blaber_pdbs/' +blaber_pdb_filename )
	header = parser.get_header()
	trailer = parser.get_trailer()
	aa_atoms = []
	
	x.append(get_temp(blaber_pdb_filename))
	print x
	
	for atom4 in list(structure.get_atoms()):
   		if (atom4.get_parent().get_resname() in amino_acids):
   			aa_atoms.append(atom4)

# get coordinates of atoms for calculation of distance between atoms
	a = []
	atoms = list(structure.get_atoms())

	at  = atoms[1].get_name()
	res = atoms[1].get_parent().get_resname()

	sum = 0
	num_atoms_of_type = [0] * 18
	num_of_surr_atoms_of_type = [0] * 18
	equals = 0
	print num_atoms_of_type
	for atom1 in aa_atoms:
		if (aa_atoms.index(atom1) > 3):
			break
		atom_type1 = get_atom_type(atom1.get_name(),atom1.get_parent().get_resname())
		if atom_type1 in types:
			index1 = types.index(atom_type1)
			num_atoms_of_type[types.index(atom_type1)] = num_atoms_of_type[types.index(atom_type1)] + 1
			for atom3 in aa_atoms:
				distance2 = atom1 - atom3
				if ( distance2 <= 6):
					num_of_surr_atoms_of_type[types.index(atom_type1)] = num_of_surr_atoms_of_type[types.index(atom_type1)]+1

		if (get_atom_type(atom1.get_name(),atom1.get_parent().get_resname()) not in types):
			print get_atom_type(atom1.get_name(),atom1.get_parent().get_resname())
	#print get_atom_type(atom1.get_name(),atom1.get_parent().get_resname())
		k_i =  atom1.get_parent().get_id()[1]
		sumlist = []
		for atom2 in aa_atoms:
			#print aa_atoms.index(atom2)
			if(aa_atoms.index(atom2) < aa_atoms.index(atom1) ):
				continue
			if (atom1 - atom2 <= 0.1) :
				#print "equals"
				equals = equals + 1
				continue
			atom_type2 = get_atom_type(atom2.get_name(),atom2.get_parent().get_resname())
			
			
			if atom_type2 in types:
				index2 = types.index(atom_type2)
   			k_j = atom2.get_parent().get_id()[1]
   			distance = atom1 - atom2
   			i = class_by_atomname(atom1.get_name())
   			j = class_by_atomname(atom2.get_name())
   			if(k_j >= k_i):
   					if not ( k_j - k_i <= Con[i][j]):
   						if ( distance <= 6):
   							if ( index1 < index2):
   								#print "==================="
   							   	#print "distance is ", distance  									
   								#print "atom2 type is : ", atom_type2
   								#print "atom1 type is : ", atom_type1
   								#print "k_j is : ", k_j
   								#print "k_i is : ", k_i
   								#print "i is : ", i
   								#print "j is : ", j
   								#print "Con[i][j] is : ", Con[i][j]
   								#print "Matrix value is : ",  Matrix[index1][index2]
								sumlist.append(Matrix[index1][index2])
   								sum += Matrix[index1][index2]
   							else:
   							   	#print "*********************"
   							   	#print "distance is ", distance  									
   								#print "atom2 type is : ", atom_type2
   								#print "atom1 type is : ", atom_type1
   								#print "k_j is : ", k_j
   								#print "k_i is : ", k_i
   								#print "i is : ", i
   								#print "j is : ", j
   								#print "Con[i][j] is : ", Con[i][j]
   								#print "Matrix value is : ",  Matrix[index2][index1]
   								sumlist.append(Matrix[index2][index1])
   								sum += Matrix[index2][index1]
			print atom1
			print sumlist
			#print aa_atoms.index(atom1)
			#print aa_atoms.index(atom2)
			#print sum
         #print atom1, atom2.get_parent()
	print sumlist
	print "types : " , types
	print "num_atoms_of_type : ", num_atoms_of_type
	print "num_of_surr_atoms_of_type : ", num_of_surr_atoms_of_type
	print "energy: ", sum
	print "number of atoms : ", len(aa_atoms)
	print equals
	y.append(sum)

plt.scatter(x, y)
plt.show()
