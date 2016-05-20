import math
import numpy as np



def class_by_atomname(atomname):
   if (atomname == 'N'):
      return 0
   if (atomname == 'CA'):
      return 1
   if (atomname == 'C' or atomname == 'O'):
      return 2
   return 3



w, h = 19, 19
Matrix = [[0 for x in range(w)] for y in range(h)]
Matrix[0][0] = -0.724
Matrix[0][1] = -0.903
Matrix[0][2] = -0.722
Matrix[0][3] = -0.322
# BITTE VERVOLLSTÄNDIGEN :)
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

from Bio.PDB.PDBParser import PDBParser

parser = PDBParser()

structure = parser.get_structure('155l', '155l.pdb')
header = parser.get_header()
trailer = parser.get_trailer()


# remove hydrogen according to  : http://pelican.rsvs.ulaval.ca/mediawiki/index.php/Manipulating_PDB_files_using_BioPython
for model in structure:
    for chain in model:
        for residue in chain:
            id = residue.id
            if id[0] == 'W':
               chain.detach_child(id)
        if len(chain) == 0:
            model.detach_child(chain.id)


# get coordinates of atoms for calculation of distance between atoms

atoms = list(structure.get_atoms())
for atom1 in atoms:
   k_i =  atom1.get_parent().get_id()[1]
   #print class_by_atomname(atom1.get_name())
   for atom2 in atoms:
      k_j = atom2.get_parent().get_id()[1]
      distance = atom1 - atom2

      if (abs(k_i - k_j > 5)):
            print "no neighbours"
      else:
         i = class_by_atomname(atom1.get_name())
         j = class_by_atomname(atom2.get_name())
         if ( k_i - k_j <= Con[i][j]):
            print "neighbours"
      if ( distance <= 6):
         atom1 - atom2
         #print atom1, atom2.get_parent()
