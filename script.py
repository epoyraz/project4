


w, h = 19, 19
Matrix = [[0 for x in range(w)] for y in range(h)]

Matrix[0][0] = -0.724
Matrix[0][1] = -0.903
Matrix[0][2] = -0.722
Matrix[0][3] = -0.322

print Matrix





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
for model in structure:
   for chain in model:
      for residue in chain:
         for atom in residue:
               print atom.get_vector(), atom.name
