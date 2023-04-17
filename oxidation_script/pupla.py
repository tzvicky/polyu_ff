#!/usr/bin/env python

"""
pupla.py.py, by L. Cwiklik, 2022, cwiklik[at]gmail.com
Randomly oxidizes polyurethne chains.
Usage:
pupla.py GRO_FILE OX_RATIO
Input parameters:
GRO_FILE - input gro file
OX_RATIO - oxidation ratio (0.0 to 1.0)
Prints the resulting gro file.
"""

import sys, copy, random

class Atom:
    def __init__(self, symbol, number, x, y, z, vx, vy, vz):
        self.symbol = symbol
        self.number = number
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz

class Molecule:
    def __init__(self):
        self.resname = ""
        self.atoms = []

class System:
    def __init__(self, filename):
        self.title = ""
        self.n_atoms = 0
        self.molecules = None
        self.box_x = 0.0
        self.box_y = 0.0
        self.box_z = 0.0

        self.__readgro__(filename)

    def __readgro__(self, filename):

        f = open(sys.argv[1])
        lines = f.readlines()
        title = lines[0].strip()
        n_atoms = int(lines[1])
        box_x, box_y, box_z = lines[-1].split()
        box_x = float(box_x)
        box_y = float(box_y)
        box_z = float(box_z)

        molecules = []

        n_residue_previous = -1
        atom_n = -1
        for i in range(2, len(lines)-1):

            n_residue = int(lines[i][:5])

            #check if a new Molecule object should be created
            if n_residue != n_residue_previous:
                molecule = Molecule()
                resname = lines[i][5:10].strip()
                molecule.resname = resname

            spl = lines[i].split()
            atom_symbol = lines[i][10:15].strip()
            atom_n = int(lines[i][15:20].strip())
            atom_x = float(lines[i][20:28].strip())
            atom_y = float(lines[i][28:36].strip())
            atom_z = float(lines[i][36:44].strip())
            if len(spl) > 6:
                atom_vx = float(lines[i][44:52].strip())
                atom_vy = float(lines[i][52:60].strip())
                atom_vz = float(lines[i][60:68].strip())
            else:
                atom_vx = atom_vy = atom_vz = 0.0

            atom = Atom(atom_symbol, atom_n, atom_x, atom_y, atom_z,
                        atom_vx, atom_vy, atom_vz)

            molecule.atoms.append(atom)

            #checking again
            if n_residue != n_residue_previous:
                molecules.append(molecule)
                n_residue_previous = n_residue

        f.close()

        self.title = title
        self.n_atoms = n_atoms
        self.molecules = molecules
        self.box_x = box_x
        self.box_y = box_y
        self.box_z = box_z
	
    def put_first(self, resname):
        self.molecules.sort(cmp=lambda x,y: int(x.resname!=resname and y.resname==resname)-1)

    def molecule_sorting_index(self, molecule):
        return self.sorting_list.index(molecule.resname) if molecule.resname in self.sorting_list else -1

    def sort_by_list(self, sorting_list):
        self.sorting_list = sorting_list #to be used by self.molecule_sorting_index(self, resname)
        self.molecules.sort(key=self.molecule_sorting_index)
        

    def order(self, order):
        new_molecules = []
        for i in range(len(order)):
            new_molecules.append(0)
        for i in range(len(self.molecules)):
            mol = self.molecules[i]
            resid = i+1
            if resid in order:
                new_molecules[order.index(resid)] = mol
            else:
                new_molecules.append(mol)
        self.molecules = new_molecules

    def oxidize(self, ox_ratio):

        new_molecules = []
        for i in range(len(self.molecules)):
            m = self.molecules[i]
            if (m.resname not in ('FRA', 'FRC')):
                new_molecules.append(m)
                continue

            if ox_ratio < random.random():
                new_molecules.append(m)
                continue

            newmol = Molecule()

            if m.resname == "FRA":
                if random.random() < 0.5:
                    new_resname = "AOX1"
                    newmol.resname = new_resname
                    newmol.atoms = copy.deepcopy(m.atoms[7:])
                    newmol.atoms[0].symbol = "O5"
                    newmol.atoms[1].symbol = "H16"
                else:
                    new_resname = "AOX2"
                    newmol.resname = new_resname
                    newmol.atoms.append(copy.deepcopy(m.atoms[19]))
                    newmol.atoms += (copy.deepcopy(m.atoms[21:]))
                    newmol.atoms[0].symbol = "O5"
                    newmol.atoms[1].symbol = "H16"

            if m.resname == "FRC":
                if random.random() < 0.5:
                    new_resname = "COX1"
                    newmol.resname = new_resname
                    newmol.atoms = copy.deepcopy(m.atoms[:29])
                    newmol.atoms[-2].symbol = "O5"
                    newmol.atoms[-1].symbol = "H16"
                else:
                    new_resname = "COX2"
                    newmol.resname = new_resname
                    newmol.atoms = copy.deepcopy(m.atoms[:16])
                    newmol.atoms[-2].symbol = "O5"
                    newmol.atoms[-1].symbol = "H16"
		
            
            
	
            new_molecules.append(newmol)


        self.molecules = new_molecules

    def printgro(self):
        print(self.title)
        n_atoms = 0
        for m in self.molecules:
            n_atoms += len(m.atoms)
        print(n_atoms)
        n_atom = 0
        n_molecule = 0
        out_line = ""
        for m in self.molecules:
            n_molecule += 1
            if n_molecule == 100000:
                n_molecule = 0 #loop back after 99,999 molecules
            
            
            for a in m.atoms:
                n_atom += 1
                if n_atom == 100000:
                    n_atom = 0 #loop back after 99,999 atoms
                out_line = str(n_molecule).rjust(5)
                out_line += m.resname.ljust(5)
                out_line += a.symbol.rjust(5)
                out_line += str(n_atom).rjust(5)
                out_line += str('%.3f'%a.x).rjust(8)
                out_line += str('%.3f'%a.y).rjust(8)
                out_line += str('%.3f'%a.z).rjust(8)
                #out_line += str('%.4f'%a.vx).rjust(8)
                #out_line += str('%.4f'%a.vy).rjust(8)
                #out_line += str('%.4f'%a.vz).rjust(8)

                print(out_line)

        print(self.box_x, self.box_y, self.box_z)

    """
    def printgro(self):
        print(self.title)
        n_atoms = 0
        for m in self.molecules:
            n_atoms += len(m.atoms)
        print(n_atoms)
        n_atom = 0
        n_molecule = 0
        out_line = ""
        for m in self.molecules:
            n_molecule += 1
            
            for a in m.atoms:
                n_atom += 1
                out_line = str(n_molecule).rjust(5)
                out_line += m.resname.ljust(5)
                out_line += a.symbol.rjust(5)
                out_line += str(n_atom).rjust(5)
                out_line += str('%.3f'%a.x).rjust(8)
                out_line += str('%.3f'%a.y).rjust(8)
                out_line += str('%.3f'%a.z).rjust(8)
                out_line += str('%.4f'%a.vx).rjust(8)
                out_line += str('%.4f'%a.vy).rjust(8)
                out_line += str('%.4f'%a.vz).rjust(8)

                print(out_line)

        print(self.box_x, self.box_y, self.box_z)

    """

    def trim(self, x_min, x_max, y_min, y_max, z_min, z_max):
    
        mol_to_remove = []
        for i in range(len(self.molecules)):
            m = self.molecules[i]
            remove = False
            for atom in m.atoms:
                if atom.x < x_min or atom.x > x_max or \
                       atom.y < y_min or atom.y > y_max or \
                       atom.z < z_min or atom.z > z_max:
                    remove = True
            if remove:
                mol_to_remove.append(i)
                
        new_molecules = []
        for i in range(len(self.molecules)):
            if i not in mol_to_remove:
                new_molecules.append(self.molecules[i])
        self.molecules = new_molecules

        # find new box sizes
        x_left = x_max # left<->max sic!
        x_right = x_min
        y_left = y_max
        y_right = y_min
        z_left = z_max
        z_right = z_min
        for m in self.molecules:
            for a in m.atoms:
                if a.x < x_left:
                    x_left = a.x
                if a.x > x_right:
                    x_right = a.x
                if a.y < y_left:
                    y_left = a.y
                if a.y > y_right:
                    y_right = a.y
                if a.z < z_left:
                    z_left = a.z
                if a.z > z_right:
                    z_right = a.z
            
        self.box_x = x_right - x_left
        self.box_y = y_right - y_left
        self.box_z = z_right - z_left
        

def main():

    if len(sys.argv) < 3:
        print(__doc__)
        return -1

    filename = sys.argv[1]
    ox_ratio = float(sys.argv[2])

    system = System(filename)
    system.oxidize(ox_ratio)

    #system.sort_by_list(('POPC', 'POXN', 'SOL'))

    system.printgro()

if __name__ == '__main__':
    main()
