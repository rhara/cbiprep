from rdkit import Chem
import numpy as np


class AtomTyper:
    """
    1 2 3 4 5 6 7  8  9
    C N O F P S Cl Br other
    """
    MAP = {6:1, 7:2, 8:3, 9:4, 15:5, 16:6, 17:7, 35:8}
    MAX = 9

    def get_atomtype(self, atom):
        an = atom.GetAtomicNum()
        return self.MAP.get(an, self.MAX)
    
    def __call__(self, mol):
        return np.array([self.get_atomtype(atom) for atom in mol.GetAtoms()])


class HybAtomTyper(AtomTyper):
    """
      1   2   3  4    5   6   7   8    9   10  11 12 13    14  15 16  17 18
    C sp3 sp2 sp ar N sp3 sp2 sp ar  O sp3 sp2 ar F  P  S  sp3 sp2 ar Cl other
    """
    MAP = {(6,3):1, (6,2):2, (6,1):3, (6,5):4, (7,3):5, (7,2):6, (7,1):7, (7,5):8, (8,3):9, (8,2):10, (8,5):11,
           (9,1):12, (15,0):13, (16,3):14, (16,2):15, (16,5):16, (17,1):17}
    MAX = 18

    def get_atomtype(self, atom):
        an = atom.GetAtomicNum()
        hyb = atom.GetHybridization()
        typ = {Chem.HybridizationType.SP:1, Chem.HybridizationType.SP2:2, Chem.HybridizationType.SP3:3}.get(hyb)
        arom = atom.GetIsAromatic()
        if arom:
            typ = 5
        if an == 15:
           typ = 0
        return self.MAP.get((an, typ), self.MAX)


def GetAtomVector(ligand_mol, pocket_atoms, ligand_thres=50, pocket_thres=400, atomtyper=HybAtomTyper):
    atomtyper = atomtyper()
    ligand_types = atomtyper(ligand_mol)
    pocket_mol = Chem.MolFromPDBBlock(str(pocket_atoms))
    pocket_types = atomtyper(pocket_mol)
    types_vec = np.zeros(ligand_thres + pocket_thres, dtype=int)
    types_vec[0:ligand_mol.GetNumAtoms()] = ligand_types
    types_vec[ligand_thres:ligand_thres + len(pocket_atoms)] = pocket_types
    return types_vec
