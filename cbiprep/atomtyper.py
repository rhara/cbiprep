from rdkit import Chem, RDConfig
from rdkit.Chem import ChemicalFeatures
import numpy as np
import os

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


def GetAtomTypes(mol):
    """
    Family Acceptor => Acc
    Family Donor = Don
    Family Hydrophobe => Hph
    Family LumpedHydrophobe = LHp
    Family NegIonizable => Neg
    Family PosIonizable => Pos
    Family ZnBinder => Znb
    Family Aromatic => Aro

    C, N, O, F, P, S, Cl, Br, X,
    Aromatic, sp3, sp2, sp, Ring,
    Acceptor, Donor, Hydrophobe, LumpedHydrophobe, NegIonizable, PosIonizable, ZnBinder
    """
    fdef_name = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
    factory = ChemicalFeatures.BuildFeatureFactory(fdef_name)
    feats = factory.GetFeaturesForMol(mol)
    features = {}
    for sym in ['Acc', 'Don', 'Hyd', 'Lum', 'Neg', 'Pos', 'ZnB', 'Aro']:
        features[sym] = set()
    for feat in feats:
        sym = feat.GetFamily()[:3]
        [features[sym].add(i) for i in feat.GetAtomIds()]
    features_lb = {}
    for k, v in features.items():
        for i in v:
            if i not in features_lb:
                features_lb[i] = []
            features_lb[i].append(k)

    vecs = []
    for atom in mol.GetAtoms():
        atom_feat = np.zeros(21, dtype=int)
        an = atom.GetAtomicNum()
        if an == 6:
            atom_feat[0] = 1
        elif an == 7:
            atom_feat[1] = 1
        elif an == 8:
            atom_feat[2] = 1
        elif an == 9:
            atom_feat[3] = 1
        elif an == 15:
            atom_feat[4] = 1
        elif an == 16:
            atom_feat[5] = 1
        elif an == 17:
            atom_feat[6] = 1
        elif an == 35:
            atom_feat[7] = 1
        else:
            atom_feat[8] = 1
        hyb = atom.GetHybridization()
        if hyb == Chem.HybridizationType.SP3:
            atom_feat[10] = 1
        elif hyb == Chem.HybridizationType.SP2:
            atom_feat[11] = 1
        elif hyb == Chem.HybridizationType.SP:
            atom_feat[12] = 1
        if atom.IsInRing():
            atom_feat[13] = 1
        idx = atom.GetIdx()
        if idx in features_lb:
            ff = features_lb[idx]
            if 'Aro' in ff:
                atom_feat[9] = 1
            if 'Acc' in ff:
                atom_feat[14] = 1
            if 'Don' in ff:
                atom_feat[15] = 1
            if 'Hyd' in ff:
                atom_feat[16] = 1
            if 'Lum' in ff:
                atom_feat[17] = 1
            if 'Neg' in ff:
                atom_feat[18] = 1
            if 'Pos' in ff:
                atom_feat[19] = 1
            if 'ZnB' in ff:
                atom_feat[20] = 1
        vecs.append(atom_feat)
    vecs = np.array(vecs)
    return vecs
