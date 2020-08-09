from rdkit import Chem
from rdkit.Chem import AllChem
import os


class LigandExpo(dict):
    def __init__(self):
        dataname = os.path.dirname(__file__) + '/../data/Components-smiles-oe.smi'
        if not os.path.exists(dataname):
            dataname = os.path.dirname(__file__) + '/../../../../cbiprep/data/Components-smiles-oe.smi'
        for lineno, line in enumerate(open(dataname, 'rt'), start=1):
            line = line.rstrip()
            its = line.split('\t')
            if len(its) < 2:
                continue
            if its[0] == '':
                its[0] = _save
            smiles = its[0]
            name = its[1]
            self[name] = smiles
            _save = its[0]

    def assign(self, ligand_atoms, name):
        smi = self[name]
        rdkit_mol = Chem.MolFromPDBBlock(str(ligand_atoms))
        template = Chem.MolFromSmiles(smi)
        mol = AllChem.AssignBondOrdersFromTemplate(template, rdkit_mol)
        return mol
