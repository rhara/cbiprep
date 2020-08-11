from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG
import subprocess as sp


def draw_mol(mol, w=300, h=150):
    d = rdMolDraw2D.MolDraw2DSVG(w, h)
    d.DrawMolecule(mol)
    d.FinishDrawing()
    display(SVG(d.GetDrawingText()))


def check_dir(pdb_code, workdir='work'):
    print(sp.check_output(f'ls -l {workdir}/{pdb_code}', shell=True, universal_newlines=True))

