from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG


def draw_mol(mol, w=300, h=150):
    d = rdMolDraw2D.MolDraw2DSVG(w, h)
    d.DrawMolecule(mol)
    d.FinishDrawing()
    display(SVG(d.GetDrawingText()))

