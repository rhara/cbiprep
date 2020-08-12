import numpy as np


def GetTopologicalMatrix(ligand_atoms, pocket_atoms, ligand_thres=50, pocket_thres=400, distance_thres=4.0):
    if ligand_thres < len(ligand_atoms):
        return None
    if pocket_thres < len(pocket_atoms):
        return None
    ligand_adjmat = ligand_atoms.get_distance_based_adjacency_matrix(ligand_atoms, thres=distance_thres, diagzero=True)
    pocket_adjmat = pocket_atoms.get_distance_based_adjacency_matrix(pocket_atoms, thres=distance_thres, diagzero=True)
    complex_adjmat = ligand_atoms.get_distance_based_adjacency_matrix(pocket_atoms, thres=distance_thres)
    mat = np.zeros((ligand_thres + pocket_thres, ligand_thres + pocket_thres), dtype=int)
    mat[0:len(ligand_atoms),0:len(ligand_atoms)] = ligand_adjmat
    mat[ligand_thres:ligand_thres + len(pocket_atoms), ligand_thres:ligand_thres+len(pocket_atoms)] = pocket_adjmat
    mat[0:len(ligand_atoms), ligand_thres:ligand_thres+len(pocket_atoms)] = complex_adjmat
    mat[ligand_thres:ligand_thres + len(pocket_atoms), 0:len(ligand_atoms)] = complex_adjmat.T
    return mat
