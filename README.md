### TODO

For the time being, please use Jupyter notebooks.

### LOG

```
ver 0.1.31  1_index_and_pdb_download.ipynb - refactored
|           2_ligand_and_pocket.ipynb      - refactored
ver 0.1.30                                 - subtle changes
ver 0.1.29                                 - subtle changes
ver 0.1.28  Dataset.ipynb                  - fixed
ver 0.1.27  cbipre/atomtyper.py            - added GetAtomTypes function
ver 0.1.26  Maserati.ipynb                 - output is compact tarball
ver 0.1.25  notebooks                      - slightly modified
|           cbiprep/pdbatoms.py            - PDBAtoms::GetResFragment method added
ver 0.1.24  notebooks                      - modified for Hasegawa selection
ver 0.1.23  Maserati.ipynb                 - essentials added
ver 0.1.22  Maserati.ipynb                 - newly added, compile dataset from all the calculations
ver 0.1.21  cbiprep/pdbatoms.py            - fixed for in-chain duplicated ligands
ver 0.1.20  SminaSelfDocking.ipynb         - small changes
ver 0.1.19  cbiprep/smina.py               - small changes
|                                          - test_smina.ipynb no longer works, so deleted
ver 0.1.18  papers/                        - collection of publications
ver 0.1.17  SminaSelfDocking.ipynb         - newly added
ver 0.1.16  cbiprep/matrix.py              - newly added
|           DataSet.ipynb                  - newly added
|           cbiprep/pdbatoms.py            - remove NMR Hs and altLoc
ver 0.1.15                                 - fixed in general
ver 0.1.14  test_smina.ipynb               - smina function to go through all PDBs
|                                          - RMSD is written out
|                                          - Cleaned
|           cbiprep/smina.py               - newly added
```

### update the local repository

```
git clone
```

### install/update by using pip

```
python setup.py sdist
pip install dist/cbiprep-X.Y.Z.tar.gz
```

or

```
pip install .
```

### run prep

```
python -m cbiprep.pdbbind_index
python -m cbiprep.pdbdl
python -m cbiprep.ligpocket
```
