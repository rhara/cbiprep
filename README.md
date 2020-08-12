### TODO

There is a problem in ligand extraction where a chain has multiple ligand with separate sequence numbers.<br/>
Ligand must be split by (chain, sequence number) identity.

For the time being, please use Jupyter notebooks.

### LOG

```
ver 0.1.20  SminaSelfDocking.ipynb      - small changes
ver 0.1.19  cbiprep/smina.py            - small changes
                                        - test_smina.ipynb no longer works, so deleted
ver 0.1.18  papers/                     - collection of publications
ver 0.1.17  SminaSelfDocking.ipynb      - newly added
ver 0.1.16  cbiprep/matrix.py           - newly added
            DataSet.ipynb               - newly added
            cbiprep/pdbatoms.py         - remove NMR Hs and altLoc
ver 0.1.15                              - fixed in general
ver 0.1.14  test_smina.ipynb            - smina function to go through all PDBs
|                                       - RMSD is written out
|                                       - Cleaned
|           cbiprep/smina.py            - newly added
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
