### update the local repository

```
ver 0.1.14 | test_smina.ipynb | smina function to go through all PDBs
           |                  | RMSD is written out
           |                  | Cleaned
           | cbiprep/smina.py | newly added
```

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
