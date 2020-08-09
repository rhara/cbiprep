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
