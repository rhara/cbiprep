from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from cbiprep.ligand_expo import LigandExpo
from cbiprep import pdbatoms
import pandas as pd
import os, sys
from multiprocessing import Pool, freeze_support


def __gen(df):
    expo_dic = LigandExpo()
    count = 0
    for idx in df.index:
        count += 1
        r = df.loc[idx]
        lig_name = r['lig']
        pdb_code = r['pdb']
        smi = expo_dic.get(lig_name)
        par = dict(idx=idx, count=count, lig_name=lig_name, smi=smi, pdb_code=pdb_code)
        yield par


def __worker(par):
    idx = par['idx']
    count = par['count']
    lig_name = par['lig_name']
    pdb_code = par['pdb_code']
    smi = par['smi']
    par['error'] = None

    lig_oname = f'work/{pdb_code}/{pdb_code}_ligand_{lig_name}.sdf'
    pocket_oname = f'work/{pdb_code}/{pdb_code}_pocket.pdb'

    if os.path.exists(lig_oname) and os.path.exists(pocket_oname):
        par['error'] = 'already exists'
        return par

    if not par['smi']:
        par['error'] = 'no entry in ligand expo'
        return par

    try:
        pdb = pdbatoms.PDBAtoms(f'pdb/{pdb_code}.pdb.gz')
    except Exception as e:
        par['error'] = 'failed read pdb'
        return par

    RDLogger.DisableLog('rdApp.warning')
    try:
        ligand_atoms = pdb.get_ligand(lig_name)
        rdkit_mol = Chem.MolFromPDBBlock(str(ligand_atoms))
        template = Chem.MolFromSmiles(smi)
        lig_mol = AllChem.AssignBondOrdersFromTemplate(template, rdkit_mol)
    except Exception as e:
        par['error'] = 'failed create ligand'
        RDLogger.EnableLog('rdApp.warning')
        return par
    RDLogger.EnableLog('rdApp.warning')

    try:
        pocket_atoms = pdb.get_pocket(ligand_atoms, thres=5.0)
    except Exception as e:
        par['error'] = 'fail create pocket'
        return par

    os.makedirs(f'work/{pdb_code}', exist_ok=True)
    lig_ofs = Chem.SDWriter(lig_oname)
    lig_ofs.write(lig_mol)
    lig_ofs.close()
    open(pocket_oname, 'wt').write(str(pocket_atoms)+'\n')

    return par


def process_ligpocket(df, nprocs=0):
    if nprocs == 0:
        nprocs = os.cpu_count()

    expo_dic = LigandExpo()
    __EXPO_DIC = expo_dic

    pool = Pool(nprocs)

    for par in pool.imap_unordered(__worker, __gen(df)):
        print(par, file=sys.stderr)


if __name__ == '__main__':
    freeze_support()

    df = pd.read_pickle('index_2019.pkl.gz')
    DF = df[(df['type'] == 'Kd') & (df['lig_ok'] == True) & (df['refined'] == True)]

    expo_dic = LigandExpo()

    pool = Pool(os.cpu_count())

    for par in pool.imap_unordered(__worker, __gen(DF)):
        print(par)
