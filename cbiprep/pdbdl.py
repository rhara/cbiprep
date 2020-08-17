import pandas as pd
from urllib import request
from multiprocessing import Pool, freeze_support
import os, sys


def getPDB(pdb_code, dest='.'):
    os.makedirs(dest, exist_ok=True)
    local_name = f'{dest}/{pdb_code}.pdb.gz'
    remote_url = f'http://files.rcsb.org/download/{pdb_code}.pdb.gz'
    if os.path.exists(local_name):
        return False
    request.urlretrieve(remote_url, local_name)
    return True

def __gen(df, dest):
    count = 0
    for i in df.index:
        count += 1
        r = df.loc[i]
        yield dict(count=count, pdb=r['pdb'], dest=dest)

def __worker(par):
    downloaded = getPDB(par['pdb'], dest=par['dest'])
    par['downloaded'] = downloaded
    return par

def download_df_PDBs(df, dest='.', nprocs=0):
    if nprocs == 0:
        nprocs = os.cpu_count()

    pool = Pool(nprocs)
    ndownloaded = 0
    for par in pool.imap_unordered(__worker, __gen(df, dest)):
        if par['downloaded']:
            ndownloaded += 1
            print(par)
    return ndownloaded


if __name__ == '__main__':
    freeze_support()

    df = pd.read_pickle('index_2019.pkl.gz')
    DF = df[(df['type'] == 'Kd') & (df['lig_ok'] == True) & (df['refined'] == True)]
    print(DF)

    download_df_PDBs(DF, 'pdb')
