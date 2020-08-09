import pandas as pd
import os, re

class PDBBindIndex:
    AMINO_ACIDS = 'ALA ARG ASN ASP ASX CYS GLU GLN GLX GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split()
    NUCLEOTIDES = 'AMP ADP ATP GMP GDP GTP TMP TDP TTP CMP CDP CTP UMP UDP UTP'.split()

    COMMENT_PAT = re.compile('^#.*$')
    TYPE_PAT = re.compile('^(Ki|Kd|IC50)[~><=].*$')
    LIGNAME_PAT = re.compile('^[A-Z0-9]{3}$')

    def __init__(self, fname):
        self.fname = fname
        self.df = self.__read()

    def is_ligand_ok(self, ligname):
        m = self.LIGNAME_PAT.match(ligname)
        if not m:
            return False
        if ligname in self.AMINO_ACIDS:
            return False
        if ligname in self.NUCLEOTIDES:
            return False
        return True

    def __read(self):
        data = []
        for line in open(self.fname):
            line = line.rstrip()
            if self.COMMENT_PAT.match(line):
                continue
            if line == '':
                continue
            items = line.split(maxsplit=7)
            pdb = items[0]
            year = int(items[2])
            pval = float(items[3])
            type = items[4]
            m = self.TYPE_PAT.match(type)
            if not m:
                raise Exception(type)
            type = m.group(1)
            ligname = items[7][1:-1]
            lig_ok = self.is_ligand_ok(ligname)
            data.append(dict(pdb=pdb, year=year, pval=pval, type=type, lig=ligname, lig_ok=lig_ok))
        df = pd.DataFrame(data)
        return df

    @staticmethod
    def read(general_index, refined_index):
        data = PDBBindIndex(general_index)
        data_refined = PDBBindIndex(refined_index)
        pdbs_refined = set(data_refined.df['pdb'])
        in_refined = []
        for i in data.df.index:
            r = data.df.iloc[i]
            in_refined.append(r['pdb'] in pdbs_refined)
        data.df['refined'] = in_refined
        return data.df


if __name__ == '__main__':
    general_index = os.path.dirname(__file__) + '/../data/INDEX_general_PL_data.2019'
    if not os.path.exists(general_index):
        general_index = os.path.dirname(__file__) + '/../../../../cbiprep/data/INDEX_general_PL_data.2019'
    refined_index = os.path.dirname(__file__) + '/../data/INDEX_refined_data.2019'
    if not os.path.exists(refined_index):
        refined_index = os.path.dirname(__file__) + '/../../../../cbiprep/data/INDEX_refined_data.2019'
    df = PDBBindIndex.read(general_index=general_index, refined_index=refined_index)
    print(df)
    df.to_pickle('index_2019.pkl.gz', compression='gzip')

