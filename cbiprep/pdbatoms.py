import numpy as np
import gzip

class PDBAtom(dict):
    def __init__(self, kwargs):
        self.update(kwargs)

    def __str__(self):
        return f"{self['record']}{self['serial']:5d}" \
               f" {self['name']:>4}{self['altLoc']}{self['resName']:<3} {self['chainID']}" \
               f"{self['resSeq']:>4}{self['iCode']}" \
               f"   {self['x']:8.3f}{self['y']:8.3f}{self['z']:8.3f}{self['occupancy']:6.2f}" \
               f"{self['tempFactor']:6.2f}" \
               f"          {self['element']:>2}{self['charge']:2}"
    
    def __repr__(self):
        return self.__str__()

class PDBAtoms(list):
    AMINO_ACIDS = 'ALA ARG ASN ASP ASX CYS GLU GLN GLX GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split()

    def __init__(self, pdb=None):
        self.pdb_file = None
        list.__init__(self)
        if pdb:
            self.read_pdb(pdb)

    def read_pdb(self, pdb):
        self.pdb_file = pdb
        openf = gzip.open if pdb.endswith('.gz') else open
        for line in openf(pdb, 'rt'):
            line = line.rstrip()
            if line[:6] in ['ATOM  ', 'HETATM']:
                kw = dict()
                kw['record'] = line[:6]
                kw['serial'] = int(line[6:11])
                kw['name'] = line[12:16]
                kw['altLoc'] = line[16]
                kw['resName'] = line[17:20]
                kw['chainID'] = line[21]
                kw['resSeq'] = int(line[22:26])
                kw['iCode'] = line[26]
                kw['x'] = float(line[30:38])
                kw['y'] = float(line[38:46])
                kw['z'] = float(line[46:54])
                kw['occupancy'] = float(line[54:60])
                kw['tempFactor'] = float(line[60:66])
                kw['element'] = line[76:78]
                kw['charge'] = line[78:80]
                if kw['altLoc'] not in [' ', 'A']:
                    continue
                if kw['element'].strip() == 'H':
                    continue
                self.append(PDBAtom(kw))
    
    def __str__(self):
        s = []
        for atom in self:
            s.append(str(atom))
        return '\n'.join(s)
    
    def __repr__(self):
        return self.__str__()
    
    def __add__(self, other):
        this = PDBAtoms()
        this += self
        this += other
        return this

    def __iadd__(self, other):
        self.extend(other)
        return self

    def GetResFragment(self, seqs, fix_carboxylate=False):
        """
        written by Iino
        # Needs review
        """
        chain_id = seqs[0][0]
        res_seqs = []
        for seq in seqs:
            res_seqs.append(seq[1])
        next_res_ns = []
        last_res = []
        frags = PDBAtoms()
        for a in self:
            if a['chainID'] == chain_id:
                if a['resSeq'] in res_seqs:
                    frags.append(a)
                if a['resSeq'] == max(res_seqs):
                    last_res.append(a)
                if a['element'].strip() == 'N' and a['resSeq'] == max(res_seqs) + 1:
                    next_res_ns.append(a)
        if fix_carboxylate and len(next_res_ns) > 0:
            for n in next_res_ns:
                n_cord = np.array([n['x'],n['y'],n['z']])
                for a in last_res:
                    a_cord = np.array([a['x'],a['y'],a['z']])
                    if np.abs(np.linalg.norm(n_cord - a_cord)) < 1.5:
                        n['name'] = '  O'
                        n['element'] = ' O'
                        n['resName'] = a['resName']
                        n['resSeq'] = a['resSeq']
                        frags.append(n)
                        break
        return frags

    def get_ligand_names(self):
        names = set()
        for atom in self:
            resName = atom['resName']
            if resName == 'HOH':
                continue
            if resName in self.AMINO_ACIDS:
                continue
            names.add(resName)
        return sorted(names)
    
    def get_ligand(self, name):
        info = set()
        for atom in self:
            if atom['resName'] == name:
                info.add((atom['chainID'], atom['resSeq']))
        if len(info) == 0:
            return PDBAtoms()
        info = sorted(info)
        chainID, resSeq = info[0][0], info[0][1]
        ligand_atoms = PDBAtoms()
        for atom in self:
            if atom['resName'] == name and atom['chainID'] == chainID and atom['resSeq'] == resSeq:
                ligand_atoms.append(atom)
        return ligand_atoms
    
    def get_protein(self):
        protein_atoms = PDBAtoms()
        for atom in self:
            if atom['resName'] in self.AMINO_ACIDS:
                protein_atoms.append(atom)
        return protein_atoms

    def get_relevant_protein(self, ligand_atoms, thres=5.0):
        protein_atoms = self.get_protein()
        distance_matrix = protein_atoms.get_distance_matrix(ligand_atoms)
        distance_mins = distance_matrix.min(axis=1)
        relevant_chain = set()
        for i in range(len(distance_mins)):
            if distance_mins[i] < thres:
                p = protein_atoms[i]
                relevant_chain.add(p['chainID'])
        relevant_atoms = PDBAtoms()
        for atom in protein_atoms:
            if atom['chainID'] in relevant_chain:
                relevant_atoms.append(atom)
        return relevant_atoms

    def get_pocket(self, ligand_atoms, thres=5.0):
        protein_atoms = self.get_protein()
        distance_matrix = protein_atoms.get_distance_matrix(ligand_atoms)
        distance_mins = distance_matrix.min(axis=1)
        pocket_res = set()
        for i in range(len(distance_mins)):
            if distance_mins[i] < thres:
                p = protein_atoms[i]
                pocket_res.add((p['chainID'], p['resSeq']))
        pocket_atoms = PDBAtoms()
        for atom in protein_atoms:
            if (atom['chainID'], atom['resSeq']) in pocket_res:
                pocket_atoms.append(atom)
        return pocket_atoms

    def get_distance_matrix(self, other):
        coords1 = np.array([(a['x'], a['y'], a['z']) for a in self])
        coords2 = np.array([(a['x'], a['y'], a['z']) for a in other])
        n1 = coords1.shape[0]
        n2 = coords2.shape[0]
        c1 = np.tile(coords1, n2).reshape(n1, n2, 3)
        c2 = np.tile(coords2.T, n1).T.reshape(n1, n2, 3)
        return np.linalg.norm(c1-c2, axis=2)

    def get_distance_based_adjacency_matrix(self, other, thres=5.0, diagzero=False):
        dmat = self.get_distance_matrix(other)
        x = (dmat < thres).astype(int)
        if diagzero:
            x -= np.eye(x.shape[0]).astype(int)
        return x
