from rdkit import Chem
from rdkit.Chem.rdMolAlign import GetBestRMS
from cbiprep.pdbatoms import PDBAtoms
import subprocess as sp
import os

def RunSmina(pdb_code, ligand_name, ncpu=0, workdir='work', num_modes=4, seed=0):
    if ncpu == 0:
        ncpu = os.cpu_count()
    ligand_sdf = f'{workdir}/{pdb_code}/{ligand_name}.sdf'
    for ligand_mol in Chem.SDMolSupplier(ligand_sdf):
        break
    pdb_atoms = PDBAtoms(f'pdb/{pdb_code}.pdb.gz')
    lig_atoms = pdb_atoms.get_ligand(ligand_name)
    protein_atoms = pdb_atoms.get_relevant_protein(lig_atoms, thres=4)
    
    protein_pdb = f'{workdir}/{pdb_code}/{pdb_code}_apo.pdb'
    open(protein_pdb, 'wt').write(str(protein_atoms))
    
    tleaprc_template = os.path.abspath(os.path.dirname(Chem.__file__) + '/../../../../../cbiprep/data/tleaprc.template')

    tleap_log = f'{workdir}/{pdb_code}/{pdb_code}_tleap.log'
    test = f'{workdir}/{pdb_code}/{pdb_code}_apo.pdb'
    save_mol2 = f'{workdir}/{pdb_code}/{pdb_code}_apo.mol2'
    tleaprc = f'{workdir}/{pdb_code}/{pdb_code}_tleaprc'
    cont = open(tleaprc_template, 'rt').read()
    cont = cont.replace('{{logFile}}', tleap_log).replace('{{test}}', test).replace('{{saveMol2}}', save_mol2)
    open(tleaprc, 'wt').write(cont)

    command = f'tleap -s -f {tleaprc} > /dev/null 2>&1'
    sp.call(command, shell=True)
    
    fixed_receptor = f'{workdir}/{pdb_code}/{pdb_code}_charged.mol2'
    command = f'obabel {save_mol2} -O {fixed_receptor} > /dev/null 2>&1'
    sp.call(command, shell=True)
    
    boxpars = sp.check_output(f'ligand_center {ligand_sdf}', shell=True).decode().strip()
    
    docked_sdf = f'{workdir}/{pdb_code}/{pdb_code}_ligand_{ligand_name}_docked.sdf'
    smina_log = f'work/{pdb_code}/{pdb_code}_ligand_{ligand_name}_smina.log'
    command = f'smina -r {fixed_receptor} -l {ligand_sdf} {boxpars} --cpu {ncpu} --num_modes {num_modes} --seed {seed} -o {docked_sdf} --log {smina_log}'
    print(command)
    sp.call(command, shell=True)

    return docked_sdf

def CheckRMSD(pdb_code, ligand_name, docked_sdf, workdir='work'):
    for xray_ligand in Chem.SDMolSupplier(f'{workdir}/{pdb_code}/{ligand_name}.sdf'):
        break
    rmsds = []
    for mol in Chem.SDMolSupplier(docked_sdf):
        rmsd = GetBestRMS(mol, xray_ligand)
        rmsds.append(round(rmsd, 3))
    return rmsds
        
