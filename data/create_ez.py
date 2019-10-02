import os, sys
import rdkit.Chem as Chem 
from rdkit.Chem.rdchem import BondStereo

from tqdm import tqdm 

'''
Prepare a dataset of molecules where there is a single cis/trans double bond
and the objective is to predict E (0) versus Z (1) for the overall molecule at
that center according to CIP rules
'''

froot = os.path.dirname(__file__)
ez_fpath = os.path.join(froot, 'ez.csv')

mols = {}

# def atom_could_be_tetra(a):
#   '''Doesn't account for symmetry of side chains'''
#     return not (a.GetDegree() < 3 or (a.GetDegree() == 3 and 'H' not in a.GetSmarts()))

for name in os.listdir(froot):
    fpath = os.path.join(froot, name)
    if fpath == ez_fpath:
        continue
    if '.py' in fpath:
        continue
    print(f'processing {name}')

    with open(fpath, 'r') as fid:
        for i, line in tqdm(enumerate(fid)):
            if i == 0:
                continue
            smi = line.split(',')[0]
            if '.' in smi:
                continue
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            Chem.AssignStereochemistry(mol)
            num_cistrans = 0
            for b in mol.GetBonds():
                if b.GetStereo() in [BondStereo.STEREOE, BondStereo.STEREOZ]:
                    num_cistrans += 1
                    mol_ez = b.GetStereo()
            # Only keep if exactly one tetrahedral center
            if num_cistrans == 1:
                mols[smi] = mol_ez

print('Found {} E/Z molecules'.format(len(mols)))
mols_aug = {}
for smi, ez in mols.items():
    # Convert EZ tag to 01
    ez_01 = 0 if ez == BondStereo.STEREOE else 1
    # Record
    mols_aug[smi] = ez_01

with open(ez_fpath, 'w') as fid:
    print('smiles,EZ_as_01', file=fid)
    for smi, ez_01 in mols_aug.items():
        print(f'{smi},{ez_01}', file=fid)
