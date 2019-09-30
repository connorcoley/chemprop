import os, sys
import rdkit.Chem as Chem 
from rdkit.Chem.rdchem import ChiralType
from tqdm import tqdm 

'''
Prepare a dataset of molecules where there is a single tetrahedral center
and the objective is to predict R (0) versus S (1) for the overall molecule at
that center according to CIP rules
'''

froot = os.path.dirname(__file__)
rs_fpath = os.path.join(froot, 'rs.csv')

mols = {}

# def atom_could_be_tetra(a):
#   '''Doesn't account for symmetry of side chains'''
#     return not (a.GetDegree() < 3 or (a.GetDegree() == 3 and 'H' not in a.GetSmarts()))

for name in os.listdir(froot):
    fpath = os.path.join(froot, name)
    if fpath == rs_fpath:
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
            num_tetra = 0
            for a in mol.GetAtoms():
                if a.HasProp('_CIPCode'):
                    num_tetra += 1
                    mol_rs = a.GetProp('_CIPCode')
            # Only keep if exactly one tetrahedral center
            if num_tetra == 1:
                mols[smi] = mol_rs

print('Found {} chiral molecules'.format(len(mols)))
mols_aug = {}
for smi, rs in mols.items():
    # Convert RS tag to 01
    rs_01 = 0 if rs == 'R' else 1
    # Add enantiomer for good measure
    smi_flip = smi.replace('@@', '@') if '@@' in smi else smi.replace('@', '@@')
    rs_01_flip = 1 - rs_01
    # Record both
    mols_aug[smi] = rs_01
    mols_aug[smi_flip] = rs_01_flip
print('Added enantiomers to make {} total'.format(len(mols_aug)))

with open(rs_fpath, 'w') as fid:
    print('smiles,RS_as_01', file=fid)
    for smi, rs_01 in mols_aug.items():
        print(f'{smi},{rs_01}', file=fid)
