import os, sys
import rdkit.Chem as Chem 
from tqdm import tqdm 

froot = os.path.dirname(__file__)
rs_fpath = os.path.join(froot, 'rs.csv')

with open(rs_fpath, 'r') as fid:
	for i, line in tqdm(enumerate(fid)):
		if i == 0:
			continue
		smi = line.split(',')[0]
		mol = Chem.MolFromSmiles(smi)
		if any(len(a.GetBonds()) > 4 for a in mol.GetAtoms()):
			print(f'### WARNING: NBS > 4 IN {smi}')