__author__ = "Cristian Orgaz Portero"
__copyright__ = "Copyright (C) 2023 Cristian Orgaz Portero"
__license__ = "Public Domain"
__version__ = "1.0"
import pandas as pd
from rdkit.Chem import Draw
from rdkit import Chem
import os
path = "../results/images"
isExist = os.path.exists(path)
if not isExist:
   os.makedirs(path)
   print("Directory results/images has been created!")
df_sc=pd.read_csv('../results/guts_scaffolds_smiles.csv',sep=";")
dd={"smiles_generic_scaffolds":df_sc.smiles_from_generic_scaffold.value_counts().index.tolist(),"count":df_sc.smiles_from_generic_scaffold.value_counts().tolist()}
dat_count_sc=pd.DataFrame(data=dd)
dat_count_sc_subset=dat_count_sc[0:10]
for i in dat_count_sc_subset.smiles_generic_scaffolds:
    m=Chem.MolFromSmiles(i)
    img = Draw.MolToImage(m)
    img.save(f'../results/images/{i}.png')