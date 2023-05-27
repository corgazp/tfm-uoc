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
   
df_sc_guts=pd.read_csv('../results/guts_scaffolds_smiles.csv',sep=";")
df_sc_drugs=pd.read_csv('../results/drugs_scaffolds_smiles.csv',sep=";")

data_guts_generic_scaffold={"smiles_generic_scaffolds":df_sc_guts.smiles_from_generic_scaffold.value_counts().index.tolist(),"count":df_sc_guts.smiles_from_generic_scaffold.value_counts().tolist()}
data_drugs_generic_scaffold={"smiles_generic_scaffolds":df_sc_drugs.smiles_from_generic_scaffold.value_counts().index.tolist(),"count":df_sc_drugs.smiles_from_generic_scaffold.value_counts().tolist()}
data_guts_scaffold={"smiles_scaffolds":df_sc_guts.smiles_from_scaffold.value_counts().index.tolist(),"count":df_sc_guts.smiles_from_scaffold.value_counts().tolist()}
data_drugs_scaffold={"smiles_scaffolds":df_sc_drugs.smiles_from_scaffold.value_counts().index.tolist(),"count":df_sc_drugs.smiles_from_scaffold.value_counts().tolist()}

dat_count_sc_generic_guts=pd.DataFrame(data=data_guts_generic_scaffold)
dat_count_sc_generic_guts_subset=dat_count_sc_generic_guts[0:10]
dat_count_sc_generic_drugs=pd.DataFrame(data=data_drugs_generic_scaffold)
dat_count_sc_generic_drugs_subset=dat_count_sc_generic_drugs[0:10]

dat_count_sc_guts=pd.DataFrame(data=data_guts_scaffold)
dat_count_sc_guts_subset=dat_count_sc_guts[0:10]
dat_count_sc_drugs=pd.DataFrame(data=data_drugs_scaffold)
dat_count_sc_drugs_subset=dat_count_sc_drugs[0:10]

def get_molecules_from_scaffolds(array_scaffolds, image_prefix):
    for i in range(0,len(array_scaffolds)):
        smiles=array_scaffolds[i]
        m=Chem.MolFromSmiles(smiles)
        img = Draw.MolToImage(m)
        img.save(f'../results/images/{image_prefix}_{smiles}_{i}.png')

get_molecules_from_scaffolds(dat_count_sc_generic_guts_subset.smiles_generic_scaffolds, "gut_generic")
get_molecules_from_scaffolds(dat_count_sc_generic_drugs_subset.smiles_generic_scaffolds, "drugs_generic")
get_molecules_from_scaffolds(dat_count_sc_guts_subset.smiles_scaffolds, "gut")
get_molecules_from_scaffolds(dat_count_sc_drugs_subset.smiles_scaffolds, "drugs")