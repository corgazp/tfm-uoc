__author__ = "Cristian Orgaz Portero"
__copyright__ = "Copyright (C) 2023 Cristian Orgaz Portero"
__license__ = "Public Domain"
__version__ = "1.0"
import pandas as pd
from rdkit.Chem.Scaffolds.MurckoScaffold import GetScaffoldForMol, MakeScaffoldGeneric
from rdkit.Chem import MolFromSmiles, MolToSmiles
df_guts=pd.read_csv('../results/guts_smiles.csv', sep=";", encoding= 'unicode_escape')
df_drugs=pd.read_csv('../results/drugs_smiles.csv', sep=";" , encoding= 'unicode_escape')
def get_smiles_from_scaffold(SMILES):
    molecule=MolFromSmiles(SMILES)
    return MolToSmiles(GetScaffoldForMol(molecule))
    
def get_smiles_from_generic_scaffold(SMILES):
    molecule=MolFromSmiles(SMILES)
    return MolToSmiles(MakeScaffoldGeneric(molecule))

df_guts["smiles_from_scaffold"]=df_guts.apply(lambda row: get_smiles_from_scaffold(row["canonical_smiles"]), axis=1)
df_guts["smiles_from_generic_scaffold"]=df_guts.apply(lambda row: get_smiles_from_generic_scaffold(row["canonical_smiles"]), axis=1)
df_guts.to_csv('../results/guts_scaffolds_smiles.csv', sep=";" , index=False)

df_drugs["smiles_from_scaffold"]=df_drugs.apply(lambda row: get_smiles_from_scaffold(row["canonical_smiles"]), axis=1)
df_drugs["smiles_from_generic_scaffold"]=df_drugs.apply(lambda row: get_smiles_from_generic_scaffold(row["canonical_smiles"]), axis=1)
df_drugs.to_csv('../results/drugs_scaffolds_smiles.csv', sep=";" , index=False)