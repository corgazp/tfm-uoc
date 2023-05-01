__author__ = "Cristian Orgaz Portero"
__copyright__ = "Copyright (C) 2023 Cristian Orgaz Portero"
__license__ = "Public Domain"
__version__ = "1.0"
import pandas as pd
from rdkit.Chem.Scaffolds.MurckoScaffold import GetScaffoldForMol, MakeScaffoldGeneric
from rdkit.Chem import MolFromSmiles, MolToSmiles
df=pd.read_csv('file_with_smiles.csv', sep=";")

def get_smiles_from_scaffold(SMILES):
    molecule=MolFromSmiles(SMILES)
    return MolToSmiles(GetScaffoldForMol(molecule))
    
def get_smiles_from_generic_scaffold(SMILES):
    molecule=MolFromSmiles(SMILES)
    return MolToSmiles(MakeScaffoldGeneric(molecule))

df["smiles_from_scaffold"]=df.apply(lambda row: get_smiles_from_scaffold(row["canonical_smiles"]), axis=1)
df["smiles_from_generic_scaffold"]=df.apply(lambda row: get_smiles_from_generic_scaffold(row["canonical_smiles"]), axis=1)
df.to_csv('file_with_scaffolds_smiles.csv', sep=";" , index=False)