__author__ = "Cristian Orgaz Portero"
__copyright__ = "Copyright (C) 2023 Cristian Orgaz Portero"
__license__ = "Public Domain"
__version__ = "1.0"
import pandas as pd
import requests
import numpy as np
from datetime import datetime
# leemos los ficheros y los guardamos en las variables df_guts y df_drugs
df_guts=pd.read_csv('../results/guts_cids.csv',sep=";")
df_drugs=pd.read_csv('../results/drugs_cids.csv',sep=";")
baseURL="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES/JSON"

def get_smiles_by_cid(cid):
    response=None
    try:
        response = requests.get(url = baseURL.replace("{cid}",str(cid).replace(".0",'')))
        if(response.status_code==200):
            return response.json()["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
        else:
            return np.nan
    except:
        if(response):
            print(f'API Resquest ERROR {response.status_code}')
            response = requests.get(url = baseURL.replace("{cid}",str(cid).replace(".0",'')))
            if(response.status_code==200):
                return response.json()["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
            else:
                return np.nan
        else:
            response="Error while execution"
            print(f'API Resquest ERROR: {response}')
            return np.nan

def check_data_na(SMILES, cid):
    if(SMILES==np.nan):
        return get_smiles_by_cid(cid)
    else:
        return SMILES

def generate_SEA_file(SMILES, CID, file):
    file.write(f"{SMILES} {CID}\n")

start_time = datetime.now()
df_guts["canonical_smiles"]=df_guts.apply(lambda row:get_smiles_by_cid(row["cid"]),axis=1)
df_guts["canonical_smiles"]=df_guts.apply(lambda row:check_data_na(row["canonical_smiles"], row["cid"]),axis=1)
df_guts.to_csv("../results/guts_smiles.csv", sep=';',index=False)
file_SEA=open("../results/input_guts.txt","w")
file_SEA.write("SMILES CID\n")
df_guts.apply(lambda row:generate_SEA_file(row["canonical_smiles"], row["cid"], file_SEA),axis=1)
file_SEA.close()
end_time=datetime.now()
print('Duration guts: {}'.format(end_time - start_time))
start_time = datetime.now()
df_drugs["canonical_smiles"]=df_drugs.apply(lambda row:get_smiles_by_cid(row["cid"]),axis=1)
df_drugs["canonical_smiles"]=df_drugs.apply(lambda row:check_data_na(row["canonical_smiles"], row["cid"]),axis=1)
df_drugs.to_csv("../results/drugs_smiles.csv", sep=';',index=False)
file_SEA=open("../results/input_drugs.txt","w")
file_SEA.write("SMILES CID\n")
df_drugs.apply(lambda row:generate_SEA_file(row["canonical_smiles"], row["cid"], file_SEA),axis=1)
file_SEA.close()
end_time=datetime.now()
print('Duration drugs: {}'.format(end_time - start_time))