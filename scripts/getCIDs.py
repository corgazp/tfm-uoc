__author__ = "Cristian Orgaz Portero"
__copyright__ = "Copyright (C) 2023 Cristian Orgaz Portero"
__license__ = "Public Domain"
__version__ = "1.0"
# Usamos las librerias panda, requests y numpy
import pandas as pd
import requests
import numpy as np
from datetime import datetime

import os
path = "../results"
isExist = os.path.exists(path)
if not isExist:
   os.makedirs(path)
   print("Directory results has been created!")
# leemos el fichero y lo guardamos en la variable df
df_guts=pd.read_csv('../data/guts.csv',sep=";")
df_drugs=pd.read_csv('../data/drugs.csv',sep=";")
baseURL="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/cids/JSON"
# URL para conseguir CID a partir de HMDBID
altURL="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/xref/RegistryID/{registryID}/cids/JSON"
# Creamos la función get_cid_by_inchi
def get_cid_by_inchi(inchi):
    PARAMS = {'inchi':inchi}
    response=None
    try:
        response = requests.get(url = baseURL, params = PARAMS)
        if(response.status_code==200):
            return response.json()["IdentifierList"]["CID"][0]
        else:
            return np.nan
    except:
        if(response):
            print(f'API Resquest ERROR {response.status_code}')
            response = requests.get(url = baseURL, params = PARAMS)
            if(response.status_code==200):
                return response.json()["IdentifierList"]["CID"][0]
            else:
                return np.nan
        else:
            response="Error while execution"
            print(f'API Resquest ERROR: {response}')
            return np.nan
# Función para conseguir CID a partir de HMDBID        
def get_cid_by_hmdbid(hmdbid):
    response=None
    try:
        response = requests.get(url = altURL.replace("{registryID}",hmdbid))
        if(response.status_code==200):
            return response.json()["IdentifierList"]["CID"][0]
        else:
            return np.nan
    except:
        if(response):
            print(f'API Resquest ERROR {response.status_code}')
            response = requests.get(url = altURL.replace("{registryID}", hmdbid))
            if(response.status_code==200):
                return response.json()["IdentifierList"]["CID"][0]
            else:
                return np.nan
        else:
            response="Error while execution"
            print(f'API Resquest ERROR: {response}')
            return np.nan
# Comprobamos la existencia de valores NA o valores igual a cero
def check_data_na(cid, hmdbid):
    if(np.isnan(cid) or cid == 0):
        return get_cid_by_hmdbid(hmdbid.split("|")[1])
    else:
        return cid


start_time = datetime.now()
# Recorremos la columna inchi y realizamos a tanta peticiones a la api como valores hay en la columna
df_guts['cid'] = df_guts.apply(lambda row:get_cid_by_inchi(row['inchi']), axis = 1)
df_guts['cid'] = df_guts.apply(lambda row:check_data_na(row['cid'], row['hmdb_id']), axis = 1)
df_guts.to_csv("../results/guts_cids.csv", sep=';',index=False)
end_time=datetime.now()
print('Duration guts: {}'.format(end_time - start_time))

start_time = datetime.now()
# Recorremos la columna inchi y realizamos a tanta peticiones a la api como valores hay en la columna
df_drugs['cid'] = df_drugs.apply(lambda row:get_cid_by_inchi(row['inchi']), axis = 1)
df_drugs['cid'] = df_drugs.apply(lambda row:check_data_na(row['cid'], row['hmdb_id']), axis = 1)
df_drugs.to_csv("../results/drugs_cids.csv", sep=';',index=False)
end_time=datetime.now()
print('Duration drugs: {}'.format(end_time - start_time))