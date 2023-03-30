__author__ = "Cristian Orgaz Portero"
__copyright__ = "Copyright (C) 2023 Cristian Orgaz Portero"
__license__ = "Public Domain"
__version__ = "1.0"
# Usamos las librerias panda, requests y numpy
import pandas as pd
import requests
import numpy as np
from datetime import datetime
# leemos el fichero y lo guardamos en la variable df
df=pd.read_csv('gut_comps.csv',sep=";")
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
            return np.NA
    except:
        if(response):
            print(f'API Resquest ERROR {response.status_code}')
            response = requests.get(url = baseURL, params = PARAMS)
            if(response.status_code==200):
                return response.json()["IdentifierList"]["CID"][0]
            else:
                return np.NA
        else:
            response="Error while execution"
            print(f'API Resquest ERROR: {response}')
            return np.NA
# Función para conseguir CID a partir de HMDBID        
def get_cid_by_hmdbid(hmdbid):
    response=None
    try:
        response = requests.get(url = altURL.replace("{registryID}",hmdbid))
        if(response.status_code==200):
            return response.json()["IdentifierList"]["CID"][0]
        else:
            return np.NA
    except:
        if(response):
            print(f'API Resquest ERROR {response.status_code}')
            response = requests.get(url = altURL.replace("{registryID}", hmdbid))
            if(response.status_code==200):
                return response.json()["IdentifierList"]["CID"][0]
            else:
                return np.NA
        else:
            response="Error while execution"
            print(f'API Resquest ERROR: {response}')
            return np.NA
# Comprobamos la existencia de valores NA o valores igual a cero
def check_data_na(cid,inchi, hmdbid):
    if(np.isnan(cid)):
        return get_cid_by_inchi(inchi)
    elif(cid==0):
        return get_cid_by_hmdbid(hmdbid)
    else:
        return cid

start_time = datetime.now()
# Recorremos la columna inchi y realizamos a tanta peticiones a la api como valores hay en la columna
df['cid'] = df.apply(lambda row:get_cid_by_inchi(row['inchi']), axis = 1)
df['cid'] = df.apply(lambda row:check_data_na(row['cid'], row['inchi'], row['hmdb_id']), axis = 1)
df.to_csv("file_with_cids.csv", sep=';',index=False)
end_time=datetime.now()
print('Duration: {}'.format(end_time - start_time))