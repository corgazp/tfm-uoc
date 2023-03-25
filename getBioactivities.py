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
df=pd.read_csv('file_with_cids.csv',sep=";")
df_bioactivity=[]
baseURL="https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=json&query={%22select%22:%22*%22,%22collection%22:%22bioactivity%22,%22where%22:{%22ands%22:[{%22cid%22:%22{cid}%22}]},%22start%22:{start},%22limit%22:10000}"

def getBioactivitiesByCID(CID,start,limit,df_result):
    response=None
    totalCount=0
    try:
        response = requests.get(url = baseURL.replace("{cid}",str(CID)).replace("{start}",str(start)))
        if(response.status_code==200):
            totalCount=response.json()["SDQOutputSet"][0]["totalCount"]
            start=start+limit
            if(totalCount>0):
                df_result.append(pd.DataFrame(response.json()["SDQOutputSet"][0]["rows"]))
            if(totalCount>start):
                getBioactivitiesByCID(CID,start,limit,df_result)
        else:
            return np.NA
    except:
        if(response):
            print(f'API Resquest ERROR {response.status_code}')
            getBioactivitiesByCID(CID,start,limit,df_result)
        else:
            response="Error while execution"
            print(f'API Resquest ERROR: {response}')
            return np.NA

def getBioactivities(CIDs,start,limit,df_result):
    for i in CIDs:
        getBioactivitiesByCID(i,start,limit,df_result)
    pd.concat(df_result).dropna(subset="targetname").to_csv("bioactivities_by_cid.csv", sep=';',index=False)

start_time = datetime.now()
getBioactivities(df["cid"],1,1000,df_bioactivity)
end_time=datetime.now()
print('Duration: {}'.format(end_time - start_time))
