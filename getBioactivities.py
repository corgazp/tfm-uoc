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
df_bioactivity=None
arr_bioactivity=[]
cids_with_error=[]
cids_without_bioactivity=[]
baseURL="https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=json&query={%22select%22:%22*%22,%22collection%22:%22bioactivity%22,%22where%22:{%22ands%22:[{%22cid%22:%22{cid}%22}]},%22start%22:{start},%22limit%22:10000}"

def getBioactivitiesByCID(CID,start,limit,arr_results,arr_errors, arr_no_bioactivity):
    response=None
    totalCount=0
    try:
        response = requests.get(url = baseURL.replace("{cid}",str(CID)).replace("{start}",str(start)))
        if(response.status_code==200):
            totalCount=response.json()["SDQOutputSet"][0]["totalCount"]
            start=start+limit
            if(totalCount>0):
                arr_results.append(pd.DataFrame(response.json()["SDQOutputSet"][0]["rows"]))
                if(totalCount>start):
                    getBioactivitiesByCID(CID,start,limit,arr_results, arr_errors, arr_no_bioactivity)
            else:
                arr_no_bioactivity.append(CID)
        else:
            arr_errors.append(CID)
    except:
        arr_errors.append(CID)
        if(response):
            print(f'API Resquest ERROR {response.status_code}')
        else:
            response="Error while execution"
            print(f'API Resquest ERROR: {response}')

def getBioactivities(CIDs,start,limit,arr_results,arr_errors, arr_no_bioactivity):
    for i in CIDs:
        getBioactivitiesByCID(i,start,limit,arr_results,arr_errors, arr_no_bioactivity)
    resultsToCSV(arr_results)
    if(len(cids_with_error)>0):
        getBioactivities(cids_with_error,start,limit,arr_results, arr_errors, arr_no_bioactivity)

def checkbioactivity(CID, list_to_check):
    if(CID in list_to_check):
        return "yes"
    else:
        return "no"

def resultsToCSV(arr_results):
    pd.concat(arr_results).dropna(subset="targetname").to_csv("bioactivities_by_cid.csv", sep=';',index=False)

def resultsToDataframe(arr_results):
    return pd.concat(arr_results).dropna(subset="targetname")

start_time = datetime.now()
getBioactivities(df["cid"],1,1000,arr_bioactivity,cids_with_error,cids_without_bioactivity)
df_bioactivity=resultsToDataframe(arr_bioactivity)
df["bioactivity"]=df.apply(lambda row:checkbioactivity(row["cid"], pd.unique(df_bioactivity["cid"])), axis = 1)
df.to_csv("gut_comps_cids_bioactivity.csv", sep=';',index=False)
end_time=datetime.now()
print('Duration: {}'.format(end_time - start_time))
