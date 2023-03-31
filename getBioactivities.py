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
df_bioactivity_filtered=None
arr_bioactivity=[]
cids_with_error=[]
cids_without_bioactivity=[]
baseURL="https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=json&query={%22select%22:%22*%22,%22collection%22:%22bioactivity%22,%22where%22:{%22ands%22:[{%22cid%22:%22{cid}%22}]},%22start%22:{start},%22limit%22:10000}"

def get_bioactivities_by_cid(cid,start,limit, arr_results, arr_errors, arr_no_bioactivity):
    response=None
    total_count=0
    try:
        response = requests.get(url = baseURL.replace("{cid}",str(cid)).replace("{start}",str(start)))
        if(response.status_code==200):
            total_count=response.json()["SDQOutputSet"][0]["totalCount"]
            start=start+limit
            if(total_count>0):
                arr_results.append(pd.DataFrame(response.json()["SDQOutputSet"][0]["rows"]))
                if(total_count>start):
                    get_bioactivities_by_cid(cid, start, limit, arr_results, arr_errors, arr_no_bioactivity)
            else:
                arr_no_bioactivity.append(cid)
        else:
            arr_errors.append(cid)
    except:
        arr_errors.append(cid)
        if(response):
            print(f'API Resquest ERROR {response.status_code}')
        else:
            response="Error while execution"
            print(f'API Resquest ERROR: {response}')

def get_bioactivities(cids, start, limit, arr_results, arr_errors, arr_no_bioactivity):
    for i in cids:
        get_bioactivities_by_cid(i, start, limit, arr_results, arr_errors, arr_no_bioactivity)
    results_to_csv(arr_results)
    results_filtered_to_csv(arr_results)
    if(len(arr_errors)>0):
        arr_errors_copied=arr_errors
        arr_errors=[]
        get_bioactivities(arr_errors_copied, start, limit, arr_results, arr_errors, arr_no_bioactivity)

def check_bioactivity(cid, list_to_check, list_to_check_filtered):
    if(cid in list_to_check and cid in list_to_check_filtered):
        return "yes"
    elif(cid in list_to_check and cid not in list_to_check_filtered):
        return "no repacxn"
    else:
        return "no"

def results_to_csv(arr_results):
    pd.concat(arr_results).to_csv("bioactivities_by_cid.csv", sep=';',index=False)
def results_filtered_to_csv(arr_results):
    pd.concat(arr_results).dropna(subset="repacxn").to_csv("bioactivities_by_cid_filtered.csv", sep=';',index=False)

def results_to_dataframe(arr_results):
    return pd.concat(arr_results)
def results_filtered_to_dataframe(arr_results):
    return pd.concat(arr_results).dropna(subset="repacxn")

start_time = datetime.now()
get_bioactivities(df["cid"], 1, 1000, arr_bioactivity, cids_with_error, cids_without_bioactivity)
df_bioactivity=results_to_dataframe(arr_bioactivity)
df_bioactivity_filtered=results_filtered_to_dataframe(arr_bioactivity)
df["bioactivity"]=df.apply(lambda row:check_bioactivity(row["cid"], pd.unique(df_bioactivity["cid"]), pd.unique(df_bioactivity_filtered["cid"])), axis = 1)
df.to_csv("gut_comps_cids_bioactivity.csv", sep=';',index=False)
end_time=datetime.now()
print('Duration: {}'.format(end_time - start_time))
