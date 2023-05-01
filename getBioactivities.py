__author__ = "Cristian Orgaz Portero"
__copyright__ = "Copyright (C) 2023 Cristian Orgaz Portero"
__license__ = "Public Domain"
__version__ = "1.0"
# Usamos las librerias panda, requests y numpy
import pandas as pd
import requests
from datetime import datetime
# leemos el fichero y lo guardamos en la variable df
df=pd.read_csv('file_with_cids.csv',sep=";")
arr_bioactivity=[]
cids_with_error=[]
cids_without_bioactivity=[]
baseURL="https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=json&query={%22select%22:%22*%22,%22collection%22:%22bioactivity%22,%22where%22:{%22ands%22:[{%22cid%22:%22{cid}%22}]},%22start%22:{start},%22limit%22:{limit}}"

def get_bioactivities_by_cid(cid,start,limit, arr_results, arr_errors, arr_no_bioactivity):
    response=None
    total_count=0
    try:
        response = requests.get(url = baseURL.replace("{cid}",str(cid)).replace("{start}",str(start)).replace("{limit}",str(limit)))
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
    if(len(arr_errors)>0):
        arr_errors_copied=arr_errors
        arr_errors=[]
        get_bioactivities(arr_errors_copied, start, limit, arr_results, arr_errors, arr_no_bioactivity)

def results_to_csv(arr_results):
    pd.concat(arr_results).to_csv("bioactivities_by_cid.csv", sep=';',index=False)

start_time = datetime.now()
get_bioactivities(df["cid"], 1, 10000, arr_bioactivity, cids_with_error, cids_without_bioactivity)
end_time=datetime.now()
print('Duration: {}'.format(end_time - start_time))
