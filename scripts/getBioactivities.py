__author__ = "Cristian Orgaz Portero"
__copyright__ = "Copyright (C) 2023 Cristian Orgaz Portero"
__license__ = "Public Domain"
__version__ = "1.0"
# Usamos las librerias panda, requests y numpy
import pandas as pd
import requests
from datetime import datetime
import os
# leemos el fichero y lo guardamos en la variable df_guts
df_guts=pd.read_csv('../results/guts_cids.csv',sep=";", encoding= 'unicode_escape')
df_drugs=pd.read_csv('../results/drugs_cids.csv',sep=";", encoding= 'unicode_escape')
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

def get_bioactivities(cids, start, limit, arr_results, arr_errors, arr_no_bioactivity, filename):
    for i in cids:
        get_bioactivities_by_cid(i, start, limit, arr_results, arr_errors, arr_no_bioactivity)
    results_to_csv(arr_results, filename)
    if(len(arr_errors)>0):
        arr_errors_copied=arr_errors
        arr_errors=[]
        get_bioactivities(arr_errors_copied, start, limit, arr_results, arr_errors, arr_no_bioactivity)

def results_to_csv(arr_results, filename):
    path = "../results/bioactivity"
    isExist = os.path.exists(path)
    if not isExist:
       os.makedirs(path)
    print("Directory results has been created!")
    pd.concat(arr_results).to_csv(f'../results/bioactivity/{filename}.csv', sep=';',index=False)

start_time = datetime.now()
arr_bioactivity_guts=[]
cids_with_error_guts=[]
cids_without_bioactivity_guts=[]
get_bioactivities(df_guts["cid"], 1, 10000, arr_bioactivity_guts, cids_with_error_guts, cids_without_bioactivity_guts, 'bioactivities_by_cid_guts')
end_time=datetime.now()
print('Duration guts: {}'.format(end_time - start_time))
start_time = datetime.now()
arr_bioactivity_drugs=[]
cids_with_error_drugs=[]
cids_without_bioactivity_drugs=[]
get_bioactivities(df_drugs["cid"], 1, 10000, arr_bioactivity_drugs, cids_with_error_drugs, cids_without_bioactivity_drugs, 'bioactivities_by_cid_drugs')
end_time=datetime.now()
print('Duration drugs: {}'.format(end_time - start_time))
