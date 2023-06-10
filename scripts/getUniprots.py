__author__ = "Cristian Orgaz Portero"
__copyright__ = "Copyright (C) 2023 Cristian Orgaz Portero"
__license__ = "Public Domain"
__version__ = "1.0"
# Usamos las librerias panda, requests y numpy
import pandas as pd
import requests
from datetime import datetime
import numpy as np
baseURL="https://pubchem.ncbi.nlm.nih.gov/assay/pcget.cgi?task=protein_similarseq&protacxn={protacxn}&start={start}&limit=10000&infmt=json&outfmt=json"
arr_protacxn_active_drugs=[]
arr_protacxn_inactive_drugs=[]
protacxns_with_error_active_drugs=[]
protacxns_with_error_inactive_drugs=[]
protacxn_without_uniprot_id_active_drugs=[]
protacxn_without_uniprot_id_inactive_drugs=[]
df_drugs=pd.read_csv('../results/bioactivity/filtered_bioactivity_result_drugs.csv', sep=";")
df_drugs_with_underscore=df_drugs[df_drugs["protacxn"].str.contains("_")]
df_drugs_with_underscore_active=df_drugs_with_underscore[df_drugs_with_underscore["my_activity"]=="Active"]
df_drugs_with_underscore_inactive=df_drugs_with_underscore[df_drugs_with_underscore["my_activity"]=="Inactive"]
df_drugs_len_upper_six=df_drugs[df_drugs["protacxn"].str.len()>6]
df_drugs_len_upper_six_active=df_drugs_len_upper_six[df_drugs_len_upper_six["my_activity"]=="Active"]
df_drugs_len_upper_six_inactive=df_drugs_len_upper_six[df_drugs_len_upper_six["my_activity"]=="Inactive"]
df_drugs_active_no_uniprot_id=pd.concat([df_drugs_len_upper_six_active,df_drugs_with_underscore_active]).drop_duplicates()
df_drugs_inactive_no_uniprot_id=pd.concat([df_drugs_len_upper_six_inactive,df_drugs_with_underscore_inactive]).drop_duplicates()

arr_protacxn_active_guts=[]
arr_protacxn_inactive_guts=[]
protacxns_with_error_active_guts=[]
protacxns_with_error_inactive_guts=[]
protacxn_without_uniprot_id_active_guts=[]
protacxn_without_uniprot_id_inactive_guts=[]
df_guts=pd.read_csv('../results/bioactivity/filtered_bioactivity_result_guts.csv', sep=";")
df_guts_with_underscore=df_guts[df_guts["protacxn"].str.contains("_")]
df_guts_with_underscore_active=df_guts_with_underscore[df_guts_with_underscore["my_activity"]=="Active"]
df_guts_with_underscore_inactive=df_guts_with_underscore[df_guts_with_underscore["my_activity"]=="Inactive"]
df_guts_len_upper_six=df_guts[df_guts["protacxn"].str.len()>6]
df_guts_len_upper_six_active=df_guts_len_upper_six[df_guts_len_upper_six["my_activity"]=="Active"]
df_guts_len_upper_six_inactive=df_guts_len_upper_six[df_guts_len_upper_six["my_activity"]=="Inactive"]
df_guts_active_no_uniprot_id=pd.concat([df_guts_len_upper_six_active,df_guts_with_underscore_active]).drop_duplicates()
df_guts_inactive_no_uniprot_id=pd.concat([df_guts_len_upper_six_inactive,df_drugs_with_underscore_inactive]).drop_duplicates()
def get_uniprot_id_by_protacxn(protacxn, start, limit, arr_results, arr_errors, arr_no_uniprot):
    response=None
    total_count=0
    try:
        response = requests.get(url = baseURL.replace("{protacxn}",str(protacxn)).replace("{start}",str(start)))
        if(response.status_code==200):
            if(response.json()["SDQOutputSet"][0]["status"]["code"]==0):
                total_count=response.json()["SDQOutputSet"][0]["totalCount"]
                start=start+limit
                if(total_count>0):
                    tmp=pd.DataFrame(response.json()["SDQOutputSet"][0]["rows"])
                    tmp["protacxn"]=protacxn
                    arr_results.append(tmp)
                    
                    if(total_count>start):
                        get_uniprot_id_by_protacxn(protacxn, start, limit, arr_results, arr_errors, arr_no_uniprot)
            else:
                arr_no_uniprot.append(protacxn)
        else:
            arr_errors.append(protacxn)
    except:
        arr_errors.append(protacxn)
        if(response):
            print(f'API Resquest ERROR {response.status_code}')
        else:
            response="Error while execution"
            print(f'API Resquest ERROR: {response}')

def get_unitprots(protacxns, start, limit, arr_results, arr_errors, arr_no_uniprot, filename):
    for i in pd.unique(protacxns):
        print(f'{np.where(pd.unique(protacxns)==i)[0]+1} of {len(pd.unique(protacxns))}')
        get_uniprot_id_by_protacxn(i, start, limit, arr_results, arr_errors, arr_no_uniprot)
    results_to_csv(arr_results,filename)
    if(len(arr_errors)>0):
        arr_errors_copied=arr_errors
        arr_errors=[]
        get_unitprots(arr_errors_copied, start, limit, arr_results, arr_errors, arr_no_uniprot)

def results_to_csv(arr_results, filename):
    result=pd.concat(arr_results)
    result[result.ident>=99].to_csv(f'{filename}.csv', sep=';',index=False)

start_time = datetime.now()
get_unitprots(df_drugs_active_no_uniprot_id["protacxn"], 1, 10000, arr_protacxn_active_drugs, protacxns_with_error_active_drugs, protacxn_without_uniprot_id_active_drugs, '../results/bioactivity/uniprots_by_protacxn_active_drugs')
end_time=datetime.now()
print('Duration for Active drugs: {}'.format(end_time - start_time))
start_time = datetime.now()
get_unitprots(df_drugs_inactive_no_uniprot_id["protacxn"], 1, 10000, arr_protacxn_inactive_drugs, protacxns_with_error_inactive_drugs, protacxn_without_uniprot_id_inactive_drugs, '../results/bioactivity/uniprots_by_protacxn_inactive_drugs')
end_time=datetime.now()
print('Duration for Inactive drugs: {}'.format(end_time - start_time))

start_time = datetime.now()
get_unitprots(df_guts_active_no_uniprot_id["protacxn"], 1, 10000, arr_protacxn_active_guts, protacxns_with_error_active_guts, protacxn_without_uniprot_id_active_guts, '../results/bioactivity/uniprots_by_protacxn_active_guts')
end_time=datetime.now()
print('Duration for Active guts: {}'.format(end_time - start_time))
start_time = datetime.now()
get_unitprots(df_guts_inactive_no_uniprot_id["protacxn"], 1, 10000, arr_protacxn_inactive_guts, protacxns_with_error_inactive_guts, protacxn_without_uniprot_id_inactive_guts, '../results/bioactivity/uniprots_by_protacxn_inactive_guts')
end_time=datetime.now()
print('Duration for Inactive guts: {}'.format(end_time - start_time))