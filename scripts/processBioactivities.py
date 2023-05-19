__author__ = "Cristian Orgaz Portero"
__copyright__ = "Copyright (C) 2023 Cristian Orgaz Portero"
__license__ = "Public Domain"
__version__ = "1.0"
# Usamos las librerias panda, requests y numpy
import pandas as pd

df_bioactivity_guts=pd.read_csv('../results/bioactivities_by_cid_guts.csv', sep=";")
df_bioactivity_drugs=pd.read_csv('../results/bioactivities_by_cid_drugs.csv', sep=";")
df_bioactivity_filtered_guts=df_bioactivity_guts.copy()
df_bioactivity_filtered_drugs=df_bioactivity_drugs.copy()
df_bioactivity_filtered_guts.dropna(subset="repacxn", inplace=True)
df_bioactivity_filtered_drugs.dropna(subset="repacxn", inplace=True)
df_bioactivity_filtered_my_activity_guts=df_bioactivity_filtered_guts.copy()
df_bioactivity_filtered_my_activity_drugs=df_bioactivity_filtered_drugs.copy()
df_cids_guts=pd.read_csv('../results/guts_cids.csv', sep=";")
df_cids_drugs=pd.read_csv('../results/guts_cids.csv', sep=";")
dict_cids_guts={}
dict_cids_drugs={}
arr_of_results_guts=[]
arr_of_results_drugs=[]
def check_activities_values(activity, acvalue):
    match activity:
        case "Inactive":
            return "Inactive"
        case "Inconclusive":
            return "Inconclusive"
        case "Active":
            return "Active"
        case _:
            if(acvalue<10):
                return "Active"
            else:
                return "Inactive"

def check_bioactivity(cid, list_to_check, list_to_check_filtered_by_repacxn, list_to_check_full_filtered):
    if(cid in list_to_check and cid in list_to_check_full_filtered):
        return "yes"
    elif(cid in list_to_check and cid not in list_to_check_filtered_by_repacxn):
        return "no repacxn"
    else:
        return "no"

def filter_unique_values(cid, protacxn, my_activity, dict_resume):
    if(dict_resume.get(cid)==None):
        if(my_activity=="Active"):
            dict_resume[cid]={str(protacxn):1}
        else:
            dict_resume[cid]={str(protacxn):-1}
    elif(dict_resume[cid].get(protacxn)==None):
        if(my_activity=="Active"):
            dict_resume[cid][str(protacxn)]=1
        else:
            dict_resume[cid][str(protacxn)]=-1
    else:
        if(my_activity=="Active"):
            dict_resume[cid][protacxn]+=1
        else:
            dict_resume[cid][protacxn]+=(-1)
    return dict_resume

def obtain_filtered_results(arr, dict_cids):
    for i in dict_cids:
        for j in dict_cids[i]:
            if(dict_cids[i][j]>0):
                arr.append([i,j,"Active"])
            else:
                arr.append([i,j,"Inactive"])
    return pd.DataFrame(arr,columns=["cid","protacxn", "my_activity"])

df_bioactivity_filtered_my_activity_guts["my_activity"]=df_bioactivity_filtered_my_activity_guts.apply(lambda row: check_activities_values(row["activity"], row["acvalue"]), axis=1)
df_bioactivity_filtered_my_activity_guts.drop(df_bioactivity_filtered_my_activity_guts[df_bioactivity_filtered_my_activity_guts["my_activity"]=="Inconclusive"].index,inplace=True)
df_bioactivity_filtered_my_activity_guts.apply(lambda row:filter_unique_values(row["cid"],row["protacxn"],row["my_activity"],dict_cids_guts),axis=1)
result_df_guts=obtain_filtered_results(arr_of_results_guts, dict_cids_guts)
result_df_guts.to_csv('../results/filtered_bioactivity_result_guts.csv',sep=";",index=False)
df_cids_guts["bioactivity"]=df_cids_guts.apply(lambda row:check_bioactivity(row["cid"], pd.unique(df_bioactivity_guts["cid"]), pd.unique(df_bioactivity_filtered_guts["cid"]),pd.unique(result_df_guts[result_df_guts["my_activity"]=="Active"]["cid"])), axis = 1)
df_cids_guts.to_csv("../results/guts_comps_cids_bioactivity.csv", sep=';',index=False)

df_bioactivity_filtered_my_activity_drugs["my_activity"]=df_bioactivity_filtered_my_activity_drugs.apply(lambda row: check_activities_values(row["activity"], row["acvalue"]), axis=1)
df_bioactivity_filtered_my_activity_drugs.drop(df_bioactivity_filtered_my_activity_drugs[df_bioactivity_filtered_my_activity_drugs["my_activity"]=="Inconclusive"].index,inplace=True)
df_bioactivity_filtered_my_activity_drugs.apply(lambda row:filter_unique_values(row["cid"],row["protacxn"],row["my_activity"],dict_cids_drugs),axis=1)
result_df_drugs=obtain_filtered_results(arr_of_results_drugs, dict_cids_drugs)
result_df_drugs.to_csv('../results/filtered_bioactivity_result_drugs.csv',sep=";",index=False)
result_df_drugs["bioactivity"]=result_df_drugs.apply(lambda row:check_bioactivity(row["cid"], pd.unique(df_bioactivity_drugs["cid"]), pd.unique(df_bioactivity_filtered_drugs["cid"]),pd.unique(result_df_drugs[result_df_drugs["my_activity"]=="Active"]["cid"])), axis = 1)
result_df_drugs.to_csv("../results/drugs_comps_cids_bioactivity.csv", sep=';',index=False)