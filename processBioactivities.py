__author__ = "Cristian Orgaz Portero"
__copyright__ = "Copyright (C) 2023 Cristian Orgaz Portero"
__license__ = "Public Domain"
__version__ = "1.0"
# Usamos las librerias panda, requests y numpy
import pandas as pd

df_bioactivity=pd.read_csv('bioactivities_by_cid.csv', sep=";")
df_bioactivity_filtered=df_bioactivity.copy().dropna(subset="repacxn")
df_bioactivity_filtered_my_activity=df_bioactivity_filtered.copy()
df_cids=pd.read_csv('file_with_cids.csv', sep=";")
dict_cids={}
arr_of_results=[]
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

def filter_unique_values(cid, repacxn, my_activity, dict_resume):
    if(dict_resume.get(cid)==None):
        if(my_activity=="Active"):
            dict_resume[cid]={str(repacxn):1}
        else:
            dict_resume[cid]={str(repacxn):-1}
    elif(dict_resume[cid].get(repacxn)==None):
        if(my_activity=="Active"):
            dict_resume[cid][str(repacxn)]=1
        else:
            dict_resume[cid][str(repacxn)]=-1
    else:
        if(my_activity=="Active"):
            dict_resume[cid][repacxn]+=1
        else:
            dict_resume[cid][repacxn]+=(-1)
    return dict_resume

def obtain_filtered_results(arr):
    for i in dict_cids:
        for j in dict_cids[i]:
            if(dict_cids[i][j]>0):
                arr.append([i,j,"Active"])
    return pd.DataFrame(arr,columns=["cid","repacxn", "my_activity"])

df_bioactivity_filtered_my_activity["my_activity"]=df_bioactivity_filtered_my_activity.apply(lambda row: check_activities_values(row["activity"], row["acvalue"]), axis=1)
df_bioactivity_filtered_my_activity.drop(df_bioactivity_filtered_my_activity[df_bioactivity_filtered_my_activity["my_activity"]=="Inconclusive"].index,inplace=True)
df_bioactivity_filtered_my_activity.apply(lambda row:filter_unique_values(row["cid"],row["repacxn"],row["my_activity"],dict_cids),axis=1)
result_df=obtain_filtered_results(arr_of_results)
result_df.to_csv('filtered_bioactivity_result.csv',sep=";",index=False)
df_cids["bioactivity"]=df_cids.apply(lambda row:check_bioactivity(row["cid"], pd.unique(df_bioactivity["cid"]), pd.unique(df_bioactivity_filtered["cid"]),pd.unique(result_df["cid"])), axis = 1)
df_cids.to_csv("gut_comps_cids_bioactivity.csv", sep=';',index=False)