import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm

df_drugs=pd.read_csv('../results/drugs_cids.csv', sep=";", encoding='windows-1252')
df_guts=pd.read_csv('../results/guts_cids.csv', sep=";")
df_drugs_active=pd.read_csv('../results/bioactivity/filtered_bioactivity_result_drugs.csv', sep=";")
df_drugs_active=df_drugs_active[df_drugs_active["my_activity"]=="Active"]
df_guts_active=pd.read_csv('../results/bioactivity/filtered_bioactivity_result_guts.csv', sep=";")
df_guts_active=df_guts_active[df_guts_active["my_activity"]=="Active"]
df_drugs_target_class=pd.read_csv('../results/bioactivity/target_class_pubchem_drugs.csv', sep=";")
df_guts_target_class=pd.read_csv('../results/bioactivity/target_class_pubchem_guts.csv', sep=";")
df_drugs_uniprots=pd.read_csv('../results/bioactivity/pubchem_targets_uniprot_id_filtered_drugs.csv', sep=";")
df_guts_uniprots=pd.read_csv('../results/bioactivity/pubchem_targets_uniprot_id_filtered_guts.csv', sep=";")



df_drugs_sea=pd.read_csv('../results/drugs_smiles.csv', sep=";", encoding='windows-1252')
df_guts_sea=pd.read_csv('../results/guts_smiles.csv', sep=";", encoding='windows-1252')
df_drugs_sea_predicted=pd.read_csv('../results/SEA/SEA_predicted_activities_drugs_filtered_duplicates.csv', sep=";")
df_guts_sea_predicted=pd.read_csv('../results/SEA/SEA_predicted_activities_guts_filtered_duplicates.csv', sep=";")
df_drugs_sea_uniprots=pd.read_csv('../results/SEA/sea_target_class_drugs.csv', sep=";")
df_guts_sea_uniprots=pd.read_csv('../results/SEA/sea_target_class_guts.csv', sep=";")


comp_order = ['Organoheterocyclic compounds',
'Glycerolipids',
'Benzenoids',
'Organic acids and derivatives',
'Organic oxygen compounds',
'Other',
'Steroids and steroid derivatives',
'Fatty Acyls',
'Phenylpropanoids and polyketides',
'Prenol lipids',
'Glycerophospholipids',
'Organic nitrogen compounds',
'Nucleosides, nucleotides, and analogues',
'Organosulfur compounds',
'Hydrocarbons',
'Sphingolipids',
'Endocannabinoids']
santos2_order=[
    '7TM1',
    'Nuclear receptor',
    'VGIC',
    'Oxidoreductase',
    'LGIC',
    'Electrochemical transporter',
    'Kinase',
    'Protease',
    'Hydrolase',
    'Transferase',
    '7TM2',
    'Lyase',
    'Isomerase',
    'Phosphodiesterase',
    'Cytochrome P450',
    'Epigenetic regulator',
    '7TM3',
    'Phosphotase',
    'Other'
    ]
def add_ccl_to_active_compounds(cid, df_ccl):
    return df_ccl[df_ccl["cid"]==str(cid)].ccl.values[0] if df_ccl[df_ccl["cid"]==str(cid)].ccl.values else df_ccl[df_ccl["cid"]==cid].ccl.values[0]

def add_protacxn_to_target(acc,df_protacxn):
    if len(df_protacxn[df_protacxn["acc"]==acc].protacxn.values)>0:
        return list(df_protacxn[df_protacxn["acc"]==acc].protacxn.values)
    else:
        return ""

def add_target_ccl_to_df(val, df_protacxn):
    target_class_found=False
    for i in range(0,len(df_protacxn["protacxn"])):
        if len(df_protacxn["protacxn"][i])>0:
            if val in df_protacxn["protacxn"][i]:
                target_class_found==True
                return df_protacxn["tcl"][i]
    if target_class_found==False:
        if len(list(df_protacxn[df_protacxn["accession"]==val].tcl.values))>0:
            return list(df_protacxn[df_protacxn["accession"]==val].tcl.values)[0]
        else:
            return np.nan

def add_ccl_to_predicted_compounds(SMILES, df_ccl):
    return df_ccl[df_ccl["canonical_smiles"]==SMILES].ccl.values[0] if len(df_ccl[df_ccl["canonical_smiles"]==SMILES].ccl.values)>0 else "pepe"


def add_target_ccl_to_df_sea(val, df_uniprots):
    target_class_found=False
    for i in range(0,len(df_uniprots["component_synonym"])):
        if len(df_uniprots["component_synonym"])>0:
            if val==df_uniprots["component_synonym"][i]:
                target_class_found==True
                return df_uniprots["tcl"][i]
    if target_class_found==False:
        if len(list(df_uniprots[df_uniprots["component_synonym"]==val].tcl.values))>0:
            return list(df_uniprots[df_uniprots["component_synonym"]==val].tcl.values)[0]
        else:
            return np.nan

def get_matrix_interaction(df, ccl_values, tcl_values,result):
    for i in range(0,len(ccl_values)):
        row=[]
        for j in range(0,len(tcl_values)):
           row.append(len(list(df.query(f'tcl== "{tcl_values[j]}" and ccl == "{ccl_values[i]}"').tcl)))
        result.append(row)
    return result

arr_int_matrix_drugs=[]
arr_int_matrix_guts=[]
arr_int_matrix_drugs_sea=[]
arr_int_matrix_guts_sea=[]

df_drugs_active["ccl"]=df_drugs_active.apply(lambda row:add_ccl_to_active_compounds(row["cid"],df_drugs),axis=1)
df_guts_active["ccl"]=df_guts_active.apply(lambda row:add_ccl_to_active_compounds(row["cid"],df_guts),axis=1)
df_drugs_target_class["protacxn"]=df_drugs_target_class.apply(lambda row: add_protacxn_to_target(row["accession"], df_drugs_uniprots),axis=1)
df_guts_target_class["protacxn"]=df_guts_target_class.apply(lambda row: add_protacxn_to_target(row["accession"], df_guts_uniprots),axis=1)
df_drugs_active["tcl"]=df_drugs_active.apply(lambda row: add_target_ccl_to_df(row["protacxn"], df_drugs_target_class),axis=1)
df_guts_active["tcl"]=df_guts_active.apply(lambda row: add_target_ccl_to_df(row["protacxn"], df_guts_target_class),axis=1)
get_matrix_interaction(df_drugs_active, comp_order, santos2_order, arr_int_matrix_drugs)
get_matrix_interaction(df_guts_active, comp_order, santos2_order, arr_int_matrix_guts)



df_drugs_sea_predicted["ccl"]=df_drugs_sea_predicted.apply(lambda row: add_ccl_to_predicted_compounds(row["Query Smiles"],df_drugs_sea),axis=1)
df_guts_sea_predicted["ccl"]=df_guts_sea_predicted.apply(lambda row: add_ccl_to_predicted_compounds(row["Query Smiles"],df_guts_sea),axis=1)
df_drugs_sea_predicted["tcl"]=df_drugs_sea_predicted.apply(lambda row: add_target_ccl_to_df_sea(row["Target ID"], df_drugs_sea_uniprots),axis=1)
df_guts_sea_predicted["tcl"]=df_guts_sea_predicted.apply(lambda row: add_target_ccl_to_df_sea(row["Target ID"], df_guts_sea_uniprots),axis=1)
get_matrix_interaction(df_drugs_sea_predicted, comp_order, santos2_order, arr_int_matrix_drugs_sea)
get_matrix_interaction(df_guts_sea_predicted, comp_order, santos2_order, arr_int_matrix_guts_sea)

# ## Heatmap de matriz de interacciones por clase química y clase de diana

def generate_heatmap(df,title, v_min,v_max):
    f = plt.figure(figsize=(14, 14))                                                                                                                                    
    plt.matshow(df, fignum=f.number, cmap=plt.cm.get_cmap("Blues"), norm=LogNorm(vmin = v_min, vmax =  v_max))
    plt.yticks(range(df.shape[0]), comp_order, fontsize=19)
    plt.xticks(range(df.shape[1]), santos2_order, fontsize=19, rotation=90)
    plt.ylim(df.shape[0]-0.5, -0.5) # Stupid fix to set correct label positions in y axis
    cb = plt.colorbar(fraction=0.038, pad=0.14)
    cb.ax.tick_params(labelsize=14)
    #plt.rcParams['axes.titley'] = 1.5
    plt.title(title, fontsize=27, y = 1.6)
    cc = df_guts_active[["cid","ccl"]].drop_duplicates().ccl.value_counts().reindex(index = comp_order, fill_value = 0)
    df_cc = df_guts_active[pd.notna(df_guts_active.tcl)][["cid","ccl"]].drop_duplicates()["ccl"].value_counts().reindex(index = comp_order, fill_value = 0)
    # for i in range(len(df_cc)):
    #     tx = str(df_cc.iloc[i]) + "(" + '{:0.2f}'.format(df_cc.iloc[i] / cc.iloc[i]*100) + ")"
    #     plt.text(23-0.3, i, tx, ha='left', va='center',fontsize = 15)
    # df_tc = df_guts_active[pd.notna(df_guts_active.tcl)][["protacxn","tcl"]].drop_duplicates()["tcl"].value_counts().reindex(index = santos2_order, fill_value = 0)
    # for i in range(len(df_tc)):
    #     plt.text(i, 18-0.3, df_tc.iloc[i], ha='center', va='top',fontsize = 15, rotation = "vertical")
    #plt.title("FooDB Class vs Target Class: Adjusted Residuals (ChEMBL+SEA)", fontsize=27)
    for (i, j), z in np.ndenumerate(df):
        if z > 500:
            col = "white"
            fs = 12
            if z > 999:
                fs = 11
                if z > 9999:
                    fs = 10
        else:
            col = "black"
            fs = 15
        plt.text(j, i, '{:0.0f}'.format(z), ha='center', va='center', fontsize = fs, color = col)
    plt.show()
 
generate_heatmap(pd.DataFrame(arr_int_matrix_drugs), "Drugbank Compound Class vs Target Class (PubChem)",1,1782)
generate_heatmap(pd.DataFrame(arr_int_matrix_guts), "Gut Compound Class vs Target Class (PubChem)",1,29)
generate_heatmap(pd.DataFrame(arr_int_matrix_drugs_sea), "Drugbank Compound Class vs Target Class (SEA+TC)",1,152)
generate_heatmap(pd.DataFrame(arr_int_matrix_guts_sea), "Gut Compound Class vs Target Class (SEA+TC)",1,2258)