#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  4 07:54:06 2022

@author: gonzalo
"""


import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools as pt
import matplotlib.pyplot as plt
from rdkit import RDLogger
from datetime import datetime
from rdkit.Chem.Draw import IPythonConsole
from scipy import stats as ss
from rdkit.Chem import Descriptors
import random
import math
import seaborn as sns

def sdf2df(fsdf, idname, encoding = 'utf-8', pattern_start = '$$$$', pattern_end = "M  END", idfield = "> <DATABASE_ID>", getmol = False, added_id = False):
    """ Generates a dataframe by reading line by line an sdf file """
    list_molblocks = []
    ids = []
    if added_id:
        ids2 = []
        idname2 = idname + "2"
    parsing = False
    nextfdid = False
    molblock = ""
    flist = open(fsdf, encoding = encoding).readlines()
    i = 0
    for line in flist:
        if parsing:
            if line.startswith(pattern_end):
                molblock = molblock + line.rstrip()
                list_molblocks.append(molblock)
                parsing = False
            else:
                if added_id and molblock == "":
                    ids2.append(line.rstrip())
                molblock += line
        if nextfdid:
            ids.append(line.rstrip())
            nextfdid = False
        if line.startswith(idfield):
            nextfdid = True 
        if (line.startswith(pattern_start)) or (i == 0):
            molblock = ""
            parsing = True
        i = i + 1
    if added_id:
        df = pd.DataFrame({idname: ids, "molblock": list_molblocks, idname2: ids2})
    else:
        df = pd.DataFrame({idname: ids, "molblock": list_molblocks})
    print("Read ", df.shape[0], " molecules")
    df['mol'] = df.molblock.apply(lambda x: Chem.MolFromMolBlock(x))
    nfail = df.mol.isnull().sum()
    if nfail > 0:
        print("Unable to generate molecule for " + str(df.mol.isnull().sum()) + " molblocks")
    df = df[~df.mol.isnull()]
    print("Imported ", df.shape[0], " molecules")
    if ~getmol:
        df.drop("mol", axis = 1, inplace = True)
    return df



def df2sdf(fname, mblist, idlist, idfield, addid = True):
    nmols = len(mblist)
    f = open(fname,"w+")
    for i in range(nmols):
        if i > 0:
            f.write("$$$$\n")
        if addid is True:
            f.write("\t" + idlist[i] + "\n" + mblist[i][1:] + "\n") 
        else:
            f.write(mblist[i] + "\n")
        f.write(idfield + "\n")
        f.write(idlist[i]+"\n\n")
    f.write("$$$$\n")
    f.close()


def df2sdf2(fname, mblist, idlist, idfield, addid = True):
    nmols = len(mblist)
    f = open(fname,"w+")
    for i in range(nmols):
        if i > 0:
            f.write("$$$$\n")
        if addid is True:
            mb = mblist[i]
            mb = "\n".join(mb.split("\n")[1:])
            f.write("\t" + idlist[i] + "\n" + mb[1:] + "\n") 
        else:
            f.write(mblist[i] + "\n")
        f.write(idfield + "\n")
        f.write(idlist[i]+"\n\n")
    f.write("$$$$\n")
    f.close()



    
def sdf2energy(fsdf, encoding = 'utf-8', pattern_end = '$$$$', pattern_end_mb = "M  END", enfield = "> <minimizedAffinity>"):
    """ Generates a dataframe by reading line by line an sdf file """
    ens = []
    ids = []
    parsing = False
    molblock = ""
    nextid = False
    nexten = False
    flist = open(fsdf, encoding = encoding).readlines()
    i = 0
    for line in flist:
        line = flist[i]
        if (i == 0):
            ids.append(line.rstrip())
            parsing = True
        if (line.startswith(pattern_end)):
            molblock = ""
            nextid = True
            parsing = True
        if parsing:
            if (nextid == True and molblock != ""):
                ids.append(line.rstrip())
                molblock = molblock + line
                nextid = False
            else: 
                if line.startswith(pattern_end_mb):
                    parsing = False
                    molblock = molblock + line.rstrip()
                else:
                    molblock = molblock + line
        else:
            if(line.startswith(enfield)):
                nexten = True
            if (nexten) and (not(line.startswith(enfield))):
                ens.append(line.rstrip())
                nexten = False
        i = i + 1
    df = pd.DataFrame({"id": ids, "energy": ens})
    print("Read ", df.shape[0], " molecules")
    return df


def is_pareto_efficient_simple(costs):
    """
    Find the pareto-efficient points
    :param costs: An (n_points, n_costs) array
    :return: A (n_points, ) boolean array, indicating whether each point is Pareto efficient
    """
    is_efficient = np.ones(costs.shape[0], dtype = bool)
    for i, c in enumerate(costs):
        if is_efficient[i]:
            is_efficient[is_efficient] = np.any(costs[is_efficient]>c, axis=1)  # Keep any point with a lower cost
            is_efficient[i] = True  # And keep self
    return is_efficient

def pareto_rank(df, id, cost_ids):
    """ 
    Pareto rank
    """
    nrows = df.shape[0]
    df["prank"] = [np.nan]*nrows 
    df_aux = df
    score = 1
    
    while df.prank.isnull().sum() > 0:
        costs = df_aux[cost_ids].to_numpy()
        effs = df_aux.id[is_pareto_efficient_simple(costs)]
        df.loc[df.id.isin(effs),"prank"] = score
        df_aux = df_aux[~df_aux.id.isin(effs)]
        score = score + 1
        df_aux
    return df

def sdf_extraction(filename, pattern_start='$', pattern_end='>'):
    molblock = """"""
    list_molblocks = []
    flist = open(filename, encoding = 'utf-8').readlines()
    parsing = True
    for line in flist:
        if line.startswith(pattern_start):
            parsing = True
            continue
        if line.startswith(pattern_end):
            parsing = False
            if len(molblock) > 0: 
                list_molblocks.append(molblock)
                molblock = """"""
                continue
        if parsing:
            molblock += line
    return list_molblocks
    

######################################
###### GENERATE INPUT SDF FILE
######################################

# path = "/results"

# fdf = sdf2df(fsdf = path + "/foodb_clean.sdf", idname = "fid", encoding = 'utf-8', pattern_start = '$$$$', 
#             pattern_end = "M  END", idfield = "> <foodb_id>", getmol = False)
# fdf.shape
# fdf.columns
# fdf["mol"] = fdf.molblock.apply(Chem.MolFromMolBlock)
# fdf["rb"] = fdf.mol.apply(Chem.Lipinski.NumRotatableBonds)
# fdf09 = fdf[fdf.rb <= 9]
# fdf09.shape
# df2sdf2(path + "/f.sdf", mblist = fdf09.molblock.to_list(), idlist = fdf09.fid.to_list(), idfield = "> <foodb_id>")


######################################
###### SELECT TOP HITS
######################################

path = '/results'
df1 = sdf2energy("../results/df_full_multi01.sdf")
df2 = sdf2energy("../results/df_full_multi02.sdf")
df3 = sdf2energy("../results/df_full_multi03.sdf")

df1["energy"] = df1.energy.apply(lambda x: float(x))
df2["energy"] = df2.energy.apply(lambda x: float(x))
df3["energy"] = df3.energy.apply(lambda x: float(x))
df3.dtypes
df3.columns

df1 = df1.groupby(["id"]).apply(lambda x: min(x.energy)).reset_index()
df2 = df2.groupby(["id"]).apply(lambda x: min(x.energy)).reset_index()
df3 = df3.groupby(["id"]).apply(lambda x: min(x.energy)).reset_index()
# df1=df1[df1[0]<=12]
# df2=df2[df2[0]<=12]
# df3=df3[df3[0]<=12]
df1.shape
df2.shape
df3.shape

df = pd.merge(df1, df2, on = ["id"], how = "outer")
df = pd.merge(df, df3, on = ["id"], how = "outer")
df.columns = ["id","en1","en2","en3"]
df.shape
df.columns
arr_result=pd.plotting.scatter_matrix(df[["en1","en2","en3"]], alpha=0.9, hist_kwds={'bins':300}, s = 0.1)
for i in range(3):
    for j in range(3):
        arr_result[i,j].set_xlim(-14,14)
        if i!=j:
            arr_result[i,j].set_ylim(-14,14)

df["en1m"] = df.en1.apply(lambda x: -1*x)
df["en2m"] = df.en2.apply(lambda x: -1*x) 
df["en3m"] = df.en3.apply(lambda x: -1*x)
df = pareto_rank(df, "id", ["en1m","en2m","en3m"])
df["rankcent"] = df.prank * df.prank.min() / df.prank.max()*100
df.prank.max()
df.prank.min()

df[["id","prank","rankcent"]].to_csv("../results/metabos-sel.csv", sep = ";", index = False)





