import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from rdkit.Chem import Descriptors
from rdkit import Chem
df_guts=pd.read_csv('../results/guts_smiles.csv', sep=";", encoding= 'unicode_escape')
df_drugs=pd.read_csv('../results/drugs_smiles.csv', sep=";", encoding= 'unicode_escape')
# Violin plots y boxplots de propiedades fisicoquímicas
df_drugs['mol'] = df_drugs.canonical_smiles.apply(lambda x: Chem.MolFromSmiles(x))
df_drugs["tpsa"] = df_drugs.mol.apply(Descriptors.TPSA)
df_drugs["logp"] = df_drugs.mol.apply(Descriptors.MolLogP)
df_drugs["rb"] = df_drugs.mol.apply(Chem.Lipinski.NumRotatableBonds) 
df_drugs["hbd"] = df_drugs.mol.apply(Chem.Lipinski.NumHDonors)
df_drugs["hba"] = df_drugs.mol.apply(Chem.Lipinski.NumHAcceptors)
df_drugs["mw"] = df_drugs.mol.apply(Descriptors.ExactMolWt)
df_drugs["qed"] = df_drugs.mol.apply(Descriptors.qed)
df_drugs["nring"] = df_drugs.mol.apply(Chem.Lipinski.RingCount)
df_drugs["naring"] = df_drugs.mol.apply(Chem.Lipinski.NumAromaticRings)
df_drugs["fsp3"] = df_drugs.mol.apply(Chem.Lipinski.FractionCSP3)
df_drugs["set"]="DB"
df_guts['mol'] = df_guts.canonical_smiles.apply(lambda x: Chem.MolFromSmiles(x))
df_guts["tpsa"] = df_guts.mol.apply(Descriptors.TPSA)
df_guts["logp"] = df_guts.mol.apply(Descriptors.MolLogP)
df_guts["rb"] = df_guts.mol.apply(Chem.Lipinski.NumRotatableBonds) 
df_guts["hbd"] = df_guts.mol.apply(Chem.Lipinski.NumHDonors)
df_guts["hba"] = df_guts.mol.apply(Chem.Lipinski.NumHAcceptors)
df_guts["mw"] = df_guts.mol.apply(Descriptors.ExactMolWt)
df_guts["qed"] = df_guts.mol.apply(Descriptors.qed)
df_guts["nring"] = df_guts.mol.apply(Chem.Lipinski.RingCount)
df_guts["naring"] = df_guts.mol.apply(Chem.Lipinski.NumAromaticRings)
df_guts["fsp3"] = df_guts.mol.apply(Chem.Lipinski.FractionCSP3)
df_guts["set"]=df_guts.ccl.apply(lambda x: "G (GL)" if x=="Glycerolipids" else "G (NoGL)" )
pchpros = ["tpsa","logp","rb","hbd","hba","mw","qed","nring","naring","fsp3"]

df_full=pd.concat([df_drugs,df_guts])
for i in range(len(pchpros)):
    var = pchpros[i]
    sz = 15
    bw = 0.8
    ax = sns.boxplot(x= "set", y= var, data = df_full, showfliers = False)
    #ax = sns.violinplot(x= "group", y= var, data = alldf_drugs, bw = 0.8)
    ax.set_xlabel("", fontsize = sz)
    ax.set_ylabel("", fontsize = sz)
    ax.set_title(var.upper(), fontsize = sz*1.2)
    ax.tick_params(labelsize=10)
    plt.show()
    
for i in range(len(pchpros)):
    var = pchpros[i]
    sz = 15
    bw = 0.8
    ax = sns.violinplot(x= "set", y= var, data = df_full, bw = bw)
    #ax = sns.violinplot(x= "group", y= var, data = alldf, bw = 0.8)
    ax.set_xlabel("", fontsize = sz)
    ax.set_ylabel("", fontsize = sz)
    ax.set_title(var.upper(), fontsize = sz*1.2)
    ax.tick_params(labelsize=10)
    plt.show()