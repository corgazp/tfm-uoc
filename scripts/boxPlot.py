import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from rdkit.Chem import Descriptors
from rdkit import Chem
df=pd.read_csv('file_with_smiles.csv', sep=";")
ccl=pd.read_csv('guts_ccl.csv', sep=";").ccl
df["ccl"]=ccl
# Violin plots y boxplots de propiedades fisicoquímicas
df['mol'] = df.canonical_smiles.apply(lambda x: Chem.MolFromSmiles(x))
df["tpsa"] = df.mol.apply(Descriptors.TPSA)
df["logp"] = df.mol.apply(Descriptors.MolLogP)
df["rb"] = df.mol.apply(Chem.Lipinski.NumRotatableBonds) 
df["hbd"] = df.mol.apply(Chem.Lipinski.NumHDonors)
df["hba"] = df.mol.apply(Chem.Lipinski.NumHAcceptors)
df["mw"] = df.mol.apply(Descriptors.ExactMolWt)
df["qed"] = df.mol.apply(Descriptors.qed)
df["nring"] = df.mol.apply(Chem.Lipinski.RingCount)
df["naring"] = df.mol.apply(Chem.Lipinski.NumAromaticRings)
df["fsp3"] = df.mol.apply(Chem.Lipinski.FractionCSP3)
df["set"]=df.ccl.apply(lambda x: "GL" if x=="Glycerolipids" else "NoGL" )
pchpros = ["tpsa","logp","rb","hbd","hba","mw","qed","nring","naring","fsp3"]
for i in range(len(pchpros)):
    var = pchpros[i]
    sz = 15
    bw = 0.8
    ax = sns.boxplot(x= "set", y= var, data = df)
    ax.set_xlabel("", fontsize = sz)
    ax.set_ylabel("", fontsize = sz)
    ax.set_title(var.upper(), fontsize = sz*1.2)
    ax.tick_params(labelsize=10)
    plt.show()
