import pandas as pd
from rdkit.Chem import Descriptors
from rdkit import Chem
df_guts=pd.read_csv('../results/guts_smiles.csv', sep=";", encoding= 'unicode_escape')
df_drugs=pd.read_csv('../results/drugs_smiles.csv', sep=";", encoding= 'unicode_escape')

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

df_full=pd.concat([df_drugs,df_guts]).drop_duplicates()

def df2sdf(fname, mblist, idlist, idfield, addid = True):
    nmols = len(mblist)
    f = open(fname,"w+")
    for i in range(nmols):
        if i > 0:
            f.write("$$$$\n")
        if addid is True:
            f.write(f'\t{idlist[i]}\n{mblist[i][1:]}') 
        else:
            f.write(mblist[i])
        f.write(f'{idfield}\n')
        f.write(f'{idlist[i]}\n\n')
    f.write("$$$$\n")
    f.close()

df2sdf("df_full.sdf", [Chem.MolToMolBlock(x) for x in df_full[df_full.rb<9].mol], df_full[df_full.rb<9].name.values.tolist(), '> <full_dataset >', addid = True)