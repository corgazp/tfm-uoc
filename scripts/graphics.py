import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from rdkit.Chem import Descriptors
df=pd.read_csv('guts_ccl.csv', sep=";")
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
'Saccharolipids',
'Endocannabinoids']

d={"ccl":df.ccl.value_counts().index.tolist(),"count":df.ccl.value_counts().tolist()}
dat_count=pd.DataFrame(data=d)
## Barplot de cuentas de compuestos por clase química
fig = plt.figure(figsize=(16,5))  
ax = sns.barplot(x = "ccl", y = "count", data = dat_count, orient = "v", order = comp_order)
ax.tick_params(axis='x', rotation=90)
delta = 130
sz = 16
ax.set_xlabel("Clase de compuesto", fontsize = sz*1.1)
ax.set_ylabel("Número de compuestos", fontsize = sz*1.1)
ax.set_title("Compuestos por clase química", fontsize = sz*1.3)
ax.set_yticklabels(ax.get_yticks(), size = sz*0.9)
ax.set_xticklabels([x.get_text().capitalize() for x in ax.get_xticklabels()], size = sz*0.9)
ax.tick_params(axis='x', rotation=90)
plt.legend(loc = "upper right", fontsize=sz)
plt.margins(y = 0.15, x = 0.01)
for p in ax.patches:
  # get the height of each bar
  height = p.get_height()
  # adding text to each bar
  ax.text(x = p.get_x()+(p.get_width()/2), # x-coordinate position of data label, padded to be in the middle of the bar
  y = height + delta, # y-coordinate position of data label, padded 100 above bar
  s = "{:.0f}".format(height), # data label, formatted to ignore decimals
  fontsize = sz*0.65,
  ha = "center",
  rotation = 90) # sets horizonta
plt.show()  
 
# ## Heatmap de matriz de interacciones por clase química y clase de diana
# f = plt.figure(figsize=(14, 14))                                                                                                                                    
# plt.matshow(df, fignum=f.number, cmap=plt.cm.get_cmap("Blues"), norm= LogNorm(vmin = 1, vmax = 27289))
# plt.yticks(range(df.shape[0]), df.index, fontsize=19)
# plt.xticks(range(df.shape[1]), df.columns, fontsize=19, rotation=90)
# plt.ylim(df.shape[0]-0.5, -0.5) # Stupid fix to set correct label positions in y axis
# cb = plt.colorbar(fraction=0.038, pad=0.14)
# cb.ax.tick_params(labelsize=14)
# #plt.rcParams['axes.titley'] = 1.5
# plt.title("Gut Compound Class vs Target Class (SEA+TC)", fontsize=27, y = 1.6)
# cc = fullchem[fullchem.set == "gut"][["hmdb_id","ccl"]].drop_duplicates().ccl.value_counts().reindex(index = comp_order, fill_value = 0)
# df_cc = fullsea[pd.notna(fullsea.tchid) & 
#                 (fullsea.set == "gut")][["hmdb_id","ccl"]].drop_duplicates()["ccl"].value_counts().reindex(index = comp_order, fill_value = 0)
# for i in range(len(df_cc)):
#     tx = str(df_cc.iloc[i]) + "(" + '{:0.2f}'.format(df_cc.iloc[i] / cc.iloc[i]*100) + ")"
#     plt.text(23-0.3, i, tx, ha='left', va='center',fontsize = 15)
# df_tc = fullsea[pd.notna(fullsea.tchid) & 
#                 (fullsea.set == "gut")][["tchid","tcl"]].drop_duplicates()["tcl"].value_counts().reindex(index = santos2_order, fill_value = 0)
# for i in range(len(df_tc)):
#     plt.text(i, 18-0.3, df_tc.iloc[i], ha='center', va='top',fontsize = 15, rotation = "vertical")
# #plt.title("FooDB Class vs Target Class: Adjusted Residuals (ChEMBL+SEA)", fontsize=27)
# for (i, j), z in np.ndenumerate(df):
#     if z > 500:
#         col = "white"
#         fs = 12
#         if z > 999:
#             fs = 11
#             if z > 9999:
#                 fs = 10
#     else:
#         col = "black"
#         fs = 15
#     plt.text(j, i, '{:0.0f}'.format(z), ha='center', va='center', fontsize = fs, color = col)
# plt.show()
 
