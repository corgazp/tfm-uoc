__author__ = "Cristian Orgaz Portero"
__copyright__ = "Copyright (C) 2023 Cristian Orgaz Portero"
__license__ = "Public Domain"
__version__ = "1.0"

import pandas as pd
from rdkit.Chem import Draw
from rdkit import Chem
import seaborn as sns
import matplotlib.pyplot as plt
import math
import os
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from scipy import ndimage
path = "../results/images"
isExist = os.path.exists(path)
if not isExist:
   os.makedirs(path)
   print("Directory results/images has been created!")
   
df_sc_guts=pd.read_csv('../results/guts_scaffolds_smiles.csv',sep=";")
df_sc_drugs=pd.read_csv('../results/drugs_scaffolds_smiles.csv',sep=";")

data_guts_generic_scaffold={"smiles_generic_scaffolds":df_sc_guts.smiles_from_generic_scaffold.value_counts().index.tolist(),"count":df_sc_guts.smiles_from_generic_scaffold.value_counts().tolist()}
data_drugs_generic_scaffold={"smiles_generic_scaffolds":df_sc_drugs.smiles_from_generic_scaffold.value_counts().index.tolist(),"count":df_sc_drugs.smiles_from_generic_scaffold.value_counts().tolist()}
data_guts_scaffold={"smiles_scaffolds":df_sc_guts.smiles_from_scaffold.value_counts().index.tolist(),"count":df_sc_guts.smiles_from_scaffold.value_counts().tolist()}
data_drugs_scaffold={"smiles_scaffolds":df_sc_drugs.smiles_from_scaffold.value_counts().index.tolist(),"count":df_sc_drugs.smiles_from_scaffold.value_counts().tolist()}

dat_count_sc_generic_guts=pd.DataFrame(data=data_guts_generic_scaffold)
dat_count_sc_generic_guts_subset=dat_count_sc_generic_guts[0:10]
dat_count_sc_generic_drugs=pd.DataFrame(data=data_drugs_generic_scaffold)
dat_count_sc_generic_drugs_subset=dat_count_sc_generic_drugs[0:10]

dat_count_sc_guts=pd.DataFrame(data=data_guts_scaffold)
dat_count_sc_guts_subset=dat_count_sc_guts[0:10]
dat_count_sc_drugs=pd.DataFrame(data=data_drugs_scaffold)
dat_count_sc_drugs_subset=dat_count_sc_drugs[0:10]

arr_generic_scaffolds_guts=[]
arr_scaffolds_guts=[]
arr_generic_scaffolds_drugs=[]
arr_scaffolds_drugs=[]

def get_molecules_from_scaffolds(array_scaffolds, image_prefix, arr_img):
    for i in range(0,len(array_scaffolds)):
        smiles=array_scaffolds[i]
        m=Chem.MolFromSmiles(smiles)
        img = Draw.MolToImage(m)
        arr_img.append(img)
        img.save(f'../results/images/{image_prefix}_{smiles}_{i}.png')

get_molecules_from_scaffolds(dat_count_sc_generic_guts_subset.smiles_generic_scaffolds, "gut_generic", arr_generic_scaffolds_guts)
get_molecules_from_scaffolds(dat_count_sc_generic_drugs_subset.smiles_generic_scaffolds, "drugs_generic", arr_generic_scaffolds_drugs)
get_molecules_from_scaffolds(dat_count_sc_guts_subset.smiles_scaffolds, "gut", arr_scaffolds_guts)
get_molecules_from_scaffolds(dat_count_sc_drugs_subset.smiles_scaffolds, "drugs", arr_scaffolds_drugs)

def get_bar_plot_with_scaffolds(arr_scaffolds, var, arr_scaffolds_img, title, x_axis_title, y_axis_title,crr_1,crr_2):
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(25,20))
    
    fig.subplots_adjust(hspace=0.05)
    sns.barplot(x = var, y = "count", data = arr_scaffolds, orient = "v", color="blue", ax=ax1)
    sns.barplot(x = var, y = "count", data = arr_scaffolds, orient = "v", color="blue", ax=ax2)
    ax1.set_ylim(arr_scaffolds["count"][0] - 20,arr_scaffolds["count"][0] + 20)  # outliers only
    ax2.set_ylim(0, arr_scaffolds["count"][1]+20)
    ax1.tick_params(axis='x', rotation=45)
    ax2.tick_params(axis='x', rotation=45)
    delta=1
    sz = 16
    ax2.set_xlabel(x_axis_title, fontsize = sz*1.1)
    ax1.set_xlabel("", fontsize = sz*1.1)
    ax1.set_ylabel(y_axis_title, fontsize = sz*1.1)
    ax1.yaxis.set_label_coords(-0.05, 0)
    ax2.set_ylabel("", fontsize = sz*1.1)
    ax1.set_title(title, fontsize = sz*1.3)
    ax1.set_yticklabels(ax1.get_yticks(), size = sz*0.9)
    ax2.set_yticklabels(ax2.get_yticks(), size = sz*0.9)
    ax1.set_xticklabels("", size = sz*0.9)
    ax1.tick_params(axis='x', rotation=45)
    ax1.spines.bottom.set_visible(False)
    ax2.spines.top.set_visible(False)
    ax1.xaxis.set_visible(False)
    ax1.tick_params(labeltop=False)  # don't put tick labels at the top
    ax2.xaxis.tick_bottom()
    plt.margins(y = 0.15, x = 0.025)
    for i in range(0,len(ax1.patches)):
      # get the height of each bar
      height = ax1.patches[i].get_height()
      # adding text to each bar
      if(ax1.patches[i].get_height()>arr_scaffolds["count"][1]):
          ax1.text(x = ax1.patches[i].get_x()+(ax1.patches[i].get_width()/2), # x-coordinate position of data label, padded to be in the middle of the bar
          y = height + delta, # y-coordinate position of data label, padded 100 above bar
          s ="{:.0f}".format(0) if math.isnan(height) else "{:.0f}".format(height), # data label, formatted to ignore decimals
          fontsize = sz*0.9,
          ha = "center")
          rotated_img = ndimage.rotate(arr_scaffolds_img[i], 90)
          ib = OffsetImage(rotated_img, zoom=.4)
          ab = AnnotationBbox(ib,
           (0,0),
           xycoords='data',
           frameon=False,
           box_alignment=(0.5, -(height+delta)/crr_1))
          ax2.add_artist(ab) 
    delta = 2
    for i in range(0,len(ax2.patches)):
      # get the height of each bar
      height = ax2.patches[i].get_height()
      # adding text to each bar
      if(ax2.patches[i].get_height()<=arr_scaffolds["count"][1]):
          ax2.text(x = ax1.patches[i].get_x()+(ax1.patches[i].get_width()/2), # x-coordinate position of data label, padded to be in the middle of the bar
          y = height + delta, # y-coordinate position of data label, padded 100 above bar
          s ="{:.0f}".format(0) if math.isnan(height) else "{:.0f}".format(height), # data label, formatted to ignore decimals
          fontsize = sz*0.9,
          ha = "center")

          rotated_img = ndimage.rotate(arr_scaffolds_img[i], 90)
          ib = OffsetImage(rotated_img, zoom=.4)
          ab = AnnotationBbox(ib,
            (i,height+delta+crr_2),
            xycoords='data',
            frameon=False)
          ax2.add_artist(ab) 
    plt.show()
get_bar_plot_with_scaffolds(dat_count_sc_drugs_subset, "smiles_scaffolds" , arr_scaffolds_drugs, "Scaffolds más frecuentes DB", "Scaffolds", "Frecuencias",16.7,30/6)
get_bar_plot_with_scaffolds(dat_count_sc_generic_drugs_subset, "smiles_generic_scaffolds" , arr_generic_scaffolds_drugs, "Scaffolds genéricos más frecuentes DB", "Scaffolds genéricos", "Frecuencias",22.7,60/6)
get_bar_plot_with_scaffolds(dat_count_sc_guts_subset, "smiles_scaffolds" , arr_scaffolds_guts, "Scaffolds más frecuentes Intestino", "Scaffolds", "Frecuencias",14.8,45/7)
get_bar_plot_with_scaffolds(dat_count_sc_generic_guts_subset, "smiles_generic_scaffolds" , arr_generic_scaffolds_guts, "Scaffolds genéricos más frecuentes Intestino", "Scaffolds genéricos", "Frecuencias",24,70/7)
