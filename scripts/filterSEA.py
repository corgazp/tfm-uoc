__author__ = "Cristian Orgaz Portero"
__copyright__ = "Copyright (C) 2023 Cristian Orgaz Portero"
__license__ = "Public Domain"
__version__ = "1.0"
import pandas as pd
import numpy as np
SEA_predicted_drugs=pd.read_csv('../results/SEA/SEA_predicted_activities_drugbank.csv', sep=",", names=["Query ID","Target ID","Affinity Threshold (nM)","P-Value","Max Tc","Cut Sum","Z-Score","Name","Description","Query Smiles","InChi Key"])
SEA_predicted_drugs_filtered=SEA_predicted_drugs[(SEA_predicted_drugs["Affinity Threshold (nM)"]>=6) & (SEA_predicted_drugs["Max Tc"]>=0.4) & (-np.log10(SEA_predicted_drugs["P-Value"])>=40)].drop_duplicates()
SEA_predicted_drugs_filtered.to_csv('../results/SEA/SEA_predicted_activities_drugs_filtered.csv', sep=";" , index=False)
SEA_predicted_drugs_filtered_duplicates=SEA_predicted_drugs_filtered.drop_duplicates(subset=["Query ID", "Target ID"])
SEA_predicted_drugs_filtered_duplicates.to_csv('../results/SEA/SEA_predicted_activities_drugs_filtered_duplicates.csv', sep=";" , index=False)

SEA_predicted_guts=pd.read_csv('../results/SEA/SEA_predicted_activities_guts.csv', sep=",")
SEA_predicted_guts_filtered=SEA_predicted_guts[(SEA_predicted_guts["Affinity Threshold (nM)"]>=6) & (SEA_predicted_guts["Max Tc"]>=0.4) & (-np.log10(SEA_predicted_guts["P-Value"])>=40)].drop_duplicates()
SEA_predicted_guts_filtered.to_csv('../results/SEA/SEA_predicted_activities_guts_filtered.csv', sep=";" , index=False)
SEA_predicted_guts_filtered_duplicates=SEA_predicted_guts_filtered.drop_duplicates(subset=["Query ID", "Target ID"])
SEA_predicted_guts_filtered_duplicates.to_csv('../results/SEA/SEA_predicted_activities_guts_filtered_duplicates.csv', sep=";" , index=False)