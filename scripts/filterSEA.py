__author__ = "Cristian Orgaz Portero"
__copyright__ = "Copyright (C) 2023 Cristian Orgaz Portero"
__license__ = "Public Domain"
__version__ = "1.0"
import pandas as pd
import numpy as np
SEA_predicted=pd.read_csv('SEA_predicted_activities.csv', sep=",", names=["Query ID","Target ID","Affinity Threshold (nM)","P-Value","Max Tc","Cut Sum","Z-Score","Name","Description","Query Smiles","InChi Key"])
SEA_predicted_filtered=SEA_predicted[(SEA_predicted["Affinity Threshold (nM)"]>=6) & (SEA_predicted["Max Tc"]>=0.4) & (-np.log10(SEA_predicted["P-Value"])>=40)].drop_duplicates()
SEA_predicted_filtered.to_csv('SEA_predicted_activities_filtered.csv', sep=";" , index=False)
SEA_predicted_filtered_duplicates=SEA_predicted_filtered.drop_duplicates(subset=["Query ID", "Target ID"])
SEA_predicted_filtered_duplicates.to_csv('SEA_predicted_activities_filtered_duplicates.csv', sep=";" , index=False)