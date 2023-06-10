__author__ = "Cristian Orgaz Portero"
__copyright__ = "Copyright (C) 2023 Cristian Orgaz Portero"
__license__ = "Public Domain"
__version__ = "1.0"

import mysql.connector as conn
from mysql.connector import Error
import pandas as pd
import numpy as np

SEA_predicted_guts_filtered=pd.read_csv('../results/SEA/SEA_predicted_activities_guts_filtered.csv', sep=";")
target_ids_SEA_guts="".join(SEA_predicted_guts_filtered["Target ID"].drop_duplicates().apply(lambda x: "'" + x + "',").to_list())[:-1]
query_SEA_guts="select distinct chembl_31.component_sequences.component_id, chembl_31.component_sequences.accession, chembl_31.component_sequences.organism, chembl_31.component_sequences.tax_id, chembl_31.protein_family_classification.*, chembl_31.organism_class.l1 as 'domain', chembl_31.component_synonyms.component_synonym  from chembl_31.component_sequences"\
+" join chembl_31.component_class on chembl_31.component_class.component_id = chembl_31.component_sequences.component_id"\
+" join chembl_31.protein_family_classification on chembl_31.protein_family_classification.protein_class_id=chembl_31.component_class.protein_class_id"\
+" join chembl_31.organism_class on chembl_31.organism_class.tax_id = chembl_31.component_sequences.tax_id"\
+" join chembl_31.component_synonyms on chembl_31.component_synonyms.component_id = chembl_31.component_sequences.component_id"\
+" where chembl_31.component_synonyms.component_synonym in ("+ target_ids_SEA_guts+")"\
+" and (binary chembl_31.organism_class.l1 like 'Bacteria' or binary chembl_31.component_sequences.organism like 'Homo sapiens');"

SEA_predicted_drugs_filtered=pd.read_csv('../results/SEA/SEA_predicted_activities_drugs_filtered.csv', sep=";")
target_ids_SEA_drugs="".join(SEA_predicted_drugs_filtered["Target ID"].drop_duplicates().apply(lambda x: f'\'{x}\',').to_list())[:-1]
query_SEA_drugs="select distinct chembl_31.component_sequences.component_id, chembl_31.component_sequences.accession, chembl_31.component_sequences.organism, chembl_31.component_sequences.tax_id, chembl_31.protein_family_classification.*, chembl_31.organism_class.l1 as 'domain', chembl_31.component_synonyms.component_synonym  from chembl_31.component_sequences"\
+" join chembl_31.component_class on chembl_31.component_class.component_id = chembl_31.component_sequences.component_id"\
+" join chembl_31.protein_family_classification on chembl_31.protein_family_classification.protein_class_id=chembl_31.component_class.protein_class_id"\
+" join chembl_31.organism_class on chembl_31.organism_class.tax_id = chembl_31.component_sequences.tax_id"\
+" join chembl_31.component_synonyms on chembl_31.component_synonyms.component_id = chembl_31.component_sequences.component_id"\
+" where chembl_31.component_synonyms.component_synonym in ("+ target_ids_SEA_drugs+")"\
+" and (binary chembl_31.organism_class.l1 like 'Bacteria' or binary chembl_31.component_sequences.organism like 'Homo sapiens');"

try:
    connection = conn.connect(host='localhost',
                                         database='chembl_31',
                                         user='root',
                                         password='admin')
    if connection.is_connected():
        db_Info = connection.get_server_info()
        print("Connected to MySQL Server version ", db_Info)
        unitProtId_guts = pd.read_sql(query_SEA_guts, con=connection)
        unitProtId_guts.to_csv('../results/SEA/target_uniprot_ids_sea_guts.csv', sep=";", index=False)
        unitProtId_drugs = pd.read_sql(query_SEA_drugs, con=connection)
        unitProtId_drugs.to_csv('../results/SEA/target_uniprot_ids_sea_drugs.csv', sep=";", index=False)

except Error as e:
    print("Error while connecting to MySQL", e)
finally:
    if connection.is_connected():
        connection.close()
        print("MySQL connection is closed")
sea_target_class_guts=pd.read_csv('../results/SEA/target_uniprot_ids_sea_guts.csv', sep=";")
sea_target_class_guts["tcl"] = sea_target_class_guts.apply(lambda x: x.l2 if x.l2 in ["Transferase","Oxidoreductase","Protease",
                                                       "Hydrolase","Nuclear receptor","Electrochemical transporter",
                                                       "Phosphodiesterase","Lyase","Phosphatase",
                                                       "Isomerase","Ligase"] else 
         "Protease" if x.l2 in ["Protease","Other ion channel|Protease"] else
         "Cytochrome P450" if x.l2 in ["Cytochrome P450"] else   
         "Oxidoreductase" if x.l2 in ["Oxidoreductase","Isomerase|Oxidoreductase"] else            
         "7TM1" if x.l2 in ["Family A G protein-coupled receptor"] else 
         "7TM2" if x.l2 in ["Family B G protein-coupled receptor"] else
         "Kinase" if x.l2 in ["Kinase","Kinase|Transferase"] else
         "7TM3" if x.l2 == "Family C G protein-coupled receptor" else
         "LGIC" if x.l2 == "Ligand-gated ion channel" else
         "VGIC" if x.l2 in ["Voltage-gated ion channel",
                            "Voltage-gated ion channel|Slow voltage-gated potassium channel accessory protein family"] else 
         "Primary active transporter" if x.l2 in ["Primary active transporter","Primary active transporter|Hydrolase",
                                                  "Voltage-gated ion channel|Primary active transporter",
                                                   "Hydrolase|Primary active transporter",
                                                   "Hydrolase|Other ion channel|Primary active transporter"] else
         "Epigenetic regulator" if x.l2 in(["Eraser","Writer","Reader",
                                            "Reader|Writer",
                                            "Eraser|Reader",
                                            "Eraser|Writer",
                                            "Protease|Reader",
                                            "Writer|Reader"]) else  
         "Enzyme other" if x.l2 in ([None,"Aminoacyltransferase",np.nan]) and x.l1 in (["Enzyme","Enzyme|Structural protein",
                                                                                        "Unclassified protein|Enzyme"]) else
         "Other cytosolic protein" if x.l2 is None and x.l1 in(["Other cytosolic protein",
                                                                "Other cytosolic protein|Unclassified protein",
                                                              "Other nuclear protein|Other cytosolic protein"]) else 
         "Secreted protein" if x.l2 is None and x.l1 == "Secreted protein" else
         "Transcription factor other" if x.l1 in (["Transcription factor","Other cytosolic protein|Transcription factor"]) and x.l2 is None else 
         "Membrane receptor other" if x.l1 in ["Membrane receptor","Surface antigen|Membrane receptor",
                                                "Adhesion|Surface antigen|Membrane receptor",
                                                "Adhesion|Membrane receptor",
                                                "Enzyme|Membrane receptor"] and x.l2 in [None,
                                                                                           "Taste family G protein-coupled receptor",
                                                                                           "Toll-like and Il-1 receptors",
                                                                                           "Hydrolase|Toll-like and Il-1 receptors"] 
                                                                                           else
         "Other", axis = 1)
                                                                                         
sea_target_class_drugs=pd.read_csv('../results/SEA/target_uniprot_ids_sea_drugs.csv', sep=";")
sea_target_class_drugs["tcl"] = sea_target_class_drugs.apply(lambda x: x.l2 if x.l2 in ["Transferase","Oxidoreductase","Protease",
                                                       "Hydrolase","Nuclear receptor","Electrochemical transporter",
                                                       "Phosphodiesterase","Lyase","Phosphatase",
                                                       "Isomerase","Ligase"] else 
         "Protease" if x.l2 in ["Protease","Other ion channel|Protease"] else
         "Cytochrome P450" if x.l2 in ["Cytochrome P450"] else   
         "Oxidoreductase" if x.l2 in ["Oxidoreductase","Isomerase|Oxidoreductase"] else            
         "7TM1" if x.l2 in ["Family A G protein-coupled receptor"] else 
         "7TM2" if x.l2 in ["Family B G protein-coupled receptor"] else
         "Kinase" if x.l2 in ["Kinase","Kinase|Transferase"] else
         "7TM3" if x.l2 == "Family C G protein-coupled receptor" else
         "LGIC" if x.l2 == "Ligand-gated ion channel" else
         "VGIC" if x.l2 in ["Voltage-gated ion channel",
                            "Voltage-gated ion channel|Slow voltage-gated potassium channel accessory protein family"] else 
         "Primary active transporter" if x.l2 in ["Primary active transporter","Primary active transporter|Hydrolase",
                                                  "Voltage-gated ion channel|Primary active transporter",
                                                   "Hydrolase|Primary active transporter",
                                                   "Hydrolase|Other ion channel|Primary active transporter"] else
         "Epigenetic regulator" if x.l2 in(["Eraser","Writer","Reader",
                                            "Reader|Writer",
                                            "Eraser|Reader",
                                            "Eraser|Writer",
                                            "Protease|Reader",
                                            "Writer|Reader"]) else  
         "Enzyme other" if x.l2 in ([None,"Aminoacyltransferase",np.nan]) and x.l1 in (["Enzyme","Enzyme|Structural protein",
                                                                                        "Unclassified protein|Enzyme"]) else
         "Other cytosolic protein" if x.l2 is None and x.l1 in(["Other cytosolic protein",
                                                                "Other cytosolic protein|Unclassified protein",
                                                              "Other nuclear protein|Other cytosolic protein"]) else 
         "Secreted protein" if x.l2 is None and x.l1 == "Secreted protein" else
         "Transcription factor other" if x.l1 in (["Transcription factor","Other cytosolic protein|Transcription factor"]) and x.l2 is None else 
         "Membrane receptor other" if x.l1 in ["Membrane receptor","Surface antigen|Membrane receptor",
                                                "Adhesion|Surface antigen|Membrane receptor",
                                                "Adhesion|Membrane receptor",
                                                "Enzyme|Membrane receptor"] and x.l2 in [None,
                                                                                           "Taste family G protein-coupled receptor",
                                                                                           "Toll-like and Il-1 receptors",
                                                                                           "Hydrolase|Toll-like and Il-1 receptors"] 
                                                                                           else
         "Other", axis = 1)
 
sea_target_class_drugs.to_csv('../results/SEA/sea_target_class_drugs.csv', sep=";", index=False)
sea_target_class_guts.to_csv('../results/SEA/sea_target_class_guts.csv',sep=";",index=False)