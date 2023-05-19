En este repositorio se encuentran los archivos de código y los resultados obtenidos tras su ejecución para la elaboración del trabajo de fin de master Análisis y modelado de actividades biológicas de metabolitos intestinales como fuente de diseño de fármacos.

La versión de python necesaria es 3.11

El orden de ejecución y función se resume a continuación:
* getCIDs.py (Obtención de component id de los metabolitos en PubChem)
* getBioactivities.py (Obtención de bioactividades de los metabolitos en PubChem)
* processBioactivities.py (Procesado de bioactividades obtenidas a partir de getBioactivities.py)
* getSmiles.py (Obtención de Canonical Smiles en PubChem para crear el archivo que usaremos en tldr input.txt)
* filterSEA.py (Para ejecutar este archivo es necesario es renombrar el archivo output.csv obtenido a partir de tldr a SEA_predicted_activities.csv )
* getUniprots.py (Pubchem no siempre devuelve el uniprot id de la diana encontrada por lo que es necesario buscar el uniprot)
* getPubChemTargetsFromChembl.py (Obtención de la clase de diana de los resultados obtenidos en pubchem)
* getSEATargetsFromChembl.py (Obtención de la clase de diana de los resultados obtenidos en tldr)