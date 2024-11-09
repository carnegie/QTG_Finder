# QTG_Finder (version 2.0)

QTG-Finder is a machine-learning pipeline to prioritize causal genes for QTLs identified by linkage mapping. We trained QTG-Finder models for Arabidopsis, rice, sorghum, Setaria viridis based on known causal genes and orthologs of known causal genes, respectively. By utilizing additional information like polymorphisms, function annotation, co-function network, paralog copy number, the models can prioritize causal genes for QTLs identified by QTL mapping.

Authors: Fan Lin, February 2020<br />
Environment: 
Python=3.7.12
pandas=0.25.3
numpy=1.19.5
scikit-learn=0.21.2

### For prediction
The source code and input files can be found in the <mark>'QTG2_prediction'</mark> folder. Running the 'QTG_Finder_predict.py' will require a QTL gene list provided by the user.

1. Users can prepare the QTL gene list as a single column table (.csv). See "SV_height_QTL_example.csv" or "AT_Seedsize_QTL_example.csv" for a example.

||
|:-| 
|//|
|QTL1 name|
|Gene1 in QTL1|
|Gene2 in QTL1|
|Gene3 in QTL1|
|…| 
|//|
|QTL2 name|
|Gene1 in QTL2|
|Gene2 in QTL2|
|Gene3 in QTL2|
|…|

2. The pre-calculated models can be downloaded from the following links: <br />
Arabidopsis: https://carnegiedpb.s3.amazonaws.com/software/QTG2_prediction/AT_model.dat.zip<br />
Rice: https://carnegiedpb.s3.amazonaws.com/software/QTG2_prediction/OS_model.dat.zip<br />
Sorghum: https://carnegiedpb.s3.amazonaws.com/software/QTG2_prediction/SB_model.dat.zip<br />
Setaria: https://carnegiedpb.s3.amazonaws.com/software/QTG2_prediction/SV_model.dat.zip<br />

3. Unzip the pre-calculated models in working directory: ./QTG2_prediction<br />
Example:
```java
jar xvf AT_model.dat.zip
```

4. Usage: “QTG_Finder_predict.py -gl QTL_gene_list -sp species_abbreviation" <br />
QTL_gene_list: this is the list of QTL genes to be ranked. See "SV_height_QTL_example.csv" for a example <br />
species_abbreviation: "AT" for Arabidopsis; "OS" for rice; "SB" for sorghum;"SV" for Setaria viridis  <br />
As a example,
```python3
python QTG_Finder_predict.py -gl SV_height_QTL_example.csv -sp 'SV'
```

For help,
```python3
python QTG_Finder_predict.py -h
```
5. “QTL_gene_rank.csv” will be the output file. 

### For analyses and replications

The source code and input files for cross-validation, feature importance analysis, literature validation and category analysis can be found in the <mark>'QTG2_analysis'</mark> folder. The usage of each scripts (.py) is described at the beginning of them.

