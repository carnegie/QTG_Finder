# QTG_Finder (version 1.1)

QTG-Finder is a novel machine-learning pipeline to prioritize causal genes for QTLs identified by linkage mapping. We trained QTG-Finder models for Arabidopsis and rice based on the known causal genes from each species, respectively. By utilizing additional information like poly-morphisms, function annotation, co-function network, and paralog copy number, the models can rank QTL genes to prioritize causal genes.


Authors: Fan Lin, March 2018
         Jue Fan, March 2018

### For prediction
The source code and input files can be found in the <mark>'prediction'</mark> folder. Running the 'QTG_Finder_predict.py' will require a QTL gene list provided by the user.

1. Users can prepare the QTL gene list as a single column table (.csv). See "SSQ_batch_QTL_genes.csv” for a example.

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

2. Make sure you have unzipped pre-calculated models "AT_model.dat" or "OS_model.dat" to your working directory. "AT_model.dat" is the Arabidopsis model. "OS_model.dat”is the rice model.

3. Usage ="QTG_Finder_predict.py -gl QTL_gene_list -sp species_abbreviation" <br />
QTL_gene_list: this is the list of QTL genes to be ranked. See 'SSQ_batch_QTL_genes.csv' for a example <br />
species_abbreviation: "AT" for Arabidopsis; "OS" for rice <br />
As a example,
```python3
python QTG_Finder_predict.py -gl SSQ_batch_QTL_genes.csv -sp 'AT'
```

For help,
```python3
python QTG_Finder_predict.py -h
```

### For analyses and replications

The source code and input files for cross-validation, feature importance analysis, literature validation and category analysis can be found in the <mark>‘tests'</mark> folder. The usage of each scripts (.py) is described in the beginnings of the file.

