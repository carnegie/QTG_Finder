# QTG_Finder (version 1.0)

QTG-Finder is a novel machine-learning pipeline to prioritize causal genes for QTLs identified by linkage mapping. We trained QTG-Finder models for Arabidopsis and rice based on the known causal genes from each species, respectively. By utilizing additional information like poly-morphisms, function annotation, co-function network, and paralog copy number, the models can rank QTL genes to prioritize causal genes.


Authors: Fan Lin, December 2018
         Jue Fan, December 2018

### For prediction

The source code and input files can be found in the <mark>'prediction'</mark> folder. Runing the 'QTG_Finder.py' will require a QTL gene list provided by the user.

1. Users can prepare the QTL gene list as a single column table (.csv). See “SSQ_batch_QTL_genes.csv” for a example.

// 

QTL1 name

Gene1 in QTL1

Gene2 in QTL1

Gene3 in QTL1

… 

// 

QTL2 name 

Gene1 in QTL2

Gene2 in QTL2

Gene3 in QTL2 

…

2. Make sure you have the feature list "Arabidopsis_features-v3.05_n.txt" or "rice_features_v1.3.11_n.txt" in the same directory.

3. Usage = "QTG_Finder.py input_feature_list QTL_gene_list species_abbreviation" <br />
QTL_gene_list: this is the list of QTL genes to be ranked. See 'SSQ_batch_QTL_genes.csv' for a example
feature list: use Arabidopsis_features-v3.05.txt for Arabidopsis; use rice_features_v1.3.11.txt for rice
species_abbreviation: "AT" for Arabidopsis; "OS" for rice <br />
As a example, ```python3 QTG_Finder.py Arabidopsis_features-v3.05.txt SSQ_batch_QTL_genes.csv 'AT'```


As an alternative, users can also run the batch file 'QTG_batch_Prediction.sh' on a server.



### For analyses and replications

The source code and input files for cross-validation, feature importance analysis, literature validation and category analysis can be found in the <mark>‘tests'</mark> folder. The usage of each scripts (.py) is described in the beginnings of the file.

