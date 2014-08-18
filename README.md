
# Patient phenotype and genotype similarity #
This package is research code designed to measure the similarity between patients, using both phenotype (specified as Human Phenotype Ontology terms) and genotype (specified as the Exomiser-processed results of whole-exome VCF files). 


## Usage example ##
```bash
# Run on example data (expected output in test/test.out)
./patient_similarity.py --log=INFO -s simgic test/test.hpo \
  data/hp.obo data/phenotype_annotation.tab \
  data/en_product1.xml data/en_product2.xml
```