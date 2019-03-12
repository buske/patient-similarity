
# Patient phenotype and genotype similarity
This package is research code designed to measure the similarity between patients, using phenotype (specified as [Human Phenotype Ontology (HPO)](human-phenotype-ontology.github.io) terms) and/or genotype (specified as the [Exomiser](http://www.sanger.ac.uk/science/tools/exomiser)-processed results of whole-exome VCF files).

## Dependencies

- Python 3.3+

## Phenotypic similarity

### Input file formats

#### JSON (default, `--patient-file-format phenotips`)

By default, phenotype data is expected in [PhenoTips](https://phenotips.org) JSON export format, e.g.: `phenotips_2017-02-01_00-01.json`.

Here is an [example JSON file](test/test.json).

#### CSV (`--patient-file-format csv`)

A simple csv file format is supported, with a patient per line and columns:
1. The patient's identifier (required)
2. The patient's first present HPO term (required)
3. The patient's second present HPO term (optional)
4. The patient's third present HPO term (optional), ...

For example :
```
Patient1,HP:0000001,HP:0000002,HP:0000003,HP:0000004
Patient2,HP:0000001,HP:0000002
```

Here is an [example CSV file](test/test.csv).


### Pair-wise phenotypic similarity

Pair-wise phenotypic similarity can be computed using a number of different similarity metrics using the `patient_similarity.py` script. For example, to compute just the simGIC score:
```bash
python -m patient_similarity --log=INFO -s simgic test/test.json \
  data/hp.obo data/phenotype_annotation.tab
```

This will print to stdout the pairwise similarity scores, e.g.:
```
A	B	simgic
P0000001	P0000002	0.146613
P0000001	P0000003	0.191716
P0000001	P0000004	0.170512
P0000002	P0000003	0.124032
P0000002	P0000004	0.167785
P0000003	P0000004	0.291074
```

Multiple scores can be added by specifying `-s` multiple times, or all scores will be computed if `-s` is not specified. Supported phenotypic similarity scores include:
- jaccard
- resnik
- lin
- jc
- owlsim
- ui
- simgic
- icca
- TODO: add [ebosimgic](http://rucs.ca/computational-biology/exponential-back-off-simgic)

_See the [PhenomeCentral paper](http://dx.doi.org/10.1002/humu.22851) for a comparison of many of these_

Many of these similarity scores use the information content of the terms in the HPO to compute a similarity score. The information content of a term is defined to be `IC(t) = -log_2(p(t))`, where `p(t)` is the probability of the term. The probability of the term can be estimated in many ways, such as the fraction of OMIM diseases that have the term associated ([10.1016/j.ajhg.2008.09.017](https://dx.doi.org/10.1016%2Fj.ajhg.2008.09.017)).

A number of options have been added to support different variants of the IC computation:
- `--use-disease-prevalence`: instead of weighting each disease uniformly, weight them by their estimated prevalence from Orphanet
- `--use-phenotype-frequency`: instead of weighting each phenotype-disease association uniformly, weight them by the frequency of the association where available
- `--use-patient-phenotyes`: count each patient as an additional entry in the corpus, alongside diseases, in the frequency estimation
- `--distribute-ic-to-leaves`: evenly divide the observed frequency of each term amongst its children, so that all non-leaf nodes have zero frequency
- `--use-aoo`: include an age-of-onset similarity penalty in the similarity scoring

## Updating the data files

This package includes data files from HPO and Orphanet sources, which should be updated occasionally.

- data/hp.obo - See http://human-phenotype-ontology.github.io/downloads.html
- data/phenotype_annotations.tab - See http://human-phenotype-ontology.github.io/downloads.html
- data/en_product1.xml - See http://www.orphadata.org/cgi-bin/inc/product1.inc.php
- data/en_product2.xml - See http://www.orphadata.org/cgi-bin/inc/product2.inc.php
