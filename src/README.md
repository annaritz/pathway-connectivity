## Directory Contents

The `run.py` script runs all subsequent experiments, and requires that the datasets and Reactome representations have already been generated. Datasets are described in the top-level `README.md` file and in the `data/` directory.

## Reactome pathway representations
1. Directed Graph -- SIF-formatted graph is in `data/`, see the `SIF/` directory for extracting conversion types and the file used to determine which types are directed, undirected, and ignored.

2. Compound Graph -- This is the BioPAXSTREAM code from [Using Biological Pathway Data with Paxtools
](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003194) (Demir et al., PLOS Computational Biology, 2013).  See the `BioPAXSTREAM/` directory for more information.


## Experiments

### Compute Statistics

```
python3 run.py --stats
python3 run.py --stats --small_molecule_filter
python3 run.py --stats --blacklist_filter
```

### Generate Connectivity Survey Histograms

### Run Permutation Tests

### 