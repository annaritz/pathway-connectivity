## Directory Contents

The `run.py` script runs all subsequent experiments, and requires that the datasets and Reactome representations have already been generated. Datasets are described in the top-level `README.md` file and in the `data/` directory.

## Reactome pathway representations
1. Directed Graph -- SIF-formatted graph is in `data/`, see the `SIF/` directory for extracting conversion types and the file used to determine which types are directed, undirected, and ignored.

2. Compound Graph -- This is the BioPAXSTREAM code from [Algorithms for effective querying of compound graph-based pathway databases](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-376) (Dogrusoz et al., 2009), which is also described in [Using Biological Pathway Data with Paxtools](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003194) (Demir et al., PLOS Computational Biology, 2013).  See the `BioPAXSTREAM/` directory for more information about running this method (which is part of `paxtools`).

3. Hypergraph and Bipartite Graph - **this has to be updated**

For historical reasons, the hyperedges and hypernodes files are listed under `pathway-connectivity/hypergraph/reactome_full_hypergraph/`.  Many intermediate files are outputted when using the BioPAX parser.  

## Experiments

There are two useful arguments: 
- `--force` will overwrite *any* file, but the default behavior is if the file exists then it is not re-generated.  In that case, a warning will be printed to say remove the file or use the `--force` option to overwrite files.
- `--printonly` will print the commands to the console, but *not* run them. This is useful if there's an OS call that causes and error and you want to re-run just that command.

If you have any issues running these functions, send Anna an email.

### Compute Statistics

```
python3 run.py --stats
python3 run.py --stats --small_molecule_filter
python3 run.py --stats --blacklist_filter
```

### Generate Connectivity Survey Histograms

```
python3 run.py --histograms
python3 run.py --histograms --small_molecule_filter
python3 run.py --histograms --blacklist_filter
```

### Run Permutation Tests 
This and all following experiments were run on the small molecule filter only; however you can easily run these on the full hypergraph or the blacklist-filtered representations.

```
python3 run.py --perm_test 1000 --small_molecule_filter
```

Uncomment the `viz()` and `write()` functions to generate the images and/or write the Excel files. (TODO this should be an option.)

### Run Case Studies
These are hard-coded in the script.

```
python3 run.py --case_studies --small_molecule_filter
```

### Run STRING channel assessment

```
python3 run.py --string_channels --small_molecule_filter
```