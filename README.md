## Hypergraph-based connectivity measures for signaling pathway topologies

Data and source code for (nearly) all figures and experiments the following paper:

Nick Franzese, Adam Groce, T.M. Murali, and Anna Ritz. **Hypergraph-based connectivity measures for signaling pathway topologies**.  _PLOS Computational Biology_ 2019; 15(10): e1007384. [(publisher link)](https://doi.org/10.1371/journal.pcbi.1007384)

Before publication in _PLOS Computational Biology_, this work was presented at GLBio 2019 (as an accepted paper, see versions on [bioRxiv](https://doi.org/10.1101/593913) and at ISMB 2019 (as an [accepted poster](https://www.reed.edu/biology/ritz/files/posters/2019-franzese-groce-murali-ritz.pdf)).

Email Anna (aritz@reed.edu) with any questions about this code base.

### Dependencies

Most of the experiment options require these tools:
- Python (tested with v3.5.1 and v3.8.2)
- [Hypergraph Algorithms Library in Python (HALP)](http://murali-group.github.io/halp/)
  - Check out `annabranch` (this passes all `pytest` tests; there are currently some issues with Travis CI)
- [NetworkX](https://networkx.github.io/) (v1.11, should work with v2.1)
- [xlsxwriter](https://xlsxwriter.readthedocs.io/) to write permutation test information to an Excel file.
- [matplotlib_venn](https://pypi.org/project/matplotlib-venn/) to visualize venn diagrams.

If you want to implement shortest B-hyperpath on directed graphs, it is implemented as a `cplex` integer linear program.
CPLEX will be imported if `--histograms` is specified, since it is used in the connectivity survey.
- [IBM ILOG CPLEX optimizer Python module](https://www.ibm.com/analytics/cplex-optimizer) (academic license available)

### Datasets

The data is generated from the Reactome Pathway, available via [PathwayCommons](http://www.pathwaycommons.org/).  PC v10 files available [here](http://www.pathwaycommons.org/archives/PC2/v10/).  See details in the `data/` directory for more information.

To run all analyses, you must generate the following datasets:
- Simple Interaction Format (SIF, instructions in `data/`)
- Hypergraph, Bipartite Graph, and Compound Graph (OWL, instructions in `data/`)

The four data representations (hypergraph, bipartite graph, directed graph, and compound graph) are available in this repo.  To re-generate all representations, you can run the scripts outlined in the `src/` directory.  

### Usage Information

The main function is `run.py` in the `src/` directory:

```
usage: run.py [-h] [--force] [--printonly] [--keep_singletons] [--small_molecule_filter]
              [--ubiquitous_filter] [--stats] [--histograms] [--perm_test #PERMS] [--set_seed]
              [--case_studies] [--string_channels]

Run experiments for pathway connectivity. The four representations are "SIF-Graph","Compound
Graph","Bipartite Graph", and "Hypergraph"

optional arguments:
  -h, --help            show this help message and exit

General Arguments:
  --force               force existing files to be overwritten (default=False)
  --printonly           print the commands instead of running them (default=False)

Dataset Arguments:
  --keep_singletons     Keep singleton nodes. Default False.
  --small_molecule_filter
                        Filter by small molecule nodes
  --ubiquitous_filter   Filter by PathwayCommons ubiquitous nodes

Experiment Arguments:
  --stats               print statistics about each representation.
  --histograms          print histograms and heatmaps for each representation.
  --perm_test #PERMS    pathway influence permutation test with #PERMS number of permutations.
  --set_seed            Set the seed of the random number generator to 123456. Default False.
  --case_studies        pathway influence case studies (hard-coded)
  --string_channels     run STRING channel assessment.
  ```

Refer to the `README.md` files for more commands to generate data, run experiments, etc.

### Developer Notes

- We need to put this in a virtual environment or dockerize it.  
- TravisCI testing fails for the hypergraph algorithms library (halp). However, the `annabranch` branch (which includes tests for B-Relaxation Distance) have all unit tests passing.
- There are two places where we describe PaxTools and downloading the JAR (src/BioPAXSTREAM and src/hypergraph_code/hypergraph_parser).Further, there are two versions of the JAR file needed.  Combine this (or better yet wrap it into an easy install script).
- `src/BioPAXSTREAM/` code still has some hard-coded lines in it. Currently users are instructed to replace these lines; they should be put in a command-line-argument.
- The `hypergraph/` directory at the top level is out of place (for historical reasons).  It needs to be movable and other scripts shouldn't have it hard-coded.

#### Developer Log

_Jun 2020:_ Renamed "blacklisted" files to "ubiquitous" files in response to the racist connotations with using "blacklist" and "whitelist."  Read more about [GitHub's phrase changes](https://www.zdnet.com/article/github-to-replace-master-with-alternative-term-to-avoid-slavery-references/), and look through your own repos to remove racist terminology.

_Jun 2020:_ Tested code on Python v3.8 and changed `scipy.stats` to `scipy.special` package in visualization script.

_Sep 2019:_ Code to parse the Reactome OWL files into hypergraphs are now documented in `src/hypergraph_code/hypergraph_parser/`.  

_Aug 2019:_ Full dependencies, code for generating _almost_ all pathway representations (hypergraphs, bipartite graphs, compound graphs, and directed graphs) are now here.  Third-party datasets are also documented.  Email Anna with any questions.

_May 2019:_ This repo is still a little rough, Anna plans to refactor code and update instructions for data acquisition and processing in the next few weeks. Email her if you have questions.
