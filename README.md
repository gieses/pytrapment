[![GitHub](https://flat.badgen.net/github/license/gieses/pytrapment)](https://www.apache.org/licenses/LICENSE-2.0)
[![Twitter](https://flat.badgen.net/twitter/follow/SvenHGiese?icon=twitter)](https://twitter.com/SvenHGiese)
[![Python 3.8.3](https://img.shields.io/badge/python-3.8.3-blue.svg)](https://www.python.org/downloads/release/python-370/)

A python package for creating entrapment databases.

---
- [Overview](#overview)
- [Installation](#Installation)
---

## overview

pytrapment allows the convenient creation of entrapment databases for proteomic mass spectrometry analysis.
Entrapment databases are build by sampling for each protein in the host fasta file an fitting
entrapment protein.

pytrapment performs the following steps to minimize the differences between the host and entrapment
database.

1. remove all entrapment proteins that share a peptide with the host (replace I with L amino acids)
2. compute the amino acid composition for each protein in the host and entrapment database
3. for each host protein find the nearest neighbor (Euclidean distance) in the composition space
4. perform some quality control metrics on the peptides
5. save the entrapment fasta file (host + entrapment proteins)   


## Installation
pytrapment is available on pypi and can be installed via ```pip install pytrapment```

## Usage

To use pytrapment, simply call the main program via a command line:

```
pytrapment -i host.fasta -t trap.fasta -o entrapment_db
```

Make sure to have the correct paths to the fasta files. The out dir will contain the entrapment
fasta and two qc plots for peptide and protein features. The repository contains example files
which can be used as follows:

```
pytrapment -i sample_data/host.fasta -t sample_data/trap.fasta -o sample_data/
```

The results (qc_plot.png and entrapment_data.fasta) can also be found in the sample_data folder.
## Contributors
- Sven Giese