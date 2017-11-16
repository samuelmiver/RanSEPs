
<center><img src="./figures/ranseps.png"></center>

RanSEPs provides a framework for genome re-annotation and novel small proteins detection adjusting the search to different genomic features that govern protein-coding capabilities.

Original publication of RanSEPs can be found [here](XXXXX).

# Preparation

RanSEPs requires:
  - Python: version 2.7 or higher. We have not tested it in version 3.
  - Propy: tool to compute protein features. Instructions for downloading [here](https://www.researchgate.net/publication/235922761_UserGuide_for_propy).
  - ViennaRNA package: framework for RNA structures predictions. Download it [here](https://www.tbi.univie.ac.at/RNA/).
  - Blast: database generation requires to run this program locally. Find information for downloadind and installation [here](https://www.ncbi.nlm.nih.gov/guide/howto/run-blast-local/).

# Installation

Specific libraries are required by RanSEPs to compute certain processes in their predictions. We provide a [requirements](./requirements.txt) file to install everything at once. To do so, you will need first to have [pip](https://pip.pypa.io/en/stable/installing/) installed and then run:

```bash
sudo apt-get install python-pip    # if you need to install pip
pip install -r requirements.txt
```




# License

RanSEPs is under a common GNU GENERAL PUBLIC LICENSE. Plese, check [LICENSE](./LICENSE) for further information. 
