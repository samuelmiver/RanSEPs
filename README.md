<div style="text-align:center"><img src ="./figures/logo_ranseps.png" /></div>

RanSEPs provides a framework for bacterial genome re-annotation and novel small proteins detection adjusting the search to different genomic features that govern protein-coding capabilities.

# How does RanSEPs work?

<center><img src="./figures/RanSEPs_functioning.png"></center>

Original publication with full description of methods can be found [here](XXXXX).

# Preparation

RanSEPs requires:
  - Python: version 2.7 or higher. We have not tested it in version 3.
  - Propy: tool to compute protein features. Instructions for downloading [here](https://www.researchgate.net/publication/235922761_UserGuide_for_propy).
  - Blast: database generation requires to run this program locally. Find information for downloadind and installation [here](https://www.ncbi.nlm.nih.gov/guide/howto/run-blast-local/).

# Installation

Specific libraries are required by RanSEPs to compute certain processes in their predictions. We provide a [requirements](./requirements.txt) file to install everything at once. To do so, you will need first to have [pip](https://pip.pypa.io/en/stable/installing/) installed and then run:

```bash
sudo apt-get install python-pip    # if you need to install pip
pip install -r requirements.txt
```

Then, move to RanSEPs directory and install the program typing:

```bash
sudo python setup.py install
```

After this you will have access to the program just typing `ranseps` in your command line.

# Usage

In order to run a prediction you will only a pair of files:
  - **Genome of reference** in fasta file.
  - The annotated **CDS** in **nucleotidic sequences** annotated in that genome in fasta.

Then just:

```bash
ranseps -g <path/to/your/genome> -c <path/to/your/cds/file>
```

This will run a simple search for proteins with size higher than 10 amino acids. However, RanSEPs allows multiple sets of parameters to explore and find the best set for your organism of interest. To check them execute:

```bash
ranseps -h
```

# RanSEPs as a python package


# Contact

[Miravet-Verde, Samuel](samuel.miravet@crg.eu)    
[Lluch-Senar, Maria](maria.lluch@crg.eu)    
[Serrano, Luis](luis.serrano@crg.eu)    

# License

RanSEPs is under a common GNU GENERAL PUBLIC LICENSE. Plese, check [LICENSE](./LICENSE) for further information.

