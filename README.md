# PRODIGY-LIG / Binding Affinity Prediction for Protein-Small Ligands

This is a collection of scripts to predict the binding affinity
values for protein-small ligands complexes, starting from their
atominc 3D-structure (in the PDB of mmCIF format)

The online version of PRODIGY-LIG predictor can be found here:
* [PRODIGY-LIG](http://milou.science.uu.nl/services/PRODIGY/)

Details of the binding affinity predictor implemented in PRODIGY can be found here:
* [The server] (Submitted: Vangone et al. "Large-scale prediction of binding affinity in protein-small ligand complexes: the PRODIGY-LIG web server")
* [Results of PRODIGY-LIG in D3F](https://www.ncbi.nlm.nih.gov/pubmed/?term=Performance+of+HADDOCK+and+a+simple+contact-based+protein%E2%80%93ligand+binding+affinity+predictor+in+the+D3R+Grand+Challenge+2)

# Quick & Dirty Installation
```bash

git clone http://github.com/haddocking/prodigy-lig

# ADD STUFF

# Have fun!
```

# Usage

```bash
python prodigy-lig.py <pdb file> [OTHER]
```

Type --help to get a list of all the possible options of the script.
# Installation & Dependencies
The scripts rely on [Biopython](www.biopython.org) to validate the PDB structures and calculate
interatomic distances. [freesasa](https://github.com/mittinatten/freesasa), with the parameter
set used in NACCESS ([Chothia, 1976](http://www.ncbi.nlm.nih.gov/pubmed/994183)), is also
required for calculating the buried surface area.


To install and use the scripts, just clone the git repository or download the tarball zip
archive

# License
These utilities are open-source and licensed under the Apache License 2.0. For more information
read the LICENSE file.

# Citing us
If our predictive model or any scripts are useful to you, consider citing them in your
publications:

**Vangone A, Schaarschmidt J, Koukos P, Geng C, Citro N, Trellet M, Xue L, Bonvin A.**: Large-scale prediction of binding affinity in protein-small ligand complexes: the PRODIGY-LIG web server. *Bioinformatics submitted*

**Kurkcuoglu Z, Koukos P, Citro N, Trellet M, Rodrigues J, Moreira I, Roel-Touris J, Melquiond A, Geng C, Schaarschmidt J,
Xue L, Vangone A, Bonvin AMJJ.**: Performance of HADDOCK and a simple contact-based protein-ligand binding affinity
predictor in the D3R Grand Challenge 2. *J Comput Aided Mol Des* 32(1):175-185 (2017). ([link](https://www.ncbi.nlm.nih.gov/pubmed/?term=Performance+of+HADDOCK+and+a+simple+contact-based+protein%E2%80%93ligand+binding+affinity+predictor+in+the+D3R+Grand+Challenge+2))

# Contact
For questions about PRODIGY usage, please contact the team at: prodigy.bonvinlab@gmail.com
