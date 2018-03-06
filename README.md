# Introduction

# Quick start

Assuming that you know what you are doing you can fetch a structure from the
PDB and start using the code immediately. The following two lines should do
the trick:

```
wget https://files.rcsb.org/download/1A0Q.pdb
prodigy_lig.py -c H,L H:HEP -i 1A0Q.pdb
```

If all goes well this will produce the following two lines of output.

```
Job name        DGprediction (low refinement) (Kcal/mol)
1a0q    -8.49
```

You can read more about the usage of prodigy_lig in the next section

# Usage

Running `prodigy_lig.py -h` will list some basic information about the code.

```
usage: prodigy_lig.py [-h] -c CHAINS CHAINS -i INPUT_FILE [-e ELECTROSTATICS]
                      [-d DISTANCE_CUTOFF] [-o] [-v] [-V]

Calculate the Binding Affinity score using the PRODIGY-LIG model

prodigy_lig dependes on biopython for the structure manipulations
and only requires a single structure file (in mmCIF or PDB format)
as input.

prodigy_lig is licensed under the Apache License 2.0 included in
the LICENSE file of this repository or at the following URL

https://github.com/haddocking/prodigy-lig/blob/master/LICENSE

If you use prodigy_lig in your research please cite the following
papers:

1. to be submitted
2. https://doi.org/10.1007/s10822-017-0049-y

optional arguments:
  -h, --help            show this help message and exit
  -c CHAINS CHAINS, --chains CHAINS CHAINS
                        Which chains to use. Expects two sets of arguments.
                        The first set refers to the protein selection and
                        multiple chains can be specified by separating the
                        chain identifiers with commas. The second set refers
                        to the ligand and requires one chain and the residue
                        identifier of the ligand. A typical use case could be
                        the following: prodigy_lig.py -c A,B A:LIG
  -i INPUT_FILE, --input_file INPUT_FILE
                        This is the PDB/mmcif file for which the score will be
                        calculated.
  -e ELECTROSTATICS, --electrostatics ELECTROSTATICS
                        This is the electrostatics energy as calculated during
                        the water refinement stage of HADDOCK.
  -d DISTANCE_CUTOFF, --distance_cutoff DISTANCE_CUTOFF
                        This is the distance cutoff for the Atomic Contacts
                        (def = 10.5Å).
  -o, --output_file     Store the processed file. The filename will be the
                        name of the input with -processed appended just before
                        the file ending.
  -v, --verbose         Include the calculated contact counts in the output.
  -V, --version         Print the version and exit.

Authors: Panagiotis Koukos, Anna Vangone, Joerg Schaarschmidt
```

## Command line arguments

The command line arguments belong to one of two categories; required and optional.
The two required arguments are the input structure file and the specification of
the chain and residue identifiers of the interactors.

### Required arguments

#### Input file

The input structure file can be specified with the `-i` flag (or `--input_file`)
and it should be the path to the input PDB/mmCIF file.

#### Chain specification

The `-c` flag (or `--chains`) must be used to specify the chain identifiers for
the interactors, and for the ligand the residue identifier as well. The first
argument allows to specify the protein chains to be used in the analysis. Multiple
chains can be specified by comma separating them. The second argument allows us to
specify the chain and residue identifier of the ligand of interest. The chain and
residue identifier must be separated by a colon (`:`).

For the next examples we will be using the structure from the quickstart [section](#quick-start), 1A0Q.

 If you didn't do so before you can fetch the PDb file with the following line of
 code from the command line.

 ```
 wget https://files.rcsb.org/download/1A0Q.pdb
 ```

Or by simply pointing your browser [here](https://files.rcsb.org/download/1A0Q.pdb)

After examining the file with your favourite molecular viewer you will note that
this is a Fab fragment with a small molecule ligand embedded between the heavy and
light chains. The chain identifiers of the heavy and light chains are H and L, and
the small molecule of interest has the residue identifier `HEP` and is part of chain
H.

The simplest analysis we can do is to include both chains

```
prodigy_lig.py -c H,L H:HEP -i 1A0Q.pdb
```

which produces the following output.

```
Job name        DGprediction (low refinement) (Kcal/mol)
1a0q    -8.49
```

In this case, atomic distances from the ligand to residues of both chains have
been calculated. **This only includes atoms of the protein.** Any cofactors,
ions or solvent molecules that might be present, are not included in the distances that are calculated.

If we wanted to only include the heavy chain (`H`) in the analysis we would use
this command instead.


```
prodigy_lig.py -c H H:HEP -i 1A0Q.pdb
```

which would produce this output

```
Job name        DGprediction (low refinement) (Kcal/mol)
1a0q    -6.69
```

### Optional arguments

The first two optional arguments will affect your results, whereas the latter ones
will only impact the formatting and content of your output.

#### Electrostatics Energy

In addition to the contact information prodigy_lig can make use of the electrostatic
energy of the interaction as well. This refers to the intermolecular electrostatic
energy as calculated by [HADDOCK](http://milou.science.uu.nl/services/HADDOCK2.2).

If you know this energy because you have refined your complex through HADDOCK then
you can specify this using the `-e` flag (`--electrostatics`). Additionaly, if your
PDB file is coming from HADDOCK, `prodigy_lig` will automatically extract the relevant
information from the input file and make use of it.

#### Distance cutoff

This is the cutoff used when calculating the atomic contacts. By default it has a
value of 10.5 Angstrom, which was identified when training the model. You can modify
this with the `-d` flag (`--distance_cutoff`).

#### Processed output file

If you would like a copy of the structure that was used to calculate the atomic contacts you can use the `-o` flag (`--output_file`). This might be useful if you
have chosen to exclude some chains from the analysis. The file will be created in
the current working directory and its filename will be `input_name-processed.pdb`.

#### Verbosity

If you specify the `-v` flag (`--verbose`), in addition to the DG values / scores,
the output will include the contact counts and if defined electrostatics as well.
Once again using the same example as above.

```
prodigy_lig.py -c H,L H:HEP -i 1A0Q.pdb -v
```

Will produce the following output

| Job name | DGprediction (low refinement) (Kcal/mol) | CC | CN | CO | CX | NN | NO | NX | OO | OX | XX
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1A0Q | -8.49 | 1838 | 589 | 1010 | 124 | 27 | 175 | 30 | 132 | 24 | 0

#### Version

Finally you can use the `-V` flag (`--version`) to get the version of the code.
