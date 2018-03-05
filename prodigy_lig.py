#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Calculate the Binding Affinity score using the PRODIGY-LIG model

This script only requires one PDB file as input and expects that
the all-atom contact script lives somewhere in the PATH. Failing
that the user can provide the path to the executable.

Authors: Panagiotis Koukos, Anna Vangone, Joerg Schaarschmidt
"""

from __future__ import print_function
from os.path import basename, splitext
import sys
import argparse
import string
from subprocess import Popen, PIPE
from StringIO import StringIO
from collections import namedtuple

from Bio.PDB import PDBParser, FastMMCIFParser, PDBIO


class ProdigyLig(object):
    """Run the prodigy-lig calculations and store all the relevant output."""
    def __init__(self, structure, chains, electrostatics, cpp_contacts, cutoff=10.5):
        """Initialise the Prodigy-lig instance."""
        self.chains = self._parse_chains(chains)
        self.structure = self._clean_structure(structure)
        self.electrostatics = electrostatics
        self.cpp_contacts = cpp_contacts
        self.cutoff = cutoff
        self.dg_score = None
        self.dg_elec = None
        self.dg = None
        self.contact_counts = None

    def predict(self):
        """
        API method used by the webserver
        """
        if self.cpp_contacts is not None:
            self.cpp_contacts = self.cpp_contacts[0]
            atomic_contacts = calc_atomic_contacts_cpp(self.cpp_contacts, self.structure, self.chains, self.cutoff)
        else:
            atomic_contacts = calc_atomic_contacts_python(self.structure, self.chains, self.cutoff)

        if len(atomic_contacts) == 0:
            raise RuntimeWarning(
                "There are no contacts between the specified chains."
            )

        self.contact_counts = calculate_contact_counts(atomic_contacts)

        if self.electrostatics is not None:
            self.dg_score = calculate_score(self.contact_counts, self.electrostatics)
            self.dg_elec = calculate_DG_electrostatics(self.contact_counts, self.electrostatics)
        self.dg = calculate_DG(self.contact_counts)

    def as_dict(self):
        """Return the data of the class as a dictionary for the server."""
        data = {
            'structure': self.structure.id,
            'chains': self.chains,
            'electrostatics': self.electrostatics,
            'cutoff': self.cutoff,
            'dg_score': self.dg_score,
            'dg_elec': self.dg_elec,
            'dg': self.dg
        }
        
        data.update(self.contact_counts)

        return data 

    def _clean_structure(self, structure):
        """
        Remove all the unnecessary elements from the structure.

        Water molecules, ions and cofactors need to be removed from the structure
        because otherwise their presence might affect the algorithm. In the case
        of multi-model structures we are only keeping the first model.
        """
        if len(structure) > 1:
            for i in xrange(1, len(structure)):
                structure.detach_child(structure[i].id)

        specified_chains = [chain for group in self.chains for chain in group]
        structure_chains = [chain.id for chain in list(structure.get_chains())]

        for chain in specified_chains:
            if chain not in structure_chains:
                raise RuntimeWarning(
                    "Chain {} specified during runtime wasn't found in "
                    "the structure".format(chain)
                )
        
        for chain in structure_chains:
            if chain not in specified_chains:
                structure[0].detach_child(chain)

        return structure

    @staticmethod
    def _parse_chains(chains):
        """
        Parse the chain and return a list of lists.

        The chains specification allows for one chain per interactor or more. If more
        than one chains per interactor are specified split on ',' and return the chains
        as a list.
        """
        def validate_chain_string(chain_string):
            """
            Check the chain string for any character other than letters and commas.
            """
            chain_string = chain_string.upper()
            try:
                chain_string.decode("ascii")
            except UnicodeDecodeError:
                raise RuntimeError(
                    "Please use uppercase ASCII characters [ A-Z ]."
                )

            if len(chain_string) == 1:
                if chain_string not in string.uppercase:
                    raise RuntimeError(
                        "Please use standard chain identifiers [ A-Z ]."
                    )
            else:
                for char in chain_string:
                    if char not in string.uppercase and char != ",":
                        raise RuntimeError(
                            "Use uppercase ASCII characters [ A-Z ] to speciy the"
                            " chains and , to separate them."
                        )

                # Make sure that "A," or "A,B," didn't slip through
                comma_count = chain_string.count(",")
                chain_count = len(set(chain_string).intersection(string.uppercase))

                if comma_count != chain_count -1:
                    raise RuntimeError(
                        "Specify multiple chains like this: prodigy_lig.py -c A,B C"
                    )

            return chain_string

        parsed_chains = []
        for chain in chains:
            try:
                chain = validate_chain_string(chain)
                parsed_chains.append(chain.split(","))
            except RuntimeError:
                raise
        return parsed_chains

    def print_prediction(self, outfile=''):
        if outfile:
            handle = open(outfile, 'w')
        else:
            handle = sys.stdout
        """Print to the File or STDOUT if no filename is specified."""
        if self.electrostatics is not None:
            handle.write("{}\t{}\t{}\n".format("Job name", "DGprediction (Kcal/mol)", "DGscore"))
            handle.write("{0}\t{1:.2f}\t{2:.2f}\n".format(self.structure.id, self.dg_elec, self.dg_score))
        else:
            handle.write("{}\t{}\n".format("Job name", "DGprediction (low refinement) (Kcal/mol)"))
            handle.write("{0}\t{1:.2f}\n".format(self.structure.id, self.dg))
        if handle is not sys.stdout:
            handle.close()


def extract_electrostatics(pdb_file):
    """
    Extracts the electrostatics energy from a HADDOCK PDB file.

    :param pdb_file: The input PDB file.
    :return: Electrostatics energy
    """
    electrostatics = None
    for line in pdb_file:
        # scan for Haddock energies line and assign electrostatics
        if line.startswith('REMARK energies'):
            line = line.rstrip()
            line = line.replace('REMARK energies: ', '')
            electrostatics = float(line.split(',')[6])
            break
        # stop on first ATOM Line as remarks should be beforehand
        elif line.startswith('ATOM'):
            break
    # try to reset file handle for further processing
    try:
        pdb_file.seek(0)
    except Exception:
        pass
    return electrostatics


def calc_atomic_contacts_cpp(contact_executable, pdb_file, chains, cutoff=10.5):
    """
    Calculate atomic contacts.

    This will call out to the executable defined during startup time and
    collect its output. After processing it will return a list of the
    atomic contacts.

    :param contact_executable: Path to the all-atom contact script
    :type contact_executable: str or unicode
    :param pdb_file: The structure object
    :type pdb_file: Bio.PDB structure object
    :param cutoff: The cutoff to use for the AC calculation
    :type cutoff: float
    :return: Str of atomic contacts
    """
    def _filter_contacts_by_chain(contacts, chains):
        """
        Filter the contacts using only the chains specified during runtime.
        """
        filtered_contacts = []

        for contact in contacts:
            words = contact.split()
            chain1 = words[1].upper()
            chain2 = words[5].upper()

            chains_are_acceptable = (
                (chain1 in chains[0] and chain2 in chains[1]) or
                (chain1 in chains[1] and chain2 in chains[0])
            )

            if chains_are_acceptable:
                filtered_contacts.append(contact)

        return filtered_contacts

    io = PDBIO()
    io.set_structure(pdb_file)
    io_stream = StringIO()
    io.save(io_stream)
    io_stream.seek(0)

    p = Popen([contact_executable, str(cutoff)], stdin=PIPE, stdout=PIPE, stderr=PIPE)
    for line in io_stream:
        p.stdin.write(line)
    p.stdin.close()
    atomic_contacts = p.stdout.readlines()

    del atomic_contacts[-1]

    return _filter_contacts_by_chain(atomic_contacts, chains)


def calc_atomic_contacts_python(structure, chains, cutoff=10.5):
    """
    Calculate the contacts without calling out to the CPP code.

    :param structure: Biopython structure object of the input file
    :return: List of contacts
    """
    # Ignore multi model structures
    interactors = [[], []]
    structure = structure[0]

    for group in chains:
        for chain in group:
            chain_index = [chain in group for group in chains].index(True)
            interactors[chain_index].append(structure[chain])

    contacts = []
    for protein_chain in interactors[0]:
        for ligand_chain in interactors[1]:
            for protein_atom in protein_chain.get_atoms():
                for ligand_atom in ligand_chain.get_atoms():
                    dist = protein_atom - ligand_atom
                    if dist <= cutoff:
                        contacts.append("\t".join([
                            protein_atom.parent.resname,
                            protein_atom.parent.parent.id,
                            str(protein_atom.parent.id[1]),
                            protein_atom.element,
                            ligand_atom.parent.resname,
                            ligand_atom.parent.parent.id,
                            str(ligand_atom.parent.id[1]),
                            ligand_atom.element,
                            str(dist)
                        ]))

    return contacts


def calculate_contact_counts(contacts):
    """
    Calculate the counts of the various atomic contact types based on the
    types of atoms that are in contact. The categories are:

    CC: Carbon-Carbon
    NN: Nitrogen-Nitrogen
    OO: Oxygen-Oxygen
    XX: Other-Other

    and the combinations: CN, CO, CX, NO, NX, OX

    :param contacts: The output of the calc_atomic_contacts functions
    :return: dict of the counts of each category defined above
    """
    def _classify_atom(atom):
        """
        Classify the atom involved in the interaction in one of the categories
        laid out in calculate_atomic_contacts.

        :param atom: The atom involved in the interaction
        :return: Atom type. One of C, N, O, X
        """
        if atom == 'C' or atom == 'N' or atom == 'O':
            return atom
        else:
            return 'X'

        return None

    def _classify_contact(atom_classes):
        """
        Classify the contact in one of the categories defined in the function
        calculate_atomic_contact_counts.

        :param atom_classes: Class of the atoms involved in the interaction
        :type atom_classes: List of length 2
        :return: One of CC, NN, OO, XX, CN, CO, CX, NO, NX, OX
        """
        atom_1, atom_2 = atom_classes
        if atom_1 == 'C' and atom_2 == 'C':
            return 'CC'
        elif (atom_1 == 'C' and atom_2 == 'N') or (atom_1 == 'N' and atom_2 == 'C'):
            return 'CN'
        elif (atom_1 == 'C' and atom_2 == 'O') or (atom_1 == 'O' and atom_2 == 'C'):
            return 'CO'
        elif (atom_1 == 'C' and atom_2 == 'X') or (atom_1 == 'X' and atom_2 == 'C'):
            return 'CX'
        elif (atom_1 == 'N' and atom_2 == 'N') or (atom_1 == 'N' and atom_2 == 'N'):
            return 'NN'
        elif (atom_1 == 'N' and atom_2 == 'O') or (atom_1 == 'O' and atom_2 == 'N'):
            return 'NO'
        elif (atom_1 == 'N' and atom_2 == 'X') or (atom_1 == 'X' and atom_2 == 'N'):
            return 'NX'
        elif (atom_1 == 'O' and atom_2 == 'O') or (atom_1 == 'O' and atom_2 == 'O'):
            return 'OO'
        elif (atom_1 == 'O' and atom_2 == 'X') or (atom_1 == 'X' and atom_2 == 'O'):
            return 'OX'
        elif (atom_1 == 'X' and atom_2 == 'X') or (atom_1 == 'X' and atom_2 == 'X'):
            return 'XX'
        else:
            return None

    counts = {
        'CC': 0,
        'NN': 0,
        'OO': 0,
        'XX': 0,
        'CN': 0,
        'CO': 0,
        'CX': 0,
        'NO': 0,
        'NX': 0,
        'OX': 0
    }

    allowed_atoms = set(['C', 'N', 'O', 'F', 'CL', 'BR', 'S', 'P'])

    for line in contacts:
        if len(line) == 0:
            continue
        words = line.split()

        atom_name_1 = words[3]
        atom_name_2 = words[7]

        if not (atom_name_1 in allowed_atoms and atom_name_2 in allowed_atoms):
            continue

        atom_class_1 = _classify_atom(atom_name_1)
        atom_class_2 = _classify_atom(atom_name_2)

        contact_class = _classify_contact([atom_class_1, atom_class_2])
        counts[contact_class] += 1

    return counts


def calculate_score(contact_counts, electrostatics_energy):
    """
    Calculates the PRODIGY-lig score based on the contact counts and the
    electrostatics energy.

    :param contact_counts: Counts of the CC, NN, OO, XX contacts
    :type contact_counts: dict
    :param electrostatics_energy: Electrostatics energy calculated by HADDOCK
    :type electrostatics_energy: float
    :return: The PRODIGY-lig score
    """
    elec_weight = 0.343794
    cc_weight = -0.037597
    nn_weight = 0.138738
    oo_weight = 0.160043
    xx_weight = -3.088861
    intercept = 187.011384

    return (
        (elec_weight * electrostatics_energy) +
        (cc_weight * contact_counts['CC']) +
        (nn_weight * contact_counts['NN']) +
        (oo_weight * contact_counts['OO']) +
        (xx_weight * contact_counts['XX']) +
        intercept
    )


def calculate_DG(contact_counts):
    """
    Calculates the PRODIGY-lig binding affinity using the weigths that
    have been trained for the prediction without electrostatics.
    """
    nn_weight = 0.0354707
    xx_weight = -0.1277895
    cn_weight = -0.0072166
    intercept = -5.1923181

    return (
        (nn_weight * contact_counts['NN']) +
        (xx_weight * contact_counts['XX']) +
        (cn_weight * contact_counts['CN']) +
        intercept
    )


def calculate_DG_electrostatics(contact_counts, electrostatics_energy):
    """
    Calculates the PRODIGY-lig binding affinity using the weights that
    have been optimised for the prediction with electrostatics.
    """
    elec_weight = 0.0115148
    cc_weight = -0.0014852
    nn_weight = 0.0057097
    xx_weight = -0.1301806
    intercept = -5.1002233

    return (
        (elec_weight * electrostatics_energy) +
        (cc_weight * contact_counts['CC']) +
        (nn_weight * contact_counts['NN']) +
        (xx_weight * contact_counts['XX']) +
        intercept
    )


def _parse_arguments():
    """Parse the command line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        '--cpp_contacts',
        required=False,
        nargs=1,
        help=(
            'Path to the all-atom contact script. By default prodigy_lig.py'
            ' will calculate the distances using a built-in implementation. In'
            ' some cases there might be a speed-up when using external compiled'
            ' code.'
        )
    )
    parser.add_argument(
        '-c',
        '--chains',
        required=True,
        nargs=2,
        help=(
            'Which chains to use. You can specify multi-chain selections'
            ' by comma separating the chain identifiers (e.g. -c A,B C,D).'
            ' In that case only contacts between chains A - C+D and B - C+D'
            ' will be considered.'
        )
    )
    parser.add_argument(
        '-i',
        '--input_file',
        required=True,
        help='This is the PDB/mmcif file for which the score will be calculated.'
    )
    parser.add_argument(
        '-e',
        '--electrostatics',
        required=False,
        type=float,
        help=u'This is the electrostatics energy as calculated during the'
             u' water refinement stage of HADDOCK.'
    )
    parser.add_argument(
        '-d',
        '--distance_cutoff',
        required=False,
        default=10.5,
        help=u'This is the distance cutoff for the Atomic Contacts '
             u' (def = 10.5Å).'
    )

    return parser.parse_args()


def main():
    """Run it."""
    args = _parse_arguments()
    fname, s_ext = splitext(basename(args.input_file))
    parser = None
    if s_ext in {'.pdb', '.ent'}:
        parser = PDBParser(QUIET=1)
    elif s_ext == ".cif":
        parser = FastMMCIFParser(QUIET=1)

    with open(args.input_file) as in_file:
        # try to set electrostatics from input file if not provided by user
        electrostatics = args.electrostatics \
            if args.electrostatics or s_ext == '.cif' \
            else extract_electrostatics(in_file)
        prodigy_lig = ProdigyLig(
            parser.get_structure(fname, in_file),
            chains=args.chains,
            electrostatics=electrostatics,
            cpp_contacts=args.cpp_contacts,
            cutoff=args.distance_cutoff
        )

    prodigy_lig.predict()
    prodigy_lig.print_prediction()


if __name__ == "__main__":
    main()
