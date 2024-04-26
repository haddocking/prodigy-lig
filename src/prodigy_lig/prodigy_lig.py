"""
Calculate the Binding Affinity score using the PRODIGY-LIG model
prodigy_lig dependes on biopython for the structure manipulations
and only requires a single structure file (in mmCIF or PDB format)
as input.

prodigy_lig is licensed under the Apache License 2.0 included in
the LICENSE file of this repository or at the following URL
https://github.com/haddocking/prodigy-lig/blob/master/LICENSE

If you use prodigy_lig in your research please cite the following
papers:

1. https://doi.org/10.1093/bioinformatics/bty816
2. https://doi.org/10.1007/s10822-017-0049-y
"""

from os.path import basename, splitext
import sys
import argparse

try:
    from Bio.PDB import PDBParser, FastMMCIFParser, PDBIO
except ImportError:
    print(
        "prodigy_lig depends on Biopython. Please install "
        "Biopython along with its dependencies and restart."
    )
    sys.exit(1)


class ProdigyLig(object):
    """Run the prodigy-lig calculations and store all the relevant output."""
    def __init__(self, structure, chains, electrostatics, cutoff=10.5):
        """Initialise the Prodigy-lig instance."""
        self.chains = self._parse_chains(chains)
        self.structure = self._clean_structure(structure)
        self.electrostatics = electrostatics
        self.cutoff = cutoff
        self.dg_score = None
        self.dg_elec = None
        self.dg = None
        self.atomic_contacts = None
        self.contact_counts = None

    def predict(self):
        """
        API method used by the webserver
        """
        self.atomic_contacts = calc_atomic_contacts(self.structure, self.chains, self.cutoff)

        if len(self.atomic_contacts) == 0:
            raise RuntimeWarning(
                "There are no contacts between the specified chains."
            )

        self.contact_counts = calculate_contact_counts(self.atomic_contacts)

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
        def _is_it_a_residue(residue):
            """
            Check for the presence of backbone atoms.

            Pretty often PDB files contain modified residues that are part of
            a protein but are classified as HETATM, such as a selenomethionine
            (MSE). We want to keep these atoms instead of discarding them.
            """
            backbone_atoms = {"C", "CA", "N", "O"}
            residue_atoms = set([_.id for _ in residue.child_list])

            if len(residue_atoms.intersection(backbone_atoms)) >= 3:
                return True
            else:
                return False

        if len(structure) > 1:
            for i in range(1, len(structure)):
                structure.detach_child(structure[i].id)

        chains = self.chains[0] + list(self.chains[1][0])
        specified_chains = {chain for chain in chains}
        structure_chains = [chain.id for chain in list(structure.get_chains())]

        for chain in specified_chains:
            if chain not in structure_chains:
                raise RuntimeWarning(
                    ('Chain {} specified during runtime wasn\'t found '
                     'in the structure').format(chain)
                )

        ligand_chain, ligand_residue = self.chains[1]
        ligand_residues = [
            _.resname for _ in structure[0][ligand_chain].child_list
        ]

        if ligand_residue not in ligand_residues:
            raise RuntimeError(
                ('Ligand identifier {} not found '
                 'in chain {} of the input file.').format(
                    ligand_residue,
                    ligand_chain
                )
            )

        if ligand_residues.count(ligand_residue) > 1:
            raise RuntimeError(
                ('Ligand identifier {} found {} times in '
                 'chain {} of the input file. Please '
                 'remove the redundant molecule(s).').format(
                    ligand_residue,
                    ligand_residues.count(ligand_residue),
                    ligand_chain
                )
            )

        for chain in structure_chains:
            if chain not in specified_chains:
                structure[0].detach_child(chain)
            else:
                for res in list(structure[0][chain].child_list):
                    if res.id[0] == "W":
                        structure[0][chain].detach_child(res.id)
                        continue
                    if chain not in self.chains[0]:
                        # We are only interested in the ligand here.
                        # Everything else can go.
                        if res.resname != ligand_residue:
                            structure[0][chain].detach_child(res.id)
                    else:
                        if res.resname != ligand_residue:
                            if _is_it_a_residue(res) or res.id[0] == " ":
                                continue
                            else:
                                structure[0][chain].detach_child(res.id)

        return structure

    @staticmethod
    def _parse_chains(chains):
        """
        Parse the chain and return a list of lists.

        The chains specification requires one chain per interactor. The first
        argument following the -c flag is the specfication for the protein
        selection and more than one chain ids can be specified by comma separating
        them. The second command line argument following -c corresponds to the
        ligand specification and requires one chain id followed by the reside
        identifier of the ligand (e.g. -c A B:LIG). The parsed chains are returned
        as a list.
        """
        parsed_chains = []
        chains = [_.upper() for _ in chains]
        flattened_chains_string = "".join(chains)
        protein_chain_string, ligand_chain_string = [_.upper() for _ in chains]

        try:
            flattened_chains_string.encode("ascii")
        except UnicodeDecodeError:
            raise RuntimeError(
                'Please use uppercase ASCII characters [ A-Z ].'
            )

        for char in protein_chain_string:
            if not char.isalpha() and char != ",":
                raise RuntimeError(
                    'Use uppercase ASCII characters [ A-Z ] to specify the'
                    ' chains and , to separate them.'
                )

        # Make sure that "A," or "A,B," didn't slip through
        comma_count = protein_chain_string.count(",")

        if 2 * comma_count + 1 != len(protein_chain_string):
            raise RuntimeError(
                'Specify multiple chains like this: prodigy_lig.py -c A,B C'
            )

        parsed_chains.append(protein_chain_string.split(","))

        ligand_specification = (
            'Use uppercase ASCII characters [ A-Z ] to speciy the '
            'chain and ligand identifiers. Separate the ligand '
            'identifier from its chain by colon. prodigy_lig.py '
            '-c A B:LIG'
        )
        if len(ligand_chain_string) != 5:
            raise RuntimeError(ligand_specification)

        for char in ligand_chain_string:
            if not char.isalnum() and char != ":":
                raise RuntimeError(ligand_specification)

        if ligand_chain_string.count(":") != 1:
            raise RuntimeError(ligand_specification)

        chain, ligand = ligand_chain_string.split(":")
        if len(chain) != 1:
            raise RuntimeError(ligand_specification)

        if len(ligand) != 3:
            raise RuntimeError(ligand_specification)

        parsed_chains.append(ligand_chain_string.split(":"))

        return parsed_chains

    def print_contacts(self, outfile=''):
        """Print the atomic contacts to STDOUT or the specified handle."""
        if outfile:
            handle = open(outfile, 'w')
        else:
            handle = sys.stdout

        for line in self.atomic_contacts:
            handle.write("{}\n".format(line))

        if handle is not sys.stdout:
            handle.close()

    def print_prediction(self, outfile='', verbose=False):
        """Print to the File or STDOUT if no filename is specified."""
        if outfile:
            handle = open(outfile, 'w')
        else:
            handle = sys.stdout
        if self.electrostatics is not None:
            header = ["Job name", "DGprediction (Kcal/mol)", "DGscore"]
            values = [self.structure.id, self.dg_elec, self.dg_score]
            if verbose:
                header.append("Electrostatics Energy")
                values.append(self.electrostatics)
                for key, value in sorted(self.contact_counts.items()):
                    header.append(key)
                    values.append(value)
                handle.write("\t".join(header) + "\n")
                handle.write(("{}" + "\t{:.2f}" * 3 + "\t{:4d}" * 10 + "\n").format(*values))
            else:
                handle.write("\t".join(header) + "\n")
                handle.write("{}\t{:.2f}\t{:.2f}\n".format(*values))
        else:
            header = ["Job name", "DGprediction (low refinement) (Kcal/mol)"]
            values = [self.structure.id, self.dg]
            if verbose:
                for key, value in sorted(self.contact_counts.items()):
                    header.append(key)
                    values.append(value)
                handle.write("\t".join(header) + "\n")
                handle.write(("{}\t{:.2f}" + "\t{:4d}" * 10 + "\n").format(*values))
            else:
                handle.write("\t".join(header) + "\n")
                handle.write("{}\t{:.2f}\n".format(*values))
        if handle is not sys.stdout:
            handle.close()

    def print_structure(self, outfile=''):
        """Store the processed structure in a file."""
        io = PDBIO()
        io.set_structure(self.structure)
        io.save(outfile)


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


def calc_atomic_contacts(structure, chains, cutoff=10.5):
    """
    Calculate the contacts without calling out to the CPP code.

    :param structure: Biopython structure object of the input file
    :param chains: the interactor selection as array (e.g [ ['A', 'B'], ['A','LIG'] ])
    :param cutoff: the distance cutoff to determine interactions in Angstrom
    :return: List of contacts
    """
    contacts = []
    structure = structure[0]
    prot_chains = {_ for _ in chains[0] + list(chains[1][0])}
    lig_chain, lig_res = chains[1]

    for lig_atom in structure[lig_chain].get_atoms():
        if lig_atom.parent.resname == lig_res:
            for prot_chain in prot_chains:
                for prot_atom in structure[prot_chain].get_atoms():
                    if prot_atom.parent.resname != lig_atom.parent.resname:
                        dist = lig_atom - prot_atom
                        if dist <= cutoff:
                            contacts.append("\t".join([
                                prot_atom.parent.resname,
                                prot_atom.parent.parent.id,
                                str(prot_atom.parent.id[1]),
                                prot_atom.element,
                                lig_atom.parent.resname,
                                lig_atom.parent.parent.id,
                                str(lig_atom.parent.id[1]),
                                lig_atom.element,
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

    allowed_atoms = {'C', 'N', 'O', 'F', 'CL', 'BR', 'S', 'P'}

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
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='Authors: Panagiotis Koukos, Anna Vangone, Joerg Schaarschmidt'
    )

    parser.add_argument(
        '-c',
        '--chains',
        required=True,
        nargs=2,
        help=(
            'Which chains to use. Expects two sets of arguments. The first '
            'set refers to the protein selection and multiple chains can be'
            ' specified by separating the chain identifiers with commas. The'
            ' second set refers to the ligand and requires one chain and the'
            ' residue identifier of the ligand. A typical use case could be '
            'the following: prodigy_lig.py -c A,B A:LIG'
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
        help='This is the electrostatics energy as calculated during the'
             ' water refinement stage of HADDOCK.'
    )
    parser.add_argument(
        '-d',
        '--distance_cutoff',
        required=False,
        default=10.5,
        type=float,
        help='This is the distance cutoff for the Atomic Contacts '
             ' (def = 10.5Å).'
    )
    parser.add_argument(
        '-o',
        '--output_file',
        default=None,
        action='store_true',
        required=False,
        help='Store the processed file. The filename will be the name '
             'of the input with -processed appended just before the '
             'file ending.'
    )
    parser.add_argument(
        '-v',
        '--verbose',
        default=False,
        action='store_true',
        required=False,
        help='Include the calculated contact counts in the output.'
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
            cutoff=args.distance_cutoff
        )

    prodigy_lig.predict()
    prodigy_lig.print_prediction('', args.verbose)

    if args.output_file is not None:
        output_file_name = splitext(prodigy_lig.structure.id)[0]
        output_file_name += "-processed.pdb"
        prodigy_lig.print_structure(output_file_name)


if __name__ == "__main__":
    main()
