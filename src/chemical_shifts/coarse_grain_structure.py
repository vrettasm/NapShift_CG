# General imports
from collections import defaultdict, OrderedDict

# Private import(s).
from src.chemical_shifts.auxiliaries import ACCEPTED_RES
from src.chemical_shifts.coarse_grain_entity import CoarseGrainEntity
from Bio.PDB.vectors import Vector, calc_dihedral, calc_angle


class CoarseGrainStructure(object):
    """
    This class creates a structure for the coarse-grain protein model
    (MARTINI). It reads a MARTINI  PDB file and computes the dihedral
    angles for each residue.
    """

    # Object variables.
    __slots__ = ("_item_dict", "residue_dict", "sequence_dict")

    # Constructor.
    def __init__(self):
        """
        Initializes the member dictionaries with default values.
        """

        # The items are hidden because
        # we only process them internally.
        self._item_dict = defaultdict(list)

        # The residue dictionary can
        # be returned upon request.
        self.residue_dict = defaultdict(OrderedDict)

        # This is the amino-acid sequence of
        # the structure in three-letter code.
        self.sequence_dict = defaultdict(list)
    # _end_def_

    @property
    def get_structure(self):
        """
        Accessor (getter) of the residue structure.

        :return: a dictionary with the residues.
        """
        return self.residue_dict
    # _end_def_

    @property
    def get_sequence(self):
        """
        Accessor (getter) of the residue sequence.

        :return: a dictionary with the residues (three-letter code).
        """
        return self.sequence_dict
    # _end_def_

    # Auxiliary method.
    def parse_martini_pdb_file(self, f_path, f_type="MRT", all_chains=False, allow_overwrite=False):
        """
        PDB file record format.

        COLUMNS        DATA  TYPE    FIELD        DEFINITION
        -------------------------------------------------------------------------------------
         1 -  6        Record name   "ATOM  "
         7 - 11        Integer       Serial       Atom  serial number.
        13 - 16        Atom          Name         The atom name.
        17             Character     AltLoc       Alternate location indicator.
        18 - 20        Residue name  ResName      Residue name.
        22             Character     ChainID      Chain identifier.
        23 - 26        Integer       ResSeq       Residue sequence number.
        27             AChar         iCode        Code for insertion of residues.
        31 - 38        Real(8.3)     X            Orthogonal coordinates for X in Angstroms.
        39 - 46        Real(8.3)     Y            Orthogonal coordinates for Y in Angstroms.
        47 - 54        Real(8.3)     Z            Orthogonal coordinates for Z in Angstroms.
        55 - 60        Real(6.2)     Occupancy    The occupancy.
        61 - 66        Real(6.2)     TempFactor   Temperature  factor.
        77 - 78        LString(2)    Element      The element symbol, right-justified.
        79 - 80        LString(2)    Charge       The charge on the atom.
        -------------------------------------------------------------------------------------

        NOTE: The numbering of the columns starts from '1', not '0'.

        :param f_path: (Path) with the MARTINI pdb file.

        :param f_type: (string) extension for the file name.

        :param all_chains: (bool) if True it will process all chains in the protein.
        By default, it processes only the first (main) chain.

        :param allow_overwrite: (bool) flag that determines to overwrite protection.
        By default, it does not allow overwriting.

        :return: None.
        """

        # Check the file type.
        if f_type not in {"MRT", "PDB"}:
            raise ValueError(f"{self.__class__.__name__}: "
                             f"The acceptable file types are 'MRT' or 'PDB': {f_type}")
        # _end_if_

        # Check for overwrite protection.
        if not allow_overwrite:

            # Check if the list contains elements.
            if self._item_dict:

                # Exit with a message.
                raise RuntimeError(f"{self.__class__.__name__}: "
                                   f"Overwrite protection is: {allow_overwrite}. "
                                   f"Items list is not empty: {len(self._item_dict)}!")
            # _end_if_

        # _end_if_

        # Break list keywords.
        break_list = ["END", "ENDMDL"]

        # If we don't want all chains we update the list.
        if not all_chains:
            break_list.append("TER")
        # _end_if_

        # Make sure the item list is empty.
        self._item_dict.clear()

        # Get the file id.
        f_id = f_path.stem

        # Read the data from the file.
        with open(f_path, "r") as f_in:

            # Process line-by-line.
            for row in f_in.readlines():

                # Make sure we capture any exceptions.
                try:

                    # Get the first entry in the line.
                    entry = str(row[0:6]).strip().upper()

                    # Break when you reach the end of the model.
                    if entry in break_list:
                        # Exit the loop.
                        break
                    else:

                        # Process only 'atom' entries.
                        if entry == "ATOM":

                            # Get the chain id separately.
                            chain_id = str(row[21]).upper()

                            # Create a new row for the list.
                            # NOTE: The "key" names should be identical
                            # to the CoarseGrainEntity input fields __init__.
                            record = {"num_id": int(row[6:11]),
                                      "atom_type": str(row[12:16]).strip().upper(),
                                      "res_name": str(row[17:20]).strip().upper(),
                                      "chain_id": str(row[21]).upper(),
                                      "res_id": int(row[22:26]),
                                      "res_x": float(row[30:38]),
                                      "res_y": float(row[38:46]),
                                      "res_z": float(row[46:54])}

                            # Update the list.
                            self._item_dict[chain_id].append(CoarseGrainEntity(**record))
                        else:

                            # Skip other entries.
                            continue
                        # _end_if_

                    # _end_if_
                except ValueError as val_err:

                    # Show a warning message if some conversion fails.
                    print(f"{self.__class__.__name__}: "
                          f"WARNING: Method failed in {f_id}: {val_err}")

                    # If methods 'int' and 'float' can't convert the strings
                    # they will throw a ValueError. In this case skip and go
                    # to the next row.
                    continue
                # _end_try_

            # _end_for_

        # _end_with_

        # Iterate through all the list items and group them by "res_id".
        # This way inside each dictionary entry will have all the reduced
        # entities related to that residue.
        for chain_id in self._item_dict:

            # Create a temporary dictionary with dictionary values.
            tmp_dict = defaultdict(dict)

            for it in self._item_dict[chain_id]:

                # Use the "res_id" as key.
                tmp_dict[it.res_id][it.atom_type] = it
            # _end_for_

            # Sort them by "key" value. The sorted method returns a
            # sorted copy of the input dictionary to avoid leaks.
            self.residue_dict[chain_id] = OrderedDict(sorted(tmp_dict.items()))
        # _end_if_

    # _end_def_

    # Auxiliary method.
    def compute_dihedral_angles_and_more(self):
        """
        Computes the dihedral angles, along with the other angles,
        for all the chains and for each residue in the item lists.

        :return: None.
        """

        # Make sure the dictionary is not empty.
        if not self.residue_dict:

            # Display a message to the user.
            print(f"{self.__class__.__name__}: "
                  f"Residue dict is empty, nothing to compute.")

            # Exit here.
            return None
        # _end_if_

        # Process all chain in the protein.
        for chain_id in self.residue_dict:

            # Make a temporary list with the residues.
            res_list = list(self.residue_dict[chain_id].values())

            # NOTE: Since Python allows negative indexes we have to make
            # sure that we will not accidentally take the wrong atoms in
            # the calculations.

            # The length of the list is used only for bound checking.
            L = len(res_list)

            # Make sure the list is clear before processing.
            self.sequence_dict[chain_id].clear()

            # Calculate the dihedral angles for all residues
            # in the dictionary.
            for i, res_i in enumerate(res_list, start=0):

                # Residue name (three-letter code).
                RES_NAME = str(res_i["BB"].res_name).strip()

                # Accepted list of residues.
                if RES_NAME in ACCEPTED_RES:
                    self.sequence_dict[chain_id].append(RES_NAME)
                else:
                    # Show a warning. This is temporary until
                    # we decide what to do if there are gaps.
                    print(f" Unknown residue name: {RES_NAME}")

                    # Skip to the next one.
                    continue
                # _end_if_

                # Temporary 'BB' atoms for the dihedral.
                tmp_BB = []

                # Left dihedral angle.
                for left_i in (i - 2, i - 1, i, i + 1):

                    # Bounds checking.
                    if left_i < 0 or left_i > (L - 1):
                        # Cancel the calculation since we
                        # don't have all the atoms needed.
                        tmp_BB.clear()

                        # Break the local loop.
                        break
                    # _end_if_

                    # Copy the 'BB' atom coordinates (as Vector object).
                    tmp_BB.append(Vector(res_list[left_i]["BB"].get_coord))
                # _end_for_

                # We need exactly four points to calculate
                # the left dihedral angle.
                if len(tmp_BB) == 4:
                    res_i["BB"].set_angle("left_dihedral", calc_dihedral(*tmp_BB))
                # _end_if_

                # Clear the list for the right dihedral.
                tmp_BB.clear()

                # Right dihedral angle.
                for right_i in (i - 1, i, i + 1, i + 2):

                    # Bounds checking.
                    if right_i < 0 or right_i > (L - 1):
                        # Cancel the calculation since we
                        # don't have all the atoms needed.
                        tmp_BB.clear()

                        # Break the local loop.
                        break
                    # _end_if_

                    # Copy the 'BB' atom coordinates (as Vector object).
                    tmp_BB.append(Vector(res_list[right_i]["BB"].get_coord))
                # _end_for_

                # We need exactly four points to calculate
                # the right dihedral angle.
                if len(tmp_BB) == 4:
                    res_i["BB"].set_angle("right_dihedral", calc_dihedral(*tmp_BB))
                # _end_if_

                # Temporary atoms for the angles.
                tmp_atoms = []

                # Angle 'alpha'.
                for a_i in (i - 1, i):

                    # Bounds checking.
                    if a_i < 0 or a_i > (L - 1):
                        # Cancel the calculation since we
                        # don't have all the atoms needed.
                        tmp_atoms.clear()

                        # Break the local loop.
                        break
                    # _end_if_

                    # Copy the 'BB' coordinates (as Vector object).
                    tmp_atoms.append(Vector(res_list[a_i]["BB"].get_coord))

                    # Copy the 'SC1' in the i-th position, if it exists.
                    if a_i == i and ("SC1" in res_list[a_i]):
                        tmp_atoms.append(Vector(res_list[a_i]["SC1"].get_coord))
                    # _end_if_

                # _end_for_

                # If the temp list is not empty calculate
                # the 'alpha' angle. We need three points.
                if len(tmp_atoms) == 3:
                    res_i["BB"].set_angle("alpha", calc_angle(*tmp_atoms))
                # _end_if_

                # Clear the list for the next angle.
                tmp_atoms.clear()

                # Angle 'beta'.
                for b_i in (i - 1, i, i + 1):

                    # Bounds checking.
                    if b_i < 0 or b_i > (L - 1):
                        # Cancel the calculation since we
                        # don't have all the atoms needed.
                        tmp_atoms.clear()

                        # Break the local loop.
                        break
                    # _end_if_

                    # Copy the 'BB' coordinates (as Vector object).
                    tmp_atoms.append(Vector(res_list[b_i]["BB"].get_coord))
                # _end_for_

                # If the temp list is not empty calculate
                # the 'beta' angle. We need three points.
                if len(tmp_atoms) == 3:
                    res_i["BB"].set_angle("beta", calc_angle(*tmp_atoms))
                # _end_if_

                # Clear the list for the next angle.
                tmp_atoms.clear()

                # Angle 'gamma'.
                for c_i in (i, i + 1):

                    # Bounds checking.
                    if c_i < 0 or c_i > (L - 1):
                        # Cancel the calculation since we
                        # don't have all the atoms needed.
                        tmp_atoms.clear()

                        # Break the local loop.
                        break
                    # _end_if_

                    # Copy the 'SC1' in the i-th position, if it exists.
                    if c_i == i and ("SC1" in res_list[c_i]):
                        tmp_atoms.append(Vector(res_list[c_i]["SC1"].get_coord))
                    # _end_if_

                    # Copy the 'BB' coordinates (as Vector object).
                    tmp_atoms.append(Vector(res_list[c_i]["BB"].get_coord))
                # _end_for_

                # If the temp list is not empty calculate
                # the 'gamma' angle. We need three points.
                if len(tmp_atoms) == 3:
                    res_i["BB"].set_angle("gamma", calc_angle(*tmp_atoms))
                # _end_if_

            # _end_for_

        # _end_for_

    # _end_def_

# _end_class_
