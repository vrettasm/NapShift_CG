"""
This module includes the main class that handles
the data generation from the (input) PDB files.
"""

# Python import(s).
import numpy as np
from numba import njit
from pathlib import Path
from pandas import DataFrame

# Private import(s).
from src.chemical_shifts.coarse_grain_structure import CoarseGrainStructure
from src.chemical_shifts.auxiliaries import (RES_3_TO_1, TARGET_ATOMS, ACCEPTED_RES)


class CoarseGrainInputVector(object):
    """
    This class creates an input vector for the artificial neural network.
    """

    # Object variables.
    __slots__ = ("_blosum",)

    # Constructor.
    def __init__(self):
        """
        Constructs an object that will create the input vectors from
        a given PDB (protein) file. The blosum field is kept hidden
        since it is only used internally.
        """
        self._blosum = {"ALA": (+4, -1, -2, -2, +0, -1, -1, +0, -2, -1, -1, -1, -1, -2, -1, +1, +0, -3, -2, +0),
                        "ARG": (-1, +5, +0, -2, -3, +1, +0, -2, +0, -3, -2, +2, -1, -3, -2, -1, -1, -3, -2, -3),
                        "ASN": (-2, +0, +6, +1, -3, +0, +0, +0, +1, -3, -3, +0, -2, -3, -2, +1, +0, -4, -2, -3),
                        "ASP": (-2, -2, +1, +6, -3, +0, +2, -1, -1, -3, -4, -1, -3, -3, -1, +0, -1, -4, -3, -3),
                        "CYS": (+0, -3, -3, -3, +9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1),
                        "GLN": (-1, +1, +0, +0, -3, +5, +2, -2, +0, -3, -2, +1, +0, -3, -1, +0, -1, -2, -1, -2),
                        "GLU": (-1, +0, +0, +2, -4, +2, +5, -2, +0, -3, -3, +1, -2, -3, -1, +0, -1, -3, -2, -2),
                        "GLY": (+0, -2, +0, -1, -3, -2, -2, +6, -2, -4, -4, -2, -3, -3, -2, +0, -2, -2, -3, -3),
                        "HIS": (-2, +0, +1, -1, -3, +0, +0, -2, +8, -3, -3, -1, -2, -1, -2, -1, -2, -2, +2, -3),
                        "ILE": (-1, -3, -3, -3, -1, -3, -3, -4, -3, +4, +2, -3, +1, +0, -3, -2, -1, -3, -1, +3),
                        "LEU": (-1, -2, -3, -4, -1, -2, -3, -4, -3, +2, +4, -2, +2, +0, -3, -2, -1, -2, -1, +1),
                        "LYS": (-1, +2, +0, -1, -3, +1, +1, -2, -1, -3, -2, +5, -1, -3, -1, +0, -1, -3, -2, -2),
                        "MET": (-1, -1, -2, -3, -1, +0, -2, -3, -2, +1, +2, -1, +5, +0, -2, -1, -1, -1, -1, +1),
                        "PHE": (-2, -3, -3, -3, -2, -3, -3, -3, -1, +0, +0, -3, +0, +6, -4, -2, -2, +1, +3, -1),
                        "PRO": (-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, +7, -1, -1, -4, -3, -2),
                        "SER": (+1, -1, +1, +0, -1, +0, +0, +0, -1, -2, -2, +0, -1, -2, -1, +4, +1, -3, -2, -2),
                        "THR": (+0, -1, +0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, +1, +5, -2, -2, +0),
                        "TRP": (-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, +1, -4, -3, -2, 11, +2, -3),
                        "TYR": (-2, -2, -2, -3, -2, -1, -2, -3, +2, -1, -1, -2, -1, +3, -3, -2, -2, +2, +7, -1),
                        "VAL": (+0, -3, -3, -3, -1, -2, -2, -3, -3, +3, +1, -2, +1, -1, -2, -2, +0, -3, -1, +4)}
    # _end_def_

    @staticmethod
    @njit(fastmath=True)
    def sine_cosine(angle: float) -> (float, float):
        """
        Auxiliary function that returns the sin(x) / cos(x) for
        a given input angle 'x'. This method is using Numba for
        speed up. Upon testing the performance we found that we
        average a  ~23x speed up comparing to using only numpy.

        NOTE:  We don't make any checks on the input angle as we
        assume that it has already been checked before this call.
        The whole point is to speed up the repeated calls to the
        trigonometric functions.

        :param angle: the angle (in radians) that we want to get
        the sine / cosine.

        :return: sin(angle), cos(angle)
        """
        # Return the results.
        return np.sin(angle), np.cos(angle)
    # _end_def_

    @staticmethod
    def save_auxiliary(f_id, rec_data, kind=None, output=None):
        """
        This auxiliary (static)  function will save the auxiliary
        bi-products of the input vector construction, such as the
        hydrogen bonds, torsion angles and aromatic rings.

        :param f_id: File id. Usually the PDB-ID is ok to be used
        as a filename to identify the contents.

        :param rec_data: The data we want to save. Usually it is
        a list of things (tuples, dicts, etc.)

        :param kind: This is the type / kind of data that we are
        saving. It can only be of four types: 1) "t_peptides" and
        2) "t_angles". Anything else will force the method to raise
        an exception.

        :param output: This is the main (parent) output directory
        where the data will be saved.

        :return: None.
        """

        # Set the local variables according to the
        # kind of data that we have to save.
        if kind == "t_peptides":
            f_ext, folder_name = "n_pept", "Poly_Peptides"

        elif kind == "t_angles":
            f_ext, folder_name = "angles", "Torsion_Angles"

        else:
            raise ValueError(f" Unknown type of data: {kind}.")
        # _end_if_

        # Check for the output.
        if output is None:
            output = Path.cwd()
        # _end_if_

        # Create an output filename.
        output_path = Path(output / folder_name)

        # This will be true only once (for each kind).
        if not output_path.is_dir():
            output_path.mkdir(parents=True)
        # _end_if_

        # Convert to DataFrame.
        df = DataFrame(data=rec_data)

        # Save to csv file.
        df.to_csv(Path(output_path / f"{f_id}.{f_ext}"), header=False)
    # _end_def_

    # Auxiliary method.
    def _zero_padding(self, x_in, n_peptides, terminal="front"):
        """
        This auxiliary method is used to add the front and back
        terminal vectors, using zero padding. The method can be
        used for any (odd) number of peptides (3-, 5-, 7-, etc.)

        :param x_in: This is the input that will be used to make
        the front or back terminal poly-peptide entries.

        :param n_peptides: This is the number (int) of the poly-
        peptides that we are constructing. It should be strictly
        an odd number.

        :param terminal: This is the type of the terminal. It
        accepts only two (string) values "front" / "back".

        :return: a list with the front or back terminal entries.
        """

        # Quick sanity check.
        if terminal not in {"front", "back"}:
            raise ValueError(f"{self.__class__.__name__}: "
                             f"Zero padding type can be 'front' | 'back': {terminal}.")
        # _end_if_

        # Extract the data from the input.
        x_peptides = x_in["poly-peptides"]
        x_targets = x_in["targets"]
        x_vector = x_in["vector"]

        # Total length of the vector.
        K = x_vector.size

        # Number of elements (per peptide).
        L = K // int(n_peptides)

        # Declare the return list.
        list_out = []

        # Number of required padded elements.
        number_of_pads = int(n_peptides) // 2

        # Iterate to get all the zero padding
        # elements (front or back).
        for i in range(number_of_pads):
            # Create a temporary vector.
            t_vector = np.zeros_like(x_vector)

            # Make a temporary peptides copy.
            t_peptides = x_peptides.copy()

            # Compute the index where we want
            # to start copying elements from.
            kappa = L * (i + 1)

            # Switch according to the front/back side.
            if terminal == "front":
                # Copy the vector elements.
                t_vector[kappa:] = x_vector[:K-kappa].copy()

                # Rotate the elements forwards.
                t_peptides = t_peptides[-(i+1):] + t_peptides[:-(i+1)]

                # Mark the empty entries with None.
                for j in range(i+1):
                    t_peptides[j] = (None, '-', None)
                # _end_for_
            else:
                # Copy the vector elements.
                t_vector[:K-kappa] = x_vector[kappa:].copy()

                # Rotate the elements backwards.
                t_peptides = t_peptides[(i+1):] + t_peptides[:(i+1)]

                # Mark the empty entries with None.
                for j in range(i+1):
                    t_peptides[-(j+1)] = (None, '-', None)
                # _end_for_
            # _end_if_

            # Add it to the final return list.
            list_out.append({"poly-peptides": t_peptides,
                             "targets": x_targets,
                             "vector": t_vector})
        # _end_for_

        # Return the list.
        return list_out

    # _end_def_

    # Main functionality.
    def get_data(self, f_path, n_peptides=3, all_chains=False, save_output_path=None):
        """
        It accepts as input a 'Path(/path/to/the/file)' and returns a list with all
        the information that is extracted. This can be used as input to the ANN for
        predicting the chemical shift values of the target backbone atoms (e.g. "N",
        "C", "CA", "CB", "H", "HA").

        :param f_path: Path of the (protein) PDB file.

        :param n_peptides: (integer) Number of peptides to consider for the input
        vectors. By default, it considers tri-peptides.

        :param all_chains: (bool) if True it will process all chains in the protein.
        By default, it processes only the first (main) chain.

        :param save_output_path: (Path/String) If given it will be the output path
        where all the auxiliary data will be saved.

        :return: A dictionary with (data + sequence) for each processed model from
        the PDB file. The data + sequence are given as:

            1) a dictionary with the following information:

                1.1) poly-peptides that have been generated
                     (index + three-letter code)

                1.2) vector with the input values

                1.3) a list with the target atoms that are
                     available for each vector.

            2) the amino-acid sequence in string (one-letter code).
        """

        # Make sure peptides is int.
        n_peptides = int(n_peptides)

        # Make sure the input peptide is a positive number.
        if n_peptides < 1:
            raise ValueError(f"{self.__class__.__name__}: "
                             f"Input peptide should be greater than zero: {n_peptides}.")
        # _end_if_

        # Make sure the input peptide is an odd number.
        if np.mod(n_peptides, 2) == 0:
            raise ValueError(f"{self.__class__.__name__}: "
                             f"Input peptide should be odd number: {n_peptides}.")
        # _end_if_

        # Ensure the f_path is Path object.
        f_path = Path(f_path)

        # Get the filename.
        f_name = f_path.stem

        # STEP:1 - Create a reduced structure object.
        r_struct = CoarseGrainStructure()

        # STEP:2 - Parse the input file.
        r_struct.parse_martini_pdb_file(f_path, f_type="MRT", all_chains=all_chains)

        # STEP:3 - Compute the internal coordinates.
        r_struct.compute_dihedral_angles_and_more()

        # Get the structure(s).
        structure = r_struct.get_structure

        # Sanity check.
        if not structure:
            raise ValueError(f"{self.__class__.__name__}: "
                             f"PDB file {f_path} is empty.")
        # _end_if_

        # Localize the static method.
        get_sin_cos = self.sine_cosine

        # Output dictionary.
        output_data = {}

        # Will hold the extracted data.
        x_out, amino_acid_seq = {}, {}

        # Process all chains in the structure.
        for chain_id in structure:

            # Make a temporary list with the residues.
            residue_list = list(structure[chain_id].values())

            # Stores all the poly-peptides / torsion angles.
            poly_peptides, torsion_angles = [], []

            # Localize the append method(s).
            poly_peptides_append = poly_peptides.append
            torsion_angles_append = torsion_angles.append

            # Get the length of the chain.
            chain_length = len(residue_list)

            # Amino-acid counter.
            aa_counter = 0

            # Create a new output entry,  for the current chain.
            # This will leave the room to expand to other chains
            # in the future.
            x_out[chain_id] = []

            # Amino-acid sequence.
            seq = []

            # Localize the append method.
            seq_append = seq.append

            # Iterate through the list of amino-acids.
            for r_name in r_struct.get_sequence[chain_id]:

                # Add the one letter residue in the list.
                seq_append(RES_3_TO_1[r_name])
            # _end_for_

            # Get the sequence in one-letter code.
            amino_acid_seq[chain_id] = "".join(seq)

            # Iterate through the list of residues.
            for i, res_i in enumerate(residue_list, start=0):

                # Check the index 'n_peptides' ahead to
                # ensure that we don't go out of bounds.
                if (i + n_peptides - 1) >= chain_length:
                    break
                # _end_if_

                # Increase by one.
                aa_counter += 1

                # Auxiliary poly-peptide level lists.
                vec_i, t_peptide, t_angles = [], [], []

                # Localize the 'append' method(s).
                t_angles_append = t_angles.append
                t_peptide_append = t_peptide.append

                # Localize the 'extend' method.
                vec_i_extend = vec_i.extend

                # Boolean flag: Reset to "True"
                # for every single poly-peptide.
                all_good = True

                # Get the poly-peptide information.
                for k, res_k in enumerate(residue_list[i:i + n_peptides], start=0):

                    # Index (one-based).
                    INDEX = aa_counter + k

                    # Residue name (three-letter code).
                    RES_NAME = res_k["BB"].res_name

                    # Sanity check: This should not happen!
                    if RES_NAME not in ACCEPTED_RES:
                        # Change the flag value.
                        all_good = False

                        # Move on to the next poly-peptide.
                        break
                    # _end_if_

                    # Get the residue-id directly
                    # from the PDB/MRT file.
                    RES_ID = res_k["BB"].res_id

                    # Append the residue info to a separate list.
                    # This builds the poly-peptides information for
                    # later storage.
                    t_peptide_append((INDEX, RES_NAME, RES_ID))

                    # Add the BLOSUM score vector.
                    vec_i_extend(self._blosum[RES_NAME])

                    # Get the 'left', 'right' dihedral angles.
                    # Units: (in radians).
                    dihedral_left_k = res_k["BB"].get_angle("left_dihedral")
                    dihedral_right_k = res_k["BB"].get_angle("right_dihedral")

                    # Get the 'alpha', 'beta' and 'gamma' angles.
                    # Units: (in radians).
                    alpha_k = res_k["BB"].get_angle("alpha")
                    beta_k = res_k["BB"].get_angle("beta")
                    gamma_k = res_k["BB"].get_angle("gamma")

                    # Add all the angles (if they exist).
                    for angle in [dihedral_left_k, dihedral_right_k,
                                  alpha_k, beta_k, gamma_k]:

                        # Check for "None" angles.
                        if angle is None:
                            # Since there is no angle in [0, 2pi]
                            # that satisfies both sine and cosine
                            # to be zero this will imply that the
                            # angle is missing (i.e. is None).
                            vec_i_extend([0.0, 0.0])
                        else:
                            # Add in the vector: [sin(x), cos(x)].
                            vec_i_extend(get_sin_cos(angle))
                        # _end_if_

                    # _end_for_

                    # Append the torsion angles (here in radians).
                    t_angles_append((INDEX, RES_3_TO_1[RES_NAME], dihedral_left_k,
                                     dihedral_right_k, alpha_k, beta_k, gamma_k))
                # _end_for_

                # The total number of entries will vary according
                # to what we will include in the data list vec_i.
                if all_good:

                    # Append all information to a final (return) list.
                    # The entries here are combined into a dictionary.
                    # This way we have easier access using the "keys",
                    # to access the right data for each record.
                    x_out[chain_id].append({"poly-peptides": t_peptide,
                                            "targets": TARGET_ATOMS,
                                            "vector": np.array(vec_i, dtype=float)})

                    # Update the list with the poly-peptides.
                    poly_peptides_append(t_peptide)

                    # Update the list with the torsion angles.
                    torsion_angles_append(t_angles)
                # _end_if_

            # _end_for_

            # Check the list if is empty.
            if len(x_out[chain_id]) > 1:
                # Copy the first entry.
                x_front = x_out[chain_id][0].copy()

                # Get the front end points.
                front_end = self._zero_padding(x_front, n_peptides, terminal="front")

                # Copy the last entry.
                x_back = x_out[chain_id][-1].copy()

                # Get the back end points.
                back_end = self._zero_padding(x_back, n_peptides, terminal="back")

                # Insert the front terminal points (in reverse order).
                # This will place them in the right position in x_out.
                for f, vec_fwd in enumerate(reversed(front_end)):
                    x_out[chain_id].insert(f, vec_fwd)
                # _end_for_

                # Append the back terminal points.
                for vec_bwd in back_end:
                    x_out[chain_id].append(vec_bwd)
                # _end_for_

            # _end_if_

            # Check if we want to save the data.
            if save_output_path:

                # Make sure the save path exists.
                if not Path(save_output_path).is_dir():
                    # This should run only the first time.
                    Path(save_output_path).mkdir(parents=True)
                # _end_if_

                # File ID includes:
                # - 0) the filename "f_name"
                # - 1) the model number = '1'
                # - 2) and chain number = chain_id
                file_id = f"{f_name}_1_{chain_id}"

                # Poly-peptides.
                if poly_peptides:
                    self.save_auxiliary(file_id, poly_peptides,
                                        "t_peptides", save_output_path)
                # _end_if_

                # Torsion angles.
                if torsion_angles:
                    self.save_auxiliary(file_id, torsion_angles,
                                        "t_angles", save_output_path)
                # _end_if_

            # _end_if_

        # _end_for_

        # Add the {"data" + "sequence"} to the output dictionary.
        output_data["model-1"] = {"data": x_out,
                                  "sequence": amino_acid_seq}
        # Output data.
        return output_data
    # _end_def_

    # Auxiliary.
    def __call__(self, *args, **kwargs):
        """
        This is only a "wrapper" method
        of the "get_data" method.
        """
        return self.get_data(*args, **kwargs)
    # _end_def_

# _end_class_
