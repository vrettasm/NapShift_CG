# Python import(s).
import h5py
import joblib
import numpy as np
from pathlib import Path
from pandas import DataFrame, read_csv

# Other custom imports.
from src.chemical_shifts.coarse_grain_input_vector import CoarseGrainInputVector
from src.chemical_shifts.auxiliaries import (RES_3_TO_1, TARGET_ATOMS,
                                             RANDOM_COIL_TBL)
# CamCoil random coil prediction engine.
from src.random_coil.camcoil import CamCoil


# Predictor Class.
class CoarseGrainCSPredictor(object):

    # Default (input) directory.
    dir_default = Path.cwd()
    """
    This will be used as the default directory.
    """

    # Object to create the input vectors.
    r_vec_in = CoarseGrainInputVector()
    """
    This will be used to load the MARTINI PDB file(s)
    before we make the prediction with the ANN model.
    """

    # Convert the random coil (average) values to DataFrame.
    df_avg = DataFrame(RANDOM_COIL_TBL, index=TARGET_ATOMS)
    """
    This dataframe will be used in case we do not provide a
    random coil file.
    """

    # Camcoil predictor of random coil chemical shifts.
    random_coil = CamCoil(pH=7.0)
    """
    This object will be used to predict the random coil
    chemical shifts. The default value for the "pH" is
    set to '7.0'.
    """

    # Object variables.
    __slots__ = ("dir_model", "dir_output", "overwrite", "nn_weights_cg",
                 "input_scaler_cg")

    # Constructor.
    def __init__(self, dir_model=None, dir_output=None, overwrite=True):
        """
        Constructs an object that will perform the chemical shift prediction.
        This is done by first constructing the necessary input values from the
        PDB (martini) file and subsequently using the pre-trained nn-models to
        perform the predictions on all backbone atoms.

        :param dir_model: (Path) Directory where the trained ANN weights exist
        (one for each target atom).

        :param dir_output: (Path) Directory where the output of the predictions
        will be saved.

        :param overwrite: (bool) Overwrite file protection. If "True", then the
        output file WILL overwrite any pre-existing output.
        """

        # Input data scaler. Set initial value to "None".
        self.input_scaler_cg = {atom: None for atom in TARGET_ATOMS}

        # Neural network weights. Set initial value to "None".
        self.nn_weights_cg = {atom: {"Hidden": None,
                                     "Output": None} for atom in TARGET_ATOMS}

        # Check if we have given explicitly a new location
        # for the input files.
        if dir_model is None:

            # This is the default input location.
            self.dir_model = self.dir_default
        else:

            # This will be the new location.
            self.dir_model = Path(dir_model)
        # _end_if_

        # Make sure the input directory ALWAYS exists.
        if not self.dir_model.is_dir():
            raise ValueError(f"{self.__class__.__name__}: "
                             f"Input directory doesn't exist: {self.dir_model}.")
        # _end_if_

        # Check if we have given explicitly a new location
        # for the output results.
        if dir_output is None:

            # This is the default output location.
            self.dir_output = self.dir_default
        else:
            # This will be the new location.
            self.dir_output = Path(dir_output)

            # If the output directory doesn't exist
            # create it, along with all its parents.
            if not self.dir_output.is_dir():
                self.dir_output.mkdir(parents=True)
            # _end_if_

        # _end_if_

        # Boolean flag. If "True" the class will allow
        # the new results to overwrite old ones (if exist).
        if isinstance(overwrite, bool):

            # Copy overwrite flag value.
            self.overwrite = overwrite
        else:
            raise TypeError(f"{self.__class__.__name__}: "
                            f"Overwrite protection flag should be bool: {type(overwrite)}.")
        # _end_if_

        # Load all the trained models (one for each atom).
        for atom in TARGET_ATOMS:

            # Make the file path.
            file_path = Path(self.dir_model / f"CG_ann_model_{atom}.h5")

            # Check if the file exists.
            if not file_path.is_file():
                # Show a message instead of raising an exception.
                print(f" {self.__class__.__name__}. WARNING:"
                      f" Model for target '{atom}', doesn't exist. Skipping ...")

                # Skip to the next atom.
                continue
            # _end_if_

            # If everything is OK load the model.
            with h5py.File(file_path, "r") as file_h5:

                # Extract the weights.
                _weights = file_h5["model_weights"]

                # Extract kernel/bias for the all layers.
                for _layer in ["Hidden", "Output"]:
                    self.nn_weights_cg[atom][_layer] = {"bias": np.array(
                        _weights[f"{_layer}_1/{_layer}_1/bias:0"]),
                                                        "kernel": np.array(
                                                            _weights[f"{_layer}_1/{_layer}_1/kernel:0"])}
                # _end_for_

            # _end_with_

            # Create the Scaler filename.
            scaler_path = Path(self.dir_model / f"CG_data_scaler_{atom}.gz")

            # Check if the scaler exists.
            if scaler_path.is_file():
                # Load the Scaler from the file.
                self.input_scaler_cg[atom] = joblib.load(scaler_path)
            else:
                # Display what went wrong.
                print(f" {self.__class__.__name__}."
                      f" WARNING: File {scaler_path} not found.")
            # _end_if_

        # _end_for_

    # _end_def_

    # Auxiliary method.
    def save_to_file(self, pdb_id, aa_sequence, predictions, target_peptides, ref_peptides,
                     random_coil=None, model_id=None, chain_id=None, talos_format=True):
        """
        This method accepts the output of the __call__() method and writes all the
        information in a text file. Optionally we can save using the TALOS format.

        :param pdb_id: (string) This is the PDB-ID from the input file. It is used
        to provide information in the final output file.

        :param aa_sequence: (string) This is the sequence of the amino acids in the
        input file. It is used for information to the final output file.

        :param predictions: These are the predictions (numerical values) of the ANN.
        It is a dictionary where each entry (key) corresponds to an atom-target.

        :param target_peptides:  These are the poly-peptides that were predicted by
        the artificial neural network. Because of the "aromatic-rings effect" there
        could be poly-peptides that were not predicted for all the targets.

        :param ref_peptides: These are ALL the poly-peptides, as constructed by
        the InputVector class. They are used as reference regarding the list of
        "poly_peptides".

        :param random_coil: This is a DataFrame with the random coil values. If it
        isn't given (default=None) we will use average values from a default table.

        :param model_id: This is the id of the model in the PDB file. We use it to
        distinguish the output results.

        :param chain_id: This is the id of the chain in the protein model.
        We use it to distinguish the output results.

        :param talos_format: (bool) Flag that defines the file format. If is set to
        "True" we will use the TALOS format. If it is set to "False" the file will
        be saved with a default tabular format.

        :return: None.
        """

        try:
            # Construct the model-id for the file name.
            model_id = "1" if model_id is None else model_id

            # Construct the chain-id for the file name.
            chain_id = "A" if chain_id is None else chain_id

            # Construct the output file name.
            f_name_out = Path(self.dir_output / f"CG_prediction_{pdb_id}_"
                                                f"model_{model_id}_chain_{chain_id}.tab")

            # Check if we have enabled overwrite protection.
            if (not self.overwrite) and f_name_out.is_file():
                raise FileExistsError(f" Output: {f_name_out} already exists.")
            # _end_if_

            # Check there is a random coil dataframe.
            if random_coil is not None:
                # This will optimize the searches.
                random_coil.set_index(["ID", "RES"], inplace=True)
            # _end_if_

            # Size of chunks.
            n = 20

            # Split the amino-acid sequence to chunk of size 'n'.
            chunks = [aa_sequence[i:i + n] for i in range(0, len(aa_sequence), n)]

            # Write the prediction data to a text file.
            with open(f_name_out, "w") as f_out:

                # Localize the write function.
                file_write = f_out.write

                # In case we need to add comments. This isn't mandatory, but
                # it will let us know what the original file was coming from.
                file_write(f"REMARK Chemical Shift predictions for {pdb_id}. \n")

                # Model/Chain information
                file_write(f"REMARK Model {model_id} / Chain {chain_id}. \n")

                # Empty line.
                file_write("\n")

                # Default value is set to N/A (optional).
                file_write("DATA FIRST_RESID N/A \n")

                # Empty line.
                file_write("\n")

                # Write the whole sequence in chunks of size "n".
                for sub_k in chunks:
                    file_write("DATA SEQUENCE " + sub_k + "\n")
                # _end_for_

                # Empty line.
                file_write("\n")

                # Check the file format.
                if talos_format:
                    # Write the (TALOS) variable names.
                    file_write("VARS RESID RESNAME ATOMNAME SHIFT \n")

                    # Write the (TALOS) file format.
                    file_write("FORMAT %4d %1s %4s %8.3f \n")
                else:
                    # Declare a dictionary to group the data
                    # values according to their atom values.
                    record = {atom: np.nan for atom in TARGET_ATOMS}

                    # Tabular text format. This will preserve the same order
                    # (of atoms) as in the TARGET_ATOMS (tuple) declaration.
                    file_write("{:>4} {:>4} {:>8} {:>8} {:>8} {:>8} {:>8} "
                               "{:>8} \n".format("ID", "RES", *record.keys()))
                # _end_if_

                # Empty line.
                file_write("\n")

                # Extract the data and write them to the file.
                for n, peptide in enumerate(ref_peptides, start=1):

                    # Extract the information.
                    index, res_name, res_id = peptide

                    # Convert the name from 3 to 1 letters.
                    res_name_1 = RES_3_TO_1[res_name]

                    # Search link.
                    search_link = (index, res_name_1)

                    # Extract the predicted chemical shifts.
                    for atom in TARGET_ATOMS:

                        # Setting to NaN will indicate that we don't
                        # have a predicted "ss" value for this atom.
                        ss_value, rc_value = np.nan, np.nan

                        # Create a search-peptide tuple.
                        search_peptide = tuple(peptide)

                        # If the peptide is in the target list.
                        if search_peptide in target_peptides[atom]:

                            # Get its index.
                            idx = target_peptides[atom].index(search_peptide)

                            # Get the predicted (secondary structure)
                            # value that comes directly from the ANN.
                            ss_value = predictions[atom][idx].item()

                            # Get the random coil chemical shift.
                            if random_coil is not None:
                                # Get the value from the random coil file.
                                rc_value = random_coil.loc[search_link, atom]
                            else:
                                # Get the average value from a Table.
                                rc_value = CoarseGrainCSPredictor.df_avg.loc[atom, res_name]
                            # _end_if_

                        # _end_if_

                        # Add the random coil value to
                        # the secondary structure value.
                        value = ss_value + rc_value

                        # Check the file format.
                        if talos_format:
                            # Rename the hydrogen from "H" to "HN".
                            atom = "HN" if atom == "H" else atom

                            # Put all the information together in one record.
                            file_write(f"{res_id:>4} {res_name_1:>3} {atom:>6} {value:>8.3f} \n")
                        else:
                            # Store it to the dictionary.
                            record[atom] = value
                        # _end_if_

                    # _end_for_

                    # Check the file format.
                    if not talos_format:
                        # NOTE: From Python "3.6" onwards, the standard dict type maintains
                        # insertion order by default!
                        file_write("{:>4} {:>4} {:>8.3f} {:>8.3f} {:>8.3f} {:>8.3f} {:>8.3f} "
                                   "{:>8.3f} \n".format(res_id, res_name_1, *record.values()))
                    # _end_if_

                # _end_for_

            # _end_with_
        except FileExistsError as e0:

            # Show the error message.
            print(f" {self.__class__.__name__}. ERROR: {e0}")

        # _end_try_

    # _end_def_

    @staticmethod
    def f_hidden_elu(x, alpha=1.0):
        """
        Exponential Linear Unit.

        Auxiliary function that returns the elu(x)
        for a given input vector 'x'.

        :param x: the input values for the NN.

        :param alpha: Scale parameter that controls the
        value to which an ELU saturates for negative net
        inputs. Default is 1.0.

        :return: Elu(x).
        """
        # Copy all the input values to the output.
        f_out = np.array(x, copy=True)

        # Find the negative values (indexes).
        neg_ = (x < 0.0)

        # Update only the negative values.
        f_out[neg_] = alpha * np.expm1(x[neg_])

        # NOTE (from the documentation):
        # Function "expm1(x)" provides greater precision
        # than "exp(x) - 1" for small values of input x.

        # Return the ELU.
        return f_out

    # _end_def_

    # Neural Network output function.
    def f_output_linear(self, x=None, atom=None):
        """
        This method implements the output (activation)
        function of the neural network. The output is
        linear (for regression problems).

        :param x: This is the input vector to the NN.

        :param atom: This is the atom to load the right
        weights (from the corresponding network).

        :return: The output of the NN. This is the "secondary CS"
        value defined as (in LaTex): $\delta{X} - \delta{X}_{rc}$
        """

        # Get the numpy weights for the "hidden" layer.
        w_h = self.nn_weights_cg[atom]["Hidden"]["kernel"]
        b_h = self.nn_weights_cg[atom]["Hidden"]["bias"]

        # Get the numpy weights for the "output" layer.
        w_o = self.nn_weights_cg[atom]["Output"]["kernel"]
        b_o = self.nn_weights_cg[atom]["Output"]["bias"]

        # Return the output function.
        return self.f_hidden_elu(x.dot(w_h) + b_h).dot(w_o) + b_o
    # _end_def_

    # Main method.
    def predict(self, f_path, n_peptides=3, random_coil_path=None, verbose=False, talos_fmt=True):
        """
        It accepts a PDB file as input, constructs the input to the trained NN
        and puts the results (predicted chemical shifts) in a new text file.

        :param f_path: (string) PDB file with the residue / atom coordinates.

        :param n_peptides: (int) Number of peptides to consider for the input
        vectors. By default, it considers tri-peptides.

        :param random_coil_path: (string) file with the random coil chemical
        shift values.

        :param verbose: (bool) If 'True' it will display more info during the
        prediction. The default setting is 'False' to avoid cluttering the
        screen with information.

        :param talos_fmt: (bool) If "True" (default) it will use the TALOS format
        to save the results. If it is set to 'False' the output format will be tabular.

        :return: It will call the save method to write the results in a TALOS
        file format.
        """

        # Make sure the input file is a Path.
        f_path = Path(f_path)

        # Sanity check.
        if not f_path.is_file():
            raise FileNotFoundError(f"{self.__class__.__name__} : "
                                    f"File {f_path} doesn't exist.")
        # _end_if_

        # Call the 'vector' method to create the input values.
        model_vectors = CoarseGrainCSPredictor.r_vec_in(f_path,
                                                        n_peptides=n_peptides)
        # Index of the middle element.
        mid_idx = int(n_peptides) // 2

        # Early exit if the input data is empty.
        # This shouldn't happen very frequently.
        if not model_vectors:
            # Display a warning message.
            print(f" {self.__class__.__name__}. WARNING:"
                  f" File: {f_path} is empty.")

            # Exit from here.
            return None
        # _end_if_

        # Random coil shift values
        # (from an input file).
        random_coil_shifts = None

        # Check if a random coil file is given.
        if random_coil_path:
            # Make sure the input is a Path.
            rc_path = Path(random_coil_path)

            # Sanity check.
            if rc_path.is_file():
                # Extract the random coil chem-shifts.
                # The first row (id->0) is the header.
                random_coil_shifts = read_csv(rc_path, header=0)
            else:
                # Display a message.
                print(f" {self.__class__.__name__}. WARNING:"
                      f" Random coil file {rc_path} doesn't exist.")
            # _end_if_

        # _end_if_

        # Predict chemical shifts of all models.
        for model_id, model_value in enumerate(model_vectors.values(), start=0):

            # Unpack the contents.
            input_data = model_value["data"]
            amino_acid_seq = model_value["sequence"]

            # Sanity check.
            if len(input_data) != len(amino_acid_seq):
                raise RuntimeError(f" {self.__class__.__name__} :"
                                   f" Data / Sequence length mismatch.")
            # _end_if_

            # Process all chains in the model.
            for chain_id, chain_seq in amino_acid_seq.items():

                # Sanity check.
                if not input_data[chain_id]:

                    # Display a warning message.
                    if verbose:
                        # Display a warning message.
                        print(f" {self.__class__.__name__}. WARNING:"
                              f" Model: {model_id} -"
                              f" Chain: {chain_id} didn't produce any input data.")
                    # _end_if_

                    # Go to the next model.
                    continue
                # _end_if_

                # Check if a random coil file is given.
                if random_coil_shifts is not None:
                    # Assign the file shift values.
                    df_random_coil = random_coil_shifts
                else:
                    # Use Camcoil algorithm to predict the values.
                    df_random_coil = CoarseGrainCSPredictor.random_coil(chain_seq)
                # _end_if_

                # Declare data dictionary.
                data = {atom: [] for atom in TARGET_ATOMS}

                # Declare peptides dictionary.
                y_peptide = {atom: [] for atom in TARGET_ATOMS}

                # Reference poly-peptide list.
                ref_peptide = []

                # Localize append method.
                ref_peptide_append = ref_peptide.append

                # Separate all the input data.
                for entry in input_data[chain_id]:

                    # Extract only the middle one.
                    peptide = entry["poly-peptides"][mid_idx]

                    # Check all target atoms.
                    for atom in TARGET_ATOMS:

                        # Check for membership.
                        if atom in entry["targets"]:
                            # Add the data (numpy) vector.
                            data[atom].append(entry["vector"])

                            # Add the peptide information.
                            y_peptide[atom].append(peptide)
                        # _end_if_

                    # _end_for_

                    # Here we add all the peptides, because
                    # these will be the "reference" peptides.
                    ref_peptide_append(peptide)
                # _end_for_

                # Model predictions: Initialize them with "None".
                y_predict = {atom: None for atom in TARGET_ATOMS}

                # Run through all atom-models.
                for atom in TARGET_ATOMS:

                    # Convert the data to a DataFrame.
                    df_data = DataFrame(data[atom])

                    # Check if we have to scale the data.
                    if self.input_scaler_cg[atom]:
                        df_data = DataFrame(self.input_scaler_cg[atom].transform(df_data))
                    # _end_if_

                    # Get the model predictions (secondary structure values).
                    y_predict[atom] = self.f_output_linear(x=df_data.to_numpy(),
                                                           atom=atom)
                # _end_for_

                # Send the predictions to the "save method".
                self.save_to_file(f_path.stem, chain_seq, y_predict, y_peptide, ref_peptide,
                                  df_random_coil, model_id, chain_id, talos_format=talos_fmt)
            # _end_for_
        # _end_for_

    # _end_def_

    # Auxiliary.
    def __call__(self, *args, **kwargs):
        """
        This is only a "wrapper" method
        of the "predict" method.
        """
        return self.predict(*args, **kwargs)

    # _end_def_

    # Auxiliary.
    def __str__(self):
        """
        Override to print a readable string presentation of the object.
        This will include its id(), along with its field values.

            NOTE: Overwrite protection is the opposite of to overwrite flag,
            so in the printed version we show the "not overwrite_flag"!

        :return: a string representation of a CoarseGrainCSPredictor object.
        """

        # Local import of new line.
        from os import linesep as new_line

        # Return the f-string.
        return f" CoarseGrainCSPredictor Id({id(self)}): {new_line}" \
               f" Models dir={self.dir_model} {new_line}" \
               f" Output dir={self.dir_output} {new_line}" \
               f" Overwrite protection={(not self.overwrite)}"
    # _end_def_

# _end_class_
