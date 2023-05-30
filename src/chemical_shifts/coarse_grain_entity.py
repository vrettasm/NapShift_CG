# General import(s).
import numpy as np

# Private import(s).
from src.chemical_shifts.auxiliaries import ACCEPTED_RES


class CoarseGrainEntity(object):
    """
    This class creates a data structure to hold a "coarse grain entity" object.
    """

    # Object variables. These are declared here to minimize the memory footprint
    # of the object. By definition, the class will occupy the specific amount of
    # memory for these objects only.
    __slots__ = ("_num_id", "_atom_type", "_res_name", "_res_id", "_chain_id",
                 "_res_x", "_res_y", "_res_z", "_angle")

    # Constructor.
    def __init__(self, num_id: int, atom_type: str, res_name: str, chain_id: str,
                 res_id: float, res_x: float, res_y: float, res_z: float):
        """
        Default constructor for the "CG" entity structure. This will hold the
        necessary information, such as residue id, atom type, etc.

        :param num_id: (integer) Atom serial number, as read from the input file.

        :param atom_type: (string) Atom type name (BB or SCx).

        :param res_name: (string) Residue name.

        :param chain_id: (string) Chain sequence identifier.

        :param res_id: (integer) Residue sequence number.

        :param res_x: (float) Orthogonal coordinates for X in Angstroms.

        :param res_y: (float) Orthogonal coordinates for Y in Angstroms.

        :param res_z: (float) Orthogonal coordinates for Z in Angstroms.
        """

        # Check for the correct type.
        if isinstance(num_id, int):

            # Assign the number id.
            self._num_id = num_id
        else:
            raise TypeError(f"{self.__class__.__name__}: "
                            f"Number ID should be integer: {type(num_id)}.")
        # _end_if_

        # Check for the correct type.
        if isinstance(atom_type, str):

            # Accept only back-bone (BB), or side-chains (SC) atoms.
            if atom_type == "BB" or atom_type.startswith("SC"):

                # Assign the atom type.
                self._atom_type = atom_type
            else:
                raise ValueError(f"{self.__class__.__name__}: "
                                 f"Atom type should be a 'BB' or 'SC1', 'SC2', etc.")
            # _end_if_
        else:
            raise TypeError(f"{self.__class__.__name__}: "
                            f"Atom type should be string: {type(atom_type)}.")
        # _end_if_

        # Check for the correct type.
        if isinstance(res_name, str):

            # Accept only the standard residues
            # here in are (three-letter codes).
            if res_name in ACCEPTED_RES:

                # Assign the residue type.
                self._res_name = res_name
            else:
                raise ValueError(f"{self.__class__.__name__}: "
                                 f"Residue name is not in the accepted list of residues.")
            # _end_if_
        else:
            raise TypeError(f"{self.__class__.__name__}: "
                            f"Residue name should be string: {type(res_name)}.")
        # _end_if_

        # Check for the correct type.
        if isinstance(res_id, int):

            # Assign the residue id.
            self._res_id = res_id
        else:
            raise TypeError(f"{self.__class__.__name__}: "
                            f"Residue ID should be integer: {type(res_id)}.")
        # _end_if_

        # Check for the correct type.
        if isinstance(chain_id, str):

            # Assign the chain id.
            self._chain_id = chain_id
        else:
            raise TypeError(f"{self.__class__.__name__}: "
                            f"Chain ID type should be string: {type(chain_id)}.")
        # _end_if_

        # Check for the correct type.
        if isinstance(res_x, float):

            # Accept only real valued coordinates.
            if np.isreal(res_x):

                # Assign the x-coordinate.
                self._res_x = res_x
            else:
                raise ValueError(f"{self.__class__.__name__}: "
                                 f"Atom 'X' coordinate should be real: {res_x}")
            # _end_if_
        else:
            raise TypeError(f"{self.__class__.__name__}: "
                            f"Atom 'X' coordinate should be float: {type(res_x)}.")
        # _end_if_

        # Check for the correct type.
        if isinstance(res_y, float):

            # Accept only real valued coordinates.
            if np.isreal(res_y):

                # Assign the y-coordinate.
                self._res_y = res_y
            else:
                raise ValueError(f"{self.__class__.__name__}: "
                                 f"Atom 'Y' coordinate should be real: {res_y}")
            # _end_if_
        else:
            raise TypeError(f"{self.__class__.__name__}: "
                            f"Atom 'Y' coordinate should be float: {type(res_y)}.")
        # _end_if_

        # Check for the correct type.
        if isinstance(res_z, float):

            # Accept only real valued coordinates.
            if np.isreal(res_z):

                # Assign the z-coordinate.
                self._res_z = res_z
            else:
                raise ValueError(f"{self.__class__.__name__}: "
                                 f"Atom 'Z' coordinate should be real: {res_z}")
            # _end_if_
        else:
            raise TypeError(f"{self.__class__.__name__}: "
                            f"Atom 'Z' coordinate should be float: {type(res_z)}.")
        # _end_if_

        # Initialize all angles to None.
        self._angle = {"LEFT_DIHEDRAL": None, "RIGHT_DIHEDRAL": None,
                       "ALPHA": None, "BETA": None, "GAMMA": None}
    # _end_def_

    @property
    def num_id(self):
        """
        Accessor (getter) of the number ID value.

        :return: integer value.
        """
        return self._num_id
    # _end_def_

    @property
    def chain_id(self):
        """
        Accessor (getter) of the chain ID value.

        :return: string value.
        """
        return self._chain_id
    # _end_def_

    @property
    def atom_type(self):
        """
        Accessor (getter) of the atom type.

        :return: string value (BB or SCx).
        """
        return self._atom_type
    # _end_def_

    @property
    def res_name(self):
        """
        Accessor (getter) of the residue name.

        :return: string value (three-letter code).
        """
        return self._res_name
    # _end_def_

    @property
    def res_id(self):
        """
        Accessor (getter) of the residue sequence ID.

        :return: integer value.
        """
        return self._res_id
    # _end_def_

    @property
    def res_x(self):
        """
        Accessor (getter) of the atom 'X' coordinate.

        :return: float value.
        """
        return self._res_x
    # _end_def_

    @property
    def res_y(self):
        """
        Accessor (getter) of the atom 'Y' coordinate.

        :return: float value.
        """
        return self._res_y
    # _end_def_

    @property
    def res_z(self):
        """
        Accessor (getter) of the atom 'Z' coordinate.

        :return: float value.
        """
        return self._res_z
    # _end_def_

    @property
    def get_coord(self):
        """
        Accessor (getter) of the atom coordinates (x, y, z),
        all together in a single tuple.

        :return: float values.
        """
        return self._res_x, self._res_y, self._res_z
    # _end_def_

    # Auxiliary method.
    def get_angle(self, key_value: str):
        """
        Get the angle (specified by key_value).

        :param key_value: (string) key value for angle.

        :return: the angle value in radians.
        """

        # Make the key word uppercase.
        key_value = key_value.upper()

        # Check if the angle exists.
        if key_value in self._angle:

            return self._angle[key_value]
        else:
            raise KeyError(f"{self.__class__.__name__}: "
                           f"Angle {key_value} does not exist.")
        # _end_if_

    # _end_def_

    # Auxiliary method.
    def set_angle(self, key_angle: str, value: float):
        """
        Set the angle (specified by key_angle) with a
        proper value.

        :param key_angle: (string) key value for angle.

        :param value: (float) for the angle in radians.

        """

        # Make sure the value is float.
        value = float(value)

        # Make the key word uppercase.
        key_angle = key_angle.upper()

        # Check if the angle exists.
        if key_angle in self._angle:

            # Allow different range to the dihedral angles.
            if key_angle in ["LEFT_DIHEDRAL", "RIGHT_DIHEDRAL"]:

                # Check for the correct range.
                if -np.pi <= value <= +np.pi:

                    # Assign the dihedral angle.
                    self._angle[key_angle] = value
                else:
                    raise ValueError(f"{self.__class__.__name__}: "
                                     f"Angle {key_angle} should be in [-pi, +pi]: {value}.")
                # _end_if_

            else:

                # Check for the correct range.
                if 0.0 <= value <= +np.pi:

                    # Assign the angle value.
                    self._angle[key_angle] = value
                else:
                    raise ValueError(f"{self.__class__.__name__}: "
                                     f"Angle {key_angle} should be in [0, +pi]: {value}.")
                # _end_if_

            # _end_if_
        else:
            raise KeyError(f"{self.__class__.__name__}: "
                           f"Angle {key_angle} does not exist.")
        # _end_if_

    # _end_def_

    # Auxiliary (used for sorting).
    def __lt__(self, other):
        """
        For this class we use the 'atom type' to define what is
        'less than' operator. 'BB' < 'SC1' < 'SC2' < ... < etc.

        :param other: object to compare with (of the same class).

        :return: true if the 'self' object is less than the 'other'.
        """
        # Check for the same type.
        if isinstance(other, type(self)):

            # Both objects should belong in the same chain.
            if self.chain_id == other.chain_id:

                return self.atom_type < other.atom_type
            else:
                raise RuntimeError(f"{self.__class__.__name__}: "
                                   f"Both objects should belong in the same chain: "
                                   f"{self.chain_id} != {other.chain_id}.")
            # _end_if_
        else:
            return NotImplemented
        # _end_if_
    # _end_def_

    # Auxiliary (used for sorting).
    def __gt__(self, other):
        """
        For this class we use the 'atom type' to define what is
        'greater than' operator. 'SCx' > ... > 'SC2' > 'SC1' > 'BB'.

        :param other: object to compare with (of the same class).

        :return: true if the 'self' object is greater than the 'other'.
        """
        # Check for the same type.
        if isinstance(other, type(self)):

            # Both objects should belong in the same chain.
            if self.chain_id == other.chain_id:

                return self.atom_type > other.atom_type
            else:
                raise RuntimeError(f"{self.__class__.__name__}: "
                                   f"Both objects should belong in the same chain: "
                                   f"{self.chain_id} != {other.chain_id}.")
            # _end_if_
        else:
            return NotImplemented
    # _end_def_

    # Auxiliary (used for sorting).
    def __le__(self, other):
        """
        Same as '__lt__' only with the equality included.

        :param other: object to compare with (of the same class).

        :return: true if the 'self' object is less than or equal
        the 'other'.
        """
        # Check for the same type.
        if isinstance(other, type(self)):

            # Both objects should belong in the same chain.
            if self.chain_id == other.chain_id:

                return self.atom_type <= other.atom_type
            else:
                raise RuntimeError(f"{self.__class__.__name__}: "
                                   f"Both objects should belong in the same chain: "
                                   f"{self.chain_id} != {other.chain_id}.")
            # _end_if_
        else:
            return NotImplemented
    # _end_def_

    # Auxiliary (used for sorting).
    def __ge__(self, other):
        """
        Same as '__gt__' only with the equality included.

        :param other: object to compare with (of the same class).

        :return: true if the 'self' object is greater than or equal
        the 'other'.
        """
        # Check for the same type.
        if isinstance(other, type(self)):

            # Both objects should belong in the same chain.
            if self.chain_id == other.chain_id:

                return self.atom_type >= other.atom_type
            else:
                raise RuntimeError(f"{self.__class__.__name__}: "
                                   f"Both objects should belong in the same chain: "
                                   f"{self.chain_id} != {other.chain_id}.")
            # _end_if_
        else:
            return NotImplemented
    # _end_def_

    # Auxiliary (used for sorting).
    def __eq__(self, other):
        """
        Check if two objects are assumed equal.

        :param other: object to compare with (of the same class).

        :return: true if the 'self' object is equal than the 'other'.
        """
        # Check for the same type.
        if isinstance(other, type(self)):

            # Both objects should belong in the same chain.
            if self.chain_id == other.chain_id:

                return self.atom_type == other.atom_type
            else:
                raise RuntimeError(f"{self.__class__.__name__}: "
                                   f"Both objects should belong in the same chain: "
                                   f"{self.chain_id} != {other.chain_id}.")
            # _end_if_
        else:
            return NotImplemented
    # _end_def_

    # Auxiliary.
    def __str__(self):
        """
        Override to print a readable string presentation of the object.
        This will include its id(), along with its values of its flags.

        :return: a string representation of a CoarseGrainEntity object.
        """

        # Return the f-string.
        return f" CoarseGrainEntity Id({id(self)}): \n" \
               f" ID={self.num_id} \n" \
               f" Name={self.res_name} \n" \
               f" Type={self.atom_type} \n" \
               f" Res-Id={self.res_id} \n" \
               f" Chain-Id={self.chain_id} \n" \
               f" Coordinates={self.get_coord}"
    # _end_def_

# _end_class_
