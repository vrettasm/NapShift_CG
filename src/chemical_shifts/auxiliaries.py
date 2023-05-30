"""
Constants that are used by all code are placed here for easy access.
These include:

    1. ACCEPTED_RES
    2. RES_3_TO_1
    3. TARGET_ATOMS
    4. RANDOM_COIL_TBL
    5. RANDOM_COIL_AVG
    6. RANDOM_COIL_STD

"""
from collections import namedtuple

# Add documentation to the NamedTuple.
__pdoc__ = {}

# Residue (amino-acid) names (three-letters code).
# This set is updated with: CYO /CYR and PRT /PRC.
ACCEPTED_RES = {"ALA", "ARG", "ASN", "ASP",
                "CYS", "GLN", "GLU", "GLY",
                "HIS", "ILE", "LEU", "LYS",
                "MET", "PHE", "PRO", "SER",
                "THR", "TRP", "TYR", "VAL",
                "CYO", "CYR", "PRT", "PRC"}

# Protein (amino-acid) 3-to-1 letters mapping.
# NOTE in the modified map:
# The CYO --> X, CYR --> C.
# The PRC --> O, PRT --> P.
RES_3_TO_1 = {"ALA": 'A', "CYS": 'C', "CYO": 'X', "CYR": 'C',
              "ASP": 'D', "GLU": 'E', "PHE": 'F', "GLY": 'G',
              "HIS": 'H', "ILE": 'I', "LYS": 'K', "LEU": 'L',
              "MET": 'M', "ASN": 'N', "PRO": 'P', "PRC": 'O',
              "PRT": 'P', "GLN": 'Q', "ARG": 'R', "SER": 'S',
              "THR": 'T', "VAL": 'V', "TRP": 'W', "TYR": 'Y'}

# These are the target atoms that are predicted.
# NOTE: The order of the entries matters, so we
#       can't change this to a set()!!
TARGET_ATOMS = ('N', 'C', 'CA', 'CB', 'H', 'HA')

# Module level declaration.
ChemShifts = namedtuple("ChemShifts", TARGET_ATOMS)

# Add documentation for the fields.
__pdoc__["ChemShifts.N"] = "Random coil chemical shift value for 'N'."
__pdoc__["ChemShifts.C"] = "Random coil chemical shift value for 'C'."
__pdoc__["ChemShifts.H"] = "Random coil chemical shift value for 'H'."

__pdoc__["ChemShifts.CA"] = "Random coil chemical shift value for 'CA'."
__pdoc__["ChemShifts.CB"] = "Random coil chemical shift value for 'CB'."
__pdoc__["ChemShifts.HA"] = "Random coil chemical shift value for 'HA'."

# Random-coil (corrections) chemical shift values.
# Source: https://pubs.acs.org/doi/10.1021/ja105656t
# -----------------------------------------------------------------------------------
RANDOM_COIL_TBL = {"ALA": ChemShifts(123.960, 178.418, 52.599, 19.102, 8.158, 4.224),
                   "CYO": ChemShifts(120.300, 173.500, 55.300, 40.500, 8.500, 4.850),
                   "CYR": ChemShifts(119.068, 174.927, 58.327, 28.085, 8.410, 4.447),
                   "ASP": ChemShifts(120.207, 176.987, 54.331, 41.089, 8.217, 4.537),
                   "GLU": ChemShifts(120.769, 177.125, 56.650, 30.225, 8.304, 4.222),
                   "PHE": ChemShifts(120.138, 176.368, 57.934, 39.660, 8.107, 4.573),
                   "GLY": ChemShifts(108.783, 174.630, 45.236, 00.000, 8.370, 3.980),
                   "HIS": ChemShifts(118.930, 175.349, 55.964, 29.719, 8.310, 4.585),
                   "ILE": ChemShifts(120.512, 176.897, 61.247, 38.563, 7.963, 4.076),
                   "LYS": ChemShifts(121.353, 177.224, 56.412, 32.921, 8.221, 4.237),
                   "LEU": ChemShifts(121.877, 178.037, 55.260, 42.212, 8.088, 4.260),
                   "MET": ChemShifts(120.002, 176.953, 55.591, 32.690, 8.209, 4.425),
                   "ASN": ChemShifts(118.668, 175.825, 53.231, 38.790, 8.366, 4.632),
                   "PRC": ChemShifts(136.612, 177.542, 63.180, 32.072, 0.000, 4.339),
                   "PRT": ChemShifts(136.612, 177.542, 63.180, 32.072, 0.000, 4.339),
                   "GLN": ChemShifts(120.224, 176.510, 55.840, 29.509, 8.258, 4.254),
                   "ARG": ChemShifts(121.288, 176.821, 56.088, 30.691, 8.232, 4.239),
                   "SER": ChemShifts(115.935, 175.236, 58.352, 63.766, 8.215, 4.392),
                   "THR": ChemShifts(114.024, 175.122, 61.926, 69.794, 8.047, 4.252),
                   "VAL": ChemShifts(120.403, 176.772, 62.347, 32.674, 8.037, 4.009),
                   "TRP": ChemShifts(120.733, 174.549, 57.500, 29.380, 7.725, 4.567),
                   "TYR": ChemShifts(120.228, 176.284, 57.761, 38.750, 8.026, 4.504)}

# Random-coil (avg) chemical shift values.
# Source: https://bmrb.io/published/Ikura_cs_study/part2_rc_aa_cs_stats.pdf
# -----------------------------------------------------------------------------------
RANDOM_COIL_AVG = {"ALA": ChemShifts(124.200, 176.800, 52.100, 19.300, 8.240, 4.390),
                   "CYO": ChemShifts(120.300, 173.500, 55.300, 40.500, 8.500, 4.850),
                   "CYR": ChemShifts(120.100, 174.700, 58.200, 29.400, 8.200, 4.960),
                   "ASP": ChemShifts(121.500, 175.700, 53.800, 41.200, 8.310, 4.710),
                   "GLU": ChemShifts(121.400, 175.900, 56.300, 30.300, 8.360, 4.390),
                   "PHE": ChemShifts(120.100, 175.000, 57.200, 40.200, 8.270, 4.650),
                   "GLY": ChemShifts(109.800, 173.900, 45.200, 00.000, 8.310, 4.120),
                   "HIS": ChemShifts(119.700, 174.400, 55.300, 30.100, 8.290, 4.730),
                   "ILE": ChemShifts(122.000, 174.900, 60.400, 38.700, 8.300, 4.310),
                   "LYS": ChemShifts(121.000, 176.000, 56.200, 32.800, 8.240, 4.360),
                   "LEU": ChemShifts(122.300, 176.400, 54.500, 42.500, 8.180, 4.470),
                   "MET": ChemShifts(121.200, 174.600, 55.400, 33.700, 8.350, 4.410),
                   "ASN": ChemShifts(119.400, 174.900, 53.000, 38.900, 8.420, 4.750),
                   "PRC": ChemShifts(135.355, 176.100, 62.600, 31.900, 8.756, 4.440),
                   "PRT": ChemShifts(135.355, 176.100, 62.600, 31.900, 8.756, 4.440),
                   "GLN": ChemShifts(120.200, 175.700, 55.500, 29.400, 8.210, 4.430),
                   "ARG": ChemShifts(121.700, 175.300, 55.900, 31.000, 8.240, 4.470),
                   "SER": ChemShifts(116.800, 174.200, 58.100, 64.100, 8.360, 4.550),
                   "THR": ChemShifts(114.600, 174.500, 60.900, 69.700, 8.270, 4.550),
                   "VAL": ChemShifts(121.800, 175.100, 61.400, 32.800, 8.320, 4.300),
                   "TRP": ChemShifts(121.700, 175.500, 57.300, 30.400, 8.190, 4.800),
                   "TYR": ChemShifts(120.000, 174.800, 57.600, 39.400, 8.240, 4.730)}

# Random-coil (std) chemical shift values.
# Source: https://bmrb.io/published/Ikura_cs_study/part2_rc_aa_cs_stats.pdf
# -----------------------------------------------------------------------------------
RANDOM_COIL_STD = {"ALA": ChemShifts(6.900, 2.100, 1.900, 2.000, 0.630, 0.430),
                   "CYO": ChemShifts(6.200, 1.700, 2.600, 2.100, 0.830, 0.840),
                   "CYR": ChemShifts(4.900, 1.700, 2.200, 4.000, 0.700, 0.340),
                   "ASP": ChemShifts(4.400, 1.600, 2.000, 1.600, 0.590, 0.360),
                   "GLU": ChemShifts(4.200, 1.800, 1.900, 1.800, 0.650, 0.440),
                   "PHE": ChemShifts(4.600, 2.100, 2.100, 1.700, 0.840, 0.470),
                   "GLY": ChemShifts(4.200, 1.600, 1.500, 1.000, 0.880, 2.290),
                   "HIS": ChemShifts(5.400, 1.900, 2.500, 2.300, 1.001, 0.810),
                   "ILE": ChemShifts(5.300, 1.800, 2.200, 1.700, 0.780, 0.450),
                   "LYS": ChemShifts(4.300, 1.700, 1.800, 1.700, 0.600, 0.420),
                   "LEU": ChemShifts(4.400, 1.900, 2.000, 2.400, 0.750, 0.430),
                   "MET": ChemShifts(4.600, 2.100, 1.700, 2.200, 0.620, 0.480),
                   "ASN": ChemShifts(4.400, 1.600, 1.900, 2.300, 0.720, 0.360),
                   "PRC": ChemShifts(5.450, 1.700, 1.300, 1.100, 0.710, 0.410),
                   "PRT": ChemShifts(5.450, 1.700, 1.300, 1.100, 0.710, 0.410),
                   "GLN": ChemShifts(4.400, 1.900, 2.000, 2.200, 0.730, 0.460),
                   "ARG": ChemShifts(4.500, 1.900, 2.000, 1.900, 0.700, 0.410),
                   "SER": ChemShifts(4.300, 1.700, 1.900, 1.600, 0.720, 0.410),
                   "THR": ChemShifts(5.000, 1.600, 2.200, 4.500, 0.660, 0.440),
                   "VAL": ChemShifts(5.100, 1.600, 2.200, 2.000, 0.700, 0.470),
                   "TRP": ChemShifts(4.300, 1.700, 2.300, 1.600, 0.740, 0.500),
                   "TYR": ChemShifts(4.500, 1.800, 2.400, 3.200, 0.750, 0.520)}


# _end_module_
