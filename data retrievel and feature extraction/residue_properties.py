VOLUME_BY_AA = {"G": 60.1, "A": 88.6, "V": 140, "L": 166.7, "I": 166.7, "M": 162.9, "F": 189.9, "Y": 193.6, "W": 227.8,
                "S": 89, "T": 116.1, "C": 108.5, "P": 112.7, "N": 114.1, "Q": 143.8, "D": 111.1, "E": 138.4,
                "K": 168.6, "R": 173.4, "H": 153.2}

HYDROPHOBICITY_BY_AA = {"G": -0.4, "A": 1.8, "V": 4.2, "L": 3.8, "I": 4.5, "M": 1.9, "F": 2.8, "Y": -1.3, "W": -0.9,
                        "S": -0.8, "T": -0.7, "C": 2.5, "P": -1.6, "N": -3.5, "Q": -3.5, "D": -3.5, "E": -3.5,
                        "K": -3.9, "R": -4.5, "H": -3.2}


# This file contains functions that return properties of amino acids.

def get_aa_volume(aa: str) -> float:
    """Returns the volume of the given amino acid in Angstroms^3."""
    return VOLUME_BY_AA[aa]


def get_aa_hydrophobicity_by_kd_scale(aa: str) -> float:
    """Returns the hydrophobicity of the given amino acid, according to the Kyte-Doolittle scale."""
    return HYDROPHOBICITY_BY_AA[aa]
