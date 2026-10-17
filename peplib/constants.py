"""
Références physico-chimiques des 20 acides aminés standards.
"""

# Les 20 acides aminés standards, en code à une lettre.
AMINO_ACIDS = tuple("ACDEFGHIKLMNPQRSTVWY")
AA_LIST = AMINO_ACIDS

# Masse monoisotopique du résidu (Da), sans H2O.
# Source : ExPASy FindMod, Amino acids — Monoisotopic masses.
# https://proteinsandproteomics.org/content/free/tables_1/table07.pdf
MONOISOTOPIC_RESIDUE_MASS = {
    "A": 71.03711, "R": 156.10111, "N": 114.04293, "D": 115.02694,
    "C": 103.00919, "E": 129.04259, "Q": 128.05858, "G": 57.02146,
    "H": 137.05891, "I": 113.08406, "L": 113.08406, "K": 128.09496,
    "M": 131.04049, "F": 147.06841, "P": 97.05276, "S": 87.03203,
    "T": 101.04768, "W": 186.07931, "Y": 163.06333, "V": 99.06841,
}
MASS = MONOISOTOPIC_RESIDUE_MASS

# pKa des chaînes latérales ionisables.
# Valeurs conformes au cahier des charges du projet.
# Référence complémentaire : article comparatif sur les valeurs de pKa.
# https://pmc.ncbi.nlm.nih.gov/articles/PMC5075173/  (table 4 , ligne wikepidea)
# Les pKa peuvent varier selon la source et l'environnement moléculaire.
SIDE_CHAIN_PKA = {
    "C": 8.18,
    "D": 3.9,
    "E": 4.07,
    "H": 6.04,
    "K": 10.54,
    "R": 12.48,
    "Y": 10.46,
}

# pKa des extrémités d'un peptide, selon les valeurs imposées par le projet.
# Valeurs de référence : N-terminal = 9.69 ; C-terminal = 2.34.
# Référence : https://www.bachem.com/knowledge-center/peptide-calculator/
N_TERMINUS_PKA = 9.69
C_TERMINUS_PKA = 2.34

# Échelle d'hydrophobicité de Kyte-Doolittle.
# Kyte J. & Doolittle R. F. (1982), Journal of Molecular Biology,
# 157, 105-132. DOI: 10.1016/0022-2836(82)90515-0
# https://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/650/Hydrophobicity_scales.html
KYTE_DOOLITTLE = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# Échelle d'hydrophobicité d'Eisenberg (consensus scale).
# Eisenberg D., Schwarz E., Komaromy M. & Wall R. (1984),
# Journal of Molecular Biology, 179, 125-142.
# https://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/650/Hydrophobicity_scales.html
EISENBERG = {
    "A": 0.62, "R": -2.53, "N": -0.78, "D": -0.90, "C": 0.29,
    "Q": -0.85, "E": -0.74, "G": 0.48, "H": -0.40, "I": 1.38,
    "L": 1.06, "K": -1.50, "M": 0.64, "F": 1.19, "P": 0.12,
    "S": -0.18, "T": -0.05, "W": 0.81, "Y": 0.26, "V": 1.08,
}
