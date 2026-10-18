"""
Descripteurs physico-chimiques des peptides.
"""

import math

from .constants import (
    AMINO_ACIDS,
    BOMAN,
    MONOISOTOPIC_RESIDUE_MASS,
    SIDE_CHAIN_PKA,
    N_TERMINUS_PKA,
    C_TERMINUS_PKA,
    KYTE_DOOLITTLE,
    EISENBERG,
)


def calculate_charge(sequence: str, pH: float = 7.0) -> float:
    """
    Calcule la charge électrique théorique d'un peptide à un pH donné.

    Args:
        sequence: Séquence du peptide, par exemple "KKD".
        pH: pH du milieu, égal à 7.0 par défaut.

    Returns:
        Charge totale théorique du peptide.
    """
    if not isinstance(sequence, str) or not sequence:
        raise ValueError("La séquence doit être une chaîne non vide.")

    sequence = sequence.upper()

    if any(aa not in AMINO_ACIDS for aa in sequence):
        raise ValueError("La séquence contient un acide aminé inconnu.")

    if isinstance(pH, bool) or not isinstance(pH, (int, float)):
        raise ValueError("Le pH doit être un nombre.")

    if not math.isfinite(pH):
        raise ValueError("Le pH doit être un nombre fini.")

    charge = 0.0

    # Extrémité N-terminale : contribution positive.
    charge += 1.0 / (1.0 + 10.0 ** (pH - N_TERMINUS_PKA))

    # Extrémité C-terminale : contribution négative.
    charge -= 1.0 / (1.0 + 10.0 ** (C_TERMINUS_PKA - pH))

    # Contributions des chaînes latérales ionisables.
    for aa in sequence:
        if aa not in SIDE_CHAIN_PKA:
            continue

        pka = SIDE_CHAIN_PKA[aa]

        if aa in ("K", "R", "H"):
            # Groupes positifs.
            charge += 1.0 / (1.0 + 10.0 ** (pH - pka))
        else:
            # Groupes négatifs : D, E, C et Y.
            charge -= 1.0 / (1.0 + 10.0 ** (pka - pH))

    return charge


def calculate_pi(
    sequence: str,
    tolerance: float = 1e-6,
    max_iterations: int = 200,
) -> float:
    """
    Calcule le point isoélectrique théorique d'un peptide
    par dichotomie.

    Args:
        sequence: Séquence du peptide.
        tolerance: Précision souhaitée sur le pH.
        max_iterations: Nombre maximal d'itérations.

    Returns:
        pH auquel la charge totale est proche de zéro.
    """
    if not isinstance(sequence, str) or not sequence:
        raise ValueError("La séquence doit être une chaîne non vide.")

    sequence = sequence.upper()

    if any(aa not in AMINO_ACIDS for aa in sequence):
        raise ValueError("La séquence contient un acide aminé inconnu.")

    if not isinstance(tolerance, (int, float)) or isinstance(tolerance, bool):
        raise ValueError("La tolérance doit être un nombre positif.")

    if not math.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("La tolérance doit être un nombre positif.")

    if isinstance(max_iterations, bool) or not isinstance(max_iterations, int):
        raise ValueError("max_iterations doit être un entier positif.")

    if max_iterations <= 0:
        raise ValueError("max_iterations doit être un entier positif.")

    low = 0.0
    high = 14.0

    charge_low = calculate_charge(sequence, low)
    charge_high = calculate_charge(sequence, high)

    if charge_low < 0 or charge_high > 0:
        raise ValueError(
            "Impossible d'encadrer le point isoélectrique entre pH 0 et 14."
        )

    for _ in range(max_iterations):
        middle = (low + high) / 2.0
        charge_middle = calculate_charge(sequence, middle)

        if abs(charge_middle) < 1e-12 or high - low <= tolerance:
            return middle

        if charge_middle > 0:
            low = middle
        else:
            high = middle

    return (low + high) / 2.0

def calculate_gravy(sequence: str) -> float:
    """
    Calcule le GRAVY, c'est-à-dire l'hydrophobicité moyenne
    d'un peptide selon l'échelle de Kyte-Doolittle.
    """
    if not isinstance(sequence, str) or not sequence:
        raise ValueError("La séquence doit être une chaîne non vide.")

    sequence = sequence.upper()

    if any(aa not in AMINO_ACIDS for aa in sequence):
        raise ValueError("La séquence contient un acide aminé inconnu.")

    total = sum(KYTE_DOOLITTLE[aa] for aa in sequence)

    return total / len(sequence)



def calculate_hydrophobic_moment(
    sequence: str,
    angle: float = 100.0,
) -> float:
    """
    Calcule le moment hydrophobe théorique d'un peptide
    en supposant une hélice alpha.
    """
    if not isinstance(sequence, str) or not sequence:
        raise ValueError("La séquence doit être une chaîne non vide.")

    sequence = sequence.upper()

    if any(aa not in AMINO_ACIDS for aa in sequence):
        raise ValueError("La séquence contient un acide aminé inconnu.")

    if isinstance(angle, bool) or not isinstance(angle, (int, float)):
        raise ValueError("L'angle doit être un nombre.")

    if not math.isfinite(angle):
        raise ValueError("L'angle doit être un nombre fini.")

    angle_rad = math.radians(angle)

    somme_x = 0.0
    somme_y = 0.0

    for i, aa in enumerate(sequence):
        h = EISENBERG[aa]
        theta = i * angle_rad

        somme_x += h * math.cos(theta)
        somme_y += h * math.sin(theta)

    moment = math.sqrt(somme_x**2 + somme_y**2)

    return moment / len(sequence)



def calculate_boman_index(sequence: str) -> float:
    """
    Calcule l'indice de Boman moyen d'un peptide.
    """
    if not isinstance(sequence, str) or not sequence:
        raise ValueError("La séquence doit être une chaîne non vide.")

    sequence = sequence.upper()

    if any(aa not in AMINO_ACIDS for aa in sequence):
        raise ValueError("La séquence contient un acide aminé inconnu.")

    total = sum(BOMAN[aa] for aa in sequence)
    return total / len(sequence)




def calculate_molecular_mass(sequence: str) -> float:
    """
    Calcule la masse monoisotopique d'un peptide en Da.
    """
    if not isinstance(sequence, str) or not sequence:
        raise ValueError("La séquence doit être une chaîne non vide.")

    sequence = sequence.upper()

    if any(aa not in AMINO_ACIDS for aa in sequence):
        raise ValueError("La séquence contient un acide aminé inconnu.")

    residue_mass = sum(
        MONOISOTOPIC_RESIDUE_MASS[aa] for aa in sequence
    )

    water_mass = 18.01056

    return residue_mass + water_mass
