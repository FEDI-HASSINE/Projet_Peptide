"""Vocabulaire S4 indépendant des identifiants historiques du modèle causal."""
from peplib.constants import AA_LIST

PAD_ID, BOS_ID, EOS_ID, MASK_ID = range(20, 24)
V_DIFF = 24
MAX_RESIDUES = 50
VOCAB = {token: i for i, token in enumerate((*AA_LIST, "<PAD>", "<BOS>", "<EOS>", "[MASK]"))}


def encode(sequence: str) -> list[int]:
    """Encode une séquence propre entre BOS et EOS."""
    if len(sequence) > MAX_RESIDUES:
        raise ValueError("Séquence trop longue")
    try:
        return [BOS_ID, *(VOCAB[letter] for letter in sequence), EOS_ID]
    except KeyError as exc:
        raise ValueError("Résidu hors alphabet") from exc


def decode(ids) -> str:
    """Décode une séquence complète, éventuellement paddée à droite."""
    ids = list(ids)
    while ids and ids[-1] == PAD_ID:
        ids.pop()
    if len(ids) < 2 or ids[0] != BOS_ID or ids[-1] != EOS_ID:
        raise ValueError("BOS et EOS attendus")
    residues = ids[1:-1]
    if len(residues) > MAX_RESIDUES or any(i not in range(20) for i in residues):
        raise ValueError("Séquence de résidus invalide")
    return "".join(AA_LIST[i] for i in residues)
