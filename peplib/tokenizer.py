"""
Tokenisation des peptides (T40) : sequence de lettres <-> liste d'entiers.

Vocabulaire, V = 23 (ids FIXES : ne jamais changer l'ordre) :
    0 = <PAD>   remplissage pour aligner des sequences de longueurs differentes
    1 = <BOS>   debut de sequence
    2 = <EOS>   fin de sequence
    3..22       les 20 acides amines dans l'ordre alphabetique ACDEFGHIKLMNPQRSTVWY

encode("AG")  -> [1, 3, 8, 2]        # BOS, A, G, EOS
decode([1, 3, 8, 2]) -> "AG"

Les lettres hors alphabet (X, U, B, minuscules...) sont REFUSEES (ValueError) :
elles sont censees etre rejetees en amont, au filtrage du corpus.

Regle pour la semaine 3 : un modele ne doit jamais tirer <PAD> ni <BOS> en sortie ;
FORBIDDEN_OUTPUT_IDS donne les ids dont les logits seront mis a -inf a l'echantillonnage.
"""
from peplib.constants import AA_LIST

PAD_TOKEN, BOS_TOKEN, EOS_TOKEN = "<PAD>", "<BOS>", "<EOS>"
PAD_ID, BOS_ID, EOS_ID = 0, 1, 2

_TOKENS = [PAD_TOKEN, BOS_TOKEN, EOS_TOKEN] + list(AA_LIST)
VOCAB = {token: i for i, token in enumerate(_TOKENS)}          # token -> id
ID_TO_TOKEN = {i: token for token, i in VOCAB.items()}          # id -> token
AA_TO_ID = {aa: VOCAB[aa] for aa in AA_LIST}
FORBIDDEN_OUTPUT_IDS = (PAD_ID, BOS_ID)


def encode(sequence: str) -> list:
    """BOS + ids des acides amines + EOS."""
    bad = sorted({c for c in sequence if c not in AA_TO_ID})
    if bad:
        raise ValueError(f"lettres hors alphabet : {bad}")
    return [BOS_ID] + [AA_TO_ID[c] for c in sequence] + [EOS_ID]


def pad(ids: list, length: int) -> list:
    """Complete avec PAD a droite jusqu'a `length` (refuse de tronquer)."""
    if len(ids) > length:
        raise ValueError(f"sequence plus longue que {length} : {len(ids)}")
    return list(ids) + [PAD_ID] * (length - len(ids))


def decode(ids: list) -> str:
    """Inverse de encode. Format strict : BOS, acides amines, EOS, puis PAD eventuels."""
    ids = list(ids)
    if not ids or ids[0] != BOS_ID:
        raise ValueError("la sequence doit commencer par BOS")
    letters, i = [], 1
    while i < len(ids) and ids[i] != EOS_ID:
        token = ID_TO_TOKEN.get(ids[i])
        if token is None or ids[i] in (PAD_ID, BOS_ID):
            raise ValueError(f"id invalide a la position {i} : {ids[i]}")
        letters.append(token)
        i += 1
    if i == len(ids):
        raise ValueError("EOS manquant")
    if any(x != PAD_ID for x in ids[i + 1:]):
        raise ValueError("contenu autre que PAD apres EOS")
    return "".join(letters)