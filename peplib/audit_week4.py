"""Audit T83 local. Seuls les agrégats du groupe admissible sont retournés."""
from collections import Counter

from scipy.stats import binomtest

from peplib.align import smith_waterman
from peplib.dedup import sequence_identity
from peplib.tokenizer_diffusion import encode


def memorization_report(records, train, flags):
    """Joint le registre v02 par séquence. Toute absence ou ambiguïté exclut.

    Calcule le maximum exhaustif de l'identité de Levenshtein sur les
    fragments locaux alignés (gaps retirés). Les bornes des déciles
    sont fermées à gauche, avec une catégorie distincte pour 1.
    """
    statuses = {}
    for sequence, status in flags[["sequence", "toxicity_status"]].itertuples(index=False, name=None):
        eligible = isinstance(status, str) and status == "not_flagged_by_uniprot"
        statuses[sequence] = statuses.get(sequence, True) and eligible
    # Le signalement d'un exemplaire prévaut pour tous ses doublons.
    sequences = []
    for sequence, flag in records[["sequence", "toxin_flag"]].itertuples(index=False, name=None):
        sequences.append(sequence)
        if flag is not False:
            statuses[sequence] = False
    references = set(train)
    if not references or "" in references:
        raise ValueError("Références train non vides requises")
    for reference in references:
        encode(reference)
    histogram = Counter()
    n = hits = excluded = 0
    for sequence in sequences:
        if not statuses.get(sequence, False):
            excluded += 1
            continue
        if not sequence:
            raise ValueError("Séquence vide")
        encode(sequence)
        best = 0.0
        if sequence in references:
            best = 1.0
        else:
            for reference in references:
                alignment = smith_waterman(sequence, reference)
                if alignment.score > 0:
                    best = max(best, sequence_identity(alignment.aligned_x.replace("-", ""),
                                                       alignment.aligned_y.replace("-", "")))
                if best == 1:
                    break
        histogram[int(best * 10)] += 1
        hits += int(best >= 0.90)
        n += 1
    interval = binomtest(hits, n).proportion_ci(method="wilson") if n else None
    return {"n_eligible": n, "n_excluded": excluded, "n_ge_90": hits,
            "fraction_ge_90": hits / n if n else None,
            "wilson_95": tuple(interval) if interval else None,
            "identity_deciles": [histogram[i] for i in range(11)]}
