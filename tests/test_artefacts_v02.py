"""Verrous des artefacts du corpus v0.2.

Les comptes et l'empreinte ci-dessous verrouillent les artefacts committés :
toute re-génération qui change silencieusement une de ces valeurs fait
échouer la suite. Toute variation doit être consignée dans PROVENANCE.md
et les tests mis à jour en même temps (relecture PR #19, point W3).
"""

import hashlib
from pathlib import Path

import pandas as pd

PROCESSED = Path("data/processed")


def _sha256(filename):
    h = hashlib.sha256()
    with open(PROCESSED / filename, "rb") as fh:
        for block in iter(lambda: fh.read(65536), b""):
            h.update(block)
    return h.hexdigest()


def test_lock_corpus_v02():
    corpus = pd.read_csv(PROCESSED / "corpus_v0.2.csv", dtype=str)

    assert list(corpus.columns)[:3] == ["sequence", "sources", "ids"]
    assert corpus["sequence"].notna().all()
    assert corpus["sequence"].is_unique
    assert len(corpus) == 42_575


def test_lock_source_views_v02():
    expected = {
        "view_uniprot_v02.csv": 9_021,
        "view_apd3_v02.csv": 2_796,
        "view_ampchallenge_v02.csv": 34_474,
    }

    total = 0
    for name, count in expected.items():
        view = pd.read_csv(PROCESSED / name, dtype=str)
        assert len(view) == count, name
        assert view["sequence"].is_unique, name
        total += count

    # 46 291 appartenances pour 42 575 uniques : les 3 716 recouvrements
    # multi-sources sont documentés dans PROVENANCE.md.
    assert total == 46_291


def test_lock_toxin_flags_v02_counts():
    flags = pd.read_csv(PROCESSED / "toxin_flags_v02.csv", dtype=str)

    assert list(flags.columns) == [
        "sequence",
        "sources",
        "ids",
        "has_uniprot_annotation",
        "is_toxin_annotation",
        "toxicity_status",
    ]
    assert len(flags) == 42_575

    counts = flags["toxicity_status"].value_counts()
    assert counts["flagged"] == 1_786
    assert counts["not_flagged_by_uniprot"] == 7_125
    assert counts["unknown"] == 33_664
    assert int(counts.sum()) == 42_575


def test_lock_toxin_flags_v02_hash():
    # L'empreinte porte sur le fichier committé tel quel (octets, LF, sans BOM).
    # Toute re-génération ou normalisation modifie le hash : elle doit être
    # consignée dans PROVENANCE.md avec une note de relecture.
    assert (
        _sha256("toxin_flags_v02.csv")
        == "171f150b7a41871c2d10b762756d8cf7ed1aae4fbbd69981557ec6056180d314"
    )
