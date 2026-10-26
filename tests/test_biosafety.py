import json

import pandas as pd
import pytest

from peplib.biosafety import (
    TOXIN_PATTERN,
    detect_id_column,
    flag_toxins,
    load_raw_tsv,
    run,
    write_toxin_flags,
)


@pytest.fixture
def sample_df():
    return pd.DataFrame(
        {
            "Entry": ["A1", "B2", "C3"],
            "Protein.name": [
                "Alpha-conotoxin OmIA",
                "Unrelated kinase",
                "Bacteriocin plantaricin",
            ],
        }
    )


def test_detect_id_column_entry():
    df = pd.DataFrame({"Entry": ["x"], "Accession": ["y"]})
    assert detect_id_column(df) == "Entry"


def test_detect_id_column_accession_is_fallback():
    df = pd.DataFrame({"Accession": ["y"], "Protein.name": ["z"]})
    assert detect_id_column(df) == "Accession"


def test_detect_id_column_first_column_is_last_resort():
    df = pd.DataFrame({"Foo": ["y"], "Protein.name": ["z"]})
    assert detect_id_column(df) == "Foo"


def test_flag_toxins_detects_conotoxin(sample_df):
    flagged = flag_toxins(sample_df)
    assert flagged["is_toxin"].tolist() == [True, False, False]


def test_flag_toxins_is_case_insensitive():
    df = pd.DataFrame({"Entry": ["X"], "Protein.name": ["DELTA TOXIN" ]})
    assert flag_toxins(df)["is_toxin"].tolist() == [True]


def test_flag_toxins_does_not_mutate_input(sample_df):
    columns_before = list(sample_df.columns)
    flag_toxins(sample_df)
    assert list(sample_df.columns) == columns_before
    assert "is_toxin" not in sample_df.columns


def test_flag_toxins_handles_missing_values():
    df = pd.DataFrame({"Entry": ["X"], "Protein.name": [None]})
    # fillna n'est pas appliqué ici : aucun crash attendu, juste False
    assert flag_toxins(df)["is_toxin"].tolist() == [False]


def test_toxin_pattern_matches_conotoxin_alone():
    # "Toxin" couvre "Conotoxin" : le motif ne doit pas rater une conotoxine.
    assert pd.Series(["conotoxin"]).str.contains(
        TOXIN_PATTERN, case=False
    ).tolist() == [True]


def test_flag_toxins_detects_venom():
    # Relecture PR #19 (B2) : les annotations de venin doivent être signalées.
    df = pd.DataFrame(
        {"Entry": ["X"], "Protein.name": ["Cysteine-rich venom protein"]}
    )
    assert flag_toxins(df)["is_toxin"].tolist() == [True]


def test_flag_toxins_detects_toxic():
    df = pd.DataFrame({"Entry": ["X"], "Protein.name": ["Toxic peptide"]})
    assert flag_toxins(df)["is_toxin"].tolist() == [True]


def test_flag_toxins_over_inclusive_is_fail_safe():
    # Mode d'échec sûr : même "non-toxic" est signalé (doute = flagué).
    df = pd.DataFrame({"Entry": ["X"], "Protein.name": ["Non-toxic marker"]})
    assert flag_toxins(df)["is_toxin"].tolist() == [True]


def test_write_toxin_flags_produces_expected_records(tmp_path, sample_df):
    out = tmp_path / "flags.json"
    write_toxin_flags(flag_toxins(sample_df), str(out))

    records = json.loads(out.read_text(encoding="utf-8"))
    assert records == [
        {"Entry": "A1", "is_toxin": True},
        {"Entry": "B2", "is_toxin": False},
        {"Entry": "C3", "is_toxin": False},
    ]


def test_run_end_to_end(tmp_path, sample_df):
    raw = tmp_path / "raw.tsv"
    sample_df.to_csv(raw, sep="\t", index=False)
    out = tmp_path / "processed" / "flags.json"

    count = run(str(raw), str(out))

    assert count == 1
    assert out.exists()
    assert len(json.loads(out.read_text(encoding="utf-8"))) == 3


def test_load_raw_tsv_fills_na_with_empty_string(tmp_path):
    raw = tmp_path / "raw.tsv"
    raw.write_text("Entry\tProtein.name\nA1\t\nB2\tSome toxin\n", encoding="utf-8")
    df = load_raw_tsv(str(raw))
    assert (df["Protein.name"] == "").tolist() == [True, False]
