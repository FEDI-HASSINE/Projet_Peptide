import pandas as pd

from peplib.download_apd3 import (
    SOURCE,
    load_apd3,
    prepare_apd3,
    write_processed,
)


FASTA_REFERENCE = (
    ">Your search led to 3 peptides\n"
    ">AP00001\nACDEFGHIKLMN\n"
    "\n"
    ">AP00002\nACDE\n"
    ">AP00003\n" + "ACDEFGHIKLMNPQRSTVWY" * 3 + "\n"
)


def _write(tmp_path, content, name="apd.fasta"):
    fasta = tmp_path / name
    fasta.write_text(content, encoding="utf-8")
    return fasta


def test_load_apd3(tmp_path):
    fasta = _write(tmp_path, FASTA_REFERENCE)
    df = load_apd3(fasta)
    assert len(df) == 3
    assert list(df.columns) == ["sequence", "source", "id"]
    assert df.iloc[0]["id"] == "AP00001"
    assert df.iloc[0]["sequence"] == "ACDEFGHIKLMN"
    assert set(df["source"]) == {SOURCE}


def test_annonce_site_ignoree(tmp_path):
    fasta = _write(tmp_path, FASTA_REFERENCE)
    df = load_apd3(fasta)
    assert "Your search led to 3 peptides" not in set(df["id"])
    assert df["id"].notna().all()


def test_prepare_apd3(tmp_path):
    fasta = _write(tmp_path, FASTA_REFERENCE)
    cleaned, rejections = prepare_apd3(fasta)
    assert len(cleaned) == 1
    assert rejections["trop_courte"] == 1
    assert rejections["trop_longue"] == 1
    assert len(cleaned) + sum(rejections.values()) == 3


def test_modified_residue_rejected(tmp_path):
    fasta = _write(
        tmp_path,
        ">AP00001\nACDEFGHIKLMN\n"
        ">AP00002\nACDEFGHIKLMX\n",
    )
    cleaned, rejections = prepare_apd3(fasta)
    assert len(cleaned) == 1
    assert rejections["lettres_ambigues"] == 1


def test_entete_non_ap_ignore(tmp_path):
    fasta = _write(
        tmp_path,
        ">NotAnApIdentifier\nACDEFGHIKLMN\n"
        ">AP00001\nACDEFGHIKLMN\n",
    )
    df = load_apd3(fasta)
    assert len(df) == 1
    assert df.iloc[0]["id"] == "AP00001"


def test_sequence_avant_entete_ignoree(tmp_path):
    fasta = _write(
        tmp_path,
        "ACDEFGHIKLMN\n"
        ">AP00001\nACDEFGHIKLMN\n",
    )
    df = load_apd3(fasta)
    assert len(df) == 1
    assert df.iloc[0]["sequence"] == "ACDEFGHIKLMN"


def test_write_processed(tmp_path):
    df = pd.DataFrame({
        "sequence": ["ACDEFGHIKLMN"],
        "source": [SOURCE],
        "id": ["AP00001"],
    })
    out = write_processed(df, tmp_path / "sous_dossier" / "amp_apd3.csv")
    assert out.exists()
    relu = pd.read_csv(out, dtype=str)
    assert list(relu.columns) == ["sequence", "source", "id"]
    assert relu.iloc[0]["id"] == "AP00001"
