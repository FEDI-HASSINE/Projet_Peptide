import pandas as pd
import pytest

from peplib.download_challenge import (
    SOURCE,
    load_challenge,
    parse_header,
    prepare_challenge,
    write_processed,
)


FASTA_REFERENCE = (
    ">AMP001 len=12 charge=2.0\nACDEFGHIKLMN\n"
    ">AMP002 len=4 charge=0\nACDE\n"
)


def _write(tmp_path, content, name="test.fasta"):
    fasta = tmp_path / name
    fasta.write_text(content, encoding="utf-8")
    return fasta


def test_parse_header():
    result = parse_header(
        ">MLAMP0003076 len=15 charge=4.73 disulfide=0 "
        "dbs=AMPDB|APD activity=antibacterial|antimicrobial"
    )
    assert result["id"] == "MLAMP0003076"
    assert result["len"] == "15"
    assert result["charge"] == "4.73"
    assert result["dbs"] == "AMPDB|APD"
    assert result["activity"] == "antibacterial|antimicrobial"


def test_parse_header_vide():
    with pytest.raises(ValueError):
        parse_header(">")


def test_load_and_filter(tmp_path):
    fasta = _write(tmp_path, FASTA_REFERENCE)
    raw = load_challenge(fasta)
    cleaned, rejections = prepare_challenge(fasta)

    assert len(raw) == 2
    assert len(cleaned) == 1
    assert rejections["trop_courte"] == 1
    assert len(cleaned) + sum(rejections.values()) == len(raw)
    assert set(raw["source"]) == {SOURCE}
    assert "dbs_amont" in raw.columns


def test_dbs_conserve_en_colonne_dbs_amont(tmp_path):
    fasta = _write(
        tmp_path,
        ">AMP001 len=12 charge=2.0 dbs=APD|CAMP\nACDEFGHIKLMN\n",
    )
    df = load_challenge(fasta)
    assert df.iloc[0]["dbs_amont"] == "APD|CAMP"
    assert "dbs" not in df.columns


def test_entete_incomplet_metadonnees_nan(tmp_path):
    fasta = _write(tmp_path, ">AMP003\nACDEFGHIKLMN\n")
    df = load_challenge(fasta)
    assert len(df) == 1
    assert df.iloc[0]["id"] == "AMP003"
    assert pd.isna(df.iloc[0]["len"])
    assert pd.isna(df.iloc[0]["dbs_amont"])


def test_texte_avant_entete_leve(tmp_path):
    fasta = _write(tmp_path, "ACDEFGHIKLMN\n>AMP001\nACDEFGHIKLMN\n")
    with pytest.raises(ValueError):
        load_challenge(fasta)


def test_commentaire_avant_entete_tolere(tmp_path):
    fasta = _write(tmp_path, "; commentaire\n>AMP001 len=12\nACDEFGHIKLMN\n")
    df = load_challenge(fasta)
    assert len(df) == 1


def test_write_processed(tmp_path):
    df = pd.DataFrame({
        "sequence": ["ACDEFGHIKLMN"],
        "source": [SOURCE],
        "id": ["AMP001"],
        "len": ["12"],
        "charge": ["2.0"],
        "disulfide": ["0"],
        "activity": ["antibacterial"],
        "dbs_amont": ["APD|CAMP"],
    })
    out = write_processed(df, tmp_path / "sous_dossier" / "challenge.csv")
    assert out.exists()
    relu = pd.read_csv(out, dtype=str)
    assert list(relu.columns) == [
        "sequence", "source", "id", "len", "charge",
        "disulfide", "activity", "dbs_amont",
    ]
    assert relu.iloc[0]["dbs_amont"] == "APD|CAMP"
