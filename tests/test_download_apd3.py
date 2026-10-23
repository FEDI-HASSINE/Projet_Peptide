from peplib.download_apd3 import load_apd3, prepare_apd3

def test_load_apd3(tmp_path):
    fasta = tmp_path / "apd.fasta"
    fasta.write_text(
        ">Your search led to 3 peptides\n"
        ">AP00001\nACDEFGHIKLMN\n\n"
        ">AP00002\nACDE\n"
        ">AP00003\nACDEFGHIKLMNPQRSTVWY\n",
        encoding="utf-8",
    )
    df = load_apd3(fasta)
    assert len(df) == 3
    assert df.iloc[0]["id"] == "AP00001"
    assert df.iloc[0]["sequence"] == "ACDEFGHIKLMN"

def test_prepare_apd3(tmp_path):
    fasta = tmp_path / "apd.fasta"
    fasta.write_text(
        ">Your search led to 3 peptides\n"
        ">AP00001\nACDEFGHIKLMN\n"
        ">AP00002\nACDE\n"
        ">AP00003\nACDEFGHIKLMNPQRSTVWY\n",
        encoding="utf-8",
    )
    cleaned, rejections = prepare_apd3(fasta)
    assert len(cleaned) == 2
    assert rejections["trop_courte"] == 1
    assert len(cleaned) + sum(rejections.values()) == 3

def test_modified_residue_rejected(tmp_path):
    fasta = tmp_path / "apd.fasta"
    fasta.write_text(
        ">AP00001\nACDEFGHIKLMN\n"
        ">AP00002\nACDEFGHIKLMX\n",
        encoding="utf-8",
    )
    cleaned, rejections = prepare_apd3(fasta)
    assert len(cleaned) == 1
    assert rejections["lettres_ambigues"] == 1
