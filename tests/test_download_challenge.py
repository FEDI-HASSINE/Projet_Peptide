from peplib.download_challenge import parse_header, load_challenge, prepare_challenge

def test_parse_header():
    result = parse_header(
        ">MLAMP0003076 len=15 charge=4.73 disulfide=0 "
        "dbs=AMPDB|APD activity=antibacterial|antimicrobial"
    )
    assert result["id"] == "MLAMP0003076"
    assert result["len"] == "15"
    assert result["charge"] == "4.73"
    assert result["dbs"] == "AMPDB|APD"

def test_load_and_filter(tmp_path):
    fasta = tmp_path / "test.fasta"
    fasta.write_text(
        ">AMP001 len=12 charge=2.0\nACDEFGHIKLMN\n"
        ">AMP002 len=4 charge=0\nACDE\n",
        encoding="utf-8"
    )
    raw = load_challenge(fasta)
    cleaned, rejections = prepare_challenge(fasta)

    assert len(raw) == 2
    assert len(cleaned) == 1
    assert rejections["trop_courte"] == 1
    assert len(cleaned) + sum(rejections.values()) == len(raw)
