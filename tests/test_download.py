import hashlib
import os
import tempfile

from download_t8_t9 import compute_sha256, count_lines


def test_compute_sha256():
    """Test unitaire de l'empreinte SHA-256 avec un contenu connu."""
    with tempfile.NamedTemporaryFile(
        delete=False, mode="w", encoding="utf-8"
    ) as tmp:
        tmp.write("test peptide sequence")
        tmp_path = tmp.name

    try:
        sha = compute_sha256(tmp_path)
        # Hash SHA-256 connu de la chaîne "test peptide sequence"
        expected_sha = hashlib.sha256(
            "test peptide sequence".encode("utf-8")
        ).hexdigest()
        assert sha == expected_sha
        assert len(sha) == 64
    finally:
        os.remove(tmp_path)


def test_count_lines():
    """Test unitaire du comptage de lignes d'un fichier."""
    with tempfile.NamedTemporaryFile(
        delete=False, mode="w", encoding="utf-8"
    ) as tmp:
        tmp.write("header\nseq1\nseq2\nseq3\n")
        tmp_path = tmp.name

    try:
        lines = count_lines(tmp_path)
        assert lines == 4
    finally:
        os.remove(tmp_path)


def test_functions_come_from_script():
    """Garantit que les tests portent sur le vrai module, pas sur une copie."""
    assert compute_sha256.__module__ == "download_t8_t9"
    assert count_lines.__module__ == "download_t8_t9"
