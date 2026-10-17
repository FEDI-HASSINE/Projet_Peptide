import hashlib
import os
import tempfile

import pytest

from peplib.download import (
    count_lines,
    compute_sha256,
    download_paginated_tsv,
    next_page_url,
    page_body_lines,
    raw_path,
)


def test_compute_sha256():
    """Test unitaire de l'empreinte SHA-256 avec un contenu connu."""
    with tempfile.NamedTemporaryFile(
        delete=False, mode="w", encoding="utf-8"
    ) as tmp:
        tmp.write("test peptide sequence")
        tmp_path = tmp.name

    try:
        sha = compute_sha256(tmp_path)
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
        assert count_lines(tmp_path) == 4
    finally:
        os.remove(tmp_path)


def test_functions_come_from_library():
    """Garantit que les tests portent sur le vrai module, pas sur une copie."""
    assert compute_sha256.__module__ == "peplib.download"
    assert count_lines.__module__ == "peplib.download"


@pytest.mark.parametrize(
    "header, expected",
    [
        (None, None),
        ("", None),
        ('<https://x/?page=2>; rel="next"', "https://x/?page=2"),
        ('<https://x/?page=1>; rel="prev"', None),
        (
            '<https://x/?page=1>; rel="prev", <https://x/?page=3>; rel="next"',
            "https://x/?page=3",
        ),
    ],
)
def test_next_page_url(header, expected):
    assert next_page_url(header) == expected


def test_page_body_lines_first_page_keeps_header():
    assert page_body_lines(["h", "a", "b"], is_first=True) == ["h", "a", "b"]


def test_page_body_lines_later_page_drops_header():
    assert page_body_lines(["h", "a", "b"], is_first=False) == ["a", "b"]


def test_page_body_lines_later_page_without_body():
    assert page_body_lines(["h"], is_first=False) == []


def test_raw_path_uses_prefix_and_date():
    path = raw_path("uniprot_reviewed", "20261017")
    assert path.endswith(os.path.join("data", "raw", "uniprot_reviewed_20261017.tsv"))


class _FakeResponse:
    def __init__(self, text, link=None):
        self.text = text
        self.headers = {"Link": link} if link else {}

    def raise_for_status(self):
        pass


class _FakeSession:
    """Session factice rejouant une liste de pages dans l'ordre."""

    def __init__(self, pages):
        self.pages = list(pages)
        self.calls = 0

    def get(self, url, params=None, timeout=None):
        page = self.pages[self.calls]
        self.calls += 1
        return page


def test_download_paginated_tsv_writes_header_once(tmp_path):
    p1 = _FakeResponse("h\nA1\nA2", '<https://x/?page=2>; rel="next"')
    p2 = _FakeResponse("h\nA3")
    session = _FakeSession([p1, p2])
    out = tmp_path / "raw.tsv"

    total = download_paginated_tsv(
        "q", str(out), "accession,sequence", session=session
    )

    assert total == 3
    assert out.read_text(encoding="utf-8").splitlines() == ["h", "A1", "A2", "A3"]
    assert session.calls == 2


def test_download_paginated_tsv_single_page(tmp_path):
    session = _FakeSession([_FakeResponse("h\nA1")])
    out = tmp_path / "raw.tsv"

    assert download_paginated_tsv("q", str(out), "f", session=session) == 1
    assert out.read_text(encoding="utf-8").splitlines() == ["h", "A1"]
