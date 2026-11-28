import csv
from pathlib import Path

import pytest

from peplib import tokenizer as s3
from peplib import tokenizer_diffusion as s4


def test_vocabularies():
    assert (len(s3.VOCAB), s3.PAD_ID, s3.BOS_ID, s3.EOS_ID) == (23, 0, 1, 2)
    assert (s4.V_DIFF, s4.PAD_ID, s4.BOS_ID, s4.EOS_ID, s4.MASK_ID) == (24, 20, 21, 22, 23)
    assert list(s4.VOCAB.values()) == list(range(24))


def test_round_trip_full_train():
    path = Path(__file__).resolve().parents[1] / "data/splits/train_v02.csv"
    with path.open(encoding="utf-8", newline="") as stream:
        count = 0
        for record in csv.DictReader(stream):
            sequence = record["sequence"]
            ids = s4.encode(sequence)
            assert s4.MASK_ID not in ids
            assert s4.decode([*ids, s4.PAD_ID]) == s3.decode(s3.encode(sequence)) == sequence
            count += 1
        assert count == 38336


@pytest.mark.parametrize("ids", [[21, 23, 22], [21, 22, 1], [21], [21, 20, 22]])
def test_reject_incomplete(ids):
    with pytest.raises(ValueError):
        s4.decode(ids)
