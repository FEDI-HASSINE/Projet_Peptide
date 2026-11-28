import pytest
import torch

from peplib.datasets import make_loader
from peplib.diffusion_data import DiffusionDataset, mask_batch
from peplib.tokenizer_diffusion import BOS_ID, EOS_ID, MASK_ID, PAD_ID
from scripts.check_masking_s4 import measurements


def test_frequency():
    assert all(row["passed"] for row in measurements())


def test_special_tokens_and_time():
    generator = torch.Generator().manual_seed(42)
    x = torch.tensor([[BOS_ID, 0, EOS_ID, PAD_ID]])
    times = []
    for _ in range(10_000):
        xt, mask, t = mask_batch(x, generator=generator)
        times.append(t)
        assert xt.dtype == torch.long and mask.dtype == torch.bool
        assert not mask[0, 0] and not mask[0, -1]
        assert torch.equal(xt, x.masked_fill(mask, MASK_ID))
    assert set(times) == set(range(1, 51))
    assert mask_batch(x, 50)[1].tolist() == [[False, True, True, False]]
    assert torch.equal(mask_batch(x, 0)[0], x)


def test_padding_and_seed():
    x = torch.tensor([[4, EOS_ID], [0, EOS_ID]])
    padded = torch.nn.functional.pad(x, (0, 5), value=PAD_ID)
    a, m, t = mask_batch(x, generator=torch.Generator().manual_seed(7))
    b, n, u = mask_batch(padded, generator=torch.Generator().manual_seed(7))
    assert t == u and torch.equal(a, b[:, :2]) and torch.equal(m, n[:, :2])
    assert not n[:, 2:].any()


@pytest.mark.parametrize("x,t", [(torch.ones(2), 1), (torch.ones((2, 2)), 1),
                                (torch.tensor([[23]]), 2), (torch.tensor([[0]]), 51)])
def test_invalid(x, t):
    with pytest.raises(ValueError):
        mask_batch(x, t)


def test_loader(tmp_path):
    (tmp_path / "val_v02.csv").write_text("sequence\nA\nAC\n", encoding="utf-8")
    batch = next(iter(make_loader(DiffusionDataset("val", tmp_path), batch_size=2)))
    assert batch.shape == (2, 51) and batch.dtype == torch.long
    assert batch[0, :3].tolist() == [0, EOS_ID, PAD_ID]
    with pytest.raises(ValueError):
        DiffusionDataset("test", tmp_path)
