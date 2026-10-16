import random
import numpy as np
import torch
from peplib.seed import set_seed


def _draw():
    return random.random(), np.random.rand(3).tolist(), torch.rand(3).tolist()


def test_set_seed_reproducible():
    set_seed(42)
    a = _draw()
    set_seed(42)
    b = _draw()
    assert a == b


def test_different_seed_differs():
    set_seed(42)
    a = _draw()
    set_seed(43)
    assert a != _draw()