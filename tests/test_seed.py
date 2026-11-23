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
def test_deterministic_algorithms_enabled():
    set_seed(42)
    assert torch.are_deterministic_algorithms_enabled()


def test_cudnn_deterministic_flags():
    set_seed(42)
    assert torch.backends.cudnn.deterministic is True
    assert torch.backends.cudnn.benchmark is False


def test_cublas_workspace_config_set_before_cuda():
    import os
    import peplib.seed  # noqa: F401

    assert os.environ.get("CUBLAS_WORKSPACE_CONFIG") in (":4096:8", ":16:8")
