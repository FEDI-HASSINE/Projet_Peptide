"""Initialisation reproductible des generateurs aleatoires.

La variable CUBLAS_WORKSPACE_CONFIG doit etre posee AVANT l'initialisation
du contexte CUDA : elle est donc definie au niveau du module, avant
l'import de torch. Sans elle, torch.use_deterministic_algorithms(True)
refuse les operations cuBLAS.
"""
import os

os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")

import random  # noqa: E402

import numpy as np  # noqa: E402
import torch  # noqa: E402


def set_seed(seed: int = 42) -> None:
    """Fixe toutes les graines et impose le mode deterministe."""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)

    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)

    # Pas de noyau non deterministe (matmul, convolutions, cudnn).
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    torch.use_deterministic_algorithms(True)
