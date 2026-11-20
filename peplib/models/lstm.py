"""
peplib/models/lstm.py — LSTM causal pour la génération de peptides (T52).

Architecture
------------
  Embedding(vocab_size, embed_dim, padding_idx=PAD_ID)
  Dropout(dropout)
  LSTM(embed_dim, hidden_dim, n_layers, batch_first=True, dropout=dropout)
  Linear(hidden_dim, vocab_size)

Entrée  : (B, L) LongTensor d'ids de tokens
Sortie  : (B, L, V) FloatTensor de logits

La causalité est garantie par la nature séquentielle des LSTM : le logit à la
position i ne dépend que des tokens 0..i (fenêtre causale gauche).

Aucun poids pré-entraîné — initialisation PyTorch par défaut (Xavier/Kaiming).

Méthode generate_step : avance d'un seul token pour l'échantillonnage autoregressif.
"""
from __future__ import annotations

import torch
import torch.nn as nn

from peplib.tokenizer import PAD_ID, VOCAB


class CausalLSTM(nn.Module):
    """LSTM causal pour la modélisation du langage des peptides (T52).

    Parameters
    ----------
    vocab_size : int
        Taille du vocabulaire (défaut 23 = PAD + BOS + EOS + 20 AA).
    embed_dim : int
        Dimension de l'embedding.
    hidden_dim : int
        Dimension de l'état caché du LSTM.
    n_layers : int
        Nombre de couches LSTM empilées.
    dropout : float
        Taux de dropout (appliqué entre les couches LSTM et après l'embedding).
    """

    def __init__(
        self,
        vocab_size: int = len(VOCAB),
        embed_dim: int = 128,
        hidden_dim: int = 256,
        n_layers: int = 2,
        dropout: float = 0.2,
    ) -> None:
        super().__init__()
        self.vocab_size = vocab_size
        self.embed_dim  = embed_dim
        self.hidden_dim = hidden_dim
        self.n_layers   = n_layers
        self.dropout_p  = dropout

        self.embedding = nn.Embedding(vocab_size, embed_dim, padding_idx=PAD_ID)
        self.drop      = nn.Dropout(dropout)
        self.lstm      = nn.LSTM(
            input_size=embed_dim,
            hidden_size=hidden_dim,
            num_layers=n_layers,
            batch_first=True,
            dropout=dropout if n_layers > 1 else 0.0,
        )
        self.head = nn.Linear(hidden_dim, vocab_size)

        self._init_weights()

    # ------------------------------------------------------------------
    # Initialisation
    # ------------------------------------------------------------------

    def _init_weights(self) -> None:
        """Initialisation légèrement améliorée : embedding uniforme, LSTM
        orthogonal, biais à zéro, head avec init Xavier."""
        nn.init.uniform_(self.embedding.weight, -0.1, 0.1)
        if PAD_ID is not None:
            with torch.no_grad():
                self.embedding.weight[PAD_ID].fill_(0.0)

        for name, param in self.lstm.named_parameters():
            if "weight_ih" in name:
                nn.init.xavier_uniform_(param)
            elif "weight_hh" in name:
                nn.init.orthogonal_(param)
            elif "bias" in name:
                nn.init.zeros_(param)
                # Bias forget gate à 1 (technique standard pour LSTM)
                n = param.shape[0]
                param.data[n // 4 : n // 2].fill_(1.0)

        nn.init.xavier_uniform_(self.head.weight)
        nn.init.zeros_(self.head.bias)

    # ------------------------------------------------------------------
    # Forward
    # ------------------------------------------------------------------

    def forward(
        self,
        input_ids: torch.Tensor,
        hidden: tuple[torch.Tensor, torch.Tensor] | None = None,
    ) -> tuple[torch.Tensor, tuple[torch.Tensor, torch.Tensor]]:
        """
        Parameters
        ----------
        input_ids : LongTensor (B, L)
        hidden    : état caché initial (h, c) ou None (→ zéros).

        Returns
        -------
        logits : FloatTensor (B, L, vocab_size)
        hidden : tuple (h_n, c_n)  — état final, utile pour generate_step.
        """
        x = self.drop(self.embedding(input_ids))   # (B, L, embed_dim)
        out, hidden = self.lstm(x, hidden)          # (B, L, hidden_dim)
        out = self.drop(out)
        logits = self.head(out)                     # (B, L, vocab_size)
        return logits, hidden

    # ------------------------------------------------------------------
    # Génération token-à-token
    # ------------------------------------------------------------------

    def generate_step(
        self,
        token_id: torch.Tensor,
        hidden: tuple[torch.Tensor, torch.Tensor] | None = None,
    ) -> tuple[torch.Tensor, tuple[torch.Tensor, torch.Tensor]]:
        """Avance d'un seul token (pour l'échantillonneur autoregressif).

        Parameters
        ----------
        token_id : LongTensor (B, 1)  — token courant.
        hidden   : état caché courant ou None.

        Returns
        -------
        logits : FloatTensor (B, 1, vocab_size)
        hidden : nouvel état caché.
        """
        return self.forward(token_id, hidden)

    # ------------------------------------------------------------------
    # Utilitaires
    # ------------------------------------------------------------------

    def count_parameters(self) -> int:
        """Nombre de paramètres entraînables."""
        return sum(p.numel() for p in self.parameters() if p.requires_grad)

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"CausalLSTM("
            f"vocab={self.vocab_size}, "
            f"embed={self.embed_dim}, "
            f"hidden={self.hidden_dim}, "
            f"layers={self.n_layers}, "
            f"dropout={self.dropout_p}, "
            f"params={self.count_parameters():,})"
        )
