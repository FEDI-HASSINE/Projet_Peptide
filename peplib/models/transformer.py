"""
peplib/models/transformer.py — Transformer causal pour la génération de peptides (T55).

Architecture
------------
  TokenEmbedding  : Embedding(vocab_size, embed_dim, padding_idx=PAD_ID)
  PosEmbedding    : Embedding(max_len, embed_dim)  — positions apprises
  Dropout
  N × TransformerEncoderLayer(embed_dim, n_heads, ffn_dim, dropout, batch_first=True)
       utilisé comme décodeur causal via attn_mask triangulaire supérieure
  LayerNorm(embed_dim)
  Linear(embed_dim, vocab_size)

Entrée  : (B, L) LongTensor d'ids de tokens
Sortie  : (B, L, V) FloatTensor de logits

Masque causal : positions futures masquées avec -inf →
  le logit à la position i ne dépend que des tokens 0..i.

Masque de padding : positions PAD ignorées via key_padding_mask.

Aucun poids pré-entraîné — initialisation Xavier/Kaiming par défaut.
"""
from __future__ import annotations

import math

import torch
import torch.nn as nn

from peplib.tokenizer import PAD_ID, VOCAB


class CausalTransformer(nn.Module):
    """Transformer causal décodeur pour la modélisation du langage des peptides (T55).

    Parameters
    ----------
    vocab_size : int
        Taille du vocabulaire (défaut 23).
    embed_dim : int
        Dimension des embeddings (d_model).
    n_heads : int
        Nombre de têtes d'attention (embed_dim % n_heads == 0).
    n_layers : int
        Nombre de couches TransformerEncoderLayer.
    ffn_dim : int
        Dimension de la couche feed-forward interne.
    dropout : float
        Taux de dropout.
    max_len : int
        Longueur maximale de séquence pour les embeddings de position.
    """

    def __init__(
        self,
        vocab_size: int = len(VOCAB),
        embed_dim: int = 128,
        n_heads: int = 4,
        n_layers: int = 4,
        ffn_dim: int = 256,
        dropout: float = 0.1,
        max_len: int = 52,
    ) -> None:
        super().__init__()

        if embed_dim % n_heads != 0:
            raise ValueError(
                f"embed_dim ({embed_dim}) doit être divisible par n_heads ({n_heads})."
            )

        self.vocab_size = vocab_size
        self.embed_dim  = embed_dim
        self.n_heads    = n_heads
        self.n_layers   = n_layers
        self.ffn_dim    = ffn_dim
        self.dropout_p  = dropout
        self.max_len    = max_len

        # Embeddings
        self.token_emb = nn.Embedding(vocab_size, embed_dim, padding_idx=PAD_ID)
        self.pos_emb   = nn.Embedding(max_len, embed_dim)
        self.drop      = nn.Dropout(dropout)

        # Transformer layers (utilisés en mode causal via attn_mask)
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=embed_dim,
            nhead=n_heads,
            dim_feedforward=ffn_dim,
            dropout=dropout,
            batch_first=True,
            # norm_first désactivé : incompatible avec enable_nested_tensor
        )
        self.transformer = nn.TransformerEncoder(
            encoder_layer,
            num_layers=n_layers,
            norm=nn.LayerNorm(embed_dim),
        )

        self.head = nn.Linear(embed_dim, vocab_size)

        # Partage des poids token_emb ↔ head (weight tying, standard LM)
        self.head.weight = self.token_emb.weight

        self._init_weights()

    # ------------------------------------------------------------------
    # Initialisation
    # ------------------------------------------------------------------

    def _init_weights(self) -> None:
        """Initialisation Xavier pour les couches linéaires.
        
        Attention : head.weight partage le tenseur avec token_emb.weight (weight tying).
        On ne doit PAS réinitialiser head.weight, sinon on écrase token_emb.
        """
        nn.init.normal_(self.token_emb.weight, std=0.02)
        if PAD_ID is not None:
            with torch.no_grad():
                self.token_emb.weight[PAD_ID].fill_(0.0)
        nn.init.normal_(self.pos_emb.weight, std=0.02)

        for module in self.modules():
            if isinstance(module, nn.Linear):
                # head.weight partage son tenseur avec token_emb.weight
                # (weight tying) : le reinitialiser ecraserait l'embedding.
                if module is self.head:
                    if module.bias is not None:
                        nn.init.zeros_(module.bias)
                    continue
                nn.init.xavier_uniform_(module.weight)
                if module.bias is not None:
                    nn.init.zeros_(module.bias)
            elif isinstance(module, nn.LayerNorm):
                nn.init.ones_(module.weight)
                nn.init.zeros_(module.bias)

    # ------------------------------------------------------------------
    # Masque causal
    # ------------------------------------------------------------------

    @staticmethod
    def _causal_mask(seq_len: int, device: torch.device) -> torch.Tensor:
        """Retourne un masque (seq_len, seq_len) booléen True aux positions futures.

        Conventions PyTorch TransformerEncoder :
          attn_mask est additionné aux scores → on met -inf aux positions True.
          On utilise float('-inf') directement.
        """
        mask = torch.triu(
            torch.full((seq_len, seq_len), float("-inf"), device=device),
            diagonal=1,
        )
        return mask  # (L, L), -inf au-dessus de la diagonale

    # ------------------------------------------------------------------
    # Forward
    # ------------------------------------------------------------------

    def forward(self, input_ids: torch.Tensor) -> torch.Tensor:
        """
        Parameters
        ----------
        input_ids : LongTensor (B, L)

        Returns
        -------
        logits : FloatTensor (B, L, vocab_size)
        """
        B, L = input_ids.shape
        device = input_ids.device

        # Positions 0..L-1
        positions = torch.arange(L, device=device).unsqueeze(0).expand(B, -1)  # (B, L)

        # Embeddings token + position
        x = self.drop(
            self.token_emb(input_ids) + self.pos_emb(positions)
        )  # (B, L, embed_dim)

        # Masque causal
        causal_mask = self._causal_mask(L, device)  # (L, L)

        # Masque de padding (True = ignorer) — converti en float pour matcher attn_mask
        pad_mask = (input_ids == PAD_ID).float()  # (B, L) — 0.0 valide, 1.0 ignoré

        # Passage dans le Transformer
        out = self.transformer(x, mask=causal_mask, src_key_padding_mask=pad_mask)
        # (B, L, embed_dim)

        logits = self.head(out)  # (B, L, vocab_size)
        return logits

    # ------------------------------------------------------------------
    # Utilitaires
    # ------------------------------------------------------------------

    def count_parameters(self) -> int:
        """Nombre de paramètres entraînables (weight-tying déduit une fois)."""
        return sum(p.numel() for p in self.parameters() if p.requires_grad)

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"CausalTransformer("
            f"vocab={self.vocab_size}, "
            f"embed={self.embed_dim}, "
            f"heads={self.n_heads}, "
            f"layers={self.n_layers}, "
            f"ffn={self.ffn_dim}, "
            f"dropout={self.dropout_p}, "
            f"params={self.count_parameters():,})"
        )
