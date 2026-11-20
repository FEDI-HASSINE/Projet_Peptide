"""peplib/models — Modèles de génération de peptides (Semaine 3)."""
from peplib.models.lstm import CausalLSTM
from peplib.models.transformer import CausalTransformer

__all__ = ["CausalLSTM", "CausalTransformer"]
