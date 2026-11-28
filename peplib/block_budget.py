"""Budget T79. Mesure les appels réels à un module PyTorch sur un lot."""
from contextlib import contextmanager
from dataclasses import dataclass


@dataclass
class BlockBudget:
    waves: int
    forward_calls: int = 0

    def result(self):
        return {"method": "blocks_approximate", "K": self.waves,
                "forward_calls_per_batch": self.forward_calls}


@contextmanager
def count_block_passes(model, waves, *, T=50):
    """Entourer l'appel T79 complet. Le hook est retiré même après exception."""
    if type(waves) is not int or not 1 <= waves <= T:
        raise ValueError("K doit être compris entre 1 et T")
    budget = BlockBudget(waves)

    def count(module, inputs):
        budget.forward_calls += 1

    hook = model.register_forward_pre_hook(count)
    try:
        yield budget
    finally:
        hook.remove()
