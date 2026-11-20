"""
peplib/generate.py -- Generation de peptides par echantillonnage autoregressif (T53/T57).

Protocole :
  1. Partir de <BOS>.
  2. A chaque etape, appeler le modele pour obtenir les logits du prochain token.
  3. Forcer PAD et BOS a -inf (peplib.sampling garantit cela).
  4. Echantillonner avec temperature / top_k / top_p via peplib.sampling.
  5. Arreter a <EOS> ou quand la longueur max est atteinte.
  6. Retourner uniquement les sequences valides (passe peplib.filters).

Usage
-----
    from peplib.generate import generate_sequences

    seqs = generate_sequences(
        model, n=1000,
        temperature=1.0, top_k=0, top_p=1.0,
        seed=42,
    )
"""
from __future__ import annotations

import torch

from peplib.tokenizer import BOS_ID, EOS_ID, PAD_ID, ID_TO_TOKEN, MAX_LEN
from peplib.sampling import sample_next_token
from peplib.filters import classify_sequence


def generate_sequences(
    model: torch.nn.Module,
    n: int = 1000,
    temperature: float = 1.0,
    top_k: int | None = None,
    top_p: float = 1.0,
    max_len: int = MAX_LEN,
    seed: int | None = None,
    device: torch.device | None = None,
) -> list[str]:
    """Genere n sequences de peptides par echantillonnage autoregressif.

    Seules les sequences passant les filtres de peplib.filters (alphabet
    canonique, longueur 10-50 residus) sont retournees. La boucle tourne
    jusqu'a avoir n sequences valides ou 10*n tentatives.

    Parameters
    ----------
    model       : CausalLSTM ou CausalTransformer en mode eval.
    n           : nombre de sequences valides a produire.
    temperature : diversite de l'echantillonnage (> 0).
    top_k       : filtre top-k (None = desactive).
    top_p       : filtre nucleus (1.0 = desactive).
    max_len     : longueur max incluant BOS et EOS (defaut : 52).
    seed        : graine pour torch.Generator (reproductibilite).
    device      : torch.device ; None = utilise le meme que le modele.

    Returns
    -------
    Liste de chaines d'acides amines valides (sans BOS/EOS/PAD).
    """
    if device is None:
        try:
            device = next(model.parameters()).device
        except StopIteration:
            device = torch.device("cpu")

    model.eval()

    # Generateur seede (reproductibilite)
    gen = torch.Generator(device=device)
    if seed is not None:
        gen.manual_seed(seed)

    # Detecter le type de modele (LSTM garde un etat cache)
    is_lstm = hasattr(model, "lstm")

    valid_seqs: list[str] = []
    max_attempts = max(10 * n, 5000)
    attempts = 0

    with torch.no_grad():
        while len(valid_seqs) < n and attempts < max_attempts:
            attempts += 1
            tokens = [BOS_ID]
            hidden = None  # etat LSTM initial

            for _ in range(max_len - 2):  # max_len - BOS - EOS = nb residus max
                input_t = torch.tensor(
                    [[tokens[-1]]], dtype=torch.long, device=device
                )  # (1, 1)

                if is_lstm:
                    logits_t, hidden = model(input_t, hidden)
                else:
                    # Transformer : on repasse toute la sequence depuis BOS
                    input_full = torch.tensor(
                        [tokens], dtype=torch.long, device=device
                    )  # (1, L)
                    logits_full = model(input_full)  # (1, L, V)
                    logits_t = logits_full[:, -1:, :]  # (1, 1, V)

                logits_1d = logits_t[0, 0, :]  # (V,)

                next_tok = sample_next_token(
                    logits_1d,
                    temperature=temperature,
                    top_k=top_k if top_k and top_k > 0 else None,
                    top_p=top_p,
                    generator=gen if device.type == "cpu" else None,
                )

                if next_tok == EOS_ID:
                    break
                tokens.append(next_tok)

            # Decoder la sequence (exclure BOS et tout ce qui suit EOS)
            residues = []
            for tok in tokens[1:]:  # sauter BOS
                if tok in (EOS_ID, PAD_ID):
                    break
                aa = ID_TO_TOKEN.get(tok, "")
                if not aa or aa.startswith("<"):
                    continue
                residues.append(aa)

            seq = "".join(residues)

            # Valider avec peplib.filters (alphabet canonique, 10-50 aa)
            if seq and classify_sequence(seq) is None:
                valid_seqs.append(seq)

    return valid_seqs


def generate_all_configs(
    model: torch.nn.Module,
    sampling_grid: dict,
    n_per_config: int = 1000,
    seed: int = 42,
    device: torch.device | None = None,
) -> dict[str, list[str]]:
    """Genere des sequences pour chaque combinaison de la grille d'echantillonnage.

    Parameters
    ----------
    model         : modele entraine.
    sampling_grid : dict avec cles 'temperature', 'top_k', 'top_p'.
    n_per_config  : sequences a generer par config.
    seed          : graine de base (incrementee par config pour diversite).

    Returns
    -------
    dict keye par 'T{temp}_k{top_k}_p{top_p}' -> liste de sequences.
    """
    temperatures = sampling_grid.get("temperature", [1.0])
    top_ks       = sampling_grid.get("top_k", [0])
    top_ps       = sampling_grid.get("top_p", [1.0])

    results: dict[str, list[str]] = {}
    config_seed = seed

    for temp in temperatures:
        for k in top_ks:
            for p in top_ps:
                key = f"T{temp}_k{k}_p{p}"
                seqs = generate_sequences(
                    model,
                    n=n_per_config,
                    temperature=temp,
                    top_k=k if k > 0 else None,
                    top_p=p,
                    seed=config_seed,
                    device=device,
                )
                results[key] = seqs
                config_seed += 1
                print(
                    f"[generate] {key} : {len(seqs)}/{n_per_config} sequences valides"
                )

    return results
