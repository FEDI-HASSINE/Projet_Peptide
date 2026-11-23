"""
peplib/generate.py -- Generation de peptides par echantillonnage autoregressif (T53/T57).

Protocole :
  1. Partir de <BOS>.
  2. A chaque etape, appeler le modele pour obtenir les logits du prochain token.
  3. Forcer PAD et BOS a -inf (peplib.sampling garantit cela).
  4. Echantillonner avec temperature / top_k / top_p via peplib.sampling.
  5. Arreter a <EOS> ou quand la longueur max est atteinte.
  6. Retourner uniquement les sequences valides (passe peplib.filters).

Deux points d'entree :
  - ``generate_sequences`` : liste des sequences valides (API historique).
  - ``generate_with_stats`` : memes sequences + compteurs bruts (tentatives,
    rejets par cause, arret sur epuisement du budget). Necessaire pour
    mesurer un taux de validite reel au lieu de le postuler a 100 %.

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

from peplib.filters import REJECTION_CAUSES, classify_sequence
from peplib.sampling import sample_next_token
from peplib.tokenizer import BOS_ID, EOS_ID, ID_TO_TOKEN, MAX_LEN, PAD_ID

# Budget de tentatives : 10 * n, avec un plancher de 5000.
ATTEMPT_FACTOR = 10
MIN_ATTEMPTS = 5000


def _resolve_device(model: torch.nn.Module, device: torch.device | None) -> torch.device:
    """Device explicite, sinon celui des parametres du modele."""
    if device is not None:
        return device
    try:
        return next(model.parameters()).device
    except StopIteration:
        return torch.device("cpu")


def _decode_tokens(tokens: list[int]) -> str:
    """Convertit des ids en sequence d'acides amines (BOS/EOS/PAD retires)."""
    residues: list[str] = []
    for token in tokens[1:]:  # sauter BOS
        if token in (EOS_ID, PAD_ID):
            break
        residue = ID_TO_TOKEN.get(token, "")
        if not residue or residue.startswith("<"):
            continue
        residues.append(residue)
    return "".join(residues)


def _sample_tokens(
    model: torch.nn.Module,
    is_lstm: bool,
    max_len: int,
    temperature: float,
    top_k: int | None,
    top_p: float,
    generator: torch.Generator,
    device: torch.device,
) -> list[int]:
    """Echantillonne une sequence complete de tokens, BOS inclus, EOS exclu."""
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

        next_token = sample_next_token(
            logits_t[0, 0, :],
            temperature=temperature,
            top_k=top_k if top_k and top_k > 0 else None,
            top_p=top_p,
            generator=generator,
        )

        if next_token == EOS_ID:
            break
        tokens.append(next_token)

    return tokens


def generate_with_stats(
    model: torch.nn.Module,
    n: int = 1000,
    temperature: float = 1.0,
    top_k: int | None = None,
    top_p: float = 1.0,
    max_len: int = MAX_LEN,
    seed: int | None = None,
    device: torch.device | None = None,
    max_attempts: int | None = None,
    keep_raw: bool = False,
) -> dict:
    """Genere jusqu'a n sequences valides et retourne les compteurs bruts.

    Parameters
    ----------
    model        : CausalLSTM ou CausalTransformer (passe en mode eval).
    n            : nombre de sequences valides demandees.
    temperature  : diversite de l'echantillonnage (> 0).
    top_k        : filtre top-k (None ou 0 = desactive).
    top_p        : filtre nucleus (1.0 = desactive).
    max_len      : longueur max incluant BOS et EOS (defaut : 52).
    seed         : graine du torch.Generator (reproductibilite, CPU et CUDA).
    device       : torch.device ; None = celui du modele.
    max_attempts : budget de tentatives ; None = max(10 * n, 5000).
    keep_raw     : si True, conserve toutes les sequences decodees.

    Returns
    -------
    dict avec les cles :
        sequences      : list[str] des sequences valides (au plus n) ;
        n_requested    : n ;
        n_attempts     : nombre de tentatives effectuees ;
        n_raw          : nombre de sequences decodees (= n_attempts) ;
        n_valid        : nombre de sequences valides ;
        n_invalid      : n_raw - n_valid ;
        rejections     : dict cause -> nombre (somme = n_invalid) ;
        exhausted      : True si le budget a ete epuise avant n valides ;
        max_attempts   : budget effectif ;
        raw_sequences  : toutes les sequences decodees si keep_raw, sinon None.
    """
    if n < 0:
        raise ValueError("n doit etre >= 0")
    if max_attempts is not None and max_attempts < 0:
        raise ValueError("max_attempts doit etre >= 0")

    device = _resolve_device(model, device)
    model.eval()

    generator = torch.Generator(device=device)
    if seed is not None:
        generator.manual_seed(seed)

    is_lstm = hasattr(model, "lstm")
    budget = max(ATTEMPT_FACTOR * n, MIN_ATTEMPTS) if max_attempts is None else max_attempts

    valid_sequences: list[str] = []
    raw_sequences: list[str] = []
    rejections = {cause: 0 for cause in REJECTION_CAUSES}
    attempts = 0

    with torch.no_grad():
        while len(valid_sequences) < n and attempts < budget:
            attempts += 1
            tokens = _sample_tokens(
                model=model,
                is_lstm=is_lstm,
                max_len=max_len,
                temperature=temperature,
                top_k=top_k,
                top_p=top_p,
                generator=generator,
                device=device,
            )
            sequence = _decode_tokens(tokens)
            if keep_raw:
                raw_sequences.append(sequence)

            cause = classify_sequence(sequence)
            if cause is None:
                valid_sequences.append(sequence)
            else:
                rejections[cause] = rejections.get(cause, 0) + 1

    n_valid = len(valid_sequences)
    return {
        "sequences": valid_sequences,
        "n_requested": n,
        "n_attempts": attempts,
        "n_raw": attempts,
        "n_valid": n_valid,
        "n_invalid": attempts - n_valid,
        "rejections": rejections,
        "exhausted": n_valid < n,
        "max_attempts": budget,
        "raw_sequences": raw_sequences if keep_raw else None,
    }


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
    """Genere n sequences de peptides valides par echantillonnage autoregressif.

    Seules les sequences passant les filtres de peplib.filters (alphabet
    canonique, longueur 10-50 residus) sont retournees. La boucle tourne
    jusqu'a avoir n sequences valides ou 10*n tentatives (plancher 5000).

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
    return generate_with_stats(
        model=model,
        n=n,
        temperature=temperature,
        top_k=top_k,
        top_p=top_p,
        max_len=max_len,
        seed=seed,
        device=device,
    )["sequences"]


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
