import math

import numpy as np
import pandas as pd
import pytest

from peplib.pwm import (
    LOG2_S,
    LN_2,
    S_ALPHABET,
    build_counts,
    pwm_matrix,
    position_entropy,
    position_information,
    small_sample_correction,
)


def test_pwm_single_sequence_with_background_quarter():
    """Ancre 1 : PWM(0, A) = log2(1 / 0.25) = 2.0 avec un fond A = 0.25."""
    pwm = pwm_matrix(["AAAA"], terminus="N", background={"A": 0.25})
    assert pwm.iloc[0]["A"] == pytest.approx(2.0)


def test_entropy_deterministic_position_is_zero():
    """Ancre 2 : position déterministe (100% A) -> H = 0.0."""
    h = position_entropy(["AAAAAAAAAA"], terminus="N")
    assert h.iloc[0] == pytest.approx(0.0)


def test_information_deterministic_n10():
    """Ancre 3 : I = log2(20) - 0 - 19/(2·ln2·10) = 2.951349..."""
    seqs = ["AAAAAAAAAA"] * 10
    info = position_information(seqs, terminus="N")
    attendu = math.log2(20) - (S_ALPHABET - 1) / (2 * math.log(2) * 10)
    assert info.iloc[0] == pytest.approx(attendu)


def test_zero_masking_no_inf():
    """Ancre 4 : composition 100% A -> H = 0.0, pas d'erreur, pas de -inf."""
    h = position_entropy(["AAAA", "AAAA"], terminus="N")
    assert h.iloc[0] == pytest.approx(0.0)
    pwm = pwm_matrix(["AAAA", "AAAA"], terminus="N")
    col = pwm.iloc[0]
    assert not np.isinf(col.to_numpy()).any()  # pas de -inf
    assert np.isnan(col["C"])  # entrée nulle masquée


def test_no_truncation_n_per_position():
    """Ancre 5 : ['AAAAAAAAAA', 'AA'] (N) -> n = [2,2,1,1,1,1,1,1,1,1]."""
    _, n = build_counts(["AAAAAAAAAA", "AA"], terminus="N")
    attendu = np.array([2, 2, 1, 1, 1, 1, 1, 1, 1, 1])
    assert n.tolist() == attendu.tolist()


def test_c_terminus_alignment():
    """Ancre 6 : ['AA', 'AAA'] (C) -> n = [2, 2, 1] depuis le dernier résidu."""
    _, n = build_counts(["AA", "AAA"], terminus="C")
    assert n.tolist() == [2, 2, 1]


def test_small_sample_correction_values():
    """Ancre 7 : e(1) = 19/(2·ln2) = 13.7057... ; e(20) = 19/(40·ln2) = 0.68529..."""
    assert small_sample_correction(1) == pytest.approx(
        (S_ALPHABET - 1) / (2 * math.log(2))
    )
    assert small_sample_correction(20) == pytest.approx(
        (S_ALPHABET - 1) / (40 * math.log(2))
    )
    # valeurs décimales commentées (calcul main)
    # e(1) = 19 / 1.3862943611198906 = 13.7056028884...
    # e(20) = 0.6852801444...
    assert small_sample_correction(np.array([1, 20])) == pytest.approx(
        np.array([13.7056028884, 0.6852801444]), rel=1e-9
    )


def test_uniform_position_entropy_is_log2_20():
    """Position uniforme sur 20 AA : H = log2(20)."""
    seqs = [aa * 2 for aa in "ACDEFGHIKLMNPQRSTVWY"]  # position 0 couvre les 20
    h = position_entropy(seqs, terminus="N")
    assert h.iloc[0] == pytest.approx(LOG2_S)


def test_background_default_is_global_composition():
    """Fond par défaut = composition globale du corpus fourni."""
    pwm = pwm_matrix(["AAAC"], terminus="N")
    # p_bg(A) = 3/4 -> PWM(0, A) = log2(1 / (3/4)) = log2(4/3)
    assert pwm.iloc[0]["A"] == pytest.approx(math.log2(4 / 3))


def test_rejects_non_alphabet_and_bad_terminus():
    with pytest.raises(ValueError):
        build_counts(["ACDX"], terminus="N")
    with pytest.raises(ValueError):
        build_counts(["ACDE"], terminus="X")


def test_module_is_real_peplib():
    """Les tests portent sur le vrai module, pas sur une copie."""
    import peplib.pwm as m

    assert build_counts.__module__ == "peplib.pwm"
    assert m.LN_2 == LN_2
