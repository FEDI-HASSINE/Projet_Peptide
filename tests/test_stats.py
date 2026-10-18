import math
import pandas as pd
import pytest
from peplib.stats import (
    AA_ORDER,
    compute_composition,
    compute_entropy,
    compute_js_matrix,
    jensen_shannon_divergence,
    kl_divergence,
    plot_length_histograms,
)


def test_composition_sum_to_one():
    comp = compute_composition(["ACDEFGHIKL", "MNPQRSTVWY"])
    assert abs(sum(comp.values()) - 1.0) < 1e-12
    assert all(abs(v - 0.05) < 1e-12 for v in comp.values())


def test_composition_empty_raises():
    with pytest.raises(ValueError):
        compute_composition([])


def test_jensen_shannon_properties():
    comp_a = {aa: (1.0 if aa == "A" else 0.0) for aa in AA_ORDER}
    comp_c = {aa: (1.0 if aa == "C" else 0.0) for aa in AA_ORDER}
    assert abs(jensen_shannon_divergence(comp_a, comp_a)) < 1e-12
    assert (
        abs(jensen_shannon_divergence(comp_a, comp_c) - math.log(2.0)) < 1e-12
    )


def test_js_matrix_and_entropy(tmp_path):
    c1 = compute_composition(["AAAAAAAAAA"])
    c2 = compute_composition(["CCCCCCCCCC"])
    mat = compute_js_matrix({"s1": c1, "s2": c2})
    assert mat.shape == (2, 2)
    assert abs(mat.loc["s1", "s1"]) < 1e-12
    assert abs(mat.loc["s1", "s2"] - math.log(2.0)) < 1e-12
    assert abs(compute_entropy(c1)) < 1e-12
    df_test = pd.DataFrame({"sequence": ["ACDEFGHIKL", "ACDEFGHIKLMN"]})
    out_png = plot_length_histograms({"s1": df_test}, out_path=str(tmp_path / "hist.png"))
    assert (tmp_path / "hist.png").exists()


def test_kl_divergence_valeur():
    import numpy as np

    p = np.array([0.5, 0.5])
    q = np.array([0.9, 0.1])
    attendu = 0.5 * math.log(0.5 / 0.9) + 0.5 * math.log(0.5 / 0.1)
    assert kl_divergence(p, q) == pytest.approx(attendu)


def test_kl_divergence_asymetrie():
    import numpy as np

    p = np.array([0.5, 0.5])
    q = np.array([0.9, 0.1])
    assert kl_divergence(p, q) != kl_divergence(q, p)
    assert kl_divergence(p, p) == 0.0


def test_kl_q_zero_donne_inf():
    import numpy as np

    p = np.array([0.5, 0.5])
    q = np.array([1.0, 0.0])
    assert kl_divergence(p, q) == float("inf")


def test_entropy_ancre():
    c = compute_composition(["AAAAAAAAAA"])
    assert abs(compute_entropy(c)) < 1e-12
    c2 = compute_composition(["ACDEFGHIKL", "MNPQRSTVWY"])
    assert compute_entropy(c2) >= 0.0
    c3 = compute_composition(["A" * 10, "C" * 10])
    assert compute_entropy(c3) >= 0.0