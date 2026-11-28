import pandas as pd
import pytest
import torch
import yaml

from peplib.audit_week4 import memorization_report
from peplib.block_budget import count_block_passes
from peplib.traceability import checkpoint_manifest, config_sha256
from scripts.check_week4_artifacts import verify_checkpoints, verify_runs


def test_audit_aggregate_only():
    flags = pd.DataFrame({"sequence": ["A", "C", "C", "D"],
                          "toxicity_status": ["not_flagged_by_uniprot", "not_flagged_by_uniprot",
                                              "flagged", "unknown"]})
    records = pd.DataFrame({"sequence": ["A", "C", "D", "E"], "toxin_flag": False})
    result = memorization_report(records, ["A"], flags)
    assert (result["n_eligible"], result["n_excluded"], result["fraction_ge_90"]) == (1, 3, 1)
    assert result["identity_deciles"] == [0] * 10 + [1]
    assert result["wilson_95"][0] < 1
    assert memorization_report(records.iloc[1:], ["A"], flags)["fraction_ge_90"] is None
    assert memorization_report(records.iloc[:1], ["CA"], flags)["identity_deciles"][-1] == 1
    records.loc[1] = ["A", True]
    assert memorization_report(records, ["A"], flags)["n_eligible"] == 0


@pytest.mark.parametrize("k", [5, 10, 25, 50])
def test_measured_calls(k):
    toy = torch.nn.Identity()
    with count_block_passes(toy, k) as budget:
        for _ in range(3):
            toy(torch.tensor([100]))
    toy(torch.tensor([100]))
    assert budget.result() == {"method": "blocks_approximate", "K": k, "forward_calls_per_batch": 3}


def test_hashes(tmp_path):
    config = {"toy": 1}
    config_path = tmp_path / "toy.yaml"
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    runs = pd.DataFrame({"seed": [42, 43, 44], "epoch": 1, "T": 50, "schedule": "linear",
                         "nelbo_val": 1, "n_params": 2, "elapsed_s": 1, "git_commit": "toy",
                         "config_sha256": config_sha256(config)})
    verify_runs(runs, config_path)
    runs.loc[0, "config_sha256"] = "wrong"
    with pytest.raises(ValueError):
        verify_runs(runs, config_path)
    checkpoint = tmp_path / "toy.pt"
    checkpoint.write_bytes(b"toy")
    manifest = tmp_path / "manifest.csv"
    checkpoint_manifest(tmp_path).to_csv(manifest, index=False)
    verify_checkpoints(manifest, tmp_path)
    checkpoint.write_bytes(b"bad")
    with pytest.raises(AssertionError):
        verify_checkpoints(manifest, tmp_path)
