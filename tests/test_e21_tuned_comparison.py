from eval.e21_tuned_comparison import rank_results


def test_model_must_improve_s1_clearly_and_preserve_external_auc() -> None:
    rows = [
        {"model": "golden", "S1_auc": 0.9, "KM_auc": 0.86},
        {"model": "internal_only", "S1_auc": 0.92, "KM_auc": 0.85},
        {"model": "small_gain", "S1_auc": 0.904, "KM_auc": 0.89},
        {"model": "qualifies", "S1_auc": 0.91, "KM_auc": 0.86},
    ]

    ranked, qualifying = rank_results(rows, golden_name="golden")
    by_model = ranked.set_index("model")

    assert qualifying == ["qualifies"]
    assert list(ranked["model"]) == ["internal_only", "qualifies", "small_gain", "golden"]
    assert by_model.loc["small_gain", "KM_rank"] == 1
    assert by_model.loc["golden", "KM_rank"] == by_model.loc["qualifies", "KM_rank"] == 2
