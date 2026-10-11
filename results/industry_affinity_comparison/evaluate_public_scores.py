"""Evaluate frozen AeroFold test predictions against pinned public AbBiBench scores.

Run from any directory with the project's Python environment. No model fitting,
pipeline changes, downloads, or target-label selection are performed here.
"""

from __future__ import annotations

import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
COMMIT = "7d9c73e23c3c535fbb73ef073cf2fb9633e1fda8"
TASKS = {
    "cr9114_full_h1_binding": ("4fqi_h1", "CR9114 / H1"),
    "cr9114_full_h3_binding": ("4fqi_h3", "CR9114 / H3"),
    "cr6261_full_h1_binding": ("3gbn_h1", "CR6261 / H1"),
    "cr6261_full_h9_binding": ("3gbn_h9", "CR6261 / H9"),
}
# Score directions are specified before evaluating the test sets. Likelihood is
# higher-better; physical binding energy is lower-better (use minus dg).
MODELS = [
    ("ESM-IF1", "ESMIF1", "log-likelihood", 1),
    ("ProteinMPNN", "ProteinMPNN", "log-likelihood", 1),
    ("AntiFold", "Antifold", "log-likelihood", 1),
    ("SaProt", "SaProt", "log-likelihood", 1),
    ("ESM-2", "ESM2", "log-likelihood", 1),
    ("FoldX", "FoldX", "dg", -1),
    ("DiffAb", "diffab", "log-likelihood", 1),
    ("DiffAb fixed backbone", "diffab", "log-likelihood (fixed backbone)", 1),
]


def read_csv(path):
    return pd.read_csv(path, float_precision="round_trip")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def protected_files_check():
    before = json.loads((HERE / "pipeline_hashes_before.json").read_text())
    changed = [name for name, value in before.items() if sha(ROOT / name) != value]
    assert not changed, f"Protected pipeline files changed: {changed}"
    return len(before)


def metrics(y, prediction, latent):
    y, prediction, latent = map(np.asarray, (y, prediction, latent))
    assert np.isfinite(y).all() and np.isfinite(prediction).all() and np.isfinite(latent).all()
    n, k = len(y), min(24, len(y))
    cutoff = np.sort(y)[-max(1, math.ceil(0.05 * n))]
    positive = y >= cutoff
    selected = np.argsort(-prediction, kind="stable")[:k]
    hits = int(positive[selected].sum())
    count = int(positive.sum())

    def rho(score):
        return float(spearmanr(y, score).statistic) if np.ptp(y) and np.ptp(score) else float("nan")

    result = {
        "spearman": rho(prediction),
        "latent_spearman": rho(latent),
        "top_recall": hits / count,
        "top_enrichment": (hits / k) / (count / n),
        "top_hits": hits,
        "top_population": count,
        "top_k": k,
        "test_count": n,
        "random_recall": k / n,
        "recall_ceiling": min(k, count) / count,
        "top5_threshold": float(cutoff),
        "top5_tie_expanded": count > math.ceil(0.05 * n),
    }
    # Independent calculation: scalar sorted indices and percentile order,
    # Pearson on average ranks, and direct confusion-count arithmetic.
    cutoff2 = sorted(map(float, y), reverse=True)[math.ceil(n / 20) - 1]
    positives2 = {i for i, value in enumerate(y) if value >= cutoff2}
    selected2 = sorted(range(n), key=lambda i: (-float(prediction[i]), i))[:k]
    hits2 = sum(i in positives2 for i in selected2)
    assert hits == hits2 and count == len(positives2) and cutoff == cutoff2
    assert abs(result["top_recall"] - hits2 / len(positives2)) < 1e-14
    assert abs(result["top_enrichment"] - hits2 * n / (k * len(positives2))) < 1e-14
    for column, score in [("spearman", prediction), ("latent_spearman", latent)]:
        check = float(np.corrcoef(rankdata(y), rankdata(score))[0, 1]) if np.ptp(y) and np.ptp(score) else float("nan")
        assert np.isclose(result[column], check, atol=1e-12, equal_nan=True)
    return result


def sequence_matrix(series):
    assert series.notna().all() and series.str.len().nunique() == 1
    return np.array([list(s) for s in series], dtype="U1")


def signatures(matrix, positions):
    return pd.Index(["".join(row) for row in matrix[:, positions]])


def expected_readme_values():
    lines = (HERE / "inputs/AbBiBench_README.md").read_text(encoding="utf-8").splitlines()
    header = next(line for line in lines if line.startswith("| Rank"))
    columns = [s.strip() for s in header.strip("|").split("|")]
    values = {}
    for line in lines:
        if not line.startswith("|"):
            continue
        items = [s.strip() for s in line.strip("|").split("|")]
        if len(items) != len(columns):
            continue
        row = dict(zip(columns, items))
        try:
            values[row["Model"]] = {assay: float(row[assay]) for assay, _ in TASKS.values()}
        except (ValueError, KeyError):
            pass
    return values


def main():
    protected_count = protected_files_check()
    manifest = json.loads((HERE / "source_manifest.json").read_text())
    for entry in manifest["files"]:
        path = ROOT / entry["file"]
        assert path.stat().st_size == entry["bytes"] and sha(path) == entry["sha256"]
    frozen = read_csv(ROOT / "results/multifamily/heldout_predictions.csv")
    frozen = frozen[frozen.model.eq("AeroFold")].copy()
    old_runs = read_csv(ROOT / "results/multifamily/benchmark_runs.csv")
    old_runs = old_runs[old_runs.model.eq("AeroFold")].set_index(["task", "split", "seed"])
    protocol = json.loads((ROOT / "results/multifamily/evaluation_protocol.json").read_text())
    assert protocol["top_k"] == 24 and protocol["seeds"] == [314159, 271828]
    expected = expected_readme_values()
    run_rows, predictions, mapping_audit, published_checks = [], [], [], []

    for task, (assay, display) in TASKS.items():
        local = read_csv(ROOT / f"data/{task}.csv")
        local_matrix = sequence_matrix(local.heavy_sequence)
        positions = np.flatnonzero((local_matrix != local_matrix[0]).any(axis=0))
        local_sig = signatures(local_matrix, positions)
        assert local_sig.is_unique and local.sequence_sha256.is_unique and local.heavy_sequence.is_unique
        local = local.assign(genotype=local_sig).set_index("genotype")
        reference = read_csv(HERE / f"inputs/{assay}_benchmarking_data_ESMIF1_scores.csv")
        reference_matrix = sequence_matrix(reference.mut_heavy_chain_seq)
        assert reference_matrix.shape[1] == local_matrix.shape[1]
        assert np.array_equal(positions, np.flatnonzero((reference_matrix != reference_matrix[0]).any(axis=0)))
        reference_sig = signatures(reference_matrix, positions)
        assert reference_sig.is_unique and set(reference_sig) == set(local.index)
        reference = reference.assign(genotype=reference_sig).set_index("genotype").loc[local.index]
        constant = [i for i in range(local_matrix.shape[1]) if i not in positions]
        fixed_differences = [
            {"position_1_based": i + 1, "local_residue": str(local_matrix[0, i]), "public_residue": str(reference_matrix[0, i])}
            for i in constant if local_matrix[0, i] != reference_matrix[0, i]
        ]
        assert fixed_differences == ([] if task.startswith("cr9114") else [
            {"position_1_based": 1, "local_residue": "Q", "public_residue": "E"},
            {"position_1_based": 113, "local_residue": "Q", "public_residue": "K"},
        ])
        ref_label_error = float(np.abs(local.score.to_numpy() - reference.binding_score.to_numpy()).max())
        assert ref_label_error < 1e-10
        source_scores = {}
        for model, filename_model, column, sign in MODELS:
            path = HERE / f"inputs/{assay}_benchmarking_data_{filename_model}_scores.csv"
            public = read_csv(path)
            original_count = len(public)
            missing_labels = int(public.binding_score.isna().sum())
            # Public releases include extra unmeasured genotypes. They are not
            # members of AeroFold's measured library and are excluded explicitly.
            public = public.dropna(subset=["binding_score"]).copy()
            pmatrix = sequence_matrix(public.mut_heavy_chain_seq)
            assert pmatrix.shape[1] == local_matrix.shape[1]
            assert (pmatrix[:, constant] == reference_matrix[0, constant]).all()
            p_sig = signatures(pmatrix, positions)
            assert p_sig.is_unique and set(p_sig) == set(local.index)
            public = public.assign(genotype=p_sig).set_index("genotype").loc[local.index]
            # Sequence and genotype identity, never experimental labels, select
            # the join. Labels are used only to validate the completed join.
            assert public.mut_heavy_chain_seq.equals(reference.mut_heavy_chain_seq)
            label_error = float(np.abs(local.score.to_numpy() - public.binding_score.to_numpy()).max())
            assert label_error < 1e-10
            score = public[column] * sign
            assert np.isfinite(score).all()
            source_scores[model] = score
            full_rho = float(spearmanr(local.score, score).statistic)
            readme_name = {"ESM-IF1": "ESMIF1", "AntiFold": "Antifold", "ESM-2": "ESM2", "DiffAb": "diffab", "DiffAb fixed backbone": "diffab_fixbb"}.get(model, model)
            reported = expected[readme_name][assay]
            assert abs(full_rho - reported) <= 0.00501, (task, model, full_rho, reported)
            published_checks.append({"task": task, "model": model, "recomputed_whole_library_spearman": full_rho, "published_rounded_spearman": reported})
            mapping_audit.append({
                "task": task, "model": model, "source_file": path.name,
                "public_rows": original_count, "unmeasured_public_rows_excluded": missing_labels,
                "measured_variants": len(local), "matched_variants": len(public),
                "variable_positions_1_based": (positions + 1).tolist(),
                "fixed_background_differences": fixed_differences,
                "maximum_label_difference": label_error,
                "public_score_column": column, "score_direction_multiplier": sign,
            })
        local_by_key = local.reset_index().set_index("sequence_sha256")
        tests = frozen[frozen.task.eq(task)]
        for (split, seed), test in tests.groupby(["split", "seed"], sort=False):
            assert test.key.is_unique and len(test) <= protocol["test_limit"]
            target = local_by_key.loc[test.key.to_list()]
            assert np.array_equal(target.score.to_numpy(), test.truth.to_numpy())
            genotype = target.genotype.to_list()
            y = test.truth.to_numpy()
            eval_scores = {"AeroFold": (test.rank_score.to_numpy(), test.latent_rank_score.to_numpy())}
            eval_scores.update({model: (score.loc[genotype].to_numpy(), score.loc[genotype].to_numpy()) for model, score in source_scores.items()})
            for model, (score, continuous) in eval_scores.items():
                result = metrics(y, score, continuous)
                if model == "AeroFold":
                    previous = old_runs.loc[(task, split, seed)]
                    for field in ["spearman", "latent_spearman", "top_recall", "top_enrichment", "top_hits", "top_population", "test_count"]:
                        assert np.isclose(result[field], previous[field], atol=1e-12, equal_nan=True), (task, split, seed, field)
                run_rows.append({"task": task, "assay": display, "split": split, "seed": int(seed), "model": model, "target_supervision": model == "AeroFold", "training_labels": int(old_runs.loc[(task, split, seed), "training_labels"]) if model == "AeroFold" else 0, **result})
                for order, (key, variant, true, pred, cont, gt) in enumerate(zip(test.key, test.variant_id, y, score, continuous, genotype)):
                    predictions.append({"task": task, "split": split, "seed": int(seed), "model": model, "test_order": order, "key": key, "variant_id": variant, "genotype_signature": gt, "truth": float(true), "rank_score": float(pred), "continuous_score": float(cont)})

    runs = pd.DataFrame(run_rows)
    assert len(runs) == len(TASKS) * 2 * 2 * (len(MODELS) + 1)
    aggregate_fields = ["spearman", "latent_spearman", "top_recall", "top_enrichment", "top_hits", "top_population", "test_count", "random_recall", "recall_ceiling"]
    summary = runs.groupby(["task", "assay", "split", "model"], sort=False)[aggregate_fields].agg(["mean", "min", "max"])
    summary.columns = ["_".join(c) for c in summary.columns]
    summary = summary.reset_index()
    winners = []
    for (task, split), group in summary.groupby(["task", "split"], sort=False):
        own = group[group.model.eq("AeroFold")].iloc[0]
        external = group[~group.model.eq("AeroFold")]
        rho_best = external.loc[external.spearman_mean.idxmax()]
        recall_max = external.top_recall_mean.max()
        discovery_best = external[np.isclose(external.top_recall_mean, recall_max)]
        enrichment_max = external.top_enrichment_mean.max()
        winners.append({
            "task": task, "assay": own.assay, "split": split,
            "aerofold_spearman": own.spearman_mean,
            "aerofold_continuous_spearman": own.latent_spearman_mean,
            "best_public_spearman": rho_best.spearman_mean,
            "best_public_spearman_model": rho_best.model,
            "aerofold_recall": own.top_recall_mean,
            "best_public_recall": recall_max,
            "best_public_recall_models": "; ".join(discovery_best.model),
            "aerofold_enrichment": own.top_enrichment_mean,
            "best_public_enrichment": enrichment_max,
            "aerofold_continuous_spearman_higher": bool(own.latent_spearman_mean > rho_best.spearman_mean),
            "aerofold_discovery_higher": bool(own.top_recall_mean > recall_max + 1e-12),
            "aerofold_discovery_equal": bool(np.isclose(own.top_recall_mean, recall_max)),
        })
    winner_frame = pd.DataFrame(winners)
    runs.to_csv(HERE / "matched_test_runs.csv", index=False)
    summary.to_csv(HERE / "matched_test_summary.csv", index=False)
    winner_frame.to_csv(HERE / "comparison_overview.csv", index=False)
    pd.DataFrame(predictions).to_csv(HERE / "matched_test_predictions.csv", index=False)
    pd.DataFrame(published_checks).to_csv(HERE / "published_correlation_checks.csv", index=False)
    pd.DataFrame(mapping_audit).to_json(HERE / "variant_mapping_audit.json", orient="records", indent=2)
    coverage = []
    for task in protocol["tasks"]:
        coverage.append({"task": task, "matched_public_scores_available": task in TASKS, "reason": "Pinned AbBiBench scores with sequence identifiers" if task in TASKS else "No compatible per-variant score release located in the checked AbBiBench and FLAb releases"})
    pd.DataFrame(coverage).to_csv(HERE / "assay_coverage.csv", index=False)
    assert protected_files_check() == protected_count
    audit = {
        "status": "passed", "evaluated_utc": datetime.now(timezone.utc).isoformat(),
        "public_source_commit": COMMIT, "pipeline_version": protocol["version"],
        "pipeline_code_sha256": protocol["code_sha256"],
        "protected_files_unchanged": protected_count,
        "verified_source_file_hashes": len(manifest["files"]),
        "matched_assays": len(TASKS), "external_comparators": len(MODELS),
        "assay_split_seed_model_rows": len(runs), "paired_prediction_rows": len(predictions),
        "independent_metric_calculation": "All runs passed: average-rank Pearson and scalar top-k confusion-count arithmetic",
        "frozen_aerofold_metrics_reproduced": True,
        "published_whole_library_correlations_reproduced": True,
        "test_identity_and_label_validation": "Unique genotype joins; full measured-library coverage; experimental labels agree to 1e-10",
        "no_model_training_or_selection_performed": True,
        "cr6261_framework_difference": "Local Q1/Q113; public E1/K113 (one-based heavy-chain positions)",
        "excluded_public_models": {"Boltz-2": "No sequence or canonical variant identifiers in downloaded confidence-score files", "AlphaFold 3": "Ordinal mutant identifiers; no verified genotype crosswalk"},
        "important_limits": [
            "AeroFold uses target-assay supervision; published comparators use pretrained likelihood or energy proxies without target-label fitting in this evaluation.",
            "CR6261 genotypes and experimental labels match, but two fixed framework residues differ between model inputs.",
            "Two seeds only; means and ranges are descriptive, not confidence intervals or proof of general industry superiority.",
            "Public libraries and previously inspected retrospective test partitions; no independent prospective affinity validation.",
            "Pretraining exposure cannot be audited from score CSVs. Light-chain inputs are not recorded in these public CSVs.",
            "Best public comparator is a post-hoc per-metric reference, not a selected or trained ensemble.",
        ],
        "output_sha256": {p.name: sha(p) for p in HERE.glob("*.csv")},
    }
    (HERE / "validation.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    write_report()
    print(winner_frame.to_string(index=False))
    print(json.dumps({key: audit[key] for key in ["status", "matched_assays", "external_comparators", "assay_split_seed_model_rows", "paired_prediction_rows", "protected_files_unchanged"]}))


def write_report():
    """Build the human-readable comparison from the validated saved outputs."""
    overview = read_csv(HERE / "comparison_overview.csv")
    summary = read_csv(HERE / "matched_test_summary.csv")
    audit = json.loads((HERE / "validation.json").read_text())
    original = read_csv(ROOT / "results/multifamily/benchmark_summary.csv")
    original = original[original.model.eq("AeroFold")]
    assert audit["status"] == "passed"
    for name, value in audit["output_sha256"].items():
        assert sha(HERE / name) == value
    protected_files_check()
    rho_wins = int(overview.aerofold_continuous_spearman_higher.sum())
    discovery_wins = int(overview.aerofold_discovery_higher.sum())
    discovery_ties = int(overview.aerofold_discovery_equal.sum())
    release = f"https://github.com/MSBMI-SAFE/AbBiBench/tree/{COMMIT}/notebooks/scoring_outputs"
    leaderboard = f"https://github.com/MSBMI-SAFE/AbBiBench/blob/{COMMIT}/README.md"
    split_names = {"random": "Interpolation", "mutation_extrapolation": "Higher-mutation extrapolation"}

    def table(header, rows):
        return "\n".join(["| " + " | ".join(header) + " |", "| " + " | ".join(["---"] * len(header)) + " |", *["| " + " | ".join(map(str, row)) + " |" for row in rows]])

    def file_link(name):
        return f"[{name}](<{(HERE / name).as_posix()}>)"

    rho_rows, discovery_rows = [], []
    for row in overview.itertuples():
        split = split_names[row.split]
        rho_rows.append([row.assay, split, f"{row.aerofold_spearman:.3f}", f"{row.aerofold_continuous_spearman:.3f}", f"{row.best_public_spearman:.3f}", row.best_public_spearman_model])
        discovery_rows.append([row.assay, split, f"{row.aerofold_recall:.1%}", f"{row.best_public_recall:.1%}", f"{row.aerofold_enrichment:.2f}x", f"{row.best_public_enrichment:.2f}x", row.best_public_recall_models.replace(";", ",")])

    lines = [
        "# AeroFold affinity comparison with published models", "",
        f"AeroFold's continuous Spearman exceeds the best of eight public comparators on **{rho_wins} of {len(overview)} assay/split combinations**. Its top-5% recall is higher on **{discovery_wins}**, tied on **{discovery_ties}**, and lower on **{len(overview) - discovery_wins - discovery_ties}**. These counts describe this release and these retrospective partitions; they are not evidence of industry-wide superiority.", "",
        "The evaluation reuses the frozen AeroFold GPU-trained model's saved predictions. No fitting, tuning, retraining, or changes to the working notebook were performed. All models are scored on the same experimental test genotypes within each assay/split/seed. CR9114 sequences match exactly. CR6261 has a disclosed fixed-background mismatch, described below.", "",
        "Eight comparators: **ESM-IF1, ProteinMPNN, AntiFold, SaProt, ESM-2, FoldX, DiffAb, and DiffAb with a fixed backbone**. Their per-variant predictions come from the pinned [AbBiBench release](" + release + "). These are likelihood/energy proxies for measured affinity, not identically trained supervised affinity predictors. AeroFold uses the target assay's development labels; the external scores have no target-label fitting in this evaluation.", "",
        "## Spearman correlation", "",
        "Means of the two saved seeds. Floor-aware is AeroFold's production score; continuous is its saved latent score before floor ties. Both are shown because floor handling substantially changes CR9114/H3 correlation. The public maximum is calculated separately for each row after evaluation; it is a reference maximum, not an ensemble or a model selected for deployment.", "",
        table(["Assay", "Test", "AeroFold floor-aware", "AeroFold continuous", "Best public", "Public model"], rho_rows), "",
        "## Discovery at the same screening budget", "",
        "Select the 24 highest predicted variants. A true positive is any variant at or above the measured 95th-percentile boundary (including boundary ties). Recall = hits / all positives; enrichment = selected-positive fraction / test-positive fraction. Each number is the mean across two seeds. For a fixed assay/split, recall and enrichment have the same ordering because all models share the test pool and budget.", "",
        table(["Assay", "Test", "AeroFold recall", "Best public recall", "AeroFold enrichment", "Best public enrichment", "Public discovery model(s)"], discovery_rows), "",
        "Smaller CR6261 test pools contain 11 or 19 positives, so 24 selections can recover every positive. CR9114 tests contain 512 variants and 26 positives; their maximum possible recall is 24/26 = 92.3%. Absolute enrichment differs across pool sizes. Compare models within rows, rather than treating enrichment across antibodies as one common scale.", "",
        "## Scope and interpretation", "",
        "- **Different training information:** AeroFold uses 1,473-65,023 target-development measurements across the matched settings, plus family-excluded transfer learning. The released comparators are evaluated as frozen pretrained likelihood or energy scorers. This is an operational comparison of available predictions, not a fair contest between algorithms with equal labels or modalities.",
        "- **CR6261 framework mismatch:** the local heavy chain has Q at positions 1 and 113; every included public score uses E and K at those positions. The 11 varying residues identify each genotype uniquely, and every measured label agrees, but the physical sequence inputs differ at two fixed residues. Treat CR6261 conclusions as qualified. CR9114 has exact heavy-chain matches across all 16 varying positions and the fixed background. Public CSVs do not record the models' light-chain inputs.",
        "- **Frozen retrospective evaluation:** two seeds provide descriptive means and ranges, not confidence intervals. Higher-mutation tests may overlap between seeds. Public-library performance does not establish affinity improvement in new laboratory experiments or generalization to a new antigen.",
        "- **Coverage:** matched public scores were located for CR9114/H1, CR9114/H3, CR6261/H1 and CR6261/H9. Compatible scores for 4-4-20, the three CH65 assays, and CR9114/FluB were not located in the checked AbBiBench and FLAb releases. Those assays retain their existing AeroFold results and have no external rank here.",
        "- **AlphaFold 3 and Boltz-2:** published files were downloaded, but they lack a verified sequence/genotype crosswalk (ordinal IDs or no IDs), so they are excluded from this matched test. Their released columns are structural-confidence metrics; they are not calibrated antibody-affinity predictions. The pinned [published leaderboard](" + leaderboard + ") provides their whole-library context, using a different evaluation population.",
        "- **Newer supervised work:** [SimBinder-IF](https://arxiv.org/html/2512.17815v1) reports mean Spearman 0.410 over 11 assays and 0.294 over four held-out complexes. Its partitions, training setting, and 10-fold-improvement precision@10 differ from this protocol; compatible per-variant predictions were not located, so its aggregate values are not used in a head-to-head rank. Other affinity/ddG benchmarks on different complexes likewise do not establish a result on these tests.",
        "- **Affinity units:** rank metrics can compare differently scaled scores. Likelihood/energy scores cannot be compared directly with pKd using RMSE without separately trained calibration. This report evaluates affinity ranking and screening, not accurate absolute Kd prediction.", "",
        "## Existing AeroFold results without a matched external score", "",
    ]
    uncovered_rows = []
    names = {"titeseq_4420_binding": "4-4-20 / fluorescein", "ch65_full_ma90_binding": "CH65 / MA90", "ch65_full_si06_binding": "CH65 / SI06", "ch65_full_g189e_binding": "CH65 / G189E", "cr9114_full_flub_binding": "CR9114 / FluB"}
    for row in original[~original.task.isin(TASKS)].itertuples():
        uncovered_rows.append([names[row.task], split_names[row.split], f"{row.spearman_mean:.3f}", f"{row.latent_spearman_mean:.3f}", f"{row.top_recall_mean:.1%}", f"{row.top_enrichment_mean:.2f}x"])
    lines += [
        table(["Assay", "Test", "Floor-aware Spearman", "Continuous Spearman", "Recall", "Enrichment"], uncovered_rows), "",
        "CR9114/FluB is overwhelmingly assay-floor limited; the tie-inclusive top-5% set expands to the entire 512-variant test pool. Its 4.7% recall and 1.00x enrichment therefore do not measure useful top-binder discovery.", "",
        "## Verification and reproducibility", "",
        f"Validation passed for **{audit['assay_split_seed_model_rows']} model/assay/split/seed runs** and **{audit['paired_prediction_rows']:,} aligned prediction rows**. Sequence/genotype joins are unique, measured-library coverage is complete, and every experimental label agrees within 1e-10. All 32 whole-library public correlations reproduce the pinned leaderboard's two-decimal values. Every test metric was independently checked using average-rank Pearson correlation and scalar selection/count arithmetic. AeroFold's previously saved metrics reproduce exactly within numerical tolerance.", "",
        f"All **{audit['protected_files_unchanged']} protected notebook, data, model and result files** are unchanged. All **{audit['verified_source_file_hashes']} downloaded inputs** match their recorded hashes. AbBiBench commit: `{COMMIT}`. Frozen pipeline: `{audit['pipeline_version']}`. Source retrieval: {json.loads((HERE / 'source_manifest.json').read_text())['downloaded_utc']}.", "",
        f"Detailed model means and seed ranges: {file_link('matched_test_summary.csv')}. Individual runs: {file_link('matched_test_runs.csv')}. Aligned scores and identities: {file_link('matched_test_predictions.csv')}. Mapping and framework checks: {file_link('variant_mapping_audit.json')}. Source URLs and hashes: {file_link('source_manifest.json')}. Validation: {file_link('validation.json')}.", "",
        "To reproduce the evaluation with the existing project environment, run from the project folder:", "",
        "```powershell", ".\\.venv\\Scripts\\python.exe .\\results\\industry_affinity_comparison\\evaluate_public_scores.py", "```", "",
        "The next discriminating comparison would fine-tune or train public affinity models on exactly the same development labels, then evaluate their frozen scores on these unchanged test partitions. That work has not been performed in this evaluation.", "",
    ]
    (HERE / "comparison.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
