# AeroFold affinity comparison with published models

AeroFold's continuous Spearman exceeds the best of eight public comparators on **8 of 8 assay/split combinations**. Its top-5% recall is higher on **8**, tied on **0**, and lower on **0**. These counts describe this release and these retrospective partitions; they are not evidence of industry-wide superiority.

The evaluation reuses the frozen AeroFold GPU-trained model's saved predictions. No fitting, tuning, retraining, or changes to the working notebook were performed. All models are scored on the same experimental test genotypes within each assay/split/seed. CR9114 sequences match exactly. CR6261 has a disclosed fixed-background mismatch, described below.

Eight comparators: **ESM-IF1, ProteinMPNN, AntiFold, SaProt, ESM-2, FoldX, DiffAb, and DiffAb with a fixed backbone**. Their per-variant predictions come from the pinned [AbBiBench release](https://github.com/MSBMI-SAFE/AbBiBench/tree/7d9c73e23c3c535fbb73ef073cf2fb9633e1fda8/notebooks/scoring_outputs). These are likelihood/energy proxies for measured affinity, not identically trained supervised affinity predictors. AeroFold uses the target assay's development labels; the external scores have no target-label fitting in this evaluation.

## Spearman correlation

Means of the two saved seeds. Floor-aware is AeroFold's production score; continuous is its saved latent score before floor ties. Both are shown because floor handling substantially changes CR9114/H3 correlation. The public maximum is calculated separately for each row after evaluation; it is a reference maximum, not an ensemble or a model selected for deployment.

| Assay | Test | AeroFold floor-aware | AeroFold continuous | Best public | Public model |
| --- | --- | --- | --- | --- | --- |
| CR9114 / H1 | Interpolation | 0.968 | 0.968 | 0.636 | FoldX |
| CR9114 / H1 | Higher-mutation extrapolation | 0.913 | 0.913 | 0.531 | SaProt |
| CR9114 / H3 | Interpolation | 0.805 | 0.505 | 0.457 | ESM-IF1 |
| CR9114 / H3 | Higher-mutation extrapolation | 0.913 | 0.847 | 0.705 | ESM-IF1 |
| CR6261 / H1 | Interpolation | 0.963 | 0.963 | 0.663 | DiffAb |
| CR6261 / H1 | Higher-mutation extrapolation | 0.931 | 0.931 | 0.670 | DiffAb |
| CR6261 / H9 | Interpolation | 0.962 | 0.962 | 0.774 | DiffAb fixed backbone |
| CR6261 / H9 | Higher-mutation extrapolation | 0.960 | 0.960 | 0.721 | DiffAb fixed backbone |

## Discovery at the same screening budget

Select the 24 highest predicted variants. A true positive is any variant at or above the measured 95th-percentile boundary (including boundary ties). Recall = hits / all positives; enrichment = selected-positive fraction / test-positive fraction. Each number is the mean across two seeds. For a fixed assay/split, recall and enrichment have the same ordering because all models share the test pool and budget.

| Assay | Test | AeroFold recall | Best public recall | AeroFold enrichment | Best public enrichment | Public discovery model(s) |
| --- | --- | --- | --- | --- | --- | --- |
| CR9114 / H1 | Interpolation | 48.1% | 21.2% | 10.26x | 4.51x | ESM-IF1, SaProt |
| CR9114 / H1 | Higher-mutation extrapolation | 30.8% | 7.7% | 6.56x | 1.64x | AntiFold |
| CR9114 / H3 | Interpolation | 78.8% | 55.8% | 16.82x | 11.90x | ESM-IF1 |
| CR9114 / H3 | Higher-mutation extrapolation | 76.9% | 42.3% | 16.41x | 9.03x | ESM-IF1 |
| CR6261 / H1 | Interpolation | 89.5% | 36.8% | 14.09x | 5.80x | DiffAb |
| CR6261 / H1 | Higher-mutation extrapolation | 100.0% | 27.3% | 8.71x | 2.38x | ESM-IF1, DiffAb |
| CR6261 / H9 | Interpolation | 86.8% | 52.6% | 13.35x | 8.09x | ProteinMPNN |
| CR6261 / H9 | Higher-mutation extrapolation | 86.4% | 54.5% | 7.48x | 4.73x | SaProt |

Smaller CR6261 test pools contain 11 or 19 positives, so 24 selections can recover every positive. CR9114 tests contain 512 variants and 26 positives; their maximum possible recall is 24/26 = 92.3%. Absolute enrichment differs across pool sizes. Compare models within rows, rather than treating enrichment across antibodies as one common scale.

## Scope and interpretation

- **Different training information:** AeroFold uses 1,473-65,023 target-development measurements across the matched settings, plus family-excluded transfer learning. The released comparators are evaluated as frozen pretrained likelihood or energy scorers. This is an operational comparison of available predictions, not a fair contest between algorithms with equal labels or modalities.
- **CR6261 framework mismatch:** the local heavy chain has Q at positions 1 and 113; every included public score uses E and K at those positions. The 11 varying residues identify each genotype uniquely, and every measured label agrees, but the physical sequence inputs differ at two fixed residues. Treat CR6261 conclusions as qualified. CR9114 has exact heavy-chain matches across all 16 varying positions and the fixed background. Public CSVs do not record the models' light-chain inputs.
- **Frozen retrospective evaluation:** two seeds provide descriptive means and ranges, not confidence intervals. Higher-mutation tests may overlap between seeds. Public-library performance does not establish affinity improvement in new laboratory experiments or generalization to a new antigen.
- **Coverage:** matched public scores were located for CR9114/H1, CR9114/H3, CR6261/H1 and CR6261/H9. Compatible scores for 4-4-20, the three CH65 assays, and CR9114/FluB were not located in the checked AbBiBench and FLAb releases. Those assays retain their existing AeroFold results and have no external rank here.
- **AlphaFold 3 and Boltz-2:** published files were downloaded, but they lack a verified sequence/genotype crosswalk (ordinal IDs or no IDs), so they are excluded from this matched test. Their released columns are structural-confidence metrics; they are not calibrated antibody-affinity predictions. The pinned [published leaderboard](https://github.com/MSBMI-SAFE/AbBiBench/blob/7d9c73e23c3c535fbb73ef073cf2fb9633e1fda8/README.md) provides their whole-library context, using a different evaluation population.
- **Newer supervised work:** [SimBinder-IF](https://arxiv.org/html/2512.17815v1) reports mean Spearman 0.410 over 11 assays and 0.294 over four held-out complexes. Its partitions, training setting, and 10-fold-improvement precision@10 differ from this protocol; compatible per-variant predictions were not located, so its aggregate values are not used in a head-to-head rank. Other affinity/ddG benchmarks on different complexes likewise do not establish a result on these tests.
- **Affinity units:** rank metrics can compare differently scaled scores. Likelihood/energy scores cannot be compared directly with pKd using RMSE without separately trained calibration. This report evaluates affinity ranking and screening, not accurate absolute Kd prediction.

## Existing AeroFold results without a matched external score

| Assay | Test | Floor-aware Spearman | Continuous Spearman | Recall | Enrichment |
| --- | --- | --- | --- | --- | --- |
| 4-4-20 / fluorescein | Interpolation | 0.723 | 0.643 | 61.5% | 13.13x |
| 4-4-20 / fluorescein | Higher-mutation extrapolation | 0.412 | 0.309 | 56.2% | 7.08x |
| CH65 / MA90 | Interpolation | 0.989 | 0.989 | 75.0% | 16.00x |
| CH65 / MA90 | Higher-mutation extrapolation | 0.948 | 0.948 | 57.7% | 12.31x |
| CH65 / SI06 | Interpolation | 0.949 | 0.930 | 75.0% | 16.00x |
| CH65 / SI06 | Higher-mutation extrapolation | 0.954 | 0.954 | 55.8% | 11.90x |
| CH65 / G189E | Interpolation | 0.985 | 0.985 | 73.1% | 15.59x |
| CH65 / G189E | Higher-mutation extrapolation | 0.944 | 0.944 | 48.1% | 10.26x |
| CR9114 / FluB | Interpolation | 0.645 | 0.146 | 4.7% | 1.00x |
| CR9114 / FluB | Higher-mutation extrapolation | -0.028 | -0.028 | 4.7% | 1.00x |

CR9114/FluB is overwhelmingly assay-floor limited; the tie-inclusive top-5% set expands to the entire 512-variant test pool. Its 4.7% recall and 1.00x enrichment therefore do not measure useful top-binder discovery.

## Verification and reproducibility

Validation passed for **144 model/assay/split/seed runs** and **57,816 aligned prediction rows**. Sequence/genotype joins are unique, measured-library coverage is complete, and every experimental label agrees within 1e-10. All 32 whole-library public correlations reproduce the pinned leaderboard's two-decimal values. Every test metric was independently checked using average-rank Pearson correlation and scalar selection/count arithmetic. AeroFold's previously saved metrics reproduce exactly within numerical tolerance.

All **78 protected notebook, data, model and result files** are unchanged. All **37 downloaded inputs** match their recorded hashes. AbBiBench commit: `7d9c73e23c3c535fbb73ef073cf2fb9633e1fda8`. Frozen pipeline: `aerofold-attention-gpu-5`. Source retrieval: 2026-10-03T02:46:10.125653+00:00.

Detailed model means and seed ranges: [matched_test_summary.csv](<C:/Users/alokd/Downloads/for training aerofold model/results/industry_affinity_comparison/matched_test_summary.csv>). Individual runs: [matched_test_runs.csv](<C:/Users/alokd/Downloads/for training aerofold model/results/industry_affinity_comparison/matched_test_runs.csv>). Aligned scores and identities: [matched_test_predictions.csv](<C:/Users/alokd/Downloads/for training aerofold model/results/industry_affinity_comparison/matched_test_predictions.csv>). Mapping and framework checks: [variant_mapping_audit.json](<C:/Users/alokd/Downloads/for training aerofold model/results/industry_affinity_comparison/variant_mapping_audit.json>). Source URLs and hashes: [source_manifest.json](<C:/Users/alokd/Downloads/for training aerofold model/results/industry_affinity_comparison/source_manifest.json>). Validation: [validation.json](<C:/Users/alokd/Downloads/for training aerofold model/results/industry_affinity_comparison/validation.json>).

To reproduce the evaluation with the existing project environment, run from the project folder:

```powershell
.\.venv\Scripts\python.exe .\results\industry_affinity_comparison\evaluate_public_scores.py
```

The next discriminating comparison would fine-tune or train public affinity models on exactly the same development labels, then evaluate their frozen scores on these unchanged test partitions. That work has not been performed in this evaluation.
