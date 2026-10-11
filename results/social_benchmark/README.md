# AeroFold social benchmark figures

Both charts use higher-mutation holdouts and means of two saved seeds. No model was retrained.

- `aerofold_vs_additive_all_antibodies.png`: all four antibody families and all nine assays; mutation-additive regression is a target-supervised internal baseline.
- `aerofold_vs_public_baselines.png`: the four assays with matched public scores, across CR9114 and CR6261; the best public scorer is selected separately after evaluation for each assay and metric.
- PNGs are intended for sharing; SVGs keep labels editable. Exact plotted rows are in the companion CSVs.

## Interpretation

Continuous Spearman avoids the change in correlation produced by predicted assay-floor ties. Top-5% recovery is tie-inclusive and uses 24 picks. FluB recall is not informative because every test row qualifies; the raw 4.6875% is preserved in the CSV, while the figure explicitly marks that row as not informative.

AeroFold has target-assay supervision; external likelihood/energy scorers do not receive target-label fitting here. CR6261 public and local inputs differ at two fixed heavy-chain residues (Q1/Q113 vs E1/K113). These are previously inspected retrospective public-library tests, with overlapping high-mutation test identities across two seeds. Seed means are descriptive, not confidence intervals, and there is no prospective laboratory confirmation. CR9114 heavy chains match exactly; public score CSVs do not record light-chain inputs.

The public source is the pinned AbBiBench release at commit `7d9c73e23c3c535fbb73ef073cf2fb9633e1fda8`; detailed provenance is preserved in `../industry_affinity_comparison/source_manifest.json` and `comparison.md`.

Reproduce from the project root: `.venv/Scripts/python.exe results/social_benchmark/create_charts.py`. Checks and input hashes are in `chart_validation.json`.
