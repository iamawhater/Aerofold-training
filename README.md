# AeroFold antibody discovery notebook

Open **AeroFold_lab.ipynb** with the project `.venv` Python kernel. This single notebook reconstructs public libraries, learns across antibody families, adapts to each antigen, evaluates held-out variants, trains production models and ranks candidate sequences.

The corpus contains **396,688 antigen-specific measurement records** across nine tasks. Different antigens can share antibody sequences; this is not a count of distinct antibodies or independent experiments.

| Antibody | Antigen tasks | Finite measured variants per task |
| --- | --- | --- |
| 4-4-20 | Fluorescein | 2,807 |
| CH65 | MA90, SI06, G189E | 65,530; 64,619; 63,840 |
| CR9114 | H1, H3, influenza B | 65,094; 65,535; 65,534 |
| CR6261 | H1, H9 | 1,887; 1,842 |

## Learning and evaluation

Supervised source training learns a shared neural mutation representation with separate antigen/assay heads. Inputs describe substitutions, biochemical changes, position bins and local sequence context. Each antibody family has a source checkpoint trained on **every other family's measurements**. All assays of the evaluated family are excluded, preventing related-antigen measurements from revealing test sequences. Loss weighting gives each source assay equal total weight; all source rows enter every epoch.

Target experts include additive regression, nonlinear global epistasis, ordinal probabilities, pairwise ranking, pair and triple interactions, boosting, screening-focused LambdaRank trees, transferred-feature Ridge, adapted neural models, and a new two-layer, four-head mutation-site transformer with an additive skip path. The verified existing AbLang-2 cache adds a frozen foundation representation for 4-4-20. Source pretraining, target neural adaptation and attention train on the NVIDIA RTX 3050 Ti with CUDA 12.8; device and allocation checks are recorded. The new attention model trains from scratch on target-assay labels; it is not a fine-tuned pretrained protein language model.

Attention trains on rank regression, tied-pair ranking with extra top-binder weights, and auxiliary measured top-5% classification. All eligible development rows contribute; epoch counts adapt to the fixed update target, with at least three full passes. Ranking trees train shuffled same-assay groups of up to 512 variants and optimize the top 24. Different antigens are never ranked against one another.

Inner selection assigns 40% of its objective to Spearman and 60% to attainable top-24 hit efficiency. Attention enters only if its augmented mixture beats a matched attention-free ensemble by more than 0.005 on internal validation; the baseline has the same screening trees, folds and objective. Source transfer uses the same gain gate against a target-only ensemble. The completed run selects attention in 2 of nine final models, screening trees in 0, and source transfer in 0. None of these experts is assumed to improve every held-out assay.

Two random validation folds and one strict higher-mutation fold select small ensembles and assay-floor cutoffs. Internal training is capped at 4,000 examples per fold and validation at 768 for CPU practicality. Selected experts then fit **all eligible development measurements**. Final models fit **all finite measurements** after benchmarking; there is no 5,000-label restriction.

The benchmark covers nine tasks, two seeds, and two splits:

- Random interpolation: held-out sequence pairs from the same measured library.
- Higher-mutation extrapolation: test variants have strictly more mutations from the published ancestor than every development variant. The fixed 4-4-20 reference is its sequence-only library consensus.

Floor-aware and continuous Spearman are both reported. Measured top-5% recall includes boundary ties, with a screening budget of 24 variants; enrichment uses the same hits. Large assay-floor ties can expand the positive set far beyond 5%. Recovery of uncensored binders and exact positive counts are recorded separately. Recall ceilings are reported: 26 positives and 24 slots permit at most 92.3% recall. Constant-label tests have undefined correlation. For exact Spearman reproduction from saved predictions, use pandas `read_csv(..., float_precision='round_trip')` to preserve numerical ties. Seed ranges are descriptive, not confidence intervals; high-mutation test identities overlap across seeds.

## Run and rank candidates

1. Select `.venv/Scripts/python.exe` in your notebook editor. To recreate the environment, use Python 3.12 and install `requirements.txt`, which pins the CUDA 12.8 PyTorch build. Training defaults to `REQUIRE_GPU=True` and stops if CUDA is unavailable; saved-model inference remains portable on CPU.
2. Run the notebook top to bottom. Raw inputs are included and verified by recorded hashes. Source checkpoints are reused only when signatures match.
3. Inspect the executed recall/enrichment benchmark, attention contribution and transfer contribution figures. Legacy or screenshot targets do not determine selection.
4. Nine final artifacts are saved in `results/multifamily/`.
5. Set `CANDIDATE_TASK` and `CANDIDATES_CSV` to rank a CSV with `variant_id,heavy_sequence,light_sequence`. Use complete canonical domains matching the chosen family reference layout.

For prediction without retraining, set `RUN_BENCHMARK=False` and `TRAIN_FINAL_MODELS=False`. Run Setup, data reconstruction, model definitions and embedding cells, then candidate ranking; the saved model loads automatically. Novel AbLang-2 embeddings, if required by the model, use installed local weights.

`rank_score` is dimensionless and higher-is-better. `latent_rank_score` retains continuous ordering before floor clipping; `predicted_floor_tie` identifies predicted assay-floor ties. `estimated_assay_score` is an approximate mapping from inner validation. Known measurements and membership in training are flagged. These models do not supply calibrated uncertainty intervals or predict arbitrary new antigens.

## Sources and artifacts

CH65 comes from the [original repository](https://github.com/amphilli/CH65-comblib) and [paper](https://doi.org/10.7554/eLife.83628). CR9114 and CR6261 use [publication source tables and supplementary plasmids](https://elifesciences.org/articles/71393/figures) and the [study repository](https://github.com/klawrence26/bnab-landscapes). The notebook translates original constructs, removes tags/linkers, reconstructs every genotype and checks published bits and mutation counts. Original local CSVs are preserved; the full original CR tables feed this run.

- `data/raw/download_manifest.json`: source URLs, sizes and hashes.
- `data/*_full_*_binding.csv`: reconstructed measured libraries.
- `cache/transfer_*.pt`: supervised source encoders, exclusions, training history and signatures.
- `results/multifamily/benchmark_summary.csv`, `benchmark_runs.csv`, `heldout_predictions.csv`: results and predictions.
- `results/multifamily/attention_contribution.csv`: matched attention-free comparison.
- `results/multifamily/attention_training.csv`: trained transformer diagnostics.
- `results/multifamily/transfer_contribution.csv`: separately selected target-only control; attention admission may also differ.
- `results/multifamily/inner_selection.csv`, `inner_validation.csv`: model-selection evidence.
- `results/multifamily/*_model.pkl`, `*_model_card.json`: final models and provenance. Load only trusted pickle files.
- `results/multifamily/evaluation_protocol.json`, `benchmark_audit.json`, `validation.json`: protocol, partitions and checks.
- `results/multifamily/notebook_preview.html`: executed preview.

Retrospective public-library results cannot establish universal or prospective superiority. Predicting assay-floor ties does not recover affinity ordering below the detection limit. The separate web-app project is unchanged. The pre-attention notebook and models are recoverable in `archive/before_attention_experiment.zip`. The outer tests were viewed in that earlier run; this is a retrospective rerun, not an independent prospective confirmation.
