# Executed multi-family model findings
Separate antigen-specific models; supervised transfer excludes the complete benchmark antibody family.
- titeseq_4420_binding / Interpolation: Spearman 0.723, continuous Spearman 0.643, measured top-5% recall 61.5%, enrichment 13.13x.
- titeseq_4420_binding / Higher-mutation extrapolation: Spearman 0.412, continuous Spearman 0.309, measured top-5% recall 56.2%, enrichment 7.08x.
- ch65_full_ma90_binding / Interpolation: Spearman 0.989, continuous Spearman 0.989, measured top-5% recall 75.0%, enrichment 16.00x.
- ch65_full_ma90_binding / Higher-mutation extrapolation: Spearman 0.948, continuous Spearman 0.948, measured top-5% recall 57.7%, enrichment 12.31x.
- ch65_full_si06_binding / Interpolation: Spearman 0.949, continuous Spearman 0.930, measured top-5% recall 75.0%, enrichment 16.00x.
- ch65_full_si06_binding / Higher-mutation extrapolation: Spearman 0.954, continuous Spearman 0.954, measured top-5% recall 55.8%, enrichment 11.90x.
- ch65_full_g189e_binding / Interpolation: Spearman 0.985, continuous Spearman 0.985, measured top-5% recall 73.1%, enrichment 15.59x.
- ch65_full_g189e_binding / Higher-mutation extrapolation: Spearman 0.944, continuous Spearman 0.944, measured top-5% recall 48.1%, enrichment 10.26x.
- cr9114_full_h1_binding / Interpolation: Spearman 0.968, continuous Spearman 0.968, measured top-5% recall 48.1%, enrichment 10.26x.
- cr9114_full_h1_binding / Higher-mutation extrapolation: Spearman 0.913, continuous Spearman 0.913, measured top-5% recall 30.8%, enrichment 6.56x.
- cr9114_full_h3_binding / Interpolation: Spearman 0.805, continuous Spearman 0.505, measured top-5% recall 78.8%, enrichment 16.82x.
- cr9114_full_h3_binding / Higher-mutation extrapolation: Spearman 0.913, continuous Spearman 0.847, measured top-5% recall 76.9%, enrichment 16.41x.
- cr9114_full_flub_binding / Interpolation: Spearman 0.645, continuous Spearman 0.146, measured top-5% recall 4.7%, enrichment 1.00x.
- cr9114_full_flub_binding / Higher-mutation extrapolation: Spearman -0.028, continuous Spearman -0.028, measured top-5% recall 4.7%, enrichment 1.00x.
- cr6261_full_h1_binding / Interpolation: Spearman 0.963, continuous Spearman 0.963, measured top-5% recall 89.5%, enrichment 14.09x.
- cr6261_full_h1_binding / Higher-mutation extrapolation: Spearman 0.931, continuous Spearman 0.931, measured top-5% recall 100.0%, enrichment 8.71x.
- cr6261_full_h9_binding / Interpolation: Spearman 0.962, continuous Spearman 0.962, measured top-5% recall 86.8%, enrichment 13.35x.
- cr6261_full_h9_binding / Higher-mutation extrapolation: Spearman 0.960, continuous Spearman 0.960, measured top-5% recall 86.4%, enrichment 7.48x.

Attention and source transfer can help, receive zero weight, or reduce individual metrics. Inspect the matched attention-free and target-only comparisons.
Top-5% labels are tie-inclusive. Floor-dominated assays may expand this set far beyond 5%; uncensored-binder recall is reported separately.
Undefined correlations identify constant test labels. Public-library results do not establish prospective or universal superiority.
Final production models fit all finite measured variants after evaluation. The frozen source checkpoints contain only other antibody families.