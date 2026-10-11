"""Render source-backed AeroFold figures for a social post; no model fitting."""

from pathlib import Path
import hashlib
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
PUBLIC = ROOT / "results" / "industry_affinity_comparison"
LOCAL = ROOT / "results" / "multifamily"
SPLIT = "mutation_extrapolation"
COLORS = {"AeroFold": "#2366B1", "baseline": "#B88322", "ink": "#202A35", "muted": "#5F6872", "grid": "#E2E6EB"}

TASKS = [
    ("titeseq_4420_binding", "4-4-20 / Fluorescein"),
    ("ch65_full_ma90_binding", "CH65 / MA90"),
    ("ch65_full_si06_binding", "CH65 / SI06"),
    ("ch65_full_g189e_binding", "CH65 / G189E"),
    ("cr9114_full_h1_binding", "CR9114 / H1"),
    ("cr9114_full_h3_binding", "CR9114 / H3"),
    ("cr9114_full_flub_binding", "CR9114 / FluB"),
    ("cr6261_full_h1_binding", "CR6261 / H1"),
    ("cr6261_full_h9_binding", "CR6261 / H9"),
]

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 14,
    "text.color": COLORS["ink"], "axes.labelcolor": COLORS["ink"],
    "xtick.color": COLORS["muted"], "ytick.color": COLORS["ink"],
    "svg.fonttype": "none", "figure.facecolor": "white", "axes.facecolor": "white",
})


def read(path):
    return pd.read_csv(path, float_precision="round_trip")


def require_equal(actual, expected, context):
    if not np.isclose(float(actual), float(expected), atol=1e-12, rtol=1e-10):
        raise ValueError(f"Mismatch in {context}: {actual} != {expected}")


def verify_saved_means(summary, runs):
    selected = runs[runs["split"].eq(SPLIT)]
    if not (selected.groupby(["task", "model"]).size() == 2).all():
        raise ValueError("Expected exactly two seed runs per plotted model/task")
    grouped = selected.groupby(["task", "model"])
    for _, row in summary[summary["split"].eq(SPLIT)].iterrows():
        group = grouped.get_group((row["task"], row["model"]))
        require_equal(row["latent_spearman_mean"], group["latent_spearman"].mean(), "continuous Spearman seed mean")
        require_equal(row["top_recall_mean"], group["top_recall"].mean(), "recall seed mean")
        if not (group["top_k"] == 24).all():
            raise ValueError("Screening budget differs from chart caption")
        for _, run in group.iterrows():
            require_equal(run["top_recall"], run["top_hits"] / run["top_population"], "hit/positive recall")


def counts_for(runs, task):
    subset = runs[runs["task"].eq(task) & runs["split"].eq(SPLIT) & runs["model"].eq("AeroFold")]
    tested = sorted(subset["test_count"].unique())
    positives = sorted(subset["top_population"].unique())
    if len(tested) != 1 or len(positives) != 1:
        raise ValueError("Denominators vary; adapt the figure rather than displaying one count")
    return int(tested[0]), int(positives[0])


def prepare_internal():
    summary = read(LOCAL / "benchmark_summary.csv")
    runs = read(LOCAL / "benchmark_runs.csv")
    verify_saved_means(summary, runs)
    table = summary[summary["split"].eq(SPLIT)].set_index(["task", "model"])
    rows = []
    for task, label in TASKS:
        aero, baseline = table.loc[(task, "AeroFold")], table.loc[(task, "Additive")]
        tested, positives = counts_for(runs, task)
        rows.append({
            "task": task, "assay": label, "split": SPLIT,
            "aerofold_continuous_spearman": aero["latent_spearman_mean"],
            "baseline_continuous_spearman": baseline["latent_spearman_mean"],
            "aerofold_recall": aero["top_recall_mean"], "baseline_recall": baseline["top_recall_mean"],
            "baseline_spearman_model": "Mutation-additive regression",
            "baseline_recall_model": "Mutation-additive regression",
            "test_count": tested, "positive_count": positives,
            "recall_informative": positives != tested,
            "framework_mismatch": False,
        })
    return pd.DataFrame(rows)


def prepare_public():
    overview = read(PUBLIC / "comparison_overview.csv")
    summary = read(PUBLIC / "matched_test_summary.csv")
    runs = read(PUBLIC / "matched_test_runs.csv")
    verify_saved_means(summary, runs)
    rows = []
    for _, row in overview[overview["split"].eq(SPLIT)].iterrows():
        task = row["task"]
        group = summary[summary["task"].eq(task) & summary["split"].eq(SPLIT)]
        aero = group[group["model"].eq("AeroFold")].iloc[0]
        comparators = group[~group["model"].eq("AeroFold")]
        if len(comparators) != 8:
            raise ValueError("Expected eight public scorers")
        require_equal(row["aerofold_continuous_spearman"], aero["latent_spearman_mean"], "public overview AeroFold score")
        require_equal(row["aerofold_recall"], aero["top_recall_mean"], "public overview AeroFold recall")
        require_equal(row["best_public_spearman"], comparators["spearman_mean"].max(), "best public correlation")
        require_equal(row["best_public_recall"], comparators["top_recall_mean"].max(), "best public recall")
        tested, positives = counts_for(runs, task)
        rows.append({
            "task": task, "assay": row["assay"], "split": SPLIT,
            "aerofold_continuous_spearman": row["aerofold_continuous_spearman"],
            "baseline_continuous_spearman": row["best_public_spearman"],
            "aerofold_recall": row["aerofold_recall"], "baseline_recall": row["best_public_recall"],
            "baseline_spearman_model": row["best_public_spearman_model"],
            "baseline_recall_model": row["best_public_recall_models"],
            "test_count": tested, "positive_count": positives,
            "recall_informative": True,
            "framework_mismatch": task.startswith("cr6261"),
        })
    return pd.DataFrame(rows)


def model_label(value):
    return value.replace("DiffAb fixed backbone", "DiffAb (fixed)").replace("; ", " / ")


def add_bar(ax, y, value, aero, is_percent=False, model=None):
    color = COLORS["AeroFold"] if aero else COLORS["baseline"]
    if abs(value) > 1e-12:
        ax.barh(y, value, height=.24, color=color if aero else "white", edgecolor=color, linewidth=1.4, zorder=3)
    else:
        ax.plot(0, y, "|", color=color, markersize=15, markeredgewidth=1.6, zorder=4)
    label = f"{value:.1f}%" if is_percent else f"{value:.3f}"
    if model:
        label += f"  {model_label(model)}"
    pad = 1.5 if is_percent else .015
    if value < 0:
        # Keep a short negative bar's label clear of the category text.
        x, align = .025, "left"
    else:
        x, align = value + pad, "left"
    ax.text(x, y, label, va="center", ha=align, fontsize=12.3 if model else 14.2,
            color=COLORS["ink"] if aero else COLORS["muted"], zorder=5)


def render(rows, public=False):
    fig = plt.figure(figsize=(16, 9.8 if public else 13), dpi=140)
    bottom, top = (.31, .70) if public else (.245, .76)
    axes = [fig.add_axes([.245, bottom, .285, top-bottom]), fig.add_axes([.61, bottom, .34, top-bottom])]
    fig.text(.055, .946, "AeroFold vs public scorers" if public else "AeroFold across four antibody families", fontsize=28, weight="bold")
    fig.text(.055, .898 if public else .907, "Higher-mutation holdouts | 4 matched assays" if public else "Higher-mutation holdouts | all 9 measured assays", fontsize=18, color=COLORS["muted"])
    handles = [Patch(facecolor=COLORS["AeroFold"], edgecolor=COLORS["AeroFold"], label="AeroFold"),
               Patch(facecolor="white", edgecolor=COLORS["baseline"], linewidth=1.5,
                     label="Best of 8 public scorers (per metric)" if public else "Mutation-additive regression")]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(.055, .868 if public else .871),
               frameon=False, ncol=2, fontsize=15.2, handlelength=1.5, columnspacing=2.1, borderaxespad=0)
    headings = [("Continuous Spearman", "Affinity ranking; higher is better"),
                ("Top-5% recall", "24 picks; higher is better")]
    n = len(rows)
    for panel, ax in enumerate(axes):
        pos = ax.get_position()
        fig.text(pos.x0, top+.059 if public else top+.042, headings[panel][0], fontsize=20, weight="bold")
        fig.text(pos.x0, top+.025 if public else top+.018, headings[panel][1], fontsize=12.7, color=COLORS["muted"])
        ax.set_ylim(n-.49, -.58)
        ax.set_yticks(range(n))
        labels = []
        for _, row in rows.iterrows():
            labels.append(row["assay"] + (" *" if row["framework_mismatch"] else "")
                          + f"\n{row['positive_count']} positives / {row['test_count']} tested")
        ax.set_yticklabels(labels if panel == 0 else [""] * n, fontsize=13.7, linespacing=1.35)
        ax.tick_params(axis="y", length=0, pad=12)
        ax.tick_params(axis="x", length=0, labelsize=12.2, pad=8)
        ax.set_axisbelow(True)
        ax.grid(axis="x", color=COLORS["grid"], linewidth=.8)
        for side in ["top", "right", "left"]:
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_color("#AEB6BF")
        ax.spines["bottom"].set_linewidth(.8)
        if panel == 0:
            ax.set_xlim((0, 1.25) if public else (-.17, 1.15))
            ax.set_xticks([0, .25, .5, .75, 1] if public else [-.1, 0, .25, .5, .75, 1])
            if not public:
                ax.set_xticklabels(["−0.1", "0", "0.25", "0.50", "0.75", "1.00"])
            ax.axvline(0, color="#949DA7", lw=.8, zorder=2)
            ax.set_xlabel("Spearman correlation", labelpad=10, fontsize=13)
        else:
            ax.set_xlim(0, 117 if public else 116)
            ax.set_xticks([0, 25, 50, 75, 100])
            ax.xaxis.set_major_formatter(PercentFormatter(100, decimals=0))
            ax.set_xlabel("Measured top-5% positives recovered", labelpad=10, fontsize=13)
        for i, row in rows.iterrows():
            if panel == 1 and not row["recall_informative"]:
                ax.text(3, i, "Not informative: all 512 labels qualify", fontsize=12.5,
                        va="center", color=COLORS["muted"])
                continue
            a = row["aerofold_recall"] * 100 if panel == 1 else row["aerofold_continuous_spearman"]
            b = row["baseline_recall"] * 100 if panel == 1 else row["baseline_continuous_spearman"]
            model = row["baseline_recall_model"] if panel == 1 else row["baseline_spearman_model"]
            add_bar(ax, i-.155, a, True, is_percent=panel == 1)
            add_bar(ax, i+.155, b, False, is_percent=panel == 1, model=model if public else None)
        # Separate families without relying on color.
        for i in range(n-1):
            if rows.iloc[i]["task"].split("_")[0] != rows.iloc[i+1]["task"].split("_")[0]:
                ax.axhline(i+.5, color=COLORS["grid"], lw=.8)

    if public:
        notes = [
            "AeroFold uses target-assay labels; public scores are frozen likelihood/energy proxies without target-label fitting.",
            "* CR6261: two fixed heavy-chain residues differ between local and public inputs (Q1/Q113 vs E1/K113).",
            "Baselines are post-hoc maxima for each assay and metric; names label the gold bars. This is not one combined model.",
            "Mean of 2 seeds; descriptive, with overlapping test identities. Previously inspected retrospective tests; no prospective validation.",
            "Top-5% positives include boundary ties. CR9114 recall ceiling: 92.3%; CR6261: 100%. Source: saved AbBiBench comparison.",
        ]
        note_top, step = .197, .031
        stem = "aerofold_vs_public_baselines"
    else:
        notes = [
            "Both models use target-assay development labels. Baseline: saved mutation-additive regression on the same test partitions.",
            "Mean of 2 seeds; descriptive, with overlapping test identities. Previously inspected retrospective tests; no prospective validation.",
            "Top-5% positives include boundary ties. 24 picks can recover every positive only when the positive set has 24 or fewer variants.",
            "CR9114/FluB: every test variant qualifies because of assay-floor ties; both recalls are 4.7% and enrichment is 1.0x.",
            "Sources: multifamily/benchmark_summary.csv and benchmark_runs.csv. Continuous Spearman uses scores before floor clipping.",
        ]
        note_top, step = .161, .024
        stem = "aerofold_vs_additive_all_antibodies"
    for j, note in enumerate(notes):
        fig.text(.055, note_top-j*step, note, fontsize=10.4, color=COLORS["muted"])
    fig.text(.055, .036 if public else .023, "AEROFOLD  /  RETROSPECTIVE PUBLIC-LIBRARY BENCHMARK", fontsize=9.5,
             color=COLORS["muted"], weight="bold")
    # Check text bounds against the actual renderer before export.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    width, height = fig.canvas.get_width_height()
    for ax in axes:
        tick_bounds = [tick.get_window_extent(renderer) for tick in ax.get_xticklabels() if tick.get_visible()]
        if any(left.x1 + 2 > right.x0 for left, right in zip(tick_bounds, tick_bounds[1:])):
            raise ValueError("Overlapping quantitative-axis tick labels")
    for text in fig.findobj(match=matplotlib.text.Text):
        if not text.get_visible() or not text.get_text():
            continue
        bounds = text.get_window_extent(renderer)
        if bounds.x0 < -1 or bounds.y0 < -1 or bounds.x1 > width+1 or bounds.y1 > height+1:
            raise ValueError(f"Text outside canvas: {text.get_text()!r}, {bounds.bounds}")
    fig.savefig(OUT / f"{stem}.png", dpi=140, facecolor="white")
    fig.savefig(OUT / f"{stem}.svg", facecolor="white")
    plt.close(fig)
    return stem


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    internal = prepare_internal()
    public = prepare_public()
    internal.to_csv(OUT / "additive_chart_data.csv", index=False)
    public.to_csv(OUT / "public_chart_data.csv", index=False)
    files = [render(internal), render(public, public=True)]
    tweet = (
        "AeroFold: supervised antibody affinity ranking across 4-4-20, CH65, CR9114 & CR6261.\n\n"
        "Charts: retrospective higher-mutation holdouts, Spearman + top-5% recall (24 picks).\n\n"
        "ML/bio folks: how do these results look? Which baselines or split designs would you add?"
    )
    if len(tweet) > 280:
        raise ValueError(f"Tweet exceeds 280 characters: {len(tweet)}")
    (OUT / "tweet.txt").write_text(tweet + "\n", encoding="utf-8")
    inputs = [LOCAL / "benchmark_summary.csv", LOCAL / "benchmark_runs.csv", PUBLIC / "comparison_overview.csv",
              PUBLIC / "matched_test_summary.csv", PUBLIC / "matched_test_runs.csv"]
    validation = {
        "status": "passed", "plot_rows": {"additive": len(internal), "public": len(public)},
        "plotted_split": SPLIT, "tweet_characters": len(tweet),
        "checks": ["saved means agree with two seed runs", "recall recomputed as hits / positives",
                   "24-pick budget verified", "public baselines equal eight-comparator maxima",
                   "all nine internal assays included", "rendered text inside canvas"],
        "input_sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs},
    }
    (OUT / "chart_validation.json").write_text(json.dumps(validation, indent=2) + "\n", encoding="utf-8")
    (OUT / "README.md").write_text(
        "# AeroFold social benchmark figures\n\n"
        "Both charts use higher-mutation holdouts and means of two saved seeds. No model was retrained.\n\n"
        "- `aerofold_vs_additive_all_antibodies.png`: all four antibody families and all nine assays; "
        "mutation-additive regression is a target-supervised internal baseline.\n"
        "- `aerofold_vs_public_baselines.png`: the four assays with matched public scores, across CR9114 and CR6261; "
        "the best public scorer is selected separately after evaluation for each assay and metric.\n"
        "- PNGs are intended for sharing; SVGs keep labels editable. Exact plotted rows are in the companion CSVs.\n\n"
        "## Interpretation\n\n"
        "Continuous Spearman avoids the change in correlation produced by predicted assay-floor ties. "
        "Top-5% recovery is tie-inclusive and uses 24 picks. FluB recall is not informative because every test row qualifies; "
        "the raw 4.6875% is preserved in the CSV, while the figure explicitly marks that row as not informative.\n\n"
        "AeroFold has target-assay supervision; external likelihood/energy scorers do not receive target-label fitting here. "
        "CR6261 public and local inputs differ at two fixed heavy-chain residues (Q1/Q113 vs E1/K113). "
        "These are previously inspected retrospective public-library tests, with overlapping high-mutation test identities across "
        "two seeds. Seed means are descriptive, not confidence intervals, and there is no prospective laboratory confirmation. "
        "CR9114 heavy chains match exactly; public score CSVs do not record light-chain inputs.\n\n"
        "The public source is the pinned AbBiBench release at commit "
        "`7d9c73e23c3c535fbb73ef073cf2fb9633e1fda8`; detailed provenance is preserved in "
        "`../industry_affinity_comparison/source_manifest.json` and `comparison.md`.\n\n"
        "Reproduce from the project root: `.venv/Scripts/python.exe results/social_benchmark/create_charts.py`. "
        "Checks and input hashes are in `chart_validation.json`.\n",
        encoding="utf-8",
    )
    print(json.dumps({"figures": files, "tweet_characters": len(tweet), "validation": "passed"}))


if __name__ == "__main__":
    main()
