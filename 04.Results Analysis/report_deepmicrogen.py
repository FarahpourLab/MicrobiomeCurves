"""Summarise and visualise DeepMicroGen evaluation results.

This module also holds the shared reporting machinery (it is imported by
``report_multiway_imputation.py``); running it directly reports DeepMicroGen.

The evaluator writes one CSV per dataset and mechanism under
``results/<tool>/<dataset>_<tool>_<mechanism>_performance.csv``; each row is one
imputation run -- a (parameter-set, condition, ratio, replicate) cell -- with the
full metric panel. ``run(tool)`` loads, per dataset, all mechanisms together and
produces a summary table plus several complementary views. Each view can be
toggled with the ``MAKE_*`` flags below.

Views produced (under results/<tool>/plot/<dataset>/)
-----------------------------------------------------
* **summary table** -- per (mechanism, parameter-set, ratio) mean/std/median of
  every metric (``summary_by_ratio.csv``).
* **per-ratio distributions** (core) -- for each parameter-set and metric, a
  **box** plot and a **violin** plot across mask ratios, plus a **degradation
  curve** (mean +/- 95% CI vs ratio).
* **mechanism overlay** -- the same metric across ratios with MCAR/MAR/MNAR drawn
  together (hue = mechanism): box, violin and line.
* **condition overlay** -- degradation curves with conditions/countries drawn
  together (hue = condition); only when a dataset has >1 condition.
* **diversity scatter** -- true vs imputed richness / evenness / Shannon
  effective, coloured by ratio, against the y = x line.

Tool comparison (hue = tool / win-rate ranking) is intentionally NOT produced
here -- it belongs in a separate cross-tool report.
"""
import glob
import os
import re

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Anchor to the script location so it works from any working directory.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(SCRIPT_DIR, "results")
# plots live alongside each tool's CSVs: results/<tool>/plot/<dataset>/...

# which view families to generate
MAKE_CORE = True                 # box / violin / degradation per parameter-set
MAKE_MECHANISM_OVERLAY = True    # MCAR vs MAR vs MNAR on one axis
MAKE_CONDITION_OVERLAY = True    # conditions/countries on one axis
MAKE_DIVERSITY_SCATTER = True    # true vs imputed alpha diversity

# Restrict which parameter-set gets plotted, per dataset. Only the tuned dataset
# (Helminth + DeepMicroGen) has multiple sets, so we keep just one; untuned
# datasets have a single fixed set and are left unrestricted (no entry here).
# For BAMITA the NN columns are absent, so this is a no-op there.
PARAM_FILTER = {
    "Helminth infection and colon cancer": {
        "Learning Rate": 0.001, "Dropout Rate": 0.3, "Number of Epochs": 3000,
    },
}

# columns that identify a run rather than measure it -- everything else is a metric
KEY_COLS = ["Taxonomic Level", "Learning Rate", "Dropout Rate", "Number of Epochs",
            "Condition", "Mechanisms of missingness", "Mask Ratio", "Iteration",
            "Pseudo Count", "Seed", "Rarefaction depth"]

# a parameter-set = one model/condition configuration; ratios vary WITHIN it.
# Tool-specific columns (e.g. the NN hyperparameters) are absent for BAMITA, so
# the effective set is intersected with each frame's columns via ``param_cols``.
PARAM_COLS = ["Taxonomic Level", "Learning Rate", "Dropout Rate",
              "Number of Epochs", "Condition"]

# (true column, imputed column, label) for the alpha-diversity scatter
DIVERSITY_PAIRS = [
    ("Richness true", "Richness imputed", "Richness"),
    ("Evenness true", "Evenness imputed", "Evenness"),
    ("Shannon eff true", "Shannon eff imputed", "Shannon effective"),
]

MECH_COL = "Mechanisms of missingness"
COND_COL = "Condition"


def param_cols(df):
    """The parameter-set columns actually present in `df` (BAMITA lacks the NN ones)."""
    return [c for c in PARAM_COLS if c in df.columns]


def apply_param_filter(df, dataset):
    """Keep only the chosen parameter-set for tuned datasets; no-op otherwise."""
    spec = PARAM_FILTER.get(dataset)
    if not spec:
        return df
    mask = pd.Series(True, index=df.index)
    for col, val in spec.items():
        if col in df.columns:
            mask &= np.isclose(df[col].astype(float), float(val))
    return df[mask]


# --------------------------------------------------------------------------
# discovery / loading
# --------------------------------------------------------------------------
def safe_name(text):
    """Return a filesystem-safe token (characters awkward in paths replaced)."""
    return re.sub(r"[\\/:*?\"<>|]+", "-", str(text)).strip()


def plot_base(tool, dataset):
    """Per-tool plot directory: results/<tool>/plot/<dataset>/."""
    return os.path.join(RESULTS_DIR, tool, "plot", safe_name(dataset))


def discover_inputs(tool):
    """Group one tool's performance CSVs by dataset.

    Reads ``results/<tool>/*_performance.csv`` and returns
    ``{dataset: [(mechanism, path), ...]}`` so every mechanism of a dataset is
    loaded together (needed for the mechanism-overlay views)."""
    groups = {}
    pattern = os.path.join(RESULTS_DIR, tool, "*_performance.csv")
    for path in sorted(glob.glob(pattern)):
        base = os.path.basename(path)[:-len("_performance.csv")]
        dataset, _tool, mechanism = base.rsplit("_", 2)   # dataset may contain spaces
        groups.setdefault(dataset, []).append((mechanism, path))
    return groups


def load_dataset(items):
    """Concatenate all mechanism CSVs of one dataset into a single tidy frame."""
    frames = []
    for _mechanism, path in items:
        d = pd.read_csv(path)
        d["Mask Ratio"] = pd.to_numeric(d["Mask Ratio"], errors="coerce")
        frames.append(d.dropna(subset=["Mask Ratio"]))
    return pd.concat(frames, ignore_index=True)


def metric_columns(df):
    """All numeric metric columns (everything that is not a key column)."""
    return [c for c in df.columns if c not in KEY_COLS]


def has_data(df, *cols):
    """True if every named column has at least one non-NaN value (plottable)."""
    return all(c in df and df[c].notna().any() for c in cols)


# --------------------------------------------------------------------------
# summary table
# --------------------------------------------------------------------------
def write_summary(df, metrics, out_csv):
    """Per (mechanism, parameter-set, ratio) mean/std/median of each metric."""
    grouped = df.groupby([MECH_COL] + param_cols(df) + ["Mask Ratio"])[metrics]
    summary = grouped.agg(["mean", "std", "median"]).reset_index()
    summary.columns = ["_".join(c).strip("_") if isinstance(c, tuple) else c
                       for c in summary.columns]
    summary.to_csv(out_csv, index=False, encoding="utf-8")


# --------------------------------------------------------------------------
# core per-ratio distributions
# --------------------------------------------------------------------------
def plot_box(df_group, metric, order, palette, title, out_png):
    """Box plot of `metric` across mask ratios, raw points overlaid."""
    plt.figure(figsize=(8, 5))
    sns.boxplot(data=df_group, x="Mask Ratio", y=metric, order=order,
                hue="Mask Ratio", palette=palette, dodge=False,
                showfliers=False, legend=False)
    sns.stripplot(data=df_group, x="Mask Ratio", y=metric, order=order,
                  color="black", jitter=True, size=6, alpha=0.6)
    plt.title(title); plt.xlabel("Mask ratio"); plt.tight_layout()
    plt.savefig(out_png, dpi=150); plt.close()


def plot_violin(df_group, metric, order, palette, title, out_png):
    """Violin plot of `metric` across mask ratios, raw points overlaid."""
    plt.figure(figsize=(8, 5))
    sns.violinplot(data=df_group, x="Mask Ratio", y=metric, order=order,
                   hue="Mask Ratio", palette=palette, dodge=False,
                   cut=0, inner="box", legend=False)
    sns.stripplot(data=df_group, x="Mask Ratio", y=metric, order=order,
                  color="black", jitter=True, size=6, alpha=0.6)
    plt.title(title); plt.xlabel("Mask ratio"); plt.tight_layout()
    plt.savefig(out_png, dpi=150); plt.close()


def plot_degradation(df_group, metric, order, title, out_png, hue=None):
    """Mean +/- 95% CI of `metric` vs mask ratio (optionally split by `hue`)."""
    plt.figure(figsize=(8, 5))
    sns.pointplot(data=df_group, x="Mask Ratio", y=metric, order=order,
                  hue=hue, errorbar=("ci", 95), dodge=0.3 if hue else False)
    plt.title(title); plt.xlabel("Mask ratio"); plt.tight_layout()
    plt.savefig(out_png, dpi=150); plt.close()


# --------------------------------------------------------------------------
# overlay views (several series on one axis)
# --------------------------------------------------------------------------
def plot_overlay_box(df, metric, order, hue, title, out_png):
    """Box plot across ratios with series split by `hue` (mechanism/condition)."""
    plt.figure(figsize=(9, 5))
    sns.boxplot(data=df, x="Mask Ratio", y=metric, order=order,
                hue=hue, showfliers=False)
    plt.title(title); plt.xlabel("Mask ratio")
    plt.legend(title=hue, bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout(); plt.savefig(out_png, dpi=150); plt.close()


def plot_overlay_violin(df, metric, order, hue, title, out_png):
    """Violin plot across ratios with series split by `hue`."""
    plt.figure(figsize=(9, 5))
    sns.violinplot(data=df, x="Mask Ratio", y=metric, order=order,
                   hue=hue, cut=0, inner="quartile")
    plt.title(title); plt.xlabel("Mask ratio")
    plt.legend(title=hue, bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout(); plt.savefig(out_png, dpi=150); plt.close()


# --------------------------------------------------------------------------
# diversity scatter
# --------------------------------------------------------------------------
def plot_diversity_scatter(df_group, title, out_png):
    """True vs imputed alpha diversity (richness/evenness/Shannon eff), coloured
    by mask ratio, against the y = x identity line -- reveals systematic
    over-/under-estimation that the abs-diff columns hide."""
    present = [(t, i, lab) for t, i, lab in DIVERSITY_PAIRS if has_data(df_group, t, i)]
    if not present:
        return
    fig, axes = plt.subplots(1, len(present), figsize=(5 * len(present), 4.5),
                             squeeze=False)
    for ax, (t, i, lab) in zip(axes[0], present):
        sub = df_group[[t, i, "Mask Ratio"]].dropna()
        if sub.empty:
            ax.set_visible(False); continue
        sc = ax.scatter(sub[t], sub[i], c=sub["Mask Ratio"], cmap="viridis",
                        s=80, alpha=0.75, edgecolors="white", linewidths=0.3)
        lo = float(min(sub[t].min(), sub[i].min()))
        hi = float(max(sub[t].max(), sub[i].max()))
        ax.plot([lo, hi], [lo, hi], "k--", lw=1)
        ax.set_xlabel(f"{lab} true"); ax.set_ylabel(f"{lab} imputed"); ax.set_title(lab)
        fig.colorbar(sc, ax=ax, label="Mask ratio")
    fig.suptitle(title)
    fig.tight_layout(rect=[0, 0, 1, 0.90])      # leave room for the 2-line title
    fig.savefig(out_png, dpi=150, bbox_inches="tight"); plt.close(fig)


# --------------------------------------------------------------------------
# orchestration
# --------------------------------------------------------------------------
def core_views(df, dataset, tool, metrics, order, palette):
    """Box / violin / degradation per (mechanism, parameter-set, metric)."""
    pcols = param_cols(df)
    for mechanism, df_mech in df.groupby(MECH_COL):
        base = os.path.join(plot_base(tool, dataset), mechanism)
        for param_values, df_group in df_mech.groupby(pcols):
            param_name = safe_name("_".join(str(v) for v in param_values))
            for metric in metrics:
                if not has_data(df_group, metric):
                    continue
                out_dir = os.path.join(base, param_name, safe_name(metric))
                os.makedirs(out_dir, exist_ok=True)
                title = f"{metric}\n{dataset} | {mechanism} | {param_name}"
                stem = os.path.join(out_dir, safe_name(metric))
                plot_box(df_group, metric, order, palette, title, f"{stem}_box.png")
                plot_violin(df_group, metric, order, palette, title, f"{stem}_violin.png")
                plot_degradation(df_group, metric, order, title, f"{stem}_degradation.png")


def mechanism_overlay_views(df, dataset, tool, metrics, order):
    """MCAR/MAR/MNAR drawn together per (parameter-set, metric)."""
    if df[MECH_COL].nunique() < 2:
        return
    for param_values, df_group in df.groupby(param_cols(df)):
        param_name = safe_name("_".join(str(v) for v in param_values))
        base = os.path.join(plot_base(tool, dataset), "mechanism_comparison", param_name)
        for metric in metrics:
            if not has_data(df_group, metric):
                continue
            out_dir = os.path.join(base, safe_name(metric))
            os.makedirs(out_dir, exist_ok=True)
            title = f"{metric}\n{dataset} | mechanisms | {param_name}"
            stem = os.path.join(out_dir, safe_name(metric))
            plot_overlay_box(df_group, metric, order, MECH_COL, title, f"{stem}_box.png")
            plot_overlay_violin(df_group, metric, order, MECH_COL, title, f"{stem}_violin.png")
            plot_degradation(df_group, metric, order, title, f"{stem}_line.png", hue=MECH_COL)


def condition_views(df, dataset, tool, metrics, order):
    """Per mechanism: condition-overlay degradation curves (>1 condition only)."""
    if df[COND_COL].nunique() < 2:
        return
    for mechanism, df_mech in df.groupby(MECH_COL):
        base = os.path.join(plot_base(tool, dataset), mechanism)
        for metric in metrics:
            if not has_data(df_mech, metric) or df_mech[COND_COL].nunique() < 2:
                continue
            title = f"{metric}\n{dataset} | {mechanism}"
            d = os.path.join(base, "condition_comparison"); os.makedirs(d, exist_ok=True)
            plot_degradation(df_mech, metric, order, title,
                             os.path.join(d, f"{safe_name(metric)}_line.png"), hue=COND_COL)


def diversity_scatter_views(df, dataset, tool):
    """True-vs-imputed alpha-diversity scatter per (mechanism, parameter-set)."""
    pcols = param_cols(df)
    for mechanism, df_mech in df.groupby(MECH_COL):
        base = os.path.join(plot_base(tool, dataset), mechanism)
        for param_values, df_group in df_mech.groupby(pcols):
            param_name = safe_name("_".join(str(v) for v in param_values))
            out_dir = os.path.join(base, param_name); os.makedirs(out_dir, exist_ok=True)
            title = f"Alpha diversity true vs imputed\n{dataset} | {mechanism} | {param_name}"
            plot_diversity_scatter(df_group, title,
                                   os.path.join(out_dir, "diversity_scatter.png"))


def process_dataset(dataset, tool, items):
    """Build the summary table and every enabled view for one (dataset, tool)."""
    df = load_dataset(items)
    df = apply_param_filter(df, dataset)        # tuned datasets -> single parameter-set
    if df.empty:
        print(f"{dataset} [{tool}]: no rows after parameter filter -- skipped")
        return
    metrics = metric_columns(df)
    order = sorted(df["Mask Ratio"].unique())
    palette = dict(zip(order, sns.color_palette("viridis", len(order))))

    base = plot_base(tool, dataset)
    os.makedirs(base, exist_ok=True)
    write_summary(df, metrics, os.path.join(base, "summary_by_ratio.csv"))

    if MAKE_CORE:
        core_views(df, dataset, tool, metrics, order, palette)
    if MAKE_MECHANISM_OVERLAY:
        mechanism_overlay_views(df, dataset, tool, metrics, order)
    if MAKE_CONDITION_OVERLAY:
        condition_views(df, dataset, tool, metrics, order)
    if MAKE_DIVERSITY_SCATTER:
        diversity_scatter_views(df, dataset, tool)

    print(f"{dataset} [{tool}]: {len(metrics)} metrics, "
          f"{df[MECH_COL].nunique()} mechanism(s), "
          f"{df.groupby(param_cols(df)).ngroups} parameter-set(s)")


def run(tool):
    """Generate the summary tables and plots for every dataset of one `tool`."""
    os.makedirs(RESULTS_DIR, exist_ok=True)
    groups = discover_inputs(tool)
    if not groups:
        raise SystemExit(f"no {tool} *_performance.csv found in "
                         f"{os.path.join(RESULTS_DIR, tool)}/")
    for dataset, items in groups.items():
        process_dataset(dataset, tool, items)
    print(f"done [{tool}]: {len(groups)} dataset(s) -> results/{tool}/plot/")


if __name__ == "__main__":
    run("DeepMicroGen")
