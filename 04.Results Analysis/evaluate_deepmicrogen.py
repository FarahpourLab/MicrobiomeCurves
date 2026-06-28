"""Evaluate imputation results (DeepMicroGen) against the ground truth.

For every masked timepoint the true vs imputed composition is compared with:

* value metrics, computed consistently in two clearly-labelled spaces --
  CLR/Aitchison space and relative-abundance (RA) space (overall and non-zero);
* compositional distances per masked sample -- Aitchison, Bray-Curtis, Jensen-
  Shannon -- and Lin's CCC + NRMSE;
* zero-detection as a classifier -- precision / recall / F1;
* alpha diversity (richness, Pielou evenness, Shannon effective number) and beta
  diversity (Bray-Curtis, Jaccard index) computed on integer counts obtained from
  the relative abundance at a chosen depth.

Counts from abundance (no real count table exists -- only RA / CLR):
CLR is back-transformed to relative abundance, which is then turned into integer
counts at depth ``D`` STOCHASTICALLY -- ``multinomial(D, RA)`` -- i.e. drawing
``D`` reads from the composition, mimicking sequencing sampling noise.
``D`` is the fixed ``DEPTH`` if set, otherwise derived per run as
``round(1 / min_nonzero_RA)`` (the library size at which the rarest observed
taxon equals one read).
"""
import glob
import os
import numpy as np
import pandas as pd
from sklearn.metrics import r2_score
from scipy.stats import pearsonr, spearmanr
from scipy.spatial.distance import braycurtis, jensenshannon, euclidean
from skbio.stats.composition import clr_inv

DEPTH = None                   # fixed sampling depth; None -> derive 1/min_nonzero_RA per run
RAREFY_SEED = 12345            # seed for the stochastic (multinomial) counts
MAX_RAREFY_DEPTH = 100_000     # cap when DEPTH is derived from a tiny min-abundance

# Anchor all paths to the script location so the script works from ANY working
# directory (relative "../" paths only resolve when run from this folder).
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)                             # .../CANDi
BENCH_DIR = os.path.join(PROJECT_ROOT, "04.Benchmarking")             # imputation outputs
DATASETS_DIR = os.path.join(PROJECT_ROOT, "01.Benchmarking Datasets")  # ground truth + masks
RESULTS_DIR = os.path.join(SCRIPT_DIR, "results")                     # where this script writes


# --------------------------------------------------------------------------
# transforms / counts
# --------------------------------------------------------------------------
def clr_to_ra(clr_vals, pseudo_count):
    """Back-transform CLR -> relative abundance: clr_inv, zero out values below
    the pseudo-count, then redistribute the lost mass over the non-zero taxa."""
    ra = np.asarray(clr_inv(np.asarray(clr_vals, dtype=float)), dtype=float)
    ra[ra < pseudo_count] = 0.0
    nz = ra > 0
    n_nz = int(nz.sum())
    if n_nz > 0:
        ra[nz] += (1.0 - ra.sum()) / n_nz
    return ra


def derive_depth(ra_matrix):
    """Implied library size: 1 / smallest non-zero relative abundance, capped."""
    pos = ra_matrix[ra_matrix > 0]
    if pos.size == 0:
        return 0
    return int(min(MAX_RAREFY_DEPTH, round(1.0 / pos.min())))


def ra_to_counts_sto(ra, depth, rng):
    """Stochastic counts at `depth`: multinomial(depth, RA)."""
    p = np.clip(np.asarray(ra, dtype=float), 0, None)
    s = p.sum()
    if s <= 0 or depth <= 0:
        return np.zeros(len(p), dtype=int)
    return rng.multinomial(depth, p / s)


# --------------------------------------------------------------------------
# diversity
# --------------------------------------------------------------------------
def alpha_diversity(counts):
    """(richness, Pielou evenness, Shannon effective number) of a count vector."""
    counts = np.asarray(counts, dtype=float)
    total = counts.sum()
    present = counts > 0
    richness = int(present.sum())
    if total <= 0 or richness == 0:
        return 0, np.nan, np.nan
    p = counts[present] / total
    shannon = float(-(p * np.log(p)).sum())            # natural log -> nats
    evenness = shannon / np.log(richness) if richness > 1 else np.nan   # Pielou's J
    shannon_eff = float(np.exp(shannon))               # Hill number q=1
    return richness, evenness, shannon_eff


def jaccard_index(u, v):
    """Presence/absence Jaccard SIMILARITY index of two count vectors."""
    a = np.asarray(u) > 0
    b = np.asarray(v) > 0
    union = np.logical_or(a, b).sum()
    return float(np.logical_and(a, b).sum() / union) if union else np.nan


def diversity_block(a_ra_list, b_ra_list, to_counts):
    """Aggregate alpha + beta diversity for one count-conversion `to_counts(ra)`."""
    rich_t, even_t, eff_t = [], [], []
    rich_i, even_i, eff_i = [], [], []
    bc, jac = [], []
    for a_ra, b_ra in zip(a_ra_list, b_ra_list):
        ct, ci = to_counts(a_ra), to_counts(b_ra)
        rt, et, ft = alpha_diversity(ct)
        ri, ei, fi = alpha_diversity(ci)
        rich_t.append(rt); even_t.append(et); eff_t.append(ft)
        rich_i.append(ri); even_i.append(ei); eff_i.append(fi)
        bc.append(braycurtis(ct, ci) if (ct.sum() and ci.sum()) else np.nan)
        jac.append(jaccard_index(ct, ci))

    nm = lambda v: float(np.nanmean(v)) if len(v) else np.nan
    df = lambda t, i: float(np.nanmean(np.abs(np.array(t, float) - np.array(i, float))))
    return [nm(rich_t), nm(rich_i), df(rich_t, rich_i),
            nm(even_t), nm(even_i), df(even_t, even_i),
            nm(eff_t), nm(eff_i), df(eff_t, eff_i),
            nm(bc), nm(jac)]


# --------------------------------------------------------------------------
# value-metric helpers (all guarded against degenerate inputs)
# --------------------------------------------------------------------------
def safe_pearson(x, y):
    """Pearson correlation, returning NaN for <2 points or a constant input
    (where the coefficient is undefined) instead of raising/warning."""
    x, y = np.asarray(x), np.asarray(y)
    if len(x) < 2 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return np.nan
    return pearsonr(x, y)[0]


def safe_spearman(x, y):
    """Spearman rank correlation, guarded like ``safe_pearson`` (NaN on <2
    points or a constant input)."""
    x, y = np.asarray(x), np.asarray(y)
    if len(x) < 2 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return np.nan
    return spearmanr(x, y)[0]


def lin_ccc(x, y):
    """Lin's concordance correlation coefficient (agreement, not just linearity)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    denom = x.var() + y.var() + (x.mean() - y.mean()) ** 2
    if denom == 0:
        return np.nan
    return float(2 * np.mean((x - x.mean()) * (y - y.mean())) / denom)


def value_metrics(true_vals, imp_vals):
    """MSE/RMSE/MAE/R2/Pearson/Spearman/CCC/NRMSE over a 1-D set of values."""
    true_vals, imp_vals = np.asarray(true_vals, float), np.asarray(imp_vals, float)
    diffs = true_vals - imp_vals
    mse = float(np.mean(diffs ** 2))
    rmse = float(np.sqrt(mse))
    mae = float(np.mean(np.abs(diffs)))
    rng = np.ptp(true_vals)
    nrmse = float(rmse / rng) if rng > 0 else np.nan
    r2 = r2_score(true_vals, imp_vals) if len(true_vals) > 1 else np.nan
    return [mse, rmse, mae, r2, safe_pearson(true_vals, imp_vals),
            safe_spearman(true_vals, imp_vals), lin_ccc(true_vals, imp_vals), nrmse]


def zero_detection(tp, fp, fn):
    """precision / recall / F1 for predicting the structural zeros."""
    precision = tp / (tp + fp) if (tp + fp) else np.nan
    recall = tp / (tp + fn) if (tp + fn) else np.nan
    if precision and recall and not np.isnan(precision) and not np.isnan(recall):
        f1 = 2 * precision * recall / (precision + recall)
    else:
        f1 = np.nan
    return precision, recall, f1


# --------------------------------------------------------------------------
# per-run aggregation
# --------------------------------------------------------------------------
def calculate_total_metrics(merged_df, samples, pseudo_count):
    """Aggregate every metric for one run over its masked samples.

    ``merged_df`` holds, for each masked sample, the true CLR column and its
    ``_imputed`` counterpart. For each sample the CLR is back-transformed to
    relative abundance, then collected to compute: value metrics in CLR, RA and
    RA-non-zero spaces; per-sample compositional distances (Aitchison, Bray-
    Curtis, Jensen-Shannon); structural-zero detection; and stochastic-count
    alpha/beta diversity at depth ``D``. Returns the flat list of values that
    becomes one result row (order matches ``COLUMNS`` after the key fields)."""
    clr_true, clr_imp = [], []
    ra_true, ra_imp = [], []
    ra_true_nz, ra_imp_nz = [], []
    tp_zero = fp_zero = fn_zero = 0
    per_sample_ait, per_sample_bc, per_sample_jsd = [], [], []
    a_ra_list, b_ra_list = [], []

    for sample in samples:
        a_clr = merged_df[sample].to_numpy(float)
        b_clr = merged_df[sample + "_imputed"].to_numpy(float)
        a_ra = clr_to_ra(a_clr, pseudo_count)
        b_ra = clr_to_ra(b_clr, pseudo_count)

        clr_true.append(a_clr); clr_imp.append(b_clr)
        ra_true.append(a_ra); ra_imp.append(b_ra)
        nz = a_ra > 0
        ra_true_nz.append(a_ra[nz]); ra_imp_nz.append(b_ra[nz])

        tz, pz = (a_ra == 0), (b_ra == 0)
        tp_zero += int(np.logical_and(tz, pz).sum())
        fp_zero += int(np.logical_and(~tz, pz).sum())
        fn_zero += int(np.logical_and(tz, ~pz).sum())

        per_sample_ait.append(euclidean(a_clr, b_clr))
        per_sample_bc.append(braycurtis(a_ra, b_ra))
        if a_ra.sum() > 0 and b_ra.sum() > 0:
            per_sample_jsd.append(jensenshannon(a_ra, b_ra, base=2))
        a_ra_list.append(a_ra); b_ra_list.append(b_ra)

    clr_m = value_metrics(np.concatenate(clr_true), np.concatenate(clr_imp))
    ra_m = value_metrics(np.concatenate(ra_true), np.concatenate(ra_imp))
    ra_nz_m = value_metrics(np.concatenate(ra_true_nz), np.concatenate(ra_imp_nz))
    precision, recall, f1 = zero_detection(tp_zero, fp_zero, fn_zero)

    nanmean = lambda v: float(np.nanmean(v)) if len(v) else np.nan

    # diversity on stochastic counts at depth D (multinomial draw)
    depth = DEPTH if DEPTH else derive_depth(np.vstack(a_ra_list))
    rng = np.random.default_rng(RAREFY_SEED)
    sto_block = diversity_block(a_ra_list, b_ra_list, lambda ra: ra_to_counts_sto(ra, depth, rng))

    return (clr_m + ra_m + ra_nz_m
            + [nanmean(per_sample_ait), nanmean(per_sample_bc), nanmean(per_sample_jsd)]
            + [precision, recall, f1]
            + [depth] + sto_block)


def calc_performance(sel_mask_df, tax_ground_truth_df, tax_imputed_df, sel_pseudo_count):
    """Select the masked samples for one run and score them.

    From the mask matrix it builds the ``subject.timepoint`` labels of the
    masked (``0``) cells, pulls those columns from the ground-truth and imputed
    CLR tables, aligns them, and hands the pair to ``calculate_total_metrics``.
    Subjects with every timepoint masked (nothing observed) are dropped first.
    Returns ``None`` when no scorable masked samples remain (the run is skipped)."""
    # drop subjects whose every timepoint is masked (nothing observed)
    sel_mask_df = sel_mask_df.loc[~(sel_mask_df.drop(columns="subject") == 0).all(axis=1)]

    long = sel_mask_df.melt(id_vars="subject", var_name="timepoint", value_name="value")
    zeros = long[long["value"] == 0].copy()
    zeros["label"] = zeros["subject"] + "." + zeros["timepoint"].astype(str)
    masked_samples = zeros["label"].tolist()

    # keep only masked cells present in BOTH ground-truth and imputed tables
    masked_samples = [s for s in masked_samples
                      if s in tax_ground_truth_df.columns and s in tax_imputed_df.columns]
    if not masked_samples:
        return None

    sel_gt_df = tax_ground_truth_df[masked_samples]
    sel_imp_df = tax_imputed_df[masked_samples]
    merged_df = sel_gt_df.merge(sel_imp_df, left_index=True, right_index=True,
                                suffixes=("", "_imputed"))
    return calculate_total_metrics(merged_df, masked_samples, sel_pseudo_count)


def find_mask_file(sel_level, sel_condition, selected_iter, selected_seed,
                   sel_dataset, missing_method, has_group):
    """Locate the single mask file matching a run and return it with its ratio.

    Globs ``03.processed/<mechanism>/mask/<level>/<ratio>/`` for the mask whose
    name encodes this level, condition (omitted when ``has_group`` is False),
    replicate and seed. Returns ``(mask_df, ratio)`` on a unique match, else
    ``(None, None)`` after warning (0 or >1 matches)."""
    cond = f"_{sel_condition}" if has_group else ""      # BONUS has no group token
    pattern = os.path.join(
        DATASETS_DIR, sel_dataset, "03.processed", missing_method, "mask", sel_level,
        "*", f"masked_metadata_{sel_level}{cond}.rep_{selected_iter}.seed_{selected_seed}.csv")
    matches = glob.glob(pattern)
    if len(matches) != 1:
        print(f"WARNING: expected 1 mask file, found {len(matches)} for pattern:\n  {pattern}")
        return None, None
    return pd.read_csv(matches[0]), os.path.basename(os.path.dirname(matches[0]))


def _value_cols(space):
    """Column names for one value-metric block, tagged with its ``space``
    (e.g. CLR / RA / RA non-zero) so the three blocks stay distinguishable."""
    return [f"MSE ({space})", f"RMSE ({space})", f"MAE ({space})", f"R2 ({space})",
            f"Pearson ({space})", f"Spearman ({space})", f"CCC ({space})", f"NRMSE ({space})"]


def _div_cols():
    """Column names for the diversity block (alpha: richness/evenness/Shannon
    effective as true/imputed/abs-diff; beta: Bray-Curtis, Jaccard index)."""
    return ["Richness true", "Richness imputed", "Richness abs diff",
            "Evenness true", "Evenness imputed", "Evenness abs diff",
            "Shannon eff true", "Shannon eff imputed", "Shannon eff abs diff",
            "Bray-Curtis rarefied", "Jaccard index rarefied"]


SEL_TOOL = "DeepMicroGen"
FIXED_PARAMS = (0.001, 0.7, 3000)      # lr, dr, ep used for the untuned (flat-layout) runs

# Per-dataset layout differences:
#   has_group -- profiles carry a condition/country token (Helminth, DIABIMMUNE) or not (BONUS)
#   layout    -- "tuned": results under <MECH>/<_lr_dr_ep>/ (params parsed from dir name)
#                "flat" : results under <MECH>/<ratio>/      (fixed params, no tuning)
DATASETS = [
    {"name": "Helminth infection and colon cancer", "has_group": True,  "layout": "tuned"},
    {"name": "DIABIMMUNE three country",            "has_group": True,  "layout": "flat"},
    {"name": "BONUS",                               "has_group": False, "layout": "flat"},
]

COLUMNS = (
    ["Taxonomic Level", "Learning Rate", "Dropout Rate", "Number of Epochs",
     "Condition", "Mechanisms of missingness", "Mask Ratio", "Iteration",
     "Pseudo Count", "Seed"]
    + _value_cols("CLR")
    + _value_cols("RA")
    + _value_cols("RA non-zero")
    + ["Aitchison distance (mean)", "Bray-Curtis (mean)", "Jensen-Shannon (mean)"]
    + ["Zero precision", "Zero recall", "Zero F1"]
    + ["Rarefaction depth"]
    + _div_cols()
)


def run_dataset(cfg):
    """Evaluate every DeepMicroGen result for one dataset; return a list of records."""
    sel_dataset = cfg["name"]
    has_group = cfg["has_group"]
    records = []

    for dm_res in glob.glob(os.path.join(BENCH_DIR, sel_dataset, SEL_TOOL, "*", "*", "*.csv")):
        filename = os.path.basename(dm_res).split('.csv')[0]
        if filename.endswith('_scaled'):          # evaluate the raw (unscaled) output only
            continue

        seed_val = filename.split('seed_')[1].split('_')[0]
        rep_num = filename.split('_rep_')[1].split('.')[0]
        after = filename.split('relative_abundance_')[1].split('_')
        tax_level = '_'.join(after[:2])
        if has_group:
            group_val = after[2]
            pseudo_str = after[3]                  # keep exact string for the filename
        else:
            group_val = "all"
            pseudo_str = after[2]
        sel_pseudo_count = float(pseudo_str)

        missing_method = os.path.basename(os.path.dirname(os.path.dirname(dm_res)))
        if cfg["layout"] == "tuned":
            params = [float(v) for v in os.path.basename(os.path.dirname(dm_res)).split('_')[1:]]
        else:
            params = list(FIXED_PARAMS)

        dm_red_df = pd.read_csv(dm_res, index_col=0)
        gt_group = f"{group_val}_" if has_group else ""
        ground_truth_df = pd.read_csv(
            os.path.join(DATASETS_DIR, sel_dataset, "03.processed", "04.sorted",
                         f"sorted_clr_transformation_relative_abundance_"
                         f"{tax_level}_{gt_group}{pseudo_str}.csv"),
            index_col=0)

        print(dm_res)
        mask_df, mask_ratio = find_mask_file(tax_level, group_val, rep_num, seed_val,
                                             sel_dataset, missing_method, has_group)
        if mask_df is None:
            continue

        performance = calc_performance(mask_df, ground_truth_df, dm_red_df, sel_pseudo_count)
        if performance is None:                    # no scorable masked samples
            print(f"  skipped (no masked samples): {os.path.basename(dm_res)}")
            continue

        record = [tax_level, params[0], params[1], params[2],
                  group_val, missing_method, mask_ratio, rep_num, sel_pseudo_count, seed_val]
        record.extend(performance)
        records.append(record)

    return records


if __name__ == '__main__':
    tool_dir = os.path.join(RESULTS_DIR, SEL_TOOL)      # results grouped per tool
    os.makedirs(tool_dir, exist_ok=True)

    for cfg in DATASETS:
        records = run_dataset(cfg)
        df = pd.DataFrame(records, columns=COLUMNS)

        # one CSV per dataset AND per missingness mechanism (MCAR / MAR / MNAR)
        for mechanism, sub in df.groupby("Mechanisms of missingness"):
            out = os.path.join(tool_dir, f'{cfg["name"]}_{SEL_TOOL}_{mechanism}_performance.csv')
            sub.to_csv(out, index=False)
            print(f"wrote {out} ({len(sub)} rows)")
