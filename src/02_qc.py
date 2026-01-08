from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA


def ensure_dirs(*paths: Path) -> None:
    for p in paths:
        p.mkdir(parents=True, exist_ok=True)


def strip_bom(s: str) -> str:
    # Handle UTF-8 BOM in headers if present
    return s.replace("\ufeff", "").strip()


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [strip_bom(str(c)) for c in df.columns]
    return df


def ensure_sample_group_columns(pheno_raw: pd.DataFrame) -> pd.DataFrame:
    """
    Return phenotype with guaranteed columns: sample_id, group.
    Handles cases where sample_id is:
      - a column with weird casing/spaces/BOM
      - the index
      - stored as 'index' after reset_index
    """
    ph = normalize_columns(pheno_raw)

    # If sample id might be in index, bring it out
    if ph.index.name:
        idx_name = strip_bom(str(ph.index.name)).lower()
        if ("sample" in idx_name and "id" in idx_name) or idx_name in ("gsm", "geo_accession"):
            ph = ph.reset_index()

    ph = normalize_columns(ph)

    # Build mapping of lowercase->original
    lower_map = {strip_bom(c).lower(): c for c in ph.columns}

    # Candidate keys
    sample_keys = ["sample_id", "geo_accession", "gsm", "sampleid", "sample"]
    group_keys = ["group", "condition", "status", "phenotype", "class", "label"]

    sample_col = None
    for k in sample_keys:
        if k in lower_map:
            sample_col = lower_map[k]
            break

    # sometimes reset_index creates 'index'
    if sample_col is None and "index" in lower_map:
        sample_col = lower_map["index"]

    group_col = None
    for k in group_keys:
        if k in lower_map:
            group_col = lower_map[k]
            break

    if sample_col is None or group_col is None:
        raise ValueError(
            "Could not identify sample/group columns in phenotype.tsv.\n"
            f"Columns found: {list(ph.columns)}\n"
            "Expected something like sample_id + group (any casing/spaces)."
        )

    ph = ph.rename(columns={sample_col: "sample_id", group_col: "group"}).copy()
    ph["sample_id"] = ph["sample_id"].astype(str).str.strip()
    ph["group"] = ph["group"].astype(str).str.strip()

    return ph[["sample_id", "group"]].copy()


def log1p_if_needed(expr: pd.DataFrame) -> tuple[pd.DataFrame, str]:
    vals = expr.to_numpy().ravel()
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return expr, "none"
    q95 = float(np.quantile(vals, 0.95))
    if q95 > 100:
        return np.log1p(expr), "log1p"
    return expr, "none"


def plot_distributions(expr: pd.DataFrame, out: Path, title: str) -> None:
    plt.figure(figsize=(10, 5))
    data = [expr[c].dropna().values for c in expr.columns]
    plt.boxplot(data, tick_labels=list(expr.columns), showfliers=False)
    plt.xticks(rotation=90)
    plt.ylabel("Expression")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()


def plot_missingness(expr: pd.DataFrame, out: Path) -> pd.DataFrame:
    miss = pd.DataFrame({
        "sample_id": list(expr.columns),
        "missing_fraction": expr.isna().mean(axis=0).values,
    }).sort_values("missing_fraction", ascending=False)

    plt.figure(figsize=(8, 4))
    plt.bar(miss["sample_id"], miss["missing_fraction"])
    plt.xticks(rotation=90)
    plt.ylabel("Missing fraction")
    plt.title("Per-sample missingness")
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return miss


def plot_sample_correlation(expr: pd.DataFrame, out: Path) -> None:
    corr = expr.corr(method="pearson", min_periods=100)

    plt.figure(figsize=(6.8, 6))
    plt.imshow(corr.values, aspect="auto")
    plt.colorbar(label="Pearson r")
    plt.xticks(range(len(corr.columns)), corr.columns, rotation=90, fontsize=7)
    plt.yticks(range(len(corr.index)), corr.index, fontsize=7)
    plt.title("Sample–sample correlation")
    plt.tight_layout()
    plt.savefig(out, dpi=220)
    plt.close()


def run_pca(expr: pd.DataFrame, pheno: pd.DataFrame, out: Path) -> pd.DataFrame:
    """
    PCA on samples with safe merge to phenotype.
    expr: features x samples
    pheno: must have columns sample_id, group (we enforce upstream)
    """
    # samples x features
    X = expr.T.copy()
    X = X.apply(pd.to_numeric, errors="coerce")
    X = X.fillna(X.median(axis=0))

    pca = PCA(n_components=2, random_state=0)
    pcs = pca.fit_transform(X.values)

    pca_df = pd.DataFrame({
        "sample_id": X.index.astype(str),
        "PC1": pcs[:, 0],
        "PC2": pcs[:, 1],
    })

    ph = pheno.copy()
    # If sample_id somehow became index again, fix it:
    if "sample_id" not in ph.columns:
        if ph.index.name:
            ph = ph.reset_index()
        if "index" in ph.columns and "sample_id" not in ph.columns:
            ph = ph.rename(columns={"index": "sample_id"})
    if "sample_id" not in ph.columns:
        raise RuntimeError(f"Phenotype still missing 'sample_id'. Columns now: {list(ph.columns)}")

    ph["sample_id"] = ph["sample_id"].astype(str)
    pca_df = pca_df.merge(ph[["sample_id", "group"]], on="sample_id", how="left")

    plt.figure(figsize=(6.5, 5))
    for g, sub in pca_df.groupby("group", dropna=False):
        plt.scatter(sub["PC1"], sub["PC2"], label=str(g), s=60)

    for _, r in pca_df.iterrows():
        plt.text(r["PC1"], r["PC2"], r["sample_id"], fontsize=7, alpha=0.85)

    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    plt.title("PCA of samples (QC)")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(out, dpi=220)
    plt.close()

    return pca_df


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--expr", default="data/raw/expression.tsv")
    ap.add_argument("--pheno", default="data/raw/phenotype.tsv")
    ap.add_argument("--min_nonmissing_frac", type=float, default=0.90)
    args = ap.parse_args()

    expr_path = Path(args.expr)
    pheno_path = Path(args.pheno)

    results_fig = Path("results/figures")
    results_tab = Path("results/tables")
    interim = Path("data/interim")
    ensure_dirs(results_fig, results_tab, interim)

    print(f">> Loading: {expr_path}")
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    expr.columns = [str(c) for c in expr.columns]
    print(f">> Expression shape: {expr.shape} (features x samples)")

    print(f">> Loading: {pheno_path}")
    pheno_raw = pd.read_csv(pheno_path, sep="\t")
    pheno = ensure_sample_group_columns(pheno_raw)

    print(f">> Phenotype columns used: {list(pheno.columns)}")
    # Align phenotype to expression columns
    missing_in_expr = sorted(set(pheno["sample_id"]) - set(expr.columns))
    if missing_in_expr:
        raise RuntimeError(f"These phenotype samples are missing in expression.tsv columns: {missing_in_expr}")

    # reorder pheno to match expr column order
    pheno = pheno.set_index("sample_id").loc[expr.columns].reset_index()
    # sometimes reset_index gives 'index' if index name gets lost; fix it:
    if "sample_id" not in pheno.columns and "index" in pheno.columns:
        pheno = pheno.rename(columns={"index": "sample_id"})

    # QC plots
    plot_distributions(expr, results_fig / "qc_boxplot_raw.png", "Raw expression distribution (boxplot)")
    miss = plot_missingness(expr, results_fig / "qc_missingness.png")

    expr2, transform = log1p_if_needed(expr)
    print(f">> Transform applied: {transform}")
    if transform != "none":
        plot_distributions(expr2, results_fig / "qc_boxplot_transformed.png",
                           f"Transformed expression distribution ({transform})")

    nonmiss_frac = expr2.notna().mean(axis=1)
    keep = nonmiss_frac >= args.min_nonmissing_frac
    expr_filt = expr2.loc[keep].copy()

    print(f">> Feature filter: kept {expr_filt.shape[0]} / {expr2.shape[0]} features "
          f"(min_nonmissing_frac={args.min_nonmissing_frac})")

    plot_sample_correlation(expr_filt, results_fig / "qc_sample_correlation.png")
    pca_df = run_pca(expr_filt, pheno, results_fig / "qc_pca.png")
    pca_df.to_csv(results_tab / "qc_pca_scores.tsv", sep="\t", index=False)

    expr_out = interim / "expression_qc.tsv"
    expr_filt.to_csv(expr_out, sep="\t", index=True)

    summary = miss.merge(pheno, on="sample_id", how="left").merge(
        pca_df[["sample_id", "PC1", "PC2"]],
        on="sample_id",
        how="left",
    )
    summary.to_csv(results_tab / "qc_summary.tsv", sep="\t", index=False)

    print(f">> Saved: {expr_out}")
    print(f">> Saved: {results_tab / 'qc_summary.tsv'}")
    print(f">> Figures saved in: {results_fig}")
    print("\n✅ Step 2 complete. Next: Differential Expression (src/03_differential_expression.py).")


if __name__ == "__main__":
    main()
