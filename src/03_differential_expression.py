from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests


def ensure_dirs(*paths: Path) -> None:
    for p in paths:
        p.mkdir(parents=True, exist_ok=True)


def load_pheno(path: Path) -> pd.DataFrame:
    ph = pd.read_csv(path, sep="\t")
    # robust normalize headers
    ph.columns = [str(c).replace("\ufeff", "").strip() for c in ph.columns]
    if "sample_id" not in ph.columns:
        # fallback: if index was saved oddly
        if "index" in ph.columns:
            ph = ph.rename(columns={"index": "sample_id"})
    if "group" not in ph.columns:
        raise ValueError(f"'group' column missing in {path}. Columns: {list(ph.columns)}")
    ph["sample_id"] = ph["sample_id"].astype(str).str.strip()
    ph["group"] = ph["group"].astype(str).str.strip()
    return ph[["sample_id", "group"]].copy()


def try_load_platform(platform_path: Path) -> Optional[pd.DataFrame]:
    if not platform_path.exists():
        return None
    plat = pd.read_csv(platform_path, sep="\t")
    plat.columns = [str(c).replace("\ufeff", "").strip() for c in plat.columns]
    # Most GPL tables have "ID" as probe id column; keep everything and decide later.
    return plat


def guess_probe_and_gene_columns(platform_df: pd.DataFrame) -> Tuple[str, Optional[str]]:
    """
    Try to infer probe id column and gene symbol column from platform table.
    Return (probe_col, gene_col_or_None)
    """
    cols = platform_df.columns.tolist()
    lower = {c.lower(): c for c in cols}

    # probe id
    probe_col = None
    for key in ["id", "id_ref", "probe", "probe_id", "identifier"]:
        if key in lower:
            probe_col = lower[key]
            break
    if probe_col is None:
        probe_col = cols[0]  # fallback

    # gene symbol
    gene_col = None
    for key in ["gene symbol", "gene_symbol", "symbol", "genesymbol", "gene"]:
        if key in lower:
            gene_col = lower[key]
            break

    # Another common GEO platform column name is "Gene Symbol" with space
    if gene_col is None:
        for c in cols:
            if "gene" in c.lower() and "symbol" in c.lower():
                gene_col = c
                break

    return probe_col, gene_col


def compute_contrast_de(expr: pd.DataFrame, ph: pd.DataFrame, group_a: str, group_b: str) -> pd.DataFrame:
    """
    For each feature/probe:
      y ~ 1 + I(group==A)   using samples from A and B only
    Returns a table with logFC (meanA-meanB), pval, adj_pval.
    """
    samples = ph[ph["group"].isin([group_a, group_b])].copy()
    samples = samples.set_index("sample_id")

    # Align expression columns to phenotype
    expr = expr.copy()
    expr.columns = [str(c) for c in expr.columns]
    samples = samples.loc[[c for c in expr.columns if c in samples.index]]

    a_mask = (samples["group"] == group_a).astype(int).values
    X = sm.add_constant(a_mask)  # intercept + group indicator

    # Precompute means for logFC
    a_samples = samples.index[samples["group"] == group_a].tolist()
    b_samples = samples.index[samples["group"] == group_b].tolist()

    mean_a = expr[a_samples].mean(axis=1, skipna=True)
    mean_b = expr[b_samples].mean(axis=1, skipna=True)
    logfc = (mean_a - mean_b)  # on log-scale if data already log; else "difference in expression"

    # Fit per feature
    pvals = np.empty(expr.shape[0], dtype=float)
    betas = np.empty(expr.shape[0], dtype=float)

    Y = expr[samples.index].apply(pd.to_numeric, errors="coerce").values  # features x samples
    for i in range(Y.shape[0]):
        y = Y[i, :]
        # handle missing
        ok = np.isfinite(y)
        if ok.sum() < 4:  # too few points
            pvals[i] = np.nan
            betas[i] = np.nan
            continue
        model = sm.OLS(y[ok], X[ok, :])
        res = model.fit()
        betas[i] = res.params[1] if len(res.params) > 1 else np.nan
        pvals[i] = res.pvalues[1] if len(res.pvalues) > 1 else np.nan

    df = pd.DataFrame({
        "feature_id": expr.index.astype(str),
        "group_a": group_a,
        "group_b": group_b,
        "logFC": logfc.values,
        "beta": betas,
        "pval": pvals,
    })

    # FDR
    valid = df["pval"].notna().values
    df["padj"] = np.nan
    if valid.sum() > 0:
        df.loc[valid, "padj"] = multipletests(df.loc[valid, "pval"].values, method="fdr_bh")[1]

    # add ranking helper
    df["neglog10_padj"] = -np.log10(df["padj"].astype(float))
    df.replace([np.inf, -np.inf], np.nan, inplace=True)

    # sort
    df = df.sort_values(["padj", "pval"], ascending=[True, True])
    return df


def volcano_plot(de: pd.DataFrame, out: Path, title: str, fc_thr: float = 1.0, fdr_thr: float = 0.05) -> None:
    """
    Volcano plot: logFC vs -log10(padj)
    """
    x = de["logFC"].astype(float)
    y = de["neglog10_padj"].astype(float)

    sig = (de["padj"].astype(float) < fdr_thr) & (x.abs() >= fc_thr)

    plt.figure(figsize=(6.5, 5))
    plt.scatter(x[~sig], y[~sig], s=8, alpha=0.6)
    plt.scatter(x[sig], y[sig], s=10, alpha=0.9)

    plt.axvline(fc_thr, linestyle="--")
    plt.axvline(-fc_thr, linestyle="--")
    plt.axhline(-np.log10(fdr_thr), linestyle="--")

    plt.xlabel("logFC (mean A - mean B)")
    plt.ylabel("-log10(FDR)")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out, dpi=220)
    plt.close()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--expr", default="data/interim/expression_qc.tsv")
    ap.add_argument("--pheno", default="data/raw/phenotype.tsv")
    ap.add_argument("--platform", default="data/raw/platform.tsv")
    args = ap.parse_args()

    expr_path = Path(args.expr)
    pheno_path = Path(args.pheno)
    platform_path = Path(args.platform)

    results_fig = Path("results/figures")
    results_tab = Path("results/tables")
    ensure_dirs(results_fig, results_tab)

    print(f">> Loading expr: {expr_path}")
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    print(f">> expr shape: {expr.shape} (features x samples)")

    print(f">> Loading pheno: {pheno_path}")
    ph = load_pheno(pheno_path)
    print(">> group counts:")
    print(ph["group"].value_counts().to_string())

    platform = try_load_platform(platform_path)
    probe_col, gene_col = (None, None)
    probe_to_gene = None

    if platform is not None:
        probe_col, gene_col = guess_probe_and_gene_columns(platform)
        print(f">> platform loaded: {platform.shape}, probe_col='{probe_col}', gene_col='{gene_col}'")
        if gene_col is not None:
            probe_to_gene = (
                platform[[probe_col, gene_col]]
                .dropna()
                .drop_duplicates(subset=[probe_col])
                .rename(columns={probe_col: "feature_id", gene_col: "gene_symbol"})
            )
        else:
            print(">> No gene symbol column detected in platform.tsv; will keep feature_id only.")

    contrasts = [
        ("active_fvm", "control_retina"),
        ("inactive_fvm", "control_retina"),
        ("active_fvm", "inactive_fvm"),
    ]

    for a, b in contrasts:
        name = f"{a}_vs_{b}"
        print(f"\n>> Running DE: {name}")
        de = compute_contrast_de(expr, ph, a, b)

        if probe_to_gene is not None:
            de = de.merge(probe_to_gene, on="feature_id", how="left")

        out_tsv = results_tab / f"deg_{name}.tsv"
        de.to_csv(out_tsv, sep="\t", index=False)
        print(f">> Saved: {out_tsv}  (rows={len(de)})")

        out_png = results_fig / f"volcano_{name}.png"
        volcano_plot(de, out_png, title=f"Volcano: {name}")
        print(f">> Saved: {out_png}")

    print("\n✅ Step 3 complete. Next: Enrichment (src/04_enrichment.py).")


if __name__ == "__main__":
    main()
