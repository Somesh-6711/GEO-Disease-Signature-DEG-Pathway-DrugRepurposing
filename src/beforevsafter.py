from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

RAW_EXPR = Path("data/raw/expression.tsv")
QC_EXPR = Path("data/interim/expression_qc.tsv")
PHENO = Path("data/raw/phenotype.tsv")
OUT = Path("results/figures")
OUT.mkdir(parents=True, exist_ok=True)

def load_expr(p: Path) -> pd.DataFrame:
    df = pd.read_csv(p, sep="\t", index_col=0)
    return df

def log2_transform(expr: pd.DataFrame) -> pd.DataFrame:
    x = expr.copy()
    # avoid log on negatives if present (microarray often already log)
    # if any negative values exist, shift minimally
    minv = np.nanmin(x.values)
    if minv <= 0:
        x = x - minv + 1.0
    return np.log2(x + 1.0)

def boxplot(expr: pd.DataFrame, out_png: Path, title: str):
    plt.figure(figsize=(12, 4))
    plt.boxplot([expr[c].dropna().values for c in expr.columns], tick_labels=expr.columns, showfliers=False)
    plt.xticks(rotation=45, ha="right")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=220)
    plt.close()

def pca_plot(expr: pd.DataFrame, out_png: Path, title: str):
    ph = pd.read_csv(PHENO, sep="\t")
    # infer sample_id column
    if "sample_id" not in ph.columns:
        ph = ph.rename(columns={ph.columns[0]: "sample_id"})
    if "group" not in ph.columns and len(ph.columns) > 1:
        ph = ph.rename(columns={ph.columns[1]: "group"})

    X = expr.T.values
    X = np.nan_to_num(X, nan=np.nanmedian(X))
    X = StandardScaler().fit_transform(X)
    pcs = PCA(n_components=2).fit_transform(X)

    dfp = pd.DataFrame({"PC1": pcs[:,0], "PC2": pcs[:,1], "sample_id": expr.columns})
    dfp = dfp.merge(ph[["sample_id","group"]], on="sample_id", how="left")

    plt.figure(figsize=(6.5, 5))
    for g in sorted(dfp["group"].dropna().unique()):
        sub = dfp[dfp["group"] == g]
        plt.scatter(sub["PC1"], sub["PC2"], label=str(g))
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=220)
    plt.close()

def ma_plot(deg_path: Path, out_png: Path, title: str):
    df = pd.read_csv(deg_path, sep="\t")
    # expected columns: logFC, avgExpr (if exists), padj
    if "avgExpr" in df.columns:
        A = pd.to_numeric(df["avgExpr"], errors="coerce")
    else:
        # fallback: use abs(logFC) as proxy if avgExpr not present
        A = pd.Series(np.nan, index=df.index)

    M = pd.to_numeric(df["logFC"], errors="coerce")
    padj = pd.to_numeric(df.get("padj", np.nan), errors="coerce")

    plt.figure(figsize=(7, 5))
    plt.scatter(A, M, s=6, alpha=0.35)
    if A.notna().any():
        sig = (padj < 0.05) & (M.abs() >= 1.0)
        plt.scatter(A[sig], M[sig], s=7, alpha=0.8)
    plt.axhline(0, linewidth=1)
    plt.xlabel("A (avg expression)")
    plt.ylabel("M (logFC)")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=220)
    plt.close()

def main():
    raw = load_expr(RAW_EXPR)
    qc = load_expr(QC_EXPR)
    raw_log = log2_transform(raw)
    qc_log = log2_transform(qc)

    boxplot(raw, OUT / "before_after_boxplot_raw.png", "Before: Raw expression distribution")
    boxplot(raw_log, OUT / "before_after_boxplot_log2.png", "After: log2(x+1) expression distribution")

    pca_plot(raw_log, OUT / "before_after_pca_rawlog.png", "Before-ish: PCA on log2(raw)")
    pca_plot(qc_log, OUT / "before_after_pca_qclog.png", "After: PCA on log2(QC-filtered)")

    # MA plot for one contrast if available
    deg = Path("results/tables/deg_active_fvm_vs_control_retina.tsv")
    if deg.exists():
        ma_plot(deg, OUT / "ma_active_fvm_vs_control_retina.png", "MA Plot: active_fvm_vs_control_retina")

    print(">> Saved before/after plots to results/figures/")

if __name__ == "__main__":
    main()
