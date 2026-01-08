from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gseapy as gp


# Libraries commonly available in Enrichr that are useful for drug repurposing.
# We will also "discover" which are available at runtime and skip if missing.
DRUG_LIBRARIES = [
    "DSigDB",
    "Drug_Perturbations_from_GEO_up",
    "Drug_Perturbations_from_GEO_down",
    "LINCS_L1000_Chem_Pert_up",
    "LINCS_L1000_Chem_Pert_down",
]


def ensure_dirs(*paths: Path) -> None:
    for p in paths:
        p.mkdir(parents=True, exist_ok=True)


def collapse_probes_to_genes(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    d["gene_symbol"] = d["gene_symbol"].astype(str).str.strip()
    d = d[d["gene_symbol"].notna()]
    d = d[d["gene_symbol"] != ""]
    d = d[d["gene_symbol"].str.lower() != "nan"]

    d["padj_num"] = pd.to_numeric(d.get("padj", pd.Series([np.nan] * len(d))), errors="coerce")
    d["logFC_num"] = pd.to_numeric(d.get("logFC", pd.Series([np.nan] * len(d))), errors="coerce")
    d["abs_logFC"] = d["logFC_num"].abs()

    d = d.sort_values(["padj_num", "abs_logFC"], ascending=[True, False])
    d = d.drop_duplicates(subset=["gene_symbol"], keep="first")
    d = d.drop(columns=["padj_num", "abs_logFC"], errors="ignore")
    return d


def split_up_down(df_gene: pd.DataFrame, fdr_thr: float, fc_thr: float) -> Tuple[List[str], List[str]]:
    d = df_gene.copy()
    d["padj"] = pd.to_numeric(d["padj"], errors="coerce")
    d["logFC"] = pd.to_numeric(d["logFC"], errors="coerce")

    sig = d[(d["padj"] < fdr_thr) & (d["logFC"].abs() >= fc_thr)].copy()
    up = sig[sig["logFC"] > 0]["gene_symbol"].dropna().astype(str).tolist()
    down = sig[sig["logFC"] < 0]["gene_symbol"].dropna().astype(str).tolist()

    def dedup(x: List[str]) -> List[str]:
        seen = set()
        out = []
        for g in x:
            if g not in seen:
                out.append(g)
                seen.add(g)
        return out

    return dedup(up), dedup(down)


def enrichr_results(genes: List[str], library: str) -> Optional[pd.DataFrame]:
    if len(genes) < 5:
        return None
    enr = gp.enrichr(
        gene_list=genes,
        gene_sets=[library],
        organism="Human",
        outdir=None,
        no_plot=True,
    )
    res = enr.results.copy()
    res.columns = [str(c).strip().replace(" ", "_") for c in res.columns]
    return res


def get_term_col(res: pd.DataFrame) -> str:
    return "Term" if "Term" in res.columns else res.columns[0]


def get_fdr_col(res: pd.DataFrame) -> str:
    # gseapy uses Adjusted_P-value typically
    if "Adjusted_P-value" in res.columns:
        return "Adjusted_P-value"
    if "Adjusted_P_value" in res.columns:
        return "Adjusted_P_value"
    if "FDR" in res.columns:
        return "FDR"
    # fallback
    if "P-value" in res.columns:
        return "P-value"
    return res.columns[1] if len(res.columns) > 1 else res.columns[0]


def safe_neglog10(x: pd.Series) -> pd.Series:
    v = pd.to_numeric(x, errors="coerce")
    v = v.replace([0, np.inf, -np.inf], np.nan)
    return -np.log10(v)


def compute_reversal_table(res_up: Optional[pd.DataFrame], res_down: Optional[pd.DataFrame]) -> pd.DataFrame:
    """
    Reversal Score idea:
      score = (-log10(FDR_down)) - (-log10(FDR_up))
    Higher score => more enriched in DOWN than UP (potential reversal).
    """
    if res_up is None:
        res_up = pd.DataFrame()
    if res_down is None:
        res_down = pd.DataFrame()

    if len(res_up) == 0 and len(res_down) == 0:
        return pd.DataFrame()

    # normalize
    def prep(res: pd.DataFrame, side: str) -> pd.DataFrame:
        if len(res) == 0:
            return pd.DataFrame(columns=["Term", f"FDR_{side}", f"score_{side}"])
        term_col = get_term_col(res)
        fdr_col = get_fdr_col(res)
        out = res.copy()
        out = out.rename(columns={term_col: "Term", fdr_col: f"FDR_{side}"})
        out[f"score_{side}"] = safe_neglog10(out[f"FDR_{side}"])
        return out[["Term", f"FDR_{side}", f"score_{side}"]]

    up = prep(res_up, "UP")
    dn = prep(res_down, "DOWN")

    merged = pd.merge(dn, up, on="Term", how="outer")
    merged["score_DOWN"] = merged["score_DOWN"].fillna(0.0)
    merged["score_UP"] = merged["score_UP"].fillna(0.0)
    merged["ReversalScore"] = merged["score_DOWN"] - merged["score_UP"]

    # also keep min FDR as a sanity indicator
    merged["minFDR"] = pd.to_numeric(
        merged[["FDR_DOWN", "FDR_UP"]].min(axis=1),
        errors="coerce"
    )

    merged = merged.sort_values(["ReversalScore", "score_DOWN", "minFDR"], ascending=[False, False, True])
    return merged


def plot_top_drugs(df: pd.DataFrame, out_png: Path, title: str, top_n: int = 20) -> None:
    if df is None or len(df) == 0:
        return
    d = df.head(top_n).copy()
    d = d.iloc[::-1]

    plt.figure(figsize=(10.5, 6))
    plt.barh(d["Term"].astype(str), pd.to_numeric(d["ReversalScore"], errors="coerce"))
    plt.xlabel("ReversalScore  (higher = more DOWN-enriched than UP-enriched)")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=220)
    plt.close()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--deg_dir", default="results/tables")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--fc", type=float, default=1.0)
    ap.add_argument("--top_n", type=int, default=25)
    args = ap.parse_args()

    deg_dir = Path(args.deg_dir)
    out_tables = Path("results/tables")
    out_figs = Path("results/figures")
    ensure_dirs(out_tables, out_figs)

    deg_files = sorted(deg_dir.glob("deg_*.tsv"))
    if not deg_files:
        raise FileNotFoundError(f"No deg_*.tsv files found in {deg_dir}")

    # Try to fetch available enrichr libs and filter our list.
    try:
        avail = set(gp.get_library_name(organism="Human"))
        libs = [l for l in DRUG_LIBRARIES if l in avail]
        missing = [l for l in DRUG_LIBRARIES if l not in avail]
        print(">> Enrichr drug libraries available:")
        for l in libs:
            print("  -", l)
        if missing:
            print(">> Skipping missing libraries (not present in Enrichr for Human):")
            for l in missing:
                print("  -", l)
    except Exception as e:
        print(f">> Warning: could not fetch Enrichr library list ({e}). Will try libraries anyway.")
        libs = DRUG_LIBRARIES

    for f in deg_files:
        contrast = f.stem.replace("deg_", "")
        print(f"\n>> Processing contrast: {contrast}")
        df = pd.read_csv(f, sep="\t")

        if "gene_symbol" not in df.columns:
            print(f"   - Skip: no gene_symbol column in {f.name}")
            continue

        df_gene = collapse_probes_to_genes(df)
        up, down = split_up_down(df_gene, fdr_thr=args.fdr, fc_thr=args.fc)

        print(f"   - Significant UP: {len(up)}  DOWN: {len(down)}  (FDR<{args.fdr}, |logFC|>={args.fc})")
        if len(up) < 5 and len(down) < 5:
            print("   - Skip: too few genes for drug repurposing.")
            continue

        for lib in libs:
            print(f"   >> Library: {lib}")
            res_up = enrichr_results(up, lib)
            res_dn = enrichr_results(down, lib)

            # save raw tables too (optional but useful)
            if res_up is not None:
                raw_up = out_tables / f"drugrep_raw_{contrast}_{lib}_UP.tsv"
                res_up.to_csv(raw_up, sep="\t", index=False)
            if res_dn is not None:
                raw_dn = out_tables / f"drugrep_raw_{contrast}_{lib}_DOWN.tsv"
                res_dn.to_csv(raw_dn, sep="\t", index=False)

            merged = compute_reversal_table(res_up, res_dn)
            if len(merged) == 0:
                print("      - No results.")
                continue

            out_all = out_tables / f"drugrep_{contrast}_{lib}.tsv"
            merged.to_csv(out_all, sep="\t", index=False)
            print(f"      - Saved: {out_all} (n_terms={len(merged)})")

            out_top = out_tables / f"drugrep_{contrast}_{lib}_TOP.tsv"
            merged.head(args.top_n).to_csv(out_top, sep="\t", index=False)
            print(f"      - Saved: {out_top} (top_n={args.top_n})")

            out_png = out_figs / f"drugrep_{contrast}_{lib}_top.png"
            plot_top_drugs(merged, out_png, title=f"{contrast} | {lib} | Top candidates", top_n=min(args.top_n, 25))
            print(f"      - Saved: {out_png}")

    print("\n✅ Step 5 complete. Next: Final report + README polishing.")


if __name__ == "__main__":
    main()
