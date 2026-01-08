from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gseapy as gp


LIBRARIES = [
    "GO_Biological_Process_2021",
    "GO_Molecular_Function_2021",
    "GO_Cellular_Component_2021",
    "KEGG_2021_Human",
    "Reactome_2022",
]


def ensure_dirs(*paths: Path) -> None:
    for p in paths:
        p.mkdir(parents=True, exist_ok=True)


def collapse_probes_to_genes(df: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse multiple probes to a single row per gene_symbol.
    Keep the row with smallest padj; tie-break by largest |logFC|.
    """
    d = df.copy()

    d["gene_symbol"] = d["gene_symbol"].astype(str).str.strip()
    d = d[d["gene_symbol"].notna()]
    d = d[d["gene_symbol"] != ""]
    d = d[d["gene_symbol"].str.lower() != "nan"]

    d["padj_num"] = pd.to_numeric(d.get("padj", pd.Series([np.nan] * len(d))), errors="coerce")
    d["logFC_num"] = pd.to_numeric(d.get("logFC", pd.Series([np.nan] * len(d))), errors="coerce")
    d["abs_logFC"] = d["logFC_num"].abs()

    # sort so the "best" probe per gene comes first
    d = d.sort_values(["padj_num", "abs_logFC"], ascending=[True, False])

    # keep first row per gene
    d = d.drop_duplicates(subset=["gene_symbol"], keep="first")

    # clean helper cols
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


def run_enrichr(
    genes: List[str],
    library: str,
    outdir_tables: Path,
    tag: str,
    contrast: str
) -> Optional[pd.DataFrame]:
    if len(genes) < 5:
        print(f"   - Skip {library} ({tag}): only {len(genes)} genes")
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

    out_tsv = outdir_tables / f"enrich_{contrast}_{library}_{tag}.tsv"
    res.to_csv(out_tsv, sep="\t", index=False)
    print(f"   - Saved: {out_tsv} (n_terms={len(res)})")
    return res


def plot_top_terms(res: pd.DataFrame, out_png: Path, title: str, top_n: int = 15) -> None:
    if res is None or len(res) == 0:
        return

    cols = res.columns.tolist()

    term_col = "Term" if "Term" in cols else cols[0]

    if "Combined_Score" in cols:
        score_col = "Combined_Score"
        d = res.head(top_n).copy()
        d["_score"] = pd.to_numeric(d[score_col], errors="coerce")
        xlab = "Combined Score"
    elif "Adjusted_P-value" in cols:
        score_col = "Adjusted_P-value"
        d = res.head(top_n).copy()
        d["_score"] = -np.log10(pd.to_numeric(d[score_col], errors="coerce"))
        xlab = "-log10(Adjusted P-value)"
    elif "P-value" in cols:
        score_col = "P-value"
        d = res.head(top_n).copy()
        d["_score"] = -np.log10(pd.to_numeric(d[score_col], errors="coerce"))
        xlab = "-log10(P-value)"
    else:
        d = res.head(top_n).copy()
        d["_score"] = np.arange(len(d), 0, -1)
        xlab = "Rank"

    d = d.iloc[::-1]

    plt.figure(figsize=(10, 5.5))
    plt.barh(d[term_col].astype(str), d["_score"].astype(float))
    plt.xlabel(xlab)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=220)
    plt.close()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--deg_dir", default="results/tables")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--fc", type=float, default=1.0)
    ap.add_argument("--top_terms", type=int, default=15)
    args = ap.parse_args()

    deg_dir = Path(args.deg_dir)
    out_tables = Path("results/tables")
    out_figs = Path("results/figures")
    ensure_dirs(out_tables, out_figs)

    deg_files = sorted(deg_dir.glob("deg_*.tsv"))
    if not deg_files:
        raise FileNotFoundError(f"No deg_*.tsv files found in {deg_dir}")

    print(">> Enrichment libraries:")
    for lib in LIBRARIES:
        print("  -", lib)

    for f in deg_files:
        contrast = f.stem.replace("deg_", "")
        print(f"\n>> Processing contrast: {contrast}")

        df = pd.read_csv(f, sep="\t")

        if "gene_symbol" not in df.columns:
            raise ValueError(f"{f} missing 'gene_symbol' column. Columns: {list(df.columns)}")

        df_gene = collapse_probes_to_genes(df)
        up, down = split_up_down(df_gene, fdr_thr=args.fdr, fc_thr=args.fc)

        print(f"   - Genes after probe→gene collapse: {df_gene['gene_symbol'].nunique()}")
        print(f"   - Significant UP: {len(up)}  DOWN: {len(down)}  (FDR<{args.fdr}, |logFC|>={args.fc})")

        for lib in LIBRARIES:
            res_up = run_enrichr(up, lib, out_tables, tag="UP", contrast=contrast)
            if res_up is not None and len(res_up) > 0:
                out_png = out_figs / f"enrich_{contrast}_{lib}_UP.png"
                plot_top_terms(res_up, out_png, f"{contrast} | {lib} | UP genes", top_n=args.top_terms)
                print(f"   - Saved: {out_png}")

            res_dn = run_enrichr(down, lib, out_tables, tag="DOWN", contrast=contrast)
            if res_dn is not None and len(res_dn) > 0:
                out_png = out_figs / f"enrich_{contrast}_{lib}_DOWN.png"
                plot_top_terms(res_dn, out_png, f"{contrast} | {lib} | DOWN genes", top_n=args.top_terms)
                print(f"   - Saved: {out_png}")

    print("\n✅ Step 4 complete. Next: Drug repurposing (src/05_drug_repurposing.py).")


if __name__ == "__main__":
    main()
