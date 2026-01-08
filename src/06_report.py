from __future__ import annotations

from pathlib import Path
import pandas as pd


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def head_table_md(df: pd.DataFrame, cols: list[str], n: int = 10) -> str:
    d = df.copy()
    d = d[cols].head(n)
    return d.to_markdown(index=False)


def load_top_enrich(contrast: str, library: str, tag: str, n: int = 10) -> str:
    p = Path("results/tables") / f"enrich_{contrast}_{library}_{tag}.tsv"
    if not p.exists():
        return "_(not found)_"
    df = read_tsv(p)
    # columns vary slightly; these are common in gseapy
    cols = [c for c in ["Term", "Adjusted_P-value", "P-value", "Combined_Score", "Genes"] if c in df.columns]
    if not cols:
        cols = df.columns[:5].tolist()
    return head_table_md(df, cols, n=n)


def load_top_drugs(contrast: str, library: str, n: int = 10) -> str:
    p = Path("results/tables") / f"drugrep_{contrast}_{library}_TOP.tsv"
    if not p.exists():
        return "_(not found)_"
    df = read_tsv(p)
    cols = [c for c in ["Term", "ReversalScore", "FDR_DOWN", "FDR_UP"] if c in df.columns]
    if not cols:
        cols = df.columns[:5].tolist()
    return head_table_md(df, cols, n=n)


def main():
    out = Path("results/summary.md")
    out.parent.mkdir(parents=True, exist_ok=True)

    contrasts = [
        "active_fvm_vs_control_retina",
        "inactive_fvm_vs_control_retina",
    ]

    lines = []
    lines.append("# GEO Disease Signature → Pathways → Drug Repurposing Report\n")
    lines.append("This report summarizes QC, differential expression, enrichment, and drug-repurposing candidates.\n")
    lines.append("**Disclaimer:** Research/portfolio use only. Not medical advice.\n")

    # QC summary
    qc_sum = Path("results/tables/qc_summary.tsv")
    if qc_sum.exists():
        dfqc = read_tsv(qc_sum)
        cols = dfqc.columns.tolist()[:10]
        lines.append("## QC Summary\n")
        lines.append(dfqc.head(20).to_markdown(index=False))
        lines.append("\n")

    for c in contrasts:
        lines.append(f"## Contrast: {c}\n")

        # DEG counts
        deg = Path("results/tables") / f"deg_{c}.tsv"
        if deg.exists():
            dfd = read_tsv(deg)
            # collapse to gene level is done in later steps; here give raw counts
            lines.append(f"- DEG rows (probe-level): **{len(dfd)}**\n")

        lines.append("### Top Enriched Pathways (Reactome)\n")
        lines.append("**UP genes**\n")
        lines.append(load_top_enrich(c, "Reactome_2022", "UP", n=10))
        lines.append("\n**DOWN genes**\n")
        lines.append(load_top_enrich(c, "Reactome_2022", "DOWN", n=10))
        lines.append("\n")

        lines.append("### Drug Repurposing Candidates\n")
        lines.append("**DSigDB (Top)**\n")
        lines.append(load_top_drugs(c, "DSigDB", n=10))
        lines.append("\n**LINCS L1000 Chem Pert (down, Top)**\n")
        lines.append(load_top_drugs(c, "LINCS_L1000_Chem_Pert_down", n=10))
        lines.append("\n")

        lines.append("### Key Figures\n")
        lines.append(f"- Volcano: `results/figures/volcano_{c}.png`\n")
        lines.append(f"- Enrichment (Reactome UP): `results/figures/enrich_{c}_Reactome_2022_UP.png`\n")
        lines.append(f"- Enrichment (Reactome DOWN): `results/figures/enrich_{c}_Reactome_2022_DOWN.png`\n")
        lines.append(f"- Drug candidates (DSigDB): `results/figures/drugrep_{c}_DSigDB_top.png`\n")
        lines.append("\n---\n")

    out.write_text("\n".join(lines), encoding="utf-8")
    print(f">> Saved report: {out}")


if __name__ == "__main__":
    main()
