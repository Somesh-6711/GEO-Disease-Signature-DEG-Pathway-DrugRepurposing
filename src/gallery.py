from __future__ import annotations
from pathlib import Path
import pandas as pd

FIG = Path("results/figures")
TAB = Path("results/tables")
OUT = Path("results/gallery.md")

CONTRASTS = [
    "active_fvm_vs_control_retina",
    "inactive_fvm_vs_control_retina",
]

def md_img(path: Path, alt: str) -> str:
    # GitHub-friendly relative path
    return f"![{alt}]({path.as_posix()})"

def md_table_tsv(tsv: Path, cols: list[str] | None = None, n: int = 10) -> str:
    if not tsv.exists():
        return "_(missing)_"
    df = pd.read_csv(tsv, sep="\t")
    if cols:
        cols = [c for c in cols if c in df.columns]
        if cols:
            df = df[cols]
    return df.head(n).to_markdown(index=False)

def main():
    lines = []
    lines.append("# Results Gallery\n")
    lines.append("Auto-generated gallery of QC, DEG, enrichment, and drug repurposing plots.\n")

    # QC
    lines.append("## QC\n")
    for fname, alt in [
        ("qc_boxplot_raw.png", "QC boxplot (raw)"),
        ("qc_missingness.png", "QC missingness"),
        ("qc_pca.png", "QC PCA"),
        ("qc_sample_correlation.png", "QC sample correlation"),
    ]:
        p = FIG / fname
        if p.exists():
            lines.append(md_img(p, alt))
            lines.append("")
    lines.append("---\n")

    # Per contrast
    for c in CONTRASTS:
        lines.append(f"## Contrast: {c}\n")

        # Volcano
        volc = FIG / f"volcano_{c}.png"
        if volc.exists():
            lines.append("### Volcano\n")
            lines.append(md_img(volc, f"Volcano {c}"))
            lines.append("")

        # Enrichment (Reactome)
        up = FIG / f"enrich_{c}_Reactome_2022_UP.png"
        dn = FIG / f"enrich_{c}_Reactome_2022_DOWN.png"
        if up.exists() or dn.exists():
            lines.append("### Reactome Enrichment\n")
            if up.exists():
                lines.append(md_img(up, f"Reactome UP {c}"))
                lines.append("")
            if dn.exists():
                lines.append(md_img(dn, f"Reactome DOWN {c}"))
                lines.append("")

        # Drug rep
        dsi = FIG / f"drugrep_{c}_DSigDB_top.png"
        lincsd = FIG / f"drugrep_{c}_LINCS_L1000_Chem_Pert_down_top.png"
        if dsi.exists() or lincsd.exists():
            lines.append("### Drug Repurposing\n")
            if dsi.exists():
                lines.append(md_img(dsi, f"DSigDB top {c}"))
                lines.append("")
            if lincsd.exists():
                lines.append(md_img(lincsd, f"LINCS down top {c}"))
                lines.append("")

        # Top tables
        lines.append("### Top Drug Candidates (tables)\n")
        dsi_t = TAB / f"drugrep_{c}_DSigDB_TOP.tsv"
        lin_t = TAB / f"drugrep_{c}_LINCS_L1000_Chem_Pert_down_TOP.tsv"
        if dsi_t.exists():
            lines.append("**DSigDB (Top 10)**\n")
            lines.append(md_table_tsv(dsi_t, cols=["Term", "ReversalScore", "FDR_DOWN", "FDR_UP"], n=10))
            lines.append("")
        if lin_t.exists():
            lines.append("**LINCS Chem Pert DOWN (Top 10)**\n")
            lines.append(md_table_tsv(lin_t, cols=["Term", "ReversalScore", "FDR_DOWN", "FDR_UP"], n=10))
            lines.append("")

        lines.append("---\n")

    OUT.write_text("\n".join(lines), encoding="utf-8")
    print(f">> Saved: {OUT}")

if __name__ == "__main__":
    main()