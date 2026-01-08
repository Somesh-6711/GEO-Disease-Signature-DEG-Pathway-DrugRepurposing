"""
Step 1: Download a GEO Series (GSE) and export:
- data/raw/expression.tsv   (features x samples)
- data/raw/phenotype.tsv    (sample metadata + inferred group labels)
- data/raw/platform.tsv     (probe annotation table, if available)

Default dataset: GSE60436 (PDR fibrovascular membrane: active/inactive + normal retina controls)
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
import GEOparse


DEFAULT_GSE = "GSE60436"


def parse_characteristics(chars: List[str]) -> Dict[str, str]:
    """
    Convert GEO 'characteristics_ch1' strings like:
      ["tissue: retina", "status: active"]
    into dict columns: {"tissue": "retina", "status": "active"}.
    """
    out: Dict[str, str] = {}
    for c in chars or []:
        if not isinstance(c, str):
            continue
        if ":" in c:
            k, v = c.split(":", 1)
            k = k.strip().lower()
            v = v.strip()
            # Keep first occurrence; don't overwrite unless empty
            if k and (k not in out or not out[k]):
                out[k] = v
    return out


def infer_group(text: str) -> str:
    """
    Heuristic group label inference from combined metadata text.
    You can later override/edit phenotype.tsv if needed.
    """
    t = (text or "").lower()

    # Common case/control cues
    if re.search(r"\binactive\b", t):
        return "inactive_fvm"
    if re.search(r"\bactive\b", t):
        return "active_fvm"
    if re.search(r"\bcontrol\b|\bhealthy\b|\bnormal\b", t) and "fvm" not in t:
        return "control_retina"
    if re.search(r"\bretina\b", t) and ("active" not in t and "inactive" not in t):
        return "control_retina"

    return "unknown"


def detect_value_column(gse: GEOparse.GSE) -> Tuple[str, str]:
    """
    Choose a column to pivot across samples.
    Prefer VALUE (common for microarray Series Matrix).
    Fallback to the first shared numeric-like column across GSM tables.
    Returns (values_col, index_col).
    """
    index_col = "ID_REF"

    # Fast path
    try:
        # Check at least one GSM has VALUE
        any_gsm = next(iter(gse.gsms.values()))
        if "VALUE" in any_gsm.table.columns:
            return "VALUE", index_col
    except Exception:
        pass

    # Fallback: find a column present in all GSM tables (besides ID_REF)
    common_cols = None
    for gsm in gse.gsms.values():
        cols = set(map(str, gsm.table.columns))
        if common_cols is None:
            common_cols = cols
        else:
            common_cols &= cols

    if not common_cols:
        raise RuntimeError("Could not find common columns across GSM tables to pivot.")

    candidates = [c for c in common_cols if c != index_col]
    # Prefer typical names
    for pref in ["VALUE", "COUNT", "COUNTS", "RPKM", "TPM", "FPKM", "signal"]:
        for c in candidates:
            if c.lower() == pref.lower():
                return c, index_col

    # Otherwise just pick first
    return sorted(candidates)[0], index_col


def export_platform(gse: GEOparse.GSE, out_path: Path) -> None:
    """
    Save platform annotation table (probe metadata) if available.
    """
    if not getattr(gse, "gpls", None):
        return
    # Usually one platform; if multiple, concatenate with a platform_id column
    frames = []
    for gpl_id, gpl in gse.gpls.items():
        tbl = gpl.table.copy()
        tbl.insert(0, "platform_id", gpl_id)
        frames.append(tbl)
    platform_df = pd.concat(frames, ignore_index=True)
    platform_df.to_csv(out_path, sep="\t", index=False)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gse", default=DEFAULT_GSE, help="GEO Series accession, e.g., GSE60436")
    ap.add_argument("--outdir", default="data/raw", help="Output directory (default: data/raw)")
    args = ap.parse_args()

    gse_id = args.gse.strip()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f">> Downloading {gse_id} into {outdir.resolve()}")
    gse = GEOparse.get_GEO(geo=gse_id, destdir=str(outdir), how="full")

    # ---- Expression matrix ----
    values_col, index_col = detect_value_column(gse)
    print(f">> Extracting expression matrix with values='{values_col}' index='{index_col}'")

    expr = gse.pivot_samples(values=values_col, index=index_col)
    # Ensure numeric where possible
    expr = expr.apply(pd.to_numeric, errors="coerce")

    # Save expression (features x samples)
    expr_path = outdir / "expression.tsv"
    expr.to_csv(expr_path, sep="\t", index=True)
    print(f">> Saved: {expr_path}  shape={expr.shape}")

    # ---- Phenotype table ----
    rows = []
    for gsm_id, gsm in gse.gsms.items():
        md = gsm.metadata or {}
        title = (md.get("title") or [""])[0]
        source = (md.get("source_name_ch1") or [""])[0]
        chars = md.get("characteristics_ch1") or []
        chars_dict = parse_characteristics(chars)
        combined_text = " ".join([title, source, " ".join(chars)])
        group = infer_group(combined_text)

        row = {
            "sample_id": gsm_id,
            "group": group,
            "title": title,
            "source_name_ch1": source,
            "raw_characteristics_ch1": " | ".join(chars),
        }
        # add parsed characteristic columns
        for k, v in chars_dict.items():
            row[f"ch1_{k}"] = v

        rows.append(row)

    pheno = pd.DataFrame(rows).sort_values("sample_id").reset_index(drop=True)

    pheno_path = outdir / "phenotype.tsv"
    pheno.to_csv(pheno_path, sep="\t", index=False)
    print(f">> Saved: {pheno_path}  n_samples={len(pheno)}")
    print(">> Group counts:")
    print(pheno["group"].value_counts(dropna=False).to_string())

    # ---- Platform annotations (probe metadata) ----
    platform_path = outdir / "platform.tsv"
    export_platform(gse, platform_path)
    if platform_path.exists():
        print(f">> Saved: {platform_path}")

    print("\n✅ Step 1 complete. Next: QC + normalization checks (src/02_qc.py).")


if __name__ == "__main__":
    main()
