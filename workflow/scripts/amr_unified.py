"""
amr_unified.py

Snakemake script (called via script: directive).

Merges AMR gene detections from two evidence streams using cpg normalization:

  1. Short-reads  – KMA alignment against AMRFinderPlus (high sensitivity)
  2. Contigs      – AMRFinderPlus protein annotation on assembled contigs
                    (provides molecule_type, taxonomy, genomic context)

Join key: gene_symbol
  = 'allele' column in short_reads_output.csv
  = 'gene'   column in contig_summary.tsv (feature_type == "AMR")

Abundance priority: contig cpg preferred when available; else short-reads cpg.

Metadata priority:
  1. contig_amr.tsv  (direct AMRFinderPlus output on contigs — most accurate)
  2. short_reads_output.csv  (already catalog-annotated by R script)
  3. ReferenceGeneCatalog.txt  (fallback)

Outputs:
  AMR_unified.csv             – one row per sample × gene_symbol
  AMR_abundance_summary.csv   – per-sample totals: AMR cpg + pBI143 + crAss001
"""

import pandas as pd
import os
import glob as _glob

# ---------------------------------------------------------------------------
# Inputs / outputs
# ---------------------------------------------------------------------------
sr_file        = str(snakemake.input.short_reads_amr) if snakemake.input.get("short_reads_amr") else None
contigs_file   = str(snakemake.input.contig_summary)  if snakemake.input.get("contig_summary")  else None
catalog_file   = snakemake.input.catalog               # always present when amr_unified runs
contig_amr_dir = snakemake.params.contig_amr_dir       # dir with *_contig_amr.tsv
markers_file   = str(snakemake.input.markers_cpg) if snakemake.input.get("markers_cpg") else None
out_unified    = snakemake.output.csv
out_summary    = snakemake.output.amr_abundance_summary

META_COLS = ["gene_family", "product_name", "scope", "type", "subtype", "class", "subclass"]

# ---------------------------------------------------------------------------
# 1. Metadata priority 1: contig_amr.tsv files
#    Columns: Element symbol, Element name, Scope, Type, Subtype, Class, Subclass
# ---------------------------------------------------------------------------
contig_amr_files = _glob.glob(os.path.join(contig_amr_dir, "*_contig_amr.tsv"))
contig_meta = pd.DataFrame()
if contig_amr_files:
    dfs = []
    for f in contig_amr_files:
        try:
            df = pd.read_csv(f, sep="\t", dtype=str)
            if not df.empty:
                dfs.append(df)
        except Exception:
            pass
    if dfs:
        raw = pd.concat(dfs, ignore_index=True)
        # normalise column names
        raw.columns = [c.lower().replace(" ", "_") for c in raw.columns]
        gene_col = next((c for c in ["element_symbol", "gene_symbol"] if c in raw.columns), None)
        if gene_col:
            raw = raw.rename(columns={gene_col: "gene_symbol",
                                      "element_name": "product_name"})
            contig_meta = (
                raw[["gene_symbol", "product_name", "scope",
                     "type", "subtype", "class", "subclass"]]
                .drop_duplicates("gene_symbol")
                .set_index("gene_symbol")
            )
            # gene_family: not in AMRFinder direct output; derive from gene_symbol
            # (same symbol used as fallback — catalog join fills it for known genes)
            contig_meta["gene_family"] = contig_meta.index

# ---------------------------------------------------------------------------
# 2. Metadata priority 2: catalog
# ---------------------------------------------------------------------------
catalog = pd.read_csv(catalog_file, sep="\t", dtype=str)
catalog_meta = (
    catalog[["allele", "gene_family", "product_name",
             "scope", "type", "subtype", "class", "subclass"]]
    .drop_duplicates("allele")
    .rename(columns={"allele": "gene_symbol"})
    .set_index("gene_symbol")
)

def get_meta(gene_sym):
    """Return metadata Series for a gene_symbol, checking sources in priority order."""
    for src in [contig_meta, catalog_meta]:
        if not src.empty and gene_sym in src.index:
            return src.loc[gene_sym]
    return pd.Series({c: "" for c in ["gene_family"] + META_COLS})

# ---------------------------------------------------------------------------
# 3. Short-reads: sum cpg per sample × gene_symbol
#    Also capture per-row metadata (from R script catalog join)
# ---------------------------------------------------------------------------
if sr_file and os.path.exists(sr_file):
    sr = pd.read_csv(sr_file)
else:
    sr = pd.DataFrame()

if not sr.empty:
    sr = sr.rename(columns={"allele": "gene_symbol"})
    sr_agg = (
        sr.groupby(["sample", "gene_symbol"])
        .agg(short_reads_cpg=("cpg", "sum"))
        .reset_index()
    )
    # Build supplementary metadata from short-reads output (already annotated)
    sr_meta_cols = [c for c in META_COLS if c in sr.columns]
    if sr_meta_cols:
        sr_meta = (
            sr[["gene_symbol"] + sr_meta_cols]
            .drop_duplicates("gene_symbol")
            .set_index("gene_symbol")
        )
    else:
        sr_meta = pd.DataFrame()
else:
    sr_agg = pd.DataFrame(columns=["sample", "gene_symbol", "short_reads_cpg"])
    sr_meta = pd.DataFrame()

# ---------------------------------------------------------------------------
# 4. Contigs: aggregate AMR features per sample × gene_symbol
# ---------------------------------------------------------------------------
if contigs_file and os.path.exists(contigs_file):
    cs = pd.read_csv(contigs_file, sep="\t", dtype=str)
else:
    cs = pd.DataFrame()
amr_cs = cs[cs["feature_type"] == "AMR"].copy() if not cs.empty else pd.DataFrame()

if not amr_cs.empty:
    amr_cs = amr_cs.rename(columns={"gene": "gene_symbol"})
    amr_cs["abundance_cpg"] = pd.to_numeric(amr_cs["abundance_cpg"], errors="coerce").fillna(0.0)

    def _agg(grp):
        mol = sorted(grp["molecule_type"].dropna().unique())
        mol_str = ",".join(m for m in mol if m.strip()) or "chromosome"
        tax_vals = grp["taxonomy"].dropna()
        tax_vals = tax_vals[tax_vals.str.strip().ne("")]
        taxonomy = tax_vals.mode().iloc[0] if not tax_vals.empty else ""
        return pd.Series({
            "contig_cpg":    grp["abundance_cpg"].sum(),
            "molecule_type": mol_str,
            "taxonomy":      taxonomy,
            "n_contigs":     int(grp["contig_id"].nunique()),
        })

    contig_agg = (
        amr_cs.groupby(["sample", "gene_symbol"])
        [["molecule_type", "taxonomy", "abundance_cpg", "contig_id"]]
        .apply(_agg)
        .reset_index()
    )
else:
    contig_agg = pd.DataFrame(
        columns=["sample", "gene_symbol", "contig_cpg",
                 "molecule_type", "taxonomy", "n_contigs"]
    )

# ---------------------------------------------------------------------------
# 5. Outer join
# ---------------------------------------------------------------------------
unified = sr_agg.merge(contig_agg, on=["sample", "gene_symbol"], how="outer")

# ---------------------------------------------------------------------------
# 6. Attach metadata
# ---------------------------------------------------------------------------
all_genes = unified["gene_symbol"].unique()
meta_rows = []
for g in all_genes:
    row = {"gene_symbol": g}
    # Priority: contig_meta > sr_meta > catalog_meta
    for src in [contig_meta, sr_meta if not sr_meta.empty else pd.DataFrame(), catalog_meta]:
        if not src.empty and g in src.index:
            for col in ["gene_family", "product_name", "scope", "type", "subtype", "class", "subclass"]:
                if col in src.columns and col not in row:
                    row[col] = src.loc[g, col]
    meta_rows.append(row)

meta_df = pd.DataFrame(meta_rows).set_index("gene_symbol")
unified = unified.join(meta_df, on="gene_symbol", how="left")

# ---------------------------------------------------------------------------
# 7. Evidence + preferred cpg
# ---------------------------------------------------------------------------
has_sr = unified["short_reads_cpg"].notna() & (unified["short_reads_cpg"] > 0)
has_ct = unified["contig_cpg"].notna()       & (unified["contig_cpg"] > 0)

unified["evidence"] = "short_reads_only"
unified.loc[has_ct & ~has_sr, "evidence"] = "contigs_only"
unified.loc[has_sr &  has_ct, "evidence"] = "both"

unified["cpg"] = unified["contig_cpg"].where(has_ct, unified["short_reads_cpg"])

unified["short_reads_cpg"] = unified["short_reads_cpg"].fillna(0.0)
unified["contig_cpg"]      = unified["contig_cpg"].fillna(0.0)
unified["n_contigs"]       = unified["n_contigs"].fillna(0).astype(int)
unified["molecule_type"]   = unified["molecule_type"].fillna("")
unified["taxonomy"]        = unified["taxonomy"].fillna("")

# ---------------------------------------------------------------------------
# 8. Sort and write AMR_unified.csv
# ---------------------------------------------------------------------------
unified = unified.sort_values(["sample", "cpg"], ascending=[True, False])

COLS = [
    "sample", "gene_symbol", "gene_family", "product_name",
    "class", "subclass", "type", "subtype",
    "cpg", "evidence",
    "short_reads_cpg", "contig_cpg",
    "molecule_type", "taxonomy", "n_contigs",
]
out_cols = [c for c in COLS if c in unified.columns]
os.makedirs(os.path.dirname(out_unified), exist_ok=True)
unified[out_cols].to_csv(out_unified, index=False)

# ---------------------------------------------------------------------------
# 9. AMR_abundance_summary.csv
#    AMR_total (cpg) = sum of unified cpg per sample (uses contig cpg where
#    available, else short-reads cpg — same logic as the unified table)
#    pBI143 and crAss001 come from markers_cpg.csv (produced by R script)
# ---------------------------------------------------------------------------
markers = pd.read_csv(markers_file) if markers_file and os.path.exists(markers_file) else pd.DataFrame()

amr_total = (
    unified.groupby("sample")
    .agg(**{"AMR_total (cpg)": ("cpg", "sum")})
    .reset_index()
)

summary = amr_total.merge(markers, on="sample", how="left") if not markers.empty else amr_total
# fill missing marker columns
for mc in ["pBI143 (cpg)", "crAss001 (cpg)"]:
    if mc not in summary.columns:
        summary[mc] = 0.0
    else:
        summary[mc] = summary[mc].fillna(0.0)

summary["AMR_total (cpg)"] = summary["AMR_total (cpg)"].fillna(0.0)
os.makedirs(os.path.dirname(out_summary), exist_ok=True)
summary.to_csv(out_summary, index=False)
