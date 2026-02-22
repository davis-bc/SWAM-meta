"""
amr_unified.py

Snakemake script (called via script: directive).

Merges AMR gene detections from two evidence streams into a single
unified output, using copies-per-genome (cpg) as the shared abundance unit:

  1. Short-reads  – KMA alignment against AMRFinderPlus (high sensitivity,
                    no assembly context)
  2. Contigs      – AMRFinderPlus protein+nucleotide annotation on assembled
                    contigs (assembly context: molecule_type, taxonomy,
                    genomic coordinates)

Join key: gene_symbol
  = 'allele' column in short_reads_output.csv
  = 'gene'   column in contig_summary.tsv (feature_type == "AMR")

Abundance priority:
  When a gene is detected by both sources, contig cpg is preferred —
  it is derived from actual assembled sequence rather than read depth,
  giving a more accurate copy count.

Output columns (AMR_unified.csv):
  sample | gene_symbol | gene_family | product_name | class | subclass |
  type | subtype | cpg | evidence |
  short_reads_cpg | contig_cpg | molecule_type | taxonomy | n_contigs
"""

import pandas as pd
import os

# ---------------------------------------------------------------------------
# Inputs / outputs
# ---------------------------------------------------------------------------
sr_file      = snakemake.input.short_reads_amr   # short_reads_output.csv
contigs_file = snakemake.input.contig_summary     # contig_summary.tsv
catalog_file = snakemake.input.catalog            # ReferenceGeneCatalog.txt
out_file     = snakemake.output.csv

# ---------------------------------------------------------------------------
# 1. Load AMRFinderPlus catalog  (allele → gene_family, class, etc.)
# ---------------------------------------------------------------------------
catalog = pd.read_csv(catalog_file, sep="\t", dtype=str)

# Canonical metadata per allele (first occurrence wins)
catalog_meta = (
    catalog[["allele", "gene_family", "product_name",
             "scope", "type", "subtype", "class", "subclass"]]
    .drop_duplicates("allele")
    .rename(columns={"allele": "gene_symbol"})
    .set_index("gene_symbol")
)

# ---------------------------------------------------------------------------
# 2. Short-reads: sum cpg per sample × gene_symbol
# ---------------------------------------------------------------------------
sr = pd.read_csv(sr_file)

if not sr.empty:
    # 'allele' is the specific gene symbol (e.g. blaTEM-1)
    sr = sr.rename(columns={"allele": "gene_symbol"})
    sr_agg = (
        sr.groupby(["sample", "gene_symbol"])
        .agg(short_reads_cpg=("cpg", "sum"))
        .reset_index()
    )
else:
    sr_agg = pd.DataFrame(columns=["sample", "gene_symbol", "short_reads_cpg"])

# ---------------------------------------------------------------------------
# 3. Contigs: aggregate AMR features per sample × gene_symbol
#    contig_summary.tsv columns:
#      sample, contig_id, abundance_cpg, molecule_type, taxonomy,
#      feature_type, gene, start, stop, strand
# ---------------------------------------------------------------------------
cs = pd.read_csv(contigs_file, sep="\t", dtype=str)
amr_cs = cs[cs["feature_type"] == "AMR"].copy() if not cs.empty else pd.DataFrame()

if not amr_cs.empty:
    amr_cs = amr_cs.rename(columns={"gene": "gene_symbol"})
    amr_cs["abundance_cpg"] = pd.to_numeric(amr_cs["abundance_cpg"], errors="coerce").fillna(0.0)

    def _agg_contig_group(grp):
        # Molecule types: comma-sep sorted unique (chromosomal default)
        mol = sorted(grp["molecule_type"].dropna().unique())
        mol_str = ",".join(mol) if mol else "chromosome"

        # Taxonomy: most common non-empty hit
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
        .apply(_agg_contig_group)
        .reset_index()
    )
else:
    contig_agg = pd.DataFrame(
        columns=["sample", "gene_symbol", "contig_cpg",
                 "molecule_type", "taxonomy", "n_contigs"]
    )

# ---------------------------------------------------------------------------
# 4. Outer join on sample × gene_symbol
# ---------------------------------------------------------------------------
unified = sr_agg.merge(contig_agg, on=["sample", "gene_symbol"], how="outer")

# ---------------------------------------------------------------------------
# 5. Annotate with catalog metadata (gene_family, class, etc.)
#    Left-join so unrecognised gene_symbols still appear with NA metadata.
# ---------------------------------------------------------------------------
unified = unified.join(catalog_meta, on="gene_symbol", how="left")

# ---------------------------------------------------------------------------
# 6. Evidence column
# ---------------------------------------------------------------------------
has_sr = unified["short_reads_cpg"].notna() & (unified["short_reads_cpg"] > 0)
has_ct = unified["contig_cpg"].notna()       & (unified["contig_cpg"] > 0)

unified["evidence"] = "short_reads_only"
unified.loc[has_ct & ~has_sr, "evidence"] = "contigs_only"
unified.loc[has_sr &  has_ct, "evidence"] = "both"

# ---------------------------------------------------------------------------
# 7. Preferred cpg (contig when available, else short-reads)
# ---------------------------------------------------------------------------
unified["cpg"] = unified["contig_cpg"].where(has_ct, unified["short_reads_cpg"])

# Tidy numeric fill / types
unified["short_reads_cpg"] = unified["short_reads_cpg"].fillna(0.0)
unified["contig_cpg"]      = unified["contig_cpg"].fillna(0.0)
unified["n_contigs"]       = unified["n_contigs"].fillna(0).astype(int)
unified["molecule_type"]   = unified["molecule_type"].fillna("")
unified["taxonomy"]        = unified["taxonomy"].fillna("")

# ---------------------------------------------------------------------------
# 8. Sort and write
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

os.makedirs(os.path.dirname(out_file), exist_ok=True)
unified[out_cols].to_csv(out_file, index=False)
