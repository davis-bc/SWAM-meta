#!/usr/bin/env python3
"""
generate_mock_data.py

Generates all mock test data for SWAM-meta:
  - 2 synthetic metagenome sample FASTQ pairs (mock1, mock2) in test/data/
  - Mini AMRFinderPlus KMA database (test/dbs/amrfinder/)
  - Mini SCG DIAMOND database (test/dbs/scg/)
  - Mini human genome for host filtering (test/dbs/human/)
  - Mini UniRef50 MMseqs2 taxonomy database (test/dbs/uniref50/)

Genomes are seeded-random (repeat-free, ~50% GC) with the real blaTEM-1
beta-lactamase gene embedded so that:
  - megahit assembles long contigs (no mutation in reads → clean overlaps)
  - MetaBAT2 can bin the chromosome contigs (300 kb > 200 kb min bin size)
  - AMRFinderPlus detects blaTEM-1 in contigs and MAGs
  - KMA detects blaTEM-1 reads in the short-reads summary

Run from the repository root:

    python test/scripts/generate_mock_data.py

Prerequisites (installed into PATH or via conda):
    kma index, diamond makedb, minimap2 (index only),
    mmseqs createdb / createtaxdb / createindex

These are the same tools used by the pipeline envs so they will be available
once --use-conda has created the environments, or if run inside them.
"""

import random
import gzip
import os
import subprocess
import sys
import textwrap

SEED = 42
random.seed(SEED)

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TEST_DATA = os.path.join(ROOT, "test", "data")
TEST_DBS  = os.path.join(ROOT, "test", "dbs")

# ---------------------------------------------------------------------------
# blaTEM-1 coding sequence (861 bp, NG_050145.1:101-961)
# Source: AMRFinderPlus database, accession WP_000027057.1
# Encodes broad-spectrum class A beta-lactamase TEM-1; detected by
# AMRFinderPlus and used for short-read AMR quantification via KMA.
# ---------------------------------------------------------------------------
BLA_TEM1 = (
    "atgagtattcaacattttcgtgtcgcccttattcccttttttgcggcattttgccttcctgtttttgct"
    "cacccagaaacgctggtgaaagtaaaagatgctgaagatcagttgggtgcacgagtgggttacatcgaa"
    "ctggatctcaacagcggtaagatccttgagagttttcgccccgaagaacgttttccaatgatgagcact"
    "tttaaagttctgctatgtggtgcggtattatcccgtgttgacgccgggcaagagcaactcggtcgccgc"
    "atacactattctcagaatgacttggttgagtactcaccagtcacagaaaagcatcttacggatggcatg"
    "acagtaagagaattatgcagtgctgccataaccatgagtgataacactgctgccaacttacttctgaca"
    "acgatcggaggaccgaaggagctaaccgcttttttgcacaacatgggggatcatgtaactcgccttgat"
    "cgttgggaaccggagctgaatgaagccataccaaacgacgagcgtgacaccacgatgcctgcagcaatg"
    "gcaacaacgttgcgcaaactattaactggcgaactacttactctagcttcccggcaacaattaataga"
    "ctggatggaggcggataaagttgcaggaccacttctgcgctcggcccttccggctggctggtttattgc"
    "tgataaatctggagccggtgagcgtgggtctcgcggtatcattgcagcactggggccagatggtaagcc"
    "ctcccgtatcgtagttatctacacgacggggagtcaggcaactatggatgaacgaaatagacagatcgc"
    "tgagataggtgcctcactgattaagcattggtaa"
).upper()

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

BASES = "ACGT"
# Weights for ~50 % GC (like E. coli)
_WEIGHTS = [0.25, 0.25, 0.25, 0.25]


def rand_seq(length, rng=None):
    r = rng or random
    return "".join(r.choices(BASES, weights=_WEIGHTS, k=length))


def rev_comp(seq):
    comp = str.maketrans("ACGT", "TGCA")
    return seq.translate(comp)[::-1]


def phred_str(length, qual=35):
    return chr(qual + 33) * length


def write_fastq_gz(path, records):
    """records: list of (name, seq, qual_str)"""
    with gzip.open(path, "wt") as fh:
        for name, seq, qual in records:
            fh.write(f"@{name}\n{seq}\n+\n{qual}\n")


def write_fasta(path, records):
    """records: list of (name, seq)"""
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for chunk in textwrap.wrap(seq, 60):
                fh.write(chunk + "\n")


def run(cmd, **kwargs):
    print(f"  $ {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, **kwargs)
    if result.returncode != 0:
        print(f"  STDOUT: {result.stdout[-500:]}")
        print(f"  STDERR: {result.stderr[-500:]}")
        sys.exit(f"Command failed: {cmd}")


def embed(genome_list, seq, position):
    """Overwrite genome_list[position:position+len(seq)] with seq."""
    for i, b in enumerate(seq):
        genome_list[position + i] = b

# ---------------------------------------------------------------------------
# 1. Reference genomes
# ---------------------------------------------------------------------------

def make_references():
    """
    Create three synthetic 'genomes' using seeded-random sequences:

      bac_chrom  300 kbp  – blaTEM-1 embedded at position 100,000
      plasmid     20 kbp  – blaTEM-1 embedded at position 1,000
      phage       20 kbp  – no AMR (viral mobile element background)

    Using a seeded RNG guarantees the sequences are deterministic AND
    repeat-free (random 300 kbp sequences have vanishingly few k-mer
    repeats), so megahit assembles them into very long contigs.

    Returns dict of name -> uppercase sequence string.
    """
    print("[1] Generating synthetic reference sequences...")
    rng = random.Random(SEED)

    chrom = list(rand_seq(300_000, rng))
    embed(chrom, BLA_TEM1, 100_000)

    plasmid = list(rand_seq(20_000, rng))
    embed(plasmid, BLA_TEM1, 1_000)

    phage = list(rand_seq(20_000, rng))

    refs = {
        "bac_chrom": "".join(chrom),
        "plasmid":   "".join(plasmid),
        "phage":     "".join(phage),
    }
    for name, seq in refs.items():
        print(f"  {name}: {len(seq):,} bp")
    return refs

# ---------------------------------------------------------------------------
# 2. FASTQs — paired-end reads tiled across each genome
# ---------------------------------------------------------------------------

def simulate_reads_from_genome(genome, coverage, read_len, insert_mean,
                                insert_std, rng, error_rate=0.001):
    """
    Simulate paired-end reads from *genome* at *coverage*×.
    Reads are tiled at even spacing with Gaussian insert-size jitter.
    A very low error rate (0.1 %) makes reads essentially error-free while
    retaining a tiny amount of natural-looking variation.
    """
    n_pairs = max(1, int(len(genome) * coverage / (2 * read_len)))
    step    = max(1, len(genome) // n_pairs)
    pairs   = []
    for i in range(n_pairs):
        ins = max(read_len * 2, int(rng.gauss(insert_mean, insert_std)))
        ins = min(ins, len(genome) - 1)
        # Start positions wrap around the circular genome
        start = (i * step + rng.randint(-step // 4, step // 4)) % (len(genome) - ins)
        start = max(0, start)
        r1 = list(genome[start: start + read_len])
        r2 = list(rev_comp(genome[start + ins - read_len: start + ins]))
        # Introduce rare sequencing errors
        for r in (r1, r2):
            for j in range(len(r)):
                if rng.random() < error_rate:
                    r[j] = rng.choice(BASES)
        pairs.append(("".join(r1), "".join(r2)))
    return pairs


def make_fastqs(refs):
    """
    Generate two samples:
      mock1: chromosome 50×, plasmid 80×, phage 30×
      mock2: same genome, different RNG seed (mimics a second environmental sample)

    Coverage levels are chosen so that:
      - chromosome (300 kbp × 50×) → ~50 k read pairs → excellent megahit assembly
      - plasmid (20 kbp × 80×) → ~5.3 k read pairs → high-copy plasmid signature
      - phage (20 kbp × 30×) → ~2 k read pairs → background viral reads
    """
    READ_LEN     = 150
    INSERT_MEAN  = 500
    INSERT_STD   = 50

    sample_cfg = {
        "mock1": {
            "bac_chrom": (50, random.Random(42)),
            "plasmid":   (80, random.Random(43)),
            "phage":     (30, random.Random(44)),
        },
        "mock2": {
            "bac_chrom": (40, random.Random(142)),
            "plasmid":   (60, random.Random(143)),
            "phage":     (25, random.Random(144)),
        },
    }

    print("[2] Generating mock FASTQ files...")
    for sample, cfg in sample_cfg.items():
        r1_recs, r2_recs = [], []
        idx = 1
        for region, (cov, rng) in cfg.items():
            pairs = simulate_reads_from_genome(
                refs[region], coverage=cov, read_len=READ_LEN,
                insert_mean=INSERT_MEAN, insert_std=INSERT_STD, rng=rng)
            for r1, r2 in pairs:
                q  = phred_str(READ_LEN)
                nm = f"{sample}.{idx}"
                r1_recs.append((f"{nm}/1", r1, q))
                r2_recs.append((f"{nm}/2", r2, q))
                idx += 1

        write_fastq_gz(os.path.join(TEST_DATA, f"{sample}_R1.fastq.gz"), r1_recs)
        write_fastq_gz(os.path.join(TEST_DATA, f"{sample}_R2.fastq.gz"), r2_recs)
        print(f"  {sample}: {idx-1:,} read pairs written")

# ---------------------------------------------------------------------------
# 3. Mini AMRFinderPlus FASTA + KMA index
# ---------------------------------------------------------------------------

def make_amrfinder_db():
    """
    Build the KMA database used for short-read AMR quantification.
    Includes blaTEM-1 so reads originating from the embedded gene are counted.
    Header format mirrors the real AMRFinderPlus AMR_CDS.fa:
        WP_accession|nucleotide_acc|1|1|allele|gene_family|product nt_acc:start-end
    """
    print("[3] Building mini AMRFinderPlus KMA database...")
    db_dir    = os.path.join(TEST_DBS, "amrfinder")
    fasta_path = os.path.join(db_dir, "AMR_CDS.fa")
    rng = random.Random(SEED + 100)

    records = []
    for i in range(10):
        seq  = rand_seq(900, rng)
        name = (f"WP_FAKE{i:06d}.1|NG_FAKE{i:06d}.1|1|1|"
                f"fakeGene_{i}|fakeFamily|synthetic_AMR_gene_{i} "
                f"NG_FAKE{i:06d}.1:101-1000")
        records.append((name, seq))

    # blaTEM-1 entry matching the real AMRFinderPlus header format
    records.append((
        "WP_000027057.1|NG_050145.1|1|1|blaTEM-1|blaTEM|"
        "broad-spectrum_class_A_beta-lactamase_TEM-1 NG_050145.1:101-961",
        BLA_TEM1
    ))

    write_fasta(fasta_path, records)

    meta_path = os.path.join(db_dir, "ReferenceGeneCatalog.txt")
    with open(meta_path, "w") as fh:
        fh.write("allele\tgene_family\tproduct_name\tscope\ttype\tsubtype\tclass\tsubclass\t"
                 "refseq_protein_accession\trefseq_nucleotide_accession\t"
                 "genbank_nucleotide_accession\n")
        for i in range(10):
            fh.write(f"fakeGene_{i}\tfakeFamily\tsynthetic_AMR_gene_{i}\tcore\tAMR\tAMR\t"
                     f"BETA-LACTAM\tCLASS_A\tWP_FAKE{i:06d}.1\t"
                     f"NG_FAKE{i:06d}.1\tNG_FAKE{i:06d}.1\n")
        fh.write("blaTEM-1\tblaTEM\tbroad-spectrum_class_A_beta-lactamase_TEM-1\t"
                 "core\tAMR\tAMR\tBETA-LACTAM\tCLASS_A\t"
                 "WP_000027057.1\tNG_050145.1\tNG_050145.1\n")

    run(f"kma index -i {fasta_path} -o {db_dir}/afp_db")
    print(f"  AMRFinder KMA db written to {db_dir}/afp_db")

# ---------------------------------------------------------------------------
# 4. Mini SCG DIAMOND database
# ---------------------------------------------------------------------------

def make_scg_db():
    print("[4] Building mini SCG DIAMOND database...")
    db_dir = os.path.join(TEST_DBS, "scg")
    fasta_path = os.path.join(db_dir, "SCGs_40.fasta")
    rng = random.Random(SEED + 200)

    records = []
    for i in range(40):
        aa_seq = "".join(rng.choices("ACDEFGHIKLMNPQRSTVWY", k=300))
        records.append((f"SCG_{i:02d}", aa_seq))

    write_fasta(fasta_path, records)
    run(f"diamond makedb --in {fasta_path} -d {db_dir}/scg_db --quiet")
    print(f"  SCG DIAMOND db written to {db_dir}/scg_db.dmnd")

# ---------------------------------------------------------------------------
# 5. Mini human genome for host filtering
# ---------------------------------------------------------------------------

def make_human_genome():
    print("[5] Writing mini human reference genome...")
    db_dir = os.path.join(TEST_DBS, "human")
    # A tiny chromosome stub — minimap2 will index it; no real human reads in
    # mock data so nothing will be filtered out.
    fasta_path = os.path.join(db_dir, "human_mini.fna")
    rng = random.Random(SEED + 300)
    records = [("chr1_stub", rand_seq(10_000, rng))]
    write_fasta(fasta_path, records)
    run(f"gzip -kf {fasta_path}")
    print(f"  Human genome stub written to {fasta_path}.gz")

# ---------------------------------------------------------------------------
# 6. Mini UniRef50 MMseqs2 taxonomy database
# ---------------------------------------------------------------------------

def make_uniref50_mmseqs_db():
    """
    Build a minimal MMseqs2 taxonomy database that mmseqs easy-taxonomy
    can run against.

    Layout required by createtaxdb:
      - protein FASTA with headers: >unirefXX_ACCESSION taxid=<int> <description>
      - mapping file: acc -> taxid (nodes.dmp + names.dmp style via --ncbi-tax-dump)
      OR simpler: use --tax-mapping-file  acc<TAB>taxid

    We use the --tax-mapping-file approach (supported in mmseqs >= 13).
    """
    print("[6] Building mini UniRef50 MMseqs2 taxonomy database...")
    db_dir = os.path.join(TEST_DBS, "uniref50")

    # NCBI-style taxonomy: 3 taxa (bacteria genus level)
    # taxid 1 = root, 131567 = cellular organisms, 2 = Bacteria,
    # 1224 = Proteobacteria, 1236 = Gammaproteobacteria,
    # 91347 = Enterobacterales, 543 = Enterobacteriaceae, 561 = Escherichia
    # We'll use fake taxids 1001–1005

    fake_taxa = {
        1: ("root", 1),
        131567: ("cellular organisms", 1),
        2: ("Bacteria", 131567),
        1224: ("Proteobacteria", 2),
        1236: ("Gammaproteobacteria", 1224),
        91347: ("Enterobacterales", 1236),
        543: ("Enterobacteriaceae", 91347),
        561: ("Escherichia", 543),
        562: ("Escherichia coli", 561),
    }

    # Write NCBI-style nodes.dmp and names.dmp
    nodes_path = os.path.join(db_dir, "nodes.dmp")
    names_path = os.path.join(db_dir, "names.dmp")
    with open(nodes_path, "w") as fh:
        for taxid, (name, parent) in fake_taxa.items():
            rank = "species" if taxid == 562 else "genus" if taxid == 561 else "no rank"
            fh.write(f"{taxid}\t|\t{parent}\t|\t{rank}\t|\t\t|\n")
    with open(names_path, "w") as fh:
        for taxid, (name, _) in fake_taxa.items():
            fh.write(f"{taxid}\t|\t{name}\t|\t\t|\tscientific name\t|\n")

    # Write empty merged.dmp (required by mmseqs createtaxdb but can be empty)
    with open(os.path.join(db_dir, "merged.dmp"), "w") as fh:
        pass

    # Protein FASTA: 50 sequences assigned to taxid 562
    fasta_path = os.path.join(db_dir, "uniref50_mini.fasta")
    acc_taxid = []
    records = []
    rng = random.Random(SEED + 400)
    for i in range(50):
        acc = f"UniRef50_TEST{i:04d}"
        aa_seq = "".join(rng.choices("ACDEFGHIKLMNPQRSTVWY", k=200))
        records.append((f"{acc} n=1 Tax=Escherichia coli TaxID=562 RepID=TEST{i:04d}_ECOLI", aa_seq))
        acc_taxid.append((acc, 562))
    write_fasta(fasta_path, records)

    # Tax mapping file
    mapping_path = os.path.join(db_dir, "tax_mapping.tsv")
    with open(mapping_path, "w") as fh:
        for acc, taxid in acc_taxid:
            fh.write(f"{acc}\t{taxid}\n")

    # Build MMseqs2 DB
    seqdb = os.path.join(db_dir, "uniref50_mmseqs")
    run(f"mmseqs createdb {fasta_path} {seqdb}")
    run(f"mmseqs createtaxdb {seqdb} {db_dir}/tmp "
        f"--ncbi-tax-dump {db_dir} --tax-mapping-file {mapping_path}")
    run(f"mmseqs createindex {seqdb} {db_dir}/tmp --search-type 2")
    print(f"  MMseqs2 UniRef50 db written to {seqdb}")

# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    refs = make_references()
    make_fastqs(refs)
    # DB-building steps require tools — skip gracefully if not found
    for fn, label in [
        (make_amrfinder_db, "AMRFinderPlus KMA db"),
        (make_scg_db,       "SCG DIAMOND db"),
        (make_human_genome, "human genome stub"),
        (make_uniref50_mmseqs_db, "UniRef50 MMseqs2 db"),
    ]:
        try:
            fn()
        except SystemExit as e:
            print(f"  WARNING: skipped {label}: {e}")

    print("\nDone. Test data ready in test/")
