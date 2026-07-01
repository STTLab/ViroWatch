#!/usr/bin/env python3
# SPDX-License-Identifier: MPL-2.0
"""ViroWatch metagenomics (Kraken2) knowledge-graph CSV exporter.

Reads one sample's Kraken2 report (produced by the ``kraken2`` classification step in
``bin/virowatch_run.sh``) and writes a single flat CSV for manual bulk import into Neo4j
via the LOAD DATA template in ``assets/virowatch_cypher_templates.csv``.

Unlike NosoGraph's exporter (which materialises an ``(:Organism)`` per taxon and links
each with an ``IDENTIFIED`` edge), ViroWatch treats Kraken2 as a *quick read-set QC* — not
a searchable taxonomy graph. So the per-taxon records are folded onto the
``TaxonomicClassification`` node itself as a single JSON string property (``taxa_json``),
sorted by abundance descending, after an adaptive z-score filter collapses the long tail of
trace taxa into one ``"Other"`` bucket.

Target subgraph:

  (Sample)-[:CLASSIFIED_IN]->(:ProcessRun:TaxonomicClassification {taxa_json, ...})
          -[:CLASSIFIED_FROM]->(:BioDataFile{FASTQ})   # reuses the node kg_export.py emits

Emitted (under ``<outdir>/kg/``):

  taxonomic_classification.csv  — one row: the TaxonomicClassification ProcessRun

Only species (rank ``S``) and genus (rank ``G``) rows are considered. Per STTLab
conventions: ``taxid`` is a STRING; abundance is the Kraken2 clade fraction
(clade reads / classified reads).

Adaptive z-score filter (``--z-min``, ``--min-taxa``): abundances are heavily right-skewed
(one dominant taxon + a long trace tail), so the score is computed in log10 space using a
**robust modified z-score** — ``z = (log_ab - median) / (1.4826 * MAD)``. Median/MAD are used
rather than mean/std because the single dominant taxon (HIV routinely >90 %) would otherwise
inflate the standard deviation and pull the whole trace tail inside 1 SD, so the "Other"
bucket would almost never fire. Taxa with ``z < z_min`` fold into a single ``"Other"`` row
(summed read_count/abundance, plus ``n_grouped``). Guards: if fewer than ``min_taxa`` taxa,
or the MAD is zero (no spread in the bulk), filtering is skipped and all taxa are kept.
"""
import argparse
import json
import math
import statistics
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

# Standard Kraken2 report: 6 tab-separated columns, no header.
KRAKEN2_COLS = ["pct", "clade_reads", "taxon_reads", "rank_code", "taxid", "name"]
# Consider species + genus only.
KEEP_RANKS = {"S", "G"}


def _read_kraken2(report: Path) -> pd.DataFrame:
    df = pd.read_csv(
        report, sep="\t", header=None, names=KRAKEN2_COLS, dtype={"taxid": str}
    )
    # Names are indented by depth in the lineage; rank codes can carry trailing spaces.
    df["rank_code"] = df["rank_code"].astype(str).str.strip()
    df["name"] = df["name"].astype(str).str.strip()
    return df


def _clade_reads(df: pd.DataFrame, rank_code: str) -> int:
    """Clade-read count of the first row with this exact rank code, or 0."""
    sub = df[df["rank_code"] == rank_code]
    return int(sub["clade_reads"].iloc[0]) if not sub.empty else 0


def _adaptive_filter(taxa, z_min, min_taxa):
    """Fold unusually-small taxa (robust log-abundance z-score < z_min) into one "Other" row.

    Uses a modified z-score ``(log_ab - median) / (1.4826 * MAD)`` so the dominant taxon does
    not desensitise the filter. ``taxa`` is a list of dicts with
    taxid/sciname/rank/read_count/abundance. Returns a new list, sorted by abundance
    descending, with the grouped tail (if any) as a trailing ``"Other"`` entry carrying
    summed read_count/abundance and ``n_grouped``.
    """
    if len(taxa) < min_taxa:
        return sorted(taxa, key=lambda t: t["abundance"], reverse=True)

    logs = [math.log10(t["abundance"]) for t in taxa]
    median = statistics.median(logs)
    mad = statistics.median([abs(x - median) for x in logs])
    if mad == 0:
        return sorted(taxa, key=lambda t: t["abundance"], reverse=True)

    kept, grouped = [], []
    for t, log_ab in zip(taxa, logs):
        z = (log_ab - median) / (1.4826 * mad)
        (kept if z >= z_min else grouped).append(t)

    if grouped:
        kept.append({
            "taxid":      "other",
            "sciname":    "Other",
            "rank":       "X",
            "read_count": sum(t["read_count"] for t in grouped),
            "abundance":  round(sum(t["abundance"] for t in grouped), 6),
            "n_grouped":  len(grouped),
        })
    return sorted(kept, key=lambda t: t["abundance"], reverse=True)


def export(kg_dir, sample_id, kraken2_report, reads, tool, z_min, min_taxa):
    process_run_id = f"{sample_id}_taxclass"
    created_at = datetime.now(timezone.utc).date().isoformat()

    df = _read_kraken2(kraken2_report)
    classified = _clade_reads(df, "R")    # root clade = all classified reads
    unclassified = _clade_reads(df, "U")

    # --- Identified taxa (species + genus) ---
    taxa = []
    for _, r in df[df["rank_code"].isin(KEEP_RANKS)].iterrows():
        clade = int(r["clade_reads"])
        abundance = round(clade / classified, 6) if classified else 0.0
        if abundance <= 0:
            continue
        taxa.append({
            "taxid":      str(r["taxid"]),
            "sciname":    r["name"],
            "rank":       r["rank_code"],
            "read_count": clade,
            "abundance":  abundance,
        })

    taxa = _adaptive_filter(taxa, z_min, min_taxa)

    reads_path = Path(reads).resolve() if reads else None
    fastq_uri = str(reads_path) if reads_path else ""

    pd.DataFrame(
        [{
            "process_run_id":     process_run_id,
            "sample_id":          sample_id,
            "tool":               tool,
            "created_at":         created_at,
            "classified_reads":   classified,
            "unclassified_reads": unclassified,
            "fastq_uri":          fastq_uri,
            "taxa_json":          json.dumps(taxa, separators=(",", ":")),
        }],
        columns=["process_run_id", "sample_id", "tool", "created_at",
                 "classified_reads", "unclassified_reads", "fastq_uri", "taxa_json"],
    ).to_csv(kg_dir / "taxonomic_classification.csv", index=False)


def main():
    p = argparse.ArgumentParser(
        description="Export ViroWatch per-sample Kraken2 QC to a Neo4j CSV"
    )
    p.add_argument("--sample", required=True, help="Sample ID")
    p.add_argument("--kraken2-report", required=True, help="Kraken2 report (.kraken2.report.txt)")
    p.add_argument("--reads", default="", help="Input FASTQ used for classification (optional)")
    p.add_argument("--tool", default="kraken2", help="Classifier tool label (default: kraken2)")
    p.add_argument("--outdir", default=".", help="Output dir; CSV is written to <outdir>/kg/")
    p.add_argument("--z-min", type=float, default=-1.0,
                   help="Log-abundance z-score cutoff; taxa below fold into 'Other' (default: -1.0)")
    p.add_argument("--min-taxa", type=int, default=3,
                   help="Below this taxa count the z-filter is skipped (default: 3)")
    args = p.parse_args()

    kg_dir = Path(args.outdir) / "kg"
    kg_dir.mkdir(parents=True, exist_ok=True)

    export(
        kg_dir,
        args.sample,
        Path(args.kraken2_report),
        args.reads or None,
        args.tool,
        args.z_min,
        args.min_taxa,
    )

    print(f"[meta_kg_export] {args.sample}: taxonomic_classification.csv -> {kg_dir}")


if __name__ == "__main__":
    main()
