#!/usr/bin/env python3
"""
ViroWatch knowledge-graph CSV exporter.

Reads per-sample pipeline outputs and writes flat CSVs for bulk import
into a Neo4j knowledge graph via BULK_MERGE_* Cypher templates:

  kg/sample.csv                  — Sample node
  kg/assembly.csv                — Assembly node (linked to Sample)
  kg/biodata_files.csv           — BioDataFile nodes (FASTQ input + consensus FASTA)
  kg/contigs.csv                 — Contig nodes (Flye info joined with Medaka FASTA)
  kg/stanford_alignments.csv     — StanfordHIVDRAlignment nodes (one row per contig×gene)
  kg/stanford_predictions.csv    — StanfordHIVDRPrediction + Drug + DrugClass rows
  kg/mutations.csv               — Mutation nodes (full flag set)
  kg/blast_hits.csv              — ReferenceGenome + Organism rows from BLAST hits
  kg/variant_calling_run.csv     — VariantCallingRun node (+ CALLED_FROM the FASTQ)
  kg/variants.csv                — Variant nodes + CALLED edge quality (one row per
                                   VCF row × SnpEff ANN entry; AGAINST the reference)

Contig IDs are namespaced as  {sample_id}:{flye_contig_name}  so they are
globally unique across samples (Flye always starts from contig_1).
"""
import argparse
import csv
import gzip
import hashlib
import os
import sys
from datetime import datetime, timezone
from pathlib import Path


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()


def _write_csv(path, fieldnames, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def _bool_str(v):
    """Normalise Python / SierraPy booleans to lowercase string."""
    if isinstance(v, bool):
        return str(v).lower()
    if isinstance(v, str):
        return v.lower() if v.lower() in ("true", "false") else v
    return ""


# ── exporters ─────────────────────────────────────────────────────────────────

def export_sample(kg_dir, sample_id):
    _write_csv(kg_dir / "sample.csv", ["sample_id"], [{"sample_id": sample_id}])


def export_assembly_and_contigs(kg_dir, sample_id, sample_dir, report_dir):
    consensus_fa = sample_dir / "medaka_consensus" / "consensus.fasta"
    flye_info    = sample_dir / "flye" / "assembly_info.txt"

    if not (consensus_fa.exists() and flye_info.exists()):
        return {}

    sys.path.insert(0, str(report_dir))
    from ContextBuilder.assembly_info import get_length_with_flye_info  # noqa: PLC0415

    assembly_id = f"{sample_id}_assembly"
    created_at  = datetime.now(timezone.utc).date().isoformat()

    _write_csv(
        kg_dir / "assembly.csv",
        ["assembly_id", "sample_id", "assembler", "created_at"],
        [{"assembly_id": assembly_id, "sample_id": sample_id,
          "assembler": "Flye", "created_at": created_at}],
    )

    raw = get_length_with_flye_info(str(consensus_fa), str(flye_info))

    contig_rows = []
    contig_id_map = {}  # flye_name → namespaced contig_id
    for item in raw:
        flye_name = item["contig_id"]
        contig_id = f"{sample_id}:{flye_name}"
        contig_id_map[flye_name] = contig_id
        contig_rows.append({
            "contig_id":          contig_id,
            "assembly_id":        assembly_id,
            "contig_name":        flye_name,
            "length":             item["length"],
            "coverage":           item["coverage"],
            "is_circular":        _bool_str(item["is_circular"]),
            "is_repeated_region": _bool_str(item["is_repeated_region"]),
            "multiplicity":       item["intracontig_copy_number"],
            "alt_group":          item["alternative_contigs_in_group"],
            "graph_path":         item["graph_path"],
            "sequence":           item["sequence"],
            "sequence_hash":      item["sequence_hash"],
            "hash_algorithm":     item["hash_algorithm"],
        })

    _write_csv(
        kg_dir / "contigs.csv",
        ["contig_id", "assembly_id", "contig_name", "length", "coverage",
         "is_circular", "is_repeated_region", "multiplicity", "alt_group",
         "graph_path", "sequence", "sequence_hash", "hash_algorithm"],
        contig_rows,
    )

    return contig_id_map, assembly_id, str(consensus_fa.resolve())


def export_biodata_files(kg_dir, sample_id, fastq_path, consensus_fa_path, assembly_id):
    fastq = Path(fastq_path)
    rows = [
        {
            "uri":         str(fastq.resolve()),
            "file_type":   "FASTQ",
            "compressed":  "true" if fastq.suffix in (".gz", ".bz2") else "false",
            "sha256":      "",   # skipped — FASTQ can be tens of GB
            "assembly_id": "",   # FASTQ is not the product of assembly
            "sample_id":   sample_id,
        },
        {
            "uri":         consensus_fa_path,
            "file_type":   "FASTA",
            "compressed":  "false",
            "sha256":      _sha256(consensus_fa_path),
            "assembly_id": assembly_id,
            "sample_id":   sample_id,
        },
    ]
    _write_csv(
        kg_dir / "biodata_files.csv",
        ["uri", "file_type", "compressed", "sha256", "assembly_id", "sample_id"],
        rows,
    )


def export_sierrapy(kg_dir, sample_id, sample_dir, report_dir, contig_id_map):
    sierrapy_json = sample_dir / "sierrapy.json"
    if not sierrapy_json.exists():
        return

    sys.path.insert(0, str(report_dir))
    from ContextBuilder.sierrapy import SierraPyResult  # noqa: PLC0415

    result      = SierraPyResult.read_sierrapy_json(sierrapy_json)
    result_sha  = result.get_result_file_hash()
    file_ts     = datetime.fromtimestamp(
        sierrapy_json.stat().st_ctime, tz=timezone.utc
    ).isoformat()

    alignment_rows   = []
    prediction_rows  = []
    mutation_rows    = []

    # Drug resistance panel + alignments
    for contig_name, dr_reports in result.get_drug_resistance_panal():
        if not dr_reports:
            continue
        contig_id = contig_id_map.get(contig_name, f"{sample_id}:{contig_name}")

        for report in dr_reports:
            gene            = report["gene"]["name"]
            db_version      = report["version"]["text"]
            db_published    = report["version"]["publishDate"]

            alignment_rows.append({
                "result_sha256":          result_sha,
                "contig_id_prefix":       contig_id,
                "gene":                   gene,
                "timestamp":              file_ts,
                "database_version":       db_version,
                "database_published_date": db_published,
            })

            for ds in report["drugScores"]:
                drug       = ds["drug"]
                drug_class = ds["drugClass"]
                drug_name  = drug["name"]
                prediction_id = f"{contig_id}:{gene}:{drug_name}"
                prediction_rows.append({
                    "prediction_id":          prediction_id,
                    "contig_id":              contig_id,
                    "sample_id":              sample_id,
                    "gene":                   gene,
                    "drug_name":              drug_name,
                    "drug_full_name":         drug.get("fullName", ""),
                    "drug_abbr":              drug.get("displayAbbr", ""),
                    "drug_class":             drug_class["name"],
                    "score":                  ds["score"],
                    "level":                  ds["level"],
                    "interpretation":         ds["text"],
                    "database_version":       db_version,
                    "database_published_date": db_published,
                })

    # Mutations (full flag set, deduplicated by gene+text+sha)
    seen_mutations = set()
    for contig_name, genes in result.get_mutations(simple_return=False).items():
        for gene, muts in genes.items():
            for m in muts:
                key = (result_sha, gene, m.get("text", ""))
                if key in seen_mutations:
                    continue
                seen_mutations.add(key)
                text = m.get("text", "")
                mutation_rows.append({
                    "gene":              gene,
                    "text":              text,
                    "mutation_id":       f"{gene}:{text}",
                    "primary_type":      m.get("primaryType", ""),
                    "is_sdrm":           _bool_str(m.get("isSDRM", False)),
                    "position":          m.get("position", ""),
                    "has_stop":          _bool_str(m.get("hasStop", False)),
                    "is_apobec_mutation": _bool_str(m.get("isApobecMutation", False)),
                    "is_apobec_drm":     _bool_str(m.get("isApobecDRM", False)),
                    "is_insertion":      _bool_str(m.get("isInsertion", False)),
                    "is_deletion":       _bool_str(m.get("isDeletion", False)),
                    "is_unusual":        _bool_str(m.get("isUnusual", False)),
                    "result_sha256":     result_sha,
                })

    _write_csv(
        kg_dir / "stanford_alignments.csv",
        ["result_sha256", "contig_id_prefix", "gene", "timestamp",
         "database_version", "database_published_date"],
        alignment_rows,
    )
    _write_csv(
        kg_dir / "stanford_predictions.csv",
        ["prediction_id", "contig_id", "sample_id", "gene",
         "drug_name", "drug_full_name", "drug_abbr", "drug_class",
         "score", "level", "interpretation",
         "database_version", "database_published_date"],
        prediction_rows,
    )
    _write_csv(
        kg_dir / "mutations.csv",
        ["gene", "text", "mutation_id", "primary_type", "is_sdrm", "position",
         "has_stop", "is_apobec_mutation", "is_apobec_drm",
         "is_insertion", "is_deletion", "is_unusual", "result_sha256"],
        mutation_rows,
    )


def export_blast_hits(kg_dir, sample_id, sample_dir, report_dir, contig_id_map):
    blast_files = [
        (sample_dir / "blast" / "los_alamos.blast.json", "LosAlamos_HIV"),
        (sample_dir / "blast" / "core_nt.blast.json",    "NCBI_core_nt"),
    ]
    existing = [(p, db) for p, db in blast_files if p.exists()]
    if not existing:
        return

    sys.path.insert(0, str(report_dir))
    from ContextBuilder.core import BLASTResult  # noqa: PLC0415

    rows = []
    seen_accessions = set()

    for blast_path, source_db in existing:
        blast = BLASTResult.read_json(blast_path)
        for contig_name in blast.get_query_names():
            contig_id = contig_id_map.get(contig_name, f"{sample_id}:{contig_name}")
            for hit in blast.get_hits(contig_name):
                acc = hit["accession"]
                if acc in seen_accessions:
                    continue
                seen_accessions.add(acc)
                taxid = hit.get("taxid")
                rows.append({
                    "accession_no":    acc,
                    "name":            hit.get("subject_title", ""),
                    "source_database": source_db,
                    "taxid":           taxid if taxid and taxid != 0 else "",
                    "sciname":         hit.get("sciname", ""),
                    "query_contig_id": contig_id,
                    "sample_id":       sample_id,
                })

    if rows:
        _write_csv(
            kg_dir / "blast_hits.csv",
            ["accession_no", "name", "source_database", "taxid", "sciname",
             "query_contig_id", "sample_id"],
            rows,
        )


# ── variants (medaka + SnpEff) ─────────────────────────────────────────────────
# VCF parsing mirrors nosograph-py's nosograph/utils/vcf.py::parse_medaka_vcf so
# the emitted CSVs match the lib's Variant / VariantCallingRun schema and the
# BULK_MERGE_Variants Cypher contract — WITHOUT importing the private lib.

def _open_vcf(vcf_path):
    """Open a VCF for text reading, transparently handling gzip (.gz)."""
    if str(vcf_path).endswith(".gz"):
        return gzip.open(vcf_path, "rt", encoding="utf-8")
    return open(vcf_path, encoding="utf-8")


def _vcf_int(val):
    if val is None or val in (".", ""):
        return None
    try:
        return int(val)
    except ValueError:
        return None


def _parse_info(info_str):
    result = {}
    for token in info_str.split(";"):
        if "=" in token:
            k, v = token.split("=", 1)
            result[k] = v
        else:
            result[token] = "true"
    return result


def _parse_format(format_str, sample_str):
    keys = format_str.split(":")
    values = sample_str.split(":")
    return dict(zip(keys, values))


def _parse_ann(ann_str):
    """SnpEff ANN INFO field → one dict per annotation entry.

    ANN pipe fields: ALLELE|EFFECT|IMPACT|GENE_NAME|GENE_ID|FEATURE_TYPE|
    TRANSCRIPT_ID|BIOTYPE|RANK|HGVS.c|HGVS.p|...
    """
    results = []
    for entry in ann_str.split(","):
        fields = entry.split("|")
        if len(fields) < 11:
            continue
        results.append({
            "effect":    fields[1],
            "impact":    fields[2],
            "gene_name": fields[3] or None,
            "hgvs_c":    fields[9],
            "hgvs_p":    fields[10],
        })
    return results


def _ref_accession(ref_fa_path):
    """First FASTA header token (the accession, e.g. 'AF164485.1'), or None."""
    if not ref_fa_path:
        return None
    try:
        with open(ref_fa_path, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith(">"):
                    return line[1:].split()[0]
    except OSError:
        return None
    return None


def _nn(val):
    """None → empty string (CSV null convention); everything else str()."""
    return "" if val is None else str(val)


def export_variants(kg_dir, sample_id, sample_dir, fastq_path, ref_fa_path):
    """Export medaka+SnpEff variants to variant_calling_run.csv + variants.csv.

    Runs independently of the assembly branch — variants are called against the
    reference from the filtered reads, so they exist even without a Flye assembly.
    Prefers the SnpEff-annotated bgzip VCF; falls back to the un-annotated tagged
    VCF (EFFECT/IMPACT/gene_name then empty).
    """
    variants_dir = sample_dir / "variants"
    vcf_path = None
    for cand in ("medaka.annotated.vcf.gz", "tagged.vcf"):
        p = variants_dir / cand
        if p.exists():
            vcf_path = p
            break
    if vcf_path is None:
        return

    ref_acc = _ref_accession(ref_fa_path) or ""
    process_run_id = f"{sample_id}:variant_calling:medaka:{ref_acc}"
    fastq_uri = str(Path(fastq_path).resolve())

    _write_csv(
        kg_dir / "variant_calling_run.csv",
        ["process_run_id", "process", "tool", "reference", "run_id",
         "fastq_uri", "sample_id"],
        [{
            "process_run_id": process_run_id,
            "process":        "variant_calling",
            "tool":           "medaka",
            "reference":      ref_acc,
            "run_id":         "",
            "fastq_uri":      fastq_uri,
            "sample_id":      sample_id,
        }],
    )

    rows = []
    with _open_vcf(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom, pos_str, _id, ref, alt, qual_str, flt, info_str = cols[:8]
            format_str = cols[8] if len(cols) > 8 else ""
            sample_str = cols[9] if len(cols) > 9 else ""

            info = _parse_info(info_str)
            fmt = _parse_format(format_str, sample_str) if format_str else {}
            annotations = _parse_ann(info.get("ANN", ""))

            base = {
                "process_run_id": process_run_id,
                "REF_ACC": ref_acc or chrom,
                "POS":     pos_str,
                "REF":     ref,
                "ALT":     alt,
                "CHROM":   chrom,
                "TYPE":    _nn(info.get("TYPE")),
                "DP":      _nn(_vcf_int(info.get("DP"))),
                "GT":      _nn(fmt.get("GT")),
                "QUAL":    "" if qual_str in (".", "") else qual_str,
                "GQ":      _nn(_vcf_int(fmt.get("GQ"))),
                "AO":      "",   # medaka does not emit AO/RO
                "RO":      "",
                "FILTER":  "" if flt == "." else flt,
            }

            if not annotations:
                rows.append({**base, "hgvs_c": "", "hgvs_p": "",
                             "EFFECT": "", "IMPACT": "", "gene_name": ""})
                continue
            for ann in annotations:
                rows.append({
                    **base,
                    "hgvs_c":    ann["hgvs_c"],
                    "hgvs_p":    ann["hgvs_p"],
                    "EFFECT":    _nn(ann["effect"]),
                    "IMPACT":    _nn(ann["impact"]),
                    "gene_name": _nn(ann["gene_name"]),
                })

    _write_csv(
        kg_dir / "variants.csv",
        ["process_run_id", "REF_ACC", "POS", "REF", "ALT", "hgvs_c", "hgvs_p",
         "CHROM", "TYPE", "EFFECT", "IMPACT", "gene_name",
         "DP", "GT", "QUAL", "GQ", "AO", "RO", "FILTER"],
        rows,
    )


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Export ViroWatch per-sample outputs to Neo4j-compatible CSVs"
    )
    parser.add_argument("--sample",     required=True, help="Sample ID")
    parser.add_argument("--sample_dir", required=True, help="Path to sample output directory")
    parser.add_argument("--report_dir", required=True, help="Path to the report/ directory")
    parser.add_argument("--fastq",      required=True, help="Path to the input FASTQ file")
    parser.add_argument("--ref_fa",     default="",    help="Reference FASTA (for the variant REF_ACC/accession)")
    args = parser.parse_args()

    sample_id  = args.sample
    sample_dir = Path(args.sample_dir)
    report_dir = Path(args.report_dir)
    kg_dir     = sample_dir / "kg"
    kg_dir.mkdir(exist_ok=True)

    export_sample(kg_dir, sample_id)

    result = export_assembly_and_contigs(kg_dir, sample_id, sample_dir, report_dir)
    if result:
        contig_id_map, assembly_id, consensus_fa_path = result
        export_biodata_files(
            kg_dir, sample_id, args.fastq, consensus_fa_path, assembly_id
        )
        export_sierrapy(kg_dir, sample_id, sample_dir, report_dir, contig_id_map)
        export_blast_hits(kg_dir, sample_id, sample_dir, report_dir, contig_id_map)
    else:
        contig_id_map = {}
        export_sierrapy(kg_dir, sample_id, sample_dir, report_dir, contig_id_map)
        export_blast_hits(kg_dir, sample_id, sample_dir, report_dir, contig_id_map)

    # Variants are independent of the assembly branch (called against the
    # reference from the filtered reads), so export them unconditionally.
    export_variants(kg_dir, sample_id, sample_dir, args.fastq, args.ref_fa)

    written = sorted(p.name for p in kg_dir.glob("*.csv"))
    print(f"[kg_export] {sample_id}: {len(written)} CSV(s) → {kg_dir}")
    for name in written:
        print(f"  {name}")


if __name__ == "__main__":
    main()
