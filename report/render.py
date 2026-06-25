#!/usr/bin/env python
import argparse
import os
import sys
import pandas as pd
from jinja2 import Environment, FileSystemLoader

GENOME_LENGTH = 9719
DRAW_WIDTH = 280
PX_PER_BP = DRAW_WIDTH / GENOME_LENGTH

ROWS = {
    1: (9.5,20.5),
    2: (24.5,32.5),
    3: (38.5,46.5),
    4: (52.5,60.5),
    5: (66.5,74.5),
}

HIV_GENES = [
    dict(name="gag",     start=790,  end=2292, row=2, strand="+"),
    dict(name="gag-pol", start=2085, end=5096, row=3, strand="+"),
    dict(name="vif",     start=5041, end=5619, row=2, strand="+"),
    dict(name="vpr",     start=5559, end=5850, row=3, strand="+"),
    dict(name="tat",     start=5831, end=6045, row=2, strand="+"),
    dict(name="tat",     start=8379, end=8469, row=2, strand="+"),
    dict(name="rev",     start=5970, end=6045, row=3, strand="+"),
    dict(name="rev",     start=8379, end=8653, row=3, strand="+"),
    dict(name="vpu",     start=6062, end=6310, row=4, strand="+"),
    dict(name="env",     start=6225, end=8795, row=5, strand="+"),
    dict(name="nef",     start=8797, end=9417, row=2, strand="+"),
]

INTERPRETATION_SIR_MAP = {
    'susceptible': 'S', 'potential low-level resistance': 'S',
    'low-level resistance': 'I', 'intermediate resistance': 'I',
    'high-level resistance': 'R'
}

DRUG_DISPLAY_TO_FULLNAME = {
    "ABC": "Abacavir",
    "AZT": "Zidovudine",
    "D4T": "Stavudine",
    "DDI": "Didanosine",
    "FTC": "Emtricitabine",
    "3TC": "Lamivudine",
    "TDF": "Tenofovir Disoproxil Fumarate",
    "DOR": "Doravirine",
    "DPV": "Dapivirine",
    "EFV": "Efavirenz",
    "ETR": "Etravirine",
    "NVP": "Nevirapine",
    "RPV": "Rilpivirine",
    "BIC": "Bictegravir",
    "CAB": "Cabotegravir",
    "DTG": "Dolutegravir",
    "EVG": "Elvitegravir",
    "RAL": "Raltegravir",
    "NFV": "Nelfinavir",
    "ATV/r": "Atazanavir/ritonavir",
    "DRV/r": "Darunavir/ritonavir",
    "LPV/r": "Lopinavir/ritonavir",
    "FPV/r": "Fosamprenavir/ritonavir",
    "IDV/r": "Indinavir/ritonavir",
    "SQV/r": "Saquinavir/ritonavir",
    "TPV/r": "Tipranavir/ritonavir",
}


def overlap(a1, a2, b1, b2):
    start = max(a1, b1)
    end = min(a2, b2)
    if start < end:
        return start, end
    return None


def bp_to_px(bp):
    return bp * PX_PER_BP


def arrow_forward(x1, x2, y1, y2):
    mid = (y1 + y2) / 2
    tip = 2
    return f"M{x1},{y1} L{x2-tip},{y1} L{x2},{mid} L{x2-tip},{y2} L{x1},{y2} L{x1+tip},{mid} Z"


def arrow_reverse(x1, x2, y1, y2):
    mid = (y1 + y2) / 2
    tip = 2
    return f"M{x1+tip},{y1} L{x2},{y1} L{x2-tip},{mid} L{x2},{y2} L{x1+tip},{y2} L{x1},{mid} Z"


def splice_ticks(start, end, row, step=300):
    y1, y2 = ROWS[row]
    ticks = []
    for bp in range(start, end, step):
        x = bp_to_px(bp)
        ticks.append({"x": x, "y1": y1+1, "y2": y2-1})
    return ticks


def build_backbone_segments(highlight_regions):
    segments = []
    for start, end, strand in highlight_regions:
        if start > end:
            start, end = end, start
        segments.append({"x": bp_to_px(start), "width": bp_to_px(end - start), "strand": strand})
    return segments


def build_genes(highlight_regions):
    out = []
    for i, g in enumerate(HIV_GENES):
        x1 = bp_to_px(g["start"])
        x2 = bp_to_px(g["end"])
        y1, y2 = ROWS[g["row"]]
        path = arrow_forward(x1, x2, y1, y2) if g["strand"] == "+" else arrow_reverse(x1, x2, y1, y2)
        segments = []
        for start, end in highlight_regions:
            if start > end:
                start, end = end, start
            ov = overlap(g["start"], g["end"], start, end)
            if ov:
                s, e = ov
                segments.append({"x": bp_to_px(s), "width": bp_to_px(e - s), "y": y1, "height": y2 - y1})
        out.append({
            "id": f"gene{i}",
            "name": g["name"],
            "path": path,
            "segments": segments,
            "ticks": splice_ticks(g["start"], g["end"], g["row"]),
        })
    return out


def build_hiv_image_from_hit(blast_record):
    hsps_subject = [(i[1], i[2]) for i in blast_record['hsps_subject']]
    hsps_query_strand = ['+' if i[0] == 'Plus' else '-' for i in blast_record['hsps_query']]
    hsps_subject_with_strand = [(*pair, n) for pair, n in zip(hsps_subject, hsps_query_strand)]
    return {
        'backbone_raw': hsps_subject_with_strand,
        'image_genes': build_genes(hsps_subject),
        'image_backbone_segments': build_backbone_segments(hsps_subject_with_strand),
        'image_genome_length': blast_record['subject_length'],
    }


def get_vl(vl_csv, sample_id):
    df = pd.read_csv(vl_csv)
    filtered = df.loc[df['sample_id'] == sample_id]
    return filtered[['date', 'vl']].to_dict('list')


def get_cd4(cd4_csv, sample_id):
    df = pd.read_csv(cd4_csv)
    filtered = df.loc[df['sample_id'] == sample_id]
    return filtered[['date', 'cd4_pct', 'cd4_count']].to_dict('list')


def build_drug_resistance_panel(sierrapy_result, sample_id):
    data = sierrapy_result.get_drug_resistance_panal()
    rows = []
    for contig, reports in data:
        for report in reports:
            version = report['version']['text']
            publish_date = report['version']['publishDate']
            gene = report['gene']['name']
            for drug_score in report['drugScores']:
                rows.append({
                    'result_sha256': sierrapy_result.get_result_file_hash(),
                    'sample_id': sample_id,
                    'contig': contig,
                    'gene': gene,
                    'version': version,
                    'publish_date': publish_date,
                    'drug_class': drug_score['drugClass']['name'],
                    'drug': drug_score['drug']['name'],
                    'drug_abbr': drug_score['drug']['displayAbbr'],
                    'score': drug_score['score'],
                    'level': drug_score['level'],
                    'interpretation': drug_score['text'],
                    'SIR': INTERPRETATION_SIR_MAP[drug_score['text'].lower()],
                })

    df = pd.DataFrame(rows)
    max_idx = df.groupby("drug_abbr")["level"].idxmax()
    summary = df.groupby("drug_abbr").agg(
        score=("score", "max"),
        level=("level", "max"),
        drug_class=("drug_class", "first"),
    )
    summary["interpretation"] = df.loc[max_idx].set_index("drug_abbr")["interpretation"]

    support_pos = df[df["level"] > 0].groupby("drug_abbr")["contig"].unique()
    support_all = df.groupby("drug_abbr")["contig"].unique()
    supported = support_pos.reindex(summary.index).combine_first(support_all)

    summary["supported_contig"] = supported.apply(list)
    summary["drug_full_name"] = summary.index.map(DRUG_DISPLAY_TO_FULLNAME)
    summary["SIR"] = summary["interpretation"].map(lambda x: INTERPRETATION_SIR_MAP[x.lower()])
    return summary.reset_index().to_dict("records")


def main():
    parser = argparse.ArgumentParser(description="Render per-sample ViroWatch HTML report")
    parser.add_argument("--sample",     required=True,  help="Sample ID")
    parser.add_argument("--sample_dir", required=True,  help="Path to sample pipeline output directory")
    parser.add_argument("--report_dir", required=True,  help="Path to the report/ directory (contains templates/)")
    parser.add_argument("--output",     required=True,  help="Output HTML path")
    parser.add_argument("--vl_csv",     default=None,   help="Optional viral load CSV (columns: sample_id, date, vl)")
    parser.add_argument("--cd4_csv",    default=None,   help="Optional CD4 CSV (columns: sample_id, date, cd4_pct, cd4_count)")
    args = parser.parse_args()

    sys.path.insert(0, args.report_dir)
    from ContextBuilder.core import BLASTResult
    from ContextBuilder.assembly_info import get_length_with_flye_info
    from ContextBuilder.sierrapy import SierraPyResult

    sample_id  = args.sample
    sample_dir = args.sample_dir

    assemblies = get_length_with_flye_info(
        os.path.join(sample_dir, "medaka_consensus", "consensus.fasta"),
        os.path.join(sample_dir, "flye", "assembly_info.txt"),
    )

    core_nt_json    = os.path.join(sample_dir, "blast", "core_nt.blast.json")
    los_alamos_json = os.path.join(sample_dir, "blast", "los_alamos.blast.json")

    core_nt_blast_result = BLASTResult.read_json(core_nt_json)    if os.path.exists(core_nt_json)    else None
    hiv_blast_result     = BLASTResult.read_json(los_alamos_json) if os.path.exists(los_alamos_json) else None

    for i, item in enumerate(assemblies):
        core_nt_hits = core_nt_blast_result.get_hits(item["contig_id"], 6) if core_nt_blast_result else []
        hiv_hits     = hiv_blast_result.get_hits(item["contig_id"], 6)     if hiv_blast_result     else []
        hits_with_image = [dict(build_hiv_image_from_hit(h), **h) for h in hiv_hits]
        assemblies[i]["blast"] = {"core_nt": core_nt_hits, "los_alamos": hits_with_image}

    sierrapy_result = SierraPyResult.read_sierrapy_json(os.path.join(sample_dir, "sierrapy.json"))
    drug_resistance_panel = build_drug_resistance_panel(sierrapy_result, sample_id)

    vl_data  = get_vl(args.vl_csv, sample_id)   if args.vl_csv  else {"date": [], "vl": []}
    cd4_data = get_cd4(args.cd4_csv, sample_id) if args.cd4_csv else {"date": [], "cd4_pct": [], "cd4_count": []}

    jinja_data = {
        "sequencing_sample_id": sample_id,
        "assemblies": assemblies,
        "genome_length": GENOME_LENGTH,
        "drug_resistance_panel": drug_resistance_panel,
        "viral_load_history": {
            "timestamp": vl_data.get("date", []),
            "data":      vl_data.get("vl", []),
        },
        "cd4_count_history": {
            "timestamp": cd4_data.get("date", []),
            "percent":   {"data": cd4_data.get("cd4_pct", [])},
            "count":     {"data": cd4_data.get("cd4_count", [])},
        },
    }

    env = Environment(loader=FileSystemLoader(os.path.join(args.report_dir, "templates")))
    template = env.get_template("report_template.jinja2")

    with open(args.output, "w", encoding="utf-8") as f:
        f.write(template.render(jinja_data))

    print(f"Report written to {args.output}")


if __name__ == "__main__":
    main()
