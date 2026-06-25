#!/usr/bin/env python
import pandas as pd
from jinja2 import Environment, FileSystemLoader
from ContextBuilder.core import BLASTResult
from ContextBuilder.assembly_info import get_length_with_flye_info
from ContextBuilder.sierrapy import SierraPyResult
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
    dict(name="gag", start=790, end=2292, row=2, strand="+"),
    dict(name="gag-pol", start=2085, end=5096, row=3, strand="+"),
    dict(name="vif", start=5041, end=5619, row=2, strand="+"),
    dict(name="vpr", start=5559, end=5850, row=3, strand="+"),
    dict(name="tat", start=5831, end=6045, row=2, strand="+"),
    dict(name="tat", start=8379, end=8469, row=2, strand="+"),
    dict(name="rev", start=5970, end=6045, row=3, strand="+"),
    dict(name="rev", start=8379, end=8653, row=3, strand="+"),
    dict(name="vpu", start=6062, end=6310, row=4, strand="+"),
    dict(name="env", start=6225, end=8795, row=5, strand="+"),
    dict(name="nef", start=8797, end=9417, row=2, strand="+"),
]

def overlap(a1,a2,b1,b2):
    start = max(a1,b1)
    end = min(a2,b2)
    if start < end:
        return start,end
    return None


def bp_to_px(bp):
    return bp * PX_PER_BP

def arrow_forward(x1, x2, y1, y2):

    mid = (y1 + y2) / 2
    tip = 2

    return f"""
M{x1},{y1}
L{x2-tip},{y1}
L{x2},{mid}
L{x2-tip},{y2}
L{x1},{y2}
L{x1+tip},{mid}
Z
"""


def arrow_reverse(x1, x2, y1, y2):

    mid = (y1 + y2) / 2
    tip = 2

    return f"""
M{x1+tip},{y1}
L{x2},{y1}
L{x2-tip},{mid}
L{x2},{y2}
L{x1+tip},{y2}
L{x1},{mid}
Z
"""

def splice_ticks(start, end, row, step=300):

    y1,y2 = ROWS[row]
    ticks = []

    for bp in range(start, end, step):

        x = bp_to_px(bp)

        ticks.append({
            "x":x,
            "y1":y1+1,
            "y2":y2-1
        })

    return ticks

def build_backbone_segments(highlight_regions):
    '''
    highlight_regions use start,end,strand format, strand should be from query
    '''
    segments = []
    for start,end,strand in highlight_regions:
        if start > end:
            start, end = end, start
        x = bp_to_px(start)
        width = bp_to_px(end - start)
        segments.append({
            "x": x,
            "width": width,
            "strand": strand
        })

    return segments

def build_genes(highlight_regions):
    out = []
    for i,g in enumerate(HIV_GENES):
        x1 = bp_to_px(g["start"])
        x2 = bp_to_px(g["end"])
        y1,y2 = ROWS[g["row"]]

        if g["strand"] == "+":
            path = arrow_forward(x1,x2,y1,y2)
        else:
            path = arrow_reverse(x1,x2,y1,y2)

        segments = []
        for start,end in highlight_regions:
            if start > end:
                start, end = end, start
            ov = overlap(g["start"],g["end"],start,end)
            if ov:
                s,e = ov
                segments.append({
                    "x": bp_to_px(s),
                    "width": bp_to_px(e-s),
                    "y": y1,
                    "height": y2-y1
                })
        out.append({
            "id": f"gene{i}",
            "name":g["name"],
            "path":path,
            "segments":segments,
            "ticks":splice_ticks(g["start"],g["end"],g["row"])
        })

    return out

INTERPRETATION_SIR_MAP = {
    'susceptible': 'S', 'potential low-level resistance': 'S',
    'low-level resistance': 'I', 'intermediate resistance': 'I',
    'high-level resistance': 'R'
}
drug_display_to_fullname = {
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
    "TPV/r": "Tipranavir/ritonavir"
    
}

def build_hiv_image_from_hit(blast_record):
    hsps_subject = [(i[1], i[2]) for i in blast_record['hsps_subject']]
    hsps_query_strand = ['+' if i[0] == 'Plue' else '-' for i in blast_record['hsps_query']]
    hsps_subject_with_strand = [(*pair, n) for pair, n in zip(hsps_subject, hsps_query_strand)]
    return {
        'backbone_raw': hsps_subject_with_strand,
        'image_genes': build_genes(hsps_subject),
        'image_backbone_segments': build_backbone_segments(hsps_subject_with_strand),
        'image_genome_length': blast_record['subject_length']
    }

def get_vl(vl_csv, sample_id):
    _df = pd.read_csv(vl_csv)
    filtered = _df.loc[_df['sample_id'] == sample_id]
    return filtered[['date', 'vl']].to_dict('list')
def get_cd4(cd4_csv, sample_id):
    _df = pd.read_csv(cd4_csv)
    filtered = _df.loc[_df['sample_id'] == sample_id]
    return filtered[['date', 'cd4_pct', 'cd4_count']].to_dict('list')

SAMPLE_IDS = [
    '03A5', '03B4', '03B8', '03C0', '03C5', 
    '03D1', '0113', '0134', '0144', '0145',
    '0151', '0228', '0241', '0242', '0250',
    '0256', '0259', '0261', '0270', '0271',
    '0277', '0278', '0279', '0264', '0125',
]

for sample_id in SAMPLE_IDS:
    output_html_path=f'/mnt/sdb/WORKSPACE/ONT_HIV_Sequencing/workdir/2k_assembly/Merged_analysis/HIV64148_report/{sample_id}_HIV64148_report.html'
    print(f'Writing report to {output_html_path}')
    assemblies = get_length_with_flye_info(
        f'/mnt/sdb/WORKSPACE/ONT_HIV_Sequencing/workdir/2k_assembly/{sample_id}/medaka_consensus/consensus.fasta',
        f'/mnt/sdb/WORKSPACE/ONT_HIV_Sequencing/workdir/2k_assembly/{sample_id}/flye/assembly_info.txt'
    )
    core_nt_blast_result = BLASTResult.read_json(f'/mnt/sdb/WORKSPACE/ONT_HIV_Sequencing/workdir/2k_assembly/{sample_id}/blast/medaka_consensus.fa.core_nt.blast.json')
    hiv_blast_result = BLASTResult.read_json(f'/mnt/sdb/WORKSPACE/ONT_HIV_Sequencing/workdir/2k_assembly/{sample_id}/blast/medaka_consensus.fa.los_alamos.json')
    for i,item in enumerate(assemblies):
        core_nt_hits = core_nt_blast_result.get_hits(item['contig_id'], 6)
        hiv_hits = hiv_blast_result.get_hits(item['contig_id'], 6)
        hits_with_hiv_genome_image = [ dict(build_hiv_image_from_hit(item), **item) for item in hiv_hits ]

        assemblies[i]['blast'] = {
            'core_nt': core_nt_hits,
            'los_alamos': hits_with_hiv_genome_image,
        }

    sirerrapy_result = SierraPyResult.read_sierrapy_json(f'/mnt/sdb/WORKSPACE/ONT_HIV_Sequencing/workdir/2k_assembly/{sample_id}/sierrapy_result.0.json')
    data = sirerrapy_result.get_drug_resistance_panal()
    records=[]
    rows = []
    for contig, reports in data:           # data = your list
        for report in reports:
            version = report['version']['text']
            publish_date = report['version']['publishDate']
            gene = report['gene']['name']

            for drug_score in report['drugScores']:
                rows.append({
                    'result_sha256': sirerrapy_result.get_result_file_hash(),
                    "sample_id": sample_id,
                    "contig": contig,
                    "gene": gene,
                    "version": version,
                    "publish_date": publish_date,
                    "drug_class": drug_score['drugClass']['name'],
                    "drug": drug_score['drug']['name'],
                    "drug_abbr": drug_score['drug']['displayAbbr'],
                    "score": drug_score['score'],
                    "level": drug_score['level'],
                    "interpretation": drug_score['text'],
                    "SIR": INTERPRETATION_SIR_MAP[drug_score['text'].lower()]
                })
    records.extend(rows)

    df = pd.DataFrame(records)

    # index of row with max level per drug
    max_idx = df.groupby("drug_abbr")["level"].idxmax()

    summary = (
        df.groupby("drug_abbr")
        .agg(score=("score", "max"),
            level=("level", "max"),
            drug_class=("drug_class", "first"))
    )

    # interpretation corresponding to max level
    summary["interpretation"] = (
        df.loc[max_idx]
        .set_index("drug_abbr")["interpretation"]
    )

    # contigs where level > 0
    support_pos = (
        df[df["level"] > 0]
        .groupby("drug_abbr")["contig"]
        .unique()
    )

    # all contigs
    support_all = (
        df.groupby("drug_abbr")["contig"]
        .unique()
    )

    supported = support_pos.reindex(summary.index).combine_first(support_all)

    summary["supported_contig"] = supported.apply(list)
    summary["drug_full_name"] = summary.index.map(drug_display_to_fullname)
    summary['SIR'] = summary['interpretation'].map(lambda x: INTERPRETATION_SIR_MAP[x.lower()])
    drug_resistance_panel = summary.reset_index().to_dict("records")

    vl_data = get_vl('/mnt/sdb/WORKSPACE/ONT_HIV_Sequencing/report/vl.csv', sample_id)
    cd4_data = get_cd4('/mnt/sdb/WORKSPACE/ONT_HIV_Sequencing/report/cd4.csv', sample_id)

    jinja_data = {
        'hn': 'xxxxxxx',
        'firstname': 'xxxxx',
        'lastname': 'xxxxxxxxxxxxxx',
        'sex': 'x',
        'age': 'xx',
        "viral_load_history":{
            'timestamp': vl_data.get('date', []),
            'data': vl_data.get('vl', [])
        },
        "cd4_count_history": {
            'timestamp': cd4_data.get('date', []),
            'percent': { 'data': cd4_data.get('cd4_pct', []) },
            'count': { 'data': cd4_data.get('cd4_count', []) }
        },
        'drug_resistance_panel': drug_resistance_panel,
        "sequencing_sample_id": sample_id,
        'assemblies': assemblies,
        "genome_length":GENOME_LENGTH,
    }

    with open(output_html_path, "w", encoding="utf-8") as output_file:
        env = Environment(
            loader=FileSystemLoader("templates")
        )
        template = env.get_template("report_template.jinja2")
        output_file.write(template.render(jinja_data))
