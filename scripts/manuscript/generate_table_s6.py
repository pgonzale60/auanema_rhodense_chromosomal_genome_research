#!/usr/bin/env python3
import csv
import gzip
import sys
from pathlib import Path

def generate_table_s6():
    """
    Generates Table S6: Coordinates of Sequences for Elimination (SFEs) at the 28 main breakage sites
    and their accompanying redundant SFE motifs in eliminated DNA (GRS).
    Uses the FIMO callset (analyses/diminution/fimo_out/fimo.tsv.gz)
    for Motif-1 (TTTTGYTYAGGCTAAGCCTRWRCWAAAATW, 30 bp) filtered programmatically at the minimum canonical
    break-site score threshold.
    """
    project_root = Path(__file__).resolve().parent.parent.parent
    break_sites_path = project_root / "analyses/genome_features/elim_coords/nxAuaRhod1_1.break_sites.tsv"
    fimo_path = project_root / "analyses/diminution/fimo_out/fimo.tsv.gz"
    out_tsv_path = project_root / "report/supplementary tables/Table_S6_SFE_coordinates.tsv"

    target_motif = "WTTTWGYTHAGGCTKAGCCWAAACAAAAA"
    canonical_chrs = {f"SUPER_{i}" for i in range(1, 7)} | {"SUPER_X"}

    # 1. Load 28 canonical break sites in chromosomal order
    break_sites = []
    with open(break_sites_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            break_sites.append((row['chrom'], int(row['coordinate'])))

    # 2. Load positive-scoring FIMO hits for Motif-1 on canonical chromosomes
    m1_raw = []
    with gzip.open(fimo_path, 'rt') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for r in reader:
            if not r or not r.get('start') or r['motif_id'] != target_motif:
                continue
            chrom = r['sequence_name']
            if chrom not in canonical_chrs:
                continue
            score = float(r['score'])
            if score < 0:
                continue
            m1_raw.append({
                'chrom': chrom,
                'start': int(r['start']),
                'end': int(r['stop']),
                'strand': r['strand'],
                'score': score,
                'p_value': float(r['p-value']),
                'sequence': r['matched_sequence']
            })

    # 3. Deduplicate overlapping +/- strand hits at the same palindromic locus (keep highest score)
    m1_raw.sort(key=lambda x: -x['score'])
    m1_dedup = []
    for h in m1_raw:
        if not any(
            h['chrom'] == d['chrom'] and not (h['end'] < d['start'] or h['start'] > d['end'])
            for d in m1_dedup
        ):
            m1_dedup.append(h)

    # 4. Identify the main breakage site SFE at each of the 28 break sites
    bs_hits = {}
    for bc, bpos in break_sites:
        candidates = [
            h for h in m1_dedup
            if h['chrom'] == bc and (abs(h['start'] - bpos) <= 25 or abs(h['end'] - bpos) <= 25 or (h['start'] <= bpos <= h['end']))
        ]
        candidates.sort(key=lambda x: -x['score'])
        bs_hits[(bc, bpos)] = candidates[0]

    # Programmatically obtain the minimum break-site score threshold
    min_bs_score = min(h['score'] for h in bs_hits.values())
    print(f"Programmatic break-site threshold: {min_bs_score:.4f} bits", file=sys.stderr)

    # 5. Identify redundant motifs with score >= min_bs_score outside the 28 main break sites
    bs_hit_coords = {(h['chrom'], h['start'], h['end']) for h in bs_hits.values()}
    redundant_hits = [
        h for h in m1_dedup
        if h['score'] >= min_bs_score and (h['chrom'], h['start'], h['end']) not in bs_hit_coords
    ]
    redundant_hits.sort(key=lambda x: (x['chrom'], x['start']))

    # 6. Pair each redundant motif with its nearest main break site on the same chromosome
    paired_redundant = {k: [] for k in break_sites}
    for rh in redundant_hits:
        c = rh['chrom']
        chr_bs = [(bc, bpos) for bc, bpos in break_sites if bc == c]
        nearest_bs = min(chr_bs, key=lambda b: min(abs(rh['start'] - b[1]), abs(rh['end'] - b[1])))
        paired_redundant[nearest_bs].append(rh)

    # 7. Write TSV in the original side-by-side border & redundant motif format
    out_tsv_path.parent.mkdir(parents=True, exist_ok=True)
    row_count = 0
    with open(out_tsv_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow([
            "scaffold",
            "motif at main breakage site", "", "",
            "redundant motif", "", "",
            "distance between redundant motif and SFE at break site"
        ])
        writer.writerow([
            "",
            "start", "end", "score",
            "start", "end", "score",
            ""
        ])

        for bs_key in break_sites:
            bh = bs_hits[bs_key]
            reds = paired_redundant[bs_key]
            if not reds:
                writer.writerow([
                    bh['chrom'],
                    bh['start'], bh['end'], f"{bh['score']:.4f}",
                    "", "", "",
                    ""
                ])
                row_count += 1
            else:
                for rh in reds:
                    dist = rh['start'] - bh['end'] if rh['start'] > bh['end'] else bh['start'] - rh['end']
                    writer.writerow([
                        bh['chrom'],
                        bh['start'], bh['end'], f"{bh['score']:.4f}",
                        rh['start'], rh['end'], f"{rh['score']:.4f}",
                        dist
                    ])
                    row_count += 1

    print(f"Generated Table S6 at {out_tsv_path} ({row_count} rows across 28 break sites).", file=sys.stderr)

if __name__ == "__main__":
    generate_table_s6()
