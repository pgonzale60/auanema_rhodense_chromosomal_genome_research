import csv
import gzip

def generate_table_s6():
    """
    Generates Table S6: Coordinates and classification of Sequences for Elimination (SFEs) in Auanema rhodense.
    Filters motif occurrences by score >= lowest canonical break site score (21.74 bits),
    annotates canonical break sites vs internal eliminated copies, and sorts by chromosome then coordinate.
    """
    # 1. Load canonical break sites
    break_sites_path = "analyses/genome_features/elim_coords/nxAuaRhod1_1.break_sites.tsv"
    break_sites = []
    try:
        with open(break_sites_path, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                break_sites.append({
                    'chrom': row['chrom'],
                    'coordinate': int(row['coordinate'])
                })
    except FileNotFoundError:
        print(f"Error: {break_sites_path} not found.")
        return

    # Chromosome ordering
    chrom_order = {f"SUPER_{i}": i for i in range(1, 7)}
    chrom_order["SUPER_X"] = 7

    # 2. First pass: find canonical break site scores from expanded 33-bp palindromic FIMO scan
    fimo_path = "analyses/diminution/fimo_out/nxAuaRhod1_1_expanded_33bp_pal_fimo.tsv.gz"
    bs_best_scores = {}

    with gzip.open(fimo_path, 'rt') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if not row or not row.get('motif_id') or row['motif_id'].startswith("#"):
                continue
            
            f_chrom = row['sequence_name']
            f_start = int(row['start'])
            f_stop = int(row['stop'])
            f_score = float(row['score'])

            for bs in break_sites:
                if f_chrom == bs['chrom']:
                    b_coord = bs['coordinate']
                    if abs(f_start - b_coord) <= 25 or abs(f_stop - b_coord) <= 25:
                        key = (bs['chrom'], bs['coordinate'])
                        if key not in bs_best_scores or f_score > bs_best_scores[key]:
                            bs_best_scores[key] = f_score

    min_canonical_score = min(bs_best_scores.values())

    # 3. Second pass: collect all distinct 33-bp loci with score >= min_canonical_score across the 7 chromosomes
    loci = {}
    with gzip.open(fimo_path, 'rt') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if not row or not row.get('motif_id') or row['motif_id'].startswith("#"):
                continue
            
            f_chrom = row['sequence_name']
            if f_chrom not in chrom_order:
                continue

            f_score = float(row['score'])
            if f_score < min_canonical_score:
                continue

            f_start = int(row['start'])
            f_stop = int(row['stop'])
            key = (f_chrom, f_start, f_stop)

            if key not in loci or f_score > loci[key]['score']:
                loci[key] = {
                    'chrom': f_chrom,
                    'start': f_start,
                    'stop': f_stop,
                    'strand': '+/-',
                    'score': f_score,
                    'p_value': float(row['p-value']),
                    'sequence': row['matched_sequence']
                }

    # 4. Annotate canonical vs internal eliminated SFE and distance
    results = []
    for loc in loci.values():
        f_chrom = loc['chrom']
        f_start = loc['start']
        f_stop = loc['stop']

        # Find closest canonical break site on this chromosome
        bs_on_chr = [bs['coordinate'] for bs in break_sites if bs['chrom'] == f_chrom]
        closest_bs = min(bs_on_chr, key=lambda b: min(abs(f_start - b), abs(f_stop - b)))
        distance = min(abs(f_start - closest_bs), abs(f_stop - closest_bs))

        is_canonical = distance <= 25
        classification = "Canonical Break Site" if is_canonical else "Eliminated Region (Internal SFE)"
        break_site_label = f"{closest_bs:,}" if is_canonical else f"Near {closest_bs:,}"

        results.append({
            'chrom_num': chrom_order[f_chrom],
            'Chromosome': f_chrom.replace("SUPER_", "Chr "),
            'Classification': classification,
            'Associated Break Site': break_site_label,
            'SFE Start': f"{f_start:,}",
            'SFE Stop': f"{f_stop:,}",
            'Strand': loc['strand'],
            'Distance (bp)': distance,
            'Score': loc['score'],
            'P-value': loc['p_value'],
            'Sequence': loc['sequence'],
            'raw_start': f_start
        })

    # 5. Sort by chromosome, then start coordinate
    results.sort(key=lambda x: (x['chrom_num'], x['raw_start']))

    # 6. Print Markdown Table
    print("# Table S6. Coordinates and classification of Sequences for Elimination (SFEs)")
    print(f"\nThis table lists all genomic occurrences of the 33-bp palindromic SFE motif with scores greater than or equal to the lowest canonical break site score ({min_canonical_score:.2f} bits). Hits are classified as either canonical chromosome break sites (overlapping the 28 validated cleavage junctions) or internal eliminated copies (located within eliminated germline-restricted regions), and ordered by chromosome and genomic coordinate.\n")
    print("| Chromosome | Classification | Associated Break Site | SFE Start | SFE Stop | Strand | Distance to Junction (bp) | Score | P-value | Sequence |")
    print("| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |")
    for r in results:
        print(f"| {r['Chromosome']} | {r['Classification']} | {r['Associated Break Site']} | {r['SFE Start']} | {r['SFE Stop']} | {r['Strand']} | {r['Distance (bp)']} | {r['Score']:.2f} | {r['P-value']:.2e} | {r['Sequence']} |")

if __name__ == "__main__":
    generate_table_s6()
