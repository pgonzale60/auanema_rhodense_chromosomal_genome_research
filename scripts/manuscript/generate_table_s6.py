import csv
import gzip
 
def generate_table_s6():
    """
    Generates Table S6: Coordinates of Sequences for Elimination (SFEs) at the breakage sites.
    """
    # 1. Load break sites
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
 
    # 2. Load and filter FIMO results
    fimo_path = "analyses/diminution/fimo_out/fimo.tsv.gz"
    results = []
    window = 500
    
    with gzip.open(fimo_path, 'rt') as f:
        # FIMO header: motif_id motif_alt_id sequence_name start stop strand score p-value q-value matched_sequence
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if not row or row['motif_id'].startswith("#"):
                continue
            
            f_chrom = row['sequence_name']
            f_start = int(row['start'])
            f_stop = int(row['stop'])
            
            # Check proximity to any break site
            for bs in break_sites:
                if f_chrom == bs['chrom']:
                    b_coord = bs['coordinate']
                    if (f_start >= b_coord - window and f_start <= b_coord + window) or \
                       (f_stop >= b_coord - window and f_stop <= b_coord + window):
                        
                        distance = min(abs(f_start - b_coord), abs(f_stop - b_coord))
                        results.append({
                            'Chromosome': f_chrom.replace("SUPER_", "Chr "),
                            'Break Site': f"{b_coord:,}",
                            'SFE Start': f"{f_start:,}",
                            'SFE Stop': f"{f_stop:,}",
                            'Strand': row['strand'],
                            'Distance': distance,
                            'Score': float(row['score']),
                            'P-value': float(row['p-value']),
                            'Sequence': row['matched_sequence']
                        })
 
    # Sort results
    results.sort(key=lambda x: (x['Chromosome'], x['Distance']))
 
    # Print as Markdown Table
    print("# Table S6. Coordinates of Sequences for Elimination (SFEs) at the breakage sites")
    print("\nThis table lists the motif occurrences (SFEs) found within 500 bp of the 28 internal chromosome breakage sites.\n")
    print("| Chromosome | Break Site | SFE Start | SFE Stop | Strand | Distance (bp) | Score | P-value | Sequence |")
    print("| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |")
    for r in results:
        print(f"| {r['Chromosome']} | {r['Break Site']} | {r['SFE Start']} | {r['SFE Stop']} | {r['Strand']} | {r['Distance']} | {r['Score']:.2f} | {r['P-value']:.2e} | {r['Sequence']} |")
 
if __name__ == "__main__":
    generate_table_s6()
