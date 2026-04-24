import gzip
import re
 
def load_grs():
    grs = {}
    with open("analyses/genome_features/elim_coords/nxAuaRhod1_1.GRS.bed", 'r') as f:
        for line in f:
            if not line.strip(): continue
            parts = line.strip().split('\t')
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            if chrom not in grs:
                grs[chrom] = []
            grs[chrom].append((start, end))
    return grs
 
def is_in_grs(chrom, start, end, grs):
    if chrom not in grs:
        return False
    for g_start, g_end in grs[chrom]:
        if start <= g_end and end >= g_start:
            return True
    return False
 
def generate_table_s4():
    """
    Generates Table S4: Summary of ncRNA genes present in eliminated DNA.
    """
    grs = load_grs()
    main_chroms = ["SUPER_1", "SUPER_2", "SUPER_3", "SUPER_4", "SUPER_5", "SUPER_6", "SUPER_X"]
    
    stats = {chrom: {"tRNA": 0, "rRNA": 0, "snRNA/Other": 0} for chrom in main_chroms}
    stats["Other Scaffolds"] = {"tRNA": 0, "rRNA": 0, "snRNA/Other": 0}
 
    # 1. tRNA
    with gzip.open("analyses/genome_features/repeats/tRNAscan/nxAuaRhod1_1.trnas.gff.gz", 'rt') as f:
        seen_genes = set()
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split('\t')
            match = re.search(r"Parent=[^0-9]+(\d+)_gen", parts[8])
            if match:
                gene_id = match.group(1)
                if gene_id in seen_genes: continue
                seen_genes.add(gene_id)
            
            if is_in_grs(parts[0], int(parts[3]), int(parts[4]), grs):
                target = parts[0] if parts[0] in main_chroms else "Other Scaffolds"
                stats[target]["tRNA"] += 1
 
    # 2. rRNA
    with gzip.open("analyses/genome_features/repeats/rnammer/nxAuaRhod1_1.rnammer.gff3.gz", 'rt') as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split('\t')
            if is_in_grs(parts[0], int(parts[3]), int(parts[4]), grs):
                target = parts[0] if parts[0] in main_chroms else "Other Scaffolds"
                stats[target]["rRNA"] += 1
 
    # 3. miscRNA
    with gzip.open("analyses/genome_features/repeats/infernal/nxAuaRhod1_1.infernal.gff.gz", 'rt') as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split('\t')
            rtype = parts[2]
            if "rRNA" in rtype or "tRNA" in rtype or "Protozoa" in rtype: continue
            if rtype in ["K_chan_RES", "RNase_MRP", "Fluoride", "GlsR7", "RAGATH-21", "snosnR60_Z15", "Histone3"]: continue
            
            if is_in_grs(parts[0], int(parts[3]), int(parts[4]), grs):
                target = parts[0] if parts[0] in main_chroms else "Other Scaffolds"
                stats[target]["snRNA/Other"] += 1
 
    print("# Table S4. ncRNA genes present in eliminated DNA in Auanema rhodense")
    print("\nSummary of non-coding RNA (ncRNA) features identified within the eliminated DNA regions (GRS and unassigned scaffolds).\n")
    print("| Chromosome | tRNA | rRNA | snRNA/Other | Total |")
    print("| :--- | :---: | :---: | :---: | :---: |")
    
    grand_total = 0
    for chrom in main_chroms + ["Other Scaffolds"]:
        s = stats[chrom]
        row_total = sum(s.values())
        grand_total += row_total
        label = chrom.replace("SUPER_", "Chr ") if "SUPER" in chrom else chrom
        print(f"| {label} | {s['tRNA']:,} | {s['rRNA']:,} | {s['snRNA/Other']:,} | **{row_total:,}** |")
    
    print(f"| **Total** | | | | **{grand_total:,}** |")
 
if __name__ == "__main__":
    generate_table_s4()
