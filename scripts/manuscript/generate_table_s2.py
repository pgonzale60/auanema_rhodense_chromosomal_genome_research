import csv

def generate_table_s2():
    # 1. Load FAI
    fai_path = "nxAuaRhod1_1.primary.fa.gz.fai"
    fai = {}
    with open(fai_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            fai[parts[0]] = int(parts[1])
    
    # 2. Load core bed
    core_bed_path = "analyses/diminution/nxAuaRhod1_1.core.bed"
    retained_per_chr = {}
    with open(core_bed_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            retained_per_chr[chrom] = retained_per_chr.get(chrom, 0) + (end - start)
    
    # 3. Define main chromosomes
    main_chroms = ["SUPER_1", "SUPER_2", "SUPER_3", "SUPER_4", "SUPER_5", "SUPER_6", "SUPER_X"]
    
    results = []
    
    # 4. Process main chromosomes and their unloc scaffolds
    processed_seqs = set(main_chroms)
    for chrom in main_chroms:
        # Length of the main scaffold
        main_len = fai[chrom]
        
        # Length of unloc scaffolds
        unloc_prefix = f"{chrom}_unloc_"
        unloc_len = 0
        for seqid, length in fai.items():
            if seqid.startswith(unloc_prefix):
                unloc_len += length
                processed_seqs.add(seqid)
        
        germline_len = main_len + unloc_len
        retained_len = retained_per_chr.get(chrom, 0)
        eliminated_len = (main_len - retained_len) + unloc_len
        
        results.append({
            "Chromosome": chrom.replace("SUPER_", "Chr "),
            "Germline Length (bp)": germline_len,
            "Retained Somatic Length (bp)": retained_len,
            "Eliminated Length (bp)": eliminated_len,
            "Percentage Eliminated": (eliminated_len / germline_len) * 100
        })
    
    # 5. Process Shrapnel
    shrapnel_germline_len = 0
    for seqid, length in fai.items():
        if seqid in processed_seqs:
            continue
        if "MT" in seqid:
            continue
        shrapnel_germline_len += length
    
    results.append({
        "Chromosome": "Shrapnel (unassigned)",
        "Germline Length (bp)": shrapnel_germline_len,
        "Retained Somatic Length (bp)": 0,
        "Eliminated Length (bp)": shrapnel_germline_len,
        "Percentage Eliminated": 100.0
    })
    
    # 6. Totals
    total_germline = sum(r["Germline Length (bp)"] for r in results)
    total_retained = sum(r["Retained Somatic Length (bp)"] for r in results)
    total_eliminated = sum(r["Eliminated Length (bp)"] for r in results)
    
    results.append({
        "Chromosome": "Total",
        "Germline Length (bp)": total_germline,
        "Retained Somatic Length (bp)": total_retained,
        "Eliminated Length (bp)": total_eliminated,
        "Percentage Eliminated": (total_eliminated / total_germline) * 100
    })
    
    # Print as Markdown Table
    print("| Chromosome | Germline Length (bp) | Retained Somatic Length (bp) | Eliminated Length (bp) | Percentage Eliminated |")
    print("| :--- | :---: | :---: | :---: | :---: |")
    for r in results:
        print(f"| {r['Chromosome']} | {r['Germline Length (bp)']:,.0f} | {r['Retained Somatic Length (bp)']:,.0f} | {r['Eliminated Length (bp)']:,.0f} | {r['Percentage Eliminated']:.2f}% |")

if __name__ == "__main__":
    generate_table_s2()
