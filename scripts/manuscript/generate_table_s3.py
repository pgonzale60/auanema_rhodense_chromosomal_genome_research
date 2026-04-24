def generate_table_s3():
    # 1. Define the mapping of TR Names to Family IDs
    tr_mapping = {
        "AT-1": "P351_F001",
        "AT-2": "P176_F001",
        "AC-I": "P348_F001",
        "AC-II": "P399_F001",
        "AC-III": "P167_F002",
        "AC-IV": "P332_F002",
        "AC-V": "P231_F004",
        "AC-VI": "P412_F002",
        "XT-1": "P347_F001",
        "XT-2": "P342_F001",
        "XR-1": "P348_F003",
        "XR-2": "P291_F001"
    }
    named_family_ids = set(tr_mapping.values())
    
    # 2. Load the TRF family summary
    summary_path = "analyses/genome_features/repeats/TRF/trf_family_summary.tsv"
    families = {}
    total_genome_tr_span = 0
    total_genome_tr_elim = 0
    
    with open(summary_path, 'r') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            data = dict(zip(header, parts))
            f_id = data['Family_ID']
            families[f_id] = data
            
            # Aggregate totals for ALL TRs
            total_genome_tr_span += int(data['Total_span_bp'])
            total_genome_tr_elim += int(data['GRS_span_bp'])
            
    # 3. Process each TR family in the requested order
    requested_order = ["AT-1", "AT-2", "AC-I", "AC-II", "AC-III", "AC-IV", "AC-V", "AC-VI", "XT-1", "XT-2", "XR-1", "XR-2"]
    
    results = []
    named_total_span = 0
    named_elim_span = 0
    
    for tr_name in requested_order:
        family_id = tr_mapping.get(tr_name)
        if family_id in families:
            f = families[family_id]
            total_span = int(f['Total_span_bp'])
            elim_span = int(f['GRS_span_bp'])
            ret_span = total_span - elim_span
            
            named_total_span += total_span
            named_elim_span += elim_span
            
            results.append({
                "TR Family": tr_name,
                "Consensus Monomer Length (bp)": f['Period'],
                "Total Span (bp)": total_span,
                "Eliminated Span (bp)": elim_span,
                "Retained Span (bp)": ret_span,
                "Percentage Eliminated": (elim_span / total_span) * 100 if total_span > 0 else 0
            })
    
    # 4. Calculate "Other"
    other_total_span = total_genome_tr_span - named_total_span
    other_elim_span = total_genome_tr_elim - named_elim_span
    other_ret_span = other_total_span - other_elim_span
    
    results.append({
        "TR Family": "Other families",
        "Consensus Monomer Length (bp)": "-",
        "Total Span (bp)": other_total_span,
        "Eliminated Span (bp)": other_elim_span,
        "Retained Span (bp)": other_ret_span,
        "Percentage Eliminated": (other_elim_span / other_total_span) * 100 if other_total_span > 0 else 0
    })
    
    # 5. Calculate "Total"
    total_ret_span = total_genome_tr_span - total_genome_tr_elim
    results.append({
        "TR Family": "**Total**",
        "Consensus Monomer Length (bp)": "-",
        "Total Span (bp)": total_genome_tr_span,
        "Eliminated Span (bp)": total_genome_tr_elim,
        "Retained Span (bp)": total_ret_span,
        "Percentage Eliminated": (total_genome_tr_elim / total_genome_tr_span) * 100 if total_genome_tr_span > 0 else 0
    })
    
    # 6. Print as Markdown Table
    print("| TR Family | Consensus Monomer Length (bp) | Total Span (bp) | Eliminated Span (bp) | Retained Span (bp) | Percentage Eliminated |")
    print("| :--- | :---: | :---: | :---: | :---: | :---: |")
    for r in results:
        label = r['TR Family']
        if label == "**Total**":
             print(f"| {label} | {r['Consensus Monomer Length (bp)']} | **{r['Total Span (bp)']:,.0f}** | **{r['Eliminated Span (bp)']:,.0f}** | **{r['Retained Span (bp)']:,.0f}** | **{r['Percentage Eliminated']:.2f}%** |")
        else:
             print(f"| {label} | {r['Consensus Monomer Length (bp)']} | {r['Total Span (bp)']:,.0f} | {r['Eliminated Span (bp)']:,.0f} | {r['Retained Span (bp)']:,.0f} | {r['Percentage Eliminated']:.2f}% |")

if __name__ == "__main__":
    generate_table_s3()
