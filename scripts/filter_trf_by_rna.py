#!/usr/bin/env python3
import gzip
import sys
from collections import defaultdict
from typing import Tuple, Dict, List

def parse_window_id(window_id: str) -> Tuple[str, int, int]:
    """Parse window ID into chromosome and coordinates."""
    chrom, coords = window_id.rsplit(':', 1)
    start, end = coords.split('-')
    return chrom, int(start), int(end)

def load_rnas(bed_file: str) -> Dict[str, List[Tuple[int, int]]]:
    """Load RNAs from BED file into a dictionary of intervals per chromosome."""
    rnas = defaultdict(list)
    with open(bed_file, 'r') as f:
        for line in f:
            chrom, start, end = line.strip().split('\t')
            rnas[chrom].append((int(start), int(end)))
    
    # Sort for easier searching later
    for chrom in rnas:
        rnas[chrom].sort()
    
    return rnas

def get_overlap(s1: int, e1: int, s2: int, e2: int) -> int:
    """Calculate overlap between two intervals."""
    return max(0, min(e1, e2) - max(s1, s2))

def is_tr_filtered(chrom: str, abs_start: int, abs_end: int, rnas: Dict[str, List[Tuple[int, int]]], threshold: float = 0.5) -> bool:
    """Check if TR overlaps with any RNA above the threshold."""
    if chrom not in rnas:
        return False
    
    tr_len = abs_end - abs_start
    if tr_len <= 0:
        return False
        
    # Find relevant RNAs (we could use binary search, but list is likely small per chromosome or fine for one pass)
    for rna_start, rna_end in rnas[chrom]:
        # Since rnas are sorted by start, we can optimize slightly
        if rna_start >= abs_end:
            break
        if rna_end <= abs_start:
            continue
            
        overlap = get_overlap(abs_start, abs_end, rna_start, rna_end)
        if overlap / tr_len > threshold:
            return True
            
    return False

def main():
    if len(sys.argv) < 4:
        print("Usage: python filter_trf_by_rna.py <input.trf.tsv.gz> <rnas.bed> <output.trf.tsv.gz>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    rna_file = sys.argv[2]
    output_file = sys.argv[3]
    
    print(f"Loading RNAs from {rna_file}...")
    rnas = load_rnas(rna_file)
    
    filtered_count = 0
    total_count = 0
    
    print(f"Filtering {input_file}...")
    with gzip.open(input_file, 'rt') as f_in, gzip.open(output_file, 'wt') as f_out:
        for line in f_in:
            fields = line.strip().split('\t')
            if len(fields) < 10:
                continue
            
            total_count += 1
            window_id = fields[0]
            chrom, window_start, _ = parse_window_id(window_id)
            
            rel_start = int(fields[1])
            rel_end = int(fields[2])
            abs_start = window_start + rel_start - 1
            abs_end = window_start + rel_end - 1
            
            if is_tr_filtered(chrom, abs_start, abs_end, rnas):
                filtered_count += 1
                continue
            
            f_out.write(line)
            
    print(f"Total TRs: {total_count}")
    print(f"Filtered out: {filtered_count}")
    print(f"Remaining: {total_count - filtered_count}")

if __name__ == "__main__":
    main()
