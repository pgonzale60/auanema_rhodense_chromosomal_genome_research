#!/usr/bin/env python3
"""
Annotate TRF tandem repeats into families using period + consensus sequence clustering.

This script:
1. Loads TRF data and converts window-relative coordinates to absolute coordinates
2. Groups TRs by exact period size
3. Clusters TRs within each period group by consensus sequence similarity
4. Calculates total genomic span for each family (merging overlaps)
5. Resolves overlaps between families (larger families take precedence)
6. Outputs family annotations and summary statistics

Author: Claude Code
Date: 2026-01-21
"""

import gzip
import pandas as pd
import numpy as np
from collections import defaultdict, Counter
from dataclasses import dataclass, field
from typing import List, Dict, Set, Tuple
import sys
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform


@dataclass
class TR:
    """Tandem repeat instance."""
    window_id: str
    chrom: str
    window_start: int
    window_end: int
    rel_start: int
    rel_end: int
    abs_start: int
    abs_end: int
    period: int
    copies: float
    percent_match: int
    percent_indels: int
    score: int
    consensus: str
    family_id: str = None

    @property
    def length(self):
        return self.abs_end - self.abs_start

    def __repr__(self):
        return f"TR({self.chrom}:{self.abs_start}-{self.abs_end}, period={self.period}, family={self.family_id})"


@dataclass
class Family:
    """Tandem repeat family."""
    family_id: str
    period: int
    sub_family_num: int
    members: List[TR] = field(default_factory=list)
    total_span: int = 0
    grs_span: int = 0
    chromosomes: Set[str] = field(default_factory=set)

    @property
    def n_members(self):
        return len(self.members)

    @property
    def n_chromosomes(self):
        return len(self.chromosomes)

    @property
    def avg_copies(self):
        if not self.members:
            return 0
        return sum(tr.copies for tr in self.members) / len(self.members)

    @property
    def avg_identity(self):
        if not self.members:
            return 0
        return sum(tr.percent_match for tr in self.members) / len(self.members)


def parse_window_id(window_id: str) -> Tuple[str, int, int]:
    """Parse window ID into chromosome and coordinates."""
    chrom, coords = window_id.rsplit(':', 1)
    start, end = coords.split('-')
    return chrom, int(start), int(end)


def load_trf_data(filepath: str, test_chrom: str = None) -> List[TR]:
    """
    Load TRF data and convert to absolute coordinates.

    Args:
        filepath: Path to gzipped TRF TSV file
        test_chrom: Optional chromosome to filter for testing (e.g., 'SUPER_1')

    Returns:
        List of TR objects
    """
    print(f"Loading TRF data from {filepath}...")

    trs = []
    with gzip.open(filepath, 'rt') as f:
        for line_num, line in enumerate(f, 1):
            if line_num % 5000 == 0:
                print(f"  Processed {line_num} lines...", end='\r')

            fields = line.strip().split('\t')
            if len(fields) < 10:
                continue

            window_id = fields[0]
            chrom, window_start, window_end = parse_window_id(window_id)

            # Filter by chromosome if specified
            if test_chrom and chrom != test_chrom:
                continue

            rel_start = int(fields[1])
            rel_end = int(fields[2])
            abs_start = window_start + rel_start - 1  # Convert to 0-based
            abs_end = window_start + rel_end - 1

            tr = TR(
                window_id=window_id,
                chrom=chrom,
                window_start=window_start,
                window_end=window_end,
                rel_start=rel_start,
                rel_end=rel_end,
                abs_start=abs_start,
                abs_end=abs_end,
                period=int(fields[3]),
                copies=float(fields[4]),
                percent_match=int(fields[5]),
                percent_indels=int(fields[6]),
                score=int(fields[7]),
                consensus=get_canonical_lsr(fields[9]) # Strand-invariant LSR
            )
            trs.append(tr)

    print(f"\nLoaded {len(trs)} tandem repeats")
    return trs


def load_grs_data(filepath: str) -> List[Tuple[str, int, int]]:
    """Load GRS regions from BED file."""
    if not filepath:
        return []
    
    print(f"Loading GRS data from {filepath}...")
    grs_regions = []
    try:
        with open(filepath, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 3:
                    grs_regions.append((fields[0], int(fields[1]), int(fields[2])))
    except FileNotFoundError:
        print(f"Warning: GRS file {filepath} not found.")
    
    return grs_regions


def group_by_period(trs: List[TR]) -> Dict[int, List[TR]]:
    """Group TRs by exact period size."""
    print("\nGrouping TRs by period...")

    period_groups = defaultdict(list)
    for tr in trs:
        period_groups[tr.period].append(tr)

    print(f"Found {len(period_groups)} distinct period sizes")

    # Show top 10 periods by count
    period_counts = [(period, len(trs)) for period, trs in period_groups.items()]
    period_counts.sort(key=lambda x: x[1], reverse=True)
    print("\nTop 10 periods by count:")
    for period, count in period_counts[:10]:
        print(f"  Period {period}: {count} TRs")

    return period_groups


def reverse_complement(seq: str) -> str:
    """Get the reverse complement of a DNA sequence."""
    if not seq:
        return ""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'U': 'A'}
    return "".join(complement.get(base, base) for base in reversed(seq.upper()))


def get_lsr(seq: str) -> str:
    """Get the Lexicographical Smallest Rotation (LSR) of a sequence."""
    if not seq:
        return ""
    s = seq.upper()
    n = len(s)
    double_s = s + s
    best = s
    for i in range(1, n):
        current = double_s[i:i+n]
        if current < best:
            best = current
    return best


def get_canonical_lsr(seq: str) -> str:
    """Get the strand-invariant Lexicographical Smallest Rotation (LSR)."""
    if not seq:
        return ""
    s1 = get_lsr(seq)
    s2 = get_lsr(reverse_complement(seq))
    return s1 if s1 < s2 else s2


def get_kmer_profile(seq: str, k: int = 5) -> Counter:
    """Generate circular, RC-invariant k-mer profile for a sequence."""
    kmers = Counter()
    if not seq:
        return kmers
    seq_upper = seq.upper()
    n = len(seq_upper)
    if n < k:
        k_rc = reverse_complement(seq_upper)
        kmers[min(seq_upper, k_rc)] += 1
        return kmers

    extended_seq = seq_upper + seq_upper[:k-1]
    for i in range(n):
        kmer = extended_seq[i:i+k]
        kmer_rc = reverse_complement(kmer)
        kmers[min(kmer, kmer_rc)] += 1
    return kmers


def jaccard_similarity(profile1: Counter, profile2: Counter) -> float:
    """Calculate Jaccard similarity between two k-mer profiles, normalized for length."""
    if not profile1 or not profile2:
        return 0.0

    # Handle multimers by normalizing counts
    # To be strand-invariant, we could also combine p1 with its RC p1,
    # but get_canonical_lsr already picked a strand. 
    # However, k-mer profiles can still be slightly different.
    # We'll just use the canonical profiles.

    sum1 = sum(profile1.values())
    sum2 = sum(profile2.values())
    
    p1 = {k: v/sum1 for k, v in profile1.items()}
    p2 = {k: v/sum2 for k, v in profile2.items()}
    
    all_kmers = set(p1.keys()) | set(p2.keys())

    intersection = sum(min(p1.get(k, 0), p2.get(k, 0)) for k in all_kmers)
    union = sum(max(p1.get(k, 0), p2.get(k, 0)) for k in all_kmers)

    if union == 0:
        return 0.0

    return intersection / union


def cluster_by_consensus(trs: List[TR], period: int, similarity_threshold: float = 0.80,
                        k_mer_size: int = 5, min_family_size: int = 2) -> List[List[TR]]:
    """
    Cluster TRs by consensus sequence similarity using k-mer profiles.

    Args:
        trs: List of TRs with the same period
        period: Period size for these TRs
        similarity_threshold: Minimum similarity for clustering
        k_mer_size: Size of k-mers for profiling
        min_family_size: Minimum number of TRs per family

    Returns:
        List of TR clusters (families)
    """
    if len(trs) == 0:
        return []

    if len(trs) == 1:
        if min_family_size <= 1:
            return [trs]
        else:
            return []

    # Generate k-mer profiles
    profiles = [get_kmer_profile(tr.consensus, k_mer_size) for tr in trs]

    # Build similarity matrix
    n = len(trs)
    similarity_matrix = np.zeros((n, n))

    for i in range(n):
        similarity_matrix[i, i] = 1.0
        for j in range(i + 1, n):
            sim = jaccard_similarity(profiles[i], profiles[j])
            similarity_matrix[i, j] = sim
            similarity_matrix[j, i] = sim

    # Convert similarity to distance
    distance_matrix = 1 - similarity_matrix

    # Hierarchical clustering
    # Convert to condensed distance matrix for scipy
    condensed_dist = squareform(distance_matrix, checks=False)

    if n == 2:
        # Special case for 2 TRs
        if similarity_matrix[0, 1] >= similarity_threshold:
            return [trs] if len(trs) >= min_family_size else []
        else:
            return []

    # Perform hierarchical clustering
    linkage_matrix = linkage(condensed_dist, method='average')

    # Cut tree at distance threshold
    distance_threshold = 1 - similarity_threshold
    cluster_labels = fcluster(linkage_matrix, distance_threshold, criterion='distance')

    # Group TRs by cluster
    clusters = defaultdict(list)
    for tr, label in zip(trs, cluster_labels):
        clusters[label].append(tr)

    # Filter by minimum family size
    families = [trs for trs in clusters.values() if len(trs) >= min_family_size]

    return families


def merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Merge overlapping intervals."""
    if not intervals:
        return []

    # Sort by start position
    sorted_intervals = sorted(intervals)
    merged = [sorted_intervals[0]]

    for start, end in sorted_intervals[1:]:
        last_start, last_end = merged[-1]

        if start <= last_end:
            # Overlapping, merge
            merged[-1] = (last_start, max(last_end, end))
        else:
            # Non-overlapping, add new interval
            merged.append((start, end))

    return merged


def calculate_family_span(family: Family) -> int:
    """Calculate total genomic span of a family (merging overlaps)."""
    # Group by chromosome
    chrom_intervals = defaultdict(list)
    for tr in family.members:
        chrom_intervals[tr.chrom].append((tr.abs_start, tr.abs_end))
        family.chromosomes.add(tr.chrom)

    # Merge intervals per chromosome and sum
    total_span = 0
    for chrom, intervals in chrom_intervals.items():
        merged = merge_intervals(intervals)
        total_span += sum(end - start for start, end in merged)

    family.total_span = total_span
    return total_span


def calculate_grs_overlap(family: Family, grs_regions: List[Tuple[str, int, int]]):
    """Calculate how much of a family overlaps with GRS regions."""
    if not grs_regions:
        return 0

    # Simplify GRS into a dict of intervals per chrom
    grs_dict = defaultdict(list)
    for chrom, start, end in grs_regions:
        grs_dict[chrom].append((start, end))

    # Merge family intervals per chromosome
    chrom_intervals = defaultdict(list)
    for tr in family.members:
        chrom_intervals[tr.chrom].append((tr.abs_start, tr.abs_end))

    total_grs_overlap = 0
    for chrom, intervals in chrom_intervals.items():
        if chrom not in grs_dict:
            continue
            
        merged_fam = merge_intervals(intervals)
        merged_grs = merge_intervals(grs_dict[chrom])
        
        # Calculate intersection
        for f_start, f_end in merged_fam:
            for g_start, g_end in merged_grs:
                overlap_start = max(f_start, g_start)
                overlap_end = min(f_end, g_end)
                if overlap_start < overlap_end:
                    total_grs_overlap += (overlap_end - overlap_start)
    
    family.grs_span = total_grs_overlap
    return total_grs_overlap


def create_families(period_groups: Dict[int, List[TR]],
                   similarity_threshold: float = 0.80,
                   k_mer_size: int = 5,
                   min_family_size: int = 2) -> List[Family]:
    """Create families from period groups."""
    print(f"\nClustering by consensus sequence (threshold={similarity_threshold}, k={k_mer_size})...")

    all_families = []
    family_counter = 0

    # Process each period group
    for period in sorted(period_groups.keys()):
        trs = period_groups[period]
        print(f"\n  Processing period {period} ({len(trs)} TRs)...", end='')

        # Cluster by consensus
        clusters = cluster_by_consensus(trs, period, similarity_threshold,
                                       k_mer_size, min_family_size)

        print(f" -> {len(clusters)} families")

        # Create Family objects
        for sub_family_num, cluster in enumerate(clusters, 1):
            family_id = f"P{period}_F{sub_family_num:03d}"
            family = Family(
                family_id=family_id,
                period=period,
                sub_family_num=sub_family_num,
                members=cluster
            )

            # Assign family_id to TRs
            for tr in cluster:
                tr.family_id = family_id

            # Calculate span
            calculate_family_span(family)

            all_families.append(family)
            family_counter += 1

    print(f"\nCreated {len(all_families)} initial families")
    return all_families


def merge_similar_families(families: List[Family], similarity_threshold: float = 0.80, k_mer_size: int = 5) -> List[Family]:
    """Merge families across different periods if they are similar."""
    if not families:
        return []

    print(f"\nMerging similar families across periods (threshold={similarity_threshold})...")
    
    # Get representative profiles for each family
    family_profiles = []
    for fam in families:
        seq_counts = Counter(tr.consensus for tr in fam.members if tr.consensus)
        if not seq_counts:
            family_profiles.append(None)
            continue
        rep_seq = seq_counts.most_common(1)[0][0]
        family_profiles.append(get_kmer_profile(rep_seq, k_mer_size))

    # Build similarity matrix for families
    n = len(families)
    if n < 2:
        return families
        
    similarity_matrix = np.eye(n)
    for i in range(n):
        if family_profiles[i] is None: continue
        for j in range(i + 1, n):
            if family_profiles[j] is None: continue
            sim = jaccard_similarity(family_profiles[i], family_profiles[j])
            if sim > 0.5:
                print(f"    Similarity {families[i].family_id} - {families[j].family_id}: {sim:.4f}")
            similarity_matrix[i, j] = sim
            similarity_matrix[j, i] = sim

    # Cluster families (use a more relaxed threshold for cross-period)
    # 0.70 is often necessary for noisy TR data or multimers
    merge_threshold = similarity_threshold * 0.85 # Relaxed threshold
    distance_matrix = 1 - similarity_matrix
    condensed_dist = squareform(distance_matrix, checks=False)
    linkage_matrix = linkage(condensed_dist, method='average')
    cluster_labels = fcluster(linkage_matrix, 1 - merge_threshold, criterion='distance')
    
    # Merge
    merged_groups = defaultdict(list)
    for fam, label in zip(families, cluster_labels):
        merged_groups[label].append(fam)
        
    new_families = []
    for label, group in merged_groups.items():
        if len(group) == 1:
            new_families.append(group[0])
            continue
            
        # Create merged family: ID comes from the largest family
        group.sort(key=lambda f: f.total_span, reverse=True)
        primary = group[0]
        
        merged_members = []
        for fam in group:
            merged_members.extend(fam.members)
            
        merged_fam = Family(
            family_id=primary.family_id,
            period=primary.period,
            sub_family_num=primary.sub_family_num,
            members=merged_members
        )
        # Re-assign family_id to all members
        for tr in merged_members:
            tr.family_id = primary.family_id
            
        calculate_family_span(merged_fam)
        new_families.append(merged_fam)
        
    print(f"  Merged {n} families into {len(new_families)} families")
    return new_families


def resolve_overlaps_thorough(families: List[Family], chrom_sizes: Dict[str, int]) -> List[Family]:
    """
    Resolve overlaps at the base level across the entire genome.
    Each base is assigned to exactly one (or zero) TR family based on priority.
    
    Priority: Family Total Span (larger families first) > TR Score.
    """
    print("\nResolving overlaps thoroughly (base-level winner-take-all)...")

    # Collect all TRs across all families
    all_trs_to_resolve = []
    for fam in families:
        for tr in fam.members:
            # Attach family_id and priority score
            all_trs_to_resolve.append({
                'tr': tr,
                'family_id': fam.family_id,
                'priority': (fam.total_span, tr.score)
            })

    # Sort all TRs by priority (descending)
    all_trs_to_resolve.sort(key=lambda x: x['priority'], reverse=True)

    # Dictionary to store resolved blocks per chromosome
    # Chrom -> List of {family_id, start, end, ...metrics}
    resolved_blocks = defaultdict(list)

    # Process each chromosome
    for chrom in sorted(chrom_sizes.keys()):
        size = chrom_sizes[chrom]
        print(f"  Processing {chrom} ({size/1e6:.2f} Mb)...")
        
        # Initialize winner array (0 means no family)
        # We use a mapping for family IDs to fit in uint16 if possible
        # but for safety and speed let's just use a list of families for this chrom
        chrom_trs = [x for x in all_trs_to_resolve if x['tr'].chrom == chrom]
        if not chrom_trs:
            continue

        # Map family_id to a small integer for the array
        unique_fams = sorted(list(set(x['family_id'] for x in chrom_trs)))
        fam_to_idx = {fam_id: i + 1 for i, fam_id in enumerate(unique_fams)}
        idx_to_fam = {i + 1: fam_id for i, fam_id in enumerate(unique_fams)}

        # Create winner array
        # Use uint16 if possible (up to 65k families)
        winner_arr = np.zeros(size + 1, dtype=np.uint16)

        # Assign winners (in reverse priority so the highest priority overwrites)
        # Wait, if we process from highest priority, we only fill if 0.
        # Let's process from LOWEST priority and overwrite.
        for tr_data in reversed(chrom_trs):
            tr = tr_data['tr']
            idx = fam_to_idx[tr_data['family_id']]
            # Clip to chrom size
            start = max(0, tr.abs_start)
            end = min(size, tr.abs_end)
            if start < end:
                winner_arr[start:end] = idx

        # Convert back to contiguous blocks
        if not np.any(winner_arr > 0):
            continue

        # Find transitions
        # Find where winner_arr changes
        diffs = np.diff(winner_arr)
        trans_indices = np.where(diffs != 0)[0] + 1
        
        # Include start and end
        boundaries = [0] + trans_indices.tolist() + [len(winner_arr) - 1]
        
        for i in range(len(boundaries) - 1):
            start = boundaries[i]
            end = boundaries[i+1]
            val = winner_arr[start]
            
            if val > 0:
                family_id = idx_to_fam[val]
                
                # We need metrics for this block. Ideally, we average from the 
                # overlapping original TRs of the SAME family.
                # To keep it simple, we'll find the TRs of this family that cover this range.
                block_trs = [x['tr'] for x in chrom_trs if x['family_id'] == family_id and 
                             x['tr'].abs_start < end and x['tr'].abs_end > start]
                
                if block_trs:
                    avg_period = sum(t.period for t in block_trs) / len(block_trs)
                    avg_copies = sum(t.copies for t in block_trs) / len(block_trs)
                    avg_identity = sum(t.percent_match for t in block_trs) / len(block_trs)
                    consensus = block_trs[0].consensus # Use any for now
                    
                    resolved_blocks[chrom].append({
                        'Family_ID': family_id,
                        'Chrom': chrom,
                        'Start': start,
                        'End': end,
                        'Length': end - start,
                        'Period': round(avg_period),
                        'Copies': round(avg_copies, 1),
                        'Identity': round(avg_identity),
                        'Consensus_seq': consensus
                    })

    # Recreate Family objects from resolved blocks
    new_families = []
    final_family_ids = sorted(list(set(b['Family_ID'] for blocks in resolved_blocks.values() for b in blocks)))
    
    # Track which family members are now these blocks
    for fam_id in final_family_ids:
        # Re-parse ID for construction
        period_str, sub_family_str = fam_id.split('_')
        period = int(period_str[1:])
        sub_family_num = int(sub_family_str[1:])
        
        fam_blocks = [b for chrom_blocks in resolved_blocks.values() for b in chrom_blocks if b['Family_ID'] == fam_id]
        
        # Convert blocks to TR objects (mocking members)
        mock_members = []
        for b in fam_blocks:
            mock_members.append(TR(
                window_id=f"{b['Chrom']}:{b['Start']}-{b['End']}",
                chrom=b['Chrom'],
                window_start=b['Start'],
                window_end=b['End'],
                rel_start=1,
                rel_end=b['Length'],
                abs_start=b['Start'],
                abs_end=b['End'],
                period=b['Period'],
                copies=b['Copies'],
                percent_match=b['Identity'],
                percent_indels=0,
                score=0, # Not strictly needed for summary
                consensus=b['Consensus_seq'],
                family_id=fam_id
            ))
            
        new_fam = Family(
            family_id=fam_id,
            period=period,
            sub_family_num=sub_family_num,
            members=mock_members
        )
        calculate_family_span(new_fam)
        new_families.append(new_fam)

    print(f"  Resolved into {sum(len(b) for b in resolved_blocks.values())} non-overlapping blocks")
    return new_families


def generate_outputs(families: List[Family], conflicts: List, output_dir: str):
    """Generate all output files."""
    print("\nGenerating output files...")

    # Sort families by total span
    sorted_families = sorted(families, key=lambda f: f.total_span, reverse=True)

    # Mark top 30
    for i, family in enumerate(sorted_families):
        family.top30 = i < 30

    # Output 1: Family summary
    summary_file = f"{output_dir}/trf_family_summary.tsv"
    with open(summary_file, 'w') as f:
        f.write("Family_ID\tPeriod\tSub_family\tN_members\tN_chromosomes\tChromosomes\t"
                "Total_span_bp\tTotal_span_Mb\tGRS_span_bp\tAvg_copies\tAvg_identity\tTop30\n")

        for family in sorted_families:
            chroms = ','.join(sorted(family.chromosomes))
            f.write(f"{family.family_id}\t{family.period}\t{family.sub_family_num}\t"
                   f"{family.n_members}\t{family.n_chromosomes}\t{chroms}\t"
                   f"{family.total_span}\t{family.total_span/1e6:.2f}\t"
                   f"{family.grs_span}\t"
                   f"{family.avg_copies:.1f}\t{family.avg_identity:.1f}\t"
                   f"{family.top30}\n")

    print(f"  Wrote {summary_file}")

    # Output 2: Family members
    members_file = f"{output_dir}/trf_family_members.tsv"
    with open(members_file, 'w') as f:
        f.write("Family_ID\tChrom\tStart\tEnd\tLength\tPeriod\t"
                "Copies\tIdentity\tConsensus_seq\n")

        for family in sorted_families:
            for tr in sorted(family.members, key=lambda t: (t.chrom, t.abs_start)):
                f.write(f"{family.family_id}\t{tr.chrom}\t{tr.abs_start}\t{tr.abs_end}\t"
                       f"{tr.length}\t{tr.period}\t{tr.copies:.1f}\t"
                       f"{tr.percent_match}\t{tr.consensus}\n")

    print(f"  Wrote {members_file}")

    # Output 3: Overlap conflicts
    conflicts_file = f"{output_dir}/trf_overlap_conflicts.tsv"
    with open(conflicts_file, 'w') as f:
        f.write("Chrom\tPosition\tTR_count\tAssigned_family\t"
                "Competing_families\tResolution\n")

        for conflict in conflicts:
            f.write(f"{conflict['chrom']}\t{conflict['position']}\t"
                   f"{conflict['tr_count']}\t{conflict['assigned_family']}\t"
                   f"{conflict['competing_families']}\t{conflict['resolution']}\n")

    print(f"  Wrote {conflicts_file}")

    # Output 4: Top 30 families summary
    top30_file = f"{output_dir}/trf_top30_families_summary.tsv"
    with open(top30_file, 'w') as f:
        f.write("Rank\tFamily_ID\tPeriod\tSub_family\tN_members\tN_chromosomes\t"
                "Chromosomes\tTotal_span_bp\tTotal_span_Mb\tAvg_copies\tAvg_identity\n")

        for rank, family in enumerate(sorted_families[:30], 1):
            chroms = ','.join(sorted(family.chromosomes))
            f.write(f"{rank}\t{family.family_id}\t{family.period}\t{family.sub_family_num}\t"
                   f"{family.n_members}\t{family.n_chromosomes}\t{chroms}\t"
                   f"{family.total_span}\t{family.total_span/1e6:.2f}\t"
                   f"{family.avg_copies:.1f}\t{family.avg_identity:.1f}\n")

    print(f"  Wrote {top30_file}")

    # Print summary statistics
    print("\n" + "="*60)
    print("SUMMARY STATISTICS")
    print("="*60)

    total_trs = sum(f.n_members for f in families)
    total_span = sum(f.total_span for f in families)
    top30_span = sum(f.total_span for f in sorted_families[:30])

    print(f"Total TRs: {total_trs:,}")
    print(f"Total families: {len(families)}")
    print(f"Total genomic span: {total_span:,} bp ({total_span/1e6:.2f} Mb)")
    print(f"Top 30 families span: {top30_span:,} bp ({top30_span/1e6:.2f} Mb)")
    print(f"Top 30 coverage: {100*top30_span/total_span:.1f}%")
    print(f"Overlapping regions: {len(conflicts):,}")

    print("\nTop 10 families by span:")
    for i, family in enumerate(sorted_families[:10], 1):
        print(f"  {i}. {family.family_id}: {family.total_span/1e6:.2f} Mb "
              f"({family.n_members} members, {family.n_chromosomes} chroms)")

    print("="*60)


def main():
    """Main execution."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Annotate TRF tandem repeats into families"
    )
    parser.add_argument(
        '--input',
        default='analyses/genome_features/repeats/TRF/windowed100K_nxAuaRhod1_1.trf.tsv.gz',
        help='Input TRF file (gzipped TSV)'
    )
    parser.add_argument(
        '--output-dir',
        default='analyses/genome_features/repeats/TRF',
        help='Output directory for results'
    )
    parser.add_argument(
        '--similarity-threshold',
        type=float,
        default=0.50,
        help='Similarity threshold for consensus clustering (default: 0.50)'
    )
    parser.add_argument(
        '--k-mer-size',
        type=int,
        default=5,
        help='K-mer size for sequence profiling (default: 5)'
    )
    parser.add_argument(
        '--overlap-threshold',
        type=float,
        default=0.50,
        help='Overlap threshold for reassignment (default: 0.50)'
    )
    parser.add_argument(
        '--min-family-size',
        type=int,
        default=2,
        help='Minimum TRs per family (default: 2)'
    )
    parser.add_argument(
        '--test-chrom',
        default=None,
        help='Test on single chromosome (e.g., SUPER_1)'
    )
    parser.add_argument(
        '--grs-bed',
        default='analyses/genome_features/elim_coords/nxAuaRhod1_1.GRS.bed',
        help='Input GRS BED file'
    )

    args = parser.parse_args()

    print("="*60)
    print("TRF FAMILY ANNOTATION")
    print("="*60)
    print(f"Input: {args.input}")
    print(f"Output directory: {args.output_dir}")
    print(f"Similarity threshold: {args.similarity_threshold}")
    print(f"K-mer size: {args.k_mer_size}")
    print(f"Overlap threshold: {args.overlap_threshold}")
    print(f"Min family size: {args.min_family_size}")
    if args.test_chrom:
        print(f"Test chromosome: {args.test_chrom}")
    print("="*60)

    # Phase 1: Load data
    trs = load_trf_data(args.input, args.test_chrom)

    # Phase 2: Group by period
    period_groups = group_by_period(trs)

    # Phase 3: Create families (includes consensus clustering)
    families = create_families(
        period_groups,
        similarity_threshold=args.similarity_threshold,
        k_mer_size=args.k_mer_size,
        min_family_size=args.min_family_size
    )
    
    # NEW: Merge similar families across periods
    families = merge_similar_families(
        families, 
        similarity_threshold=args.similarity_threshold,
        k_mer_size=args.k_mer_size
    )

    # Phase 4: Resolve overlaps thoroughly
    # Get chromosome sizes for base-level resolution
    chrom_sizes = {}
    for chromosome in sorted(list(set(tr.chrom for tr in trs))):
        chrom_trs = [t for t in trs if t.chrom == chromosome]
        if chrom_trs:
            chrom_sizes[chromosome] = max(t.abs_end for t in chrom_trs)

    final_families = resolve_overlaps_thorough(families, chrom_sizes)
    
    # Dummy conflicts for output compatibility
    conflicts = []

    # NEW: Calculate GRS overlap
    grs_regions = load_grs_data(args.grs_bed)
    for family in final_families:
        calculate_grs_overlap(family, grs_regions)

    # Phase 6: Generate outputs
    generate_outputs(final_families, conflicts, args.output_dir)

    print("\nDone!")


if __name__ == '__main__':
    main()
