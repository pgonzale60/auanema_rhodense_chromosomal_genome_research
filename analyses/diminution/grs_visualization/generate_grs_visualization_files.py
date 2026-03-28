#!/usr/bin/env python3
"""
Generate visualization files for GRS (germline-restricted sequence) coverage analysis.

Creates three complementary data files:
1. Windowed coverage across GRS region with flanking context (BED format)
2. Per-base coverage at GRS boundaries (BED format)
3. Telomeric breakpoint positions from MILTEL analysis (TSV format with detailed attributes)

The telomeric breakpoints file includes rich annotation from MILTEL:
- breakage_score: ratio of telomere-bearing reads to total coverage
- num_telomere_reads: calculated number of telomere-bearing reads
- average_coverage: total read coverage at position
- orientation: left or right of junction
- telomere_sense: + (sense) or - (reverse complement)
- telomere_gap_avg: mean gap from telomere start to clipping junction
- telomere_senses_dist: distribution of senses across reads
- telomere_gaps_dist: distribution of gap distances across reads

Author: Claude Code
Date: 2026-02-06
"""

import gzip
import argparse
import statistics
from typing import List, Tuple, Dict
from pathlib import Path


# GRS Region Constants
GRS_START = 10706804
GRS_END = 10913775
GRS_LENGTH = GRS_END - GRS_START
FLANK_SIZE = int(GRS_LENGTH * 0.2)  # 20% flanking regions (41,394 bp)
WINDOW_SIZE = 100
BORDER_BP = 200

# Calculated region boundaries
REGION_START = GRS_START - FLANK_SIZE  # 10665410
REGION_END = GRS_END + FLANK_SIZE      # 10955169
LEFT_BORDER_START = GRS_START - BORDER_BP  # 10706784
LEFT_BORDER_END = GRS_START + BORDER_BP    # 10706824
RIGHT_BORDER_START = GRS_END - BORDER_BP   # 10913755
RIGHT_BORDER_END = GRS_END + BORDER_BP     # 10913795


def read_mosdepth_bed(bed_path: str, chromosome: str, start: int, end: int) -> List[Tuple[str, int, int, float]]:
    """
    Read mosdepth per-base BED file and filter for specific region.

    Args:
        bed_path: Path to gzipped BED file
        chromosome: Chromosome name to filter
        start: Start position (inclusive)
        end: End position (exclusive)

    Returns:
        List of tuples (chr, start, end, coverage)
    """
    entries = []

    with gzip.open(bed_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 4:
                continue

            chrom = fields[0]
            pos_start = int(fields[1])
            pos_end = int(fields[2])
            coverage = float(fields[3])

            # Filter for target chromosome and region
            if chrom == chromosome and pos_start >= start and pos_end <= end:
                entries.append((chrom, pos_start, pos_end, coverage))

    return entries


def generate_windowed_coverage(coverage_path: str, output_path: str, chromosome: str = "SUPER_X"):
    """
    Generate windowed coverage file (100bp windows) for GRS region + 20% flanks.

    Args:
        coverage_path: Path to mosdepth per-base BED file
        output_path: Output BED file path
        chromosome: Chromosome name
    """
    print(f"Generating windowed coverage for {chromosome}:{REGION_START}-{REGION_END}")
    print(f"  GRS region: {GRS_START}-{GRS_END} ({GRS_LENGTH} bp)")
    print(f"  Flanking: {FLANK_SIZE} bp each side")
    print(f"  Window size: {WINDOW_SIZE} bp")

    # Read per-base coverage
    coverage_data = read_mosdepth_bed(coverage_path, chromosome, REGION_START, REGION_END)

    if not coverage_data:
        print(f"WARNING: No coverage data found for {chromosome}:{REGION_START}-{REGION_END}")
        return

    print(f"  Read {len(coverage_data)} per-base entries")

    # Create coverage dictionary for fast lookup
    coverage_dict = {}
    for chrom, start, end, cov in coverage_data:
        for pos in range(start, end):
            coverage_dict[pos] = cov

    # Generate windowed coverage
    windows = []
    for window_start in range(REGION_START, REGION_END, WINDOW_SIZE):
        window_end = min(window_start + WINDOW_SIZE, REGION_END)

        # Collect coverage values in this window
        window_coverages = []
        for pos in range(window_start, window_end):
            if pos in coverage_dict:
                window_coverages.append(coverage_dict[pos])

        # Calculate mean coverage
        if window_coverages:
            mean_cov = statistics.mean(window_coverages)
            windows.append((chromosome, window_start, window_end, mean_cov))

    # Write output
    with open(output_path, 'w') as out:
        for chrom, start, end, cov in windows:
            out.write(f"{chrom}\t{start}\t{end}\t{cov:.2f}\n")

    print(f"  Wrote {len(windows)} windows to {output_path}")


def generate_border_coverage(coverage_path: str, output_path: str, chromosome: str = "SUPER_X"):
    """
    Generate per-base coverage at GRS boundaries (±20bp).

    Args:
        coverage_path: Path to mosdepth per-base BED file
        output_path: Output BED file path
        chromosome: Chromosome name
    """
    print(f"\nGenerating border coverage for GRS boundaries")
    print(f"  Left border: {LEFT_BORDER_START}-{LEFT_BORDER_END}")
    print(f"  Right border: {RIGHT_BORDER_START}-{RIGHT_BORDER_END}")

    # Read left border
    left_data = read_mosdepth_bed(coverage_path, chromosome, LEFT_BORDER_START, LEFT_BORDER_END)

    # Read right border
    right_data = read_mosdepth_bed(coverage_path, chromosome, RIGHT_BORDER_START, RIGHT_BORDER_END)

    # Combine and sort
    all_borders = sorted(left_data + right_data, key=lambda x: x[1])

    print(f"  Left border: {len(left_data)} entries")
    print(f"  Right border: {len(right_data)} entries")

    # Write output
    with open(output_path, 'w') as out:
        for chrom, start, end, cov in all_borders:
            out.write(f"{chrom}\t{start}\t{end}\t{cov:.2f}\n")

    print(f"  Wrote {len(all_borders)} entries to {output_path}")


def parse_gff3_attributes(attr_string: str) -> Dict[str, str]:
    """
    Parse GFF3 attributes field into dictionary.

    Args:
        attr_string: GFF3 attributes string (key=value pairs)

    Returns:
        Dictionary of attributes
    """
    attrs = {}
    for pair in attr_string.split(';'):
        if '=' in pair:
            key, value = pair.split('=', 1)
            attrs[key.strip()] = value.strip()
    return attrs


def generate_telomeric_breakpoints(gff_path: str, output_path: str, chromosome: str = "SUPER_X"):
    """
    Extract telomeric breakpoint positions from MILTEL GFF3 file.

    MILTEL produces chromosome_breakage_site annotations with rich attributes:
    - score: breakage score = (telomere-bearing reads) / (total coverage)
    - orientation: left or right of junction
    - telomere_sense: + (sense) or - (reverse complement)
    - telomere_gap_average: mean gap from telomere start to clipping junction
    - telomere_senses: distribution of senses (e.g., "-*326" = 326 reads with -)
    - telomere_gaps: distribution of gap distances (e.g., "0*317 15*2")
    - average_coverage: total read coverage at position

    Args:
        gff_path: Path to MILTEL GFF3 file
        output_path: Output TSV file path
        chromosome: Chromosome name
    """
    print(f"\nGenerating telomeric breakpoint positions from {gff_path}")

    breakpoints = []

    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom = fields[0]
            feature_type = fields[2]
            position = int(fields[3])  # start = end for point features
            score = fields[5]
            attributes = parse_gff3_attributes(fields[8])

            # Filter for chromosome_breakage_site features on target chromosome
            if chrom == chromosome and feature_type == "chromosome_breakage_site":
                # Extract all meaningful attributes
                orientation = attributes.get('orientation', 'unknown')
                telomere_sense = attributes.get('telomere_sense', 'unknown')
                telomere_gap_avg = attributes.get('telomere_gap_average', '0')
                average_coverage = attributes.get('average_coverage', '0')
                telomere_senses = attributes.get('telomere_senses', '')
                telomere_gaps = attributes.get('telomere_gaps', '')

                # Convert score to float (breakage score)
                try:
                    breakage_score = float(score)
                except ValueError:
                    breakage_score = 0.0

                # Calculate number of telomere-bearing reads
                try:
                    avg_cov = float(average_coverage)
                    num_telomere_reads = int(breakage_score * avg_cov)
                except (ValueError, ZeroDivisionError):
                    avg_cov = 0.0
                    num_telomere_reads = 0

                # Parse telomere gap average
                try:
                    gap_avg = float(telomere_gap_avg)
                except ValueError:
                    gap_avg = 0.0

                breakpoints.append({
                    'chrom': chrom,
                    'position': position,
                    'breakage_score': breakage_score,
                    'orientation': orientation,
                    'telomere_sense': telomere_sense,
                    'telomere_gap_avg': gap_avg,
                    'average_coverage': avg_cov,
                    'num_telomere_reads': num_telomere_reads,
                    'telomere_senses_dist': telomere_senses,
                    'telomere_gaps_dist': telomere_gaps
                })

    # Sort by position
    breakpoints.sort(key=lambda x: x['position'])

    print(f"  Found {len(breakpoints)} chromosome_breakage_site features")

    # Write TSV output with header
    with open(output_path, 'w') as out:
        # Write header
        header = [
            'chromosome',
            'position',
            'breakage_score',
            'num_telomere_reads',
            'average_coverage',
            'orientation',
            'telomere_sense',
            'telomere_gap_avg',
            'telomere_senses_dist',
            'telomere_gaps_dist'
        ]
        out.write('\t'.join(header) + '\n')

        # Write data
        for bp in breakpoints:
            out.write(
                f"{bp['chrom']}\t"
                f"{bp['position']}\t"
                f"{bp['breakage_score']:.4f}\t"
                f"{bp['num_telomere_reads']}\t"
                f"{bp['average_coverage']:.1f}\t"
                f"{bp['orientation']}\t"
                f"{bp['telomere_sense']}\t"
                f"{bp['telomere_gap_avg']:.2f}\t"
                f"{bp['telomere_senses_dist']}\t"
                f"{bp['telomere_gaps_dist']}\n"
            )

    print(f"  Wrote {len(breakpoints)} breakpoints to {output_path}")

    # Report key breakpoints near GRS boundaries
    grs_breakpoints = [bp for bp in breakpoints if GRS_START - 100 <= bp['position'] <= GRS_END + 100]
    if grs_breakpoints:
        print(f"  Breakpoints near GRS boundaries ({GRS_START}-{GRS_END}):")
        for bp in grs_breakpoints:
            print(
                f"    {bp['chrom']}:{bp['position']} "
                f"(breakage_score={bp['breakage_score']:.4f}, "
                f"telomere_reads={bp['num_telomere_reads']}, "
                f"coverage={bp['average_coverage']:.1f}, "
                f"orientation={bp['orientation']}, "
                f"sense={bp['telomere_sense']}, "
                f"gap_avg={bp['telomere_gap_avg']:.2f})"
            )


def main():
    parser = argparse.ArgumentParser(
        description="Generate visualization files for GRS coverage analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        '--coverage',
        required=True,
        help='Path to mosdepth per-base BED file (gzipped)'
    )

    parser.add_argument(
        '--telomeres',
        required=True,
        help='Path to MILTEL GFF3 file'
    )

    parser.add_argument(
        '--output-dir',
        default='.',
        help='Output directory for generated files (default: current directory)'
    )

    parser.add_argument(
        '--chromosome',
        default='SUPER_X',
        help='Chromosome name (default: SUPER_X)'
    )

    args = parser.parse_args()

    # Create output directory if needed
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Define output files
    windowed_output = output_dir / f"{args.chromosome}_GRS_windowed_coverage.bed"
    borders_output = output_dir / f"{args.chromosome}_GRS_borders.bed"
    telomeric_output = output_dir / f"{args.chromosome}_telomeric_breakpoints.tsv"

    print("="*70)
    print("GRS Visualization File Generator")
    print("="*70)

    # Generate windowed coverage
    generate_windowed_coverage(args.coverage, windowed_output, args.chromosome)

    # Generate border coverage
    generate_border_coverage(args.coverage, borders_output, args.chromosome)

    # Generate telomeric breakpoints
    generate_telomeric_breakpoints(args.telomeres, telomeric_output, args.chromosome)

    print("\n" + "="*70)
    print("Generation complete!")
    print("="*70)
    print(f"\nOutput files:")
    print(f"  1. Windowed coverage (BED):      {windowed_output}")
    print(f"  2. Border coverage (BED):        {borders_output}")
    print(f"  3. Telomeric breakpoints (TSV):  {telomeric_output}")
    print(f"\nThe TSV file includes:")
    print(f"  - breakage_score (telomere reads / total coverage)")
    print(f"  - num_telomere_reads (calculated from breakage_score × coverage)")
    print(f"  - orientation (left/right of junction)")
    print(f"  - telomere_sense (+/- strand)")
    print(f"  - telomere_gap_avg (bp from telomere to junction)")
    print(f"  - Full distribution data for senses and gaps")


if __name__ == "__main__":
    main()
