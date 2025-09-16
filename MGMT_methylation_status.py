#! /usr/bin/env python
# This script determines the methylation status of the MGMT gene in a given coverage file.

import gzip
import argparse

def load_coverage(cov_file, reference="hg38"):
    """Load coverage data from a bismark coverage file.
       Skip over anything not in the region of interest (MGMT gene).
    """
    coverage_DMR1 = {}
    coverage_DMR2 = {}
    assert reference in ["hg19", "hg38"], "Reference must be either 'hg19' or 'hg38'"
    if reference == "hg19":
        # DMR1: chr10:131264949-131265710
        # DMR2: chr10:131265496-131265627
        dmr1_start, dmr1_end = 131264949, 131265710
        dmr2_start, dmr2_end = 131265496, 131265627
    else:  # hg38
        # DMR1: chr10:129466685-129467446
        # DMR2: chr10:129467232-129467363
        dmr1_start, dmr1_end = 129466685, 129467446
        dmr2_start, dmr2_end = 129467232, 129467363
    if cov_file.endswith('.gz'):
        with gzip.open(cov_file, 'rt') as f:
            lines = f.readlines()
    else:
        with open(cov_file, 'r') as f:
            lines = f.readlines()
    
    for line in lines:
        if line.startswith('#') or not line.strip():
            continue
        parts = line.strip().split('\t')
        if len(parts) < 6:
            continue  # Skip malformed lines
        chrom, pos, _, percent_methylation, meth_count, unmeth_count = parts[:6]
        if (chrom == 'chr10' or chrom == 10) and (dmr1_start <= int(pos) <= dmr1_end):
            key = (chrom, int(pos))
            coverage_DMR1[key] = {
                'percent_methylation': float(percent_methylation),
                'meth_count': int(meth_count),
                'unmeth_count': int(unmeth_count),
                'total_coverage': int(meth_count) + int(unmeth_count)
            }
        elif (chrom == 'chr10' or chrom == 10) and (dmr2_start < int(pos) < dmr2_end):
            key = (chrom, int(pos))
            coverage_DMR2[key] = {
                'percent_methylation': float(percent_methylation),
                'meth_count': int(meth_count),
                'unmeth_count': int(unmeth_count),
                'total_coverage': int(meth_count) + int(unmeth_count)
            }
    return coverage_DMR1, coverage_DMR2

def determine_mgmt_status(coverage_DMR1, coverage_DMR2, coverage_threshold=3, 
                          meth_threshold_upper_DMR1=0.172, meth_threshold_lower_DMR1=0.076,
                          meth_threshold_upper_DMR2=0.256, meth_threshold_lower_DMR2=0.041):
    """Determine the methylation status of the MGMT gene"""
    used_positions = {}
    discarded_positions = {}
    total_meth_DMR1 = 0
    total_unmeth_DMR1 = 0
    total_meth_DMR2 = 0
    total_unmeth_DMR2 = 0

    # Get usable positions for DMR1
    for key, data in coverage_DMR1.items():
        if data['total_coverage'] >= coverage_threshold:
            used_positions[key] = data
            total_meth_DMR1 += data['meth_count']
            total_unmeth_DMR1 += data['unmeth_count']
        else:
            discarded_positions[key] = data

    # Calculate methylation percentage for DMR1
    # Avoid division by zero
    assert (total_meth_DMR1 + total_unmeth_DMR1) > 0, "No usable positions in DMR1"
    percent_methylation_DMR1 = total_meth_DMR1 / (total_meth_DMR1 + total_unmeth_DMR1)
    if percent_methylation_DMR1 >= meth_threshold_upper_DMR1:
        status_DMR1 = "methylated"
    elif percent_methylation_DMR1 <= meth_threshold_lower_DMR1:
        status_DMR1 = "unmethylated"
    else:
        status_DMR1 = "greyzone"
    # If the first DMR is inonclusive, check the second DMR
    if status_DMR1 == "greyzone":
        # Get usable positions for DMR2
        for key, data in coverage_DMR2.items():
            if data['total_coverage'] >= coverage_threshold:
                used_positions[key] = data
                total_meth_DMR2 += data['meth_count']
                total_unmeth_DMR2 += data['unmeth_count']
            else:
                discarded_positions[key] = data
        # Calculate methylation percentage for DMR2
        # Avoid division by zero
        assert (total_meth_DMR2 + total_unmeth_DMR2) > 0, "No usable positions in DMR2"
        percent_methylation_DMR2 = total_meth_DMR2 / (total_meth_DMR2 + total_unmeth_DMR2)
        if percent_methylation_DMR2 >= meth_threshold_upper_DMR2:
            status_DMR2 = "methylated"
        elif percent_methylation_DMR2 <= meth_threshold_lower_DMR2:
            status_DMR2 = "unmethylated"
        else:
            status_DMR2 = "greyzone"
        return status_DMR1, percent_methylation_DMR1, status_DMR2, percent_methylation_DMR2, used_positions, discarded_positions
    else:
        return status_DMR1, percent_methylation_DMR1, None, None, used_positions, discarded_positions


def write_report(report_file, status_DMR1, percent_methylation_DMR1, status_DMR2, percent_methylation_DMR2, used_positions, discarded_positions, input_file):
    """Write the methylation status report to a file"""
    with open(report_file, 'w') as f:
        f.write(f"MGMT Methylation Status Report\n")
        f.write(f"=============================\n\n")
        f.write(f"Input file used for analysis: {input_file}\n\n")

        f.write(f"DMR1 Status: {status_DMR1}\nDMR1 methylation percentage: {percent_methylation_DMR1:.3}\n\n")
        if status_DMR2:
            f.write(f"DMR2 Status: {status_DMR2}\nDMR2 methylation percentage: {percent_methylation_DMR2:.3}\n\n")
        f.write(f"Based on {len(used_positions)} used positions and {len(discarded_positions)} discarded positions due to low coverage.\n\n")
        f.write(f"\nUsed Positions:\n")
        for key, data in used_positions.items():
            f.write(f"{key}: {data}\n")
        f.write(f"\nDiscarded Positions (due to low coverage):\n")
        for key, data in discarded_positions.items():
            f.write(f"{key}: {data}\n")

def main():
    parser = argparse.ArgumentParser(description="Determine MGMT methylation status from bismark coverage file.")
    parser.add_argument("input_file", help="Input bismark coverage file (can be gzipped)")
    parser.add_argument("output_report", help="Output report file")
    parser.add_argument("-r", "--reference", choices=["hg19", "hg38"], default="hg38", help="Reference genome version (default: hg38)")
    parser.add_argument("-c", "--coverage_threshold", type=int, default=3, help="Minimum coverage to consider a position (default: 3)")
    parser.add_argument("--meth_threshold_upper_DMR1", type=float, default=0.172, help="Upper methylation threshold for DMR1 (default: 0.172)")
    parser.add_argument("--meth_threshold_lower_DMR1", type=float, default=0.076, help="Lower methylation threshold for DMR1 (default: 0.076)")
    parser.add_argument("--meth_threshold_upper_DMR2", type=float, default=0.256, help="Upper methylation threshold for DMR2 (default: 0.256)")
    parser.add_argument("--meth_threshold_lower_DMR2", type=float, default=0.041, help="Lower methylation threshold for DMR2 (default: 0.041)")
    args = parser.parse_args()

    coverage_DMR1, coverage_DMR2 = load_coverage(args.input_file, args.reference)
    status_DMR1, percent_methylation_DMR1, status_DMR2, percent_methylation_DMR2, used_positions, discarded_positions = determine_mgmt_status(coverage_DMR1, coverage_DMR2,
                                                                                          args.coverage_threshold, 
                                                                                          args.meth_threshold_upper_DMR1, args.meth_threshold_lower_DMR1,
                                                                                          args.meth_threshold_upper_DMR2, args.meth_threshold_lower_DMR2)
    write_report(args.output_report, status_DMR1, percent_methylation_DMR1, status_DMR2, percent_methylation_DMR2, used_positions, discarded_positions, args.input_file)
    print(f"MGMT methylation status report written to {args.output_report}")
    print(f"DMR1 Status: {status_DMR1}, Methylation Percentage: {percent_methylation_DMR1:.3}")
    if status_DMR2:
        print(f"DMR2 Status: {status_DMR2}, Methylation Percentage: {percent_methylation_DMR2:.3}")

if __name__ == "__main__":
    main()