#!/usr/bin/env python3
# Script to create bedMethyl files from bismark coverage files
import gzip
import argparse
import re
import pandas as pd

def parse_arguments():
    parser = argparse.ArgumentParser(description="Convert bismark coverage files to bedMethyl format.")
    parser.add_argument("input_file", help="Input bismark coverage file (can be gzipped)")
    parser.add_argument("output_file", help="Output bedMethyl file")
    parser.add_argument("--filter", type=int, default=0, help="Filter threshold for total coverage (default: 0)")
    parser.add_argument("--epic_annotation", type=str, help="Path to EPIC annotation file to add Probe IDs")
    return parser.parse_args()  

def read_bismark_coverage(file_path):
    # Load coverage file into a DataFrame
    df = pd.read_csv(file_path, sep='\t', header=None, dtype={0: str, 1: int, 2: int, 3: float, 4: int, 5: int}, compression='infer')
    # Expect columns: chrom, start, end, percent_methylation, meth_count, unmeth_count
    df.columns = ['chrom', 'start', 'end', 'percent_methylation', 'meth_count', 'unmeth_count']
    return df

def add_probe_ids(df, epic_annotation):
    # Load EPIC annotation
    annotation_df = pd.read_csv(epic_annotation, sep='\t')
    annotation_df = annotation_df[['CpG_chrm', 'CpG_beg', 'Probe_ID']]
    # Remove 'chr' prefix if present
    annotation_df['CpG_chrm'] = annotation_df['CpG_chrm'].apply(lambda x: re.sub(r'^chr', '', str(x)))
    # Make start pos integer and add 1 to match bismark coverage format
    annotation_df['CpG_beg'] = pd.to_numeric(annotation_df['CpG_beg'], errors='coerce').fillna(-1).astype(int) + 1
    # Remove characters after _ from probe ids
    annotation_df['Probe_ID'] = annotation_df['Probe_ID'].apply(lambda x: x.split('_')[0] if isinstance(x, str) else x)
    # Remove 'chr' prefix from df if present
    df['chrom'] = df['chrom'].apply(lambda x: re.sub(r'^chr', '', str(x)))
    df['start'] = pd.to_numeric(df['start'], errors='coerce').fillna(-1).astype(int)
    # Merge with annotation to add Probe_IDs
    merged_df = pd.merge(
        df, 
        annotation_df, 
        left_on=['chrom', 'start'], 
        right_on=['CpG_chrm', 'CpG_beg'], 
        how='left'
    )
    merged_df = merged_df.drop(columns=['CpG_chrm', 'CpG_beg'])
    return merged_df

def convert_to_bedmethyl(df, filter_threshold=0):
    # Calculate total coverage
    df['meth_count'] = pd.to_numeric(df['meth_count'], errors='coerce').fillna(0).astype(int)
    df['unmeth_count'] = pd.to_numeric(df['unmeth_count'], errors='coerce').fillna(0).astype(int)
    df['total_coverage'] = df['meth_count'] + df['unmeth_count']
    # Filter by coverage
    df = df[df['total_coverage'] >= filter_threshold]
    # Use Probe_ID if available, else "m"
    df['name'] = df['Probe_ID'].combine_first(pd.Series(["m"] * len(df)))
    # Format bedMethyl entry
    bedmethyl_df = pd.DataFrame({
        'chrom': df['chrom'],
        'chromStart': df['start'],
        'chromEnd': df['end'],
        'name': df['name'],
        'score': df['total_coverage'],
        'strand': ['.'] * len(df),
        'thickStart': df['start'],
        'thickEnd': df['end'],
        'itemRgb': ['255,0,0'] * len(df),
        'coverage': df['total_coverage'],
        'MAF': df['percent_methylation']
    })
    # Remove duplicate rows based on chrom, chromStart, chromEnd, keeping the first occurrence
    bedmethyl_df = bedmethyl_df.drop_duplicates(subset=['chrom', 'chromStart', 'chromEnd'])
    # Drop rows where 'name' is still 'm' (optional, if you want only annotated probes)
    # bedmethyl_df = bedmethyl_df[~bedmethyl_df['name'].str.contains('m', na=False)]
    return bedmethyl_df

def write_bedmethyl_file(bedmethyl_df, output_file):
    bedmethyl_df.to_csv(output_file, sep='\t', header=False, index=False)

def main():
    args = parse_arguments()
    df = read_bismark_coverage(args.input_file)
    if args.epic_annotation:
        df = add_probe_ids(df, args.epic_annotation)
    bedmethyl_df = convert_to_bedmethyl(df, filter_threshold=args.filter)
    write_bedmethyl_file(bedmethyl_df, args.output_file)

if __name__ == "__main__":
    main()


