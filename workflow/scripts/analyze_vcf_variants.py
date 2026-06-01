#!/usr/bin/env python3
"""
VCF Structural Variant Analysis Script
Analyzes NanoVar VCF files to count and categorize structural variants
"""

import os
import sys
import glob
import argparse
from collections import defaultdict, Counter
import pandas as pd
import re

def parse_vcf_file(vcf_path):
    """
    Parse a single VCF file and extract structural variant information
    
    Args:
        vcf_path (str): Path to the VCF file
    
    Returns:
        dict: Dictionary containing variant counts and details
    """
    variant_counts = Counter()
    variant_details = []
    total_variants = 0
    
    print(f"Processing: {os.path.basename(vcf_path)}")
    
    try:
        with open(vcf_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                # Skip header lines
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 8:
                    continue
                
                chrom = fields[0]
                pos = int(fields[1])
                var_id = fields[2]
                ref = fields[3]
                alt = fields[4]
                qual = fields[5]
                filter_field = fields[6]
                info = fields[7]
                
                # Parse INFO field
                info_dict = {}
                for item in info.split(';'):
                    if '=' in item:
                        key, value = item.split('=', 1)
                        info_dict[key] = value
                
                # Extract variant type
                svtype = info_dict.get('SVTYPE', 'UNKNOWN')
                end_pos = int(info_dict.get('END', pos))
                svlen = int(info_dict.get('SVLEN', 0)) if info_dict.get('SVLEN', '0').lstrip('-').isdigit() else 0
                support_reads = int(info_dict.get('SR', 0)) if info_dict.get('SR', '0').isdigit() else 0
                nn_confidence = float(info_dict.get('NN', 0.0)) if info_dict.get('NN', '0.0').replace('.', '').isdigit() else 0.0
                
                # Count variants
                variant_counts[svtype] += 1
                total_variants += 1
                
                # Store detailed information
                variant_details.append({
                    'chromosome': chrom,
                    'position': pos,
                    'variant_id': var_id,
                    'svtype': svtype,
                    'length': abs(svlen),
                    'end_position': end_pos,
                    'quality': qual,
                    'filter': filter_field,
                    'support_reads': support_reads,
                    'nn_confidence': nn_confidence,
                    'line_number': line_num
                })
    
    except FileNotFoundError:
        print(f"Error: File {vcf_path} not found")
        return None
    except Exception as e:
        print(f"Error processing {vcf_path}: {e}")
        return None
    
    return {
        'file_path': vcf_path,
        'variant_counts': variant_counts,
        'variant_details': variant_details,
        'total_variants': total_variants
    }

def analyze_directory(directory_path, pattern="*.vcf"):
    """
    Analyze all VCF files in a directory - each file separately
    
    Args:
        directory_path (str): Path to directory containing VCF files
        pattern (str): File pattern to match (default: "*.vcf")
    
    Returns:
        dict: Analysis results for all files
    """
    # Find all VCF files
    search_pattern = os.path.join(directory_path, pattern)
    vcf_files = glob.glob(search_pattern)
    
    if not vcf_files:
        print(f"No VCF files found matching pattern: {search_pattern}")
        return None
    
    print(f"Found {len(vcf_files)} VCF files to process")
    print("=" * 80)
    
    all_results = {}
    
    for vcf_file in sorted(vcf_files):
        result = parse_vcf_file(vcf_file)
        if result:
            filename = os.path.basename(vcf_file)
            all_results[filename] = result
    
    return {
        'individual_results': all_results,
        'total_files': len(all_results)
    }

def generate_individual_report(filename, result):
    """
    Generate a report for a single VCF file
    
    Args:
        filename (str): Name of the VCF file
        result (dict): Analysis result for the file
    
    Returns:
        str: Report text for the file
    """
    variant_counts = result['variant_counts']
    variant_details = result['variant_details']
    total_variants = result['total_variants']
    
    report = []
    report.append("=" * 80)
    report.append(f"STRUCTURAL VARIANT ANALYSIS: {filename}")
    report.append("=" * 80)
    report.append(f"Total variants found: {total_variants}")
    report.append("")
    
    if total_variants == 0:
        report.append("No structural variants found in this file.")
        return "\n".join(report)
    
    # Variant type summary
    report.append("VARIANT TYPE SUMMARY:")
    report.append("-" * 30)
    for svtype in sorted(variant_counts.keys()):
        count = variant_counts[svtype]
        percentage = (count / total_variants) * 100
        report.append(f"{svtype:>10}: {count:>6} ({percentage:>5.1f}%)")
    report.append("")
    
    # Length and quality statistics
    length_stats = defaultdict(list)
    quality_stats = defaultdict(list)
    
    for variant in variant_details:
        svtype = variant['svtype']
        if variant['length'] > 0:
            length_stats[svtype].append(variant['length'])
        if variant['quality'] != '.' and variant['quality'].replace('.', '').isdigit():
            try:
                quality_stats[svtype].append(float(variant['quality']))
            except ValueError:
                pass
    
    report.append("DETAILED STATISTICS:")
    report.append("-" * 25)
    
    for svtype in sorted(variant_counts.keys()):
        report.append(f"\n{svtype} variants:")
        report.append(f"  Count: {variant_counts[svtype]}")
        
        if length_stats[svtype]:
            lengths = length_stats[svtype]
            report.append(f"  Length - Min: {min(lengths)}, Max: {max(lengths)}, "
                         f"Median: {sorted(lengths)[len(lengths)//2]}")
        
        if quality_stats[svtype]:
            qualities = quality_stats[svtype]
            report.append(f"  Quality - Min: {min(qualities):.1f}, Max: {max(qualities):.1f}, "
                         f"Mean: {sum(qualities)/len(qualities):.1f}")
    
    # INV and DUP size analysis
    if 'INV' in variant_counts or 'DUP' in variant_counts:
        report.append("\n" + "="*60)
        report.append("INV AND DUP SIZE DISTRIBUTION ANALYSIS")
        report.append("="*60)
        
        # Analyze INV variants
        if 'INV' in variant_counts and length_stats['INV']:
            inv_lengths = length_stats['INV']
            inv_short = [l for l in inv_lengths if l < 1000]
            inv_long = [l for l in inv_lengths if l >= 1000]
            
            report.append(f"\nINVERSION (INV) SIZE ANALYSIS:")
            report.append(f"  Total INV variants: {len(inv_lengths)}")
            report.append(f"  INV < 1000bp: {len(inv_short)} ({len(inv_short)/len(inv_lengths)*100:.1f}%)")
            report.append(f"  INV >= 1000bp: {len(inv_long)} ({len(inv_long)/len(inv_lengths)*100:.1f}%)")
            
            if inv_short:
                report.append(f"  Short INV stats - Min: {min(inv_short)}bp, Max: {max(inv_short)}bp, "
                             f"Mean: {sum(inv_short)/len(inv_short):.0f}bp")
            if inv_long:
                report.append(f"  Long INV stats - Min: {min(inv_long)}bp, Max: {max(inv_long)}bp, "
                             f"Mean: {sum(inv_long)/len(inv_long):.0f}bp")
        elif 'INV' in variant_counts:
            report.append(f"\nINVERSION (INV) SIZE ANALYSIS:")
            report.append(f"  Total INV variants: {variant_counts['INV']}")
            report.append(f"  (No length data available for size distribution)")
        
        # Analyze DUP variants
        if 'DUP' in variant_counts and length_stats['DUP']:
            dup_lengths = length_stats['DUP']
            dup_short = [l for l in dup_lengths if l < 1000]
            dup_long = [l for l in dup_lengths if l >= 1000]
            
            report.append(f"\nDUPLICATION (DUP) SIZE ANALYSIS:")
            report.append(f"  Total DUP variants: {len(dup_lengths)}")
            report.append(f"  DUP < 1000bp: {len(dup_short)} ({len(dup_short)/len(dup_lengths)*100:.1f}%)")
            report.append(f"  DUP >= 1000bp: {len(dup_long)} ({len(dup_long)/len(dup_lengths)*100:.1f}%)")
            
            if dup_short:
                report.append(f"  Short DUP stats - Min: {min(dup_short)}bp, Max: {max(dup_short)}bp, "
                             f"Mean: {sum(dup_short)/len(dup_short):.0f}bp")
            if dup_long:
                report.append(f"  Long DUP stats - Min: {min(dup_long)}bp, Max: {max(dup_long)}bp, "
                             f"Mean: {sum(dup_long)/len(dup_long):.0f}bp")
        elif 'DUP' in variant_counts:
            report.append(f"\nDUPLICATION (DUP) SIZE ANALYSIS:")
            report.append(f"  Total DUP variants: {variant_counts['DUP']}")
            report.append(f"  (No length data available for size distribution)")
        
        # Combined summary
        total_inv = variant_counts.get('INV', 0)
        total_dup = variant_counts.get('DUP', 0)
        
        if total_inv > 0 or total_dup > 0:
            report.append(f"\nCOMBINED INV + DUP SUMMARY:")
            report.append(f"  Total structural variants (INV + DUP): {total_inv + total_dup}")
            report.append(f"  INV variants: {total_inv}")
            report.append(f"  DUP variants: {total_dup}")
            
            if length_stats['INV'] or length_stats['DUP']:
                all_inv_dup_lengths = length_stats['INV'] + length_stats['DUP']
                if all_inv_dup_lengths:
                    short_total = len([l for l in all_inv_dup_lengths if l < 1000])
                    long_total = len([l for l in all_inv_dup_lengths if l >= 1000])
                    report.append(f"  Combined < 1000bp: {short_total} ({short_total/len(all_inv_dup_lengths)*100:.1f}%)")
                    report.append(f"  Combined >= 1000bp: {long_total} ({long_total/len(all_inv_dup_lengths)*100:.1f}%)")
    
    return "\n".join(report)

def generate_summary_table(analysis_results):
    """
    Generate a summary table with all samples and their variant counts
    
    Args:
        analysis_results (dict): Results from analyze_directory
        
    Returns:
        str: Summary table text
    """
    individual_results = analysis_results['individual_results']
    
    # Collect all variant types
    all_svtypes = set()
    for result in individual_results.values():
        all_svtypes.update(result['variant_counts'].keys())
    all_svtypes = sorted(all_svtypes)
    
    # Prepare summary data
    summary_data = []
    
    for filename, result in sorted(individual_results.items()):
        variant_counts = result['variant_counts']
        variant_details = result['variant_details']
        
        # Calculate length statistics for INV, DUP, and INS
        inv_short = inv_long = dup_short = dup_long = ins_short = ins_long = 0
        
        for variant in variant_details:
            svtype = variant['svtype']
            length = variant['length']
            
            if svtype == 'INV' and length > 0:
                if length < 1000:
                    inv_short += 1
                else:
                    inv_long += 1
            elif svtype == 'DUP' and length > 0:
                if length < 1000:
                    dup_short += 1
                else:
                    dup_long += 1
            elif svtype == 'INS' and length > 0:
                if length < 1000:
                    ins_short += 1
                else:
                    ins_long += 1
        
        # Create row data
        row_data = {
            'filename': filename,
            'total': result['total_variants'],
            'inv_short': inv_short,
            'inv_long': inv_long,
            'dup_short': dup_short,
            'dup_long': dup_long,
            'ins_short': ins_short,
            'ins_long': ins_long
        }
        
        # Add counts for all variant types
        for svtype in all_svtypes:
            row_data[svtype.lower()] = variant_counts.get(svtype, 0)
        
        summary_data.append(row_data)
    
    # Generate table
    table = []
    table.append("=" * 120)
    table.append("SUMMARY TABLE - ALL SAMPLES")
    table.append("=" * 120)
    
    # Header
    header = f"{'Sample':<25}"
    for svtype in all_svtypes:
        header += f"{svtype:>8}"
    header += f"{'INV<1kb':>8}{'INV>=1kb':>9}{'DUP<1kb':>8}{'DUP>=1kb':>9}{'INS<1kb':>8}{'INS>=1kb':>9}{'Total':>8}"
    table.append(header)
    table.append("-" * len(header))
    
    # Data rows
    totals = {'total': 0, 'inv_short': 0, 'inv_long': 0, 'dup_short': 0, 'dup_long': 0, 'ins_short': 0, 'ins_long': 0}
    for svtype in all_svtypes:
        totals[svtype.lower()] = 0
    
    for row_data in summary_data:
        # Truncate filename if too long
        sample_name = row_data['filename']
        if len(sample_name) > 24:
            sample_name = sample_name[:21] + "..."
        
        row = f"{sample_name:<25}"
        
        # Add variant type counts
        for svtype in all_svtypes:
            count = row_data[svtype.lower()]
            row += f"{count:>8}"
            totals[svtype.lower()] += count
        
        # Add INV, DUP, and INS size categories
        row += f"{row_data['inv_short']:>8}{row_data['inv_long']:>9}"
        row += f"{row_data['dup_short']:>8}{row_data['dup_long']:>9}"
        row += f"{row_data['ins_short']:>8}{row_data['ins_long']:>9}"
        row += f"{row_data['total']:>8}"
        
        table.append(row)
        
        # Update totals
        totals['total'] += row_data['total']
        totals['inv_short'] += row_data['inv_short']
        totals['inv_long'] += row_data['inv_long']
        totals['dup_short'] += row_data['dup_short']
        totals['dup_long'] += row_data['dup_long']
        totals['ins_short'] += row_data['ins_short']
        totals['ins_long'] += row_data['ins_long']
    
    # Totals row
    table.append("-" * len(header))
    totals_row = f"{'TOTALS':<25}"
    for svtype in all_svtypes:
        totals_row += f"{totals[svtype.lower()]:>8}"
    totals_row += f"{totals['inv_short']:>8}{totals['inv_long']:>9}"
    totals_row += f"{totals['dup_short']:>8}{totals['dup_long']:>9}"
    totals_row += f"{totals['ins_short']:>8}{totals['ins_long']:>9}"
    totals_row += f"{totals['total']:>8}"
    table.append(totals_row)
    
    table.append("")
    table.append("Legend:")
    table.append(f"{'INV<1kb':<10}: Inversions < 1000bp")
    table.append(f"{'INV>=1kb':<10}: Inversions >= 1000bp")
    table.append(f"{'DUP<1kb':<10}: Duplications < 1000bp")
    table.append(f"{'DUP>=1kb':<10}: Duplications >= 1000bp")
    table.append(f"{'INS<1kb':<10}: Insertions < 1000bp")
    table.append(f"{'INS>=1kb':<10}: Insertions >= 1000bp")
    
    return "\n".join(table)

def generate_report(analysis_results, output_file=None):
    """
    Generate comprehensive reports for each VCF file separately
    
    Args:
        analysis_results (dict): Results from analyze_directory
        output_file (str): Optional output file path for the report
    """
    if not analysis_results:
        print("No analysis results to report")
        return
    
    individual_results = analysis_results['individual_results']
    total_files = analysis_results['total_files']
    
    all_reports = []
    
    print(f"Generating individual reports for {total_files} files...")
    print("\n")
    
    # Generate report for each file
    for filename, result in sorted(individual_results.items()):
        report_text = generate_individual_report(filename, result)
        all_reports.append(report_text)
        print(report_text)
        print("\n" + "="*80 + "\n")
    
    # Generate and display summary table
    summary_table = generate_summary_table(analysis_results)
    print(summary_table)
    all_reports.append(summary_table)
    
    # Save to file if requested
    if output_file:
        try:
            full_report = "\n\n".join(all_reports)
            with open(output_file, 'w') as f:
                f.write(full_report)
            print(f"\nAll individual reports and summary table saved to: {output_file}")
        except Exception as e:
            print(f"Error saving report to {output_file}: {e}")

def create_summary_csv(analysis_results, output_file):
    """
    Create a CSV file with detailed variant information
    
    Args:
        analysis_results (dict): Results from analyze_directory
        output_file (str): Output CSV file path
    """
    if not analysis_results:
        print("No analysis results to export")
        return
    
    all_variants = []
    
    for filename, result in analysis_results['individual_results'].items():
        for variant in result['variant_details']:
            variant_copy = variant.copy()
            variant_copy['source_file'] = filename
            all_variants.append(variant_copy)
    
    if all_variants:
        df = pd.DataFrame(all_variants)
        df.to_csv(output_file, index=False)
        print(f"Detailed variant data saved to: {output_file}")
    else:
        print("No variant data to export")

def main():
    parser = argparse.ArgumentParser(
        description="Analyze structural variants in NanoVar VCF files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze all VCF files in current directory
  python analyze_vcf_variants.py .
  
  # Analyze specific pattern
  python analyze_vcf_variants.py /path/to/vcf/files -p "*_combined.nanovar.pass.vcf"
  
  # Generate report and CSV output
  python analyze_vcf_variants.py /path/to/vcf/files -r report.txt -c variants.csv
        """
    )
    
    parser.add_argument('directory', help='Directory containing VCF files')
    parser.add_argument('-p', '--pattern', default='*.vcf', 
                        help='File pattern to match (default: *.vcf)')
    parser.add_argument('-r', '--report', 
                        help='Output file for detailed report')
    parser.add_argument('-c', '--csv', 
                        help='Output CSV file for detailed variant data')
    
    args = parser.parse_args()
    
    # Check if directory exists
    if not os.path.isdir(args.directory):
        print(f"Error: Directory '{args.directory}' not found")
        sys.exit(1)
    
    # Analyze the directory
    print(f"Analyzing VCF files in: {args.directory}")
    print(f"Pattern: {args.pattern}")
    print("=" * 50)
    
    results = analyze_directory(args.directory, args.pattern)
    
    if results:
        # Generate individual reports for each file
        generate_report(results, args.report)
        
        # Create CSV if requested (still contains all files' data)
        if args.csv:
            create_summary_csv(results, args.csv)
    else:
        print("No results to display")
        sys.exit(1)

if __name__ == "__main__":
    main()