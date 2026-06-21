#!/usr/bin/env python3
"""
FASTA to Karyotype TSV Converter (with Header Parsing)

Acknowledgments:
This script was developed with the assistance of an AI language model 
to facilitate the automated extraction of sequence lengths from FASTA 
files, including dynamic parsing of organism and chromosome metadata 
from FASTA headers, formatted for standard karyotype visualization.
"""

import argparse
import sys
import os
import re

def parse_header_for_label(header_line):
    """
    Parses a FASTA header to extract a short label.
    Example Input:  >OX637739.1 Avena sativa genome assembly, chromosome: 1A
    Example Output: A.sativa_1A
    """
    clean_header = header_line[1:].strip()
    tokens = clean_header.split()
    
    if len(tokens) < 2:
        return None
        
    # 1. Guess the organism name (assuming it's the 2nd and 3rd words in the header)
    org_short = ""
    if len(tokens) >= 3:
        genus = tokens[1]
        species = tokens[2]
        # Make sure they look like words to avoid parsing weird numbers/symbols
        if genus.isalpha() and species.rstrip(',;').isalpha():
            # Takes first letter of genus + "." + species (removing any trailing commas)
            org_short = f"{genus[0]}.{species.rstrip(',;')}"

    # 2. Look for the "chromosome:" tag using Regular Expressions
    # This matches "chromosome: 1A", "chromosome 1A", "Chromosome: 1A", etc.
    chrom_match = re.search(r'chromosome[:\s]+([a-zA-Z0-9_]+)', clean_header, re.IGNORECASE)
    
    # 3. Combine them if found
    if chrom_match and org_short:
        chrom_id = chrom_match.group(1)
        return f"{org_short}_{chrom_id}"
    elif chrom_match:
        # If we found a chromosome but no organism name
        chrom_id = chrom_match.group(1)
        return f"Chr_{chrom_id}"
        
    # If no chromosome is specified in this header, return None to trigger the fallback
    return None

def get_fasta_data(fasta_path):
    """
    Reads a FASTA file, calculates length, and attempts to parse the header label.
    Returns a dictionary: {seq_id: {'length': int, 'label': str or None}}
    """
    data = {}
    if not os.path.exists(fasta_path):
        sys.exit(f"Error: The file {fasta_path} does not exist.")

    with open(fasta_path, 'r') as f:
        curr_id = None
        curr_len = 0
        curr_label = None
        
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # Save previous sequence data
                if curr_id is not None:
                    data[curr_id] = {'length': curr_len, 'label': curr_label}
                
                # Start new sequence
                curr_id = line[1:].split()[0] 
                curr_len = 0
                curr_label = parse_header_for_label(line)
            else:
                curr_len += len(line)
                
        # Add the final sequence
        if curr_id is not None:
            data[curr_id] = {'length': curr_len, 'label': curr_label}
            
    return data

def write_karyotype_tsv(data, output_path, color, prefix, suffix):
    """
    Writes the parsed data into the standard 7-column karyotype format.
    """
    unplaced_counter = 1
    
    with open(output_path, 'w') as out_f:
        for seq_id, info in data.items():
            seq_len = info['length']
            parsed_label = info['label']
            
            # If the script successfully parsed "A.sativa_1A" from the header, use it.
            # If it didn't find "chromosome:" in the header, fall back to the prefix/suffix logic.
            if parsed_label is not None:
                final_label = parsed_label
            else:
                final_label = f"{prefix}{unplaced_counter}{suffix}"
                unplaced_counter += 1
            
            # Construct columns
            col1_chr    = "chr"                           
            col2_dash   = "-"                             
            col3_id     = seq_id                          
            col4_label  = final_label          
            col5_start  = "0"                             
            col6_end    = str(seq_len)                    
            col7_color  = color                           
            
            row = f"{col1_chr}\t{col2_dash}\t{col3_id}\t{col4_label}\t{col5_start}\t{col6_end}\t{col7_color}\n"
            out_f.write(row)
            
    print(f"Successfully wrote {len(data)} sequences to {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Convert a FASTA file to a Karyotype TSV file.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Path to the output TSV file")
    parser.add_argument("-c", "--color", default="green", help="Color for the chromosomes (default: green)")
    # Prefix/Suffix are now used as fallbacks for contigs/scaffolds that don't have "chromosome:" in the header
    parser.add_argument("-p", "--prefix", default="scaf.", help="Fallback prefix for unplaced scaffolds (default: scaf.)")
    parser.add_argument("-s", "--suffix", default="", help="Fallback suffix for unplaced scaffolds")

    # FIXED: This uses parse_args() now!
    args = parser.parse_args()

    # 1. Parse Data
    seq_data = get_fasta_data(args.input)
    
    if not seq_data:
        sys.exit("Error: No sequences found in the FASTA file.")

    # 2. Write TSV
    write_karyotype_tsv(
        data=seq_data, 
        output_path=args.output, 
        color=args.color, 
        prefix=args.prefix, 
        suffix=args.suffix
    )

if __name__ == "__main__":
    main()