#!/usr/bin/env python3
import os
from glob import glob
import pandas as pd
from pathlib import Path
import argparse
from tqdm import tqdm

# All taxonomic levels 
taxonomic_levels = {
    'R': 'Realm', 
    'R1': 'Subrealm', 'R2': 'Sub-subrealm',
    'D':'Domain', 'K': 'Kingdom', 
    'P': 'Phylum', 'C': 'Class', 
    'O': 'Order', 'F': 'Family', 
    'G': 'Genus', 'S': 'Species'
}
VALID_LEVELS = list(taxonomic_levels.keys())

if __name__ == "__main__":


    #####################################
    ####            ARGS             ####
    #####################################
    
    parser = argparse.ArgumentParser(description='Create phyloseq-format tables from Kraken2 reports.')
    parser.add_argument(
        '-i', 
        '--inputFolder', 
        help='Folder with Kraken2 reports.', 
        required=True, 
        type=Path)
    parser.add_argument(
        '-o', 
        '--otu', 
        help='Path to save OTU table. Default: "./otu_table.csv"', 
        default='otu_table.csv', 
        type=Path)
    parser.add_argument(
        '-t', 
        '--taxonomy', 
        help='Path to save taxonomy table. Default: "./taxonomy_table.csv"', 
        default='taxonomy_table.csv', 
        type=Path)
    parser.add_argument(
        '-f', 
        '--format', 
        help='Format of created Kraken2 report. Possible options: "default", "minimizer" (--report-minimizer-data option), "mpa" (--use-mpa-style option).', 
        choices=['default', 'minimizer', 'mpa'], 
        default='default', 
        type=str)
    parser.add_argument(
        '-e', 
        '--extension', 
        help='Extension of Kraken2 reports. Default = "txt".', 
        default='txt', 
        type=str)
    args = parser.parse_args()

    folder_with_reports = args.inputFolder
    reports = sorted(list(folder_with_reports.glob(f'*.{args.extension}')))

    if not reports:
        print(f'No Kraken2 reports found in {folder_with_reports}')
        exit()

    # Global dictionaries to collect data across all samples
    all_terminal_otus = set()                   # For lookups
    sample_counts = dict()                      # Dictionary: {sample_name: {otu_id: count}}
    taxonomy_cache = dict()                     # Dictionary: {otu_id: {Domain: x, Kingdom: y, ... Species: z}}

    # Define column indices based on format
    if args.format == 'default':
        otuID_col, otuCount_col, name_col, taxLevel_col = 4, 1, 5, 3
    elif args.format == 'minimizer':
        otuID_col, otuCount_col, name_col, taxLevel_col = 6, 1, 7, 5

    
    
    ################################################
    ####          PROCESS ALL REPORTS           ####
    ################################################

    for report_path in tqdm(reports, desc="Processing Kraken2 reports"):
        sample_name = report_path.stem
        sample_counts[sample_name] = {}
        
        if args.format in ['default', 'minimizer']:
            # 1. Read Kraken2 report with Pandas
            df = pd.read_csv(report_path, sep='\t', header=None)
            
            # Filter only valid taxonomic levels
            df = df[df.iloc[:, taxLevel_col].isin(VALID_LEVELS)].copy()
            
            # 2. Search for terminal node of identified taxonomy detection
            # Level is terminal if index of current taxonomic level is greater or equal to index of next taxonomic level in VALID_LEVELS list
            # Create dictionary with indices of taxonomic levels from VALID_LEVELS list
            tax_level_indices = {level: idx for idx, level in enumerate(VALID_LEVELS)}

            # Create a column with True/False statements indicating whether taxonomic level is terminal (True) or not (False)
            df['is_terminal'] = df.iloc[:,taxLevel_col].map(tax_level_indices) >= df.iloc[:,taxLevel_col].map(tax_level_indices).shift(-1)
            # Last row is always terminal
            df.iloc[-1, df.columns.get_loc('is_terminal')] = True
            
            # Create new dataframe with terminal nodes only
            terminal_df = df[df.loc[:, 'is_terminal'] == True]
            
            # 3. Extract OTU IDs and counts
            terminal_otus = terminal_df.iloc[:, otuID_col].astype(str).values
            terminal_counts = terminal_df.iloc[:, otuCount_col].values
            
            # Add new terminal OTU IDs if any are present
            all_terminal_otus.update(terminal_otus)
            # Add OTU IDs and corresponding counts to according column with current sample name in sample_counts dataframe
            sample_counts[sample_name] = dict(zip(terminal_otus, terminal_counts))
            
            # 4. Extract taxonomy of terminal OTU IDs
            # Dictionary with taxonomic lineage for current terminal OTU
            current_lineage = {level: None for level in taxonomic_levels.values()}
            
            # Parse Kraken2 report file
            for row in df.itertuples(index=False):
                # One-letter code for taxonomic level
                tax_level = row[taxLevel_col]
                # ORU ID
                otu_id = str(row[otuID_col])
                # Taxon's name
                taxonomic_name = row[name_col].strip()
                # Full name of taxonomic level
                level_name = taxonomic_levels[tax_level]
                # Add current taxon's name to current taxonomic level in dictionary of current taxonomic lineage
                current_lineage[level_name] = taxonomic_name
                
                # Clear lower levels
                level_idx = VALID_LEVELS.index(tax_level)
                for key in VALID_LEVELS[level_idx+1:]:
                    current_lineage[taxonomic_levels[key]] = None
                    
                # Cache the current state for this OTU ID
                taxonomy_cache[otu_id] = current_lineage.copy()


        elif args.format == 'mpa':
            # 1. 1. Read Kraken2 report in MPA format with Pandas 
            df = pd.read_csv(report_path, sep='\t', header=None, names=['taxonomy', 'count'])
            
            # 2. Search for terminal node of identified taxonomy detection
            # Depth in MPA-format report is determined by the number of '|' separators
            df['depth'] = df['taxonomy'].str.count(r'\|')
            next_depth = df['depth'].shift(-1)
            is_terminal = (next_depth <= df['depth'])
            # The last row is always terminal
            is_terminal.iloc[-1] = True
            
            # Create new dataframe with terminal nodes only
            terminal_df = df[is_terminal]
            
            # 3. Extract OTU IDs and counts
            terminal_otus = terminal_df['taxonomy'].str.split('|').str[-1].str.split('__').str[-1].values
            terminal_counts = terminal_df['count'].values
            
            # Add new OTU IDs if any are present
            all_terminal_otus.update(terminal_otus)
            # Add OTU IDs and corresponding counts to according column with current sample name in sample_counts dataframe
            sample_counts[sample_name] = dict(zip(terminal_otus, terminal_counts))
            
            # 4. Extract taxonomy of terminal OTU IDs
            # Parse Kraken2 report file
            for row in df.itertuples(index=False):
                tax_line = row.taxonomy
                # Extract the lowest level name for the OTU ID key from defined taxonomic levels
                otu_id = tax_line.split('|')[-1].split('__')[-1].strip()
                
                # Dictionary with taxonomic lineage for current terminal OTU
                lineage = {level: None for level in taxonomic_levels.values()}
                # Split current taxonomy0lineage line by '|'
                # Parts are taxonomic names
                for part in tax_line.split('|'):
                    level_prefix = part.split('__')[0].upper()
                    if level_prefix in taxonomic_levels:
                        # Add current taxon's name to current taxonomic level in dictionary of current taxonomic lineage
                        lineage[taxonomic_levels[level_prefix]] = part.split('__')[-1]

                # Cache the current state for this OTU ID      
                taxonomy_cache[otu_id] = lineage



    ################################################
    ####        BUILD FINAL DATAFRAMES          ####
    ################################################

    print("Building final tables...")
    
    # 1. Build OTU Table
    otu_data = []
    for otu in tqdm(all_terminal_otus, desc="Building OTU table"):
        row_data = {'otu': otu}
        for sample in sample_counts:
            # Map OTU counts 
            row_data[sample] = sample_counts[sample].get(otu, 0)
        otu_data.append(row_data)
        
    otu_counts_df = pd.DataFrame(otu_data)
    # Sort sample columns
    cols = ['otu'] + sorted([col for col in otu_counts_df.columns if col != 'otu'])
    otu_counts_df = otu_counts_df[cols]

    # 2. Build Taxonomy Table
    tax_data = []
    for otu in tqdm(all_terminal_otus, desc="Building taxonomy table"):
        # Map taxonomies
        if otu in taxonomy_cache:
            row_data = {'otu': otu}
            row_data.update(taxonomy_cache[otu])
            tax_data.append(row_data)
            
    taxonomy_table_df = pd.DataFrame(tax_data)
    # Order taxonomy columns
    tax_cols = ['otu'] + list(taxonomic_levels.values())
    taxonomy_table_df = taxonomy_table_df[tax_cols]

    
    
    #############################################
    ####             SAVE RESULTS            ####
    #############################################
    
    otu_counts_df.to_csv(args.otu, index=False)
    print(f"OTU table saved as {args.otu}")
    
    taxonomy_table_df.to_csv(args.taxonomy, index=False)
    print(f"Taxonomy table saved as {args.taxonomy}")