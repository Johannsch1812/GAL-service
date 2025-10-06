import pandas as pd
import json
import os
import sys
import glob
from pathlib import Path
from typing import List, Optional

def read_database(database_file: str) -> Optional[pd.DataFrame]:
    """
    Read the ALMAGAL database.xlsx with 1017 targets and MOUS mapping
    -----------
    database_file : str
        Path to database.xlsx file
        
    Returns:
    --------
    pandas.DataFrame
        Database with target-MOUS mapping
    """
    try:
        import warnings # warnings from openpyxl is suppressed

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning, module="openpyxl")
            df = pd.read_excel(database_file, engine='openpyxl')

        print(f"Successfully read database with {len(df)} targets")
        print(f"Database columns: {df.columns.tolist()}")
        return df
    except Exception as e:
        print(f"Error reading database file: {e}")
        return None

def find_fits_files(fits_directory: str) -> List:
    """
    Find all FITS files in the specified directory
    -----------
    fits_directory : str
        Directory containing FITS files
        
    Returns:
    --------
    list
        List of FITS file paths
    """
    # Look for common FITS file patterns
    patterns = [
        "*.pbcor.fits", 
        # "*cont*.pbcor.fits", 
        # "*cube*.pbcor.fits", 
        # "*spw*.pbcor.fits",
    ]
    
    fits_files = []
    for pattern in patterns:
        search_path = os.path.join(fits_directory, "**", pattern)
        found_files = glob.glob(search_path, recursive=True)
        fits_files.extend(found_files)
    
    # Remove duplicates
    fits_files = list(set(fits_files))
    
    print(f"Found {len(fits_files)} FITS files")
    return fits_files

def check_fits_contains_terms(fits_filename: str, source_name: str, mous_id: str) -> bool:
    """
    Check if a FITS filename contains both source name and MOUS ID
    -----------
    fits_filename : str
        FITS filename (basename)
    source_name : str
        Source name to search for
    mous_id : str
        MOUS ID to search for
        
    Returns:
    --------
    bool
        True if both terms are found in filename
    """
    if pd.isna(source_name) or pd.isna(mous_id):
        return False
    
    source = str(source_name).strip()
    # if there is '+' in the Source name, '+' is replaced with 'p' in fits name
    source_convent = source.replace('+', 'p') 
    mous = str(mous_id).strip()
    
    # Check if both Source value and MOUS are in the filename
    source_in_filename = (source in fits_filename) or (source_convent in fits_filename)
    mous_in_filename = mous in fits_filename
    
    return source_in_filename and mous_in_filename

def find_matching_fits(source_name: str, mous_id: str, fits_files: List) -> List:
    """
    Find FITS files that contain both source name and MOUS ID
    -----------
    source_name : str
        Source name
    mous_id : str
        MOUS ID
    fits_files : list
        List of all FITS file paths
        
    Returns:
    --------
    list
        List of matching FITS file paths
    """
    matching_files = []
    
    for fits_file in fits_files:
        filename = os.path.basename(fits_file)
        if check_fits_contains_terms(filename, source_name, mous_id):
            matching_files.append(fits_file)
    
    return matching_files

def create_source_mapping(database_df: pd.DataFrame, fits_files: List) -> List:
    """
    Check whether the fits image of each source exists in every configuration
    -----------
    database_df : pandas.DataFrame
        Database with Source-MOUS mapping
    fits_files : list
        List of all FITS file paths
        
    Returns:
    --------
    list
        List of dictionaries with checking information
    """
    check_list = []
    
    print("Processing sources...")
    
    for index, row in database_df.iterrows():
        source_name = row['Source']
        
        # Create dictionary for this source
        source_dict = {
            'Source': str(source_name),
            'SGOUS': row.get('SGOUS', None),
            'GOUS': row.get('GOUS', None),
        }
        
        # Check each configuration
        configurations = ['7M', 'TM1', 'TM2']
        
        for config in configurations:
            mous_col = f'MOUS_{config}'
            mous_id = row.get(mous_col, None)

            # Create dictionary for configuration
            config_dict = {
                'CONFIG': str(config),
                'MOUS': mous_id,
            }
            
            if not mous_id:
                continue
            else:
                # Find matching FITS files  
                matching_fits = find_matching_fits(str(source_name), mous_id, fits_files)
            
            if matching_fits:
                # Only one file is expected in 'matching_fits' list
                config_dict['cont_pbcor_fits'] = matching_fits[0]  # absolute path
                config_dict['cont_check'] = 1
                if len(matching_fits) > 1: # if there is more than one fits in the list
                    config_dict[f'{config}_all_matches'] = matching_fits
            else:
                config_dict['cont_pbcor_fits'] = None
                config_dict['cont_check'] = 0
            source_dict[config] = config_dict
        
        check_list.append(source_dict)
        
    return check_list

def print_statistics(check_list: List) -> None:
    """
    Print statistics about the mapping
    -----------
    mapping_list : list
        List of mapping dictionaries
    """
    total_sources = len(check_list)
    
    stats = {
        '7M': {'found': 0, 'missing': 0},
        'TM1': {'found': 0, 'missing': 0},
        'TM2': {'found': 0, 'missing': 0}
    }
    
    for source_dict in check_list:
        for config in ['7M', 'TM1', 'TM2']:
            if source_dict[config].get('cont_check', 0) == 1:
                stats[config]['found'] += 1
            else:
                stats[config]['missing'] += 1
    
    print("\nChecking Statistics:")
    print(f"Total sources: {total_sources}")
    print("-" * 40)
    
    for config in ['7M', 'TM1', 'TM2']:
        found = stats[config]['found']
        missing = stats[config]['missing']
        percentage = (found / total_sources) * 100 if total_sources > 0 else 0
        print(f"{config:4}: {found:4} found, {missing:4} missing ({percentage:.1f}%)")
    
    # Find sources with all configurations
    all_configs = sum(1 for s in check_list if all(s[c].get('cont_check', 0) == 1 for c in ['7M', 'TM1', 'TM2']))
    print(f"\nSources with all 3 configurations: {all_configs}/{total_sources} ({all_configs/total_sources*100:.1f}%)")

def save_to_json(check_list: List, output_file: str) -> None:
    """
    Save mapping list to JSON file
    -----------
    mapping_list : list
        List of mapping dictionaries
    output_file : str
        Output JSON file path
    """
    try:
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(check_list, f, indent=2, ensure_ascii=False)
        print(f"\nMapping saved to {output_file}")
        print(f"Total entries: {len(check_list)}")
    except Exception as e:
        print(f"Error saving JSON file: {e}")

def save_summary_csv(check_list: List, output_file: str) -> None:
    """
    Save a summary CSV for easy viewing
    -----------
    mapping_list : list
        List of mapping dictionaries
    output_file : str
        Output CSV file path
    """
    try:
        # Create a simplified version for CSV
        summary_data = []
        
        for source_dict in check_list:
            summary_row = {
                'Source': source_dict['Source'],
                'SGOUS': source_dict.get('SGOUS'),
                'GOUS': source_dict.get('GOUS'),
                '7M_check': source_dict['7M'].get('cont_check', 0),
                'TM1_check': source_dict['TM1'].get('cont_check', 0),
                'TM2_check': source_dict['TM2'].get('cont_check', 0),
                'MOUS_7M': source_dict['7M'].get('MOUS'),
                'MOUS_TM1': source_dict['TM1'].get('MOUS'),
                'MOUS_TM2': source_dict['TM2'].get('MOUS'),
                '7M_fits': os.path.basename(source_dict['7M']['cont_pbcor_fits']) if source_dict['7M'].get('cont_pbcor_fits') else None,
                'TM1_fits': os.path.basename(source_dict['TM1']['cont_pbcor_fits']) if source_dict['TM1'].get('cont_pbcor_fits') else None,
                'TM2_fits': os.path.basename(source_dict['TM2']['cont_pbcor_fits']) if source_dict['TM2'].get('cont_pbcor_fits') else None,
            }
            summary_data.append(summary_row)
        
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv(output_file, index=False)
        print(f"Summary saved to {output_file}")
        
    except Exception as e:
        print(f"Error saving CSV file: {e}")

if __name__ == "__main__":

    database_file = "tables/database.xlsx"  # Path to database.xlsx
    fits_directory = "../data16/"  # Directory containing FITS files
    output_json = "almagal_source_mous_checking.json"
    output_csv = "almagal_source_mous_checking.csv"
    
    # Read database
    print("Reading database.xlsx")
    database_df = read_database(database_file)
    if database_df is None:
        print("Failed to read database.xlsx Exiting.")
        sys.exit(1)
    
    # Find FITS files
    print(f"\n Finding FITS files in {fits_directory}...")
    fits_files = find_fits_files(fits_directory)
    if not fits_files:
        print("No FITS files found. Exiting.")
        sys.exit(1)
    
    # Show sample FITS filenames for verification
    # print("\nSample FITS filenames:")
    # for i, fits_file in enumerate(fits_files[:5]):
    #     print(f"  {i+1}. {os.path.basename(fits_file)}")
    # if len(fits_files) > 5:
    #     print(f"  ... and {len(fits_files) - 5} more")
    
    # Checking
    print(f"\nStart to check source-MOUS-FITS ")
    check_list = create_source_mapping(database_df, fits_files)
    
    # Result Statistics
    print_statistics(check_list)
    
    # Save results
    save_to_json(check_list, output_json)
    save_summary_csv(check_list, output_csv)
    print(f"\nResults has been stored in JSON file {output_json} and CSV file {output_csv}")
    
    print("\nChekcing completed successfully!")
    
