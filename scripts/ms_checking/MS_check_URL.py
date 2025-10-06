import pandas as pd
import json
import os
import sys
import glob
from pathlib import Path
from typing import List, Optional
import requests
from bs4 import BeautifulSoup
import time
import re

def make_request_with_retry(url: str, max_retries: int = 3, timeout: int = 120) -> Optional[requests.Response]:
    """
    Make HTTP request with retry logic
    -----------
    url : str
        URL to request
    max_retries : int
        Maximum number of retries (default: 3)
    timeout : int
        Request timeout in seconds (default: 120)
        
    Returns:
    --------
    requests.Response or None
        Response object if successful, None if all retries failed
    """
    for attempt in range(max_retries):
        try:
            response = requests.get(url, timeout=timeout)
            response.raise_for_status()
            return response
        except requests.RequestException as e:
            print(f"    Attempt {attempt + 1}/{max_retries} failed for {url}: {e}")
            if attempt < max_retries - 1:
                wait_time = 2 ** attempt  # Exponential backoff
                print(f"    Retrying in {wait_time} seconds...")
                time.sleep(wait_time)
            else:
                print(f"    All {max_retries} attempts failed for {url}")
                return None
    return None

def read_database(database_file: str) -> Optional[pd.DataFrame]:
    """
    Read the ALMAGAL database.xlsx and prepare for MS data checking
    -----------
    database_file : str
        Path to database.xlsx file
        
    Returns:
    --------
    pandas.DataFrame
        Database with MOUS mapping only (Source column removed, duplicates removed)
    """
    try:
        import warnings # warnings from openpyxl is suppressed

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning, module="openpyxl")
            df = pd.read_excel(database_file, engine='openpyxl')

        print(f"Successfully read database with {len(df)} targets")
        print(f"Original database columns: {df.columns.tolist()}")
        
        # Remove Source column as it's not needed for MS data checking
        if 'Source' in df.columns:
            df = df.drop('Source', axis=1)
            print("Removed 'Source' column from database")
        
        # Remove duplicate rows based on MOUS combinations
        initial_count = len(df)
        df = df.drop_duplicates()
        final_count = len(df)
        print(f"Removed {initial_count - final_count} duplicate rows")
        print(f"Final database has {final_count} unique MOUS combinations")
        
        return df
    except Exception as e:
        print(f"Error reading database file: {e}")
        return None

def find_ms_data(base_url: str) -> tuple[List, List]:
    """
    Find all MS data on the remote server
    -----------
    base_url : str
        Base URL to scan for MS data
        
    Returns:
    --------
    tuple
        (ms_data_list, tar_files_list) - MS data info and tar.gz files found
    """
    ms_data = []
    tar_files = []
    
    try:
        # Get the main page content
        response = make_request_with_retry(base_url)
        if response is None:
            print(f"Failed to access {base_url} after retries")
            return []
        soup = BeautifulSoup(response.content, 'html.parser')
        
        # Find all links on the main page
        links = soup.find_all('a', href=True)
        
        print(f"Scanning {base_url} for MS data directories...")
        
        # Process each link to find uid___<MOUS> patterns
        for link in links:
            href = link['href']
            if href.startswith('uid___') and href.endswith('_calibrated_final.ms/'):
                # This is a MOUS MS directory
                mous_url = base_url.rstrip('/') + '/' + href
                print(f"  Found MOUS MS directory: {href}")
                
                # Extract MOUS from href (e.g., "uid___A001_X1467_X1d3_calibrated_final.ms/" -> "uid___A001_X1467_X1d3")
                mous = href.replace('_calibrated_final.ms/', '').rstrip('/')
                
                # Check if this MS directory contains table.info
                if check_for_table_info(mous_url):
                    ms_info = {
                        'MOUS': mous,
                        'ms_url': mous_url
                    }
                    ms_data.append(ms_info)
                    print(f"    Found MS data with table.info: {href}")
                else:
                    print(f"    No table.info found in: {href}")
            elif href.startswith('uid___') and href.endswith('/'):
                # This is a regular MOUS directory, check for subdirectories
                mous_url = base_url.rstrip('/') + '/' + href
                print(f"  Found MOUS directory: {href}")
                
                # Extract MOUS from href
                mous = href.rstrip('/')
                
                # Scan this MOUS directory for MS data
                mous_ms_data, mous_tar_files = scan_mous_directory(mous_url, mous)
                ms_data.extend(mous_ms_data)
                tar_files.extend(mous_tar_files)
        
        # Remove duplicates based on ms_url
        seen_urls = set()
        unique_ms_data = []
        for ms_info in ms_data:
            if ms_info['ms_url'] not in seen_urls:
                seen_urls.add(ms_info['ms_url'])
                unique_ms_data.append(ms_info)
        
        print(f"Found {len(unique_ms_data)} MS data entries")
        print(f"Found {len(tar_files)} tar.gz files")
        return unique_ms_data, tar_files
        
    except requests.RequestException as e:
        print(f"Error accessing {base_url}: {e}")
        return []
    except Exception as e:
        print(f"Unexpected error: {e}")
        return []

def scan_mous_directory(mous_url: str, mous: str) -> tuple[List, List]:
    """
    Scan a MOUS directory for MS data and tar.gz files
    -----------
    mous_url : str
        URL of the MOUS directory
    mous : str
        MOUS ID
        
    Returns:
    --------
    tuple
        (ms_data_list, tar_files_list) - MS data info and tar.gz files found
    """
    ms_data = []
    tar_files = []
    
    try:
        response = make_request_with_retry(mous_url)
        if response is None:
            print(f"    Failed to access {mous_url} after retries")
            return []
        soup = BeautifulSoup(response.content, 'html.parser')
        
        # Find all links in the MOUS directory
        links = soup.find_all('a', href=True)
        
        # First, look for .ms directories
        for link in links:
            href = link['href']
            if href.endswith('.ms/') or href.endswith('.ms'):
                ms_dir_url = mous_url.rstrip('/') + '/' + href.rstrip('/')
                print(f"    Found MS directory: {href}")
                
                # Check if this MS directory contains table.info
                if check_for_table_info(ms_dir_url):
                    ms_info = {
                        'MOUS': mous,
                        'ms_url': ms_dir_url
                    }
                    ms_data.append(ms_info)
                    print(f"      Found MS data with table.info: {href}")
                    # Break here since we found table.info - no need to check other subdirectories
                    return ms_data, tar_files
                else:
                    print(f"      No table.info found in: {href}")
        
        # Look for tar.gz files
        for link in links:
            href = link['href']
            if href.endswith('.tar.gz'):
                tar_url = mous_url.rstrip('/') + '/' + href
                tar_info = {
                    'MOUS': mous,
                    'tar_url': tar_url,
                    'filename': href
                }
                tar_files.append(tar_info)
                print(f"    Found tar.gz file: {href}")
        
        # Only check subdirectories if no .ms directory with table.info was found
        for link in links:
            href = link['href']
            if href.endswith('/') and not href.startswith('..'):
                subdir_url = mous_url.rstrip('/') + '/' + href
                print(f"    Checking subdirectory: {href}")
                subdir_ms_data, subdir_tar_files = scan_mous_directory(subdir_url, mous)
                ms_data.extend(subdir_ms_data)
                tar_files.extend(subdir_tar_files)
                # If we found MS data in subdirectory, break to avoid further scanning
                if subdir_ms_data:
                    break
        
                
    except requests.RequestException as e:
        print(f"    Error accessing {mous_url}: {e}")
    except Exception as e:
        print(f"    Unexpected error in MOUS {mous_url}: {e}")
    
    return ms_data, tar_files

def check_for_table_info(ms_dir_url: str) -> bool:
    """
    Check if a MS directory contains table.info file
    -----------
    ms_dir_url : str
        URL of the MS directory
        
    Returns:
    --------
    bool
        True if table.info is found, False otherwise
    """
    try:
        response = make_request_with_retry(ms_dir_url)
        if response is None:
            return False
        soup = BeautifulSoup(response.content, 'html.parser')
        
        # Look for table.info file
        links = soup.find_all('a', href=True)
        for link in links:
            href = link['href']
            if href == 'table.info':
                print(f"        Found table.info in {ms_dir_url}")
                return True
        
        return False
        
    except Exception as e:
        print(f"        Error checking table.info in {ms_dir_url}: {e}")
        return False

def save_ms_data_to_csv(ms_data: List, output_file: str) -> None:
    """
    Save MS data information to CSV file
    -----------
    ms_data : list
        List of MS data information dictionaries
    output_file : str
        Output CSV file path
    """
    try:
        if not ms_data:
            print("No MS data to save")
            return
            
        # Create DataFrame from the list of dictionaries
        df = pd.DataFrame(ms_data)
        
        # Save to CSV
        df.to_csv(output_file, index=False)
        print(f"MS data information saved to {output_file}")
        print(f"Total MS data entries: {len(ms_data)}")
        
        # Show sample of the data
        print("\nSample of saved data:")
        print(df.head())
        
    except Exception as e:
        print(f"Error saving MS data CSV: {e}")

def save_tar_files_to_csv(tar_files: List, output_file: str) -> None:
    """
    Save tar.gz files information to CSV file
    -----------
    tar_files : list
        List of tar.gz file information dictionaries
    output_file : str
        Output CSV file path
    """
    try:
        if not tar_files:
            print("No tar.gz files to save")
            return
            
        # Create DataFrame from the list of dictionaries
        df = pd.DataFrame(tar_files)
        
        # Save to CSV
        df.to_csv(output_file, index=False)
        print(f"Tar.gz files information saved to {output_file}")
        print(f"Total tar.gz files: {len(tar_files)}")
        
        # Show sample of the data
        print("\nSample of tar.gz files found:")
        print(df.head())
        
    except Exception as e:
        print(f"Error saving tar.gz files CSV: {e}")

def load_ms_data_from_csv(csv_file: str) -> List:
    """
    Load MS data information from CSV file
    -----------
    csv_file : str
        Path to CSV file containing MS data information
        
    Returns:
    --------
    list
        List of MS data information dictionaries, or empty list if file doesn't exist
    """
    try:
        if not os.path.exists(csv_file):
            print(f"CSV file {csv_file} does not exist")
            return []
            
        df = pd.read_csv(csv_file)
        ms_data = df.to_dict('records')
        print(f"Loaded {len(ms_data)} MS data entries from {csv_file}")
        return ms_data
        
    except Exception as e:
        print(f"Error loading MS data from CSV: {e}")
        return []

def find_matching_ms(mous_id: str, ms_data: List) -> List:
    """
    Find MS data that matches the MOUS ID
    -----------
    mous_id : str
        MOUS ID to search for
    ms_data : list
        List of all MS data information dictionaries
        
    Returns:
    --------
    list
        List of matching MS data information dictionaries
    """
    matching_ms = []
    
    if pd.isna(mous_id):
        return matching_ms
    
    mous = str(mous_id).strip()
    
    for ms_info in ms_data:
        if ms_info['MOUS'] == mous:
            matching_ms.append(ms_info)
    
    return matching_ms

def create_source_mapping_ms(database_df: pd.DataFrame, ms_data: List) -> List:
    """
    Check whether the MS data exists for each MOUS configuration
    -----------
    database_df : pandas.DataFrame
        Database with MOUS mapping (Source column removed)
    ms_data : list
        List of all MS data information dictionaries
        
    Returns:
    --------
    list
        List of dictionaries with checking information
    """
    check_list = []
    
    print("Processing MOUS configurations for MS data...")
    
    for index, row in database_df.iterrows():
        # Create dictionary for this MOUS combination
        mous_dict = {
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
                config_dict['ms_data_url'] = None
                config_dict['ms_check'] = 0
            else:
                # Find matching MS data  
                matching_ms = find_matching_ms(mous_id, ms_data)
            
                if matching_ms:
                    # Only one MS data is expected in 'matching_ms' list
                    config_dict['ms_data_url'] = matching_ms[0]['ms_url']  # URL
                    config_dict['ms_check'] = 1
                    if len(matching_ms) > 1: # if there is more than one MS data in the list
                        config_dict[f'{config}_all_matches'] = [m['ms_url'] for m in matching_ms]
                else:
                    config_dict['ms_data_url'] = None
                    config_dict['ms_check'] = 0
            mous_dict[config] = config_dict
        
        check_list.append(mous_dict)
        
    return check_list

def print_statistics_ms(check_list: List) -> None:
    """
    Print statistics about the MS data mapping
    -----------
    check_list : list
        List of mapping dictionaries
    """
    total_mous_combinations = len(check_list)
    
    stats = {
        '7M': {'found': 0, 'missing': 0},
        'TM1': {'found': 0, 'missing': 0},
        'TM2': {'found': 0, 'missing': 0}
    }
    
    for mous_dict in check_list:
        for config in ['7M', 'TM1', 'TM2']:
            if mous_dict[config].get('ms_check', 0) == 1:
                stats[config]['found'] += 1
            else:
                stats[config]['missing'] += 1
    
    print("\nMS Data Checking Statistics:")
    print(f"Total MOUS combinations: {total_mous_combinations}")
    print("-" * 40)
    
    for config in ['7M', 'TM1', 'TM2']:
        found = stats[config]['found']
        missing = stats[config]['missing']
        percentage = (found / total_mous_combinations) * 100 if total_mous_combinations > 0 else 0
        print(f"{config:4}: {found:4} found, {missing:4} missing ({percentage:.1f}%)")
    
    # Find MOUS combinations with all configurations
    all_configs = sum(1 for s in check_list if all(s[c].get('ms_check', 0) == 1 for c in ['7M', 'TM1', 'TM2']))
    print(f"\nMOUS combinations with all 3 MS configurations: {all_configs}/{total_mous_combinations} ({all_configs/total_mous_combinations*100:.1f}%)")

def save_to_json_ms(check_list: List, output_file: str) -> None:
    """
    Save mapping list to JSON file
    -----------
    check_list : list
        List of mapping dictionaries
    output_file : str
        Output JSON file path
    """
    try:
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(check_list, f, indent=2, ensure_ascii=False)
        print(f"\nMS mapping saved to {output_file}")
        print(f"Total entries: {len(check_list)}")
    except Exception as e:
        print(f"Error saving JSON file: {e}")

def save_summary_csv_ms(check_list: List, output_file: str) -> None:
    """
    Save a summary CSV for easy viewing
    -----------
    check_list : list
        List of mapping dictionaries
    output_file : str
        Output CSV file path
    """
    try:
        # Create a simplified version for CSV
        summary_data = []
        
        for mous_dict in check_list:
            summary_row = {
                'SGOUS': mous_dict.get('SGOUS'),
                'GOUS': mous_dict.get('GOUS'),
                '7M_ms_check': mous_dict['7M'].get('ms_check', 0),
                'TM1_ms_check': mous_dict['TM1'].get('ms_check', 0),
                'TM2_ms_check': mous_dict['TM2'].get('ms_check', 0),
                'MOUS_7M': mous_dict['7M'].get('MOUS'),
                'MOUS_TM1': mous_dict['TM1'].get('MOUS'),
                'MOUS_TM2': mous_dict['TM2'].get('MOUS'),
                '7M_ms_url': mous_dict['7M'].get('ms_data_url'),
                'TM1_ms_url': mous_dict['TM1'].get('ms_data_url'),
                'TM2_ms_url': mous_dict['TM2'].get('ms_data_url'),
            }
            summary_data.append(summary_row)
        
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv(output_file, index=False)
        print(f"MS summary saved to {output_file}")
        
    except Exception as e:
        print(f"Error saving CSV file: {e}")

def report_tar_files(tar_files: List) -> None:
    """
    Report tar.gz files found during scanning
    -----------
    tar_files : list
        List of tar.gz file information dictionaries
    """
    if not tar_files:
        print("\nNo tar.gz files found.")
        return
    
    print(f"\n{'='*60}")
    print(f"TAR.GZ FILES FOUND ({len(tar_files)} files)")
    print(f"{'='*60}")
    
    # Group by MOUS for better organization
    mous_tar_files = {}
    for tar_info in tar_files:
        mous = tar_info['MOUS']
        if mous not in mous_tar_files:
            mous_tar_files[mous] = []
        mous_tar_files[mous].append(tar_info)
    
    for mous, files in mous_tar_files.items():
        print(f"\nMOUS: {mous}")
        print("-" * 40)
        for tar_info in files:
            print(f"  File: {tar_info['filename']}")
            print(f"  URL:  {tar_info['tar_url']}")
            print()
    
    print(f"\nTotal tar.gz files found: {len(tar_files)}")
    print("These files may contain MS data that needs to be extracted.")
    print(f"{'='*60}")

if __name__ == "__main__":

    # Configuration variables
    database_file = "tables/database.xlsx"  # Path to database.xlsx
    base_url = "https://casdc.china-vo.org/mirror/ALMAGAL/2019.1.00195.L/MS"  # Remote server URL for MS data
    output_json = "almagal_source_mous_ms_checking.json"
    output_csv = "almagal_source_mous_ms_checking.csv"
    ms_data_csv = "ms_data_found.csv"  # CSV file for MS data found
    tar_files_csv = "tar_files_found.csv"  # CSV file for tar.gz files found
    
    # Read database
    print("Reading database.xlsx")
    database_df = read_database(database_file)
    if database_df is None:
        print("Failed to read database.xlsx Exiting.")
        sys.exit(1)
    
    # Check if MS data CSV already exists
    print(f"\nChecking for existing MS data CSV: {ms_data_csv}")
    ms_data = load_ms_data_from_csv(ms_data_csv)
    
    if not ms_data:
        # CSV doesn't exist or is empty, scan the remote server
        print(f"Scanning remote server {base_url} for MS data and tar.gz files...")
        ms_data, tar_files = find_ms_data(base_url)
        if not ms_data and not tar_files:
            print("No MS data or tar.gz files found. Exiting.")
            sys.exit(1)
        
        # Save MS data information to CSV for future use
        print(f"\nSaving MS data information to {ms_data_csv}...")
        save_ms_data_to_csv(ms_data, ms_data_csv)
        
        # Save tar.gz files information to CSV
        print(f"\nSaving tar.gz files information to {tar_files_csv}...")
        save_tar_files_to_csv(tar_files, tar_files_csv)
    else:
        print(f"Using existing MS data from {ms_data_csv}")
        print(f"Found {len(ms_data)} MS data entries from CSV")
        
        # Load tar files if available
        tar_files = []
        if os.path.exists(tar_files_csv):
            try:
                df = pd.read_csv(tar_files_csv)
                tar_files = df.to_dict('records')
                print(f"Found {len(tar_files)} tar.gz files from CSV")
            except Exception as e:
                print(f"Error loading tar files CSV: {e}")
                tar_files = []
    
    # Checking
    print(f"\nStart to check MOUS-MS data")
    check_list = create_source_mapping_ms(database_df, ms_data)
    
    # Result Statistics
    print_statistics_ms(check_list)
    
    # Save results
    save_to_json_ms(check_list, output_json)
    save_summary_csv_ms(check_list, output_csv)
    print(f"\nResults has been stored in JSON file {output_json} and CSV file {output_csv}")
    print(f"MS data information has been stored in CSV file {ms_data_csv}")
    
    # Report tar.gz files found
    report_tar_files(tar_files)
    
    print("\nMS data checking completed successfully!")
