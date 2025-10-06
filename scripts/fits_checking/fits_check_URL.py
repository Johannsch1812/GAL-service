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

def find_fits_files(base_url: str, fits_pattern: str = "cont.I.pbcor.fits") -> List:
    """
    Find all FITS files ending with specified pattern on the remote server
    -----------
    base_url : str
        Base URL to scan for FITS files
    fits_pattern : str
        FITS file pattern to match (default: "cont.I.pbcor.fits")
        
    Returns:
    --------
    list
        List of dictionaries with FITS file information (SOUS, GOUS, MOUS, fits_url)
    """
    fits_files = []
    
    try:
        # Get the main page content
        response = make_request_with_retry(base_url)
        if response is None:
            print(f"Failed to access {base_url} after retries")
            return []
        soup = BeautifulSoup(response.content, 'html.parser')
        
        # Find all links on the main page
        links = soup.find_all('a', href=True)
        
        print(f"Scanning {base_url} for science goal directories...")
        
        # Process each link to find science goal directories
        for link in links:
            href = link['href']
            if href.endswith('/') and 'science_goal' in href:
                # This is a science goal directory
                science_goal_url = base_url.rstrip('/') + '/' + href
                print(f"  Found science goal: {href}")
                
                # Extract SOUS from href (e.g., "science_goal.uid___A001_X1467_X1d3/" -> "uid___A001_X1467_X1d3")
                sous = href.replace('science_goal.', '').rstrip('/')
                
                # Scan this science goal directory
                science_goal_fits = scan_science_goal_directory(science_goal_url, sous, fits_pattern)
                fits_files.extend(science_goal_fits)
        
        # Remove duplicates based on fits_url
        seen_urls = set()
        unique_fits_files = []
        for fits_info in fits_files:
            if fits_info['fits_url'] not in seen_urls:
                seen_urls.add(fits_info['fits_url'])
                unique_fits_files.append(fits_info)
        
        print(f"Found {len(unique_fits_files)} FITS files ending with '{fits_pattern}'")
        return unique_fits_files
        
    except requests.RequestException as e:
        print(f"Error accessing {base_url}: {e}")
        return []
    except Exception as e:
        print(f"Unexpected error: {e}")
        return []

def save_fits_files_to_csv(fits_files: List, output_file: str) -> None:
    """
    Save FITS files information to CSV file
    -----------
    fits_files : list
        List of FITS file information dictionaries
    output_file : str
        Output CSV file path
    """
    try:
        if not fits_files:
            print("No FITS files to save")
            return
            
        # Create DataFrame from the list of dictionaries
        df = pd.DataFrame(fits_files)
        
        # Save to CSV
        df.to_csv(output_file, index=False)
        print(f"FITS files information saved to {output_file}")
        print(f"Total FITS files: {len(fits_files)}")
        
        # Show sample of the data
        print("\nSample of saved data:")
        print(df.head())
        
    except Exception as e:
        print(f"Error saving FITS files CSV: {e}")

def load_fits_files_from_csv(csv_file: str) -> List:
    """
    Load FITS files information from CSV file
    -----------
    csv_file : str
        Path to CSV file containing FITS files information
        
    Returns:
    --------
    list
        List of FITS file information dictionaries, or empty list if file doesn't exist
    """
    try:
        if not os.path.exists(csv_file):
            print(f"CSV file {csv_file} does not exist")
            return []
            
        df = pd.read_csv(csv_file)
        fits_files = df.to_dict('records')
        print(f"Loaded {len(fits_files)} FITS files from {csv_file}")
        return fits_files
        
    except Exception as e:
        print(f"Error loading FITS files from CSV: {e}")
        return []

def scan_science_goal_directory(science_goal_url: str, sous: str, fits_pattern: str) -> List:
    """
    Scan a science goal directory for group and member subdirectories
    -----------
    science_goal_url : str
        URL of the science goal directory
    sous : str
        Science goal UID (SOUS)
    fits_pattern : str
        FITS file pattern to match
        
    Returns:
    --------
    list
        List of FITS file information dictionaries found in this science goal
    """
    fits_files = []
    
    try:
        response = make_request_with_retry(science_goal_url)
        if response is None:
            print(f"    Failed to access {science_goal_url} after retries")
            return []
        soup = BeautifulSoup(response.content, 'html.parser')
        
        # Find group directories
        links = soup.find_all('a', href=True)
        
        for link in links:
            href = link['href']
            if href.endswith('/') and 'group.' in href:
                group_url = science_goal_url.rstrip('/') + '/' + href
                print(f"    Found group: {href}")
                
                # Extract GOUS from href (e.g., "group.uid___A001_X1467_X1d4/" -> "uid___A001_X1467_X1d4")
                gous = href.replace('group.', '').rstrip('/')
                
                # Scan this group directory
                group_fits = scan_group_directory(group_url, sous, gous, fits_pattern)
                fits_files.extend(group_fits)
                
    except requests.RequestException as e:
        print(f"    Error accessing {science_goal_url}: {e}")
    except Exception as e:
        print(f"    Unexpected error in science goal {science_goal_url}: {e}")
    
    return fits_files

def scan_group_directory(group_url: str, sous: str, gous: str, fits_pattern: str) -> List:
    """
    Scan a group directory for member subdirectories
    -----------
    group_url : str
        URL of the group directory
    sous : str
        Science goal UID (SOUS)
    gous : str
        Group UID (GOUS)
    fits_pattern : str
        FITS file pattern to match
        
    Returns:
    --------
    list
        List of FITS file information dictionaries found in this group
    """
    fits_files = []
    
    try:
        response = make_request_with_retry(group_url)
        if response is None:
            print(f"      Failed to access {group_url} after retries")
            return []
        soup = BeautifulSoup(response.content, 'html.parser')
        
        # Find member directories
        links = soup.find_all('a', href=True)
        
        for link in links:
            href = link['href']
            if href.endswith('/') and 'member.' in href:
                member_url = group_url.rstrip('/') + '/' + href
                print(f"      Found member: {href}")
                
                # Extract MOUS from href (e.g., "member.uid___A001_X1467_X1d5/" -> "uid___A001_X1467_X1d5")
                mous = href.replace('member.', '').rstrip('/')
                
                # Scan this member directory
                member_fits = scan_member_directory(member_url, sous, gous, mous, fits_pattern)
                fits_files.extend(member_fits)
                
    except requests.RequestException as e:
        print(f"      Error accessing {group_url}: {e}")
    except Exception as e:
        print(f"      Unexpected error in group {group_url}: {e}")
    
    return fits_files

def scan_member_directory(member_url: str, sous: str, gous: str, mous: str, fits_pattern: str) -> List:
    """
    Scan a member directory for product subdirectory and FITS files
    -----------
    member_url : str
        URL of the member directory
    sous : str
        Science goal UID (SOUS)
    gous : str
        Group UID (GOUS)
    mous : str
        Member UID (MOUS)
    fits_pattern : str
        FITS file pattern to match
        
    Returns:
    --------
    list
        List of FITS file information dictionaries found in this member
    """
    fits_files = []
    
    try:
        response = make_request_with_retry(member_url)
        if response is None:
            print(f"        Failed to access {member_url} after retries")
            return []
        soup = BeautifulSoup(response.content, 'html.parser')
        
        # Look for product directory
        links = soup.find_all('a', href=True)
        
        for link in links:
            href = link['href']
            if href == 'product/':
                product_url = member_url.rstrip('/') + '/' + href
                print(f"        Found product directory")
                
                # Scan product directory for FITS files
                product_fits = scan_product_directory(product_url, sous, gous, mous, fits_pattern)
                fits_files.extend(product_fits)
                break
                # The 'break' is used to stop searching after the first 'product/' directory is found,
                # since only one product directory is expected per member. 
                
    except requests.RequestException as e:
        print(f"        Error accessing {member_url}: {e}")
    except Exception as e:
        print(f"        Unexpected error in member {member_url}: {e}")
    
    return fits_files

def scan_product_directory(product_url: str, sous: str, gous: str, mous: str, fits_pattern: str) -> List:
    """
    Scan a product directory for FITS files ending with specified pattern
    -----------
    product_url : str
        URL of the product directory
    sous : str
        Science goal UID (SOUS)
    gous : str
        Group UID (GOUS)
    mous : str
        Member UID (MOUS)
    fits_pattern : str
        FITS file pattern to match
        
    Returns:
    --------
    list
        List of FITS file information dictionaries found in this product directory
    """
    fits_files = []
    
    try:
        response = make_request_with_retry(product_url)
        if response is None:
            print(f"          Failed to access {product_url} after retries")
            return []
        soup = BeautifulSoup(response.content, 'html.parser')
        
        # Find FITS files ending with specified pattern
        links = soup.find_all('a', href=True)
        
        for link in links:
            href = link['href']
            if href.endswith(fits_pattern):
                fits_url = product_url.rstrip('/') + '/' + href
                fits_info = {
                    'SOUS': sous,
                    'GOUS': gous,
                    'MOUS': mous,
                    'fits_url': fits_url
                }
                fits_files.append(fits_info)
                print(f"          Found FITS: {href}")
                
    except requests.RequestException as e:
        print(f"          Error accessing {product_url}: {e}")
    except Exception as e:
        print(f"          Unexpected error in product {product_url}: {e}")
    
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
        List of all FITS file information dictionaries
        
    Returns:
    --------
    list
        List of matching FITS file information dictionaries
    """
    matching_files = []
    
    for fits_info in fits_files:
        # Extract filename from URL
        filename = os.path.basename(fits_info['fits_url'])
        if check_fits_contains_terms(filename, source_name, mous_id):
            matching_files.append(fits_info)
    
    return matching_files

def create_source_mapping(database_df: pd.DataFrame, fits_files: List) -> List:
    """
    Check whether the fits image of each source exists in every configuration
    -----------
    database_df : pandas.DataFrame
        Database with Source-MOUS mapping
    fits_files : list
        List of all FITS file URLs
        
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
                config_dict['cont_pbcor_fits'] = matching_fits[0]['fits_url']  # URL
                config_dict['cont_check'] = 1
                if len(matching_fits) > 1: # if there is more than one fits in the list
                    config_dict[f'{config}_all_matches'] = [f['fits_url'] for f in matching_fits]
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

    # Configuration variables
    database_file = "tables/database.xlsx"  # Path to database.xlsx
    base_url = "https://casdc.china-vo.org/mirror/ALMAGAL/2019.1.00195.L/QA2"  # Remote server URL
    output_json = "almagal_QA2_cube_checking.json"
    output_csv = "almagal_QA2_cube_checking.csv"

    fits_pattern = "cube.I.pbcor.fits"  # FITS file pattern to match
    fits_files_csv = "fits_cube_files_found.csv"  # CSV file for FITS files found
    
    # Read database
    print("Reading database.xlsx")
    database_df = read_database(database_file)
    if database_df is None:
        print("Failed to read database.xlsx Exiting.")
        sys.exit(1)
    
    # Check if FITS files CSV already exists
    print(f"\nChecking for existing FITS files CSV: {fits_files_csv}")
    fits_files = load_fits_files_from_csv(fits_files_csv)
    
    if not fits_files:
        # CSV doesn't exist or is empty, scan the remote server
        print(f"Scanning remote server {base_url} for FITS files with pattern '{fits_pattern}'...")
        fits_files = find_fits_files(base_url, fits_pattern)
        if not fits_files:
            print("No FITS files found. Exiting.")
            sys.exit(1)
        
        # Save FITS files information to CSV for future use
        print(f"\nSaving FITS files information to {fits_files_csv}...")
        save_fits_files_to_csv(fits_files, fits_files_csv)
    else:
        print(f"Using existing FITS files from {fits_files_csv}")
        print(f"Found {len(fits_files)} FITS files from CSV")
    
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
    print(f"FITS files information has been stored in CSV file {fits_files_csv}")
    
    print("\nChecking completed successfully!")
    
