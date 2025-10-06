import csv
import os
import sys
from time import sleep
import pyvo

def read_missing_mous(csv_file: str) -> list:
    """
    Read missing MOUS IDs from CSV file
    -----------
    csv_file : str
        Path to CSV file containing missing MOUS IDs
        
    Returns:
    --------
    list
        List of MOUS IDs with 'uid___' prefix preserved
    """
    mous_ids = []
    
    try:
        with open(csv_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                mous_id = row['MOUS_ID'].strip()
                # Keep the 'uid___' prefix
                mous_ids.append(mous_id)
        
        print(f"Read {len(mous_ids)} MOUS IDs from {csv_file}")
        return mous_ids
        
    except Exception as e:
        print(f"Error reading {csv_file}: {e}")
        return []

def query_mous_data(mous_id: str, index: int) -> dict:
    """
    Query data size and download links for a single MOUS ID using pyvo
    -----------
    mous_id : str
        MOUS ID with 'uid___' prefix
    index : int
        Index for server selection (not used in this implementation)
        
    Returns:
    --------
    dict
        Dictionary containing query results
    """
    result = {
        'MOUS_ID': mous_id,
        'Main_Data_URL': None,
        'Main_Data_Size_GB': None,
        'Auxiliary_Data_URL': None,
        'Auxiliary_Data_Size_GB': None,
        'README_URL': None,
        'Status': 'Failed',
        'Error': None
    }
    
    try:
        print(f"\nQuerying MOUS: {mous_id}")
        
        # Query ALMA datalink service using pyvo
        datalink_url = f"https://almascience.eso.org/datalink/sync?ID={mous_id}"
        datalink = pyvo.dal.adhoc.DatalinkResults.from_result_url(datalink_url)
        
        if datalink and len(datalink) > 0:
            # Extract access URLs and content lengths
            access_urls = datalink['access_url'].tolist()
            content_lengths = datalink['content_length'].tolist()
            
            # Find specific URLs and their sizes based on patterns
            main_data_url = None
            main_data_size_gb = None
            auxiliary_data_url = None
            auxiliary_data_size_gb = None
            readme_url = None
            
            for i, url in enumerate(access_urls):
                if url and isinstance(url, str):
                    content_length = content_lengths[i] if i < len(content_lengths) else None
                    
                    if url.endswith('of_001.tar'):
                        main_data_url = url
                        if content_length is not None:
                            main_data_size_gb = round(content_length / (1024**3), 2)
                    elif url.endswith('_auxiliary.tar'):
                        auxiliary_data_url = url
                        if content_length is not None:
                            auxiliary_data_size_gb = round(content_length / (1024**3), 2)
                    elif url.endswith('README.txt'):
                        readme_url = url
            
            result['Main_Data_URL'] = main_data_url
            result['Main_Data_Size_GB'] = main_data_size_gb
            result['Auxiliary_Data_URL'] = auxiliary_data_url
            result['Auxiliary_Data_Size_GB'] = auxiliary_data_size_gb
            result['README_URL'] = readme_url
            result['Status'] = 'Success'
            
            print(f"  ✓ Found {len(access_urls)} URLs")
            if main_data_url:
                size_str = f" ({main_data_size_gb} GB)" if main_data_size_gb else ""
                print(f"    Main data: {main_data_url}{size_str}")
            if auxiliary_data_url:
                size_str = f" ({auxiliary_data_size_gb} GB)" if auxiliary_data_size_gb else ""
                print(f"    Auxiliary data: {auxiliary_data_url}{size_str}")
            if readme_url:
                print(f"    README: {readme_url}")
            
        else:
            result['Status'] = 'No Data Found'
            print(f"  ✗ No data found for {mous_id}")
            
    except Exception as e:
        result['Error'] = str(e)
        result['Status'] = 'Error'
        print(f"  ✗ Error querying {mous_id}: {e}")
    
    return result

def query_all_missing_mous(csv_file: str, output_file: str) -> None:
    """
    Query all missing MOUS IDs and save results
    -----------
    csv_file : str
        Path to CSV file containing missing MOUS IDs
    output_file : str
        Path to output CSV file
    """
    # Read missing MOUS IDs
    missing_mous = read_missing_mous(csv_file)
    
    if not missing_mous:
        print("No MOUS IDs found to query")
        return
    
    print(f"\n{'='*60}")
    print(f"QUERYING ALMA DATA FOR {len(missing_mous)} MISSING MOUS IDs")
    print(f"{'='*60}")
    
    results = []
    successful_queries = 0
    main_data_found = 0
    auxiliary_data_found = 0
    readme_found = 0
    total_main_data_size = 0.0
    total_auxiliary_data_size = 0.0
    
    for i, mous_id in enumerate(missing_mous):
        print(f"\nProgress: {i+1}/{len(missing_mous)}")
        
        # Query data for this MOUS
        result = query_mous_data(mous_id, i)
        results.append(result)
        
        if result['Status'] == 'Success':
            successful_queries += 1
            if result['Main_Data_URL']:
                main_data_found += 1
                if result['Main_Data_Size_GB']:
                    total_main_data_size += result['Main_Data_Size_GB']
            if result['Auxiliary_Data_URL']:
                auxiliary_data_found += 1
                if result['Auxiliary_Data_Size_GB']:
                    total_auxiliary_data_size += result['Auxiliary_Data_Size_GB']
            if result['README_URL']:
                readme_found += 1
        
        # Add delay between queries to avoid overwhelming the server
        sleep(2)
    
    # Save results to CSV
    try:
        with open(output_file, 'w', newline='', encoding='utf-8') as f:
            if results:
                # Write header
                fieldnames = results[0].keys()
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                
                # Write data
                for result in results:
                    writer.writerow(result)
        
        print(f"\nResults saved to {output_file}")
        
        # Print summary
        print(f"\n{'='*60}")
        print("QUERY SUMMARY")
        print(f"{'='*60}")
        print(f"Total MOUS IDs queried: {len(missing_mous)}")
        print(f"Successful queries: {successful_queries}")
        print(f"Failed queries: {len(missing_mous) - successful_queries}")
        print(f"Main data URLs found: {main_data_found}")
        print(f"Auxiliary data URLs found: {auxiliary_data_found}")
        print(f"README URLs found: {readme_found}")
        print(f"Total main data size: {total_main_data_size:.2f} GB")
        print(f"Total auxiliary data size: {total_auxiliary_data_size:.2f} GB")
        print(f"Total data size: {total_main_data_size + total_auxiliary_data_size:.2f} GB")
        
        # Show first few results
        print(f"\nFirst 5 results:")
        for i, result in enumerate(results[:5]):
            status_icon = "✓" if result['Status'] == 'Success' else "✗"
            urls_found = []
            size_info = []
            if result['Main_Data_URL']:
                urls_found.append("Main")
                if result['Main_Data_Size_GB']:
                    size_info.append(f"Main: {result['Main_Data_Size_GB']} GB")
            if result['Auxiliary_Data_URL']:
                urls_found.append("Aux")
                if result['Auxiliary_Data_Size_GB']:
                    size_info.append(f"Aux: {result['Auxiliary_Data_Size_GB']} GB")
            if result['README_URL']:
                urls_found.append("README")
            urls_str = ", ".join(urls_found) if urls_found else "None"
            size_str = f" - Sizes: {', '.join(size_info)}" if size_info else ""
            print(f"  {i+1}. {result['MOUS_ID']} - {status_icon} {result['Status']} - URLs: {urls_str}{size_str}")
        
        if len(results) > 5:
            print(f"  ... and {len(results) - 5} more")
            
    except Exception as e:
        print(f"Error saving results: {e}")

def main():
    """Main function"""
    # File paths
    input_csv = "missing_mous_list.csv"
    output_csv = "missing_mous_query_results.csv"
    
    # Check if input file exists
    if not os.path.exists(input_csv):
        print(f"Error: {input_csv} not found")
        print("Please run extract_missing_mous.py first to generate the missing MOUS list")
        sys.exit(1)
    
    print("ALMA Missing MOUS Data Query Tool")
    print("="*40)
    print(f"Input file: {input_csv}")
    print(f"Output file: {output_csv}")
    
    # Query all missing MOUS
    query_all_missing_mous(input_csv, output_csv)
    
    print(f"\nQuery completed! Results saved to {output_csv}")

if __name__ == "__main__":
    main()
