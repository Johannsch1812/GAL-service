import csv
import os
import sys
from time import sleep
import urllib.request
import urllib.parse
import json

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

def query_alma_datalink(mous_id: str, location: str = 'US') -> dict:
    """
    Query ALMA datalink for a MOUS ID using direct HTTP requests
    -----------
    mous_id : str
        MOUS ID with 'uid___' prefix
    location : str
        Server location ('US' or 'EU')
        
    Returns:
    --------
    dict
        Dictionary containing query results
    """
    result = {
        'MOUS_ID': mous_id,  # Keep original MOUS ID with prefix
        'MOUS_ID_Query': mous_id,
        'Data_Size_GB': None,
        'Download_Links': None,
        'Status': 'Failed',
        'Error': None,
        'Server': location
    }
    
    try:
        # ALMA datalink endpoints
        endpoints = {
            'US': 'https://almascience.nrao.edu/datalink/sync',
            'EU': 'https://almascience.eso.org/datalink/sync'
        }
        
        if location not in endpoints:
            result['Error'] = f"Invalid location: {location}"
            return result
        
        # Construct query URL
        query_url = f"{endpoints[location]}?ID={mous_id}"
        print(f"  Querying: {query_url}")
        
        # Make HTTP request
        with urllib.request.urlopen(query_url, timeout=30) as response:
            if response.status == 200:
                # Parse JSON response
                data = json.loads(response.read().decode('utf-8'))
                
                if 'links' in data and data['links']:
                    # Extract download links and calculate size
                    download_links = []
                    total_size_bytes = 0
                    
                    for link in data['links']:
                        if 'access_url' in link:
                            download_links.append(link['access_url'])
                        
                        if 'content_length' in link and link['content_length']:
                            try:
                                total_size_bytes += int(link['content_length'])
                            except (ValueError, TypeError):
                                pass
                    
                    if download_links:
                        result['Download_Links'] = '; '.join(download_links)
                        result['Status'] = 'Success'
                    
                    if total_size_bytes > 0:
                        total_size_gb = total_size_bytes / (1024**3)
                        result['Data_Size_GB'] = round(total_size_gb, 2)
                        print(f"  ✓ Found data: {result['Data_Size_GB']} GB")
                    else:
                        print(f"  ✓ Found data but no size information")
                else:
                    result['Status'] = 'No Data Found'
                    print(f"  ✗ No data found for {mous_id}")
            else:
                result['Error'] = f"HTTP {response.status}"
                result['Status'] = 'HTTP Error'
                print(f"  ✗ HTTP Error {response.status}")
                
    except urllib.error.URLError as e:
        result['Error'] = f"URL Error: {e}"
        result['Status'] = 'Connection Error'
        print(f"  ✗ Connection Error: {e}")
    except json.JSONDecodeError as e:
        result['Error'] = f"JSON Error: {e}"
        result['Status'] = 'Parse Error'
        print(f"  ✗ Parse Error: {e}")
    except Exception as e:
        result['Error'] = str(e)
        result['Status'] = 'Error'
        print(f"  ✗ Error querying {mous_id}: {e}")
    
    return result

def query_mous_data(mous_id: str, index: int) -> dict:
    """
    Query data size and download links for a single MOUS ID
    -----------
    mous_id : str
        MOUS ID with 'uid___' prefix
    index : int
        Index for server selection
        
    Returns:
    --------
    dict
        Dictionary containing query results
    """
    # Try US server first, then EU if needed
    servers = ['US', 'EU']
    server = servers[index % 2]
    
    result = query_alma_datalink(mous_id, server)
    
    # If US server fails, try EU server
    if result['Status'] not in ['Success', 'No Data Found'] and server == 'US':
        print(f"  Retrying with EU server...")
        sleep(1)
        result = query_alma_datalink(mous_id, 'EU')
    
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
    total_data_size = 0
    
    for i, mous_id in enumerate(missing_mous):
        print(f"\nProgress: {i+1}/{len(missing_mous)}")
        
        # Query data for this MOUS
        result = query_mous_data(mous_id, i)
        results.append(result)
        
        if result['Status'] == 'Success':
            successful_queries += 1
            if result['Data_Size_GB']:
                total_data_size += result['Data_Size_GB']
        
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
        print(f"Total data size found: {total_data_size:.2f} GB")
        
        # Show first few results
        print(f"\nFirst 5 results:")
        for i, result in enumerate(results[:5]):
            status_icon = "✓" if result['Status'] == 'Success' else "✗"
            size_str = f"{result['Data_Size_GB']} GB" if result['Data_Size_GB'] else "N/A"
            print(f"  {i+1}. {result['MOUS_ID']} - {status_icon} {result['Status']} - {size_str}")
        
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
