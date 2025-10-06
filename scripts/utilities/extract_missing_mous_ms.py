import csv
import os
from typing import List, Set

def extract_missing_mous_ms(ms_csv: str) -> None:
    """
    Extract MOUS IDs for missing MS data from MS checking CSV file
    -----------
    ms_csv : str
        Path to MS checking CSV file
    """
    
    def read_csv_and_extract_missing(csv_file: str) -> Set[str]:
        """
        Read CSV file and extract missing MOUS IDs
        -----------
        csv_file : str
            Path to CSV file
            
        Returns:
        --------
        set
            Set of missing MOUS IDs
        """
        missing_mous = set()
        
        try:
            if not os.path.exists(csv_file):
                print(f"Warning: {csv_file} not found")
                return missing_mous
                
            with open(csv_file, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                rows = list(reader)
                
            print(f"\nReading MS data from {csv_file}")
            print(f"Total rows: {len(rows)}")
            
            # Check each configuration
            configurations = ['7M', 'TM1', 'TM2']
            
            for config in configurations:
                check_col = f'{config}_ms_check'
                mous_col = f'MOUS_{config}'
                
                missing_count = 0
                for row in rows:
                    if check_col in row and mous_col in row:
                        # Check if the value is 0 (missing data)
                        if row[check_col] == '0' and row[mous_col].strip():
                            missing_mous.add(row[mous_col].strip())
                            missing_count += 1
                
                print(f"  {config}: {missing_count} missing MOUS IDs")
            
            return missing_mous
            
        except Exception as e:
            print(f"Error reading {csv_file}: {e}")
            return missing_mous
    
    # Extract missing MOUS from MS file
    print("="*60)
    print("EXTRACTING MISSING MOUS IDs FOR MS DATA")
    print("="*60)
    
    ms_missing = read_csv_and_extract_missing(ms_csv)
    
    print(f"\nSummary:")
    print(f"  MS missing MOUS: {len(ms_missing)}")
    
    # Save results
    if ms_missing:
        # Save to text file
        output_file = "missing_mous_ms_list.txt"
        with open(output_file, 'w') as f:
            f.write("Missing MOUS IDs for MS Data\n")
            f.write("="*50 + "\n\n")
            f.write("MS missing MOUS:\n")
            for mous in sorted(ms_missing):
                f.write(f"  {mous}\n")
        
        print(f"\nMissing MS MOUS IDs saved to {output_file}")
        
        # Save to CSV file
        csv_output_file = "missing_mous_ms_list.csv"
        with open(csv_output_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(['MOUS_ID', 'Data_Type'])
            for mous in sorted(ms_missing):
                writer.writerow([mous, 'MS'])
        
        print(f"Missing MS MOUS IDs also saved to {csv_output_file}")
        
        # Display first 10 missing MOUS
        print(f"\nFirst 10 missing MS MOUS IDs:")
        for i, mous in enumerate(sorted(ms_missing)[:10]):
            print(f"  {i+1:2d}. {mous}")
        
        if len(ms_missing) > 10:
            print(f"  ... and {len(ms_missing) - 10} more")
    else:
        print("\nNo missing MS MOUS IDs found!")

if __name__ == "__main__":
    # File path
    ms_csv = "almagal_source_mous_ms_checking.csv"
    
    # Check if file exists
    if not os.path.exists(ms_csv):
        print(f"Error: {ms_csv} not found")
        exit(1)
    
    # Extract missing MOUS
    extract_missing_mous_ms(ms_csv)
    
    print("\nExtraction completed!")
