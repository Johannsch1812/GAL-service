import pandas as pd
import os
from typing import List, Set

def extract_missing_mous(cont_csv: str, cube_csv: str) -> None:
    """
    Extract MOUS IDs for missing data from cont and cube checking CSV files
    -----------
    cont_csv : str
        Path to cont checking CSV file
    cube_csv : str
        Path to cube checking CSV file
    """
    
    def read_csv_and_extract_missing(csv_file: str, data_type: str) -> Set[str]:
        """
        Read CSV file and extract missing MOUS IDs
        -----------
        csv_file : str
            Path to CSV file
        data_type : str
            Type of data (cont or cube)
            
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
                
            df = pd.read_csv(csv_file)
            print(f"\nReading {data_type} data from {csv_file}")
            print(f"Total rows: {len(df)}")
            
            # Check each configuration
            configurations = ['7M', 'TM1', 'TM2']
            
            for config in configurations:
                check_col = f'{config}_check'
                mous_col = f'MOUS_{config}'
                
                if check_col in df.columns and mous_col in df.columns:
                    # Find rows where check is 0 (missing data)
                    missing_mask = df[check_col] == 0
                    missing_mous_config = df[missing_mask][mous_col].dropna().tolist()
                    
                    print(f"  {config}: {len(missing_mous_config)} missing MOUS IDs")
                    missing_mous.update(missing_mous_config)
                else:
                    print(f"  Warning: Columns {check_col} or {mous_col} not found in {csv_file}")
            
            return missing_mous
            
        except Exception as e:
            print(f"Error reading {csv_file}: {e}")
            return missing_mous
    
    # Extract missing MOUS from both files
    print("="*60)
    print("EXTRACTING MISSING MOUS IDs")
    print("="*60)
    
    cont_missing = read_csv_and_extract_missing(cont_csv, "cont")
    cube_missing = read_csv_and_extract_missing(cube_csv, "cube")
    
    # Combine and get unique missing MOUS
    all_missing = cont_missing.union(cube_missing)
    
    print(f"\nSummary:")
    print(f"  Cont missing MOUS: {len(cont_missing)}")
    print(f"  Cube missing MOUS: {len(cube_missing)}")
    print(f"  Total unique missing MOUS: {len(all_missing)}")
    
    # Save results
    if all_missing:
        # Save to text file
        output_file = "missing_mous_list.txt"
        with open(output_file, 'w') as f:
            f.write("Missing MOUS IDs\n")
            f.write("="*50 + "\n\n")
            f.write("Cont missing MOUS:\n")
            for mous in sorted(cont_missing):
                f.write(f"  {mous}\n")
            f.write(f"\nCube missing MOUS:\n")
            for mous in sorted(cube_missing):
                f.write(f"  {mous}\n")
            f.write(f"\nAll unique missing MOUS:\n")
            for mous in sorted(all_missing):
                f.write(f"  {mous}\n")
        
        print(f"\nMissing MOUS IDs saved to {output_file}")
        
        # Save to CSV file
        csv_output_file = "missing_mous_list.csv"
        missing_df = pd.DataFrame({
            'MOUS_ID': sorted(all_missing),
            'Missing_in_Cont': [mous in cont_missing for mous in sorted(all_missing)],
            'Missing_in_Cube': [mous in cube_missing for mous in sorted(all_missing)]
        })
        missing_df.to_csv(csv_output_file, index=False)
        print(f"Missing MOUS IDs also saved to {csv_output_file}")
        
        # Display first 10 missing MOUS
        print(f"\nFirst 10 missing MOUS IDs:")
        for i, mous in enumerate(sorted(all_missing)[:10]):
            cont_status = "✓" if mous in cont_missing else "✗"
            cube_status = "✓" if mous in cube_missing else "✗"
            print(f"  {i+1:2d}. {mous} (Cont:{cont_status} Cube:{cube_status})")
        
        if len(all_missing) > 10:
            print(f"  ... and {len(all_missing) - 10} more")
    else:
        print("\nNo missing MOUS IDs found!")

if __name__ == "__main__":
    # File paths
    cont_csv = "almagal_QA2_cont_checking.csv"
    cube_csv = "almagal_QA2_cube_checking.csv"
    
    # Check if files exist
    if not os.path.exists(cont_csv):
        print(f"Error: {cont_csv} not found")
        exit(1)
    
    if not os.path.exists(cube_csv):
        print(f"Error: {cube_csv} not found")
        exit(1)
    
    # Extract missing MOUS
    extract_missing_mous(cont_csv, cube_csv)
    
    print("\nExtraction completed!")
