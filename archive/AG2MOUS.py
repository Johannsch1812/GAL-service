import json
import os
import sys
import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS

from shapely.geometry import Point, Polygon, shape
from shapely.prepared import prep
from rtree import index

from typing import List, Dict, Optional, Tuple

def read_catalog_dat(dat_file: str) -> Optional[pd.DataFrame]:
    """
    Read the ALMAGAL catalog from table3.dat file
    -----------
    dat_file : str
        Path to the table3.dat file
        
    Returns:
    --------
    pandas.DataFrame
        Catalog data with 6348 sources
    """
    # read_fwf uses 0-base indexing. For example, 1-17 indexes becomes 0-16 for read_fwf.
    colspecs = [
        (0, 17),   # Name, it reads bytes 0-16. Byte 17 is not taken.
        (18, 20),  # Seq
        (21, 30),  # RAdeg
        (31, 40),  # DEdeg
        (41, 50),  # GLON
        (51, 59),  # GLAT
        (60, 67),  # Fpeak
        (68, 74),  # Fbkg
        (75, 79),  # Diam1
        (80, 84),  # Diam2
        (85, 89),  # PAs
        (90, 97),  # Fint
        (98, 102), # MajDiam
        (103, 107),# MinDiam
        (108, 111),# PAb
        (112, 113),# Flag1
        (114, 115) # Flag2
    ]
    
    column_names = ['Name', 'Seq', 'RAdeg', 'DEdeg', 'GLON', 'GLAT', 
                   'Fpeak', 'Fbkg', 'Diam1', 'Diam2', 'PAs', 'Fint',
                   'MajDiam', 'MinDiam', 'PAb', 'Flag1', 'Flag2']
    
    try:
        df = pd.read_fwf(dat_file, colspecs=colspecs, names=column_names)
        print(f"Successfully read {len(df)} sources from catalog")
        return df
    except Exception as e:
        print(f"Error reading catalog file: {e}")
        return None

def extract_field_coverage_from_fits(fits_file: str) -> Optional[Dict]:
    """
    Extract field coverage (WCS + image dimensions) from FITS header
    -----------
    fits_file : str
        Path to FITS file
        
    Returns:
    --------
    dict or None
        Dictionary with field coverage information
    """
    try:
        with fits.open(fits_file) as hdul:
            header = hdul[0].header
            
            naxis1 = header.get('NAXIS1', 0)
            naxis2 = header.get('NAXIS2', 0)
            
            if naxis1 == 0 or naxis2 == 0:
                print(f"Warning: Could not determine image dimensions for {fits_file}")
                return None
            
            # Create WCS object
            try:
                wcs = WCS(header).celestial
            except Exception as e:
                print(f"Warning: Could not create WCS for {fits_file}: {e}")
                return None
            
            # Get corner coordinates of the image
            corners_pix = [
                [0, 0],                    # Bottom-left
                [naxis1-1, 0],             # Bottom-right  
                [naxis1-1, naxis2-1],      # Top-right
                [0, naxis2-1]              # Top-left
            ]
            
            # Convert pixel coordinates to world coordinates
            corners_world = wcs.pixel_to_world_values(corners_pix)
            
            # Get center coordinates
            center_pix = [(naxis1-1)/2, (naxis2-1)/2]
            center_world = wcs.pixel_to_world_values([center_pix])
            
            return {
                'center_ra': center_world[0][0],
                'center_dec': center_world[0][1],
                'corners_ra': [c[0] for c in corners_world],
                'corners_dec': [c[1] for c in corners_world],
                'naxis1': naxis1,
                'naxis2': naxis2,
                'wcs': wcs,
                'polygon': None  # A place holder. It will be set later
            }
            
    except Exception as e:
        print(f"Error extracting coverage from {fits_file}: {e}")
        return None

def create_coverage_polygon(coverage: Dict) -> Optional[Polygon]:
    """
    Create a Shapely polygon from field coverage
    -----------
    coverage : dict
        Coverage information from extract_field_coverage_from_fits
        
    Returns:
    --------
    shapely.geometry.Polygon or None
        Polygon representing the field coverage
    """
    try:
        # Handle RA wrap-around at 0/360 degrees
        corners_ra = np.array(coverage['corners_ra'])
        corners_dec = np.array(coverage['corners_dec'])
        
        # Check for RA wrap-around (difference > 180 degrees)
        if np.ptp(corners_ra) > 180:
            # Adjust RA coordinates for wrap-around
            corners_ra = np.where(corners_ra > 180, corners_ra - 360, corners_ra)
        
        # Create polygon from corners
        coords = list(zip(corners_ra, corners_dec))
        polygon = Polygon(coords)
        
        # Validate polygon
        if not polygon.is_valid:
            print(f"Warning: Invalid polygon created for field coverage")
            return None
        
        return polygon
        
    except Exception as e:
        print(f"Error creating coverage polygon: {e}")
        return None

def get_target_coverages_from_json(fits_check_json: str, config: str = '7M') -> Optional[pd.DataFrame]:
    """
    Extract target field coverages from the FITS checking JSON results
    -----------
    fits_check_json : str
        Path to the JSON file from fits checking script
    config : str
        Configuration to use ('7M', 'TM1', or 'TM2')
        
    Returns:
    --------
    pandas.DataFrame or None
        DataFrame with target field coverages, or None if failed
    """
    try:
        with open(fits_check_json, 'r') as f:
            data = json.load(f)
        
        print(f"Using {config} configuration for field coverage extraction...")
        
        target_coverages = []
        extracted_count = 0
        failed_count = 0
        
        for entry in data:
            target_name = entry['Source']
            
            # Check if the specified configuration has a FITS file
            if entry[config].get('cont_check') == 1:
                fits_file = entry[config].get('cont_pbcor_fits')
                if fits_file and os.path.exists(fits_file):
                    print(f"Processing {target_name}...")
                    coverage = extract_field_coverage_from_fits(fits_file)
                    if coverage:
                        polygon = create_coverage_polygon(coverage)
                        if polygon:
                            coverage['polygon'] = polygon
                            coverage['prepared_polygon'] = prep(polygon)  # For faster point-in-polygon tests
                            
                            target_coverages.append({
                                'Target': target_name,
                                'RA_center': coverage['center_ra'],
                                'Dec_center': coverage['center_dec'],
                                'FITS_file': fits_file,
                                'coverage': coverage
                            })
                            extracted_count += 1
                        else:
                            failed_count += 1
                    else:
                        failed_count += 1
                else:
                    if fits_file:
                        print(f"Warning: FITS file not found: {fits_file}")
                    failed_count += 1
        
        print(f"Successfully extracted coverage from {extracted_count} {config} FITS files")
        print(f"Failed to extract coverage from {failed_count} files")
        
        return pd.DataFrame(target_coverages)
    
    except Exception as e:
        print(f"Error reading FITS check results: {e}")
        return None

def create_spatial_coverage_index(targets_df: pd.DataFrame) -> Tuple[index.Index, Dict]:
    """
    Create an R-tree spatial index for target field coverages
    -----------
    targets_df : pandas.DataFrame
        DataFrame with target field coverages
        
    Returns:
    --------
    tuple
        (R-tree index, dictionary mapping index_id to target info)
    """
    print("Creating spatial coverage index (R-tree) for fast polygon searches...")
    
    # Create R-tree index
    idx = index.Index()
    target_lookup = {}
    
    for i, target_info in targets_df.iterrows():
        coverage = target_info['coverage']
        polygon: Polygon = coverage['polygon']
        
        # Get bounding box of the polygon
        minx, miny, maxx, maxy = polygon.bounds
        
        # Handle RA wrap-around in bounding box
        corners_ra = np.array(coverage['corners_ra'])
        if np.ptp(corners_ra) > 180:  # Field crosses RA=0
            # Create two bounding boxes if necessary
            if minx < 0:  # Already adjusted coordinates
                # Insert the polygon with adjusted bounds
                idx.insert(i, (minx, miny, maxx, maxy))
            else:
                # Insert with original bounds but mark as wrapped
                idx.insert(i, (minx, miny, maxx, maxy))
        else:
            # Normal case - no RA wrap
            idx.insert(i, (minx, miny, maxx, maxy))
        
        # Store target information for lookup
        target_lookup[i] = {
            'target_name': target_info['Target'],
            'ra_center': target_info['RA_center'],
            'dec_center': target_info['Dec_center'],
            'fits_file': target_info['FITS_file'],
            'coverage': coverage,
            'polygon': polygon,
            'prepared_polygon': coverage['prepared_polygon']
        }
    
    print(f"Spatial index created for {len(target_lookup)} target coverages")
    return idx, target_lookup

def map_sources_to_coverage_efficient(catalog_df: pd.DataFrame, targets_df: pd.DataFrame) -> List[Dict]:
    """
    Map catalog sources to targets using spatial indexing of coverage polygons
    -----------
    catalog_df : pandas.DataFrame
        Catalog with 6348 sources
    targets_df : pandas.DataFrame
        Targets with field coverage information
        
    Returns:
    --------
    list
        List of mapping dictionaries
    """
    print("Starting coverage-based mapping with spatial indexing...")
    
    # Create spatial index for target coverages
    spatial_idx, target_lookup = create_spatial_coverage_index(targets_df)
    
    mapping_list = []
    
    # Process sources in batches for better memory management
    batch_size = 1000
    total_sources = len(catalog_df)
    
    for batch_start in range(0, total_sources, batch_size):
        batch_end = min(batch_start + batch_size, total_sources)
        print(f"Processing sources {batch_start+1} to {batch_end}...")
        
        batch_sources = catalog_df.iloc[batch_start:batch_end]
        
        for i, (_, source_row) in enumerate(batch_sources.iterrows()):
            source_ra, source_dec = source_row['RAdeg'], source_row['DEdeg']
            
            # Create source mapping entry
            source_entry = {
                'Source_Name': source_row['Name'].strip(),
                'Source_Seq': int(source_row['Seq']),
                'Source_RA': source_ra,
                'Source_Dec': source_dec,
                'Source_GLON': source_row['GLON'],
                'Source_GLAT': source_row['GLAT'],
                'Source_Fpeak': source_row['Fpeak'],
                'Source_Fint': source_row['Fint'],
            }
            
            matching_targets = []
            
            # Handle RA wrap-around for source coordinate
            test_ra = source_ra
            
            # Query spatial index for potential matches
            # Use a small buffer around the source point for the bounding box query
            buffer = 0.001  # Small buffer in degrees
            potential_matches = list(spatial_idx.intersection((test_ra - buffer, source_dec - buffer, 
                                                             test_ra + buffer, source_dec + buffer)))
            
            # Also check for RA wrap-around case
            if source_ra > 180:
                test_ra_wrapped = source_ra - 360
                potential_matches_wrapped = list(spatial_idx.intersection((test_ra_wrapped - buffer, source_dec - buffer,
                                                                         test_ra_wrapped + buffer, source_dec + buffer)))
                potential_matches.extend(potential_matches_wrapped)
            elif source_ra < 180:
                test_ra_wrapped = source_ra + 360
                potential_matches_wrapped = list(spatial_idx.intersection((test_ra_wrapped - buffer, source_dec - buffer,
                                                                         test_ra_wrapped + buffer, source_dec + buffer)))
                potential_matches.extend(potential_matches_wrapped)
            
            # Remove duplicates
            potential_matches = list(set(potential_matches))
            
            # Test actual polygon containment for potential matches
            for target_idx in potential_matches:
                if target_idx in target_lookup:
                    target_info = target_lookup[target_idx]
                    
                    # Determine which RA coordinate to test
                    coverage = target_info['coverage']
                    corners_ra = np.array(coverage['corners_ra'])
                    
                    test_ra_final = source_ra
                    if np.ptp(corners_ra) > 180:  # Field crosses RA=0
                        if source_ra > 180:
                            test_ra_final = source_ra - 360
                    
                    # Test if source is within field coverage
                    source_point = Point(test_ra_final, source_dec)
                    
                    try:
                        if target_info['prepared_polygon'].contains(source_point):
                            # Calculate distance to field center for ranking
                            center_coord = SkyCoord(ra=target_info['ra_center']*u.deg, 
                                                  dec=target_info['dec_center']*u.deg, 
                                                  frame='icrs')
                            source_coord = SkyCoord(ra=source_ra*u.deg, 
                                                  dec=source_dec*u.deg, 
                                                  frame='icrs')
                            distance = center_coord.separation(source_coord)
                            
                            target_match = {
                                'Target_Name': target_info['target_name'],
                                'Target_RA_center': target_info['ra_center'],
                                'Target_Dec_center': target_info['dec_center'],
                                'Distance_to_center_arcsec': distance.to(u.arcsec).value,
                                'FITS_file': target_info['fits_file']
                            }
                            matching_targets.append(target_match)
                    except Exception as e:
                        print(f"Warning: Error testing point in polygon for target {target_info['target_name']}: {e}")
                        continue
            
            # Sort matching targets by distance to center
            if matching_targets:
                matching_targets.sort(key=lambda x: x['Distance_to_center_arcsec'])
                
                # Add distance ranks
                for k, target in enumerate(matching_targets):
                    target['Distance_rank'] = k + 1
                
                source_entry['Matching_Targets'] = matching_targets
                source_entry['Primary_Target'] = matching_targets[0]['Target_Name']
                source_entry['Primary_Distance_arcsec'] = matching_targets[0]['Distance_to_center_arcsec']
                source_entry['Primary_FITS'] = matching_targets[0]['FITS_file']
                source_entry['Num_Matching_Targets'] = len(matching_targets)
                source_entry['Spatial_Candidates'] = len(potential_matches)  # For debugging
            else:
                # No matching targets
                source_entry['Matching_Targets'] = []
                source_entry['Primary_Target'] = None
                source_entry['Primary_Distance_arcsec'] = None
                source_entry['Primary_FITS'] = None
                source_entry['Num_Matching_Targets'] = 0
                source_entry['Spatial_Candidates'] = len(potential_matches)  # For debugging
            
            mapping_list.append(source_entry)
    
    return mapping_list

def print_mapping_statistics(mapping_list: List[Dict]) -> None:
    """
    Print statistics about the source-target mapping
    -----------
    mapping_list : list
        List of mapping dictionaries
    """
    total_sources = len(mapping_list)
    mapped_sources = sum(1 for entry in mapping_list if entry['Primary_Target'] is not None)
    unmapped_sources = total_sources - mapped_sources
    
    # Count sources with multiple targets
    multiple_targets = sum(1 for entry in mapping_list if entry['Num_Matching_Targets'] > 1)
    
    # Distance statistics for mapped sources
    distances = [entry['Primary_Distance_arcsec'] for entry in mapping_list 
                if entry['Primary_Distance_arcsec'] is not None]
    
    print("\nSource-Target Coverage Mapping Statistics:")
    print(f"Total sources: {total_sources}")
    print(f"Mapped sources: {mapped_sources} ({mapped_sources/total_sources*100:.1f}%)")
    print(f"Unmapped sources: {unmapped_sources} ({unmapped_sources/total_sources*100:.1f}%)")
    print(f"Sources with multiple target coverage: {multiple_targets} ({multiple_targets/total_sources*100:.1f}%)")
    
    if distances:
        print(f"\nDistance to Center Statistics (for mapped sources):")
        print(f"Mean distance: {np.mean(distances):.1f} arcsec")
        print(f"Median distance: {np.median(distances):.1f} arcsec")
        print(f"Max distance: {np.max(distances):.1f} arcsec")
        print(f"Min distance: {np.min(distances):.1f} arcsec")

def save_mapping_results(mapping_list: List[Dict], output_json: str, output_csv: str) -> None:
    """
    Save mapping results to JSON and CSV files
    -----------
    mapping_list : list
        List of mapping dictionaries
    output_json : str
        Output JSON file path
    output_csv : str
        Output CSV file path
    """
    # Save full mapping to JSON (exclude non-serializable objects)
    try:
        # Create a clean version for JSON (remove WCS and polygon objects)
        json_mapping = []
        for entry in mapping_list:
            clean_entry = entry.copy()
            # Remove any non-serializable objects if they exist
            json_mapping.append(clean_entry)
        
        with open(output_json, 'w', encoding='utf-8') as f:
            json.dump(json_mapping, f, indent=2, ensure_ascii=False)
        print(f"Full mapping saved to {output_json}")
    except Exception as e:
        print(f"Error saving JSON: {e}")
    
    # Create simplified CSV
    try:
        csv_data = []
        for entry in mapping_list:
            csv_row = {
                'Source_Name': entry['Source_Name'],
                'Source_Seq': entry['Source_Seq'],
                'Source_RA': entry['Source_RA'],
                'Source_Dec': entry['Source_Dec'],
                'Source_GLON': entry['Source_GLON'],
                'Source_GLAT': entry['Source_GLAT'],
                'Source_Fpeak': entry['Source_Fpeak'],
                'Source_Fint': entry['Source_Fint'],
                'Primary_Target': entry['Primary_Target'],
                'Primary_Distance_arcsec': entry['Primary_Distance_arcsec'],
                'Num_Matching_Targets': entry['Num_Matching_Targets'],
                'All_Targets': ', '.join([t['Target_Name'] for t in entry['Matching_Targets']]) if entry['Matching_Targets'] else '',
                'Primary_FITS': os.path.basename(entry['Primary_FITS']) if entry['Primary_FITS'] else None
            }
            csv_data.append(csv_row)
        
        csv_df = pd.DataFrame(csv_data)
        csv_df.to_csv(output_csv, index=False)
        print(f"Summary saved to {output_csv}")
    except Exception as e:
        print(f"Error saving CSV: {e}")

if __name__ == "__main__":
    """
    Main function
    """
    # Configuration
    catalog_file = "tables/table3.dat"  # Path to catalog file
    fits_check_json = "almagal_source_mous_checking.json"  # From previous script
    config = "7M"  # Configuration to use: '7M', 'TM1', or 'TM2'
    output_json = "source_target_coverage_mapping.json"
    output_csv = "source_target_coverage_mapping.csv"
    
    print("ALMAGAL Source to Target Coverage Mapping Script")
    print("=" * 55)
    print(f"Using {config} configuration for field coverage extraction")
    
    # Step 1: Read catalog
    print("Step 1: Reading catalog...")
    catalog_df = read_catalog_dat(catalog_file)
    if catalog_df is None:
        print("Failed to read catalog. Exiting.")
        sys.exit(1)
    
    # Step 2: Get target field coverages from JSON
    print("Step 2: Extracting field coverages from FITS files...")
    if not os.path.exists(fits_check_json):
        print(f"FITS checking JSON file not found: {fits_check_json}")
        sys.exit(1)
    
    targets_df = get_target_coverages_from_json(fits_check_json, config)
    if targets_df is None or len(targets_df) == 0:
        print("Could not obtain target field coverages. Exiting.")
        sys.exit(1)
    
    print(f"Extracted field coverage for {len(targets_df)} targets")
    
    # Step 3: Perform efficient coverage-based mapping
    print("Step 3: Mapping sources to target field coverage using spatial indexing...")
    mapping_list = map_sources_to_coverage_efficient(catalog_df, targets_df)
    
    # Step 4: Print statistics
    print_mapping_statistics(mapping_list)
    
    # Step 5: Save results
    print("\nStep 4: Saving results...")
    save_mapping_results(mapping_list, output_json, output_csv)
    
    print("\nCoverage mapping completed successfully!")

