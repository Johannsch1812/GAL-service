#!/usr/bin/env python3
"""
Compare missing MOUS lists between FITS data and MS data
"""

import csv
import pandas as pd

def read_missing_mous_fits(csv_file):
    """Read missing MOUS IDs from FITS data CSV"""
    mous_ids = set()
    try:
        with open(csv_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                mous_id = row['MOUS_ID'].strip()
                mous_ids.add(mous_id)
        return mous_ids
    except Exception as e:
        print(f"Error reading {csv_file}: {e}")
        return set()

def read_missing_mous_ms(csv_file):
    """Read missing MOUS IDs from MS data CSV"""
    mous_ids = set()
    try:
        with open(csv_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                mous_id = row['MOUS_ID'].strip()
                mous_ids.add(mous_id)
        return mous_ids
    except Exception as e:
        print(f"Error reading {csv_file}: {e}")
        return set()

def compare_missing_mous():
    """Compare missing MOUS lists between FITS and MS data"""
    
    # Read the missing MOUS lists
    fits_missing = read_missing_mous_fits('missing_mous_list.csv')
    ms_missing = read_missing_mous_ms('missing_mous_ms_list.csv')
    
    print("="*80)
    print("MISSING MOUS COMPARISON ANALYSIS")
    print("="*80)
    
    print(f"\nFITS Data Missing MOUS IDs: {len(fits_missing)}")
    print(f"MS Data Missing MOUS IDs: {len(ms_missing)}")
    
    # Find common missing MOUS IDs
    common_missing = fits_missing.intersection(ms_missing)
    print(f"Common Missing MOUS IDs: {len(common_missing)}")
    
    # Find FITS-only missing MOUS IDs
    fits_only = fits_missing - ms_missing
    print(f"FITS-only Missing MOUS IDs: {len(fits_only)}")
    
    # Find MS-only missing MOUS IDs
    ms_only = ms_missing - fits_missing
    print(f"MS-only Missing MOUS IDs: {len(ms_only)}")
    
    # Calculate percentages
    if len(fits_missing) > 0:
        common_percentage = (len(common_missing) / len(fits_missing)) * 100
        print(f"Common missing as % of FITS missing: {common_percentage:.1f}%")
    
    if len(ms_missing) > 0:
        common_percentage_ms = (len(common_missing) / len(ms_missing)) * 100
        print(f"Common missing as % of MS missing: {common_percentage_ms:.1f}%")
    
    # Show detailed breakdown
    print(f"\n{'='*80}")
    print("DETAILED BREAKDOWN")
    print(f"{'='*80}")
    
    if common_missing:
        print(f"\nCommon Missing MOUS IDs ({len(common_missing)}):")
        for mous_id in sorted(common_missing):
            print(f"  {mous_id}")
    
    if fits_only:
        print(f"\nFITS-only Missing MOUS IDs ({len(fits_only)}):")
        for mous_id in sorted(fits_only):
            print(f"  {mous_id}")
    
    if ms_only:
        print(f"\nMS-only Missing MOUS IDs ({len(ms_only)}):")
        for mous_id in sorted(ms_only):
            print(f"  {mous_id}")
    
    # Analysis summary
    print(f"\n{'='*80}")
    print("ANALYSIS SUMMARY")
    print(f"{'='*80}")
    
    if len(common_missing) > 0:
        print(f"✓ {len(common_missing)} MOUS IDs are missing from BOTH FITS and MS data")
        print("  This suggests systematic data availability issues for these MOUS IDs")
    
    if len(fits_only) > 0:
        print(f"⚠ {len(fits_only)} MOUS IDs are missing from FITS data but present in MS data")
        print("  This suggests FITS data processing issues or different data availability")
    
    if len(ms_only) > 0:
        print(f"⚠ {len(ms_only)} MOUS IDs are missing from MS data but present in FITS data")
        print("  This suggests MS data processing issues or different data availability")
    
    # Consistency check
    total_unique = len(fits_missing.union(ms_missing))
    if total_unique > 0:
        consistency = (len(common_missing) / total_unique) * 100
        print(f"\nConsistency Score: {consistency:.1f}%")
        if consistency > 80:
            print("✓ High consistency between FITS and MS missing data")
        elif consistency > 50:
            print("⚠ Moderate consistency between FITS and MS missing data")
        else:
            print("✗ Low consistency between FITS and MS missing data")
    
    return {
        'fits_missing': fits_missing,
        'ms_missing': ms_missing,
        'common_missing': common_missing,
        'fits_only': fits_only,
        'ms_only': ms_only
    }

if __name__ == "__main__":
    results = compare_missing_mous()
