import argparse
import os
import sys
from typing import List, Optional

import pandas as pd


def read_database(database_file: str) -> Optional[pd.DataFrame]:
    """
    Read the ALMAGAL database.xlsx with target and MOUS mapping.
    """
    try:
        import warnings

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning, module="openpyxl")
            df = pd.read_excel(database_file, engine="openpyxl")

        print(f"Successfully read database with {len(df)} targets")
        print(f"Database columns: {df.columns.tolist()}")
        return df
    except Exception as e:
        print(f"Error reading database file: {e}")
        return None


def load_scan_csv(csv_files: List[str]) -> pd.DataFrame:
    """
    Load one or multiple CSV files produced by scan_fits.py and concatenate them.
    """
    frames: List[pd.DataFrame] = []
    for path in csv_files:
        df = pd.read_csv(path)
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    combined = pd.concat(frames, ignore_index=True)
    # Normalize helpful columns
    if "filename" in combined.columns:
        combined["filename_lower"] = combined["filename"].astype(str).str.lower()
    else:
        combined["filename_lower"] = combined.get("path", pd.Series(dtype=str)).astype(str).str.lower()
    if "image_type" not in combined.columns:
        combined["image_type"] = None
    return combined


def check_fits_contains_terms(filename: str, source_name: str, mous_id: str) -> bool:
    if pd.isna(source_name) or pd.isna(mous_id):
        return False
    source = str(source_name).strip()
    source_convent = source.replace('+', 'p')
    mous = str(mous_id).strip()
    name = str(filename)
    return ((source in name) or (source_convent in name)) and (mous in name)


def find_matching_rows(
    source_name: str,
    mous_id: str,
    scan_df: pd.DataFrame,
    allowed_types: List[str],
    require_pbcor: bool,
) -> List[pd.Series]:
    matches: List[pd.Series] = []
    allowed_set = set(t.lower() for t in allowed_types)
    for _, row in scan_df.iterrows():
        file_ok = True
        if allowed_set:
            file_type = str(row.get("image_type", "")).lower()
            if file_type not in allowed_set:
                file_ok = False
        if file_ok and require_pbcor:
            filename_lower = str(row.get("filename_lower", "")).lower()
            if "pbcor" not in filename_lower:
                file_ok = False
        if not file_ok:
            continue
        filename_to_check = row.get("filename") or row.get("path")
        if filename_to_check is None:
            continue
        if check_fits_contains_terms(str(filename_to_check), source_name, mous_id):
            matches.append(row)
    return matches


def create_source_mapping_from_csv(
    database_df: pd.DataFrame,
    scan_df: pd.DataFrame,
    include_cubes: bool,
    require_pbcor: bool,
) -> List[dict]:
    check_list = []
    allowed_types = ["cont"] if not include_cubes else ["cont", "cube"]
    for _, row in database_df.iterrows():
        source_name = row['Source']
        source_dict = {
            'Source': str(source_name),
            'SGOUS': row.get('SGOUS', None),
            'GOUS': row.get('GOUS', None),
        }
        for config in ['7M', 'TM1', 'TM2']:
            mous_col = f'MOUS_{config}'
            mous_id = row.get(mous_col, None)
            config_dict = {
                'CONFIG': str(config),
                'MOUS': mous_id,
            }
            if not mous_id:
                source_dict[config] = config_dict
                continue
            matches = find_matching_rows(str(source_name), mous_id, scan_df, allowed_types, require_pbcor)
            if matches:
                best = matches[0]
                config_dict['cont_pbcor_fits'] = best.get('path') or best.get('filename')
                config_dict['cont_check'] = 1
                if len(matches) > 1:
                    config_dict[f'{config}_all_matches'] = [m.get('path') or m.get('filename') for m in matches]
            else:
                config_dict['cont_pbcor_fits'] = None
                config_dict['cont_check'] = 0
            source_dict[config] = config_dict
        check_list.append(source_dict)
    return check_list


def print_statistics(check_list: List[dict]) -> None:
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
    all_configs = sum(1 for s in check_list if all(s[c].get('cont_check', 0) == 1 for c in ['7M', 'TM1', 'TM2']))
    if total_sources > 0:
        print(f"\nSources with all 3 configurations: {all_configs}/{total_sources} ({all_configs/total_sources*100:.1f}%)")


def save_to_json(check_list: List[dict], output_file: str) -> None:
    import json
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(check_list, f, indent=2, ensure_ascii=False)
    print(f"\nMapping saved to {output_file}")
    print(f"Total entries: {len(check_list)}")


def save_summary_csv(check_list: List[dict], output_file: str) -> None:
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


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Read server scan CSV(s) and perform ALMAGAL source-MOUS-FITS checks without directory scanning."
    )
    parser.add_argument("--database", default="tables/database.xlsx", help="Path to database.xlsx")
    parser.add_argument("--scan-csv", nargs="+", required=True, help="One or more CSV files generated by scan_fits.py")
    parser.add_argument("--output-json", default="almagal_source_mous_checking.json", help="Output JSON path")
    parser.add_argument("--output-csv", default="almagal_source_mous_checking.csv", help="Output CSV path")
    parser.add_argument("--include-cubes", action="store_true", help="Allow matching 3D cubes in addition to 2D cont images")
    parser.add_argument("--allow-non-pbcor", action="store_true", help="Do not require 'pbcor' in filename")
    return parser.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    print("Reading database.xlsx")
    database_df = read_database(args.database)
    if database_df is None:
        print("Failed to read database.xlsx Exiting.")
        return 1
    print("\nLoading scan CSV(s)...")
    scan_df = load_scan_csv(args.scan_csv)
    if scan_df.empty:
        print("No rows found in scan CSV(s). Exiting.")
        return 1
    print(f"Loaded {len(scan_df)} scanned files")
    print("\nStart to check source-MOUS-FITS from CSV")
    check_list = create_source_mapping_from_csv(
        database_df,
        scan_df,
        include_cubes=args.include_cubes,
        require_pbcor=(not args.allow_non_pbcor),
    )
    print_statistics(check_list)
    save_to_json(check_list, args.output_json)
    save_summary_csv(check_list, args.output_csv)
    print(f"\nResults have been stored in JSON file {args.output_json} and CSV file {args.output_csv}")
    print("\nChecking completed successfully!")
    return 0


if __name__ == "__main__":
    sys.exit(main())


