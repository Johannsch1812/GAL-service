# ALMA Data Checking Project - File Structure

## Overview

This document describes the organized file structure for the ALMA data checking project. The project has been restructured to improve maintainability, clarity, and ease of navigation.

## Directory Structure

```
GAL-service/
├── scripts/                          # All Python scripts organized by function
│   ├── fits_checking/               # FITS data checking scripts
│   │   ├── fits_check_URL.py       # Main FITS file scanner
│   │   ├── fits_check.py           # Original FITS checker
│   │   ├── fits_check_from_csv.py  # CSV-based FITS checker
│   │   └── scan_fits.py            # FITS file scanner utility
│   ├── ms_checking/                # MS data checking scripts
│   │   └── MS_check_URL.py         # MS data scanner
│   ├── querying/                   # ALMA archive query scripts
│   │   ├── query_missing_mous.py   # Main query script with pyvo
│   │   ├── query_missing_mous_standalone.py  # Standalone query script
│   │   └── query.py                # Query utilities library
│   └── utilities/                  # Utility and extraction scripts
│       ├── extract_missing_mous.py      # Extract missing FITS MOUS
│       ├── extract_missing_mous_ms.py   # Extract missing MS MOUS
│       ├── compare_missing_mous.py      # Compare missing MOUS lists
│       └── download_*.sh               # Download scripts
├── data/                           # All data files
│   ├── input/                      # Input data files
│   │   └── tables/                 # Database and table files
│   │       ├── database.xlsx       # Main database
│   │       ├── table3.dat          # Table 3 data
│   │       └── table7.dat          # Table 7 data
│   ├── output/                     # Generated output files
│   │   ├── *.csv                   # All CSV output files
│   │   └── *.json                  # All JSON output files
│   └── intermediate/               # Intermediate processing files
├── results/                        # Analysis results and statistics
│   └── statistics/                 # Statistical output files
│       ├── statistics_cont.txt     # Continuum statistics
│       ├── statistics_cube.txt     # Cube statistics
│       ├── missing_mous_list.txt   # Missing MOUS list (text)
│       └── missing_mous_ms_list.txt # Missing MS MOUS list (text)
├── docs/                          # Documentation
│   ├── reports/                   # Analysis reports and summaries
│   │   ├── ALMA_Data_Checking_Summary.md
│   │   ├── Missing_MOUS_Consistency_Analysis.md
│   │   ├── ALAMGAL_Timing_Analysis.md
│   │   ├── ALMAGAL_Pipeline_Task.md
│   │   ├── ALMAGAL_QA2_Data_Task.md
│   │   ├── ALMAGAL_source_split.md
│   │   ├── conversion_flow.md
│   │   └── summary_query_missing.md
│   ├── analysis/                  # Analysis documentation
│   └── guides/                    # User guides and instructions
├── analysis/                      # Analysis tools and reports
│   ├── comparisons/               # Comparison analyses
│   └── reports/                   # Analysis reports
├── archive/                       # Archived and test files
│   └── (test files and old scripts)
├── ALMAGAL/                      # Original ALMAGAL project files
│   ├── data/                     # ALMAGAL data directory
│   ├── software/                 # ALMAGAL software
│   └── README.md                 # ALMAGAL documentation
└── README.md                     # Main project README
```

## File Categories

### Scripts (`scripts/`)

#### FITS Checking (`scripts/fits_checking/`)
- **`fits_check_URL.py`**: Main script for scanning remote ALMA QA2 server for FITS files
- **`fits_check.py`**: Original FITS file checker
- **`fits_check_from_csv.py`**: CSV-based FITS file checker
- **`scan_fits.py`**: Utility for scanning FITS files

#### MS Checking (`scripts/ms_checking/`)
- **`MS_check_URL.py`**: Script for scanning remote ALMA MS server for Measurement Set data

#### Querying (`scripts/querying/`)
- **`query_missing_mous.py`**: Main script for querying ALMA archive using pyvo
- **`query_missing_mous_standalone.py`**: Standalone version without external dependencies
- **`query.py`**: Library of query utilities and functions

#### Utilities (`scripts/utilities/`)
- **`extract_missing_mous.py`**: Extract missing MOUS IDs from FITS checking results
- **`extract_missing_mous_ms.py`**: Extract missing MOUS IDs from MS checking results
- **`compare_missing_mous.py`**: Compare missing MOUS lists between FITS and MS data
- **`download_*.sh`**: Shell scripts for downloading data

### Data Files (`data/`)

#### Input Data (`data/input/`)
- **`tables/database.xlsx`**: Main ALMAGAL database with source and MOUS information
- **`tables/table3.dat`**: Table 3 data
- **`tables/table7.dat`**: Table 7 data

#### Output Data (`data/output/`)
- **FITS Files**: `fits_cont_files_found.csv`, `fits_cube_files_found.csv`
- **MS Files**: `ms_data_found.csv`, `tar_files_found.csv`
- **Missing MOUS**: `missing_mous_list.csv`, `missing_mous_ms_list.csv`
- **Query Results**: `missing_mous_query_results.csv`
- **Checking Results**: `almagal_QA2_*_checking.csv`, `almagal_source_mous_*_checking.csv`
- **JSON Files**: All corresponding JSON output files

### Results (`results/`)

#### Statistics (`results/statistics/`)
- **`statistics_cont.txt`**: Continuum data statistics
- **`statistics_cube.txt`**: Cube data statistics
- **`missing_mous_list.txt`**: Human-readable missing MOUS list
- **`missing_mous_ms_list.txt`**: Human-readable missing MS MOUS list

### Documentation (`docs/`)

#### Reports (`docs/reports/`)
- **`ALMA_Data_Checking_Summary.md`**: Comprehensive project summary
- **`Missing_MOUS_Consistency_Analysis.md`**: Consistency analysis between FITS and MS data
- **`ALAMGAL_Timing_Analysis.md`**: Timing analysis documentation
- **`ALMAGAL_Pipeline_Task.md`**: Pipeline task documentation
- **`ALMAGAL_QA2_Data_Task.md`**: QA2 data task documentation
- **`ALMAGAL_source_split.md`**: Source splitting documentation
- **`conversion_flow.md`**: Data conversion flow documentation
- **`summary_query_missing.md`**: Query missing data summary

## Usage Guidelines

### Running Scripts
1. **FITS Checking**: Run scripts from `scripts/fits_checking/`
2. **MS Checking**: Run scripts from `scripts/ms_checking/`
3. **Querying**: Run scripts from `scripts/querying/`
4. **Utilities**: Run scripts from `scripts/utilities/`

### Data Access
- **Input Data**: Place input files in `data/input/`
- **Output Data**: Generated files will be saved in `data/output/`
- **Results**: Analysis results are stored in `results/`

### Documentation
- **Project Overview**: See `docs/reports/ALMA_Data_Checking_Summary.md`
- **Analysis Reports**: Check `docs/reports/` for detailed analyses
- **File Structure**: This document (`README_File_Structure.md`)

## Benefits of This Structure

1. **Clear Separation**: Scripts, data, and documentation are clearly separated
2. **Functional Organization**: Scripts are grouped by their primary function
3. **Easy Navigation**: Logical directory structure makes finding files intuitive
4. **Maintainability**: Related files are grouped together for easier maintenance
5. **Scalability**: Structure can easily accommodate new scripts and data files
6. **Documentation**: All documentation is centralized and organized

## Migration Notes

- All original functionality is preserved
- File paths in scripts may need updating if they reference other files
- The `ALMAGAL/` directory remains unchanged for compatibility
- Archive directory contains test files and old scripts that are no longer needed

This structure provides a clean, organized, and maintainable framework for the ALMA data checking project.
