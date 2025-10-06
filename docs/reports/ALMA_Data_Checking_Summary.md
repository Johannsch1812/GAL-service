# ALMA Data Checking Work Summary

## Overview

This document summarizes the comprehensive ALMA data checking work performed on the ALMAGAL project data. The work involved checking the availability of FITS files, MS (Measurement Set) data, and querying missing data from the ALMA archive.

## Project Context

- **Project**: ALMAGAL (ALMA Large Program)
- **Proposal ID**: 2019.1.00195.L
- **Total Sources**: 1,017 targets
- **Data Types Checked**: FITS files (continuum and cube), MS data
- **Configurations**: 7M, TM1, TM2

## Data Sources

### Remote Data Servers
- **FITS Files**: `https://casdc.china-vo.org/mirror/ALMAGAL/2019.1.00195.L/QA2`
- **MS Data**: `https://casdc.china-vo.org/mirror/ALMAGAL/2019.1.00195.L/MS`
- **ALMA Archive**: `https://almascience.eso.org/datalink/sync`

### Local Database
- **File**: `database.xlsx`
- **Columns**: Source, SGOUS, GOUS, MOUS_TM1, MOUS_TM2, MOUS_7M

## Scripts Developed

### 1. FITS File Checking (`fits_check_URL.py`)
- **Purpose**: Scan remote ALMA QA2 server for FITS files
- **Pattern**: `cont.I.pbcor.fits` and `cube.I.pbcor.fits`
- **Features**:
  - Recursive URL scanning
  - Retry mechanism with exponential backoff
  - CSV caching system
  - Configurable file patterns

### 2. MS Data Checking (`MS_check_URL.py`)
- **Purpose**: Scan remote ALMA MS server for Measurement Set data
- **Pattern**: `uid___<MOUS>` directories containing `table.info`
- **Features**:
  - Prioritizes `_calibrated_final.ms/` directories
  - Detects tar.gz files (untarred MS data)
  - Database preprocessing (removes Source column, deduplicates)

### 3. Missing MOUS Extraction
- **`extract_missing_mous.py`**: Extracts missing MOUS IDs from FITS checking results
- **`extract_missing_mous_simple.py`**: Extracts missing MOUS IDs from MS checking results

### 4. ALMA Archive Querying (`query_missing_mous.py`)
- **Purpose**: Query ALMA archive for missing data
- **Method**: Uses `pyvo` to access ALMA datalink service
- **Output**: Download URLs and file sizes for main data, auxiliary data, and README files

## Results Summary

### FITS Files Availability

#### Continuum FITS Files
- **Total Sources**: 1,017
- **7M Configuration**: 1,009 found, 8 missing (99.2%)
- **TM1 Configuration**: 89 found, 928 missing (8.8%)
- **TM2 Configuration**: 718 found, 299 missing (70.6%)
- **Sources with all 3 configurations**: 88/1,017 (8.7%)
- **Total FITS files found**: 1,816

#### Cube FITS Files
- **Total Sources**: 1,017
- **7M Configuration**: 1,009 found, 8 missing (99.2%)
- **TM1 Configuration**: 89 found, 928 missing (8.8%)
- **TM2 Configuration**: 718 found, 299 missing (70.6%)
- **Sources with all 3 configurations**: 88/1,017 (8.7%)
- **Total FITS files found**: 7,240

### MS Data Availability
- **Missing MOUS IDs**: 26 unique MOUS IDs
- **Tar.gz files detected**: Multiple compressed MS data files found

### Missing Data Query Results

#### Archive Query Statistics
- **Total MOUS IDs queried**: 37
- **Successful queries**: 37 (100% success rate)
- **Failed queries**: 0
- **Main data URLs found**: 37
- **Auxiliary data URLs found**: 37
- **README URLs found**: 37

#### Data Size Summary
- **Total main data size**: 5,772.50 GB (5.77 TB)
- **Total auxiliary data size**: 25.06 GB
- **Total data size**: 5,797.56 GB (5.80 TB)

## Key Findings

### 1. Data Availability Patterns
- **7M Configuration**: Excellent coverage (99.2% for both continuum and cube)
- **TM1 Configuration**: Very low coverage (8.8% for both continuum and cube)
- **TM2 Configuration**: Moderate coverage (70.6% for both continuum and cube)

### 2. Missing Data Characteristics
- **Primary missing data**: TM1 and TM2 configurations
- **Missing MOUS IDs**: 37 unique identifiers
- **Data size range**: Individual files range from 3.38 GB to 302.28 GB
- **Average main data size**: ~156 GB per MOUS
- **Average auxiliary data size**: ~0.68 GB per MOUS

### 3. Data Distribution
- **Largest single file**: 302.28 GB (uid___A001_X146c_Xe6)
- **Smallest single file**: 3.38 GB (uid___A001_X1467_X218)
- **Most common size range**: 130-220 GB for main data files

## Technical Achievements

### 1. Robust Data Access
- Implemented retry mechanisms for network requests
- Added timeout handling (120 seconds)
- Created caching systems to avoid redundant scans

### 2. Comprehensive Coverage
- Checked all three ALMA configurations
- Verified both continuum and cube data
- Included MS data verification
- Queried ALMA archive for missing data

### 3. Data Quality Assurance
- Validated file patterns and naming conventions
- Cross-referenced with ALMA archive
- Provided detailed size information for planning

## Output Files Generated

### CSV Files
- `fits_cont_files_found.csv`: Continuum FITS files inventory
- `fits_cube_files_found.csv`: Cube FITS files inventory
- `ms_data_found.csv`: MS data inventory
- `tar_files_found.csv`: Compressed MS data inventory
- `missing_mous_list.csv`: Missing MOUS IDs for FITS data
- `missing_mous_ms_list.csv`: Missing MOUS IDs for MS data
- `missing_mous_query_results.csv`: ALMA archive query results

### JSON Files
- `almagal_source_mous_checking.json`: Source-MOUS mapping for FITS data
- `almagal_QA2_cont_checking.json`: Continuum checking results
- `almagal_QA2_cube_checking.json`: Cube checking results

### Text Files
- `statistics_cont.txt`: Continuum data statistics
- `statistics_cube.txt`: Cube data statistics
- `missing_mous_list.txt`: Human-readable missing MOUS list
- `missing_mous_ms_list.txt`: Human-readable missing MS MOUS list

## Recommendations

### 1. Data Recovery Priority
1. **High Priority**: TM1 configuration data (8.8% coverage)
2. **Medium Priority**: TM2 configuration data (70.6% coverage)
3. **Low Priority**: 7M configuration data (99.2% coverage)

### 2. Storage Planning
- **Total missing data**: ~5.8 TB
- **Recommended storage**: At least 7-8 TB for safe margins
- **Download strategy**: Prioritize by configuration and file size

### 3. Data Management
- Implement automated monitoring for new data availability
- Set up alerts for data recovery completion
- Maintain regular synchronization with ALMA archive

## Conclusion

The ALMA data checking work has successfully identified and quantified missing data across all three ALMA configurations. The comprehensive analysis provides a clear roadmap for data recovery efforts, with detailed information about file locations, sizes, and download URLs. The 100% success rate in querying the ALMA archive ensures that all missing data can be recovered when available.

The work demonstrates the importance of systematic data verification in large astronomical projects and provides a robust framework for ongoing data management and monitoring.
