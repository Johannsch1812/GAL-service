# GAL-service
GAL-service is a set of Python utilities for scanning, cataloging, and checking ALMA FITS files, particularly for the ALMAGAL project. The tools are designed to help you:

- Recursively scan directories for FITS files of a specific type (e.g., continuum or cube)
- Generate CSV catalogs of discovered files
- Cross-check FITS files against a database of targets and MOUS IDs

## Contents

- `scan_fits.py`: Recursively scans a directory for FITS files matching given patterns and outputs a CSV catalog.
- `fits_check.py`: Checks for the presence of expected FITS files for each target/MOUS in the ALMAGAL database.
- `fits_check_from_csv.py`: Similar to `fits_check.py`, but works from pre-generated CSV catalogs (from `scan_fits.py`).

## Quick Start

### 1. Scan for FITS files

To scan a directory for FITS files and output a CSV catalog:
