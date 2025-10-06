# Quick Reference Guide

## Most Important Files

### Main Scripts
- **FITS Checking**: `scripts/fits_checking/fits_check_URL.py`
- **MS Checking**: `scripts/ms_checking/MS_check_URL.py`
- **Archive Querying**: `scripts/querying/query_missing_mous.py`
- **Compare Results**: `scripts/utilities/compare_missing_mous.py`

### Key Data Files
- **Database**: `data/input/tables/database.xlsx`
- **Missing MOUS (FITS)**: `data/output/missing_mous_list.csv`
- **Missing MOUS (MS)**: `data/output/missing_mous_ms_list.csv`
- **Query Results**: `data/output/missing_mous_query_results.csv`

### Important Reports
- **Project Summary**: `docs/reports/ALMA_Data_Checking_Summary.md`
- **Consistency Analysis**: `docs/reports/Missing_MOUS_Consistency_Analysis.md`

## Quick Commands

### Run FITS Checking
```bash
cd scripts/fits_checking
python fits_check_URL.py
```

### Run MS Checking
```bash
cd scripts/ms_checking
python MS_check_URL.py
```

### Query Missing Data
```bash
cd scripts/querying
python query_missing_mous.py
```

### Compare Missing Lists
```bash
cd scripts/utilities
python compare_missing_mous.py
```

## File Locations

| What You Need | Where to Find It |
|---------------|------------------|
| FITS checking scripts | `scripts/fits_checking/` |
| MS checking scripts | `scripts/ms_checking/` |
| Query scripts | `scripts/querying/` |
| Utility scripts | `scripts/utilities/` |
| Input data | `data/input/` |
| Output CSV/JSON | `data/output/` |
| Statistics | `results/statistics/` |
| Reports | `docs/reports/` |
| Database | `data/input/tables/database.xlsx` |
