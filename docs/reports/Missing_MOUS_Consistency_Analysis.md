# Missing MOUS Consistency Analysis

## Overview

This analysis compares the missing MOUS IDs between FITS data and MS data to identify consistency patterns and potential data availability issues.

## Summary Statistics

| Metric | Count | Percentage |
|--------|-------|------------|
| **FITS Data Missing MOUS IDs** | 37 | 100% |
| **MS Data Missing MOUS IDs** | 25 | 100% |
| **Common Missing MOUS IDs** | 19 | 51.4% of FITS, 76.0% of MS |
| **FITS-only Missing** | 18 | 48.6% of FITS |
| **MS-only Missing** | 6 | 24.0% of MS |

## Consistency Analysis

### Consistency Score: 44.2% (Low Consistency)

The low consistency score indicates that there are significant differences between what's missing in FITS data versus MS data.

## Detailed Breakdown

### 1. Common Missing MOUS IDs (19)
These MOUS IDs are missing from BOTH FITS and MS data, suggesting systematic data availability issues:

```
uid___A001_X1467_X1e5
uid___A001_X1467_X1f1
uid___A001_X1467_X1f8
uid___A001_X1467_X1ff
uid___A001_X1467_X206
uid___A001_X1467_X20d
uid___A001_X1467_X214
uid___A001_X1467_X21b
uid___A001_X1467_X222
uid___A001_X1467_X229
uid___A001_X1467_X230
uid___A001_X1467_X237
uid___A001_X146c_Xac
uid___A001_X146c_Xcf
uid___A001_X146c_Xd6
uid___A001_X146c_Xdd
uid___A001_X146c_Xe4
uid___A001_X146c_Xeb
uid___A001_X146c_Xf9
```

### 2. FITS-only Missing MOUS IDs (18)
These MOUS IDs are missing from FITS data but present in MS data:

```
uid___A001_X1467_X1d5
uid___A001_X1467_X1dc
uid___A001_X1467_X1de
uid___A001_X1467_X1e3
uid___A001_X1467_X1ea
uid___A001_X1467_X1ec
uid___A001_X1467_X1fa
uid___A001_X1467_X216
uid___A001_X1467_X218
uid___A001_X1467_X22b
uid___A001_X1467_X232
uid___A001_X1467_X23e
uid___A001_X1467_X240
uid___A001_X146c_X97
uid___A001_X146c_Xc1
uid___A001_X146c_Xc8
uid___A001_X146c_Xe6
uid___A001_X146c_Xf2
```

### 3. MS-only Missing MOUS IDs (6)
These MOUS IDs are missing from MS data but present in FITS data:

```
uid___A001_X146c_X9e
uid___A001_X146c_Xa5
uid___A001_X146c_Xb7
uid___A001_X146c_Xba
uid___A001_X146c_Xca
uid___A001_X146c_Xdf
```

## Key Findings

### 1. Systematic Data Issues
- **19 MOUS IDs** (51.4% of FITS missing) are missing from both data types
- This suggests fundamental data availability problems for these specific MOUS IDs
- These should be the highest priority for data recovery

### 2. Data Processing Differences
- **18 MOUS IDs** are missing from FITS but present in MS data
- **6 MOUS IDs** are missing from MS but present in FITS data
- This suggests different data processing pipelines or availability windows

### 3. Data Type Specificity
- MS data appears to have better coverage (25 missing vs 37 for FITS)
- Some MOUS IDs may have been processed for MS but not for FITS, or vice versa

## Recommendations

### 1. Priority Actions
1. **Immediate**: Focus on the 19 common missing MOUS IDs for data recovery
2. **Investigate**: Check why 18 MOUS IDs are missing from FITS but present in MS
3. **Verify**: Confirm why 6 MOUS IDs are missing from MS but present in FITS

### 2. Data Quality Checks
- Verify the accuracy of both missing data lists
- Check if there are any data processing errors
- Confirm the completeness of the MS data scanning

### 3. Recovery Strategy
- Prioritize recovery of common missing MOUS IDs
- Investigate the discrepancies in data type coverage
- Implement cross-validation between FITS and MS data availability

## Conclusion

The analysis reveals significant inconsistencies between FITS and MS data availability, with only 44.2% consistency in missing data patterns. This suggests that:

1. **Systematic issues** affect 19 MOUS IDs across both data types
2. **Data processing differences** result in different coverage patterns
3. **Further investigation** is needed to understand the root causes of these discrepancies

The findings highlight the importance of cross-validating data availability across different data types and implementing robust data quality assurance processes.
