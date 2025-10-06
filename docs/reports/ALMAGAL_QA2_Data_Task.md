# ✅ Task: ALMAGAL ALMA QA2 Data Verification & Recovery

## 📁 Project: ALMAGAL – ALMA Large Program  
**Proposal ID:** 2019.1.00195.L  
**Target Count:** 1,017 sources  
**Data Types:** FITS (continuum + cube), MS (Measurement Sets)  
**Configurations:** 7M, TM1, TM2  

---

## 🧠 Objective  
Verify availability and completeness of QA2-calibrated ALMA data (FITS + MS) across three configurations, and recover missing datasets via ALMA archive queries.

---

## 📌 Work Summary

### ✅ Developed Scripts
- `fits_check_URL.py`: Checks for continuum/cube FITS files on CASDC server
- `MS_check_URL.py`: Checks for MS directories and untarred data
- `extract_missing_mous.py`: Identifies missing MOUS from FITS results
- `query_missing_mous.py`: Queries ALMA archive via `pyvo` to recover data

### 🔎 Data Sources
- **CASDC QA2 Server:**  
  `https://casdc.china-vo.org/mirror/ALMAGAL/2019.1.00195.L/QA2`
- **CASDC MS Server:**  
  `https://casdc.china-vo.org/mirror/ALMAGAL/2019.1.00195.L/MS`
- **ALMA Archive (VO):**  
  `https://almascience.eso.org/datalink/sync`

---

## 📊 Results Summary

### 🔬 FITS Availability
- 7M: 99.2% found  
- TM1: 8.8% found  
- TM2: 70.6% found  
- **Complete 3-config sources:** 88 / 1,017  
- **Total FITS files (cont + cube):** 9,056

### 📂 MS Availability
- 26 missing MOUS IDs
- Several `.tar.gz` archives found

### 🌐 Archive Recovery
- 37 MOUS queried, 100% success
- ~5.8 TB recoverable
- Main data: 5.77 TB | Auxiliary: 25 GB

---

## ⚙️ Technical Achievements
- Recursive remote directory scanning
- Retry + timeout (120s) + backoff mechanisms
- CSV/JSON/Text-based result outputs
- Full ALMA datalink integration via `pyvo`

---

## 📁 Output Files

- CSV: `fits_cont_files_found.csv`, `ms_data_found.csv`, `missing_mous_list.csv`, etc.  
- JSON: `almagal_QA2_cont_checking.json`, etc.  
- TXT: `statistics_cube.txt`, `missing_mous_list.txt`, etc.

---

## 📌 Recommendations

### 🔧 Recovery Priority
1. TM1 (critical missing)
2. TM2 (partial)
3. 7M (mostly complete)

### 💾 Storage Planning
- Missing data ~5.8 TB
- Reserve ≥7–8 TB for future growth

### 🔁 Data Monitoring
- Automate QA2 sync
- Archive alerts for new availability

---

## ✅ Status: Complete  
- All sources scanned and logged  
- Archive access confirmed and working  
- Download list generated  
- Ready for batch recovery
