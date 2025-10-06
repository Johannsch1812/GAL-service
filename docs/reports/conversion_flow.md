Conversion Workflow

convert_msv1_to_msv2.sh checks an MS directory, moves the MSv1 copy into calibrated_v1, and launches CASA with convert_msv1_to_msv2_batch.py to rewrite the Measurement Set in MSv2 format (data/mpi-runs/convert_msv1_to_msv2.sh:1-101). It now logs run-times per UID to process_ms_data_timings.log.
convert_msv1_to_msv2_batch.py (data/mpi-runs/convert_msv1_to_msv2_batch.py:1-70) is the CASA-side script that runs mstransform and reports success/failure.
convert_msv1_to_msv2_simple.py (data/mpi-runs/convert_msv1_to_msv2_simple.py:1-114) offers the same conversion logic without the wrapper’s safety steps; useful for manual one-offs.
configALMAGAL.py (data/mpi-runs/configALMAGAL.py:107-111) provides the path definitions (my_mainPath, my_storagePath, etc.) consumed by every script; keeping these accurate underpins the whole pipeline.
process_ms_data.sh orchestrates batch conversions/splits from a CSV, calling convert_msv1_to_msv2.sh, split_sources.sh, and update_table_info.sh as needed; it records per-entry timing in process_ms_data_timings.log and process_ms_data_split_timings.log (data/mpi-runs/process_ms_data.sh:1-138).
update_table_info.sh (not opened but included in the zip) normalises metadata ahead of splits so CASA sees the MS as MSv2.
split_sources.sh (data/mpi-runs/split_sources.sh:1-86) validates MSv2 input and runs scriptToSplitSources.py.
scriptToSplitSources.py (data/mpi-runs/scriptToSplitSources.py:1-360) is the CASA driver that executes split, tars each per-source MS, and copies it to the storage tree.
Key Observations

Conversion scripts run reliably; the long runtimes reflect full MS rewrites. Batch operation is safe because each UID is processed sequentially, with a rollback path (calibrated_v1) if CASA fails.
Source splitting succeeds, but we saw “FAILED” entries in the timing log because the copy step tries to scp a tar into an identical path when my_storagePath == my_mainPath. The split completes and tar files are valid, but the guard should be tightened (skip copy when source and destination match or provide a distinct storage location).
Timing logs now capture start/end timestamps and durations, which helps track slow datasets.
Temporary scripts (tmpExecute_*.py) are generated per source by the splitter; we haven’t modified those directly.
Overall the conversion pipeline is in good shape—MSv1→MSv2 works, batch tooling is logging properly, and the splitter produces tar archives. The outstanding clean-up is to prevent the redundant copy in scriptToSplitSources.py from aborting the shell run so the timing log reflects actual success.
