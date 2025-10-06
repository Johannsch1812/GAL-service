# Summary of Query
## New Features Added:
1. Content Length Extraction: Now extracts content_length from the datalink results alongside access_url
2. Size Conversion: Converts bytes to GB using the formula: content_length / (1024**3) and rounds to 2 decimal places
3. New CSV Columns:
* Main_Data_Size_GB: Size of the main data file (ending with "of_001.tar")
* Auxiliary_Data_Size_GB: Size of the auxiliary data file (ending with "auxiliary.tar")
4. Enhanced Output:
* Shows file sizes in the console output during processing
* Displays individual file sizes in the summary
* Calculates and reports total data sizes
## Results Summary:
✅ Successfully queried all 37 missing MOUS IDs
✅ 100% success rate - all queries found data
✅ All three URL types found for every MOUS
✅ Size information extracted for all files
## Total Data Sizes:
* Main Data: 5,772.50 GB (5.77 TB)
* Auxiliary Data: 25.06 GB
* Total Data: 5,797.56 GB (5.80 TB)
