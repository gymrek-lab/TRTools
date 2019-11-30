# MergeSTR 

MergeSTR is a tool for merging vcf files from GangSTR. 

## Usage 
To run mergeSTR use the following command: 
```
./mergeSTR.py \
  --vcfs <vcf file1, vcf file2> \
  --out test
  [additional options]
```

Required Parameters: 
* **--vcf <GangSTR VCF>** VCF file output by GangSTR. Must be bgzipped, sorted, and indexed. 
* **--out <string>** prefix to name output files

Special Merge Options: 
* **`--update-sample-from-file`**: Uses file names, rather than sample header names, when merging
* **`--merge-ggl`**: Merge GGL fields (Under construction) 

Optional Additional Parameters: 
* **`--verbose`**: Prints out extra information 
* **`--quiet`**: Doesn't print out anything 

## Basic MergeSTR command 

```
./mergeSTR.py \
   --vcfs ../test_files/test_file.vcf.gz,../test_files/test_file2.vcf.gz
   --out test_run_original
```

Check test_run_original in Example_Files for example output. 

