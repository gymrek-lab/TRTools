# MergeSTR 

MergeSTR is a tool for merging vcf files from GangSTR. 

## Usage 
To run mergeSTR use the following command: 
```
./mergeSTR.py \
  --vcfs ERR1955393.17.sorted.vcf.gz,ERR1955394.17.sorted.vcf.gz \
  --out test
  [additional options]
```

Required Parameters: 
* **--vcf <GangSTR VCF>** VCF file output by GangSTR. Can be unzipped or bgzipped.
* **--out <string>** prefix to name output files

Special Merge Options: 
* **`--update-sample-from-file`**: Uses file names, rather than sample header names, when merging
* **`--merge-ggl`**: Merge GGL fields 

Optional Additional Parameters: 
* **`--verbose`**: Prints out extra information 
* **`--quiet`**: Doesn't print out anything 
