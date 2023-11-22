#!/usr/bin/env bash
set -e

if (( $# != 4 )) ; then
	echo "Usage: trtools_prep_beagle_vcf.sh  <vcftype> <ref VCF> <imputed VCF> <output VCF>" 1>&2
	echo "Expects exactly four arguments:" 1>&2
	echo "* original genotyping tool - one of 'advntr' 'eh' 'gangstr' or 'hipstr'" 1>&2
	echo "* Imputation reference panel VCF" 1>&2
	echo "* Imputed genotypes VCF produced by Beagle v5.4" 1>&2
	echo "* Output VCF location" 1>&2
	exit 1
fi

genotyper="$1"
ref_panel="$2"
imputed="$3"
output="$4"

if ! [[ "$genotyper" =~ advntr|eh|gangstr|hipstr ]] ; then
	echo " genotyping tool (first arg) should be one of 'advntr' 'eh' 'gangstr' or 'hipstr'" 1>&2
	exit 2
fi

if ! [[ "$ref_panel" =~ .*\.vcf.gz ]] ; then
	echo "Reference panel VCF (second arg) should be a .vcf.gz file" 1>&2
	exit 3
fi

if ! [[ -f "$ref_panel" ]]  ; then
	echo "Reference panel file not found: ${ref_panel}" 1>&2
	exit 4
fi

if ! [[ "$imputed" =~ .*\.vcf\.gz ]] ; then
	echo "Imputed file (third arg) should be a .vcf.gz file" 1>&2
	exit 5
fi

if ! [[ -f "$imputed" ]]  ; then
	echo "Imputed file not found: ${imputed}" 1>&2
	exit 6
fi

if ! [[ "$output" =~ .*\.vcf\.gz ]] ; then
	echo "Output VCF (fourth arg) should be a .vcf.gz file" 1>&2
	exit 7
fi

if ! tabix --help > /dev/null ; then
	echo "Could not find the tabix command. Exiting"
	exit  8
fi

if ! bcftools --help > /dev/null ; then
	echo "Could not find the bcftools command. Exiting"
	exit  8
fi

echo -n "Creating temporary file ... "
temp_file=$(mktemp)
mv "$temp_file" "$temp_file".vcf
temp_file="$temp_file".vcf
echo "$temp_file"

echo "Copying over meta header lines from the reference panel, then copying the imputed file contents"
line_num=1
# copy contents of $imputed over, inserting a bunch of meta lines in the middle
while read -r line ; do
	echo "$line" >> "$temp_file"
	# this assumes the output meta line ordering of Beagle 5.4 VCFs
	if (( line_num == 3 )) ; then
		for meta_line in '##source' '##command' ; do
			zcat < "$ref_panel" | awk '/^'"$meta_line"'/ {print $0} /^#CHROM/ {exit}' | sed -e 's/##/##preimuptation_/' >> "$temp_file"
		done
		for meta_line in '##contig' '##ALT' '##INFO=<ID=END' ; do
			zcat < "$ref_panel" | awk '/^'"$meta_line"'/ {print $0} /^#CHROM/ {exit}' >> "$temp_file"
		done
	fi
	line_num=$(( line_num+1 ))
done < <(zcat < "$imputed")


echo "bgzipping and tabix indexing"
bgzip "$temp_file"
temp_file="$temp_file".gz
tabix "$temp_file"

# Remove loci that don't have the necessary info fields
if [[ "$genotyper" == advntr ]] ; then
	INFO_fields="RU VID"
elif [[ "$genotyper" == eh ]] ; then
	INFO_fields="RU VARID RL"
elif [[ "$genotyper" == gangstr ]] ; then
	INFO_fields="RU"
elif [[ "$genotyper" == hipstr ]] ; then
	INFO_fields="START END PERIOD"
else 
	echo 'Unexpected genotyper!' 1>&2
	rm "$temp_file"
	rm "$temp_file".tbi
	exit 9
fi

echo "Adding INFO fields and values from the ref_panel, then removing loci which are missing required INFO fields"
# copy INFO fields over
bcftools annotate \
	-a "$ref_panel" \
	-c "$(for field in $INFO_fields ; do echo -n "INFO/$field", ; done | sed -e 's/,$//')" \
	"$temp_file" | \
bcftools filter \
	-i "$( for field in $INFO_fields ; do echo -n '(INFO/'"$field"'!=".") && ' ; done | sed -e 's/ && $//' )" \
	-o "$output" \
	-O z \
	-

tabix "$output"

echo "Removing temp file"
rm "$temp_file"
rm "$temp_file".tbi

echo "Done"
