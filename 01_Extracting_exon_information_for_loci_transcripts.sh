#!/bin/bash


data_dir="path/to/working/directory"

### Extracting loci transcript IDs ###
######################################

input_file="$data_dir/loci_file"
output_file="$input_file.2"
id_counter=1

# Check if input file exists
if [ ! -f "$input_file.gff3" ]; then
    echo "Input file not found: $input_file.gff3"
    exit 1
fi

# Read the input file line by line, skip lines starting with "No results for",
# add IDs to the rest, and write to output file

while IFS= read -r line; do
    # Skip lines starting with "No results for". This is a feature of TargetFinder output. 
    if [[ $line == "No results for"* ]]; then
        continue
    fi

    # Skip comment lines and empty lines
    if [[ $line == \#* || -z $line ]]; then
        echo "$line" >> "$output_file.gff3"
    else
        echo "$id_counter	$line" >> "$output_file.gff3"
        ((id_counter++))
    fi
done < "$input_file.gff3"

echo "Modified GFF3 file with IDs written to: $output_file.gff3"

input_file="$output_file"
output_file="$input_file.3"

# Calculating the length of the loci
awk -F'\t' 'BEGIN {OFS = FS} {difference = $6 - $5 + 1; print $0, difference}' "$input_file.gff3" > "$output_file.gff3"

echo "Modified GFF3 file with calculated values written to: $output_file.gff3"

input_file="$output_file"
output_file="$input_file.4"

# Extracting transcript IDs from the GFF3 file descriptions field
awk -F '\t' '{split($2, arr, " "); print $0 "\t" arr[1]}' $input_file.gff3 > $output_file.gff3

echo "Modified GFF3 file with transcript IDs written to: $output_file.gff3"

input_file="$output_file"
loci_id_file="$input_file.non-redundant_transcript_ids.txt"

# Extracting transcript IDs into a non-redundant list
awk -F '\t' '{split($2, arr, " "); print arr[1]}' $input_file.gff3 | sort | uniq > $loci_id_file

echo "Transcript IDs extracted into non-redundant list of transcripts $loci_id_file"

### Extracting list of exon positions for transcripts ###
#########################################################

input_file="$data_dir/annotation_file"
output_file="$input_file.exons"

# Extracting exons from genome annotation GFF3 file
awk -F '\t' '$3 == "exon"' $input_file.gff3 > $output_file.gff3

echo "Exons extracted from genome annotation file into $output_file.gff3"

# Extracting only the relevant exons from the bigger GFF3 file
input_file="$output_file"
output_file="$output_file.targets"

grep -F -f $loci_id_file $input_file.gff3 > $output_file.gff3

echo "Transcript features extracted from genome annotation file into $output_file.gff3"

# Separating positive and negative strand genes
input_file="$output_file"
output_file_negative="$input_file.negative"
output_file_positive="$input_file.positive"

awk -F '\t' '{ if($7=="-") print $0}' $input_file.gff3 | tr -d '\r' > $output_file_negative.gff3

echo "Negative strand features extracted"

awk -F '\t' '{ if($7=="+") print $0}' $input_file.gff3 | tr -d '\r' > $output_file_positive.gff3

echo "Positive strand features extracted"

# Ordering negative strand exons by exon start location
input_file="$output_file_negative"
output_file="$input_file.2"

sort -k9,9 -k4,4nr "$input_file.gff3" > "$output_file.gff3" #sorting exons by start location, then transcript ID

# Numbering negative strand exons according to their order in this sorted file
input_file="$output_file"
output_file="$input_file.3"

awk 'BEGIN {FS = "\t"; exon_count = 0; last_tid = "" } $9 != last_tid { exon_count = 1 } $9 == last_tid { exon_count++ } { print $0, "\t" "exon_" exon_count; last_tid = $9 }' "$input_file.gff3" > "$output_file.gff3"

# Calculating negative strand exon lengths
input_file="$output_file"
output_file="$input_file.4"

awk -F'\t' '{ diff = $5 - $4 + 1; print $0 "\t" diff }' "$input_file.gff3" > "$output_file.gff3" 

# Ordering positive strand exons by exon start location
input_file="$output_file_positive"
output_file="$input_file.2"

sort -k9,9 -k4,4n "$input_file.gff3" > "$output_file.gff3" #sorting exons by start location, then transcript ID

# Numbering positive strand exons according to their order in this sorted file
input_file="$output_file"
output_file="$input_file.3"

awk 'BEGIN {FS = "\t"; exon_count = 0; last_tid = "" } $9 != last_tid { exon_count = 1 } $9 == last_tid { exon_count++ } { print $0, "\t" "exon_" exon_count; last_tid = $9 }' "$input_file.gff3" > "$output_file.gff3"

# Calculating positive strand exon lengths
input_file="$output_file"
output_file="$input_file.4"

awk -F'\t' '{ diff = $5 - $4 + 1; print $0 "\t" diff }' "$input_file.gff3" > "$output_file.gff3" 
