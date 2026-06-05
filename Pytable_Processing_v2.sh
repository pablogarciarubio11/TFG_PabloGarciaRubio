#!/bin/bash

#The directories used were coding and noncoding, this is just an example
WORK_DIR="D:\TFG\Filtrado_noncoding"
OUTPUT_DIR="D:\TFG\Pytables_noncoding"

#Create the output directory if it does not exist
mkdir -p "$OUTPUT_DIR"

if [[ ! -d "$WORK_DIR" ]]; then
    echo "Error: The directory '$WORK_DIR' does not exist"
    exit 1
fi

cd "$WORK_DIR"

count=0

echo "Starting the processing of Pytables files..."
echo "Origin directory: $(pwd)"
echo "Destination directory: $OUTPUT_DIR"
echo "----------------------------------------"

for file in GSM*_out_filtered.h5; do
    if [[ ! -f "$file" ]]; then
        echo "Files were not found"
        continue
    fi

    GSM_NUM=""
    output_file=""

    # Extract GSM
    if [[ "$file" =~ (GSM[0-9]+) ]]; then
        GSM_NUM=${BASH_REMATCH[1]}
    else
        echo "The GSM could not be extracted: $file"
        continue
    fi

    output_file="$OUTPUT_DIR/${GSM_NUM}.h5"

    echo "Processing: $file"
    echo "  → $output_file"

    rm -f "$output_file"

    if ptrepack --complevel 5 --overwrite-nodes "${file}:/matrix" "${output_file}:/matrix"; then
        echo "OK"
        ((count++))
    else
        echo "Error in $file"
    fi

    echo "----------------------------------------"
done

echo "Processing completed: $count files"
echo "Files saved in: $OUTPUT_DIR"


