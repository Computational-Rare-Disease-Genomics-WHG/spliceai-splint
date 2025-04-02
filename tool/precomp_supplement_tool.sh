#!/bin/bash

# Exit on errors
set -e

# Function to show usage
usage() {
    echo "Usage: $0 -i <input.vcf> [-o <output.vcf>]"
    echo "  -i, --input    Input VCF file"
    echo "  -o, --output   Output VCF file (default: <input_vcf>_spliceai_supplement.vcf.gz)"
    exit 1
}

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)
            VCF_FILE=$2
            shift 2
            ;;
        -o|--output)
            FINAL_OUTPUT=$2
            shift 2
            ;;
        *)
            usage
            ;;
    esac
done

# Check if input VCF file is provided
if [[ -z "$VCF_FILE" ]]; then
    echo "Error: Input VCF file is required."
    usage
fi
# Check if the input VCF file ends with .vcf.gz
if [[ ! "$VCF_FILE" =~ \.vcf\.gz$ ]]; then
    echo "Error: Input must be a bgzipped VCF file (.vcf.gz)"
    exit 1
fi

# Set default output if not specified
FINAL_OUTPUT=${FINAL_OUTPUT:-"${VCF_FILE/.vcf.gz/}.spliceai_supplement.vcf.gz"}

# Check the chromosome naming style in the VCF file
echo "Checking chromosome naming style in the input VCF..."
FIRST_CHROM=$(gunzip -c "$VCF_FILE" | grep -v "^#" | awk '{print $1; exit}')
if [[ -z "$FIRST_CHROM" ]]; then
    echo "Error: Could not determine chromosome style from the input VCF."
    exit 1
elif [[ "$FIRST_CHROM" =~ ^chr ]]; then
    echo "Detected 'chr'-style chromosomes."
    INDEL_ANNOTATION_FILE="tool_annotation_files/indel_annotation_chr.sorted.txt.gz"
    SMALLVAR_ANNOTATION_FILE="tool_annotation_files/precomp_supplement_annotation_chr.sorted.txt.gz"
else
    echo "Detected non-'chr'-style chromosomes."
    INDEL_ANNOTATION_FILE="tool_annotation_files/indel_annotation_nochr.sorted.txt.gz"
    SMALLVAR_ANNOTATION_FILE="tool_annotation_files/precomp_supplement_annotation_nochr.sorted.txt.gz"
fi

HEADER_FILE="${SMALLVAR_ANNOTATION_FILE/_chr.sorted.txt.gz/.hdr}"
HEADER_FILE="${HEADER_FILE/_nochr.sorted.txt.gz/.hdr}"

if [[ ! -f "$INDEL_ANNOTATION_FILE" ]]; then
    echo "Error: Annotation file '$INDEL_ANNOTATION_FILE' not found!"
    exit 1
fi
if [[ ! -f "$SMALLVAR_ANNOTATION_FILE" ]]; then
    echo "Error: Annotation file '$SMALLVAR_ANNOTATION_FILE' not found!"
    exit 1
fi
if [[ ! -f "$HEADER_FILE" ]]; then
    echo "Error: Annotation file '$HEADER_FILE' not found!"
    exit 1
fi


# Output file names for intermediate files
OUTPUT_LARGE_INDEL="large_indel_annotated.vcf.gz"
OUTPUT_SMALL_VARS="small_vars_annotated.vcf.gz"

# Step 2: Annotate large indels (using bcftools filter directly for length-based filtering)
echo "Annotating large indels..."
bcftools filter -i '(ILEN < -4 | ILEN > 1)' $VCF_FILE | \
    bcftools annotate -a $INDEL_ANNOTATION_FILE -h $HEADER_FILE \
    -c CHROM,FROM,TO,SPLICEAI_RECOMMENDATION,REASON \
    -Oz -o "$OUTPUT_LARGE_INDEL"

# Create index for large indels output
echo "Indexing large indels VCF..."
tabix -p vcf "$OUTPUT_LARGE_INDEL"

# Step 3: Annotate remaining variants with precomputed annotations
echo "Annotating remaining variants with precomputed annotations..."
bcftools filter -e '(ILEN < -4 | ILEN > 1)' $VCF_FILE | \
bcftools annotate --force -a $SMALLVAR_ANNOTATION_FILE -h $HEADER_FILE \
    -c CHROM,FROM,TO,SPLICEAI_RECOMMENDATION,REASON \
    -Oz -o "$OUTPUT_SMALL_VARS"

# Create index for small variants output
echo "Indexing small variants VCF..."
tabix -p vcf "$OUTPUT_SMALL_VARS"

# Step 4: Combine and sort the final output VCF
echo "Creating and indexing the final annotated VCF..."
bcftools concat -a -O z $OUTPUT_LARGE_INDEL $OUTPUT_SMALL_VARS | \
    bcftools sort -Oz -o $FINAL_OUTPUT

# Index the final combined output VCF
echo "Indexing the final output VCF..."
tabix -p vcf "$FINAL_OUTPUT"

# Clean up intermediate files
echo "Cleaning up intermediate files..."
rm -f "$OUTPUT_LARGE_INDEL" "${OUTPUT_LARGE_INDEL}.tbi" "$OUTPUT_SMALL_VARS" "${OUTPUT_SMALL_VARS}.tbi" "temp_final.vcf.gz"

echo "Annotation complete. Final output saved to $FINAL_OUTPUT"