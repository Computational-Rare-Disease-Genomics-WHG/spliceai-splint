#!/bin/bash

# Exit on errors
set -e

# Function to show usage
usage() {
    echo "Usage: $0 -i <input.vcf> [-o <output.vcf>] --anno_files <annotation1,annotation2,...>"
    echo "  -i, --input    Input VCF file"
    echo "  -o, --output   Output VCF file (default: <input_vcf>_spliceai_supplement.vcf.gz)"
    echo "  --anno_files   Comma-separated list of annotation types to use (options: large_indel, annotation_change, distance, or 'all')"
    exit 1
}

# Parse command-line arguments
ANNOTATIONS=("large_indel" "annotation_change" "distance")  # Default to 'all' if no argument is given
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
        --anno_files)
            if [[ -n "$2" && "$2" != "all" ]]; then
                IFS=',' read -ra ANNOTATIONS <<< "$2"
            fi
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

# Check if the input file is bgzipped, if not, compress it
if [[ "$VCF_FILE" =~ \.vcf$ ]]; then
    echo "Input is not gzipped. Compressing..."
    bgzip -c "$VCF_FILE" > "${VCF_FILE}.gz"
    tabix -p vcf "${VCF_FILE}.gz"
    VCF_FILE="${VCF_FILE}.gz"
    TEMP_GZIPPED=true
elif [[ ! "$VCF_FILE" =~ \.vcf\.gz$ ]]; then
    echo "Error: Input must be a VCF file (.vcf or .vcf.gz)"
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
    DISTANCE_ANNOTATION_FILE="tool_annotation_files/distance_annotation_chr.sorted.txt.gz"
else
    echo "Detected non-'chr'-style chromosomes."
    INDEL_ANNOTATION_FILE="tool_annotation_files/indel_annotation_nochr.sorted.txt.gz"
    SMALLVAR_ANNOTATION_FILE="tool_annotation_files/precomp_supplement_annotation_nochr.sorted.txt.gz"
    DISTANCE_ANNOTATION_FILE="tool_annotation_files/distance_annotation_nochr.sorted.txt.gz"
fi

HEADER_FILE="${SMALLVAR_ANNOTATION_FILE/_chr.sorted.txt.gz/.hdr}"
HEADER_FILE="${HEADER_FILE/_nochr.sorted.txt.gz/.hdr}"

# Output file names for intermediate files
OUTPUT_LARGE_INDEL="large_indel_annotated.vcf.gz"
OUTPUT_SMALL_VARS="small_vars_annotated.vcf.gz"
OUTPUT_DISTANCE="distance_annotated.vcf.gz"

# Annotation steps based on user selection
ANNOTATED_FILES=()
if [[ " ${ANNOTATIONS[@]} " =~ " large_indel " ]]; then
    if [[ ! -f "$INDEL_ANNOTATION_FILE" ]]; then
        echo "Error: Annotation file '$INDEL_ANNOTATION_FILE' not found!"
        exit 1
    fi
    echo "Annotating large indels..."
    bcftools filter -i '(ILEN < -4 | ILEN > 1)' $VCF_FILE | \
        bcftools annotate -a $INDEL_ANNOTATION_FILE -h $HEADER_FILE \
        -c CHROM,FROM,TO,SPLICEAI_RECOMMENDATION,REASON \
        -Oz -o "$OUTPUT_LARGE_INDEL"
    tabix -p vcf "$OUTPUT_LARGE_INDEL"
    ANNOTATED_FILES+=("$OUTPUT_LARGE_INDEL")
fi

if [[ " ${ANNOTATIONS[@]} " =~ " annotation_change " ]]; then
    if [[ ! -f "$SMALLVAR_ANNOTATION_FILE" ]]; then
        echo "Error: Annotation file '$SMALLVAR_ANNOTATION_FILE' not found!"
        exit 1
    fi
    echo "Annotating remaining variants with precomputed annotations..."
    bcftools filter -e '(ILEN < -4 | ILEN > 1)' "$VCF_FILE" | \
    bcftools annotate --force -a "$SMALLVAR_ANNOTATION_FILE" -h "$HEADER_FILE" \
        -c CHROM,FROM,TO,SPLICEAI_RECOMMENDATION,REASON \
        -Oz -o "$OUTPUT_SMALL_VARS"
    tabix -p vcf "$OUTPUT_SMALL_VARS"
    
    # Only add to ANNOTATED_FILES if distance annotation is not being applied later
    if [[ ! " ${ANNOTATIONS[@]} " =~ " distance " ]]; then
        ANNOTATED_FILES+=("$OUTPUT_SMALL_VARS")
    fi
fi

if [[ " ${ANNOTATIONS[@]} " =~ " distance " ]]; then
    if [[ ! -f "$DISTANCE_ANNOTATION_FILE" ]]; then
        echo "Error: Annotation file '$DISTANCE_ANNOTATION_FILE' not found!"
        exit 1
    fi
    echo "Annotating with distance annotations..."
    
    # Use OUTPUT_SMALL_VARS if annotation_change was applied, otherwise use VCF_FILE
    INPUT_FILE="$OUTPUT_SMALL_VARS"
    if [[ -z "$INPUT_FILE" ]]; then
        INPUT_FILE="$VCF_FILE"
    fi
    
    bcftools annotate --force -a "$DISTANCE_ANNOTATION_FILE" -h "$HEADER_FILE" \
        -c CHROM,FROM,TO,SPLICEAI_RECOMMENDATION,REASON \
        -Oz -o "$OUTPUT_DISTANCE" "$INPUT_FILE"
    tabix -p vcf "$OUTPUT_DISTANCE"
    ANNOTATED_FILES+=("$OUTPUT_DISTANCE")  # Only store the final file
fi

# Combine and sort the final output VCF
echo "Creating and indexing the final annotated VCF..."
bcftools concat -a -O z ${ANNOTATED_FILES[@]} | \
    bcftools sort -Oz -o $FINAL_OUTPUT

tabix -p vcf "$FINAL_OUTPUT"

# Clean up intermediate files
echo "Cleaning up intermediate files..."
rm -f ${ANNOTATED_FILES[@]} ${ANNOTATED_FILES[@]/%.vcf.gz/.vcf.gz.tbi}
if [[ "$TEMP_GZIPPED" == true ]]; then
    echo "Cleaning up temporary gzipped file..."
    rm -f "$VCF_FILE" "$VCF_FILE.tbi"
fi
echo "Annotation complete. Final output saved to $FINAL_OUTPUT"