#!/bin/bash

# Exit on errors
set -e

# Function to show usage
usage() {
    echo "Usage: $0 -i <input.vcf> [-o <output.vcf>] --anno_files <annotation1,annotation2,...>"
    echo "  -i, --input    Input VCF file"
    echo "  -o, --output   Output VCF file (default: <input_vcf>_spliceai_supplement.vcf.gz)"
    exit 1
}

# Parse command-line arguments
ANNOTATIONS=("large_indel" "annotation_change" "distance", "wrong_ref")  # Default to 'all' if no argument is given
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

# Check that the VCF header indicates GRCh38
if ! zgrep -q "^##reference=GRCh38" "$VCF_FILE"; then
    echo "Error: Input VCF does not declare reference genome as GRCh38 in the header. SpliceAI-splint can only be run on GRCh38 VCFs"
    exit 1
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
FINAL_OUTPUT=${FINAL_OUTPUT:-"${VCF_FILE/.vcf.gz/}.spliceai_splint.vcf.gz"}

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
    WRONGREF_ANNOTATION_FILE="tool_annotation_files/wrongref_annotation_chr.sorted.txt.gz"
else
    echo "Detected non-'chr'-style chromosomes."
    INDEL_ANNOTATION_FILE="tool_annotation_files/indel_annotation_nochr.sorted.txt.gz"
    SMALLVAR_ANNOTATION_FILE="tool_annotation_files/precomp_supplement_annotation_nochr.sorted.txt.gz"
    DISTANCE_ANNOTATION_FILE="tool_annotation_files/distance_annotation_nochr.sorted.txt.gz"
    WRONGREF_ANNOTATION_FILE="tool_annotation_files/wrongref_annotation_nochr.sorted.txt.gz"
fi

HEADER_FILE="${SMALLVAR_ANNOTATION_FILE/_chr.sorted.txt.gz/.hdr}"
HEADER_FILE="${HEADER_FILE/_nochr.sorted.txt.gz/.hdr}"

# Output file names for intermediate files
OUTPUT_LARGE_INDEL="large_indel_annotated.vcf.gz"
OUTPUT_SMALL_VARS="small_vars_annotated.vcf.gz"
OUTPUT_DISTANCE="distance_annotated.vcf.gz"
OUTPUT_WRONGREF="wrong_ref_annotated.vcf.gz"

# add indel 
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

# add transcript changes
if [[ ! -f "$SMALLVAR_ANNOTATION_FILE" ]]; then
    echo "Error: Annotation file '$SMALLVAR_ANNOTATION_FILE' not found!"
    exit 1
fi
echo "Annotating remaining variants if in transcript change regions..."
bcftools filter -e '(ILEN < -4 | ILEN > 1)' "$VCF_FILE" | \
bcftools annotate --force -a "$SMALLVAR_ANNOTATION_FILE" -h "$HEADER_FILE" \
    -c CHROM,FROM,TO,SPLICEAI_RECOMMENDATION,REASON \
    -Oz -o "$OUTPUT_SMALL_VARS"
tabix -p vcf "$OUTPUT_SMALL_VARS"

# add distance
if [[ ! -f "$DISTANCE_ANNOTATION_FILE" ]]; then
    echo "Error: Annotation file '$DISTANCE_ANNOTATION_FILE' not found!"
    exit 1
fi
echo "Annotating variants if we recommend running with increased distance parameter..."

bcftools annotate --force -a "$DISTANCE_ANNOTATION_FILE" -h "$HEADER_FILE" \
    -c CHROM,FROM,TO,SPLICEAI_RECOMMENDATION,REASON \
    -Oz -o "$OUTPUT_DISTANCE" "$OUTPUT_SMALL_VARS"
tabix -p vcf "$OUTPUT_DISTANCE"

if [[ ! -f "$WRONGREF_ANNOTATION_FILE" ]]; then
    echo "Error: Annotation file '$WRONGREF_ANNOTATION_FILE' not found!"
    exit 1
fi
echo "Annotating variants with wrong reference allele..."


bcftools annotate --force -a "$WRONGREF_ANNOTATION_FILE" -h "$HEADER_FILE"  \
    -c CHROM,FROM,TO,SPLICEAI_RECOMMENDATION,REASON \
    -Oz -o "$OUTPUT_WRONGREF" "$OUTPUT_DISTANCE"
tabix -p vcf "$OUTPUT_WRONGREF"

ANNOTATED_FILES+=("$OUTPUT_WRONGREF")

# Combine and sort the final output VCF
echo "Creating and indexing the final annotated VCF..."
echo $ANNOTATED_FILES
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