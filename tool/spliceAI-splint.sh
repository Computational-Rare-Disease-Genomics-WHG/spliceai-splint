#!/bin/bash

# Exit on errors
set -e

# Function to show usage
usage() {
    cat <<'USAGE'
Usage:   spliceAI-splint.sh -i <input.vcf> [-o <output.vcf>] [-g <genome>] [-sc]
Options:
  -i, --input            Input VCF file (.vcf or .vcf.gz) (required)
  -o, --output           Output VCF file (default: <input_vcf>.spliceai_splint.vcf.gz)
  -g, --genome-version   User-supplied genome version (e.g. GRCh38) -- skips auto GRCh38 header check
  -sc, --skip-checks     Auto-insert a compatible ##fileformat=VCFvX.X header if missing
USAGE
    exit 1
}

# -------------------------
# Argument parsing
# -------------------------
ANNOTATIONS=("large_indel" "annotation_change" "distance" "wrong_ref")
GENOME_VERSION=""
SKIP_CHECKS=false
VCF_FILE=""
FINAL_OUTPUT=""

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
        -g|--genome-version)
            GENOME_VERSION=$2
            shift 2
            ;;
        -sc|--skip-checks)
            SKIP_CHECKS=true
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown argument: $1"
            usage
            ;;
    esac
done

# Check if input VCF file is provided
if [[ -z "$VCF_FILE" ]]; then
    echo "Error: Input VCF file is required."
    usage
fi

if [[ ! -f "$VCF_FILE" ]]; then
    echo "Error: Input VCF file '$VCF_FILE' not found."
    echo "Please check the file path and try again."
    exit 1
fi

# -------------------------
# Tool presence & version checks
# -------------------------
# helper: compare semantic versions (returns 0 if $1 >= $2)
ver_gte() {
    # safe numeric compare that handles missing parts
    IFS=. read -r a1 a2 a3 <<<"$1"
    IFS=. read -r b1 b2 b3 <<<"$2"
    a1=${a1:-0}; a2=${a2:-0}; a3=${a3:-0}
    b1=${b1:-0}; b2=${b2:-0}; b3=${b3:-0}
    if ((10#$a1 > 10#$b1)); then return 0; fi
    if ((10#$a1 < 10#$b1)); then return 1; fi
    if ((10#$a2 > 10#$b2)); then return 0; fi
    if ((10#$a2 < 10#$b2)); then return 1; fi
    if ((10#$a3 >= 10#$b3)); then return 0; else return 1; fi
}

missing_tools=()
for tool in bcftools tabix bgzip; do
    if ! command -v "$tool" &>/dev/null; then
        missing_tools+=("$tool")
    fi
done

if (( ${#missing_tools[@]} > 0 )); then
    echo "Error: Required tools not found: ${missing_tools[*]}"
    echo "Installation requirements:"
    echo "  - bcftools (v1.10+)"
    echo "  - tabix/bgzip (part of htslib)"
    echo "  - Bash (Linux/macOS environment)"
    exit 1
fi

# Check bcftools version
BCFVER_RAW=$(bcftools --version 2>/dev/null | head -n1 | awk '{print $2}')
if [[ -z "$BCFVER_RAW" ]]; then
    echo "Warning: Could not parse bcftools version. Ensure bcftools >= 1.10 is installed."
else
    if ! ver_gte "$BCFVER_RAW" "1.10"; then
        echo "Error: Detected bcftools version $BCFVER_RAW. This script requires bcftools >= 1.10."
        echo "Please upgrade bcftools (see https://samtools.github.io/bcftools/)."
        exit 1
    fi
fi

# -------------------------
# Helper: pick compatible VCF fileformat by probing bcftools
# -------------------------
pick_compatible_vcf_version() {
    local candidates=("##fileformat=VCFv4.5" "##fileformat=VCFv4.3" "##fileformat=VCFv4.2")
    for candidate in "${candidates[@]}"; do
        tmp=$(mktemp)
        # minimal VCF with that fileformat and a tiny record
        printf "%s\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t1\t.\tA\tT\t.\t.\t.\n" "$candidate" > "$tmp"
        # bcftools view -h will exit 0 for acceptable file; suppress output
        if bcftools view -h "$tmp" &>/dev/null || bcftools view "$tmp" &>/dev/null; then
            rm -f "$tmp"
            echo "$candidate"
            return 0
        fi
        rm -f "$tmp"
    done
    return 1
}

# -------------------------
# VCF fileformat header check & optional auto-insert
# -------------------------
check_fileformat_header() {
    local file=$1
    local desired_line

    # read existing fileformat header if present
    if [[ "$file" =~ \.vcf\.gz$ ]]; then
        header_line=$(zgrep -m1 '^##fileformat=' "$file" || true)
    else
        header_line=$(grep -m1 '^##fileformat=' "$file" || true)
    fi

    if [[ -n "$header_line" ]]; then
    echo "Found VCF fileformat header: $header_line"

    # parse declared version (e.g. "VCFv5.2" -> "5.2")
    declared_vcf_version=$(echo "$header_line" | grep -Eo 'VCFv[0-9]+(\.[0-9]+)?' | sed 's/VCFv//')
    if [[ -n "$declared_vcf_version" ]]; then
        declared_major=$(echo "$declared_vcf_version" | cut -d. -f1)
        declared_minor=$(echo "$declared_vcf_version" | cut -s -d. -f2)
        declared_minor=${declared_minor:-0}
        echo "Declared VCF version: v${declared_major}.${declared_minor}"

        # conservative rule: treat any VCF major version >= 5 as unsupported
        if (( declared_major >= 5 )); then
            if [[ "$SKIP_CHECKS" == true ]]; then
                echo "Warning: Declared VCF version v${declared_major}.${declared_minor} may be newer than your installed bcftools (v${BCFVER_RAW}); proceeding because --skip-checks was set."
            else
                echo "Error: Detected VCF version v${declared_major}.${declared_minor}, which may not be supported by your installed bcftools (v${BCFVER_RAW})."
                echo "If you understand the risks and wish to proceed, re-run with -sc or --skip-checks to bypass this check."
                echo "SpliceAI-splint cannot guarantee correct parsing of unrecognized or future VCF versions."
                exit 1
            fi
        fi
    else
        echo "Warning: Could not parse VCF version number from header: $header_line"
    fi
     # finished header-present compatibility checks — continue the script
        return 0
    fi

    
    # -------------------------
    # If header missing
    # -------------------------
    if [[ "$SKIP_CHECKS" != true ]]; then
        echo "Warning: VCF does not specify fileformat (e.g. '##fileformat=VCFv4.2')."
        echo "Ensure VCF version is compatible with installed bcftools, or re-create the VCF with a valid header."
        echo "To auto-add a compatible fileformat header, re-run with -sc or --skip-checks."
        echo "Error: Missing '##fileformat' header. Exiting."
        exit 1
    fi

    # -------------------------
    # If skip-checks enabled, silently add compatible header
    # -------------------------
    if ! desired_line=$(pick_compatible_vcf_version); then
        echo "Error: Could not find a VCF version that your installed bcftools accepts."
        echo "Please ensure bcftools >= 1.10 is installed, or manually add a '##fileformat=' header."
        exit 1
    fi

    echo "Auto-adding compatible header: $desired_line"

    if [[ "$file" =~ \.vcf$ ]]; then
        tmp="$(mktemp)"
        { printf "%s\n" "$desired_line"; cat "$file"; } > "$tmp"
        mv "$tmp" "$file"
    elif [[ "$file" =~ \.vcf\.gz$ ]]; then
        tmp="$(mktemp)"
        zcat "$file" > "${tmp}.orig"
        { printf "%s\n" "$desired_line"; cat "${tmp}.orig"; } > "$tmp"
        bgzip -c "$tmp" > "${file}.new"
        mv "${file}.new" "$file"
        if [[ -f "${file}.tbi" ]]; then rm -f "${file}.tbi"; fi
        tabix -p vcf "$file"
        rm -f "${tmp}" "${tmp}.orig"
    else
        echo "Unrecognized VCF filename: $file"
        exit 1
    fi
    echo "VCF header updated with: $desired_line"
}

# Run fileformat header check now (before we compress .vcf -> .vcf.gz if needed)
check_fileformat_header "$VCF_FILE"

# -------------------------
# Check or skip GRCh38 validation
# -------------------------
if [[ -n "$GENOME_VERSION" ]]; then
    echo "User input genome version = $GENOME_VERSION."
else
    # choose grep vs zgrep based on filename
    if [[ "$VCF_FILE" =~ \.vcf\.gz$ ]]; then
        if ! zgrep -q "GRCh38" "$VCF_FILE"; then
            echo "Error: Input VCF does not declare reference genome as GRCh38 in the header. SpliceAI-splint can only be run on GRCh38 VCFs."
            echo "If this is incorrect, rerun with '--genome-version GRCh38' to override this check."
            exit 1
        fi
    else
        if ! grep -q "GRCh38" "$VCF_FILE"; then
            echo "Error: Input VCF does not declare reference genome as GRCh38 in the header. SpliceAI-splint can only be run on GRCh38 VCFs."
            echo "If this is incorrect, rerun with '--genome-version GRCh38' to override this check."
            exit 1
        fi
    fi
    echo "Detected genome version = GRCh38 (from VCF header)."
fi

# -------------------------
# If input is plain .vcf compress it (bgzip + tabix)
# -------------------------
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

# -------------------------
# Check chromosome naming style
# -------------------------
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

# initialize annotated files array
ANNOTATED_FILES=()

# -------------------------
# add indel 
# -------------------------
if [[ ! -f "$INDEL_ANNOTATION_FILE" ]]; then
    echo "Error: Annotation file '$INDEL_ANNOTATION_FILE' not found!"
    exit 1
fi
echo "Annotating large indels..."
bcftools filter -i '(ILEN < -4 | ILEN > 1)' "$VCF_FILE" | \
    bcftools annotate -a "$INDEL_ANNOTATION_FILE" -h "$HEADER_FILE" \
    -c CHROM,FROM,TO,SPLICEAI_RECOMMENDATION,REASON \
    -Oz -o "$OUTPUT_LARGE_INDEL"
tabix -p vcf "$OUTPUT_LARGE_INDEL"
ANNOTATED_FILES+=("$OUTPUT_LARGE_INDEL")

# -------------------------
# add wrong ref
# -------------------------
if [[ ! -f "$WRONGREF_ANNOTATION_FILE" ]]; then
    echo "Error: Annotation file '$WRONGREF_ANNOTATION_FILE' not found!"
    exit 1
fi
echo "Annotating variants with wrong reference allele..."
bcftools filter -e '(ILEN < -4 | ILEN > 1)' "$VCF_FILE" | \
    bcftools annotate --force -a "$WRONGREF_ANNOTATION_FILE" -h "$HEADER_FILE"  \
    -c CHROM,FROM,TO,SPLICEAI_RECOMMENDATION,REASON \
    -Oz -o "$OUTPUT_WRONGREF"
tabix -p vcf "$OUTPUT_WRONGREF"

# -------------------------
# add transcript changes
# -------------------------
if [[ ! -f "$SMALLVAR_ANNOTATION_FILE" ]]; then
    echo "Error: Annotation file '$SMALLVAR_ANNOTATION_FILE' not found!"
    exit 1
fi

echo "Annotating variants if in transcript change regions..."
bcftools annotate --force -a "$SMALLVAR_ANNOTATION_FILE" -h "$HEADER_FILE" \
    -c CHROM,FROM,TO,SPLICEAI_RECOMMENDATION,REASON \
    -Oz -o "$OUTPUT_SMALL_VARS" "$OUTPUT_WRONGREF"
tabix -p vcf "$OUTPUT_SMALL_VARS"

# -------------------------
# add distance
# -------------------------
if [[ ! -f "$DISTANCE_ANNOTATION_FILE" ]]; then
    echo "Error: Annotation file '$DISTANCE_ANNOTATION_FILE' not found!"
    exit 1
fi
echo "Annotating variants if we recommend running with increased distance parameter..."
bcftools annotate --force -a "$DISTANCE_ANNOTATION_FILE" -h "$HEADER_FILE" \
    -c CHROM,FROM,TO,SPLICEAI_RECOMMENDATION,REASON \
    -Oz -o "$OUTPUT_DISTANCE" "$OUTPUT_SMALL_VARS"
tabix -p vcf "$OUTPUT_DISTANCE"

ANNOTATED_FILES+=("$OUTPUT_DISTANCE")

# -------------------------
# Combine and sort the final output VCF
# -------------------------
echo "Creating and indexing the final annotated VCF..."
echo "${ANNOTATED_FILES[@]}"
bcftools concat -a -O z "${ANNOTATED_FILES[@]}" | \
    bcftools sort -Oz -o "$FINAL_OUTPUT"

tabix -p vcf "$FINAL_OUTPUT"

# -------------------------
# Clean up intermediate files
# -------------------------
echo "Cleaning up intermediate files..."
rm -f "${ANNOTATED_FILES[@]}" "${ANNOTATED_FILES[@]/%.vcf.gz/.vcf.gz.tbi}"
rm -f "$OUTPUT_WRONGREF" "$OUTPUT_WRONGREF.tbi" "$OUTPUT_SMALL_VARS" "$OUTPUT_SMALL_VARS.tbi"
if [[ "$TEMP_GZIPPED" == true ]]; then
    echo "Cleaning up temporary gzipped file..."
    rm -f "$VCF_FILE" "$VCF_FILE.tbi"
fi

echo "Annotation complete. Final output saved to $FINAL_OUTPUT"






