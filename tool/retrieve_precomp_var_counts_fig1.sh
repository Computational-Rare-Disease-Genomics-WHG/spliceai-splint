#!/bin/bash

# Configuration
OUTPUT_FILE="annotation_summary.txt"
R_OUTPUT_FILE="annotation_summary.csv"
VCF_FILE="../spliceai_precomp/spliceai_scores.masked.snv.hg38.vcf.gz"
ANNOTATION_FILE="tool_annotation_files/precomp_supplement_annotation_nochr.sorted.txt.gz"
WRONGREF_FILE="tool_annotation_files/wrongref_annotation_nochr.sorted.txt.gz"
HEADER_FILE="tool_annotation_files/precomp_supplement_annotation.hdr"

# Create output files
{
    echo "ANNOTATION SUMMARY REPORT"
    echo "========================"
    echo "Generated: $(date)"
    echo ""
    echo "REASON                          COUNT"
    echo "------                          -----"
} > "$OUTPUT_FILE"

echo "reason,count,category,percentage" > "$R_OUTPUT_FILE"

echo "Starting annotation analysis..."

# 1. Get total variants in original file
echo "Counting total variants in original file..."
total_variants=$(bcftools view -H "$VCF_FILE" | wc -l)
echo "TOTAL_ORIGINAL                  $total_variants" >> "$OUTPUT_FILE"
echo "TOTAL_ORIGINAL,$total_variants,total,100.00" >> "$R_OUTPUT_FILE"

# 2. Count categories from first annotation (single pass)
echo "Counting small variants annotation categories..."

bcftools annotate --force \
    -a "$ANNOTATION_FILE" \
    -h "$HEADER_FILE" \
    "$VCF_FILE" \
    -c CHROM,FROM,TO,SPLICEAI_RECOMMENDATION,REASON | \
bcftools query -f '%INFO/REASON\n' | \
awk '
{
    if ($0 == "missing gene") missing_gene++
    else if ($0 == "transcript boundary change") transcript_boundary++
    else if ($0 == "missing exon") missing_exon++
    else if ($0 == ".") small_vars_missing++
}
END {
    print missing_gene+0, transcript_boundary+0, missing_exon+0, small_vars_missing+0
}' > temp_counts1.txt

# Read the counts
read missing_gene transcript_boundary missing_exon small_vars_missing < temp_counts1.txt

# Calculate percentages for first annotation
missing_gene_pct=$(echo "scale=2; $missing_gene * 100 / $total_variants" | bc -l)
transcript_boundary_pct=$(echo "scale=2; $transcript_boundary * 100 / $total_variants" | bc -l)
missing_exon_pct=$(echo "scale=2; $missing_exon * 100 / $total_variants" | bc -l)
small_vars_missing_pct=$(echo "scale=2; $small_vars_missing * 100 / $total_variants" | bc -l)

# Write to output files
echo "missing gene                    $missing_gene" >> "$OUTPUT_FILE"
echo "transcript boundary change      $transcript_boundary" >> "$OUTPUT_FILE"
echo "missing exon                    $missing_exon" >> "$OUTPUT_FILE"
echo "small_vars_missing (.)          $small_vars_missing" >> "$OUTPUT_FILE"

echo "missing_gene,$missing_gene,annotated,$missing_gene_pct" >> "$R_OUTPUT_FILE"
echo "transcript_boundary_change,$transcript_boundary,annotated,$transcript_boundary_pct" >> "$R_OUTPUT_FILE"
echo "missing_exon,$missing_exon,annotated,$missing_exon_pct" >> "$R_OUTPUT_FILE"
echo "small_vars_missing,$small_vars_missing,missing,$small_vars_missing_pct" >> "$R_OUTPUT_FILE"

# 3. Count categories from second annotation (single pass)
echo "Counting wrong REF annotation categories..."

bcftools annotate --force \
    -a "$WRONGREF_FILE" \
    -h "$HEADER_FILE" \
    "$VCF_FILE" \
    -c CHROM,FROM,TO,SPLICEAI_RECOMMENDATION,REASON | \
bcftools query -f '%INFO/REASON\n' | \
awk '
{
    if ($0 == "wrong REF") wrong_ref++
    else if ($0 == ".") wrong_ref_missing++
}
END {
    print wrong_ref+0, wrong_ref_missing+0
}' > temp_counts2.txt

# Read the counts
read wrong_ref wrong_ref_missing < temp_counts2.txt

# Calculate percentages for second annotation
wrong_ref_pct=$(echo "scale=2; $wrong_ref * 100 / $total_variants" | bc -l)
wrong_ref_missing_pct=$(echo "scale=2; $wrong_ref_missing * 100 / $total_variants" | bc -l)

# Write to output files
echo "wrong REF                       $wrong_ref" >> "$OUTPUT_FILE"
echo "wrong_ref_missing (.)           $wrong_ref_missing" >> "$OUTPUT_FILE"

echo "wrong_REF,$wrong_ref,annotated,$wrong_ref_pct" >> "$R_OUTPUT_FILE"
echo "wrong_ref_missing,$wrong_ref_missing,missing,$wrong_ref_missing_pct" >> "$R_OUTPUT_FILE"

# 4. Add summary calculations
{
    echo ""
    echo "SUMMARY CALCULATIONS"
    echo "===================="
    
    small_vars_annotated=$((missing_gene + transcript_boundary + missing_exon))
    wrong_ref_annotated=$wrong_ref
    total_annotated=$((small_vars_annotated + wrong_ref_annotated))
    total_missing=$((small_vars_missing + wrong_ref_missing))
    
    echo "Total annotated (small_vars):   $small_vars_annotated"
    echo "Total annotated (wrong_ref):    $wrong_ref_annotated"
    echo "Total annotated (combined):     $total_annotated"
    echo "Total missing annotations:      $total_missing"
    echo "Total processed:                $((total_annotated + total_missing))"
    
    # Calculate percentages
    if [[ $total_variants -gt 0 ]]; then
        annotated_pct=$(echo "scale=2; $total_annotated * 100 / $total_variants" | bc -l)
        missing_pct=$(echo "scale=2; $total_missing * 100 / $total_variants" | bc -l)
        echo "Annotation rate:                ${annotated_pct}%"
        echo "Missing rate:                   ${missing_pct}%"
    fi
    
} >> "$OUTPUT_FILE"

# Add summary rows to CSV
if [[ $total_variants -gt 0 ]]; then
    annotated_pct=$(echo "scale=2; $total_annotated * 100 / $total_variants" | bc -l)
    missing_pct=$(echo "scale=2; $total_missing * 100 / $total_variants" | bc -l)
    echo "total_annotated,$total_annotated,summary,$annotated_pct" >> "$R_OUTPUT_FILE"
    echo "total_missing,$total_missing,summary,$missing_pct" >> "$R_OUTPUT_FILE"
fi

# Clean up temporary files
rm -f temp_counts1.txt temp_counts2.txt

echo ""
echo "Analysis completed!"
echo "Human-readable results: $OUTPUT_FILE"
echo "R-friendly results: $R_OUTPUT_FILE"
echo "No intermediate VCF files created - saved ~60GB of disk space!"

# Display summary
echo ""
echo "QUICK SUMMARY:"
echo "=============="
echo "Total variants:     $total_variants"
echo "Total annotated:    $total_annotated"
echo "Total missing:      $total_missing"
if [[ $total_variants -gt 0 ]]; then
    echo "Annotation rate:    ${annotated_pct}%"
fi
