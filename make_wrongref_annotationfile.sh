#!/bin/bash

# remove gene
# add last two columns
cut -f1-3 REF_issues/allele_check_results_20250423_221108_merged_padding_5.bed | awk -v OFS='\t' '{print $1, $2, $3, "update spliceAI scores", "wrong REF"}' > tool/tool_annotation_files/wrongref_annotation_nochr.txt
# make no chr and chr version
cut -f1-3 REF_issues/allele_check_results_20250423_221108_merged_padding_5.bed | \
awk -v OFS='\t' '{print "chr"$1, $2, $3, "update spliceAI scores", "wrong REF"}'  > tool/tool_annotation_files/wrongref_annotation_chr.txt

# sort (already done)
sort -k1,1 -k2,2n tool/tool_annotation_files/wrongref_annotation_chr.txt > tool/tool_annotation_files/wrongref_annotation_chr.sorted.txt
sort -k1,1 -k2,2n tool/tool_annotation_files/wrongref_annotation_nochr.txt > tool/tool_annotation_files/wrongref_annotation_nochr.sorted.txt
# bgzip
bgzip tool/tool_annotation_files/wrongref_annotation_nochr.sorted.txt
bgzip tool/tool_annotation_files/wrongref_annotation_chr.sorted.txt

# index
tabix -p bed -f tool/tool_annotation_files/wrongref_annotation_nochr.sorted.txt.gz
tabix -p bed -f tool/tool_annotation_files/wrongref_annotation_chr.sorted.txt.gz

# remove intermediate files

rm tool/tool_annotation_files/wrongref_annotation_nochr.txt
rm tool/tool_annotation_files/wrongref_annotation_chr.txt
