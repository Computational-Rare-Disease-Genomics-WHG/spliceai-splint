# SpliceAI-splint

This is a repository of scripts, documents, and data associated with the `SpliceAI-splint` tool and spliceAI reanalysis work performed in the Genomics England 100k Genomes Project research environment.

The `spliceAI-splint` tool provides a lightweight, command-line method to rapidly annotate variants for which the precomputed SpliceAI scores may not be adequate. This means users can re-run SpliceAI on only the subset of variants which are most likely to require updated scores.

---

## Table of Contents

- [Overview](#overview)
- [Installation Requirements](#installation-requirements)
- [Usage](#usage)
- [Annotation files](#annotation-files)
- [Cite](#cite)

---

## Overview

The `spliceAI-splint` script allows users to:
- Annotate large indels excluded from precomputed scores.
- Annotate variants in regions where there have been transcript changes since the precomputed scores were calculated
- Recommend when variants may benefit from increased spliceAI distance parameter.
- Identify variants where reference mismatches occur in the precomputed scores

Input is any .vcf 4.0 file, which may or may not be compressed. Output is a `.vcf.gz` file with additional `INFO` fields `SPLICEAI_RECOMMENDATION` and `REASON`.

---

## Installation Requirements

- [`bcftools`](https://samtools.github.io/bcftools/) (v1.10+)
- [`tabix`](http://www.htslib.org/doc/tabix.html)
- Bash (Linux/macOS environment)

---

## Usage

```bash
./spliceai_supplement.sh -i <input.vcf.gz> [-o <output.vcf.gz>]
```

The 'test_vcfs' folder contains two, minimal test .vcf files, one 'toy' dataset, and one containing the variants identified in our assocated analysis. These data may be used as an initial test to ensure the splice-splint tool is running as expected.

## Annotation files
The 'tool_annotation_files' folder contains transcript annotation coordinates used by the SpliceAI-splint tool to identify variants which require recomputed spliceAI scores. These are generated using the `create_tool_annotation_files.Rmd` script.

## What to do with output
Any variants which are annotated by spliceAI-splint with `SPLICEAI_RECOMMENDATION=update spliceAI scores` should be annotated using the SpliceAI command-line tool, available here: [`https://github.com/Illumina/SpliceAI`](https://github.com/Illumina/SpliceAI)

## Cite
If you use SpliceAI-splint, please cite [insert paper here]
