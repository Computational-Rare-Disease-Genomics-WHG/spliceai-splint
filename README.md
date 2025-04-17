# SpliceAI-splint
This is a repository of scripts, documents and data associated with the spliceAI-splint tool, and splice reanalysis work performed in the Genomics England 100k genomes research environment.
<br />
## The 'test_vcfs' folder : 
This folder contains two, minimal test .vcf files, one 'toy' dataset, and one containing the variants identified in our assocated analysis. <br />
These data may be used as an initial test to ensure the tool is running as expected. <br />
<br />
## The 'tool_annotation_files' folder : 
This folder contains transcript annotation coordinates used in the calculation of spliceAI precomputed scores, and new data generated using MANE v1.0.<br />
Scripts used to generate the required files for updated transcript annotations can be found here: 
https://github.com/RubyDawes/spliceai_custom_annotation <br />
<br />
## The 'spliceAI-splint' script : 
This command line script allows users to rapidly annotate large numbers of variants at a fraction of the computational cost of running SpliceAI directly. <br />
The input is any .vcf 4.0 file, which may or may not be compressed <br />
Output will be a compressed .vcf with recommendations added to the INFO column. <br />
This script requires prior installations of tabic and BCFtools: <br />
<br />
tabix: http://www.htslib.org/doc/tabix.html <br />
bcftools: https://samtools.github.io/bcftools/howtos/query.html <br />
<br />
