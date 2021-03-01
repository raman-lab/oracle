These scripts are used in sequence to analyze DMS data. 

Paired-end Illumina sequencing reads are first merged with FLASH (Fast Length Adjustment of SHort reads) using the default software parameters.

Sequences are 'cleaned' by deleting the header, the third empty line, and putting the nucleotide sequence on the same line as the quality scores. 
This is done with a bash script, and example of which is clean_merged_fastq.

Poor quality reads are then removed with an R script, an example of which is shown as SampleRemoveBadReads.R

Samples are translated and binned using an internal barcode in another R script, an example of which is shown in SampleDNAtoAA_WithBarcodes.R

Single variants are then counted with SampleAnalysis.R
