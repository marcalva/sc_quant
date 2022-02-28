
# SC Quant

A tool for counting reads in single-cell data. SC Quant can generate 
gene counts and variant allele counts. Gene counts are split by 
spliced, unspliced, or ambiguous.
Currently only works with UMI-based data.

## Installation

To install SC Quant, clone the github repository recursively to get the 
latest version.
```bash
git clone --recursive git@github.com:marcalva/sc_quant.git
```
Cloning recursively makes sure that htslib is downloaded in a subdirectory.

To generate the binary file, type in Make. You can then copy or move the 
binary to another directory.

## Running SC Quant

There are two modes to generate counts: gene-based and allele-based. 
Gene-based counts will assign UMIs to genes as spliced, unspliced, or 
ambiguous. Allele-based counts assign UMIs to variants and generate counts 
consistent with the reference, alternate, or other allele. To run 
gene-based counts, use `sc_quant gc`. Similarly, to run allele-based
counts, use `sc_quant ac`. Type in `sc_quant gc --help` and 
`sc_quant ac --help` to see a list of options.

The required input file for SC Quant in both modes is a BAM file. The 
SAM/BAM tag for both the barcode and UMI must be present. These are 
by default `CB` and `UB`, respectively, to use the error-corrected 
oligo sequences from Cell Ranger or STARSolo.

For allele counts, a VCF or BCF file is required. The format fields containing 
the sample genotype data are not required. However, the file must be 
properly formatted and indexed. For gene counts, a GTF file is required. The 
`gene_id` field is used to uniquely determine each gene. The `gene_id`, 
`transcript_id`, `exon_id`, `gene_name`, and `gene_type` attributes are 
required.

A list of barcodes can be provided so that counts will be output for all 
those listed. This is typically a whitelist, for example the list of valid 
barcodes from [10X](https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist-).
If not provided, only barcodes that are detected and have valid counts 
are output.

When considering multi-mapping reads by setting `--max-nh` greater than 
1, only the primary alignment is considered. 
It is recommended to keep the default of disregarding multi-mapping 
reads by setting `--max-nh 1`.

SC Quant can also be run for a specific region, for example by chromosome. 
This is useful to speed up runs in parallel.
Valid regions are those accepted by htslib, but only chromosome runs have 
been tested. If run by chromosome, the output file for gene counts will 
still include all genes in the GTF file. However, for allele counts, only 
variants in the region are output. When combining output from multiple runs, 
the allele count matrices should be concatenated and the gene count 
matrices should be added. To avoid this, the GTF file can be subsetted 
before running SC Quant.

## Counting procedure

For gene counts, only UMIs that overlap a single gene are considered. Even 
if a UMI is uniquely aligned to the genome, it can still overlap two or more 
genes. These are disregarded. Similarly, isoforms where the UMI is not 
fully contained within it are disregarded. Isoforms of a gene which fully 
contain the UMI are considered compatible.

A UMI is considered spliced if it crosses an annotated splice junction from 
a compatible isoform. It is considered intronic if at least one base pair 
overlaps an intron from all compatible transcript isoforms.
Otherwise, it is ambiguous. For example, if read overlaps an exon or UTR 
fully, it is ambiguous. If there are no compatible isoforms, the UMI is also 
ambiguous.

For allele counts, if a UMI overlaps two or more variants, each variant 
gets one count. This means that reads can be counted more than once.

## Output

The output for both gene counts and allele counts are a series of sparse 
matrices in [matrix market format](https://math.nist.gov/MatrixMarket/formats.html).
For gene counts, three matrices are produced with the suffixes
- `gc.spl.mtx.gz` for spliced counts
- `gc.uns.mtx.gz` for unspliced counts
- `gc.amb.mtx.gz` for ambiguous counts

The rows (genes) and columns (barcodes) are given in the files ending 
`gc.gene.txt.gz` and `gc.barcodes.txt.gz`, respectively.

For allele counts, three matrices are produced with the suffixes
- `ac.ref.mtx.gz` for reference allele counts
- `ac.alt.mtx.gz` for alternate allele counts
- `ac.oth.mtx.gz` for other allele counts

Similarly, the rows (variants) and columns (barcodes) are given in the files 
ending `ac.var.txt.gz` and `ac.barcodes.txt.gz`, respectively.

The files `gc.gene.txt.gz` and `ac.var.txt.gz` contain additonal information 
on the genes and variants from the GTF file and the variant file.

