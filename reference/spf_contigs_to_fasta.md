# Export the contig consensus sequences in a Sequencher project as FASTA

Writes one FASTA record per contig, named after the contig, holding its
saved consensus sequence without gaps.

## Usage

``` r
spf_contigs_to_fasta(path, file)
```

## Arguments

- path:

  Path to a Sequencher project file.

- file:

  Path of the FASTA file to write.

## Value

`file`, invisibly.

## Details

The consensus is the one stored when the project was last saved, so it
includes edits made up to then and none made later. Ambiguity codes are
written as stored. The gapped consensus is not written, because each
contig is gapped against its own alignment columns and the sequences do
not form a shared alignment.

## See also

[`read_spf_contigs()`](https://ethanbass.github.io/seqsalvage/reference/read_spf_contigs.md),
[`spf_to_abif()`](https://ethanbass.github.io/seqsalvage/reference/spf_to_abif.md)

## Examples

``` r
spf <- "project.SPF"
if (file.exists(spf)) spf_contigs_to_fasta(spf, tempfile(fileext = ".fasta"))
```
