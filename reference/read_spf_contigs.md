# Read the contigs in a Sequencher project

Lists every contig in a Sequencher project file (`.SPF`) with its name,
consensus sequence and the reads assembled into it.

## Usage

``` r
read_spf_contigs(path, deleted = FALSE)
```

## Arguments

- path:

  Path to a Sequencher project file.

- deleted:

  If `TRUE`, also return contigs deleted from the project whose records
  are still in the file.

## Value

A list with one element per contig, named after the contigs. Each
element is a list with

- `name`: the contig's name in the project;

- `consensus`: the saved consensus sequence as IUPAC codes, without
  gaps;

- `consensus_gapped`: the same consensus with `-` for alignment columns
  where it has no base;

- `reads`: the names of the reads in the contig, matching
  [`names()`](https://rdrr.io/r/base/names.html) of the list
  [`read_spf()`](https://ethanbass.github.io/seqsalvage/reference/read_spf.md)
  returns, so each copy of a read imported twice has its own name;

- `references`: the names of any reference sequences in the contig;

- `unresolved`: the number of members that could not be traced to a read
  or reference, for example reads deleted from the project;

- `deleted`: whether the contig was deleted from the project, as for
  reads in
  [`read_spf()`](https://ethanbass.github.io/seqsalvage/reference/read_spf.md).

## Details

The consensus is the one Sequencher saved with the project, so it
reflects the reads and edits at the time the project was last saved; it
is not recomputed from the reads.
[spf_format](https://ethanbass.github.io/seqsalvage/reference/spf_format.md)
describes how contigs are stored. The position and orientation of each
read within its contig are not read yet.

## See also

[`read_spf()`](https://ethanbass.github.io/seqsalvage/reference/read_spf.md),
[spf_format](https://ethanbass.github.io/seqsalvage/reference/spf_format.md)

## Examples

``` r
spf <- "project.SPF"
if (file.exists(spf)) {
  contigs <- read_spf_contigs(spf)
  names(contigs)
  lengths(lapply(contigs, `[[`, "reads"))
}
```
