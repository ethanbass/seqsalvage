# Read the trace files embedded in a Sequencher project

Extracts every chromatogram stored in a Sequencher project file
(`.SPF`), with the base calls, quality scores and peak positions that go
with it.

## Usage

``` r
read_spf(path, deleted = FALSE)
```

## Arguments

- path:

  Path to a Sequencher project file.

- deleted:

  If `TRUE`, also return reads deleted from the project whose data is
  still in the file.

## Value

A list with one element per read, named after the reads. Where names
repeat, ignoring case, later copies get a suffix (`_2`, `_3`, ...);
reads that are not deleted are numbered first, so their names are the
same whatever `deleted` is. These are the names
[`spf_to_abif()`](https://ethanbass.github.io/seqsalvage/reference/spf_to_abif.md)
gives its files and
[`read_spf_contigs()`](https://ethanbass.github.io/seqsalvage/reference/read_spf_contigs.md)
gives contig members. Each element is a list with

- `name`: the file name the read was imported under, as shown in the
  project. Names can repeat: a file imported twice gives two reads with
  one name. A deleted read whose name has been overwritten is named
  `unnamed_read_<n>`, where `n` is its position among all reads in the
  file;

- `bases`: base calls as IUPAC codes, with `-` for a gap;

- `quality`: Phred quality scores, one per base, or an empty vector if
  the project stores none for the read;

- `peaks`: the scan position of each base call, counted from 0 as in the
  ABIF `PLOC` tags, so base `i` sits at row `peaks[i] + 1` of `traces`;

- `traces`: an integer matrix of processed trace intensities with
  columns `A`, `C`, `G` and `T`, one row per scan;

- `deleted`: `TRUE` for a read deleted from the project, `FALSE`
  otherwise, and `NA` if the file has no object table to tell them
  apart.

## Details

The traces are the processed (analysed) channels, the data an `.ab1`
file holds in its `DATA9`–`DATA12` tags. The sequencer's software has
already subtracted the baseline, separated the dye signals and corrected
for mobility, so these are the curves a trace viewer shows. The raw
detector signal (`DATA1`–`DATA4`) is not stored in the project.

The bases, quality scores and peak positions are the ones stored in the
project, so they include any edits made in Sequencher. The project keeps
no copy of the original calls and no flag marking which reads were
edited, so edited and unedited reads cannot be told apart. In a project
where 150 reads could be checked against their original files, two had
been edited, each by deleting one base.

Sequencher can leave the data of deleted reads in the file. A read
counts as deleted when no object table in the file lists it;
[spf_format](https://ethanbass.github.io/seqsalvage/reference/spf_format.md)
describes the tables. If the file has none, every read is returned, with
a warning.

The `.SPF` format is undocumented.
[spf_format](https://ethanbass.github.io/seqsalvage/reference/spf_format.md)
describes how it was decoded and checked, and which parts are not read
yet.

## Examples

``` r
spf <- "project.SPF"
if (file.exists(spf)) {
  reads <- read_spf(spf)
  names(reads)
  matplot(reads[[1]]$traces[1:500, ], type = "l", lty = 1)
}
```
