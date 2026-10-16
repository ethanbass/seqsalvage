# Write a read as an ABIF (.ab1) trace file

Writes a minimal ABIF file with the four trace channels
(`DATA9`–`DATA12`, in the order given by `FWO_1`), the sample name
(`SMPL1`), and the base calls, quality scores and peak positions.

## Usage

``` r
write_abif(read, file, basecaller = FALSE)
```

## Arguments

- read:

  A list with elements `name`, `bases`, `quality`, `peaks` and `traces`,
  as returned by
  [`read_spf()`](https://ethanbass.github.io/seqsalvage/reference/read_spf.md).

- file:

  Path of the file to write.

- basecaller:

  If `TRUE`, also write the base calls, quality scores and peak
  positions to the base-caller tags `PBAS2`, `PCON2` and `PLOC2`.

## Value

`file`, invisibly.

## Details

ABIF keeps two sets of base calls: `PBAS1`, `PCON1` and `PLOC1` hold the
calls, quality scores and peak positions as edited by the user, and
`PBAS2`, `PCON2` and `PLOC2` hold them as called by the base caller. A
Sequencher project keeps only the edited set, and it does not record
which reads were edited, so by default `write_abif()` fills only the
edited tags.

Some readers need the base-caller tags. Biopython's `SeqIO.read()`
raises an error on a file without `PBAS2`, and
[`sangerseqR::readsangerseq()`](https://rdrr.io/pkg/sangerseqR/man/readsangerseq.html)
fails without `PBAS2` and `PLOC2`. Set `basecaller = TRUE` to copy the
values into those tags as well, so that these tools can read the file.
The copies are then labelled as base-caller output even where the read
was edited in Sequencher.
