# Export every read in a Sequencher project as an .ab1 file

Export every read in a Sequencher project as an .ab1 file

## Usage

``` r
spf_to_abif(path, out_dir, basecaller = FALSE, deleted = FALSE)
```

## Arguments

- path:

  Path to a Sequencher project file.

- out_dir:

  Directory to write the `.ab1` files to. It is created if needed.

- basecaller:

  If `TRUE`, also write the base calls, quality scores and peak
  positions to the base-caller tags `PBAS2`, `PCON2` and `PLOC2`.

- deleted:

  If `TRUE`, also return reads deleted from the project whose data is
  still in the file.

## Value

The paths of the files written, invisibly. Each file takes the read's
name in
[`read_spf()`](https://ethanbass.github.io/seqsalvage/reference/read_spf.md),
so a read whose name repeats another gets a suffix (`_2`, `_3`, ...) and
a warning, and no file overwrites another.

## See also

[`read_spf()`](https://ethanbass.github.io/seqsalvage/reference/read_spf.md),
[`write_abif()`](https://ethanbass.github.io/seqsalvage/reference/write_abif.md),
[spf_format](https://ethanbass.github.io/seqsalvage/reference/spf_format.md)

## Examples

``` r
spf <- "project.SPF"
if (file.exists(spf)) spf_to_abif(spf, tempdir())
```
