#' Write a read as an ABIF (.ab1) trace file
#'
#' Writes a minimal ABIF file with the four trace channels (`DATA9`–`DATA12`,
#' in the order given by `FWO_1`), the sample name (`SMPL1`), and the base
#' calls, quality scores and peak positions.
#'
#' ABIF keeps two sets of base calls: `PBAS1`, `PCON1` and `PLOC1` hold the
#' calls, quality scores and peak positions as edited by the user, and `PBAS2`,
#' `PCON2` and `PLOC2` hold them as called by the base caller. A Sequencher
#' project keeps only the edited set, and it does not record which reads were
#' edited, so by default `write_abif()` fills only the edited tags.
#'
#' Some readers need the base-caller tags. Biopython's `SeqIO.read()` raises an
#' error on a file without `PBAS2`, and [sangerseqR::readsangerseq()] fails
#' without `PBAS2` and `PLOC2`. Set `basecaller = TRUE` to copy the values into
#' those tags as well, so that these tools can read the file. The copies are
#' then labelled as base-caller output even where the read was edited in
#' Sequencher.
#'
#' @param read A list with elements `name`, `bases`, `quality`, `peaks` and
#'   `traces`, as returned by [read_spf()].
#' @param file Path of the file to write.
#' @param basecaller If `TRUE`, also write the base calls, quality scores and
#'   peak positions to the base-caller tags `PBAS2`, `PCON2` and `PLOC2`.
#'
#' @return `file`, invisibly.
#' @export
write_abif <- function(read, file, basecaller = FALSE) {
  smpl <- charToRaw(sub("\\.ab1$", "", read$name))
  calls <- list(bases = charToRaw(read$bases), qual = as.raw(pmin(read$quality, 255L)), peaks = abif_shorts(read$peaks))
  entries <- list(
    list("DATA", 9L, 4L, 2L, abif_shorts(read$traces[, "A"])),
    list("DATA", 10L, 4L, 2L, abif_shorts(read$traces[, "C"])),
    list("DATA", 11L, 4L, 2L, abif_shorts(read$traces[, "G"])),
    list("DATA", 12L, 4L, 2L, abif_shorts(read$traces[, "T"])),
    list("FWO_", 1L, 2L, 1L, charToRaw("ACGT")),
    list("SMPL", 1L, 18L, 1L, c(as.raw(length(smpl)), smpl)))
  for (n in if (basecaller) 1:2 else 1L) {
    entries <- c(entries, list(list("PBAS", n, 2L, 1L, calls$bases), list("PCON", n, 2L, 1L, calls$qual),
                               list("PLOC", n, 4L, 2L, calls$peaks)))
  }
  entries <- Filter(function(e) length(e[[5]]) > 0L, entries)
  sizes <- vapply(entries, function(e) length(e[[5]]), integer(1))
  external <- ifelse(sizes > 4L, sizes, 0L)
  offsets <- 128L + c(0L, cumsum(external))[seq_along(sizes)]
  con <- file(file, "wb")
  on.exit(close(con))
  writeChar("ABIF", con, eos = NULL)
  writeBin(101L, con, size = 2, endian = "big")
  writeChar("tdir", con, eos = NULL)
  writeBin(1L, con, size = 4, endian = "big")
  writeBin(c(1023L, 28L), con, size = 2, endian = "big")
  writeBin(c(length(entries), 28L * length(entries), 128L + sum(external), 0L), con, size = 4, endian = "big")
  writeBin(raw(128L - 34L), con)
  for (i in which(sizes > 4L)) writeBin(entries[[i]][[5]], con)
  for (i in seq_along(entries)) {
    e <- entries[[i]]
    writeChar(e[[1]], con, eos = NULL)
    writeBin(e[[2]], con, size = 4, endian = "big")
    writeBin(c(e[[3]], e[[4]]), con, size = 2, endian = "big")
    writeBin(c(sizes[i] %/% e[[4]], sizes[i]), con, size = 4, endian = "big")
    if (sizes[i] > 4L) writeBin(offsets[i], con, size = 4, endian = "big") else writeBin(c(e[[5]], raw(4L - sizes[i])), con)
    writeBin(0L, con, size = 4, endian = "big")
  }
  invisible(file)
}

abif_shorts <- function(v) {
  con <- rawConnection(raw(0), "wb")
  on.exit(close(con))
  writeBin(as.integer(v), con, size = 2, endian = "big")
  rawConnectionValue(con)
}
