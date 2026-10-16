#' Read the trace files embedded in a Sequencher project
#'
#' Extracts every chromatogram stored in a Sequencher project file (`.SPF`),
#' with the base calls, quality scores and peak positions that go with it.
#'
#' @param path Path to a Sequencher project file.
#' @param deleted If `TRUE`, also return reads deleted from the project whose
#'   data is still in the file.
#'
#' @return A list with one element per read, named after the reads. Where
#'   names repeat, ignoring case, later copies get a suffix (`_2`, `_3`, ...);
#'   reads that are not deleted are numbered first, so their names are the same
#'   whatever `deleted` is. These are the names [spf_to_abif()] gives its files
#'   and [read_spf_contigs()] gives contig members. Each element is a list with
#'   * `name`: the file name the read was imported under, as shown in the project.
#'     Names can repeat: a file imported twice gives two reads with one name. A
#'     deleted read whose name has been overwritten is named `unnamed_read_<n>`,
#'     where `n` is its position among all reads in the file;
#'   * `bases`: base calls as IUPAC codes, with `-` for a gap;
#'   * `quality`: Phred quality scores, one per base, or an empty vector if the
#'     project stores none for the read;
#'   * `peaks`: the scan position of each base call, counted from 0 as in the
#'     ABIF `PLOC` tags, so base `i` sits at row `peaks[i] + 1` of `traces`;
#'   * `traces`: an integer matrix of processed trace intensities with columns
#'     `A`, `C`, `G` and `T`, one row per scan;
#'   * `deleted`: `TRUE` for a read deleted from the project, `FALSE` otherwise,
#'     and `NA` if the file has no object table to tell them apart.
#'
#' @details
#' The traces are the processed (analysed) channels, the data an `.ab1` file
#' holds in its `DATA9`–`DATA12` tags. The sequencer's software has already
#' subtracted the baseline, separated the dye signals and corrected for
#' mobility, so these are the curves a trace viewer shows. The raw detector
#' signal (`DATA1`–`DATA4`) is not stored in the project.
#'
#' The bases, quality scores and peak positions are the ones stored in the
#' project, so they include any edits made in Sequencher. The project keeps no
#' copy of the original calls and no flag marking which reads were edited, so
#' edited and unedited reads cannot be told apart. In a project where 150 reads
#' could be checked against their original files, two had been edited, each
#' by deleting one base.
#'
#' Sequencher can leave the data of deleted reads in the file. A read counts as
#' deleted when no object table in the file lists it; [spf_format] describes
#' the tables. If the file has none, every read is returned, with a warning.
#'
#' The `.SPF` format is undocumented. [spf_format] describes how it was decoded
#' and checked, and which parts are not read yet.
#'
#' @examples
#' spf <- "project.SPF"
#' if (file.exists(spf)) {
#'   reads <- read_spf(spf)
#'   names(reads)
#'   matplot(reads[[1]]$traces[1:500, ], type = "l", lty = 1)
#' }
#' @export
read_spf <- function(path, deleted = FALSE) {
  x <- readBin(path, "raw", file.size(path))
  hits <- spf_read_hits(x)
  if (!length(hits)) stop("no reads found in ", path, call. = FALSE)
  del <- spf_deleted(x, hits - 1L)
  keys <- spf_read_keys(vapply(seq_along(hits), function(i) spf_read_name(x, hits[i], i), ""), del)
  keep <- if (deleted) seq_along(hits) else which(!del %in% TRUE)
  reads <- lapply(keep, function(i) c(read_spf_record(x, hits[i], i), deleted = del[i]))
  stats::setNames(reads, keys[keep])
}

spf_read_keys <- function(names, deleted) {
  ab1 <- grepl("\\.ab1$", names, ignore.case = TRUE)
  stems <- ifelse(ab1, substr(names, 1L, nchar(names) - 4L), names)
  ord <- order(deleted %in% TRUE)
  n <- integer(length(names))
  n[ord] <- stats::ave(seq_along(ord), tolower(stems[ord]), FUN = seq_along)
  ifelse(n > 1L, paste0(stems, "_", n, substring(names, nchar(stems) + 1L)), names)
}

abif_file_names <- function(keys) {
  ifelse(grepl("\\.ab1$", keys, ignore.case = TRUE), keys, paste0(keys, ".ab1"))
}

#' Export every read in a Sequencher project as an .ab1 file
#'
#' @param path Path to a Sequencher project file.
#' @param out_dir Directory to write the `.ab1` files to. It is created if needed.
#' @inheritParams write_abif
#' @inheritParams read_spf
#'
#' @return The paths of the files written, invisibly. Each file takes the
#'   read's name in [read_spf()], so a read whose name repeats another gets a
#'   suffix (`_2`, `_3`, ...) and a warning, and no file overwrites another.
#'
#' @seealso [read_spf()], [write_abif()], [spf_format]
#' @examples
#' spf <- "project.SPF"
#' if (file.exists(spf)) spf_to_abif(spf, tempdir())
#' @export
spf_to_abif <- function(path, out_dir, basecaller = FALSE, deleted = FALSE) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  reads <- read_spf(path, deleted = deleted)
  renamed <- sum(names(reads) != vapply(reads, `[[`, "", "name"))
  if (renamed) warning(renamed, " reads share a name with another read; their files get a numbered suffix", call. = FALSE)
  files <- file.path(out_dir, abif_file_names(names(reads)))
  for (i in seq_along(reads)) write_abif(reads[[i]], files[i], basecaller = basecaller)
  invisible(files)
}

iupac_codes <- c("-", "A", "C", "M", "G", "R", "S", "V", "T", "W", "Y", "H", "K", "D", "B", "N")

spf_read_tags <- lapply(c(0x31, 0x32), function(v) as.raw(c(0xc0, v, 0, 0, 0, 1, 0x63, 0x43, 0x72, 0x6d)))

spf_u32 <- function(x, at) readBin(x[at:(at + 3L)], "integer", size = 4, endian = "big")

# Apple's Mac OS Roman table (ROMAN.TXT), bytes 0x80-0xFF
mac_roman <- c(
  0x00C4, 0x00C5, 0x00C7, 0x00C9, 0x00D1, 0x00D6, 0x00DC, 0x00E1,
  0x00E0, 0x00E2, 0x00E4, 0x00E3, 0x00E5, 0x00E7, 0x00E9, 0x00E8,
  0x00EA, 0x00EB, 0x00ED, 0x00EC, 0x00EE, 0x00EF, 0x00F1, 0x00F3,
  0x00F2, 0x00F4, 0x00F6, 0x00F5, 0x00FA, 0x00F9, 0x00FB, 0x00FC,
  0x2020, 0x00B0, 0x00A2, 0x00A3, 0x00A7, 0x2022, 0x00B6, 0x00DF,
  0x00AE, 0x00A9, 0x2122, 0x00B4, 0x00A8, 0x2260, 0x00C6, 0x00D8,
  0x221E, 0x00B1, 0x2264, 0x2265, 0x00A5, 0x00B5, 0x2202, 0x2211,
  0x220F, 0x03C0, 0x222B, 0x00AA, 0x00BA, 0x03A9, 0x00E6, 0x00F8,
  0x00BF, 0x00A1, 0x00AC, 0x221A, 0x0192, 0x2248, 0x2206, 0x00AB,
  0x00BB, 0x2026, 0x00A0, 0x00C0, 0x00C3, 0x00D5, 0x0152, 0x0153,
  0x2013, 0x2014, 0x201C, 0x201D, 0x2018, 0x2019, 0x00F7, 0x25CA,
  0x00FF, 0x0178, 0x2044, 0x20AC, 0x2039, 0x203A, 0xFB01, 0xFB02,
  0x2021, 0x00B7, 0x201A, 0x201E, 0x2030, 0x00C2, 0x00CA, 0x00C1,
  0x00CB, 0x00C8, 0x00CD, 0x00CE, 0x00CF, 0x00CC, 0x00D3, 0x00D4,
  0xF8FF, 0x00D2, 0x00DA, 0x00DB, 0x00D9, 0x0131, 0x02C6, 0x02DC,
  0x00AF, 0x02D8, 0x02D9, 0x02DA, 0x00B8, 0x02DD, 0x02DB, 0x02C7
)

spf_text <- function(b) {
  i <- as.integer(b)
  high <- i >= 128L
  if (!any(high)) return(rawToChar(b))
  i[high] <- mac_roman[i[high] - 127L]
  intToUtf8(i)
}

dnfx_residuals <- function(nib, n) {
  out <- integer(n)
  k <- 0L
  pos <- 1L
  eights <- which(nib == 8L)
  while (k < n) {
    nxt <- eights[findInterval(pos - 1L, eights) + 1L]
    if (is.na(nxt)) nxt <- length(nib) + 1L
    m <- min(nxt - pos, n - k)
    if (m > 0L) {
      v <- nib[pos:(pos + m - 1L)]
      out[k + seq_len(m)] <- v - 16L * (v >= 8L)
      k <- k + m
      pos <- pos + m
    }
    if (k >= n) break
    pos <- pos + 1L
    repeat {
      byte <- nib[pos] * 16L + nib[pos + 1L]
      pos <- pos + 2L
      s <- byte %% 128L
      if (s >= 64L) s <- s - 128L
      if (s == 0L) {
        s <- sum(nib[pos:(pos + 3L)] * c(4096L, 256L, 16L, 1L))
        pos <- pos + 4L
        if (s >= 32768L) s <- s - 65536L
      }
      k <- k + 1L
      out[k] <- if (s > 0L) s + 7L else s - 7L
      if (byte >= 128L || k >= n) break
    }
  }
  out
}

read_dnfx <- function(x, ptr) {
  p <- ptr + 1L
  if (rawToChar(x[p:(p + 3L)]) != "DnFx") stop("no DnFx block at offset ", ptr, call. = FALSE)
  usize <- spf_u32(x, p + 4L)
  csize <- spf_u32(x, p + 8L)
  order <- as.integer(x[p + 21L])
  pay <- as.integer(x[(p + 24L):(p + 23L + csize)])
  nchunk <- pay[1] * 256L + pay[2]
  b <- pay[(4L + 14L * nchunk + 1L):length(pay)]
  v <- as.numeric(dnfx_residuals(as.vector(rbind(b %/% 16L, b %% 16L)), usize %/% 2L))
  for (i in seq_len(order)) v <- cumsum(v)
  as.integer(v)
}

spf_read_hits <- function(x) sort(unlist(lapply(spf_read_tags, grepRaw, x, fixed = TRUE, all = TRUE)))

spf_read_name <- function(x, h, i) {
  name <- spf_string(x, h + 44L)
  if (is.na(name) || !nzchar(name) || grepl("[[:cntrl:]]", name)) paste0("unnamed_read_", i) else name
}

spf_deleted <- function(x, offsets, table = spf_object_table(x), warn = TRUE) {
  if (length(table)) return(!offsets %in% table)
  if (warn) warning("no object table found, so deleted records cannot be told apart and are all returned", call. = FALSE)
  rep(NA, length(offsets))
}

read_spf_record <- function(x, h, i) {
  tab <- vapply(0:7, function(i) c(spf_u32(x, h + 259L + 8L * i), spf_u32(x, h + 263L + 8L * i)), numeric(2))
  bases <- as.integer(x[tab[1, 1] + seq_len(tab[2, 1])])
  n_scan <- spf_u32(x, tab[1, 5] + 5L) %/% 2L
  traces <- vapply(5:8, function(i) read_dnfx(x, tab[1, i]), integer(n_scan))
  colnames(traces) <- c("A", "C", "G", "T")
  list(name = spf_read_name(x, h, i),
       bases = paste(iupac_codes[bases + 1L], collapse = ""),
       quality = as.integer(x[tab[1, 2] + seq_len(tab[2, 2])]),
       peaks = read_dnfx(x, tab[1, 4]),
       traces = traces)
}
