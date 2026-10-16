spf_string <- function(x, at) {
  p <- spf_u32(x, at)
  l <- spf_u32(x, at + 4L)
  if (p <= 0 || l < 2 || l > 1000 || p + l > length(x)) return(NA_character_)
  b <- x[(p + 1):(p + l - 1)]
  if (any(b == as.raw(0))) return(NA_character_)
  spf_text(b)
}

spf_u16 <- function(x, at) as.integer(x[at]) * 256L + as.integer(x[at + 1L])

spf_is_record <- function(x, o) {
  !is.na(o) && o >= 0 && o + 6 <= length(x) && x[o + 1] >= as.raw(0x80) && x[o + 1] <= as.raw(0xcf) &&
    x[o + 2] %in% as.raw(c(0x30, 0x31, 0x32)) && all(x[(o + 3):(o + 6)] == as.raw(c(0, 0, 0, 1)))
}

spf_blocks <- function(x, classes) {
  hits <- lapply(classes, function(a) {
    lapply(0x30:0x32, function(v) grepRaw(as.raw(c(a, v, 0, 0, 0, 1)), x, fixed = TRUE, all = TRUE))
  })
  sort(unlist(hits)) - 1L
}

spf_block_values <- function(x, o, n) {
  if (!n || o + 8L + 4L * n > length(x)) return(integer(0))
  readBin(x[(o + 9L):(o + 8L + 4L * n)], "integer", n = n, size = 4, endian = "big")
}

spf_object_table <- function(x) {
  v <- unlist(lapply(spf_blocks(x, c(0x90, 0x92)), function(o) {
    p <- spf_block_values(x, o, 2L * spf_u16(x, o + 7L))
    if (length(p) && all(vapply(p[c(TRUE, FALSE)], function(q) spf_is_record(x, q), TRUE))) p
  }))
  if (is.null(v)) return(stats::setNames(integer(0), character(0)))
  offs <- v[c(TRUE, FALSE)]
  ids <- v[c(FALSE, TRUE)]
  keep <- !duplicated(ids)
  stats::setNames(offs[keep], ids[keep])
}

spf_members <- function(x, lp) {
  if (!(x[lp + 2L] %in% as.raw(0x30:0x32)) || !all(x[lp + 3:6] == as.raw(c(0, 0, 0, 1)))) return(NULL)
  cls <- bitwAnd(as.integer(x[lp + 1L]), 0xfd)
  n <- spf_u16(x, lp + 7L)
  if (cls == 0xd0) return(spf_block_values(x, lp, n))
  if (cls == 0x50) {
    pages <- spf_block_values(x, lp, n)
    return(unlist(lapply(pages, function(p) spf_block_values(x, p, spf_u16(x, p + 7L)))))
  }
  NULL
}

spf_object <- function(x, table, id) {
  o <- table[as.character(id)]
  if (!is.na(o)) return(unname(o))
  hits <- grepRaw(writeBin(as.integer(id), raw(), size = 4, endian = "big"), x, fixed = TRUE, all = TRUE)
  for (h in hits[hits > 4]) {
    o <- spf_u32(x, h - 4L)
    if (spf_is_record(x, o)) return(o)
  }
  NA_real_
}

#' Read the contigs in a Sequencher project
#'
#' Lists every contig in a Sequencher project file (`.SPF`) with its name,
#' consensus sequence and the reads assembled into it.
#'
#' @param path Path to a Sequencher project file.
#' @param deleted If `TRUE`, also return contigs deleted from the project whose
#'   records are still in the file.
#'
#' @return A list with one element per contig, named after the contigs. Each
#'   element is a list with
#'   * `name`: the contig's name in the project;
#'   * `consensus`: the saved consensus sequence as IUPAC codes, without gaps;
#'   * `consensus_gapped`: the same consensus with `-` for alignment columns
#'     where it has no base;
#'   * `reads`: the names of the reads in the contig, matching `names()` of the
#'     list [read_spf()] returns, so each copy of a read imported twice has its
#'     own name;
#'   * `references`: the names of any reference sequences in the contig;
#'   * `unresolved`: the number of members that could not be traced to a read
#'     or reference, for example reads deleted from the project;
#'   * `deleted`: whether the contig was deleted from the project, as for reads
#'     in [read_spf()].
#'
#' @details
#' The consensus is the one Sequencher saved with the project, so it reflects
#' the reads and edits at the time the project was last saved; it is not
#' recomputed from the reads. [spf_format] describes how contigs are stored.
#' The position and orientation of each read within its contig are not read yet.
#'
#' @seealso [read_spf()], [spf_format]
#' @examples
#' spf <- "project.SPF"
#' if (file.exists(spf)) {
#'   contigs <- read_spf_contigs(spf)
#'   names(contigs)
#'   lengths(lapply(contigs, `[[`, "reads"))
#' }
#' @export
read_spf_contigs <- function(path, deleted = FALSE) {
  x <- readBin(path, "raw", file.size(path))
  table <- spf_object_table(x)
  read_hits <- spf_read_hits(x)
  read_keys <- spf_read_keys(vapply(seq_along(read_hits), function(i) spf_read_name(x, read_hits[i], i), ""),
                             spf_deleted(x, read_hits - 1L, table, warn = FALSE))
  read_names <- stats::setNames(read_keys, read_hits - 1L)
  starts <- sort(unlist(lapply(c(0x31, 0x32), function(v) {
    grepRaw(as.raw(c(0xc0, v, 0, 0, 0, 1, 0, 0, 0, 0)), x, fixed = TRUE, all = TRUE)
  })))
  starts <- starts[x[starts + 28L] == as.raw(2)]
  del <- spf_deleted(x, starts - 1L, table)
  keep <- if (deleted) seq_along(starts) else which(!del %in% TRUE)
  contigs <- lapply(keep, function(i) {
    h <- starts[i]
    lp <- spf_u32(x, h + 49L)
    if (lp <= 0 || lp + 8 > length(x)) return(NULL)
    ids <- spf_members(x, lp)
    if (is.null(ids)) return(NULL)
    reads <- references <- character(0)
    unresolved <- 0L
    for (id in ids) {
      o <- spf_object(x, table, id)
      if (is.na(o)) {
        unresolved <- unresolved + 1L
        next
      }
      target <- spf_object(x, table, spf_u32(x, o + 80L))
      if (!is.na(target) && as.character(target) %in% names(read_names)) {
        reads <- c(reads, read_names[[as.character(target)]])
      } else if (!is.na(target) && any(vapply(c("cBAS", "BASE"), function(s) identical(x[target + 7:10], charToRaw(s)), TRUE))) {
        ref <- spf_string(x, o + 206L)
        references <- c(references, if (is.na(ref)) "reference" else ref)
      } else {
        unresolved <- unresolved + 1L
      }
    }
    cp <- spf_u32(x, h + 71L)
    cl <- spf_u32(x, h + 75L)
    gapped <- paste(iupac_codes[bitwAnd(as.integer(x[cp + seq_len(cl)]), 15L) + 1L], collapse = "")
    list(name = spf_string(x, h + 205L), consensus = gsub("-", "", gapped, fixed = TRUE),
         consensus_gapped = gapped, reads = reads, references = references, unresolved = unresolved,
         deleted = del[i])
  })
  contigs <- Filter(Negate(is.null), contigs)
  stats::setNames(contigs, vapply(contigs, `[[`, "", "name"))
}

#' Export the contig consensus sequences in a Sequencher project as FASTA
#'
#' Writes one FASTA record per contig, named after the contig, holding its
#' saved consensus sequence without gaps.
#'
#' @param path Path to a Sequencher project file.
#' @param file Path of the FASTA file to write.
#'
#' @return `file`, invisibly.
#'
#' @details
#' The consensus is the one stored when the project was last saved, so it
#' includes edits made up to then and none made later. Ambiguity codes are
#' written as stored. The gapped consensus is not written, because each contig
#' is gapped against its own alignment columns and the sequences do not form a
#' shared alignment.
#'
#' @seealso [read_spf_contigs()], [spf_to_abif()]
#' @examples
#' spf <- "project.SPF"
#' if (file.exists(spf)) spf_contigs_to_fasta(spf, tempfile(fileext = ".fasta"))
#' @export
spf_contigs_to_fasta <- function(path, file) {
  contigs <- read_spf_contigs(path)
  writeLines(paste0(">", vapply(contigs, `[[`, "", "name"), "\n", vapply(contigs, `[[`, "", "consensus")), file)
  invisible(file)
}
