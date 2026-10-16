# Builds the trimmed test fixtures in tests/testthat/fixtures from the published
# projects. Run from the package root, with SEQSALVAGE_TEST_DATA set to a
# directory holding the S1 and S8 projects (URLs in tests/testthat/helper-published.R).

pkgload::load_all()
invisible(testthat::source_test_helpers("tests/testthat", env = environment()))

spf_used_bytes <- function(x, h) {
  name_ptr <- spf_u32(x, h + 44L)
  tab <- vapply(0:7, function(i) c(spf_u32(x, h + 259L + 8L * i), spf_u32(x, h + 263L + 8L * i)), numeric(2))
  dnfx <- lapply(4:8, function(i) tab[1, i] + seq_len(24L + spf_u32(x, tab[1, i] + 9L)))
  c(h + 0:322, name_ptr + seq_len(spf_u32(x, h + 48L)),
    tab[1, 1] + seq_len(tab[2, 1]), tab[1, 2] + seq_len(tab[2, 2]), unlist(dnfx))
}

trim_spf <- function(src, keep, keep_deleted, dest) {
  x <- readBin(src, "raw", file.size(src))
  hits <- spf_read_hits(x)
  names <- vapply(seq_along(hits), function(i) spf_read_name(x, hits[i], i), "")
  del <- spf_deleted(x, hits - 1L)
  sel <- which((names %in% keep & !del) | (names %in% keep_deleted & del))
  stopifnot(length(sel) == length(keep) + length(keep_deleted))
  idx <- unlist(lapply(hits[sel], spf_used_bytes, x = x))
  for (o in spf_blocks(x, c(0x90, 0x92))) {
    p <- spf_block_values(x, o, 2L * spf_u16(x, o + 7L))
    offs <- p[c(TRUE, FALSE)]
    valid <- length(offs) && all(vapply(offs, function(q) spf_is_record(x, q), TRUE))
    if (valid && any(offs %in% (hits[sel] - 1L))) idx <- c(idx, o + seq_len(8L + 4L * length(p)), unlist(lapply(offs, `+`, 1:6)))
  }
  out <- raw(max(idx))
  out[idx] <- x[idx]
  writeBin(memCompress(out, "xz"), dest)
  lapply(read_spf(src, deleted = TRUE)[sel], summarise_read)
}

data_dir <- Sys.getenv("SEQSALVAGE_TEST_DATA")
if (data_dir == "") stop("set SEQSALVAGE_TEST_DATA to the directory holding the published projects")
fixtures <- list(
  "S1-trimmed.spf.xz" = list(file = "S1 Sequencher File.SPF", keep = "5C-est2-6-M13-FP.ab1", keep_deleted = character(0)),
  "S8-trimmed.spf.xz" = list(file = "S8 Sequencher File.SPF",
                             keep = c("MC-1A-2-M13F(-47).ab1", "MC-1A-37-M13F(-47)_R.ab1"),
                             keep_deleted = "MC-1A-14-M13F(-47).ab1")
)

expected <- lapply(names(fixtures), function(f) {
  trim_spf(published_path(fixtures[[f]]$file, data_dir), fixtures[[f]]$keep, fixtures[[f]]$keep_deleted,
           file.path("tests/testthat/fixtures", f))
})
names(expected) <- names(fixtures)
saveRDS(expected, "tests/testthat/fixtures/expected.rds", compress = "xz")
