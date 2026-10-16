nibbles <- function(hex) strtoi(strsplit(hex, "")[[1]], 16L)

test_that("residual codes decode", {
  expect_equal(seqsalvage:::dnfx_residuals(nibbles("3d07"), 4), c(3L, -3L, 0L, 7L))
  expect_equal(seqsalvage:::dnfx_residuals(nibbles("881"), 1), 8L)
  expect_equal(seqsalvage:::dnfx_residuals(nibbles("8ff"), 1), -8L)
  expect_equal(seqsalvage:::dnfx_residuals(nibbles("87e06fb2"), 4), c(-9L, 13L, -12L, 2L))
  expect_equal(seqsalvage:::dnfx_residuals(nibbles("800042a00f37e000c7c80fbdddca"), 7),
               c(1073L, -3209L, 3203L, -1066L, -3L, -4L, -6L))
})

test_that("a DnFx block integrates residuals by its predictor order", {
  pay <- c(0x00, 0x00, 0x01, 0x00, 0x13, 0x10)
  block <- as.raw(c(charToRaw("DnFx"), 0, 0, 0, 6, 0, 0, 0, length(pay), rep(0, 9), 2, 0, 0, pay))
  expect_equal(seqsalvage:::read_dnfx(block, 0L), c(1L, 5L, 10L))
})

abif_tags <- function(f) {
  x <- readBin(f, "raw", file.size(f))
  n <- readBin(x[19:22], "integer", size = 4, endian = "big")
  dir <- readBin(x[27:30], "integer", size = 4, endian = "big")
  list(x = x, dir = dir, tags = vapply(seq_len(n) - 1L, function(i) {
    paste0(rawToChar(x[dir + 28L * i + 1:4]), readBin(x[dir + 28L * i + 5:8], "integer", size = 4, endian = "big"))
  }, ""))
}
read <- list(name = "x.ab1", bases = "ACG", quality = c(10L, 20L, 30L), peaks = c(2L, 5L, 8L),
             traces = cbind(A = 1:10, C = 11:20, G = 21:30, T = 31:40))

test_that("write_abif writes the edited tags by default", {
  a <- abif_tags(write_abif(read, tempfile(fileext = ".ab1")))
  expect_equal(rawToChar(a$x[1:4]), "ABIF")
  expect_equal(a$tags, c("DATA9", "DATA10", "DATA11", "DATA12", "FWO_1", "SMPL1", "PBAS1", "PCON1", "PLOC1"))
  off <- readBin(a$x[a$dir + 28L * 2L + 21:24], "integer", size = 4, endian = "big")
  expect_equal(readBin(a$x[off + 1:20], "integer", n = 10, size = 2, endian = "big"), 21:30)
})

test_that("basecaller = TRUE adds the base-caller tags", {
  a <- abif_tags(write_abif(read, tempfile(fileext = ".ab1"), basecaller = TRUE))
  expect_true(all(c("PBAS2", "PCON2", "PLOC2") %in% a$tags))
})

test_that("a read without quality scores gets no PCON tags", {
  a <- abif_tags(write_abif(modifyList(read, list(quality = integer(0))), tempfile(fileext = ".ab1"), basecaller = TRUE))
  expect_false(any(c("PCON1", "PCON2") %in% a$tags))
  expect_true(all(c("PBAS1", "PBAS2", "PLOC1", "PLOC2") %in% a$tags))
})

test_that("repeated read names get numbered keys and file names", {
  keys <- spf_read_keys(c("a.ab1", "b", "a.ab1", "A.ab1", "b.ab1"), rep(FALSE, 5))
  expect_equal(keys, c("a.ab1", "b", "a_2.ab1", "A_3.ab1", "b_2.ab1"))
  expect_equal(abif_file_names(keys), c("a.ab1", "b.ab1", "a_2.ab1", "A_3.ab1", "b_2.ab1"))
  expect_equal(spf_read_keys(c("a.ab1", "b.ab1"), c(FALSE, FALSE)), c("a.ab1", "b.ab1"))
})

test_that("reads that are not deleted keep the plain name", {
  expect_equal(spf_read_keys(c("x.ab1", "x.ab1"), c(TRUE, FALSE)), c("x_2.ab1", "x.ab1"))
  expect_equal(spf_read_keys(c("x.ab1", "x.ab1"), c(NA, NA)), c("x.ab1", "x_2.ab1"))
})
