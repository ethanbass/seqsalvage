put_u32 <- function(x, at, v) {
  x[at + 1:4] <- writeBin(as.integer(v), raw(), size = 4, endian = "big")
  x
}
put_raw <- function(x, at, v) {
  x[at + seq_along(v)] <- v
  x
}
record <- function(tag) as.raw(c(0xc0, 0x31, 0, 0, 0, 1, tag))

contig_project <- function() {
  x <- raw(1300)
  x <- put_raw(x, 100, record(charToRaw("cCrm")))
  x <- put_u32(put_u32(x, 144, 900), 148, 7)
  x <- put_raw(x, 900, c(charToRaw("r1.ab1"), as.raw(0)))
  x <- put_raw(x, 300, record(raw(4)))
  x <- put_u32(x, 379, 7)
  x <- put_raw(x, 400, record(charToRaw("cBAS")))
  x <- put_raw(x, 500, record(raw(4)))
  x <- put_u32(put_u32(put_u32(x, 579, 9), 705, 920), 709, 4)
  x <- put_raw(x, 920, c(charToRaw("ref"), as.raw(0)))
  x <- put_raw(x, 600, record(raw(4)))
  x <- put_raw(x, 628, as.raw(2))
  x <- put_u32(put_u32(x, 649, 1000), 653, 3)
  x <- put_u32(put_u32(x, 671, 1100), 675, 6)
  x <- put_u32(put_u32(x, 805, 940), 809, 5)
  x <- put_raw(x, 940, c(charToRaw("ctg1"), as.raw(0)))
  x <- put_raw(x, 1000, as.raw(c(0xd2, 0x30, 0, 0, 0, 1, 0, 3)))
  x <- put_u32(put_u32(put_u32(x, 1008, 5), 1012, 8), 1016, 0x7777)
  x <- put_raw(x, 1100, as.raw(c(0x01, 0x12, 0x04, 0x00, 0x08, 0x02)))
  x <- put_raw(x, 1200, as.raw(c(0x92, 0x30, 0, 0, 0, 1, 0, 5)))
  for (i in 0:4) x <- put_u32(put_u32(x, 1208 + 8 * i, c(100, 300, 400, 500, 600)[i + 1]), 1212 + 8 * i, c(7, 5, 9, 8, 11)[i + 1])
  x
}

read_contigs_from <- function(x) {
  f <- tempfile(fileext = ".SPF")
  on.exit(unlink(f))
  writeBin(x, f)
  read_spf_contigs(f)
}

expect_members <- function(contigs) {
  expect_length(contigs, 1)
  expect_equal(contigs[[1]]$reads, "r1.ab1")
  expect_equal(contigs[[1]]$references, "ref")
  expect_equal(contigs[[1]]$unresolved, 1L)
}

test_that("contigs are read with their members and consensus", {
  f <- tempfile(fileext = ".SPF")
  writeBin(contig_project(), f)
  contigs <- read_spf_contigs(f)
  expect_members(contigs)
  expect_equal(contigs[[1]]$name, "ctg1")
  expect_named(contigs, "ctg1")
  expect_equal(contigs[[1]]$consensus_gapped, "ACG-TC")
  expect_equal(contigs[[1]]$consensus, "ACGTC")
  fa <- spf_contigs_to_fasta(f, tempfile(fileext = ".fasta"))
  expect_equal(readLines(fa), c(">ctg1", "ACGTC"))
})

test_that("member lists and object tables are found with later version bytes", {
  x <- put_raw(contig_project(), 1000, as.raw(c(0xd0, 0x32, 0, 0, 0, 1)))
  x <- put_raw(x, 1200, as.raw(c(0x90, 0x32, 0, 0, 0, 1)))
  expect_members(read_contigs_from(x))
})

test_that("paged member lists are followed to their pages", {
  x <- put_raw(contig_project(), 1000, as.raw(c(0x50, 0x32, 0, 0, 0, 1, 0, 1)))
  x <- put_u32(x, 1008, 1260)
  x <- put_raw(x, 1260, as.raw(c(0x90, 0x32, 0, 0, 0, 1, 0, 3)))
  x <- put_u32(put_u32(put_u32(x, 1268, 5), 1272, 8), 1276, 0x7777)
  expect_members(read_contigs_from(x))
})

test_that("members resolve without an object table", {
  expect_warning(contigs <- read_contigs_from(put_raw(contig_project(), 1200, raw(6))), "no object table")
  expect_members(contigs)
  expect_identical(contigs[[1]]$deleted, NA)
})

test_that("contigs missing from the object table are returned only on request", {
  x <- put_raw(contig_project(), 1206, as.raw(c(0, 4)))
  expect_length(read_contigs_from(x), 0)
  f <- tempfile(fileext = ".SPF")
  writeBin(x, f)
  contigs <- read_spf_contigs(f, deleted = TRUE)
  expect_members(contigs)
  expect_true(contigs[[1]]$deleted)
})

test_that("records of other types are not read as contigs", {
  expect_length(read_contigs_from(put_raw(contig_project(), 628, as.raw(0))), 0)
})

test_that("names are decoded from Mac Roman", {
  x <- put_raw(contig_project(), 940, c(charToRaw("ctg"), as.raw(c(0xc6, 0))))
  expect_equal(read_contigs_from(x)[[1]]$name, "ctg\u2206")
})
