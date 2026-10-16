expected <- readRDS(test_path("fixtures", "expected.rds"))

test_that("trimmed projects read the same as the published originals", {
  for (fixture in names(expected)) {
    expect_equal(lapply(read_trimmed(fixture, deleted = TRUE), summarise_read), expected[[fixture]], label = fixture)
  }
})

test_that("reads are named after their names in the project", {
  reads <- read_trimmed("S8-trimmed.spf.xz", deleted = TRUE)
  expect_equal(names(reads), vapply(reads, `[[`, "", "name", USE.NAMES = FALSE))
})

test_that("records with either tag byte are found", {
  expect_length(read_trimmed("S1-trimmed.spf.xz"), 1)
  expect_length(read_trimmed("S8-trimmed.spf.xz"), 2)
})

test_that("a read without quality scores returns an empty quality vector", {
  reads <- read_trimmed("S8-trimmed.spf.xz")
  expect_identical(reads[["MC-1A-2-M13F(-47).ab1"]]$quality, integer(0))
  expect_length(reads[["MC-1A-37-M13F(-47)_R.ab1"]]$quality, nchar(reads[["MC-1A-37-M13F(-47)_R.ab1"]]$bases))
})

test_that("deleted reads are returned only on request", {
  expect_false(any(vapply(read_trimmed("S8-trimmed.spf.xz"), `[[`, NA, "deleted")))
  reads <- read_trimmed("S8-trimmed.spf.xz", deleted = TRUE)
  deleted <- vapply(reads, `[[`, NA, "deleted")
  expect_equal(sum(deleted), 1)
  expect_equal(reads[deleted][[1]]$name, "MC-1A-14-M13F(-47).ab1")
})

test_that("published projects read consistently", {
  skip_on_cran()
  skip_on_ci()
  dir <- Sys.getenv("SEQSALVAGE_TEST_DATA")
  skip_if(dir == "", "set SEQSALVAGE_TEST_DATA to a directory holding the published projects")
  present <- file.exists(file.path(dir, published_spf$file))
  skip_if(!any(present), "no published projects found in SEQSALVAGE_TEST_DATA")
  for (i in which(present)) {
    p <- published_spf[i, ]
    path <- published_path(p$file, dir)
    all_reads <- read_spf(path, deleted = TRUE)
    reads <- all_reads[!vapply(all_reads, `[[`, NA, "deleted")]
    expect_length(all_reads, p$reads + p$deleted)
    expect_length(reads, p$reads)
    n_qual <- vapply(reads, function(r) length(r$quality), 1L)
    n_base <- vapply(reads, function(r) nchar(r$bases), 1L)
    expect_equal(sum(n_qual == 0), p$no_quality, label = p$file)
    expect_true(all(n_qual == 0 | n_qual == n_base), label = p$file)
    expect_true(all(vapply(reads, function(r) {
      length(r$peaks) == nchar(r$bases) && all(diff(r$peaks) >= 0) && max(r$peaks) < nrow(r$traces)
    }, TRUE)), label = p$file)
    contigs <- read_spf_contigs(path)
    expect_length(contigs, p$contigs)
    members <- unlist(lapply(contigs, `[[`, "reads"))
    expect_equal(length(members), p$placed, label = p$file)
    expect_true(all(members %in% names(reads)) && !anyDuplicated(members), label = p$file)
    expect_false(anyDuplicated(tolower(names(reads))) > 0, label = p$file)
    expect_equal(sum(vapply(contigs, `[[`, 1L, "unresolved")), 0, label = p$file)
  }
})
