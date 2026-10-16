figshare <- function(id) paste0("https://ndownloader.figshare.com/files/", id)
zenodo <- function(file) paste0("https://zenodo.org/records/5557681/files/", file, "?download=1")

published_spf <- data.frame(
  file = c(sprintf("S%d Sequencher File.SPF", 1:8),
           paste0("Dataset_S3-", c("AcroCR", "ASNA", "DOPR", "FZD", "PMCA", "TDH"), ".SPF")),
  url = c(figshare(c(6105210, 6105219, 6105225, 6105234, 6105249, 6105264, 6105267, 6105276)),
          zenodo(paste0("Dataset_S3-", c("AcroCR", "ASNA", "DOPR", "FZD", "PMCA", "TDH"), ".SPF"))),
  md5 = c("8d5d3da32212842461b57e3a694f9f39", "34c185d669f11f5ee367315c62a5c31e",
          "27bd10325cc0be14773610a702e26248", "13fc94d06208670ba9acc4f620872e15",
          "3b8118640f1132c7c7059db56ce9d45c", "738a084f92536836f9df25ac53d60aa5",
          "871488cefbc424ff6e0f74ca108bbb3a", "5e911ec130cfe873bed524ca473c0cc2",
          "499fc39a52e71eb802751f641d582a95", "bbbeb27bd2b8c2eac8444bc75a836816",
          "300be778d3052dc07fcc2209a25f820f", "48cd2b3de6c98ce55a4f55b6f0c47d65",
          "56a83a3ee14c7a2e47bac197096b71ad", "0362563b57118eb2652390f7f19970a0"),
  reads = c(191, 324, 152, 133, 265, 175, 169, 137, 80, 108, 99, 222, 72, 86),
  deleted = c(0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0),
  no_quality = c(0, 0, 0, 0, 0, 0, 0, 59, 72, 72, 69, 170, 72, 72),
  contigs = c(2, 2, 1, 2, 2, 1, 1, 1, 36, 36, 36, 36, 48, 36),
  placed = c(191, 324, 152, 133, 265, 175, 169, 137, 80, 108, 99, 222, 72, 86)
)

published_path <- function(file, dir) {
  row <- published_spf[published_spf$file == file, ]
  path <- file.path(dir, file)
  if (!file.exists(path)) stop(path, " not found; download it from ", row$url, call. = FALSE)
  if (unname(tools::md5sum(path)) != row$md5) stop(path, " does not match its published checksum", call. = FALSE)
  path
}

summarise_read <- function(r) {
  list(name = r$name, bases = r$bases, quality = r$quality, peaks = r$peaks,
       trace_dim = dim(r$traces), trace_head = r$traces[1:200, ], trace_sums = colSums(r$traces),
       deleted = r$deleted)
}

read_trimmed <- function(fixture, ...) {
  src <- test_path("fixtures", fixture)
  tmp <- tempfile(fileext = ".SPF")
  on.exit(unlink(tmp))
  writeBin(memDecompress(readBin(src, "raw", file.size(src)), "xz"), tmp)
  read_spf(tmp, ...)
}
