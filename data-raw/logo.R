pkgload::load_all()

src <- "tests/testthat/fixtures/S8-trimmed.spf.xz"
tmp <- tempfile(fileext = ".SPF")
writeBin(memDecompress(readBin(src, "raw", file.size(src)), "xz"), tmp)
read <- read_spf(tmp, deleted = TRUE)[[2]]

window <- 1060:1200
tr <- read$traces[window, ]
scale <- max(tr)
df <- data.frame(
  scan = rep(window, 4),
  signal = as.vector(tr) / scale,
  base = rep(colnames(tr), each = length(window))
)
called <- which((read$peaks + 1) %in% window[3:(length(window) - 3)])
calls <- data.frame(
  scan = read$peaks[called] + 1,
  base = strsplit(read$bases, "")[[1]][called]
)
menlo_bold <- file.path(tempdir(), "menlo-bold.ttf")
system2("python3", c("-I", "-c", shQuote(sprintf(
  "from fontTools.ttLib import TTCollection; TTCollection('/System/Library/Fonts/Menlo.ttc').fonts[1].save('%s')",
  menlo_bold
))))
sysfonts::font_add("menlo_bold", menlo_bold)

cols <- c(A = "#3ddc84", C = "#4dabf7", G = "#ffd43b", T = "#ff5c5c")

p <- ggplot2::ggplot(df, ggplot2::aes(scan, signal, colour = base)) +
  ggplot2::geom_line(linewidth = 0.56) +
  ggplot2::geom_text(
    ggplot2::aes(y = 1.15, label = base), data = calls,
    size = 7, fontface = "bold", family = "mono"
  ) +
  ggplot2::scale_colour_manual(values = cols, guide = "none") +
  ggplot2::coord_cartesian(ylim = c(0, 1.25), clip = "off") +
  ggplot2::theme_void()

hexSticker::sticker(
  p,
  package = "seqsalvage", p_size = 15, p_y = 0.58, p_color = "white", p_family = "menlo_bold",
  s_x = 1, s_y = 1.08, s_width = 1.35, s_height = 0.8,
  h_fill = "#252a37", h_color = "#000000", h_size = 2,
  filename = "man/figures/logo.png", dpi = 300
)
