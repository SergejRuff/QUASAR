rm(list = ls())

library(hexSticker)
library(magick)
library(devtools)
library(usethis)

# ---------------------------------------------------------------------------
# 1. Clean the source artwork
#    Crops a 12% inset off each edge to drop the stray red dot (top left) and
#    the handwritten signature (bottom right), then trims residual transparency
#    so the star itself defines the bounding box.
# ---------------------------------------------------------------------------
src <- magick::image_read("data-raw/quasar_logo.png")
info <- magick::image_info(src)

inset <- 0.12
geom <- sprintf(
  "%dx%d+%d+%d",
  round(info$width  * (1 - 2 * inset)),
  round(info$height * (1 - 2 * inset)),
  round(info$width  * inset),
  round(info$height * inset)
)

src |>
  magick::image_crop(geom) |>
  magick::image_trim() |>
  magick::image_write("data-raw/quasar_logo_clean.png")

# Inspect before continuing — if the star is clipped, lower `inset`;
# if the artefacts survive, raise it.
magick::image_read("data-raw/quasar_logo_clean.png")

# ---------------------------------------------------------------------------
# 2. Build the sticker
# ---------------------------------------------------------------------------
hexSticker::sticker(
  "data-raw/quasar_logo_clean.png",
  package  = "QUASAR",
  p_color  = "transparent",     # wordmark lives in the artwork
  p_size   = 18,
  p_y      = 0.5,
  p_family = "Aller_Rg",        # bundled with hexSticker; avoids font warnings
  h_size   = 1.5,
  h_color  = "#663399",
  h_fill   = "white",
  s_x      = 1,
  s_y      = 1.0,
  s_width  = 0.62,              # 0.8 let the crosshairs run past the border
  s_height = 0.62,
  url      = "https://github.com/SergejRuff/QUASAR",
  u_size   = 2,
  u_y      = 0.08,
  u_color  = "#663399",
  filename = "tools/quasar_logo.png"
)

# ---------------------------------------------------------------------------
# 3. Scaled copy (186px = 150 x 1.24)
# ---------------------------------------------------------------------------
magick::image_read("tools/quasar_logo.png") |>
  magick::image_scale("186") |>
  magick::image_write("tools/logo.png", format = "png")

# ---------------------------------------------------------------------------
# 4. Install into man/figures/ so the README reference resolves
# ---------------------------------------------------------------------------
usethis::use_logo("tools/quasar_logo.png")