rm(list = ls())

library(hexSticker)
library(magick)
library(devtools)
library(usethis)

# 1. Trim transparent padding so the artwork's visible ink is what gets centred
magick::image_read("data-raw/quasar_logo.png") |>
  magick::image_trim() |>
  magick::image_write("data-raw/quasar_logo_trim.png")

# 2. Build the sticker from the TRIMMED file
hexSticker::sticker(
  "data-raw/quasar_logo_trim.png",   # <- was the untrimmed file
  package  = "QUASAR",
  p_color  = "transparent",          # wordmark lives in the artwork
  p_size   = 18,
  p_y      = 0.5,
  p_family = "Aller_Rg",             # bundled with hexSticker; no font warnings
  h_size   = 1.5,
  h_color  = "#663399",
  h_fill   = "white",
  s_x      = 1,
  s_y      = 1.0,                    # true vertical centre of the hex
  s_width  = 0.8,
  s_height = 0.8,
  url      = "https://github.com/SergejRuff/QUASAR",
  u_size   = 2,
  u_y      = 0.08,                   # 0.05 sat too close to the bottom vertex
  u_color  = "#663399",
  filename = "tools/quasar_logo.png"
)

# 3. Scaled copy for the README
magick::image_read("tools/quasar_logo.png") |>
  magick::image_scale("150") |>
  magick::image_write("tools/logo.png", format = "png")

# 4. Install into man/figures/ so the README's <img src="man/figures/logo.png"> resolves
usethis::use_logo("tools/quasar_logo.png")