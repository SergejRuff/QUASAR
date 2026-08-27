rmarkdown::render("README.Rmd", output_format = "github_document")


# vignette
usethis::use_vignette("quasar-deconvolution", "Cell-type deconvolution with QUASAR")


# render vignette
rmarkdown::render("vignettes/quasar-deconvolution.Rmd")
#devtools::install(build_vignettes = TRUE)
vignette("quasar-deconvolution", package = "quasar")



# built site
unlink("docs", recursive = TRUE)
pkgdown::build_site()


pkgdown::preview_site()


# Create github page
#usethis::use_pkgdown_github_pages()


usethis::use_github_pages()
usethis::use_github_action("pkgdown")
usethis::use_build_ignore("docs")



dir.create("vignettes/figures", showWarnings = FALSE)

magick::image_read("data-raw/r_python_comparison.tif") |>
  magick::image_convert(format = "png") |>
  magick::image_resize("1600x") |>          # cap the width; TIFFs are often huge
  magick::image_write("vignettes/figures/r_python_comparison.png")