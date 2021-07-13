library(renv)

Sys.setenv(RENV_DOWNLOAD_FILE_METHOD = "libcurl")
options(repos = c(CRAN = "https://cloud.r-project.org"))

renv::restore(
  # clean = TRUE
  # , 
  confirm = FALSE, 
  repos = "https://cloud.r-project.org"
)