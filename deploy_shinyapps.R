.libPaths(c(".r-lib", .libPaths()))

if (!requireNamespace("rsconnect", quietly = TRUE)) {
  stop("Install rsconnect into .r-lib before deploying.", call. = FALSE)
}

dir.create(".renv-cache", showWarnings = FALSE)
Sys.setenv(RENV_PATHS_CACHE = normalizePath(".renv-cache"))

required_vars <- c("SHINYAPPS_NAME", "SHINYAPPS_TOKEN", "SHINYAPPS_SECRET")
missing_vars <- required_vars[!nzchar(Sys.getenv(required_vars, unset = ""))]

if (length(missing_vars) > 0) {
  stop(
    sprintf(
      "Set these environment variables before deploying: %s",
      paste(missing_vars, collapse = ", ")
    ),
    call. = FALSE
  )
}

rsconnect::setAccountInfo(
  name = Sys.getenv("SHINYAPPS_NAME"),
  token = Sys.getenv("SHINYAPPS_TOKEN"),
  secret = Sys.getenv("SHINYAPPS_SECRET")
)

app_files <- c(
  "app.R",
  "shinny.R",
  sort(list.files("PandemicsR-main/R", pattern = "\\.R$", full.names = TRUE))
)

rsconnect::deployApp(
  appDir = ".",
  appFiles = app_files,
  appPrimaryDoc = "app.R",
  appName = Sys.getenv("SHINYAPPS_APP_NAME", unset = "pandemicsr"),
  appVisibility = "public",
  logLevel = "verbose",
  forceUpdate = TRUE
)
