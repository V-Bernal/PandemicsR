.libPaths(c(".r-lib", .libPaths()))

if (!requireNamespace("rsconnect", quietly = TRUE)) {
  stop("Install rsconnect into .r-lib before deploying.", call. = FALSE)
}

dir.create(".renv-cache", showWarnings = FALSE)
Sys.setenv(RENV_PATHS_CACHE = normalizePath(".renv-cache"))

account_name <- Sys.getenv("SHINYAPPS_NAME", unset = "")
account_token <- Sys.getenv("SHINYAPPS_TOKEN", unset = "")
account_secret <- Sys.getenv("SHINYAPPS_SECRET", unset = "")

if (nzchar(account_name) && nzchar(account_token) && nzchar(account_secret)) {
  rsconnect::setAccountInfo(
    name = account_name,
    token = account_token,
    secret = account_secret
  )
} else {
  registered_accounts <- rsconnect::accounts()
  shinyapps_accounts <- registered_accounts[registered_accounts$server == "shinyapps.io", , drop = FALSE]
  if (nzchar(account_name)) {
    shinyapps_accounts <- shinyapps_accounts[shinyapps_accounts$name == account_name, , drop = FALSE]
  }
  if (nrow(shinyapps_accounts) != 1L) {
    stop(
      "Set SHINYAPPS_NAME, SHINYAPPS_TOKEN, and SHINYAPPS_SECRET, or register exactly one shinyapps.io account locally.",
      call. = FALSE
    )
  }
  account_name <- shinyapps_accounts$name[[1]]
}

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
  account = account_name,
  server = "shinyapps.io",
  appVisibility = "public",
  logLevel = "verbose",
  forceUpdate = TRUE
)
