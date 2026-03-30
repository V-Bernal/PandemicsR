source("shinny.R", local = TRUE)

if (!exists("app", inherits = FALSE)) {
  app <- shinyApp(ui, server)
}

app
