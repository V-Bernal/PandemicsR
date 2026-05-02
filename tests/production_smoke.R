local({
  script_file <- tryCatch(sys.frame(1)$ofile, error = function(err) NULL)
  repo_root <- if (!is.null(script_file)) {
    normalizePath(file.path(dirname(script_file), ".."), mustWork = FALSE)
  } else {
    getwd()
  }
  if (!file.exists(file.path(repo_root, "app.R"))) {
    repo_root <- getwd()
  }
  setwd(repo_root)

  .libPaths(c(".r-lib", .libPaths()))
  library(Matrix)
  library(igraph)

  for (r_file in sort(list.files("PandemicsR-main/R", pattern = "[.]R$", full.names = TRUE))) {
    source(r_file)
  }

  assert <- function(condition, message) {
    if (!isTRUE(condition)) {
      stop(message, call. = FALSE)
    }
  }

  assert_close <- function(x, y, tolerance = 1e-6, message = "Values differ") {
    if (any(abs(x - y) > tolerance, na.rm = TRUE)) {
      stop(message, call. = FALSE)
    }
  }

  assert_result_integrity <- function(result, expected_mode) {
    assert(identical(result$simulation_mode, expected_mode), "Unexpected simulation mode.")
    assert(length(result$time_history) == nrow(result$frac_mat), "Opinion history/time mismatch.")
    assert(length(result$time_history) == nrow(result$sir_mat), "SIR history/time mismatch.")
    assert(length(result$camp_sir_history) == length(result$time_history), "Camp SIR history/time mismatch.")
    assert(!any(result$sir_mat < -1e-6), "Negative SIR counts detected.")
    assert(!any(result$frac_mat < -1e-6 | result$frac_mat > 1 + 1e-6), "Opinion fractions outside [0, 1].")
    assert_close(rowSums(result$sir_mat), rep(result$population_size, nrow(result$sir_mat)), 1e-4, "Population not conserved in SIR history.")
    assert_close(rowSums(result$frac_mat), rep(1, nrow(result$frac_mat)), 1e-6, "Opinion fractions do not sum to one.")
    assert(all(is.na(result$camp_attack_rate) | (result$camp_attack_rate >= -1e-6 & result$camp_attack_rate <= 1 + 1e-6)), "Attack rates outside [0, 1].")
    assert_close(sum(result$final_camp_sir), result$population_size, 1e-4, "Final camp SIR counts do not sum to population.")
    assert(all(diff(result$time_history) >= -1e-9), "Time history is not monotone.")
  }

  render_core_plots <- function(result) {
    png(tempfile(fileext = ".png"), width = 900, height = 600)
    visual_step_time(result$frac_mat, result$num_opinions, result$time_history)
    dev.off()

    png(tempfile(fileext = ".png"), width = 900, height = 600)
    if (identical(result$simulation_mode, "aggregate")) {
      visual_opinion_shares(result$frac_mat, result$num_opinions)
    } else {
      visual_histo(result$opinion_history, result$num_opinions)
    }
    dev.off()

    png(tempfile(fileext = ".png"), width = 900, height = 600)
    visual_sir_time(result$sir_mat, result$time_history)
    dev.off()

    png(tempfile(fileext = ".png"), width = 900, height = 800)
    visual_sir_camp_time(result$camp_sir_history, result$time_history)
    dev.off()

    png(tempfile(fileext = ".png"), width = 900, height = 600)
    visual_cumulative_infections(result$infection_events, result$final_time, result$camp_labels)
    dev.off()

    png(tempfile(fileext = ".png"), width = 900, height = 600)
    visual_epidemic_camps(result$final_camp_sir)
    dev.off()
  }

  default_args <- list(
    m = 3,
    t_max = 100,
    lambda = 0.5,
    c_param = 0.4,
    gamma_light = 0.1,
    gamma_dark = 0.02,
    infected_dark_multiplier = 1.5,
    beta_plus = 0.5,
    beta_minus = 0.2,
    T_threshold = 0.3,
    num_opinions = 4,
    beta_red = 1.6,
    beta_blue = 0.05,
    gamma_sir = 0.30,
    initial_infected_fraction = 0.05
  )
  with_args <- function(...) utils::modifyList(default_args, list(...))

  set.seed(42)
  exact_result <- do.call(simulate_hybrid_model, with_args(n = 80, t_max = 30))
  assert_result_integrity(exact_result, "exact")
  assert(isTRUE(exact_result$graph_available), "Small exact run should keep graph outputs available.")
  assert(exact_result$final_time <= 30 + 1e-9, "Exact run exceeded the requested time horizon.")
  render_core_plots(exact_result)

  aggregate_timing <- system.time({
    aggregate_result <- do.call(simulate_hybrid_model, with_args(n = 200000, m = 20))
  })
  assert_result_integrity(aggregate_result, "aggregate")
  assert(!isTRUE(aggregate_result$graph_available), "Aggregate run should not materialize graph outputs.")
  assert(aggregate_timing[["elapsed"]] < 20, "Aggregate 200k run exceeded 20 seconds.")
  render_core_plots(aggregate_result)

  default_result <- do.call(simulate_hybrid_model, with_args(n = 120))
  assert_result_integrity(default_result, "exact")
  assert(isTRUE(default_result$graph_available), "Default run should keep graph outputs available.")
  assert(default_result$camp_attack_rate[["Red camp"]] > default_result$camp_attack_rate[["Blue camp"]], "Default parameters should make the red camp sicker.")

  no_epidemic_result <- simulate_hybrid_model(
    n = 40,
    m = 4,
    t_max = 20,
    lambda = 0,
    c_param = 0,
    gamma_light = 0,
    gamma_dark = 0,
    infected_dark_multiplier = 1,
    beta_plus = 0,
    beta_minus = 0,
    T_threshold = 0,
    num_opinions = 2,
    run_voter = FALSE,
    run_schelling = FALSE,
    beta_red = 0,
    beta_blue = 0,
    gamma_sir = 0,
    initial_infected_fraction = 0
  )
  assert_result_integrity(no_epidemic_result, "exact")
  assert_close(no_epidemic_result$overall_attack_rate, 0, 1e-9, "Zero-epidemic run should have zero attack rate.")
  render_core_plots(no_epidemic_result)

  invalid_result <- tryCatch(
    simulate_hybrid_model(
      n = 10,
      m = 11,
      t_max = 1,
      lambda = 0,
      c_param = 0,
      gamma_light = 0,
      gamma_dark = 0,
      infected_dark_multiplier = 1,
      beta_plus = 0,
      beta_minus = 0,
      T_threshold = 0,
      num_opinions = 2,
      beta_red = 0,
      beta_blue = 0,
      gamma_sir = 0,
      initial_infected_fraction = 0
    ),
    error = identity
  )
  assert(inherits(invalid_result, "error"), "Invalid m > n parameters should fail fast.")

  app_env <- new.env(parent = globalenv())
  source("app.R", local = app_env)
  assert(exists("app", envir = app_env), "app.R did not create a Shiny app object.")
  if (requireNamespace("htmltools", quietly = TRUE)) {
    rendered_shell <- htmltools::renderTags(app_env$ui)
    rendered_dashboard <- htmltools::renderTags(app_env$dashboard_ui)
    assert(grepl("<title>PandemicsR Sandbox</title>", rendered_shell$head, fixed = TRUE), "Outer page browser title was not set cleanly.")
    assert(grepl("<title>PandemicsR Sandbox</title>", rendered_dashboard$head, fixed = TRUE), "Browser title was not set cleanly.")
    assert(!grepl("<title>.*h2", rendered_dashboard$head), "Browser title appears to contain an h2 tag.")
  }

  message("production_smoke.R: all checks passed")
})
