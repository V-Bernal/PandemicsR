#' Get parameters for predefined scenario
#' @export
get_scenario_params <- function(scenario) {

  switch(
    scenario,

    resilient = list(
      gamma = 0.50,
      c_param = 0.30,
      beta_plus = 0.05,
      beta_minus = 0.01,
      T_threshold = 0.50,

      alpha0_rad = 0.01,
      alpha = 0.05,
      alpha0_derad = 0.05,
      alpha_deradicalization = 0.10,

      beta_red_red = 0.20,
      beta_red_blue = 0.10,
      beta_blue_red = 0.10,
      beta_blue_blue = 0.20,
      gamma_epi = 0.20
    ),

    polarized = list(
      gamma = 1.00,
      c_param = 0.10,
      beta_plus = 0.10,
      beta_minus = 0.005,
      T_threshold = 0.70,

      alpha0_rad = 0.02,
      alpha = 0.15,
      alpha0_derad = 0.01,
      alpha_deradicalization = 0.02,

      beta_red_red = 0.30,
      beta_red_blue = 0.05,
      beta_blue_red = 0.05,
      beta_blue_blue = 0.30,
      gamma_epi = 0.15
    ),

    radicalization = list(
      gamma = 0.80,
      c_param = 0.20,
      beta_plus = 0.08,
      beta_minus = 0.01,
      T_threshold = 0.60,

      alpha0_rad = 0.05,
      alpha = 0.30,
      alpha0_derad = 0.005,
      alpha_deradicalization = 0.01,

      beta_red_red = 0.35,
      beta_red_blue = 0.15,
      beta_blue_red = 0.15,
      beta_blue_blue = 0.35,
      gamma_epi = 0.08
    ),

    epidemic = list(
      gamma = 0.40,
      c_param = 0.50,
      beta_plus = 0.05,
      beta_minus = 0.02,
      T_threshold = 0.50,

      alpha0_rad = 0.03,
      alpha = 0.10,
      alpha0_derad = 0.02,
      alpha_deradicalization = 0.05,

      beta_red_red = 0.60,
      beta_red_blue = 0.40,
      beta_blue_red = 0.40,
      beta_blue_blue = 0.60,
      gamma_epi = 0.05
    ),

    NULL
  )
}
