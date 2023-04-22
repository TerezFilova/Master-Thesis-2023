# Description:
# Estimation of yield curve for given time period using Svensson model. Function
# returns vector of mean squared error, an interactive 3D plot and creates two 
# .csv files which contain estimated coefficients and estimated values for each 
# time point. Optionally function returns dataframe of estimated yield values 
# for selected maturities only. This can serve as a basis for further 
# multivariate time series modeling.
#
# Arguments:
# tl - character vector, timeline of the input data;
# mat - integer vector, maturities (in years) of the input data;
# y - yield matrix of dimensions length(tl) x length(mat), input data for
#     Svensson model series; typically market historical data
# x_by - integer, by what value is sequenced the maturity axis, default value is
#        0.05;
# T_min - integer, initial maturity of the resulting yield curve figure, default 
#         value is 1;
# T_max - integer, final maturity of the resulting yield curve figure, default 
#         value is 30;
# pd -  integer, affects the plot density, frequency of plotting the estimated
#       curves in the resulting 3d graph, value 1 means every estimated curve is
#       being plotted, default value is 63;
# mts - integer vector or null (default value), if not null function returns
#       dataframe of estimated yield values for maturities specified in this
#       argument.
#
# R packages used:
# YieldCurve, xts, plotly, progress.
#
# Author(s):
# Terezia Filova, Gabor Szucs.
#
# References:
# o Filova, T.: Matematické modelovanie vynosovych kriviek, bakalarska praca,
#   FMFI UK, Bratislava, 2021.
# o Filova, T.: Matematické modelovanie zavislosti vynosovych kriviek pomocou
#   kopula funkcii, diplomova praca, FMFI UK, Bratislava, 2023.
# o Guirreri, S. S.: Modelling and Estimation of the Yield Curve, R package,
#   (2022), <https://CRAN.R-project.org/package=YieldCurve>.
################################################################################
yc_3d_estim <- function(tl, mat, y, T_min = 1, T_max = 30, x_by = 0.05, 
                        pd = 63, mts = NULL) {
  # logical checks  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (length(mat) != ncol(y)) {
    stop("Dimension mismatch - arguments y & mat.")}
  if (max(mat) > T_max) {
    warning("The resulting yield curve is not going to be extrapolated.")}
  if (max(y) > 1) {
    warning("Svensson's model needs input to be decimal - please enter e.g. 20% 
            as 0.2.")}
  x <- seq(T_min, T_max, by = x_by) # x axis of resulting yield curve
  if (sum(is.na(match(mts, x))) >= 1) {
    stop("Each mts component must be a multiple of x_by and in T_min and T_max 
          boundary.")}
  # local variables - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  t_i <- length(tl) # time length
  m <- length(mat)
  coeff_memo <- matrix(NA, nrow = t_i, ncol = 6) # memory of Svensson parameters
    rownames(coeff_memo) <- tl
    colnames(coeff_memo) <- c("beta_0", "beta_1", "beta_2", "beta_3",
                            "alpha_1", "alpha_2")
  yc_estim <- matrix(NA, nrow = t_i, ncol = length(x)) # estimated yield curve values
    rownames(yc_estim) <- tl
    colnames(yc_estim) <- x
  yc_plot <- matrix(NA, nrow = floor(t_i / pd), ncol = length(x)); j <- 1
  yc_plot_time <- c()
  r_tilde <- c()
  mse <- matrix(NA, nrow = t_i, ncol = 1) # memory of mean squared errors
    rownames(mse) <- tl
  pb <- progress::progress_bar$new(total = t_i) # progress bar initialization
  # yield curves estimation for given time points - - - - - - - - - - - - - - - 
  for (i in 1:t_i) {
    coeff_i <- YieldCurve::Svensson(rate = y[i, ], maturity = mat)
    yc_i <- YieldCurve::Srates(xts(coeff_i, order.by = Sys.Date()),
                   maturity = x, whichRate = "Spot")
    coeff_memo[i, ] <- coeff_i; yc_estim[i, ] <- yc_i
    if (i %% pd == 0) { 
      yc_plot[j, ] <- yc_i; yc_plot_time[j] <- tl[i]
      j <- j + 1
    } 
    # mean squared error - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    if (min(mat) >= T_min && max(mat) <= T_max) {
      for (l in 1:m) {r_tilde[l] <- yc_estim[i, which(x == mat[l])]}
    } else {
      # extrapolated values - beginning
      m_lower <- 0
      if  (min(mat) < T_min) {
        x_plus <- mat[mat < T_min]; m_lower <- length(x_plus);
        for (ll in 1:m_lower) {
          r_tilde[ll] <- YieldCurve::Srates(xts(coeff_i, order.by = Sys.Date()),
                          maturity = x_plus[ll], whichRate = "Spot")}
      }
      # interpolated values
      m_ok <- sum(mat <= T_max) - m_lower
      for (l in 1:m_ok) {
        r_tilde <- c(r_tilde, yc_estim[i, which(x == mat[l])])}
      # extrapolated values - tail
      if (max(mat) > T_max) {
        x_plus <- mat[mat > T_max]; m_higher <- length(x_plus)
        for (ll in 1:m_higher) { 
          r_tilde <- c(r_tilde, 
                       YieldCurve::Srates(xts(coeff_i, order.by = Sys.Date()),
                        maturity = x_plus[ll], whichRate = "Spot"))}
      }
    }
    mse[i] <- (1 / m) * (sum(r_tilde - y[i, ])^2); length(r_tilde)
    pb$tick()
  }
  # 3D plot - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fig_yc <- plotly::plot_ly(x = x, y = yc_plot_time, z = yc_plot)
  fig_yc <- fig_yc %>% add_surface() %>% layout(scene = list(
                xaxis = list(title = "Maturity"),            
                yaxis = list(title = ""),
                zaxis = list(title = "Yield")))
  # results of estimation - - - - - - - - - - - - - - - - - - - - - - - - - - -
  write.csv(coeff_memo, file = "Coefficients.csv", row.names = TRUE)
  write.csv(yc_estim, file = "YieldCurves.csv", row.names = TRUE)
  # output for multivariate time series - - - - - - - - - - - - - - - - - - - -
  if (is.null(mts) == FALSE) {
    n <- length(mts)
    mtsb <- matrix(NA, nrow = t_i, ncol = n)
    for (i in 1:n) {mtsb[, i] <- yc_estim[, which(x == mts[i])]}
    rownames(mtsb) <- tl; colnames(mtsb) <- mts
    return(list("YC_selection" = as.data.frame(mtsb), 
                "MSE" = mse, 
                "Plot 3d" = fig_yc))
  } else {return(list("MSE" = mse, "Plot 3d" = fig_yc))}
}