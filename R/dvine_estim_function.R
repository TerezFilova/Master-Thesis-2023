# Description:
# Automated estimation of D-vine multivariate copula structure. Function returns
# list of selected bivariate copulas from estimated D-vine and their parameters.
# Objects in the list are from BiCop class. For better orientation in the list, 
# function also prints the D-vine tree using copula names from the list.
#
# Arguments:
# u - matrix of pseudo observations, columns ordered hierarchically;
# indt - logical, whether a hypothesis test for the independence of u1 and u2 is
#        performed before bivariate copula selection, TRUE is default 
#
# R packages used:
# VineCopula.
# 
# Author(s):
# Terezia Filova, Gabor Szucs.
#
# References:
# o Filova, T.: Matematické modelovanie zavislosti vynosovych kriviek pomocou 
#   kopula funkcii, diplomova praca, FMFI UK, Bratislava, 2023.
# o Nagler, T., Schepsmeier, U., Stoeber, J., Brechmann, E., Graeler, B., 
#   Erhardt, T.: Statistical Inference of Vine Copulas, R package, (2022)
#   <https://CRAN.R-project.org/package=VineCopula>.
# o Righi, M. B., Schlender, S. G.,  Ceretta, P. S.,: Pair copula constructions 
#   to determine the dependence structure of Treasury bond yields, IIMB Manage-
#   ment Review (2015) 27, p. 216–227.
################################################################################
dvine_estim <- function(u, indt = TRUE) {
  # local variables 
  L <- ncol(u) - 1 # dimension of the vine
  col_nam <- c()
  if (is.null(colnames(u))) {
    col_nam <- c(col_nam, paste0("u", 1:ncol(u)))
    colnames(u) <- col_nam
  } else {col_nam <- colnames(u)}
  attach(as.data.frame(u))
  var_nam <- matrix(NA, ncol = L, nrow = L)
  cdf_memory <- matrix(NA, ncol = L + 1, nrow = L)
  cdf_memory[1, ] <- col_nam
  #
  # MAIN LOOP FOR D-VINE ESTIMATION ============================================
  for (i in 1:L) {
    print(paste("Calculating: D-vine - Level", i))
    # level 1 of the vine
    if (i == 1) {
      for (j in 1:L) {
        value <- VineCopula::BiCopSelect(u[, j], u[, j + 1], indeptest = indt)
        var_nam[i, j] <- paste0(col_nam[j], ".", col_nam[j + 1])
        assign(var_nam[i, j], value)
      }
    } else { # other levels of the vine
      obj1 <- get(var_nam[i - 1, 1])
      arg1 <- get(cdf_memory[i - 1, 1])
      arg2 <- get(cdf_memory[i - 1, 2])
      value <- VineCopula::BiCopHfunc2(u1 = arg1, u2 = arg2, obj = obj1)
      cdf_memory[i, 1] <- paste0("F", col_nam[1],
                                 "_", paste0(col_nam[2:i], collapse = "."))
      assign(cdf_memory[i, 1], value)
      for (j in 1:(L + 1 - i)) {
        # getting values from variable names
        obj2 <- get(var_nam[i - 1, j + 1])
        arg2 <- get(cdf_memory[i - 1, j + 1])
        arg3 <- get(cdf_memory[i - 1, j + 2])
        # constructing conditional distribution
        cdf1 <- get(cdf_memory[i, j])
        cdf2 <- VineCopula::BiCopHfunc1(u1 = arg2, u2 = arg3, obj = obj2)
        cdf_memory[i, j + 1] <- paste0("F", col_nam[j + i],
                                       "_", paste0(col_nam[(j + 1):(j + i - 1)],
                                                   collapse = "."))
        assign(cdf_memory[i, j + 1], cdf2)
        # conditional copula selection
        value <- VineCopula::BiCopSelect(cdf1, cdf2, indeptest = indt)
        var_nam[i, j] <- paste0(col_nam[j], ".", col_nam[j + i],
                                "_", paste0(col_nam[(j + 1):(j + i - 1)],
                                            collapse = "."))
        assign(var_nam[i, j], value)
      }
    }
  } # ==========================================================================
  detach()
  # results of the estimation
  copula_nam <- as.vector(t(var_nam)); copula_nam <- copula_nam[!is.na(copula_nam)]
  bicop_list <- c()
  for (i in 1:length(copula_nam)) {
    bicop_list <- rbind(get(copula_nam[i]), bicop_list)
  }
  rownames(bicop_list) <- rev(copula_nam)
  print("D-vine")
  print(var_nam)
  print("Bivariate Copulas")
  return(bicop_list)
}