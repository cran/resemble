## -----------------------------------------------------------------------------
#| echo: false
options(cli.num_colors = 1)
Sys.setenv(OMP_NUM_THREADS = 2)


## -----------------------------------------------------------------------------
#| message: false
#| echo: false
library(resemble)


## -----------------------------------------------------------------------------
#| label: fit-methods
#| eval: true

# Creates an object with instructions to build PLS models
my_pls <- fit_pls(ncomp = 15)
my_pls

# Creates an object with instructions to build WAPLS models
my_wapls <- fit_wapls(min_ncomp = 3, max_ncomp = 20)
my_wapls

# Creates an object with instructions to build GPR models
my_gpr <- fit_gpr()
my_gpr


## -----------------------------------------------------------------------------
#| label: validation-control
#| eval: true
# Create an object with instructions to conduct both validation types
# "NNv" and "local_cv"
two_val_control <- mbl_control(
  validation_type = c("NNv", "local_cv"),
  number = 10,
  p = 0.75
)


## -----------------------------------------------------------------------------
#| message: false
library(resemble)
library(prospectr)

# obtain a numeric vector of the wavelengths at which spectra is recorded 
wavs <- as.numeric(colnames(NIRsoil$spc))

# pre-process the spectra:
# - use detrend
# - use first order derivative
diff_order <- 1
poly_order <- 1
window <- 7

# Preprocess spectra
NIRsoil$spc_pr <- savitzkyGolay(
  detrend(NIRsoil$spc, wav = wavs),
  m = diff_order, p = poly_order, w = window
)
train_x <- NIRsoil$spc_pr[NIRsoil$train == 1, ]
train_y <- NIRsoil$Ciso[NIRsoil$train == 1]

test_x  <- NIRsoil$spc_pr[NIRsoil$train == 0, ]
test_y  <- NIRsoil$Ciso[NIRsoil$train == 0]


## -----------------------------------------------------------------------------
#| label: mbl-local
#| eval: true
#| results: hide
# Define the neighborhood sizes to test
my_ks <- seq(80, 160, by = 40)

# Define how to use the dissimilarity information (ignore it)
ignore_diss <- "none"

# Define the regression method to be used at each neighborhood
my_wapls <- fit_wapls(min_ncomp = 3, max_ncomp = 25, scale = TRUE)

# For the moment use only "NNv" validation (faster)
nnv_val_control <- mbl_control()

# Predict Total Carbon
# (remove missing values)
local_ciso <- mbl(
  Xr = train_x[!is.na(train_y), ],
  Yr = train_y[!is.na(train_y)],
  Xu = test_x,
  neighbors = neighbors_k(my_ks),
  diss_method = diss_correlation(center = FALSE),
  diss_usage = ignore_diss,
  fit_method = my_wapls,
  gh = TRUE,
  control = nnv_val_control
)


## -----------------------------------------------------------------------------
#| label: fig-localresultsciso
#| fig-cap: "MBL results for Total Carbon predictions using the LOCAL algorithm. NNv: nearest-neighbor cross-validation."
#| fig-align: "center"
#| fig-width: 7.5
#| fig-height: 4
plot(local_ciso, main = "")
local_ciso


## -----------------------------------------------------------------------------
#| label: bestk
#| eval: true
#| echo: false
#| results: hide
bestk <- which.min(local_ciso$validation_results$nearest_neighbor_validation$rmse)
bestk <- local_ciso$validation_results$nearest_neighbor_validation$k[bestk]


## -----------------------------------------------------------------------------
#| label: get-predictions
#| eval: true
#| results: hide
#| fig-show: hide
bki <- which.min(local_ciso$validation_results$nearest_neighbor_validation$rmse)
bk <- local_ciso$validation_results$nearest_neighbor_validation$k[bki]

# All the prediction results are stored in:
local_ciso$results

# The get_predictions function makes easier to retrieve the
# predictions from the previous object
ciso_hat <- as.matrix(get_predictions(local_ciso))[, bki]


## -----------------------------------------------------------------------------
#| label: fig-plot-predictions
#| fig-cap: "Predicted vs reference values for Total Carbon using the LOCAL algorithm with the best neighborhood size according to nearest-neighbor validation."
#| eval: true
#| fig-width: 5
#| fig-height: 5
# Plot predicted vs reference
rng <- range(ciso_hat, test_y, na.rm = TRUE)
plot(ciso_hat, test_y,
     xlim = rng,
     ylim = rng,
     xlab = "Predicted Total Carbon, %",
     ylab = "Total Carbon, %",
     main = "LOCAL using a fixed k", 
     cex = 1.5,
     pch = 16, col = rgb(0.5, 0.5, 0.5, 0.6))
grid(lty = 1)
abline(0, 1, col = "red")


## -----------------------------------------------------------------------------
#| eval: true
# Prediction RMSE:
sqrt(mean((ciso_hat - test_y)^2, na.rm = TRUE))

# Squared R
cor(ciso_hat, test_y, use = "complete.obs")^2


## -----------------------------------------------------------------------------
#| label: mbl-diss-threshold
#| eval: true
#| results: hide
#| fig-show: hide
# Create a vector of dissimilarity thresholds to evaluate
# since the correlation dissimilarity will be used
# these thresholds need to be > 0 and <= 1
dths <- seq(0.025, 0.15, by = 0.025)

# Indicate the minimum and maximum sizes allowed for the neighborhood
k_min <- 30
k_max <- 150

local_ciso_diss <- mbl(
  Xr = train_x[!is.na(train_y), ],
  Yr = train_y[!is.na(train_y)],
  Xu = test_x,
  neighbors = neighbors_diss(threshold = dths, k_min = k_min, k_max = k_max),
  diss_method = diss_correlation(center = FALSE),
  diss_usage = ignore_diss,
  fit_method = my_wapls,
  control = nnv_val_control
)


## -----------------------------------------------------------------------------
#| label: plot-diss-results
#| eval: false
# plot(local_ciso_diss)


## -----------------------------------------------------------------------------
#| label: bestd
#| eval: true
#| echo: false
#| results: hide
bestd <- which.min(local_ciso_diss$validation_results$nearest_neighbor_validation$rmse)
bestd <- local_ciso_diss$validation_results$nearest_neighbor_validation$k_diss[bestd]


## -----------------------------------------------------------------------------
#| label: show-diss-results
#| eval: true
local_ciso_diss


## -----------------------------------------------------------------------------
#| label: diss-predictions
#| eval: true
#| results: hide
#| fig-show: hide
# Best distance threshold
bdi <- which.min(local_ciso_diss$validation_results$nearest_neighbor_validation$rmse)
bd <- local_ciso_diss$validation_results$nearest_neighbor_validation$k_diss[bdi]

# Predictions for the best distance
ciso_diss_hat <- as.matrix(get_predictions(local_ciso_diss))[, bdi]


## -----------------------------------------------------------------------------
#| label: fig-diss-predictions
#| fig-cap: "Predicted vs reference values for Total Carbon using the LOCAL algorithm with the best neighborhood size according to nearest-neighbor validation and distance thresholds."
#| eval: true
# Plot predicted vs reference
plot(ciso_diss_hat, test_y,
     xlim = rng,
     ylim = rng,
     xlab = "Predicted Total Carbon, %",
     ylab = "Total Carbon, %",
     main = "LOCAL using a distance threshold \nfor neighbor retrieval", 
     cex = 1.5,
     pch = 16, col = rgb(0.5, 0.5, 0.5, 0.6))
grid()
abline(0, 1, col = "red")

## -----------------------------------------------------------------------------
# RMSE
sqrt(mean((ciso_diss_hat - test_y)^2, na.rm = TRUE))

# Squared R
cor(ciso_diss_hat, test_y, use = "complete.obs")^2


## -----------------------------------------------------------------------------
train_x <- NIRsoil$spc_pr[NIRsoil$train == 1, ]
train_cec <- NIRsoil$CEC[NIRsoil$train == 1]

test_x  <- NIRsoil$spc_pr[NIRsoil$train == 0, ]
test_cec  <- NIRsoil$CEC[NIRsoil$train == 0]


## -----------------------------------------------------------------------------
#| label: cec-examples
#| results: hide
#| eval: true
#| fig-show: hide
# Define the WAPLS fitting method
my_wapls <- fit_wapls(min_ncomp = 2, max_ncomp = 25, scale = FALSE)

# mbl_cor: LOCAL algorithm with correlation dissimilarity
dth_cor <- seq(0.01, 0.3, by = 0.08)
mbl_cor <- mbl(
  Xr = train_x[!is.na(train_cec), ],
  Yr = train_cec[!is.na(train_cec)],
  Xu = test_x,
  neighbors = neighbors_diss(threshold = dth_cor, k_min = 80, k_max = 150),
  diss_method = diss_correlation(),
  diss_usage = "none",
  fit_method = my_wapls,
  control = nnv_val_control
)

# mbl_pc: PCA dissimilarity with dissimilarity matrix as predictors
dth_pc <- seq(0.05, 1, by = 0.4)
mbl_pc <- mbl(
  Xr = train_x[!is.na(train_cec), ],
  Yr = train_cec[!is.na(train_cec)],
  Xu = test_x,
  neighbors = neighbors_diss(threshold = dth_pc, k_min = 80, k_max = 150),
  diss_method = diss_pca(ncomp = ncomp_by_opc(), scale = TRUE),
  diss_usage = "predictors",
  fit_method = my_wapls,
  control = nnv_val_control
)

# mbl_pls: PLS dissimilarity without dissimilarity predictors
mbl_pls <- mbl(
  Xr = train_x[!is.na(train_cec), ],
  Yr = train_cec[!is.na(train_cec)],
  Xu = test_x,
  Yu = test_cec,
  neighbors = neighbors_diss(threshold = dth_pc, k_min = 80, k_max = 150),
  diss_method = diss_pls(ncomp = ncomp_by_opc(), scale = TRUE),
  diss_usage = "none",
  fit_method = my_wapls,
  control = nnv_val_control
)

# mbl_gpr: Gaussian process regression with PCA dissimilarity as predictors
mbl_gpr <- mbl(
  Xr = train_x[!is.na(train_cec), ],
  Yr = train_cec[!is.na(train_cec)],
  Xu = test_x,
  neighbors = neighbors_diss(threshold = dth_pc, k_min = 80, k_max = 150),
  diss_method = diss_pca(ncomp = ncomp_by_opc(), scale = TRUE),
  diss_usage = "predictors",
  fit_method = fit_gpr(),
  control = nnv_val_control
)


## -----------------------------------------------------------------------------
#| label: collect-predictions
#| eval: true
#| fig-show: hide
# Get the indices of the best results according to
# nearest neighbor validation statistics
c_val_name <- "validation_results"
c_nn_val_name <- "nearest_neighbor_validation"

bi_cor <- which.min(mbl_cor[[c_val_name]][[c_nn_val_name]]$rmse)
bi_pc  <- which.min(mbl_pc[[c_val_name]][[c_nn_val_name]]$rmse)
bi_pls <- which.min(mbl_pls[[c_val_name]][[c_nn_val_name]]$rmse)
bi_gpr <- which.min(mbl_gpr[[c_val_name]][[c_nn_val_name]]$rmse)

preds <- cbind(
  get_predictions(mbl_cor)[, bi_cor],
  get_predictions(mbl_pc)[, bi_pc],
  get_predictions(mbl_pls)[, bi_pls],
  get_predictions(mbl_gpr)[, bi_gpr]
)

colnames(preds) <- c("mbl_cor", "mbl_pc", "mbl_pls", "mbl_gpr")
preds <- as.matrix(preds)

# R2s
cor(test_cec, preds, use = "complete.obs")^2

# RMSEs
colMeans((preds - test_cec)^2, na.rm = TRUE)^0.5


## -----------------------------------------------------------------------------
#| label: cec-plots-code
#| eval: false
#| fig-show: hide
# old_par <- par("mfrow", "mar")
# par(mfrow = c(2, 2))
# 
# plot(test_cec, preds[, "mbl_cor"],
#      xlab = "CEC, meq/100g",
#      ylab = "Predicted CEC, meq/100g",
#      main = "mbl_cor (LOCAL)")
# abline(0, 1, col = "red")
# 
# plot(test_cec, preds[, "mbl_pc"],
#      xlab = "CEC, meq/100g",
#      ylab = "Predicted CEC, meq/100g",
#      main = "mbl_pc")
# abline(0, 1, col = "red")
# 
# plot(test_cec, preds[, "mbl_pls"],
#      xlab = "CEC, meq/100g",
#      ylab = "Predicted CEC, meq/100g",
#      main = "mbl_pls")
# abline(0, 1, col = "red")
# 
# plot(test_cec, preds[, "mbl_gpr"],
#      xlab = "CEC, meq/100g",
#      ylab = "Predicted CEC, meq/100g",
#      main = "mbl_gpr")
# abline(0, 1, col = "red")
# 
# par(old_par)


## -----------------------------------------------------------------------------
#| label: fig-mblcomparisons
#| fig-cap: "CEC prediction results for the different MBL configurations tested"
#| fig-align: center
#| fig-width: 8
#| fig-height: 8
#| echo: false
old_par <- par("mfrow", "mar")
par(mfrow = c(2, 2), pch = 16, mar = c(4, 4, 4, 4))

my_cols <- c(
  "mbl_cor" = "#D42B08CC",
  "mbl_pc"  = "#750E3380",
  "mbl_pls" = "#EFBF4780",
  "mbl_gpr" = "#5177A180"
)

# R2s
r2s <- drop(round(cor(test_cec, preds, use = "complete.obs")^2, 2))

# RMSEs
rmses <- round(colMeans((preds - test_cec)^2, na.rm = TRUE)^0.5, 2)

plot_titles <- c(
  "mbl_cor" = "mbl_cor (LOCAL)",
  "mbl_pc"  = "mbl_pc",
  "mbl_pls" = "mbl_pls",
  "mbl_gpr" = "mbl_gpr"
)

p <- sapply(
  colnames(preds),
  FUN = function(label, y, yhats, cols, rsq, rmse, titles) {
    plot(
      x = y,
      y = yhats[, label],
      xlim = range(y, na.rm = TRUE),
      ylim = range(y, na.rm = TRUE),
      xlab = "CEC, meq/100g",
      ylab = "Predicted CEC, meq/100g",
      col = cols[label]
    )
    title(titles[label])
    title(paste0("\n\n\n RMSE: ", rmse[label], "; R²: ", rsq[label]), cex.main = 1)
    grid(col = "#80808080", lty = 1)
    abline(0, 1, col = "#FF1A0080")
  },
  y = test_cec,
  yhats = preds,
  cols = my_cols,
  rsq = r2s,
  rmse = rmses,
  titles = plot_titles
)

par(old_par)


## -----------------------------------------------------------------------------
#| label: yu-argument
#| eval: true
# Use Yu argument to validate the predictions
pc_pred_cec_yu <- mbl(
  Xr = train_x[!is.na(train_cec), ],
  Yr = train_cec[!is.na(train_cec)],
  Xu = test_x,
  Yu = test_cec,
  neighbors = neighbors_k(80),
  diss_method = diss_pca(scale = TRUE),
  diss_usage = "none",
  verbose = FALSE,
  control = mbl_control()
)

pc_pred_cec_yu


## -----------------------------------------------------------------------------
#| label: parallel-example
#| eval: false
# # Running the mbl function using multiple cores
# 
# # Execute with two cores, if available, ...
# n_cores <- 2
# 
# # ... if not then go with 1 core
# if (parallel::detectCores() < 2) {
#   n_cores <- 1
# }
# 
# # Set the number of cores
# library(doParallel)
# clust <- makeCluster(n_cores)
# registerDoParallel(clust)
# 
# # Alternatively:
# # library(doSNOW)
# # clust <- makeCluster(n_cores, type = "SOCK")
# # registerDoSNOW(clust)
# # getDoParWorkers()
# 
# pc_pred_cec <- mbl(
#   Xr = train_x[!is.na(train_cec), ],
#   Yr = train_cec[!is.na(train_cec)],
#   Xu = test_x,
#   neighbors = neighbors_k(seq(40, 100, by = 10)),
#   diss_method = diss_pca(ncomp = ncomp_by_opc(), scale = TRUE),
#   diss_usage = "none",
#   control = mbl_control()
# )
# 
# # Go back to sequential processing
# registerDoSEQ()
# try(stopCluster(clust))
# 
# pc_pred_cec

