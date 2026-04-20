## -----------------------------------------------------------------------------
#| message: false
#| echo: false
library(resemble)
library(prospectr)
Sys.setenv(OMP_NUM_THREADS = 2)


## -----------------------------------------------------------------------------
#| label: data-prep
#| message: false
library(resemble)
library(prospectr)

data(NIRsoil)

# Wavelengths
wavs <- as.numeric(colnames(NIRsoil$spc))

# Preprocess: detrend + first derivative
NIRsoil$spc_pr <- savitzkyGolay(
  detrend(NIRsoil$spc, wav = wavs),
  m = 1, p = 1, w = 7
)

# Split into training and test sets
# Note: missing values in the response are allowed in liblex
train_x <- NIRsoil$spc_pr[NIRsoil$train == 1, ]
train_y <- NIRsoil$Ciso[NIRsoil$train == 1]
test_x  <- NIRsoil$spc_pr[NIRsoil$train == 0, ]
test_y  <- NIRsoil$Ciso[NIRsoil$train == 0]

cat("Training set:", nrow(train_x), "observations\n")
cat("Test set:", nrow(test_x), "observations\n")


## -----------------------------------------------------------------------------
#| label: neighbors-example1
# Fixed neighborhood sizes to evaluate
neighbors_k(k = seq(40, 120, by = 20))


## -----------------------------------------------------------------------------
#| label: neighbors-example2
# Threshold-based selection with bounds
neighbors_diss(
  threshold = seq(0.05, 0.5, length.out = 10), 
  k_min = 40, 
  k_max = 150
)


## -----------------------------------------------------------------------------
#| label: diss-example
# Moving-window correlation dissimilarity (recommended for spectra)
diss_correlation(ws = 37, scale = TRUE)

# PCA-based Mahalanobis distance with optimized components
diss_pca()


## -----------------------------------------------------------------------------
#| label: fit-method-example
# waPLS with modified PLS algorithm and scaling
fit_wapls(
  min_ncomp = 3, 
  max_ncomp = 15, 
  scale = TRUE, 
  method = "mpls"
)


## -----------------------------------------------------------------------------
#| label: control-example
#| eval: false
# # Build library with hyperparameter tuning
# liblex_control(mode = "build", tune = TRUE)


## -----------------------------------------------------------------------------
#| eval: false
# ciso_lib_k <- liblex(
#   Xr = train_x,
#   Yr = train_y,
#   neighbors = neighbors_k(seq(40, 80, by = 20)),
#   diss_method = diss_correlation(ws = 37, scale = TRUE),
#   fit_method = fit_wapls(
#     min_ncomp = 3,
#     max_ncomp = 15,
#     scale = TRUE,
#     method = "mpls"
#   ),
#   control = liblex_control(tune = TRUE)
# )
# 
# ciso_lib_k


## -----------------------------------------------------------------------------
#| label: liblex-fixed-k
#| cache: true
#| echo: false
#| message: false
#| warning: false
ciso_lib_k <- liblex(
  Xr = train_x, 
  Yr = train_y, 
  neighbors = neighbors_k(seq(40, 80, by = 20)),
  diss_method = diss_correlation(ws = 37, scale = TRUE),
  fit_method = fit_wapls(
    min_ncomp = 3,
    max_ncomp = 15,
    scale = TRUE,
    method = "mpls"
  ),
  control = liblex_control(tune = TRUE),
  verbose = FALSE
)

ciso_lib_k


## -----------------------------------------------------------------------------
#| label: fig-liblex
#| fig-cap: "Top: Best model obtained for each neighborhood size (based on the RMSE). Bottom: Centroids of the neighborhoods, usled later in prediction to selet the models to be used."
#| fig-align: "center"
#| fig-width: 7.5
#| fig-height: 6.5
plot(ciso_lib_k)


## -----------------------------------------------------------------------------
#| label: liblex-threshold
#| cache: true
#| echo: false
ciso_lib_thr <- liblex(
  Xr = train_x, 
  Yr = train_y, 
  neighbors = neighbors_diss(
    threshold = seq(0.05, 0.5, length.out = 5), 
    k_min = 40, 
    k_max = 150
  ),
  diss_method = diss_correlation(ws = 37, scale = TRUE),
  fit_method = fit_wapls(
    min_ncomp = 3,
    max_ncomp = 15,
    scale = TRUE,
    method = "mpls"
  ),
  control = liblex_control(tune = TRUE),
  verbose = FALSE
)

ciso_lib_thr


## -----------------------------------------------------------------------------
#| cache: true
#| eval: false
# ciso_lib_thr <- liblex(
#   Xr = train_x,
#   Yr = train_y,
#   neighbors = neighbors_diss(
#     threshold = seq(0.05, 0.5, length.out = 5),
#     k_min = 40,
#     k_max = 150
#   ),
#   diss_method = diss_correlation(ws = 37, scale = TRUE),
#   fit_method = fit_wapls(
#     min_ncomp = 3,
#     max_ncomp = 15,
#     scale = TRUE,
#     method = "mpls"
#   ),
#   control = liblex_control(tune = TRUE)
# )
# 
# ciso_lib_thr


## -----------------------------------------------------------------------------
#| label: fig-centroids
#| fig-cap: "Neighborhood centroids (blue) compared to test spectra (red)"
#| fig-width: 8
#| fig-height: 5
#| fig-align: center
wavs_pr <- as.numeric(colnames(ciso_lib_k$scaling$local_x_center))

matplot(
  wavs_pr,
  t(test_x),
  col = rgb(1, 0, 0, 0.3),
  lty = 1,
  type = "l",
  xlab = "Wavelength (nm)",
  ylab = "First derivative detrended absorbance"
)

matlines(
  wavs_pr,
  t(ciso_lib_k$scaling$local_x_center),
  col = rgb(0, 0, 1, 0.3),
  lty = 1
)

grid(lty = 1)

legend(
  "topright",
  legend = c("Samples to predict", "Neighborhood centroids"),
  col = c(rgb(1, 0, 0, 0.8), rgb(0, 0, 1, 0.8)),
  lty = 1,
  lwd = 2,
  bty = "n"
)


## -----------------------------------------------------------------------------
#| label: fig-regression-coefficients
#| fig-cap: "Regression coefficients (top) and intercept distribution (bottom) of the local experts"
#| fig-width: 8
#| fig-height: 7
#| fig-align: center

par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))

# Regression coefficients across wavelengths
models <- ciso_lib_k$coefficients

matplot(
  x = wavs_pr,
  y = t(models$B),
  type = "l",
  lty = 1,
  col = rgb(0.23, 0.51, 0.96, 0.15),
  xlab = "Wavelength (nm)",
  ylab = "Regression coefficient",
  main = paste0("Regression coefficients (", nrow(models$B), " experts)")
)

# Add mean coefficient profile
lines(wavs_pr, colMeans(models$B), col = "red", lwd = 2)
grid(lty = 1)
legend(
  "topright",
  legend = c("Individual experts", "Mean"),
  col = c(rgb(0.23, 0.51, 0.96, 0.6), "red"),
  lty = 1,
  lwd = c(1, 2),
  bty = "n"
)

# Distribution of intercepts
plot(
  density(models$B0, na.rm = TRUE),
  col = rgb(0.23, 0.51, 0.96),
  lwd = 2,
  xlab = "Intercept value",
  ylab = "Density",
  main = "Distribution of intercepts"
)
polygon(
  density(models$B0, na.rm = TRUE),
  col = rgb(0.23, 0.51, 0.96, 0.2),
  border = NA
)
abline(v = mean(models$B0, na.rm = TRUE), col = "red", lty = 2, lwd = 2)
grid(lty = 1)
legend(
  "topright",
  legend = c("Density", "Mean"),
  col = c(rgb(0.23, 0.51, 0.96), "red"),
  lty = c(1, 2),
  lwd = 2,
  bty = "n"
)

par(mfrow = c(1, 1))


## -----------------------------------------------------------------------------
#| label: kmeans-anchors
# Select 350 representative anchors using k-means sampling
# on the first 20 principal components
set.seed(1124)
kms <- naes(
  train_x, 
  k = 350, 
  pc = 20, 
  iter.max = 100,
  .center = TRUE, 
  .scale = TRUE
)

anchor_km <- kms$model

cat("Selected", length(anchor_km), "anchors via k-means\n")


## -----------------------------------------------------------------------------
#| label: liblex-anchored
#| cache: true
#| echo: false
ciso_lib_anchored <- liblex(
  Xr = train_x, 
  Yr = train_y, 
  neighbors = neighbors_k(seq(40, 80, by = 20)),
  diss_method = diss_correlation(ws = 37, scale = TRUE),
  fit_method = fit_wapls(
    min_ncomp = 3,
    max_ncomp = 15,
    scale = TRUE,
    method = "mpls"
  ),
  anchor_indices = anchor_km,
  control = liblex_control(tune = TRUE),
  verbose = FALSE
)

ciso_lib_anchored


## -----------------------------------------------------------------------------
#| cache: true
#| eval: false
# ciso_lib_anchored <- liblex(
#   Xr = train_x,
#   Yr = train_y,
#   neighbors = neighbors_k(seq(40, 80, by = 20)),
#   diss_method = diss_correlation(ws = 37, scale = TRUE),
#   fit_method = fit_wapls(
#     min_ncomp = 3,
#     max_ncomp = 15,
#     scale = TRUE,
#     method = "mpls"
#   ),
#   anchor_indices = anchor_km,
#   control = liblex_control(tune = TRUE)
# )
# 
# ciso_lib_anchored


## -----------------------------------------------------------------------------
# Number of experts in the anchored library
n_experts <- nrow(ciso_lib_anchored$coefficients$B)
cat("Number of experts in the anchored library:", n_experts, "\n")


## -----------------------------------------------------------------------------
#| label: fig-centroids-anchored
#| fig-cap: "Neighborhood centroids from anchored library compared to test spectra"
#| fig-width: 8
#| fig-height: 5
#| fig-align: center

wavs_pr <- as.numeric(colnames(ciso_lib_anchored$scaling$local_x_center))
n_experts <- nrow(ciso_lib_anchored$scaling$local_x_center)
n_test <- nrow(test_x)

# Plot test spectra first (background)
matplot(
  wavs_pr,
  t(test_x),
  col = rgb(0.8, 0.2, 0.2, 0.2),
  lty = 1,
  type = "l",
  xlab = "Wavelength (nm)",
  ylab = "First derivative detrended absorbance",
  main = paste0("Coverage: ", n_experts, " experts vs ", n_test, " test samples")
)

# Overlay centroids
matlines(
  wavs_pr,
  t(ciso_lib_anchored$scaling$local_x_center),
  col = rgb(0.23, 0.51, 0.96, 0.3),
  lty = 1
)
grid(lty = 1)

legend(
  "topright",
  legend = c(
    paste0("Test samples (n = ", n_test, ")"),
    paste0("Centroids (n = ", n_experts, ")")
  ),
  col = c(
    rgb(0.8, 0.2, 0.2, 0.5),
    rgb(0.23, 0.51, 0.96, 0.5)
  ),
  lty = 1,
  lwd = c(1, 1, 2, 2),
  bty = "n"
)


## -----------------------------------------------------------------------------
#| label: predict-basic
ciso_pred <- predict(ciso_lib_k, newdata = test_x, verbose = FALSE)

# Prediction output structure
names(ciso_pred)

# Main predictions with uncertainty
head(ciso_pred$predictions)


## -----------------------------------------------------------------------------
#| label: weighting-options
#| eval: false
# # Gaussian kernel (default)
# predict(ciso_lib_k, newdata = test_x, weighting = "gaussian")
# 
# # Tricubic kernel (similar to LOCAL algorithm)
# predict(ciso_lib_k, newdata = test_x, weighting = "tricube")
# 
# # Equal weights
# predict(ciso_lib_k, newdata = test_x, weighting = "none")


## -----------------------------------------------------------------------------
#| label: enforce-indices
#| eval: false
# # Always include specific experts in predictions
# predict(ciso_lib_k, newdata = test_x, enforce_indices = c(1, 5, 10))


## -----------------------------------------------------------------------------
#| label: evaluation
pred_values <- ciso_pred$predictions$pred
pred_sd <- ciso_pred$predictions$pred_sd

# Performance metrics (excluding missing values)
complete <- !is.na(test_y)
r2 <- cor(pred_values[complete], test_y[complete])^2
rmse <- sqrt(mean((pred_values[complete] - test_y[complete])^2))

cat("R²:  ", round(r2, 3), "\n")
cat("RMSE:", round(rmse, 3), "\n")


## -----------------------------------------------------------------------------
#| label: uncertainty-filter
# Filter predictions by uncertainty threshold
unc_threshold <- quantile(pred_sd[complete], 0.75)
reliable <- complete & (pred_sd < unc_threshold)

r2_filtered <- cor(pred_values[reliable], test_y[reliable])^2
rmse_filtered <- sqrt(mean((pred_values[reliable] - test_y[reliable])^2))

cat("After filtering high-uncertainty predictions:\n")
cat("  Retained:", sum(reliable), "/", sum(complete), "observations\n")
cat("  R²:      ", round(r2_filtered, 3), "\n")
cat("  RMSE:    ", round(rmse_filtered, 3), "\n")


## -----------------------------------------------------------------------------
#| label: fig-validation
#| fig-cap: "Predicted versus observed values"
#| fig-width: 6
#| fig-height: 5
#| fig-align: center
lims <- range(pred_values[complete], test_y[complete], na.rm = TRUE)

plot(
  pred_values[complete],
  test_y[complete],
  pch = 16,
  col = rgb(0, 0, 0, 0.5),
  xlab = "Predicted",
  ylab = "Observed",
  xlim = lims,
  ylim = lims
)

abline(0, 1, col = "red", lwd = 2)
grid(lty = 1)


## -----------------------------------------------------------------------------
#| label: parallel-example
#| eval: false
# library(doParallel)
# 
# # Register parallel backend
# n_cores <- min(parallel::detectCores() - 1, 4)
# cl <- makeCluster(n_cores)
# registerDoParallel(cl)
# 
# # Build library with parallel processing
# ciso_lib_parallel <- liblex(
#   Xr = train_x,
#   Yr = train_y,
#   neighbors = neighbors_k(seq(40, 120, by = 20)),
#   diss_method = diss_correlation(ws = 37, scale = TRUE),
#   fit_method = fit_wapls(
#     min_ncomp = 3,
#     max_ncomp = 15,
#     scale = TRUE,
#     method = "mpls"
#   ),
#   control = liblex_control(
#     tune = TRUE,
#     allow_parallel = TRUE,
#     chunk_size = 10
#   )
# )
# 
# # Clean up
# stopCluster(cl)
# registerDoSEQ()

