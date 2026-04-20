## -----------------------------------------------------------------------------
#| echo: false
Sys.setenv(OMP_NUM_THREADS = 2)


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
#| eval: true
# Retain components that individually explain at least 1% of variance
proj_var <- ortho_projection(train_x, ncomp = ncomp_by_var(0.01))
proj_var


## -----------------------------------------------------------------------------
#| eval: true
# Retain enough components to explain at least 99% of variance
proj_cumvar <- ortho_projection(train_x, ncomp = ncomp_by_cumvar(0.99))
proj_cumvar


## -----------------------------------------------------------------------------
#| eval: true
# Optimal component selection using side information
# using default max_ncomp which is 40
proj_opc <- ortho_projection(train_x, Yr = train_y, ncomp = ncomp_by_opc())
proj_opc


## -----------------------------------------------------------------------------
# this needs to be passed to the ncomp argument in ortho_projection()
ncomp_by_opc(max_ncomp = 15)


## -----------------------------------------------------------------------------
#| eval: false
# # Using ncomp_fixed()
# proj_fixed <- ortho_projection(train_x, ncomp = ncomp_fixed(10))
# 
# # Equivalent: passing an integer directly
# proj_fixed <- ortho_projection(train_x, ncomp = 10)


## -----------------------------------------------------------------------------
#| results: hide
# PCA with default component selection (ncomp_by_var(0.01))
pca_train <- ortho_projection(train_x, method = "pca")
pca_train


## -----------------------------------------------------------------------------
#| label: fig-pca-variance
#| fig-cap: "Individual contribution to the explained variance for each component (left) and cumulative variance explained by the principal components (right)."
#| fig-width: 7
#| fig-height: 3
plot(pca_train)


## -----------------------------------------------------------------------------
#| results: hide
# PCA using NIPALS
pca_nipals_train <- ortho_projection(train_x, method = "pca_nipals")


## -----------------------------------------------------------------------------
#| results: hide
# PLS using SIMPLS with 10 components
pls_train <- ortho_projection(
  train_x,
  Yr = train_y,
  method = "simpls",
  ncomp = 10
)
pls_train


## -----------------------------------------------------------------------------
#| results: hide
# PCA with optimal component selection
pca_opc <- ortho_projection(
  train_x,
  Yr = train_y,
  method = "pca",
  ncomp = ncomp_by_opc(40)
)
pca_opc


## -----------------------------------------------------------------------------
#| label: fig-pca-opc
#| fig-cap: "RMSD between observations and their nearest neighbors as a function of the number of principal components."
#| fig-width: 5
#| fig-height: 4
plot(pca_opc, col = "#F59E0B")


## -----------------------------------------------------------------------------
#| results: hide
# Fit PLS projection model
pls_model <- ortho_projection(
  train_x,
  Yr = train_y,
  method = "simpls",
  ncomp = ncomp_by_opc(40),
  scale = TRUE
)

# Project test data
pls_test_scores <- predict(pls_model, newdata = test_x)


## -----------------------------------------------------------------------------
#| label: fig-pc-proj
#| fig-width: 4.5
#| fig-height: 5
#| fig-cap: "Projection of the training observations in the first two PLS components (blue) and projection of the test observations onto the same latent space (red)."
plot(
  pls_model$scores[, 1:2], 
  col = rgb(0.231, 0.51, 0.965, 0.3), 
  pch = 16, cex = 1.5,
  xlab = "PLS1", ylab = "PLS2"
)
grid(lty = 1)
points(
  pls_test_scores[, 1:2], 
  col = rgb(0.961, 0.620, 0.043, 0.7), 
  pch = 16, cex = 1.5
)
legend(
  "topright",
  legend = c("Training projection", "Predicted projection"),
  col = c(
    rgb(0.231, 0.51, 0.965, 0.3),
    rgb(0.961, 0.620, 0.043, 0.7)
  ),
  pch = 16,
  pt.cex = 1.5,
  box.col = "black"
)


## -----------------------------------------------------------------------------
#| results: hide
# Project training and test data simultaneously
pca_both <- ortho_projection(
  train_x,
  Xu = test_x,
  Yr = train_y,
  method = "pca",
  ncomp = ncomp_by_opc(40),
  scale = TRUE
)


## -----------------------------------------------------------------------------
#| results: hide
#| eval: false
# train_y2 <- NIRsoil$Nt[NIRsoil$train == 1]
# train_y3 <- NIRsoil$CEC[NIRsoil$train == 1]
# 
# # PLS with multivariate side information
# pls_multi <- ortho_projection(
#   train_x,
#   Xu = test_x,
#   Yr = cbind(train_y, train_y2, train_y3),
#   method = "pca",
#   ncomp = ncomp_by_opc(40),
#   scale = TRUE
# )
# 
# plot(pls_multi)

