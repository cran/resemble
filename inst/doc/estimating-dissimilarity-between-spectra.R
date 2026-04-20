## -----------------------------------------------------------------------------
#| echo: false
Sys.setenv(OMP_NUM_THREADS = 2)


## -----------------------------------------------------------------------------
#| message: false
#| echo: false
library(resemble)


## -----------------------------------------------------------------------------
diss_correlation()


## -----------------------------------------------------------------------------
diss_correlation(center = TRUE, scale = TRUE)


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
#| label: pca-diss
#| eval: TRUE

# Default: variance-based component selection (ncomp_by_var(0.01))
d_pca <- dissimilarity(train_x, diss_method = diss_pca())

# With OPC-based component selection (requires Yr)
d_pca_opc <- dissimilarity(
  train_x,
  Yr = train_y,
  diss_method = diss_pca(ncomp = ncomp_by_opc(30))
)

# Between training and test sets
d_pca_tr_ts <- dissimilarity(
 train_x,
 Xu = test_x,
 Yr = train_y,
 diss_method = diss_pca(
   ncomp = ncomp_by_opc(30),
   return_projection = TRUE
 )
)

d_pca_tr_ts


## -----------------------------------------------------------------------------
# Example of the first 5 rows and columns of the dissimilarity 
# matrix between training and test sets
# Rows: training observations; 
# Columns: test observations
d_pca_tr_ts$dissimilarity[1:5, 1:5]


## -----------------------------------------------------------------------------
d_pca_tr_ts$projection


## ----pls-diss, eval = FALSE---------------------------------------------------
# # Default: OPC-based component selection
# d_pls <- dissimilarity(
#   train_x,
#   Xu = test_x,
#   Yr = train_y,
#   diss_method = diss_pls()
# )
# 
# # Fixed number of components
# d_pls_fixed <- dissimilarity(
#   train_x,
#   Xu = test_x,
#   Yr = train_y,
#   diss_method = diss_pls(ncomp = 15)
# )


## ----cor-diss, eval = FALSE---------------------------------------------------
# # Standard correlation dissimilarity
# d_cor <- dissimilarity(train_x, diss_method = diss_correlation())
# d_cor


## -----------------------------------------------------------------------------
correlation_baser <- 0.5 * (1 - cor(t(train_x)))
corrrelation_resemble <- dissimilarity(
  train_x, diss_method = diss_correlation(center = FALSE)
)
# the maximum discreepancy between the two implementations should 
# be very close to zero
max(abs(corrrelation_resemble$dissimilarity - correlation_baser))


## -----------------------------------------------------------------------------
#| eval: false
# # Compare computational speed: resemble vs base R correlation dissimilarity
# n_iter <- 50
# 
# # Transpose once before timing (fair comparison)
# train_x_t <- t(train_x)
# 
# # resemble implementation
# time_resemble <- system.time({
#   for (i in seq_len(n_iter)) {
#     cor_resemble <- dissimilarity(train_x, diss_method = diss_correlation(center = FALSE))
#   }
# })
# 
# # base R implementation
# time_base <- system.time({
#   for (i in seq_len(n_iter)) {
#     cor_base <- 0.5 * (1 - cor(train_x_t))
#   }
# })
# 
# # Results
# data.frame(
#   method = c("resemble", "base R"),
#   elapsed_sec = c(time_resemble["elapsed"], time_base["elapsed"]),
#   per_call_ms = c(time_resemble["elapsed"], time_base["elapsed"]) / n_iter * 1000
# )


## ----cor-diss-mw, eval = FALSE------------------------------------------------
# # Moving window correlation (window size must be odd)
# d_cor_mw <- dissimilarity(
#   train_x,
#   Xu = test_x,
#   diss_method = diss_correlation(ws = 51)
# )
# d_cor_mw


## ----omp-threads-use, eval = FALSE--------------------------------------------
# # Use 4 threads
# Sys.setenv(OMP_NUM_THREADS = 4)
# 
# # Or use all available cores
# Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())


## ----omp-threads-restrict, eval = FALSE---------------------------------------
# # Limit OpenMP threads
# Sys.setenv(OMP_NUM_THREADS = 1)
# 
# # If using OpenBLAS
# Sys.setenv(OPENBLAS_NUM_THREADS = 1)
# 
# # If using MKL
# Sys.setenv(MKL_NUM_THREADS = 1)


## ----blas-threads, eval = FALSE-----------------------------------------------
# RhpcBLASctl::blas_set_num_threads(1)
# RhpcBLASctl::omp_set_num_threads(1)


## ----euclid-diss, eval = FALSE------------------------------------------------
# d_euclid <- dissimilarity(train_x, Xu = test_x, diss_method = diss_euclidean())
# d_euclid


## -----------------------------------------------------------------------------
ed_resemble <- dissimilarity(train_x, diss_method = diss_euclidean(center = FALSE))
ed_baser <- as.matrix(dist(train_x, method = "euclidean"))

# scaling the base R Euclidean distance by the number of variables:
ed_baser_scaled <- sqrt((ed_baser^2) / ncol(train_x))

# differences between the two implementations
max(abs(ed_baser_scaled - ed_resemble$dissimilarity))


## -----------------------------------------------------------------------------
#| label: euclid-speed
#| eval: false
# # Compare computational speed: resemble vs base R
# n_iter <- 50
# 
# # resemble implementation
# time_resemble <- system.time({
#   for (i in seq_len(n_iter)) {
#     ed_resemble <- dissimilarity(train_x, diss_method = diss_euclidean(center = FALSE))
#   }
# })
# 
# # base R implementation
# time_base <- system.time({
#   for (i in seq_len(n_iter)) {
#     ed_base <- as.matrix(dist(train_x, method = "euclidean"))
#   }
# })
# 
# # Results
# data.frame(
#   method = c("resemble", "base R"),
#   elapsed_sec = c(time_resemble["elapsed"], time_base["elapsed"]),
#   per_call_ms = c(time_resemble["elapsed"], time_base["elapsed"]) / n_iter * 1000
# )


## ----cosine-diss, eval = FALSE------------------------------------------------
# d_cosine <- dissimilarity(train_x, Xu = test_x, diss_method = diss_cosine())
# d_cosine


## -----------------------------------------------------------------------------
# Compute dissimilarity matrices using different methods
d_pca_var <- dissimilarity(train_x, diss_method = diss_pca(scale = TRUE))

d_pca_opc <- dissimilarity(
  train_x,
  Yr = train_y,
  diss_method = diss_pca(ncomp = ncomp_by_opc(30), scale = TRUE)
)

d_pls_opc <- dissimilarity(
  train_x,
  Yr = train_y,
  diss_method = diss_pls(ncomp = ncomp_by_opc(), scale = TRUE)
)

d_cor <- dissimilarity(
  train_x,
  diss_method = diss_correlation(ws = 51, scale = TRUE)
)

d_euclid <- dissimilarity(
  train_x,
  diss_method = diss_euclidean(scale = TRUE)
)

d_cosine <- dissimilarity(
  train_x,
  diss_method = diss_cosine(scale = TRUE)
)


## -----------------------------------------------------------------------------
# Evaluate each method using Ciso as side information
side_info <- as.matrix(train_y)

eval_results <- list(
  "PCA (var)" = diss_evaluate(
    d_pca_var$dissimilarity, side_info
  ),
  "PCA (opc)" = diss_evaluate(
    d_pca_opc$dissimilarity, side_info
  ),
  "PLS (opc)" = diss_evaluate(
    d_pls_opc$dissimilarity, side_info
  ),
  correlation = diss_evaluate(
    d_cor$dissimilarity, side_info
  ),
  Euclidean = diss_evaluate(
    d_euclid$dissimilarity, side_info
  ),
  cosine = diss_evaluate(
    d_cosine$dissimilarity, side_info
  )
)

# Extract RMSD and correlation for each method
comparison <- do.call(rbind, lapply(names(eval_results), function(nm) {
  data.frame(
    method = nm,
    rmsd   = eval_results[[nm]]$eval[, "rmsd"],
    r      = eval_results[[nm]]$eval[, "r"]
  )
}))
comparison


## -----------------------------------------------------------------------------
#| label: diss-eval-plot
#| fig-cap: "Comparison of dissimilarity methods based on first nearest-neighbor evaluation using `Ciso` as side information. In each panel, the x-axis shows the `Ciso` value of the nearest neighbor identified from the spectral dissimilarity matrix, and the y-axis shows the observed `Ciso` value of the corresponding sample. The dashed line represents perfect agreement (y = x)."
#| fig-align: "center"
#| fig-width: 8
#| fig-height: 6
#| echo: true
#| fig-show: hold
blue  <- "#3B82F6"
amber <- "#F59E0B"
slate_grid <- "#33415540"

old_par <- par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

for (nm in names(eval_results)) {
  xy <- eval_results[[nm]]$first_nn
  plot(
    xy[, 2], xy[, 1],
    xlab = "Ciso (1-NN)",
    ylab = "Ciso",
    col = adjustcolor(blue, alpha.f = 0.5),
    pch = 16,
    main = nm
  )
  grid(lty = 1, col = slate_grid)
  abline(0, 1, col = amber, lwd = 2)
}

par(old_par)

