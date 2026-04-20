## -----------------------------------------------------------------------------
#| echo: false
Sys.setenv(OMP_NUM_THREADS = 2)


## ----neighbors-k-example, eval = FALSE----------------------------------------
# neighbors_k(k = 50)
# neighbors_k(k = c(40, 60, 80, 100))


## ----neighbors-diss-example, eval = FALSE-------------------------------------
# neighbors_diss(threshold = 0.3)
# neighbors_diss(threshold = c(0.1, 0.2, 0.3), k_min = 10, k_max = 150)


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
set.seed(8011)
rnd_idc <- sample(nrow(test_x), 3)

k_fixed <- search_neighbors(
  Xr = train_x,
  Xu = test_x[rnd_idc, ],
  diss_method = diss_pca(ncomp = 2, return_projection = TRUE),
  neighbors = neighbors_k(30)
)


## -----------------------------------------------------------------------------
k_diss <- search_neighbors(
  Xr = train_x,
  Xu = test_x[rnd_idc, ],
  diss_method = diss_pca(ncomp = 2),
  neighbors = neighbors_diss(threshold = 0.25)
)


## -----------------------------------------------------------------------------
test_scores_indices <- grep("^Xu_", rownames(k_fixed$projection$scores))


## -----------------------------------------------------------------------------
#| label: fig-neighbors
#| fig-cap: "Nearest neighbors identified using the same dissimilarity metric but retained using two different selection methods."
#| fig-width: 8.5
#| fig-height: 4.5

old_par <- par(no.readonly = TRUE)
on.exit(par(old_par), add = TRUE)

par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

plot(
  k_fixed$projection$scores,
  pch = 16,
  col = rgb(0.5, 0.5, 0.5, 0.3),
  main = "Fixed k-nearest neighbors"
)
grid(lty = 1)
points(
  k_fixed$projection$scores[k_fixed$unique_neighbors, ],
  pch = 16,
  col = "dodgerblue"
)
points(
  k_fixed$projection$scores[test_scores_indices, ],
  pch = 16,
  cex = 1.5,
  col = "red"
)

plot(
  k_fixed$projection$scores,
  pch = 16,
  col = rgb(0.5, 0.5, 0.5, 0.3),
  main = "Neighbors selected by a \nthreshold dissimilarity"
)
grid(lty = 1)
points(
  k_fixed$projection$scores[k_diss$unique_neighbors, ],
  pch = 16,
  col = "dodgerblue"
)
points(
  k_fixed$projection$scores[test_scores_indices, ],
  pch = 16,
  cex = 1.5,
  col = "red"
)


## ----knn-pca, eval = FALSE----------------------------------------------------
# # matrix of neighbors
# k_fixed$neighbors
# 
# # matrix of neighbor distances (dissimilarity scores)
# k_fixed$neighbors_diss
# 
# # the index (in the training set) of the first two closest neighbors found in
# # training for the first observation in testing:
# k_fixed$neighbors[1:2, 1, drop = FALSE]
# 
# # the distances of the two closest neighbors found in
# # training for the first observation in testing:
# k_fixed$neighbors_diss[1:2, 1, drop = FALSE]
# 
# # the indices in training that fall in any of the
# # neighborhoods of testing
# k_fixed$unique_neighbors


## -----------------------------------------------------------------------------
#| eval: TRUE
#| label: knn-other-methods
# using PC dissimilarity with optimal selection of components
knn_opc <- search_neighbors(
  Xr = train_x,
  Xu = test_x,
  diss_method = diss_pca(
    ncomp = ncomp_by_opc(),
    scale = TRUE,
    return_projection = TRUE
  ),
  Yr = train_y,
  neighbors = neighbors_k(50)
)

# using PLS dissimilarity with optimal selection of components
knn_pls <- search_neighbors(
  Xr = train_x,
  Xu = test_x,
  diss_method = diss_pls(
    ncomp = ncomp_by_opc(),
    scale = TRUE
  ),
  Yr = train_y,
  neighbors = neighbors_k(50)
)

# using correlation dissimilarity
knn_cor <- search_neighbors(
  Xr = train_x,
  Xu = test_x,
  diss_method = diss_correlation(),
  neighbors = neighbors_k(50)
)

# using moving window correlation dissimilarity
knn_mw <- search_neighbors(
  Xr = train_x,
  Xu = test_x,
  diss_method = diss_correlation(ws = 51),
  neighbors = neighbors_k(50)
)


## -----------------------------------------------------------------------------
#| eval: TRUE
#| label: knn-threshold-example
#| results: hide
# a dissimilarity threshold
d_th <- 1

# the minimum number of observations required in each neighborhood
k_min <- 20

# the maximum number of observations allowed in each neighborhood
k_max <- 300

dnn_pca <- search_neighbors(
  Xr = train_x,
  Xu = test_x,
  diss_method = diss_pca(scale = TRUE),
  neighbors = neighbors_diss(threshold = d_th, k_min = k_min, k_max = k_max)
)

# matrix of neighbors. The minimum number of indices is 20 (given by k_min)
# and the maximum number of indices is 300 (given by k_max).
# NAs indicate "not a neighbor"
dnn_pca$neighbors

# this reports how many neighbors were found for each observation in 
# testing using the input distance threshold (column n_k) and how 
# many were finally selected (column final_n_k)
dnn_pca$k_diss_info

# matrix of neighbor distances
dnn_pca$neighbors_diss

# the indices in training that fall in any of the 
# neighborhoods of testing
dnn_pca$unique_neighbors


## -----------------------------------------------------------------------------
#| eval: TRUE
#| label: fig-knn-hist
#| fig-cap: "Histogram of the final neighborhood sizes after applying the dissimilarity threshold and the minimum and maximum neighborhood size constraints."
#| fig-width: 5
#| fig-height: 4
hist(
  dnn_pca$k_diss_info$final_n_k,
  breaks = 20,
  xlab = "Final neighborhood size",
  main = "",
  col = "dodgerblue"
)


## -----------------------------------------------------------------------------
#| eval: TRUE
#| label: knn-spike-example
# the indices of the observations that we want to "invite" to every neighborhood
forced_guests <- c(1, 5, 8, 9)

# using PC dissimilarity with optimal selection of components
knn_spiked <- search_neighbors(
  Xr = train_x,
  Xu = test_x,
  diss_method = diss_pca(
    ncomp = ncomp_by_opc(20)
  ),
  Yr = train_y,
  neighbors = neighbors_k(50),
  spike = forced_guests
)

# check the first 8 neighbors found in training for the 
# first 2 observations in testing
knn_spiked$neighbors[1:8, 1:2]

