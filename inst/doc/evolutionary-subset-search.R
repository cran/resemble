## -----------------------------------------------------------------------------
#| message: false
#| echo: false
library(resemble)
Sys.setenv(OMP_NUM_THREADS = 2)


## -----------------------------------------------------------------------------
#| message: false
library(prospectr)
data(NIRsoil)

# Preprocess spectra
NIRsoil$spc_pr <- savitzkyGolay(
  detrend(NIRsoil$spc, wav = as.numeric(colnames(NIRsoil$spc))),
  m = 1, p = 1, w = 7
)

# Missing values in the response are allowed
train_x <- NIRsoil$spc_pr[NIRsoil$train == 1, ]
train_y <- NIRsoil$Ciso[NIRsoil$train == 1]
test_x  <- NIRsoil$spc_pr[NIRsoil$train == 0, ]
test_y  <- NIRsoil$Ciso[NIRsoil$train == 0]

md <- dissimilarity(
  test_x,
  apply(test_x, 2, median),
  diss_method = diss_correlation(ws = 101, center = T, scale = T)
)

threshold <- 0.4


## -----------------------------------------------------------------------------
#| label: fig-diss-threshold
#| fig-cap: "Left: Dissimilarity of the test spectra to the median test spectrum. Samples above the threshold (dashed line and red) were excluded from the target set as a simple attempt to reduce spectral heterogeneity and potential non-linearities in the spectra-response relationship. Right: Spectra of the samples above (red, to be excluded) and below (blue, to be kept) the dissimilarity threshold."
#| fig-align: "center"
#| fig-width: 7.5
#| fig-height: 4.5
op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2, 1))

plot(
  seq_along(md$dissimilarity),
  md$dissimilarity, 
  ylim = c(0, 1), 
  col = "steelblue",
  pch = 16,
  cex = 1.5,
  xlab = "Sample index", 
  ylab = "Dissimilarity to median spectrum"
)
points(
  which(md$dissimilarity >= threshold),
  md$dissimilarity[md$dissimilarity >= threshold],
  col = "tomato",
  pch = 16,
  cex = 1.5
)
abline(h = threshold, col = "tomato", lty = 2, lwd = 2)
grid(lty = 1)

matplot(
  as.numeric(colnames(train_x)),
  t(test_x[md$dissimilarity >= threshold, ]),
  col = "tomato",
  ylab = "Preprocessed spectra",
  xlab = "Wavelengths, nm",
  lty = 1,
  type = "l"
)
grid(lty = 1)
matlines(
  as.numeric(colnames(train_x)),
  t(test_x[md$dissimilarity < threshold, ]),
  col = "steelblue",
  lty = 1,
  type = "l"
)

par(op)


## -----------------------------------------------------------------------------
keep <- md$dissimilarity < threshold

cat("Number of samples retained:", sum(keep))

test_x_test <- test_x[keep, ]
test_y_test <- test_y[keep]

# sample a very small subset to guide the search, 
# and use the rest for testing the final model
set.seed(1124)
ref_ind <- sample(sum(keep), 8)
test_x_ref <- test_x_test[ref_ind, ]
test_y_ref <- test_y_test[ref_ind]

test_x_test <- test_x_test[-ref_ind, ]
test_y_test <- test_y_test[-ref_ind]


## -----------------------------------------------------------------------------
set.seed(1124)

pls_model <- model(
  Xr = rbind(
    train_x[!is.na(train_y), ],
    test_x_ref[!is.na(test_y_ref), ]
  ), 
  Yr = c(
    train_y[!is.na(train_y)],
    test_y_ref[!is.na(test_y_ref)]
  ),
  fit_method = fit_pls(ncomp = 15, method = "simpls", scale = TRUE),
  control = model_control(validation_type = "lgo", number = 10)
)

best_ncomp <- which(pls_model$cv_results$optimal)

pred <- predict(pls_model, test_x_test, ncomp = best_ncomp)


## -----------------------------------------------------------------------------
reg_metrics <- function(obs, pred, na.rm = TRUE) {
  if (na.rm) {
    ok <- complete.cases(obs, pred)
    obs <- obs[ok]
    pred <- pred[ok]
  }
  c(RMSE = sqrt(mean((obs - pred)^2)), R2 = cor(obs, pred)^2)
}


## -----------------------------------------------------------------------------
reg_metrics(test_y_test, pred)


## -----------------------------------------------------------------------------
my_control <- gesearch_control(
  retain_by = "probability",
  tune = FALSE,
  number = 10L,
  p = 0.75
)


## -----------------------------------------------------------------------------
#| eval: true
#| echo: false 
gs <- gesearch(
  Xr = train_x[!is.na(train_y), ], # the spectra of the reference "set/library"
  Yr = train_y[!is.na(train_y)],   # the response of the reference "set/library"
  Xu = test_x_ref, # the available target domain spectra
  Yu = test_y_ref, # the available target domain response (accepts NAs)
  k = 50, b = 100, retain = 0.95,
  target_size = 120,
  fit_method = fit_pls(ncomp = 20, method = "simpls", scale = TRUE),
  optimization = c("reconstruction", "similarity", "response"),
  control = my_control,
  seed = 1124, 
  verbose = FALSE
)


## -----------------------------------------------------------------------------
#| eval: false
# gs <- gesearch(
#   Xr = train_x[!is.na(train_y), ], # the spectra of the reference "set/library"
#   Yr = train_y[!is.na(train_y)],   # the response of the reference "set/library"
#   Xu = test_x_ref, # the available target domain spectra
#   Yu = test_y_ref, # the available target domain response (accepts NAs)
#   k = 50, b = 100, retain = 0.95,
#   target_size = 120,
#   fit_method = fit_pls(ncomp = 20, method = "simpls", scale = TRUE),
#   optimization = c("reconstruction", "similarity", "response"),
#   control = my_control,
#   seed = 1124
# )


## -----------------------------------------------------------------------------
#| eval: false
# plot(gs)


## -----------------------------------------------------------------------------
#| label: fig-evolution-gs
#| echo: false 
#| fig-cap: "Evolution of the three weakness metrics (response, reconstruction and dissimilarity) used in the search."
#| fig-align: "center"
#| fig-width: 8.5
#| fig-height: 5.5
plot(gs, main = NA)


## -----------------------------------------------------------------------------
# the internal leave-group out CV results for the data found
gs$validation_results[[1]]$results$train


## -----------------------------------------------------------------------------
# the CV on the target samples that drove the search (note that this is not a proper test set, as these samples were used to guide the search, but it can provide some insight into the performance of the final model on the target domain)
gs$validation_results[[1]]$results$test


## -----------------------------------------------------------------------------
best_ncomp <- which.min(gs$validation_results[[1]]$results$train$rmse_cv)


## -----------------------------------------------------------------------------
pred_gs <- predict(gs, test_x_test)[[1]]

reg_metrics(test_y_test, pred_gs[, best_ncomp])


## -----------------------------------------------------------------------------
#| label: fig-pred-comparison
#| fig-cap: "Comparison of observed versus predicted total carbon values. Left: model fitted using all available training samples. Right: model fitted using the subset of samples selected by gesearch()."
#| fig-align: "center"
#| fig-width: 7.5
#| fig-height: 4.5
op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))

rng <- range(pred, pred_gs, test_y_test, na.rm = TRUE)

plot(
  pred, test_y_test,
  xlim = rng,
  ylim = rng,
  xlab = "Predicted Total Carbon, %",
  ylab = "Observed Total Carbon, %",
  main = "Model with all training \nsamples",
  cex = 1.5,
  pch = 16,
  col = rgb(0.5, 0.5, 0.5, 0.6)
)
grid(lty = 1)
abline(0, 1, col = "red")

plot(
  pred_gs[, best_ncomp], test_y_test,
  xlim = rng,
  ylim = rng,
  xlab = "Predicted Total Carbon, %",
  ylab = "Observed Total Carbon, %",
  main = "Model with gesearch-selected \nsamples",
  cex = 1.5,
  pch = 16,
  col = rgb(0.5, 0.5, 0.5, 0.6)
)
grid(lty = 1)
abline(0, 1, col = "red")

par(op)


## -----------------------------------------------------------------------------
#| eval: true
#| echo: false 
gs_nr <- gesearch(
  Xr = train_x[!is.na(train_y), ], # the spectra of the reference "set/library"
  Yr = train_y[!is.na(train_y)],   # the response of the reference "set/library"
  Xu = test_x_ref, # the available target domain spectra
  Yu = NULL, # NO available target domain response (accepts NAs)
  k = 50, b = 80, retain = 0.98,
  target_size = 160,
  fit_method = fit_pls(ncomp = 20, method = "simpls", scale = TRUE),
  optimization = c("reconstruction", "similarity"),
  control = my_control,
  seed = 1124, 
  verbose = FALSE
)


## -----------------------------------------------------------------------------
#| eval: false
# gs_nr <- gesearch(
#   Xr = train_x[!is.na(train_y), ], # the spectra of the reference "set/library"
#   Yr = train_y[!is.na(train_y)],   # the response of the reference "set/library"
#   Xu = test_x_ref, # the available target domain spectra
#   Yu = NULL, # NO available target domain response (accepts NAs)
#   k = 50, b = 80, retain = 0.98,
#   target_size = 160,
#   fit_method = fit_pls(ncomp = 20, method = "simpls", scale = TRUE),
#   optimization = c("reconstruction", "similarity"),
#   control = my_control,
#   seed = 1124
# )


## -----------------------------------------------------------------------------
best_ncomp_nr <- which.min(gs_nr$validation_results[[1]]$results$train$rmse_cv)
best_ncomp_nr
pred_gs_nr <- predict(gs_nr, test_x_test)[[1]]

reg_metrics(test_y_test, pred_gs_nr[, best_ncomp_nr])


## -----------------------------------------------------------------------------
#| eval: true
#| echo: false 
gs_nr2 <- gesearch(
  Xr = train_x[!is.na(train_y), ], # the spectra of the reference "set/library"
  Yr = train_y[!is.na(train_y)],   # the response of the reference "set/library"
  Xu = test_x_ref, # the available target domain spectra
  Yu = NULL, # NO available target domain response (accepts NAs)
  Yu_lims = c(0, 5), # the response range in the target domain
  k = 50, b = 80, retain = 0.98,
  target_size = 160,
  fit_method = fit_pls(ncomp = 20, method = "simpls", scale = TRUE),
  optimization = c("reconstruction", "similarity", "range"),
  control = my_control,
  seed = 1124, 
  verbose = FALSE
)


## -----------------------------------------------------------------------------
#| eval: false
# gs_nr2 <- gesearch(
#   Xr = train_x[!is.na(train_y), ], # the spectra of the reference "set/library"
#   Yr = train_y[!is.na(train_y)],   # the response of the reference "set/library"
#   Xu = test_x_ref, # the available target domain spectra
#   Yu = NULL, # NO available target domain response (accepts NAs)
#   Yu_lims = c(0, 5), # the response range in the target domain
#   k = 50, b = 80, retain = 0.98,
#   target_size = 160,
#   fit_method = fit_pls(ncomp = 20, method = "simpls", scale = TRUE),
#   optimization = c("reconstruction", "similarity", "range"),
#   control = my_control,
#   seed = 1124
# )


## -----------------------------------------------------------------------------
best_ncomp_nr2 <- which.min(gs_nr2$validation_results[[1]]$results$train$rmse_cv)
best_ncomp_nr2
pred_gs_nr2 <- predict(gs_nr2, test_x_test)[[1]]

reg_metrics(test_y_test, pred_gs_nr2[, best_ncomp_nr2])


## -----------------------------------------------------------------------------
#| label: parallel-gesearch
#| eval: false
# # Running gesearch() using multiple cores
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
# my_control <- gesearch_control(
#   retain_by = "probability",
#   tune = FALSE,
#   number = 10L,
#   p = 0.75
# )
# 
# gs_p <- gesearch(
#   Xr = train_x[!is.na(train_y), ],
#   Yr = train_y[!is.na(train_y)],
#   Xu = test_x_ref,
#   Yu = test_y_ref,
#   k = 50,
#   b = 100,
#   retain = 0.98,
#   target_size = 80,
#   fit_method = fit_pls(ncomp = 20, method = "simpls", scale = TRUE),
#   optimization = c("reconstruction", "similarity", "response"),
#   control = my_control,
#   seed = 1124,
#   pchunks = 2L
# )
# 
# # Go back to sequential processing
# registerDoSEQ()
# try(stopCluster(clust))
# 
# gs_p

