## ----setup, include = FALSE---------------------------------------------------
library(formatR)
knitr::opts_chunk$set(
  collapse = TRUE, eval.after = "fig.cap"
)

## ----classdiagram, echo = FALSE, out.width = '20%', fig.align = 'right'-------
knitr::include_graphics("logo.jpg")

## ----eval = TRUE--------------------------------------------------------------
citation(package = "resemble")

## ----libraries, tidy = TRUE, message = FALSE----------------------------------
library(resemble)
library(prospectr)
library(magrittr)

## ---- tidy = FALSE, message = FALSE, results = 'hide'-------------------------
data(NIRsoil)
dim(NIRsoil)
str(NIRsoil)

## ----NIRsoil, tidy = FALSE, message = FALSE-----------------------------------
# obtain a numeric vector of the wavelengths at which spectra is recorded 
wavs <- NIRsoil$spc %>% colnames() %>% as.numeric()

# pre-process the spectra:
# - resample it to a resolution of 6 nm
# - use first order derivative
new_res <- 5
poly_order <- 1
window <- 5
diff_order <- 1

NIRsoil$spc_p <- NIRsoil$spc %>% 
  resample(wav = wavs, new.wav = seq(min(wavs), max(wavs), by = new_res)) %>% 
  savitzkyGolay(p = poly_order, w = window, m = diff_order)

## ----plotspectra, fig.cap = "Raw spectral absorbance data (top) and first derivative of the absorbance spectra (bottom).",  fig.cap.style = "Image Caption", fig.align = "center", fig.width = 7, fig.height = 7, echo = FALSE, fig.retina = 0.85----
old_par <- par("mfrow", "mar")
par(mfrow = c(2, 1), mar = c(4, 4, 1, 4))

new_wavs <- as.matrix(as.numeric(colnames(NIRsoil$spc_p)))
plot(range(wavs), range(NIRsoil$spc), col = NA,
     xlab = "",
     ylab = "Absorbance")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#EFBF4780")
grid(lty = 1, col = "#E47C4E80")
matlines(x = wavs, y = t(NIRsoil$spc), 
         lty = 1, col = "#5177A133")

plot(range(new_wavs), range(NIRsoil$spc_p), col = NA,
     xlab = "Wavelengths, nm",
     ylab = "1st derivative")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#EFBF4780")
grid(lty = 1, col = "#E47C4E80")
matlines(x = new_wavs, y = t(NIRsoil$spc_p), 
        lty = 1, col = "#5177A133")
par(old_par)

## ----eval = FALSE-------------------------------------------------------------
#  new_wavs <- as.matrix(as.numeric(colnames(NIRsoil$spc_p)))
#  
#  matplot(x = wavs, y = t(NIRsoil$spc),
#          xlab = "Wavelengths, nm",
#          ylab = "Absorbance",
#          type = "l", lty = 1, col = "#5177A133")
#  
#  matplot(x = new_wavs, y = t(NIRsoil$spc_p),
#          xlab = "Wavelengths, nm",
#          ylab = "1st derivative",
#          type = "l", lty = 1, col = "#5177A133")

## -----------------------------------------------------------------------------
# training dataset
training  <- NIRsoil[NIRsoil$train == 1, ]
# testing dataset
testing  <- NIRsoil[NIRsoil$train == 0, ]

## ---- results = 'hide'--------------------------------------------------------
# principal component (pc) analysis with the default 
# method (singular value decomposition) 
pca_tr <- ortho_projection(Xr = training$spc_p, method = "pca")

pca_tr

## ----plotpcsvariance, fig.cap = "Individual contribution to the explained variance for each component (left) and cumulative variance explained by the principal components (right).",  fig.cap.style = "Image Caption", fig.align = "center", fig.width = 7, fig.height = 3, fig.retina = 0.85----
plot(pca_tr, col = "#D42B08CC")

## ---- results = 'hide'--------------------------------------------------------
# principal component (pc) analysis with the default 
# NIPALS algorithm
pca_nipals_tr <- ortho_projection(Xr = training$spc_p,
                                  method = "pca.nipals")

pca_nipals_tr

## ---- results = 'hide', eval = FALSE------------------------------------------
#  # Partial Least Squares decomposition using
#  # Total carbon as side information
#  # (this might take some seconds)
#  pls_tr <- ortho_projection(Xr = training$spc_p,
#                             Yr = training$Ciso,
#                             method = "pls")
#  pls_tr

## ---- results = 'hide', eval = FALSE------------------------------------------
#  # This retains components that alone explain at least 5% of the original
#  # variation in training$spc_p
#  var_sel <-  list(method = "var", value = 0.05)
#  pca_tr_minvar5 <- ortho_projection(Xr = training$spc_p,
#                                     method = "pca",
#                                     pc_selection = var_sel)
#  
#  pca_tr_minvar5

## ---- results = 'hide', eval = FALSE------------------------------------------
#  # This retains components that together explain at least 90% of the original
#  # variation in training$spc_p
#  cumvar_sel <-  list(method = "cumvar", value = 0.90)
#  
#  pca_tr_cumvar90 <- ortho_projection(Xr = training$spc_p,
#                                      method = "pca",
#                                      pc_selection = cumvar_sel)
#  
#  pca_tr_cumvar90

## ---- results = 'hide'--------------------------------------------------------
# This uses optimal component selection
# variation in training$spc_p
optimal_sel <-  list(method = "opc", value = 40)
pca_tr_opc <- ortho_projection(Xr = training$spc_p,
                               Yr = training$Ciso,
                               method = "pca", 
                               pc_selection = optimal_sel)
pca_tr_opc

## ----pcrmsd, fig.cap = "Root mean squared difference between the samples and their corresponding nearest neighbors (for Total Carbon as side finormation) found by using dissimilarity matrices computed with different number of PCs.", fig.id = "plot_pcs_opc", fig.cap.style = "Image Caption", fig.align = "center", fig.width = 5, fig.height = 4, fig.retina = 0.85----
plot(pca_tr_opc, col = "#FF1A00CC")

## ----rmsdscatter, fig.cap = paste("Comparison between each sample and its corresponding nearest neighbor (in terms of  Total Carbon) when ", pca_tr_opc$n_components, "are used for dissimilarity matrix computations."), fig.id = "plot_pcs_opc2", fig.cap.style = "Image Caption", fig.align = "center", fig.width = 4, fig.height = 4, fig.retina = 0.85----
# compute the dissimilarity matrix using all the retained scores
pc_diss <- f_diss(pca_tr_opc$scores, diss_method = "mahalanobis")
# get the nearest neighbor for each sample
nearest_n <- apply(pc_diss, MARGIN = 1, FUN = function(x) order(x)[2])
# compute the RMSD
rmsd <- sqrt(mean((training$Ciso - training$Ciso[nearest_n])^2, na.rm = TRUE))
rmsd
# the RSMD for all the components is already available in 
# ...$opc_evaluation
pca_tr_opc$opc_evaluation[pca_tr_opc$n_components, , drop = FALSE]
plot(training$Ciso[nearest_n], 
     training$Ciso, 
     ylab = "Ciso of the nearest neighbor, %", xlab = "Ciso, %",
     col = "#D19C17CC", pch = 16)
grid()

## ---- results = 'hide'--------------------------------------------------------
# This uses manual component selection 
manual_sel <-  list(method = "manual", value = 9)
# PC
pca_tr_manual <- ortho_projection(Xr = training$spc_p,
                                  method = "pca", 
                                  pc_selection = manual_sel)
pca_tr_manual

# PLS
pls_tr_manual <- ortho_projection(Xr = training$spc_p,
                                  Yr = training$Ciso,
                                  method = "pls", 
                                  pc_selection = manual_sel)
pls_tr_manual

## ---- results = 'hide'--------------------------------------------------------
optimal_sel <-  list(method = "opc", value = 40)
# PLS
pls_tr_opc <- ortho_projection(Xr = training$spc_p,
                               Yr = training$Ciso,
                               method = "pls", 
                               pc_selection = optimal_sel,
                               scale = TRUE)
# the pls projection matrix
pls_tr_opc$projection_mat

pls_projected <- predict(pls_tr_opc, newdata = testing$spc_p)

# PC
pca_tr_opc <- ortho_projection(Xr = training$spc_p,
                               Yr = training$Ciso,
                               method = "pca", 
                               pc_selection = optimal_sel,
                               scale = TRUE)
# the pca projection matrix
t(pca_tr_opc$X_loadings)

pca_projected <- predict(pca_tr_opc, newdata = testing$spc_p)

## ---- results = 'hide', eval = FALSE------------------------------------------
#  optimal_sel <-  list(method = "opc", value = 40)
#  pca_tr_ts <- ortho_projection(Xr = training$spc_p,
#                                Xu = testing$spc_p,
#                                Yr = training$Ciso,
#                                method = "pca",
#                                pc_selection = optimal_sel,
#                                scale = TRUE)
#  plot(pca_tr_ts)

## ---- results = 'hide', eval = FALSE------------------------------------------
#  optimal_sel <-  list(method = "opc", value = 40)
#  pls_tr_ts <- ortho_projection(Xr = training$spc_p,
#                                Xu = testing$spc_p,
#                                Yr = training$Ciso,
#                                method = "pls",
#                                pc_selection = optimal_sel,
#                                scale = TRUE)
#  
#  # the same PLS projection model can be obtained with:
#  pls_tr_ts2 <- ortho_projection(Xr = training$spc_p[!is.na(training$Ciso),],
#                                 Yr = training$Ciso[!is.na(training$Ciso)],
#                                 method = "pls",
#                                 pc_selection = optimal_sel,
#                                 scale = TRUE)
#  
#  identical(pls_tr_ts$projection_mat, pls_tr_ts2$projection_mat)

## ---- results = 'hide', eval = FALSE------------------------------------------
#  optimal_sel <-  list(method = "opc", value = 40)
#  pls_multi_yr <- ortho_projection(Xr = training$spc_p,
#                                   Xu = testing$spc_p,
#                                   Yr = training[, c("Ciso", "Nt", "CEC")],
#                                   method = "pls",
#                                   pc_selection = optimal_sel,
#                                   scale = TRUE)
#  plot(pls_multi_yr)

## ---- results = 'hide', eval = FALSE------------------------------------------
#  pls_multi_yr$opc_evaluation

## ---- results = 'hide', eval = FALSE------------------------------------------
#  # for PC dissimilarity using the default settings
#  pcd <- dissimilarity(Xr = training$spc_p,
#                       diss_method = "pca")
#  dim(pcd$dissimilarity)
#  
#  # for PC dissimilarity using the optimized component selection method
#  pcd2 <- dissimilarity(Xr = training$spc_p,
#                        diss_method = "pca.nipals",
#                        Yr = training$Ciso,
#                        pc_selection = list("opc", 20),
#                        return_projection = TRUE)
#  dim(pcd2$dissimilarity)
#  pcd2$dissimilarity
#  pcd2$projection # the projection used to compute the dissimilarity matrix
#  
#  # for PLS dissimilarity
#  plsd <- dissimilarity(Xr = training$spc_p,
#                        diss_method = "pls",
#                        Yr = training$Ciso,
#                        pc_selection = list("opc", 20),
#                        return_projection = TRUE)
#  dim(plsd$dissimilarity)
#  plsd$dissimilarity
#  plsd$projection # the projection used to compute the dissimilarity matrix

## ---- results = 'hide', eval = TRUE-------------------------------------------
# For PC dissimilarity using the optimized component selection method
pcd_tr_ts <- dissimilarity(Xr = training$spc_p,
                           Xu = testing$spc_p,
                           diss_method = "pca.nipals",
                           Yr = training$Ciso,
                           pc_selection = list("opc", 20))
dim(pcd_tr_ts$dissimilarity)

# For PLS dissimilarity
plsd_tr_ts <- dissimilarity(Xr = training$spc_p,
                            Xu = testing$spc_p,
                            diss_method = "pls",
                            Yr = training$Ciso,
                            pc_selection = list("opc", 20))
dim(plsd_tr_ts$dissimilarity)

## ----localdiss, eval = TRUE---------------------------------------------------
# for localized PC dissimilarity using the optimized component selection method
# set the number of neighbors to retain
knn <- 200
local_pcd_tr_ts <- dissimilarity(Xr = training$spc_p,
                                 Xu = testing$spc_p,
                                 diss_method = "pca",
                                 Yr = training$Ciso,
                                 pc_selection = list("opc", 20),
                                 .local = TRUE, 
                                 pre_k = knn)
dim(local_pcd_tr_ts$dissimilarity)

# For PLS dissimilarity
local_plsd_tr_ts <- dissimilarity(Xr = training$spc_p,
                                  Xu = testing$spc_p,
                                  diss_method = "pls",
                                  Yr = training$Ciso,
                                  pc_selection = list("opc", 20),
                                  .local = TRUE, 
                                  pre_k = knn)
dim(local_plsd_tr_ts$dissimilarity)

# check the dissimilarity scores between the first two 
# observations in the testing dataset and the first 10 
# observations in the training dataset
local_plsd_tr_ts$dissimilarity[1:10, 1:2]


## ---- results = 'hide', eval = FALSE------------------------------------------
#  cd_tr <- dissimilarity(Xr = training$spc_p, diss_method = "cor")
#  dim(cd_tr$dissimilarity)
#  cd_tr$dissimilarity

## ---- results = 'hide', eval = TRUE-------------------------------------------
cd_tr_ts <- dissimilarity(Xr = training$spc_p,
                          Xu = testing$spc_p,
                          diss_method = "cor")
dim(cd_tr_ts$dissimilarity)
cd_tr_ts$dissimilarity

## ---- results = 'hide', eval = FALSE------------------------------------------
#  # a moving window correlation dissimilarity between training and testing
#  # using a window size of 19 spectral data points (equivalent to 95 nm)
#  cd_mw <- dissimilarity(Xr = training$spc_p,
#                         Xu = testing$spc_p,
#                         diss_method = "cor",
#                         ws = 19)
#  cd_mw$dissimilarity

## ---- results = 'hide', eval = FALSE------------------------------------------
#  # compute the dissimilarity between all the training observations
#  ed <- dissimilarity(Xr = training$spc_p, diss_method = "euclid")
#  ed$dissimilarity

## ---- results = 'hide', eval = FALSE------------------------------------------
#  # compute the dissimilarity between all the training observations
#  pre_time_resemble <- proc.time()
#  ed_resemble <- dissimilarity(Xr = training$spc_p, diss_method = "euclid")
#  post_time_resemble <- proc.time()
#  post_time_resemble - pre_time_resemble
#  
#  pre_time_stats <- proc.time()
#  ed_stats <- dist(training$spc_p, method = "euclid")
#  post_time_stats <- proc.time()
#  post_time_stats - pre_time_stats
#  
#  # scale the results of dist() based on the number of input columns
#  ed_stats_tr <- sqrt((as.matrix(ed_stats)^2)/ncol(training$spc_p))
#  ed_stats_tr[1:2, 1:3]
#  
#  # compare resemble and R stats results of Euclidean distances
#  ed_resemble$dissimilarity[1:2, 1:3]

## ---- results = 'hide', eval = FALSE------------------------------------------
#  # compute the dissimilarity between the training and testing observations
#  ed_tr_ts <- dissimilarity(Xr = training$spc_p,
#                            Xu = testing$spc_p,
#                            diss_method = "euclid")

## ---- results = 'hide', eval = FALSE------------------------------------------
#  # compute the dissimilarity between the training and testing observations
#  cosine_tr_ts <- dissimilarity(Xr = training$spc_p,
#                                Xu = testing$spc_p,
#                                diss_method = "cosine")
#  dim(cosine_tr_ts$dissimilarity)
#  cosine_tr_ts$dissimilarity

## ---- results = 'hide', eval = FALSE------------------------------------------
#  sid_tr_ts <- dissimilarity(Xr = training$spc_p,
#                             Xu = testing$spc_p,
#                             diss_method = "sid")
#  dim(sid_tr_ts$dissimilarity)
#  sid_tr_ts$dissimilarity

## ---- results = 'hide', eval = TRUE-------------------------------------------
# PC dissimilarity with default settings (variance-based 
# of components)
pcad <- dissimilarity(training$spc_p, diss_method = "pca", scale = TRUE)

# PLS dissimilarity with default settings (variance-based 
# of components)
plsd <- dissimilarity(training$spc_p, diss_method = "pls", Yr = training$Ciso,
                      scale = TRUE)

# PC dissimilarity with optimal selection of components
opc_sel <- list("opc", 30)
o_pcad <- dissimilarity(training$spc_p,
                        diss_method = "pca",
                        Yr = training$Ciso,
                        pc_selection = opc_sel, 
                        scale = TRUE)

# PLS dissimilarity with optimal selection of components
o_plsd <- dissimilarity(training$spc_p,
                        diss_method = "pls",
                        Yr = training$Ciso,
                        pc_selection = opc_sel, 
                        scale = TRUE)

# Correlation dissimilarity 
cd <- dissimilarity(training$spc_p, diss_method = "cor", scale = TRUE)

# Moving window correlation dissimilarity 
mcd <- dissimilarity(training$spc_p, diss_method = "cor", ws = 51, scale = TRUE)

# Euclidean dissimilarity 
ed <- dissimilarity(training$spc_p, diss_method = "euclid", scale = TRUE)

# Cosine dissimilarity 
cosd <- dissimilarity(training$spc_p, diss_method = "cosine", scale = TRUE)

# Spectral information divergence/dissimilarity 
sinfd <- dissimilarity(training$spc_p, diss_method = "sid", scale = TRUE)

## ---- results = 'hide', eval = TRUE-------------------------------------------
Ciso <- as.matrix(training$Ciso)
ev <- NULL
ev[["pcad"]] <- sim_eval(pcad$dissimilarity, side_info = Ciso)
ev[["plsd"]] <- sim_eval(plsd$dissimilarity, side_info = Ciso)
ev[["o_pcad"]] <- sim_eval(o_pcad$dissimilarity, side_info = Ciso)
ev[["o_plsd"]] <- sim_eval(o_plsd$dissimilarity, side_info = Ciso)
ev[["cd"]] <- sim_eval(cd$dissimilarity, side_info = Ciso)
ev[["mcd"]] <- sim_eval(mcd$dissimilarity, side_info = Ciso)
ev[["ed"]] <- sim_eval(ed$dissimilarity, side_info = Ciso)
ev[["cosd"]] <- sim_eval(cosd$dissimilarity, side_info = Ciso)
ev[["sinfd"]] <- sim_eval(sinfd$dissimilarity, side_info = Ciso)

## ---- results = 'hide', eval = TRUE, echo = FALSE-----------------------------
fig_cap <- paste("Comparison between observations and their corresponding 
                 nearest neighbor (1-NN)","observation in terms of Tocal Carbon (Ciso). The 1-NNs",
                 "are retrieved with the", 
                 "following dissimilarity merics:",
                 "pcad: PC dissimilarity with default settings (variance-based of components);",
                 "plsd: PLS dissimilarity with default settings (variance-based of components);",
                 "o_pcad: PC dissimilarity with optimal selection of components;",  
                 "o_plsd: PLS dissimilarity with optimal selection of components;",  
                 "cd: Correlation dissimilarity;",  
                 "mcd: Moving window correlation dissimilarity;",  
                 "ed: Euclidean dissimilarity;",  
                 "sinfd: Spectral information divergence/dissimilarity.")

## ----results = 'hide', eval = TRUE, echo = FALSE------------------------------
comparisons <- lapply(names(ev), 
                      FUN = function(x, label) {
                        irmsd <- round(x[[label]]$eval[1], 2)
                        ir <- round(x[[label]]$eval[2], 2)
                        if (label %in% c("o_pcad", "o_plsd")) {
                          label <- paste0("**", label, "**")
                          irmsd <- paste0("**", irmsd, "**")
                          ir <- paste0("**", ir, "**")
                        }
                        
                        data.frame(Measure = label, 
                                   RMSD = irmsd, 
                                   r =  ir)
                      },
                      x = ev)
comparisons

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  comparisons <- lapply(names(ev),
#                        FUN = function(x, label) {
#                          irmsd <- x[[label]]$eval[1]
#                          ir <- x[[label]]$eval[2]
#                          data.frame(Measure = label,
#                                     RMSD = irmsd,
#                                     r =  ir)
#                        },
#                        x = ev)
#  comparisons

## ----tcomparisons, eval = TRUE, echo = FALSE----------------------------------
knitr::kable(do.call("rbind", comparisons), 
             caption = paste("Root mean squared difference (RMSD)",
                             "and correlation coefficients",
                             "for between the observations and their", 
                             "corrresponding closest observations", 
                             "retrieved with the different dissimilarity 
                             methods."), 
             format = "simple", digits = 2, align = "l", padding = 2)

## ----eval = FALSE-------------------------------------------------------------
#  old_par <- par("mfrow")
#  par(mfrow = c(3, 3))
#  p <- sapply(names(ev),
#              FUN = function(x, label, labs = c("Ciso (1-NN), %", "Ciso, %")) {
#                xy <- x[[label]]$first_nn[,2:1]
#                plot(xy, xlab = labs[1], ylab = labs[2], col = "red")
#                title(label)
#                grid()
#                abline(0, 1)
#  
#              },
#              x = ev)
#  par(old_par)

## ----pcomparisons, fig.cap = paste(fig_cap), fig.cap.style = "Image Caption", fig.align = "center", fig.width = 8, fig.height = 8, echo = FALSE, fig.retina = 0.85----
old_par <- par("mfrow", "mar")

par(mfrow = c(3, 3), pch = 16, mar = c(4,4,4,4))
my_cols <- c("#750E3380", 
             "#C3BC6180", 
             "#FFE64F80", 
             "#EFBF4780", 
             "#E47C4E80", 
             "#F1A07380", 
             "#A1CFC480", 
             "#6B8EB580", 
             "#5177A180")
names(my_cols) <- sample(names(ev))
p <- sapply(names(ev), 
            FUN = function(x, 
                           label, 
                           labs = c("Ciso (1-NN), %", "Ciso, %"),
                           cols) {
              xy <- x[[label]]$first_nn[,2:1]
              plot(xy, xlab = labs[1], ylab = labs[2], col = cols[label])
              title(label)
              grid(col= "#80808080", lty = 1)
              abline(0, 1, col = "#FF1A0080")
            },
            x = ev,
            cols = my_cols)
par(old_par)

## ---- results = 'hide', eval = TRUE-------------------------------------------
knn_pc <- search_neighbors(Xr = training$spc_p, 
                           Xu = testing$spc_p, 
                           diss_method = "pca.nipals",
                           k = 50)

# matrix of neighbors
knn_pc$neighbors

# matrix of neighbor distances (dissimilarity scores)
knn_pc$neighbors_diss

# the index (in the training set) of the first two closest neighbors found in 
# training for the first observation in testing:
knn_pc$neighbors[1:2, 1, drop = FALSE]

# the distances of the two closest neighbors found in 
# training for the first observation in testing:
knn_pc$neighbors_diss[1:2, 1, drop = FALSE]

# the indices in training that fall in any of the 
# neighborhoods of testing
knn_pc$unique_neighbors

## ---- results = 'hide', eval = FALSE------------------------------------------
#  # using PC dissimilarity with optimal selection of components
#  knn_opc <- search_neighbors(Xr = training$spc_p,
#                              Xu = testing$spc_p,
#                              diss_method = "pca.nipals",
#                              Yr = training$Ciso,
#                              k = 50,
#                              pc_selection = list("opc", 20),
#                              scale = TRUE)
#  
#  # using PLS dissimilarity with optimal selection of components
#  knn_pls <- search_neighbors(Xr = training$spc_p,
#                              Xu = testing$spc_p,
#                              diss_method = "pls",
#                              Yr = training$Ciso,
#                              k = 50,
#                              pc_selection = list("opc", 20),
#                              scale = TRUE)
#  
#  # using correlation dissimilarity
#  knn_c <- search_neighbors(Xr = training$spc_p,
#                            Xu = testing$spc_p,
#                            diss_method = "cor",
#                            k = 50, scale = TRUE)
#  
#  # using moving window correlation dissimilarity
#  knn_mwc <- search_neighbors(Xr = training$spc_p,
#                              Xu = testing$spc_p,
#                              diss_method = "cor",
#                              k = 50,
#                              ws = 51, scale = TRUE)

## ---- results = 'hide', eval = FALSE------------------------------------------
#  # using localized PC dissimilarity with optimal selection of components
#  knn_local_opc <- search_neighbors(Xr = training$spc_p,
#                                    Xu = testing$spc_p,
#                                    diss_method = "pca.nipals",
#                                    Yr = training$Ciso,
#                                    k = 50,
#                                    pc_selection = list("opc", 20),
#                                    scale = TRUE,
#                                    .local = TRUE,
#                                    pre_k = 250)
#  
#  # using localized PLS dissimilarity with optimal selection of components
#  knn_local_opc <- search_neighbors(Xr = training$spc_p,
#                                    Xu = testing$spc_p,
#                                    diss_method = "pls",
#                                    Yr = training$Ciso,
#                                    k = 50,
#                                    pc_selection = list("opc", 20),
#                                    scale = TRUE,
#                                    .local = TRUE,
#                                    pre_k = 250)

## ---- results = 'hide', eval = TRUE-------------------------------------------
# a dissimilarity threshold
d_th <- 1

# the minimum number of observations required in each neighborhood
k_min <- 20

# the maximum number of observations allowed in each neighborhood
k_max <- 300

dnn_pc <- search_neighbors(Xr = training$spc_p, 
                           Xu = testing$spc_p, 
                           diss_method = "pca.nipals",
                           k_diss = d_th,
                           k_range = c(k_min, k_max), 
                           scale = TRUE)

# matrix of neighbors. The minimum number of indices is 20 (given by k_min)
# and the maximum number of indices is 300 (given by k_max).
# NAs indicate "not a neighbor"
dnn_pc$neighbors

# this reports how many neighbors were found for each observation in 
# testing using the input distance threshold (column n_k) and how 
# many were finally selected (column final_n_k)
dnn_pc$k_diss_info

# matrix of neighbor distances
dnn_pc$neighbors_diss

# the indices in training that fall in any of the 
# neighborhoods of testing
dnn_pc$unique_neighbors

## ---- knnhist, fig.cap = "Histogram of the original neighborhood sizes", fig.cap.style = "Image Caption", fig.align = "center", fig.width = 5, fig.height = 4, fig.retina = 0.80----
hist(dnn_pc$k_diss_info$final_n_k, 
     breaks = k_min,
     xlab = "Final neighborhood size", 
     main = "", col = "#EFBF47CC")

## ---- results = 'hide', eval = FALSE------------------------------------------
#  # the indices of the observations that we want to "invite" to every neighborhood
#  forced_guests <- c(1, 5, 8, 9)
#  
#  # using PC dissimilarity with optimal selection of components
#  knn_spiked <- search_neighbors(Xr = training$spc_p,
#                                 Xu = testing$spc_p,
#                                 diss_method = "pca.nipals",
#                                 Yr = training$Ciso,
#                                 k = 50,
#                                 spike = forced_guests,
#                                 pc_selection = list("opc", 20))
#  
#  # check the first 8 neighbors found in training for the
#  # first 2 observations in testing
#  knn_spiked$neighbors[1:8, 1:2]

## ----mblgif, fig.cap =  "Example of the main steps in memory-based learning for predicting a response variable in five different observations based on set of p-dimensional space.", echo = FALSE, out.width = '65%', fig.align = 'center', fig.retina = 0.85----
knitr::include_graphics("MBL.gif")

## ---- eval = TRUE-------------------------------------------------------------
# creates an object with instructions to build PLS models
my_plsr <- local_fit_pls(pls_c = 15)
my_plsr

# creates an object with instructions to build WAPLS models
my_waplsr <- local_fit_wapls(min_pls_c = 3, max_pls_c = 20)
my_waplsr
  
# creates an object with instructions to build GPR models
my_gpr <- local_fit_gpr()
my_gpr

## ---- eval = TRUE-------------------------------------------------------------
# create an object with instructions to conduct both validation types 
# "NNv" "local_cv"
two_val_control <- mbl_control(validation_type = c("NNv", "local_cv"),
                               number = 10,
                               p = 0.75)

## ---- results = 'hide', eval = TRUE-------------------------------------------
# define the dissimilarity method 
my_diss <- "cor"

# define the neighborhood sizes to test
my_ks <- seq(80, 200, by = 40)
  
# define how to use the dissimilarity information (ignore it)
ignore_diss <- "none"

# define the regression method to be used at each neighborhood 
my_waplsr <- local_fit_wapls(min_pls_c = 3, max_pls_c = 20)

# for the moment use only "NNv" validation (it will be faster)
nnv_val_control <- mbl_control(validation_type = "NNv")
  
# predict Total Carbon
# (remove missing values)
local_ciso <- mbl(
  Xr = training$spc_p[!is.na(training$Ciso),],
  Yr = training$Ciso[!is.na(training$Ciso)],
  Xu = testing$spc_p,
  k = my_ks,
  method = my_waplsr,
  diss_method = my_diss,
  diss_usage = ignore_diss,
  control = nnv_val_control,
  scale = TRUE
)

## ----localresultsciso, fig.cap = "MBL results for Total Carbon predictions using the LOCAL algorithm. NNv: nearest-neighbor cross-validation.", fig.align = "center", fig.width = 7.5, fig.height = 4, fig.retina = 0.85----
plot(local_ciso, main = "")
local_ciso 

## ---- results = 'hide', eval = TRUE, echo = FALSE-----------------------------
bestk <- which.min(local_ciso$validation_results$nearest_neighbor_validation$rmse)
bestk <- local_ciso$validation_results$nearest_neighbor_validation$k[bestk]

## ---- results = 'hide', eval = TRUE,  fig.show = 'hide'-----------------------
bki <- which.min(local_ciso$validation_results$nearest_neighbor_validation$rmse)
bk <- local_ciso$validation_results$nearest_neighbor_validation$k[bki]

# all the prediction results are stored in:
local_ciso$results

# the get_predictions function makes easier to retrieve the
# predictions from the previous object
ciso_hat <- as.matrix(get_predictions(local_ciso))[, bki]

## ---- eval = F----------------------------------------------------------------
#  # Plot predicted vs reference
#  plot(ciso_hat, testing$Ciso,
#       xlim = c(0, 14),
#       ylim = c(0, 14),
#       xlab = "Predicted Total Carbon, %",
#       ylab = "Total Carbon, %",
#       main = "LOCAL using argument k")
#  grid()
#  abline(0, 1, col = "red")

## ---- eval = TRUE-------------------------------------------------------------
# prediction RMSE:
sqrt(mean((ciso_hat - testing$Ciso)^2, na.rm = TRUE))

# squared R
cor(ciso_hat, testing$Ciso, use = "complete.obs")^2

## ---- results = 'hide', eval = TRUE,  fig.show = 'hide'-----------------------
# create a vector of dissimilarity thresholds to evaluate
# since the correlation dissimilarity will be used
# these thresholds need to be > 0 and <= 1 
dths <- seq(0.025, 0.3, by = 0.025)

# indicate the minimum and maximum sizes allowed for the neighborhood
k_min <- 30 
k_max <- 200 

local_ciso_diss <- mbl(
  Xr = training$spc_p[!is.na(training$Ciso),],
  Yr = training$Ciso[!is.na(training$Ciso)],
  Xu = testing$spc_p,
  k_diss = dths,
  k_range = c(k_min, k_max),
  method = my_waplsr,
  diss_method = my_diss,
  diss_usage = ignore_diss,
  control = nnv_val_control,
  scale = TRUE
)

## ---- eval = FALSE------------------------------------------------------------
#  plot(local_ciso_diss)

## ---- results = 'hide', eval = TRUE, echo = FALSE-----------------------------
bestd <- which.min(local_ciso_diss$validation_results$nearest_neighbor_validation$rmse)
bestd <- local_ciso_diss$validation_results$nearest_neighbor_validation$k[bestd]

## ---- eval = TRUE-------------------------------------------------------------
local_ciso_diss

## ---- results = 'hide', eval = TRUE,  fig.show = 'hide'-----------------------
# best distance threshold
bdi <- which.min(local_ciso_diss$validation_results$nearest_neighbor_validation$rmse)
bd <- local_ciso_diss$validation_results$nearest_neighbor_validation$k[bdi]

# predictions for the best distance
ciso_diss_hat <- as.matrix(get_predictions(local_ciso_diss))[, bdi]

## ---- eval = FALSE------------------------------------------------------------
#  # Plot predicted vs reference
#  plot(ciso_diss_hat, testing$Ciso,
#       xlim = c(0, 14),
#       ylim = c(0, 14),
#       xlab = "Predicted Total Carbon, %",
#       ylab = "Total Carbon, %",
#       main = "LOCAL using argument k_diss")
#  grid()
#  abline(0, 1, col = "red")

## ----addexamples, eval = TRUE, echo = FALSE-----------------------------------
abbr <- c("`local_cec`",
          "`pc_pred_cec`",
          "`pls_pred_cec`",
          "`local_gpr_cec`")

diss <- c("Correlation",
          "optimized PC",
          "optimized PLS",
          "optimized PC")

d_usage <- c("None",
             "Source of predictors",
             "None",
             "Source of predictors")

reg_m <- c("Weighted average PLS",
           "Weighted average PLS",
           "Weighted average PLS",
           "Gaussian process")

my_tab <-  as.data.frame(cbind(abbr, diss, d_usage, reg_m))
colnames(my_tab) <- c("Abreviation",
                      "Dissimilarity method",
                      "Dissimilarity usage",
                      "Local regression")

knitr::kable(my_tab, 
             caption = paste("Basic description of the different MBL", 
                             "configurations in the examples to predict",
                             "Cation Exhange Capacity (CEC)."), 
             format = "simple", align = "l", padding = 2)

## ---- results = 'hide', eval = TRUE,  fig.show = 'hide'-----------------------
# Lets define some methods:
my_wapls <- local_fit_wapls(2, 25)
k_min_max <- c(80, 200)

# use the LOCAL algorithm
# specific thresholds for cor dissimilarity
dth_cor <- seq(0.01, 0.3, by = 0.03)
local_cec <- mbl(
  Xr = training$spc_p[!is.na(training$CEC),],
  Yr = training$CEC[!is.na(training$CEC)],
  Xu = testing$spc_p,
  k_diss = dth_cor, 
  k_range = k_min_max,
  method = my_wapls,
  diss_method = "cor",
  diss_usage = "none",
  control = nnv_val_control,
  scale = TRUE
)

# use one where pca dissmilarity is used and the dissimilarity matrix  
# is used as source of additional predictors
# lets define first a an appropriate vector of dissimilarity thresholds 
# for the PC dissimilarity method
dth_pc <- seq(0.05, 1, by = 0.1)
pc_pred_cec <- mbl(
  Xr = training$spc_p[!is.na(training$CEC),],
  Yr = training$CEC[!is.na(training$CEC)],
  Xu = testing$spc_p,
  k_diss = dth_pc,
  k_range = k_min_max,
  method = my_wapls,
  diss_method = "pca",
  diss_usage = "predictors",
  control = nnv_val_control,
  scale = TRUE
)

# use one where PLS dissmilarity is used and the dissimilarity matrix  
# is used as source of additional predictors
pls_pred_cec <- mbl(
  Xr = training$spc_p[!is.na(training$CEC),],
  Yr = training$CEC[!is.na(training$CEC)],
  Xu = testing$spc_p,
  Yu = testing$CEC,
  k_diss = dth_pc,
  k_range = k_min_max,
  method = my_wapls,
  diss_method = "pls",
  diss_usage = "none",
  control = nnv_val_control,
  scale = TRUE
)

# use one where Gaussian process regression and pca dissmilarity are used 
# and the dissimilarity matrix  is used as source of additional predictors
local_gpr_cec <- mbl(
  Xr = training$spc_p[!is.na(training$CEC),],
  Yr = training$CEC[!is.na(training$CEC)],
  Xu = testing$spc_p,
  k_diss = dth_pc,
  k_range = k_min_max,
  method = local_fit_gpr(),
  diss_method = "pca",
  diss_usage = "predictors",
  control = nnv_val_control,
  scale = TRUE
)

## ---- eval = TRUE,  fig.show = 'hide'-----------------------------------------
# get the indices of the best results according to 
# nearest neighbor validation statistics
c_val_name <- "validation_results"
c_nn_val_name <- "nearest_neighbor_validation"

bi_local <- which.min(local_cec[[c_val_name]][[c_nn_val_name]]$rmse)
bi_pc_pred <- which.min(pc_pred_cec[[c_val_name]][[c_nn_val_name]]$rmse)
bi_pls_pred <- which.min(pls_pred_cec[[c_val_name]][[c_nn_val_name]]$rmse)
bi_local_gpr  <- which.min(local_gpr_cec[[c_val_name]][[c_nn_val_name]]$rmse)

preds <- cbind(get_predictions(local_cec)[, ..bi_local],
               get_predictions(pc_pred_cec)[, ..bi_pc_pred],
               get_predictions(pls_pred_cec)[, ..bi_pls_pred],
               get_predictions(local_gpr_cec)[, ..bi_local_gpr])

colnames(preds) <- c("local_cec", 
                     "pc_pred_cec", 
                     "pls_pred_cec", 
                     "local_gpr_cec")
preds <- as.matrix(preds)

# R2s
cor(testing$CEC, preds, use = "complete.obs")^2

#RMSEs
colMeans((preds - testing$CEC)^2, na.rm = TRUE)^0.5

## ---- eval = FALSE,  fig.show = 'hide'----------------------------------------
#  old_par <- par("mfrow", "mar")
#  
#  par(mfrow = c(2, 2))
#  plot(testing$CEC, preds[, 2],
#       xlab = "Predicted CEC, meq/100g",
#       ylab = "CEC, meq/100g", main = colnames(preds)[2])
#  abline(0, 1, col = "red")
#  
#  plot(testing$CEC, preds[, 3],
#       xlab = "Predicted CEC, meq/100g",
#       ylab = "CEC, meq/100g", main = colnames(preds)[3])
#  abline(0, 1, col = "red")
#  
#  plot(testing$CEC, preds[, 4],
#       xlab = "Predicted CEC, meq/100g",
#       ylab = "CEC, meq/100g", main = colnames(preds)[4])
#  abline(0, 1, col = "red")
#  par(old_par)

## ----mblcomparisons, fig.cap = "CEC prediction results for the different MBL configurations tested" , fig.cap.style = "Image Caption", fig.align = "center", fig.width = 8, fig.height = 8, echo = FALSE, fig.retina = 0.85----
old_par <- par("mfrow", "mar")

par(mfrow = c(2, 2), pch = 16, mar = c(4,4,4,4))
my_cols <- c("#D42B08CC",
             "#750E3380",
             "#EFBF4780", 
             "#5177A180")
names(my_cols) <- colnames(preds)

# R2s
r2s <- round(cor(testing$CEC, preds, use = "complete.obs")^2, 2)

#RMSEs
rmses <- round(colMeans((preds - testing$CEC)^2, na.rm = TRUE)^0.5,2)

names(r2s) <- colnames(preds)
names(rmses) <- colnames(preds)

p <- sapply(colnames(preds), 
            FUN = function(y,
                           yhats,
                           label, 
                           labs = c("Predicted CEC, meq/100g", "CEC, meq/100g"),
                           rsq,
                           rmse,
                           cols) {
              plot(x = yhats[,label], 
                   y = y, 
                   ylim = range(y, na.rm = T), xlim = range(y, na.rm = T),
                   xlab = labs[1], ylab = labs[2], 
                   col = cols[label])
              title(label)
              title(paste0("\n\n\n RMSE: ", rmse[label], "; Rsq: ", rsq[label]), cex.main=1)
              grid(col= "#80808080", lty = 1)
              abline(0, 1, col = "#FF1A0080")
            },
            yhats = preds,
            y = testing$CEC,
            cols = my_cols,
            rsq = r2s,
            rmse = rmses)
par(old_par)

## ---- eval = FALSE------------------------------------------------------------
#  # use Yu argument to validate the predictions
#  pc_pred_nt_yu <- mbl(
#    Xr = training$spc_p[!is.na(training$Nt),],
#    Yr = training$Nt[!is.na(training$Nt)],
#    Xu = testing$spc_p,
#    Yu = testing$Nt,
#    k = seq(40, 100, by = 10),
#    diss_usage = "none",
#    control = mbl_control(validation_type = "NNv"),
#    scale = TRUE
#  )
#  
#  pc_pred_nt_yu

## ---- eval = FALSE------------------------------------------------------------
#  # Running the mbl function using multiple cores
#  
#  # Execute with two cores, if available, ...
#  n_cores <- 2
#  
#  # ... if not then go with 1 core
#  if (parallel::detectCores() < 2) {
#    n_cores <- 1
#  }
#  
#  # Set the number of cores
#  library(doParallel)
#  clust <- makeCluster(n_cores)
#  registerDoParallel(clust)
#  
#  # Alternatively:
#  # library(doSNOW)
#  # clust <- makeCluster(n_cores, type = "SOCK")
#  # registerDoSNOW(clust)
#  # getDoParWorkers()
#  
#  pc_pred_nt <- mbl(
#    Xr = training$spc_p[!is.na(training$Nt),],
#    Yr = training$Nt[!is.na(training$Nt)],
#    Xu = testing$spc_p,
#    k = seq(40, 100, by = 10),
#    diss_usage = "none",
#    control = mbl_control(validation_type = "NNv"),
#    scale = TRUE
#  )
#  
#  # go back to sequential processing
#  registerDoSEQ()
#  try(stopCluster(clust))
#  
#  pc_pred_nt

