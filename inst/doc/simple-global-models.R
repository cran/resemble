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
# PLS with 15 components
fit_pls(ncomp = 15)

# SIMPLS with scaling
fit_pls(ncomp = 15, method = "simpls", scale = TRUE)

# mPLS with scaling
fit_pls(ncomp = 15, method = "mpls")

# GPR with default noise variance
fit_gpr()


## -----------------------------------------------------------------------------
#| label: data-prep
#| message: false
library(prospectr)

wavs <- as.numeric(colnames(NIRsoil$spc))

# Preprocess: detrend + first derivative
NIRsoil$spc_pr <- savitzkyGolay(
  detrend(NIRsoil$spc, wav = wavs),
  m = 1, p = 1, w = 7
)

# Split data
train_x <- NIRsoil$spc_pr[NIRsoil$train == 1, ]
train_y <- NIRsoil$Ciso[NIRsoil$train == 1]
test_x  <- NIRsoil$spc_pr[NIRsoil$train == 0, ]
test_y  <- NIRsoil$Ciso[NIRsoil$train == 0]

# Remove missing values
ok_train <- !is.na(train_y)
ok_test  <- !is.na(test_y)
train_x <- train_x[ok_train, ]
train_y <- train_y[ok_train]
test_x  <- test_x[ok_test, ]
test_y  <- test_y[ok_test]


## -----------------------------------------------------------------------------
#| label: fit-pls

set.seed(1124) # guarantee same CV splits for all methods
pls_mod <- model(
 Xr = train_x,
 Yr = train_y,
 fit_method = fit_pls(ncomp = 12, method = "mpls", scale = TRUE),
 control = model_control(validation_type = "lgo", number = 10)
)

pls_mod


## -----------------------------------------------------------------------------
#| label: fit-gpr
set.seed(1124) # guarantee same CV splits for all methods
gpr_mod <- model(
 Xr = train_x,
 Yr = train_y,
 fit_method = fit_gpr(noise_variance = 0.5),
 control = model_control(validation_type = "lgo", number = 10)
)

gpr_mod


## -----------------------------------------------------------------------------
#| label: predict
# Predict using optimal number of components (from CV)
pls_pred <- predict(pls_mod, newdata = test_x)

# Predict with GPR
gpr_pred <- predict(gpr_mod, newdata = test_x)

# Compare predictions
data.frame(
  observed = test_y,
  pls = pls_pred[, which(pls_mod$cv_results$optimal)],
  gpr = as.vector(gpr_pred)
) |> head(10)


## -----------------------------------------------------------------------------
#| label: validation
# Function to compute stats
eval_pred <- function(obs, pred) {
  data.frame(
    rmse = sqrt(mean((obs - pred)^2)),
    r2 = cor(obs, pred)^2
  )
}


## -----------------------------------------------------------------------------
## PLS evaluation
eval_pred(test_y, pls_pred[, which(pls_mod$cv_results$optimal)])


## -----------------------------------------------------------------------------
## GPR evaluation
eval_pred(test_y, gpr_pred)

