## -----------------------------------------------------------------------------
#| echo: false
Sys.setenv(OMP_NUM_THREADS = 2)


## -----------------------------------------------------------------------------
#| eval: true
citation(package = "resemble")


## -----------------------------------------------------------------------------
#| label: libraries
#| message: false
library(resemble)
library(prospectr)


## -----------------------------------------------------------------------------
#| message: false
#| results: hide
data(NIRsoil)
dim(NIRsoil)
str(NIRsoil)


## -----------------------------------------------------------------------------
#| label: NIRsoil
#| message: false
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


## -----------------------------------------------------------------------------
#| label: fig-plotspectra
#| fig-cap: "Raw spectral absorbance data (top) and first derivative of the absorbance spectra (bottom)."
#| fig-align: center
#| fig-width: 7
#| fig-height: 7
#| echo: false
# Hex sticker palette#
bg_dark <- NA #"#0F172A"      # background
blue <- "#3B82F6"         # border / spiral high
amber <- "#F59E0B"        # spiral low
amber_light <- "#FBBF24"  # text accent
slate <- "#64748B"        # muted elements

old_par <- par("mfrow", "mar", "bg")
par(mfrow = c(2, 1), mar = c(4, 4, 1, 4), bg = bg_dark)

new_wavs <- as.matrix(as.numeric(colnames(NIRsoil$spc_p)))

text_col <- "black"#"white"

# Plot 1: Raw spectra
plot(range(wavs), range(NIRsoil$spc), col = NA,
     xlab = "",
     ylab = "Absorbance",
     col.lab = text_col, col.axis = text_col)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#1E293B")
grid(lty = 1, col = "#334155")
matlines(x = wavs, y = t(NIRsoil$spc), 
         lty = 1, col = paste0(blue, "33"))

# Plot 2: First derivative
plot(range(new_wavs), range(NIRsoil$spc_p), col = NA,
     xlab = "Wavelengths, nm",
     ylab = expression(d(detrended~A)/d*lambda),
     col.lab = text_col, col.axis = text_col)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#1E293B")
grid(lty = 1, col = "#334155")
matlines(x = new_wavs, y = t(NIRsoil$spc_p), 
         lty = 1, col = paste0(amber, "33"))

par(old_par)


## -----------------------------------------------------------------------------
train_x <- NIRsoil$spc_pr[NIRsoil$train == 1, ]
train_y <- NIRsoil$Ciso[NIRsoil$train == 1]

test_x  <- NIRsoil$spc_pr[NIRsoil$train == 0, ]
test_y  <- NIRsoil$Ciso[NIRsoil$train == 0]

