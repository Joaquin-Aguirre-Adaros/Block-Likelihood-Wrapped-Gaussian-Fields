# ============================================================
# Wrapped-Gaussian field: Pairwise composite likelihood + bootstrap
# - Selecciona mejor covarianza entre eleccion=1:3 por log-CL
# - ds = 0.2
# - Genera 200 datasets sintéticos y re-estima en paralelo
# - Guarda: est_bootstrap_pair.csv y simul_bootstrap_pair.csv
# ============================================================

library(parallel)
library(mvtnorm)

set.seed(42)

# ---- 0) Cargar datos observados ----
datos <- read.csv("df.csv")

stopifnot(all(c("lon", "lat") %in% names(datos)))

# Ángulo observado (radianes en [0, 2*pi))
ang_col <- if ("ang_dir" %in% names(datos)) "ang_dir" else if ("value" %in% names(datos)) "value" else NA
if (is.na(ang_col)) stop("No encontré columna de ángulo. Usa 'ang_dir' o 'value' en df.csv.")

w_obs <- (as.numeric(datos[[ang_col]]) %% (2 * pi))

# ---- 1) Grid y distancias (plano) ----
grid <- cbind(lon = as.numeric(datos$lon), lat = as.numeric(datos$lat))
n <- nrow(grid)
dd <- as.matrix(dist(grid))

# todas las parejas (triángulo superior)
index_all <- which(upper.tri(dd, diag = FALSE), arr.ind = TRUE)

# umbral ds
ds <- 0.2
d_pairs <- dd[cbind(index_all[, 1], index_all[, 2])]
index_ds <- index_all[d_pairs < ds, , drop = FALSE]

cat("============================================================\n")
cat("Datos:", n, "puntos\n")
cat("Total parejas:", nrow(index_all), "\n")
cat("Parejas con ds <", ds, ":", nrow(index_ds), "\n")
cat("============================================================\n\n")

# ---- 2) Check cov_f ----
if (!exists("cov_f")) {
  stop("No existe cov_f en tu entorno. Debes tener definida cov_f(distmat, s2, alfa, eleccion).")
}

# ---- 3) Parámetros iniciales (optim en escala transformada) ----
# par = (log(s2), log(alfa), p3) con mu = 2*atan(p3) + pi
inic <- c(log(1), log(0.1), 1)

# grilla de envoltura kvec (9 combinaciones)
kvec <- as.matrix(expand.grid(list(-1:1, -1:1)))

# ---- 4) Función negativa log composite (parejas) ----
neglog_cl_pairs <- function(par, w, eleccion, index_ds, dd, kvec) {
  s2   <- exp(par[1])
  alfa <- exp(par[2])
  mu   <- 2 * atan(par[3]) + pi  # en (0, 2*pi)

  suma <- 0

  for (ind in seq_len(nrow(index_ds))) {
    ii <- index_ds[ind, 1]
    jj <- index_ds[ind, 2]

    d_ij <- dd[ii, jj]
    dd_ij <- matrix(c(0, d_ij, d_ij, 0), 2, 2)

    cov_ij <- cov_f(dd_ij, s2, alfa, eleccion)

    ww <- c(w[ii], w[jj])

    # mezcla de envolturas
    w_aug <- ww - 2 * pi * kvec - rep(mu, 2)  # 9x2
    dens <- mvtnorm::dmvnorm(w_aug, mean = c(0, 0), sigma = cov_ij)

    suma <- suma - log(max(sum(dens), 1e-300))
  }

  suma
}

# ---- 5) Ajuste por parejas para una covarianza ----
fit_pair_model <- function(w, eleccion, inic, index_ds, dd, kvec) {
  timing <- system.time(
    opt <- optim(
      par = inic,
      fn = neglog_cl_pairs,
      w = w,
      eleccion = eleccion,
      index_ds = index_ds,
      dd = dd,
      kvec = kvec,
      control = list(reltol = 1e-16, maxit = 1e4, trace = FALSE)
    )
  )

  s2_hat   <- exp(opt$par[1])
  alfa_hat <- exp(opt$par[2])
  mu_hat   <- 2 * atan(opt$par[3]) + pi

  list(
    eleccion = eleccion,
    s2 = s2_hat,
    alfa = alfa_hat,
    mu = mu_hat,
    neglogcl = opt$value,
    logcl = -opt$value,                 # función objetivo (mayor es mejor)
    conv = opt$convergence,
    time = timing[["elapsed"]]
  )
}

# ---- 6) Estimar para las 3 covarianzas y elegir mejor ----
fits <- lapply(1:3, function(m) fit_pair_model(w_obs, m, inic, index_ds, dd, kvec))

for (ft in fits) {
  cat(">>> Modelo eleccion =", ft$eleccion, "\n")
  cat(sprintf("    s2   = %.6f\n", ft$s2))
  cat(sprintf("    alfa = %.6f\n", ft$alfa))
  cat(sprintf("    mu   = %.6f rad\n", ft$mu))
  cat(sprintf("    neglogCL = %.6f\n", ft$neglogcl))
  cat(sprintf("    logCL    = %.6f  (OBJ)\n", ft$logcl))
  cat(sprintf("    conv=%d | time=%.2fs\n\n", ft$conv, ft$time))
}

best_idx <- which.max(sapply(fits, `[[`, "logcl"))
best_fit <- fits[[best_idx]]

cat("============================================================\n")
cat("MEJOR MODELO (por mayor logCL): eleccion =", best_fit$eleccion, "\n")
cat(sprintf("s2=%.6f | alfa=%.6f | mu=%.6f | logCL=%.6f\n",
            best_fit$s2, best_fit$alfa, best_fit$mu, best_fit$logcl))
cat("============================================================\n\n")

# ---- 7) Preparar simulación (mismo grid) con mejor modelo ----
# Covarianza completa + Cholesky (con jitter por si hace falta)
covmat <- cov_f(dd, best_fit$s2, best_fit$alfa, best_fit$eleccion)

# jitter robusto si chol falla
eps <- 1e-10
ok <- FALSE
for (k in 1:8) {
  try_chol <- try(chol(covmat + diag(eps, n)), silent = TRUE)
  if (!inherits(try_chol, "try-error")) {
    chol_cov <- try_chol
    ok <- TRUE
    break
  }
  eps <- eps * 10
}
if (!ok) stop("No pude hacer Cholesky de la covarianza (ni con jitter). Revisa cov_f/parametrización.")

B <- 200L

# inicial para bootstrap: partir desde el estimado (más estable)
p3_init <- tan((best_fit$mu - pi) / 2)
# cap por si mu cae muy cerca de pi ± pi
p3_init <- max(min(p3_init, 1e6), -1e6)
inic_boot <- c(log(best_fit$s2), log(best_fit$alfa), p3_init)

# ---- 8) Bootstrap en paralelo: simular y re-estimar ----
n_cores <- min(max(1L, parallel::detectCores() - 1L), B)
cat("Usando", n_cores, "núcleos para bootstrap.\n\n")

cl <- makeCluster(n_cores)
on.exit(stopCluster(cl), add = TRUE)

clusterEvalQ(cl, library(mvtnorm))

clusterExport(
  cl,
  varlist = c(
    "n", "grid", "dd", "index_ds", "kvec",
    "cov_f", "neglog_cl_pairs", "fit_pair_model",
    "best_fit", "chol_cov", "inic_boot"
  ),
  envir = environment()
)

clusterSetRNGStream(cl, iseed = 42)

worker_boot <- function(b) {
  # Simular campo gaussiano latente y envolver
  z <- as.numeric(t(chol_cov) %*% rnorm(n) + best_fit$mu)
  w <- z %% (2 * pi)

  # Re-estimar (mismo modelo ganador)
  ft <- fit_pair_model(w, best_fit$eleccion, inic_boot, index_ds, dd, kvec)

  # devolver estimación + simulación (en formato largo)
  list(
    est = data.frame(
      sim_id = b,
      eleccion = ft$eleccion,
      s2 = ft$s2,
      alfa = ft$alfa,
      mu = ft$mu,
      neglogcl = ft$neglogcl,
      logcl = ft$logcl,
      conv = ft$conv,
      time = ft$time
    ),
    sim = data.frame(
      sim_id = b,
      lon = grid[, "lon"],
      lat = grid[, "lat"],
      value = w
    )
  )
}

cat("Corriendo bootstrap (B=200)...\n")
boot_list <- parLapplyLB(cl, X = seq_len(B), fun = worker_boot)

# ---- 9) Armar data.frames finales ----
est_bootstrap_pair <- do.call(rbind, lapply(boot_list, `[[`, "est"))
simul_bootstrap_pair <- do.call(rbind, lapply(boot_list, `[[`, "sim"))

cat("\nBootstrap terminado.\n")
cat("Resumen estimaciones (primeras 6 filas):\n")
print(head(est_bootstrap_pair))
cat("\nMejores/peores logCL (bootstrap):\n")
print(summary(est_bootstrap_pair$logcl))

# ---- 10) Guardar a CSV ----
write.csv(est_bootstrap_pair, "est_bootstrap_pair.csv", row.names = FALSE)
write.csv(simul_bootstrap_pair, "simul_bootstrap_pair.csv", row.names = FALSE)

cat("\nGuardado:\n")
cat(" - est_bootstrap_pair.csv\n")
cat(" - simul_bootstrap_pair.csv\n")
