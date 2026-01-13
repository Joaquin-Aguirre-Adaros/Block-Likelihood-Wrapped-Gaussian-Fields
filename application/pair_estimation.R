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

datos <- read.csv("df.csv")

ang_col <- "ang_dir"
w_obs <- (as.numeric(datos[[ang_col]]) %% (2 * pi))

grid <- cbind(lon = as.numeric(datos$lon), lat = as.numeric(datos$lat))
n <- nrow(grid)
dd <- as.matrix(dist(grid))

# Todas las parejas (triángulo superior)
index_all <- which(upper.tri(dd, diag = FALSE), arr.ind = TRUE)

# Umbral ds
ds <- 0.2
d_pairs <- dd[cbind(index_all[, 1], index_all[, 2])]
index_ds <- index_all[d_pairs < ds, , drop = FALSE]

cat("============================================================\n")
cat("Datos:", n, "puntos\n")
cat("Total parejas:", nrow(index_all), "\n")
cat("Parejas con ds <", ds, ":", nrow(index_ds), "\n")
cat("============================================================\n\n")



# ---- Parámetros iniciales (optim en escala transformada) ----

inic <- c(log(1), log(0.1), 1)

kvec <- as.matrix(expand.grid(k1 = -1:1, k2 = -1:1))

# ---- Función negativa log composite (parejas) ----
cl_pairs <- function(par, w, eleccion, index_ds, dd, kvec) {
  s2 <- exp(par[1])
  alfa <- exp(par[2])
  mu <- 2 * atan(par[3]) + pi # en (0, 2*pi)

  suma <- 0

  for (ind in seq_len(nrow(index_ds))) {
    ii <- index_ds[ind, 1]
    jj <- index_ds[ind, 2]

    d_ij <- dd[ii, jj]
    dd_ij <- matrix(c(0, d_ij, d_ij, 0), 2, 2)

    cov_ij <- cov_f(dd_ij, s2, alfa, eleccion)

    # Jitter por estabilidad numérica (por si no es SPD)
    cov_ij <- cov_ij + diag(1e-10, 2)

    ww <- c(w[ii], w[jj])

    dens_total <- 0
    for (kk in seq_len(nrow(kvec))) {
      k1 <- kvec[kk, 1]
      k2 <- kvec[kk, 2]

      w_wrapped <- c(
        ww[1] - 2 * pi * k1 - mu,
        ww[2] - 2 * pi * k2 - mu
      )

      dens_total <- dens_total +
        mvtnorm::dmvnorm(w_wrapped, mean = c(0, 0), sigma = cov_ij)
    }

    suma <- suma - log(max(dens_total, 1e-300))
  }

  suma
}

# ---- Ajuste por parejas para una covarianza ----
fit_pair_model <- function(w, eleccion, inic, index_ds, dd, kvec) {
  timing <- system.time(
    opt <- optim(
      par = inic,
      fn = cl_pairs,
      w = w,
      eleccion = eleccion,
      index_ds = index_ds,
      dd = dd,
      kvec = kvec,
      control = list(reltol = 1e-16, maxit = 1e4, trace = FALSE)
    )
  )

  s2_hat <- exp(opt$par[1])
  alfa_hat <- exp(opt$par[2])
  mu_hat <- 2 * atan(opt$par[3]) + pi

  list(
    eleccion = eleccion,
    s2 = s2_hat,
    alfa = alfa_hat,
    mu = mu_hat,
    neglogcl = opt$value,
    logcl = -opt$value,
    conv = opt$convergence,
    time = timing[["elapsed"]]
  )
}

# ---- Estimar para las 3 covarianzas y elegir mejor ----
cat("Estimando modelos...\n\n")
fits <- lapply(1:3, function(m) {
  cat("Ajustando modelo", m, "...\n")
  fit_pair_model(w_obs, m, inic, index_ds, dd, kvec)
})

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
cat(sprintf(
  "s2=%.6f | alfa=%.6f | mu=%.6f | logCL=%.6f\n",
  best_fit$s2, best_fit$alfa, best_fit$mu, best_fit$logcl
))
cat("============================================================\n\n")

# ---- Preparar para Bootstrap: Calcular Cholesky ----
Sigma_full <- cov_f(dd, best_fit$s2, best_fit$alfa, best_fit$eleccion)
Sigma_full <- Sigma_full + diag(1e-10, n)

chol_cov <- tryCatch(
  chol(Sigma_full),
  error = function(e) {
    cat("Error en Cholesky, usando eigen-descomposición...\n")
    eig <- eigen(Sigma_full, symmetric = TRUE)
    eig$vectors %*% diag(sqrt(pmax(eig$values, 0))) %*% t(eig$vectors)
  }
)

cat("Matriz de covarianza calculada (", n, "x", n, ")\n")

# ---- Bootstrap en paralelo : simular y re-estimar ----
B <- 200L

RNGkind("L'Ecuyer-CMRG")
set.seed(42)

n_cores <- parallel::detectCores(logical = TRUE) - 1
cat("Usando", n_cores, "núcleos para bootstrap (mclapply/fork).\n\n")

worker_boot <- function(b) {
  z <- as.numeric(t(chol_cov) %*% rnorm(n) + best_fit$mu)
  w <- z %% (2 * pi)

  ft <- tryCatch(
    fit_pair_model(w, best_fit$eleccion, inic, index_ds, dd, kvec),
    error = function(e) {
      list(
        eleccion = best_fit$eleccion,
        s2 = NA, alfa = NA, mu = NA,
        neglogcl = NA, logcl = NA,
        conv = 999, time = 0
      )
    }
  )

  if (b %% 25 == 0) {
    cat(sprintf(
      "Bootstrap %d/%d | logCL=%.3f | conv=%d\n",
      b, B, ft$logcl, ft$conv
    ))
  }

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
boot_list <- parallel::mclapply(
  X = seq_len(B),
  FUN = worker_boot,
  mc.cores = n_cores,
  mc.preschedule = FALSE
)

# ---- Armar data.frames finales ----
est_bootstrap_pair <- do.call(rbind, lapply(boot_list, `[[`, "est"))
simul_bootstrap_pair <- do.call(rbind, lapply(boot_list, `[[`, "sim"))

cat("\nBootstrap terminado.\n")
cat("Resumen estimaciones (primeras 6 filas):\n")
print(head(est_bootstrap_pair))
cat("\nResumen logCL (bootstrap):\n")
print(summary(est_bootstrap_pair$logcl))
cat("\nConvergencias exitosas:", sum(est_bootstrap_pair$conv == 0, na.rm = TRUE), "/", B, "\n")

# ---- Guardar a CSV ----
write.csv(est_bootstrap_pair, "est_bootstrap_pair.csv", row.names = FALSE)
write.csv(simul_bootstrap_pair, "simul_bootstrap_pair.csv", row.names = FALSE)

cat("\nGuardado:\n")
cat(" - est_bootstrap_pair.csv\n")
cat(" - simul_bootstrap_pair.csv\n")
