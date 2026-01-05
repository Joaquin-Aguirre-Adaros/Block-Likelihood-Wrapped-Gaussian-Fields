# Paquetes (mínimos para este script)
library(parallel)
library(mvtnorm)

set.seed(42)

# Configuración
nucleos <- 115
simul_total <- readRDS("simul.rds")

# Grid y distancias (comunes para todas las simulaciones)
grid <- cbind(simul_total[[1]]$x, simul_total[[1]]$y)
dd <- as.matrix(dist(grid))
n <- nrow(grid)

# Índices de todas las parejas posibles (triángulo superior)
index_all <- which(upper.tri(matrix(nrow = n, ncol = n), diag = FALSE), arr.ind = TRUE)

# Parámetros iniciales
inic <- c(log(1), log(0.1), 1)

# Precalcular grilla de envoltura (kvec) una vez
kvec <- as.matrix(expand.grid(list(-1:1, -1:1)))

# Función de covarianza
cov_f <- function(h, s2, alfa, eleccion) {
  if (eleccion == 1) {
    return(s2 * exp(-3 * h / alfa))
  }
  if (eleccion == 2) {
    return(s2 * exp(-4.7619 * h / alfa) * (1 + 4.7619 * h / alfa))
  }
  if (eleccion == 3) {
    return(s2 * exp(-5.91865 * h / alfa) *
      (1 + 5.91865 * h / alfa + (5.91865^2 / 3) * h^2 / alfa^2))
  }
  stop("`eleccion` debe ser 1, 2 o 3.")
}

# Worker: estima para una simulación i dado un escenario (ds, eleccion, index_ds)
est_pair <- function(i, ds, eleccion, index_ds) {
  w <- simul_total[[i]]$value

  pcl_pairs <- function(p) {
    p1 <- exp(p[1])
    p2 <- exp(p[2])
    p3 <- 2 * atan(p[3]) + pi

    suma <- 0

    for (ind in seq_len(nrow(index_ds))) {
      idx_i <- index_ds[ind, 1]
      idx_j <- index_ds[ind, 2]

      d_ij <- dd[idx_i, idx_j]
      dd_ij <- matrix(c(0, d_ij, d_ij, 0), 2, 2)
      cov_ij <- cov_f(dd_ij, p1, p2, eleccion)

      ww <- c(w[idx_i], w[idx_j])

      w_augmented <- ww - 2 * pi * kvec - rep(p3, 2)
      aux <- mvtnorm::dmvnorm(w_augmented, mean = c(0, 0), sigma = cov_ij)

      suma <- suma - log(max(sum(aux), 1e-300))
    }

    suma
  }

  timing <- system.time(
    estimacion <- optim(
      par = inic,
      fn = pcl_pairs,
      control = list(reltol = 1e-16, maxit = 1e4, trace = FALSE)
    )
  )

  # <<< AQUÍ estaba el HTML roto en tu script >>>
  c(
    s2 = exp(estimacion$par[1]),
    alfa = exp(estimacion$par[2]),
    mu = 2 * atan(estimacion$par[3]) + pi,
    time = timing[["elapsed"]]
  )
}

# ===== Paralelización: cluster una vez =====
cl <- makeCluster(nucleos)

clusterEvalQ(cl, library(mvtnorm))
clusterExport(
  cl,
  varlist = c("simul_total", "dd", "inic", "kvec", "cov_f", "est_pair"),
  envir = environment()
)

# (Opcional pero recomendado para reproducibilidad)
clusterSetRNGStream(cl, iseed = 42)

run_scenario <- function(ds, eleccion, sim_ids, out_file) {
  # Filtrar parejas por ds una vez por escenario
  d_pairs <- dd[cbind(index_all[, 1], index_all[, 2])]
  index_ds <- index_all[d_pairs < ds, , drop = FALSE]

  # Exportar index_ds (cambia por escenario)
  clusterExport(cl, "index_ds", envir = environment())

  res <- parLapplyLB(
    cl,
    X = sim_ids,
    fun = function(i, ds, eleccion) est_pair(i, ds = ds, eleccion = eleccion, index_ds = index_ds),
    ds = ds,
    eleccion = eleccion
  )

  res_df <- do.call(rbind, res)
  saveRDS(res_df, file = out_file)
  invisible(res_df)
}

# ===== Escenarios (mismo diseño que tu script) =====
# ds = 0.1
run_scenario(0.10, 1, 1:600,     "est_pair_1_01.rds")
run_scenario(0.10, 2, 601:1200,  "est_pair_2_01.rds")
run_scenario(0.10, 3, 1201:1800, "est_pair_3_01.rds")

# ds = 0.15
run_scenario(0.15, 1, 1:600,     "est_pair_1_015.rds")
run_scenario(0.15, 2, 601:1200,  "est_pair_2_015.rds")
run_scenario(0.15, 3, 1201:1800, "est_pair_3_015.rds")

# ds = 0.2
run_scenario(0.20, 1, 1:600,     "est_pair_1_02.rds")
run_scenario(0.20, 2, 601:1200,  "est_pair_2_02.rds")
run_scenario(0.20, 3, 1201:1800, "est_pair_3_02.rds")

stopCluster(cl)
