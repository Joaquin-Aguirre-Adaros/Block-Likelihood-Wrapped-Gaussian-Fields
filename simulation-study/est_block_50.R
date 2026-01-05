library(parallel)
library(mvtnorm)

set.seed(42)

# ---- Inputs ----
nucleos <- 115  # ideal: detectCores() - 1
simul_total <- simul  

# ---- Grid / distancias ----
grid <- cbind(simul_total[[1]]$x, simul_total[[1]]$y)
dd <- as.matrix(dist(grid))

# ---- K-means blocks ----
nb <- 50
km <- kmeans(grid, centers = nb)
centers <- km$centers
d_c <- as.matrix(dist(centers))
Index <- which(upper.tri(matrix(nrow = nb, ncol = nb), diag = FALSE), arr.ind = TRUE)

# ---- Parámetros iniciales ----
inic <- c(log(1), log(0.1), 1)

kvec_cache <- new.env(parent = emptyenv())
get_kvec <- function(n) {
  key <- as.character(n)
  if (exists(key, envir = kvec_cache, inherits = FALSE)) {
    return(get(key, envir = kvec_cache, inherits = FALSE))
  }
  lista <- rep(list(-1:1), n)
  kv <- as.matrix(expand.grid(lista))
  assign(key, kv, envir = kvec_cache)
  kv
}


est_block <- function(sim_id, ds, eleccion) {
  w <- simul_total[[sim_id]]$value

  pcl_block <- function(p) {
    p1 <- exp(p[1])                  # varianza
    p2 <- exp(p[2])                  # rango
    p3 <- 2 * atan(p[3]) + pi        # media

    suma <- 0

    for (ind in seq_len(nrow(Index))) {
      b1 <- Index[ind, 1]
      b2 <- Index[ind, 2]

      if (d_c[b1, b2] >= ds) next

      sqi <- which(km$cluster == b1)
      sqj <- which(km$cluster == b2)

      idx <- c(sqi, sqj)
      ww <- w[idx]
      nloc <- length(ww)

      dd_ij <- dd[idx, idx, drop = FALSE]
      cov_ij <- cov_f(dd_ij, p1, p2, eleccion)

      kvec <- get_kvec(nloc)

      # Matriz (#kvec x nloc): cada fila es un "wrap"
      w_augmented <- sweep(2 * pi * kvec, 2, ww, FUN = function(k, wv) wv - k)
      w_augmented <- sweep(w_augmented, 2, rep(p3, nloc), FUN = "-")

      aux <- mvtnorm::dmvnorm(w_augmented, mean = rep(0, nloc), sigma = cov_ij)

      suma <- suma - log(max(sum(aux), 1e-300))
    }

    suma
  }

  timing <- system.time({
    estimacion <- optim(
      par = inic,
      fn = pcl_block,
      control = list(reltol = 1e-16, maxit = 1e4, trace = FALSE)
    )
  })

  c(
    s2   = exp(estimacion$par[1]),
    alfa = exp(estimacion$par[2]),
    mu   = 2 * atan(estimacion$par[3]) + pi,
    time = timing[["elapsed"]]
  )
}


cl <- makeCluster(nucleos, type = "SOCK")
clusterEvalQ(cl, library(mvtnorm))

clusterExport(
  cl,
  varlist = c(
    "simul_total", "dd", "nb", "km", "d_c", "Index", "inic",
    "cov_f", "est_block", "get_kvec", "kvec_cache"
  ),
  envir = environment()
)

clusterSetRNGStream(cl, iseed = 42)

run_scenario <- function(ds, eleccion, sim_ids, out_file) {
  res <- parLapplyLB(
    cl,
    X = sim_ids,
    fun = function(id) est_block(sim_id = id, ds = ds, eleccion = eleccion)
  )
  res_mat <- do.call(rbind, res)
  saveRDS(res_mat, file = out_file)
  invisible(res_mat)
}

# ds = 0.1
run_scenario(0.10, 1, 1:600,     "est_block50_1_01.rds")
run_scenario(0.10, 2, 601:1200,  "est_block50_2_01.rds")
run_scenario(0.10, 3, 1201:1800, "est_block50_3_01.rds")

# ds = 0.15
run_scenario(0.15, 1, 1:600,     "est_block50_1_015.rds")
run_scenario(0.15, 2, 601:1200,  "est_block50_2_015.rds")
run_scenario(0.15, 3, 1201:1800, "est_block50_3_015.rds")

# ds = 0.2
run_scenario(0.20, 1, 1:600,     "est_block50_1_02.rds")
run_scenario(0.20, 2, 601:1200,  "est_block50_2_02.rds")
run_scenario(0.20, 3, 1201:1800, "est_block50_3_02.rds")

stopCluster(cl)
