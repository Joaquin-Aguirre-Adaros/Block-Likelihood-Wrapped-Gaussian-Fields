set.seed(42)

alpha <- c(0.1, 0.15,0.20)
n <- 100
n_sim <- 200


wrapped_simulations <- function(n_sim, n, alpha) {
  lista <- list()
  s2 <- 1
  mu <- pi
  point_x <- runif(n)
  point_y <- runif(n)
  grid <- cbind(point_x, point_y)
  dd <- as.matrix(dist(grid))
  for (k in 1:3) {
    for (j in 1:length(alpha)) {
      alfa <- alpha[j]
      for (i in 1:n_sim) {
        covmat <- cov_f(dd, s2, alfa,k) 
        cc <- chol(covmat)
        z <- t(cc) %*% rnorm(n) + mu
        w <- z %% (2 * pi)
        data <- data.frame(x = grid[, 1], y = grid[, 2], value = w)
        lista[[i + (k-1)*(n_sim*length(alpha)) +n_sim * (j - 1)]] <- data
    }
    }
  }
  return(lista)
}

simul <- wrapped_simulations(n_sim, n, alpha)

nombre_archivo <- "simul.rds"

# 2. Guardamos en la instancia (Linux)
saveRDS(simul, file = nombre_archivo)

# 3. Confirmación
print(paste("Archivo guardado exitosamente en la nube en:", getwd(), "/", nombre_archivo))
