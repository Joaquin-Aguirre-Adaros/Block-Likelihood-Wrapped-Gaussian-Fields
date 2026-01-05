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
