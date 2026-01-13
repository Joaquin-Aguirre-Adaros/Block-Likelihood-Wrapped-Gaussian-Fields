library(rWind)
library(dplyr)
library(circular)
library(CircStats)
library(ggplot2)
library(dplyr)


data(wind.data)
datos <- wind.data

# Si wind.data viene como lista con 1 data.frame, descomenta esto:
# if (is.list(datos) && !is.data.frame(datos)) datos <- datos[[1]]

# 1) Quedarse con lon, lat y ang_dir (radianes en [0, 2*pi))
df <- datos %>%
  transmute(
    lon = as.numeric(lon),
    lat = as.numeric(lat),
    ang_dir = ((as.numeric(dir) %% 360) * 2 * pi / 360) %% (2 * pi)
  ) %>%
  filter(is.finite(lon), is.finite(lat), is.finite(ang_dir)) %>%
  filter(lat <= 41 & 35 <= lat)%>%
  filter(lon <= 0 & -6<=lon)



# 2) Largo de flecha automático en UNIDADES lon/lat
span <- max(diff(range(df$lon)), diff(range(df$lat)))
arrow_len <- 0.03 * span   

# Componentes
dx <- arrow_len * cos(df$ang_dir)
dy <- arrow_len * sin(df$ang_dir)

# 3) Plot base “bonito”
op <- par(mar = c(4, 4, 3, 1))
plot(df$lon, df$lat,
     type = "n",
     asp = 1,
     xlab = "Longitud",
     ylab = "Latitud",
     main = "Campo direccional (wind.data)")

# puntos suaves (da contexto espacial)
points(df$lon, df$lat, pch = 16, cex = 0.6, col = gray(0.6, alpha = 0.4))

# flechas: punta más chica (length en pulgadas) + palo más grueso (lwd)
arrows(
  x0 = df$lon,
  y0 = df$lat,
  x1 = df$lon + dx,
  y1 = df$lat + dy,
  length = 0.04,   # <-- esta es la clave: punta más pequeña (0.02–0.06)
  angle  = 25,
  lwd    = 1.1,    # palo más visible
  col    = "black"
)

box()
par(op)

# 4) Rosa (convención matemática: 0=Este, CCW)
ang_circ <- circular::circular(
  df$ang_dir,
  units = "radians",
  modulo = "2pi",
  zero = 0,
  rotation = "counter"
)

circular::rose.diag(
  ang_circ,
  bins = 36,
  col = "lightblue",
  prop = 2,
  main = "Rosa de direcciones (0=Este, CCW)"
)


ang_circ <- circular::circular(
  df$ang_dir,
  units = "radians",
  modulo = "2pi",
  zero = 0,              # 0 = Este
  rotation = "counter"   # antihorario (convención matemática)
)

mu <- circular::mean.circular(ang_circ)
sd_c <- circular::sd.circular(ang_circ)
rho <- circular::rho.circular(ang_circ)  # longitud media resultante (concentración)

cat("\n--- MÉTRICAS DIRECCIONALES ---\n")
cat(sprintf("Media direccional (rad): %.4f\n", mu))
cat(sprintf("Media direccional (deg): %.2f\n", (as.numeric(mu) * 180/pi) %% 360))
cat(sprintf("Desviación direccional (rad): %.4f\n", sd_c))
cat(sprintf("Concentración rho (0-1): %.4f\n", rho))

# (Opcional) estimación tipo Von Mises (kappa), útil como medida adicional de concentración
vm <- circular::mle.vonmises(ang_circ)
cat(sprintf("Kappa (Von Mises, opcional): %.4f\n", vm$kappa))

df$lon <- (df$lon +6)/6
df$lat <- (df$lat -35)/6


# --- vector de ángulos (radianes en [0, 2*pi)) ---
w <- df$ang_dir %% (2 * pi)

# --- distancias euclídeas en el plano lon/lat (en "grados") ---
coords <- as.matrix(df[, c("lon", "lat")])
dd <- as.matrix(dist(coords))

# --- parámetros de binning ---
tol_dist <- 0.025

# si quieres que el grid se adapte automáticamente al tamaño del dominio:
max_d <- max(dd[upper.tri(dd)], na.rm = TRUE)
h_values <- seq(0, 1, by = 0.05)

# --- correlación circular empírica por bin de distancia ---
corr_bins <- lapply(h_values, function(h) {
  pares <- which(dd >= (h - tol_dist) & dd <= (h + tol_dist), arr.ind = TRUE)
  
  if (is.null(pares) || nrow(pares) == 0) return(NULL)
  
  # eliminar duplicados (quedarse con triángulo superior)
  pares <- pares[pares[, 1] < pares[, 2], , drop = FALSE]
  if (nrow(pares) == 0) return(NULL)
  
  # circ.cor entrega 1 valor (correlación) para los vectores emparejados
  corr_h <- CircStats::circ.cor(w[pares[, 1]], w[pares[, 2]])
  
  data.frame(
    h = h,
    dist_prom = mean(dd[pares], na.rm = TRUE),
    n_pares = nrow(pares),
    corr = as.numeric(corr_h)
  )
})

corr_df <- bind_rows(corr_bins) %>%
  filter(is.finite(corr))

corr_df <- bind_rows(
  data.frame(h = 0, dist_prom = 0, n_pares = NA_integer_, corr = 1),
  corr_df
) %>%
  distinct(dist_prom, .keep_all = TRUE) %>%  # por si ya apareció ~0 por tolerancia
  arrange(dist_prom)

# --- plot: correlación circular empírica vs distancia ---
ggplot(corr_df, aes(x = dist_prom, y = corr)) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "gray60") +
  geom_hline(yintercept = 0.05, linewidth = 0.6, color = "red") +  
  geom_line(linewidth = 0.9, color = "#1f77b4") +
  geom_point(aes(size = n_pares), color = "#1f77b4", alpha = 0.85, na.rm = TRUE) +
  geom_point( 
    data = data.frame(dist_prom = 0, corr = 1),
    aes(x = dist_prom, y = corr),
    inherit.aes = FALSE,
    size = 2.8,
    color = "#1f77b4"
  ) +
  scale_size(range = c(1.5, 5), guide = "none") +
  coord_cartesian(ylim = c(-1, 1)) +
  theme_minimal() +
  labs(
    title = "Correlación circular empírica ",
    subtitle = paste0(
      "Tolerancia = ±", tol_dist,
      " (distancia en plano transformado). Línea roja: 0.05"
    ),
    x = "Distancia",
    y = "Correlación circular"
  )


# 2) Largo de flecha automático en UNIDADES lon/lat
span <- max(diff(range(df$lon)), diff(range(df$lat)))
arrow_len <- 0.03 * span   

# Componentes
dx <- arrow_len * cos(df$ang_dir)
dy <- arrow_len * sin(df$ang_dir)

# 3) Plot base “bonito”
op <- par(mar = c(4, 4, 3, 1))
plot(df$lon, df$lat,
     type = "n",
     asp = 1,
     xlab = "X",
     ylab = "Y",
     main = "Campo direccional Transformado")

# puntos suaves (da contexto espacial)
points(df$lon, df$lat, pch = 16, cex = 0.6, col = gray(0.6, alpha = 0.4))

# flechas: punta más chica (length en pulgadas) + palo más grueso (lwd)
arrows(
  x0 = df$lon,
  y0 = df$lat,
  x1 = df$lon + dx,
  y1 = df$lat + dy,
  length = 0.04,   # <-- esta es la clave: punta más pequeña (0.02–0.06)
  angle  = 25,
  lwd    = 1.1,    # palo más visible
  col    = "black"
)

box()
par(op)

readr::write_csv(df, "C:/Users/cacoa/Desktop/Aplicacion_tesis/df.csv")






