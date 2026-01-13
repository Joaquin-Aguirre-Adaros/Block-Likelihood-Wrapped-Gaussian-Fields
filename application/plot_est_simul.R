# ---- 1) Cargar simulaciones ----
simul <- readRDS("C:/Users/cacoa/Desktop/Tesis/simul.rds")

# ---- 2) Librerías ----
library(ggplot2)
library(circular)

# ---- 3) Función para visualizar un campo direccional ----
plot_directional_field <- function(simul,
                                   i,
                                   x_col = "x",
                                   y_col = "y",
                                   angle_col = "value",
                                   arrow_length = 0.05,
                                   bins = 36,
                                   bw = 10,
                                   title_prefix = "Campo direccional") {
  # Extraer data
  df <- simul[[i]]
  
  # Validación mínima
  needed <- c(x_col, y_col, angle_col)
  if (!all(needed %in% names(df))) {
    stop(
      "Faltan columnas. Se esperaban: ",
      paste(needed, collapse = ", "),
      "\nColumnas disponibles: ",
      paste(names(df), collapse = ", ")
    )
  }
  
  # Ángulo en radianes (asumido en [0, 2*pi))
  angle_rad <- df[[angle_col]] %% (2 * pi)
  
  # Componentes para flechas (convención matemática: 0 = +x, CCW)
  df$dx <- arrow_length * cos(angle_rad)
  df$dy <- arrow_length * sin(angle_rad)
  
  # Objeto circular (útil para rosa/densidad)
  # template = "geographics" se usa para los gráficos circulares (rosa),
  # pero las flechas siguen cos/sin "matemático".
  df$angle_circular <- circular::circular(
    angle_rad,
    units = "radians",
    template = "geographics"
  )
  
  # --- 3.1) Mapa con flechas ---
  p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(alpha = 0.25, size = 1, color = "gray50") +
    geom_segment(
      aes(
        xend = .data[[x_col]] + dx,
        yend = .data[[y_col]] + dy
      ),
      arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
      color = "black",
      linewidth = 0.6,
      alpha = 0.85
    ) +
    coord_fixed() +
    theme_minimal() +
    labs(
      title = paste0(title_prefix, " - Simulación ", i),
      x = "Coordenada X",
      y = "Coordenada Y"
    )
  
  print(p)
  
  # --- 3.2) Rosa de direcciones ---
  circular::rose.diag(
    df$angle_circular,
    bins = bins,
    col = "lightblue",
    prop = 2,
    main = paste0("Rosa de direcciones - ", title_prefix, " ", i)
  )
  
  # --- 3.3) Densidad circular ---
  densidad <- circular::density.circular(df$angle_circular, bw = bw)
  plot(
    densidad,
    ylim = c(-1.1, 1.5),
    main = paste0("Densidad circular - ", title_prefix, " ", i),
    xlab = "Dirección (radianes)",
    ylab = "Densidad",
    col = "blue"
  )
  
  invisible(df$angle_circular)
}

# ---- 4) Función para calcular estadísticas circulares ----
calculate_circular_stats <- function(circular_data, title = "Simulación") {
  cat("--- ESTADÍSTICAS CIRCULARES (", title, ") ---\n", sep = "")
  cat("(Todos los valores en radianes)\n\n")
  
  media_circular <- circular::mean.circular(circular_data)
  cat(sprintf("Media circular: %.3f rad\n", media_circular))
  
  desviacion_circular <- circular::sd.circular(circular_data)
  cat(sprintf("Desviación estándar circular: %.3f rad\n", desviacion_circular))
  
  concentracion <- circular::rho.circular(circular_data)
  cat(sprintf("Concentración circular (ρ): %.3f\n", concentracion))
  
  rayleigh_test <- circular::rayleigh.test(circular_data)
  cat(sprintf("Prueba de Rayleigh (p-value): %.4f\n", rayleigh_test$p.value))
  
  if (rayleigh_test$p.value < 0.05) {
    cat("→ La distribución NO es uniforme (existe dirección predominante)\n")
  } else {
    cat("→ La distribución es uniforme (no existe dirección predominante)\n")
  }
  
  cat("\n")
  invisible(list(
    mean = media_circular,
    sd = desviacion_circular,
    rho = concentracion,
    rayleigh = rayleigh_test
  ))
}

# ---- 5) Ejemplo de uso ----
i <- 1800

ang_circ <- plot_directional_field(
  simul = simul,
  i = i,
  x_col = "x",
  y_col = "y",
  angle_col = "value",
  arrow_length = 0.05,
  bins = 36,
  bw = 10,
  title_prefix = "Campo direccional"
)

calculate_circular_stats(ang_circ, title = paste("Simulación", i))