#' Interactive Plotting for MGWRSAR Models
#'
#' @description
#' Visualizes the results of an MGWRSAR model (coefficients, t-statistics, residuals, or fitted values)
#' using interactive plots via the \code{plotly} package.
#'
#' If a Coordinate Reference System (CRS) is provided (either via the \code{crs} argument or stored in the model object),
#' the function generates an interactive map. Otherwise, it generates a standard scatter plot of the values
#' against the coordinates.
#'
#' @param x An object of class \code{mgwrsar}.
#' @param type Character. The type of output to plot. Options are:
#' \itemize{
#'   \item \code{'coef'}: Spatially varying coefficients (default).
#'   \item \code{'t_coef'}: t-statistics of the coefficients (visualized with a significance threshold of 1.96).
#'   \item \code{'residuals'}: Model residuals.
#'   \item \code{'fitted'}: Fitted values.
#' }
#' @param var Character. The name of the variable (covariate) to plot. Required if \code{type} is \code{'coef'} or \code{'t_coef'}.
#' @param crs Numeric or character. The Coordinate Reference System (e.g., EPSG code) of the coordinates.
#' If \code{NULL}, the function attempts to retrieve the CRS from \code{x@my_crs}.
#' @param mypalette Character. The color palette to use for the points (e.g., "RdYlGn", "Viridis"). Default is "RdYlGn".
#' @param size Numeric. The size of the markers. Default is 5.
#' @param opacity Numeric. The opacity of the markers, between 0 and 1. Default is 0.8.
#' @param title Character. A custom title for the plot. If \code{NULL}, a default title is automatically generated.
#' @param show_legend Logical. Whether to display the legend. Default is \code{TRUE}.
#' @param ... Additional arguments passed to the internal plot function.
#'
#' @return A \code{plotly} object representing the interactive plot (or map).
#'
#' @seealso \code{\link{MGWRSAR}}
#'
#' @importFrom plotly plot_ly add_trace layout config
#' @importFrom sf st_as_sf st_crs st_transform st_coordinates
#' @method plot mgwrsar
#' @export
plot.mgwrsar <- function(x,
                         type = 'coef',
                         var = NULL,
                         crs = NULL,
                         mypalette = "RdYlGn",
                         size = 5,
                         opacity = 0.8,
                         title = NULL,
                         show_legend = TRUE,
                         n_time_steps = 10,
                         ...) {

  # ============================================================
  # 0. CHECKS & BACKWARD COMPATIBILITY
  # ============================================================
  dots <- list(...)

  deprecated_args <- c("fopacity", "nbins", "radius", "mytile", "myzoom",
                       "myresolution", "LayersControl", "myzoomControl",
                       "mytile2", "ScaleBar", "ScaleBarOptions")
  used_deprecated <- intersect(names(dots), deprecated_args)

  if (length(used_deprecated) > 0) {
    warning(sprintf("Deprecated arguments ignored: %s", paste(used_deprecated, collapse=", ")), call.=FALSE)
  }

  model <- x
  if (!inherits(model, "mgwrsar")) stop("A mgwrsar class object is needed.")

  is_gdt <- (model@Type == 'GDT' && !is.null(model@Z))

  if (!is.null(n_time_steps) && !is_gdt) {
    warning("n_time_steps ignored (Model is not Spatio-Temporal GDT).")
    n_time_steps <- NULL
  }

  # ============================================================
  # 1. DATA PREPARATION
  # ============================================================

  # --- CAS A : ANIMATION FLUIDE (Prédiction sur grille N x T) ---
  if (is_gdt && !is.null(n_time_steps) && type == 'B_coef') {
    # ============================================================
    # INSERTION : GESTION DU CYCLE TEMPOREL (MODULO)
    # ============================================================
      # 1. Récupération du kernel temporel (le 2ème si vecteur de 2, sinon le 1er)
      k_t <- if(length(model@kernels) > 1) model@kernels[2] else model@kernels[1]

      # 2. Parsing du nom (format attendu : "nom_type_periode", ex: "gauss_modulo_365")
      parts <- unlist(strsplit(k_t, "_"))

      # 3. Application de la transformation si "modulo" est détecté
      if (length(parts) >= 3) {
        period <- as.numeric(parts[3])

        if (!is.na(period) && period > 0) {
          # Transformation Modulo standard
          z_transformed <- model@Z %% period

          # Gestion des indices : Si le modulo donne 0, on le remplace souvent par la période
          # (ex: jour 365 %% 365 = 0 -> on remet 365 si les données sont en base 1)
          # On suppose ici que si le min > 0, c'est du 1-based index.
          if (min(model@Z, na.rm = TRUE) > 0) {
            z_transformed[z_transformed == 0] <- period
          }

          # Mise à jour locale de Z pour la suite du plot
          model@Z <- z_transformed

          # message(sprintf("Cycle temporel appliqué pour l'affichage : Modulo %s", period))
        }
      }


    if (is.null(var)) {
      if(is.matrix(model@Betav)) var <- colnames(model@Betav)[1]
      else stop("Argument 'var' is required.")
    }

    message(sprintf("Generating space-time grid predictions for %d time steps...", n_time_steps))

    # 1. Lieux UNIQUES
    coords_mat <- as.matrix(model@coords)
    u_coords <- unique(coords_mat)
    n_loc <- nrow(u_coords)

    # 2. Séquence temporelle
    t_min <- min(model@Z)
    t_max <- max(model@Z)
    t_seq <- seq(t_min, t_max, length.out = n_time_steps)

    # 3. Construction de la Grille (Expand Grid manuel)
    # Ordre strict : T1(Loc1..N), T2(Loc1..N)...

    # Temps répété (Blocs)
    new_coords_t <- rep(t_seq, each = n_loc)

    # Lieux répétés (Cycles)
    new_coords_s <- u_coords[rep(1:n_loc, times = n_time_steps), , drop = FALSE]

    # Matrice coords pour predict
    newdata_coords_st <- cbind(new_coords_s, new_coords_t)

    # 4. Dummy data
    newdata_dummy <- model@data[rep(1, nrow(newdata_coords_st)), , drop = FALSE]

    # 5. Prédiction
    #browser()
    B_pred <- predict(model,
                      newdata = newdata_dummy,
                      newdata_coords = newdata_coords_st,
                      type = "B_pred",method='model',beta_proj=TRUE)

    val_to_plot <- B_pred$Beta_proj_out[, var]

    # 6. DataFrame pour Plotly
    df_plot <- data.frame(
      Value = val_to_plot,
      raw_x = newdata_coords_st[,1],
      raw_y = newdata_coords_st[,2],
      Time = newdata_coords_st[,3]
    )

    # ID constant pour info-bulle (facultatif pour le mapping maintenant)
    #browser()
    df_plot$LocationID <- rep(1:n_loc, times = n_time_steps)
    df_plot$ID <- 1:nrow(df_plot)

    var_name <- paste("Pred:", var)
    if (is.null(title)) title <- paste("Spatio-Temporal Evolution of", var)

    # IMPORTANT : On ne supprime PAS les NA ici pour garder la symétrie des frames
    # df_plot <- df_plot[!is.na(df_plot$Value), ]  <-- SUPPRIMÉ

  } else {
    # --- CAS B : AFFICHAGE CLASSIQUE (Observed Data) ---
    val_to_plot <- NULL
    var_name <- ""
    is_t_coef <- (type == 't_coef')

    if (type == 't_coef') {
      if (is.null(var)) stop("Argument 'var' is required.")
      if (is.matrix(model@Betav)) val_to_plot <- model@Betav[, var] / model@sev[, var]
      else val_to_plot <- as.vector(model@Betav) / as.vector(model@sev)
      var_name <- paste("t-stat:", var)
    } else if (type == 'residuals') {
      val_to_plot <- as.vector(model@residuals)
      var_name <- "Residuals"
    } else if (type == 'fitted') {
      val_to_plot <- as.vector(model@fit)
      var_name <- "Fitted"
    } else { # coef
      if (is.null(var)) {
        if(is.matrix(model@Betav)) var <- colnames(model@Betav)[1]
        else stop("Argument 'var' is required.")
      }
      if (is.matrix(model@Betav)) val_to_plot <- model@Betav[, var]
      else val_to_plot <- as.vector(model@Betav)
      var_name <- var
    }

    df_plot <- data.frame(ID = 1:length(val_to_plot), Value = val_to_plot)
    coords <- model@coords
    df_plot$raw_x <- coords[,1]
    df_plot$raw_y <- coords[,2]

    if (is_gdt) {
      df_plot$Time <- model@Z
      # Ici, comme les données observées ne sont pas sur une grille régulière,
      # on filtre les NA car on ne peut pas garantir la symétrie de toute façon.
      df_plot <- df_plot[!is.na(df_plot$Value), ]
      # Tri indispensable
      df_plot <- df_plot[order(df_plot$Time), ]
    } else {
      df_plot <- df_plot[!is.na(df_plot$Value), ]
    }

    if (is.null(title)) title <- var_name
  }

  # ============================================================
  # 2. CRS & SPATIAL TRANSFORMATION
  # ============================================================
  target_crs <- crs
  if (is.null(target_crs) && .hasSlot(model, "my_crs") && !is.null(model@my_crs) && !identical(model@my_crs, NA)) {
    target_crs <- model@my_crs
  }

  is_geospatial <- !is.null(target_crs)

  if (is_geospatial) {
    df_sf <- sf::st_as_sf(df_plot, coords = c("raw_x", "raw_y"), remove = FALSE)
    sf::st_crs(df_sf) <- sf::st_crs(target_crs)
    df_sf <- sf::st_transform(df_sf, 4326)
    coords_wgs84 <- sf::st_coordinates(df_sf)
    df_plot$lon <- coords_wgs84[,1]
    df_plot$lat <- coords_wgs84[,2]
  }

  if(nrow(df_plot) == 0) stop("No data to plot.")

  # Helper for Tooltip
  make_hover <- function(val, vname, time_val=NULL) {
    txt <- paste0("<b>", vname, ":</b> ", round(val, 4))
    if (!is.null(time_val)) txt <- paste0(txt, "<br><b>Time:</b> ", round(time_val, 2))
    return(txt)
  }

  # ============================================================
  # 3. PLOTLY CONSTRUCTION
  # ============================================================
  p <- plotly::plot_ly()

  trace_type <- if (is_geospatial) 'scattermapbox' else 'scatter'
  x_col <- if (is_geospatial) NULL else ~raw_x
  y_col <- if (is_geospatial) NULL else ~raw_y
  lon_col <- if (is_geospatial) ~lon else NULL
  lat_col <- if (is_geospatial) ~lat else NULL

  frame_col <- if ("Time" %in% names(df_plot)) ~Time else NULL

  # --- AJOUT DE LA TRACE ---
  # CORRECTION CRITIQUE : Suppression de 'ids' pour le mode grille régulière
  # Le tri implicite (Ligne i Frame 1 -> Ligne i Frame 2) fonctionne mieux.

  p <- plotly::add_trace(
    p,
    data = df_plot,
    type = trace_type,
    mode = 'markers',
    x = x_col, y = y_col, lon = lon_col, lat = lat_col,

    frame = frame_col,
    # ids = ~LocationID,  <-- SUPPRIME POUR EVITER LE BUG "1 POINT"

    text = make_hover(df_plot$Value, var_name, if(!is.null(frame_col)) df_plot$Time else NULL),
    hoverinfo = "text",
    name = var_name,
    marker = list(
      size = size,
      opacity = opacity,
      color = ~Value,
      colors = mypalette,
      colorbar = list(title = var_name),
      line = list(width = 0.5, color = 'black')
    ),
    showlegend = FALSE
  )

  # ============================================================
  # 4. LAYOUT & ANIMATION CONFIG
  # ============================================================
  my_legend_style <- list(x = 0.01, y = 0.99, bgcolor = "rgba(255,255,255,0.8)")


  # --- 4.1 SUBTITLE GENERATION (Bandwidth Info) ---
  subtitle_text <- ""
  border=''

  # Helper to format bandwidth
  format_bw <- function(bw) {
    if(length(bw) > 1) return(paste0("Adaptive [", round(min(bw),1), "-", round(max(bw),1), "]"))
    return(as.character(round(bw, 2)))
  }

  if (model@Type == 'GDT') {
    # Spatio-Temporal Bandwidths
    # Try to find specific bandwidth for the variable if available (multiscale)
    h_s <- model@H
    h_t <- model@Ht



    # If variable specific H exists (named vector)
    if (!is.null(var) && !is.null(names(h_s)) && var %in% names(h_s)) {
      h_s <- h_s[var]
    } else if (length(h_s) > 1 && !is.null(var)) {
      # Fallback: if H corresponds to cols of Betav but is not named
      idx <- match(var, colnames(model@Betav))
      if(!is.na(idx) && idx <= length(h_s)) h_s <- h_s[idx]
    }

    if (!is.null(var) && !is.null(names(h_t)) && var %in% names(h_t)) {
      h_t <- h_t[var]
    } else if (length(h_t) > 1 && !is.null(var)) {
      idx <- match(var, colnames(model@Betav))
      if(!is.na(idx) && idx <= length(h_t)) h_t <- h_t[idx]
    } else if (length(h_t) == 0 && length(h_s) > 1) {
      # Fallback if Ht is empty but H has 2 elements (old structure)
      h_t <- h_s[2]
      h_s <- h_s[1]
    }
    if(h_s==max(model@V)) border="\n(Border only temporal case)" else if(h_t==max(model@Vt)) border="\n(Border only spatial case)"
    subtitle_text <- paste0("Bandwidths: Spatial (Hs) = ", format_bw(h_s),
                            " | Temporal (Ht) = ", format_bw(h_t))

  } else {
    # Spatial or Temporal Only
    h_val <- model@H
    if (!is.null(var) && !is.null(names(h_val)) && var %in% names(h_val)) {
      h_val <- h_val[var]
    } else if (length(h_val) > 1 && !is.null(var)) {
      idx <- match(var, colnames(model@Betav))
      if(!is.na(idx) && idx <= length(h_val)) h_val <- h_val[idx]
    }

    lbl <- if(model@Type == 'T') "Temporal (Ht)" else "Spatial (Hs)"
    subtitle_text <- paste0("Bandwidth: ", lbl, " = ", format_bw(h_val))
  }

  # Add subtitle to title via HTML line break
  title <- paste0(title, "<br><sup>", subtitle_text, border, "</sup>")



  if (is_geospatial) {
    # Zoom auto
    lon_rng <- range(df_plot$lon, na.rm = TRUE)
    lat_rng <- range(df_plot$lat, na.rm = TRUE)
    target_span <- max(diff(lon_rng), diff(lat_rng)) * 1.1
    if(target_span == 0) target_span <- 0.01
    my_zoom <- max(0, min(22, log2(360 / target_span)))

    p <- plotly::layout(
      p,
      title = list(text = title, y = 0.98, x = 0.05, xanchor = "left"),
      mapbox = list(
        style = "carto-positron",
        zoom = my_zoom,
        center = list(lat = mean(lat_rng), lon = mean(lon_rng))
      ),
      legend = my_legend_style,
      margin = list(l = 0, r = 0, t = 40, b = 0)
    )
  } else {
    # 2D Layout
    rng_x <- range(df_plot$raw_x, na.rm = TRUE)
    rng_y <- range(df_plot$raw_y, na.rm = TRUE)
    span_x <- diff(rng_x); if(span_x==0) span_x <- 1
    span_y <- diff(rng_y); if(span_y==0) span_y <- 1
    margin_x <- span_x * 0.05
    margin_y <- span_y * 0.05

    p <- plotly::layout(
      p,
      title = list(text = title, y = 0.98),
      xaxis = list(title = "X", zeroline = FALSE,
                   range = c(rng_x[1] - margin_x, rng_x[2] + margin_x)),
      yaxis = list(title = "Y", zeroline = FALSE, scaleanchor = "x",
                   range = c(rng_y[1] - margin_y, rng_y[2] + margin_y)),
      legend = my_legend_style,
      margin = list(l = 50, r = 0, t = 50, b = 50)
    )
  }

  # Animation Options
  if (!is.null(frame_col)) {
    p <- plotly::animation_opts(
      p,
      frame = 1000,
      transition = 0, # Pas de transition floue
      redraw = TRUE   # Force le redessin complet (crucial pour éviter les fantômes)
    )

    p <- plotly::animation_slider(
      p,
      currentvalue = list(prefix = "Time: ", font = list(color = "red"))
    )
  }

  p <- plotly::config(p, scrollZoom = TRUE)
  return(p)
}
