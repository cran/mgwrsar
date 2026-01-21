#' Helper function to extract coordinates and data from spatial objects
#' Handles sf, sp, and raw dataframes with optional CRS via control list.
#' @noRd
get_SFdata <- function(env = parent.frame()) {

  # --- 1. Retrieve objects from parent environment ---
  if (!exists("data", envir = env)) stop("Object 'data' not found.")

  local_data <- get("data", envir = env)

  # Handle coords (optional if data is spatial)
  local_coords <- if (exists("coords", envir = env)) get("coords", envir = env) else NULL

  # Handle control (to retrieve manual CRS)
  local_control <- if (exists("control", envir = env)) get("control", envir = env) else list()
  user_crs <- local_control$crs # Can be NULL

  # Initialize final CRS to NULL (safe default for ANY slot)
  local_crs <- NULL

  # --- 2. Spatial Management Logic ---

  # >>> CASE 1: 'sf' Object
  if (inherits(local_data, "sf")) {
    if (!requireNamespace("sf", quietly = TRUE)) stop("Package 'sf' required.")

    # Try to retrieve CRS from object
    obj_crs <- sf::st_crs(local_data)

    # If object has valid CRS, keep it. Otherwise, check if user provided one.
    if (!is.na(obj_crs)) {
      local_crs <- obj_crs
    } else if (!is.null(user_crs)) {
      # Object is sf but has no projection, apply user's one
      local_crs <- user_crs
      # Optional: could try to convert to st_crs(user_crs) here to standardize
    }

    if (is.null(local_coords)) {
      local_coords <- sf::st_coordinates(local_data)
      if(ncol(local_coords) > 2) local_coords <- local_coords[, 1:2]
    }
    local_data <- sf::st_drop_geometry(local_data)

    # >>> CASE 2: 'sp' Object
  } else if (inherits(local_data, "Spatial")) {
    if (!requireNamespace("sp", quietly = TRUE)) stop("Package 'sp' required.")

    # Handle sp CRS
    input_proj <- sp::proj4string(local_data)

    if (!is.na(input_proj)) {
      local_crs <- sp::CRS(input_proj)
    } else if (!is.null(user_crs)) {
      local_crs <- user_crs
    }

    if (is.null(local_coords)) local_coords <- sp::coordinates(local_data)
    if (.hasSlot(local_data, "data")) local_data <- local_data@data

    # >>> CASE 3: Classical Data.frame + Coords
  } else {
    # User input via control$crs has priority here
    if (!is.null(user_crs)) {
      local_crs <- user_crs
    }
    # If user_crs is NULL, local_crs remains NULL (no known projection)
  }

  # --- 3. Final Safety Checks ---

  if (is.null(local_coords)) stop("Argument 'coords' must be provided if 'data' is not a spatial object.")
  local_coords <- as.matrix(local_coords)
  if(!is.numeric(local_coords)) stop("Coordinates must be numeric.")

  # --- 4. Re-assign to parent environment ---

  assign("data", local_data, envir = env)
  assign("coords", local_coords, envir = env)
  assign("my_crs", local_crs, envir = env)
}
