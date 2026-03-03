# Setup Functions -------------------------------------------------------------

#' Generate test data for GD and GDT cases
#' 
#' @return List with GD and GDT datasets (300 obs each)
setup_test_data_full <- function(n=300,lambda=NULL,config_beta = 'default',config_snr = 0.9) {
  # GD case (spatial only)
  simuGD <- simu_multiscale(
    n = n, 
    lambda=lambda,
    myseed = 1, 
    config_snr = 0.9,
    config_beta = config_beta, 
    config_eps = 'normal'
  )
  
  # GDT case (spatiotemporal)
  simuGDT <- simu_multiscale(
    n = n, 
    myseed = 1, 
    config_snr = 0.9,
    config_beta = 'spatiotemp', 
    config_eps = 'normal'
  )
  
  TRUE_DIST=as.matrix(dist(simuGD$coords))
  coords_3col = cbind(simuGDT$coords, simuGDT$mydata$time)
  TRUE_DIST_t=as.matrix(dist(coords_3col[,3]))
  TRUE_DIST_tc <- pmin(TRUE_DIST_t %% 365,365 - (TRUE_DIST_t %% 365))
  list(
    GD = list(
      TRUE_DIST=TRUE_DIST,
      mydata = simuGD$mydata,
      coords = simuGD$coords,
      TRUEBETA = as.matrix(simuGD$mydata[, c('Beta1', 'Beta2', 'Beta3', 'Beta4')]),
      formula = as.formula('Y ~ X1 + X2 + X3')
    ),
    GDT = list(
      TRUE_DIST=  TRUE_DIST,
      TRUE_DIST_t=TRUE_DIST_t,
      TRUE_DIST_tc=TRUE_DIST_tc,
      mydata = simuGDT$mydata,
      coords = simuGDT$coords,
      coords_3col = coords_3col,  # coords avec time en col 3
      time = simuGDT$mydata$time,
      TRUEBETA = as.matrix(simuGDT$mydata[, c('Beta1', 'Beta2', 'Beta3', 'Beta4')]),
      formula = as.formula('Y ~ X1 + X2 + X3')
    )
  )
}

#' Extract subset of 30 observations for faster testing
#' 
#' @param data_list Output from setup_test_data_full()
#' @return Same structure with 30 observations
get_subset_30obs <- function(data_list) {
  idx <- 1:30
  
  list(
    GD = list(
      TD=data_list$GD$TRUE_DIST[idx,idx ],
      mydata = data_list$GD$mydata[idx, ],
      coords = data_list$GD$coords[idx, ],
      TRUEBETA = data_list$GD$TRUEBETA[idx, ],
      formula = data_list$GD$formula
    ),
    GDT = list(
      TD=data_list$GDT$TRUE_DIST[idx,idx ],
      TD_t=data_list$GDT$TRUE_DIST_t[idx,idx ],
      TD_tc=data_list$GDT$TRUE_DIST_tc[idx,idx ],
      mydata = data_list$GDT$mydata[idx, ],
      coords = data_list$GDT$coords[idx, ],
      coords_3col = data_list$GDT$coords_3col[idx, ],
      time = data_list$GDT$time[idx],
      TRUEBETA = data_list$GDT$TRUEBETA[idx, ],
      formula = data_list$GDT$formula
    )
  )
}
