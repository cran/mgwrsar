
configs_pred_MGWR <- list(
cfgm1 = list(n=300, lambda=NULL, config_beta="default", config_snr=0.9,
            Type="GD", Model_est="tds_mgwr", kernels_est=c("gauss"),
            adaptive_est=FALSE, kt=20,
            methods=c("model")),
cfgm2 = list(n=300, lambda=NULL, config_beta="default", config_snr=0.9,
            Type="GD", Model_est="tds_mgwr", kernels_est=c("gauss"),
            adaptive_est=T, kt=20, methods=c("model")),
cfgm3 = list(n=300, lambda=NULL, config_beta="default", config_snr=0.9,
            Type="GD", Model_est="atds_mgwr", kernels_est=c("gauss"),
            adaptive_est=FALSE, kt=20,
            methods=c("model")),
cfgm4 = list(n=300, lambda=NULL, config_beta="spatiotemp", config_snr=0.9,
                        Type="GDT", Model_est="tds_mgtwr", kernels_est=c("gauss","gauss"),
                        adaptive_est=c(FALSE,FALSE), kt=20,
                        methods=c("model"))
)

