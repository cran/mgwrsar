## configs_list.R
configs_gwr<-list(
  #####
  cfg_gwr_gd_adaptT_H60_NN600_gauss_SEFALSE=list(n=600,
                                                 lambda=0.3,config_beta="default",config_snr=0.9,
                                                 fml="",H=60,Model="GWR",Type="GD",kernels="gauss",
                                                 adaptive=TRUE,NN=600,fixed_vars=NULL,get_s=FALSE,
                                                 SE=FALSE),
  
  cfg_gwr_gd_adaptT_H60_NN600_bisq_SEFALSE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=60,Model="GWR",Type="GD",kernels="bisq",
    adaptive=TRUE,NN=600,fixed_vars=NULL,get_s=FALSE,
    SE=FALSE),   
  
  cfg_mgwr_gd_adaptT_H60_NN600_gauss_SEFALSE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=60,Model="GWR",Type="GD",kernels="gauss",
    adaptive=TRUE,NN=600,fixed_vars="X1",get_s=FALSE,
    SE=FALSE),   
  
  cfg_mgwr_gd_adaptT_H60_NN600_bisq_SEFALSE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=60,Model="GWR",Type="GD",kernels="bisq",
    adaptive=TRUE,NN=600,fixed_vars="X1",get_s=FALSE,
    SE=FALSE),   
  
  cfg_mgwrsar_1_0_kv_gd_adaptT_H60_NN600_gauss_SEFALSE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=60,Model="MGWRSAR_1_0_kv",Type="GD",
    kernels="gauss",adaptive=TRUE,NN=600,fixed_vars=NULL,
    get_s=FALSE,SE=FALSE),   
  
  cfg_mgwrsar_1_0_kv_gd_adaptT_H60_NN600_bisq_SEFALSE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=60,Model="MGWRSAR_1_0_kv",Type="GD",
    kernels="bisq",adaptive=TRUE,NN=600,fixed_vars=NULL,
    get_s=FALSE,SE=FALSE),   
  
  cfg_mgwrsar_0_0_kv_gd_adaptT_H60_NN600_gauss_SEFALSE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=60,Model="MGWRSAR_0_0_kv",Type="GD",
    kernels="gauss",adaptive=TRUE,NN=600,fixed_vars=NULL,
    get_s=FALSE,SE=FALSE),   
  
  cfg_mgwrsar_0_0_kv_gd_adaptT_H60_NN600_bisq_SEFALSE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=60,Model="MGWRSAR_0_0_kv",Type="GD",
    kernels="bisq",adaptive=TRUE,NN=600,fixed_vars=NULL,
    get_s=FALSE,SE=FALSE),   
  
  cfg_gwr_gd_adaptF_H0p3_NN600_gauss_SEFALSE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=0.3,Model="GWR",Type="GD",kernels="gauss",
    adaptive=FALSE,NN=600,fixed_vars=NULL,get_s=FALSE,
    SE=FALSE),   
  
  cfg_gwr_gd_adaptF_H0p3_NN600_bisq_SEFALSE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=0.3,Model="GWR",Type="GD",kernels="bisq",
    adaptive=FALSE,NN=600,fixed_vars=NULL,get_s=FALSE,
    SE=FALSE),   
  
  cfg_mgwr_gd_adaptF_H0p3_NN600_gauss_SEFALSE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=0.3,Model="GWR",Type="GD",kernels="gauss",
    adaptive=FALSE,NN=600,fixed_vars="X1",get_s=FALSE,
    SE=FALSE),   
  
  cfg_mgwr_gd_adaptF_H0p3_NN600_bisq_SEFALSE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=0.3,Model="GWR",Type="GD",kernels="bisq",
    adaptive=FALSE,NN=600,fixed_vars="X1",get_s=FALSE,
    SE=FALSE),   
  
  cfg_mgwrsar_1_0_kv_gd_adaptF_H0p3_NN600_gauss_SEFALSE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=0.3,Model="MGWRSAR_1_0_kv",Type="GD",
    kernels="gauss",adaptive=FALSE,NN=600,fixed_vars=NULL,
    get_s=FALSE,SE=FALSE),   
  
  cfg_mgwrsar_1_0_kv_gd_adaptF_H0p3_NN600_bisq_SEFALSE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=0.3,Model="MGWRSAR_1_0_kv",Type="GD",
    kernels="bisq",adaptive=FALSE,NN=600,fixed_vars=NULL,
    get_s=FALSE,SE=FALSE),   
  
  cfg_mgwrsar_0_0_kv_gd_adaptF_H0p3_NN600_gauss_SEFALSE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=0.3,Model="MGWRSAR_0_0_kv",Type="GD",
    kernels="gauss",adaptive=FALSE,NN=600,fixed_vars=NULL,
    get_s=FALSE,SE=FALSE),   
  
  cfg_mgwrsar_0_0_kv_gd_adaptF_H0p3_NN600_bisq_SEFALSE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=0.3,Model="MGWRSAR_0_0_kv",Type="GD",
    kernels="bisq",adaptive=FALSE,NN=600,fixed_vars=NULL,
    get_s=FALSE,SE=FALSE),   
  
  cfg_gwr_gd_adaptT_H60_NN300_gauss_SEFALSE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=60,Model="GWR",Type="GD",kernels="gauss",
    adaptive=TRUE,NN=300,fixed_vars=NULL,get_s=FALSE,
    SE=FALSE),   
  
  cfg_gwr_gd_adaptT_H60_NN300_bisq_SEFALSE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=60,Model="GWR",Type="GD",kernels="bisq",
    adaptive=TRUE,NN=300,fixed_vars=NULL,get_s=FALSE,
    SE=FALSE),   
  
  cfg_mgwr_gd_adaptT_H60_NN300_gauss_SEFALSE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=60,Model="GWR",Type="GD",kernels="gauss",
    adaptive=TRUE,NN=300,fixed_vars="X1",get_s=FALSE,
    SE=FALSE),   
  
  cfg_mgwr_gd_adaptT_H60_NN300_bisq_SEFALSE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=60,Model="GWR",Type="GD",kernels="bisq",
    adaptive=TRUE,NN=300,fixed_vars="X1",get_s=FALSE,
    SE=FALSE),   
  
  cfg_mgwrsar_1_0_kv_gd_adaptT_H60_NN300_gauss_SEFALSE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=60,Model="MGWRSAR_1_0_kv",Type="GD",
    kernels="gauss",adaptive=TRUE,NN=300,fixed_vars=NULL,
    get_s=FALSE,SE=FALSE),   
  
  cfg_mgwrsar_1_0_kv_gd_adaptT_H60_NN300_bisq_SEFALSE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=60,Model="MGWRSAR_1_0_kv",Type="GD",
    kernels="bisq",adaptive=TRUE,NN=300,fixed_vars=NULL,
    get_s=FALSE,SE=FALSE),   
  
  cfg_mgwrsar_0_0_kv_gd_adaptT_H60_NN300_gauss_SEFALSE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=60,Model="MGWRSAR_0_0_kv",Type="GD",
    kernels="gauss",adaptive=TRUE,NN=300,fixed_vars=NULL,
    get_s=FALSE,SE=FALSE),   
  
  cfg_mgwrsar_0_0_kv_gd_adaptT_H60_NN300_bisq_SEFALSE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=60,Model="MGWRSAR_0_0_kv",Type="GD",
    kernels="bisq",adaptive=TRUE,NN=300,fixed_vars=NULL,
    get_s=FALSE,SE=FALSE),   
  
  cfg_gwr_gd_adaptF_H0p3_NN300_gauss_SEFALSE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=0.3,Model="GWR",Type="GD",kernels="gauss",
    adaptive=FALSE,NN=300,fixed_vars=NULL,get_s=FALSE,
    SE=FALSE),   
  
  cfg_gwr_gd_adaptF_H0p3_NN300_bisq_SEFALSE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=0.3,Model="GWR",Type="GD",kernels="bisq",
    adaptive=FALSE,NN=300,fixed_vars=NULL,get_s=FALSE,
    SE=FALSE),   
  
  cfg_mgwr_gd_adaptF_H0p3_NN300_gauss_SEFALSE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=0.3,Model="GWR",Type="GD",kernels="gauss",
    adaptive=FALSE,NN=300,fixed_vars="X1",get_s=FALSE,
    SE=FALSE),   
  
  cfg_mgwr_gd_adaptF_H0p3_NN300_bisq_SEFALSE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=0.3,Model="GWR",Type="GD",kernels="bisq",
    adaptive=FALSE,NN=300,fixed_vars="X1",get_s=FALSE,
    SE=FALSE),   
  
  cfg_mgwrsar_1_0_kv_gd_adaptF_H0p3_NN300_gauss_SEFALSE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=0.3,Model="MGWRSAR_1_0_kv",Type="GD",
    kernels="gauss",adaptive=FALSE,NN=300,fixed_vars=NULL,
    get_s=FALSE,SE=FALSE),   
  
  cfg_mgwrsar_1_0_kv_gd_adaptF_H0p3_NN300_bisq_SEFALSE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=0.3,Model="MGWRSAR_1_0_kv",Type="GD",
    kernels="bisq",adaptive=FALSE,NN=300,fixed_vars=NULL,
    get_s=FALSE,SE=FALSE),   
  
  cfg_mgwrsar_0_0_kv_gd_adaptF_H0p3_NN300_gauss_SEFALSE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=0.3,Model="MGWRSAR_0_0_kv",Type="GD",
    kernels="gauss",adaptive=FALSE,NN=300,fixed_vars=NULL,
    get_s=FALSE,SE=FALSE),   
  
  cfg_mgwrsar_0_0_kv_gd_adaptF_H0p3_NN300_bisq_SEFALSE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=0.3,Model="MGWRSAR_0_0_kv",Type="GD",
    kernels="bisq",adaptive=FALSE,NN=300,fixed_vars=NULL,
    get_s=FALSE,SE=FALSE),   
  
  cfg_gwr_gd_adaptT_H60_NN600_gauss_SETRUE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=60,Model="GWR",Type="GD",kernels="gauss",
    adaptive=TRUE,NN=600,fixed_vars=NULL,get_s=FALSE,
    SE=TRUE),   
  
  cfg_gwr_gd_adaptT_H60_NN600_bisq_SETRUE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=60,Model="GWR",Type="GD",kernels="bisq",
    adaptive=TRUE,NN=600,fixed_vars=NULL,get_s=FALSE,
    SE=TRUE),   
  
  cfg_mgwr_gd_adaptT_H60_NN600_gauss_SETRUE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=60,Model="GWR",Type="GD",kernels="gauss",
    adaptive=TRUE,NN=600,fixed_vars="X1",get_s=FALSE,
    SE=TRUE),   
  
  cfg_mgwr_gd_adaptT_H60_NN600_bisq_SETRUE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=60,Model="GWR",Type="GD",kernels="bisq",
    adaptive=TRUE,NN=600,fixed_vars="X1",get_s=FALSE,
    SE=TRUE),   
  
  cfg_mgwrsar_1_0_kv_gd_adaptT_H60_NN600_gauss_SETRUE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=60,Model="MGWRSAR_1_0_kv",Type="GD",
    kernels="gauss",adaptive=TRUE,NN=600,fixed_vars=NULL,
    get_s=FALSE,SE=TRUE),   
  
  cfg_mgwrsar_1_0_kv_gd_adaptT_H60_NN600_bisq_SETRUE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=60,Model="MGWRSAR_1_0_kv",Type="GD",
    kernels="bisq",adaptive=TRUE,NN=600,fixed_vars=NULL,
    get_s=FALSE,SE=TRUE),   
  
  cfg_mgwrsar_0_0_kv_gd_adaptT_H60_NN600_gauss_SETRUE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=60,Model="MGWRSAR_0_0_kv",Type="GD",
    kernels="gauss",adaptive=TRUE,NN=600,fixed_vars=NULL,
    get_s=FALSE,SE=TRUE),   
  
  cfg_mgwrsar_0_0_kv_gd_adaptT_H60_NN600_bisq_SETRUE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=60,Model="MGWRSAR_0_0_kv",Type="GD",
    kernels="bisq",adaptive=TRUE,NN=600,fixed_vars=NULL,
    get_s=FALSE,SE=TRUE),   
  
  cfg_gwr_gd_adaptF_H0p3_NN600_gauss_SETRUE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=0.3,Model="GWR",Type="GD",kernels="gauss",
    adaptive=FALSE,NN=600,fixed_vars=NULL,get_s=FALSE,
    SE=TRUE),   
  
  cfg_gwr_gd_adaptF_H0p3_NN600_bisq_SETRUE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=0.3,Model="GWR",Type="GD",kernels="bisq",
    adaptive=FALSE,NN=600,fixed_vars=NULL,get_s=FALSE,
    SE=TRUE),   
  
  cfg_mgwr_gd_adaptF_H0p3_NN600_gauss_SETRUE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=0.3,Model="GWR",Type="GD",kernels="gauss",
    adaptive=FALSE,NN=600,fixed_vars="X1",get_s=FALSE,
    SE=TRUE),   
  
  cfg_mgwr_gd_adaptF_H0p3_NN600_bisq_SETRUE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=0.3,Model="GWR",Type="GD",kernels="bisq",
    adaptive=FALSE,NN=600,fixed_vars="X1",get_s=FALSE,
    SE=TRUE),   
  
  cfg_mgwrsar_1_0_kv_gd_adaptF_H0p3_NN600_gauss_SETRUE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=0.3,Model="MGWRSAR_1_0_kv",Type="GD",
    kernels="gauss",adaptive=FALSE,NN=600,fixed_vars=NULL,
    get_s=FALSE,SE=TRUE),   
  
  cfg_mgwrsar_1_0_kv_gd_adaptF_H0p3_NN600_bisq_SETRUE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=0.3,Model="MGWRSAR_1_0_kv",Type="GD",
    kernels="bisq",adaptive=FALSE,NN=600,fixed_vars=NULL,
    get_s=FALSE,SE=TRUE),   
  
  cfg_mgwrsar_0_0_kv_gd_adaptF_H0p3_NN600_gauss_SETRUE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=0.3,Model="MGWRSAR_0_0_kv",Type="GD",
    kernels="gauss",adaptive=FALSE,NN=600,fixed_vars=NULL,
    get_s=FALSE,SE=TRUE),   
  
  cfg_mgwrsar_0_0_kv_gd_adaptF_H0p3_NN600_bisq_SETRUE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=0.3,Model="MGWRSAR_0_0_kv",Type="GD",
    kernels="bisq",adaptive=FALSE,NN=600,fixed_vars=NULL,
    get_s=FALSE,SE=TRUE),   
  
  cfg_gwr_gd_adaptT_H60_NN300_gauss_SETRUE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=60,Model="GWR",Type="GD",kernels="gauss",
    adaptive=TRUE,NN=300,fixed_vars=NULL,get_s=FALSE,
    SE=TRUE),   
  
  cfg_gwr_gd_adaptT_H60_NN300_bisq_SETRUE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=60,Model="GWR",Type="GD",kernels="bisq",
    adaptive=TRUE,NN=300,fixed_vars=NULL,get_s=FALSE,
    SE=TRUE),   
  
  cfg_mgwr_gd_adaptT_H60_NN300_gauss_SETRUE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=60,Model="GWR",Type="GD",kernels="gauss",
    adaptive=TRUE,NN=300,fixed_vars="X1",get_s=FALSE,
    SE=TRUE),   
  
  cfg_mgwr_gd_adaptT_H60_NN300_bisq_SETRUE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=60,Model="GWR",Type="GD",kernels="bisq",
    adaptive=TRUE,NN=300,fixed_vars="X1",get_s=FALSE,
    SE=TRUE),   
  
  cfg_mgwrsar_1_0_kv_gd_adaptT_H60_NN300_gauss_SETRUE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=60,Model="MGWRSAR_1_0_kv",Type="GD",
    kernels="gauss",adaptive=TRUE,NN=300,fixed_vars=NULL,
    get_s=FALSE,SE=TRUE),   
  
  cfg_mgwrsar_1_0_kv_gd_adaptT_H60_NN300_bisq_SETRUE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=60,Model="MGWRSAR_1_0_kv",Type="GD",
    kernels="bisq",adaptive=TRUE,NN=300,fixed_vars=NULL,
    get_s=FALSE,SE=TRUE),   
  
  cfg_mgwrsar_0_0_kv_gd_adaptT_H60_NN300_gauss_SETRUE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=60,Model="MGWRSAR_0_0_kv",Type="GD",
    kernels="gauss",adaptive=TRUE,NN=300,fixed_vars=NULL,
    get_s=FALSE,SE=TRUE),   
  
  cfg_mgwrsar_0_0_kv_gd_adaptT_H60_NN300_bisq_SETRUE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=60,Model="MGWRSAR_0_0_kv",Type="GD",
    kernels="bisq",adaptive=TRUE,NN=300,fixed_vars=NULL,
    get_s=FALSE,SE=TRUE),   
  
  cfg_gwr_gd_adaptF_H0p3_NN300_gauss_SETRUE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=0.3,Model="GWR",Type="GD",kernels="gauss",
    adaptive=FALSE,NN=300,fixed_vars=NULL,get_s=FALSE,
    SE=TRUE),   
  
  cfg_gwr_gd_adaptF_H0p3_NN300_bisq_SETRUE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=0.3,Model="GWR",Type="GD",kernels="bisq",
    adaptive=FALSE,NN=300,fixed_vars=NULL,get_s=FALSE,
    SE=TRUE),   
  
  cfg_mgwr_gd_adaptF_H0p3_NN300_gauss_SETRUE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=0.3,Model="GWR",Type="GD",kernels="gauss",
    adaptive=FALSE,NN=300,fixed_vars="X1",get_s=FALSE,
    SE=TRUE),   
  
  cfg_mgwr_gd_adaptF_H0p3_NN300_bisq_SETRUE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=0.3,Model="GWR",Type="GD",kernels="bisq",
    adaptive=FALSE,NN=300,fixed_vars="X1",get_s=FALSE,
    SE=TRUE),   
  
  cfg_mgwrsar_1_0_kv_gd_adaptF_H0p3_NN300_gauss_SETRUE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=0.3,Model="MGWRSAR_1_0_kv",Type="GD",
    kernels="gauss",adaptive=FALSE,NN=300,fixed_vars=NULL,
    get_s=FALSE,SE=TRUE),   
  
  cfg_mgwrsar_1_0_kv_gd_adaptF_H0p3_NN300_bisq_SETRUE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=0.3,Model="MGWRSAR_1_0_kv",Type="GD",
    kernels="bisq",adaptive=FALSE,NN=300,fixed_vars=NULL,
    get_s=FALSE,SE=TRUE),   
  
  cfg_mgwrsar_0_0_kv_gd_adaptF_H0p3_NN300_gauss_SETRUE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=0.3,Model="MGWRSAR_0_0_kv",Type="GD",
    kernels="gauss",adaptive=FALSE,NN=300,fixed_vars=NULL,
    get_s=FALSE,SE=TRUE),   
  
  cfg_mgwrsar_0_0_kv_gd_adaptF_H0p3_NN300_bisq_SETRUE=list(
    n=600,lambda=0.3,config_beta="default",config_snr=0.9,
    fml="",H=0.3,Model="MGWRSAR_0_0_kv",Type="GD",
    kernels="bisq",adaptive=FALSE,NN=300,fixed_vars=NULL,
    get_s=FALSE,SE=TRUE)
  #####
)

configs_mgwr <- list(
  #####
  # NN = 600 
  # --- tds_mgwr (get_AIC = FALSE) ---
  cfg_tds_mgwr_gauss_adaptF_getAIC_F_NN600_SEFALSE = list(
    n=600, lambda=NULL, config_beta="default", config_snr=0.9, fml="",
    Model="tds_mgwr", kernels="gauss", adaptive=FALSE, fixed_vars=NULL, NN=600,
    SE=FALSE,Type='GD',
    control_tds=list(nns=20, get_AIC=FALSE, verbose=FALSE, ncore=8)
  ),
  
  cfg_tds_mgwr_gauss_adaptF_getAIC_F_NN600_SETRUE = list(
    n=600, lambda=NULL, config_beta="default", config_snr=0.9, fml="",
    Model="tds_mgwr", kernels="gauss", adaptive=FALSE, fixed_vars=NULL, NN=600,
    SE=TRUE,Type='GD',
    control_tds=list(nns=20, get_AIC=FALSE, verbose=FALSE, ncore=8)
  ),
  
  # --- tds_mgwr (get_AIC = TRUE) ---
  cfg_tds_mgwr_gauss_adaptF_getAIC_T_NN600_SEFALSE = list(
    n=600, lambda=NULL, config_beta="default", config_snr=0.9, fml="",
    Model="tds_mgwr", kernels="gauss", adaptive=FALSE, fixed_vars=NULL, NN=600,
    SE=FALSE,Type='GD',
    control_tds=list(nns=20, get_AIC=TRUE, verbose=FALSE, ncore=8)
  ),
  
  cfg_tds_mgwr_gauss_adaptF_getAIC_T_NN600_SETRUE = list(
    n=600, lambda=NULL, config_beta="default", config_snr=0.9, fml="",
    Model="tds_mgwr", kernels="gauss", adaptive=FALSE, fixed_vars=NULL, NN=600,
    SE=TRUE,Type='GD',
    control_tds=list(nns=20, get_AIC=TRUE, verbose=FALSE, ncore=8)
  ),
  
  # --- atds_mgwr (get_AIC = TRUE) ---
  cfg_atds_mgwr_gauss_adaptF_getAIC_T_NN600_SEFALSE = list(
    n=600, lambda=NULL, config_beta="default", config_snr=0.9, fml="",
    Model="atds_mgwr", kernels="gauss", adaptive=FALSE, fixed_vars=NULL, NN=600,
    SE=FALSE,Type='GD',
    control_tds=list(nns=20, get_AIC=TRUE, verbose=FALSE, ncore=8)
  ),
  
  cfg_atds_mgwr_gauss_adaptF_getAIC_T_NN600_SETRUE = list(
    n=600, lambda=NULL, config_beta="default", config_snr=0.9, fml="",
    Model="atds_mgwr", kernels="gauss", adaptive=FALSE, fixed_vars=NULL, NN=600,
    SE=TRUE,Type='GD',
    control_tds=list(nns=20, get_AIC=TRUE, verbose=FALSE, ncore=8)
  ),
  
  # NN = 300
  
  # --- tds_mgwr (get_AIC = FALSE) ---
  cfg_tds_mgwr_gauss_adaptF_getAIC_F_NN300_SEFALSE = list(
    n=600, lambda=NULL, config_beta="default", config_snr=0.9, fml="",
    Model="tds_mgwr", kernels="gauss", adaptive=FALSE, fixed_vars=NULL, NN=300,
    SE=FALSE,Type='GD',
    control_tds=list(nns=20, get_AIC=FALSE, verbose=FALSE, ncore=8)
  ),
  
  cfg_tds_mgwr_gauss_adaptF_getAIC_F_NN300_SETRUE = list(
    n=600, lambda=NULL, config_beta="default", config_snr=0.9, fml="",
    Model="tds_mgwr", kernels="gauss", adaptive=FALSE, fixed_vars=NULL, NN=300,
    SE=TRUE,Type='GD',
    control_tds=list(nns=20, get_AIC=FALSE, verbose=FALSE, ncore=8)
  ),
  
  # --- tds_mgwr (get_AIC = TRUE) ---
  cfg_tds_mgwr_gauss_adaptF_getAIC_T_NN300_SEFALSE = list(
    n=600, lambda=NULL, config_beta="default", config_snr=0.9, fml="",
    Model="tds_mgwr", kernels="gauss", adaptive=FALSE, fixed_vars=NULL, NN=300,
    SE=FALSE,Type='GD',
    control_tds=list(nns=20, get_AIC=TRUE, verbose=FALSE, ncore=8)
  ),
  
  cfg_tds_mgwr_gauss_adaptF_getAIC_T_NN300_SETRUE = list(
    n=600, lambda=NULL, config_beta="default", config_snr=0.9, fml="",
    Model="tds_mgwr", kernels="gauss", adaptive=FALSE, fixed_vars=NULL, NN=300,
    SE=TRUE,Type='GD',
    control_tds=list(nns=20, get_AIC=TRUE, verbose=FALSE, ncore=8)
  ),
  
  # --- atds_mgwr (get_AIC = TRUE) ---
  cfg_atds_mgwr_gauss_adaptF_getAIC_T_NN300_SEFALSE = list(
    n=600, lambda=NULL, config_beta="default", config_snr=0.9, fml="",
    Model="atds_mgwr", kernels="gauss", adaptive=FALSE, fixed_vars=NULL, NN=300,
    SE=FALSE,Type='GD',
    control_tds=list(nns=20, get_AIC=TRUE, verbose=FALSE, ncore=8)
  ),
  
  cfg_atds_mgwr_gauss_adaptF_getAIC_T_NN300_SETRUE = list(
    n=600, lambda=NULL, config_beta="default", config_snr=0.9, fml="",
    Model="atds_mgwr", kernels="gauss", adaptive=FALSE, fixed_vars=NULL, NN=300,
    SE=TRUE,Type='GD',
    control_tds=list(nns=20, get_AIC=TRUE, verbose=FALSE, ncore=8)
  ),
  #####
  # --- tds_mgtwr (get_AIC = F) ---

  
  cfg_tds_mgtwr_gauss_adaptF_getAIC_F_NN1000_SEFALSE = list(
    n=1000, lambda=NULL, config_beta="spatiotemp", config_snr=0.9, fml="",
    Model="tds_mgtwr", kernels=c("gauss","gauss"), adaptive=c(FALSE,FALSE), fixed_vars=NULL, NN=1000,
    SE=FALSE,Type='GDT',
    control_tds=list(nns=20, get_AIC=FALSE, verbose=FALSE, ncore=8)
  ),
  
  cfg_tds_mgtwr_gauss_adaptF_getAIC_T_NN1000_SETRUE = list(
    n=1000, lambda=NULL, config_beta="spatiotemp", config_snr=0.9, fml="",
    Model="tds_mgtwr", kernels=c("gauss","gauss"), adaptive=c(FALSE,FALSE), fixed_vars=NULL, NN=1000,
    SE=TRUE,Type='GDT',
    control_tds=list(nns=20, get_AIC=TRUE, verbose=FALSE, ncore=8)
  ),
  
  cfg_tds_mgtwr_gauss_adaptF_getAIC_F_NN1000_SEFALSE = list(
    n=1000, lambda=NULL, config_beta="spatiotemp", config_snr=0.9, fml="",
    Model="tds_mgtwr", kernels=c("gauss","gauss_SYM_365"), adaptive=c(FALSE,FALSE), fixed_vars=NULL, NN=1000,
    SE=FALSE,Type='GDT',
    control_tds=list(nns=20, get_AIC=FALSE, verbose=FALSE, ncore=8)
  ),
  cfg_tds_mgtwr_gauss_adaptF_getAIC_T_NN1000_SETRUE = list(
    n=1000, lambda=NULL, config_beta="spatiotemp", config_snr=0.9, fml="",
    Model="tds_mgtwr", kernels=c("gauss","gauss_SYM_365"), adaptive=c(FALSE,FALSE), fixed_vars=NULL, NN=1000,
    SE=TRUE,Type='GDT',
    control_tds=list(nns=20, get_AIC=TRUE, verbose=FALSE, ncore=8)
  )
  #####
)

