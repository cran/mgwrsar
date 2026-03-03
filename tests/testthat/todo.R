# =============================================================================
# TODO
# =============================================================================
if(FALSE){
library(testthat)
library(mgwrsar)

#Sys.setenv(RUN_LONG_TESTS = "1") # pour faire les test MGWR
Sys.setenv(RUN_LONG_TESTS = "2")
Sys.setenv(SAVE_TESTS = "")


devtools::test()
devtools::check(args='--as-cran',remote=T)
}
# search
## good model > badly specified model
#  DGP MGWR MGWR > GWR > GTWR sur Beta et Y_test
#  DGP MGTWR MGTWR > MGWR  et  GTWR > GWR sur Beta et Y_test
#  DGP MLGWRSAR_1_0_kv >  GWR sur Beta, lambda et Y_test

