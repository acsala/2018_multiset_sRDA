#setwd("~/DATA/data koos/")

library(tester)
library(turner)
library(diagram)
library(shape)
library(amap)
library(elasticnet)
#library(sRDA)
source("/home/Attila/sRDA/sparse_RDA.R")

Generated_data <- generate_data(number_of_ksi = 1,
                                number_of_patients = 250,
                                number_of_Xs_associated_with_ksis = c(100),
                                number_of_not_associated_Xs = 1000,
                                mean_of_the_regression_weights_of_the_associated_Xs = c(0.7),
                                sd_of_the_regression_weights_of_the_associated_Xs = c(0.05),
                                Xnoise_min = -0.3, Xnoise_max = 0.3,
                                number_of_Ys_associated_with_ksis = c(100),
                                number_of_not_associated_Ys = 200,
                                mean_of_the_regression_weights_of_the_associated_Ys = c(0.7),
                                sd_of_the_regression_weights_of_the_associated_Ys = c(0.05),
                                Ynoise_min = -0.3, Ynoise_max = 0.3)

X <- Generated_data$X
Y <- Generated_data$Y

data_info1 <- Generated_data$data_info

Generated_data <- generate_data(number_of_ksi = 1,
                                number_of_patients = 250,
                                number_of_Xs_associated_with_ksis = c(100),
                                number_of_not_associated_Xs = 150,
                                mean_of_the_regression_weights_of_the_associated_Xs = c(0.7),
                                sd_of_the_regression_weights_of_the_associated_Xs = c(0.05),
                                Xnoise_min = -0.3, Xnoise_max = 0.3,
                                number_of_Ys_associated_with_ksis = c(100),
                                number_of_not_associated_Ys = 500,
                                mean_of_the_regression_weights_of_the_associated_Ys = c(0.7),
                                sd_of_the_regression_weights_of_the_associated_Ys = c(0.05),
                                Ynoise_min = -0.3, Ynoise_max = 0.3)

Z <- Generated_data$X
V <- Generated_data$Y

data_info2 <- Generated_data$data_info


data_sets <- cbind(X,Y,Z,V)

dim(X)
dim(Y)
dim(Z)
dim(V)

dim(data_sets)

# ONLY 2 DATASETS####
EXPL_X = c(0,0,0,0)
RESP_Y = c(1,0,0,0)
EXPL_Z = c(0,1,0,0)
RESP_V = c(0,0,1,0)
sat.inner = rbind(EXPL_X, RESP_Y,EXPL_Z,RESP_V)

# blocks of outer model
sat.outer = list(1:dim(X)[2], dim(X)[2]+1:dim(Y)[2],dim(Y)[2]+1:dim(Z)[2],
                 dim(Z)[2]+1:dim(V)[2])
# define vector of reflective modes
sat.mod = c("B", "A","B","A")

library(sPLSPM)

time_data <- system.time(
  s_satpls <- splspm(data_sets, sat.inner, sat.outer, sat.mod, scheme="path",
                  scaled=T, penalization = "ust", nonzero = 50, lambda = 50)
  )

# plot path coefficients
# plot(s_satpls, what="inner")

# plot loadings
# plot(s_satpls, what="loadings")
# 
# # plot outer weights
# plot(s_satpls, what="weights")

save(s_satpls,time_data, file = "sPLM_model_output_1.Rdata")

# ALL 3 DATASETS####
# 
# data_sets <- cbind(X,Y,Z)
# 
# # path matrix
# METHYL = c(0,0,0)
# EXPRES = c(1,0,0)
# CYTO = c(0,1,0)
# sat.inner = rbind(METHYL, EXPRES, CYTO)
# 
# # blocks of outer model
# sat.outer = list(1:dim(X)[2], dim(X)[2]+1:dim(Y)[2], dim(Y)[2]+1:dim(Z)[2])
# # define vector of reflective modes
# sat.mod = c("B","A", "A")
# 
# library(sPLSPM)
# 
# s_satpls = splspm(data_sets, sat.inner, sat.outer, sat.mod, scheme="path",
#                   scaled=T, penalization = "ust", nonzero = 5, lambda = 10,
#                   tol = 0.0001, maxiter = 1000)
# 
# # plot path coefficients
# # plot(s_satpls, what="inner")
# 
# # plot loadings
# # plot(s_satpls, what="loadings")
# # 
# # # plot outer weights
# # plot(s_satpls, what="weights")
# 
# save(s_satpls, file = "sPLM_model_output_3.Rdata")