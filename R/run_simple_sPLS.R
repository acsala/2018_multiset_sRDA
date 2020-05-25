
library(sPLSPM)

Generated_data <- generate_data(number_of_ksi = 1,
                                number_of_patients = 150,
                                number_of_Xs_associated_with_ksis = c(15),
                                number_of_not_associated_Xs = 100,
                                mean_of_the_regression_weights_of_the_associated_Xs = c(0.9),
                                sd_of_the_regression_weights_of_the_associated_Xs = c(0.05),
                                Xnoise_min = -0.2, Xnoise_max = 0.2,
                                number_of_Ys_associated_with_ksis = c(15),
                                number_of_not_associated_Ys = 100,
                                mean_of_the_regression_weights_of_the_associated_Ys = c(0.9),
                                sd_of_the_regression_weights_of_the_associated_Ys = c(0.05),
                                Ynoise_min = -0.2, Ynoise_max = 0.2)

X <- Generated_data$X
Y <- Generated_data$Y

data_info1 <- Generated_data$data_info


Data <- cbind(X,Y)

dim(X)
dim(Y)

dim(Data)

# Only 2 DATASETS####
EXPL_X = c(0,0)
RESP_Y = c(1,0)
path_matrix = rbind(EXPL_X, RESP_Y)

# blocks of outer model
blocks = list(1:dim(X)[2], dim(X)[2]+1:dim(Y)[2])

dim(Data[,blocks[[1]]])

# define vector of reflective modes

modes = c("A","B")



time_data <- system.time(
  s_satpls <- splspm(Data, path_matrix, blocks, modes, scheme="path",
                     scaled=T, penalization = "enet", nonzero = 40, lambda = 1, maxiter = 100)
)

s_satpls$outer_model
s_satpls$model$iter

s_satpls$nonzero
s_satpls$lambda




# ALL DATASETS####

Generated_data <- generate_data(number_of_ksi = 1,
                                number_of_patients = 150,
                                number_of_Xs_associated_with_ksis = c(10),
                                number_of_not_associated_Xs = 150,
                                mean_of_the_regression_weights_of_the_associated_Xs = c(0.9),
                                sd_of_the_regression_weights_of_the_associated_Xs = c(0.05),
                                Xnoise_min = -0.2, Xnoise_max = 0.2,
                                number_of_Ys_associated_with_ksis = c(10),
                                number_of_not_associated_Ys = 100,
                                mean_of_the_regression_weights_of_the_associated_Ys = c(0.9),
                                sd_of_the_regression_weights_of_the_associated_Ys = c(0.05),
                                Ynoise_min = -0.3, Ynoise_max = 0.3)

Z <- Generated_data$X
V <- Generated_data$Y

data_info2 <- Generated_data$data_info

Data <- cbind(X,Y,Z,V)

dim(X)
dim(Y)
dim(Z)
dim(V)

dim(Data)

EXPL_X = c(0,0,0,0)
RESP_Y = c(1,0,0,0)
EXPL_Z = c(0,1,0,0)
RESP_V = c(0,0,1,0)
path_matrix = rbind(EXPL_X, RESP_Y,EXPL_Z,RESP_V)

# blocks of outer model
blocks = list(1:dim(X)[2], dim(X)[2]+1:dim(Y)[2],dim(Y)[2]+1:dim(Z)[2],
              dim(Z)[2]+1:dim(V)[2])

dim(Data[,blocks[[1]]])

# define vector of reflective modes

modes = c("B","B","B","B")



time_data <- system.time(
  s_satpls <- splspm(Data, path_matrix, blocks, modes, scheme="path",
                     scaled=T, penalization = "enet", nonzero = c(20,40),
                     lambda = c(0.1,1), maxiter = 100, cross_validate = T, nr_subsets = 10)
)

s_satpls$outer_model[s_satpls$outer_model[,2]=="EXPL_X",]
s_satpls$outer_model[s_satpls$outer_model[,2]=="RESP_Y",]
s_satpls$outer_model[s_satpls$outer_model[,2]=="EXPL_Z",]

s_satpls$nonzero
s_satpls$lambda
