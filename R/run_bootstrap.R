load("../02_real_data_Marfan/datasetXYZ.RData")

# ALL 3 DATASETS####

data_sets <- cbind(X,Y,Z)

print(dim(data_sets))

# path matrix
METHYL = c(0,0,0)
EXPRES = c(1,0,0)
CYTO = c(0,1,0)
sat.inner = rbind(METHYL, EXPRES, CYTO)

# blocks of outer model
sat.outer = list(1:dim(X)[2], dim(X)[2]+1:dim(Y)[2], dim(Y)[2]+1:dim(Z)[2])
# define vector of reflective modes
sat.mod = c("B","B", "A")

library(sPLSPM)


#********************#Permutation***********************
#
# 
print("Permutation study with reestimation")
start_position <- 1
nr_permutations <- 10

sum_abs_crors.m <- matrix(NA,
                          nrow = max(length(start_position:nr_permutations)),
                          ncol = 4)

nr_counter <- seq(start_position,nr_permutations,1)


for (i in 1:length(nr_counter)){
  
  print(nr_counter[i])
  
  new_row_positions <- sample(nrow(Z), size = 40, replace = T)
  
  size_of_bootstrap_sample <- 37
  rows_for_bootstrap <- new_row_positions[1:size_of_bootstrap_sample]
  
  X_bstrap <- X[rows_for_bootstrap,]
  Y_bstrap <- Y[rows_for_bootstrap,]
  Z_bstrap <- Z[rows_for_bootstrap,]
  
  dim(X_bstrap)
  dim(Y_bstrap)
  dim(Z_bstrap)

  data_sets_bt <- cbind(X_bstrap,Y_bstrap,Z_bstrap)
  sat.outer_bt <- list(1:dim(X_bstrap)[2], 
                      dim(X_bstrap)[2]+1:dim(Y_bstrap)[2], 
                      dim(Y_bstrap)[2]+1:dim(Z_bstrap)[2])
  

  #pdf(paste(i, "resultsFULL.pdf", sep=""))
  print("nr of run")
  print(i)

  
  time_data <- system.time(
    s_satpls2 <- splspm(data_sets_bt, sat.inner, sat.outer_bt, sat.mod, scheme="path",
                        scaled=T, penalization = "ust", 
                        nonzero = c(150), lambda = 1,
                        cross_validate = F)
  )
  
  
  res1.outer <- s_satpls2$outer_model[which(s_satpls2$outer_model[,2] == "METHYL"),]
  res2.outer <- s_satpls2$outer_model[which(s_satpls2$outer_model[,2] == "EXPRES"),]
  res3.outer <- s_satpls2$outer_model[which(s_satpls2$outer_model[,2] == "CYTO"),]
  
  res1.inner <-  s_satpls2$crossloadings[which(s_satpls2$crossloadings[,2]=="METHYL"),]
  res2.inner <-  s_satpls2$crossloadings[which(s_satpls2$crossloadings[,2]=="EXPRES"),]
  res3.inner <-  s_satpls2$crossloadings[which(s_satpls2$crossloadings[,2]=="CYTO"),]
  
  # Sum abs weights of Y with latent variable of X, thus Beta weights for Methyl set
  #sum(abs(res2.inner[,"METHYL"]))
  #print(sum(abs(res2.inner[,"METHYL"])))
  # Sum abs correlation of latent variable X and the Y variables
  sum(abs(cor(s_satpls2$scores[,"METHYL"],Y_bstrap)))
  print(sum(abs(cor(s_satpls2$scores[,"METHYL"],Y_bstrap))))
  sum(abs(res2.inner[,"METHYL"]))
  
  #sum_abs_crors.m[i,1] <- sum(abs(cor(s_satpls2$scores[,"METHYL"],Y)))
  
  # Sum abs weights of Z with latent variable of Y
  #sum(abs(res3.inner[,"EXPRES"]))
  #print(sum(abs(res3.inner[,"EXPRES"])))
  
  
  # Sum abs correlation of latent variable X and the Y variables
  sum(abs(cor(s_satpls2$scores[,"EXPRES"],Z_bstrap)))
  print(sum(abs(cor(s_satpls2$scores[,"EXPRES"],Z_bstrap))))
  sum(abs(res3.inner[,"EXPRES"]))
  
  #sum_abs_crors.m[i,2] <-sum(abs(cor(s_satpls2$scores[,"EXPRES"],Z)))
  
  
  sum_abs_crors.m[i,1] <- sum(abs(cor(s_satpls2$scores[,"METHYL"],Y)))
  sum_abs_crors.m[i,2] <- print(sum(abs(cor(s_satpls2$scores[,"EXPRES"],Z))))
  sum_abs_crors.m[i,3] <- sum(res1.outer[,"weight"]!=0)
  sum_abs_crors.m[i,4] <- sum(res2.outer[,"weight"]!=0)
  
  # plot path coefficients
  # plot(s_satpls, what="inner")
  
  # plot loadings
  # plot(s_satpls, what="loadings")
  # 
  # # plot outer weights
  # plot(s_satpls, what="weights")
  
  #save(s_satpls2,time_data,cv_results, file = "sPLM_results.Rdata")
  
  
  #save objects in RData file
  save(time_data,
       sum_abs_crors.m,
       file = paste("Bootstrap_run_",nr_counter[i], ".RData", sep=""))
  
}

print("We're done have a good1")
