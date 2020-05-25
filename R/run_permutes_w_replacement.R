load("../datasetXYZ.RData")

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
start_position <- 401
nr_permutations <- 410

sum_abs_crors.m <- matrix(NA,
                          nrow = max(length(start_position:nr_permutations)),
                          ncol = 4)

nr_counter <- seq(start_position,nr_permutations,1)


for (i in 1:length(nr_counter)){
  
  print(nr_counter[i])
  
  Zp <- Z[sample(nrow(Z)),sample(ncol(Z))]
  
  print(rownames(Zp) == rownames(Y))
  
  data_sets <- cbind(X,Y[sample(nrow(Y)),],Zp)
  
  
  #X[sample(nrow(X)),sample(ncol(X))]
  #Y[sample(nrow(Y)),sample(ncol(Y))]
  #Z[sample(nrow(Z)),sample(ncol(Z))]
  
  #pdf(paste(i, "resultsFULL.pdf", sep=""))
  print("nr of run")
  print(i)

  
  time_data <- system.time(
    s_satpls2 <- splspm(data_sets, sat.inner, sat.outer, sat.mod, scheme="path",
                        scaled=T, penalization = "ust", 
                        nonzero = c(10,20,150), lambda = 1,
                        cross_validate = T)
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
  sum(abs(cor(s_satpls2$scores[,"METHYL"],Y)))
  print(sum(abs(cor(s_satpls2$scores[,"METHYL"],Y))))
  
  
  #sum_abs_crors.m[i,1] <- sum(abs(cor(s_satpls2$scores[,"METHYL"],Y)))
  
  # Sum abs weights of Z with latent variable of Y
  #sum(abs(res3.inner[,"EXPRES"]))
  #print(sum(abs(res3.inner[,"EXPRES"])))
  
  
  # Sum abs correlation of latent variable X and the Y variables
  sum(abs(cor(s_satpls2$scores[,"EXPRES"],Z)))
  print(sum(abs(cor(s_satpls2$scores[,"EXPRES"],Z))))
  
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
       file = paste("Permutation_run_",nr_counter[i], ".RData", sep=""))
  
}

print("We're done have a good1")