library(sPLSPM)

multiset_data_generator <- function(N = 500, k = 4, m = 2,
                                    p0=c(500,300,200,10),
                                    p1=c(10,8,8,2)){
  print("alright")

#GENERATE DATA####
# input
## size of the data
# N = 500    # number of individuals
# k = 4      # number of datasets
# m = 2      # number of latent variables (LV's) per dataset

# association parameters between the LV's of the k datasets
# e.g. ksi(dataset4) = b1*ksi(dataset3) + b2*ksi(dataset2) + ...
b=array(0,dim=c(m,k,k))
b[1,1,2]=0.8
b[1,2,3]=0.7
b[1,3,4]=0.4
b[2,1,2]=0.6
b[2,2,3]=0.6
b[2,3,4]=0.2

# # number of irrelevant variables per dataset
# p0=c(500,300,200,10)
# # number of variables associated with the LV's
# p1=c(10,8,8,2)

# specify the regression coefficients of the relevant variables per dataset on the associated LVs
# X = a1*ksi1 + a2*ksi2 + ....
a=array(0,dim=c(max(p1),m,k))
a[1:p1[1],1:2,1]=
  matrix(
    c(.5,0,
      .5,0,
      .5,0,
      .3,0,
      .3,0,
      .2,0.2,
      0,0.3,
      0,0.7,
      0,0.6,
      0,0.6),nrow=p1[1],ncol=m,byrow=TRUE)
a[1:p1[2],1:2,2]=
  matrix(
    c(.5,0,
      .5,0,
      .5,0,
      .2,0.2,
      0.2,0.3,
      0.2,0.7,
      0.2,0.6,
      0.2,0.6),nrow=p1[2],ncol=m,byrow=TRUE)
a[1:p1[3],1:2,3]=
  matrix(
    c(.5,0,
      .5,0,
      .5,0,
      .3,0,
      .3,0,
      .2,0.2,
      0,0.3,
      0,0.6),nrow=p1[3],ncol=m,byrow=TRUE)
a[1:p1[4],1:2,4]=
  matrix(
    c(.5,0,
      0,0.6),nrow=p1[4],ncol=m,byrow=TRUE)


# generate ksi's
ksi=array(NA,dim=c(N,m,k))
for (j in 1:m) {                                     # loop over the number of LV's
  ksi[,j,1] = rnorm(N,0,1)                          # generate values for the LV of the first set of variables
  for (el in 2:k) {                                 # loop over the other sets of variables
    meanx=0
    sumb2=0
    for (ell in (el-1):1) {                        # calculate per person the mean and sd of the LVs
      meanx=meanx+b[j,ell,el]*ksi[,j,ell]
      sumb2=sumb2+b[j,ell,el]
    }
    sdx=max(0.0001,sqrt(1-sumb2))
    ksi[,j,el] = rnorm(N,meanx,sdx)                # sample LV values from the normal distribution
  }
  ksi[,j,]=scale(ksi[,j,])
}

# generate manifest data
X=array(NA,dim=c(N,max(p0+p1),k))
for (el in 1:k) {
  for (j in 1:p0[el]) {                            # sample values for the irrelevant manifest variables
    X[,j,el]=rnorm(N,0,1)
  }
  meanx=ksi[,1:m,el]%*%t(a[1:p1[el],1:m,el])       # calculate per person means and sds of the relevant manifest variables
  suma2=apply(a[,,el],1,sum)
  for (j in (p0[el]+1):(p0[el]+p1[el])) {
    sdx=min(0.0001,sqrt(1-suma2[(j-p0[el])]))
    X[,j,el]=rnorm(N,meanx[,(j-p0[el])],sdx)      # sample values for the manifest variables
  }
}

X1=X[,,1]
X2=X[,,2]
X3=X[,,3]
X4=X[,,4]
if (length(which(is.na(apply(X1,2,sd,na.rm=TRUE))))>0) {X1=X1[,-which(is.na(apply(X1,2,sd,na.rm=TRUE)))]}
if (length(which(is.na(apply(X2,2,sd,na.rm=TRUE))))>0) {X2=X2[,-which(is.na(apply(X2,2,sd,na.rm=TRUE)))]}
if (length(which(is.na(apply(X3,2,sd,na.rm=TRUE))))>0) {X3=X3[,-which(is.na(apply(X3,2,sd,na.rm=TRUE)))]}
if (length(which(is.na(apply(X4,2,sd,na.rm=TRUE))))>0) {X4=X4[,-which(is.na(apply(X4,2,sd,na.rm=TRUE)))]}

c(dim(X1),dim(X2),dim(X3),dim(X4))

X1=scale(X1)
X2=scale(X2)
X3=scale(X3)
X4=scale(X4)

result <- list(X1,
               X2,
               X3,
               X4
               )

names(result) <- c("X1",
                   "X2",
                   "X3",
                   "X4"
)

result
}

#end of function **********

## data generation
N = 50    # number of individuals
k = 4      # number of datasets
m = 2      # number of latent variables (LV's) per dataset
# number of irrelevant variables per dataset
p0=c(200,150,120,110)
# number of variables associated with the LV's
p1=c(10,8,8,2)


#repeat simulation ####

nr_of_simulations <- 1

for (i in nr_of_simulations){
  
  #pdf(paste(i, "resultsFULL.pdf", sep=""))
  print("nr of simulation")
  print(i)    
  
  # for replication
  nr_seed = runif(1, 10, 10^8)
  set.seed(nr_seed)
  
  multiblockdata <- multiset_data_generator(N, k,m,
                                            p0,
                                            p1)
  X1 <- multiblockdata$X1
  X2 <- multiblockdata$X2
  X3 <- multiblockdata$X3
  X4 <- multiblockdata$X4
  
  Data <- cbind(X1,X2,X3,X4)
  EXPL_X = c(0,0,1,0)
  RESP_Y = c(1,0,0,0)
  EXPL_Z = c(1,0,0,0)
  RESP_V = c(1,0,1,0)
  path_matrix = rbind(EXPL_X, RESP_Y,EXPL_Z,RESP_V)
  
  # blocks of outer model
  blocks = list(1:dim(X1)[2], (dim(X1)[2]+1):(dim(X1)[2]+dim(X2)[2]),
                (dim(X1)[2]+dim(X2)[2]+1):(dim(X1)[2]+dim(X2)[2]+dim(X3)[2]),
                (dim(X1)[2]+dim(X2)[2]+dim(X3)[2]+1):(dim(X1)[2]+dim(X2)[2]+dim(X3)[2]+dim(X4)[2]))
  
  modes <- replace( rep("A",dim(path_matrix)[1]),
          c(apply(path_matrix, 2, sum)!=0), c("B"))
  
  
  #modes = c("B","A","B","A")
  
  time_data <- system.time(
    s_satpls <- splspm(Data, path_matrix, blocks, modes, scheme="path",
                       scaled=T, penalization = "ust", nonzero = c(5,10), 
                       lambda = 1, maxiter = 100, cross_validate = T)
  )
  
  s_satpls$cv_results
  
  #s_satpls$outer_model
  print("s_satpls$model$iter")
  s_satpls$model$iter
  
  print("s_satpls$nonzero")
  s_satpls$nonzero
  print("s_satpls$lambda")
  s_satpls$lambda
  
  # checking non zeros
  nzero_X_positiong <- which(abs(s_satpls$outer_model[s_satpls$outer_model[,2]=="EXPL_X",3])>0)
  nzero_Z_positiong <- which(abs(s_satpls$outer_model[s_satpls$outer_model[,2]=="EXPL_Z",3])>0)
  nzero_Y_positiong <- which(abs(s_satpls$outer_model[s_satpls$outer_model[,2]=="RESP_Y",3])>0)
  
  # mode A response will not have nonzeros
  #nzero_V_positiong <- which(abs(s_satpls$outer_model[s_satpls$outer_model[,2]=="RESP_V",3])>0)
  
  nzero_X_positiong
  nzero_Y_positiong
  nzero_Z_positiong
  c(dim(X1),dim(X2),dim(X3),dim(X4))
  
  print("nr of associated variables per dataset")
  p1
  
  # s_satpls$outer_model[which(abs(s_satpls$outer_model[s_satpls$outer_model[,2]=="EXPL_X",3])>0),]
  # s_satpls$outer_model[510+which(abs(s_satpls$outer_model[s_satpls$outer_model[,2]=="RESP_Y",3])>0),]
  # s_satpls$outer_model[818+which(abs(s_satpls$outer_model[s_satpls$outer_model[,2]=="EXPL_Z",3])>0),]
  # s_satpls$outer_model[1026+which(abs(s_satpls$outer_model[s_satpls$outer_model[,2]=="RESP_V",3])>0),]
  
  
  #save objects in RData file
  save(nr_seed,
       s_satpls,
       N,k,m,
       p0,
       p1,
       file = paste(i, "_nrResults.RData", sep=""))
  
  #close pdf
  #dev.off()
  
  
}


res1.outer <- s_satpls$outer_model[which(s_satpls$outer_model[,2] == "EXPL_X"),]
res2.outer <- s_satpls$outer_model[which(s_satpls$outer_model[,2] == "RESP_Y"),]
res3.outer <- s_satpls$outer_model[which(s_satpls$outer_model[,2] == "EXPL_Z"),]
res4.outer <- s_satpls$outer_model[which(s_satpls$outer_model[,2] == "RESP_V"),]

res1.inner <-  s_satpls$crossloadings[which(s_satpls$crossloadings[,2]=="EXPL_X"),]
res2.inner <-  s_satpls$crossloadings[which(s_satpls$crossloadings[,2]=="RESP_Y"),]
res3.inner <-  s_satpls$crossloadings[which(s_satpls$crossloadings[,2]=="EXPL_Z"),]
res4.inner <-  s_satpls$crossloadings[which(s_satpls$crossloadings[,2]=="RESP_V"),]

# Sum abs weights of Y with latent variable of X
sum(abs(res2.inner[,"EXPL_X"]))
# Sum abs correlation of latent variable X and the Y variables
sum(abs(cor(s_satpls$scores[,"EXPL_X"],X2)))


# Sum abs weights of Z with latent variable of Y
sum(abs(res3.inner[,"RESP_Y"]))
# Sum abs correlation of latent variable X and the Y variables
sum(abs(cor(s_satpls$scores[,"RESP_Y"],X3)))

# Sum abs weights of Z with latent variable of Y
sum(abs(res4.inner[,"RESP_V"]))
# Sum abs correlation of latent variable X and the Y variables
sum(abs(cor(s_satpls$scores[,"RESP_V"],X4)))
