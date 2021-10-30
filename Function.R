#############################################################
#####                                                   #####
#####                                                   #####
#####                    Data Generation                #####
#####                                                   #####
#####                                                   #####
#############################################################

################ Simulation Data Functions ###################
### Generate normal distribution data for in all Cases.    ###
### There are G=2,3,4 three cases.                         ###

# Packages
library(MASS)
library(dplyr)
library(sparcl)
library(cluster)
library(factoextra)
library(CEoptim)
library(mclust)
library(mvtnorm)
library(JGL)
library(glasso)
library(Matrix)

################################################################################

# G=2
PECM_simu_gene_2 <- function(K,K1){
  ###chain network
  A <- matrix(0,nrow = K,ncol = K)
  t <- 1
  while(t==1||!all(ev >= -1e-06 * abs(ev[1L]))){
    for (i in 61:79) {
      A[i,(i+1)] <- exp(-runif(1,0.5,1))
      A[(i+1),i] <- A[i,(i+1)]
    }
    right_omega <- A
    diag(A) <- 1.1
    A <- try(solve(A),silent = TRUE)
    while('try-error' %in% class(A)){
      A <- matrix(0,nrow = K,ncol = K)
      for (i in 61:79) {
        A[i,(i+1)] <- exp(-runif(1,0.5,1))
        A[(i+1),i] <- A[i,(i+1)]
      }
      right_omega <- A
      diag(A) <- 1.1
      A <- try(solve(A),silent = TRUE)        
    }
    for (i in 61:80) {
      A[i,] <- A[i,]/A[i,i]
    }
    ev <- eigen(A,symmetric = T)$values
    t <- t + 1
  }
  D <- matrix(0,nrow = K,ncol = K)
  diag(D) <- 1
  
  x_1 <- mvrnorm(n=100, mu = c(rep(0.7,K1),rep(0,K-K1)), A)
  x_2 <- mvrnorm(n=100, mu = c(rep(0,K1),rep(0,K-K1)), D)
  
  label <- rep(1, nrow(x_1))
  yang1 <- data.frame(x_1,label)
  label <- rep(2, nrow(x_2))
  yang2 <- data.frame(x_2,label)
  
  yang <- bind_rows(yang1, yang2)
  X <- yang[,1:K]
  return(list(X=X,
              label=yang$label,
              right_omega=right_omega))
}

# G=3
PECM_simu_gene_3 <- function(K,K1){
  ###chain network
  A <- matrix(0,nrow = K,ncol = K)
  t <- 1
  while(t==1||!all(ev >= -1e-06 * abs(ev[1L]))){
    for (i in 61:79) {
      A[i,(i+1)] <- exp(-runif(1,0.5,1))
      A[(i+1),i] <- A[i,(i+1)]
    }
    right_omega <- A
    diag(A) <- 1.1
    A <- try(solve(A),silent = TRUE)
    while('try-error' %in% class(A)){
      A <- matrix(0,nrow = K,ncol = K)
      for (i in 61:79) {
        A[i,(i+1)] <- exp(-runif(1,0.5,1))
        A[(i+1),i] <- A[i,(i+1)]
      }
      right_omega <- A
      diag(A) <- 1.1
      A <- try(solve(A),silent = TRUE)        
    }
    for (i in 61:80) {
      A[i,] <- A[i,]/A[i,i]
    }
    ev <- eigen(A,symmetric = T)$values
    t <- t + 1
  }
  D <- matrix(0,nrow = K,ncol = K)
  diag(D) <- 1
  
  x_1 <- mvrnorm(n=100, mu = c(rep(0.7,K1),rep(0,K-K1)), A)
  x_2 <- mvrnorm(n=100, mu = c(rep(0,K1),rep(0.7,K1),rep(0,K-2*K1)), D)
  x_3 <- mvrnorm(n=100, mu = c(rep(0,2*K1),rep(0.7,K1),rep(0,K-3*K1)), D)
  label <- rep(1, nrow(x_1))
  yang1 <- data.frame(x_1,label)
  label <- rep(2, nrow(x_2))
  yang2 <- data.frame(x_2,label)
  label <- rep(3, nrow(x_3))
  yang3 <- data.frame(x_3,label)
  yang <- bind_rows(yang1, yang2, yang3)
  X <- yang[,1:K]
  return(list(X=X,
              label=yang$label,
              right_omega=right_omega))
}


# G=4
PECM_simu_gene_4 <- function(K,K1){
  ###chain network
  A <- matrix(0,nrow = K,ncol = K)
  t <- 1
  while(t==1||!all(ev >= -1e-06 * abs(ev[1L]))){
    for (i in 61:79) {
      A[i,(i+1)] <- exp(-runif(1,0.5,1))
      A[(i+1),i] <- A[i,(i+1)]
    }
    right_omega <- A
    diag(A) <- 1.1
    A <- try(solve(A),silent = TRUE)
    while('try-error' %in% class(A)){
      A <- matrix(0,nrow = K,ncol = K)
      for (i in 61:79) {
        A[i,(i+1)] <- exp(-runif(1,0.5,1))
        A[(i+1),i] <- A[i,(i+1)]
      }
      right_omega <- A
      diag(A) <- 1.1
      A <- try(solve(A),silent = TRUE)        
    }
    for (i in 61:80) {
      A[i,] <- A[i,]/A[i,i]
    }
    ev <- eigen(A,symmetric = T)$values
    t <- t + 1
  }
  D <- matrix(0,nrow = K,ncol = K)
  diag(D) <- 1
  
  x_1 <- mvrnorm(n=100, mu = c(rep(0.7,K1),rep(0,K-K1)), A)
  x_2 <- mvrnorm(n=100, mu = c(rep(0,K1),rep(0.7,K1),rep(0,K-2*K1)), D)
  x_3 <- mvrnorm(n=100, mu = c(rep(0,2*K1),rep(0.7,K1),rep(0,K-3*K1)), D)
  x_4 <- mvrnorm(n=100, mu = c(rep(0,3*K1),rep(0.7,K1),rep(0,K-4*K1)), D)
  label <- rep(1, nrow(x_1))
  yang1 <- data.frame(x_1,label)
  label <- rep(2, nrow(x_2))
  yang2 <- data.frame(x_2,label)
  label <- rep(3, nrow(x_3))
  yang3 <- data.frame(x_3,label)
  label <- rep(4, nrow(x_4))
  yang4 <- data.frame(x_4,label)
  yang <- bind_rows(yang1, yang2, yang3, yang4)
  X <- yang[,1:K]
  return(list(X=X,
              label=yang$label,
              right_omega=right_omega))
}


################################################################################









#############################################################
#####                                                   #####
#####                                                   #####
#####                    Main Functions                 #####
#####                                                   #####
#####                                                   #####
#############################################################

################## Preparatory functions ####################

############# Prior information transformation ##############

clu_mat <- function(Y_prior){
  clu_matrix <- matrix(0, nrow = length(Y_prior), ncol = max(Y_prior))
  for (i in 1:length(Y_prior)) {
    clu_matrix[i, Y_prior[i]] <- 1
  }
  return(clu_matrix)
}

#############################################################

################## Calculation function #####################
myfun_g2_dong <- function(x){colSums(Sigma[[i]][c(marvel),c(marvel)]*as.numeric(x))}

is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
is.wholenumber(1)

#############################################################


############### Calculate CER of results ####################
cal_CER <- function(predict, real){
  P <- matrix(0,nrow = length(predict), ncol = length(predict))
  for (i in 1:(length(predict)-1)) {
    for (j in (i+1):length(predict)) {
      if(predict[i]==predict[j]) P[i,j] <- 1
    }
  }
  Q <- matrix(0,nrow = length(real), ncol = length(real))
  for (i in 1:(length(real)-1)) {
    for (j in (i+1):length(real)) {
      if(real[i]==real[j]) Q[i,j] <- 1
    }
  }
  CER <- sum(abs(P-Q))*2/(length(predict)*(length(predict)-1))
  return(CER)
}

#############################################################


########### Calculate TPR and FPR of results ################

cal_TPR_edge <- function(E,G){
  TPR <- NULL
  FPR <- NULL
  PME <- NULL
  for (i in 1:G) {
    if('ddiMatrix'%in%class(E[[i]])){
      TPR <- c(TPR,0)
      FPR <- c(FPR,0)
    } else{
      edge_m <- as.matrix(summary(E[[i]]))
      edge_selected <- edge_m[,1:2][which(edge_m[,1]!=edge_m[,2]),]
      right_edge <- matrix(c(61:79,62:80,62:80,61:79),ncol = 2)
      total_edge <- rbind(right_edge,edge_selected)
      TPR <- c(TPR,nrow(total_edge[duplicated(total_edge),])/38)
      FPR <- c(FPR,(nrow(edge_selected)-nrow(total_edge[duplicated(total_edge),]))/999000)
    }
  }
  for (i in 1:G) {
    if(i==which.max(TPR)){
      PME <- c(PME,sqrt(sum((E[[i]]-simul_data$right_omega)^2)))
    } else{
      PME <- c(PME,sqrt(sum((E[[i]]-diag(1,nrow = 1000,ncol = 1000))^2)))
    }
  }
  which.max(TPR)
  
  return(list(TPR=max(TPR),
              FPR=1-sum(FPR),
              PME=sum(PME)/G))
}






####################### Main functions #########################


###### The Method that fully trust prior information. ##########
###### G is the number of clusters to be divided. ##############
###### prior_mu is the prior information of mean. ##############
###### prior_sigma is the prior information of interactions. ###
###### lambda1 and lambda2 are penalties. ######################
###### B is the upper bound of iteration. ######################

APECM <- function(data, G, prior_mu, prior_sigma,
                  lambda1, lambda2,
                  B = 40){

##### Input: #####
#   data: Data with continuous variables.
#   G: The number of clusters.
#   prior_mu: The prior information of the variables.
#   prior_sigma: The prior information of interactions.
#   lambda1: The penalty parameter controlling variables.
#   lambda2: The penalty parameter controlling interactions.
#   B: The upper bound of iteration.
##### Output: #####
#   variable_selected: The significant variables selected by the model. 
#   net_selected: The significant interactions selected by the model.
#   label_list: Clustering label of each iteration.
#   tau_list: Probability of each iteration.
#   miu_list: Mean values of each iteration.
#   label: The clustering result.
#   Omega: The interaction estimations.
#   miu: The mean value for each cluster.

  lambda1 <- lambda1
  G <- G
  B <- B
  
  X <- data
  init_X <- X
  N <- nrow(X)
  K <- ncol(X)
  
  init_cluster <- kmeans(X,G)$cluster
  pai <- c(rep(1/G, G))
  label <- init_cluster
  
  
  Y <- list()
  for (i in 1:G) {
    Y[[i]] <- X[which(label==i),]
  }
  
  miu <- NULL
  for (i in 1:G) {
    miu <- cbind(miu, colMeans(Y[[i]]))
  }
  Sigma <- list()
  for (i in 1:G) {
    Sigma[[i]] <- diag(rep(1,K))
  }
  Omega <- list()
  for (i in 1:G) {
    Omega[[i]] <- diag(rep(1,K))
  }
  left <- as.numeric(which((rowSums(miu==0)<G)))
  b <- 1
  
  g2 <- matrix(0, nrow = K, ncol = G)
  
  lambda2_m <- matrix(lambda2, nrow = K, ncol = K)
  lambda2_m[prior_sigma] <- 0
  
  
  label_list <- list()
  tau_list <- list()
  miu_list <- list()
  connected <- NULL
  
  miu_old <- miu - 1
  
  
  while (b<=B&(sqrt(sum((miu-miu_old)^2))/sqrt(sum((miu)^2)))) {
    
    miu_old <- miu
    Omega_old <- Omega
    ###E step
    left <- as.numeric(which((rowSums(miu==0)<G)))
    if(b>1){
      for (i in 1:G) {
        connected <- c(connected,
                       c(unique(which(Omega[[i]]!=0,arr.ind = T)[,1][duplicated(which(Omega[[i]]!=0,arr.ind = T)[,1])]),
                         unique(which(Omega[[i]]!=0,arr.ind = T)[,2][duplicated(which(Omega[[i]]!=0,arr.ind = T)[,2])])))
      }
    }
    connected <- unique(connected)
    marvel <- union(left, connected)
    X <- init_X[,c(marvel)]
    N <- nrow(X)
    K <- ncol(X)
    f <- matrix(0, nrow = N, ncol = G)
    for (i in 1:G) {
      f[,i] <- as.numeric(exp(dmvnorm(X, mean = miu[c(marvel),i], sigma = Sigma[[i]][c(marvel),c(marvel)], log = T)+length(marvel)*0.8))
    }
    #f <- matrix(0, nrow = N, ncol = G)
    #for (i in 1:G) {
    #  f[,i] <- as.numeric(dmvnorm(X, mean = miu[,i], sigma = Sigma[[i]], log = F))
    #}
    tau <- t(pai*t(f))/rowSums(t(pai*t(f)))
    tau[is.nan(tau)] <- 0.5
    tau_list[[b]] <- tau
    pai <- colMeans(tau)
    ###M step
    eta_p <- tau
    #Sigma <- shu$theta
    
    for (i in 1:G) {
      myfun_g2_dong <- function(x){colSums(Omega[[i]][c(marvel),c(marvel)]*as.numeric(x))}
      g2[c(marvel),i] <- rowSums(t(t(apply(t(t(X)-miu[c(marvel),i]), 1, myfun_g2_dong)+miu[c(marvel),i]*diag(Omega[[i]][c(marvel),c(marvel)]))*eta_p[,i]))
    }
    for (i in 1:G) {
      myfun_g2_dong <- function(x){colSums(Omega[[i]][c(marvel),c(marvel)]*as.numeric(x))}
      miu[c(marvel),i] <- colSums(t(apply(X,1,myfun_g2_dong))*eta_p[,i])/(sum(eta_p[,i])*diag(Omega[[i]][c(marvel),c(marvel)]))-
        colSums(miu[c(marvel),i]*Omega[[i]][c(marvel),c(marvel)])/diag(Omega[[i]][c(marvel),c(marvel)])+miu[c(marvel),i]-
        sign(miu[c(marvel),i])/diag(Omega[[i]][c(marvel),c(marvel)])*lambda1/sum(eta_p[,i])
      miu[setdiff(as.numeric(which(abs(g2[,i])<(lambda1/abs(miu[,i])))), c(prior_mu)),i] <- 0
    }
    miu_list[[b]] <- miu
    #if(b==1) label <- apply(eta_p, 1, which.max)
    label <- apply(eta_p, 1, which.max)
    label_list[[b]] <- label
    for (i in 1:G) {
      Y[[i]] <- init_X[which(label==i),]
      gla <- glasso(var(Y[[i]]),rho = lambda2_m,penalize.diagonal = FALSE,approx = FALSE)
      Sigma[[i]] <- gla$w
      Omega[[i]] <- Matrix(gla$wi,sparse = TRUE)
    }
    #for (i in 1:G) {
    #  Sigma[[i]] <- solve(shu$theta[[i]])
    #}
    print(b)
    b <- b + 1
  }
  
  result <- list(variable_selected = left, 
                 net_selected = connected,
                 label_list = label_list,
                 tau_list = tau_list,
                 miu_list = miu_list,
                 label = label,
                 Omega = Omega,
                 miu = miu)
  return(result)
  
}

################################################################


#################### The proposed method. ######################
###### G is the number of clusters to be divided. ##############
###### prior_information is the result of APECM. ###############
###### lambda1 and lambda2 are penalties. ######################
###### eta is the parameter controlling the weight of prior. ###
###### B is the upper bound of iteration. ######################

PECM <- function(data, G, prior_information,
                 lambda1, lambda2, eta,
                 B = 40){
  
##### Input: #####
#   data: Data with continuous variables.
#   G: The number of clusters.
#   prior_information: The labels obtained by prior information.
#   lambda1: The penalty parameter controlling variables.
#   lambda2: The penalty parameter controlling interactions.
#   eta: The tuning parameter controlling the weight of the prior information.
#   B: The upper bound of iteration.
##### Output: #####
#   variable_selected: The significant variables selected by the model. 
#   net_selected: The significant interactions selected by the model.
#   label_list: Clustering label of each iteration.
#   tau_list: Probability of each iteration.
#   miu_list: Mean values of each iteration.
#   label: The clustering result.
#   Omega: The interaction estimations.
#   miu: The mean value for each cluster.
#   mBIC: The mBIC critrion of each iteration.
  
  ###parameters
  lambda1 <- lambda1
  lambda2 <- lambda2
  G <- G
  eta <- eta
  B <- B
  
  ###initialize
  X <- data
  init_X <- X
  N <- nrow(X)
  K <- ncol(X)
  
  ###initialize π
  init_cluster <- kmeans(X,G)$cluster
  pai <- c(rep(1/G, G))
  label <- init_cluster
  
  ###prior information
  z <- clu_mat(prior_information)
  
  ###initialize data
  Y <- list()
  for (i in 1:G) {
    Y[[i]] <- X[which(label==i),]
  }
  
  ###initialize mean and covariance
  miu <- NULL
  for (i in 1:G) {
    miu <- cbind(miu, colMeans(Y[[i]]))
  }
  Sigma <- list()
  for (i in 1:G) {
    Sigma[[i]] <- diag(rep(1,K))
  }
  Omega <- list()
  for (i in 1:G) {
    Omega[[i]] <- diag(rep(1,K))
  }
  g2 <- matrix(0, nrow = K, ncol = G)
  
  left <- as.numeric(which((rowSums(miu==0)<G)))
  b <- 1
  
  ###initialize record
  label_list <- list()
  tau_list <- list()
  miu_list <- list()
  connected <- NULL
  
  miu_old <- miu - 1
  
  while (b<=B&(sqrt(sum((miu-miu_old)^2))/sqrt(sum((miu)^2)))) {
    
    miu_old <- miu
    Omega_old <- Omega
    sqrt(sum((miu-miu_old)^2))
    ###E step
    left <- as.numeric(which((rowSums(miu==0)<G)))
    if(b>1){
      for (i in 1:G) {
        connected <- c(connected,
                       c(unique(which(Omega[[i]]!=0,arr.ind = T)[,1][duplicated(which(Omega[[i]]!=0,arr.ind = T)[,1])]),
                         unique(which(Omega[[i]]!=0,arr.ind = T)[,2][duplicated(which(Omega[[i]]!=0,arr.ind = T)[,2])])))
      }
    }
    connected <- unique(connected)
    marvel <- union(left, connected)
    X <- init_X[,c(marvel)]
    N <- nrow(X)
    K <- ncol(X)
    f <- matrix(0, nrow = N, ncol = G)
    for (i in 1:G) {
      f[,i] <- as.numeric(exp(dmvnorm(X, mean = miu[c(marvel),i], sigma = Sigma[[i]][c(marvel),c(marvel)], log = T)+length(marvel)*0.8))
    }
    #f <- matrix(0, nrow = N, ncol = G)
    #for (i in 1:G) {
    #  f[,i] <- as.numeric(dmvnorm(X, mean = miu[,i], sigma = Sigma[[i]], log = F))
    #}
    tau <- t(pai*t(f))/rowSums(t(pai*t(f)))
    tau[is.nan(tau)] <- 0.5
    tau_list[[b]] <- tau
    pai <- colMeans(tau)
    ###M step
    eta_p <- (1-eta)*tau+eta*z
    #Sigma <- shu$theta
    
    for (i in 1:G) {
      myfun_g2_dong <- function(x){colSums(Omega[[i]][c(marvel),c(marvel)]*as.numeric(x))}
      g2[c(marvel),i] <- rowSums(t(t(apply(t(t(X)-miu[c(marvel),i]), 1, myfun_g2_dong)+miu[c(marvel),i]*diag(Omega[[i]][c(marvel),c(marvel)]))*eta_p[,i]))
    }
    for (i in 1:G) {
      myfun_g2_dong <- function(x){colSums(Omega[[i]][c(marvel),c(marvel)]*as.numeric(x))}
      miu[c(marvel),i] <- colSums(t(apply(X,1,myfun_g2_dong))*eta_p[,i])/(sum(eta_p[,i])*diag(Omega[[i]][c(marvel),c(marvel)]))-
        colSums(miu[c(marvel),i]*Omega[[i]][c(marvel),c(marvel)])/diag(Omega[[i]][c(marvel),c(marvel)])+miu[c(marvel),i]-
        sign(miu[c(marvel),i])/diag(Omega[[i]][c(marvel),c(marvel)])*lambda1/sum(eta_p[,i])
      miu[which(abs(g2[,i])<(lambda1/abs(miu[,i]))),i] <- 0
    }
    miu_list[[b]] <- miu
    #if(b==1) label <- apply(eta_p, 1, which.max)
    if(is.wholenumber(b/4)){
      label <- kmeans(init_X[,c(left)],G)$cluster
      label_list[[b]] <- label
    } else{
      label <- apply(eta_p, 1, which.max)
      label_list[[b]] <- label
    }
    
    for (i in 1:G) {
      Y[[i]] <- init_X[which(label==i),]
      gla <- glasso(var(Y[[i]]),rho = lambda2,penalize.diagonal = FALSE,approx = FALSE)
      Sigma[[i]] <- gla$w
      Omega[[i]] <- Matrix(gla$wi,sparse = TRUE)
    }
    
    #for (i in 1:G) {
    #  Sigma[[i]] <- solve(shu$theta[[i]])
    #}
    print(b)
    b <- b + 1
  }
  
  #mBIC
  X <- data
  init_X <- X
  N <- nrow(X)
  K <- ncol(X)
  f <- matrix(0, nrow = N, ncol = G)
  for (i in 1:G) {
    f[,i] <- as.numeric(exp(dmvnorm(X, mean = miu[,i], sigma = Sigma[[i]], log = T)+K*0.8))
  }
  q <- K - sum((rowSums(miu==0)<G))
  de <- G + K + G*K - 1 - q
  eta <- abs(eta-1e-4)
  mBIC <- (-2)*(sum(log(colSums(pai*t(f))))
                +sum(log(rowSums(z*f)*eta))
  )
  +log(N)*de
  
  result <- list(variable_selected = left, 
                 net_selected = connected,
                 label_list = label_list,
                 tau_list = tau_list,
                 miu_list = miu_list,
                 label = label,
                 Omega = Omega,
                 miu = miu,
                 mBIC = mBIC)
  return(result)
}



################################################################

############# The method without interactions. #################
###### X is the data and G is the number of clusters. ##########
###### Assume that there is no interactions. ###################


###benchmark
gouji <- function(X, G, prior_information, eta, lambda, B=5){

  N <- nrow(X)
  K <- ncol(X)
  

  init_cluster <- kmeans(X,G)$cluster
  category <- c(unique(init_cluster))
  init_X <- cbind.data.frame(X, category)
  init_weight <- data.frame(table(init_cluster)/length(init_cluster))
  pai <- init_weight$Freq
  #pai <- c(rep(0.5, G))
  

  z <- clu_mat(prior_information)
  

  miu <- NULL
  for (i in 1:G) {
    miu <- cbind(miu, runif(K, -0.5, 0.5))
  }
  sigma <- data.frame(diag=diag(cov(X)))
  
  b <- 1
  
  while (b<=B) {
    ###E step
    left <- which((rowSums(miu==0)!=2))
    #left <- sample(left,20,replace = F)
    f <- matrix(0, nrow = N, ncol = G)
    for (i in 1:G) {
      f[,i] <- as.numeric(exp(dmvnorm(X[,c(left)], mean = miu[c(left),i], sigma = diag(sigma$diag)[c(left), c(left)], log = T)+length(left)*0.8))
    }
    
    tau <- t(pai*t(f))/rowSums(t(pai*t(f)))
    #tau[is.nan(tau)] <- 0.5
    
    pai <- colMeans(tau)
    #pai <- colSums(f/rowSums(f))/N
    
    ###M step
    eta_p <- (1-eta)*tau+eta*z
    
    sigma <- data.frame(diag=rep(0,K))
    for (i in 1:G) {
      sigma <- sigma + colSums(t((t(X) - miu[,i])^2)*eta_p[,i])
    }
    sigma <- sigma/N
    temp_sig <- sigma
    for (i in 1:(G-1)) {
      temp_sig <- cbind(temp_sig, sigma)
    }
    
    for (i in 1:G) {
      miu[,i] <- colSums(eta_p[,i]*X)/sum(eta_p[,i])
    }
    
    for (i in 1:G) {
      miu[,i] <- sign(miu[,i])*quzheng(abs(miu[,i]) - temp_sig[,i]*lambda/colSums(tau)[i]/abs(miu[,i]))
    }
    
    b <- b + 1
  }
  
  ###mBIC
  q <- K - sum((rowSums(miu==0)<2))
  de <- G + K + G*K - 1 - q
  eta <- abs(eta-0.0000001)
  mBIC <- (-2)*(sum(log(colSums(pai*t(f))*(1-eta)))+sum(log(rowSums(z*f)*eta)))+log(N)*de
  
  result <- list(miu = miu, sigma = sigma, pai = pai, b = b, tau = tau, mBIC = mBIC, f = f, left = left)
  return(result)
}





