#' DEEF
#'
#' This Function is to conduct DEEF method
#' @param disP          sample by grid plobability matrix, rowsums should be all 1.
#' @param ip_mat        optional. if sample by sample matrix are given, it is used as innerproduct matrix.
#' @return list that contains eigenvalue,Theta,Cx,Fx
#' @details
#' @keywords NA
#' @export
#' @examples
#'
DEEF <- function(disP, ip_mat=NULL){
  if(is.null(ip_mat)){
    ip_mat <- disP %*% t(disP)
  }
  theta_ip_est_mat <- log(ip_mat)/2

  #Eigenvalue decomposition
  eigen_out <- eigen(theta_ip_est_mat)
  eigen_value <- eigen_out[[1]]
  V <- eigen_out[[2]]
  Sigma <- diag(sqrt(abs(eigen_value)))
  Theta <- V %*% Sigma
  S <- diag(sign(eigen_value))

  #Calculate F'
  grid_num <- ncol(disP)
  sample_num <- nrow(disP)

  Psi <- matrix(NA,sample_num,grid_num)
  for(i in 1:sample_num){
    psi <- sum(sign(eigen_value) * Theta[i,]^2)
    Psi[i,] <- rep(psi,grid_num)
  }
  disP[disP==0] <- .Machine$double.xmin
  P_dash <- log(disP) + Psi
  Theta_dash <- cbind(Theta,rep(1,nrow(Theta)))
  F_dash <- MASS::ginv(Theta_dash) %*% P_dash

  #Output
  Theta <- Theta_dash[,-ncol(Theta_dash)]
  Cx <- F_dash[nrow(F_dash),]
  Fx <- F_dash[-nrow(F_dash),]
  result <- list(eigen_value,Theta,Cx,Fx)
  names(result) <- c("eigenvalue","Theta","Cx","Fx")
  return(result)

}

#' dist_repro
#'
#' This Function is to reconstruct distribution set from the output of DEEF
#' @param eigenvalue    Output of DEEF function, eigenvalue
#' @param Theta         Output of DEEF function, Theta
#' @param Cx            Output of DEEF function, Cx
#' @param Fx    　　　　Output of DEEF function, Fx
#' @param K             the number of the used top theta coordinates
#' @return samply by grid  reconstructed probability matrix
#' @details
#' @keywords NA
#' @export
#' @examples
#'
dist_repro <- function(eigenvalue,Theta,Cx,Fx, K){

  #F_dash and Theta_dash
  Theta_dash <- cbind(Theta,rep(1,nrow(Theta)))
  F_dash <- rbind(Fx,Cx)

  # Get top K coordinate' index
  use_idx <- order(abs(eigen_value),decreasing=T)[1:K]
  grid_num <- ncol(F_dash)
  bunpu_num <- nrow(Theta_dash)

  #Reconstruct Sita' and F'
  Theta_dash_sub <- Theta_dash[,c(use_idx, ncol(Theta_dash))]
  F_dash_sub <- F_dash[c(use_idx,nrow(F_dash)),]

  #Reconstruct P
  pos_eig_idx <- intersect(use_idx,which(eigen_value >= 0))
  neg_eig_idx <- intersect(use_idx,which(eigen_value < 0))
  tmp <- apply(Theta_dash[,pos_eig_idx,drop=F]^2,1,sum)-apply(Theta_dash[,neg_eig_idx,drop=F]^2,1,sum)
  Psi <- matrix(rep(tmp,grid_num),bunpu_num,grid_num)
  P_est <- exp(Theta_dash_sub %*% F_dash_sub - Psi)
  P_est <- diag(1/rowSums(P_est)) %*%  P_est

  return(P_est)
}

#' Distset2D
#'
#' This dataset of distribution set 2D
#' @format list contains the four objects
#' \describe{
#' \item{mu_sd}{900 by 2 matrix.The mean and sd of 900 member normal distributions.}
#' \item{ip_mat}{900 by 900 functional inner matrix.}
#' \item{grid}{The 10,000 grid values in discretization}
#' \item{P}{900 by 10,000 probability matrix.}
#' }
"Distset2D"
