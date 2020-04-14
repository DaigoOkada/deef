deef - Decomposition into Extended Exponential Family in R
====

## Description
This R packages conduct decomposition into Extended Exponential Family (DEEF).This statistical method is described in our research paper titled "Decomposition of arbitrary sets of distributions inextended exponential family form for distinguishing multiple expression profiles of single-cell populations and visualizing their dynamics".We will update the documentation on the application of this package.


## Requirement
MASS (>= 7.3.50)

## Install
```{r}
install.packages("devtools")
devtools::install_github("DaigoOkada/deef")
```

## Overview

### Apply DEEF to normal distribution set

First, we applied DEEF to the Normal distribution set of 900 normal distributions with the different mean and sd.
This distribution set is the sampe example of the article and available as built-in dataset Distset2D.

```{r}
library(deef)
P <- Distset2D$P
ip_mat <- Distset2D$ip_mat
result <- DEEF(Distset2D$P,Distset2D$ip_mat)
eigen_value <- result$eigenvalue
Theta <- result$Theta
Fx <- result$Fx
Cx <- result$Cx
```

### Reconstruction using top three coordinates

dist_repro fucntion reconstracte the distribution set using the output of DEEF function.

```{r}
P_est <- dist_repro(eigenvalue,Theta,Cx,Fx, K=3)
```

### Visulalization
```{r}
library(scatterplot3d)
mu_sd <- Distset2D$mu_sd
grid <- Distset2D$grid
colpal <- colorRampPalette(c("red","blue"))
col <- colpal(nrow(mu_sd))

#Visualize original parameter structure
plot(mu_sd,xlab="mean",ylab="sd",col=col,main="Original Parameters",cex.lab=2,cex.main=4,cex.axis=2,pch=16)

#Visualize image of original disribution set
matplot(grid,t(P),xlab="x",type="l",ylab="",col=col,main="Distiburion Images",cex.lab=2,cex.main=4,cex.axis=2,yaxt="n")

#Visualiza Theta coordinate system
scatterplot3d(Theta[,order(abs(eigen_value),decreasing=T)[1:3]],color=col,angle=70,xlab="Theta1",ylab="Theta2",zlab="Theta3",tick.marks=FALSE,cex.lab=2,cex.main=4,pch=16)

#Visualiza the reconstructed distribution
matplot(grid,t(P_est),xlab="x",type="l",ylab="",col=col,main="Reconstructed",cex.lab=2,cex.main=4,cex.axis=2,yaxt="n")
```

### Check F and C
```{r}
Ft <- cbind(t(Fx)[,order(abs(eigen_value),decreasing=T)[1:3]],Cx)
labels <- c("F1","F2","F3","Cx")
cols <- c("red","blue","purple","black")
matplot(seq(from=0,to=1,length=1000),Ft,col=cols,type="l",lty="solid",xlab="",main="1D",cex.main=2,xaxp=c(0, 1, 1),ylab="")
legend("topright", legend = labels, col = cols, lty = "solid",ncol = 1)
```

### Create Probability matrix P from cytometry dataset.
In application of DEEF to cytometry dataset amalysis, Probability matrix P have to be created with the deciding the grids and 
k-nearest neighbor estimation of the probability values of the grids.
The following code is the sample code to create matrix P from GvHD dataset which contain 35 cytometry samples, with the setting that the number of markers (d) is 2, the number of grid for each marker is 100 and alpha parameter (contorolling the range) is 0.15.
The detail  of the methodology is written in our article.

```{r}
#Preproseccing of the cytometry data
library(flowCore) #install from bioconducter
library(TDA) #install from CRAN
data(GvHD)
n <- length(GvHD)
expr_list <- list()
for(i in 1:n){
  samp <- GvHD[[i]]
  expr_list[[i]] <- asinh(samp@exprs[,c('FL1-H','FL2-H')])
}

#the num of marker(d)=2,the num of grids (m)=100
d <- 2
m <- 100
min_list <- lapply(1:d,function(x){NULL})
max_list <- lapply(1:d,function(x){NULL})
alpha <- 0.15
for(i in 1:n){
  expr <- expr_list[[i]]
  for(j in 1:d){
    min_list[[j]] <- c(min_list[[j]],quantile(expr[,j],alpha))
    max_list[[j]] <- c(max_list[[j]],quantile(expr[,j],1-alpha))
  }
}

seq_list <- list()
for(j in 1:d){
  total_min <- min(min_list[[j]])
  total_max <- max(max_list[[j]])
  seq_list[[j]] <- seq(from=total_min,to=total_max,length=m)
}

#Generate Grid matrix
code <- "x_grid <- expand.grid("
for(i in 1:d){
  if(i != d) code <- paste0(code,"seq_list[[",i,"]],")
  if(i == d) code <- paste0(code,"seq_list[[",i,"]])")
}
eval(parse(text=code))

#knn estimate and generate P
P <- matrix(NA,n,m^d)
for(i in 1:n){
  expr <- expr_list[[i]]
  knni <- knnDE(expr, x_grid, k=100)
  P[i,] <- knni/sum(knni)
}
```

## Licence
GPL-2

## Author
Daigo Okada <dokada@genome.med.kyoto-u.ac.jp>
