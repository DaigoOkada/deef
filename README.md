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

## Licence
GPL-2

## Author
Daigo Okada <dokada@genome.med.kyoto-u.ac.jp>
