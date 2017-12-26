genData <- function(n,p,beta)
{
X <- matrix(runif(n*p),nrow=n,ncol=p)
y <- rnorm(n,X%*%beta,2)
return(list(X=X,y=y))
}

genData.2 <- function(n,p,beta,rho=0)
{
  library("mvtnorm")
  x=c ( 1,rho^seq(1:(length( beta) - 1)))
  x= toeplitz(x)
  X <- matrix(rmvnorm(n,rep( 0 ,length( beta )),x),nrow=n , ncol=p )
  y <- rnorm(n,X%*%beta,2.2)
  return(list(X=X,y=y))
}

ridge <- function(Data,n.lam=301)
{
require(MASS)
lam <- c(0,exp(seq(log(1e-2),log(1e7),len=n.lam)))
fit <- lm.ridge(y~X,Data,lambda=lam)
return(coef(fit)[which.min(fit$GCV),])
}

bestsub <- function(Data)
{
require(leaps)
fit <- regsubsets(y~X,Data)
b <- numeric(ncol(Data$X)+1)
names(b) <- fit$xnames
bb <- coef(fit,which.min(summary(fit)[["bic"]]))
b[names(bb)] <- bb
return(b)
}

lasso <- function(Data)
{
require(glmnet)
cvfit <- cv.glmnet(Data$X,Data$y)
return(as.numeric(coef(cvfit,s=cvfit$lambda.min)))
}

alasso <- function(Data)
{
require(glmnet)
cvfit <- cv.glmnet(Data$X,Data$y)
cvfit2 <- cv.glmnet(Data$X,Data$y, penalty.factor = abs(1/as.numeric(coef(cvfit,s=cvfit$lambda.min)+(1/n))))
return(as.numeric(coef(cvfit2,s=cvfit2$lambda.min)))
}




enet <- function(Data)
{
require(glmnet)
cvfit <- cv.glmnet(Data$X,Data$y,alpha=0.5)
return(as.numeric(coef(cvfit,s=cvfit$lambda.min)))
}

cvfit <- cv.enet(Data$X,Data$y,lambda= seq(0.05,1,0.05), s=c(0,1,100), mode="fraction", trace = T, max.steps=300)


aenet <- function(Data)
{
require(glmnet)
cvfit <- cv.glmnet(Data$X,Data$y,alpha=0.5)
cvfit2 <- cv.glmnet(Data$X,Data$y, alpha=0.5, penalty.factor = abs(1/as.numeric(coef(cvfit,s=cvfit$lambda.min)+(1/n))))
return(as.numeric(coef(cvfit2,s=cvfit2$lambda.min)))
}
aenet2 <- function(Data){#New
  cvfit <- cv.glmnet(Data$X,Data$y,alpha=0)
  cva.list = list()
  min.mse = rep(0,11)
  for(i in seq(1,11,1)) {
    cv = cv.glmnet(Data$X,Data$y, alpha=(i-1)/10, penalty.factor = abs(1/as.numeric(coef(cvfit,s=cvfit$lambda.min)+(1/n))))
    cva.list[[i]] =as.numeric(coef(cv,s=cv$lambda.min))
    min.mse[i] =  cv$cvm 
  }
  min.ind = which.min(min.mse)
  cva.min = cva.list[[min.ind]]
  return(cva.min)
}

lm.after.aenet2 <- function(Data){#New
  require(glmnet)
  coef.in = as.integer(aenet2(Data)!=0)
  coef.in=coef.in[-1]   # Do not include intercept
  coef.in.index = which(coef.in != 0)
  result.vector = rep(0,ncol(Data$X))
  if(length(coef.in.index) == 0 ){
    return(result.vector)
  }
  vector = as.numeric(coef(lm(Data$y~Data$X[,coef.in.index])))[-1] 
  j = 1
  for(i in coef.in.index){
    result.vector[i] = vector[j]
    j = j +1
  }
  return(result.vector)
}

alasso2 <- function(Data){#New
  require(glmnet)
  cvfit <- cv.glmnet(Data$X,Data$y,alpha=0)
  cvfit2 <- cv.glmnet(Data$X,Data$y, alpha =1,penalty.factor = abs(1/as.numeric(coef(cvfit,s=cvfit$lambda.min)+(1/n))))
  return(as.numeric(coef(cvfit2,s=cvfit2$lambda.min)))
}  

lm.after.alasso2 <- function(Data){ #new
  require(glmnet)
  coef.in <- as.integer(alasso2(Data)!=0)
  coef.in=coef.in[-1]   # Do not include intercept
  coef.in.index = which(coef.in != 0)
  result.vector = rep(0,ncol(Data$X))
  if(length(coef.in.index) == 0 ){
    return(result.vector)
  }
  vector = as.numeric(coef(lm(Data$y~Data$X[,coef.in.index])))[-1] 
  j = 1
  for(i in coef.in.index){
    result.vector[i] = vector[j]
    j = j +1
  }
  return(result.vector)
}

 aladlasso <- function(Data){
 require(quantreg)
        tempcfQR1<- coef(rq(Data$y~Data$X,.5,method="lasso",lambda = 1))
			  lam= log(length(Data$y))/ abs(tempcfQR1)
        lam[1]=0
        tempcfQR<- coef(rq(Data$y~Data$X,.5,method="lasso",lambda = lam))
     return(tempcfQR)
 }
aladlasso.2 <- function(Data,n.lam=301){#new
  require(quantreg)
  tempcfQR1 = coef(rq(Data$y~Data$X,.5,method="lasso",lambda = .1))
  lam = c(0,exp(seq(log(1e-2),log(1e7),len=n.lam)))
  tempcfQR = list()
  mse = rep(0,n.lam)
  for(i in 1:n.lam){
    tempcfQR[[i]] =  rq(Data$y~Data$X,.5,method="lasso",lambda = lam[i])
    mse[i] = sum((tempcfQR[[i]]$residuals)^2 )
  }
  min.ind = which.min(mse)
  result = coef(tempcfQR[[min.ind]])
  return(result)
}


require(foreach)
n=100
N <- 250
beta = c(2,-2 ,1,-1,0.5,0.2,-0.3,-0.15,rep(0,12))
#beta3 = c(2,-2,1-,1,0.5,0.2,-0.3,-0.15,rep (0,12))
p <- length(beta)
rh = 
results <- array(NA,dim=c(N,p,8+6,length(rh)))
dimnames=list(1:N,1:p,c("Subset","Lasso","Ridge","El-Net","Ad. Lasso", "Ad. El-Net","Ad. LAD-Lasso","OLS","aenet2","lm.after.aenet2","alasso2","lm.after.alasso2","aladlasso.2","FIRST"))

#Parralel
require('doParallel')
require('foreach')
dyn.load("/Users/saaryalov/Desktop/Classic_FIRST_STORM_code/FIRST.so")
registerDoParallel(cores=8)
getDoParWorkers()
set.seed(1)

#start time
strt<-Sys.time()
for(j in 1:4){
foreach (i=1:N) %dopar% 
{
#Data <- genData.2(n,p,beta,rho=rh)
Data <- genData.2(n,p,beta,rho=j/5)
results[i,,1,j] <- bestsub(Data)[-1]
results[i,,2,j] <- lasso(Data)[-1]
results[i,,3,j] <- ridge(Data)[-1]
results[i,,4,j] <- enet(Data)[-1]
results[i,,5,j] <- alasso(Data)[-1]
results[i,,6,j] <- aenet(Data)[-1]
results[i,,7,j] <- aladlasso(Data)[-1]
results[i,,8,j] <- coef(lm(y~X,Data))[-1]

results[i,,9,j] <- aenet2(Data)[-1]
results[i,,10,j] <- lm.after.aenet2(Data)
results[i,,11,j] <- alasso2(Data)[-1]
results[i,,12,j] <- lm.after.alasso2(Data)
results[i,,13,j] <- aladlasso.2(Data)[-1]
results[i,,14,j] <- first(Data$y,Data$X,1,10)$estimates



#displayProgressBar(i,N)

}
}
print(Sys.time()-strt)



results.p <- results

set.seed(1)
results <- array(NA,dim=c(N,p,8+6))
strt<-Sys.time()
for(i in 1:N){
  Data <- genData.2(n,p,beta,rho=rh)
  results[i,,1] <- bestsub(Data)[-1]
  results[i,,2] <- lasso(Data)[-1]
  results[i,,3] <- ridge(Data)[-1]
  results[i,,4] <- enet(Data)[-1]
  results[i,,5] <- alasso(Data)[-1]
  results[i,,6] <- aenet(Data)[-1]
  results[i,,7] <- aladlasso(Data)[-1]
  results[i,,8] <- coef(lm(y~X,Data))[-1]
  
  results[i,,9] <- aenet2(Data)[-1]
  results[i,,10] <- lm.after.aenet2(Data)
  results[i,,11] <- alasso2(Data)[-1]
  results[i,,12] <- lm.after.alasso2(Data)
  results[i,,13] <- aladlasso.2(Data)[-1]
  results[i,,14] <- first(Data$y,Data$X,1,10)$estimates
}
print(Sys.time()-strt)
stopCluster(cl)

B <- apply(results[,,,1],2:3,mean)-beta
V <- apply(results,2:3,var)
MSE <- B^2+V
apply(MSE,2,sum)
MSE.table = apply(MSE,2,sum)
MSE.table = round(MSE.table,4)
print(rbind.data.frame(dimnames[[3]],MSE.table))



setwd("/Users/saaryalov/STP 530")
ylim <- range(results)

B <- apply(results[,,,1],2:3,mean)-beta
png("MC2.png",width = 2610, height = 3900,pointsize = 12,res=400)
myImagePlot(B)
dev.off()

B <- apply(results[,,,1],2:3,mean)-beta
png("MC3.png",width = 2610, height = 3900,pointsize = 12,res=400)
myImagePlot(B)
dev.off()

B <- apply(results[,,,2],2:3,mean)-beta
png("MC32.png",width = 2610, height = 3900,pointsize = 12,res=400)
myImagePlot(B)
dev.off()

B <- apply(results[,,,3],2:3,mean)-beta
png("MC33.png",width = 2610, height = 3900,pointsize = 12,res=400)
myImagePlot(B)
dev.off()

B <- apply(results[,,,4],2:3,mean)-beta
png("MC34.png",width = 2610, height = 3900,pointsize = 12,res=400)
myImagePlot(B)
dev.off()



V <- apply(results[,,,1],2:3,var)
png("MC4.png",width = 2610, height = 3900,pointsize = 12,res=400)
myImagePlot(V)
dev.off()

V <- apply(results[,,,2],2:3,var)
png("MC42.png",width = 2610, height = 3900,pointsize = 12,res=400)
myImagePlot(V)
dev.off()

V <- apply(results[,,,3],2:3,var)
png("MC43.png",width = 2610, height = 3900,pointsize = 12,res=400)
myImagePlot(V)
dev.off()

V <- apply(results[,,,4],2:3,var)
png("MC44.png",width = 2610, height = 3900,pointsize = 12,res=400)
myImagePlot(V)
dev.off()



png("MC5.png",width = 2610, height = 3900,pointsize = 12,res=400)
boxplot(results[,,2],col="bisque",ylim=ylim,main="Lasso (Scenario 1)")
abline(h=c(2,-2,1,-1), lty=2)
dev.off()

png("MC6.png",width = 2610, height = 3900,pointsize = 12,res=400)
boxplot(results[,,3],col="cadetblue",ylim=ylim,main="Ridge (Scenario 1)")
abline(h=c(2,-2,1,-1), lty=2)
dev.off()

png("MC7.png",width = 2610, height = 3900,pointsize = 12,res=400)
boxplot(results[,,4],col="beige",ylim=ylim,main="Elastic Net (Scenario 1)")
abline(h=c(2,-2,1,-1), lty=2)
dev.off()

png("MC8.png",width = 2610, height = 3900,pointsize = 12,res=400)
boxplot(results[,,5],col="azure",ylim=ylim,main="Adaptive Lasso (Scenario 1)")
abline(h=c(2,-2,1,-1), lty=2)
dev.off()

png("MC9.png",width = 2610, height = 3900,pointsize = 12,res=400)
boxplot(results[,,6],col="darkgoldenrod",ylim=ylim,main="Adaptive El. Net(Scenario 1)")
abline(h=c(2,-2,1,-1), lty=2)
dev.off()



png("MC10.png",width = 2610, height = 3900,pointsize = 12,res=400)
boxplot(results[,,7],col="darkgoldenrod2",ylim=ylim,main="Adaptive Lad Lasso(Scenario 1)")
abline(h=c(2,-2,1,-1), lty=2)
dev.off()



png("MC11.png",width = 2610, height = 3900,pointsize = 12,res=400)
boxplot(results[,,8],col="burlywood",ylim=ylim,main="OLS (Scenario 1)")
abline(h=c(2,-2,1,-1), lty=2)
dev.off()




myImagePlot <- function(x, ...){
     min <- min(x)
     max <- max(x)
     yLabels <- rownames(x)
     xLabels <- colnames(x)
     title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
       min <- Lst$zlim[1]
       max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
       yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
       xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
       title <- Lst$title
    }
  }
# check for null values
if( is.null(xLabels) ){
   xLabels <- c(1:ncol(x))
}
if( is.null(yLabels) ){
   yLabels <- c(1:nrow(x))
}

layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))

 # Red and green range from 0 to 1 while Blue ranges from 1 to 0
 ColorRamp <- rgb( seq(0,1,length=256),  # Red
                   seq(0,1,length=256),  # Green
                   seq(1,0,length=256))  # Blue
 ColorLevels <- seq(min, max, length=length(ColorRamp))

 # Reverse Y axis
 reverse <- nrow(x) : 1
 yLabels <- yLabels[reverse]
 x <- x[reverse,]

 # Data Map
 par(mar = c(3,5,2.5,2))
 image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
 ylab="", axes=FALSE, zlim=c(min,max))
 if( !is.null(title) ){
    title(main=title)
 }
axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
 axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
 cex.axis=0.7)

 # Color Scale
 par(mar = c(3,2.5,2.5,2))
 image(1, ColorLevels,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=ColorRamp,
      xlab="",ylab="",
      xaxt="n")

 layout(1)
}
# ----- END plot function ----- #
