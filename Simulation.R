# Author: Axel Grøn Roepstorff
# Date: Thu Mar 16 15:18:12 2023
# Initial Parameters and Necessary Packages --------------
rm(list=ls())
set.seed(69420)
bRi <- function(pack){
  if(sum(installed.packages()[,1]==pack)==1){
    library(paste(pack),character.only=TRUE)}
  else{
    install.packages(pack)
    library(paste(pack),character.only=TRUE)}  
}
#setwd('C:/Users/axelr/OneDrive - Aarhus universitet/Studiegruppe/6. Semester/Bachelor/Bachelor - Axel/Koder')
bRi('VGAM')
bRi('distributionsrd')
bRi('tidyverse')
bRi('copula')
bRi('PDFEstimator')
bRi('EnvStats')
bRi('calculus')
dStart <- 60000
dNG <- dStart/12
vProd <- rpareto(dStart,2)
fPerc <- ecdf(vProd)
vOmega <- rep(1,times=3)
dHousingShock <- menu(c('Yes','No'),TRUE, 'Simulate Housing Shock?')
dEdGrowth <- 2
# Relevant Functions ------------------------------------------------------

fPhiP <- function(dPhi){
  vPP <- sapply(lM,function(x){
    approxfun(density(x[,2],to=10))(dPhi)
  })
  return(vPP)
  }  
fPhiP1 <- function(y,x){approxfun(density(lM[[y]][,2],to=10))(x)}
fWage <- function(dPhi,vO=vOmega,vPB=vMeans,vEL=vL,vNN=vN,rho=2/3){
  t(vO*vPB%*%t(dPhi))*t(apply(t(apply(fPhiP(dPhi),1,function(x){x*vNN})),
                                     1,function(y){vEL/y}))^rho
}
fWbar <- function(Sam=lM){
  sapply(1:3,function(x){
  dW1 <- fWage(Sam[[x]][,2])[x]
  mean(dW1)
  })}
fH <- function(dPhi,rho=2/3){
  (t(t(dPhi))%*%vH)/t(vL^(rho)*(vN*t(sapply(estSh,function(x)dpareto(dPhi,1,x))))^(1-rho))
}
fTest <- function(x,y){t((x*t(dpareto(x,1,estSh[y]))*vN[y])^(2/3))}
fC <- function(dPhi,rho=2/3,dAlpha = 2/3){
  t(t((dAlpha*dPhi%*%t(vMeans*vTau))/t(((1+dMu)*(t(sapply(estSh,function(x)dpareto(dPhi,1,x)))*vN)^(1-rho)))^((dSubs-1)/dSubs))*
      (1/(vF*vM))^(1/dSubs)*vL^((dSubs-2)/dSubs))
}
fUtil <- function(dPhi,dBeg=vSh,dAlpha=2/3, dMax=FALSE){
  mA <- t(sapply(dPhi,function(x){
    c(rpareto(1,1,dBeg[1]/(x^vPsi[1])), rpareto(1,1,dBeg[2]/(x^vPsi[2])),rpareto(1,1,dBeg[3]/(x^vPsi[3])))}))
  mW <- fWage(dPhi)
  mRes <- mA*fC(dPhi)^(dAlpha)*fH(dPhi)^(1-dAlpha)
  if(dMax==TRUE){
    mMax <- matrix(c(sapply(1:nrow(mA),function(y){mA[y,which(max(mRes[y,])==mRes[y,])]}),
                     dPhi,round(fPerc(dPhi)+0.05,1)*10,
                     c(t(sapply(1:nrow(mA),function(x){
                       c(which(max(mRes[x,])==mRes[x,])[1],
                         mW[x,which(max(mRes[x,])==mRes[x,])][1],(mRes[x,which(max(mRes[x,])==mRes[x,])][1])/mA[x,which(max(mRes[x,])==mRes[x,])][1])}))))
                   ,ncol=6)
    return(mMax)
  }else{
    return(mRes)
  }
}


# Actual Simulation -------------------------------------------------------
### Remember to set vH to cover whether housing policy should be enacted or not.
  set.seed(125)
  vSh <- rep(5,3)
  vF <- 50
  dSubs <- 4
  dMu <- 1/(dSubs-1)
  vH <- rep(dStart/3,3)
  vTau <- c(1,1,1)
  vPsi <- c(1,1.5,2)
  vT <- c(1,1,1)
  vM <- c(dSubs,dSubs,dSubs)
  vA <- rpareto(dStart,1,2)
  vProd <- rpareto(dStart,1,3)
  vR <- c(sapply(1:3,function(x)rep(x,times=length(vA)/3)))
  mAll <- cbind(vA,vProd,round(fPerc(vProd)+0.050000001,1)*10,vR,rep(1,times=dStart),
                rep(2,times=dStart))
  lM <- split(mAll,mAll[,4])
  lM <- lapply(1:3,function(x){matrix(lM[[x]],ncol=6)})
  estSh <- sapply(lM,function(x){epareto(x[,2])$parameters[2]})
  vN <- sapply(lM,function(x){nrow(x)})
  vL <- sapply(1:3,function(Z){integral(fTest,list(y=Z,x=c(1,10)))$value})^(3/2)
  vMeans <- unlist(lapply(1:3,function(x){mean(lM[[x]][,2])}))
  mWage <- fWage(mAll[,2])
  mAll[,5] <- sapply(1:nrow(mAll),function(x){mWage[x,mAll[x,4]]})
  lM <- split(mAll,mAll[,4])
  lM <- lapply(1:3,function(x){matrix(lM[[x]],ncol=6)})
  vL <- sapply(1:3,function(Z){integral(fTest,list(y=Z,x=c(1,10)))$value})^(3/2)
  vN <- vapply(1:3,function(x){nrow(lM[[x]])},FUN.VALUE = 1)
  vMeans <- unlist(lapply(1:3,function(x){mean(lM[[x]][,2])}))
  mWage <- fWage(mAll[,2])
  dIter <- 0
  vCon <- 0
  mC <- matrix(dStart/3,ncol=3)
  estSh <- sapply(lM,function(x){epareto(x[,2])$parameters[2]})
  mEstSh <- estSh
  mSE1 <- matrix(0,nrow=100, ncol = 3)
  mSE2 <- matrix(0,nrow=100, ncol = 3)

  
  
  while(dIter<100){
    if(dIter>50 && dHousingShock==1){
      vH <- c(dStart/3*1.1,dStart/3,dStart/3*0.9) 
      fUtil <- function(dPhi,dBeg=vSh,dAlpha=2/3, dMax=FALSE){
        mA <- t(sapply(dPhi,function(x){
          c(rpareto(1,1,dBeg[1]/(x^vPsi[1])), rpareto(1,1,dBeg[2]/(x^vPsi[2])),rpareto(1,1,dBeg[3]/(x^vPsi[3])))}))
        mW <- fWage(dPhi)
        mRes <- mA*fC(dPhi)^(dAlpha)*fH(dPhi)^(1-dAlpha)
        if(dMax==TRUE){
          mMax <- matrix(c(sapply(1:nrow(mA),function(y){mA[y,which(max(mRes[y,])==mRes[y,])]}),
                           dPhi,round(fPerc(dPhi)+0.05,1)*10,
                           c(t(sapply(1:nrow(mA),function(x){
                             c(which(max(mRes[x,])==mRes[x,])[1],
                               mW[x,which(max(mRes[x,])==mRes[x,])][1],(mRes[x,which(max(mRes[x,])==mRes[x,])][1])/mA[x,which(max(mRes[x,])==mRes[x,])][1])}))))
                         ,ncol=6)
          return(mMax)
        }else{
          return(mRes)
        }
      }
    }
    fTest <- function(x,y){t((x*t(dpareto(x,1,estSh[y]))*vN[y])^(2/3))}
    vMeans <- unlist(lapply(1:3,function(x){mean(lM[[x]][,2])}))
    vN <- sapply(lM,nrow)
    vL <- sapply(1:3,function(Z){integral(fTest,list(y=Z,x=c(1,10)))$value})^(3/2)
    vWbar <- fWbar(lM)
    if(dIter>50 && dEdGrowth==1){
      dEdScale <- 3*(0.999^(dIter-50))
      vNG <- rpareto(dNG,1,dEdScale)
    }else{
    vNG <- rpareto(dNG,1,3)}
    vNG <- sapply(vNG,function(x){min(x,10)})
    mNG <- fUtil(vNG,dMax=TRUE)
    while(mode(mNG)!= 'numeric'){
      mNG <- fUtil(vNG,dMax=TRUE)
    }

    vKill <- round(runif(length(vNG),1,length(vA)))
    while (length(unique(vKill))<length(vKill)) {
      vKill[which(duplicated(vKill))] <- round(runif((length(vKill)-
                                                        length(unique(vKill))),
                                                     1,length(vA)))
    }
    mAllOld <- mAll
    mAll <- rbind(mAll[-vKill,],mNG)
    if(abs(nrow(mAll[which(mAll[,4]==1),])-nrow(mAllOld[which(mAllOld[,4]==1),]))<20 &&
       abs(nrow(mAll[which(mAll[,4]==2),])-nrow(mAllOld[which(mAllOld[,4]==2),]))<20 &&
       abs(nrow(mAll[which(mAll[,4]==3),])-nrow(mAllOld[which(mAllOld[,4]==3),]))<20){
      print('Convergence')
      vCon[dIter] <- 1
    } else{vCon[dIter] <- 0}
    mWage <- fWage(mAll[,2])
    mAll[,5] <- sapply(1:nrow(mAll),function(x){mWage[x,mAll[x,4]]})
    lM <- split(mAll,mAll[,4])
    lM <- sapply(1:3,function(x){matrix(lM[[x]],ncol=6)})
    dIter <- dIter+1
    mC <- rbind(mC,c(nrow(mAll[which(mAll[,4]==1),]),nrow(mAll[which(mAll[,4]==2),]),nrow(mAll[which(mAll[,4]==3),])))
    estSh <- sapply(lM,function(x){epareto(x[,2])$parameters[2]})
    mEstSh <- rbind(mEstSh,estSh)
    mSE1[dIter,] <- (vMeans-estSh/(estSh-1))^2
    mSE2[dIter,] <- (sapply(lM,function(x)var(x[,2]))-(estSh/((estSh-1)^2*(estSh-2))))^2
    cat('This is iteration number', dIter,'\n')
    
  }
  mP <- cbind(c(mC),c(rep(1,times=dIter+1),rep(2,times=dIter+1),rep(3,times=dIter+1)))


mRes <- data.frame(mAll)
colnames(mRes) <- c('Amenity','Productivity',
            'Productivity Type','Region','Wage','Utility')
v100 <- mRes %>% 
  count(Region)
v100
hadet <- mP



# Some plots --------------------------------------------------------------
pShares <- ggplot(mapping=aes(y=hadet[,1]/dStart,x=rep(1:(nrow(hadet)/3),times= 3),colour=factor(hadet[,2]))) +
  geom_line() + xlab('Period') + ylab('Population Share') + labs(color='Region')
if (dHousingShock == 1) {
  pShares <- pShares + geom_vline(xintercept=51)
}
pShares
plot(mEstSh[,1],type='l',lwd=1.5,main='Shape Parameter Estimates',xlab='Period',ylab='Estimated Shape',col='red',ylim=c(2.5,3.7))
lines(mEstSh[,2],type='l',lwd=1.5,col='green')
lines(mEstSh[,3],type='l',lwd=1.5,col='blue')
legend('topleft',legend=c('Region 1', 'Region 2', 'Region 3'),
       lty=c(1,1,1),col=c('red','green','blue'),cex=0.8)
if(dHousingShock == 1){
  abline(v=51)
}
plot(ecdf(lM[[3]][,2]))
lines(seq(1,10,length.out=nrow(lM[[3]])),
      ppareto(seq(1,10,length.out=nrow(lM[[3]])),1,estSh[3]),col='blue')

# Unconditional distributions ---------------------------------------------
dScale <- 5
dSNG <- 3
val2 <- 1.5
val3 <- 2
fun1var <- function(k,a,phi=dPhi){
  max(a,k)^(-dScale/phi)*
    distributionsrd::dpareto(a,dScale/(phi^val2))*distributionsrd::dpareto(k,dScale/(phi^(val3)))
}

fun2var <- function(k,a,phi=dPhi){
  max(a,k)^(-dScale/(phi^val2))*
    distributionsrd::dpareto(a,dScale/(phi))*distributionsrd::dpareto(k,dScale/(phi^(val3)))
}

fun3var <- function(k,a,phi=dPhi){
  max(a,k)^(-dScale/(phi^val3))*
    distributionsrd::dpareto(a,dScale/(phi^val2))*distributionsrd::dpareto(k,dScale/(phi))
}
fT <-  function(k,a,x)fun1var(k,a,x)*distributionsrd::dpareto(x,dSNG)
fT2 <- function(k,a,x)fun2var(k,a,x)*distributionsrd::dpareto(x,dSNG)
fT3 <- function(k,a,x)fun3var(k,a,x)*distributionsrd::dpareto(x,dSNG)
c(integral(fT,bounds=list(k=c(1,Inf),a=c(1,Inf),x=c(1,Inf)))$value,
  integral(fT2,bounds=list(k=c(1,Inf),a=c(1,Inf),x=c(1,Inf)))$value,
  integral(fT3,bounds=list(k=c(1,Inf),a=c(1,Inf),x=c(1,Inf)))$value)
probs <- sapply(seq(1,6,by=0.1),function(y){
  c(integral(fT,bounds=list(k=c(1,Inf),a=c(1,Inf),x=y))$value,
    integral(fT2,bounds=list(k=c(1,Inf),a=c(1,Inf),x=y))$value,
    integral(fT3,bounds=list(k=c(1,Inf),a=c(1,Inf),x=y))$value)
})
probsN <- data.frame(c(t(probs)),c(rep(seq(1,6,by=0.1),times=3)),
                     c(rep('Region 1',times=length(probs)/3),
                       rep('Region 2',times=length(probs)/3),
                       rep('Region 3',times=length(probs)/3)))
colnames(probsN) <- c('Probs','Prod','Region')
ggplot(probsN,aes(Prod,Probs,colour=Region)) + geom_line()


  
