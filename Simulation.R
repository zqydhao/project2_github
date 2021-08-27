#######################################################
#version 3.1
#change simulated sample size to be 100,1000

#simuated data Y = X\beta+w; w~N(0, exp(-thetaV*DistMat))
#Gibbs Z = X\beta+S\eta+\Ts

#fixed the error when updating Ts[-inx,i] according to
# the fiveminutes paper algorithm1
######################################################

#1#-----load library-------
library(MASS)
library(pscl)
library(mvtnorm)
library(coda)
library(matrixStats) #for rowMedians()
#
library(FRK)
library(sp)
library(ggplot2)
library(gridExtra)
#library(INLA)
library(splancs)

library(matrixsampling)
library(cIRT)

#2#-----simulated data-------

# for (N in c(1e2,1e3,1e4)){        #simulated dataset size
simuData = function(seed, N=1e3, thetaV = 3, psV = 1){
  set.seed(seed)
  #betastar = mvrnorm(1, matrix(0,3,1), diag(3), tol = 1e-3) #!reason i stucked
  betastar = matrix(rnorm(3),3,1)
  ind = seq_len(N)
  x1 = ind/N
  X = cbind(1,x1,x1^2)
  Xbstar = X%*%betastar
  
  #locations = matrix(runif(2*N),N,2)
  locations = round( matrix(runif(2*N),N,2),3)
  colnames(locations) = c("North","East")
  DistMat = as.matrix(dist(locations, method = "euclidean"))
  
  # thetaV = 3
  covV = psV*exp(-thetaV*DistMat)  #effective range= 3/phi
  V = mvrnorm(1, matrix(0,N,1),covV, tol = 1e-3)
  
  c = as.numeric(sqrt(var(V)/var(Xbstar)))
  Xbeta = c*Xbstar
  
  Y = as.vector(Xbeta+V)
  
  # #Divide the screen in 2 line and 2 columns
  # par(
  #   mfrow=c(2,2), 
  #   oma = c(0, 0, 2, 0)
  # ) 
  # #Make the margin around each graph a bit smaller
  # par(mar=c(4,4,4,4))
  # plot(ind, V)
  # plot(ind, Xbeta)
  # plot(ind,Y)
  # hist(Y)
  # mtext(paste("N=",N,"thetaV =",fractions(thetaV),"seed =",seed,"psV=",psV) ,outer = T, cex = 1.5 )
  
  backlist =  list(X=X,locations=locations,Y=Y,X=X)
  return(backlist)
}

# a1 = simuData(seed = 1234, N=1e3, thetaV = 1/3)
# a2 = simuData(seed = 1234, N=1e3, thetaV = 3)
# a3 = simuData(seed = 1234, N=1e3, thetaV = 9)
# a4 = simuData(seed =3000, N=1e3, thetaV = 1/3)

#aa1 = simuData(seed = 3000, N = 1e2, thetaV = 3)
aa2 = simuData(seed = 3000, N = 1e3, thetaV = 3)

#### 3-ploting#
# a = a5
# a_df = as.data.frame(cbind(as.matrix(a[["Y"]]),a[["locations"]]))
# # ggplot(a_df,aes(x=North,y=East))+
# #   geom_point(col="red")+
# #   coord_equal()
# library(plotly)
# plot_ly(x=a_df[,"North"], y=a_df[,"East"], z=a_df[,"V1"], type="scatter3d", mode="markers",color=a_df[,"North"])

#3#-----Gibbs part-------
# 
a = aa2
X = a[["X"]]
Y = a[["Y"]]
locations = a[["locations"]]
N = length(Y)

#  
p = dim(X)[2]
B = 10000
burnin = 2000
epsilon = 0.5


nj=100#repeat times
response  <- matrix(NA,nj*3*3,2)
fct1 = matrix(0,nj*3*3,1)
fct2 = matrix(0,nj*3*3,1)

count = 0



# SNR=3
# subpct=0.1
#nj=5

for(SNR in c(3,5,10)){
  epl =  var(Y)/SNR
  #set.seed(10) # by doing so,no matter you run together or separately for SNR you CANT get the same result(hpc cases)
  for(subpct in c(0.1,0.25,0.5)){ #percentage for holdout sample size
    for(j in 1:nj){
      # setup-----
      Z = Y+sqrt(epl)*rnorm(N) 
      ZXlmat = cbind(matrix(Z),X,locations)
      
      Zldat = data.frame(cbind(Z,locations))
      coordinates(Zldat) = ~North+East
      basis <- auto_basis(plane(),              # we are on the plane
                          data = Zldat,         # data around which to make basis
                          nres = 2, # max_basis = 100,
                          type = "bisquare",    # bisquare(radial) basis functions
                          scale_aperture = 1.25,
      )   
      #show_basis(basis)
      
      S = as.matrix(eval_basis(basis,  Zldat))
      r = dim(S)[2]
      
      #3.1#full Gibbs part----
      Etas = matrix(0,r,B) 
      betas = matrix(0,p,B)
      taus =  matrix(1,1,B) #variance for Xi
      K = riwishart(500, diag(r))
      Ts = matrix(0,N,B)
      
      StS = crossprod(S) #t(S)%*%S
      XtX = crossprod(X)
      
      for(i in 2:B){
        #update-------
        #update eta
        XbetasB = X%*%betas[,(i-1)]
        Etcovar <- solve(StS*(1/epl) +solve(K))  
        Etmean <- (1/epl)*(Etcovar%*%crossprod(S,(Z-XbetasB-Ts[,(i-1)])) )
        Etas[,i] <- as.matrix(mvrnorm(1, Etmean, Etcovar, tol = 1e-3))
        
        SEt <- as.matrix(S%*%Etas[,i])
        
        # update Ts
        Tscovar <- 1/(1/epl+1/taus[i-1])
        Tsmean <- (1/epl)*(Tscovar*(Z-XbetasB-SEt))
        Ts[,i]  <- as.matrix(Tsmean+sqrt(Tscovar)*rnorm(N)) 
        
        #update betas
        BetAcovar <- solve((1/epl)*XtX+1/10)
        BetAmean <- (1/epl)*BetAcovar%*%crossprod(X,(Z-SEt-Ts[,i]) ) #Xo.T%*%(Zo-SoEt-Tso[,(i-1)])
        betas[,i] <- mvrnorm(1, BetAmean, BetAcovar, tol = 1e-3)
        
        #update taus (variance for Ts)
        alphaTau <- 1+N/2
        betaTau <- 0.01 + 0.5*crossprod(Ts[,i]) #(t(Tso[,i])%*%Tso[,i])
        taus[i]<- rigamma(1,alphaTau,betaTau)
        
        #update K
        Kdgree <- 500+r
        Kscale <- diag(r)+tcrossprod(Etas[,i]) #Etas[,i]%*%t(Etas[,i])
        K <- riwishart(Kdgree, Kscale)
        
      }
      
      # traditional full model without truncation
      Ys = X%*%betas + S%*%Etas + Ts
      postmnY  = apply(Ys[,(burnin+1):B],1,mean)
      
      
      #3.2#method part------
      
      #initial
      CVtrack = matrix(1,1,B)
      yhat_alltrn= matrix(0,N,B)
      
      Etas = matrix(0,r,B) 
      betas = matrix(0,p,B)
      taus =  matrix(1,1,B) #variance for Xi
      K = riwishart(500, diag(r))
      Ts = matrix(0,N,B)
      
      Etas2 = matrix(0,r,B) 
      betas2 = matrix(0,p,B)
      taus2 =  matrix(1,1,B) #variance for Xi
      Ts2 = matrix(0,N,B)
      
      #cross-validation and then truncated by quantile
      #5-min with un-truncation
      for (i in 2:B){
        #subsample step
        sub_hd = floor(N*subpct)
        sub = N-sub_hd
        
        subsmp_ind =  sample(seq_len(nrow(ZXlmat)), size = sub)
        subsmp =  ZXlmat[subsmp_ind, ]
        subhd =  ZXlmat[-subsmp_ind, ]
        # inx[,i] = subsmp_ind
        inx = subsmp_ind
        
        Zo = subsmp[,1]
        Xo = subsmp[,2:4]
        
        Zm = subhd[,1]
        Xm = subhd[,2:4]
        
        
        So = as.matrix(eval_basis(basis,  Zldat[subsmp_ind,]) ) #So:to evalute S . dim(S)=(N,r)
        Sm = as.matrix(eval_basis(basis,  Zldat[-subsmp_ind,]))
        
        SotSo = crossprod(So) #t(So)%*%So
        XotXo=crossprod(Xo)
        
        #update-------
        #update eta
        XobetasB = Xo%*%betas[,(i-1)]
        
        Etcovar <- solve(SotSo*(1/epl) +solve(K))  
        Etmean <- (1/epl)*(Etcovar%*%crossprod(So,(Zo- XobetasB-Ts[inx,(i-1)])) )#So.T%*%(Zo-Xo%*%betas[,(i-1)]-Tso[,(i-1)]))
        Etas[,i] <- as.matrix(mvrnorm(1, Etmean, Etcovar, tol = 1e-3))
        Etas2[,i] = Etas[,i]
        
        SoEt <- as.matrix(So%*%Etas[,i])
        
        # update Tso
        # Tsocovar <- 1/(1/epl+1/taus[i-1])
        # Tsomean <- (1/epl)*(Tsocovar*(Zo-Xo%*%betas[,i]-SoEt))
        # Tso[,i] <- as.matrix(Tsomean+sqrt(Tsocovar)*rnorm(sub))
        
        # update Ts (using reverse jump)
        Tsocovar <- 1/(1/epl+1/taus[i-1])
        Tsomean <- (1/epl)*(Tsocovar*(Zo-XobetasB-SoEt))
        Ts[inx,i]  <- as.matrix(Tsomean+sqrt(Tsocovar)*rnorm(sub)) #Tso[,i]
       # Ts[-inx,i] <- mvrnorm(1, matrix(0,sub_hd,1), taus[i]*diag(sub_hd), tol = 1e-3)
       # Ts2[,i] = Ts[,i]  
        
        #update betas
        BetAcovar <- solve((1/epl)*XotXo+1/10)
        BetAmean <- (1/epl)*BetAcovar%*%crossprod(Xo,(Zo-SoEt-Ts[inx,i]) ) #Xo.T%*%(Zo-SoEt-Tso[,(i-1)])
        betas[,i] <- mvrnorm(1, BetAmean, BetAcovar, tol = 1e-3)
        betas2[,i] =  betas[,i]
        
        #update taus (variance for Ts)
        alphaTau <- 1+sub/2
        betaTau <- 0.01 + 0.5*crossprod(Ts[inx,i]) #(t(Tso[,i])%*%Tso[,i])
        taus[i]<- rigamma(1,alphaTau,betaTau)
        taus2[i] = taus[i]
        
        # ##update epl
        # alphaEpl <- 1+sub/2
        # betaEpl <- 1 + 0.5*(t(Zo-Xo%*%betas_cur-SoEt-Tso_cur)%*%(Zo-Xo%*%betas_cur-SoEt-Tso_cur))
        # epl_cur<- rigamma(1,alphaEpl,betaEpl)           
        
        #update K
        Kdgree <- 500+r
        Kscale <- diag(r)+tcrossprod(Etas[,i]) #Etas[,i]%*%t(Etas[,i])
        K <- riwishart(Kdgree, Kscale)
        
        # #update Tsm
        Ts[-inx,i] <- mvrnorm(1, matrix(0,sub_hd,1), taus[i]*diag(sub_hd), tol = 1e-3)
        Ts2[,i] = Ts[,i]  
        
        yhat_hd = as.matrix(Xm%*%betas[,i] + Sm%*%Etas[,i]+ Ts[-inx,i])
        CV = sqrt(mean(crossprod(Zm-yhat_hd)))
        CVtrack[i] = CV
        
      }
      # #Trace plots
      # plot(mcmc(betas[3,]))
      # plot(mcmc(Ws[100,]))
      # plot(mcmc(phis[1,(burnin+1):B]))
      # #plot(mcmc(taus[1,]))
      # plot(mcmc(sigma2s[1,]))
      
      epsil = quantile(CVtrack[(burnin+1):B],epsilon)
      
      counter2 = 0
      for (i in 2:B){
        if (CVtrack[i] > epsil) {
          Etas[,i] = Etas2[,(i-1)]
          betas[,i] = betas2[,(i-1)]
          #Tso[,i] = Tso2[,(i-1)]
          taus[i] = taus2[,(i-1)]
          #Tsm[,i] = Tsm2[,(i-1)]
          Ts[,i] = Ts2[,(i-1)]
          counter2 = counter2+1
        }
        
        yhat_alltrn[,i] = as.matrix(X%*%betas[,i] + S%*%Etas[,i]+ Ts[,i])
        
      }
      
      
      #4#-----comparision------
      
      # 5-min method with truncation(posterior median)
      #yhat_hd_trn = apply(yhat_hd_trn[,(burnin+1):B],1,median) 
      #yhat_alltrn.final2 = rowMeans(yhat_alltrn[,(burnin+1):B])
      yhat_all_Mtrn = rowMedians(yhat_alltrn[,(burnin+1):B])
      
      #5-min method with un-truncation (posterior mean)
      Ymethod = X%*%betas2 + S%*%Etas2 + Ts2  ##!!!##(mistake for version 1.0) when working on Version 2.0.should be Etas2 not Etas2[,i] 
      yhat_all_Munt  = apply(Ymethod[,(burnin+1):B],1,mean) 
      
      # # traditional full model without truncation
      # Ys = X%*%betas2 + S%*%Etas2 + Ts2
      # postmnY  = apply(Ys[,(burnin+1):B],1,mean) 
      
      count = count+1
      
      response[count,1]  = crossprod(Y-yhat_all_Mtrn) - crossprod(Y-postmnY)
      response[count,2]  = crossprod(Y-yhat_all_Munt) - crossprod(Y-postmnY)
      #print(crossprod(Y-postmnY))
      #response[count]=t(Y-yhat.final)%*%(Y-yhat.final) - t(Y-postmnY)%*%(Y-postmnY)
      fct1[count] = SNR
      fct2[count] = subpct
      
      print(c(response[count,],SNR,subpct,count,counter2))
    }
  }
}
#save.image("/gpfs/home/qz16b/prj2V3_1_Ne2.RData")


#5#------analysis -------
# mydata = read.table("C:/00_EE/research/project2_github/hpc_project2/results/simulation3_Ne2_3510.txt",fill = TRUE)
# #mydata = read.table("C:/00_EE/research/project2_github/hpc_project2/results/simulation3_Ne3_3510 .txt.txt",fill = TRUE)
# head(mydata)
# mydata1 = mydata[complete.cases(mydata),2:5]
# colnames(mydata1) = c("rspTrn","rspUntrn","SNR","esplison")
# 
# mydata1$SNR = factor(mydata1$SNR,
#                      levels=unique(mydata1$SNR))
# mydata1$esplison = factor(mydata1$esplison,
#                           levels=unique(mydata1$esplison))
# str(mydata1)
# 
# library(ggplot2)
# library(latex2exp)
# 
# #response: median+ Trn----
# #boxplot
# 
# p1 = ggplot(mydata1, aes(x=SNR , y=rspTrn, fill=SNR )) +
#   geom_boxplot(alpha=0.3) +
#   theme(legend.position="none") +
#   scale_fill_brewer(palette="Blues")+
#   geom_hline(yintercept = 0, color="red")+
#   labs(y="ResponseTrn", x = "SNR")
# print(p1)
# 
# p2 = ggplot(mydata1, aes(x=esplison, y=rspTrn, fill=esplison)) +
#   geom_boxplot(alpha=0.3) +
#   theme(legend.position="none") +
#   scale_fill_brewer(palette = "Greens")+
#   geom_hline(yintercept = 0, color="red")+
#   labs(y="ResponseTrn", x = TeX("$d$-th Percentile"))
# 
# print(p2)
# 
# 
# # interaction plot
# p3 = ggplot(mydata1, aes(x = SNR, y =rspTrn)) +
#   stat_summary(aes(group =esplison, color = esplison),
#                fun = "mean", geom = "line", size = 1)+
#   labs(y="ResponseTrn", x = "SNR",color=TeX("$d$-th Percentile"))
# 
# print(p3)
# 
# 
# 
# #response: mean Untrn----
# #boxplot
# 
# p4 = ggplot(mydata1, aes(x=SNR , y=rspUntrn, fill=SNR )) +
#   geom_boxplot(alpha=0.3) +
#   theme(legend.position="none") +
#   scale_fill_brewer(palette="Blues")+
#   geom_hline(yintercept = 0, color="red")+
#   labs(y="ResponseUntrn", x = "SNR")
# print(p4)
# 
# p5 = ggplot(mydata1, aes(x=esplison, y=rspUntrn, fill=esplison)) +
#   geom_boxplot(alpha=0.3) +
#   theme(legend.position="none") +
#   scale_fill_brewer(palette = "Greens")+
#   geom_hline(yintercept = 0, color="red")+
#   labs(y="ResponseUntrn", x = TeX("$d$-th Percentile"))
# 
# print(p5)
# 
# 
# # interaction plot
# p6 = ggplot(mydata1, aes(x = SNR, y =rspUntrn)) +
#   stat_summary(aes(group =esplison, color = esplison),
#                fun = "mean", geom = "line", size = 1)+
#   labs(y="ResponseUntrn", x = "SNR",color=TeX("$d$-th Percentile"))
# 
# print(p6)
