#######################################################
#version 2.0
#######################################################
# fixed epl:  epl = matrix(var(Y)/SNR,1,B)
# fixed epsolion = 0.5
# f1: SNR  f2: subsample size for CV, i.e. 'subpct' of all location for hold-out
#------load library-------

# list.of.packages <- c("MASS","pscl","mvtnorm","coda",
#                       "matrixStats","FRK","sp","ggplot2",
#                       "gridExtra","INLA","splancs","matrixsampling","cIRT")
# 
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages)


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
library(INLA)
library(splancs)

library(matrixsampling)
library(cIRT)
#-------load data-------

# setwd("F:/ky/code6-gibbs_sampling")
# #setwd("C:/Users/bradley/Desktop/Students/Qingying")
# mydata<-read.table("lithology.txt",na.strings="NA")
# set.seed(2000)
# N = 112
# #covariates
# X0=matrix(1,N,1)
# X1 = mydata[,3] #Surf Elevation (ft amsl)
# X2 = mydata[,5] #A-B Elevation (ft amsl)
# X = cbind(X0,X1,X2)
# 
# #Data, locations, distance matrix
# Y = log(mydata[,4])   #log of Thickness(ft)   
# locations = mydata[,1:2] #Northing and Easting
# 
# loctrans= cbind((locations[,1]-min(locations[,1]))/max(locations[,1]-min(locations[,1])),(locations[,2]-min(locations[,2]))/max(locations[,2]-min(locations[,2])))
# locations=loctrans
load("C:/00_EE/old_laptop/ky/dropbox/0_UPDATE_zqy/codes/0902model1_a_FIXED/simulation_20_1.RData")
rm(list=setdiff(ls(), c("X","locations","Y","N")))

colnames(locations) = c("North","East")

#DistMat = as.matrix(dist(locations, method = "euclidean"))
p = dim(X)[2]


#loops------
epsilon = 0.5
nj = 100 #repeat times

response  <- matrix(NA,nj*3*3,2)

fct1 = matrix(0,nj*3*3,1)
fct2 = matrix(0,nj*3*3,1)

count = 0

# #subsample size for CV, i.e. 20% of all location for hold-out 
# sub = floor(N*0.8)
# sub_hd = N-sub


n = 112
B = 10000
burnin = 2000


# # #
# SNR=3
# subpct=0.25
#  j=1

for(SNR in c(3,5,10)){
  epl =  var(Y)/SNR
  for(subpct in c(0.1,0.25,0.5)){ #percentage for holdout sample size
    for(j in 1:nj){
     # setup-----
      Z = Y+(sqrt(var(Y))/sqrt(SNR))*rnorm(N) 
      ZXlmat = cbind(matrix(Z),X,locations)
      
      Zldat = data.frame(cbind(Z,locations))
	    coordinates(Zldat) = ~North+East
	    basis <- auto_basis(plane(),          # we are on the plane
	                        data = Zldat,       # data around which to make basis
	                        nres = 2,
	                        type = "bisquare",
	                        scale_aperture = 1.25)   # bisquare(radial) basis functions
	    #show_basis(basis)
	    
	    S = as.matrix(eval_basis(basis,  Zldat))
	    r = dim(S)[2]
	    
	    #1.full Gibbs part------
	    Etas = matrix(0,r,B) 
	    betas = matrix(0,p,B)
	    taus =  matrix(1,1,B) #variance for Xi
	    K = riwishart(500, diag(r))
	    Ts = matrix(0,n,B)
	    
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
	      Ts[,i]  <- as.matrix(Tsmean+sqrt(Tscovar)*rnorm(n)) 
	      
	      #update betas
	      BetAcovar <- solve((1/epl)*XtX+1/10)
	      BetAmean <- (1/epl)*BetAcovar%*%crossprod(X,(Z-SEt-Ts[,i]) ) #Xo.T%*%(Zo-SoEt-Tso[,(i-1)])
	      betas[,i] <- mvrnorm(1, BetAmean, BetAcovar, tol = 1e-3)

	      #update taus (variance for Ts)
	      alphaTau <- 1+n/2
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
	    
	    
      #2.method part------
	    
	    #initial
	    CVtrack = matrix(1,1,B)
	    yhat_alltrn= matrix(0,n,B)

	    Etas = matrix(0,r,B) 
	    betas = matrix(0,p,B)
	    taus =  matrix(1,1,B) #variance for Xi
	    K = riwishart(500, diag(r))
	    Ts = matrix(0,n,B)
	    
	    Etas2 = matrix(0,r,B) 
	    betas2 = matrix(0,p,B)
	    taus2 =  matrix(1,1,B) #variance for Xi
	    Ts2 = matrix(0,n,B)
	    
      #cross-validation and then truncated by quantile
	    #5-min with un-truncation
      for (i in 2:B){
        #subsample step
        sub_hd = floor(n*subpct)
        sub = n-sub_hd
        
        subsmp_ind =  sample(seq_len(nrow(ZXlmat)), size = sub)
        subsmp =  ZXlmat[subsmp_ind, ]
        subhd =  ZXlmat[-subsmp_ind, ]
        # inx[,i] = subsmp_ind
        inx = subsmp_ind
        
        Zo = subsmp[,1]
        Xo = subsmp[,2:4]
       
        Zm = subhd[,1]
        Xm = subhd[,2:4]
        
       
        So = eval_basis(basis,  Zldat[subsmp_ind,]) #So:to evalute S.o at BAUs. dim(S)=(N,r)
        Sm = eval_basis(basis,  Zldat[-subsmp_ind,])
 
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
        
        # update Ts
        Tsocovar <- 1/(1/epl+1/taus[i-1])
        Tsomean <- (1/epl)*(Tsocovar*(Zo-XobetasB-SoEt))
        Ts[inx,i]  <- as.matrix(Tsomean+sqrt(Tsocovar)*rnorm(sub)) #Tso[,i]
        Ts[-inx,i] <- mvrnorm(1, matrix(0,sub_hd,1), taus[i]*diag(sub_hd), tol = 1e-3)
        Ts2[,i] = Ts[,i]  
        
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
        # Tsm[,i] <- mvrnorm(1, matrix(0,sub_hd,1), taus[i]*diag(sub_hd), tol = 1e-3)
        
        
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
          #Ts[inx[,i],] = Tso[,i]
          #Ts[-inx[,i],] = Tsm[,i]
          yhat_alltrn[,i] = as.matrix(X%*%betas[,i] + S%*%Etas[,i]+ Ts[,i])
          
      }
      
      
      #comparision------
      
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
      #response[count]=t(Y-yhat.final)%*%(Y-yhat.final) - t(Y-postmnY)%*%(Y-postmnY)
      fct1[count] = SNR
      fct2[count] = subpct
      
      print(c(response[count,],SNR,subpct,count,counter2))
      }
    }
  }
#}
save.image("/gpfs/home/qz16b/prj2V2_0.RData")


# #analysis -------
# # pro2_600 only contains SNR=3,5 quantile=0.1,0.5,0.9
# mydata = read.table("pro2_600.txt",na.strings="NA")
# mydata1 = mydata[,c(2,3,4,5)]
# colnames(mydata1) = c("rspMn","rspMn","SNR","esplison")
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
# #response: median
# #boxplot
# 
# p1 = ggplot(mydata1, aes(x=SNR , y=rspMd, fill=SNR )) + 
#   geom_boxplot(alpha=0.3) +
#   theme(legend.position="none") +
#   scale_fill_brewer(palette="Blues")+
#   geom_hline(yintercept = 0, color="red")+
#   labs(y="ResponseMd", x = "SNR")
# print(p1)
# 
# p2 = ggplot(mydata1, aes(x=esplison, y=rspMd, fill=esplison)) + 
#   geom_boxplot(alpha=0.3) +
#   theme(legend.position="none") +
#   scale_fill_brewer(palette = "Greens")+
#   geom_hline(yintercept = 0, color="red")+
#   labs(y="ResponseMd", x = TeX("$d$-th Percentile"))
# 
# print(p2)
# 
# 
# # interaction plot
# p3 = ggplot(mydata1, aes(x = SNR, y =rspMd)) +
#   stat_summary(aes(group =esplison, color = esplison),
#                fun = "mean", geom = "line", size = 1)+
#   labs(y="ResponseMd", x = "SNR",color=TeX("$d$-th Percentile"))
# 
# print(p3)


