#######################################################
#version 1.0
#Edited based on 
#http://localhost:8888/notebooks/project1/5%20simulation%20codes.ipynb
#######################################################
# fixed epl:  epl = matrix(var(Y)/SNR,1,B)
#------load library-------
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
nj = 100 #repeat times

# elm <- rep(NA,2)
# response <- rep(list(elm),nj*3*3)
response  <- matrix(NA,nj*3*3,2)
#response = matrix(0,nj,1)
#factors
#no.1:SNR=var(Y)/var(epso) is used(1,5,10).NOTE:Z = Y+ (sqrt(var(Y))/sqrt(SNR))*rnorm(N)
#no.2: quantiles(.1,.5,.9);

#SNR=c(1,5,10)
fct1 = matrix(0,nj*3*3,1)
fct2 = matrix(0,nj*3*3,1)

count = 0

#subsample size for CV, i.e. 20% of all location for hold-out 
sub = floor(N*0.8)
sub_hd = N-sub


n = 112
B = 10000
burnin = 2000


# #
# SNR=3
# epsilon=0.5
# j=1

for(SNR in c(3,5,10)){
  for(epsilon in c(0.1,0.5,0.9)){ #epsilon stands for quantile now 
    #        for(epsilon in c(206,208,210)){
    #                for (esti in 0:1){
    for(j in 1:nj){
      
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


      
      ##initializations
      
      # Ws = matrix(0,n,B)
      # betas = matrix(0,3,B)
      # taus = matrix(var(Y)/SNR,1,B)
      # sigma2s = matrix(1,1,B)
      # phis = matrix(1,1,B)
      # 
      # Ws_sub = matrix(0,sub,B)
      # betas_sub = matrix(0,3,B)
      # taus_sub = matrix(var(Y)/SNR,1,B)
      # sigma2s_sub = matrix(1,1,B)
      # phis_sub = matrix(1,1,B)
       CVtrack = matrix(1,1,B)
       yhat_alltrn= matrix(0,n,B)
       yhat_hd = matrix(0,sub_hd,B)
      # 
      # Ws_sub2 = matrix(0,sub,B)
      # betas_sub2 = matrix(0,3,B)
      # taus_sub2 = matrix(var(Y)/SNR,1,B)
      # sigma2s_sub2 = matrix(1,1,B)
      # phis_sub2 = matrix(1,1,B)
      
      #----- initial-------
       p = 3 #dim(X)[2]
       r = 90 #dim(So)[2]
       epl =  var(Y)/SNR
       
       Etas = matrix(0,r,B) 
       betas = matrix(0,p,B)
       taus =  matrix(1,1,B) #variance for Xi
       Tso = matrix(0,sub,B)   #Xi for obs
       Tsm = matrix(0,sub_hd,B)   #Xi for hd
       K = riwishart(500, diag(r))
      
       Etas2 = matrix(0,r,B) 
       betas2 = matrix(0,p,B)
       taus2 =  matrix(1,1,B) #variance for Xi
       #Tso2 = matrix(0,sub,B)   #Xi for obs
       #Tsm2 = matrix(0,sub_hd,B)   #Xi for hd
       
       inx = matrix(0,sub,B)
       Ts = matrix(0,n,B)
       Ts2 = matrix(0,n,B)
       
     
      # #Gibbs untruncated------
      # for (i in 2:B){
      #   #update W
      #   Hphi = exp(-phis[i-1]*DistMat)
      #   Wcovar<-solve((1/sigma2s[i-1])*solve(Hphi) + (1/taus[i-1])*diag(n))  
      #   Wmean<-(1/taus[i-1])*(Wcovar%*%(Z-X%*%betas[,i-1]))
      #   Ws[,i]=mvrnorm(1, Wmean, Wcovar, tol = 1e-3)
      #   #Ws2[,i] = Ws[,i]
      #   
      #   #update betas
      #   BetAcovar<-solve((1/taus[i-1])*t(X)%*%X + (1/10)*diag(p))
      #   BetAmean<-(1/taus[i-1])*BetAcovar%*%(t(X)%*%(Z-Ws[,i]))             #*******find some wrong: lost 1/taus??19/05
      #   betas[,i]=mvrnorm(1, BetAmean, BetAcovar, tol = 1e-3)
      #   #betas2[,i]=betas[,i]
      #   
      #   #update taus
      #   #alphaTau = 1+n/2
      #   #betaTau = 0.01 + 0.5*(t(Z-X%*%betas[,i] - Ws[,i])%*%(Z-X%*%betas[,i] - Ws[,i]))
      #   #taus[i]=rigamma(1,alphaTau,betaTau)
      #   #taus2[i]=taus[i]
      #   
      #   #update sigma2s
      #   alphasigma2s = 1+n/2
      #   betasigma2s = 0.01 + 0.5*(t(Ws[,i])%*%solve(Hphi)%*%Ws[,i])
      #   sigma2s[i]=rigamma(1,alphasigma2s,betasigma2s)
      #   #sigma2s2[i] = sigma2s[i]
      #   
      #   #update phis
      #   probs = matrix(0,6,1)
      #   for (b in 1:6){
      #     Hphi = exp(-phi[b]*DistMat)
      #     probs[b] = dmvnorm(Ws[,i], matrix(0,n,1), Hphi, log=FALSE)
      #   }
      #   #print(c("itern",sum(probs)))
      #   if (sum(probs)>0){
      #     probs = probs/sum(probs)
      #     
      #     phis[i]=sample(phi, 1, replace = FALSE, prob = probs) #???why update like this way
      #     #phis2[i] = phis[i]
      #   }
      #   
      #   
      #   if (sum(probs)==0){
      #     phis[i]=phis[i-1]
      #     #phis2[i] = phis[i]#2020
      #   }
      #   
      #   #print(c("iter",phis[i]))
      #   
      #   # cH=sigma2s[i]*exp(-phis[i]*DistMat)
      #   # 
      #   # cHtaus=solve(cH+taus[i]*diag(n))
      #   # cH2=cH%*%cHtaus
      #   # 
      #   # yhat[,i]=X%*%betas[,i]+cH2%*%(Z-X%*%betas[,i])
      #   
      #   #SURE=t((Z-yhat[,i]))%*%(Z-yhat[,i])+2*(var(Y)/SNR)*sum(diag(cH2))
      #   #Suretrack[i]=SURE
      #   
      #   
      # }
      
      # #Trace plots
      # plot(mcmc(betas[3,]))
      # plot(mcmc(Ws[100,]))
      # plot(mcmc(phis[1,(burnin+1):B]))
      # #plot(mcmc(taus[1,]))
      # plot(mcmc(sigma2s[1,]))
      
      
      #cross-validation part------
      #cross-validation and then truncated by quantile
      for (i in 2:B){
        #subsample step
        subsmp_ind =  sample(seq_len(nrow(ZXlmat)), size = sub)
        subsmp =  ZXlmat[subsmp_ind, ]
        subhd =  ZXlmat[-subsmp_ind, ]
        inx[,i] = subsmp_ind
        
        Zo = subsmp[,1]
        Xo = subsmp[,2:4]
       
        Zm = subhd[,1]
        Xm = subhd[,2:4]
        
        S = eval_basis(basis,  Zldat) 
        So = eval_basis(basis,  Zldat[subsmp_ind,]) #So:to evalute S.o at BAUs. dim(S)=(N,r)
        Sm = eval_basis(basis,  Zldat[-subsmp_ind,])
 
        
        # So.T = t(So)
        # Xo.T = t(Xo)
        SotSo = crossprod(So) #t(So)%*%So
        XotXo=crossprod(Xo)
        
        
        
        #update
        
        # #update W------
        # Hphi_sub = exp(-phis_sub[i-1]*DistMat_sub)
        # Wcovar_sub = solve((1/sigma2s_sub[i-1])*solve(Hphi_sub) + (1/taus_sub[i-1])*diag(sub))  
        # Wmean_sub = (1/taus_sub[i-1])*(Wcovar_sub%*%(Z_sub-X_sub%*%betas_sub[,i-1]))
        # Ws_sub[,i] = mvrnorm(1, Wmean_sub, Wcovar_sub, tol = 1e-3)
        # Ws_sub2[,i] = Ws_sub[,i]
        # 
        # #update betas
        # BetAcovar_sub = solve((1/taus_sub[i-1])*t(X_sub)%*%X_sub + (1/10)*diag(p))
        # BetAmean_sub = (1/taus_sub[i-1])*BetAcovar_sub%*%(t(X_sub)%*%(Z_sub-Ws_sub[,i]))             #*******find some wrong: lost 1/taus??19/05
        # betas_sub[,i] = mvrnorm(1, BetAmean_sub, BetAcovar_sub, tol = 1e-3)
        # betas_sub2[,i] = betas_sub[,i]
        # 
        # #update taus
        # #alphaTau = 1+n/2
        # #betaTau = 0.01 + 0.5*(t(Z-X%*%betas[,i] - Ws[,i])%*%(Z-X%*%betas[,i] - Ws[,i]))
        # #taus[i]=rigamma(1,alphaTau,betaTau)
        # #taus2[i]=taus[i]
        # 
        # #update sigma2s
        # alphasigma2s_sub = 1+n/2
        # betasigma2s_sub = 0.01 + 0.5*(t(Ws_sub[,i])%*%solve(Hphi_sub)%*%Ws_sub[,i])
        # sigma2s_sub[i] = rigamma(1,alphasigma2s_sub,betasigma2s_sub)
        # sigma2s_sub2[i] = sigma2s_sub[i]
        
        #update eta-------
        Etcovar <- solve(SotSo*(1/epl) +solve(K))  
        Etmean <- (1/epl)*(Etcovar%*%crossprod(So,(Zo-Xo%*%betas[,(i-1)]-Tso[,(i-1)])) )#So.T%*%(Zo-Xo%*%betas[,(i-1)]-Tso[,(i-1)]))
        Etas[,i] <- as.matrix(mvrnorm(1, Etmean, Etcovar, tol = 1e-3))
        Etas2[,i] = Etas[,i]
        
        SoEt <- as.matrix(So%*%Etas[,i])
        
        #update betas
        BetAcovar <- solve((1/epl)*XotXo+1/10)
        BetAmean <- (1/epl)*BetAcovar%*%crossprod(Xo,(Zo-SoEt-Tso[,(i-1)]) ) #Xo.T%*%(Zo-SoEt-Tso[,(i-1)])
        betas[,i] <- mvrnorm(1, BetAmean, BetAcovar, tol = 1e-3)
        betas2[,i] =  betas[,i]
        
        ### #update Tso
        Tsocovar <- 1/(1/epl+1/taus[i-1])
        Tsomean <- (1/epl)*(Tsocovar*(Zo-Xo%*%betas[,i]-SoEt))
        Tso[,i] <- as.matrix(Tsomean+sqrt(Tsocovar)*rnorm(sub))
        #Tso2[,i] = Tso[,i]
        
        #update taus (variance for Ts)
        alphaTau <- 1+sub/2
        betaTau <- 0.01 + 0.5*crossprod(Tso[,i]) #(t(Tso[,i])%*%Tso[,i])
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
        
        #update Tsm
        Tsm[,i] <- mvrnorm(1, matrix(0,sub_hd,1), taus[i]*diag(sub_hd), tol = 1e-3)
        #Tsm2[,i] = Tsm[,i]
        #print(c("iter",phis[i]))
        
        #---- calc on hold-out dataset----
        # cH = sigma2s_sub[i]*exp(-phis_sub[i]*DistMat_hd)
        # cHtaus = solve(cH+taus_sub[i]*diag(sub_hd))
        # cH2 = cH %*% cHtaus
        # Xb_hd = X_hd %*% betas_sub[,i]
        # yhat_hd[,i] =Xb_hd + cH2 %*% (Z_hd-Xb_hd)
        # Zc_hd = Z_hd-yhat_hd[,i]
        # CV = crossprod(Zc_hd,Zc_hd)
        
        Ts[inx[,i],i] = Tso[,i]
        Ts[-inx[,i],i] = Tsm[,i]
        Ts2[,i] = Ts[,i]
        
        
        yhat_hd[,i] = as.matrix(Xm%*%betas[,i] + Sm%*%Etas[,i]+ Tsm[,i])
        CV = sqrt(mean(crossprod(Zm-yhat_hd[,i])))
        CVtrack[i] = CV
       # print(CVtrack[i])
        
      }
      
      
      epsil = quantile(CVtrack[(burnin+1):B],epsilon)
      #fm = data.frame(yhat[,1], matrix(0,N,B-1))
     # yhat_hd_trn = cbind(matrix(1,N,1),matrix(0,N,B-1))

      
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
        
        # cH = sigma2s_sub[i]*exp(-phis_sub[i]*DistMat)
        # cHtaus = solve(cH+taus_sub[i]*diag(n))
        # cH2 = cH%*%cHtaus
        # Xb = X%*%betas_sub[,i]
        # 
        # yhat_trnALL[,i] = Xb+cH2%*%Xb   
        
       
        
    
      
      #comparision------
      
      #yhat_hd_trn = apply(yhat_hd_trn[,(burnin+1):B],1,median)
      yhat_alltrn.final = rowMedians(yhat_alltrn[,(burnin+1):B])
      yhat_alltrn.final2 = rowMeans(yhat_alltrn[,(burnin+1):B])
      
      #untruncated
      Ys = X%*%betas2 + S%*%Etas2[,i]+ Ts2
      postmnY  = apply(Ys[,(burnin+1):B],1,mean) 
      
      count = count+1
      
      # response[[count]][1]  = crossprod(Y-yhat_alltrn.final) - crossprod(Y-postmnY)
      # response[[count]][2]  = crossprod(Y-yhat_alltrn.final2) - crossprod(Y-postmnY)
      response[count,1]  = crossprod(Y-yhat_alltrn.final) - crossprod(Y-postmnY)
      response[count,2]  = crossprod(Y-yhat_alltrn.final2) - crossprod(Y-postmnY)
      #response[count]=t(Y-yhat.final)%*%(Y-yhat.final) - t(Y-postmnY)%*%(Y-postmnY)
      fct1[count] = SNR
      fct2[count] = epsilon
      
      print(c(response[count,],SNR,epsilon,count,counter2))
      }
    }
  }
#}
#save.image("/gpfs/home/qz16b/fixed_taus_codes2020.RData")


#analysis -------

mydata = read.table("pro2_600.txt",na.strings="NA")
mydata1 = mydata[,c(2,3,4,5)]
colnames(mydata1) = c("rspMn","rspMn","SNR","esplison")

mydata1$SNR = factor(mydata1$SNR,
                     levels=unique(mydata1$SNR))
mydata1$esplison = factor(mydata1$esplison,
                          levels=unique(mydata1$esplison))
str(mydata1)

library(ggplot2)
library(latex2exp)

#response: median----
#boxplot-----

p1 = ggplot(mydata1, aes(x=SNR , y=rspMd, fill=SNR )) + 
  geom_boxplot(alpha=0.3) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="Blues")+
  geom_hline(yintercept = 0, color="red")+
  labs(y="ResponseMd", x = "SNR")
print(p1)

p2 = ggplot(mydata1, aes(x=esplison, y=rspMd, fill=esplison)) + 
  geom_boxplot(alpha=0.3) +
  theme(legend.position="none") +
  scale_fill_brewer(palette = "Greens")+
  geom_hline(yintercept = 0, color="red")+
  labs(y="ResponseMd", x = TeX("$d$-th Percentile"))

print(p2)


# interaction plot----
p3 = ggplot(mydata1, aes(x = SNR, y =rspMd)) +
  stat_summary(aes(group =esplison, color = esplison),
               fun = "mean", geom = "line", size = 1)+
  labs(y="ResponseMd", x = "SNR",color=TeX("$d$-th Percentile"))

print(p3)

#response mean-----
#boxplot-----

p4 = ggplot(mydata1, aes(x=SNR , y=rspMn, fill=SNR )) + 
  geom_boxplot(alpha=0.3) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="Blues")+
  geom_hline(yintercept = 0, color="red")+
  labs(y="ResponseMn", x = "SNR")
print(p4)

p5 = ggplot(mydata1, aes(x=esplison, y=rspMn, fill=esplison)) + 
  geom_boxplot(alpha=0.3) +
  theme(legend.position="none") +
  scale_fill_brewer(palette = "Greens")+
  geom_hline(yintercept = 0, color="red")+
  labs(y="ResponseMn", x = TeX("$d$-th Percentile"))

print(p5)


# #interaction plot----

# ggplot(data, aes(x = f1, y =rps)) +
#   stat_summary(aes(group =f2, color = f2),
#                  fun.y = "mean", geom = "line", size = 1)+
#   labs(y="Response", x = "SNR")+
#   scale_color_manual(name=TeX("$d_{th}$ Quantile"),
#                     labels = c("0.1",
#                                "0.5",
#                                "0.9"),
#                     values = c("0.1"="blue",
#                                "0.5"="brown",
#                                "0.9"="orange"))
#interaction plot
p6 = ggplot(mydata1, aes(x = SNR, y =rspMn)) +
  stat_summary(aes(group =esplison, color = esplison),
               fun = "mean", geom = "line", size = 1)+
  labs(y="ResponseMn", x = "SNR",color=TeX("$d$-th Percentile"))

print(p6)
