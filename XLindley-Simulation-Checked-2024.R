library(nleqslv)
library(truncnorm)
library(coda)
######################################################################
#                          Sample                                    #
######################################################################
rxlindley = function(n,theta){
  j=1
  t=c()
  while(j<=n){
    z=runif(1,0,1)
  model = function(u){
    v = numeric(1)
    v[1] = 1- (1+theta*u[1]/(1+theta)^2)*exp(-theta*u[1])-z
    v
  }
  res = nleqslv(c(.5),model)
  if(abs(res$fvec)<0.05) {t[j] = res$x; j=j+1}
  }
  t
}
 sh1=sh2=sh3=sh4=sh5=sh6=sh7=sh8=sh9=sh10=0
sh11=sh12=sh13=sh14=sh15=sh16=sh17=sh18=sh19=sh20=0
sh21=sh22=sh23=sh24=sh25=sh26=sh27=sh28=sh29=sh30=0
##############################################################
sh = 0
theta = 2
thetaml = thetaL = thetaU = thetamlm = thetaLm = thetaUm = c()
thetaS =thetaLc1 =thetaLc2 =thetaHPDL =thetaHPDU  =c()
thetaSm=thetaLc1m=thetaLc2m=thetaHPDLm=thetaHPDUm =c()
predU  = predL  = gglinexc2  = gglinexc1  = ggsel  = c()
predUm = predLm = gglinexc2m = gglinexc1m = ggselm = c()
m = 5
c1 =  0.5
c2 = -0.5
a1 = b1 = 0.1
sspar = c(0.9*theta)
records = c()
record = c()
rs = c()
j = 1
J = 10 
while (j<=J){
  plot(j,j,xlim = c(j-1,j+1))
  mm = m+1
  ttime = rep(1,mm)
  
  records[1] = rxlindley(1,theta)
  i = 2
  while(i<=mm){
    z = rxlindley(1,theta)
    jj = i-1
    if(z>records[jj]) ttime[jj] = ttime[jj]+1
    if(z<records[jj]){
      records[i] = z
      i = i+1
    }
    if(z == records[jj])   sh = sh+1 
  }
  
  rs[j] = records[mm]
  rm = records[m]
  record = records[1:m]
  recordm = records[1:(m-1)]
  time = ttime[1:m]
  time[m] = 1
  ######################################################################
  #                   MLE with records and times                      #
  ######################################################################
  xi = function (r, theta) (1+(theta*r)/(1+theta)^2)*exp(-r*theta)
  loglike  = function(par){
    theta = par[1]
    A = 2*m*log(theta)-2*m*log(1+theta)-theta*sum(record)+sum(log(2+theta+record))+
        sum((time-1)*log(xi(record,theta)))
    -A
  }
  opt1 = optim(par = sspar, fn = loglike, method = "BFGS", control = c(maxit = 10000), hessian = TRUE)
  if(opt1$convergence!=0) sh1=sh1+1
  if(opt1$par[1] <=0) sh2=sh2+1
  if(opt1$hessian<=0) sh3=sh3+1
  if(opt1$convergence==0 & opt1$par[1]>0 & opt1$hessian >0){
  thetaml[j] = opt1$par[1]
   thetaL[j] = max(0,thetaml[j]-1.96/sqrt(opt1$hessian))
  if(thetaml[j]-1.96/sqrt(opt1$hessian)<0) sh4 = sh4+1
  thetaU[j] =thetaml[j]+1.96/sqrt(opt1$hessian)
  
  ######################################################################
  #                   MLE with records without times                   #
  ######################################################################
  loglikem  = function(par){
    theta = par[1]
    A = 2*m*log(theta)-2*m*log(1+theta)-theta*sum(record)+sum(log(2+theta+record))+
        sum(log(1-xi(recordm,theta)))
    -A
  }
  opt2 = optim(par = sspar, fn = loglikem, method = "BFGS", control = c(maxit = 10000), hessian = TRUE)
  if(opt2$convergence!=0) sh5=sh5+1
  if(opt2$par[1] <=0) sh6=sh6 +1
  if(opt2$hessian<=0) sh7=sh7 +1
  if(opt2$convergence==0 & opt2$par[1]>0 & opt2$hessian >0){
    thetamlm[j] = opt2$par[1]
     thetaLm[j] = max(0,thetamlm[j]-1.96/sqrt(opt2$hessian))
    if(thetamlm[j]-1.96/sqrt(opt2$hessian)<0) sh8 = sh8+1
    thetaUm[j] =thetamlm[j]+1.96/sqrt(opt2$hessian)
    
  ########################################################################
  #                 Posterior with records time                        #
  ######################################################################
  passtheta = function(thetao,thetan){
    A=(thetan/thetao)^(2*m+a1-1)*((thetao+1)/(thetan+1))^(2*m)*
      exp(-(b1+sum(record))*(thetan-thetao))*
      prod((2+thetan+record)/(2+thetao+record))*
      exp(sum((time-1)*log(xi(record,thetan)/xi(record,thetao))))
    A
  }
  thetamcmc=c()
  q = 0
  while (q != 1) {
    thetamc =  c()
    thetamc[1] = thetaml[j]
    N = 11001
    M = 1001
    i = 2
    while(i <= N){
      k = i-1 
      thetaold = thetamc[k]
      thetanew = rtruncnorm(1, a = 0, b = Inf, mean = thetaold, sd = 1)
      
      Ptheta = min(1,passtheta(thetaold,thetanew)*
                     dtruncnorm(thetaold, a = 0, b = Inf, mean = thetanew, sd = 1)/
                     dtruncnorm(thetanew, a = 0, b = Inf, mean = thetaold, sd = 1))
      u = runif(1)
      if(u<=Ptheta){thetamc[i] = thetanew}
      if(u> Ptheta){thetamc[i] = thetaold}
      i = i+1
    }
    thetamcmc = thetamc[(M+1):N]

    ######################################################################
    #               Geweke & Raftery & heidel                             #
    ######################################################################
    if(abs(as.numeric(geweke.diag(mcmc(thetamcmc))[1]))<1.96 & 
       raftery.diag(mcmc(thetamcmc),r = 0.01)$resmatrix[4]<5 &
       heidel.diag(mcmc(thetamcmc), eps = 0.1, pvalue = 0.05)[3]>0.05 ){q = 1}
  }
 # plot(thetamcmc, type = "l")
  length(thetamcmc)
  thetaS[j]=mean(thetamcmc)
  thetaLc1[j]=-1/c1*log(mean(exp(-c1*thetamcmc)))
  thetaLc2[j]=-1/c2*log(mean(exp(-c2*thetamcmc)))
  thetaHPDL[j]= HPDinterval(mcmc(thetamcmc))[1]
  thetaHPDU[j]= HPDinterval(mcmc(thetamcmc))[2]
  
  
  ########################################################################
  #                 Posterior without   time                        #
  ######################################################################
  passtheta = function(thetao,thetan){
    A=(thetan/thetao)^(2*m+a1-1)*((thetao+1)/(thetan+1))^(2*m)*
      exp(-(b1+sum(record))*(thetan-thetao))*
      prod((2+thetan+record)/(2+thetao+record))*
      prod((1-xi(recordm,thetao))/(1-xi(recordm,thetan)))
    A
  }
   thetamcmc2=c() 
  q = 0
  while (q != 1) {
    thetamc2 =  c()
    thetamc2[1] = thetamlm[j]
    N = 11001
    M = 1001
    i = 2
    while(i <= N){
      k = i-1 
      thetaold = thetamc2[k]
      thetanew = rtruncnorm(1, a = 0, b = Inf, mean = thetaold, sd = 1)
      
      Ptheta = min(1,passtheta(thetaold,thetanew)*
                     dtruncnorm(thetaold, a = 0, b = Inf, mean = thetanew, sd = 1)/
                     dtruncnorm(thetanew, a = 0, b = Inf, mean = thetaold, sd = 1))
      u = runif(1)
      if(u<=Ptheta){thetamc2[i] = thetanew}
      if(u> Ptheta){thetamc2[i] = thetaold}
      i = i+1
    }
    thetamcmc2 = thetamc2[(M+1):N]
    
    ######################################################################
    #               Geweke & Raftery & Heidel                            #
    ######################################################################
    if(abs(as.numeric(geweke.diag(mcmc(thetamcmc2))[1]))<1.96 & 
       raftery.diag(mcmc(thetamcmc2),r = 0.01)$resmatrix[4]<5 &
       heidel.diag(mcmc(thetamcmc2), eps = 0.1, pvalue = 0.05)[3]>0.05 ){q = 1}
  }
 # plot(thetamcmc2, type = "l")
  length(thetamcmc2)
  thetaSm[j]=mean(thetamcmc2)
  thetaLc1m[j]=-1/c1*log(mean(exp(-c1*thetamcmc2)))
  thetaLc2m[j]=-1/c2*log(mean(exp(-c2*thetamcmc2)))
  thetaHPDLm[j] = HPDinterval(mcmc(thetamcmc2))[1]
  thetaHPDUm[j] = HPDinterval(mcmc(thetamcmc2))[2]
  
  
  
  ######################################################################
  #        one-sequence prediction by SEL with records time            #
  ######################################################################
  fsel = function(rs,theta){
 A = exp(-rs*theta)*((theta+2)*(theta*rs+1)+rs*(theta*rs+2)+2/theta)/(1-xi(rm,theta))/(1 + theta)^2
    return(-A)
  }
  gsel = c()
  for (i in 1:10000) {
      A = fsel(rm,thetamcmc[i])
      B = fsel(0, thetamcmc[i])
    gsel[i] = A-B
  }
  ggsel[j] = mean(gsel)
  ######################################################################
  #    one-sequence prediction by LINEX for c1>0  with records time    #
  ######################################################################
  flinexc1 = function(rs,theta){
    A = theta^2/(1 + theta)^2*exp(-(theta+c1)*rs)*(1+(theta+rs+2)*(theta+c1))/(1-xi(rm,theta))/(theta+c1)^2
        return(-A)
  }
  glinexc1 = c()
  for (i in 1:10000) {
    A = flinexc1(rm,thetamcmc[i])
    B = flinexc1(0, thetamcmc[i])
    glinexc1[i] =A-B
  }
  gglinexc1[j] = -1/c1*log(mean(glinexc1))
  ######################################################################
  #  one-sequence prediction by LINEX for c2<0 with records time       #
  ######################################################################
  flinexc2 = function(rs,theta){
    A = theta^2/(1 + theta)^2*exp(-(theta+c2)*rs)*(1+(theta+rs+2)*(theta+c2))/(1-xi(rm,theta))/(theta+c2)^2
    return(-A)
  }
  glinexc2 = c()
  for (i in 1:10000) {
    A = flinexc2(rm,thetamcmc[i])
    B = flinexc2(0, thetamcmc[i])
    glinexc2[i] =A-B
  }
  gglinexc2 [j] = -1/c2*log(mean(glinexc2)) 
######################################################################
  #        one-sequence prediction by SEL without records time            #
  ######################################################################

gselm = c()
for (i in 1:10000) {
  A = fsel(rm,thetamcmc2[i])
  B = fsel(0, thetamcmc2[i])
  gselm[i] = A-B
}
  ggselm[j] = mean(gselm)
  ######################################################################
  #    one-sequence prediction by LINEX for c1>0  without records time    #
  ######################################################################

  glinexc1m = c()
  for (i in 1:10000) {
    A = flinexc1(rm,thetamcmc2[i])
    B = flinexc1(0, thetamcmc2[i])
    glinexc1m[i] =A-B
  }
  gglinexc1m[j] = -1/c1*log(mean(glinexc1m))
  ######################################################################
  #  one-sequence prediction by LINEX for c2<0 without records time       #
  ######################################################################

  glinexc2m = c()
  for (i in 1:10000) {
    A = flinexc2(rm,thetamcmc2[i])
    B = flinexc2(0, thetamcmc2[i])
    glinexc2m[i] =A-B
  }
  gglinexc2m [j] = -1/c2*log(mean(glinexc2m)) 
 ########################################################################
 #            Interval prediction 
#######################################################################
fpred = function(rs, theta){
  A = (1-xi(rs,theta))/(1-xi(rm,theta))
  return(A)
}
  fhat=function(rs){
    d = 0 
    for (i in 1:10000) {
      fhatt=function(rs){
        A=fpred(rs,thetamcmc[i])
        A
      }
       d=d +fhatt(rs) 
    }
    d/10000
  }

  
 modelU = function(u){
   v = numeric(1)
   v[1] = fhat(u[1]) -0.975
 } 
  t1 = nleqslv(rm, modelU,control = list(maxit = 5000))
 modelL = function(u){
   v = numeric(1)
   v[1] =  fhat(u[1]) -0.025
 } 
  t2 = nleqslv(rm, modelL,control = list(maxit = 5000))
 
 if(abs(t1$fvec)>=0.05 | t1$x >= rm | abs(t2$fvec)>=0.05 | t2$x >= rm) sh9 = sh9 + 1
 if(abs(t1$fvec)< 0.05 & t1$x <  rm & abs(t2$fvec)< 0.05 & t2$x <  rm){
predU[j]=t1$x;
predL[j]=t2$x;
 ########################################################################
 #            Interval prediction without times
#######################################################################
  fhatm = function(rs){
    d = 0 
    for (i in 1:10000) {
      fhattm=function(rs){
        A=fpred(rs,thetamcmc2[i])
        A
              }
        
        d=d +fhattm(rs) 
    }
    d/10000
  }
  
 modelUm = function(u){
   v = numeric(1)
   v[1] = fhatm(u[1]) -0.975
 } 
  t1m = nleqslv(rm, modelUm,control = list(maxit = 5000))
 modelLm = function(u){
   v = numeric(1)
   v[1] = fhatm(u[1]) -0.025
 } 
 t2m = nleqslv(rm,modelLm,control = list(maxit = 5000))
 
 if(abs(t1m$fvec)>=0.05 | t1m$x >= rm | abs(t2m$fvec)>=0.05 | t2m$x >= rm) sh10 = sh10 + 1
 if(abs(t1m$fvec)< 0.05 & t1m$x <  rm & abs(t2m$fvec)< 0.05 & t2m$x <  rm){
predUm[j]=t1m$x;
predLm[j]=t2m$x;
  
  ####################################################################
  
  j = j+1
  }
  } 
  }
 }
}

LOSS1 = function(t,s) exp(c1*(t-s))-c1*(t-s)-1
LOSS2 = function(t,s) exp(c2*(t-s))-c2*(t-s)-1

table1=matrix(0,ncol=4,nrow=8)
rownames(table1)=c(
"MLE          Times","MLE          Without Time", 
"SEL          Times","SEL          Without Time", 
"LINEX;c=+0.5 Times","LINEX;c=+0.5 Without Time", 
"LINEX;c=-0.5 Times","LINEX;c=-0.5 Without Time")
colnames(table1)=c("    BIAS", "       ER_S","    ER_L;c=0.5","   ER_L;c=-0.5")
############################################   BIAS
table1[1,1]=mean(thetaml  -theta)
table1[2,1]=mean(thetamlm -theta)
table1[3,1]=mean(thetaS   -theta)
table1[4,1]=mean(thetaSm  -theta)
table1[5,1]=mean(thetaLc1 -theta)
table1[6,1]=mean(thetaLc1m-theta)
table1[7,1]=mean(thetaLc2 -theta)
table1[8,1]=mean(thetaLc2m-theta)
############################################   RISK_SEL
table1[1,2]=mean((thetaml  -theta)^2)
table1[2,2]=mean((thetamlm -theta)^2)
table1[3,2]=mean((thetaS   -theta)^2)
table1[4,2]=mean((thetaSm  -theta)^2)
table1[5,2]=mean((thetaLc1 -theta)^2)
table1[6,2]=mean((thetaLc1m-theta)^2)
table1[7,2]=mean((thetaLc2 -theta)^2)
table1[8,2]=mean((thetaLc2m-theta)^2)
############################################   RISK_LINEX_c1
table1[1,3]=mean(LOSS1(thetaml,theta))
table1[2,3]=mean(LOSS1(thetamlm,theta))
table1[3,3]=mean(LOSS1(thetaS,theta))
table1[4,3]=mean(LOSS1(thetaSm,theta))
table1[5,3]=mean(LOSS1(thetaLc1,theta))
table1[6,3]=mean(LOSS1(thetaLc1m,theta))
table1[7,3]=mean(LOSS1(thetaLc2,theta))
table1[8,3]=mean(LOSS1(thetaLc2m,theta))
############################################   RISK_LINEX_c2
table1[1,4]=mean(LOSS2(thetaml,theta))
table1[2,4]=mean(LOSS2(thetamlm,theta))
table1[3,4]=mean(LOSS2(thetaS,theta))
table1[4,4]=mean(LOSS2(thetaSm,theta))
table1[5,4]=mean(LOSS2(thetaLc1,theta))
table1[6,4]=mean(LOSS2(thetaLc1m,theta))
table1[7,4]=mean(LOSS2(thetaLc2,theta))
table1[8,4]=mean(LOSS2(thetaLc2m,theta))



table2=matrix(0,ncol=2,nrow=6)
rownames(table2)=c(
"Confidence    Times","Confidence   Without Time", 
"Bayes         Times","Bayes        Without Time",
"Prediction    Times","Prediction   Without Time")
colnames(table2)=c("    AW", "       CP")

############################################   AW
table2[1,1]=mean(thetaU -thetaL )
table2[2,1]=mean(thetaUm-thetaLm)
table2[3,1]=mean(thetaHPDU -thetaHPDL )
table2[4,1]=mean(thetaHPDUm-thetaHPDLm)
table2[5,1]=mean(predU -predL )
table2[6,1]=mean(predUm-predLm)

############################################   CP
table2[1,2]=mean((thetaL<theta)*(theta<thetaU))
table2[2,2]=mean((thetaLm<theta)*(theta<thetaUm))
table2[3,2]=mean((thetaHPDL<theta)*(theta<thetaHPDU))
table2[4,2]=mean((thetaHPDLm<theta)*(theta<thetaHPDUm))
table2[5,2]=mean((predL<rs)*(rs<predU))
table2[6,2]=mean((predLm<rs)*(rs<predUm))


############################################################# Prediction

table3=matrix(0,ncol=4,nrow=6)
rownames(table3)=c(
"SEL          Times","SEL          Without Time", 
"LINEX;c=+0.5 Times","LINEX;c=+0.5 Without Time", 
"LINEX;c=-0.5 Times","LINEX;c=-0.5 Without Time")
colnames(table3)=c("    BIAS", "       ER_S","    ER_L;c=0.5","   ER_L;c=-0.5")
############################################   BIAS
table3[1,1]=mean(ggsel     -rs)
table3[2,1]=mean(ggselm    -rs)
table3[3,1]=mean(gglinexc1 -rs)
table3[4,1]=mean(gglinexc1m-rs)
table3[5,1]=mean(gglinexc2 -rs)
table3[6,1]=mean(gglinexc2m-rs)
############################################   RISK_SEL
table3[1,2]=mean((ggsel     -rs)^2)
table3[2,2]=mean((ggselm    -rs)^2)
table3[3,2]=mean((gglinexc1 -rs)^2)
table3[4,2]=mean((gglinexc1m-rs)^2)
table3[5,2]=mean((gglinexc2 -rs)^2)
table3[6,2]=mean((gglinexc2m-rs)^2)
############################################   RISK_LINEX_c1
table3[1,3]=mean(LOSS1(ggsel,rs))
table3[2,3]=mean(LOSS1(ggselm,rs))
table3[3,3]=mean(LOSS1(gglinexc1,rs))
table3[4,3]=mean(LOSS1(gglinexc1m,rs))
table3[5,3]=mean(LOSS1(gglinexc2,rs))
table3[6,3]=mean(LOSS1(gglinexc2m,rs))
############################################   RISK_LINEX_c2
table3[1,4]=mean(LOSS2(ggsel,rs))
table3[2,4]=mean(LOSS2(ggselm,rs))
table3[3,4]=mean(LOSS2(gglinexc1,rs))
table3[4,4]=mean(LOSS2(gglinexc1m,rs))
table3[5,4]=mean(LOSS2(gglinexc2,rs))
table3[6,4]=mean(LOSS2(gglinexc2m,rs))

################################################################

length(thetaml);length(thetaL   );length( thetaU   );length(thetamlm  );length( thetaLm  );length( thetaUm);
length(thetaS );length(thetaLc1 );length(thetaLc2  );length(thetaHPDL );length(thetaHPDU )
length(thetaSm);length(thetaLc1m);length(thetaLc2m );length(thetaHPDLm);length(thetaHPDUm);
length(predU  );length( predL   );length(gglinexc2 );length(gglinexc1 );length( ggsel    );
length(predUm );length( predLm  );length(gglinexc2m);length(gglinexc1m);length( ggselm   )


sum(thetaml<0);sum(thetaL   <0);sum( thetaU   <0);sum(thetamlm  <0);sum( thetaLm  <0);sum( thetaUm<0);
sum(thetaS <0);sum(thetaLc1 <0);sum(thetaLc2  <0);sum(thetaHPDL <0);sum(thetaHPDU <0)
sum(thetaSm<0);sum(thetaLc1m<0);sum(thetaLc2m <0);sum(thetaHPDLm<0);sum(thetaHPDUm<0);
sum(predU  <0);sum( predL   <0);sum(gglinexc2 <0);sum(gglinexc1 <0);sum( ggsel    <0);
sum(predUm <0);sum( predLm  <0);sum(gglinexc2m<0);sum(gglinexc1m<0);sum( ggselm   <0)

sum(thetaml==0);sum(thetaL   ==0);sum( thetaU   ==0);sum(thetamlm  ==0);sum( thetaLm  ==0);sum( thetaUm==0);
sum(thetaS ==0);sum(thetaLc1 ==0);sum(thetaLc2  ==0);sum(thetaHPDL ==0);sum(thetaHPDU ==0)
sum(thetaSm==0);sum(thetaLc1m==0);sum(thetaLc2m ==0);sum(thetaHPDLm==0);sum(thetaHPDUm==0);
sum(predU  ==0);sum( predL   ==0);sum(gglinexc2 ==0);sum(gglinexc1 ==0);sum( ggsel    ==0);
sum(predUm ==0);sum( predLm  ==0);sum(gglinexc2m==0);sum(gglinexc1m==0);sum( ggselm   ==0)

round(table1, 4)   ####### RISKS AND BIAS FOR ESTIMATORS
round(table3, 6)   ####### RISKS AND BIAS FOR PREDICTORS
round(table2, 5)   ####### AWs AND CPs



shM=c(sh1,sh2,sh3,sh4,sh5,sh6,sh7,sh8,sh9,sh10)
shM;sh
j; m
theta


