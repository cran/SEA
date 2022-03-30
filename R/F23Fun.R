F23Fun<-function(df,model,m_nf){

  data<-sapply(df,as.character)

  dF23<-data[-1,which(data[1,]=="F23")];F23<-as.numeric(dF23[which(is.na(as.numeric(dF23))==FALSE)]);df<-as.data.frame(F23)



##################################################
F3colname <- c("Model","Log_Max_likelihood_Value","AIC","mean[1]","mean[2]","mean[3]","mean[4]","mean[5]","mean[6]","mean[7]","mean[8]","mean[9]",
               "Var[1]","Var[2]","Var[3]","Var[4]","Var[5]","Var[6]","Var[7]","Var[8]","Var[9]","Proportion[1]","Proportion[2]","Proportion[3]","Proportion[4]","Proportion[5]","Proportion[6]",
               "Proportion[7]","Proportion[8]","Proportion[9]","m","da(d)","db","ha(h)","hb","i","jab","jba","l","Major-Gene Var","Heritability(Major-Gene)(%)",
               "U1 square","P(U1 square)","U2 square","P(U2 square)","U3 square","P(U3 square)","nW square","P(nW square)","Dn","P(Dn)")

F3ModelFun<-list(NA)
###################define each model function##################
############################################0MG model#########################################
F3ModelFun[[1]] <-function(K1,logL,df,m_nf){
  data <- as.matrix(as.numeric(df[,1]))
  m_sam <- dim(data)[1]; mean0 <- mean(data); sigma0 <- var(data)
  #########0MG-Model#########
  m<-mean0
  sigma<-sigma0
  mix_pi <- 1
  lp<-matrix(0,m_sam,1)
  for(j in 1:m_sam){lp[j]<-log(dnorm(data[j],m,sqrt(sigma)))}
  abc<-sum(lp)
  AIC<--2*abc+2*2
  #########hypothesis testing##############
  data<-sort(data)
  w1<-1/(12*m_sam)
  bmw <- matrix(0,m_sam,1)
  gg <- (data - m)/sqrt(as.vector(sigma))
  bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
  bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(bmw)))[1]
  if(nn < m_sam){bmw <- bmw+runif(m_sam)/1e4}
  ##################################################
  dd<-c((sum(bmw)),(sum(bmw^2)),sum((bmw-0.5)^2))
  w<-w1+sum((bmw - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u<- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D<-as.numeric(ks.test(bmw,"punif")[2])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(w),D))
  D<-as.numeric(ks.test(bmw,"punif")[[1]][1])

  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  output <- data.frame("0MG",round(abc,4),round(AIC,4),round(m,4)," "," "," "," "," "," "," "," ",round(sigma,4)," "," "," "," "," "," "," "," ",round(mix_pi,4)," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(w,4),tt[4],round(D,4),tt[5])
 output<-as.matrix(output)
  OUTPUT<-list(output)
  return(OUTPUT)
}
############################################1MG-AD model#########################################
F3ModelFun[[2]]  <-function(K1,logL,df,m_nf){
  data <- as.matrix(as.numeric(df[,1]))
  ##############procedure start####################
  m_sam <- dim(data)[1];  d2<- 3; m_esp <- 0.0001
  mean0 <- mean(data);sigma0 <- var(data)
  ##########1MG-AD(A-1)#####################
  a1 <- sqrt(sigma0/m_sam)
  m <- as.matrix(c((mean0+2*a1),mean0,(mean0-2*a1)))
  mi <- as.matrix(c(0.25,0.5,0.25))
  sigma <- matrix((sigma0/2),3,1)
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,3,m_sam); swx <- matrix(0,3,1); s0 <- matrix(0,2,1)
  while(stopa > m_esp && iteration <= 1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    m <- sumwx/n0
    ######first order genetic parameter process##########
    hh1 <- matrix(c(1,1,1,1,0,-1,0,0.5,0),3,3)
    b_line1 <- m; B1 <- solve(hh1,b_line1)
    aa3<-(0.5*B1[2]^2+0.25*B1[3]^2)/m_nf
    ##########obtain variance##########
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0[1]<-swx[1]+swx[3]
    s0[2]<-n0[1]+n0[3]
    aaa0<-sigma[1]
    n_iter<-0
    aa5<-1000
    while(aa5>0.0001){
      n_iter<-n_iter+1.0
      aa4<-sigma[1]/(sigma[1]+aa3)
      sigma[1]<-(s0[1]+aa4^2*swx[2])/(s0[2]+aa4*n0[2])
      aa5<-abs(sigma[1]-aaa0)
      aaa0<-sigma[1]
      if(n_iter>20)break
    }
    sigma[2]<-sigma[1]+aa3
    sigma[3]<-sigma[1]
    if(sum(sigma < 1e-30)>=1){break}
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*5
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0
  #########hypothesis testing##############
  data <- sort(data)
  w1<-1/(12*m_sam)
  bmw <- matrix(0,m_sam,1); bmwsl <- matrix(0,m_sam,d2)
  for(i in 1:d2){
    gg <- (data - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2)))[1]
  if(nn < m_sam){P2 <- P2+runif(m_sam)/1e4}
  #########################################################
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P2 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[2])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),D))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])

  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  output <- data.frame("1MG-AD",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," "," ",round(t(sigma),4)," "," "," "," "," "," ",round(t(mix_pi),4)," "," "," "," "," "," ",round(B1[1],4),round(B1[2],4)," ",round(B1[3],4)," "," "," "," "," ",round(jj,4),round(ll*100,4),round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
 output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################1MG-A model#########################################
F3ModelFun[[3]]  <-function(K1,logL,df,m_nf)
{
  data <- as.matrix(as.numeric(as.numeric(df[,1])))
  ##############procedure start####################
  m_sam <- dim(data)[1]; d2 <- 3; m_esp <- 0.0001
  mean0 <- mean(data);sigma0 <- var(data)
  ##########1MG-A(A-2)#####################
  a1 <- sqrt(sigma0)
  m <- as.matrix(c((mean0+2*a1),mean0,(mean0-2*a1)))
  mi <- as.matrix(c(0.25,0.5,0.25))
  sigma <- matrix((sigma0/2),3,1)
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,3,m_sam); swx <- matrix(0,3,1);s0 <- matrix(0,2,1)
  while(stopa > m_esp && iteration <= 1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    aa1<-sigma[1]/n0[1]+4*sigma[2]/n0[2]+sigma[3]/n0[3]
    aa2<-sumwx[1]/n0[1]-2.0*sumwx[2]/n0[2]+sumwx[3]/n0[3]
    aa3<-aa2/aa1
    m[1]<-(sumwx[1]-sigma[1]*aa3)/n0[1]
    m[2]<-(sumwx[2]+sigma[1]*2*aa3)/n0[2]
    m[3]<-(sumwx[3]-sigma[1]*aa3)/n0[3]
    ######first order genetic parameter process##########
    hh1 <- matrix(c(1,1,1,1,0,-1),3,2)
    B1<-solve(crossprod(hh1,hh1))%*%crossprod(hh1,m)
    aa3<-0.5*B1[2]^2/m_nf
    ##########obtain variance##########
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0[1]<-swx[1]+swx[3]
    s0[2]<-n0[1]+n0[3]
    aaa0<-sigma[1]
    n_iter<-0
    aa5<-1000
    while(aa5>0.0001){
      n_iter<-n_iter+1.0
      aa4<-sigma[1]/(sigma[1]+aa3)
      sigma[1]<-(s0[1]+aa4^2*swx[2])/(s0[2]+aa4*n0[2])
      aa5<-abs(sigma[1]-aaa0)
      aaa0<-sigma[1]
      if(n_iter>20)break
    }
    sigma[2]<-sigma[1]+aa3
    sigma[3]<-sigma[1]
    if(sum(sigma < 1e-30)>=1){break}
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*3
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0
  #########hypothesis testing##############
  data <- sort(data)
  w1<-1/(12*m_sam)
  bmw <- matrix(0,m_sam,1); bmwsl <- matrix(0,m_sam,d2)
  for(i in 1:d2){
    gg <- (data - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2)))[1]
  if(nn < m_sam){P2 <- P2+runif(m_sam)/1e4}
  #####################################################################
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P2 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[2])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),D))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])

  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  output <- data.frame("1MG-A",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," "," ",round(t(sigma),4)," "," "," "," "," "," ",round(t(mix_pi),4)," "," "," "," "," "," ",round(B1[1],4),round(B1[2],4)," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4),round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
 output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################1MG-EAD model#########################################
F3ModelFun[[4]]  <-function(K1,logL,df,m_nf){
  data <- as.matrix(as.numeric(as.numeric(df[,1])))
  ##############procedure start####################
  m_sam <- dim(data)[1]; d2 <- 3; m_esp <- 0.0001
  mean0 <- mean(data);sigma0 <- var(data)
  ##########1MG-EAD(A-3)#####################
  a1 <- sqrt(sigma0/m_sam)
  m <- as.matrix(c((mean0+2*a1),mean0,(mean0-2*a1)))
  mi <- as.matrix(c(0.25,0.5,0.25))
  sigma <- matrix((sigma0/2),3,1)
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi)))
    ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,3,m_sam); swx <- matrix(0,3,1);s0 <- matrix(0,2,1)
  while(stopa > m_esp && iteration <= 1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    aa1<-3*sumwx[1]/n0[1]-4*sumwx[2]/n0[2]+sumwx[3]/n0[3]
    aa2<-9*sigma[1]/n0[1]+16*sigma[2]/n0[2]+sigma[3]/n0[3]
    aa3<-aa1/aa2
    m[1]<-(sumwx[1]-sigma[1]*3*aa3)/n0[1]
    m[2]<-(sumwx[2]+sigma[2]*4*aa3)/n0[2]
    m[3]<-(sumwx[3]-sigma[3]*aa3)/n0[3]
    ######first order genetic parameter process##########
    hh1 <- matrix(c(1,1,1,1,0.5,-1),3,2)
    B1<-solve(crossprod(hh1,hh1))%*%crossprod(hh1,m)
    aa3<-0.75*B1[2]^2/m_nf
    ##########obtain variance##########
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0[1]<-swx[1]+swx[3]
    s0[2]<-n0[1]+n0[3]
    aaa0<-sigma[1]
    n_iter<-0
    aa5<-1000
    while(aa5>0.0001){
      n_iter<-n_iter+1.0
      aa4<-sigma[1]/(sigma[1]+aa3)
      sigma[1]<-(s0[1]+aa4^2*swx[2])/(s0[2]+aa4*n0[2])
      aa5<-abs(sigma[1]-aaa0)
      aaa0<-sigma[1]
      if(n_iter>20)break
    }
    sigma[2]<-sigma[1]+aa3
    sigma[3]<-sigma[1]
    if(sum(sigma < 1e-30)>=1){break}
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*4
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0
  #########hypothesis testing##############
  data <- sort(data);
  w1<-1/(12*m_sam)
  bmw <- matrix(0,m_sam,1); bmwsl <- matrix(0,m_sam,d2)
  for(i in 1:d2){
    gg <- (data - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2)))[1]
  if(nn < m_sam){P2 <- P2+runif(m_sam)/1e4}
  #####################################################################
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P2 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[2])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),D))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])

  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  output <- data.frame("1MG-EAD",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," "," ",round(t(sigma),4)," "," "," "," "," "," ",round(t(mix_pi),4)," "," "," "," "," "," ",round(B1[1],4),round(B1[2],4)," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4),round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
 output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################1MG-NCD model#########################################
F3ModelFun[[5]]  <-function(K1,logL,df,m_nf){
  data <- as.matrix(as.numeric(as.numeric(df[,1])))
  ##############procedure start####################
  m_sam <- dim(data)[1]; d2 <- 3; m_esp <- 0.0001
  mean0 <- mean(data);sigma0 <- var(data)
  ##########1MG-NCD(A-4)#####################
  a1 <- sqrt(sigma0/m_sam)
  m <- as.matrix(c((mean0+2*a1),mean0,(mean0-2*a1)))
  mi <- as.matrix(c(0.25,0.5,0.25))
  sigma <- matrix((sigma0/2),3,1)
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,3,m_sam); swx <- matrix(0,3,1); s0 <- matrix(0,2,1)
  while(stopa > m_esp && iteration <= 1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    aa1<-sumwx[1]/n0[1]-4.0*sumwx[2]/n0[2]+3.0*sumwx[3]/n0[3]
    aa2<-sigma[1]/n0[1]+16.0*sigma[2]/n0[2]+9.0*sigma[3]/n0[3]
    aa3<-aa1/aa2
    m[1]<-(sumwx[1]-sigma[1]*aa3)/n0[1]
    m[2]<-(sumwx[2]+sigma[2]*aa3*4.0)/n0[2]
    m[3]<-(sumwx[3]-sigma[3]*aa3*3.0)/n0[3]
    ######first order genetic parameter process##########
    hh1 <- matrix(c(1,1,1,1,-0.5,-1),3,2)
    B1<-solve(crossprod(hh1,hh1))%*%crossprod(hh1,m)
    aa3<-0.75*B1[2]^2/m_nf
    ##########obtain variance##########
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0[1]<-swx[1]+swx[3]
    s0[2]<-n0[1]+n0[3]
    aaa0<-sigma[1]
    n_iter<-0
    aa5<-1000
    while(aa5>0.0001){
      n_iter<-n_iter+1.0
      aa4<-sigma[1]/(sigma[1]+aa3)
      sigma[1]<-(s0[1]+aa4^2*swx[2])/(s0[2]+aa4*n0[2])
      aa5<-abs(sigma[1]-aaa0)
      aaa0<-sigma[1]
      if(n_iter>20)break
    }
    sigma[2]<-sigma[1]+aa3
    sigma[3]<-sigma[1]
    if(sum(sigma < 1e-30)>=1){break}
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*4
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0
  #########hypothesis testing##############
  data <- sort(data)
  w1<-1/(12*m_sam)
  bmw <- matrix(0,m_sam,1); bmwsl <- matrix(0,m_sam,d2)
  for(i in 1:d2){
    gg <- (data - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2)))[1]
  if(nn < m_sam){P2 <- P2+runif(m_sam)/1e4}
  #####################################################################
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P2 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[2])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),D))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])

  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  output <- data.frame("1MG-NCD",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," "," ",round(t(sigma),4)," "," "," "," "," "," ",round(t(mix_pi),4)," "," "," "," "," "," ",round(B1[1],4),round(B1[2],4)," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4),round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
 output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################2MG-ADI model#########################################
F3ModelFun[[6]]  <-function(K1,logL,df,m_nf){
  data <- as.matrix(as.numeric(as.numeric(df[,1])))
  ##############procedure start####################
  m_sam <- dim(data)[1]; d4 <- 9; m_esp <- 0.0001
  mean0 <- mean(data);sigma0 <- var(data)
  #############2MG-ADI model(B-1)###########################
  a1 <- sqrt(sigma0)
  sigma <- matrix((sigma0/2),9,1)
  m <- as.matrix(c((mean0+3.2*a1),(mean0+2.4*a1),(mean0+1.6*a1),(mean0+0.8*a1),mean0,
                   (mean0-0.8*a1),(mean0-1.6*a1),(mean0-2.4*a1),(mean0-3.2*a1)))
  mi <- as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,9,m_sam); swx <- matrix(0,9,1);g_a<-matrix(0,5,1);s0 <- matrix(0,2,1)
  while(stopa > m_esp && iteration <= 1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d4) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    m<-sumwx/n0
    hh5 <- matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,-1,-1,-1,1,0,-1,1,0,-1,1,0,-1,0,0,0,0.5,0.5,0.5,0,0,0,
                    0,0.5,0,0,0.5,0,0,0.5,0,1,0,-1,0,0,0,-1,0,1,0,0.5,0,0,0,0,0,-0.5,0,0,0,0,0.5,0,-0.5,0,
                    0,0,0,0,0,0,0.25,0,0,0,0),9,9)
    b_line5 <- m; B5 <- solve(hh5,b_line5)
    g_a[1]<-(0.5*(B5[3]+B5[6])^2+0.25*(B5[5]+B5[7])^2)/m_nf
    g_a[2]<-(0.5*(B5[2]+B5[6])^2+0.25*(B5[4]+B5[8])^2)/m_nf
    g_a[3]<-(0.5*(B5[2]-B5[6])^2+0.25*(B5[4]-B5[8])^2)/m_nf
    g_a[4]<-(0.5*(B5[3]-B5[6])^2+0.25*(B5[5]-B5[7])^2)/m_nf
    g_a[5]<-0.25*(B5[2]^2+B5[3]^2+B5[6]^2+(B5[2]+B5[7])^2+(B5[3]+B5[8])^2+(B5[4]+0.5*B5[9])^2+(B5[5]+0.5*B5[9])^2+0.25*B5[9]^2)/m_nf
    for(i in 1:d4) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0[1]<-swx[1]+swx[3]+swx[7]+swx[9]
    s0[2]<-n0[1]+n0[3]+n0[7]+n0[9]
    aaa0<-sigma[1]
    n_iter<-0
    aa6<-1000
    while(aa6>0.0001){
      n_iter<-n_iter+1
      aa1<-sigma[1]/(sigma[1]+g_a[1]);aa2<-sigma[1]/(sigma[1]+g_a[2])
      aa3<-sigma[1]/(sigma[1]+g_a[3]);aa4<-sigma[1]/(sigma[1]+g_a[4])
      aa5<-sigma[1]/(sigma[1]+g_a[5])
      sigma[1]<-(s0[1]+aa1^2*swx[2]+aa2^2*swx[4]+aa5^2*swx[5]+aa3^2*swx[6]+aa4^2*swx[8])/(s0[2]+aa1*n0[2]+aa2*n0[4]+aa5*n0[5]+aa3*n0[6]+aa4*n0[8])
      aa6<-abs(sigma[1]-aaa0)
      aaa0<-sigma[1]
      if(n_iter>20)break
    }
    sigma[3]<-sigma[1];sigma[7]<-sigma[1];sigma[9]<-sigma[1]
    sigma[2]<-sigma[1]+g_a[1];sigma[4]<-sigma[1]+g_a[2];sigma[5]<-sigma[1]+g_a[5]
    sigma[6]<-sigma[1]+g_a[3];sigma[8]<-sigma[1]+g_a[4]
    if(sum(sigma < 1e-30)>=1){break}
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
    if(iteration>1000)break
  }
  abc<-L0
  AIC<--2*abc+2*10
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0
  #########hypothesis testing##############
  data <- sort(data)
  w1<-1/(12*m_sam)
  bmw <- matrix(0,m_sam,1); bmwsl <- matrix(0,m_sam,d4)
  for(i in 1:d4){
    gg <- (data - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2)))[1]
  if(nn < m_sam){P2 <- P2+runif(m_sam)/1e4}
  ####################################################################
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P2 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[2])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),D))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])

  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  output <- data.frame("2MG-ADI",round(abc,4),round(AIC,4),round(t(m),4),round(t(sigma),4),round(t(mix_pi),4),round(t(B5),4),round(jj,4),round(ll*100,4),round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
 output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################2MG-AD model#########################################
F3ModelFun[[7]]  <-function(K1,logL,df,m_nf){
  data <- as.matrix(as.numeric(as.numeric(df[,1])))
  ##############procedure start####################
  m_sam <- dim(data)[1]; d4 <- 9; m_esp <- 0.0001
  mean0 <- mean(data);sigma0 <- var(data)
  #########2MG-AD model(B-2)#############
  a1 <- sqrt(sigma0)
  sigma <- matrix((sigma0/2),9,1)
  m <- as.matrix(c((mean0+3.2*a1),(mean0+2.4*a1),(mean0+1.6*a1),(mean0+0.8*a1),mean0,(mean0-0.8*a1),(mean0-1.6*a1),(mean0-2.4*a1),(mean0-3.2*a1)))
  mi <- as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,9,m_sam); swx <- matrix(0,9,1); hh <- matrix(0,4,4); b_line <- matrix(0,4,1)
  g_a<-matrix(0,3,1);s0 <- matrix(0,2,1)
  while(stopa > m_esp && iteration <= 1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d4) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    ############solve the linear equation############
    hh[1,1]<-sigma[1]/n0[1]+sigma[3]/n0[3]+sigma[7]/n0[7]+sigma[9]/n0[9]
    hh[1,2]<-sigma[1]/n0[1]+sigma[7]/n0[7]
    hh[1,3]<-sigma[1]/n0[1]+sigma[3]/n0[3]
    hh[1,4]<-sigma[7]/n0[7]-sigma[9]/n0[9]
    hh[2,2]<-sigma[1]/n0[1]+sigma[2]/n0[2]+sigma[7]/n0[7]+sigma[8]/n0[8]
    hh[2,3]<-sigma[1]/n0[1]
    hh[2,4]<-sigma[7]/n0[7]+2.0*sigma[8]/n0[8]
    hh[3,3]<-sigma[1]/n0[1]+sigma[3]/n0[3]+sigma[4]/n0[4]+sigma[6]/n0[6]
    hh[3,4]<--sigma[4]/n0[4]+sigma[6]/n0[6]
    hh[4,4]<-sigma[4]/n0[4]+4.0*sigma[5]/n0[5]+sigma[6]/n0[6]+sigma[7]/n0[7]+4.0*sigma[8]/n0[8]+sigma[9]/n0[9]
    for(i in 2:4)
    {
      for(j in 1:(i-1))
      {
        hh[i,j]<-hh[j,i]
      }
    }
    #################################################
    b_line[1]<-sumwx[1]/n0[1]-sumwx[3]/n0[3]-sumwx[7]/n0[7]+sumwx[9]/n0[9]
    b_line[2]<-sumwx[1]/n0[1]-sumwx[2]/n0[2]-sumwx[7]/n0[7]+sumwx[8]/n0[8]
    b_line[3]<-sumwx[1]/n0[1]-sumwx[3]/n0[3]-sumwx[4]/n0[4]+sumwx[6]/n0[6]
    b_line[4]<-sumwx[4]/n0[4]-2.0*sumwx[5]/n0[5]+sumwx[6]/n0[6]-sumwx[7]/n0[7]+2.0*sumwx[8]/n0[8]-sumwx[9]/n0[9]
    B6<- solve(hh,b_line)
    #################################################
    m[1]<-(sumwx[1]-sigma[1]*(B6[1]+B6[2]+B6[3]))/n0[1]
    m[2]<-(sumwx[2]+sigma[2]*B6[2])/n0[2]
    m[3]<-(sumwx[3]+sigma[3]*(B6[1]+B6[3]))/n0[3]
    m[4]<-(sumwx[4]+sigma[4]*(B6[3]-B6[4]))/n0[4]
    m[5]<-(sumwx[5]+sigma[5]*2*B6[4])/n0[5]
    m[6]<-(sumwx[6]-sigma[6]*(B6[3]+B6[4]))/n0[6]
    m[7]<-(sumwx[7]+sigma[7]*(B6[1]+B6[2]+B6[4]))/n0[7]
    m[8]<-(sumwx[8]-sigma[8]*(B6[2]+2*B6[4]))/n0[8]
    m[9]<-(sumwx[9]-sigma[9]*(B6[1]-B6[4]))/n0[9]
    #########first order genetic parameter process##########
    hh61 <- matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,-1,-1,-1,1,0,-1,1,0,-1,1,0,-1,
                     0,0,0,0.5,0.5,0.5,0,0,0,0,0.5,0,0,0.5,0,0,0.5,0),9,5)
    B61<-solve(crossprod(hh61,hh61))%*%crossprod(hh61,m)
    g_a[1]<-(0.5*B61[3]^2+0.25*B61[5]^2)/m_nf
    g_a[2]<-(0.5*B61[2]^2+0.25*B61[4]^2)/m_nf
    g_a[3]<-(0.5*(B61[2]^2+B61[3]^2)+0.25*(B61[4]^2+B61[5]^2))/m_nf
    #######################obtain variance#######################################
    for(i in 1:d4) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0[1]<-swx[1]+swx[3]+swx[7]+swx[9]
    s0[2]<-n0[1]+n0[3]+n0[7]+n0[9]
    aaa0<-sigma[1]
    n_iter<-0
    aa6<-1000
    while(aa6>0.0001){
      n_iter<-n_iter+1
      aa1<-sigma[1]/(sigma[1]+g_a[1])
      aa2<-sigma[1]/(sigma[1]+g_a[2])
      aa3<-sigma[1]/(sigma[1]+g_a[3])
      sigma[1]<-(s0[1]+aa1^2*(swx[2]+swx[8])+aa2^2*(swx[4]+swx[6])+aa3^2*swx[5])/(s0[2]+aa1*(n0[2]+n0[8])+aa2*(n0[4]+n0[6])+aa3*n0[5])
      aa6<-abs(sigma[1]-aaa0)
      aaa0<-sigma[1]
      if(n_iter>20)break
    }
    sigma[3]<-sigma[1];sigma[7]<-sigma[1];sigma[9]<-sigma[1]
    sigma[2]<-sigma[1]+g_a[1];sigma[4]<-sigma[1]+g_a[2]
    sigma[5]<-sigma[1]+g_a[3];sigma[6]<-sigma[1]+g_a[2];sigma[8]<-sigma[1]+g_a[1]
    if(sum(sigma < 1e-30)>=1){break}
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*6
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0
  #########hypothesis testing##############
  data <- sort(data)
  w1<-1/(12*m_sam)
  bmw <- matrix(0,m_sam,1); bmwsl <- matrix(0,m_sam,d4)
  for(i in 1:d4){
    gg <- (data - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2)))[1]
  if(nn < m_sam){P2 <- P2+runif(m_sam)/1e4}
  #####################################################################
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P2 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[2])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),D))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])

  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  output <- data.frame("2MG-AD",round(abc,4),round(AIC,4),round(t(m),4),round(t(sigma),4),round(t(mix_pi),4),round(t(B61),4)," "," "," "," ",round(jj,4),round(ll*100,4),round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
 output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################2MG-A model#########################################
F3ModelFun[[8]]  <-function(K1,logL,df,m_nf){
  data <- as.matrix(as.numeric(as.numeric(df[,1])))
  ##############procedure start####################
  m_sam <- dim(data)[1]; d4 <- 9; m_esp <- 0.0001
  mean0 <- mean(data);sigma0 <- var(data)
  ########2MG-A model(B-3)###############
  a1 <- sqrt(sigma0)
  sigma <- matrix((sigma0/2),9,1)
  m <- as.matrix(c((mean0+3.2*a1),(mean0+2.4*a1),(mean0+1.6*a1),(mean0+0.8*a1),mean0,(mean0-0.8*a1),(mean0-1.6*a1),(mean0-2.4*a1),(mean0-3.2*a1)))
  mi <- as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,9,m_sam); swx <- matrix(0,9,1); hh <- matrix(0,6,6); b_line <- matrix(0,6,1)
  g_a<-matrix(0,3,1); s0 <- matrix(0,2,1)
  while(stopa > m_esp && iteration <= 1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d4) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    ############solve the linear equation############
    hh[1,1]<-sigma[1]/n0[1]+sigma[3]/n0[3]+sigma[7]/n0[7]+sigma[9]/n0[9]
    hh[1,2]<-sigma[1]/n0[1]+sigma[7]/n0[7]
    hh[1,3]<--sigma[3]/n0[3]-sigma[7]/n0[7]
    hh[1,4]<-2.0*sigma[3]/n0[3]
    hh[1,5]<--sigma[7]/n0[7]+sigma[9]/n0[9]
    hh[1,6]<-sigma[1]/n0[1]-sigma[7]/n0[7]
    hh[2,2]<-sigma[1]/n0[1]+sigma[2]/n0[2]+sigma[7]/n0[7]+sigma[8]/n0[8]
    hh[2,3]<--sigma[7]/n0[7]
    hh[2,4]<--2.0*sigma[2]/n0[2]
    hh[2,5]<--sigma[7]/n0[7]-2.0*sigma[8]/n0[8]
    hh[2,6]<-hh[1,6]
    hh[3,3]<-sigma[3]/n0[3]+4.0*sigma[5]/n0[5]+sigma[7]/n0[7]
    hh[3,4]<--2.0*sigma[3]/n0[3]
    hh[3,5]<-sigma[7]/n0[7]
    hh[3,6]<-hh[3,5]
    hh[4,4]<-4.0*sigma[2]/n0[2]+4.0*sigma[3]/n0[3]+sigma[4]/n0[4]+sigma[6]/n0[6]
    hh[4,5]<-0
    hh[4,6]<-2.0*sigma[4]/n0[4]
    hh[5,5]<-sigma[7]/n0[7]+4.0*sigma[8]/n0[8]+sigma[9]/n0[9]
    hh[5,6]<-hh[3,5]
    hh[6,6]<-sigma[1]/n0[1]+4.0*sigma[4]/n0[4]+sigma[7]/n0[7]
    for(i in 2:6)
    {
      for(j in 1:(i-1))
      {
        hh[i,j]<-hh[j,i]
      }
    }
    ###############################################################
    b_line[1]<-sumwx[1]/n0[1]-sumwx[3]/n0[3]-sumwx[7]/n0[7]+sumwx[9]/n0[9]
    b_line[2]<-sumwx[1]/n0[1]-sumwx[2]/n0[2]-sumwx[7]/n0[7]+sumwx[8]/n0[8]
    b_line[3]<-sumwx[3]/n0[3]-2*sumwx[5]/n0[5]+sumwx[7]/n0[7]
    b_line[4]<-2*sumwx[2]/n0[2]-2*sumwx[3]/n0[3]-sumwx[4]/n0[4]+sumwx[6]/n0[6]
    b_line[5]<-sumwx[7]/n0[7]-2*sumwx[8]/n0[8]+sumwx[9]/n0[9]
    b_line[6]<-sumwx[1]/n0[1]-2*sumwx[4]/n0[4]+sumwx[7]/n0[7]
    B7<- solve(hh,b_line)
    ############################################################################
    m[1]<-(sumwx[1]-sigma[1]*(B7[1]+B7[2]+B7[6]))/n0[1]
    m[2]<-(sumwx[2]+sigma[2]*(B7[2]-2*B7[4]))/n0[2]
    m[3]<-(sumwx[3]+sigma[3]*(B7[1]-B7[3]+2*B7[4]))/n0[3]
    m[4]<-(sumwx[4]+sigma[4]*(-B7[4]+2*B7[6]))/n0[4]
    m[5]<-(sumwx[5]+sigma[5]*(2*B7[3]))/n0[5]
    m[6]<-(sumwx[6]-sigma[6]*B7[4])/n0[6]
    m[7]<-(sumwx[7]+sigma[7]*(B7[1]+B7[2]-B7[3]-B7[5]-B7[6]))/n0[7]
    m[8]<-(sumwx[8]+sigma[8]*(-B7[2]+2*B7[5]))/n0[8]
    m[9]<-(sumwx[9]-sigma[9]*(B7[1]+B7[5]))/n0[9]
    #########first order genetic parameter process##########
    hh71 <- matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,-1,-1,-1,1,0,-1,1,0,-1,1,0,-1),9,3)
    B71<-solve(crossprod(hh71,hh71))%*%crossprod(hh71,m)
    g_a[1]<-(0.5*B71[3]^2)/m_nf
    g_a[2]<-(0.5*B71[2]^2)/m_nf
    g_a[3]<-g_a[1]+g_a[2]
    #######################obtain variance#######################################
    for(i in 1:d4) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0[1]<-swx[1]+swx[3]+swx[7]+swx[9]
    s0[2]<-n0[1]+n0[3]+n0[7]+n0[9]
    aaa0<-sigma[1]
    n_iter<-0
    aa6<-1000
    while(aa6>0.0001){
      n_iter<-n_iter+1
      aa1<-sigma[1]/(sigma[1]+g_a[1])
      aa2<-sigma[1]/(sigma[1]+g_a[2])
      aa3<-sigma[1]/(sigma[1]+g_a[3])
      sigma[1]<-(s0[1]+aa1^2*(swx[2]+swx[8])+aa2^2*(swx[4]+swx[6])+aa3^2*swx[5])/(s0[2]+aa1*(n0[2]+n0[8])+aa2*(n0[4]+n0[6])+aa3*n0[5])
      aa6<-abs(sigma[1]-aaa0)
      aaa0<-sigma[1]
      if(n_iter>20)break
    }
    sigma[3]<-sigma[1];sigma[7]<-sigma[1];sigma[9]<-sigma[1]
    sigma[2]<-sigma[1]+g_a[1];sigma[4]<-sigma[1]+g_a[2]
    sigma[5]<-sigma[1]+g_a[3];sigma[6]<-sigma[1]+g_a[2];sigma[8]<-sigma[1]+g_a[1]
    if(sum(sigma < 1e-30)>=1){break}
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
    if(iteration>60)break
  }
  abc<-L0
  AIC<--2*abc+2*4
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0
  #########hypothesis testing##############
  data <- sort(data)
  w1<-1/(12*m_sam)
  bmw <- matrix(0,m_sam,1); bmwsl <- matrix(0,m_sam,d4)
  for(i in 1:d4){
    gg <- (data - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2)))[1]
  if(nn < m_sam){P2 <- P2+runif(m_sam)/1e4}
  #####################################################################
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P2 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[2])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),D))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])


  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  output <- data.frame("2MG-A",round(abc,4),round(AIC,4),round(t(m),4),round(t(sigma),4),round(t(mix_pi),4),round(t(B71),4)," "," "," "," "," "," ",round(jj,4),round(ll*100,4),round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
 output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################2MG-EA model#########################################
F3ModelFun[[9]]  <-function(K1,logL,df,m_nf){
  data <- as.matrix(as.numeric(as.numeric(df[,1])))
  ##############procedure start####################
  m_sam <- dim(data)[1]; d5 <- 6; m_esp <- 0.0001
  mean0 <- mean(data);sigma0 <- var(data)
  ##############2MG-EA model(B-4)#############################
  a1 <- sqrt(sigma0)
  sigma <- matrix((sigma0/2),6,1)
  m <- as.matrix(c((mean0+3*a1),(mean0+2.4*a1),mean0,
                   (mean0+1.7*a1),(mean0-0.7*a1),(mean0-3*a1)))
  mi <- as.matrix(c(0.0625,0.25,0.125,0.25,0.25,0.0625))
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,6,m_sam); swx <- matrix(0,6,1); hh <- matrix(0,3,3); b_line <- matrix(0,3,1)
  g_a<-matrix(0,2,1); s0 <- matrix(0,2,1)
  while(stopa > m_esp && iteration <= 1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d5) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    ############solve the linear equation############
    hh[1,1]<-sigma[1]/n0[1]+4.0*sigma[2]/n0[2]+sigma[3]/n0[3]
    hh[1,2]<--2.0*sigma[2]/n0[2]-2.0*sigma[3]/n0[3]
    hh[1,3]<-sigma[3]/n0[3]
    hh[2,2]<-sigma[2]/n0[2]+4.0*sigma[3]/n0[3]+sigma[5]/n0[5]
    hh[2,3]<--2.0*sigma[3]/n0[3]-2.0*sigma[5]/n0[5]
    hh[3,3]<-sigma[3]/n0[3]+4.0*sigma[5]/n0[5]+sigma[6]/n0[6]
    for(i in 2:3)
    {
      for(j in 1:(i-1))
      {
        hh[i,j]<-hh[j,i]
      }
    }
    ###############################################################
    b_line[1]<-sumwx[1]/n0[1]-2.0*sumwx[2]/n0[2]+sumwx[3]/n0[3]
    b_line[2]<-sumwx[2]/n0[2]-2.0*sumwx[3]/n0[3]+sumwx[5]/n0[5]
    b_line[3]<-sumwx[3]/n0[3]-2.0*sumwx[5]/n0[5]+sumwx[6]/n0[6]
    B8<- solve(hh,b_line)
    ############################################################################
    m[1]<-(sumwx[1]-sigma[1]*B8[1])/n0[1]
    m[2]<-(sumwx[2]+sigma[2]*(2.0*B8[1]-B8[2]))/n0[2]
    m[3]<-(sumwx[3]+(sigma[3]/sigma[4])*sumwx[4]+sigma[3]*(-B8[1]+2.0*B8[2]-B8[3]))/(n0[3]+(sigma[3]/sigma[4])*n0[4])
    m[4]<-m[3]
    m[5]<-(sumwx[5]-sigma[5]*(B8[2]-2.0*B8[3]))/n0[5]
    m[6]<-(sumwx[6]-sigma[6]*B8[3])/n0[6]
    #########first order genetic parameter process##########
    hh81 <- matrix(c(1,1,1,1,1,2,1,0,-1,-2),5,2)
    B81<-solve(crossprod(hh81,hh81))%*%crossprod(hh81,m[-4])
    g_a[1]<-0.5*B81[2]^2/m_nf
    g_a[2]<-B81[2]^2/m_nf
    #######################obtain variance#######################################
    for(i in 1:d5) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0[1]<-swx[1]+swx[6]
    s0[2]<-n0[1]+n0[6]
    aaa0<-sigma[1]
    n_iter<-0
    aa6<-1000
    while(aa6>0.0001){
      n_iter<-n_iter+1
      aa1<-sigma[1]/(sigma[1]+g_a[1])
      aa2<-sigma[1]/(sigma[1]+g_a[2])
      sigma[1]<-(s0[1]+aa1^2*(swx[2]+swx[4]+swx[5])+aa2^2*swx[3])/(s0[2]+aa1*(n0[2]+n0[4]+n0[5])+aa2*n0[3])
      aa6<-abs(sigma[1]-aaa0)
      aaa0<-sigma[1]
      if(n_iter>20)break
    }
    sigma[6]<-sigma[1];sigma[2]<-sigma[1]+g_a[1];sigma[4]<-sigma[1]+g_a[1]
    sigma[5]<-sigma[1]+g_a[1];sigma[3]<-sigma[1]+g_a[2]
    if(sum(sigma < 1e-30)>=1){break}
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
    if(iteration>600)break
  }
  abc<-L0
  AIC<--2*abc+2*3
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0
  #########hypothesis testing##############
  data <- sort(data);
  w1<-1/(12*m_sam)
  bmw <- matrix(0,m_sam,1); bmwsl <- matrix(0,m_sam,d5)
  for(i in 1:d5){
    gg <- (data - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2)))[1]
  if(nn < m_sam){P2 <- P2+runif(m_sam)/1e4}
  #####################################################################
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P2 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[2])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),D))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])

  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  output <- data.frame("2MG-EA",round(abc,4),round(AIC,4),round(t(m),4)," "," "," ",round(t(sigma),4)," "," "," ",round(t(mix_pi),4)," "," "," ",round(t(B81),4)," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4),round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
 output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################2MG-CD model#########################################
F3ModelFun[[10]]  <-function(K1,logL,df,m_nf){
  data <- as.matrix(as.numeric(as.numeric(df[,1])))
  ##############procedure start####################
  m_sam <- dim(data)[1]; d4 <- 9; m_esp <- 0.0001
  mean0 <- mean(data);sigma0 <- var(data)
  ######2MG-CD model(B-5)#####################
  a1 <- sqrt(sigma0)
  sigma <- matrix((sigma0/2),9,1)
  m <- as.matrix(c((mean0+3.2*a1),(mean0+2.4*a1),(mean0+1.6*a1),(mean0+0.8*a1),mean0,(mean0-0.8*a1),(mean0-1.6*a1),(mean0-2.4*a1),(mean0-3.2*a1)))
  mi <- as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,9,m_sam); swx <- matrix(0,9,1); hh <- matrix(0,6,6); b_line <- matrix(0,6,1)
  g_a<-matrix(0,3,1);s0 <- matrix(0,2,1)
  while(stopa > m_esp && iteration <= 1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d4) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    ############solve the linear equation############
    hh[1,1]<-9.0*sigma[1]/n0[1]+16.0*sigma[2]/n0[2]+sigma[3]/n0[3]
    hh[1,2]<-0.0
    hh[1,3]<--12.0*sigma[2]/n0[2]
    hh[1,4]<-3.0*sigma[1]/n0[1]+4.0*sigma[2]/n0[2]
    hh[1,5]<-3.0*sigma[1]/n0[1]-sigma[3]/n0[3]
    hh[1,6]<-0.0
    hh[2,2]<-9.0*sigma[4]/n0[4]+16.0*sigma[5]/n0[5]+sigma[6]/n0[6]
    hh[2,3]<-16.0*sigma[5]/n0[5]
    hh[2,4]<-0.0
    hh[2,5]<--3.0*sigma[4]/n0[4]+sigma[6]/n0[6]
    hh[2,6]<--4.0*sigma[5]/n0[5]-sigma[6]/n0[6]
    hh[3,3]<-9.0*sigma[2]/n0[2]+16.0*sigma[5]/n0[5]+sigma[8]/n0[8]
    hh[3,4]<--3.0*sigma[2]/n0[2]+sigma[8]/n0[8]
    hh[3,5]<-0
    hh[3,6]<--4.0*sigma[5]/n0[5]-sigma[8]/n0[8]
    hh[4,4]<-sigma[1]/n0[1]+sigma[2]/n0[2]+sigma[7]/n0[7]+sigma[8]/n0[8]
    hh[4,5]<-sigma[1]/n0[1]
    hh[4,6]<--sigma[8]/n0[8]
    hh[5,5]<-sigma[1]/n0[1]+sigma[3]/n0[3]+sigma[4]/n0[4]+sigma[6]/n0[6]
    hh[5,6]<--sigma[6]/n0[6]
    hh[6,6]<-sigma[5]/n0[5]+sigma[6]/n0[6]+sigma[8]/n0[8]+sigma[9]/n0[9]
    for(i in 2:6)
    {
      for(j in 1:(i-1))
      {
        hh[i,j]<-hh[j,i]
      }
    }
    #######################################################################################
    b_line[1]<-3*sumwx[1]/n0[1]-4*sumwx[2]/n0[2]+sumwx[3]/n0[3]
    b_line[2]<-3*sumwx[4]/n0[4]-4*sumwx[5]/n0[5]+sumwx[6]/n0[6]
    b_line[3]<-3*sumwx[2]/n0[2]-4*sumwx[5]/n0[5]+sumwx[8]/n0[8]
    b_line[4]<-sumwx[1]/n0[1]-sumwx[2]/n0[2]-sumwx[7]/n0[7]+sumwx[8]/n0[8]
    b_line[5]<-sumwx[1]/n0[1]-sumwx[3]/n0[3]-sumwx[4]/n0[4]+sumwx[6]/n0[6]
    b_line[6]<-sumwx[5]/n0[5]-sumwx[6]/n0[6]-sumwx[8]/n0[8]+sumwx[9]/n0[9]
    B9<- solve(hh,b_line)
    #####################################################################################
    m[1]<-(sumwx[1]-sigma[1]*(3.0*B9[1]+B9[4]+B9[5]))/n0[1]
    m[2]<-(sumwx[2]+sigma[2]*(4.0*B9[1]-3.0*B9[3]+B9[4]))/n0[2]
    m[3]<-(sumwx[3]-sigma[3]*(B9[1]-B9[5]))/n0[3]
    m[4]<-(sumwx[4]+sigma[4]*(-3.0*B9[2]+B9[5]))/n0[4]
    m[5]<-(sumwx[5]+sigma[5]*(4.0*B9[2]+4.0*B9[3]-B9[6]))/n0[5]
    m[6]<-(sumwx[6]+sigma[6]*(-B9[2]-B9[5]+B9[6]))/n0[6]
    m[7]<-(sumwx[7]+sigma[7]*B9[4])/n0[7]
    m[8]<-(sumwx[8]+sigma[8]*(-B9[3]-B9[4]+B9[6]))/n0[8]
    m[9]<-(sumwx[9]-sigma[9]*B9[6])/n0[9]
    #########first order genetic parameter process##########
    hh91 <- matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,0.5,0.5,0.5,-1,-1,-1,1,0.5,-1,1,0.5,-1,1,0.5,-1),9,3)
    B91<-solve(crossprod(hh91,hh91))%*%crossprod(hh91,m)
    g_a[1]<-0.75*B91[3]^2/m_nf
    g_a[2]<-0.75*B91[2]^2/m_nf
    g_a[3]<-g_a[1]+g_a[2]
    for(i in 1:d4) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0[1]<-swx[1]+swx[3]+swx[7]+swx[9]
    s0[2]<-n0[1]+n0[3]+n0[7]+n0[9]
    aaa0<-sigma[1]
    n_iter<-0
    aa6<-1000
    while(aa6>0.0001){
      n_iter<-n_iter+1
      aa1<-sigma[1]/(sigma[1]+g_a[1])
      aa2<-sigma[1]/(sigma[1]+g_a[2])
      aa3<-sigma[1]/(sigma[1]+g_a[3])
      sigma[1]<-(s0[1]+aa1^2*(swx[2]+swx[8])+aa2^2*(swx[4]+swx[6])+aa3^2*swx[5])/(s0[2]+aa1*(n0[2]+n0[8])+aa2*(n0[4]+n0[6])+aa3*n0[5])
      aa6<-abs(sigma[1]-aaa0)
      aaa0<-sigma[1]
      if(n_iter>20)break
    }
    sigma[3]<-sigma[1];sigma[7]<-sigma[1];sigma[9]<-sigma[1]
    sigma[2]<-sigma[1]+g_a[1];sigma[4]<-sigma[1]+g_a[2]
    sigma[5]<-sigma[1]+g_a[3];sigma[6]<-sigma[1]+g_a[2];sigma[8]<-sigma[1]+g_a[1]
    if(sum(sigma < 1e-30)>=1){break}
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*4
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0
  #########hypothesis testing##############
  data <- sort(data)
  w1<-1/(12*m_sam)
  bmw <- matrix(0,m_sam,1); bmwsl <- matrix(0,m_sam,d4)
  for(i in 1:d4){
    gg <- (data - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2)))[1]
  if(nn < m_sam){P2 <- P2+runif(m_sam)/1e4}
  #####################################################################
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P2 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[2])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),D))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])

  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  output <- data.frame("2MG-CD",round(abc,4),round(AIC,4),round(t(m),4),round(t(sigma),4),round(t(mix_pi),4),round(t(B91),4)," "," "," "," "," "," ",round(jj,4),round(ll*100,4),round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################2MG-EAD model#########################################
F3ModelFun[[11]]  <-function(K1,logL,df,m_nf){
  data <- as.matrix(as.numeric(as.numeric(df[,1])))
  ##############procedure start####################
  m_sam <- dim(data)[1]; d5 <- 6; m_esp <- 0.0001
  mean0 <- mean(data);sigma0 <- var(data)
  #########2MG-EAD model(B-6)#########################
  a1 <- sqrt(sigma0/m_sam)
  sigma <- matrix((sigma0/2),6,1)
  m <- as.matrix(c((mean0+3*a1),(mean0+2.4*a1),(mean0+1.7*a1),mean0,(mean0-0.7*a1),(mean0-3*a1)))
  mi <- as.matrix(c(0.0625,0.25,0.25,0.125,0.25,0.0625))
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,6,m_sam); swx <- matrix(0,6,1); hh <- matrix(0,4,4); b_line <- matrix(0,4,1)
  g_a<-matrix(0,2,1);s0 <- matrix(0,2,1)
  while(stopa > m_esp && iteration <= 1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d5) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    ############solve the linear equation############
    hh[1,1]<-sigma[1]/n0[1]+4.0*sigma[2]/n0[2]+sigma[3]/n0[3]
    hh[1,2]<-sigma[1]/n0[1]
    hh[1,3]<--2.0*sigma[2]/n0[2]-sigma[3]/n0[3]
    hh[1,4]<-2.0*sigma[3]/n0[3]
    hh[2,2]<-sigma[1]/n0[1]+4.0*sigma[4]/n0[4]+sigma[6]/n0[6]
    hh[2,3]<-2.0*sigma[4]/n0[4]
    hh[2,4]<-6.0*sigma[4]/n0[4]+sigma[6]/n0[6]
    hh[3,3]<-sigma[2]/n0[2]+sigma[3]/n0[3]+sigma[4]/n0[4]+sigma[5]/n0[5]
    hh[3,4]<--2.0*sigma[3]/n0[3]+3.0*sigma[4]/n0[4]
    hh[4,4]<-4.0*sigma[3]/n0[3]+9.0*sigma[4]/n0[4]+sigma[6]/n0[6]
    for(i in 2:4)
    {
      for(j in 1:(i-1))
      {
        hh[i,j]<-hh[j,i]
      }
    }
    #######################################################################
    b_line[1]<-sumwx[1]/n0[1]-2.0*sumwx[2]/n0[2]+sumwx[3]/n0[3]
    b_line[2]<-sumwx[1]/n0[1]-2.0*sumwx[4]/n0[4]+sumwx[6]/n0[6]
    b_line[3]<-sumwx[2]/n0[2]-sumwx[3]/n0[3]-sumwx[4]/n0[4]+sumwx[5]/n0[5]
    b_line[4]<-2.0*sumwx[3]/n0[3]-3.0*sumwx[4]/n0[4]+sumwx[6]/n0[6]
    B10<- solve(hh,b_line)
    #####################################################################
    m[1]<-(sumwx[1]-sigma[1]*(B10[1]+B10[2]))/n0[1]
    m[2]<-(sumwx[2]+sigma[2]*(2.0*B10[1]-B10[3]))/n0[2]
    m[3]<-(sumwx[3]-sigma[3]*(B10[1]-B10[3]+2.0*B10[4]))/n0[3]
    m[4]<-(sumwx[4]+sigma[4]*(2.0*B10[2]+B10[3]+3.0*B10[4]))/n0[4]
    m[5]<-(sumwx[5]-sigma[5]*B10[3])/n0[5]
    m[6]<-(sumwx[6]-sigma[6]*(B10[2]+B10[4]))/n0[6]
    #########first order genetic parameter process##########
    hh101 <- matrix(c(1,1,1,1,1,1,2,1.5,1,0,-0.5,-2),6,2)
    B101<-solve(crossprod(hh101,hh101))%*%crossprod(hh101,m)
    g_a[1]<-0.75*B101[2]^2/m_nf
    g_a[2]<-1.5*B101[2]^2/m_nf
    for(i in 1:d5) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0[1]<-swx[1]+swx[4]+swx[6]
    s0[2]<-n0[1]+n0[4]+n0[6]
    aaa0<-sigma[1]
    n_iter<-0
    aa6<-1000
    while(aa6>0.0001){
      n_iter<-n_iter+1
      aa1<-sigma[1]/(sigma[1]+g_a[1])
      aa2<-sigma[1]/(sigma[1]+g_a[2])
      sigma[1]<-(s0[1]+aa1^2*(swx[2]+swx[5])+aa2^2*swx[3])/(s0[2]+aa1*(n0[2]+n0[5])+aa2*n0[3])
      aa6<-abs(sigma[1]-aaa0)
      aaa0<-sigma[1]
      if(n_iter>20)break
    }
    sigma[4]<-sigma[1];sigma[6]<-sigma[1]
    sigma[2]<-sigma[1]+g_a[1];sigma[5]<-sigma[1]+g_a[1]
    sigma[3]<-sigma[1]+g_a[2]
    if(sum(sigma < 1e-30)>=1){break}
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*3
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0
  #########hypothesis testing##############
  data <- sort(data)
  w1<-1/(12*m_sam)
  bmw <- matrix(0,m_sam,1); bmwsl <- matrix(0,m_sam,d5)
  for(i in 1:d5){
    gg <- (data - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2)))[1]
  if(nn < m_sam){P2 <- P2+runif(m_sam)/1e4}
  #####################################################################
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P2 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[2])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),D))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])

  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  output <- data.frame("2MG-EAD",round(abc,4),round(AIC,4),round(t(m),4)," "," "," ",round(t(sigma),4)," "," "," ",round(t(mix_pi),4)," "," "," ",round(t(B101),4)," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4),round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}

K1F23 <- function(x){
  V0 <- 0
  for(j in 0:2)
  {I1 <- 0;I2 <- 0
  for(k in 0:8)
  {I1 <- I1 + (((4*j+1)^2/(32*x))^(-0.25+2*k))/(gamma(k+1)*gamma(0.75+k))
  I2 <- I2 + ((4*j+1)^2/(32*x))^(0.25+2*k)/(gamma(k+1)*gamma(1.25+k))}
  V0 <- V0 + (gamma(j+0.5)*sqrt(4*j+1)/(gamma(0.5)*gamma(j+1)))*exp(-(4*j+1)^2/(16*x))*(I1-I2)}
  V <- (1/sqrt(2*x))*V0
  return (1-V)
}

logLF23 <- function(nm,nng,mi,mn,s,d1) { sum2 <- sum(log(dmixnorm(d1,mn,sqrt(s),mi)));return (sum2) }


if(model=="All models"){
  cl.cores <- detectCores()
  if(cl.cores<=2){
    cl.cores<-1
  }else if(cl.cores>2){
    if(cl.cores>10){
      cl.cores<-10
    }else {
      cl.cores <- detectCores()-1
    }
  }
  cl <- makeCluster(cl.cores)
  registerDoParallel(cl)
  i<-NULL
  allresult=foreach(i=1:11,.combine = 'rbind')%dopar%{
    requireNamespace("KScorrect")
    F3ModelFun[[i]](K1F23,logLF23,df,m_nf)[[1]]
  }
  stopCluster(cl)
  mi<-NULL
}else{
  allresultq<-switch(model,"0MG" = F3ModelFun[[1]](K1F23,logLF23,df,m_nf),"1MG-AD"=F3ModelFun[[2]](K1F23,logLF23,df,m_nf),"1MG-A"=F3ModelFun[[3]](K1F23,logLF23,df,m_nf),
                     "1MG-EAD"=F3ModelFun[[4]](K1F23,logLF23,df,m_nf),"1MG-NCD"=F3ModelFun[[5]](K1F23,logLF23,df,m_nf),"2MG-ADI"=F3ModelFun[[6]](K1F23,logLF23,df,m_nf),
                     "2MG-AD"=F3ModelFun[[7]](K1F23,logLF23,df,m_nf),"2MG-A"=F3ModelFun[[8]](K1F23,logLF23,df,m_nf),"2MG-EA"=F3ModelFun[[9]](K1F23,logLF23,df,m_nf),
                     "2MG-CD"=F3ModelFun[[10]](K1F23,logLF23,df,m_nf),"2MG-EAD"=F3ModelFun[[11]](K1F23,logLF23,df,m_nf))

  allresult<-allresultq[[1]]
  if(model!="0MG"){
    mi<-allresultq[[2]]
  }else{
    mi<-NULL
  }
}
  colnames(allresult) <- F3colname
  out<-list(allresult,mi)
  return(out)

}






















































