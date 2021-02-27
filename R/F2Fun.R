F2Fun<-function(df,model){

data<-sapply(df,as.character)

dF2<-data[-1,which(data[1,]=="F2")];F2<-as.numeric(dF2[which(is.na(as.numeric(dF2))==FALSE)]);df<-as.data.frame(F2)


F2colname <- c("Model","Log_Max_likelihood_Value","AIC","mean[1]","mean[2]","mean[3]","mean[4]","mean[5]","mean[6]","mean[7]","mean[8]","mean[9]",
               "Var(Residual+Polygene)","Proportion[1]","Proportion[2]","Proportion[3]","Proportion[4]","Proportion[5]","Proportion[6]",
               "Proportion[7]","Proportion[8]","Proportion[9]","m","da(d)","db","ha(h)","hb","i","jab","jba","l","Major-Gene Var","Heritability(Major-Gene)(%)",
               "U1 square","P(U1 square)","U2 square","P(U2 square)","U3 square","P(U3 square)","nW square","P(nW square)","Dn","P(Dn)")

F2ModelFun<-list(NA)
###################################0MG Model##################################
F2ModelFun[[1]] <- function(K1,logL,df){

  data<-as.matrix(as.numeric(df[,1]))
  m_sam <- dim(data)[1]; mean0 <- mean(data); sigma0 <- (sd(data))^2
  #################0MG Model##################
  d1 <- 1; m <- mean0; sigma <- sigma0; mix_pi <- 1.0
  abc <- sum(log(dnorm(data,mean0,sqrt(sigma0))))
  AICm <- -2.0*abc + 4
  ################hypothesis testing################
  data <- sort(data); P1 <- matrix(0,m_sam,1)
  gg <- (data - m)/sqrt(as.vector(sigma))
  P1[which(gg>=0)] <- pnorm(gg[gg>=0])
  P1[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
  #############deal with ties in P1##############
  nn <- dim(as.matrix(unique(P1)))[1]
  if(nn < m_sam) {P1 <- P1 + runif(m_sam)/1e4}
  ###############################################
  dd <- as.matrix(c(sum(P1),sum(P1^2),sum((P1-0.5)^2)))
  WW <- 1/(12*m_sam) + sum((P1 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P1,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW),(1-pkolm(D,m_sam))))

  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  output <- data.frame("0MG",round(abc,4),round(AICm,4),round(m,4)," "," "," "," "," "," "," "," ",round(sigma,4),round(mix_pi,4)," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW,4),tt[4],round(D,4), tt[5])

  output<-as.matrix(output)
  OUTPUT<-list(output)
  return(OUTPUT)
}

######################################1MG-AD Model##########################
F2ModelFun[[2]] <- function(K1,logL,df){
  data<-as.matrix(as.numeric(df[,1]))
  ##############procedure start####################
   m_sam <- dim(data)[1]; d2 <- 3;m_esp <- 0.0001
   mean0 <- mean(data);sigma0 <- (sd(data))^2
  #########1MG-AD Model#########
  a1 <- sqrt(sigma0)
  sigma <- matrix((sigma0/2),3,1)
  m <- as.matrix(c((mean0+2*a1),mean0,(mean0-2*a1)))
  mi <- as.matrix(c(0.25,0.5,0.25))
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,3,m_sam); swx <- matrix(0,3,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[n0<0.000001] <- 0.000001
    m <- sumwx/n0
    ##########obtain variance##########
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0 <- sum(swx)
    sigma <- matrix(c(swx/n0),3,1)
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- logL(m_sam,d2,mix_pi,m,sigma,data)
  AICm <- -2*abc + 8
  #########first order genetic parameter process##########
  hh1 <- matrix(c(1,1,1,1,0,-1,0,1,0),3,3)
  b_line1 <- m; B1 <- solve(t(hh1)%*%hh1)%*%(t(hh1)%*%b_line1)
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0
  #########hypothesis testing##############
  data <- sort(data); bmw <- matrix(0,m_sam,1); bmwsl <- matrix(0,m_sam,d2)
  for(i in 1:d2){
    gg <- (data - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  #############deal with ties in P2##############
  nn <- dim(as.matrix(unique(P2)))[1]
  if(nn < m_sam) {P2 <- P2 + runif(m_sam)/1e4}
  ###############################################
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P2 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,m_sam))))

  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  output <- data.frame("1MG-AD",round(abc,4),round(AICm,4),round(t(m),4)," "," "," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4)," "," "," "," "," "," ",round(B1[1],4),round(B1[2],4)," ",round(B1[3],4)," "," "," "," "," ",round(jj,4),round(ll*100,4),round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4), tt[5])


  output<-as.matrix(output)

  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
#################################1MG-A Model###########################################
F2ModelFun[[3]] <- function(K1,logL,df)
{
  data<-as.matrix(as.numeric(df[,1]))
  ##############procedure start####################
  m_sam <- dim(data)[1]; d2 <- 3;m_esp <- 0.0001
  mean0 <- mean(data);sigma0 <- (sd(data))^2
  #########1MG-A Model#########
  a1 <- sqrt(sigma0/m_sam)
  sigma <- matrix((sigma0/2),3,1)
  m <- as.matrix(c((mean0+2*a1),mean0,(mean0-2*a1)))
  mi <- as.matrix(c(0.25,0.5,0.25))
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,3,m_sam); swx <- matrix(0,3,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[n0<0.000001] <- 0.000001
    aa1 <- sigma[1]*(1/n0[1]+4/n0[2]+1/n0[3])
    aa2 <- sumwx[1]/n0[1]-2*sumwx[2]/n0[2]+sumwx[3]/n0[3]
    aa3 <- aa2/aa1
    m[1] <- (sumwx[1]-sigma[1]*aa3)/n0[1]
    m[2] <- (sumwx[2]+sigma[1]*2*aa3)/n0[2]
    m[3] <- (sumwx[3]-sigma[1]*aa3)/n0[3]
    ##########obtain variance##########
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0 <- sum(swx)
    sigma <- matrix(c(swx/n0),3,1)
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- logL(m_sam,d2,mix_pi,m,sigma,data)
  AICm <- -2*abc + 6
  #########first order genetic parameter process##########
  hh1 <- matrix(c(1,1,1,1,0,-1),3,2)
  B1 <- solve(t(hh1)%*%hh1)%*%(t(hh1)%*%m)
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0
  #########hypothesis testing##############
  data <- sort(data); bmw <- matrix(0,m_sam,1); bmwsl <- matrix(0,m_sam,d2)
  for(i in 1:d2){
    gg <- (data - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P3 <- rowSums(bmwsl)
  #############deal with ties in P3##############
  nn <- dim(as.matrix(unique(P3)))[1]
  if(nn < m_sam) {P3 <- P3 + runif(m_sam)/1e4}
  ###############################################
  dd <- as.matrix(c(sum(P3),sum(P3^2),sum((P3-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P3 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P3,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,m_sam))))

  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  output <- data.frame("1MG-A",round(abc,4),round(AICm,4),round(t(m),4)," "," "," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4)," "," "," "," "," "," ",round(B1[1],4),round(B1[2],4)," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4),round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4), tt[5])


  output<-as.matrix(output)

  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################1MG-EAD model#########################################
F2ModelFun[[4]] <- function(K1,logL,df){
  data<-as.matrix(as.numeric(df[,1]))
  ##############procedure start####################
  m_sam <- dim(data)[1]; d3 <- 2;m_esp <- 0.0001
  mean0 <- mean(data);sigma0 <- (sd(data))^2
  #########1MG-EAD Model#########
  a1 <- sqrt(sigma0/m_sam)
  sigma <- matrix((sigma0/2),2,1)
  m <- as.matrix(c((mean0+2*a1),(mean0-2*a1)))
  mi <- as.matrix(c(0.75,0.25))
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,2,m_sam); swx <- matrix(0,2,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d3) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[n0<0.000001] <- 0.000001
    m <- sumwx/n0
    ##########obtain variance##########
    for(i in 1:d3) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0 <- sum(swx)
    sigma <- matrix(c(swx/n0),2,1)
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- logL(m_sam,d3,mix_pi,m,sigma,data)
  AICm <- -2*abc + 8
  #########first order genetic parameter process##########
  hh1 <- matrix(c(1,1,1,-1),2,2)
  B1 <- solve(t(hh1)%*%hh1)%*%(t(hh1)%*%m)
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0
  #########hypothesis testing##############
  data <- sort(data); bmw <- matrix(0,m_sam,1); bmwsl <- matrix(0,m_sam,d3)
  for(i in 1:d3){
    gg <- (data - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P4 <- rowSums(bmwsl)
  #############deal with ties in P4##############
  nn <- dim(as.matrix(unique(P4)))[1]
  if(nn < m_sam) {P4 <- P4 + runif(m_sam)/1e4}
  ###############################################
  dd <- as.matrix(c(sum(P4),sum(P4^2),sum((P4-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P4 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P4,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,m_sam))))

  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  output <- data.frame("1MG-EAD",round(abc,4),round(AICm,4),round(t(m),4)," "," "," "," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4)," "," "," "," "," "," "," ",round(B1[1],4),round(B1[2],4)," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4),round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4), tt[5])


  output<-as.matrix(output)

  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################1MG-NCD model#########################################
F2ModelFun[[5]] <- function(K1,logL,df){
  data<-as.matrix(as.numeric(df[,1]))
  ##############procedure start####################
  m_sam <- dim(data)[1]; d3 <- 2;m_esp <- 0.0001
  mean0 <- mean(data);sigma0 <- (sd(data))^2
  #########1MG-NCD Model#########
  a1 <- sqrt(sigma0/m_sam)
  sigma <- matrix((sigma0/2),2,1)
  m <- as.matrix(c((mean0-2*a1),(mean0+2*a1)))
  mi <- as.matrix(c(0.75,0.25))
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,2,m_sam); swx <- matrix(0,2,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d3) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[n0<0.000001] <- 0.000001
    m <- sumwx/n0
    ##########obtain variance##########
    for(i in 1:d3) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0 <- sum(swx)
    sigma <- matrix(s0/m_sam,2,1)
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- logL(m_sam,d3,mix_pi,m,sigma,data)
  AICm <- -2*abc + 8
  #########first order genetic parameter process##########
  hh1 <- matrix(c(1,1,-1,1),2,2)
  B1 <- solve(t(hh1)%*%hh1)%*%(t(hh1)%*%m)
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0
  #########hypothesis testing##############
  data <- sort(data); bmw <- matrix(0,m_sam,1); bmwsl <- matrix(0,m_sam,d3)
  for(i in 1:d3){
    gg <- (data - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P4 <- rowSums(bmwsl)
  #############deal with ties in P4##############
  nn <- dim(as.matrix(unique(P4)))[1]
  if(nn < m_sam) {P4 <- P4 + runif(m_sam)/1e4}
  ###############################################
  dd <- as.matrix(c(sum(P4),sum(P4^2),sum((P4-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P4 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P4,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,m_sam))))

  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  output <- data.frame("1MG-NCD",round(abc,4),round(AICm,4),round(t(m),4)," "," "," "," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4)," "," "," "," "," "," "," ",round(B1[1],4),round(B1[2],4)," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4),round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4), tt[5])


  output<-as.matrix(output)

  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################2MG-ADI model#########################################
F2ModelFun[[6]] <- function(K1,logL,df)
{
  data<-as.matrix(as.numeric(df[,1]))
  ##############procedure start####################
  m_sam <- dim(data)[1]; d4 <- 9;m_esp <- 0.0001
  mean0 <- mean(data);sigma0 <- (sd(data))^2
  #########2MG-ADI Model#########
  a1 <- sqrt(sigma0/m_sam)
  sigma <- matrix((sigma0/2),9,1)
  m <- as.matrix(c((mean0+3.2*a1),(mean0+2.4*a1),(mean0+1.6*a1),(mean0+0.8*a1),mean0,
                   (mean0-0.8*a1),(mean0-1.6*a1),(mean0-2.4*a1),(mean0-3.2*a1)))
  mi <- as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,9,m_sam); swx <- matrix(0,9,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d4) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[n0<0.000001] <- 0.000001
    m <- sumwx/n0
    ##########obtain variance##########
    for(i in 1:d4) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0 <- sum(swx)
    sigma <- matrix(s0/m_sam,9,1)
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- logL(m_sam,d4,mix_pi,m,sigma,data)
  AICm <- -2*abc + 20
  #########first order genetic parameter process##########
  hh5 <- matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,-1,-1,-1,1,0,-1,1,0,-1,1,0,-1,
                  0,0,0,1,1,1,0,0,0,0,1,0,0,1,0,0,1,0,1,0,-1,0,0,0,-1,0,1,0,1,
                  0,0,0,0,0,-1,0,0,0,0,1,0,-1,0,0,0,0,0,0,0,1,0,0,0,0),9,9)
  b_line5 <- m; B5 <- solve(t(hh5)%*%hh5)%*%(t(hh5)%*%b_line5)
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0
  #########hypothesis testing##############
  data <- sort(data); bmw <- matrix(0,m_sam,1); bmwsl <- matrix(0,m_sam,d4)
  for(i in 1:d4){
    gg <- (data - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P6 <- rowSums(bmwsl)
  #############deal with ties in P6##############
  nn <- dim(as.matrix(unique(P6)))[1]
  if(nn < m_sam) {P6 <- P6 + runif(m_sam)/1e4}
  ###############################################
  dd <- as.matrix(c(sum(P6),sum(P6^2),sum((P6-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P6 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P6,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,m_sam))))

  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  output <- data.frame("2MG-ADI",round(abc,4),round(AICm,4),round(t(m),4),round(sigma[1],4),round(t(mix_pi),4),round(t(B5),4),round(jj,4),round(ll*100,4),round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4), tt[5])


  output<-as.matrix(output)

  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################2MG-AD model#########################################
F2ModelFun[[7]] <- function(K1,logL,df){
  data<-as.matrix(as.numeric(df[,1]))
  ##############procedure start####################
  m_sam <- dim(data)[1]; d4 <- 9;m_esp <- 0.0001
  mean0 <- mean(data);sigma0 <- (sd(data))^2
  #########2MG-AD Model#########
  a1 <- sqrt(sigma0/m_sam)
  sigma <- matrix((sigma0/2),9,1)
  m <- as.matrix(c((mean0+3.2*a1),(mean0+2.4*a1),(mean0+1.6*a1),(mean0+0.8*a1),mean0,
                   (mean0-0.8*a1),(mean0-1.6*a1),(mean0-2.4*a1),(mean0-3.2*a1)))
  mi <- as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,9,m_sam); swx <- matrix(0,9,1); hh6 <- matrix(0,4,4); b_line6 <- matrix(0,4,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d4) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[n0<0.000001] <- 0.000001
    ############solve the linear equation############
    hh6[1,1] <- sigma[1]/n0[1]+sigma[3]/n0[3]+sigma[4]/n0[4]+sigma[6]/n0[6]
    hh6[1,2] <- sigma[3]/n0[3]
    hh6[1,3] <- sigma[1]/n0[1]+sigma[3]/n0[3]
    hh6[1,4] <- -sigma[1]/n0[1]-sigma[4]/n0[4]+sigma[6]/n0[6]
    hh6[2,2] <- sigma[2]/n0[2]+sigma[3]/n0[3]+sigma[8]/n0[8]+sigma[9]/n0[9]
    hh6[2,3] <- sigma[3]/n0[3]+sigma[9]/n0[9]
    hh6[2,4] <- sigma[2]/n0[2]-sigma[8]/n0[8]-sigma[9]/n0[9]
    hh6[3,3] <- sigma[1]/n0[1]+sigma[3]/n0[3]+sigma[7]/n0[7]+sigma[9]/n0[9]
    hh6[3,4] <- -sigma[1]/n0[1]-sigma[9]/n0[9]
    hh6[4,4] <- sigma[1]/n0[1]+sigma[2]/n0[2]+sigma[4]/n0[4]+4.0*sigma[5]/n0[5]+sigma[6]/n0[6]+sigma[8]/n0[8]+sigma[9]/n0[9]
    for(i in 2:4)
    {
      for(j in 1:(i-1))
      {
        hh6[i,j]<-hh6[j,i]
      }
    }
    #################################################
    b_line6[1] <- sumwx[1]/n0[1]-sumwx[3]/n0[3]-sumwx[4]/n0[4]+sumwx[6]/n0[6]
    b_line6[2] <- sumwx[2]/n0[2]-sumwx[3]/n0[3]-sumwx[8]/n0[8]+sumwx[9]/n0[9]
    b_line6[3] <- sumwx[1]/n0[1]-sumwx[3]/n0[3]-sumwx[7]/n0[7]+sumwx[9]/n0[9]
    b_line6[4] <- sumwx[2]/n0[2]+sumwx[4]/n0[4]+sumwx[6]/n0[6]+sumwx[8]/n0[8]-sumwx[1]/n0[1]-2.0*sumwx[5]/n0[5]-sumwx[9]/n0[9]
    B6 <- solve(hh6,b_line6)
    #################################################
    m[1] <- (sumwx[1]-sigma[1]*(B6[1]+B6[3]-B6[4]))/n0[1]; m[2] <- (sumwx[2]-sigma[2]*(B6[2]+B6[4]))/n0[2]
    m[3] <- (sumwx[3]+sigma[3]*(B6[1]+B6[2]+B6[3]))/n0[3]; m[4] <- (sumwx[4]+sigma[4]*(B6[1]-B6[4]))/n0[4]
    m[5] <- (sumwx[5]+sigma[5]*2.0*B6[4])/n0[5]; m[6] <- (sumwx[6]-sigma[6]*(B6[1]+B6[4]))/n0[6]
    m[7] <- (sumwx[7]+sigma[7]*B6[3])/n0[7]; m[8] <- (sumwx[8]+sigma[8]*(B6[2]-B6[4]))/n0[8]
    m[9] <- (sumwx[9]-sigma[9]*(B6[2]+B6[3]-B6[4]))/n0[9]
    ##########obtain variance##########
    for(i in 1:d4) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0 <- sum(swx)
    sigma <- matrix(swx/n0,9,1)
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- logL(m_sam,d4,mix_pi,m,sigma,data)
  AICm <- -2*abc + 12
  #########first order genetic parameter process##########
  hh61 <- matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,-1,-1,-1,1,0,-1,1,0,
                   -1,1,0,-1,0,0,0,1,1,1,0,0,0,0,1,0,0,1,0,0,1,0),9,5)
  B61 <- solve(t(hh61)%*%hh61)%*%(t(hh61)%*%m)
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0
  #########hypothesis testing##############
  data <- sort(data); bmw <- matrix(0,m_sam,1); bmwsl <- matrix(0,m_sam,d4)
  for(i in 1:d4){
    gg <- (data - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P6 <- rowSums(bmwsl)
  #############deal with ties in P6##############
  nn <- dim(as.matrix(unique(P6)))[1]
  if(nn < m_sam) {P6 <- P6 + runif(m_sam)/1e4}
  ###############################################
  dd <- as.matrix(c(sum(P6),sum(P6^2),sum((P6-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P6 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P6,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,m_sam))))

  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  output <- data.frame("2MG-AD",round(abc,4),round(AICm,4),round(t(m),4),round(sigma[1],4),round(t(mix_pi),4),round(t(B61),4)," "," "," "," ",round(jj,4),round(ll*100,4),round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4), tt[5])


  output<-as.matrix(output)

  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################2MG-A model#########################################
F2ModelFun[[8]] <- function(K1,logL,df){
  data<-as.matrix(as.numeric(df[,1]))
  ##############procedure start####################
  m_sam <- dim(data)[1]; d4 <- 9;m_esp <- 0.0001
  mean0 <- mean(data);sigma0 <- (sd(data))^2
  #########2MG-A Model#########
  a1 <- sqrt(sigma0/m_sam)
  sigma <- matrix((sigma0/2),9,1)
  m <- as.matrix(c((mean0+3.2*a1),(mean0+2.4*a1),(mean0+1.6*a1),(mean0+0.8*a1),mean0,
                   (mean0-0.8*a1),(mean0-1.6*a1),(mean0-2.4*a1),(mean0-3.2*a1)))
  mi <- as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,9,m_sam); swx <- matrix(0,9,1); hh7 <- matrix(0,6,6); b_line7 <- matrix(0,6,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d4) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[n0<0.000001] <- 0.000001
    ############solve the linear equation############
    hh7[1,1] <- sigma[1]/n0[1]+4.0*sigma[2]/n0[2]+sigma[3]/n0[3];
    hh7[1,2] <- sigma[1]/n0[1]
    hh7[1,3] <- -2.0*sigma[2]/n0[2]-sigma[3]/n0[3];
    hh7[1,4] <- sigma[3]/n0[3];
    hh7[1,5] <- 0.0;
    hh7[1,6] <- 0.0
    hh7[2,2] <- sigma[1]/n0[1]+4.0*sigma[5]/n0[5]+sigma[9]/n0[9]
    hh7[2,3] <- sigma[9]/n0[9];
    hh7[2,4] <- 4.0*sigma[5]/n0[5];
    hh7[2,5] <- 4.0*sigma[5]/n0[5];
    hh7[2,6] <- 2.0*sigma[9]/n0[9]
    hh7[3,3] <- sigma[2]/n0[2]+sigma[3]/n0[3]+sigma[8]/n0[8]+sigma[9]/n0[9]
    hh7[3,4] <- -sigma[3]/n0[3];
    hh7[3,5] <- 0.0;
    hh7[3,6] <- 2.0*sigma[8]/n0[8]+2.0*sigma[9]/n0[9]
    hh7[4,4] <- sigma[3]/n0[3]+4.0*sigma[5]/n0[5]+sigma[7]/n0[7]
    hh7[4,5] <- 4.0*sigma[5]/n0[5]
    hh7[4,6] <- 0.0
    hh7[5,5] <- sigma[4]/n0[4]+4.0*sigma[5]/n0[5]+sigma[6]/n0[6]
    hh7[5,6] <- sigma[4]/n0[4]-sigma[6]/n0[6]
    hh7[6,6] <- sigma[4]/n0[4]+sigma[6]/n0[6]+4.0*sigma[8]/n0[8]+4.0*sigma[9]/n0[9]
    for(i in 2:6)
    {
      for(j in 1:(i-1))
      {
        hh7[i,j]<-hh7[j,i]
      }
    }
    #################################################
    b_line7[1] <- sumwx[1]/n0[1]-2.0*sumwx[2]/n0[2]+sumwx[3]/n0[3]
    b_line7[2] <- sumwx[1]/n0[1]-2.0*sumwx[5]/n0[5]+sumwx[9]/n0[9]
    b_line7[3] <- sumwx[2]/n0[2]-sumwx[3]/n0[3]-sumwx[8]/n0[8]+sumwx[9]/n0[9]
    b_line7[4] <- sumwx[3]/n0[3]-2.0*sumwx[5]/n0[5]+sumwx[7]/n0[7]
    b_line7[5] <- sumwx[4]/n0[4]-2.0*sumwx[5]/n0[5]+sumwx[6]/n0[6]
    b_line7[6] <- sumwx[4]/n0[4]-sumwx[6]/n0[6]-2.0*sumwx[8]/n0[8]+2.0*sumwx[9]/n0[9]
    B7 <- solve(t(hh7)%*%hh7)%*%(t(hh7)%*%b_line7)
    #################################################
    m[1] <- (sumwx[1]-sigma[1]*(B7[1]+B7[2]))/n0[1]; m[2] <- (sumwx[2]+sigma[2]*(2.0*B7[1]-B7[3]))/n0[2]
    m[3] <- (sumwx[3]-sigma[3]*(B7[1]-B7[3]+B7[4]))/n0[3]; m[4] <- (sumwx[4]-sigma[4]*(B7[5]+B7[6]))/n0[4]
    m[5] <- (sumwx[5]+sigma[5]*2.0*(B7[2]+B7[4]+B7[5]))/n0[5]; m[6] <- (sumwx[6]-sigma[6]*(B7[5]+B7[6]))/n0[6]
    m[7] <- (sumwx[7]-sigma[7]*B7[4])/n0[7]; m[8] <- (sumwx[8]+sigma[8]*(B7[3]+2.0*B7[6]))/n0[8]
    m[9] <- (sumwx[6]-sigma[6]*(B7[2]+B7[3]+2*B7[6]))/n0[6]
    ##########obtain variance##########
    for(i in 1:d4) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0 <- sum(swx)
    sigma <- matrix(c(swx/n0),9,1)
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- logL(m_sam,d4,mix_pi,m,sigma,data)
  AICm <- -2*abc + 8
  #########first order genetic parameter process##########
  hh71 <- matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,-1,-1,-1,1,0,-1,1,0,-1,1,0,-1),9,3)
  B71 <- solve(t(hh71)%*%hh71)%*%(t(hh71)%*%m)
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0
  #########hypothesis testing##############
  data <- sort(data); bmw <- matrix(0,m_sam,1); bmwsl <- matrix(0,m_sam,d4)
  for(i in 1:d4){
    gg <- (data - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P6 <- rowSums(bmwsl)
  #############deal with ties in P6##############
  nn <- dim(as.matrix(unique(P6)))[1]
  if(nn < m_sam) {P6 <- P6 + runif(m_sam)/1e4}
  ###############################################
  dd <- as.matrix(c(sum(P6),sum(P6^2),sum((P6-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P6 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P6,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,m_sam))))

  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  output <- data.frame("2MG-A",round(abc,4),round(AICm,4),round(t(m),4),round(sigma[1],4),round(t(mix_pi),4),round(t(B71),4)," "," "," "," "," "," ",round(jj,4),round(ll*100,4),round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4), tt[5])


  output<-as.matrix(output)

  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################2MG-EA model#########################################
F2ModelFun[[9]] <- function(K1,logL,df){
  data<-as.matrix(as.numeric(df[,1]))
  ##############procedure start####################
  m_sam <- dim(data)[1]; d5 <- 5;m_esp <- 0.0001
  mean0 <- mean(data);sigma0 <- (sd(data))^2
  #########2MG-EA Model#########
  a1 <- sqrt(sigma0/m_sam)
  sigma <- matrix((sigma0/2),5,1)
  m <- as.matrix(c((mean0+3*a1),(mean0+1.5*a1),mean0,(mean0-1.5*a1),(mean0-3*a1)))
  mi <- as.matrix(c(0.0625,0.25,0.375,0.25,0.0625))
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,5,m_sam); swx <- matrix(0,5,1); hh8 <- matrix(0,3,3); b_line8 <- matrix(0,3,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d5) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[n0<0.000001] <- 0.000001
    ############solve the linear equation############
    hh8[1,1] <- sigma[1]/n0[1]+4.0*sigma[2]/n0[2]+sigma[3]/n0[3]
    hh8[1,2] <- -2.0*sigma[2]/n0[2]-2.0*sigma[3]/n0[3]
    hh8[1,3] <- sigma[3]/n0[3]
    hh8[2,2] <- sigma[2]/n0[2]+4.0*sigma[3]/n0[3]+sigma[4]/n0[4]
    hh8[2,3] <- -2.0*sigma[3]/n0[3]-2.0*sigma[4]/n0[4]
    hh8[3,3] <- sigma[3]/n0[3]+4.0*sigma[4]/n0[4]+sigma[5]/n0[5]
    for(i in 2:3)
    {
      for(j in 1:(i-1))
      {
        hh8[i,j]<-hh8[j,i]
      }
    }
    #################################################
    b_line8[1] <- sumwx[1]/n0[1]-2.0*sumwx[2]/n0[2]+sumwx[3]/n0[3]
    b_line8[2] <- sumwx[2]/n0[2]-2.0*sumwx[3]/n0[3]+sumwx[4]/n0[4]
    b_line8[3] <- sumwx[3]/n0[3]-2.0*sumwx[4]/n0[4]+sumwx[5]/n0[5]
    B8 <- solve(hh8,b_line8)
    #################################################
    m[1] <- (sumwx[1]-sigma[1]*B8[1])/n0[1]; m[2] <- (sumwx[2]+sigma[2]*(2.0*B8[1]-B8[2]))/n0[2]
    m[3] <- (sumwx[3]-sigma[3]*(B8[1]-2.0*B8[2]+B8[3]))/n0[3]; m[4] <- (sumwx[4]-sigma[4]*(B8[2]-2.0*B8[3]))/n0[4]
    m[5] <- (sumwx[5]-sigma[5]*B8[3])/n0[5]
    ##########obtain variance##########
    for(i in 1:d5) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0 <- sum(swx)
    sigma <- matrix(c(swx/n0),5,1)
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- logL(m_sam,d5,mix_pi,m,sigma,data)
  AICm <- -2*abc + 6
  #########first order genetic parameter process##########
  hh81 <- matrix(c(1,1,1,1,1,2,1,0,-1,-2),5,2)
  B81 <- solve(t(hh81)%*%hh81)%*%(t(hh81)%*%m)
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0
  #########hypothesis testing##############
  data <- sort(data); bmw <- matrix(0,m_sam,1); bmwsl <- matrix(0,m_sam,d5)
  for(i in 1:d5){
    gg <- (data - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P6 <- rowSums(bmwsl)
  #############deal with ties in P6##############
  nn <- dim(as.matrix(unique(P6)))[1]
  if(nn < m_sam) {P6 <- P6 + runif(m_sam)/1e4}
  ###############################################
  dd <- as.matrix(c(sum(P6),sum(P6^2),sum((P6-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P6 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P6,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,m_sam))))

  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  output <- data.frame("2MG-EA",round(abc,4),round(AICm,4),round(t(m),4)," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4)," "," "," "," ",round(t(B81),4)," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4),round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4), tt[5])


  output<-as.matrix(output)

  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################2MG-CD model#########################################
F2ModelFun[[10]] <- function(K1,logL,df){
  data<-as.matrix(as.numeric(df[,1]))
  ##############procedure start####################
  m_sam <- dim(data)[1]; d6 <- 4;m_esp <- 0.0001
  mean0 <- mean(data);sigma0 <- (sd(data))^2
  #########2MG-CD Model#########
  a1 <- sqrt(sigma0/m_sam)
  sigma <- matrix((sigma0/2),4,1)
  m <- as.matrix(c((mean0+3*a1),(mean0+0.5*a1),(mean0-0.5*a1),(mean0-3*a1)))
  mi <- as.matrix(c(0.5625,0.1875,0.1875,0.0625))
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,4,m_sam); swx <- matrix(0,4,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d6) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[n0<0.000001] <- 0.000001
    ############obtain constraints############
    aa1 <- sumwx[1]/n0[1]-sumwx[2]/n0[2]-sumwx[3]/n0[3]+sumwx[4]/n0[4]
    aa2 <- sigma[1]/n0[1]+sigma[2]/n0[2]+sigma[3]/n0[3]+sigma[4]/n0[4]; aa3 <- aa1/aa2
    m[1] <- (sumwx[1]-sigma[1]*aa3)/n0[1]; m[2] <- (sumwx[2]+sigma[2]*aa3)/n0[2]
    m[3] <- (sumwx[3]+sigma[3]*aa3)/n0[3]; m[4] <- (sumwx[4]-sigma[4]*aa3)/n0[4]
    ##########obtain variance##########
    for(i in 1:d6) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0 <- sum(swx)
    sigma <- matrix(s0/m_sam,4,1)
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- logL(m_sam,d6,mix_pi,m,sigma,data)
  AICm <- -2*abc + 8
  #########first order genetic parameter process##########
  hh9 <- matrix(c(1,1,1,1,1,1,-1,-1,1,-1,1,-1),4,3)
  B9 <- solve(t(hh9)%*%hh9)%*%(t(hh9)%*%m)
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0
  #########hypothesis testing##############
  data <- sort(data); bmw <- matrix(0,m_sam,1); bmwsl <- matrix(0,m_sam,d6)
  for(i in 1:d6){
    gg <- (data - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P6 <- rowSums(bmwsl)
  #############deal with ties in P6##############
  nn <- dim(as.matrix(unique(P6)))[1]
  if(nn < m_sam) {P6 <- P6 + runif(m_sam)/1e4}
  ###############################################
  dd <- as.matrix(c(sum(P6),sum(P6^2),sum((P6-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P6 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P6,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,m_sam))))

  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  output <- data.frame("2MG-CD",round(abc,4),round(AICm,4),round(t(m),4)," "," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4)," "," "," "," "," ",round(t(B9),4)," "," "," "," "," "," ",round(jj,4),round(ll*100,4),round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4), tt[5])


  output<-as.matrix(output)

  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################2MG-EAD model#########################################
F2ModelFun[[11]] <- function(K1,logL,df){
  data <- as.matrix(as.numeric(df[,1]))
  ##############procedure start####################
  m_sam <- dim(data)[1]; d7 <- 3;m_esp <- 0.0001
  mean0 <- mean(data);sigma0 <- (sd(data))^2
  #########2MG-EAD Model#########
  a1 <- sqrt(sigma0/m_sam)
  sigma <- matrix((sigma0/2),3,1)
  m <- as.matrix(c((mean0+3*a1),mean0,(mean0-3*a1)))
  mi <- as.matrix(c(0.5625,0.375,0.0625))
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,3,m_sam); swx <- matrix(0,3,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d7) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[n0<0.000001] <- 0.000001
    ############obtain constraints############
    aa1 <- sumwx[1]/n0[1]-2.0*sumwx[2]/n0[2]+sumwx[3]/n0[3]
    aa2 <- sigma[1]/n0[1]+4.0*sigma[2]/n0[2]+sigma[3]/n0[3]; aa3 <- aa1/aa2
    m[1] <- (sumwx[1]-sigma[1]*aa3)/n0[1]; m[2] <- (sumwx[2]+sigma[2]*2.0*aa3)/n0[2]
    m[3] <- (sumwx[3]-sigma[3]*aa3)/n0[3]
    ##########obtain variance##########
    for(i in 1:d7) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0 <- sum(swx)
    sigma <- matrix(s0/m_sam,3,1)
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- logL(m_sam,d7,mix_pi,m,sigma,data)
  AICm <- -2*abc + 6
  #########first order genetic parameter process##########
  hh10 <- matrix(c(1,1,1,2,0,-2),3,2)
  B10 <- solve(t(hh10)%*%hh10)%*%(t(hh10)%*%m)
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0
  #########hypothesis testing##############
  data <- sort(data); bmw <- matrix(0,m_sam,1); bmwsl <- matrix(0,m_sam,d7)
  for(i in 1:d7){
    gg <- (data - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P6 <- rowSums(bmwsl)
  #############deal with ties in P6##############
  nn <- dim(as.matrix(unique(P6)))[1]
  if(nn < m_sam) {P6 <- P6 + runif(m_sam)/1e4}
  ###############################################
  dd <- as.matrix(c(sum(P6),sum(P6^2),sum((P6-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P6 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P6,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,m_sam))))

  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  output <- data.frame("2MG-EAD",round(abc,4),round(AICm,4),round(t(m),4)," "," "," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4)," "," "," "," "," "," ",round(t(B10),4)," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4),round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4), tt[5])


  output<-as.matrix(output)

  OUTPUT<-list(output,mi)
  return(OUTPUT)
}

K1F2 <- function(x){
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

logLF2 <- function(nm,nng,mi,mn,s,d1) { sum2 <- sum(log(dmixnorm(d1,mn,sqrt(s),mi)));return (sum2) }



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
   requireNamespace("kolmim")
    F2ModelFun[[i]](K1F2,logLF2,df)[[1]]
  }

  stopCluster(cl)
  mi<-NULL
}else{
 allresultq<-switch(model,"0MG" = F2ModelFun[[1]](K1F2,logLF2,df),"1MG-AD"=F2ModelFun[[2]](K1F2,logLF2,df),"1MG-A"=F2ModelFun[[3]](K1F2,logLF2,df),"1MG-EAD"=F2ModelFun[[4]](K1F2,logLF2,df),
                          "1MG-NCD"=F2ModelFun[[5]](K1F2,logLF2,df),"2MG-ADI"=F2ModelFun[[6]](K1F2,logLF2,df),"2MG-AD"=F2ModelFun[[7]](K1F2,logLF2,df),"2MG-A"=F2ModelFun[[8]](K1F2,logLF2,df),
                          "2MG-EA"=F2ModelFun[[9]](K1F2,logLF2,df),"2MG-CD"=F2ModelFun[[10]](K1F2,logLF2,df),"2MG-EAD"=F2ModelFun[[11]](K1F2,logLF2,df))

 allresult<-allresultq[[1]]
 if(model!="0MG"){
 mi<-allresultq[[2]]
 }else{
 mi<-NULL
 }
}
colnames(allresult) <- F2colname
out<-list(allresult,mi)
return(out)
}



































