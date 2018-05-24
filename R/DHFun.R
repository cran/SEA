DHFun<-function(df,model){
  
  
data<-sapply(df,as.character)
  
dDH<-data[-1,which(data[1,]=="DH")];DH<-as.numeric(dDH[which(is.na(as.numeric(dDH))==FALSE)]);df<-as.data.frame(DH) 
  

DHcolname <- c("Model","Log_Max_likelihood_Value","AIC","mean[1]","mean[2]","mean[3]","mean[4]","mean[5]","mean[6]","mean[7]","mean[8]","mean[9]","mean[10]","mean[11]",
               "mean[12]","mean[13]","mean[14]","mean[15]","mean[16]","Var(Residual+Polygene)","Proportion[1]","Proportion[2]","Proportion[3]","Proportion[4]","Proportion[5]",
               "Proportion[6]","Proportion[7]","Proportion[8]","Proportion[9]","Proportion[10]","Proportion[11]","Proportion[12]","Proportion[13]","Proportion[14]",
               "Proportion[15]","Proportion[16]","m","da","db","dc","dd","iab(i*)","iac","iad","ibc","ibd","icd","iabc","Major-Gene Var","Heritability(Major-Gene)(%)",	
               "U1 square","P(U1 square)","U2 square","P(U2 square)","U3 square","P(U3 square)","nW square","P(nW square)","Dn","P(Dn)")

DHModelFun<-list(NA)
###################define each model function##################
############################################0MG model#########################################  
DHModelFun[[1]] <- function(K1,logL,df){
  data <- as.matrix(as.numeric(df[,1]))
  ##############procedure start####################
  m_sam<-dim(data)[1];mean0<-mean(data);sigma0<-as.numeric(var(data))
  m_esp<-0.0001
  ###############################################################################
  ######## 0MG Model##########  (A0) ########
  mix_pi<-1.0;sigma<-sigma0;m<-mean0
  abc <-logL(m_sam,1,mix_pi,m,sigma,data) 
  AICm<- -2.0*abc+2.0*2.0   
  ##################hypothesis testing####################
  data<-sort(data);P2 <- matrix(0,m_sam,1)
  gg <- (data - m)/sqrt(as.vector(sigma))
  P2[which(gg>=0)] <- pnorm(gg[gg>=0])
  P2[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
  #############deal with ties in P1##############
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<m_sam){P2<-P2+runif(m_sam)/1e4}
  ###############################################
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P2 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)  
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,m_sam))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("0MG",round(abc,4),round(AICm,4),round(m,4)," "," "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(sigma,4),round(mix_pi,4),
                       " "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(u[1],4),
                       tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output)
  return(OUTPUT)
  
}
############################################1MG-A model######################################### 
DHModelFun[[2]] <- function(K1,logL,df){
  data <- as.matrix(as.numeric(df[,1]))
  ##############procedure start####################
  m_sam <- dim(data)[1];m_esp <- 0.0001
  mean0 <- mean(data);sigma0<-as.numeric(var(data))
  ######### 1MG-A Model########  (A1) #########
  d2 <- 2; a1 <- sqrt(sigma0/m_sam)
  mi <- as.matrix(c(0.5,0.5))
  sigma <- matrix((sigma0/2),d2,1)
  m <- as.matrix(c((mean0+1.5*a1),(mean0-1.5*a1))) 
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,m_sam); swx <- matrix(0,d2,1)
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi))) 
  while(stopa > m_esp && iteration<=1000){
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
    sigma <- matrix(s0/m_sam,d2,1)
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- logL(m_sam,d2,mix_pi,m,sigma,data)
  AICm <- -2*abc + 2*3
  
  #########first order genetic parameter process##########
  hh1 <- matrix(c(1,1,1,-1),2,2)
  b_line1 <- m; B1 <- solve(hh1,b_line1)
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
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<m_sam){P2<-P2+runif(m_sam)/1e4}
  
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P2 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,m_sam))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("1MG-A",round(abc,4),round(AICm,4),round(t(m),4)," "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," "," "," ", round(B1[1],4),round(B1[2],4)," "," "," "," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4),
                       round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################2MG-AI model######################################### 
DHModelFun[[3]] <- function(K1,logL,df){
  data <- as.matrix(as.numeric(df[,1]))
  ##############procedure start####################
  m_sam <- dim(data)[1];m_esp <- 0.0001
  mean0 <- mean(data);sigma0<-as.numeric(var(data))
  #############2MG-AI model#######  (B1) ######
  d2 <- 4;mi <- as.matrix(c(0.25,0.25,0.25,0.25))
  sigma <- matrix((sigma0/2),d2,1)
  a1 <- sqrt(sigma0/m_sam)
  m <- as.matrix(c((mean0+3*a1),(mean0+a1),(mean0-a1),(mean0-3*a1))) 
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,m_sam); swx <- matrix(0,d2,1)
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi))) 
  while(stopa > m_esp && iteration<=1000){
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
    sigma <- matrix(s0/m_sam,d2,1)
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- logL(m_sam,d2,mix_pi,m,sigma,data)
  AICm <- -2*abc + 2*5
  
  #########first order genetic parameter process##########
  hh1 <- matrix(c(1,1,1,1, 1,1,-1,-1, 1,-1,1,-1, 1,-1,-1,1),4,4)
  b_line1 <- m; B1 <- solve(hh1,b_line1)
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
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<m_sam){P2<-P2+runif(m_sam)/1e4}
  
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P2 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,m_sam))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("2MG-AI",round(abc,4),round(AICm,4),round(t(m),4)," "," "," "," "," "," "," "," "," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4),round(B1[2],4),round(B1[3],4)," "," ",round(B1[4],4)," "," "," "," "," "," ",round(jj,4),round(ll*100,4),
                       round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
} 
############################################2MG-A model######################################### 
DHModelFun[[4]] <- function(K1,logL,df){
  data <- as.matrix(as.numeric(df[,1]))
  ##############procedure start####################
  m_sam <- dim(data)[1];m_esp <- 0.0001
  mean0 <- mean(data);sigma0<-as.numeric(var(data))
  #############2MG-A model########  (B2) #####
  d2 <- 4; mi <- as.matrix(c(0.25,0.25,0.25,0.25))
  sigma <- matrix((sigma0/2),d2,1)
  a1 <- sqrt(sigma0/m_sam)
  m <- as.matrix(c((mean0+3*a1),(mean0+a1),(mean0-a1),(mean0-3*a1))) 
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,m_sam); swx <- matrix(0,d2,1)
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi))) 
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1   
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[n0<0.000001] <- 0.000001
    aa1<-sumwx[1]/n0[1]- sumwx[2]/n0[2]- sumwx[3]/n0[3]+ sumwx[4]/n0[4]
    aa2<-sum(sigma/n0)
    aa3<-aa1/aa2
    
    m[1]<-(sumwx[1]-sigma[1]*aa3)/n0[1]
    m[2]<-(sumwx[2]+sigma[2]*aa3)/n0[2]
    m[3]<-(sumwx[3]+sigma[3]*aa3)/n0[3]
    m[4]<-(sumwx[4]-sigma[4]*aa3)/n0[4]
    
    ##########obtain variance##########
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0 <- sum(swx)
    sigma <- matrix(s0/m_sam,d2,1)
    
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- logL(m_sam,d2,mix_pi,m,sigma,data)
  AICm <- -2*abc + 2*4
  
  #########first order genetic parameter process##########
  hh1 <- matrix(c(1,1,1,1, 1,1,-1,-1, 1,-1,1,-1),4,3)
  b_line1 <- m; B1 <- solve(t(hh1)%*%hh1)%*%(t(hh1)%*%m) 
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
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<m_sam){P2<-P2+runif(m_sam)/1e4}
  
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P2 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,m_sam))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("2MG-A",round(abc,4),round(AICm,4),round(t(m),4)," "," "," "," "," "," "," "," "," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4),round(B1[2],4),round(B1[3],4)," "," "," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4),
                       round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
} 
############################################2MG-EA model######################################### 
DHModelFun[[5]] <- function(K1,logL,df){
  data <- as.matrix(as.numeric(df[,1]))
  ##############procedure start####################
  m_sam <- dim(data)[1];m_esp <- 0.0001
  mean0 <- mean(data);sigma0<-as.numeric(var(data))
  #############2MG-EA model#######  (B3) ######
  d2 <- 3; mi <- as.matrix(c(0.25,0.5,0.25))
  sigma <- matrix((sigma0/2),d2,1)
  a1 <- sqrt(sigma0/m_sam)
  m <- as.matrix(c((mean0+2.5*a1),mean0,(mean0-2.5*a1))) 
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,m_sam); swx <- matrix(0,d2,1)
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi))) 
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1   
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[n0<0.000001] <- 0.000001
    aa1<-sumwx[1]/n0[1]-2*sumwx[2]/n0[2]+ sumwx[3]/n0[3]
    aa2<-sigma[1]/n0[1]+4*sigma[2]/n0[2]+ sigma[3]/n0[3]
    aa3<-aa1/aa2
    
    m[1]<-(sumwx[1]-sigma[1]*aa3)/n0[1]
    m[2]<-(sumwx[2]+2*sigma[2]*aa3)/n0[2]
    m[3]<-(sumwx[3]-sigma[3]*aa3)/n0[3]
    
    ##########obtain variance##########
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0 <- sum(swx)
    sigma <- matrix(s0/m_sam,d2,1)
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- logL(m_sam,d2,mix_pi,m,sigma,data)
  AICm <- -2*abc + 2*3
  
  #########first order genetic parameter process##########
  hh1 <- matrix(c(1,1,1,2,0,-2),3,2)
  b_line1 <- m;  B1 <- solve(t(hh1)%*%hh1)%*%(t(hh1)%*%m) 
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
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<m_sam){P2<-P2+runif(m_sam)/1e4}
  
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P2 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,m_sam))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("2MG-EA",round(abc,4),round(AICm,4),round(t(m),4)," "," "," "," "," "," "," "," "," "," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4),round(B1[2],4),round(B1[3],4)," "," "," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4),
                       round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################2MG-ED model######################################### 
DHModelFun[[6]] <- function(K1,logL,df){
  data <- as.matrix(as.numeric(df[,1]))
  ##############procedure start####################
  m_sam <- dim(data)[1];m_esp <- 0.0001
  mean0 <- mean(data);sigma0<-as.numeric(var(data))
  #############2MG-ED model######  (B4) #######
  d2 <- 3; mi <- as.matrix(c(0.5,0.25,0.25))
  sigma <- matrix((sigma0/2),d2,1)
  a1 <- sqrt(sigma0/m_sam)
  m <- as.matrix(c((mean0+2*a1),mean0,(mean0-2*a1))) 
  
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,m_sam); swx <- matrix(0,d2,1)
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi))) 
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1   
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[n0<0.000001] <- 0.000001
    m<- sumwx/n0
    
    ##########obtain variance##########
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0 <- sum(swx)
    sigma <- matrix(s0/m_sam,d2,1)
    
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- logL(m_sam,d2,mix_pi,m,sigma,data)
  AICm <- -2*abc + 2*4
  
  #########first order genetic parameter process##########
  hh1 <- matrix(c(1,1,1, 1,-1,-1, 0,1,-1),3,3)
  b_line1 <- m; B1 <- solve(hh1,b_line1)
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
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<m_sam){P2<-P2+runif(m_sam)/1e4}
  
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P2 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,m_sam))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("2MG-ED",round(abc,4),round(AICm,4),round(t(m),4)," "," "," "," "," "," "," "," "," "," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4),round(B1[2],4),round(B1[3],4)," "," "," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4),
                       round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
} 
############################################2MG-ER model######################################### 
DHModelFun[[7]] <- function(K1,logL,df){
  data <- as.matrix(as.numeric(df[,1]))
  ##############procedure start####################
  m_sam <- dim(data)[1];m_esp <- 0.0001
  mean0 <- mean(data);sigma0<-as.numeric(var(data))
  #############2MG-ER model#######  (B5) ######
  d2 <- 3; mi <- as.matrix(c(0.25,0.25,0.5))
  sigma <- matrix((sigma0/2),d2,1)
  a1 <- sqrt(sigma0/m_sam)
  m <- as.matrix(c((mean0+2*a1),mean0,(mean0-2*a1))) 
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,m_sam); swx <- matrix(0,d2,1)
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi))) 
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1   
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[n0<0.000001] <- 0.000001
    m<- sumwx/n0
    ##########obtain variance##########
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0 <- sum(swx)
    sigma <- matrix(s0/m_sam,d2,1)
    
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- logL(m_sam,d2,mix_pi,m,sigma,data)
  AICm <- -2*abc + 2*4
  
  #########first order genetic parameter process##########
  hh1 <- matrix(c(1,1,1, 1,1,-1, 1,-1,0),3,3)
  b_line1 <- m; B1 <- solve(hh1,b_line1)
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
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<m_sam){P2<-P2+runif(m_sam)/1e4}
  
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P2 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,m_sam))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("2MG-ER",round(abc,4),round(AICm,4),round(t(m),4)," "," "," "," "," "," "," "," "," "," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4),round(B1[2],4),round(B1[3],4)," "," "," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4),
                       round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
} 
############################################2MG-AE model######################################### 
DHModelFun[[8]] <- function(K1,logL,df){
  data <- as.matrix(as.numeric(df[,1]))
  ##############procedure start####################
  m_sam <- dim(data)[1];m_esp <- 0.0001
  mean0 <- mean(data);sigma0<-as.numeric(var(data))
  #############2MG-AE model######  (B6) #######
  d2 <- 3; mi <- as.matrix(c(0.25,0.5,0.25))
  sigma <- matrix((sigma0/2),d2,1)
  a1 <- sqrt(sigma0/m_sam)
  m <- as.matrix(c((mean0+2*a1),mean0,(mean0-2*a1))) 
  
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,m_sam); swx <- matrix(0,d2,1)
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi))) 
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1   
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[n0<0.000001] <- 0.000001
    m<-sumwx/n0
    
    ##########obtain variance##########
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0 <- sum(swx)
    sigma <- matrix(s0/m_sam,d2,1)
    
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- logL(m_sam,d2,mix_pi,m,sigma,data)
  AICm <- -2*abc + 2*4
  
  #########first order genetic parameter process##########
  hh1 <- matrix(c(1,1,1, 2,0,-2, 1,-1,1),3,3)
  b_line1 <- m; B1 <- solve(hh1,b_line1)
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
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<m_sam){P2<-P2+runif(m_sam)/1e4}
  
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P2 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,m_sam))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("2MG-AE",round(abc,4),round(AICm,4),round(t(m),4)," "," "," "," "," "," "," "," "," "," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4),round(B1[2],4)," "," "," ",round(B1[3],4)," "," "," "," "," "," ",round(jj,4),round(ll*100,4),
                       round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
} 
############################################2MG-CE model######################################### 
DHModelFun[[9]] <- function(K1,logL,df){
  data <- as.matrix(as.numeric(df[,1]))
  ##############procedure start####################
  m_sam <- dim(data)[1];m_esp <- 0.0001
  mean0 <- mean(data);sigma0<-as.numeric(var(data))
  #############2MG-CE model####  (B7) #########
  d2 <- 2; mi <- as.matrix(c(0.25,0.75))
  sigma <- matrix((sigma0/2),d2,1)
  a1 <- sqrt(sigma0/m_sam)
  m <- as.matrix(c((mean0+2*a1),(mean0-2*a1))) 
  
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,m_sam); swx <- matrix(0,d2,1)
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi))) 
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1   
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[n0<0.000001] <- 0.000001
    m<-sumwx/n0
    
    ##########obtain variance##########
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0 <- sum(swx)
    sigma <- matrix(s0/m_sam,d2,1)
    
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- logL(m_sam,d2,mix_pi,m,sigma,data)
  AICm <- -2*abc + 2*3
  
  #########first order genetic parameter process##########
  hh1 <- matrix(c(1,1,1,-1),2,2)
  b_line1 <- m; B1 <- solve(hh1,b_line1)
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
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<m_sam){P2<-P2+runif(m_sam)/1e4}
  
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P2 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,m_sam))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("2MG-CE",round(abc,4),round(AICm,4),round(t(m),4)," "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4)," "," "," "," ",round(B1[2],4)," "," "," "," "," "," ",round(jj,4),round(ll*100,4),
                       round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################2MG-DE model######################################### 
DHModelFun[[10]] <- function(K1,logL,df){
  data <- as.matrix(as.numeric(df[,1]))
  ##############procedure start####################
  m_sam <- dim(data)[1];m_esp <- 0.0001
  mean0 <- mean(data);sigma0<-as.numeric(var(data))
  #############2MG-DE model########  (B8) #####
  d2 <- 2; mi <- as.matrix(c(0.75,0.25))
  sigma <- matrix((sigma0/2),d2,1)
  a1 <- sqrt(sigma0/m_sam)
  m <- as.matrix(c((mean0+2*a1),(mean0-2*a1))) 
  
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,m_sam); swx <- matrix(0,d2,1)
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi))) 
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1   
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[n0<0.000001] <- 0.000001
    m<-sumwx/n0 
    
    ##########obtain variance##########
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0 <- sum(swx)
    sigma <- matrix(s0/m_sam,d2,1)
    
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- logL(m_sam,d2,mix_pi,m,sigma,data)
  AICm <- -2*abc + 2*3
  
  #########first order genetic parameter process##########
  hh1 <- matrix(c(1,1,1,-1),2,2)
  b_line1 <- m; B1 <- solve(hh1,b_line1)
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
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<m_sam){P2<-P2+runif(m_sam)/1e4}
  
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P2 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,m_sam))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("2MG-DE",round(abc,4),round(AICm,4),round(t(m),4)," "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4)," "," "," "," ",round(B1[2],4)," "," "," "," "," "," ",round(jj,4),round(ll*100,4),
                       round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
} 
############################################2MG-IE model######################################### 
DHModelFun[[11]] <- function(K1,logL,df){
  data <- as.matrix(as.numeric(df[,1]))
  ##############procedure start####################
  m_sam <- dim(data)[1];m_esp <- 0.0001
  mean0 <- mean(data);sigma0<-as.numeric(var(data))
  #############2MG-IE model########  (B9) #####
  d2 <- 2; mi <- as.matrix(c(0.75,0.25))
  sigma <- matrix((sigma0/2),d2,1)
  a1 <- sqrt(sigma0/m_sam)
  m <- as.matrix(c((mean0+2*a1),(mean0-2*a1))) 
  
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,m_sam); swx <- matrix(0,d2,1)
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi))) 
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1   
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[n0<0.000001] <- 0.000001
    m<-sumwx/n0 
    
    ##########obtain variance##########
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0 <- sum(swx)
    sigma <- matrix(s0/m_sam,d2,1)
    
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- logL(m_sam,d2,mix_pi,m,sigma,data)
  AICm <- -2*abc + 2*3
  
  #########first order genetic parameter process##########
  hh1 <- matrix(c(1,1,-1,1),2,2)
  b_line1 <- m; B1 <- solve(hh1,b_line1)
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
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<m_sam){P2<-P2+runif(m_sam)/1e4}
  
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P2 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,m_sam))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("2MG-IE",round(abc,4),round(AICm,4),round(t(m),4)," "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4)," "," "," "," ",round(B1[2],4)," "," "," "," "," "," ",round(jj,4),round(ll*100,4),
                       round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
} 
############################################3MG-AI model######################################### 
DHModelFun[[12]] <- function(K1,logL,df){
  data <- as.matrix(as.numeric(df[,1]))
  ##############procedure start####################
  m_sam <- dim(data)[1];m_esp <- 0.0001
  mean0 <- mean(data);sigma0<-as.numeric(var(data))
  #############3MG-AI model########  (F1) #####
  d2 <- 8; mi <- as.matrix(c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125))
  sigma <- matrix((sigma0/2),d2,1)
  a1 <- sqrt(sigma0/m_sam)
  m <- as.matrix(c((mean0+3*a1),(mean0+2.1*a1),(mean0+1.2*a1),(mean0+0.3*a1),(mean0+1.5*a1),(mean0+0.5*a1),(mean0-1.5*a1),(mean0-2.5*a1))) 
  
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,m_sam); swx <- matrix(0,d2,1)
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi))) 
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1   
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[n0<0.000001] <- 0.000001
    m<-sumwx/n0 
    
    ##########obtain variance##########
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0 <- sum(swx)
    sigma <- matrix(s0/m_sam,d2,1)
    
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- logL(m_sam,d2,mix_pi,m,sigma,data)
  AICm <- -2*abc + 2*9
  
  #########first order genetic parameter process##########
  hh1 <- matrix(c(1,1,1,1,1,1,1,1,
                  1,1,-1,-1,1,1,-1,-1, 
                  1,-1,1,-1,1,-1,1,-1, 
                  1,1,1,1,-1,-1,-1,-1, 
                  1,-1,-1,1,1,-1,-1,1, 
                  1,1,-1,-1,-1,-1,1,1, 
                  1,-1,1,-1,-1,1,-1,1, 
                  1,-1,-1,1,-1,1,1,-1),8,8)
  b_line1 <- m; B1 <- solve(hh1,b_line1)
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
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<m_sam){P2<-P2+runif(m_sam)/1e4}
  
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P2 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,m_sam))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("3MG-AI",round(abc,4),round(AICm,4),round(t(m),4)," "," "," "," "," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4), " "," "," "," "," "," "," "," ",
                       round(B1[1],4),round(B1[2],4),round(B1[3],4),round(B1[4],4)," ",round(B1[5],4),round(B1[6],4)," ",round(B1[7],4)," "," ",round(B1[8],4),round(jj,4),round(ll*100,4),
                       round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################3MG-A model######################################### 
DHModelFun[[13]] <- function(K1,logL,df){
  data <- as.matrix(as.numeric(df[,1]))
  ##############procedure start####################
  m_sam <- dim(data)[1];m_esp <- 0.0001
  mean0 <- mean(data);sigma0<-as.numeric(var(data))
  #############3MG-A model########  (F2) #####
  d2 <- 8; mi <- as.matrix(c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125))
  sigma <- matrix((sigma0/2),d2,1)
  a1 <- sqrt(sigma0/m_sam)
  m <- as.matrix(c((mean0+3*a1),(mean0+2.1*a1),(mean0+1.2*a1),(mean0+0.3*a1),(mean0+1.5*a1),(mean0+0.5*a1),(mean0-1.5*a1),(mean0-2.5*a1))) 
  
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,m_sam); swx <- matrix(0,d2,1)
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi))) 
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1   
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[n0<0.000001] <- 0.000001
    ############solve the linear equation############
    hh<-matrix(0,4,4)
    hh[1,1]<- sigma[1]/n0[1]+sigma[3]/n0[3]+sigma[6]/n0[6]+sigma[8]/n0[8]
    hh[2,2]<- sigma[2]/n0[2]+sigma[4]/n0[4]+sigma[5]/n0[5]+sigma[7]/n0[7]
    hh[3,3]<- sigma[1]/n0[1]+sigma[2]/n0[2]+sigma[5]/n0[5]+sigma[6]/n0[6]
    hh[4,4]<- sigma[3]/n0[3]+sigma[4]/n0[4]+sigma[7]/n0[7]+sigma[8]/n0[8]
    
    hh[1,2]<- hh[2,1]<-hh[3,4]<-hh[4,3]<-0
    hh[1,3]<- hh[3,1]<- sigma[1]/n0[1]-sigma[6]/n0[6]
    hh[1,4]<- hh[4,1]<- -sigma[3]/n0[3]+sigma[8]/n0[8]
    hh[2,3]<- hh[3,2]<- -sigma[2]/n0[2]+sigma[5]/n0[5]
    hh[2,4]<- hh[4,2]<- sigma[4]/n0[4]-sigma[7]/n0[7]
    ##################################################
    b_line<-matrix(0,4,1)
    b_line[1]<-sumwx[1]/n0[1]-sumwx[3]/n0[3]-sumwx[6]/n0[6]+sumwx[8]/n0[8]
    b_line[2]<-sumwx[2]/n0[2]-sumwx[4]/n0[4]-sumwx[5]/n0[5]+sumwx[7]/n0[7]
    b_line[3]<-sumwx[1]/n0[1]-sumwx[2]/n0[2]-sumwx[5]/n0[5]+sumwx[6]/n0[6]
    b_line[4]<-sumwx[3]/n0[3]-sumwx[4]/n0[4]-sumwx[7]/n0[7]+sumwx[8]/n0[8]
    B <- solve(hh,b_line) 
    ##################################################
    m[1]<-(sumwx[1]-sigma[1]*(B[1]+B[3]))/n0[1]
    m[2]<-(sumwx[2]-sigma[2]*(B[2]-B[3]))/n0[2]
    m[3]<-(sumwx[3]+sigma[3]*(B[1]-B[4]))/n0[3]
    m[4]<-(sumwx[4]+sigma[4]*(B[2]+B[4]))/n0[4]
    m[5]<-(sumwx[5]+sigma[5]*(B[2]+B[3]))/n0[5]
    m[6]<-(sumwx[6]+sigma[6]*(B[1]-B[3]))/n0[6]
    m[7]<-(sumwx[7]-sigma[7]*(B[2]-B[4]))/n0[7]
    m[8]<-(sumwx[8]-sigma[8]*(B[1]+B[4]))/n0[8]
    
    ##########obtain variance##########
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0 <- sum(swx)
    sigma <- matrix(s0/m_sam,d2,1)
    
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- logL(m_sam,d2,mix_pi,m,sigma,data)
  AICm <- -2*abc + 2*5
  
  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1,1, 1,1,-1,-1,1,1,-1,-1, 1,-1,1,-1,1,-1,1,-1, 1,1,1,1,-1,-1,-1,-1),8,4)
  b_line1 <- m;B1 <- ginv(t(aa)%*%aa)%*%(t(aa)%*%b_line1) 
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
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<m_sam){P2<-P2+runif(m_sam)/1e4}
  
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P2 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,m_sam))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("3MG-A",round(abc,4),round(AICm,4),round(t(m),4)," "," "," "," "," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4), " "," "," "," "," "," "," "," ",
                       round(B1[1],4),round(B1[2],4),round(B1[3],4),round(B1[4],4)," "," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4),
                       round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
} 
############################################3MG-CEA model######################################### 
DHModelFun[[14]] <- function(K1,logL,df){
  data <- as.matrix(as.numeric(df[,1]))
  ##############procedure start####################
  m_sam <- dim(data)[1];m_esp <- 0.0001
  mean0 <- mean(data);sigma0<-as.numeric(var(data))
  #############3MG-CEA model#######  (F3) ######
  d2 <- 4; mi <- as.matrix(c(0.125,0.375,0.375,0.125))
  sigma <- matrix((sigma0/2),d2,1)
  a1 <- sqrt(sigma0/m_sam)
  m <- as.matrix(c((mean0+3*a1),(mean0+a1),(mean0-a1),(mean0-3*a1))) 
  
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,m_sam); swx <- matrix(0,d2,1)
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi))) 
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1   
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[n0<0.000001] <- 0.000001
    ############solve the linear equation############
    hh<-matrix(0,2,2)
    hh[1,1]<- sigma[1]/n0[1]+sigma[2]/n0[2]+sigma[3]/n0[3]+sigma[4]/n0[4]
    hh[2,2]<- sigma[1]/n0[1]+9*sigma[2]/n0[2]+9*sigma[3]/n0[3]+sigma[4]/n0[4]
    
    hh[1,2]<- hh[2,1]<- sigma[1]/n0[1]+3*sigma[2]/n0[2]-3*sigma[3]/n0[3]-sigma[4]/n0[4]
    
    ##################################################
    b_line<-matrix(0,2,1)
    b_line[1]<- sumwx[1]/n0[1]-sumwx[2]/n0[2]-sumwx[3]/n0[3]+sumwx[4]/n0[4]
    b_line[2]<- sumwx[1]/n0[1]-3*sumwx[2]/n0[2]+3*sumwx[3]/n0[3]-sumwx[4]/n0[4]
    
    B <- solve(hh,b_line) 
    ##################################################
    m[1]<-(sumwx[1]-sigma[1]*(B[1]+B[2]))/n0[1]
    m[2]<-(sumwx[2]+sigma[2]*(B[1]+3*B[2]))/n0[2]
    m[3]<-(sumwx[3]+sigma[3]*(B[1]-3*B[2]))/n0[3]
    m[4]<-(sumwx[4]-sigma[4]*(B[1]-B[2]))/n0[4]
    
    ##########obtain variance##########
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0 <- sum(swx)
    sigma <- matrix(s0/m_sam,d2,1)
    
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- logL(m_sam,d2,mix_pi,m,sigma,data)
  AICm <- -2*abc + 2*3
  
  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1, 3,1,-1,-3),4,2)
  b_line1 <- m; B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1) 
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
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<m_sam){P2<-P2+runif(m_sam)/1e4}
  
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P2 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,m_sam))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("3MG-CEA",round(abc,4),round(AICm,4),round(t(m),4)," "," "," "," "," "," "," "," "," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4),round(B1[2],4),round(B1[3],4),round(B1[4],4)," "," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4),
                       round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
} 
############################################3MG-PEA model######################################### 
DHModelFun[[15]] <- function(K1,logL,df){
  data <- as.matrix(as.numeric(df[,1]))
  ##############procedure start####################
  m_sam <- dim(data)[1];m_esp <- 0.0001
  mean0 <- mean(data);sigma0<-as.numeric(var(data))
  #############3MG-PEA model########  (F4) #####
  d2 <- 6; mi <- as.matrix(c(0.125,0.125,0.25,0.25,0.125,0.125))
  sigma <- matrix((sigma0/2),d2,1)
  a1 <- sqrt(sigma0/m_sam)
  m <- as.matrix(c((mean0+3*a1),(mean0+2*a1),(mean0+a1),(mean0-a1),(mean0-2*a1),(mean0-3*a1))) 
  
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,m_sam); swx <- matrix(0,d2,1)
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi))) 
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1   
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[n0<0.000001] <- 0.000001
    ############solve the linear equation############
    hh<-matrix(0,3,3)
    hh[1,1]<- sigma[1]/n0[1]+sigma[2]/n0[2]+sigma[5]/n0[5]+sigma[6]/n0[6]
    hh[2,2]<- sigma[1]/n0[1]+sigma[2]/n0[2]+sigma[3]/n0[3]+sigma[4]/n0[4]
    hh[3,3]<- sigma[1]/n0[1]+4*sigma[3]/n0[3]+sigma[5]/n0[5]
    
    hh[1,2]<- hh[2,1]<- sigma[1]/n0[1]+sigma[2]/n0[2]
    hh[1,3]<- hh[3,1]<- sigma[1]/n0[1]-sigma[5]/n0[5]
    hh[2,3]<- hh[3,2]<- sigma[1]/n0[1]+2*sigma[3]/n0[3]
    
    ##################################################
    b_line<-matrix(0,3,1)
    b_line[1]<-sumwx[1]/n0[1]-sumwx[2]/n0[2]-sumwx[5]/n0[5]+sumwx[6]/n0[6]
    b_line[2]<-sumwx[1]/n0[1]-sumwx[2]/n0[2]-sumwx[3]/n0[3]+sumwx[4]/n0[4]
    b_line[3]<-sumwx[1]/n0[1]-2*sumwx[3]/n0[3]+sumwx[5]/n0[5]
    
    B <- solve(hh,b_line) 
    ##################################################
    m[1]<-(sumwx[1]-sigma[1]*(B[1]+B[2]+B[3]))/n0[1]
    m[2]<-(sumwx[2]+sigma[2]*(B[1]+B[2]))/n0[2]
    m[3]<-(sumwx[3]+sigma[3]*(B[2]+2*B[3]))/n0[3]
    m[4]<-(sumwx[4]-sigma[4]*B[2])/n0[4]
    m[5]<-(sumwx[5]+sigma[5]*(B[1]-B[3]))/n0[5]
    m[6]<-(sumwx[6]-sigma[6]*B[1])/n0[6]
    
    ##########obtain variance##########
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0 <- sum(swx)
    sigma <- matrix(s0/m_sam,d2,1)
    
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- logL(m_sam,d2,mix_pi,m,sigma,data)
  AICm <- -2*abc + 2*4
  
  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1, 2,2,0,0,-2,-2, 1,-1,1,-1,1,-1),6,3)
  b_line1 <- m;  B1 <- ginv(t(aa)%*%aa)%*%(t(aa)%*%b_line1) 
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
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<m_sam){P2<-P2+runif(m_sam)/1e4}
  
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P2 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,m_sam))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("3MG-PEA",round(abc,4),round(AICm,4),round(t(m),4)," "," "," "," "," "," "," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," ",round(B1[1],4),round(B1[2],4),round(B1[3],4),round(B1[4],4)," "," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4),
                       round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
} 
############################################4MG-AI model######################################### 
DHModelFun[[16]] <- function(K1,logL,df){
  data <- as.matrix(as.numeric(df[,1]))
  ##############procedure start####################
  m_sam <- dim(data)[1];m_esp <- 0.0001
  mean0 <- mean(data);sigma0<-as.numeric(var(data))
  #############4MG-AI model######  (H1) #######
  d2 <- 16; mi <- as.matrix(c(0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625))
  sigma <- matrix((sigma0/2),d2,1)
  a1 <- sqrt(sigma0/m_sam)
  m <- as.matrix(c((mean0+3*a1),(mean0+2.7*a1),(mean0+2.4*a1),(mean0+2.1*a1),(mean0+1.8*a1),(mean0+1.5*a1),(mean0+1.2*a1),(mean0+0.9*a1),
                   (mean0-0.9*a1),(mean0-1.2*a1),(mean0-1.5*a1),(mean0-1.8*a1),(mean0-2.1*a1),(mean0-2.4*a1),(mean0-2.7*a1),(mean0-3*a1))) 
  
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,m_sam); swx <- matrix(0,d2,1)
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi))) 
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1   
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[n0<0.000001] <- 0.000001
    ############solve the linear equation############
    hh<-matrix(0,5,5)
    hh[1,1]<- sigma[1]/n0[1]+sigma[2]/n0[2]+sigma[3]/n0[3]+sigma[4]/n0[4]+sigma[9]/n0[9]+sigma[10]/n0[10]+sigma[11]/n0[11]+sigma[12]/n0[12]
    hh[1,2]<- sigma[1]/n0[1]+sigma[2]/n0[2]+sigma[3]/n0[3]+sigma[4]/n0[4]
    hh[1,3]<- hh[1,2]
    hh[1,4]<- -(sigma[1]/n0[1]+sigma[2]/n0[2]+sigma[9]/n0[9]+sigma[10]/n0[10])
    hh[1,5]<- sigma[2]/n0[2]+sigma[3]/n0[3]+sigma[10]/n0[10]+sigma[11]/n0[11]
    
    hh[2,2]<- sigma[1]/n0[1]+sigma[2]/n0[2]+sigma[3]/n0[3]+sigma[4]/n0[4]+sigma[13]/n0[13]+sigma[14]/n0[14]+sigma[15]/n0[15]+sigma[16]/n0[16]
    hh[2,3]<- hh[1,2]
    hh[2,4]<- -(sigma[1]/n0[1]+sigma[2]/n0[2]-sigma[13]/n0[13]-sigma[14]/n0[14])
    hh[2,5]<- sigma[2]/n0[2]+sigma[3]/n0[3]+sigma[13]/n0[13]+sigma[16]/n0[16]
    
    hh[3,3]<- sigma[1]/n0[1]+sigma[2]/n0[2]+sigma[3]/n0[3]+sigma[4]/n0[4]+sigma[5]/n0[5]+sigma[6]/n0[6]+sigma[7]/n0[7]+sigma[8]/n0[8]
    hh[3,4]<- -(sigma[1]/n0[1]+sigma[2]/n0[2]+sigma[5]/n0[5]+sigma[6]/n0[6])
    hh[3,5]<- sigma[2]/n0[2]+sigma[3]/n0[3]-sigma[5]/n0[5]-sigma[8]/n0[8]
    
    hh[4,4]<- sigma[1]/n0[1]+sigma[2]/n0[2]+sigma[5]/n0[5]+sigma[6]/n0[6]+sigma[9]/n0[9]+sigma[10]/n0[10]+sigma[13]/n0[13]+sigma[14]/n0[14]
    hh[4,5]<- -sigma[2]/n0[2]+sigma[5]/n0[5]-sigma[10]/n0[10]+sigma[13]/n0[13]
    
    hh[5,5]<- sigma[2]/n0[2]+sigma[3]/n0[3]+sigma[5]/n0[5]+sigma[8]/n0[8]+sigma[10]/n0[10]+sigma[11]/n0[11]+sigma[13]/n0[13]+sigma[16]/n0[16]
    for(i in 2:5){
      for(j in 1:(i-1)){
        hh[i,j]<- hh[j,i]
      }
    }
    ##################################################
    b_line<-matrix(0,5,1)
    b_line[1]<- -sumwx[1]/n0[1]+sumwx[2]/n0[2]-sumwx[3]/n0[3]+sumwx[4]/n0[4]+sumwx[9]/n0[9]-sumwx[10]/n0[10]+sumwx[11]/n0[11]-sumwx[12]/n0[12]
    b_line[2]<- -sumwx[1]/n0[1]+sumwx[2]/n0[2]-sumwx[3]/n0[3]+sumwx[4]/n0[4]+sumwx[13]/n0[13]-sumwx[14]/n0[14]+sumwx[15]/n0[15]-sumwx[16]/n0[16]
    b_line[3]<- -sumwx[1]/n0[1]+sumwx[2]/n0[2]-sumwx[3]/n0[3]+sumwx[4]/n0[4]+sumwx[5]/n0[5]-sumwx[6]/n0[6]+sumwx[7]/n0[7]-sumwx[8]/n0[8]
    b_line[4]<- sumwx[1]/n0[1]-sumwx[2]/n0[2]-sumwx[5]/n0[5]+sumwx[6]/n0[6]-sumwx[9]/n0[9]+sumwx[10]/n0[10]+sumwx[13]/n0[13]-sumwx[14]/n0[14]
    b_line[5]<- sumwx[2]/n0[2]-sumwx[3]/n0[3]-sumwx[5]/n0[5]+sumwx[8]/n0[8]-sumwx[10]/n0[10]+sumwx[11]/n0[11]+sumwx[13]/n0[13]-sumwx[16]/n0[16]
    B <- solve(hh,b_line)
    ##################################################
    m[1]<-(sumwx[1]+(B[1]+B[2]+B[3]-B[4])*sigma[1])/n0[1]
    m[2]<-(sumwx[2]-(B[1]+B[2]+B[3]-B[4]+B[5])*sigma[2])/n0[2]
    m[3]<-(sumwx[3]+(B[1]+B[2]+B[3]+B[5])*sigma[3])/n0[3]
    m[4]<-(sumwx[4]-(B[1]+B[2]+B[3])*sigma[4])/n0[4]
    m[5]<-(sumwx[5]+(-B[3]+B[4]+B[5])*sigma[5])/n0[5]
    m[6]<-(sumwx[6]+(B[3]-B[4])*sigma[6])/n0[6]
    m[7]<-(sumwx[7]+(-B[3])*sigma[7])/n0[7]
    m[8]<-(sumwx[8]+(B[3]-B[5])*sigma[8])/n0[8]
    
    m[9]<-(sumwx[9]-(B[1]-B[4])*sigma[9])/n0[9]
    m[10]<-(sumwx[10]+(B[1]-B[4]+B[5])*sigma[10])/n0[10]
    m[11]<-(sumwx[11]-(B[1]+B[5])*sigma[11])/n0[11]
    m[12]<-(sumwx[12]+B[1]*sigma[12])/n0[12]
    m[13]<-(sumwx[13]-(B[2]+B[4]+B[5])*sigma[13])/n0[13]
    m[14]<-(sumwx[14]+(B[2]+B[4])*sigma[14])/n0[14]
    m[15]<-(sumwx[15]-B[2]*sigma[15])/n0[15]
    m[16]<-(sumwx[16]+(B[2]+B[5])*sigma[16])/n0[16]
    
    ##########obtain variance##########
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0 <- sum(swx)
    sigma <- matrix(s0/m_sam,d2,1)
    
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- logL(m_sam,d2,mix_pi,m,sigma,data)
  AICm <- -2*abc + 2*20
  
  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,
                1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,
                1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,
                
                1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,
                1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,
                1,1,-1,-1,1,1,-1,-1,-1,-1,1,1,-1,-1,1,1,
                1,-1,-1,1,1,-1,-1,1,-1,1,1,-1,-1,1,1,-1,
                
                1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,
                1,-1,-1,1,-1,1,1,-1,1,-1,-1,1,-1,1,1,-1,
                1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1),16,11)
  
  b_line1 <- m; B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1) 
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
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<m_sam){P2<-P2+runif(m_sam)/1e4}
  
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P2 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,m_sam))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("4MG-AI",round(abc,4),round(AICm,4),round(t(m),4),round(sigma[1],4),round(t(mix_pi),4),round(B1[1],4),round(B1[2],4),round(B1[3],4),round(B1[4],4),
                       round(B1[5],4),round(B1[6],4),round(B1[7],4),round(B1[8],4),round(B1[9],4),round(B1[10],4),round(B1[11],4)," ",round(jj,4),round(ll*100,4),
                       round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
} 
############################################4MG-CEA model######################################### 
DHModelFun[[17]] <- function(K1,logL,df){
  data <- as.matrix(as.numeric(df[,1]))
  ##############procedure start####################
  m_sam <- dim(data)[1];m_esp <- 0.0001
  mean0 <- mean(data);sigma0<-as.numeric(var(data))
  #############4MG-CEA model########  (H3) #####
  d2 <- 5; mi <- as.matrix(c(0.0625,0.25,0.375,0.25,0.0625))
  sigma <- matrix((sigma0/2),d2,1)
  a1 <- sqrt(sigma0/m_sam)
  m <- as.matrix(c((mean0+3*a1),(mean0+2*a1),mean0,(mean0-2*a1),(mean0-3*a1))) 
  
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,m_sam); swx <- matrix(0,d2,1)
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi))) 
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1   
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[n0<0.000001] <- 0.000001
    ############solve the linear equation############
    hh<-matrix(0,3,3)
    hh[1,1]<-sigma[1]/n0[1]+4*sigma[2]/n0[2]+sigma[3]/n0[3]
    hh[2,2]<-4*sigma[1]/n0[1]+9*sigma[2]/n0[2]+sigma[4]/n0[4]
    hh[3,3]<-sigma[1]/n0[1]+sigma[2]/n0[2]+sigma[4]/n0[4]+sigma[5]/n0[5]
    
    hh[1,2]<- hh[2,1]<- 2*sigma[1]/n0[1]+6*sigma[2]/n0[2]
    hh[1,3]<- hh[3,1]<- sigma[1]/n0[1]+2*sigma[2]/n0[2]
    hh[2,3]<- hh[3,2]<- 2*sigma[1]/n0[1]+3*sigma[2]/n0[2]-sigma[4]/n0[4]
    
    ##################################################
    b_line<-matrix(0,3,1)
    b_line[1]<- -sumwx[1]/n0[1]+2*sumwx[2]/n0[2]-sumwx[3]/n0[3]
    b_line[2]<- -2*sumwx[1]/n0[1]+3*sumwx[2]/n0[2]-sumwx[4]/n0[4]
    b_line[3]<- -sumwx[1]/n0[1]+sumwx[2]/n0[2]+sumwx[4]/n0[4]-sumwx[5]/n0[5]
    
    B <- solve(hh,b_line)
    ##################################################
    m[1]<-(sumwx[1]+sigma[1]*(B[1]+2*B[2]+B[3]))/n0[1]
    m[2]<-(sumwx[2]+sigma[2]*(-2*B[1]-3*B[2]-B[3]))/n0[2]
    m[3]<-(sumwx[3]+sigma[3]*B[1])/n0[3]
    m[4]<-(sumwx[4]+sigma[4]*(B[2]-B[3]))/n0[4]
    m[5]<-(sumwx[5]+sigma[5]*B[3])/n0[5]
    
    ##########obtain variance##########
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0 <- sum(swx)
    sigma <- matrix(s0/m_sam,d2,1)
    
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- logL(m_sam,d2,mix_pi,m,sigma,data)
  AICm <- -2*abc + 2*3
  
  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1, 4,2,0,-2,-4),5,2)
  b_line1 <- m; B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1) 
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
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<m_sam){P2<-P2+runif(m_sam)/1e4}
  
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P2 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,m_sam))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("4MG-CEA",round(abc,4),round(AICm,4),round(t(m),4)," "," "," "," "," "," "," "," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4),round(B1[2],4),round(B1[3],4),round(B1[4],4),round(B1[5],4),
                       " "," "," "," "," "," "," ",round(jj,4),round(ll*100,4),round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],
                       round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
} 
############################################4MG-EEA model######################################### 
DHModelFun[[18]] <- function(K1,logL,df){
  data <- as.matrix(as.numeric(df[,1]))
  ##############procedure start####################
  m_sam <- dim(data)[1];m_esp <- 0.0001
  mean0 <- mean(data);sigma0<-as.numeric(var(data))
  #############4MG-EEA model########  (H4) #####
  d2 <- 9; mi <- as.matrix(c(0.0625,0.0625,0.125,0.125,0.25,0.125,0.125,0.0625,0.0625))
  sigma <- matrix((sigma0/2),d2,1)
  a1 <- sqrt(sigma0/m_sam)
  m <- as.matrix(c((mean0+3*a1),(mean0+2.5*a1),(mean0+2*a1),(mean0+1.5*a1),(mean0+a1),(mean0-1.5*a1),(mean0-2*a1),(mean0-2.5*a1),(mean0-3*a1))) 
  
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,m_sam); swx <- matrix(0,d2,1)
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi))) 
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1   
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[n0<0.000001] <- 0.000001
    ############solve the linear equation############
    hh<-matrix(0,6,6)
    hh[1,1]<-sigma[1]/n0[1]+sigma[2]/n0[2]+sigma[4]/n0[4]+sigma[6]/n0[6]
    hh[1,2]<- hh[2,1]<- sigma[1]/n0[1]
    hh[1,3]<- 0
    hh[1,4]<- hh[4,1]<- sigma[4]/n0[4]-sigma[6]/n0[6]
    hh[1,5]<- hh[5,1]<- -(sigma[1]/n0[1]+sigma[4]/n0[4])
    hh[1,6]<- hh[6,1]<- -(sigma[1]/n0[1]+2*sigma[4]/n0[4])
    
    hh[2,2]<-sigma[1]/n0[1]+sigma[3]/n0[3]+sigma[7]/n0[7]+sigma[8]/n0[8]
    hh[2,3]<- hh[3,2]<- 2*sigma[7]/n0[7]+sigma[8]/n0[8]
    hh[2,4]<- 0
    hh[2,5]<- hh[5,2]<- -(sigma[1]/n0[1]+sigma[3]/n0[3])
    hh[2,6]<- hh[6,2]<- -(sigma[1]/n0[1]+2*sigma[7]/n0[7])
    
    hh[3,3]<-4*sigma[7]/n0[7]+sigma[8]/n0[8]+sigma[9]/n0[9]
    hh[3,4]<- hh[3,5]<- 0
    hh[3,6]<- hh[6,3]<- -(4*sigma[7]/n0[7]+sigma[9]/n0[9])
    
    hh[4,4]<-sigma[4]/n0[4]+4*sigma[5]/n0[5]+sigma[6]/n0[6]
    hh[4,5]<- hh[5,4]<- -(sigma[4]/n0[4]+2*sigma[5]/n0[5])
    hh[4,6]<- hh[6,4]<- -2*sigma[4]/n0[4]
    
    hh[5,5]<-sigma[1]/n0[1]+sigma[3]/n0[3]+sigma[4]/n0[4]+sigma[5]/n0[5]
    hh[5,6]<- hh[6,5]<- sigma[1]/n0[1]+2*sigma[4]/n0[4]
    
    hh[6,6]<-sigma[1]/n0[1]+4*sigma[4]/n0[4]+4*sigma[7]/n0[7]+sigma[9]/n0[9]
    for(i in 2:6){
      for(j in 1:(i-1)){
        hh[i,j]<- hh[j,i]
      }
    }
    ##################################################
    b_line<-matrix(0,6,1)
    b_line[1]<- sumwx[1]/n0[1]-sumwx[2]/n0[2]-sumwx[4]/n0[4]+sumwx[6]/n0[6]
    b_line[2]<- sumwx[1]/n0[1]-sumwx[3]/n0[3]+sumwx[7]/n0[7]-sumwx[8]/n0[8]
    b_line[3]<- 2*sumwx[7]/n0[7]-sumwx[8]/n0[8]-sumwx[9]/n0[9]
    b_line[4]<- -sumwx[4]/n0[4]+2*sumwx[5]/n0[5]-sumwx[6]/n0[6]
    b_line[5]<- -sumwx[1]/n0[1]+sumwx[3]/n0[3]+sumwx[4]/n0[4]-sumwx[5]/n0[5]
    b_line[6]<- -sumwx[1]/n0[1]+2*sumwx[4]/n0[4]+sumwx[9]/n0[9]-2*sumwx[7]/n0[7]
    B <- solve(hh,b_line)
    ##################################################
    m[1]<-(sumwx[1]-sigma[1]*(B[1]+B[2]-B[5]-B[6]))/n0[1]
    m[2]<-(sumwx[2]+sigma[2]*B[1])/n0[2]
    m[3]<-(sumwx[3]+sigma[3]*(B[2]-B[5]))/n0[3]
    m[4]<-(sumwx[4]+sigma[4]*(B[1]+B[4]-B[5]-2*B[6]))/n0[4]
    m[5]<-(sumwx[5]+sigma[5]*(-2*B[4]+B[5]))/n0[5]
    m[6]<-(sumwx[6]-sigma[6]*(B[1]-B[4]))/n0[6]
    m[7]<-(sumwx[7]-sigma[7]*(B[2]+2*B[3]-2*B[6]))/n0[7]
    m[8]<-(sumwx[8]+sigma[8]*(B[2]+B[3]))/n0[8]
    m[9]<-(sumwx[9]+sigma[9]*(B[3]-B[6]))/n0[9]
    
    ##########obtain variance##########
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0 <- sum(swx)
    sigma <- matrix(s0/m_sam,d2,1)
    
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- logL(m_sam,d2,mix_pi,m,sigma,data)
  AICm <- -2*abc + 2*4
  
  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1,1,1, 2,2,2,0,0,0,-2,-2,-2, 2,-2,0,2,0,-2,0,2,-2),9,3)
  
  b_line1 <- m; B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1) 
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
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<m_sam){P2<-P2+runif(m_sam)/1e4}
  
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P2 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,m_sam))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("4MG-EEA",round(abc,4),round(AICm,4),round(t(m),4)," "," "," "," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," ",round(B1[1],4),round(B1[2],4),round(B1[3],4),round(B1[4],4),round(B1[5],4),
                       " "," "," "," "," "," "," ",round(jj,4),round(ll*100,4),round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],
                       round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
} 
############################################4MG-EEEA model######################################### 
DHModelFun[[19]] <- function(K1,logL,df){
  data <- as.matrix(as.numeric(df[,1]))
  ##############procedure start####################
  m_sam <- dim(data)[1];m_esp <- 0.0001
  mean0 <- mean(data);sigma0<-as.numeric(var(data))
  #############4MG-EEEA model########  (H5) #####
  d2 <- 8; mi <- as.matrix(c(0.0625,0.0625,0.1875,0.1875,0.1875,0.1875,0.0625,0.0625))
  sigma <- matrix((sigma0/2),d2,1)
  a1 <- sqrt(sigma0/m_sam)
  m <- as.matrix(c((mean0+3*a1),(mean0+2.5*a1),(mean0+2*a1),(mean0+1.5*a1),(mean0-1.5*a1),(mean0-2*a1),(mean0-2.5*a1),(mean0-3*a1))) 
  
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,m_sam); swx <- matrix(0,d2,1)
  L0 <- sum(log(dmixnorm(data,m,sqrt(sigma),mi))) 
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1   
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam)
    sumwx <- WW%*%data
    n0 <- m_sam*mix_pi
    n0[n0<0.000001] <- 0.000001
    ############solve the linear equation############
    hh<-matrix(0,5,5)
    hh[1,1]<- 4*sigma[1]/n0[1]+9*sigma[4]/n0[4]+sigma[8]/n0[8]
    hh[1,2]<- -(2*sigma[1]/n0[1]+3*sigma[4]/n0[4])
    hh[1,3]<- 0
    hh[1,4]<- -2*sigma[1]/n0[1]
    hh[1,5]<- 0
    
    hh[2,2]<- sigma[1]/n0[1]+sigma[2]/n0[2]+sigma[3]/n0[3]+sigma[4]/n0[4]
    hh[2,3]<- sigma[2]/n0[2]+2*sigma[3]/n0[3]
    hh[2,4]<- sigma[1]/n0[1]+sigma[2]/n0[2]
    hh[2,5]<- 2*sigma[2]/n0[2]+3*sigma[3]/n0[3]
    
    hh[3,3]<- sigma[2]/n0[2]+4*sigma[3]/n0[3]+sigma[5]/n0[5]
    hh[3,4]<- sigma[2]/n0[2]-sigma[5]/n0[5]
    hh[3,5]<- 2*sigma[2]/n0[2]+6*sigma[3]/n0[3]
    
    hh[4,4]<- sigma[1]/n0[1]+sigma[2]/n0[2]+sigma[5]/n0[5]+sigma[6]/n0[6]
    hh[4,5]<- 2*sigma[2]/n0[2]
    
    hh[5,5]<- 4*sigma[2]/n0[2]+9*sigma[3]/n0[3]+sigma[7]/n0[7]
    for(i in 2:5){
      for(j in 1:(i-1)){
        hh[i,j]<- hh[j,i]
      }
    }
    ##################################################
    b_line<-matrix(0,5,1)
    b_line[1]<- -2*sumwx[1]/n0[1]+3*sumwx[4]/n0[4]-sumwx[8]/n0[8]
    b_line[2]<- sumwx[1]/n0[1]-sumwx[2]/n0[2]+sumwx[3]/n0[3]-sumwx[4]/n0[4]
    b_line[3]<- -sumwx[2]/n0[2]+2*sumwx[3]/n0[3]-sumwx[5]/n0[5]
    b_line[4]<- sumwx[1]/n0[1]-sumwx[2]/n0[2]+sumwx[5]/n0[5]-sumwx[6]/n0[6]
    b_line[5]<- -2*sumwx[2]/n0[2]+3*sumwx[3]/n0[3]-sumwx[7]/n0[7]
    B <- solve(hh,b_line)
    ##################################################
    m[1]<-(sumwx[1]+sigma[1]*(2*B[1]-B[2]-B[4]))/n0[1]
    m[2]<-(sumwx[2]+sigma[2]*(B[2]+B[3]+B[4]+2*B[5]))/n0[2]
    m[3]<-(sumwx[3]-sigma[3]*(B[2]+2*B[3]+3*B[5]))/n0[3]
    m[4]<-(sumwx[4]+sigma[4]*(-3*B[1]+B[2]))/n0[4]
    m[5]<-(sumwx[5]+sigma[5]*(B[3]-B[4]))/n0[5]
    m[6]<-(sumwx[6]+sigma[6]*B[4])/n0[6]
    m[7]<-(sumwx[7]+sigma[7]*B[5])/n0[7]
    m[8]<-(sumwx[8]+sigma[8]*B[1])/n0[8]
    
    ##########obtain variance##########
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(data-m[i])^2 }
    s0 <- sum(swx)
    sigma <- matrix(s0/m_sam,d2,1)
    
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(data,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- logL(m_sam,d2,mix_pi,m,sigma,data)
  AICm <- -2*abc + 2*4
  
  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1,1, 3,3,1,1,-1,-1,-3,-3, 1,-1,-1,1,-1,1,-1,1),8,3)
  b_line1 <- m; B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1) 
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
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<m_sam){P2<-P2+runif(m_sam)/1e4}
  
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*m_sam) + sum((P2 - (as.matrix(c(1:m_sam)) - 0.5)/m_sam)^2)
  u <- as.matrix(c(12*m_sam*((dd[1]/m_sam-0.5)^2),((45*m_sam)/4)*((dd[2]/m_sam-1/3)^2),180*m_sam*((dd[3]/m_sam-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,m_sam))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("4MG-EEEA",round(abc,4),round(AICm,4),round(t(m),4)," "," "," "," "," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," ",round(B1[1],4),round(B1[2],4),round(B1[3],4),round(B1[4],4),round(B1[5],4),
                       " "," "," "," "," "," "," ",round(jj,4),round(ll*100,4),round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],
                       round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
} 

K1DH <- function(x){
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

logLDH <- function(nm,nng,mi,mn,s,d1) { sum2 <- sum(log(dmixnorm(d1,mn,sqrt(s),mi)));return (sum2) } 


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
  allresult=foreach(i=1:19,.combine = 'rbind')%dopar%{
    requireNamespace("KScorrect")
    requireNamespace("kolmim")
    requireNamespace("MASS")
    DHModelFun[[i]](K1DH,logLDH,df)[[1]]
  }
  stopCluster(cl)
  mi<-NULL
}else{
  
allresultq<-switch(model,"0MG" = DHModelFun[[1]](K1DH,logLDH,df),"1MG-A"=DHModelFun[[2]](K1DH,logLDH,df),"2MG-AI"=DHModelFun[[3]](K1DH,logLDH,df),"2MG-A"=DHModelFun[[4]](K1DH,logLDH,df),
                   "2MG-EA"=DHModelFun[[5]](K1DH,logLDH,df),"2MG-ED"=DHModelFun[[6]](K1DH,logLDH,df),"2MG-ER"=DHModelFun[[7]](K1DH,logLDH,df),"2MG-AE"=DHModelFun[[8]](K1DH,logLDH,df),
                   "2MG-CE"=DHModelFun[[9]](K1DH,logLDH,df),"2MG-DE"=DHModelFun[[10]](K1DH,logLDH,df),"2MG-IE"=DHModelFun[[11]](K1DH,logLDH,df),"3MG-AI"=DHModelFun[[12]](K1DH,logLDH,df),
                   "3MG-A"=DHModelFun[[13]](K1DH,logLDH,df),"3MG-CEA"=DHModelFun[[14]](K1DH,logLDH,df),"3MG-PEA"=DHModelFun[[15]](K1DH,logLDH,df),"4MG-AI"=DHModelFun[[16]](K1DH,logLDH,df),
                   "4MG-CEA"=DHModelFun[[17]](K1DH,logLDH,df),"4MG-EEA"=DHModelFun[[18]](K1DH,logLDH,df),"4MG-EEEA"=DHModelFun[[19]](K1DH,logLDH,df))
  
allresult<-allresultq[[1]]
if(model!="0MG"){
  mi<-allresultq[[2]]  
}else{
  mi<-NULL
}
} 
colnames(allresult) <- DHcolname 
out<-list(allresult,mi)
return(out)
}







