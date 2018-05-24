BILFun<-function(df,model,BILfr){
 
data<-sapply(df,as.character)
  
dBIL<-data[-1,which(data[1,]=="BIL")];BIL<-as.numeric(dBIL[which(is.na(as.numeric(dBIL))==FALSE)]);df<-as.data.frame(BIL) 
  

BILcolname <- c("Model","Log_Max_likelihood_Value","AIC","mean[1]","mean[2]","mean[3]","mean[4]","mean[5]","mean[6]","mean[7]","mean[8]",
                "Var(Residual+Polygene)","Proportion[1]","Proportion[2]","Proportion[3]","Proportion[4]","Proportion[5]","Proportion[6]",
                "Proportion[7]","Proportion[8]","m","da(d)","db(d1)","dc(d2)","iab(i)","iac","ibc","iabc","i*","Major-Gene Var","Heritability(Major-Gene)(%)",	
                "U1 square","P(U1 square)","U2 square","P(U2 square)","U3 square","P(U3 square)","nW square","P(nW square)","Dn","P(Dn)")

BILModelFun<-list(NA)
###################define each model function##################
#######################0MG model############################## 
BILModelFun[[1]] <- function(K1,logL,df,BILfr){
  dataB <- as.matrix(as.numeric(df[,1]))
  mm <- dim(dataB)[1]; mean0 <-mean(dataB); sigma0<-var(dataB)
  ###############procedure start##############################
  m<-mean0
  sigma<-sigma0
  mix_pi<-1
  L0<-sum(log(dnorm(dataB,m,sqrt(sigma))))
  abc<-L0
  AIC<--2*abc+2*2
  ##############hypothesis testing#######################
  dataB<-sort(dataB)
  w1<-1/(12*mm)
  bmw <- matrix(0,mm,1)
  gg <- (dataB - m)/sqrt(as.vector(sigma))
  bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
  bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(bmw)))[1]
  if(nn < mm){bmw <- bmw+runif(mm)/1e4}
  #########################################################
  dd<-c((sum(bmw)),(sum(bmw^2)),sum((bmw-0.5)^2))
  w<-w1+sum((bmw - (as.matrix(c(1:mm)) - 0.5)/mm)^2)
  u<- as.matrix(c(12*mm*((dd[1]/mm-0.5)^2),((45*mm)/4)*((dd[2]/mm-1/3)^2),180*mm*((dd[3]/mm-1/12)^2)))
  D<-as.numeric(ks.test(bmw,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(w),(1-pkolm(D,mm))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)])
  
  output <- data.frame("0MG",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4)," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," ",
                       round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(w,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output)
  return(OUTPUT)
}
####################1MG-A(A-1)########################
BILModelFun[[2]] <- function(K1,logL,df,BILfr){
  dataB <- as.matrix(as.numeric(df[,1]))
  mm <- dim(dataB)[1];mean0 <-mean(dataB); sigma0<-var(dataB);m_esp <- 0.0001
  ###############procedure start##############################
  a1<-sqrt(sigma0)
  d2<-2
  if(BILfr=="BIL1(F1xP1)")
  {
    mi<-as.matrix(c(0.75,0.25))
    m<-as.matrix(c((mean0+0.25*a1),(mean0-0.75*a1))) 
  }else if(BILfr=="BIL2(F1xP2)")
  {
    mi<-as.matrix(c(0.25,0.75))
    m<-as.matrix(c((mean0-0.75*a1),(mean0+0.25*a1)))  
  }
  sigma<-matrix((sigma0/2),2,1)
  L0 <- sum(log(dmixnorm(dataB,m,sqrt(sigma),mi)))
  ##############iteration process###############
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,2,mm); swx <- matrix(0,2,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataB,m[i],sqrt(sigma[i]))/dmixnorm(dataB,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/mm)
    sumwx <- WW%*%dataB
    n0 <- mm*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    m <- sumwx/n0
    ##########obtain variance##########
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataB-m[i])^2 }
    s0 <- sum(swx)
    sigma[1]<-s0[1]/mm
    sigma[2]<-sigma[1]
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(dataB,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<- L0
  AIC<--2*abc+2*3
  ############genetic parameter ##########
  hh1<- matrix(c(1,1,1,-1),2,2)
  B1 <- solve(hh1,m)
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0
  #########hypothesis testing##############
  dataB <- sort(dataB); bmw <- matrix(0,mm,1); bmwsl <- matrix(0,mm,d2)
  for(i in 1:d2){
    gg <- (dataB - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  #############deal with ties in P2##############
  nn <- dim(as.matrix(unique(P2)))[1]
  if(nn < mm) {P2 <- P2 + runif(mm)/1e4}
  ###############################################
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*mm) + sum((P2 - (as.matrix(c(1:mm)) - 0.5)/mm)^2)
  u <- as.matrix(c(12*mm*((dd[1]/mm-0.5)^2),((45*mm)/4)*((dd[2]/mm-1/3)^2),180*mm*((dd[3]/mm-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,mm))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)])
  
  output <- data.frame("1MG-A",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4)," "," "," "," "," "," ",round(B1[1],4),round(B1[2],4)," "," "," "," "," "," "," ",round(jj*100,4),round(ll,4),
                       round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5]) 
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
########################2MG-AI(B-1)#################################
BILModelFun[[3]] <- function(K1,logL,df,BILfr){
  dataB <- as.matrix(as.numeric(df[,1]))
  mm <- dim(dataB)[1];mean0 <-mean(dataB); sigma0<-var(dataB);m_esp <- 0.0001
  ###############procedure start##############################  
  a1<-sqrt(sigma0)
  d4<-4
  if(BILfr=="BIL1(F1xP1)")
  {
    mi<-as.matrix(c(0.5625,0.1875,0.1875,0.0625))
    m<-as.matrix(c((mean0+a1),(mean0+0.25*a1),(mean0-0.5*a1),(mean0-1.5*a1))) 
  }else if(BILfr=="BIL2(F1xP2)")
  {
    mi<-as.matrix(c(0.0625,0.1875,0.1875,0.5625))
    m<-as.matrix(c((mean0-1.5*a1),(mean0+0.25*a1),(mean0-0.5*a1),(mean0+a1))) 
  }
  sigma<-matrix((sigma0/2),4,1)
  L0 <- sum(log(dmixnorm(dataB,m,sqrt(sigma),mi)))
  ##############iteration process###############
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,4,mm); swx <- matrix(0,4,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d4) { WW[i,] <- mi[i]*dnorm(dataB,m[i],sqrt(sigma[i]))/dmixnorm(dataB,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/mm)
    sumwx <- WW%*%dataB
    n0 <- mm*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    m<-sumwx/n0
    ###########obtain variance############
    for(i in 1:d4) {  swx[i] <- WW[i,]%*%(dataB-m[i])^2 }
    s0 <- sum(swx)
    sigma[1]<-s0/mm
    sigma[c(2,3,4)]<-sigma[1]
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(dataB,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<- L0
  AIC<--2*abc+2*5
  ############genetic parameter ##########
  hh2<- matrix(c(1,1,1,1,1,1,-1,-1,1,-1,1,-1,1,-1,-1,1),4,4)
  B2 <- solve(hh2,m)
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0
  #########hypothesis testing##############
  dataB <- sort(dataB); bmw <- matrix(0,mm,1); bmwsl <- matrix(0,mm,d4)
  for(i in 1:d4){
    gg <- (dataB - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  #############deal with ties in P2##############
  nn <- dim(as.matrix(unique(P2)))[1]
  if(nn < mm) {P2 <- P2 + runif(mm)/1e4}
  ###############################################
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*mm) + sum((P2 - (as.matrix(c(1:mm)) - 0.5)/mm)^2)
  u <- as.matrix(c(12*mm*((dd[1]/mm-0.5)^2),((45*mm)/4)*((dd[2]/mm-1/3)^2),180*mm*((dd[3]/mm-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,mm))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)])
  
  output <- data.frame("2MG-AI",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4)," "," "," "," ",round(B2[1],4),round(B2[2],4),round(B2[3],4)," ",round(B2[4],4)," "," "," "," ",round(jj*100,4),round(ll,4),
                       round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5]) 
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
#########################2MG-A(B-2)##########################
BILModelFun[[4]] <- function(K1,logL,df,BILfr){
  dataB <- as.matrix(as.numeric(df[,1]))
  mm <- dim(dataB)[1];mean0 <-mean(dataB); sigma0<-var(dataB);m_esp <- 0.0001
  ###############procedure start##############################  
  a1<-sqrt(sigma0)
  d4<-4
  if(BILfr=="BIL1(F1xP1)")
  {
    mi<-as.matrix(c(0.5625,0.1875,0.1875,0.0625))
    m<-as.matrix(c((mean0+a1),(mean0+0.25*a1),(mean0-0.5*a1),(mean0-1.5*a1)))
  }else if(BILfr=="BIL2(F1xP2)")
  {
    mi<-as.matrix(c(0.0625,0.1875,0.1875,0.5625))
    m<-as.matrix(c((mean0-1.5*a1),(mean0+0.25*a1),(mean0-0.5*a1),(mean0+a1)))
  }
  sigma<-matrix((sigma0/2),4,1)
  L0 <- sum(log(dmixnorm(dataB,m,sqrt(sigma),mi)))
  ##############iteration process###############
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,4,mm); swx <- matrix(0,4,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d4) { WW[i,] <- mi[i]*dnorm(dataB,m[i],sqrt(sigma[i]))/dmixnorm(dataB,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/mm)
    sumwx <- WW%*%dataB
    n0 <- mm*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    ###########restriction###########################
    aa1<-sumwx[1]/n0[1]-sumwx[2]/n0[2]-sumwx[3]/n0[3]+sumwx[4]/n0[4]
    aa2<-sigma[1]/n0[1]+sigma[2]/n0[2]+sigma[3]/n0[3]+sigma[4]/n0[4]
    aa1<-aa1/aa2
    m[1]<-(sumwx[1]-aa1*sigma[1])/n0[1]
    m[2]<-(sumwx[2]+aa1*sigma[2])/n0[2]
    m[3]<-(sumwx[3]+aa1*sigma[3])/n0[3]
    m[4]<-(sumwx[4]-aa1*sigma[4])/n0[4]
    #############obtain variance############
    for(i in 1:d4) {  swx[i] <- WW[i,]%*%(dataB-m[i])^2 }
    s0 <- sum(swx)
    sigma[1]<-s0/mm
    sigma[c(2,3,4)]<-sigma[1]
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(dataB,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<- L0
  AIC<--2*abc+2*4
  ############genetic parameter ##########
  hh3<- matrix(c(1,1,1,1,1,1,-1,-1,1,-1,1,-1),4,3)
  B3<-solve(crossprod(hh3,hh3))%*%crossprod(hh3,m)
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0  
  #########hypothesis testing##############
  dataB <- sort(dataB); bmw <- matrix(0,mm,1); bmwsl <- matrix(0,mm,d4)
  for(i in 1:d4){
    gg <- (dataB - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  #############deal with ties in P2##############
  nn <- dim(as.matrix(unique(P2)))[1]
  if(nn < mm) {P2 <- P2 + runif(mm)/1e4}
  ###############################################
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*mm) + sum((P2 - (as.matrix(c(1:mm)) - 0.5)/mm)^2)
  u <- as.matrix(c(12*mm*((dd[1]/mm-0.5)^2),((45*mm)/4)*((dd[2]/mm-1/3)^2),180*mm*((dd[3]/mm-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,mm))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)])
  
  output <- data.frame("2MG-A",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4)," "," "," "," ",round(B3[1],4),round(B3[2],4),round(B3[3],4)," "," "," "," "," "," ",round(jj*100,4),round(ll,4),
                       round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5]) 
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############2MG-EA(B-3)#############################
BILModelFun[[5]] <- function(K1,logL,df,BILfr){
  dataB <- as.matrix(as.numeric(df[,1]))
  mm <- dim(dataB)[1];mean0 <-mean(dataB); sigma0<-var(dataB);m_esp <- 0.0001
  ###############procedure start############################## 
  a1<-sqrt(sigma0)
  d3<-3
  if(BILfr=="BIL1(F1xP1)")
  {
    mi<-as.matrix(c(0.5625,0.375,0.0625))
    m<-as.matrix(c((mean0+a1),(mean0),(mean0-a1))) 
  }else if(BILfr=="BIL2(F1xP2)")
  {
    mi<-as.matrix(c(0.0625,0.375,0.5625))
    m<-as.matrix(c((mean0-a1),(mean0),(mean0+a1))) 
  }
  sigma<-matrix((sigma0/2),3,1)
  L0 <- sum(log(dmixnorm(dataB,m,sqrt(sigma),mi)))
  ##############iteration process###############
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,3,mm); swx <- matrix(0,3,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d3) { WW[i,] <- mi[i]*dnorm(dataB,m[i],sqrt(sigma[i]))/dmixnorm(dataB,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/mm)
    sumwx <- WW%*%dataB
    n0 <- mm*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001 
    #########restriction#############
    aa1<-sumwx[1]/n0[1]-2*sumwx[2]/n0[2]+sumwx[3]/n0[3]
    aa2<-sigma[1]/n0[1]+4*sigma[2]/n0[2]+sigma[3]/n0[3]
    aa1<-aa1/aa2
    m[1]<-(sumwx[1]-aa1*sigma[1])/n0[1]
    m[2]<-(sumwx[2]+2*aa1*sigma[2])/n0[2]
    m[3]<-(sumwx[3]+aa1*sigma[3])/n0[3]  
    #############obtain variance############
    for(i in 1:d3) {  swx[i] <- WW[i,]%*%(dataB-m[i])^2 }
    s0 <- sum(swx)
    sigma[1]<-s0/mm
    sigma[c(2,3)]<-sigma[1]
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(dataB,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
    if (iteration >=10)break  
  }
  abc<- L0
  AIC<--2*abc+2*3
  ############genetic parameter ##########
  hh4<- matrix(c(1,1,1,2,0,-2),3,2)
  B4<-solve(crossprod(hh4,hh4))%*%crossprod(hh4,m)
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0  
  #########hypothesis testing##############
  dataB <- sort(dataB); bmw <- matrix(0,mm,1); bmwsl <- matrix(0,mm,d3)
  for(i in 1:d3){
    gg <- (dataB - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  #############deal with ties in P2##############
  nn <- dim(as.matrix(unique(P2)))[1]
  if(nn < mm) {P2 <- P2 + runif(mm)/1e4}
  ###############################################
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*mm) + sum((P2 - (as.matrix(c(1:mm)) - 0.5)/mm)^2)
  u <- as.matrix(c(12*mm*((dd[1]/mm-0.5)^2),((45*mm)/4)*((dd[2]/mm-1/3)^2),180*mm*((dd[3]/mm-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,mm))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)])
  
  output <- data.frame("2MG-EA",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4)," "," "," "," "," ",round(B4[1],4),round(B4[2],4)," "," "," "," "," "," "," ",round(jj*100,4),round(ll,4),
                       round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5]) 
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
#############2MG-ED(B-4)###############
BILModelFun[[6]] <- function(K1,logL,df,BILfr){
  dataB <- as.matrix(as.numeric(df[,1]))
  mm <- dim(dataB)[1];mean0 <-mean(dataB); sigma0<-var(dataB);m_esp <- 0.0001
  ###############procedure start############################## 
  a1<-sqrt(sigma0)
  d3<-3
  if(BILfr=="BIL1(F1xP1)")
  {
    mi<-as.matrix(c(0.75,0.1875,0.0625))
    m<-as.matrix(c((mean0+0.5*a1),(mean0-0.25*a1),(mean0-a1)))
  }else if(BILfr=="BIL2(F1xP2)")
  {
    mi<-as.matrix(c(0.25,0.1875,0.5625))
    m<-as.matrix(c((mean0-0.25*a1),(mean0-a1),(mean0+0.5*a1))) 
  }
  sigma<-matrix((sigma0),3,1)
  L0 <- sum(log(dmixnorm(dataB,m,sqrt(sigma),mi)))
  ##############iteration process###############
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,3,mm); swx <- matrix(0,3,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d3) { WW[i,] <- mi[i]*dnorm(dataB,m[i],sqrt(sigma[i]))/dmixnorm(dataB,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/mm)
    sumwx <- WW%*%dataB
    n0 <- mm*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    m<-sumwx/n0
    #############obtain variance############
    for(i in 1:d3) {  swx[i] <- WW[i,]%*%(dataB-m[i])^2 }
    s0 <- sum(swx)
    sigma[1]<-s0/mm
    sigma[c(2,3)]<-sigma[1]
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(dataB,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<- L0
  AIC<--2*abc+2*4
  ############genetic parameter ##########
  hh5<- matrix(c(1,1,1,1,-1,-1,0,1,-1),3,3)
  B5 <- solve(hh5,m)
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0
  #########hypothesis testing##############
  dataB <- sort(dataB); bmw <- matrix(0,mm,1); bmwsl <- matrix(0,mm,d3)
  for(i in 1:d3){
    gg <- (dataB - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  #############deal with ties in P2##############
  nn <- dim(as.matrix(unique(P2)))[1]
  if(nn < mm) {P2 <- P2 + runif(mm)/1e4}
  ###############################################
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*mm) + sum((P2 - (as.matrix(c(1:mm)) - 0.5)/mm)^2)
  u <- as.matrix(c(12*mm*((dd[1]/mm-0.5)^2),((45*mm)/4)*((dd[2]/mm-1/3)^2),180*mm*((dd[3]/mm-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,mm))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)])
  
  output <- data.frame("2MG-ED",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4)," "," "," "," "," ",round(B5[1],4),round(B5[2],4),round(B5[3],4)," "," "," "," "," "," ",round(jj*100,4),round(ll,4),
                       round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5]) 
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
#############2MG-ER(B-5)##############
BILModelFun[[7]] <- function(K1,logL,df,BILfr){
  dataB <- as.matrix(as.numeric(df[,1]))
  mm <- dim(dataB)[1];mean0 <-mean(dataB); sigma0<-var(dataB);m_esp <- 0.0001
  ###############procedure start##############################  
  a1<-sqrt(sigma0)
  d3<-3
  if(BILfr=="BIL1(F1xP1)")
  {
    mi<-as.matrix(c(0.5625,0.1875,0.25))
    m<-as.matrix(c((mean0+0.5*a1),(mean0-0.25*a1),(mean0-a1))) 
  }else if(BILfr=="BIL2(F1xP2)")
  {
    mi<-as.matrix(c(0.0625,0.1875,0.75))
    m<-as.matrix(c((mean0-a1),(mean0-0.25*a1),(mean0+0.5*a1))) 
  }
  sigma<-matrix((sigma0/2),3,1)
  L0 <- sum(log(dmixnorm(dataB,m,sqrt(sigma),mi)))
  ##############iteration process###############
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,3,mm); swx <- matrix(0,3,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d3) { WW[i,] <- mi[i]*dnorm(dataB,m[i],sqrt(sigma[i]))/dmixnorm(dataB,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/mm)
    sumwx <- WW%*%dataB
    n0 <- mm*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    m<-sumwx/n0
    #############obtain variance############
    for(i in 1:d3) {  swx[i] <- WW[i,]%*%(dataB-m[i])^2 }
    s0 <- sum(swx)
    sigma[1]<-s0/mm
    sigma[c(2,3)]<-sigma[1]
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(dataB,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<- L0
  AIC<--2*abc+2*4
  ############genetic parameter ##########
  hh6<-matrix(c(1,1,1,1,1,-1,1,-1,0),3,3)
  B6<-solve(hh6,m)
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0
  #########hypothesis testing##############
  dataB <- sort(dataB); bmw <- matrix(0,mm,1); bmwsl <- matrix(0,mm,d3)
  for(i in 1:d3){
    gg <- (dataB - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  #############deal with ties in P2##############
  nn <- dim(as.matrix(unique(P2)))[1]
  if(nn < mm) {P2 <- P2 + runif(mm)/1e4}
  ###############################################
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*mm) + sum((P2 - (as.matrix(c(1:mm)) - 0.5)/mm)^2)
  u <- as.matrix(c(12*mm*((dd[1]/mm-0.5)^2),((45*mm)/4)*((dd[2]/mm-1/3)^2),180*mm*((dd[3]/mm-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,mm))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)])
  
  output <- data.frame("2MG-ER",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4)," "," "," "," "," ",round(B6[1],4),round(B6[2],4),round(B6[3],4)," "," "," "," "," "," ",round(jj*100,4),round(ll,4),
                       round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5]) 
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
#############2MG-AE(B-6)#####################
BILModelFun[[8]] <- function(K1,logL,df,BILfr){
  dataB <- as.matrix(as.numeric(df[,1]))
  mm <- dim(dataB)[1];mean0 <-mean(dataB); sigma0<-var(dataB);m_esp <- 0.0001
  ###############procedure start##############################
  a1<-sqrt(sigma0)
  d3<-3
  if(BILfr=="BIL1(F1xP1)")
  {
    mi<-as.matrix(c(0.5625,0.375,0.0625))
    m<-as.matrix(c((mean0+0.5*a1),(mean0-0.25*a1),(mean0-a1))) 
  }else if(BILfr=="BIL2(F1xP2)")
  {
    mi<-as.matrix(c(0.0625,0.375,0.5625))
    m<-as.matrix(c((mean0-a1),(mean0-0.25*a1),(mean0+0.5*a1)))  
  }
  sigma<-matrix((sigma0/2),3,1)
  L0 <- sum(log(dmixnorm(dataB,m,sqrt(sigma),mi)))
  ##############iteration process###############
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,3,mm); swx <- matrix(0,3,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d3) { WW[i,] <- mi[i]*dnorm(dataB,m[i],sqrt(sigma[i]))/dmixnorm(dataB,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/mm)
    sumwx <- WW%*%dataB
    n0 <- mm*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    m<-sumwx/n0
    #############obtain variance############
    for(i in 1:d3) {  swx[i] <- WW[i,]%*%(dataB-m[i])^2 }
    s0 <- sum(swx)
    sigma[1]<-s0/mm
    sigma[c(2,3)]<-sigma[1]
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(dataB,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<- L0
  AIC<--2*abc+2*4
  ############genetic parameter ##########
  hh7<-matrix(c(1,1,1,2,0,-2,1,-1,1),3,3)
  B7<-solve(hh7,m)
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0
  #########hypothesis testing##############
  dataB <- sort(dataB); bmw <- matrix(0,mm,1); bmwsl <- matrix(0,mm,d3)
  for(i in 1:d3){
    gg <- (dataB - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  #############deal with ties in P2##############
  nn <- dim(as.matrix(unique(P2)))[1]
  if(nn < mm) {P2 <- P2 + runif(mm)/1e4}
  ###############################################
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*mm) + sum((P2 - (as.matrix(c(1:mm)) - 0.5)/mm)^2)
  u <- as.matrix(c(12*mm*((dd[1]/mm-0.5)^2),((45*mm)/4)*((dd[2]/mm-1/3)^2),180*mm*((dd[3]/mm-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,mm))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)])
  
  output <- data.frame("2MG-AE",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4)," "," "," "," "," ",round(B7[1],4),round(B7[2],4)," "," ",round(B7[3],4)," "," "," "," ",round(jj*100,4),round(ll,4),
                       round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5]) 
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
###################2MG-CE(B-7)##########################################
BILModelFun[[9]] <- function(K1,logL,df,BILfr){
  dataB <- as.matrix(as.numeric(df[,1]))
  mm <- dim(dataB)[1];mean0 <-mean(dataB); sigma0<-var(dataB);m_esp <- 0.0001
  ###############procedure start############################## 
  a1<-sqrt(sigma0)
  d2<-2
  if(BILfr=="BIL1(F1xP1)")
  {
    mi<-as.matrix(c(0.5625,0.4375))
    m<-as.matrix(c((mean0+a1),(mean0-a1)))  
  }else if(BILfr=="BIL2(F1xP2)")
  {
    mi<-as.matrix(c(0.0625,0.9375))
    m<-as.matrix(c((mean0-a1),(mean0+a1))) 
  }
  sigma<-matrix((sigma0/2),2,1)
  L0 <- sum(log(dmixnorm(dataB,m,sqrt(sigma),mi)))
  ##############iteration process###############
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,2,mm); swx <- matrix(0,2,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataB,m[i],sqrt(sigma[i]))/dmixnorm(dataB,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/mm)
    sumwx <- WW%*%dataB
    n0 <- mm*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    m <- sumwx/n0
    ##########obtain variance##########
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataB-m[i])^2 }
    s0 <- sum(swx)
    sigma[1]<-s0/mm
    sigma[2]<-sigma[1]
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(dataB,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<- L0
  AIC<--2*abc+2*3
  ############genetic parameter ##########
  hh8<- matrix(c(1,1,1,-1),2,2)
  B8 <- solve(hh8,m)
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0
  #########hypothesis testing##############
  dataB <- sort(dataB); bmw <- matrix(0,mm,1); bmwsl <- matrix(0,mm,d2)
  for(i in 1:d2){
    gg <- (dataB - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  #############deal with ties in P2##############
  nn <- dim(as.matrix(unique(P2)))[1]
  if(nn < mm) {P2 <- P2 + runif(mm)/1e4}
  ###############################################
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*mm) + sum((P2 - (as.matrix(c(1:mm)) - 0.5)/mm)^2)
  u <- as.matrix(c(12*mm*((dd[1]/mm-0.5)^2),((45*mm)/4)*((dd[2]/mm-1/3)^2),180*mm*((dd[3]/mm-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,mm))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)])
  
  output <- data.frame("2MG-CE",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4)," "," "," "," "," "," ",round(B8[1],4)," "," "," "," "," "," "," ",round(B8[2],4),round(jj*100,4),round(ll,4),
                       round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5]) 
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############2MG-DE(B-8)##############
BILModelFun[[10]] <- function(K1,logL,df,BILfr){
  dataB <- as.matrix(as.numeric(df[,1]))
  mm <- dim(dataB)[1];mean0 <-mean(dataB); sigma0<-var(dataB);m_esp <- 0.0001
  ###############procedure start##############################  
  a1<-sqrt(sigma0)
  d2<-2
  if(BILfr=="BIL1(F1xP1)")
  {
    mi<-as.matrix(c(0.9375,0.0625))
    m<-as.matrix(c((mean0+0.5*a1),(mean0-0.5*a1))) 
  }else if(BILfr=="BIL2(F1xP2)")
  {
    mi<-as.matrix(c(0.4375,0.5625))
    m<-as.matrix(c((mean0-0.5*a1),(mean0+0.5*a1))) 
  }
  sigma<-matrix((sigma0/2),2,1)
  L0 <- sum(log(dmixnorm(dataB,m,sqrt(sigma),mi)))
  ##############iteration process###############
  iteration <- 0; stopa <- 1000
  m_esp<-0.0001
  WW <- matrix(0,2,mm); swx <- matrix(0,2,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataB,m[i],sqrt(sigma[i]))/dmixnorm(dataB,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/mm)
    sumwx <- WW%*%dataB
    n0 <- mm*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    m <- sumwx/n0
    ##########obtain variance##########
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataB-m[i])^2 }
    s0 <- sum(swx)
    sigma[1]<-s0/mm
    sigma[2]<-sigma[1]
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(dataB,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<- L0
  AIC<--2*abc+2*3
  ############genetic parameter ##########
  hh9<- matrix(c(1,1,1,-1),2,2)
  B9 <- solve(hh9,m)
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0
  #########hypothesis testing##############
  dataB <- sort(dataB); bmw <- matrix(0,mm,1); bmwsl <- matrix(0,mm,d2)
  for(i in 1:d2){
    gg <- (dataB - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  #############deal with ties in P2##############
  nn <- dim(as.matrix(unique(P2)))[1]
  if(nn < mm) {P2 <- P2 + runif(mm)/1e4}
  ###############################################
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*mm) + sum((P2 - (as.matrix(c(1:mm)) - 0.5)/mm)^2)
  u <- as.matrix(c(12*mm*((dd[1]/mm-0.5)^2),((45*mm)/4)*((dd[2]/mm-1/3)^2),180*mm*((dd[3]/mm-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,mm))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)])
  
  output <- data.frame("2MG-DE",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4)," "," "," "," "," "," ",round(B9[1],4)," "," "," "," "," "," "," ",round(B9[2],4),round(jj*100,4),round(ll,4),
                       round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5]) 
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
#############2MG-IE(B-9)##############
BILModelFun[[11]] <- function(K1,logL,df,BILfr){
  dataB <- as.matrix(as.numeric(df[,1]))
  mm <- dim(dataB)[1];mean0 <-mean(dataB); sigma0<-var(dataB);m_esp <- 0.0001
  ###############procedure start##############################
  a1<-sqrt(sigma0)
  d2<-2
  mi<-as.matrix(c(0.8125,0.1875))
  m<-as.matrix(c((mean0+a1),(mean0-a1)))
  sigma<-matrix((sigma0/2),2,1)
  L0 <- sum(log(dmixnorm(dataB,m,sqrt(sigma),mi)))
  ##############iteration process###############
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,2,mm); swx <- matrix(0,2,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataB,m[i],sqrt(sigma[i]))/dmixnorm(dataB,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/mm)
    sumwx <- WW%*%dataB
    n0 <- mm*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    m <- sumwx/n0
    ##########obtain variance##########
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataB-m[i])^2 }
    s0 <- sum(swx)
    sigma[1]<-s0/mm
    sigma[2]<-sigma[1]
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(dataB,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<- L0
  AIC<--2*abc+2*3
  ############genetic parameter ##########
  hh10<- matrix(c(1,1,-1,1),2,2)
  B10 <- solve(hh10,m)
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0
  #########hypothesis testing##############
  dataB <- sort(dataB); bmw <- matrix(0,mm,1); bmwsl <- matrix(0,mm,d2)
  for(i in 1:d2){
    gg <- (dataB - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  #############deal with ties in P2##############
  nn <- dim(as.matrix(unique(P2)))[1]
  if(nn < mm) {P2 <- P2 + runif(mm)/1e4}
  ###############################################
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*mm) + sum((P2 - (as.matrix(c(1:mm)) - 0.5)/mm)^2)
  u <- as.matrix(c(12*mm*((dd[1]/mm-0.5)^2),((45*mm)/4)*((dd[2]/mm-1/3)^2),180*mm*((dd[3]/mm-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,mm))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)])
  
  output <- data.frame("2MG-IE",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4)," "," "," "," "," "," ",round(B10[1],4)," "," "," "," "," "," "," ",round(B10[2],4),round(jj*100,4),round(ll,4),
                       round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5]) 
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
##############3MG-AI(F-1)##################
BILModelFun[[12]] <- function(K1,logL,df,BILfr){
  dataB <- as.matrix(as.numeric(df[,1]))
  mm <- dim(dataB)[1];mean0 <-mean(dataB); sigma0<-var(dataB);m_esp <- 0.0001
  ###############procedure start############################## 
  d8<-8
  a1<-sqrt(sigma0)
  if(BILfr=="BIL1(F1xP1)")
  {
    mi<-as.matrix(c(0.421875,0.140625,0.140625,0.046875,0.140625,0.046875,0.046875,0.015625))
    m<-as.matrix(c(mean0+a1,mean0+0.5*a1,mean0+0.25*a1,mean0,mean0-0.25*a1,mean0-0.5*a1,mean0-0.75*a1,mean0-a1))
  }else if(BILfr=="BIL2(F1xP2)")
  {
    mi<-as.matrix(c(0.015625,0.046875,0.046875,0.140625,0.046875,0.140625,0.140625,0.421875))
    m<-as.matrix(c(mean0-a1,mean0-0.75*a1,mean0-0.5*a1,mean0-0.25*a1,mean0,mean0+0.25*a1,mean0+0.5*a1,mean0+a1))
  }
  sigma<-matrix(c(sigma0/2),8,1)
  L0 <- sum(log(dmixnorm(dataB,m,sqrt(sigma),mi)))
  ##############iteration process###############
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,8,mm); swx <- matrix(0,8,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d8) { WW[i,] <- mi[i]*dnorm(dataB,m[i],sqrt(sigma[i]))/dmixnorm(dataB,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/mm)
    sumwx <- WW%*%dataB
    n0 <- mm*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    m<-sumwx/n0
    ########obtain variance#############
    for(i in 1:d8) {  swx[i] <- WW[i,]%*%(dataB-m[i])^2 }
    s0 <- sum(swx)
    sigma[1]<-s0/mm
    sigma[c(2:8)]<-sigma[1]
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(dataB,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<- L0
  AIC<--2*abc+2*9
  ############genetic parameter ##########
  hh11<- matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,-1,1,-1,
                  1,1,-1,-1,-1,-1,1,1,1,-1,1,-1,-1,1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,-1,1,1,-1),8,8)
  B11 <- solve(hh11,m)
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0
  #########hypothesis testing##############
  dataB <- sort(dataB); bmw <- matrix(0,mm,1); bmwsl <- matrix(0,mm,d8)
  for(i in 1:d8){
    gg <- (dataB - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  #############deal with ties in P2##############
  nn <- dim(as.matrix(unique(P2)))[1]
  if(nn < mm) {P2 <- P2 + runif(mm)/1e4}
  ###############################################
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*mm) + sum((P2 - (as.matrix(c(1:mm)) - 0.5)/mm)^2)
  u <- as.matrix(c(12*mm*((dd[1]/mm-0.5)^2),((45*mm)/4)*((dd[2]/mm-1/3)^2),180*mm*((dd[3]/mm-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,mm))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)])
  
  output <- data.frame("3MG-AI",round(abc,4),round(AIC,4),round(t(m),4),round(sigma[1],4),round(t(mix_pi),4),round(t(B11),4)," ",round(jj*100,4),round(ll,4),
                       round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5]) 
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
#####################3MG-A(F-2)###################
BILModelFun[[13]] <- function(K1,logL,df,BILfr){
  dataB <- as.matrix(as.numeric(df[,1]))
  mm <- dim(dataB)[1];mean0 <-mean(dataB); sigma0<-var(dataB);m_esp <- 0.0001
  ###############procedure start############################## 
  d8<-8
  a1<-sqrt(sigma0)
  if(BILfr=="BIL1(F1xP1)")
  {
    mi<-as.matrix(c(0.421875,0.140625,0.140625,0.046875,0.140625,0.046875,0.046875,0.015625))
    m<-as.matrix(c(mean0+a1,mean0+0.5*a1,mean0+0.25*a1,mean0,mean0-0.25*a1,mean0-0.5*a1,mean0-0.75*a1,mean0-a1))
    
  }else if(BILfr=="BIL2(F1xP2)")
  {
    mi<-as.matrix(c(0.015625,0.046875,0.046875,0.140625,0.046875,0.140625,0.140625,0.421875))
    m<-as.matrix(c(mean0-a1,mean0-0.75*a1,mean0-0.5*a1,mean0-0.25*a1,mean0,mean0+0.25*a1,mean0+0.5*a1,mean0+a1))
  }
  sigma<-matrix(c(sigma0/2),8,1)
  L0 <- sum(log(dmixnorm(dataB,m,sqrt(sigma),mi)))
  ##############iteration process###############
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,8,mm); swx <- matrix(0,8,1)
  hh<-matrix(0,4,4);b_line<-matrix(0,4,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d8) { WW[i,] <- mi[i]*dnorm(dataB,m[i],sqrt(sigma[i]))/dmixnorm(dataB,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/mm)
    sumwx <- WW%*%dataB
    n0 <- mm*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    #############restriction##################
    hh[1,1]<-sigma[1]/n0[1]+sigma[3]/n0[3]+sigma[6]/n0[6]+sigma[8]/n0[8]
    hh[1,2]<-0
    hh[1,3]<-sigma[1]/n0[1]-sigma[6]/n0[6]
    hh[1,4]<--sigma[3]/n0[3]+sigma[8]/n0[8]
    hh[2,2]<-sigma[2]/n0[2]+sigma[4]/n0[4]+sigma[5]/n0[5]+sigma[7]/n0[7]
    hh[2,3]<--sigma[2]/n0[2]+sigma[5]/n0[5]
    hh[2,4]<-sigma[4]/n0[4]-sigma[7]/n0[7] 
    hh[3,3]<-sigma[1]/n0[1]+sigma[2]/n0[2]+sigma[5]/n0[5]+sigma[6]/n0[6]
    hh[3,4]<-0
    hh[4,4]<-sigma[3]/n0[3]+sigma[4]/n0[4]+sigma[7]/n0[7]+sigma[8]/n0[8] 
    for(i in 2:4)
    {
      for(j in 1:(i-1))
      {
        hh[i,j]<-hh[j,i]
      }
    }
    ##############################################
    b_line[1]<-sumwx[1]/n0[1]-sumwx[3]/n0[3]-sumwx[6]/n0[6]+sumwx[8]/n0[8]
    b_line[2]<-sumwx[2]/n0[2]-sumwx[4]/n0[4]-sumwx[5]/n0[5]+sumwx[7]/n0[7]
    b_line[3]<-sumwx[1]/n0[1]-sumwx[2]/n0[2]-sumwx[5]/n0[5]+sumwx[6]/n0[6]
    b_line[4]<-sumwx[3]/n0[3]-sumwx[4]/n0[4]-sumwx[7]/n0[7]+sumwx[8]/n0[8]
    B121<-solve(hh,b_line)
    ############################################
    m[1]<-(sumwx[1]-sigma[1]*(B121[1]+B121[3]))/n0[1]
    m[2]<-(sumwx[2]-sigma[2]*(B121[2]-B121[3]))/n0[2]
    m[3]<-(sumwx[3]+sigma[3]*(B121[1]-B121[4]))/n0[3]
    m[4]<-(sumwx[4]+sigma[4]*(B121[2]+B121[4]))/n0[4]
    m[5]<-(sumwx[5]+sigma[5]*(B121[2]+B121[3]))/n0[5]
    m[6]<-(sumwx[6]+sigma[6]*(B121[1]-B121[3]))/n0[6]
    m[7]<-(sumwx[7]-sigma[7]*(B121[2]-B121[4]))/n0[7]
    m[8]<-(sumwx[8]-sigma[8]*(B121[1]+B121[4]))/n0[8]
    #########obtain variance######################
    for(i in 1:d8) {  swx[i] <- WW[i,]%*%(dataB-m[i])^2 }
    s0 <- sum(swx)
    sigma[1]<-s0/mm
    sigma[c(2:8)]<-sigma[1]
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(dataB,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<- L0
  AIC<--2*abc+2*5
  ############genetic parameter ##########
  hh12<- matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,-1,1,-1),8,4)
  B12<-solve(crossprod(hh12,hh12))%*%crossprod(hh12,m)
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0
  #########hypothesis testing##############
  dataB <- sort(dataB); bmw <- matrix(0,mm,1); bmwsl <- matrix(0,mm,d8)
  for(i in 1:d8){
    gg <- (dataB - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  #############deal with ties in P2##############
  nn <- dim(as.matrix(unique(P2)))[1]
  if(nn < mm) {P2 <- P2 + runif(mm)/1e4}
  ###############################################
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*mm) + sum((P2 - (as.matrix(c(1:mm)) - 0.5)/mm)^2)
  u <- as.matrix(c(12*mm*((dd[1]/mm-0.5)^2),((45*mm)/4)*((dd[2]/mm-1/3)^2),180*mm*((dd[3]/mm-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,mm))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)])
  
  output <- data.frame("3MG-A",round(abc,4),round(AIC,4),round(t(m),4),round(sigma[1],4),round(t(mix_pi),4),round(t(B12),4)," "," "," "," "," ",round(jj*100,4),round(ll,4),
                       round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5]) 
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############3MG-CEA(F-3)####################################
BILModelFun[[14]] <- function(K1,logL,df,BILfr){
  dataB <- as.matrix(as.numeric(df[,1]))
  mm <- dim(dataB)[1];mean0 <-mean(dataB); sigma0<-var(dataB);m_esp <- 0.0001
  ###############procedure start############################## 
  a1<-sqrt(sigma0)
  d4<-4
  if(BILfr=="BIL1(F1xP1)")
  {
    mi<-as.matrix(c(0.421875,0.421875,0.140625,0.015625))
    m<-as.matrix(c(mean0+0.75*a1,mean0-0.5*a1,mean0-a1,mean0-1.5*a1)) 
  }else if(BILfr=="BIL2(F1xP2)")
  {
    mi<-as.matrix(c(0.015625,0.140625,0.421875,0.421875))
    m<-as.matrix(c(mean0-1.5*a1,mean0-a1,mean0-0.5*a1,mean0+0.75*a1)) 
  }
  sigma<-matrix(c(sigma0/2),4,1)
  L0 <- sum(log(dmixnorm(dataB,m,sqrt(sigma),mi)))
  ##############iteration process###############
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,4,mm); swx <- matrix(0,4,1)
  hh<-matrix(0,2,2);b_line<-matrix(0,2,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d4) { WW[i,] <- mi[i]*dnorm(dataB,m[i],sqrt(sigma[i]))/dmixnorm(dataB,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/mm)
    sumwx <- WW%*%dataB
    n0 <- mm*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    #############restriction##################
    hh[1,1]<-sigma[1]/n0[1]+sigma[2]/n0[2]+sigma[3]/n0[3]+sigma[4]/n0[4]
    hh[1,2]<-sigma[1]/n0[1]+3*sigma[2]/n0[2]-3*sigma[3]/n0[3]-sigma[4]/n0[4]
    hh[2,1]<-hh[1,2]
    hh[2,2]<-sigma[1]/n0[1]+9*sigma[2]/n0[2]+9*sigma[3]/n0[3]+sigma[4]/n0[4]
    #####################################
    b_line[1]<-sumwx[1]/n0[1]-sumwx[2]/n0[2]-sumwx[3]/n0[3]+sumwx[4]/n0[4]
    b_line[2]<-sumwx[1]/n0[1]-3*sumwx[2]/n0[2]+3*sumwx[3]/n0[3]-sumwx[4]/n0[4]
    B131<-solve(hh,b_line)
    ####################################
    m[1]<-(sumwx[1]-sigma[1]*(B131[1]+B131[2]))/n0[1]
    m[2]<-(sumwx[2]+sigma[2]*(B131[1]+3*B131[2]))/n0[2]
    m[3]<-(sumwx[3]+sigma[3]*(B131[1]-3*B131[2]))/n0[3]
    m[4]<-(sumwx[4]-sigma[4]*(B131[1]-B131[2]))/n0[4]
    #########obtain variance######################
    for(i in 1:d4) {  swx[i] <- WW[i,]%*%(dataB-m[i])^2 }
    s0 <- sum(swx)
    sigma[1]<-s0/mm
    sigma[c(2:4)]<-sigma[1]
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(dataB,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<- L0
  AIC<--2*abc+2*4
  ############genetic parameter ##########
  hh13<- matrix(c(1,1,1,1,3,1,-1,-3),4,2)
  B13<-solve(crossprod(hh13,hh13))%*%crossprod(hh13,m)
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0
  #########hypothesis testing##############
  dataB <- sort(dataB); bmw <- matrix(0,mm,1); bmwsl <- matrix(0,mm,d4)
  for(i in 1:d4){
    gg <- (dataB - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  #############deal with ties in P2##############
  nn <- dim(as.matrix(unique(P2)))[1]
  if(nn < mm) {P2 <- P2 + runif(mm)/1e4}
  ###############################################
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*mm) + sum((P2 - (as.matrix(c(1:mm)) - 0.5)/mm)^2)
  u <- as.matrix(c(12*mm*((dd[1]/mm-0.5)^2),((45*mm)/4)*((dd[2]/mm-1/3)^2),180*mm*((dd[3]/mm-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,mm))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)])
  
  output <- data.frame("3MG-CEA",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," ",round(sigma[1],4),round(t(mix_pi),4)," "," "," "," ",round(t(B13),4)," "," "," "," "," "," "," ",round(jj*100,4),round(ll,4),
                       round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5]) 
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
#############3MG-PEA(F-4)###################
BILModelFun[[15]] <- function(K1,logL,df,BILfr){
  dataB <- as.matrix(as.numeric(df[,1]))
  mm <- dim(dataB)[1];mean0 <-mean(dataB); sigma0<-var(dataB);m_esp <- 0.0001
  ###############procedure start##############################
  d6<-6
  a1<-sqrt(sigma0)
  if(BILfr=="BIL1(F1xP1)")
  {
    mi<-as.matrix(c(0.421875,0.140625,0.28125,0.09375,0.046875,0.015625))
    m<-as.matrix(c(mean0+a1,mean0+0.5*a1,mean0,mean0-0.25*a1,mean0-0.50*a1,mean0-a1))  
  }else if(BILfr=="BIL2(F1xP2)")
  {
    mi<-as.matrix(c(0.015625,0.046875,0.09375,0.28125,0.140625,0.421875))
    m<-as.matrix(c(mean0-a1,mean0-0.5*a1,mean0-0.25*a1,mean0,mean0+0.5*a1,mean0+a1))  
  }
  sigma<-matrix(c(sigma0/2),6,1)
  L0 <- sum(log(dmixnorm(dataB,m,sqrt(sigma),mi)))
  ##############iteration process###############
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,6,mm); swx <- matrix(0,6,1)
  hh<-matrix(0,3,3);b_line<-matrix(0,3,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d6) { WW[i,] <- mi[i]*dnorm(dataB,m[i],sqrt(sigma[i]))/dmixnorm(dataB,m,sqrt(sigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/mm)
    sumwx <- WW%*%dataB
    n0 <- mm*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    #############restriction##################
    hh[1,1]<-sigma[1]/n0[1]+sigma[2]/n0[2]+sigma[5]/n0[5]+sigma[6]/n0[6]
    hh[1,2]<-sigma[1]/n0[1]+sigma[2]/n0[2]
    hh[1,3]<-sigma[1]/n0[1]-sigma[5]/n0[5]
    hh[2,2]<-sigma[1]/n0[1]+sigma[2]/n0[2]+sigma[3]/n0[3]+sigma[4]/n0[4]
    hh[2,3]<-sigma[1]/n0[1]+2*sigma[3]/n0[3]
    hh[3,3]<-sigma[1]/n0[1]+4*sigma[3]/n0[3]+sigma[5]/n0[5]
    for(i in 2:3)
    {
      for(j in 1:(i-1))
      {
        hh[i,j]<-hh[j,i]
      }
    }
    #####################################
    b_line[1]<-sumwx[1]/n0[1]-sumwx[2]/n0[2]-sumwx[5]/n0[5]+sumwx[6]/n0[6]
    b_line[2]<-sumwx[1]/n0[1]-sumwx[2]/n0[2]-sumwx[3]/n0[3]+sumwx[4]/n0[4]
    b_line[3]<-sumwx[1]/n0[1]-2*sumwx[3]/n0[3]+sumwx[5]/n0[5]
    B141<-solve(hh,b_line)
    ####################################
    m[1]<-(sumwx[1]-sigma[1]*(B141[1]+B141[2]+B141[3]))/n0[1]
    m[2]<-(sumwx[2]+sigma[2]*(B141[1]+B141[2]))/n0[2]
    m[3]<-(sumwx[3]+sigma[3]*(B141[2]+2.0*B141[3]))/n0[3]
    m[4]<-(sumwx[4]-sigma[4]*B141[2])/n0[4]
    m[5]<-(sumwx[5]+sigma[5]*(B141[1]-B141[3]))/n0[5]
    m[6]<-(sumwx[6]-sigma[6]*B141[1])/n0[6]
    #########obtain variance######################
    for(i in 1:d6) {  swx[i] <- WW[i,]%*%(dataB-m[i])^2 }
    s0 <- sum(swx)
    sigma[1]<-s0/mm
    sigma[c(2:6)]<-sigma[1]
    ########criteria for iterations to stop#######
    L1 <- sum(log(dmixnorm(dataB,m,sqrt(sigma),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<- L0
  AIC<--2*abc+2*3
  ############genetic parameter ##########
  hh14<- matrix(c(1,1,1,1,1,1,2,2,0,0,-2,-2,1,-1,1,-1,1,-1),6,3)
  B14<-solve(crossprod(hh14,hh14))%*%crossprod(hh14,m)
  jj <- sigma0 - sigma[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma0
  #########hypothesis testing##############
  dataB <- sort(dataB); bmw <- matrix(0,mm,1); bmwsl <- matrix(0,mm,d6)
  for(i in 1:d6){
    gg <- (dataB - m[i])/sqrt(sigma[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  #############deal with ties in P2##############
  nn <- dim(as.matrix(unique(P2)))[1]
  if(nn < mm) {P2 <- P2 + runif(mm)/1e4}
  ###############################################
  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*mm) + sum((P2 - (as.matrix(c(1:mm)) - 0.5)/mm)^2)
  u <- as.matrix(c(12*mm*((dd[1]/mm-0.5)^2),((45*mm)/4)*((dd[2]/mm-1/3)^2),180*mm*((dd[3]/mm-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,mm))))
  
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)])
  
  output <- data.frame("3MG-PEA",round(abc,4),round(AIC,4),round(t(m),4)," "," ",round(sigma[1],4),round(t(mix_pi),4)," "," ",round(B14[1],4)," ",round(B14[2],4),round(B14[3],4)," "," "," "," "," ",round(jj*100,4),round(ll,4),
                       round(u[1],4),tt[1],round(u[2],4),tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5]) 
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}


K1BIL <- function(x){
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

logLBIL <- function(nm,nng,mi,mn,s,d1) { sum2 <- sum(log(dmixnorm(d1,mn,sqrt(s),mi)));return (sum2) } 


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
  allresult=foreach(i=1:15,.combine = 'rbind')%dopar%{
    requireNamespace("KScorrect")
    requireNamespace("kolmim")
    BILModelFun[[i]](K1BIL,logLBIL,df,BILfr)[[1]]
  }
  stopCluster(cl)
  mi<-NULL

}else{
  
allresultq=switch(model,"0MG" = BILModelFun[[1]](K1BIL,logLBIL,df,BILfr),"1MG-A"=BILModelFun[[2]](K1BIL,logLBIL,df,BILfr),"2MG-AI"=BILModelFun[[3]](K1BIL,logLBIL,df,BILfr),"2MG-A"=BILModelFun[[4]](K1BIL,logLBIL,df,BILfr),"2MG-EA"=BILModelFun[[5]](K1BIL,logLBIL,df,BILfr),
               "2MG-ED"=BILModelFun[[6]](K1BIL,logLBIL,df,BILfr),"2MG-ER"=BILModelFun[[7]](K1BIL,logLBIL,df,BILfr),"2MG-AE"=BILModelFun[[8]](K1BIL,logLBIL,df,BILfr),"2MG-CE"=BILModelFun[[9]](K1BIL,logLBIL,df,BILfr),"2MG-DE"=BILModelFun[[10]](K1BIL,logLBIL,df,BILfr),
               "2MG-IE"=BILModelFun[[11]](K1BIL,logLBIL,df,BILfr),"3MG-AI"=BILModelFun[[12]](K1BIL,logLBIL,df,BILfr),"3MG-A"=BILModelFun[[13]](K1BIL,logLBIL,df,BILfr),"3MG-CEA"=BILModelFun[[14]](K1BIL,logLBIL,df,BILfr),"3MG-PEA"=BILModelFun[[15]](K1BIL,logLBIL,df,BILfr))
  
allresult<-allresultq[[1]]
if(model!="0MG"){
  mi<-allresultq[[2]]  
}else{
  mi<-NULL
} 
} 
colnames(allresult) <- BILcolname
out<-list(allresult,mi)
return(out) 
}


