G3DHFun<-function(df,model,G3DHtext2){

data<-sapply(df,as.character)

dP1<-data[-1,which(data[1,]=="P1")];P1<-as.numeric(dP1[which(is.na(as.numeric(dP1))==FALSE)]);df11<-as.data.frame(P1)
dP2<-data[-1,which(data[1,]=="P2")];P2<-as.numeric(dP2[which(is.na(as.numeric(dP2))==FALSE)]);df21<-as.data.frame(P2)
dDH<-data[-1,which(data[1,]=="DH")];DH<-as.numeric(dDH[which(is.na(as.numeric(dDH))==FALSE)]);df31<-as.data.frame(DH)

G3DHcolname <- c("Model","Log_Max_likelihood_Value","AIC","mean[P1]","mean[P2]","Var(P1 & P2)","mean[1]","mean[2]","mean[3]","mean[4]","mean[5]","mean[6]","mean[7]","mean[8]","mean[9]","mean[10]",
                 "mean[11]","mean[12]","mean[13]","mean[14]","mean[15]","mean[16]","Var(Residual+Polygene)","Proportion[1]","Proportion[2]","Proportion[3]","Proportion[4]","Proportion[5]",
                 "Proportion[6]","Proportion[7]","Proportion[8]","Proportion[9]","Proportion[10]","Proportion[11]","Proportion[12]","Proportion[13]","Proportion[14]","Proportion[15]","Proportion[16]",
                 "m(m1)","m2","m3","d(da)","db","dc","dd","iab(i*)","iac","iad","ibc","ibd","icd","iabc","[d]","Major-Gene Var","Heritability(Major-Gene)(%)","Polygenes Var","Heritability(Polygenes-Var)(%)",
                 "U1 square-P1","P(U1 square-P1)","U2 square-P1","P(U2 square-P1)","U3 square-P1","P(U3 square-P1)","nW square-P1","P(nW square-P1)","Dn-P1","P(Dn-P1)","U1 square-P2","P(U1 square-P2)","U2 square-P2","P(U2 square-P2)","U3 square-P2","P(U3 square-P2)","nW square-P2","P(nW square-P2)","Dn-P2","P(Dn-P2)",
                 "U1 square-DH","P(U1 square-DH)","U2 square-DH","P(U2 square-DH)","U3 square-DH","P(U3 square-DH)","nW square-DH","P(nW square-DH)","Dn-DH","P(Dn-DH)")

G3DHModelFun<-list(NA)
###################define each model function##################
##################### 0MG model##############################
G3DHModelFun[[1]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(A-0)Model##################################
  mi <- as.matrix(1); mix_pi <-1.0; meanA<-mean(dataDH); sigma <-2.0*sigma0; sigmaA <-sigma_dh/2
  abc <-logL(n_samP1,1,mix_pi,mean11,sigma,dataP1)+logL(n_samP2,1,mix_pi,mean12,sigma,dataP2)+logL(n_samDH,1,mix_pi,meanA,sigmaA,dataDH)
  AIC <- -2.0*abc+2.0*2.0
  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1)
  bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)
  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)
  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,1)

  gg <- (dataDH - meanA)/sqrt(as.vector(sigmaA))
  bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
  bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
  bmwsl[,1] <- bmw*mix_pi

  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("0MG",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanA),4)," "," "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(sigmaA,4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}

########################### 1MG-A model#########################################
G3DHModelFun[[2]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  ##############################(A-1)Model##########################################
  d2<-2 ; mi <- as.matrix(c(0.5,0.5)); meanA<-mean(dataDH); sigma <- 2.0*sigma0;
  sigmaA <- matrix((sigma_dh/2),d2,1)
  a1 <- sqrt(sigma_dh/n_samDH)
  if(mean11<mean12){a1=-a1}
  meanA <- as.matrix(c((meanA+1.5*a1),(meanA-1.5*a1)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanA,sqrt(sigmaA),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanA[i],sqrt(sigmaA[i]))/dmixnorm(dataDH,meanA,sqrt(sigmaA),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    mean11<-(sum(dataP1)+sumwx[1]*m_fam)/(n_samP1+n0[1]*m_fam)
    mean12<-(sum(dataP2)+sumwx[2]*m_fam)/(n_samP2+n0[2]*m_fam)
    meanA <- as.matrix(c(mean11,mean12))
    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanA[i])^2 }

    sigma<-(ss1+ss2+sum(swx)*m_fam)/(n_samP1+n_samP2+n_samDH)
    sigmaA<-matrix((sigma/m_fam),d2,1)
    ########criteria for iterations to stop#######

    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanA,sqrt(sigmaA),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*3
  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,-1),2,2)
  b_line1 <- meanA
  B1 <- solve(aa,b_line1)
  jj <- sigma_dh - sigmaA[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh

  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanA[i])/sqrt(sigmaA[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("1MG-A",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanA),4)," "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(sigmaA[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4)," "," ",round(B1[2],4)," "," "," "," "," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4)," "," ",round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################## 2MG-AI model##############################
G3DHModelFun[[3]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(B-1)Model##########################################
  d2<-4 ; mi <- as.matrix(c(0.25,0.25,0.25,0.25))
  meanB<-mean(dataDH); sigmaB <- matrix((sigma_dh/2),d2,1); sigma <- sigma0
  a1 <- sqrt(sigma_dh/n_samDH)
  if(mean11<mean12){a1=-a1}
  meanB <- as.matrix(c((meanB+3*a1),(meanB+a1),(meanB-a1),(meanB-3*a1)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanB,sqrt(sigmaB),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanB[i],sqrt(sigmaB[i]))/dmixnorm(dataDH,meanB,sqrt(sigmaB),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    mean11<-(sum(dataP1)+sumwx[1]*m_fam)/(n_samP1+n0[1]*m_fam)
    mean12<-(sum(dataP2)+sumwx[4]*m_fam)/(n_samP2+n0[4]*m_fam)
    meanB <- as.matrix(c(mean11,sumwx[2]/n0[2],sumwx[3]/n0[3],mean12))

    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanB[i])^2 }

    sigma<-(ss1+ss2+sum(swx)*m_fam)/(n_samP1+n_samP2+n_samDH)
    sigmaB<-matrix((sigma/m_fam),d2,1)

    ########criteria for iterations to stop#######
    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanB,sqrt(sigmaB),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*5

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1, 1,1,-1,-1, 1,-1,1,-1, 1,-1,-1,1),4,4)
  b_line1 <- meanB
  B1 <- solve(aa,b_line1)
  jj <- sigma_dh - sigmaB[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh

  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanB[i])/sqrt(sigmaB[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("2MG-AI",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanB),4)," "," "," "," "," "," "," "," "," "," "," "," ",round(sigmaB[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4)," "," ",round(B1[2],4),round(B1[3],4)," "," ",round(B1[4],4)," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4)," "," ",round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}

############################################ 2MG-A model#########################################
G3DHModelFun[[4]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(B-2)Model##########################################
  d2<-4 ; mi <- as.matrix(c(0.25,0.25,0.25,0.25))
  meanB<-mean(dataDH); sigmaB <- matrix((sigma_dh/2),d2,1); sigma <- sigma0
  a1 <- sqrt(sigma_dh/n_samDH)
  if(mean11<mean12){a1=-a1}
  meanB <- as.matrix(c((meanB+3*a1),(meanB+a1),(meanB-a1),(meanB-3*a1)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanB,sqrt(sigmaB),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanB[i],sqrt(sigmaB[i]))/dmixnorm(dataDH,meanB,sqrt(sigmaB),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    s0<-matrix(0,4,1)
    s0[1]<-sum(dataP1)+m_fam*sumwx[1]; s0[2]<-n_samP1+m_fam*n0[1]
    s0[3]<-sum(dataP2)+m_fam*sumwx[4]; s0[4]<-n_samP2+m_fam*n0[4]
    rr<-(s0[1]/s0[2]+s0[3]/s0[4]-sumwx[2]/n0[2]-sumwx[3]/n0[3])/(sigma/s0[2]+sigma/s0[4]+sigmaB[2]/n0[2]+sigmaB[3]/n0[3])

    mean11<-(s0[1]-rr*sigma)/s0[2]
    mean12<-(s0[3]-rr*sigma)/s0[4]
    meanB <- as.matrix(c(mean11,(sumwx[2]+sigmaB[2]*rr)/n0[2],(sumwx[3]+sigmaB[3]*rr)/n0[3],mean12))
    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanB[i])^2 }

    sigma<-(ss1+ss2+sum(swx)*m_fam)/(n_samP1+n_samP2+n_samDH)
    sigmaB<-matrix((sigma/m_fam),d2,1)

    ########criteria for iterations to stop#######
    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanB,sqrt(sigmaB),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*4

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1, 1,1,-1,-1, 1,-1,1,-1),4,3)
  b_line1 <- meanB
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)
  jj <- sigma_dh - sigmaB[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh

  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanB[i])/sqrt(sigmaB[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("2MG-A",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanB),4)," "," "," "," "," "," "," "," "," "," "," "," ",round(sigmaB[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4)," "," ",round(B1[2],4),round(B1[3],4)," "," "," "," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4)," "," ",round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################ 2MG-EA model#########################################
G3DHModelFun[[5]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(B-3)Model##########################################
  d2<-3
  mi <- as.matrix(c(0.25,0.5,0.25))
  meanB<-mean(dataDH)
  sigmaB <- matrix((sigma_dh/2),d2,1)
  sigma <- sigma0
  a1 <- sqrt(sigma_dh/n_samDH)
  if(mean11<mean12){a1=-a1}
  meanB <- as.matrix(c((meanB+2.5*a1),meanB,(meanB-2.5*a1)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanB,sqrt(sigmaB),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanB[i],sqrt(sigmaB[i]))/dmixnorm(dataDH,meanB,sqrt(sigmaB),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    s0<-matrix(0,4,1)
    s0[1]<-sum(dataP1)+m_fam*sumwx[1];s0[2]=n_samP1+m_fam*n0[1]
    s0[3]<-sum(dataP2)+m_fam*sumwx[3];s0[4]=n_samP2+m_fam*n0[3]
    rr<-(s0[1]/s0[2]-2.0*sumwx[2]/n0[2]+s0[3]/s0[4])/(sigma/s0[2]+sigma/s0[4]+4*sigmaB[2]/n0[2])

    mean11<-(s0[1]-rr*sigma)/s0[2]
    mean12<-(s0[3]-rr*sigma)/s0[4]
    meanB <- as.matrix(c(mean11,(sumwx[2]+2*sigmaB[2]*rr)/n0[2],mean12))
    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanB[i])^2 }

    sigma<-(ss1+ss2+sum(swx)*m_fam)/(n_samP1+n_samP2+n_samDH)
    sigmaB<-matrix((sigma/m_fam),d2,1)
    ########criteria for iterations to stop#######

    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanB,sqrt(sigmaB),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*3
  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,2,0,-2),3,2)
  b_line1 <- meanB
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)
  jj <- sigma_dh - sigmaB[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh

  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanB[i])/sqrt(sigmaB[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("2MG-EA",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanB),4)," "," "," "," "," "," "," "," "," "," "," "," "," ",round(sigmaB[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4)," "," ",round(B1[2],4),round(B1[2],4)," "," "," "," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4)," "," ",round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################ 2MG-ED model#########################################
G3DHModelFun[[6]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(B-4)Model##########################################
  d2<-3
  mi <- as.matrix(c(0.5,0.25,0.25))
  meanB<-mean(dataDH)
  sigmaB <- matrix((sigma_dh/2),d2,1)
  sigma <- sigma0
  a1 <- sqrt(sigma_dh/n_samDH)
  if(mean11<mean12){a1=-a1}

  meanB <- as.matrix(c((meanB+2*a1),meanB,(meanB-2*a1)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanB,sqrt(sigmaB),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanB[i],sqrt(sigmaB[i]))/dmixnorm(dataDH,meanB,sqrt(sigmaB),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    mean11<-(sum(dataP1)+m_fam*sumwx[1])/(n_samP1+m_fam*n0[1])
    mean12<-(sum(dataP2)+m_fam*sumwx[3])/(n_samP2+m_fam*n0[3])
    meanB <- as.matrix(c(mean11,sumwx[2]/n0[2],mean12))
    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanB[i])^2 }

    sigma<-(ss1+ss2+sum(swx)*m_fam)/(n_samP1+n_samP2+n_samDH)
    sigmaB<-matrix((sigma/m_fam),d2,1)

    ########criteria for iterations to stop#######

    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanB,sqrt(sigmaB),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*4
  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1, 1,-1,-1, 0,1,-1),3,3)
  b_line1 <- meanB
  B1 <- solve(aa,b_line1)
  jj <- sigma_dh - sigmaB[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh

  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanB[i])/sqrt(sigmaB[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("2MG-ED",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanB),4)," "," "," "," "," "," "," "," "," "," "," "," "," ",round(sigmaB[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4)," "," ",round(B1[2],4),round(B1[3],4)," "," "," "," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4)," "," ",round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################ 2MG-ER model#########################################
G3DHModelFun[[7]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(B-5)Model##########################################
  d2<-3
  mi <- as.matrix(c(0.25,0.25,0.5))
  meanB<-mean(dataDH)
  sigmaB <- matrix((sigma_dh/2),d2,1)
  sigma <- sigma0
  a1 <- sqrt(sigma_dh/n_samDH)
  if(mean11<mean12){a1=-a1}

  meanB <- as.matrix(c((meanB+2*a1),meanB,(meanB-2*a1)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanB,sqrt(sigmaB),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanB[i],sqrt(sigmaB[i]))/dmixnorm(dataDH,meanB,sqrt(sigmaB),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    mean11<-(sum(dataP1)+m_fam*sumwx[1])/(n_samP1+m_fam*n0[1])
    mean12<-(sum(dataP2)+m_fam*sumwx[3])/(n_samP2+m_fam*n0[3])
    meanB <- as.matrix(c(mean11,sumwx[2]/n0[2],mean12))

    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanB[i])^2 }

    sigma<-(ss1+ss2+sum(swx)*m_fam)/(n_samP1+n_samP2+n_samDH)
    sigmaB<-matrix((sigma/m_fam),d2,1)

    ########criteria for iterations to stop#######
    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanB,sqrt(sigmaB),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*4
  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1, 1,1,-1, 1,-1,0),3,3)
  b_line1 <- meanB
  B1 <- solve(aa,b_line1)
  jj <- sigma_dh - sigmaB[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh

  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanB[i])/sqrt(sigmaB[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("2MG-ER",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanB),4)," "," "," "," "," "," "," "," "," "," "," "," "," ",round(sigmaB[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4)," "," ",round(B1[2],4),round(B1[3],4)," "," "," "," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4)," "," ",round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################ 2MG-AE model#########################################
G3DHModelFun[[8]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(B-6)Model##########################################
  d2<-3
  mi <- as.matrix(c(0.25,0.5,0.25))
  meanB<-mean(dataDH)
  sigmaB <- matrix((sigma_dh/2),d2,1)
  sigma <- sigma0
  a1 <- sqrt(sigma_dh/n_samDH)
  if(mean11<mean12){a1=-a1}

  meanB <- as.matrix(c((meanB+2*a1),meanB,(meanB-2*a1)))

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanB,sqrt(sigmaB),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanB[i],sqrt(sigmaB[i]))/dmixnorm(dataDH,meanB,sqrt(sigmaB),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    mean11<-(sum(dataP1)+m_fam*sumwx[1])/(n_samP1+m_fam*n0[1])
    mean12<-(sum(dataP2)+m_fam*sumwx[3])/(n_samP2+m_fam*n0[3])
    meanB <- as.matrix(c(mean11,sumwx[2]/n0[2],mean12))

    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanB[i])^2 }

    sigma<-(ss1+ss2+sum(swx)*m_fam)/(n_samP1+n_samP2+n_samDH)
    sigmaB<-matrix((sigma/m_fam),d2,1)

    ########criteria for iterations to stop#######
    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanB,sqrt(sigmaB),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*4

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1, 2,0,-2, 1,-1,1),3,3)
  b_line1 <- meanB
  B1 <- solve(aa,b_line1)
  jj <- sigma_dh - sigmaB[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh

  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanB[i])/sqrt(sigmaB[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("2MG-AE",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanB),4)," "," "," "," "," "," "," "," "," "," "," "," "," ",round(sigmaB[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4)," "," ",round(B1[2],4)," "," "," ",round(B1[3],4)," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4)," "," ",round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################ 2MG-CE model#########################################
G3DHModelFun[[9]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(B-7)Model##########################################
  d2<-2
  mi <- as.matrix(c(0.25,0.75))
  meanB<-mean(dataDH)
  sigmaB <- matrix((sigma_dh/2),d2,1)
  sigma <- sigma0
  a1 <- sqrt(sigma_dh/n_samDH)
  if(mean11<mean12){a1=-a1}

  meanB <- as.matrix(c((meanB+2*a1),(meanB-2*a1)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanB,sqrt(sigmaB),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanB[i],sqrt(sigmaB[i]))/dmixnorm(dataDH,meanB,sqrt(sigmaB),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    mean11<-(sum(dataP1)+m_fam*sumwx[1])/(n_samP1+m_fam*n0[1])
    mean12<-(sum(dataP2)+m_fam*sumwx[2])/(n_samP2+m_fam*n0[2])
    meanB <- as.matrix(c(mean11,mean12))
    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanB[i])^2 }

    sigma<-(ss1+ss2+sum(swx)*m_fam)/(n_samP1+n_samP2+n_samDH)
    sigmaB<-matrix((sigma/m_fam),d2,1)

    ########criteria for iterations to stop#######

    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanB,sqrt(sigmaB),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*3

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,-1),2,2)
  b_line1 <- meanB
  B1 <- solve(aa,b_line1)
  jj <- sigma_dh - sigmaB[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh

  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanB[i])/sqrt(sigmaB[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("2MG-CE",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanB),4)," "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(sigmaB[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4)," "," "," "," "," "," ",round(B1[2],4)," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4)," "," ",round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################ 2MG-DE#########################################
G3DHModelFun[[10]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(B-8)Model##########################################
  d2<-2
  mi <- as.matrix(c(0.75,0.25))
  meanB<-mean(dataDH)
  sigmaB <- matrix((sigma_dh/2),d2,1)
  sigma <- sigma0
  a1 <- sqrt(sigma_dh/n_samDH)
  if(mean11<mean12){a1=-a1}

  meanB <- as.matrix(c((meanB+2*a1),(meanB-2*a1)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanB,sqrt(sigmaB),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanB[i],sqrt(sigmaB[i]))/dmixnorm(dataDH,meanB,sqrt(sigmaB),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    mean11<-(sum(dataP1)+m_fam*sumwx[1])/(n_samP1+m_fam*n0[1])
    mean12<-(sum(dataP2)+m_fam*sumwx[2])/(n_samP2+m_fam*n0[2])
    meanB <- as.matrix(c(mean11,mean12))

    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanB[i])^2 }

    sigma<-(ss1+ss2+sum(swx)*m_fam)/(n_samP1+n_samP2+n_samDH)
    sigmaB<-matrix((sigma/m_fam),d2,1)

    ########criteria for iterations to stop#######

    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanB,sqrt(sigmaB),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*3

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,-1),2,2)
  b_line1 <- meanB
  B1 <- solve(aa,b_line1)
  jj <- sigma_dh - sigmaB[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh

  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanB[i])/sqrt(sigmaB[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("2MG-DE",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanB),4)," "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(sigmaB[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4)," "," "," "," "," "," ",round(B1[2],4)," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4)," "," ",round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################ 2MG-IE#########################################
G3DHModelFun[[11]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(B-9)Model##########################################
  d2<-2
  mi <- as.matrix(c(0.75,0.25))
  meanB<-mean(dataDH)
  sigmaB <- matrix((sigma_dh/2),d2,1)
  sigma <- sigma0
  a1 <- sqrt(sigma_dh/n_samDH)
  if(mean11<mean12){a1=-a1}

  meanB <- as.matrix(c((meanB+2*a1),(meanB-2*a1)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanB,sqrt(sigmaB),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanB[i],sqrt(sigmaB[i]))/dmixnorm(dataDH,meanB,sqrt(sigmaB),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    mean11<-(sum(dataP1)+sum(dataP2)+m_fam*sumwx[1])/(n_samP1+n_samP2+m_fam*n0[1])
    mean12<- mean11
    meanB <- as.matrix(c(mean11,sumwx[2]/n0[2]))

    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanB[i])^2 }

    sigma<-(ss1+ss2+sum(swx)*m_fam)/(n_samP1+n_samP2+n_samDH)
    sigmaB<-matrix((sigma/m_fam),d2,1)

    ########criteria for iterations to stop#######
    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanB,sqrt(sigmaB),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*3

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,-1,1),2,2)
  b_line1 <- meanB
  B1 <- solve(aa,b_line1)
  jj <- sigma_dh - sigmaB[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh

  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanB[i])/sqrt(sigmaB[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("2MG-IE",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanB),4)," "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(sigmaB[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4)," "," "," "," "," "," ",round(B1[2],4)," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4)," "," ",round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################ PG-AI model#########################################
G3DHModelFun[[12]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(C-0)Model##########################################
  d2<- 1;mi <- as.matrix(1)
  sigma <- 2*sigma0
  meanC<- mean(dataDH)
  sigmaC <- sigma_dh/2
  mix_pi<-1.0
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+logL(n_samDH,1,1,meanC,sigmaC,dataDH)
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    swx<- sum((dataDH-meanC)^2)
    s0<-matrix(0,2,1)
    s0[1]<-(ss1+ss2);s0[2]<-swx*m_fam
    s1<-sigmaC - sigma/m_fam
    if (s1<0.0){ s1<- 0.0 }
    abc2<- sigma ; ss1<- 0 ;abc3<-1000
    while(abc3>0.0001 && ss1<1000 ){
      abc1<- (sigma/m_fam)/(sigma/m_fam+s1)
      sigma<- (s0[1]+abc1*abc1*s0[2])/(n_samP1+n_samP2+abc1*n_samDH)
      abc3<- abs(abc2-sigma)
      abc2<- sigma;ss1<-ss1+1.0
    }
    sigmaC <-s1+ sigma/m_fam
    ######## CM3 variance of polygenes ###########
    s1<- swx/n_samDH-sigma/m_fam;
    if (s1<0.0) { s1<- 0.0 }
    sigmaC<- s1+sigma/m_fam
    ######## criteria for iterations to stop #######

    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+logL(n_samDH,1,1,meanC,sigmaC,dataDH)
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) { stopa <- -stopa }
  }

  abc <- L1
  AIC <- -2*abc + 2*5

  #########first order genetic parameter process##########
  ma1<- mean11
  ma2<- mean12
  ma3<- meanC
  B1 <- as.matrix(c(ma1,ma2,ma3))
  mm <- sigma_dh-sigma
  if (mm<0) {mm<- 0}
  nnn<- mm/sigma_dh

  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)

  gg <- (dataDH - meanC)/sqrt(as.vector(sigmaC))
  bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
  bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
  bmwsl[,1] <- bmw*mix_pi
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("PG-AI",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(meanC,4)," "," "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(sigmaC[1],4),round(mix_pi[1],4),
                       " "," "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4),round(B1[2],4),round(B1[3],4)," "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(mm,4),round(nnn*100,4),round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################ PG-A model#########################################
G3DHModelFun[[13]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(C-1)Model##########################################
  d2<- 1;mi <- as.matrix(1)
  sigma <- 2*sigma0
  meanC<- mean(dataDH)
  sigmaC <- sigma_dh/2
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+logL(n_samDH,1,1,meanC,sigmaC,dataDH)
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    #######  E-step in ECM algorithm #######
    #######  CM-step in ECM algorithm #######
    #######  CM1-step for means #######
    mix_pi<-1.0
    rr<-(sum(dataP1)/n_samP1+sum(dataP2)/n_samP2-2.0*sum(dataDH)/n_samDH)/(sigma/n_samP1+sigma/n_samP2+4.0*sigmaC/n_samDH)

    mean11<- (sum(dataP1)-rr*sigma)/n_samP1
    mean12<- (sum(dataP2)-rr*sigma)/n_samP2
    meanC <- (sum(dataDH)+2.0*rr*sigmaC)/n_samDH
    #########  iteratively CM2-step for variance.
    #########  to calculate ss1, ss2
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    swx<- sum((dataDH-meanC)^2)

    s0<-matrix(0,2,1)
    s0[1]<-(ss1+ss2);s0[2]<-swx*m_fam
    s1<-sigmaC -sigma/m_fam  ########variance of polygenes.
    if (s1<0.0){ s1<- 0.0 }
    abc2<- sigma; ss1<- 0.0; abc3<-1000; n_iter<- 0
    while(abc3>0.0001){
      n_iter<- n_iter+1
      abc1<- (sigma/m_fam)/(sigma/m_fam+s1)
      sigma<- (s0[1]+abc1*abc1*s0[2])/(n_samP1+n_samP2+abc1*n_samDH)
      abc3<- abs(abc2-sigma)
      abc2<- sigma;ss1<- ss1+1.0
      if (n_iter>20) break
    }
    sigmaC <- s1+sigma/m_fam
    ########CM3 variance of polygenes #############
    s1<- swx/n_samDH-sigma/m_fam
    sigmaC <- s1+ sigma/m_fam

    ########criteria for iterations to stop#######

    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+logL(n_samDH,1,1,meanC,sigmaC,dataDH)
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*4

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1, 1,-1,0),3,2)
  b_line1 <- as.matrix(c(mean11,mean12,meanC))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)
  mm <- sigma_dh - sigma
  if(mm < 0) {mm <- 0}
  nnn <- mm/sigma_dh

  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)

  gg <- (dataDH - meanC)/sqrt(as.vector(sigmaC))
  bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
  bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
  bmwsl[,1] <- bmw*mix_pi

  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("PG-A",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(meanC,4)," "," "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(sigmaC[1],4),round(mix_pi[1],4),
                       " "," "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4)," "," "," "," "," "," "," "," "," "," "," "," "," ",round(B1[2],4)," "," ",round(mm,4),round(nnn*100,4),round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################ MX1-A-AI model#########################################
G3DHModelFun[[14]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(D-0)Model##########################################
  d2<- 2
  mi <- as.matrix(c(0.5,0.5))
  sigma <- 2*sigma0
  meanD<- mean(dataDH)
  sigmaD <- matrix((sigma_dh/2),d2,1)
  a1 <- sqrt(sigma_dh/n_samDH)
  if(mean11<mean12){a1=-a1}

  meanD <- as.matrix(c((meanD+3*a1/2),(meanD-3*a1/2)))

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanD,sqrt(sigmaD),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanD[i],sqrt(sigmaD[i]))/dmixnorm(dataDH,meanD,sqrt(sigmaD),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    mean11<- sum(dataP1)/n_samP1
    mean12<- sum(dataP2)/n_samP2
    meanD <- sumwx/n0

    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanD[i])^2 }

    s0<-matrix(0,2,1)
    s0[1]<- ss1+ss2
    s0[2]<- sum(swx)*m_fam
    s1<- sigmaD[1]-sigma/m_fam
    if (s1<0.0){s1<- 0.000001}
    abc2<- sigma

    abc3<-1000; n_iter<- 0
    while(abc3>0.0001){
      n_iter<- n_iter+1
      abc1<- (sigma/m_fam)/(sigma/m_fam+s1)
      sigma<- (s0[1]+abc1*abc1*s0[2])/(n_samP1+n_samP2+abc1*n_samDH)
      abc3<- abs(abc2-sigma)
      abc2<- sigma
      if (n_iter>20) break
    }
    if (sigma<0.1*sigma0){ sigma<- 0.1*sigma0 }
    sigmaD[1]<- s1+sigma/m_fam
    s1<- sum(swx)/n_samDH-sigma/m_fam
    if (s1<0.0){ s1<- 0.000001 }
    sigmaD<- matrix((s1+sigma/m_fam),d2,1)

    ########criteria for iterations to stop#######

    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanD,sqrt(sigmaD),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*6

  #########first order genetic parameter process##########
  aa<- matrix(c(1,0,0,0, 0,1,0,0, 0,0,1,1, 1,-1,1,-1),4,4)
  b_line1 <- as.matrix(c(mean11,mean12,meanD))
  B1 <- solve(aa,b_line1)
  jj <- sigma_dh - sigmaD[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh
  mm<- sigmaD[1]-sigma
  if (mm<0) {mm<- 0}
  nnn<- mm/sigma_dh
  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanD[i])/sqrt(sigmaD[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX1-A-AI",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanD),4)," "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(sigmaD[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4),round(B1[2],4),round(B1[3],4),round(B1[3],4)," "," "," "," "," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4),round(mm,4),round(nnn*100,4),round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################ MX1-A-A model#########################################
G3DHModelFun[[15]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(D-1)Model##########################################
  d2<- 2
  mi <- as.matrix(c(0.5,0.5))
  sigma <- 2*sigma0
  meanD<- mean(dataDH)
  sigmaD <- matrix((sigma_dh/2),d2,1)
  a1 <- sqrt(sigma_dh/n_samDH)
  if(mean11<mean12){a1=-a1}

  meanD <- as.matrix(c((meanD+3*a1/2),(meanD-3*a1/2)))

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanD,sqrt(sigmaD),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanD[i],sqrt(sigmaD[i]))/dmixnorm(dataDH,meanD,sqrt(sigmaD),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    rr<- (sum(dataP1)/n_samP1+sum(dataP2)/n_samP2-sumwx[1]/n0[1]-sumwx[2]/n0[2])/(sigma/n_samP1+sigma/n_samP2+sigmaD[1]/n0[1]+sigmaD[1]/n0[2])

    mean11<- (sum(dataP1)-rr*sigma)/n_samP1
    mean12<- (sum(dataP2)-rr*sigma)/n_samP2
    meanD <- (sumwx+rr*sigmaD[1])/n0
    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanD[i])^2 }

    s0<-matrix(0,2,1)
    s0[1]<- ss1+ss2
    s0[2]<- sum(swx)*m_fam
    s1<- sigmaD[1]-sigma/m_fam
    if (s1<0.0){s1<- 0.000001}
    abc2<- sigma

    abc3<-1000; n_iter<- 0
    while(abc3>0.0001){
      n_iter<- n_iter+1
      abc1<- (sigma/m_fam)/(sigma/m_fam+s1)
      sigma<- (s0[1]+abc1*abc1*s0[2])/(n_samP1+n_samP2+abc1*n_samDH)
      abc3<- abs(abc2-sigma)
      abc2<- sigma
      if (n_iter>20) break
    }
    if (sigma<0.1*sigma0){ sigma<- 0.1*sigma0 }
    sigmaD[1]<- s1+sigma/m_fam
    s1<- sum(swx)/n_samDH-sigma/m_fam
    if (s1<0.0){ s1<- 0 }
    sigmaD<- matrix((s1+sigma/m_fam),d2,1)

    ########criteria for iterations to stop#######

    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanD,sqrt(sigmaD),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*5

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1, 1,-1,1,-1, 1,-1,0,0),4,3)
  b_line1 <- as.matrix(c(mean11,mean12,meanD))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)
  jj <- sigma_dh - sigmaD[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh
  mm<- sigmaD[1]-sigma
  if (mm<0) {mm<- 0}
  nnn<- mm/sigma_dh

  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanD[i])/sqrt(sigmaD[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX1-A-A",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanD),4)," "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(sigmaD[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4)," "," ",round(B1[2],4)," "," "," "," "," "," "," "," "," "," ",round(B1[3],4),round(jj,4),round(ll*100,4),round(mm,4),round(nnn*100,4),round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################ MX2-AI-AI model#########################################
G3DHModelFun[[16]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(E-O)Model##########################################
  d2<-4
  mi <- as.matrix(c(0.25,0.25,0.25,0.25))
  meanE<-mean(dataDH)
  sigmaE <- matrix((sigma_dh/2*1.2222101),d2,1)
  sigma <- sigma0
  a1 <- sqrt(sigma_dh/n_samDH)
  if(mean11<mean12){a1=-a1}

  meanE <- as.matrix(c((meanE+3*a1),(meanE+1.5*a1),(meanE-1.5*a1),(meanE-3*a1)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanE,sqrt(sigmaE),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanE[i],sqrt(sigmaE[i]))/dmixnorm(dataDH,meanE,sqrt(sigmaE),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    mean11<- sum(dataP1)/n_samP1
    mean12<- sum(dataP2)/n_samP2
    meanE <- sumwx/n0

    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanE[i])^2 }

    s0<-matrix(0,2,1)
    s0[1]<- ss1+ss2
    s0[2]<- sum(swx)*m_fam
    s1<- sigmaE[1]-sigma/m_fam
    if (s1<0.0){s1<- 0.000001}
    abc2<- sigma; abc3<-1000; n_iter<- 0
    while(abc3>0.0001){
      n_iter<- n_iter+1
      abc1<- (sigma/m_fam)/(sigma/m_fam+s1)
      sigma<- (s0[1]+abc1*abc1*s0[2])/(n_samP1+n_samP2+abc1*n_samDH)
      abc3<- abs(abc2-sigma)
      abc2<- sigma
      if (n_iter>20) break
    }
    if (sigma<0.1*sigma0){ sigma<- 0.1*sigma0 }
    sigmaE[1]<- s1+sigma/m_fam
    s1<- sum(swx)/n_samDH-sigma/m_fam
    if (s1<0.0){ s1<- 0.000001 }
    sigmaE<- matrix((s1+sigma/m_fam),d2,1)

    ########criteria for iterations to stop#######

    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanE,sqrt(sigmaE),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*8

  #########first order genetic parameter process##########
  aa<- matrix(c(1,0,0,0,0,0, 0,1,0,0,0,0, 0,0,1,1,1,1, 1,-1,1,1,-1,-1, 1,-1,1,-1,1,-1, 1,1,1,-1,-1,1),6,6)
  b_line1 <- as.matrix(c(mean11,mean12,meanE))
  B1 <- solve(aa,b_line1)
  jj <- sigma_dh - sigmaE[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh
  mm<- sigmaE[1]-sigma[1]
  if(mm<0){mm<-0}
  nnn<- mm/sigma_dh

  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanE[i])/sqrt(sigmaE[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX2-AI-AI",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanE),4)," "," "," "," "," "," "," "," "," "," "," "," ",round(sigmaE[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4),round(B1[2],4),round(B1[3],4),round(B1[4],4),round(B1[5],4)," "," ",round(B1[6],4)," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4),round(mm,4),round(nnn*100,4),round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}

############################################ MX2-AI-A model#########################################
G3DHModelFun[[17]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(E-1)Model##########################################
  d2<-4
  mi <- as.matrix(c(0.25,0.25,0.25,0.25))
  meanE<-mean(dataDH)
  sigmaE <- matrix((sigma_dh/2*1.2222101),d2,1)
  sigma <- sigma0
  a1 <- sqrt(sigma_dh/n_samDH)
  if(mean11<mean12){a1=-a1}

  meanE <- as.matrix(c((meanE+3*a1),(meanE+1.5*a1),(meanE-1.5*a1),(meanE-3*a1)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanE,sqrt(sigmaE),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanE[i],sqrt(sigmaE[i]))/dmixnorm(dataDH,meanE,sqrt(sigmaE),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    rr<-(sum(dataP1)/n_samP1+sum(dataP2)/n_samP2-sumwx[1]/n0[1]-sumwx[4]/n0[4])/(sigma/n_samP1+sigma/n_samP2+sigmaE[1]/n0[1]+sigmaE[4]/n0[4])

    mean11<- (sum(dataP1)-rr*sigma)/n_samP1
    mean12<- (sum(dataP2)-rr*sigma)/n_samP2
    meanE[1]<- (sumwx[1]+rr*sigmaE[1])/n0[1]
    meanE[2]<- sumwx[2]/n0[2]
    meanE[3]<- sumwx[3]/n0[3]
    meanE[4]<- (sumwx[4]+rr*sigmaE[4])/n0[4]
    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanE[i])^2 }

    s0<-matrix(0,2,1)
    s0[1]<- ss1+ss2
    s0[2]<- sum(swx)*m_fam
    s1<- sigmaE[1]-sigma/m_fam
    if (s1<0.0){s1<- 0.000001}
    abc2<- sigma

    abc3<-1000; n_iter<- 0
    while(abc3>0.0001){
      n_iter<- n_iter+1
      abc1<- (sigma/m_fam)/(sigma/m_fam+s1)
      sigma<- (s0[1]+abc1*abc1*s0[2])/(n_samP1+n_samP2+abc1*n_samDH)
      abc3<- abs(abc2-sigma)
      abc2<- sigma
      if (n_iter>20) break
    }
    if (sigma<0.1*sigma0){ sigma<- 0.1*sigma0 }
    sigmaE[1]<- s1+sigma/m_fam
    s1<- sum(swx)/n_samDH-sigma/m_fam
    if (s1<0.0){ s1<- 0.000001 }
    sigmaE<- matrix((s1+sigma/m_fam),d2,1)

    ########criteria for iterations to stop#######
    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanE,sqrt(sigmaE),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*7

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1, 1,-1,1,1,-1,-1, 1,-1,1,-1,1,-1, 1,-1,0,0,0,0, 1,1,1,-1,-1,1),6,5)
  b_line1 <- as.matrix(c(mean11,mean12,meanE))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)
  jj <- sigma_dh - sigmaE[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh
  mm<- sigmaE[1]-sigma[1]
  if(mm<0){mm<-0}
  nnn<- mm/sigma_dh

  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanE[i])/sqrt(sigmaE[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX2-AI-A",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanE),4)," "," "," "," "," "," "," "," "," "," "," "," ",round(sigmaE[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4)," "," ",round(B1[2],4),round(B1[3],4)," "," ",round(B1[4],4)," "," "," "," "," "," ",round(B1[5],4),round(jj,4),round(ll*100,4),round(mm,4),round(nnn*100,4),round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}

############################################ MX2-A-A model#########################################
G3DHModelFun[[18]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(E-2)Model##########################################
  d2<-4
  mi <- as.matrix(c(0.25,0.25,0.25,0.25))
  meanE<-mean(dataDH)
  sigmaE <- matrix((sigma_dh/2*1.2222101),d2,1)
  sigma <- sigma0
  a1 <- sqrt(sigma_dh)
  if(mean11<mean12){a1=-a1}

  meanE <- as.matrix(c((meanE+3*a1),(meanE+1.5*a1),(meanE-1.5*a1),(meanE-3*a1)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanE,sqrt(sigmaE),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanE[i],sqrt(sigmaE[i]))/dmixnorm(dataDH,meanE,sqrt(sigmaE),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    s0<- matrix(0,6,1)
    s0[1]<- sum(dataP1)/n_samP1+sum(dataP2)/n_samP2-sumwx[1]/n0[1]-sumwx[4]/n0[4]
    s0[2]<- sumwx[1]/n0[1]-sumwx[2]/n0[2]-sumwx[3]/n0[3]+sumwx[4]/n0[4]
    s0[3]<- sigma/n_samP1+sigma/n_samP2+sigmaE[1]/n0[1]+sigmaE[4]/n0[4]
    s0[4]<- -sigmaE[1]/n0[1]-sigmaE[4]/n0[4]
    s0[5]<- sigmaE[1]/n0[1]+sigmaE[2]/n0[2]+sigmaE[3]/n0[3]+sigmaE[4]/n0[4]
    s0[6]<- s0[3]*s0[5]-s0[4]*s0[4]
    rr<- matrix(0,2,1)
    rr[1]<- (s0[1]*s0[5]-s0[2]*s0[4])/s0[6]
    rr[2]<- (s0[2]*s0[3]-s0[1]*s0[4])/s0[6]

    mean11<- (sum(dataP1)-rr[1]*sigma)/n_samP1
    mean12<- (sum(dataP2)-rr[1]*sigma)/n_samP2
    meanE[1]<- (sumwx[1]+sigmaE[1]*(rr[1]-rr[2]))/n0[1]
    meanE[2]<- (sumwx[2]+sigmaE[2]*rr[2])/n0[2]
    meanE[3]<- (sumwx[3]+sigmaE[3]*rr[2])/n0[3]
    meanE[4]<- (sumwx[4]+sigmaE[4]*(rr[1]-rr[2]))/n0[4]

    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanE[i])^2 }

    s0E<-matrix(0,2,1)
    s0E[1]<- ss1+ss2
    s0E[2]<- sum(swx)*m_fam
    s1<- sigmaE[1]-sigma/m_fam
    if (s1<0.0){s1<- 0.000001 }
    abc2<- sigma

    abc3<-1000; n_iter<- 0
    while(abc3>0.0001){
      n_iter<- n_iter+1
      abc1<- (sigma/m_fam)/(sigma/m_fam+s1)
      sigma<- (s0E[1]+abc1*abc1*s0E[2])/(n_samP1+n_samP2+abc1*n_samDH)
      abc3<- abs(abc2-sigma)
      abc2<- sigma
      if (n_iter>20) break
    }
    if (sigma<0.1*sigma0){ sigma<- 0.1*sigma0 }
    sigmaE[1]<- s1+sigma/m_fam
    s1<- sum(swx)/n_samDH-sigma/m_fam
    if (s1<0.0){ s1<- 0.000001 }
    sigmaE<- matrix((s1+sigma/m_fam),d2,1)

    ########criteria for iterations to stop#######

    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanE,sqrt(sigmaE),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*6

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1, 1,-1,1,1,-1,-1, 1,-1,1,-1,1,-1, 1,-1,0,0,0,0),6,4)
  b_line1 <- as.matrix(c(mean11,mean12,meanE))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)
  jj <- sigma_dh - sigmaE[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh
  mm<- sigmaE[1]-sigma[1]
  if(mm<0){mm<-0}
  nnn<- mm/sigma_dh

  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanE[i])/sqrt(sigmaE[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX2-A-A",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanE),4)," "," "," "," "," "," "," "," "," "," "," "," ",round(sigmaE[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4)," "," ",round(B1[2],4),round(B1[3],4)," "," "," "," "," "," "," "," "," ",round(B1[4],4),round(jj,4),round(ll*100,4),round(mm,4),round(nnn*100,4),round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}

############################################ MX2-EA-A model#########################################
G3DHModelFun[[19]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(E-3)Model##########################################
  d2<-3
  mi <- as.matrix(c(0.25,0.5,0.25))
  meanE<-mean(dataDH)
  sigmaE <- matrix((sigma_dh/2*1.2222101),d2,1)
  sigma <- sigma0
  a1 <- sqrt(sigma_dh)
  if(mean11<mean12){a1=-a1}

  #meanE <- as.matrix(c((meanE+2.5*a1),meanE,(meanE-2.5*a1)))
  meanE <- as.matrix(c(114,100,86))

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanE,sqrt(sigmaE),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanE[i],sqrt(sigmaE[i]))/dmixnorm(dataDH,meanE,sqrt(sigmaE),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    s0<- matrix(0,6,1)
    s0[1]<- sum(dataP1)/n_samP1+sum(dataP2)/n_samP2-2*sumwx[2]/n0[2]
    s0[2]<- sumwx[1]/n0[1]-2*sumwx[2]/n0[2]+sumwx[3]/n0[3]
    s0[3]<- sigma/n_samP1+sigma/n_samP2+4*sigmaE[2]/n0[2]
    s0[4]<- 4*sigmaE[2]/n0[2]
    s0[5]<- sigmaE[1]/n0[1]+4*sigmaE[2]/n0[2]+sigmaE[3]/n0[3]
    s0[6]<- s0[3]*s0[5]-s0[4]*s0[4]
    rr<- matrix(0,2,1)
    rr[1]<- (s0[1]*s0[5]-s0[2]*s0[4])/s0[6]
    rr[2]<- (s0[2]*s0[3]-s0[1]*s0[4])/s0[6]

    mean11<- (sum(dataP1)-rr[1]*sigma)/n_samP1
    mean12<- (sum(dataP2)-rr[1]*sigma)/n_samP2
    meanE[1]<- (sumwx[1]-sigmaE[1]*rr[2])/n0[1]
    meanE[2]<- (sumwx[2]+sigmaE[2]*(2.0*rr[1]+2.0*rr[2]))/n0[2]
    meanE[3]<- (sumwx[3]-sigmaE[3]*rr[2])/n0[3]

    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanE[i])^2 }

    s0E<-matrix(0,2,1)
    s0E[1]<- ss1+ss2
    s0E[2]<- sum(swx)*m_fam
    s1<- sigmaE[1]-sigma/m_fam
    if (s1<0.0){s1<- 0.000001}
    abc2<- sigma

    abc3<-1000; n_iter<- 0
    while(abc3>0.0001){
      n_iter<- n_iter+1
      abc1<- (sigma/m_fam)/(sigma/m_fam+s1)
      sigma<- (s0E[1]+abc1*abc1*s0E[2])/(n_samP1+n_samP2+abc1*n_samDH)
      abc3<- abs(abc2-sigma)
      abc2<- sigma
      if (n_iter>20) break
    }
    if (sigma<0.1*sigma0){ sigma<- 0.1*sigma0 }
    sigmaE[1]<- s1+sigma/m_fam
    s1<- sum(swx)/n_samDH-sigma/m_fam
    if (s1<0.0){ s1<- 0.000001 }
    sigmaE<- matrix((s1+sigma/m_fam),d2,1)

    ########criteria for iterations to stop#######
    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanE,sqrt(sigmaE),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*5

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1, 2,-2,2,0,-2, 1,1,0,0,0),5,3)
  b_line1 <- as.matrix(c(mean11,mean12,meanE))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)
  jj <- sigma_dh - sigmaE[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh
  mm<- sigmaE[1]-sigma[1]
  if(mm<0){mm<-0}
  nnn<- mm/sigma_dh

  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanE[i])/sqrt(sigmaE[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX2-EA-A",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanE),4)," "," "," "," "," "," "," "," "," "," "," "," "," ",round(sigmaE[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4)," "," ",round(B1[2],4),round(B1[2],4)," "," "," "," "," "," "," "," "," ",round(B1[3],4),round(jj,4),round(ll*100,4),round(mm,4),round(nnn*100,4),round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}

############################################ MX2-ED-A model#########################################
G3DHModelFun[[20]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(E-4)Model##########################################
  d2<-3
  mi <- as.matrix(c(0.5,0.25,0.25))
  meanE<-mean(dataDH)
  sigmaE <- matrix((sigma_dh/2*1.2222101),d2,1)
  sigma <- sigma0
  a1 <- sqrt(sigma_dh/n_samDH)
  if(mean11<mean12){a1=-a1}

  meanE <- as.matrix(c((meanE+2*a1),meanE,(meanE-2*a1)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanE,sqrt(sigmaE),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanE[i],sqrt(sigmaE[i]))/dmixnorm(dataDH,meanE,sqrt(sigmaE),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    rr<- (sum(dataP1)/n_samP1+sum(dataP2)/n_samP2-sumwx[1]/n0[1]-sumwx[3]/n0[3])/(sigma/n_samP1+sigma/n_samP2+sigmaE[1]/n0[1]+sigmaE[3]/n0[3])

    mean11<- (sum(dataP1)-rr[1]*sigma)/n_samP1
    mean12<- (sum(dataP2)-rr[1]*sigma)/n_samP2
    meanE[1]<- (sumwx[1]+sigmaE[1]*rr)/n0[1]
    meanE[2]<- sumwx[2]/n0[2]
    meanE[3]<- (sumwx[3]+sigmaE[3]*rr)/n0[3]

    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanE[i])^2 }

    s0E<-matrix(0,2,1)
    s0E[1]<- ss1+ss2
    s0E[2]<- sum(swx)*m_fam
    s1<- sigmaE[1]-sigma/m_fam
    if (s1<0.0){s1<- 0.000001}
    abc2<- sigma

    abc3<-1000; n_iter<- 0
    while(abc3>0.0001){
      n_iter<- n_iter+1
      abc1<- (sigma/m_fam)/(sigma/m_fam+s1)
      sigma<- (s0E[1]+abc1*abc1*s0E[2])/(n_samP1+n_samP2+abc1*n_samDH)
      abc3<- abs(abc2-sigma)
      abc2<- sigma
      if (n_iter>20) break
    }
    if (sigma<0.1*sigma0){ sigma<- 0.1*sigma0 }
    sigmaE[1]<- s1+sigma/m_fam
    s1<- sum(swx)/n_samDH-sigma/m_fam
    if (s1<0.0){ s1<- 0.000001 }
    sigmaE<- matrix((s1+sigma/m_fam),d2,1)

    ########criteria for iterations to stop#######
    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanE,sqrt(sigmaE),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*6

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1, 1,-1,1,-1,-1, 0,-1,0,1,-1, 1,-1,0,0,0),5,4)
  b_line1 <- as.matrix(c(mean11,mean12,meanE))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)
  jj <- sigma_dh - sigmaE[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh
  mm<- sigmaE[1]-sigma[1]
  if(mm<0){mm<-0}
  nnn<- mm/sigma_dh

  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanE[i])/sqrt(sigmaE[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX2-ED-A",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanE),4)," "," "," "," "," "," "," "," "," "," "," "," "," ",round(sigmaE[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4)," "," ",round(B1[2],4),round(B1[3],4)," "," "," "," "," "," "," "," "," ",round(B1[4],4),round(jj,4),round(ll*100,4),round(mm,4),round(nnn*100,4),round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################ MX2-ER-A model#########################################
G3DHModelFun[[21]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(E-5)Model##########################################
  d2<-3
  mi <- as.matrix(c(0.25,0.25,0.5))
  meanE<-mean(dataDH)
  sigmaE <- matrix((sigma_dh/2*1.2222101),d2,1)
  sigma <- sigma0
  a1 <- sqrt(sigma_dh/n_samDH)
  if(mean11<mean12){a1=-a1}

  meanE <- as.matrix(c((meanE+2*a1),meanE,(meanE-2*a1)))

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanE,sqrt(sigmaE),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanE[i],sqrt(sigmaE[i]))/dmixnorm(dataDH,meanE,sqrt(sigmaE),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    rr<- (sum(dataP1)/n_samP1+sum(dataP2)/n_samP2-sumwx[1]/n0[1]-sumwx[3]/n0[3])/(sigma/n_samP1+sigma/n_samP2+sigmaE[1]/n0[1]+sigmaE[3]/n0[3])

    mean11<- (sum(dataP1)-rr[1]*sigma)/n_samP1
    mean12<- (sum(dataP2)-rr[1]*sigma)/n_samP2
    meanE[1]<- (sumwx[1]+sigmaE[1]*rr)/n0[1]
    meanE[2]<- sumwx[2]/n0[2]
    meanE[3]<- (sumwx[3]+sigmaE[3]*rr)/n0[3]

    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanE[i])^2 }

    s0E<-matrix(0,2,1)
    s0E[1]<- ss1+ss2
    s0E[2]<- sum(swx)*m_fam
    s1<- sigmaE[1]-sigma/m_fam
    if (s1<0.0){s1<- 0.000001}
    abc2<- sigma

    abc3<-1000; n_iter<- 0
    while(abc3>0.0001){
      n_iter<- n_iter+1
      abc1<- (sigma/m_fam)/(sigma/m_fam+s1)
      sigma<- (s0E[1]+abc1*abc1*s0E[2])/(n_samP1+n_samP2+abc1*n_samDH)
      abc3<- abs(abc2-sigma)
      abc2<- sigma
      if (n_iter>20) break
    }
    if (sigma<0.1*sigma0){ sigma<- 0.1*sigma0 }
    sigmaE[1]<- s1+sigma/m_fam
    s1<- sum(swx)/n_samDH-sigma/m_fam
    if (s1<0.0){ s1<- 0.000001 }
    sigmaE<- matrix((s1+sigma/m_fam),d2,1)

    ########criteria for iterations to stop#######

    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanE,sqrt(sigmaE),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*6

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1, 1,-1,1,1,-1, 1,0,1,-1,0, 1,-1,0,0,0),5,4)
  b_line1 <- as.matrix(c(mean11,mean12,meanE))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)
  jj <- sigma_dh - sigmaE[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh
  mm<- sigmaE[1]-sigma[1]
  if(mm<0){mm<-0}
  nnn<- mm/sigma_dh

  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanE[i])/sqrt(sigmaE[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX2-ER-A",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanE),4)," "," "," "," "," "," "," "," "," "," "," "," "," ",round(sigmaE[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4)," "," ",round(B1[2],4),round(B1[3],4)," "," "," "," "," "," "," "," "," ",round(B1[4],4),round(jj,4),round(ll*100,4),round(mm,4),round(nnn*100,4),round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################ MX2-AE-A model#########################################
G3DHModelFun[[22]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(E-6)Model##########################################
  d2<-3
  mi <- as.matrix(c(0.25,0.5,0.25))
  meanE<-mean(dataDH)
  sigmaE <- matrix((sigma_dh/2*1.2222101),d2,1)
  sigma <- sigma0
  a1 <- sqrt(sigma_dh/n_samDH)
  if(mean11<mean12){a1=-a1}

  meanE <- as.matrix(c((meanE+2*a1),meanE,(meanE-2*a1)))

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanE,sqrt(sigmaE),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanE[i],sqrt(sigmaE[i]))/dmixnorm(dataDH,meanE,sqrt(sigmaE),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    rr<- (sum(dataP1)/n_samP1+sum(dataP2)/n_samP2-sumwx[1]/n0[1]-sumwx[3]/n0[3])/(sigma/n_samP1+sigma/n_samP2+sigmaE[1]/n0[1]+sigmaE[3]/n0[3])

    mean11<- (sum(dataP1)-rr*sigma)/n_samP1
    mean12<- (sum(dataP2)-rr*sigma)/n_samP2
    meanE[1]<- (sumwx[1]+sigmaE[1]*rr)/n0[1]
    meanE[2]<- sumwx[2]/n0[2]
    meanE[3]<- (sumwx[3]+sigmaE[1]*rr)/n0[3]

    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanE[i])^2 }

    s0E<-matrix(0,2,1)
    s0E[1]<- ss1+ss2
    s0E[2]<- sum(swx)*m_fam
    s1<- sigmaE[1]-sigma/m_fam
    if (s1<0.0){s1<- 0.000001}
    abc2<- sigma

    abc3<-1000; n_iter<- 0
    while(abc3>0.0001){
      n_iter<- n_iter+1
      abc1<- (sigma/m_fam)/(sigma/m_fam+s1)
      sigma<- (s0E[1]+abc1*abc1*s0E[2])/(n_samP1+n_samP2+abc1*n_samDH)
      abc3<- abs(abc2-sigma)
      abc2<- sigma
      if (n_iter>20) break
    }
    if (sigma<0.1*sigma0){ sigma<- 0.1*sigma0 }
    sigmaE[1]<- s1+sigma/m_fam
    s1<- sum(swx)/n_samDH-sigma/m_fam
    if (s1<0.0){ s1<- 0.000001 }
    sigmaE<- matrix((s1+sigma/m_fam),d2,1)

    ########criteria for iterations to stop#######

    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanE,sqrt(sigmaE),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*6

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1, 2,-2,2,0,-2, 1,1,1,-1,1, 1,-1,0,0,0),5,4)
  b_line1 <- as.matrix(c(mean11,mean12,meanE))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)
  jj <- sigma_dh - sigmaE[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh
  mm<- sigmaE[1]-sigma[1]
  if(mm<0){mm<-0}
  nnn<- mm/sigma_dh

  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanE[i])/sqrt(sigmaE[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX2-AE-A",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanE),4)," "," "," "," "," "," "," "," "," "," "," "," "," ",round(sigmaE[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4)," "," ",round(B1[2],4)," "," "," ",round(B1[3],4)," "," "," "," "," "," ",round(B1[4],4),round(jj,4),round(ll*100,4),round(mm,4),round(nnn*100,4),round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################ MX2-CE-A model#########################################
G3DHModelFun[[23]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(E-7)Model##########################################
  d2<-2
  mi <- as.matrix(c(0.25,0.75))
  meanE<-mean(dataDH)
  sigmaE <- matrix((sigma_dh/2*1.2222101),d2,1)
  sigma <- sigma0
  a1 <- sqrt(sigma_dh/n_samDH)
  if(mean11<mean12){a1=-a1}

  meanE <- as.matrix(c((meanE+1.5*a1),(meanE-1.5*a1)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanE,sqrt(sigmaE),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanE[i],sqrt(sigmaE[i]))/dmixnorm(dataDH,meanE,sqrt(sigmaE),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    rr<- (sum(dataP1)/n_samP1+sum(dataP2)/n_samP2-2*sumwx[1]/n0[1])/(sigma/n_samP1+sigma/n_samP2+sigmaE[1]/n0[1])
    mean11<- (sum(dataP1)-rr*sigma)/n_samP1
    mean12<- (sum(dataP2)-rr*sigma)/n_samP2
    meanE[1]<- (sumwx[1]+sigmaE[1]*rr)/n0[1]
    meanE[2]<- (sumwx[2]+sigmaE[2]*rr)/n0[2]

    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanE[i])^2 }

    s0E<-matrix(0,2,1)
    s0E[1]<- ss1+ss2
    s0E[2]<- sum(swx)*m_fam
    s1<- sigmaE[1]-sigma/m_fam
    if (s1<0.0){s1<- 0.000001}
    abc2<- sigma

    abc3<-1000; n_iter<- 0
    while(abc3>0.0001){
      n_iter<- n_iter+1
      abc1<- (sigma/m_fam)/(sigma/m_fam+s1)
      sigma<- (s0E[1]+abc1*abc1*s0E[2])/(n_samP1+n_samP2+abc1*n_samDH)
      abc3<- abs(abc2-sigma)
      abc2<- sigma
      if (n_iter>20) break
    }
    if (sigma<0.1*sigma0){ sigma<- 0.1*sigma0 }
    sigmaE[1]<- s1+sigma/m_fam
    s1<- sum(swx)/n_samDH-sigma/m_fam
    if (s1<0.0){ s1<- 0.000001 }
    sigmaE<- matrix((s1+sigma/m_fam),d2,1)

    ########criteria for iterations to stop#######

    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanE,sqrt(sigmaE),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*5

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1, 1,1,1,-1, 1,-1,0,0),4,3)
  b_line1 <- as.matrix(c(mean11,mean12,meanE))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)
  jj <- sigma_dh - sigmaE[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh
  mm<- sigmaE[1]-sigma[1]
  if(mm<0){mm<-0}
  nnn<- mm/sigma_dh

  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanE[i])/sqrt(sigmaE[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX2-CE-A",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanE),4)," "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(sigmaE[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4)," "," "," "," "," "," ",round(B1[2],4)," "," "," "," "," "," ",round(B1[3],4),round(jj,4),round(ll*100,4),round(mm,4),round(nnn*100,4),round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################ MX2-DE-A model#########################################
G3DHModelFun[[24]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(E-8)Model##########################################
  d2<-2
  mi <- as.matrix(c(0.75,0.25))
  meanE<-mean(dataDH)
  sigmaE <- matrix((sigma_dh/2*1.2222101),d2,1)
  sigma <- sigma0
  a1 <- sqrt(sigma_dh/n_samDH)
  if(mean11<mean12){a1=-a1}
  meanE <- as.matrix(c((meanE+2*a1),(meanE-2*a1)))

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanE,sqrt(sigmaE),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanE[i],sqrt(sigmaE[i]))/dmixnorm(dataDH,meanE,sqrt(sigmaE),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    rr<- (sum(dataP1)/n_samP1+sum(dataP2)/n_samP2-sumwx[1]/n0[1]-sumwx[2]/n0[2])/(sigma/n_samP1+sigma/n_samP2+sigmaE[1]/n0[1]+sigmaE[2]/n0[2])
    mean11<- (sum(dataP1)-rr*sigma)/n_samP1
    mean12<- (sum(dataP2)-rr*sigma)/n_samP2
    meanE[1]<- (sumwx[1]+sigmaE[1]*rr)/n0[1]
    meanE[2]<- (sumwx[2]+sigmaE[2]*rr)/n0[2]

    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanE[i])^2 }

    s0E<-matrix(0,2,1)
    s0E[1]<- ss1+ss2
    s0E[2]<- sum(swx)*m_fam
    s1<- sigmaE[1]-sigma/m_fam
    if (s1<0.0){s1<- 0.000001}
    abc2<- sigma

    abc3<-1000; n_iter<- 0
    while(abc3>0.0001){
      n_iter<- n_iter+1
      abc1<- (sigma/m_fam)/(sigma/m_fam+s1)
      sigma<- (s0E[1]+abc1*abc1*s0E[2])/(n_samP1+n_samP2+abc1*n_samDH)
      abc3<- abs(abc2-sigma)
      abc2<- sigma
      if (n_iter>20) break
    }
    if (sigma<0.1*sigma0){ sigma<- 0.1*sigma0 }
    sigmaE[1]<- s1+sigma/m_fam
    s1<- sum(swx)/n_samDH-sigma/m_fam
    if (s1<0.0){ s1<- 0.000001 }
    sigmaE<- matrix((s1+sigma/m_fam),d2,1)

    ########criteria for iterations to stop#######
    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanE,sqrt(sigmaE),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*5

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1, 1,1,1,-1, 1,-1,0,0),4,3)
  b_line1 <- as.matrix(c(mean11,mean12,meanE))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)
  jj <- sigma_dh - sigmaE[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh
  mm<- sigmaE[1]-sigma[1]
  if(mm<0){mm<-0}
  nnn<- mm/sigma_dh

  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanE[i])/sqrt(sigmaE[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX2-DE-A",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanE),4)," "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(sigmaE[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4)," "," "," "," "," "," ",round(B1[2],4)," "," "," "," "," "," ",round(B1[3],4),round(jj,4),round(ll*100,4),round(mm,4),round(nnn*100,4),round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################ MX2-IE-A model#########################################
G3DHModelFun[[25]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(E-9)Model##########################################
  d2<-2
  mi <- as.matrix(c(0.75,0.25))
  meanE<-mean(dataDH)
  sigmaE <- matrix((sigma_dh/2*1.2222101),d2,1)
  sigma <- sigma0
  a1 <- sqrt(sigma_dh/n_samDH)
  if(mean11<mean12){a1=-a1}
  meanE <- as.matrix(c((meanE+2*a1),(meanE-2*a1)))

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanE,sqrt(sigmaE),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanE[i],sqrt(sigmaE[i]))/dmixnorm(dataDH,meanE,sqrt(sigmaE),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    rr<-(sum(dataP1)/n_samP1+sum(dataP2)/n_samP2-2.0*sumwx[2]/n0[2])/(sigma/n_samP1+sigma/n_samP2+4.0*sigmaE[2]/n0[2])
    mean11<- (sum(dataP1)-rr*sigma)/n_samP1
    mean12<- (sum(dataP2)-rr*sigma)/n_samP2
    meanE[1]<- sumwx[1]/n0[1]
    meanE[2]<- (sumwx[2]+2.0*rr*sigmaE[2])/n0[2]

    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanE[i])^2 }

    s0E<-matrix(0,2,1)
    s0E[1]<- ss1+ss2
    s0E[2]<- sum(swx)*m_fam
    s1<- sigmaE[1]-sigma/m_fam
    if (s1<0.0){s1<- 0.000001}
    abc2<- sigma

    abc3<-1000; n_iter<- 0
    while(abc3>0.0001){
      n_iter<- n_iter+1
      abc1<- (sigma/m_fam)/(sigma/m_fam+s1)
      sigma<- (s0E[1]+abc1*abc1*s0E[2])/(n_samP1+n_samP2+abc1*n_samDH)
      abc3<- abs(abc2-sigma)
      abc2<- sigma
      if (n_iter>20) break
    }
    if (sigma<0.1*sigma0){ sigma<- 0.1*sigma0 }
    sigmaE[1]<- s1+sigma/m_fam
    s1<- sum(swx)/n_samDH-sigma/m_fam
    if (s1<0.0){ s1<- 0.000001 }
    sigmaE<- matrix((s1+sigma/m_fam),d2,1)

    ########criteria for iterations to stop#######

    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanE,sqrt(sigmaE),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*5

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1, 1,1,-1,1, 1,-1,0,0),4,3)
  b_line1 <- as.matrix(c(mean11,mean12,meanE))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)
  jj <- sigma_dh - sigmaE[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh
  mm<- sigmaE[1]-sigma[1]
  if(mm<0){mm<-0}
  nnn<- mm/sigma_dh

  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanE[i])/sqrt(sigmaE[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX2-IE-A",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanE),4)," "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(sigmaE[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4)," "," "," "," "," "," ",round(B1[2],4)," "," "," "," "," "," ",round(B1[3],4),round(jj,4),round(ll*100,4),round(mm,4),round(nnn*100,4),round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################ 3MG-AI model#########################################
G3DHModelFun[[26]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(F-1)Model##########################################
  d2<-8
  mi <- as.matrix(c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125))
  meanF<-mean(dataDH)
  sigmaF <- matrix((sigma_dh/(2*1.2222101)),d2,1)
  sigma <- sigma0
  a1 <- sqrt(sigma_dh/n_samDH)
  if(mean11<mean12){a1=-a1}

  meanF <- as.matrix(c((meanF+3*a1),(meanF+2.1*a1),(meanF+1.2*a1),(meanF+0.3*a1),(meanF+1.5*a1),(meanF+0.5*a1),(meanF-1.5*a1),(meanF-2.5*a1)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanF,sqrt(sigmaF),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanF[i],sqrt(sigmaF[i]))/dmixnorm(dataDH,meanF,sqrt(sigmaF),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    meanF<- sumwx/n0
    meanF[1]<- (sum(dataP1)+m_fam*sumwx[1])/(n_samP1+m_fam*n0[1])
    meanF[8]<- (sum(dataP2)+m_fam*sumwx[8])/(n_samP2+m_fam*n0[8])

    mean11<- meanF[1]
    mean12<- meanF[8]

    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanF[i])^2 }

    s0F<-matrix(0,2,1)
    s0F[1]<- ss1+ss2
    s0F[2]<- sum(swx)*m_fam
    s1<- sigmaF[1]-sigma/m_fam
    if (s1<0.0){s1<- 0.000001}
    abc2<- sigma

    abc3<-1000; n_iter<- 0
    while(abc3>0.0001){
      n_iter<- n_iter+1
      sigma<- (s0F[1]+s0F[2])/(n_samP1+n_samP2+n_samDH)
      abc3<- abs(abc2-sigma)
      abc2<- sigma
      if (n_iter>20) break
    }
    if (sigma<0.1*sigma0){ sigma<- 0.1*sigma0 }
    sigmaF<- matrix((sigma/m_fam),d2,1)

    ########criteria for iterations to stop#######

    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanF,sqrt(sigmaF),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*9

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1,1, 1,1,-1,-1,1,1,-1,-1, 1,-1,1,-1,1,-1,1,-1, 1,1,1,1,-1,-1,-1,-1,
                1,-1,-1,1,1,-1,-1,1, 1,1,-1,-1,-1,-1,1,1, 1,-1,1,-1,-1,1,-1,1, 1,-1,-1,1,-1,1,1,-1 ),8,8)
  b_line1 <- meanF
  B1 <- ginv(t(aa)%*%aa)%*%(t(aa)%*%b_line1)
  jj <- sigma_dh - sigmaF[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh
  mm<- sigmaF[1]-sigma[1]
  if(mm<0){mm<-0}
  nnn<- mm/sigma_dh

  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanF[i])/sqrt(sigmaF[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("3MG-AI",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanF),4)," "," "," "," "," "," "," "," ",round(sigmaF[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," ",round(B1[1],4)," "," ",round(B1[2],4),round(B1[3],4),round(B1[4],4)," ",round(B1[5],4),round(B1[6],4)," ",round(B1[7],4)," "," ",round(B1[8],4)," ",
                       round(jj,4),round(ll*100,4)," "," ",round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}

############################################ 3MG-A model#########################################
G3DHModelFun[[27]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(F-2)Model##########################################
  d2<-8
  mi <- as.matrix(c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125))
  meanF<-mean(dataDH)
  sigmaF <- matrix((sigma_dh/(2*1.2222101)),d2,1)
  sigma <- sigma0
  a1 <- sqrt(sigma_dh/n_samDH)
  if(mean11<mean12){a1=-a1}
  meanF <- as.matrix(c((meanF+3*a1),(meanF+2.1*a1),(meanF+1.2*a1),(meanF+0.3*a1),(meanF+1.5*a1),(meanF+0.5*a1),(meanF-1.5*a1),(meanF-2.5*a1)))

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanF,sqrt(sigmaF),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanF[i],sqrt(sigmaF[i]))/dmixnorm(dataDH,meanF,sqrt(sigmaF),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    hh<- matrix(0,4,4)
    hh[1,1]<- sigma*(1.0/(n_samP1+m_fam*n0[1])+1.0/(n_samP2+m_fam*n0[2]))+sigmaF[3]*(1.0/n0[3]+1.0/n0[6])
    hh[1,2]<- 0
    hh[1,3]<- sigma/(n_samP1+m_fam*n0[1])-sigmaF[6]/n0[6]
    hh[1,4]<- -sigmaF[3]/n0[3]+sigma/(n_samP2+m_fam*n0[8])

    hh[2,2]<- sigmaF[2]*(1.0/n0[2]+1.0/n0[4]+1.0/n0[5]+1.0/n0[7])
    hh[2,3]<- sigmaF[2]*(-1.0/n0[2]+1.0/n0[5])
    hh[2,4]<- sigmaF[4]*(1.0/n0[4]-1.0/n0[7])

    hh[3,3]<- sigma/(n_samP1+m_fam*n0[1])+sigmaF[2]*(1.0/n0[2]+1.0/n0[5]+1.0/n0[6])
    hh[3,4]<- 0

    hh[4,4]<- sigma/(n_samP2+m_fam*n0[8])+sigmaF[3]*(1.0/n0[3]+1.0/n0[4]+1.0/n0[7])
    for(i in 2:4){
      for(j in 1:(i-1)){
        hh[i,j]<- hh[j,i]
      }
    }
    b_line<- matrix(0,4,1)
    b_line[1]<-(sum(dataP1)+m_fam*sumwx[1])/(n_samP1+m_fam*n0[1])+(sum(dataP2)+m_fam*n0[8])/(n_samP2+m_fam*n0[8])-sumwx[3]/n0[3]-sumwx[6]/n0[6]
    b_line[2]<-sumwx[2]/n0[2]-sumwx[4]/n0[4]-sumwx[5]/n0[5]+sumwx[7]/n0[7]
    b_line[3]<-(sum(dataP1)+m_fam*sumwx[1])/(n_samP1+m_fam*n0[1])-sumwx[2]/n0[2]-sumwx[5]/n0[5]+sumwx[6]/n0[6]
    b_line[4]<-(sum(dataP2)+m_fam*sumwx[8])/(n_samP2+m_fam*n0[8])+sumwx[3]/n0[3]-sumwx[4]/n0[4]-sumwx[7]/n0[7]
    B <- solve(hh,b_line)

    meanF[1]<- (sum(dataP1)+m_fam*sumwx[1]-(B[1]+B[3])*sigma)/(n_samP1+m_fam*n0[1])
    meanF[2]<- (sumwx[2]+(-B[2]+B[3])*sigmaF[2])/n0[2]
    meanF[3]<- (sumwx[3]+(B[1]-B[4])*sigmaF[3])/n0[3]
    meanF[4]<- (sumwx[4]+(B[2]+B[4])*sigmaF[4])/n0[4]
    meanF[5]<- (sumwx[5]+(B[2]+B[3])*sigmaF[5])/n0[5]
    meanF[6]<- (sumwx[6]+(B[1]-B[3])*sigmaF[6])/n0[6]
    meanF[7]<- (sumwx[7]+(-B[2]+B[4])*sigmaF[7])/n0[7]
    meanF[8]<- (sum(dataP2)+m_fam*sumwx[8]-(B[1]+B[4])*sigma)/(n_samP2+m_fam*n0[8])

    mean11<- meanF[1]
    mean12<- meanF[8]

    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanF[i])^2 }

    s0F<-matrix(0,2,1)
    s0F[1]<- ss1+ss2
    s0F[2]<- sum(swx)*m_fam
    s1<- sigmaF[1]-sigma/m_fam
    if (s1<0.0){s1<- 0.000001}
    abc2<- sigma

    abc3<-1000; n_iter<- 0
    while(abc3>0.0001){
      n_iter<- n_iter+1
      sigma<- (s0F[1]+s0F[2])/(n_samP1+n_samP2+n_samDH)
      abc3<- abs(abc2-sigma)
      abc2<- sigma
      if (n_iter>20) break
    }
    if (sigma<0.1*sigma0){ sigma<- 0.1*sigma0 }
    sigmaF<- matrix((sigma/m_fam),d2,1)

    ########criteria for iterations to stop#######
    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanF,sqrt(sigmaF),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*5

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1,1, 1,1,-1,-1,1,1,-1,-1, 1,-1,1,-1,1,-1,1,-1, 1,1,1,1,-1,-1,-1,-1),8,4)
  b_line1 <- meanF
  B1 <- ginv(t(aa)%*%aa)%*%(t(aa)%*%b_line1)
  jj <- sigma_dh - sigmaF[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh
  mm<- sigmaF[1]-sigma[1]
  if(mm<0){mm<-0}
  nnn<- mm/sigma_dh

  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanF[i])/sqrt(sigmaF[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("3MG-A",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanF),4)," "," "," "," "," "," "," "," ",round(sigmaF[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," ",round(B1[1],4)," "," ",round(B1[2],4),round(B1[3],4),round(B1[4],4)," "," "," "," "," "," "," "," "," ",
                       round(jj,4),round(ll*100,4)," "," ",round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################ 3MG-CEA model#########################################
G3DHModelFun[[28]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(F-3)Model##########################################
  d2<-4
  mi <- as.matrix(c(0.125,0.375,0.375,0.125))
  meanF<-mean(dataDH)
  sigmaF <- matrix((sigma_dh/(2*1.2222101)),d2,1)
  sigma <- sigma0
  a1 <- sqrt(sigma_dh/n_samDH)
  if(mean11<mean12){a1=-a1}

  meanF <- as.matrix(c((meanF+3*a1),(meanF+a1),(meanF-a1),(meanF-3*a1)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanF,sqrt(sigmaF),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanF[i],sqrt(sigmaF[i]))/dmixnorm(dataDH,meanF,sqrt(sigmaF),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    aa1<- sigma*(1.0/n_samP1+1.0/n_samP2)+sigmaF[2]/n0[2]+sigmaF[3]/n0[3]
    aa2<- sigma*(1.0/n_samP1-1.0/n_samP2)+3.0*sigmaF[2]/n0[2]-3.0*sigmaF[3]/n0[3]
    aa3<- sigma*(1.0/n_samP1+1.0/n_samP2)+9.0*(sigmaF[2]/n0[2]+sigmaF[3]/n0[3])
    aa4<- sum(dataP1)/n_samP1+sum(dataP2)/n_samP2-sumwx[2]/n0[2]-sumwx[3]/n0[3]
    aa5<- sum(dataP1)/n_samP1-sum(dataP2)/n_samP2-3.0*sumwx[2]/n0[2]+3.0*sumwx[3]/n0[3]
    aa6<- aa1*aa3-aa2*aa2

    rr<- matrix(0,2,1)
    rr[1]<- (aa3*aa4-aa2*aa5)/aa6
    rr[2]<- (aa1*aa5-aa2*aa4)/aa6

    meanF[1]<- (sum(dataP1)+m_fam*sumwx[1]-(rr[1]+rr[2])*sigma)/(n_samP1+m_fam*n0[1])
    meanF[2]<- (sumwx[2]+(rr[1]+3.0*rr[2])*sigmaF[2])/n0[2]
    meanF[3]<- (sumwx[3]+(rr[1]-3.0*rr[2])*sigmaF[3])/n0[3]
    meanF[4]<- (sum(dataP2)+m_fam*sumwx[4]-(rr[1]-rr[2])*sigma)/(n_samP2+m_fam*n0[4])

    mean11<- meanF[1]
    mean12<- meanF[4]

    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanF[i])^2 }

    s0F<-matrix(0,2,1)
    s0F[1]<- ss1+ss2
    s0F[2]<- sum(swx)*m_fam
    s1<- sigmaF[1]-sigma/m_fam
    if (s1<0.0){s1<- 0.000001}
    abc2<- sigma

    abc3<-1000; n_iter<- 0
    while(abc3>0.0001){
      n_iter<- n_iter+1
      sigma<- (s0F[1]+s0F[2])/(n_samP1+n_samP2+n_samDH)
      abc3<- abs(abc2-sigma)
      abc2<- sigma
      if (n_iter>20) break
    }
    if (sigma<0.1*sigma0){ sigma<- 0.1*sigma0 }
    sigmaF<- matrix((sigma/m_fam),d2,1)

    ########criteria for iterations to stop#######
    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanF,sqrt(sigmaF),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*3
  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1, 3,1,-1,-3),4,2)
  b_line1 <- meanF
  B1 <- ginv(t(aa)%*%aa)%*%(t(aa)%*%b_line1)
  jj <- sigma_dh - sigmaF[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh
  mm<- sigmaF[1]-sigma[1]
  if(mm<0){mm<-0}
  nnn<- mm/sigma_dh

  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanF[i])/sqrt(sigmaF[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("3MG-CEA",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanF),4)," "," "," "," "," "," "," "," "," "," "," "," ",round(sigmaF[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4)," "," ",round(B1[2],4),round(B1[2],4),round(B1[2],4)," "," "," "," "," "," "," "," "," ",
                       round(jj,4),round(ll*100,4)," "," ",round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################ 3MG-PEA model#########################################
G3DHModelFun[[29]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(F-4)Model##########################################
  d2<-6
  mi <- as.matrix(c(0.125,0.125,0.25,0.25,0.125,0.125))
  meanF<-mean(dataDH)
  sigmaF <- matrix((sigma_dh/(2*1.2222101)),d2,1)
  sigma <- sigma0
  a1 <- sqrt(sigma_dh)
  if(mean11<mean12){a1=-a1}

  meanF <- as.matrix(c((meanF+3*a1),(meanF+2*a1),(meanF+a1),(meanF-a1),(meanF-2*a1),(meanF-3*a1)))
  #meanF <- as.matrix(c(290,110,190,10,90,-90))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanF,sqrt(sigmaF),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanF[i],sqrt(sigmaF[i]))/dmixnorm(dataDH,meanF,sqrt(sigmaF),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    hh<- matrix(0,3,3)
    hh[1,1]<- sigma*(1.0/(n_samP1+m_fam*n0[1])+1.0/(n_samP2+m_fam*n0[6]))+sigmaF[2]*(1.0/n0[2]+1.0/n0[5])
    hh[1,2]<- sigma/(n_samP1+m_fam*n0[1])+sigmaF[2]/n0[2]
    hh[1,3]<- sigma/(n_samP1+m_fam*n0[1])-sigmaF[5]/n0[5]
    hh[2,2]<- sigma/(n_samP1+m_fam*n0[1])+sigmaF[2]*(1.0/n0[2]+1.0/n0[3]+1.0/n0[4])
    hh[2,3]<- sigma/(n_samP1+m_fam*n0[1])+2.0*sigmaF[3]/n0[3]
    hh[3,3]<- sigma/(n_samP1+m_fam*n0[1])+4.0*sigmaF[3]/n0[3]+sigmaF[5]/n0[5]
    for(i in 2:3){
      for(j in 1:(i-1)){
        hh[i,j]<- hh[j,i]
      }
    }
    b_line<- matrix(0,3,1)
    b_line[1]<-(sum(dataP1)+m_fam*sumwx[1])/(n_samP1+m_fam*n0[1])+(sum(dataP2)+m_fam*sumwx[6])/(n_samP2+m_fam*n0[6])-sumwx[2]/n0[2]-sumwx[5]/n0[5]
    b_line[2]<-(sum(dataP1)+m_fam*sumwx[1])/(n_samP1+m_fam*n0[1])-sumwx[2]/n0[2]-sumwx[3]/n0[3]+sumwx[4]/n0[4]
    b_line[3]<-(sum(dataP1)+m_fam*sumwx[1])/(n_samP1+m_fam*n0[1])-2*sumwx[3]/n0[3]+sumwx[5]/n0[5]
    B <- solve(hh,b_line)

    meanF[1]<- (sum(dataP1)+m_fam*sumwx[1]-(B[1]+B[2]+B[3])*sigma)/(n_samP1+m_fam*n0[1])
    meanF[2]<- (sumwx[2]+(B[1]+B[2])*sigmaF[2])/n0[2]
    meanF[3]<- (sumwx[3]+(B[2]+2.0*B[3])*sigmaF[3])/n0[3]
    meanF[4]<- (sumwx[4]-B[2]*sigmaF[4])/n0[4]
    meanF[5]<- (sumwx[5]+(B[1]-B[3])*sigmaF[5])/n0[5]
    meanF[6]<- (sum(dataP2)+m_fam*sumwx[6]-B[1]*sigma)/(n_samP2+m_fam*n0[6])

    mean11<- meanF[1]
    mean12<- meanF[6]
    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanF[i])^2 }

    s0F<-matrix(0,2,1)
    s0F[1]<- ss1+ss2
    s0F[2]<- sum(swx)*m_fam
    s1<- sigmaF[1]-sigma/m_fam
    if (s1<0.0){s1<- 0.000001}
    abc2<- sigma

    abc3<-1000; n_iter<- 0
    while(abc3>0.0001){
      n_iter<- n_iter+1
      sigma<- (s0F[1]+s0F[2])/(n_samP1+n_samP2+n_samDH)
      abc3<- abs(abc2-sigma)
      abc2<- sigma
      if (n_iter>20) break
    }
    if (sigma<0.1*sigma0){ sigma<- 0.1*sigma0 }
    sigmaF<- matrix((sigma/m_fam),d2,1)

    ########criteria for iterations to stop#######
    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanF,sqrt(sigmaF),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*4

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1, 2,2,0,0,-2,-2, 1,-1,1,-1,1,-1),6,3)
  b_line1 <- meanF
  B1 <- ginv(t(aa)%*%aa)%*%(t(aa)%*%b_line1)
  jj <- sigma_dh - sigmaF[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh
  mm<- sigmaF[1]-sigma[1]
  if(mm<0){mm<-0}
  nnn<- mm/sigma_dh

  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanF[i])/sqrt(sigmaF[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("3MG-PEA",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanF),4)," "," "," "," "," "," "," "," "," "," ",round(sigmaF[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," ",round(B1[1],4)," "," ",round(B1[2],4),round(B1[2],4),round(B1[3],4)," "," "," "," "," "," "," "," "," ",
                       round(jj,4),round(ll*100,4)," "," ",round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}

############################################ MX3-AI-AI model#########################################
G3DHModelFun[[30]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(G-0)Model##########################################
  d2<-8
  mi <- as.matrix(c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125))
  meanG<-mean(dataDH)
  sigmaG <- matrix((sigma_dh/(2*1.2222101)),d2,1)
  sigma <- sigma0
  a1 <- sqrt(sigma_dh/n_samDH)
  if(mean11<mean12){a1=-a1}

  meanG <- as.matrix(c((meanG+3*a1),(meanG+2.1*a1),(meanG+1.2*a1),(meanG+0.3*a1),(meanG+1.5*a1),(meanG+0.5*a1),(meanG-1.5*a1),(meanG-2.5*a1)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanG,sqrt(sigmaG),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanG[i],sqrt(sigmaG[i]))/dmixnorm(dataDH,meanG,sqrt(sigmaG),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    mean11<- sum(dataP1)/n_samP1
    mean12<- sum(dataP2)/n_samP2
    meanG<- sumwx/n0
    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanG[i])^2 }

    s0G<-matrix(0,2,1)
    s0G[1]<- ss1+ss2
    s0G[2]<- sum(swx)*m_fam
    s1<- sigmaG[1]-sigma/m_fam    # variance of polygenes.
    if (s1<0.0){s1<- 0.000001}
    abc2<- sigma

    abc3<-1000; n_iter<- 0
    while(abc3>0.0001){
      n_iter<- n_iter+1
      abc1<- (sigma/m_fam)/(sigma/m_fam+s1)
      sigma<- (s0G[1]+abc1*abc1*s0G[2])/(n_samP1+n_samP2+abc1*n_samDH)
      abc3<- abs(abc2-sigma)
      abc2<- sigma
      if (n_iter>20) break
    }
    if (sigma<0.1*sigma0){ sigma<- 0.1*sigma0 }
    sigmaG[1]<- s1+sigma/m_fam

    ############CM3 variance of polygenes#############################
    s1<- sum(swx)/n_samDH-sigma/m_fam
    if (s1<0.0){ s1<- 0.000001 }
    sigmaG<- matrix((s1+sigma/m_fam),d2,1)

    ########criteria for iterations to stop#######

    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanG,sqrt(sigmaG),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*12
  #########first order genetic parameter process##########
  aa<- matrix(c(1,0,0,0,0,0,0,0,0,0, 0,1,0,0,0,0,0,0,0,0, 0,0,1,1,1,1,1,1,1,1, 1,-1,1,1,-1,-1,1,1,-1,-1, 1,-1,1,-1,1,-1,1,-1,1,-1,
                1,-1,1,1,1,1,-1,-1,-1,-1, 1,1,1,-1,-1,1,1,-1,-1,1, 1,1,1,1,-1,-1,-1,-1,1,1, 1,1,1,-1,1,-1,-1,1,-1,1, 1,-1,1,-1,-1,1,-1,1,1,-1 ),10,10)
  b_line1 <- matrix(c(mean11,mean12,meanG))
  B1 <- ginv(t(aa)%*%aa)%*%(t(aa)%*%b_line1)
  jj <- sigma_dh - sigmaG[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh
  mm<- sigmaG[1]-sigma[1]
  if(mm<0){mm<-0}
  nnn<- mm/sigma_dh

  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanG[i])/sqrt(sigmaG[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX3-AI-AI",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanG),4)," "," "," "," "," "," "," "," ",round(sigmaG[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," ",round(B1[1],4),round(B1[2],4),round(B1[3],4),round(B1[4],4),round(B1[5],4),round(B1[6],4)," ",round(B1[7],4),round(B1[8],4)," ",round(B1[9],4)," "," ",round(B1[10],4)," ",
                       round(jj,4),round(ll*100,4),round(mm,4),round(nnn*100,4),round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}

############################################ MX3-AI-A model#########################################
G3DHModelFun[[31]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(G-1)Model##########################################
  d2<-8
  mi <- as.matrix(c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125))
  meanG<-mean(dataDH)
  sigmaG <- matrix((sigma_dh/(2*1.2222101)),d2,1)
  sigma <- sigma0
  a1 <- sqrt(sigma_dh)
  if(mean11<mean12){a1=-a1}

  meanG <- as.matrix(c((meanG+3*a1),(meanG+2.1*a1),(meanG+1.2*a1),(meanG+0.3*a1),(meanG+1.5*a1),(meanG+0.5*a1),(meanG-1.5*a1),(meanG-2.5*a1)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanG,sqrt(sigmaG),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanG[i],sqrt(sigmaG[i]))/dmixnorm(dataDH,meanG,sqrt(sigmaG),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    aa1<- sum(dataP1)/n_samP1+sum(dataP2)/n_samP2-sumwx[1]/n0[1]-sumwx[8]/n0[8]
    aa2<- sigma*(1.0/n_samP1+1.0/n_samP2)+sigmaG[1]/n0[1]+sigmaG[8]/n0[8]
    rr<- aa1/aa2
    mean11<- (sum(dataP1)-rr*sigma)/n_samP1
    mean12<- (sum(dataP2)-rr*sigma)/n_samP2
    meanG<- sumwx/n0
    meanG[1]<-(sumwx[1]+rr*sigmaG[1])/n0[1]
    meanG[8]<-(sumwx[8]+rr*sigmaG[8])/n0[8]

    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanG[i])^2 }

    s0G<-matrix(0,2,1)
    s0G[1]<- ss1+ss2
    s0G[2]<- sum(swx)*m_fam
    s1<- sigmaG[1]-sigma/m_fam    # variance of polygenes.
    if (s1<0.0){s1<- 0.000001}
    abc2<- sigma; abc3<-1000; n_iter<- 0
    while(abc3>0.0001){
      n_iter<- n_iter+1
      abc1<- (sigma/m_fam)/(sigma/m_fam+s1)
      sigma<- (s0G[1]+abc1*abc1*s0G[2])/(n_samP1+n_samP2+abc1*n_samDH)
      abc3<- abs(abc2-sigma)
      abc2<- sigma
      if (n_iter>20) break
    }
    if (sigma<0.1*sigma0){ sigma<- 0.1*sigma0 }
    sigmaG[1]<- s1+sigma/m_fam

    ############CM3 variance of polygenes#############################
    s1<- sum(swx)/n_samDH-sigma/m_fam
    if (s1<0.0){ s1<- 0.000001 }
    sigmaG<- matrix((s1+sigma/m_fam),d2,1)

    ########criteria for iterations to stop#######

    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanG,sqrt(sigmaG),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*11
  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1,1,1,1, 1,-1,1,1,-1,-1,1,1,-1,-1, 1,-1,1,-1,1,-1,1,-1,1,-1, 1,-1,1,1,1,1,-1,-1,-1,-1, 1,1,1,-1,-1,1,1,-1,-1,1,
                1,1,1,1,-1,-1,-1,-1,1,1, 1,1,1,-1,1,-1,-1,1,-1,1, 1,-1,1,-1,-1,1,-1,1,1,-1, 1,-1,0,0,0,0,0,0,0,0),10,9)
  b_line1 <- matrix(c(mean11,mean12,meanG))
  B1 <- ginv(t(aa)%*%aa)%*%(t(aa)%*%b_line1)
  jj <- sigma_dh - sigmaG[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh
  mm<- sigmaG[1]-sigma[1]
  if(mm<0){mm<-0}
  nnn<- mm/sigma_dh

  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanG[i])/sqrt(sigmaG[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX3-AI-A",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanG),4)," "," "," "," "," "," "," "," ",round(sigmaG[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," ",round(B1[1],4)," "," ",round(B1[2],4),round(B1[3],4),round(B1[4],4)," ",round(B1[5],4),round(B1[6],4)," ",round(B1[7],4)," "," ",round(B1[8],4),round(B1[9],4),
                       round(jj,4),round(ll*100,4),round(mm,4),round(nnn*100,4),round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}

############################################ MX3-A-A model#########################################
G3DHModelFun[[32]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(G-2)Model##########################################
  d2<-8
  mi <- as.matrix(c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125))
  meanG<-mean(dataDH)
  sigmaG <- matrix((sigma_dh/(2*1.2222101)),d2,1)
  sigma <- sigma0
  a1 <- sqrt(sigma_dh)
  if(mean11<mean12){a1=-a1}

  meanG <- as.matrix(c((meanG+3*a1),(meanG+2.1*a1),(meanG+1.2*a1),(meanG+0.3*a1),(meanG+1.5*a1),(meanG+0.5*a1),(meanG-1.5*a1),(meanG-2.5*a1)))

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanG,sqrt(sigmaG),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanG[i],sqrt(sigmaG[i]))/dmixnorm(dataDH,meanG,sqrt(sigmaG),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    ############solve the linear equation############
    hh<-matrix(0,5,5)
    hh[1,1]<- sigma*(1.0/n_samP1+1.0/n_samP2)+sigmaG[1]*(1.0/n0[1]+1.0/n0[8])
    hh[1,2]<- -sigmaG[1]*(1.0/n0[1]+1.0/n0[8])
    hh[1,3]<- 0
    hh[1,4]<- -sigmaG[1]/n0[1]
    hh[1,5]<- -sigmaG[8]/n0[8]
    hh[2,2]<- sigmaG[1]*(1.0/n0[1]+1.0/n0[3]+1.0/n0[6]+1.0/n0[8])
    hh[2,3]<- 0
    hh[2,4]<- sigmaG[1]*(1.0/n0[1]-1.0/n0[6])
    hh[2,5]<- sigmaG[1]*(-1.0/n0[3]+1.0/n0[8])
    hh[3,3]<- sigmaG[2]*(1.0/n0[2]+1.0/n0[4]+1.0/n0[5]+1.0/n0[7])
    hh[3,4]<- sigmaG[2]*(-1.0/n0[2]+1.0/n0[5])
    hh[3,5]<- sigmaG[4]*(1.0/n0[4]-1.0/n0[7])
    hh[4,4]<- sigmaG[1]*(1.0/n0[1]+1.0/n0[2]+1.0/n0[5]+1.0/n0[6])
    hh[4,5]<- 0
    hh[5,5]<- sigmaG[3]*(1.0/n0[3]+1.0/n0[4]+1.0/n0[7]+1.0/n0[8])
    for(i in 2:5){
      for(j in 1:(i-1)){
        hh[i,j]<- hh[j,i]
      }
    }
    b_line<-matrix(0,5,1)
    b_line[1]<- sum(dataP1)/n_samP1+sum(dataP2)/n_samP2-sumwx[1]/n0[1]-sumwx[8]/n0[8]
    b_line[2]<- sumwx[1]/n0[1]-sumwx[3]/n0[3]-sumwx[6]/n0[6]+sumwx[8]/n0[8]
    b_line[3]<- sumwx[2]/n0[2]-sumwx[4]/n0[4]-sumwx[5]/n0[5]+sumwx[7]/n0[7]
    b_line[4]<- sumwx[1]/n0[1]-sumwx[2]/n0[2]-sumwx[5]/n0[5]+sumwx[6]/n0[6]
    b_line[5]<- sumwx[3]/n0[3]-sumwx[4]/n0[4]-sumwx[7]/n0[7]+sumwx[8]/n0[8]
    B <- solve(hh,b_line)

    mean11<- (sum(dataP1)-B[1]*sigma)/n_samP1
    mean12<- (sum(dataP2)-B[1]*sigma)/n_samP2
    meanG[1]<- (sumwx[1]+(B[1]-B[2]-B[4])*sigmaG[1])/n0[1]
    meanG[2]<- (sumwx[2]+(-B[3]+B[4])*sigmaG[2])/n0[2]
    meanG[3]<- (sumwx[3]+(B[2]-B[5])*sigmaG[3])/n0[3]
    meanG[4]<- (sumwx[4]+(B[3]+B[5])*sigmaG[4])/n0[4]
    meanG[5]<- (sumwx[5]+(B[3]+B[4])*sigmaG[5])/n0[5]
    meanG[6]<- (sumwx[6]+(B[2]-B[4])*sigmaG[6])/n0[6]
    meanG[7]<- (sumwx[7]-(B[3]-B[5])*sigmaG[7])/n0[7]
    meanG[8]<- (sumwx[8]+(B[1]-B[2]-B[5])*sigmaG[8])/n0[8]

    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanG[i])^2 }

    s0G<-matrix(0,2,1)
    s0G[1]<- ss1+ss2
    s0G[2]<- sum(swx)*m_fam
    s1<- sigmaG[1]-sigma/m_fam    # variance of polygenes.
    if (s1<0.0){s1<- 0.000001}
    abc2<- sigma

    abc3<-1000; n_iter<- 0
    while(abc3>0.0001){
      n_iter<- n_iter+1
      abc1<- (sigma/m_fam)/(sigma/m_fam+s1)
      sigma<- (s0G[1]+abc1*abc1*s0G[2])/(n_samP1+n_samP2+n_samDH)
      abc3<- abs(abc2-sigma)
      abc2<- sigma
      if (n_iter>20) break
    }
    if (sigma<0.1*sigma0){ sigma<- 0.1*sigma0 }
    sigmaG[1]<- s1+sigma/m_fam

    ############CM3 variance of polygenes#############################
    s1<- sum(swx)/n_samDH-sigma/m_fam
    if (s1<0.0){ s1<- 0.000001 }
    sigmaG<- matrix((s1+sigma/m_fam),d2,1)

    ########criteria for iterations to stop#######
    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanG,sqrt(sigmaG),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*6
  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1,1,1,1, 1,-1,1,1,-1,-1,1,1,-1,-1, 1,-1,1,-1,1,-1,1,-1,1,-1,
                1,-1,1,1,1,1,-1,-1,-1,-1, 1,-1,0,0,0,0,0,0,0,0),10,5)
  b_line1 <- matrix(c(mean11,mean12,meanG))
  B1 <- ginv(t(aa)%*%aa)%*%(t(aa)%*%b_line1)
  jj <- sigma_dh - sigmaG[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh
  mm<- sigmaG[1]-sigma[1]
  if(mm<0){mm<-0}
  nnn<- mm/sigma_dh

  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanG[i])/sqrt(sigmaG[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX3-A-A",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanG),4)," "," "," "," "," "," "," "," ",round(sigmaG[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," ",round(B1[1],4)," "," ",round(B1[2],4),round(B1[3],4),round(B1[4],4)," "," "," "," "," "," "," "," ",round(B1[5],4),
                       round(jj,4),round(ll*100,4),round(mm,4),round(nnn*100,4),round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################ MX3-CEA-A model#########################################
G3DHModelFun[[33]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(G-3)Model##########################################
  d2<-4
  mi <- as.matrix(c(0.125,0.375,0.375,0.125))
  meanG<-mean(dataDH)
  sigmaG <- matrix((sigma_dh/(2*1.2222101)),d2,1)
  sigma <- sigma0
  a1 <- sqrt(sigma_dh)
  if(mean11<mean12){a1=-a1}

  meanG <- as.matrix(c((meanG+3*a1),(meanG+a1),(meanG-a1),(meanG-3*a1)))

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanG,sqrt(sigmaG),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanG[i],sqrt(sigmaG[i]))/dmixnorm(dataDH,meanG,sqrt(sigmaG),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    hh<-matrix(0,3,3)
    hh[1,1]<- sigma*(1.0/n_samP1+1.0/n_samP2)+sigmaG[1]*(1.0/n0[1]+1.0/n0[4])
    hh[1,2]<- -sigmaG[1]*(1.0/n0[1]+1.0/n0[4])
    hh[1,3]<- -sigmaG[1]*(1.0/n0[1]-1.0/n0[4])
    hh[2,2]<- sigmaG[1]*(1.0/n0[1]+1.0/n0[2]+1.0/n0[3]+1.0/n0[4])
    hh[2,3]<- sigmaG[1]*(1.0/n0[1]+3.0/n0[2]-3.0/n0[3]-1.0/n0[4])
    hh[3,3]<- sigmaG[1]*(1.0/n0[1]+9.0/n0[2]+9.0/n0[3]+1.0/n0[4])
    for(i in 2:3){
      for(j in 1:(i-1)){
        hh[i,j]<- hh[j,i]
      }
    }
    b_line<-matrix(0,3,1)
    b_line[1]<- sum(dataP1)/n_samP1+sum(dataP2)/n_samP2-sumwx[1]/n0[1]-sumwx[4]/n0[4]
    b_line[2]<- sumwx[1]/n0[1]-sumwx[2]/n0[2]-sumwx[3]/n0[3]+sumwx[4]/n0[4]
    b_line[3]<- sumwx[1]/n0[1]-3.0*sumwx[2]/n0[2]+3.0*sumwx[3]/n0[3]-sumwx[4]/n0[4]
    B <- solve(hh,b_line)

    mean11<- (sum(dataP1)-B[1]*sigma)/n_samP1
    mean12<- (sum(dataP2)-B[1]*sigma)/n_samP2

    meanG[1]<-(sumwx[1]+(B[1]-B[2]-B[3])*sigmaG[1])/n0[1]
    meanG[2]<-(sumwx[2]+(B[2]+3.0*B[3])*sigmaG[2])/n0[2]
    meanG[3]<-(sumwx[3]+(B[2]-3.0*B[3])*sigmaG[3])/n0[3]
    meanG[4]<-(sumwx[4]+(B[1]-B[2]+B[3])*sigmaG[4])/n0[4]

    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanG[i])^2 }

    s0G<-matrix(0,2,1)
    s0G[1]<- ss1+ss2
    s0G[2]<- sum(swx)*m_fam
    s1<- sigmaG[1]-sigma/m_fam    # variance of polygenes.
    if (s1<0.0){s1<- 0.000001}
    abc2<- sigma

    abc3<-1000; n_iter<- 0
    while(abc3>0.0001){
      n_iter<- n_iter+1
      abc1<- (sigma/m_fam)/(sigma/m_fam+s1)
      sigma<- (s0G[1]+abc1*abc1*s0G[2])/(n_samP1+n_samP2+abc1*n_samDH)
      abc3<- abs(abc2-sigma)
      abc2<- sigma
      if (n_iter>20) break
    }
    if (sigma<0.1*sigma0){ sigma<- 0.1*sigma0 }
    sigmaG[1]<- s1+sigma/m_fam

    ############CM3 variance of polygenes#############################
    s1<- sum(swx)/n_samDH-sigma/m_fam
    if (s1<0.0){ s1<- 0.000001 }
    sigmaG<- matrix((s1+sigma/m_fam),d2,1)

    ########criteria for iterations to stop#######
    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanG,sqrt(sigmaG),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*5
  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,3,-3,3,1,-1,-3, 1,-1,0,0,0,0),6,3)
  b_line1 <- matrix(c(mean11,mean12,meanG))
  B1 <- ginv(t(aa)%*%aa)%*%(t(aa)%*%b_line1)
  jj <- sigma_dh - sigmaG[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh
  mm<- sigmaG[1]-sigma[1]
  if(mm<0){mm<-0}
  nnn<- mm/sigma_dh

  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanG[i])/sqrt(sigmaG[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX3-CEA-A",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanG),4)," "," "," "," "," "," "," "," "," "," "," "," ",round(sigmaG[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4)," "," ",round(B1[2],4),round(B1[2],4),round(B1[2],4)," "," "," "," "," "," "," "," ",round(B1[3],4),
                       round(jj,4),round(ll*100,4),round(mm,4),round(nnn*100,4),round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################ MX3-PEA-A model#########################################
G3DHModelFun[[34]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(G-4)Model##########################################
  d2<- 6
  mi <- as.matrix(c(0.125,0.125,0.25,0.25,0.125,0.125))
  meanG<-mean(dataDH)
  sigmaG <- matrix((sigma_dh/(2*1.2222101)),d2,1)
  sigma <- sigma0
  a1 <- sqrt(sigma_dh/n_samDH)
  if(mean11<mean12){a1=-a1}

  meanG <- as.matrix(c((meanG+3*a1),(meanG+2*a1),(meanG+a1),(meanG-a1),(meanG-2*a1),(meanG-3*a1)))

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanG,sqrt(sigmaG),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanG[i],sqrt(sigmaG[i]))/dmixnorm(dataDH,meanG,sqrt(sigmaG),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    hh<-matrix(0,4,4)
    hh[1,1]<- sigma*(1.0/n_samP1+1.0/n_samP2)+sigmaG[1]*(1.0/n0[1]+1.0/n0[6])
    hh[1,2]<- -sigmaG[1]*(1.0/n0[1]+1.0/n0[6])
    hh[1,3]<- -sigmaG[1]/n0[1]
    hh[1,4]<- -sigmaG[1]/n0[1]
    hh[2,2]<- sigmaG[1]*(1.0/n0[1]+1.0/n0[2]+1.0/n0[5]+1.0/n0[6])
    hh[2,3]<- sigmaG[1]*(1.0/n0[1]+1.0/n0[2])
    hh[2,4]<- sigmaG[1]*(1.0/n0[1]-1.0/n0[5])
    hh[3,3]<- sigmaG[1]*(1.0/n0[1]+1.0/n0[2]+1.0/n0[3]+1.0/n0[4])
    hh[3,4]<- sigmaG[1]*(1.0/n0[1]+2.0/n0[3])
    hh[4,4]<- sigmaG[1]*(1.0/n0[1]+4.0/n0[3]+1.0/n0[5])
    for(i in 2:4){
      for(j in 1:(i-1)){
        hh[i,j]<- hh[j,i]
      }
    }
    b_line<-matrix(0,4,1)
    b_line[1]<- sum(dataP1)/n_samP1+sum(dataP2)/n_samP2-sumwx[1]/n0[1]-sumwx[6]/n0[6]
    b_line[2]<- sumwx[1]/n0[1]-sumwx[2]/n0[2]-sumwx[5]/n0[5]+sumwx[6]/n0[6]
    b_line[3]<- sumwx[1]/n0[1]-sumwx[2]/n0[2]-sumwx[3]/n0[3]+sumwx[4]/n0[4]
    b_line[4]<- sumwx[1]/n0[1]-2.0*sumwx[3]/n0[3]+sumwx[5]/n0[5]
    B <- solve(hh,b_line)

    mean11<- (sum(dataP1)-B[1]*sigma)/n_samP1
    mean12<- (sum(dataP2)-B[1]*sigma)/n_samP2

    meanG[1]<-(sumwx[1]-(-B[1]+B[2]+B[3]+B[4])*sigmaG[1])/n0[1]
    meanG[2]<-(sumwx[2]+(B[2]+B[3])*sigmaG[2])/n0[2]
    meanG[3]<-(sumwx[3]+(B[3]+2.0*B[4])*sigmaG[3])/n0[3]
    meanG[4]<-(sumwx[4]-B[3]*sigmaG[4])/n0[4]
    meanG[5]<-(sumwx[5]+(B[2]-B[4])*sigmaG[5])/n0[5]
    meanG[6]<-(sumwx[6]+(B[1]-B[2])*sigmaG[6])/n0[6]

    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanG[i])^2 }

    s0G<-matrix(0,2,1)
    s0G[1]<- ss1+ss2
    s0G[2]<- sum(swx)*m_fam
    s1<- sigmaG[1]-sigma/m_fam    # variance of polygenes.
    if (s1<0.0){s1<- 0.000001}
    abc2<- sigma

    abc3<-1000; n_iter<- 0
    while(abc3>0.0001){
      n_iter<- n_iter+1
      abc1<- (sigma/m_fam)/(sigma/m_fam+s1)
      sigma<- (s0G[1]+abc1*abc1*s0G[2])/(n_samP1+n_samP2+n_samDH)
      abc3<- abs(abc2-sigma)
      abc2<- sigma
      if (n_iter>20) break
    }
    if (sigma<0.1*sigma0){ sigma<- 0.1*sigma0 }
    sigmaG[1]<- s1+sigma/m_fam

    ############CM3 variance of polygenes#############################
    s1<- sum(swx)/n_samDH-sigma/m_fam
    if (s1<0.0){ s1<- 0.000001 }
    sigmaG<- matrix((s1+sigma/m_fam),d2,1)

    ########criteria for iterations to stop#######
    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanG,sqrt(sigmaG),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*6
  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1,1, 2,-2,2,2,0,0,-2,-2, 1,-1,1,-1,1,-1,1,-1, 1,-1,0,0,0,0,0,0),8,4)
  b_line1 <- matrix(c(mean11,mean12,meanG))
  B1 <- ginv(t(aa)%*%aa)%*%(t(aa)%*%b_line1)
  jj <- sigma_dh - sigmaG[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh
  mm<- sigmaG[1]-sigma[1]
  if(mm<0){mm<-0}
  nnn<- mm/sigma_dh

  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanG[i])/sqrt(sigmaG[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX3-PEA-A",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanG),4)," "," "," "," "," "," "," "," "," "," ",round(sigmaG[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," ",round(B1[1],4)," "," ",round(B1[2],4),round(B1[3],4),","," "," "," "," "," "," "," "," ",round(B1[4],4),
                       round(jj,4),round(ll*100,4),round(mm,4),round(nnn*100,4),round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################ 4MG-AI model#########################################
G3DHModelFun[[35]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(H-1)Model##########################################
  d2<-16
  mi <- as.matrix(c(0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625))
  meanH<-mean(dataDH)
  sigmaH <- matrix((sigma_dh/(2*1.2222101)),d2,1)
  sigma <- sigma0
  a1 <- sqrt(sigma_dh)
  if(mean11<mean12){a1=-a1}

  meanH <- as.matrix(c(222,146,114,138,114,78,50,54,152,96,76,100,84,68,72,56))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanH,sqrt(sigmaH),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanH[i],sqrt(sigmaH[i]))/dmixnorm(dataDH,meanH,sqrt(sigmaH),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    s0<-matrix(0,4,1)
    s0[1]<- sigma;s0[2]<- n_samP1+m_fam*n0[1]
    s0[3]<- sigma;s0[4]<- n_samP2+m_fam*n0[16]

    hh<-matrix(0,5,5)
    hh[1,1]<- s0[1]/s0[2]+sigmaH[2]/n0[2]+sigmaH[3]/n0[3]+sigmaH[4]/n0[4]+sigmaH[9]/n0[9]+sigmaH[10]/n0[10]+sigmaH[11]/n0[11]+sigmaH[12]/n0[12]
    hh[1,2]<- s0[1]/s0[2]+sigmaH[2]/n0[2]+sigmaH[3]/n0[3]+sigmaH[4]/n0[4]
    hh[1,3]<- s0[1]/s0[2]+sigmaH[2]/n0[2]+sigmaH[3]/n0[3]+sigmaH[4]/n0[4]
    hh[1,4]<- -(s0[1]/s0[2]+sigmaH[2]/n0[2]+sigmaH[9]/n0[9]+sigmaH[10]/n0[10])
    hh[1,5]<- sigmaH[2]/n0[2]+sigmaH[3]/n0[3]+sigmaH[10]/n0[10]+sigmaH[11]/n0[11]

    hh[2,2]<- s0[1]/s0[2]+sigmaH[2]/n0[2]+sigmaH[3]/n0[3]+sigmaH[4]/n0[4]+sigmaH[13]/n0[13]+sigmaH[14]/n0[14]+sigmaH[15]/n0[15]+s0[3]/s0[4]
    hh[2,3]<- s0[1]/s0[2]+sigmaH[2]/n0[2]+sigmaH[3]/n0[3]+sigmaH[4]/n0[4]
    hh[2,4]<- -(s0[1]/s0[2]+sigmaH[2]/n0[2]-sigmaH[13]/n0[13]-sigmaH[14]/n0[14])
    hh[2,5]<- sigmaH[2]/n0[2]+sigmaH[3]/n0[3]+sigmaH[13]/n0[13]+s0[3]/s0[4]

    hh[3,3]<- s0[1]/s0[2]+sigmaH[2]/n0[2]+sigmaH[3]/n0[3]+sigmaH[4]/n0[4]+sigmaH[5]/n0[5]+sigmaH[6]/n0[6]+sigmaH[7]/n0[7]+sigmaH[8]/n0[8]
    hh[3,4]<- -(s0[1]/s0[2]+sigmaH[2]/n0[2]+sigmaH[5]/n0[5]+sigmaH[6]/n0[6])
    hh[3,5]<- sigmaH[2]/n0[2]+sigmaH[3]/n0[3]-sigmaH[5]/n0[5]-sigmaH[8]/n0[8]

    hh[4,4]<- s0[1]/s0[2]+sigmaH[2]/n0[2]+sigmaH[5]/n0[5]+sigmaH[6]/n0[6]+sigmaH[9]/n0[9]+sigmaH[10]/n0[10]+sigmaH[13]/n0[13]+sigmaH[14]/n0[14]
    hh[4,5]<- -sigmaH[2]/n0[2]+sigmaH[5]/n0[5]-sigmaH[10]/n0[10]+sigmaH[13]/n0[13]

    hh[5,5]<- sigmaH[2]/n0[2]+sigmaH[3]/n0[3]+sigmaH[5]/n0[5]+sigmaH[8]/n0[8]+sigmaH[10]/n0[10]+sigmaH[11]/n0[11]+sigmaH[13]/n0[13]+s0[3]/s0[4]
    for(i in 2:5){
      for(j in 1:(i-1)){
        hh[i,j]<- hh[j,i]
      }
    }
    ##############################################################################################################################
    b_line<-matrix(0,5,1)
    b_line[1]<- -(sum(dataP1)+m_fam*sumwx[1])/(n_samP1+m_fam*n0[1])+sumwx[2]/n0[2]-sumwx[3]/n0[3]+sumwx[4]/n0[4]+sumwx[9]/n0[9]-sumwx[10]/n0[10]+sumwx[11]/n0[11]-sumwx[12]/n0[12]
    b_line[2]<- -(sum(dataP1)+m_fam*sumwx[1])/(n_samP1+m_fam*n0[1])+sumwx[2]/n0[2]-sumwx[3]/n0[3]+sumwx[4]/n0[4]+sumwx[13]/n0[13]-sumwx[14]/n0[14]+sumwx[15]/n0[15]-(sum(dataP2)+m_fam*sumwx[16])/(n_samP2+m_fam*n0[16])
    b_line[3]<- -(sum(dataP1)+m_fam*sumwx[1])/(n_samP1+m_fam*n0[1])+sumwx[2]/n0[2]-sumwx[3]/n0[3]+sumwx[4]/n0[4]+sumwx[5]/n0[5]-sumwx[6]/n0[6]+sumwx[7]/n0[7]-sumwx[8]/n0[8]
    b_line[4]<- (sum(dataP1)+m_fam*sumwx[1])/(n_samP1+m_fam*n0[1])-sumwx[2]/n0[2]-sumwx[5]/n0[5]+(sum(dataP2)+m_fam*sumwx[16])/(n_samP2+m_fam*n0[16])-sumwx[9]/n0[9]+sumwx[10]/n0[10]+sumwx[13]/n0[13]-sumwx[14]/n0[14];
    b_line[5]<- sumwx[2]/n0[2]-sumwx[3]/n0[3]-sumwx[5]/n0[5]+sumwx[8]/n0[8]-sumwx[10]/n0[10]+sumwx[11]/n0[11]+sumwx[13]/n0[13]-(sum(dataP2)+m_fam*sumwx[16])/(n_samP2+m_fam*n0[16])
    B <- solve(hh,b_line)
    ###############################################################################################################################
    meanH[1]<-(sum(dataP1)+m_fam*sumwx[1]+(B[1]+B[2]+B[3]-B[4])*sigma)/(n_samP1+m_fam*n0[1])
    meanH[2]<-(sumwx[2]-(B[1]+B[2]+B[3]-B[4]+B[5])*sigmaH[2])/n0[2]
    meanH[3]<-(sumwx[3]+(B[1]+B[2]+B[3]+B[5])*sigmaH[3])/n0[3]
    meanH[4]<-(sumwx[4]-(B[1]+B[2]+B[3])*sigmaH[4])/n0[4]
    meanH[5]<-(sumwx[5]+(-B[3]+B[4]+B[5])*sigmaH[5])/n0[5]
    meanH[6]<-(sumwx[6]+(B[3]-B[4])*sigmaH[6])/n0[6]
    meanH[7]<-(sumwx[7]+(-B[3])*sigmaH[7])/n0[7]
    meanH[8]<-(sumwx[8]+(B[3]-B[5])*sigmaH[8])/n0[8]

    meanH[9]<-(sumwx[9]-(B[1]-B[4])*sigmaH[9])/n0[9]
    meanH[10]<-(sumwx[10]+(B[1]-B[4]+B[5])*sigmaH[10])/n0[10]
    meanH[11]<-(sumwx[11]-(B[1]+B[5])*sigmaH[11])/n0[11]
    meanH[12]<-(sumwx[12]+B[1]*sigmaH[12])/n0[12]
    meanH[13]<-(sumwx[13]-(B[2]+B[4]+B[5])*sigmaH[13])/n0[13]
    meanH[14]<-(sumwx[14]+(B[2]+B[4])*sigmaH[14])/n0[14]
    meanH[15]<-(sumwx[15]-B[2]*sigmaH[15])/n0[15]
    meanH[16]<-(sum(dataP2)+m_fam*sumwx[16]+(B[2]+B[5])*sigma)/(n_samP2+m_fam*n0[16])

    mean11<- meanH[1]
    mean12<- meanH[16]
    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanH[i])^2 }

    s0H<-matrix(0,2,1)
    s0H[1]<- ss1+ss2
    s0H[2]<- sum(swx)*m_fam
    s1<- sigmaH[1]-sigma/m_fam    # variance of polygenes.
    if (s1<0.0){s1<- 0.000001}
    abc2<- sigma

    abc3<-1000; n_iter<- 0
    while(abc3>0.0001){
      n_iter<- n_iter+1
      sigma<- (s0H[1]+s0H[2])/(n_samP1+n_samP2+n_samDH)
      abc3<- abs(abc2-sigma)
      abc2<- sigma
      if (n_iter>20) break
    }
    if (sigma<0.1*sigma0){ sigma<- 0.1*sigma0 }
    sigmaH<-matrix((sigma/m_fam),d2,1)

    ########criteria for iterations to stop#######
    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanH,sqrt(sigmaH),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*11
  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1, 1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,
                1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1, 1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1, 1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,
                1,1,-1,-1,1,1,-1,-1,-1,-1,1,1,-1,-1,1,1, 1,-1,-1,1,1,-1,-1,1,-1,1,1,-1,-1,1,1,-1, 1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,
                1,-1,1,1,-1,1,1,-1,1,-1,-1,1,-1,1,1,-1, 1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1),16,11)
  b_line1 <- meanH
  B1 <- ginv(t(aa)%*%aa)%*%(t(aa)%*%b_line1)
  jj <- sigma_dh - sigmaH[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh

  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanH[i])/sqrt(sigmaH[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("4MG-AI",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanH),4),round(sigmaH[1],4),round(t(mix_pi),4),
                       round(B1[1],4)," "," ",round(B1[2],4),round(B1[3],4),round(B1[4],4),round(B1[5],4),round(B1[6],4),round(B1[7],4),round(B1[8],4),round(B1[9],4),round(B1[10],4),round(B1[11],4)," "," ",
                       round(jj,4),round(ll*100,4)," "," ",round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}


############################################ 4MG-CEA model#########################################
G3DHModelFun[[36]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(H-3)Model##########################################
  d2<-5
  mi <- as.matrix(c(0.0625,0.25,0.375,0.25,0.0625))
  meanH<-mean(dataDH)
  sigmaH <- matrix((sigma_dh/(2*1.2222101)),d2,1)
  sigma <- sigma0
  a1 <- sqrt(sigma_dh/n_samDH)
  if(mean11<mean12){a1=-a1}

  meanH <- as.matrix(c((meanH+3*a1),(meanH+2*a1),meanH,(meanH-2*a1),(meanH-3*a1)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanH,sqrt(sigmaH),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanH[i],sqrt(sigmaH[i]))/dmixnorm(dataDH,meanH,sqrt(sigmaH),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    s0<-matrix(0,4,1)
    s0[1]<- sigma;s0[2]<- n_samP1+m_fam*n0[1]
    s0[3]<- sigma;s0[4]<- n_samP2+m_fam*n0[5]

    hh<-matrix(0,3,3)
    hh[1,1]<- s0[1]/s0[2]+4.0*sigmaH[2]/n0[2]+sigmaH[3]/n0[3]
    hh[1,2]<- 2.0*s0[1]/s0[2]+6.0*sigmaH[2]/n0[2]
    hh[1,3]<- s0[1]/s0[2]+2.0*sigmaH[2]/n0[2]

    hh[2,2]<- 4.0*s0[1]/s0[2]+9.0*sigmaH[2]/n0[2]+sigmaH[4]/n0[4]
    hh[2,3]<- 2.0*s0[1]/s0[2]+3.0*sigmaH[2]/n0[2]-sigmaH[4]/n0[4]

    hh[3,3]<- s0[1]/s0[2]+sigmaH[2]/n0[2]+sigmaH[4]/n0[4]+s0[3]/s0[4]
    for(i in 2:3){
      for(j in 1:(i-1)){
        hh[i,j]<- hh[j,i]
      }
    }
    ##############################################################################################################################
    b_line<-matrix(0,3,1)
    b_line[1]<- -(sum(dataP1)+m_fam*sumwx[1])/(n_samP1+m_fam*n0[1])+2*sumwx[2]/n0[2]-sumwx[3]/n0[3]
    b_line[2]<- -2*(sum(dataP1)+m_fam*sumwx[1])/(n_samP1+m_fam*n0[1])+3*sumwx[2]/n0[2]-sumwx[4]/n0[4]
    b_line[3]<- -(sum(dataP1)+m_fam*sumwx[1])/(n_samP1+m_fam*n0[1])+sumwx[2]/n0[2]+sumwx[4]/n0[4]-(sum(dataP2)+m_fam*sumwx[5])/(n_samP2+m_fam*n0[5])
    B <- solve(hh,b_line)
    ###############################################################################################################################
    meanH[1]<-(sum(dataP1)+m_fam*sumwx[1]+(B[1]+2*B[2]+B[3])*sigma)/(n_samP1+m_fam*n0[1])
    meanH[2]<-(sumwx[2]+(-2*B[1]-3*B[2]-B[3])*sigmaH[2])/n0[2]
    meanH[3]<-(sumwx[3]+B[1]*sigmaH[3])/n0[3]
    meanH[4]<-(sumwx[4]+(B[2]-B[3])*sigmaH[4])/n0[4]
    meanH[5]<-(sum(dataP2)+m_fam*sumwx[5]+B[3]*sigma)/(n_samP2+m_fam*n0[5])

    mean11<- meanH[1]
    mean12<- meanH[5]
    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanH[i])^2 }

    s0H<-matrix(0,2,1)
    s0H[1]<- ss1+ss2
    s0H[2]<- sum(swx)*m_fam
    s1<- sigmaH[1]-sigma/m_fam    # variance of polygenes.
    if (s1<0.0){s1<- 0.000001}
    abc2<- sigma

    abc3<-1000; n_iter<- 0
    while(abc3>0.0001){
      n_iter<- n_iter+1
      sigma<- (s0H[1]+s0H[2])/(n_samP1+n_samP2+n_samDH)
      abc3<- abs(abc2-sigma)
      abc2<- sigma
      if (n_iter>20) break
    }
    if (sigma<0.1*sigma0){ sigma<- 0.1*sigma0 }
    sigmaH<-matrix((sigma/m_fam),d2,1)

    ########criteria for iterations to stop#######
    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanH,sqrt(sigmaH),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*3
  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1, 4,2,0,-2,-4),5,2)
  b_line1 <- meanH
  B1 <- ginv(t(aa)%*%aa)%*%(t(aa)%*%b_line1)
  jj <- sigma_dh - sigmaH[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh
  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanH[i])/sqrt(sigmaH[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("4MG-CEA",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanH),4)," "," "," "," "," "," "," "," "," "," "," ",round(sigmaH[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," "," "," "," ",round(B1[1],4)," "," ",round(B1[2],4),round(B1[2],4),round(B1[2],4),round(B1[2],4)," "," "," "," "," "," "," "," ",
                       round(jj,4),round(ll*100,4)," "," ",round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################################ 4MG-EEA model#########################################
G3DHModelFun[[37]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(H-4)Model##########################################
  d2<-9
  mi <- as.matrix(c(0.0625,0.0625,0.125,0.125,0.25,0.125,0.125,0.0625,0.0625))
  meanH<-mean(dataDH)
  sigmaH <- matrix((sigma_dh/(2*1.2222101)),d2,1)
  sigma <- sigma0
  a1 <- sqrt(sigma_dh/n_samDH)
  if(mean11<mean12){a1=-a1}

  meanH <- as.matrix(c((meanH+3*a1),(meanH+2.5*a1),(meanH+2*a1),(meanH+1.5*a1),(meanH+a1),
                       (meanH-1.5*a1),(meanH-2*a1),(meanH-2.5*a1),(meanH-3*a1)))

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanH,sqrt(sigmaH),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanH[i],sqrt(sigmaH[i]))/dmixnorm(dataDH,meanH,sqrt(sigmaH),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    s0<-matrix(0,4,1)
    s0[1]<- sigma;s0[2]<- n_samP1+m_fam*n0[1]
    s0[3]<- sigma;s0[4]<- n_samP2+m_fam*n0[9]

    hh<-matrix(0,6,6)
    hh[1,1]<- s0[1]/s0[2]+sigmaH[2]/n0[2]+sigmaH[4]/n0[4]+sigmaH[6]/n0[6]
    hh[1,2]<- s0[1]/s0[2]
    hh[1,3]<- 0
    hh[1,4]<- sigmaH[4]/n0[4]-sigmaH[6]/n0[6]
    hh[1,5]<- -(s0[1]/s0[2]+sigmaH[4]/n0[4])
    hh[1,6]<- -(s0[1]/s0[2]+2.0*sigmaH[4]/n0[4])
    hh[2,2]<- s0[1]/s0[2]+sigmaH[3]/n0[3]+sigmaH[7]/n0[7]+sigmaH[8]/n0[8]
    hh[2,3]<- 2.0*sigmaH[7]/n0[7]+sigmaH[8]/n0[8]
    hh[2,4]<- 0
    hh[2,5]<- -(s0[1]/s0[2]+sigmaH[3]/n0[3])
    hh[2,6]<- -(s0[1]/s0[2]+2.0*sigmaH[7]/n0[7])
    hh[3,3]<- 4.0*sigmaH[7]/n0[7]+sigmaH[8]/n0[8]+s0[3]/s0[4]
    hh[3,4]<- hh[3,5]<- 0
    hh[3,6]<- -(4.0*sigmaH[7]/n0[7]+s0[3]/s0[4])
    hh[4,4]<- sigmaH[4]/n0[4]+4.0*sigmaH[5]/n0[5]+sigmaH[6]/n0[6]
    hh[4,5]<- -(sigmaH[4]/n0[4]+2.0*sigmaH[5]/n0[5])
    hh[4,6]<- -2.0*sigmaH[4]/n0[4]
    hh[5,5]<- s0[1]/s0[2]+sigmaH[3]/n0[3]+sigmaH[4]/n0[4]+sigmaH[5]/n0[5]
    hh[5,6]<- s0[1]/s0[2]+2.0*sigmaH[4]/n0[4]
    hh[6,6]<- s0[1]/s0[2]+4.0*sigmaH[4]/n0[4]+4.0*sigmaH[7]/n0[7]+s0[3]/s0[4]
    for(i in 2:6){
      for(j in 1:(i-1)){
        hh[i,j]<- hh[j,i]
      }
    }
    ##############################################################################################################################
    b_line<-matrix(0,6,1)
    b_line[1]<- (sum(dataP1)+m_fam*sumwx[1])/(n_samP1+m_fam*n0[1])-sumwx[2]/n0[2]-sumwx[4]/n0[4]+sumwx[6]/n0[6]
    b_line[2]<- (sum(dataP1)+m_fam*sumwx[1])/(n_samP1+m_fam*n0[1])-sumwx[3]/n0[3]+sumwx[7]/n0[7]-sumwx[8]/n0[8]
    b_line[3]<- 2.0*sumwx[7]/n0[7]-sumwx[8]/n0[8]-(sum(dataP2)+m_fam*sumwx[9])/(n_samP2+m_fam*n0[9])
    b_line[4]<- -sumwx[4]/n0[4]+2.0*sumwx[5]/n0[5]-sumwx[6]/n0[6]
    b_line[5]<- -(sum(dataP1)+m_fam*sumwx[1])/(n_samP1+m_fam*n0[1])+sumwx[3]/n0[3]+sumwx[4]/n0[4]-sumwx[5]/n0[5]
    b_line[6]<- -(sum(dataP1)+m_fam*sumwx[1])/(n_samP1+m_fam*n0[1])+2.0*sumwx[4]/n0[4]+(sum(dataP2)+m_fam*sumwx[9])/(n_samP2+m_fam*n0[9])-2.0*sumwx[7]/n0[7]
    B <- solve(hh,b_line)
    ###############################################################################################################################
    meanH[1]<-(sum(dataP1)+m_fam*sumwx[1]-(B[1]+B[2]-B[5]-B[6])*sigma)/(n_samP1+m_fam*n0[1])
    meanH[2]<-(sumwx[2]+B[1]*sigmaH[2])/n0[2]
    meanH[3]<-(sumwx[3]+(B[2]-B[5])*sigmaH[3])/n0[3]
    meanH[4]<-(sumwx[4]+(B[1]+B[4]-B[5]-2.0*B[6])*sigmaH[4])/n0[4]
    meanH[5]<-(sumwx[5]+(-2.0*B[4]+B[5])*sigmaH[5])/n0[5]
    meanH[6]<-(sumwx[6]-(B[1]-B[4])*sigmaH[6])/n0[6]
    meanH[7]<-(sumwx[7]-(B[2]+2.0*B[3]-2.0*B[6])*sigmaH[7])/n0[7]
    meanH[8]<-(sumwx[8]+(B[2]+B[3])*sigmaH[8])/n0[8]
    meanH[9]<-(sum(dataP2)+m_fam*sumwx[9]+(B[3]-B[6])*sigma)/(n_samP2+m_fam*n0[9])

    mean11<- meanH[1]
    mean12<- meanH[9]
    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanH[i])^2 }

    s0H<-matrix(0,2,1)
    s0H[1]<- ss1+ss2
    s0H[2]<- sum(swx)*m_fam
    s1<- sigmaH[1]-sigma/m_fam    # variance of polygenes.
    if (s1<0.0){s1<- 0.000001}
    abc2<- sigma

    abc3<-1000; n_iter<- 0
    while(abc3>0.0001){
      n_iter<- n_iter+1
      sigma<- (s0H[1]+s0H[2])/(n_samP1+n_samP2+n_samDH)
      abc3<- abs(abc2-sigma)
      abc2<- sigma
      if (n_iter>20) break
    }
    if (sigma<0.1*sigma0){ sigma<- 0.1*sigma0 }
    sigmaH<-matrix((sigma/m_fam),d2,1)

    ########criteria for iterations to stop#######

    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanH,sqrt(sigmaH),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*4
  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1,1,1, 2,2,2,0,0,0,-2,-2,-2, 2,-2,0,2,0,-2,0,2,-2),9,3)
  b_line1 <- meanH
  B1 <- ginv(t(aa)%*%aa)%*%(t(aa)%*%b_line1)
  jj <- sigma_dh - sigmaH[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh
  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanH[i])/sqrt(sigmaH[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("4MG-EEA",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanH),4)," "," "," "," "," "," "," ",round(sigmaH[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," ",round(B1[1],4)," "," ",round(B1[2],4),round(B1[2],4),round(B1[3],4),round(B1[3],4)," "," "," "," "," "," "," "," ",
                       round(jj,4),round(ll*100,4)," "," ",round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}

############################################ 4MG-EEEA model#########################################
G3DHModelFun[[38]] <- function(K1,logL,df11,df21,df31,G3DHtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataP2 <- as.matrix(as.numeric(df21[,1]));dataDH <- as.matrix(as.numeric(df31[,1]))
  n_samP1<-dim(dataP1)[1]; n_samP2<-dim(dataP2)[1];n_samDH<-dim(dataDH)[1]
  mean11<-mean(dataP1);mean12<-mean(dataP2)
  sigmaP1<- as.numeric(var(dataP1)); sigmaP2<- as.numeric(var(dataP2))
  ss1<-(n_samP1-1)*sigmaP1; ss2<-(n_samP2-1)*sigmaP2
  sigma0<-(ss1+ss2)/(n_samP1+n_samP2-2);sigma_dh<-as.numeric(var(dataDH))
  m_esp<- 0.0001 ;m_fam<- as.numeric(G3DHtext2)
  #####################(H-5)Model##########################################
  d2<-8
  mi <- as.matrix(c(0.0625,0.0625,0.1875,0.1875,0.1875,0.1875,0.0625,0.0625))
  meanH<-mean(dataDH)
  sigmaH <- matrix((sigma_dh/(2*1.2222101)),d2,1)
  sigma <- sigma0
  a1 <- sqrt(sigma_dh/n_samDH)
  if(mean11<mean12){a1=-a1}

  meanH <- as.matrix(c((meanH+3*a1),(meanH+2.5*a1),(meanH+2*a1),(meanH+1.5*a1),
                       (meanH-1.5*a1),(meanH-2*a1),(meanH-2.5*a1),(meanH-3*a1)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,d2,n_samDH); swx <- matrix(0,d2,1)
  L0 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanH,sqrt(sigmaH),mi)))
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataDH,meanH[i],sqrt(sigmaH[i]))/dmixnorm(dataDH,meanH,sqrt(sigmaH),mi) }
    mix_pi <- as.matrix(rowSums(WW)/n_samDH)
    sumwx <- WW%*%dataDH
    n0 <- n_samDH*mix_pi
    n0[n0<0.000001] <- 0.000001

    s0<-matrix(0,4,1)
    s0[1]<- sigma;s0[2]<- n_samP1+m_fam*n0[1]
    s0[3]<- sigma;s0[4]<- n_samP2+m_fam*n0[8]

    hh<-matrix(0,5,5)
    hh[1,1]<- 4.0*s0[1]/s0[2]+9.0*sigmaH[4]/n0[4]+s0[3]/s0[4]
    hh[1,2]<- -(2.0*s0[1]/s0[2]+3.0*sigmaH[4]/n0[4])
    hh[1,3]<- 0
    hh[1,4]<- -2.0*s0[1]/s0[2]
    hh[1,5]<- 0
    hh[2,2]<- s0[1]/s0[2]+sigmaH[4]/n0[4]+sigmaH[2]/n0[2]+sigmaH[3]/n0[3]
    hh[2,3]<- (sigmaH[2]/n0[2]+2.0*sigmaH[3]/n0[3])
    hh[2,4]<- s0[1]/s0[2]+sigmaH[2]/n0[2]
    hh[2,5]<- 2.0*sigmaH[2]/n0[2]+3.0*sigmaH[3]/n0[3]
    hh[3,3]<- sigmaH[2]/n0[2]+4.0*sigmaH[3]/n0[3]+sigmaH[5]/n0[5]
    hh[3,4]<- sigmaH[2]/n0[2]-sigmaH[5]/n0[5]
    hh[3,5]<- 2.0*sigmaH[2]/n0[2]+6.0*sigmaH[3]/n0[3]
    hh[4,4]<- s0[1]/s0[2]+sigmaH[2]/n0[2]+sigmaH[5]/n0[5]+sigmaH[6]/n0[6]
    hh[4,5]<- 2.0*sigmaH[2]/n0[2]
    hh[5,5]<- 4.0*sigmaH[2]/n0[2]+9.0*sigmaH[3]/n0[3]+sigmaH[7]/n0[7]
    for(i in 2:5){
      for(j in 1:(i-1)){
        hh[i,j]<- hh[j,i]
      }
    }
    ##############################################################################################################################
    b_line<-matrix(0,5,1)
    b_line[1]<- -2.0*(sum(dataP1)+m_fam*sumwx[1])/(n_samP1+m_fam*n0[1])+3.0*sumwx[4]/n0[4]-(sum(dataP2)+m_fam*sumwx[8])/(n_samP2+m_fam*n0[8])
    b_line[2]<- (sum(dataP1)+m_fam*sumwx[1])/(n_samP1+m_fam*n0[1])-sumwx[2]/n0[2]+sumwx[3]/n0[3]-sumwx[4]/n0[4]
    b_line[3]<- -sumwx[2]/n0[2]+2.0*sumwx[3]/n0[3]-sumwx[5]/n0[5]
    b_line[4]<- (sum(dataP1)+m_fam*sumwx[1])/(n_samP1+m_fam*n0[1])-sumwx[2]/n0[2]+sumwx[5]/n0[5]-sumwx[6]/n0[6]
    b_line[5]<- -2.0*sumwx[2]/n0[2]+3.0*sumwx[3]/n0[3]-sumwx[7]/n0[7]

    B <- solve(hh,b_line)
    ###############################################################################################################################
    meanH[1]<-(sum(dataP1)+m_fam*sumwx[1]+(2.0*B[1]-B[2]-B[4])*sigma)/(n_samP1+m_fam*n0[1])
    meanH[2]<-(sumwx[2]+(B[2]+B[3]+B[4]+2.0*B[5])*sigmaH[2])/n0[2]
    meanH[3]<-(sumwx[3]-(B[2]+2.0*B[3]+3.0*B[5])*sigmaH[3])/n0[3]
    meanH[4]<-(sumwx[4]+(-3.0*B[1]+B[2])*sigmaH[4])/n0[4]
    meanH[5]<-(sumwx[5]+(B[3]-B[4])*sigmaH[5])/n0[5]
    meanH[6]<-(sumwx[6]+B[4]*sigmaH[6])/n0[6]
    meanH[7]<-(sumwx[7]+B[5]*sigmaH[7])/n0[7]
    meanH[8]<-(sum(dataP2)+m_fam*sumwx[8]+B[1]*sigma)/(n_samP2+m_fam*n0[8])

    mean11<- meanH[1]
    mean12<- meanH[8]
    ##########obtain variance##########
    ss1<- sum((dataP1-mean11)^2)
    ss2<- sum((dataP2-mean12)^2)
    for(i in 1:d2) {  swx[i] <- WW[i,]%*%(dataDH-meanH[i])^2 }

    s0H<-matrix(0,2,1)
    s0H[1]<- ss1+ss2
    s0H[2]<- sum(swx)*m_fam
    s1<- sigmaH[1]-sigma/m_fam    # variance of polygenes.
    if (s1<0.0){s1<- 0.000001}
    abc2<- sigma

    abc3<-1000; n_iter<- 0
    while(abc3>0.0001){
      n_iter<- n_iter+1
      sigma<- (s0H[1]+s0H[2])/(n_samP1+n_samP2+n_samDH)
      abc3<- abs(abc2-sigma)
      abc2<- sigma
      if (n_iter>20) break
    }
    if (sigma<0.1*sigma0){ sigma<- 0.1*sigma0 }
    sigmaH<-matrix((sigma/m_fam),d2,1)

    ########criteria for iterations to stop#######

    L1 <- logL(n_samP1,1,1,mean11,sigma,dataP1)+logL(n_samP2,1,1,mean12,sigma,dataP2)+sum(log(dmixnorm(dataDH,meanH,sqrt(sigmaH),mix_pi)))

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*4
  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1,1, 3,3,1,1,-1,-1,-3,-3, 1,-1,-1,1,-1,1,-1,1),8,3)
  b_line1 <- meanH
  B1 <- ginv(t(aa)%*%aa)%*%(t(aa)%*%b_line1)
  jj <- sigma_dh - sigmaH[1]
  if(jj < 0) {jj <- 0}
  ll <- jj/sigma_dh
  #########hypothesis testing#########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - mean11)/sqrt(as.vector(sigma))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - mean12)/sqrt(as.vector(sigma))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  DH  ###################
  dataDH<-sort(dataDH);bmw <- matrix(0,n_samDH,1); bmwsl <- matrix(0,n_samDH,d2)
  for(i in 1:d2){
    gg <- (dataDH - meanH[i])/sqrt(sigmaH[i])
    bmw[which(gg>=0)] <- pnorm(gg[gg>=0])
    bmw[which(gg<0)] <- 1 - pnorm(abs(gg[gg<0]))
    bmwsl[,i] <- bmw*mix_pi[i]
  }
  P2 <- rowSums(bmwsl)
  nn<-dim(as.matrix(unique(P2)))[1]
  if(nn<n_samDH){P2<-P2+runif(n_samDH)/1e4}

  dd <- as.matrix(c(sum(P2),sum(P2^2),sum((P2-0.5)^2)))
  WW2 <- 1/(12*n_samDH) + sum((P2 - (as.matrix(c(1:n_samDH)) - 0.5)/n_samDH)^2)
  u <- as.matrix(c(12*n_samDH*((dd[1]/n_samDH-0.5)^2),((45*n_samDH)/4)*((dd[2]/n_samDH-1/3)^2),180*n_samDH*((dd[3]/n_samDH-1/12)^2)))
  D <- as.numeric(ks.test(P2,"punif")[[1]][1])
  tt <- as.matrix(c((1 - pchisq(u[1],1)),(1 - pchisq(u[2],1)),(1 - pchisq(u[3],1)),K1(WW2),(1-pkolm(D,n_samDH))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt[which(tt>=10e-4)]<-round(tt[which(tt>=10e-4)],4);tt[which(tt<10e-4)]<-format(tt[which(tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("4MG-EEEA",round(abc,4),round(AIC,4),round(mean11,4),round(mean12,4),round(sigma,4),round(t(meanH),4)," "," "," "," "," "," "," "," ",round(sigmaH[1],4),round(t(mix_pi),4),
                       " "," "," "," "," "," "," "," ",round(B1[1],4)," "," ",round(B1[2],4),round(B1[2],4),round(B1[2],4),round(B1[3],4)," "," "," "," "," "," "," "," ",
                       round(jj,4),round(ll*100,4)," "," ",round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u[1],4),tt[1],round(u[2],4),
                       tt[2],round(u[3],4),tt[3],round(WW2,4),tt[4],round(D,4),tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}



K1G3DH <- function(x){
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

logLG3DH <- function(nm,nng,mi,mn,s,d1) { sum2 <- sum(log(dmixnorm(d1,mn,sqrt(s),mi)));return (sum2) }


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
allresult=foreach(i=1:38,.combine = 'rbind')%dopar%{
  requireNamespace("KScorrect")
  requireNamespace("kolmim")
  requireNamespace("MASS")
  G3DHModelFun[[i]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2)[[1]]
}
stopCluster(cl)
mi<-NULL
}else{


  allresultq<-switch(model,"0MG"=G3DHModelFun[[1]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),"1MG-A"=G3DHModelFun[[2]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),"2MG-AI"=G3DHModelFun[[3]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),
                     "2MG-A"=G3DHModelFun[[4]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),"2MG-EA"=G3DHModelFun[[5]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),"2MG-ED"=G3DHModelFun[[6]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),
                     "2MG-ER"=G3DHModelFun[[7]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),"2MG-AE"=G3DHModelFun[[8]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),"2MG-CE"=G3DHModelFun[[9]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),
                     "2MG-DE"=G3DHModelFun[[10]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),"PG-AI"=G3DHModelFun[[12]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),
                     "PG-A"=G3DHModelFun[[13]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),"MX1-A-AI"=G3DHModelFun[[14]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),"MX1-A-A"=G3DHModelFun[[15]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),
                     "MX2-AI-AI"=G3DHModelFun[[16]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),"MX2-AI-A"=G3DHModelFun[[17]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),"MX2-A-A"=G3DHModelFun[[18]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),
                     "MX2-EA-A"=G3DHModelFun[[19]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),"MX2-ED-A"=G3DHModelFun[[20]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),"MX2-ER-A"=G3DHModelFun[[21]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),
                     "MX2-AE-A"=G3DHModelFun[[22]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),"MX2-CE-A"=G3DHModelFun[[23]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),"MX2-DE-A"=G3DHModelFun[[24]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),
                     "MX2-IE-A"=G3DHModelFun[[25]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),"3MG-AI"=G3DHModelFun[[26]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),"3MG-A"=G3DHModelFun[[27]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),
                     "3MG-CEA"=G3DHModelFun[[28]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),"3MG-PEA"=G3DHModelFun[[29]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),"MX3-AI-AI"=G3DHModelFun[[30]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),
                     "MX3-AI-A"=G3DHModelFun[[31]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),"MX3-A-A"=G3DHModelFun[[32]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),"MX3-CEA-A"=G3DHModelFun[[33]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),
                     "MX3-PEA-A"=G3DHModelFun[[34]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),"4MG-AI"=G3DHModelFun[[35]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),"4MG-CEA"=G3DHModelFun[[36]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),
                     "4MG-EEA"=G3DHModelFun[[37]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2),"4MG-EEEA"=G3DHModelFun[[38]](K1G3DH,logLG3DH,df11,df21,df31,G3DHtext2))

  allresult<-allresultq[[1]]
  if(model=="0MG"||model=="PG-A"||model=="PG-AI"){
    mi<-NULL
  }else{
    mi<-allresultq[[2]]
  }
}
colnames(allresult) <- G3DHcolname
out<-list(allresult,mi)
return(out)
}
