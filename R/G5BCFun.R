G5BCFun<-function(df,model){


data<-sapply(df,as.character)

dP1<-data[-1,which(data[1,]=="P1")];P1<-as.numeric(dP1[which(is.na(as.numeric(dP1))==FALSE)]);df11<-as.data.frame(P1)
dF1<-data[-1,which(data[1,]=="F1")];F1<-as.numeric(dF1[which(is.na(as.numeric(dF1))==FALSE)]);df21<-as.data.frame(F1)
dP2<-data[-1,which(data[1,]=="P2")];P2<-as.numeric(dP2[which(is.na(as.numeric(dP2))==FALSE)]);df31<-as.data.frame(P2)
dB1<-data[-1,which(data[1,]=="B1")];B1<-as.numeric(dB1[which(is.na(as.numeric(dB1))==FALSE)]);df41<-as.data.frame(B1)
dB2<-data[-1,which(data[1,]=="B2")];B2<-as.numeric(dB2[which(is.na(as.numeric(dB2))==FALSE)]);df51<-as.data.frame(B2)

G5BCcolname <- c("Model","Log_Max_likelihood_Value","AIC","mean[P1]","mean[F1]","mean[P2]","Var(P1 & P2 & F1)","B1-mean[1]","B1-mean[2]","B1-mean[3]","B1-mean[4]",
                 "B1-Var(Residual+Polygene)","B1-Proportion[1]","B1-Proportion[2]","B1-Proportion[3]","B1-Proportion[4]","B2-mean[1]","B2-mean[2]","B2-mean[3]","B2-mean[4]",
                 "B2-Var(Residual+Polygene)","B2-Proportion[1]","B2-Proportion[2]","B2-Proportion[3]","B2-Proportion[4]","m","da","db","ha","hb",
                 "[d]","[h]","B1-MajorGene Var","B1-Heritability(MajorGene)(%)","B1-Polygenes Var","B1-Heritability(Polygenes)(%)","B2-MajorGene Var","B2-Heritability(MajorGene)(%)",
                 "B2-Polygenes Var","B2-Heritability(Polygenes)(%)","U1 square-P1","P(U1 square-P1)","U2 square-P1","P(U2 square-P1)","U3 square-P1","P(U3 square-P1)","nW square-P1","P(nW square-P1)","Dn-P1","P(Dn-P1)",
                 "U1 square-F1","P(U1 square-F1)","U2 square-F1","P(U2 square-F1)","U3 square-F1","P(U3 square-F1)","nW square-F1","P(nW square-F1)","Dn-F1","P(Dn-F1)","U1 square-P2","P(U1 square-P2)","U2 square-P2","P(U2 square-P2)","U3 square-P2","P(U3 square-P2)","nW square-P2","P(nW square-P2)","Dn-P2","P(Dn-P2)",
                 "U1 square-B1","P(U1 square-B1)","U2 square-B1","P(U2 square-B1)","U3 square-B1","P(U3 square-B1)","nW square-B1","P(nW square-B1)","Dn-B1","P(Dn-B1)","U1 square-B2","P(U1 square-B2)","U2 square-B2","P(U2 square-B2)","U3 square-B2","P(U3 square-B2)","nW square-B2","P(nW square-B2)","Dn-B2","P(Dn-B2)")

G5BCModelFun<-list(NA)
################### define each model function ##################
########################### (A1)#####################################
G5BCModelFun[[1]] <- function(K1,logL,df11,df21,df31,df41,df51){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]))
  n_samP1<-dim(dataP1)[1];n_samP2<-dim(dataP2)[1];n_samF1<-dim(dataF1)[1];n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2);sumx4<-sum(dataB1);sumx5<-sum(dataB2)
  ss1<- (n_samP1-1)*as.numeric(var(dataP1)); meanP1<- mean(dataP1)
  ss2<- (n_samF1-1)*as.numeric(var(dataF1)); meanF1<- mean(dataF1)
  ss3<- (n_samP2-1)*as.numeric(var(dataP2)); meanP2<- mean(dataP2)
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  sigmaP1<- sigma0;sigmaF1<- sigma0;sigmaP2<- sigma0
  mean1<- mean(dataB1); sigma1<- as.numeric(var(dataB1))
  mean2<- mean(dataB2); sigma2<- as.numeric(var(dataB2))

  m_esp<-0.0001
  ####################### 1MG-AD Model ############ (A1)#####################
  d21<-2; d22<-2
  mi_1<- matrix(0.5,d21,1); a1<- sqrt(sigma1/(n_samB1-1))
  sigma1A<- matrix(sigma1/3,d21,1); mean1A<- as.matrix(c(mean1+a1,mean1-a1))

  mi_2<- matrix(0.5,d22,1); a2<- sqrt(sigma2/(n_samB2-1))
  sigma2A<- matrix(sigma2/3,d22,1); mean2A<- as.matrix(c(mean2+a2,mean2-a2))

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1A,sigma1A,dataB1)+logL(n_samB2,d22,mi_2,mean2A,sigma2A,dataB2)
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1A[i],sqrt(sigma1A[i]))/dmixnorm(dataB1,mean1A,sqrt(sigma1A),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2A[i],sqrt(sigma2A[i]))/dmixnorm(dataB2,mean2A,sqrt(sigma2A),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    ############ CM1-step for means ##############
    n01 <- n_samB1*mix_pi_1; n01[n01<0.000001] <- 0.000001
    n02 <- n_samB2*mix_pi_2; n02[n02<0.000001] <- 0.000001

    meanP1<- (sumx1)/(n_samP1)
    meanF1<- (sumx2)/(n_samF1)
    meanP2<- (sumx3)/(n_samP2)
    mean1A<- as.matrix(c(sumwx_B1[1]/n01[1],(sumwx_B1[2])/n01[2]));
    mean2A<- as.matrix(c((sumwx_B2[1])/(n02[1]),sumwx_B2[2]/n02[2]))

    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2); ss2<- sum((dataF1-meanF1)^2); ss3<- sum((dataP2-meanP2)^2)
    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1A[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2A[i])^2 }

    sigma1A[1]<- (sum(swx_B1)+sum(swx_B2)+ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2+n_samB1+n_samB2)
    sigma1A<- matrix(sigma1A[1],d21,1);sigma2A<- matrix(sigma1A[1],d22,1)
    sigmaP1<- sigma1A[1]; sigmaF1<- sigma1A[1]; sigmaP2<- sigma1A[1]
    if(sum(sigma1A < 1e-30)>=1){break}
    if(sum(sigma2A < 1e-30)>=1){break}
    ####################### the stop criterion for iteration #####################################
    L1<- logL(n_samP1,1,1,meanP1,sigmaP1,dataP1)+logL(n_samF1,1,1,meanF1,sigmaF1,dataF1)+logL(n_samP2,1,1,meanP2,sigmaP2,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1A,sigma1A,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2A,sigma2A,dataB2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*8

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1, 1,0,-1, 0,1,0),3,3)
  b_line1 <- matrix(c(mean1A,mean2A[2]))
  B1 <- solve(aa,b_line1)

  jj_1 <- sigma1 - sigma1A[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigma1

  jj_2 <- sigma2 - sigma2A[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigma2

  ######### hypothesis testing #########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - meanP1)/sqrt(as.vector(sigmaP1))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[2])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),D_P1))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])

  ###############  F1  ###################
  dataF1<-sort(dataF1);bmw_F1 <- matrix(0,n_samF1,1); bmwsl_F1 <- matrix(0,n_samF1,1)

  gg_F1 <- (dataF1 - meanF1)/sqrt(as.vector(sigmaF1))
  bmw_F1[which(gg_F1>=0)] <- pnorm(gg_F1[gg_F1>=0])
  bmw_F1[which(gg_F1<0)] <- 1 - pnorm(abs(gg_F1[gg_F1<0]))
  bmwsl_F1[,1] <- bmw_F1

  P2_F1 <- rowSums(bmwsl_F1)
  nn<-dim(as.matrix(unique(P2_F1)))[1]
  if(nn<n_samF1){P2_F1<-P2_F1+runif(n_samF1)/1e4}

  dd_F1 <- as.matrix(c(sum(P2_F1),sum(P2_F1^2),sum((P2_F1-0.5)^2)))
  WW2_F1 <- 1/(12*n_samF1) + sum((P2_F1 - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  u_F1 <- as.matrix(c(12*n_samF1*((dd_F1[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((dd_F1[2]/n_samF1-1/3)^2),180*n_samF1*((dd_F1[3]/n_samF1-1/12)^2)))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[2])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),D_F1))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - meanP2)/sqrt(as.vector(sigmaP2))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[2])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),D_P2))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])

  ###############  B1  ###################
  dataB1<-sort(dataB1);bmw_B1 <- matrix(0,n_samB1,1); bmwsl_B1 <- matrix(0,n_samB1,d21)
  for(i in 1:d21){
    gg_B1 <- (dataB1 - mean1A[i])/sqrt(sigma1A[i])
    bmw_B1[which(gg_B1>=0)] <- pnorm(gg_B1[gg_B1>=0])
    bmw_B1[which(gg_B1<0)] <- 1 - pnorm(abs(gg_B1[gg_B1<0]))
    bmwsl_B1[,i] <- bmw_B1*mix_pi_1[i]
  }
  P2_B1 <- rowSums(bmwsl_B1)
  nn<-dim(as.matrix(unique(P2_B1)))[1]
  if(nn<n_samB1){P2_B1<-P2_B1+runif(n_samB1)/1e4}

  dd_B1 <- as.matrix(c(sum(P2_B1),sum(P2_B1^2),sum((P2_B1-0.5)^2)))
  WW2_B1 <- 1/(12*n_samB1) + sum((P2_B1 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  u_B1 <- as.matrix(c(12*n_samB1*((dd_B1[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((dd_B1[2]/n_samB1-1/3)^2),180*n_samB1*((dd_B1[3]/n_samB1-1/12)^2)))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[2])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),D_B1))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])

  ###############  B2  ###################
  dataB2<-sort(dataB2);bmw_B2 <- matrix(0,n_samB2,1); bmwsl_B2 <- matrix(0,n_samB2,d22)
  for(i in 1:d22){
    gg_B2 <- (dataB2 - mean2A[i])/sqrt(sigma2A[i])
    bmw_B2[which(gg_B2>=0)] <- pnorm(gg_B2[gg_B2>=0])
    bmw_B2[which(gg_B2<0)] <- 1 - pnorm(abs(gg_B2[gg_B2<0]))
    bmwsl_B2[,i] <- bmw_B2*mix_pi_2[i]
  }
  P2_B2 <- rowSums(bmwsl_B2)
  nn<-dim(as.matrix(unique(P2_B2)))[1]
  if(nn<n_samB2){P2_B2<-P2_B2+runif(n_samB2)/1e4}

  dd_B2 <- as.matrix(c(sum(P2_B2),sum(P2_B2^2),sum((P2_B2-0.5)^2)))
  WW2_B2 <- 1/(12*n_samB2) + sum((P2_B2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  u_B2 <- as.matrix(c(12*n_samB2*((dd_B2[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((dd_B2[2]/n_samB2-1/3)^2),180*n_samB2*((dd_B2[3]/n_samB2-1/12)^2)))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[2])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),D_B2))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("1MG-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigmaP1,4),round(t(mean1A),4)," "," ",round(sigma1A[1],4),round(t(mix_pi_1),4)," "," ",
                       round(t(mean2A),4)," "," ",round(sigma2A[1],4),round(t(mix_pi_2),4)," "," ",
                       round(B1[1],4),round(B1[2],4)," ",round(B1[3],4)," "," "," ",round(jj_1,4),round(ll_1*100,4)," "," ",round(jj_2,4),round(ll_2*100,4)," "," ",
                       round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],
                       round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],
                       round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],
                       round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],
                       round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
}
########################### (A2)#####################################
G5BCModelFun[[2]] <- function(K1,logL,df11,df21,df31,df41,df51){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]))
  n_samP1<-dim(dataP1)[1];n_samP2<-dim(dataP2)[1];n_samF1<-dim(dataF1)[1];n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2);sumx4<-sum(dataB1);sumx5<-sum(dataB2)
  ss1<- (n_samP1-1)*as.numeric(var(dataP1)); meanP1<- mean(dataP1)
  ss2<- (n_samF1-1)*as.numeric(var(dataF1)); meanF1<- mean(dataF1)
  ss3<- (n_samP2-1)*as.numeric(var(dataP2)); meanP2<- mean(dataP2)
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  sigmaP1<- sigma0;sigmaF1<- sigma0;sigmaP2<- sigma0
  mean1<- mean(dataB1); sigma1<- as.numeric(var(dataB1))
  mean2<- mean(dataB2); sigma2<- as.numeric(var(dataB2))

  m_esp<-0.0001
  ####################### 1MG-A Model ############ (A2)#####################
  d21<-2; d22<-2
  mi_1<- matrix(0.5,d21,1); a1<- sqrt(sigma1/(n_samB1-1))
  sigma1A<- matrix(sigma1/3,d21,1); mean1A<- as.matrix(c(mean1+a1,mean1-a1))

  mi_2<- matrix(0.5,d22,1); a2<- sqrt(sigma2/(n_samB2-1))
  sigma2A<- matrix(sigma2/3,d22,1); mean2A<- as.matrix(c(mean2+a2,mean2-a2))

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1A,sigma1A,dataB1)+logL(n_samB2,d22,mi_2,mean2A,sigma2A,dataB2)
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1A[i],sqrt(sigma1A[i]))/dmixnorm(dataB1,mean1A,sqrt(sigma1A),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2A[i],sqrt(sigma2A[i]))/dmixnorm(dataB2,mean2A,sqrt(sigma2A),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    ############ CM1-step for means ##############
    n01 <- n_samB1*mix_pi_1; n01[n01<0.000001] <- 0.000001
    n02 <- n_samB2*mix_pi_2; n02[n02<0.000001] <- 0.000001

    s0<- matrix(0,6,1)
    s0[1]<- sumx1+sumwx_B1[1]; s0[2]<- n_samP1+n01[1]
    s0[3]<- sumx2+sumwx_B1[2]+sumwx_B2[1]; s0[4]<- n_samF1+n01[2]+n02[1]
    s0[5]<- sumx3+sumwx_B2[2]; s0[6]<- n_samP2+n02[2]

    aa1<- s0[1]/s0[2]-2.0*s0[3]/s0[4]+s0[5]/s0[6]
    aa2<- sigmaP1/s0[2]+4.0*sigmaF1/s0[4]+sigmaP2/s0[6]
    aa3<- aa1/aa2

    meanP1<- (s0[1]-sigmaP1*aa3)/s0[2]
    meanF1<- (s0[3]+sigmaF1*aa3*2.0)/s0[4]
    meanP2<- (s0[5]-sigmaP2*aa3)/s0[6]
    mean1A<- as.matrix(c(meanP1,meanF1));mean2A<- as.matrix(c(meanF1,meanP2))
    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2); ss2<- sum((dataF1-meanF1)^2); ss3<- sum((dataP2-meanP2)^2)
    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1A[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2A[i])^2 }

    sigma1A[1]<- (sum(swx_B1)+sum(swx_B2)+ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2+n_samB1+n_samB2)
    sigma1A<- matrix(sigma1A[1],d21,1);sigma2A<- matrix(sigma1A[1],d22,1)
    sigmaP1<- sigma1A[1]; sigmaF1<- sigma1A[1]; sigmaP2<- sigma1A[1]
    if(sum(sigma1A < 1e-30)>=1){break}
    if(sum(sigma2A < 1e-30)>=1){break}
    ####################### the stop criterion for iteration #####################################
    L1<- logL(n_samP1,1,1,meanP1,sigmaP1,dataP1)+logL(n_samF1,1,1,meanF1,sigmaF1,dataF1)+logL(n_samP2,1,1,meanP2,sigmaP2,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1A,sigma1A,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2A,sigma2A,dataB2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*7

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1, 1,0,-1),3,2)
  b_line1 <- matrix(c(mean1A,mean2A[2]))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigma1 - sigma1A[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigma1

  jj_2 <- sigma2 - sigma2A[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigma2

  ######### hypothesis testing #########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - meanP1)/sqrt(as.vector(sigmaP1))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[2])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),D_P1))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])

  ###############  F1  ###################
  dataF1<-sort(dataF1);bmw_F1 <- matrix(0,n_samF1,1); bmwsl_F1 <- matrix(0,n_samF1,1)

  gg_F1 <- (dataF1 - meanF1)/sqrt(as.vector(sigmaF1))
  bmw_F1[which(gg_F1>=0)] <- pnorm(gg_F1[gg_F1>=0])
  bmw_F1[which(gg_F1<0)] <- 1 - pnorm(abs(gg_F1[gg_F1<0]))
  bmwsl_F1[,1] <- bmw_F1

  P2_F1 <- rowSums(bmwsl_F1)
  nn<-dim(as.matrix(unique(P2_F1)))[1]
  if(nn<n_samF1){P2_F1<-P2_F1+runif(n_samF1)/1e4}

  dd_F1 <- as.matrix(c(sum(P2_F1),sum(P2_F1^2),sum((P2_F1-0.5)^2)))
  WW2_F1 <- 1/(12*n_samF1) + sum((P2_F1 - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  u_F1 <- as.matrix(c(12*n_samF1*((dd_F1[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((dd_F1[2]/n_samF1-1/3)^2),180*n_samF1*((dd_F1[3]/n_samF1-1/12)^2)))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[2])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),D_F1))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - meanP2)/sqrt(as.vector(sigmaP2))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[2])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),D_P2))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])

  ###############  B1  ###################
  dataB1<-sort(dataB1);bmw_B1 <- matrix(0,n_samB1,1); bmwsl_B1 <- matrix(0,n_samB1,d21)
  for(i in 1:d21){
    gg_B1 <- (dataB1 - mean1A[i])/sqrt(sigma1A[i])
    bmw_B1[which(gg_B1>=0)] <- pnorm(gg_B1[gg_B1>=0])
    bmw_B1[which(gg_B1<0)] <- 1 - pnorm(abs(gg_B1[gg_B1<0]))
    bmwsl_B1[,i] <- bmw_B1*mix_pi_1[i]
  }
  P2_B1 <- rowSums(bmwsl_B1)
  nn<-dim(as.matrix(unique(P2_B1)))[1]
  if(nn<n_samB1){P2_B1<-P2_B1+runif(n_samB1)/1e4}

  dd_B1 <- as.matrix(c(sum(P2_B1),sum(P2_B1^2),sum((P2_B1-0.5)^2)))
  WW2_B1 <- 1/(12*n_samB1) + sum((P2_B1 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  u_B1 <- as.matrix(c(12*n_samB1*((dd_B1[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((dd_B1[2]/n_samB1-1/3)^2),180*n_samB1*((dd_B1[3]/n_samB1-1/12)^2)))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[2])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),D_B1))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])

  ###############  B2  ###################
  dataB2<-sort(dataB2);bmw_B2 <- matrix(0,n_samB2,1); bmwsl_B2 <- matrix(0,n_samB2,d22)
  for(i in 1:d22){
    gg_B2 <- (dataB2 - mean2A[i])/sqrt(sigma2A[i])
    bmw_B2[which(gg_B2>=0)] <- pnorm(gg_B2[gg_B2>=0])
    bmw_B2[which(gg_B2<0)] <- 1 - pnorm(abs(gg_B2[gg_B2<0]))
    bmwsl_B2[,i] <- bmw_B2*mix_pi_2[i]
  }
  P2_B2 <- rowSums(bmwsl_B2)
  nn<-dim(as.matrix(unique(P2_B2)))[1]
  if(nn<n_samB2){P2_B2<-P2_B2+runif(n_samB2)/1e4}

  dd_B2 <- as.matrix(c(sum(P2_B2),sum(P2_B2^2),sum((P2_B2-0.5)^2)))
  WW2_B2 <- 1/(12*n_samB2) + sum((P2_B2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  u_B2 <- as.matrix(c(12*n_samB2*((dd_B2[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((dd_B2[2]/n_samB2-1/3)^2),180*n_samB2*((dd_B2[3]/n_samB2-1/12)^2)))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[2])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),D_B2))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("1MG-A",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigmaP1,4),round(t(mean1A),4)," "," ",round(sigma1A[1],4),round(t(mix_pi_1),4)," "," ",
                       round(t(mean2A),4)," "," ",round(sigma2A[1],4),round(t(mix_pi_2),4)," "," ",
                       round(B1[1],4),round(B1[2],4)," "," "," "," "," ",round(jj_1,4),round(ll_1*100,4)," "," ",round(jj_2,4),round(ll_2*100,4)," "," ",
                       round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],
                       round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],
                       round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],
                       round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],
                       round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
}
########################### (A3)#####################################
G5BCModelFun[[3]] <- function(K1,logL,df11,df21,df31,df41,df51){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]))
  n_samP1<-dim(dataP1)[1];n_samP2<-dim(dataP2)[1];n_samF1<-dim(dataF1)[1];n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2);sumx4<-sum(dataB1);sumx5<-sum(dataB2)
  ss1<- (n_samP1-1)*as.numeric(var(dataP1));meanP1<- mean(dataP1)
  ss2<- (n_samF1-1)*as.numeric(var(dataF1));meanF1<- mean(dataF1)
  ss3<- (n_samP2-1)*as.numeric(var(dataP2));meanP2<- mean(dataP2)
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  sigmaP1<- sigma0;sigmaF1<- sigma0;sigmaP2<- sigma0
  mean1<- mean(dataB1); sigma1<- as.numeric(var(dataB1))
  mean2<- mean(dataB2); sigma2<- as.numeric(var(dataB2))

  m_esp<-0.0001
  ####################### 1MG-EAD Model ############ (A3)#####################
  d21<-1; d22<-2
  mi_1<- matrix(1,d21,1)
  sigma1A<- as.matrix(sigma1/3); mean1A<- as.matrix(mean1)

  mi_2<- matrix(0.5,d22,1); a2<- sqrt(sigma2/(n_samB2-1))
  sigma2A<- matrix(sigma2/3,d22,1); mean2A<- as.matrix(c(mean2+a2,mean2-a2))

  L0 <- logL(n_samP1,1,1,meanP1,sigmaP1,dataP1)+logL(n_samF1,1,1,meanF1,sigmaF1,dataF1)+logL(n_samP2,1,1,meanP2,sigmaP2,dataP2)+logL(n_samB1,d21,mi_1,mean1A,sigma1A,dataB1)+logL(n_samB2,d22,mi_2,mean2A,sigma2A,dataB2)
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1A[i],sqrt(sigma1A[i]))/dmixnorm(dataB1,mean1A,sqrt(sigma1A),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2A[i],sqrt(sigma2A[i]))/dmixnorm(dataB2,mean2A,sqrt(sigma2A),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    ############ CM1-step for means ##############
    n01 <- n_samB1*mix_pi_1; n01[n01<0.000001] <- 0.000001
    n02 <- n_samB2*mix_pi_2; n02[n02<0.000001] <- 0.000001

    meanP1<- (sumx1)/(n_samP1);
    meanF1<- sumx2/n_samF1
    meanP2<- (sumx3)/(n_samP2)
    mean1A<- sumwx_B1[1]/n01 ; mean2A<- as.matrix(c(sumwx_B2[1]/n02[1],sumwx_B2[2]/n02[2]))

    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2); ss2<- sum((dataF1-meanF1)^2); ss3<- sum((dataP2-meanP2)^2)
    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1A[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2A[i])^2 }

    sigma1A[1]<- (sum(swx_B1)+sum(swx_B2)+ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2+n_samB1+n_samB2)
    sigma1A<- matrix(sum(swx_B1)/n_samB1,d21,1);sigma2A<- matrix(sum(swx_B2)/n_samB2,d22,1)
    sigmaP1<- ss1/n_samP1; sigmaF1<- ss2/n_samF1; sigmaP2<- ss3/n_samP2
    if(sum(sigma1A < 1e-30)>=1){break}
    if(sum(sigma2A < 1e-30)>=1){break}
    ####################### the stop criterion for iteration #####################################
    L1<- logL(n_samP1,1,1,meanP1,sigmaP1,dataP1)+logL(n_samF1,1,1,meanF1,sigmaF1,dataF1)+logL(n_samP2,1,1,meanP2,sigmaP2,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1A,sigma1A,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2A,sigma2A,dataB2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*7

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1, 1,-1),2,2)
  b_line1 <- matrix(c(mean1A,mean2A[2]))
  B1 <- solve(aa,b_line1)

  jj_1 <- sigma1 - sigma1A[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigma1

  jj_2 <- sigma2 - sigma2A[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigma2

  ######### hypothesis testing #########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - meanP1)/sqrt(as.vector(sigmaP1))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[2])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),D_P1))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])

  ###############  F1  ###################
  dataF1<-sort(dataF1);bmw_F1 <- matrix(0,n_samF1,1); bmwsl_F1 <- matrix(0,n_samF1,1)

  gg_F1 <- (dataF1 - meanF1)/sqrt(as.vector(sigmaF1))
  bmw_F1[which(gg_F1>=0)] <- pnorm(gg_F1[gg_F1>=0])
  bmw_F1[which(gg_F1<0)] <- 1 - pnorm(abs(gg_F1[gg_F1<0]))
  bmwsl_F1[,1] <- bmw_F1

  P2_F1 <- rowSums(bmwsl_F1)
  nn<-dim(as.matrix(unique(P2_F1)))[1]
  if(nn<n_samF1){P2_F1<-P2_F1+runif(n_samF1)/1e4}

  dd_F1 <- as.matrix(c(sum(P2_F1),sum(P2_F1^2),sum((P2_F1-0.5)^2)))
  WW2_F1 <- 1/(12*n_samF1) + sum((P2_F1 - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  u_F1 <- as.matrix(c(12*n_samF1*((dd_F1[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((dd_F1[2]/n_samF1-1/3)^2),180*n_samF1*((dd_F1[3]/n_samF1-1/12)^2)))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[2])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),D_F1))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - meanP2)/sqrt(as.vector(sigmaP2))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[2])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),D_P2))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])

  ###############  B1  ###################
  dataB1<-sort(dataB1);bmw_B1 <- matrix(0,n_samB1,1); bmwsl_B1 <- matrix(0,n_samB1,d21)
  for(i in 1:d21){
    gg_B1 <- (dataB1 - mean1A[i])/sqrt(as.vector(sigma1A[i]))
    bmw_B1[which(gg_B1>=0)] <- pnorm(gg_B1[gg_B1>=0])
    bmw_B1[which(gg_B1<0)] <- 1 - pnorm(abs(gg_B1[gg_B1<0]))
    bmwsl_B1[,i] <- bmw_B1*mix_pi_1[i]
  }
  P2_B1 <- rowSums(bmwsl_B1)
  nn<-dim(as.matrix(unique(P2_B1)))[1]
  if(nn<n_samB1){P2_B1<-P2_B1+runif(n_samB1)/1e4}

  dd_B1 <- as.matrix(c(sum(P2_B1),sum(P2_B1^2),sum((P2_B1-0.5)^2)))
  WW2_B1 <- 1/(12*n_samB1) + sum((P2_B1 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  u_B1 <- as.matrix(c(12*n_samB1*((dd_B1[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((dd_B1[2]/n_samB1-1/3)^2),180*n_samB1*((dd_B1[3]/n_samB1-1/12)^2)))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[2])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),D_B1))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])

  ###############  B2  ###################
  dataB2<-sort(dataB2);bmw_B2 <- matrix(0,n_samB2,1); bmwsl_B2 <- matrix(0,n_samB2,d22)
  for(i in 1:d22){
    gg_B2 <- (dataB2 - mean2A[i])/sqrt(sigma2A[i])
    bmw_B2[which(gg_B2>=0)] <- pnorm(gg_B2[gg_B2>=0])
    bmw_B2[which(gg_B2<0)] <- 1 - pnorm(abs(gg_B2[gg_B2<0]))
    bmwsl_B2[,i] <- bmw_B2*mix_pi_2[i]
  }
  P2_B2 <- rowSums(bmwsl_B2)
  nn<-dim(as.matrix(unique(P2_B2)))[1]
  if(nn<n_samB2){P2_B2<-P2_B2+runif(n_samB2)/1e4}

  dd_B2 <- as.matrix(c(sum(P2_B2),sum(P2_B2^2),sum((P2_B2-0.5)^2)))
  WW2_B2 <- 1/(12*n_samB2) + sum((P2_B2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  u_B2 <- as.matrix(c(12*n_samB2*((dd_B2[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((dd_B2[2]/n_samB2-1/3)^2),180*n_samB2*((dd_B2[3]/n_samB2-1/12)^2)))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[2])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),D_B2))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("1MG-EAD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigmaP1,4),round(t(mean1A),4)," "," "," ",round(sigma1A[1],4),round(t(mix_pi_1),4)," "," "," ",
                       round(t(mean2A),4)," "," ",round(sigma2A[1],4),round(t(mix_pi_2),4)," "," ",
                       round(B1[1],4),round(B1[2],4)," ",round(B1[2],4)," "," "," ",round(jj_1,4),round(ll_1*100,4)," "," ",round(jj_2,4),round(ll_2*100,4)," "," ",
                       round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],
                       round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],
                       round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],
                       round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],
                       round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
}
########################### (A4)#####################################
G5BCModelFun[[4]] <- function(K1,logL,df11,df21,df31,df41,df51){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]))
  n_samP1<-dim(dataP1)[1];n_samP2<-dim(dataP2)[1];n_samF1<-dim(dataF1)[1];n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2);sumx4<-sum(dataB1);sumx5<-sum(dataB2)
  ss1<- (n_samP1-1)*as.numeric(var(dataP1));meanP1<- mean(dataP1)
  ss2<- (n_samF1-1)*as.numeric(var(dataF1));meanF1<- mean(dataF1)
  ss3<- (n_samP2-1)*as.numeric(var(dataP2));meanP2<- mean(dataP2)
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  sigmaP1<- sigma0;sigmaF1<- sigma0;sigmaP2<- sigma0
  mean1<- mean(dataB1); sigma1<- as.numeric(var(dataB1))
  mean2<- mean(dataB2); sigma2<- as.numeric(var(dataB2))

  m_esp<-0.0001
  ####################### 1MG-NCD Model ############ (A4)#####################
  d21<-2; d22<-1
  mi_1<- matrix(0.5,d21,1); a1<- sqrt(sigma1/(n_samB1-1))
  sigma1A<- matrix(sigma1/3,d21,1); mean1A<- as.matrix(c(mean1+a1,mean1-a1))

  mi_2<- matrix(1,d22,1)
  sigma2A<- as.matrix(sigma2/3) ; mean2A<- as.matrix(mean2)

  L0 <- logL(n_samP1,1,1,meanP1,sigmaP1,dataP1)+logL(n_samF1,1,1,meanF1,sigmaF1,dataF1)+logL(n_samP2,1,1,meanP2,sigmaP2,dataP2)+logL(n_samB1,d21,mi_1,mean1A,sigma1A,dataB1)+logL(n_samB2,d22,mi_2,mean2A,sigma2A,dataB2)
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1A[i],sqrt(sigma1A[i]))/dmixnorm(dataB1,mean1A,sqrt(sigma1A),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2A[i],sqrt(sigma2A[i]))/dmixnorm(dataB2,mean2A,sqrt(sigma2A),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    ############ CM1-step for means ##############
    n01 <- n_samB1*mix_pi_1; n01[n01<0.000001] <- 0.000001
    n02 <- n_samB2*mix_pi_2; n02[n02<0.000001] <- 0.000001

    meanP1<- (sumx1)/(n_samP1)
    meanF1<- (sumx2)/(n_samF1)
    meanP2<- sumx3/n_samP2
    mean1A<- as.matrix(c(sumwx_B1[1]/n01[1],(sumwx_B1[2])/(n01[2])));mean2A<- (sumwx_B2[1])/(n_samB2)

    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2); ss2<- sum((dataF1-meanF1)^2); ss3<- sum((dataP2-meanP2)^2)
    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1A[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2A[i])^2 }

    sigma1A[1]<- (sum(swx_B1)+sum(swx_B2)+ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2+n_samB1+n_samB2)
    sigma1A<- matrix(sum(swx_B1)/n_samB1,d21,1);sigma2A<- matrix(sum(swx_B2)/n_samB2,d22,1)
    sigmaP1<- ss1/n_samP1; sigmaF1<- ss2/n_samF1; sigmaP2<- ss3/n_samP2
    if(sum(sigma1A < 1e-30)>=1){break}
    if(sum(sigma2A < 1e-30)>=1){break}
    ####################### the stop criterion for iteration #####################################
    L1<- logL(n_samP1,1,1,meanP1,sigmaP1,dataP1)+logL(n_samF1,1,1,meanF1,sigmaF1,dataF1)+logL(n_samP2,1,1,meanP2,sigmaP2,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1A,sigma1A,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2A,sigma2A,dataB2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*7

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,-1),2,2)
  b_line1 <- matrix(c(mean1A[1],mean2A[1]))
  B1 <- solve(aa,b_line1)

  jj_1 <- sigma1 - sigma1A[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigma1

  jj_2 <- sigma2 - sigma2A[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigma2

  ######### hypothesis testing #########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - meanP1)/sqrt(as.vector(sigmaP1))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[2])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),D_P1))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])

  ###############  F1  ###################
  dataF1<-sort(dataF1);bmw_F1 <- matrix(0,n_samF1,1); bmwsl_F1 <- matrix(0,n_samF1,1)

  gg_F1 <- (dataF1 - meanF1)/sqrt(as.vector(sigmaF1))
  bmw_F1[which(gg_F1>=0)] <- pnorm(gg_F1[gg_F1>=0])
  bmw_F1[which(gg_F1<0)] <- 1 - pnorm(abs(gg_F1[gg_F1<0]))
  bmwsl_F1[,1] <- bmw_F1

  P2_F1 <- rowSums(bmwsl_F1)
  nn<-dim(as.matrix(unique(P2_F1)))[1]
  if(nn<n_samF1){P2_F1<-P2_F1+runif(n_samF1)/1e4}

  dd_F1 <- as.matrix(c(sum(P2_F1),sum(P2_F1^2),sum((P2_F1-0.5)^2)))
  WW2_F1 <- 1/(12*n_samF1) + sum((P2_F1 - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  u_F1 <- as.matrix(c(12*n_samF1*((dd_F1[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((dd_F1[2]/n_samF1-1/3)^2),180*n_samF1*((dd_F1[3]/n_samF1-1/12)^2)))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[2])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),D_F1))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - meanP2)/sqrt(as.vector(sigmaP2))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[2])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),D_P2))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])

  ###############  B1  ###################
  dataB1<-sort(dataB1);bmw_B1 <- matrix(0,n_samB1,1); bmwsl_B1 <- matrix(0,n_samB1,d21)
  for(i in 1:d21){
    gg_B1 <- (dataB1 - mean1A[i])/sqrt(sigma1A[i])
    bmw_B1[which(gg_B1>=0)] <- pnorm(gg_B1[gg_B1>=0])
    bmw_B1[which(gg_B1<0)] <- 1 - pnorm(abs(gg_B1[gg_B1<0]))
    bmwsl_B1[,i] <- bmw_B1*mix_pi_1[i]
  }
  P2_B1 <- rowSums(bmwsl_B1)
  nn<-dim(as.matrix(unique(P2_B1)))[1]
  if(nn<n_samB1){P2_B1<-P2_B1+runif(n_samB1)/1e4}

  dd_B1 <- as.matrix(c(sum(P2_B1),sum(P2_B1^2),sum((P2_B1-0.5)^2)))
  WW2_B1 <- 1/(12*n_samB1) + sum((P2_B1 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  u_B1 <- as.matrix(c(12*n_samB1*((dd_B1[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((dd_B1[2]/n_samB1-1/3)^2),180*n_samB1*((dd_B1[3]/n_samB1-1/12)^2)))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[2])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),D_B1))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])

  ###############  B2  ###################
  dataB2<-sort(dataB2);bmw_B2 <- matrix(0,n_samB2,1); bmwsl_B2 <- matrix(0,n_samB2,d22)
  for(i in 1:d22){
    gg_B2 <- (dataB2 - mean2A[i])/sqrt(as.vector(sigma2A[i]))
    bmw_B2[which(gg_B2>=0)] <- pnorm(gg_B2[gg_B2>=0])
    bmw_B2[which(gg_B2<0)] <- 1 - pnorm(abs(gg_B2[gg_B2<0]))
    bmwsl_B2[,i] <- bmw_B2*mix_pi_2[i]
  }
  P2_B2 <- rowSums(bmwsl_B2)
  nn<-dim(as.matrix(unique(P2_B2)))[1]
  if(nn<n_samB2){P2_B2<-P2_B2+runif(n_samB2)/1e4}

  dd_B2 <- as.matrix(c(sum(P2_B2),sum(P2_B2^2),sum((P2_B2-0.5)^2)))
  WW2_B2 <- 1/(12*n_samB2) + sum((P2_B2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  u_B2 <- as.matrix(c(12*n_samB2*((dd_B2[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((dd_B2[2]/n_samB2-1/3)^2),180*n_samB2*((dd_B2[3]/n_samB2-1/12)^2)))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[2])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),D_B2))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("1MG-NCD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigmaP1,4),round(t(mean1A),4)," "," ",round(sigma1A[1],4),round(t(mix_pi_1),4)," "," ",
                       round(t(mean2A),4)," "," "," ",round(sigma2A[1],4),round(t(mix_pi_2),4)," "," "," ",
                       round(B1[1],4),round(B1[2],4)," ",round(-B1[2],4)," "," "," ",round(jj_1,4),round(ll_1*100,4)," "," ",round(jj_2,4),round(ll_2*100,4)," "," ",
                       round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],
                       round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],
                       round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],
                       round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],
                       round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
}
########################### 2MG-AD Model (B1)#####################################
G5BCModelFun[[5]] <- function(K1,logL,df11,df21,df31,df41,df51){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]))
  n_samP1<-dim(dataP1)[1];n_samP2<-dim(dataP2)[1];n_samF1<-dim(dataF1)[1];n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2);sumx4<-sum(dataB1);sumx5<-sum(dataB2)
  ss1<- (n_samP1-1)*as.numeric(var(dataP1));meanP1<- mean(dataP1)
  ss2<- (n_samF1-1)*as.numeric(var(dataF1));meanF1<- mean(dataF1)
  ss3<- (n_samP2-1)*as.numeric(var(dataP2));meanP2<- mean(dataP2)
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  sigmaP1<- sigma0;sigmaF1<- sigma0;sigmaP2<- sigma0
  mean1<- mean(dataB1); sigma1<- as.numeric(var(dataB1))
  mean2<- mean(dataB2); sigma2<- as.numeric(var(dataB2))

  m_esp<-0.0001
  ####################### 2MG-AD Model ############ (B1)#####################
  d21<-4; d22<-4
  mi_1<- matrix(0.25,d21,1); a1<- sqrt(sigma1/(n_samB1-1))
  sigma1B<- matrix(sigma1/3,d21,1); mean1B<- as.matrix(c(mean1+2*a1,mean1+0.8*a1,mean1-0.8*a1,mean1-2*a1))

  mi_2<- matrix(0.25,d22,1); a2<- sqrt(sigma2/(n_samB2-1))
  sigma2B<- matrix(sigma2/3,d22,1); mean2B<- as.matrix(c(mean2+2*a2,mean2+0.8*a2,mean2-0.8*a2,mean2-2*a2))

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1B,sigma1B,dataB1)+logL(n_samB2,d22,mi_2,mean2B,sigma2B,dataB2)
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1B[i],sqrt(sigma1B[i]))/dmixnorm(dataB1,mean1B,sqrt(sigma1B),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2B[i],sqrt(sigma2B[i]))/dmixnorm(dataB2,mean2B,sqrt(sigma2B),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    ############ CM1-step for means ##############
    n01 <- n_samB1*mix_pi_1; n01[n01<0.000001] <- 0.000001
    n02 <- n_samB2*mix_pi_2; n02[n02<0.000001] <- 0.000001

    s0<- matrix(0,8,1)
    s0[1]<- sumx1+sumwx_B1[1]; s0[2]<- n_samP1+n01[1];
    s0[3]<- sumx2+sumwx_B1[4]+sumwx_B2[1]; s0[4]<- n_samF1+n01[4]+n02[1];
    s0[5]<- sumx3+sumwx_B2[4]; s0[6]<- n_samP2+n02[4];
    aa1<- sigma1B[1]/s0[2]+sigma1B[1]/s0[4]+sigma1B[1]/n01[2]+sigma1B[1]/n01[3];
    aa2<- sigma1B[1]/s0[4];
    aa3<- sigma1B[1]/s0[4]+sigma1B[1]/s0[6]+sigma1B[1]/n02[2]+sigma1B[1]/n02[3];
    s0[7]<- s0[1]/s0[2]+s0[3]/s0[4]-sumwx_B1[2]/n01[2]-sumwx_B1[3]/n01[3];
    s0[8]<- s0[3]/s0[4]+s0[5]/s0[6]-sumwx_B2[2]/n02[2]-sumwx_B2[3]/n02[3];
    aa4<- aa1*aa3-aa2*aa2;
    aa6<- s0[7]*aa3-s0[8]*aa2;
    aa7<- aa1*s0[8]-aa2*s0[7];
    rr<- matrix(0,2,1)
    rr[1]<- aa6/aa4;rr[2]<- aa7/aa4;

    meanP1<- (s0[1]-sigma1B[1]*rr[1])/s0[2]
    meanF1<- (s0[3]-sigma1B[1]*(rr[1]+rr[2]))/s0[4]
    meanP2<- (s0[5]-sigma1B[1]*rr[2])/s0[6]
    mean1B<- as.matrix(c(meanP1,(sumwx_B1[2]+sigma1B[2]*rr[1])/n01[2],(sumwx_B1[3]+sigma1B[3]*rr[1])/n01[3],meanF1))
    mean2B<- as.matrix(c(meanF1,(sumwx_B2[2]+sigma2B[2]*rr[2])/n02[2],(sumwx_B2[3]+sigma2B[3]*rr[2])/n02[3],meanP2))

    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2); ss2<- sum((dataF1-meanF1)^2); ss3<- sum((dataP2-meanP2)^2)
    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1B[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2B[i])^2 }

    sigma1B[1]<- (sum(swx_B1)+sum(swx_B2)+ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2+n_samB1+n_samB2)
    sigma1B<- matrix(sigma1B[1],d21,1);sigma2B<- matrix(sigma1B[1],d22,1)
    sigmaP1<- sigma1B[1]; sigmaF1<- sigma1B[1]; sigmaP2<- sigma1B[1]
    if(sum(sigma1B < 1e-30)>=1){break}
    if(sum(sigma2B < 1e-30)>=1){break}
    ####################### the stop criterion for iteration #################################
    L1<- logL(n_samP1,1,1,meanP1,sigmaP1,dataP1)+logL(n_samF1,1,1,meanF1,sigmaF1,dataF1)+logL(n_samP2,1,1,meanP2,sigmaP2,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1B,sigma1B,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2B,sigma2B,dataB2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*8

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1, 1,0,-1,1,0,0,-1, 1,0,-1,0,1,-1,0, 0,1,0,0,1,1,0, 0,1,0,1,0,0,1),7,5)
  b_line1 <- matrix(c(meanP1,meanF1,meanP2,mean1B[2],mean1B[3],mean2B[2],mean2B[3]))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigma1 - sigma1B[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigma1

  jj_2 <- sigma2 - sigma2B[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigma2

  ######### hypothesis testing #########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - meanP1)/sqrt(as.vector(sigmaP1))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[2])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),D_P1))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])

  ###############  F1  ###################
  dataF1<-sort(dataF1);bmw_F1 <- matrix(0,n_samF1,1); bmwsl_F1 <- matrix(0,n_samF1,1)

  gg_F1 <- (dataF1 - meanF1)/sqrt(as.vector(sigmaF1))
  bmw_F1[which(gg_F1>=0)] <- pnorm(gg_F1[gg_F1>=0])
  bmw_F1[which(gg_F1<0)] <- 1 - pnorm(abs(gg_F1[gg_F1<0]))
  bmwsl_F1[,1] <- bmw_F1

  P2_F1 <- rowSums(bmwsl_F1)
  nn<-dim(as.matrix(unique(P2_F1)))[1]
  if(nn<n_samF1){P2_F1<-P2_F1+runif(n_samF1)/1e4}

  dd_F1 <- as.matrix(c(sum(P2_F1),sum(P2_F1^2),sum((P2_F1-0.5)^2)))
  WW2_F1 <- 1/(12*n_samF1) + sum((P2_F1 - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  u_F1 <- as.matrix(c(12*n_samF1*((dd_F1[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((dd_F1[2]/n_samF1-1/3)^2),180*n_samF1*((dd_F1[3]/n_samF1-1/12)^2)))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[2])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),D_F1))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - meanP2)/sqrt(as.vector(sigmaP2))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[2])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),D_P2))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])

  ###############  B1  ###################
  dataB1<-sort(dataB1);bmw_B1 <- matrix(0,n_samB1,1); bmwsl_B1 <- matrix(0,n_samB1,d21)
  for(i in 1:d21){
    gg_B1 <- (dataB1 - mean1B[i])/sqrt(sigma1B[i])
    bmw_B1[which(gg_B1>=0)] <- pnorm(gg_B1[gg_B1>=0])
    bmw_B1[which(gg_B1<0)] <- 1 - pnorm(abs(gg_B1[gg_B1<0]))
    bmwsl_B1[,i] <- bmw_B1*mix_pi_1[i]
  }
  P2_B1 <- rowSums(bmwsl_B1)
  nn<-dim(as.matrix(unique(P2_B1)))[1]
  if(nn<n_samB1){P2_B1<-P2_B1+runif(n_samB1)/1e4}

  dd_B1 <- as.matrix(c(sum(P2_B1),sum(P2_B1^2),sum((P2_B1-0.5)^2)))
  WW2_B1 <- 1/(12*n_samB1) + sum((P2_B1 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  u_B1 <- as.matrix(c(12*n_samB1*((dd_B1[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((dd_B1[2]/n_samB1-1/3)^2),180*n_samB1*((dd_B1[3]/n_samB1-1/12)^2)))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[2])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),D_B1))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])

  ###############  B2  ###################
  dataB2<-sort(dataB2);bmw_B2 <- matrix(0,n_samB2,1); bmwsl_B2 <- matrix(0,n_samB2,d22)
  for(i in 1:d22){
    gg_B2 <- (dataB2 - mean2B[i])/sqrt(sigma2B[i])
    bmw_B2[which(gg_B2>=0)] <- pnorm(gg_B2[gg_B2>=0])
    bmw_B2[which(gg_B2<0)] <- 1 - pnorm(abs(gg_B2[gg_B2<0]))
    bmwsl_B2[,i] <- bmw_B2*mix_pi_2[i]
  }
  P2_B2 <- rowSums(bmwsl_B2)
  nn<-dim(as.matrix(unique(P2_B2)))[1]
  if(nn<n_samB2){P2_B2<-P2_B2+runif(n_samB2)/1e4}

  dd_B2 <- as.matrix(c(sum(P2_B2),sum(P2_B2^2),sum((P2_B2-0.5)^2)))
  WW2_B2 <- 1/(12*n_samB2) + sum((P2_B2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  u_B2 <- as.matrix(c(12*n_samB2*((dd_B2[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((dd_B2[2]/n_samB2-1/3)^2),180*n_samB2*((dd_B2[3]/n_samB2-1/12)^2)))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[2])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),D_B2))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("2MG-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigmaP1,4),round(t(mean1B),4),round(sigma1B[1],4),round(t(mix_pi_1),4),
                       round(t(mean2B),4),round(sigma2B[1],4),round(t(mix_pi_2),4),
                       round(B1[1],4),round(B1[2],4),round(B1[3],4),round(B1[4],4),round(B1[5],4),round(B1[6],4)," ",round(jj_1,4),round(ll_1*100,4)," "," ",round(jj_2,4),round(ll_2*100,4)," "," ",
                       round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],
                       round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],
                       round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],
                       round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],
                       round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
}
########################### (B2)#####################################
G5BCModelFun[[6]] <- function(K1,logL,df11,df21,df31,df41,df51){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]))
  n_samP1<-dim(dataP1)[1];n_samP2<-dim(dataP2)[1];n_samF1<-dim(dataF1)[1];n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2);sumx4<-sum(dataB1);sumx5<-sum(dataB2)
  ss1<- (n_samP1-1)*as.numeric(var(dataP1));meanP1<- mean(dataP1)
  ss2<- (n_samF1-1)*as.numeric(var(dataF1));meanF1<- mean(dataF1)
  ss3<- (n_samP2-1)*as.numeric(var(dataP2));meanP2<- mean(dataP2)
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  sigmaP1<- sigma0;sigmaF1<- sigma0;sigmaP2<- sigma0
  mean1<- mean(dataB1); sigma1<- as.numeric(var(dataB1))
  mean2<- mean(dataB2); sigma2<- as.numeric(var(dataB2))

  m_esp<-0.0001
  ####################### 2MG-A Model ############ (B2)#####################
  d21<-4; d22<-4
  mi_1<- matrix(0.25,d21,1); a1<- sqrt(sigma1/(n_samB1-1))
  sigma1B<- matrix(sigma1/3,d21,1); mean1B<- as.matrix(c(mean1+2*a1,mean1+0.8*a1,mean1-0.8*a1,mean1-2*a1))

  mi_2<- matrix(0.25,d22,1); a2<- sqrt(sigma2/(n_samB2-1))
  sigma2B<- matrix(sigma2/3,d22,1); mean2B<- as.matrix(c(mean2+2*a2,mean2+0.8*a2,mean2-0.8*a2,mean2-2*a2))

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1B,sigma1B,dataB1)+logL(n_samB2,d22,mi_2,mean2B,sigma2B,dataB2)
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1B[i],sqrt(sigma1B[i]))/dmixnorm(dataB1,mean1B,sqrt(sigma1B),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2B[i],sqrt(sigma2B[i]))/dmixnorm(dataB2,mean2B,sqrt(sigma2B),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    ############ CM1-step for means ##############
    n01 <- n_samB1*mix_pi_1; n01[n01<0.000001] <- 0.000001
    n02 <- n_samB2*mix_pi_2; n02[n02<0.000001] <- 0.000001

    s0<- matrix(0,8,1)
    s0[1]<- sumx1+sumwx_B1[1]; s0[2]<- n_samP1+n01[1];
    s0[3]<- sumx2+sumwx_B1[4]+sumwx_B2[1]; s0[4]<- n_samF1+n01[4]+n02[1];
    s0[5]<- sumx3+sumwx_B2[4]; s0[6]<- n_samP2+n02[4];

    hh<- matrix(0,4,4)
    hh[1,1]<- sigma1B[1]/s0[2]+sigma1B[1]/s0[4]+sigma1B[1]/n01[2]+sigma1B[1]/n01[3]
    hh[1,2]<- sigma1B[1]/s0[4]
    hh[1,3]<- sigma1B[1]/s0[2]-2.0*sigma1B[1]/s0[6]
    hh[1,4]<- -sigma1B[1]/n01[2]+sigma1B[1]/n01[3]
    hh[2,2]<- sigma1B[1]/s0[4]+sigma1B[1]/s0[6]+sigma1B[1]/n02[2]+sigma2B[1]/n02[3]
    hh[2,3]<- -2.0*sigma1B[1]/s0[4]+sigma1B[1]/s0[6]
    hh[2,4]<- sigma1B[1]/n02[2]-sigma1B[1]/n02[3]
    hh[3,3]<- sigma1B[1]/s0[2]+4.0*sigma1B[1]/s0[4]+sigma1B[1]/s0[6]
    hh[3,4]<- 0
    hh[4,4]<- sigma1B[1]/n01[2]+sigma1B[1]/n01[3]+sigma1B[1]/n02[2]+sigma1B[1]/n02[3]
    for(i in 2:4){
      for(j in 1:(i-1)){
        hh[i,j]<- hh[j,i]
      }
    }
    b_line<- matrix(0,4,1)
    b_line[1]<- s0[1]/s0[2]+s0[3]/s0[4]-sumwx_B1[2]/n01[2]-sumwx_B1[3]/n01[3]
    b_line[2]<- s0[3]/s0[4]+s0[5]/s0[6]-sumwx_B2[2]/n02[2]-sumwx_B2[3]/n02[3]
    b_line[3]<- s0[1]/s0[2]-2.0*s0[3]/s0[4]+s0[5]/s0[6]
    b_line[4]<- sumwx_B1[2]/n01[2]-sumwx_B1[3]/n01[3]-sumwx_B2[2]/n02[2]+sumwx_B2[3]/n02[3]
    B <- solve(hh,b_line)

    meanP1<- (s0[1]-sigma1B[1]*(B[1]+B[3]))/s0[2]
    meanF1<- (s0[3]-sigma1B[1]*(B[1]+B[2]-2.0*B[3]))/s0[4]
    meanP2<- (s0[5]-sigma1B[1]*(B[2]+B[3]))/s0[6]
    mean1B<- as.matrix(c(meanP1,(sumwx_B1[2]+sigma1B[1]*(B[1]-B[4]))/n01[2],(sumwx_B1[3]+sigma1B[1]*(B[1]+B[4]))/n01[3],meanF1))
    mean2B<- as.matrix(c(meanF1,(sumwx_B2[2]+sigma1B[1]*(B[2]+B[4]))/n02[2],(sumwx_B2[3]+sigma1B[1]*(B[2]-B[4]))/n02[3],meanP2))

    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2); ss2<- sum((dataF1-meanF1)^2); ss3<- sum((dataP2-meanP2)^2)
    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1B[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2B[i])^2 }

    sigma1B[1]<- (sum(swx_B1)+sum(swx_B2)+ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2+n_samB1+n_samB2)
    sigma1B<- matrix(sigma1B[1],d21,1);sigma2B<- matrix(sigma1B[1],d22,1)
    sigmaP1<- sigma1B[1]; sigmaF1<- sigma1B[1]; sigmaP2<- sigma1B[1]
    if(sum(sigma1B < 1e-30)>=1){break}
    if(sum(sigma2B < 1e-30)>=1){break}
    ####################### the stop criterion for iteration #################################
    L1<- logL(n_samP1,1,1,meanP1,sigmaP1,dataP1)+logL(n_samF1,1,1,meanF1,sigmaF1,dataF1)+logL(n_samP2,1,1,meanP2,sigmaP2,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1B,sigma1B,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2B,sigma2B,dataB2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*7

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1, 1,0,-1,1,0,0,-1, 1,0,-1,0,1,-1,0),7,3)
  b_line1 <- matrix(c(meanP1,meanF1,meanP2,mean1B[2],mean1B[3],mean2B[2],mean2B[3]))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigma1 - sigma1B[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigma1

  jj_2 <- sigma2 - sigma2B[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigma2

  ######### hypothesis testing #########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - meanP1)/sqrt(as.vector(sigmaP1))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[2])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),D_P1))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])

  ###############  F1  ###################
  dataF1<-sort(dataF1);bmw_F1 <- matrix(0,n_samF1,1); bmwsl_F1 <- matrix(0,n_samF1,1)

  gg_F1 <- (dataF1 - meanF1)/sqrt(as.vector(sigmaF1))
  bmw_F1[which(gg_F1>=0)] <- pnorm(gg_F1[gg_F1>=0])
  bmw_F1[which(gg_F1<0)] <- 1 - pnorm(abs(gg_F1[gg_F1<0]))
  bmwsl_F1[,1] <- bmw_F1

  P2_F1 <- rowSums(bmwsl_F1)
  nn<-dim(as.matrix(unique(P2_F1)))[1]
  if(nn<n_samF1){P2_F1<-P2_F1+runif(n_samF1)/1e4}

  dd_F1 <- as.matrix(c(sum(P2_F1),sum(P2_F1^2),sum((P2_F1-0.5)^2)))
  WW2_F1 <- 1/(12*n_samF1) + sum((P2_F1 - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  u_F1 <- as.matrix(c(12*n_samF1*((dd_F1[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((dd_F1[2]/n_samF1-1/3)^2),180*n_samF1*((dd_F1[3]/n_samF1-1/12)^2)))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[2])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),D_F1))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - meanP2)/sqrt(as.vector(sigmaP2))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[2])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),D_P2))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])

  ###############  B1  ###################
  dataB1<-sort(dataB1);bmw_B1 <- matrix(0,n_samB1,1); bmwsl_B1 <- matrix(0,n_samB1,d21)
  for(i in 1:d21){
    gg_B1 <- (dataB1 - mean1B[i])/sqrt(sigma1B[i])
    bmw_B1[which(gg_B1>=0)] <- pnorm(gg_B1[gg_B1>=0])
    bmw_B1[which(gg_B1<0)] <- 1 - pnorm(abs(gg_B1[gg_B1<0]))
    bmwsl_B1[,i] <- bmw_B1*mix_pi_1[i]
  }
  P2_B1 <- rowSums(bmwsl_B1)
  nn<-dim(as.matrix(unique(P2_B1)))[1]
  if(nn<n_samB1){P2_B1<-P2_B1+runif(n_samB1)/1e4}

  dd_B1 <- as.matrix(c(sum(P2_B1),sum(P2_B1^2),sum((P2_B1-0.5)^2)))
  WW2_B1 <- 1/(12*n_samB1) + sum((P2_B1 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  u_B1 <- as.matrix(c(12*n_samB1*((dd_B1[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((dd_B1[2]/n_samB1-1/3)^2),180*n_samB1*((dd_B1[3]/n_samB1-1/12)^2)))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[2])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),D_B1))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])

  ###############  B2  ###################
  dataB2<-sort(dataB2);bmw_B2 <- matrix(0,n_samB2,1); bmwsl_B2 <- matrix(0,n_samB2,d22)
  for(i in 1:d22){
    gg_B2 <- (dataB2 - mean2B[i])/sqrt(sigma2B[i])
    bmw_B2[which(gg_B2>=0)] <- pnorm(gg_B2[gg_B2>=0])
    bmw_B2[which(gg_B2<0)] <- 1 - pnorm(abs(gg_B2[gg_B2<0]))
    bmwsl_B2[,i] <- bmw_B2*mix_pi_2[i]
  }
  P2_B2 <- rowSums(bmwsl_B2)
  nn<-dim(as.matrix(unique(P2_B2)))[1]
  if(nn<n_samB2){P2_B2<-P2_B2+runif(n_samB2)/1e4}

  dd_B2 <- as.matrix(c(sum(P2_B2),sum(P2_B2^2),sum((P2_B2-0.5)^2)))
  WW2_B2 <- 1/(12*n_samB2) + sum((P2_B2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  u_B2 <- as.matrix(c(12*n_samB2*((dd_B2[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((dd_B2[2]/n_samB2-1/3)^2),180*n_samB2*((dd_B2[3]/n_samB2-1/12)^2)))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[2])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),D_B2))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("2MG-A",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigmaP1,4),round(t(mean1B),4),round(sigma1B[1],4),round(t(mix_pi_1),4),
                       round(t(mean2B),4),round(sigma2B[1],4),round(t(mix_pi_2),4),
                       round(B1[1],4),round(B1[2],4),round(B1[3],4)," "," "," "," ",round(jj_1,4),round(ll_1*100,4)," "," ",round(jj_2,4),round(ll_2*100,4)," "," ",
                       round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],
                       round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],
                       round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],
                       round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],
                       round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
}
########################### (B3)#####################################
G5BCModelFun[[7]] <- function(K1,logL,df11,df21,df31,df41,df51){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]))
  n_samP1<-dim(dataP1)[1];n_samP2<-dim(dataP2)[1];n_samF1<-dim(dataF1)[1];n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2);sumx4<-sum(dataB1);sumx5<-sum(dataB2)
  ss1<- (n_samP1-1)*as.numeric(var(dataP1));meanP1<- mean(dataP1)
  ss2<- (n_samF1-1)*as.numeric(var(dataF1));meanF1<- mean(dataF1)
  ss3<- (n_samP2-1)*as.numeric(var(dataP2));meanP2<- mean(dataP2)
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  sigmaP1<- sigma0;sigmaF1<- sigma0;sigmaP2<- sigma0
  mean1<- mean(dataB1); sigma1<- as.numeric(var(dataB1))
  mean2<- mean(dataB2); sigma2<- as.numeric(var(dataB2))

  m_esp<-0.0001
  ####################### 2MG-EA Model ############ (B3)#####################
  d21<-3; d22<-3
  mi_1<- as.matrix(c(0.25,0.5,0.25)); a1<- sqrt(sigma1/(n_samB1-1))
  sigma1B<- matrix(sigma1/3,d21,1); mean1B<- as.matrix(c(mean1+2*a1,mean1,mean1-2*a1))

  mi_2<- as.matrix(c(0.25,0.5,0.25)); a2<- sqrt(sigma2/(n_samB2-1))
  sigma2B<- matrix(sigma2/3,d22,1); mean2B<- as.matrix(c(mean2+2*a2,mean2,mean2-2*a2))

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1B,sigma1B,dataB1)+logL(n_samB2,d22,mi_2,mean2B,sigma2B,dataB2)
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1B[i],sqrt(sigma1B[i]))/dmixnorm(dataB1,mean1B,sqrt(sigma1B),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2B[i],sqrt(sigma2B[i]))/dmixnorm(dataB2,mean2B,sqrt(sigma2B),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    ############ CM1-step for means ##############
    n01 <- n_samB1*mix_pi_1; n01[n01<0.000001] <- 0.000001
    n02 <- n_samB2*mix_pi_2; n02[n02<0.000001] <- 0.000001

    s0<- matrix(0,6,1)
    s0[1]<- sumx1+sumwx_B1[1]; s0[2]<- n_samP1+n01[1];
    s0[3]<- sumx2+sumwx_B1[3]+sumwx_B2[1]; s0[4]<- n_samF1+n01[3]+n02[1];
    s0[5]<- sumx3+sumwx_B2[3]; s0[6]<- n_samP2+n02[3];

    hh<- matrix(0,3,3)
    hh[1,1]<- sigma1B[1]/s0[2]+4.0*sigma1B[1]/s0[4]+sigma1B[1]/s0[6]
    hh[1,2]<- -4.0*sigma1B[1]/s0[4]
    hh[1,3]<- sigma1B[1]/s0[2]+2.0*sigma1B[1]/s0[4]
    hh[2,2]<- 4.0*sigma1B[1]/s0[4]+sigma1B[1]/n01[2]+sigma1B[1]/n02[2]
    hh[2,3]<- -2.0*sigma1B[1]/s0[4]+sigma1B[1]/n01[2]-sigma1B[1]/n02[2]
    hh[3,3]<- sigma1B[1]/s0[2]+sigma1B[1]/s0[4]+sigma1B[1]/n01[2]+sigma1B[1]/n02[2]
    for(i in 2:3){
      for(j in 1:(i-1)){
        hh[i,j]<- hh[j,i]
      }
    }
    b_line<- matrix(0,3,1)
    b_line[1]<- s0[1]/s0[2]-2.0*s0[3]/s0[4]+s0[5]/s0[6]
    b_line[2]<- 2.0*s0[3]/s0[4]-sumwx_B1[2]/n01[2]-sumwx_B2[2]/n02[2]
    b_line[3]<- s0[1]/s0[2]-s0[3]/s0[4]-sumwx_B1[2]/n01[2]+sumwx_B2[2]/n02[2]

    B <- solve(hh,b_line)

    meanP1<- (s0[1]-sigma1B[1]*(B[1]+B[3]))/s0[2]
    meanF1<- (s0[3]+sigma1B[1]*(2.0*B[1]-2.0*B[2]+B[3]))/s0[4]
    meanP2<- (s0[5]-sigma1B[1]*B[1])/s0[6]
    mean1B<- as.matrix(c(meanP1,(sumwx_B1[2]+sigma1B[1]*(B[2]+B[3]))/n01[2],meanF1))
    mean2B<- as.matrix(c(meanF1,(sumwx_B2[2]+sigma1B[1]*(B[2]-B[3]))/n02[2],meanP2))

    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2); ss2<- sum((dataF1-meanF1)^2); ss3<- sum((dataP2-meanP2)^2)
    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1B[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2B[i])^2 }

    sigma1B[1]<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2)
    sigma1B<- matrix(swx_B1/n01,d21,1);sigma2B<- matrix(swx_B2/n02,d22,1)
    sigmaP1<- ss1/n_samP1; sigmaF1<- ss2/n_samF1; sigmaP2<- ss3/n_samP2
    if(sum(sigma1B < 1e-30)>=1){break}
    if(sum(sigma2B < 1e-30)>=1){break}
    ####################### the stop criterion for iteration #################################
    L1<- logL(n_samP1,1,1,meanP1,sigmaP1,dataP1)+logL(n_samF1,1,1,meanF1,sigmaF1,dataF1)+logL(n_samP2,1,1,meanP2,sigmaP2,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1B,sigma1B,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2B,sigma2B,dataB2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*7

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1, 2,0,-2,1,-1),5,2)
  b_line1 <- matrix(c(meanP1,meanF1,meanP2,mean1B[2],mean2B[2]))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigma1 - sigma1B[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigma1

  jj_2 <- sigma2 - sigma2B[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigma2

  ######### hypothesis testing #########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - meanP1)/sqrt(as.vector(sigmaP1))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[2])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),D_P1))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])

  ###############  F1  ###################
  dataF1<-sort(dataF1);bmw_F1 <- matrix(0,n_samF1,1); bmwsl_F1 <- matrix(0,n_samF1,1)

  gg_F1 <- (dataF1 - meanF1)/sqrt(as.vector(sigmaF1))
  bmw_F1[which(gg_F1>=0)] <- pnorm(gg_F1[gg_F1>=0])
  bmw_F1[which(gg_F1<0)] <- 1 - pnorm(abs(gg_F1[gg_F1<0]))
  bmwsl_F1[,1] <- bmw_F1

  P2_F1 <- rowSums(bmwsl_F1)
  nn<-dim(as.matrix(unique(P2_F1)))[1]
  if(nn<n_samF1){P2_F1<-P2_F1+runif(n_samF1)/1e4}

  dd_F1 <- as.matrix(c(sum(P2_F1),sum(P2_F1^2),sum((P2_F1-0.5)^2)))
  WW2_F1 <- 1/(12*n_samF1) + sum((P2_F1 - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  u_F1 <- as.matrix(c(12*n_samF1*((dd_F1[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((dd_F1[2]/n_samF1-1/3)^2),180*n_samF1*((dd_F1[3]/n_samF1-1/12)^2)))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[2])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),D_F1))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - meanP2)/sqrt(as.vector(sigmaP2))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[2])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),D_P2))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])

  ###############  B1  ###################
  dataB1<-sort(dataB1);bmw_B1 <- matrix(0,n_samB1,1); bmwsl_B1 <- matrix(0,n_samB1,d21)
  for(i in 1:d21){
    gg_B1 <- (dataB1 - mean1B[i])/sqrt(sigma1B[i])
    bmw_B1[which(gg_B1>=0)] <- pnorm(gg_B1[gg_B1>=0])
    bmw_B1[which(gg_B1<0)] <- 1 - pnorm(abs(gg_B1[gg_B1<0]))
    bmwsl_B1[,i] <- bmw_B1*mix_pi_1[i]
  }
  P2_B1 <- rowSums(bmwsl_B1)
  nn<-dim(as.matrix(unique(P2_B1)))[1]
  if(nn<n_samB1){P2_B1<-P2_B1+runif(n_samB1)/1e4}

  dd_B1 <- as.matrix(c(sum(P2_B1),sum(P2_B1^2),sum((P2_B1-0.5)^2)))
  WW2_B1 <- 1/(12*n_samB1) + sum((P2_B1 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  u_B1 <- as.matrix(c(12*n_samB1*((dd_B1[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((dd_B1[2]/n_samB1-1/3)^2),180*n_samB1*((dd_B1[3]/n_samB1-1/12)^2)))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[2])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),D_B1))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])

  ###############  B2  ###################
  dataB2<-sort(dataB2);bmw_B2 <- matrix(0,n_samB2,1); bmwsl_B2 <- matrix(0,n_samB2,d22)
  for(i in 1:d22){
    gg_B2 <- (dataB2 - mean2B[i])/sqrt(sigma2B[i])
    bmw_B2[which(gg_B2>=0)] <- pnorm(gg_B2[gg_B2>=0])
    bmw_B2[which(gg_B2<0)] <- 1 - pnorm(abs(gg_B2[gg_B2<0]))
    bmwsl_B2[,i] <- bmw_B2*mix_pi_2[i]
  }
  P2_B2 <- rowSums(bmwsl_B2)
  nn<-dim(as.matrix(unique(P2_B2)))[1]
  if(nn<n_samB2){P2_B2<-P2_B2+runif(n_samB2)/1e4}

  dd_B2 <- as.matrix(c(sum(P2_B2),sum(P2_B2^2),sum((P2_B2-0.5)^2)))
  WW2_B2 <- 1/(12*n_samB2) + sum((P2_B2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  u_B2 <- as.matrix(c(12*n_samB2*((dd_B2[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((dd_B2[2]/n_samB2-1/3)^2),180*n_samB2*((dd_B2[3]/n_samB2-1/12)^2)))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[2])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),D_B2))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("2MG-EA",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigmaP1,4),round(t(mean1B),4)," ",round(sigma1B[1],4),round(t(mix_pi_1),4)," ",
                       round(t(mean2B),4)," ",round(sigma2B[1],4),round(t(mix_pi_2),4)," ",
                       round(B1[1],4),round(B1[2],4),round(B1[2],4)," "," "," "," ",round(jj_1,4),round(ll_1*100,4)," "," ",round(jj_2,4),round(ll_2*100,4)," "," ",
                       round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],
                       round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],
                       round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],
                       round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],
                       round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
}
########################### (B4)#####################################
G5BCModelFun[[8]] <- function(K1,logL,df11,df21,df31,df41,df51){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]))
  n_samP1<-dim(dataP1)[1];n_samP2<-dim(dataP2)[1];n_samF1<-dim(dataF1)[1];n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2);sumx4<-sum(dataB1);sumx5<-sum(dataB2)
  ss1<- (n_samP1-1)*as.numeric(var(dataP1));meanP1<- mean(dataP1)
  ss2<- (n_samF1-1)*as.numeric(var(dataF1));meanF1<- mean(dataF1)
  ss3<- (n_samP2-1)*as.numeric(var(dataP2));meanP2<- mean(dataP2)
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  sigmaP1<- sigma0;sigmaF1<- sigma0;sigmaP2<- sigma0
  mean1<- mean(dataB1); sigma1<- as.numeric(var(dataB1))
  mean2<- mean(dataB2); sigma2<- as.numeric(var(dataB2))

  m_esp<-0.0001
  ####################### 2MG-CD Model ############ (B4)#####################
  d21<-1; d22<-4
  mi_1<- as.matrix(1)
  sigma1B<- as.matrix(sigma1/3); mean1B<- as.matrix(mean1)

  mi_2<- matrix(0.25,d22,1); a2<- sqrt(sigma2/(n_samB2-1))
  sigma2B<- matrix(sigma2/3,d22,1); mean2B<- as.matrix(c(mean2+2*a2,mean2+0.5*a2,mean2-0.5*a2,mean2-2*a2))

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1B,sigma1B,dataB1)+logL(n_samB2,d22,mi_2,mean2B,sigma2B,dataB2)
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1B[i],sqrt(sigma1B[i]))/dmixnorm(dataB1,mean1B,sqrt(sigma1B),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2B[i],sqrt(sigma2B[i]))/dmixnorm(dataB2,mean2B,sqrt(sigma2B),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    ############ CM1-step for means ##############
    n01 <- n_samB1*mix_pi_1; n01[n01<0.000001] <- 0.000001
    n02 <- n_samB2*mix_pi_2; n02[n02<0.000001] <- 0.000001

    s0<- matrix(0,4,1)
    s0[1]<- sumx1+sumx2+sumwx_B1[1]+sumwx_B2[1]; s0[2]<- n_samP1+n_samF1+n_samB1+n02[1]
    s0[3]<- sumx3+sumwx_B2[4]; s0[4]<- n_samF1+n02[4]
    aa3<- (s0[1]/s0[2]+s0[3]/s0[4]-sumwx_B2[2]/n02[2]-sumwx_B2[3]/n02[3])/(sigma1B[1]/s0[2]+sigma1B[1]/s0[4]+sigma1B[1]/n02[2]+sigma1B[1]/n02[3])

    meanP1<- (s0[1]-sigma1B[1]*aa3)/s0[2]
    meanF1<- meanP1
    meanP2<- (s0[3]-sigma1B[1]*aa3)/s0[4]
    mean1B<- meanP1
    mean2B<- as.matrix(c(meanP1,(sumwx_B2[2]+sigma1B[1]*aa3)/n02[2],(sumwx_B2[3]+sigma1B[1]*aa3)/n02[3],meanP2))

    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2); ss2<- sum((dataF1-meanF1)^2); ss3<- sum((dataP2-meanP2)^2)
    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1B[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2B[i])^2 }

    sigma1B[1]<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2)
    sigma1B<- matrix(sum(swx_B1)/n_samB1,d21,1);sigma2B<- matrix(sum(swx_B2)/n_samB2,d22,1)
    sigmaP1<- sigma1B[1]; sigmaF1<- sigma1B[1]; sigmaP2<- sigma1B[1]
    if(sum(sigma1B < 1e-30)>=1){break}
    if(sum(sigma2B < 1e-30)>=1){break}
    ####################### the stop criterion for iteration #################################
    L1<- logL(n_samP1,1,1,meanP1,sigmaP1,dataP1)+logL(n_samF1,1,1,meanF1,sigmaF1,dataF1)+logL(n_samP2,1,1,meanP2,sigmaP2,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1B,sigma1B,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2B,sigma2B,dataB2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*7

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1, 1,-1,1,-1, 1,-1,-1,1),4,3)
  b_line1 <- matrix(c(meanP1,meanP2,mean2B[2],mean2B[3]))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigma1 - sigma1B[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigma1

  jj_2 <- sigma2 - sigma2B[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigma2

  ######### hypothesis testing #########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - meanP1)/sqrt(as.vector(sigmaP1))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[2])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),D_P1))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])

  ###############  F1  ###################
  dataF1<-sort(dataF1);bmw_F1 <- matrix(0,n_samF1,1); bmwsl_F1 <- matrix(0,n_samF1,1)

  gg_F1 <- (dataF1 - meanF1)/sqrt(as.vector(sigmaF1))
  bmw_F1[which(gg_F1>=0)] <- pnorm(gg_F1[gg_F1>=0])
  bmw_F1[which(gg_F1<0)] <- 1 - pnorm(abs(gg_F1[gg_F1<0]))
  bmwsl_F1[,1] <- bmw_F1

  P2_F1 <- rowSums(bmwsl_F1)
  nn<-dim(as.matrix(unique(P2_F1)))[1]
  if(nn<n_samF1){P2_F1<-P2_F1+runif(n_samF1)/1e4}

  dd_F1 <- as.matrix(c(sum(P2_F1),sum(P2_F1^2),sum((P2_F1-0.5)^2)))
  WW2_F1 <- 1/(12*n_samF1) + sum((P2_F1 - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  u_F1 <- as.matrix(c(12*n_samF1*((dd_F1[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((dd_F1[2]/n_samF1-1/3)^2),180*n_samF1*((dd_F1[3]/n_samF1-1/12)^2)))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[2])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),D_F1))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - meanP2)/sqrt(as.vector(sigmaP2))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[2])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),D_P2))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])

  ###############  B1  ###################
  dataB1<-sort(dataB1);bmw_B1 <- matrix(0,n_samB1,1); bmwsl_B1 <- matrix(0,n_samB1,d21)
  for(i in 1:d21){
    gg_B1 <- (dataB1 - mean1B[i])/sqrt(as.vector(sigma1B[i]))
    bmw_B1[which(gg_B1>=0)] <- pnorm(gg_B1[gg_B1>=0])
    bmw_B1[which(gg_B1<0)] <- 1 - pnorm(abs(gg_B1[gg_B1<0]))
    bmwsl_B1[,i] <- bmw_B1*mix_pi_1[i]
  }
  P2_B1 <- rowSums(bmwsl_B1)
  nn<-dim(as.matrix(unique(P2_B1)))[1]
  if(nn<n_samB1){P2_B1<-P2_B1+runif(n_samB1)/1e4}

  dd_B1 <- as.matrix(c(sum(P2_B1),sum(P2_B1^2),sum((P2_B1-0.5)^2)))
  WW2_B1 <- 1/(12*n_samB1) + sum((P2_B1 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  u_B1 <- as.matrix(c(12*n_samB1*((dd_B1[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((dd_B1[2]/n_samB1-1/3)^2),180*n_samB1*((dd_B1[3]/n_samB1-1/12)^2)))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[2])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),D_B1))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])

  ###############  B2  ###################
  dataB2<-sort(dataB2);bmw_B2 <- matrix(0,n_samB2,1); bmwsl_B2 <- matrix(0,n_samB2,d22)
  for(i in 1:d22){
    gg_B2 <- (dataB2 - mean2B[i])/sqrt(sigma2B[i])
    bmw_B2[which(gg_B2>=0)] <- pnorm(gg_B2[gg_B2>=0])
    bmw_B2[which(gg_B2<0)] <- 1 - pnorm(abs(gg_B2[gg_B2<0]))
    bmwsl_B2[,i] <- bmw_B2*mix_pi_2[i]
  }
  P2_B2 <- rowSums(bmwsl_B2)
  nn<-dim(as.matrix(unique(P2_B2)))[1]
  if(nn<n_samB2){P2_B2<-P2_B2+runif(n_samB2)/1e4}

  dd_B2 <- as.matrix(c(sum(P2_B2),sum(P2_B2^2),sum((P2_B2-0.5)^2)))
  WW2_B2 <- 1/(12*n_samB2) + sum((P2_B2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  u_B2 <- as.matrix(c(12*n_samB2*((dd_B2[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((dd_B2[2]/n_samB2-1/3)^2),180*n_samB2*((dd_B2[3]/n_samB2-1/12)^2)))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[2])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),D_B2))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("2MG-CD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigmaP1,4),round(t(mean1B),4)," "," "," ",round(sigma1B[1],4),round(t(mix_pi_1),4)," "," "," ",
                       round(t(mean2B),4),round(sigma2B[1],4),round(t(mix_pi_2),4),
                       round(B1[1],4),round(B1[2],4),round(B1[3],4),round(B1[2],4),round(B1[3],4)," "," ",round(jj_1,4),round(ll_1*100,4)," "," ",round(jj_2,4),round(ll_2*100,4)," "," ",
                       round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],
                       round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],
                       round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],
                       round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],
                       round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
}
########################### (B5)#####################################
G5BCModelFun[[9]] <- function(K1,logL,df11,df21,df31,df41,df51){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]))
  n_samP1<-dim(dataP1)[1];n_samP2<-dim(dataP2)[1];n_samF1<-dim(dataF1)[1];n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2);sumx4<-sum(dataB1);sumx5<-sum(dataB2)
  ss1<- (n_samP1-1)*as.numeric(var(dataP1));meanP1<- mean(dataP1)
  ss2<- (n_samF1-1)*as.numeric(var(dataF1));meanF1<- mean(dataF1)
  ss3<- (n_samP2-1)*as.numeric(var(dataP2));meanP2<- mean(dataP2)
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  sigmaP1<- sigma0;sigmaF1<- sigma0;sigmaP2<- sigma0
  mean1<- mean(dataB1); sigma1<- as.numeric(var(dataB1))
  mean2<- mean(dataB2); sigma2<- as.numeric(var(dataB2))

  m_esp<-0.0001
  ####################### 2MG-EAD Model ############ (B5)#####################
  d21<-1; d22<-3
  mi_1<- as.matrix(1)
  sigma1B<- as.matrix(sigma1/3) ; mean1B<- as.matrix(mean1)

  mi_2<- as.matrix(c(0.25,0.5,0.25)); a2<- sqrt(sigma2/(n_samB2-1))
  sigma2B<- matrix(sigma2/3,d22,1); mean2B<- as.matrix(c(mean2+2*a2,mean2,mean2-2*a2))

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1B,sigma1B,dataB1)+logL(n_samB2,d22,mi_2,mean2B,sigma2B,dataB2)
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1B[i],sqrt(sigma1B[i]))/dmixnorm(dataB1,mean1B,sqrt(sigma1B),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2B[i],sqrt(sigma2B[i]))/dmixnorm(dataB2,mean2B,sqrt(sigma2B),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    ############ CM1-step for means ##############
    n01 <- n_samB1*mix_pi_1; n01[n01<0.000001] <- 0.000001
    n02 <- n_samB2*mix_pi_2; n02[n02<0.000001] <- 0.000001

    s0<- matrix(0,4,1)
    s0[1]<- sumx1+sumx2+sumwx_B1[1]+sumwx_B2[1]; s0[2]<-n_samP1+n_samF1+n_samB1+n02[1]
    s0[3]<- sumx3+sumwx_B2[3]; s0[4]<- n_samP2+n02[3]

    aa3<- (s0[1]/s0[2]+s0[3]/s0[4]-2.0*sumwx_B2[2]/n02[2])/(sigma1B[1]/s0[2]+sigma1B[1]/s0[4]+4.0*sigma1B[1]/n02[2])

    meanP1<- (s0[1]-sigma1B[1]*aa3)/s0[2]
    meanF1<- meanP1
    meanP2<- (s0[3]-sigma1B[1]*aa3)/s0[4]
    mean1B<- meanP1
    mean2B<- as.matrix(c(meanP1,(sumwx_B2[2]+sigma1B[1]*2.0*aa3)/n02[2],meanP2))
    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2); ss2<- sum((dataF1-meanF1)^2); ss3<- sum((dataP2-meanP2)^2)
    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1B[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2B[i])^2 }

    sigma1B[1]<- (sum(swx_B1)+sum(swx_B2)+ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2+n_samB1+n_samB2)
    sigma1B<- matrix(sigma1B[1],d21,1);sigma2B<- matrix(sigma1B[1],d22,1)
    sigmaP1<- sigma1B[1]; sigmaF1<- sigma1B[1]; sigmaP2<- sigma1B[1]
    if(sum(sigma1B < 1e-30)>=1){break}
    if(sum(sigma2B < 1e-30)>=1){break}
    ####################### the stop criterion for iteration ############################################################
    L1<- logL(n_samP1,1,1,meanP1,sigmaP1,dataP1)+logL(n_samF1,1,1,meanF1,sigmaF1,dataF1)+logL(n_samP2,1,1,meanP2,sigmaP2,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1B,sigma1B,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2B,sigma2B,dataB2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*7

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1, 2,-2,0),3,2)
  b_line1 <- matrix(c(meanP1,meanP2,mean2B[2]))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigma1 - sigma1B[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigma1

  jj_2 <- sigma2 - sigma2B[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigma2

  ######### hypothesis testing #########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - meanP1)/sqrt(as.vector(sigmaP1))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[2])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),D_P1))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])

  ###############  F1  ###################
  dataF1<-sort(dataF1);bmw_F1 <- matrix(0,n_samF1,1); bmwsl_F1 <- matrix(0,n_samF1,1)

  gg_F1 <- (dataF1 - meanF1)/sqrt(as.vector(sigmaF1))
  bmw_F1[which(gg_F1>=0)] <- pnorm(gg_F1[gg_F1>=0])
  bmw_F1[which(gg_F1<0)] <- 1 - pnorm(abs(gg_F1[gg_F1<0]))
  bmwsl_F1[,1] <- bmw_F1

  P2_F1 <- rowSums(bmwsl_F1)
  nn<-dim(as.matrix(unique(P2_F1)))[1]
  if(nn<n_samF1){P2_F1<-P2_F1+runif(n_samF1)/1e4}

  dd_F1 <- as.matrix(c(sum(P2_F1),sum(P2_F1^2),sum((P2_F1-0.5)^2)))
  WW2_F1 <- 1/(12*n_samF1) + sum((P2_F1 - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  u_F1 <- as.matrix(c(12*n_samF1*((dd_F1[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((dd_F1[2]/n_samF1-1/3)^2),180*n_samF1*((dd_F1[3]/n_samF1-1/12)^2)))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[2])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),D_F1))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - meanP2)/sqrt(as.vector(sigmaP2))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[2])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),D_P2))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])

  ###############  B1  ###################
  dataB1<-sort(dataB1);bmw_B1 <- matrix(0,n_samB1,1); bmwsl_B1 <- matrix(0,n_samB1,d21)
  for(i in 1:d21){
    gg_B1 <- (dataB1 - mean1B[i])/sqrt(as.vector(sigma1B[i]))
    bmw_B1[which(gg_B1>=0)] <- pnorm(gg_B1[gg_B1>=0])
    bmw_B1[which(gg_B1<0)] <- 1 - pnorm(abs(gg_B1[gg_B1<0]))
    bmwsl_B1[,i] <- bmw_B1*mix_pi_1[i]
  }
  P2_B1 <- rowSums(bmwsl_B1)
  nn<-dim(as.matrix(unique(P2_B1)))[1]
  if(nn<n_samB1){P2_B1<-P2_B1+runif(n_samB1)/1e4}

  dd_B1 <- as.matrix(c(sum(P2_B1),sum(P2_B1^2),sum((P2_B1-0.5)^2)))
  WW2_B1 <- 1/(12*n_samB1) + sum((P2_B1 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  u_B1 <- as.matrix(c(12*n_samB1*((dd_B1[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((dd_B1[2]/n_samB1-1/3)^2),180*n_samB1*((dd_B1[3]/n_samB1-1/12)^2)))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[2])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),D_B1))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])

  ###############  B2  ###################
  dataB2<-sort(dataB2);bmw_B2 <- matrix(0,n_samB2,1); bmwsl_B2 <- matrix(0,n_samB2,d22)
  for(i in 1:d22){
    gg_B2 <- (dataB2 - mean2B[i])/sqrt(sigma2B[i])
    bmw_B2[which(gg_B2>=0)] <- pnorm(gg_B2[gg_B2>=0])
    bmw_B2[which(gg_B2<0)] <- 1 - pnorm(abs(gg_B2[gg_B2<0]))
    bmwsl_B2[,i] <- bmw_B2*mix_pi_2[i]
  }
  P2_B2 <- rowSums(bmwsl_B2)
  nn<-dim(as.matrix(unique(P2_B2)))[1]
  if(nn<n_samB2){P2_B2<-P2_B2+runif(n_samB2)/1e4}

  dd_B2 <- as.matrix(c(sum(P2_B2),sum(P2_B2^2),sum((P2_B2-0.5)^2)))
  WW2_B2 <- 1/(12*n_samB2) + sum((P2_B2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  u_B2 <- as.matrix(c(12*n_samB2*((dd_B2[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((dd_B2[2]/n_samB2-1/3)^2),180*n_samB2*((dd_B2[3]/n_samB2-1/12)^2)))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[2])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),D_B2))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("2MG-EAD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigmaP1,4),round(t(mean1B),4)," "," "," ",round(sigma1B[1],4),round(t(mix_pi_1),4)," "," "," ",
                       round(t(mean2B),4)," ",round(sigma2B[1],4),round(t(mix_pi_2),4)," ",
                       round(B1[1],4),round(B1[2],4),round(B1[2],4),round(B1[2],4),round(B1[2],4)," "," ",round(jj_1,4),round(ll_1*100,4)," "," ",round(jj_2,4),round(ll_2*100,4)," "," ",
                       round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],
                       round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],
                       round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],
                       round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],
                       round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
}
########################## (D1)#####################################
G5BCModelFun[[10]] <- function(K1,logL,df11,df21,df31,df41,df51){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]))
  n_samP1<-dim(dataP1)[1];n_samP2<-dim(dataP2)[1];n_samF1<-dim(dataF1)[1];n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2)
  sumx4<-sum(dataB1);sumx5<-sum(dataB2)
  ss1<- (n_samP1-1)*as.numeric(var(dataP1));meanP1<- mean(dataP1)
  ss2<- (n_samF1-1)*as.numeric(var(dataF1));meanF1<- mean(dataF1)
  ss3<- (n_samP2-1)*as.numeric(var(dataP2));meanP2<- mean(dataP2)
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  sigmaP1<- sigma0;sigmaF1<- sigma0;sigmaP2<- sigma0
  mean1<- mean(dataB1); sigma1<- as.numeric(var(dataB1))
  mean2<- mean(dataB2); sigma2<- as.numeric(var(dataB2))

  m_esp<-0.0001
  ####################### MX1-AD-AD Model ############ (D1)#####################
  d21<-2; d22<-2
  mi_1<- matrix(0.5,d21,1); a1<- sqrt(sigma1/(n_samB1-1))
  sigma1D<- matrix(sigma1/3,d21,1); mean1D<- as.matrix(c(mean1+a1,mean1-a1))

  mi_2<- matrix(0.5,d22,1); a2<- sqrt(sigma2/(n_samB2-1))
  sigma2D<- matrix(sigma2/3,d22,1); mean2D<- as.matrix(c(mean2+a2,mean2-a2))

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1D,sigma1D,dataB1)+logL(n_samB2,d22,mi_2,mean2D,sigma2D,dataB2)
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1D[i],sqrt(sigma1D[i]))/dmixnorm(dataB1,mean1D,sqrt(sigma1D),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2D[i],sqrt(sigma2D[i]))/dmixnorm(dataB2,mean2D,sqrt(sigma2D),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    ############ CM1-step for means ##############
    n01 <- n_samB1*mix_pi_1; n01[n01<0.000001] <- 0.000001
    n02 <- n_samB2*mix_pi_2; n02[n02<0.000001] <- 0.000001

    aa1<- sigmaP1/n_samP1+sigmaP1/n_samF1+sigma1D[1]/n01[1]+sigma1D[2]/n01[2]
    aa2<- sigmaP1/n_samF1
    aa3<- sigmaP1/n_samF1+sigmaP1/n_samP2+sigma2D[1]/n02[1]+sigma2D[2]/n02[2]
    aa4<- aa1*aa3-aa2*aa2
    s0<- matrix(0,2,1)
    s0[1]<- sumx1/n_samP1+sumx2/n_samF1-sumwx_B1[1]/n01[1]-sumwx_B1[2]/n01[2]
    s0[2]<- sumx2/n_samF1+sumx3/n_samP2-sumwx_B2[1]/n02[1]-sumwx_B2[2]/n02[2]
    aa5<- s0[1]*aa3-s0[2]*aa2
    aa6<- aa1*s0[2]-s0[1]*aa2
    rr<- matrix(0,2,1)
    rr[1]<- aa5/aa4;rr[2]<- aa6/aa4

    meanP1<- (sumx1-sigmaP1*rr[1])/n_samP1
    meanF1<- (sumx2-sigmaP1*(rr[1]+rr[2]))/n_samF1
    meanP2<- (sumx3-sigmaP1*rr[2])/n_samP2
    mean1D[1]<- (sumwx_B1[1]+sigma1D[1]*rr[1])/n01[1]; mean1D[2]<- (sumwx_B1[2]+sigma1D[1]*rr[1])/n01[2]
    mean2D[1]<- (sumwx_B2[1]+sigma2D[1]*rr[2])/n02[1]; mean2D[2]<- (sumwx_B2[2]+sigma2D[1]*rr[2])/n02[2]

    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2); ss2<- sum((dataF1-meanF1)^2); ss3<- sum((dataP2-meanP2)^2)
    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1D[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2D[i])^2 }

    sigma1D[1]<- sum(swx_B1)/n_samB1; sigma1D<- matrix(sigma1D[1],d21,1)
    sigma2D[1]<- sum(swx_B2)/n_samB1; sigma2D<- matrix(sigma2D[1],d22,1)
    #################  CM3 to estimate the variance #################
    aa1<- (sigma1D[1]-sigmaP1);aa2<- (sigma2D[1]-sigmaP1)
    aa1[aa1<0.0] <- 0.0; aa2[aa2<0.0] <- 0.0;
    n_iter<- 0.0; aaa0<- sigmaP1 ;aa5<- 1000
    while(aa5>0.0001){
      n_iter<- n_iter+1
      aa3<- sigmaP1/(sigmaP1+aa1); aa4<- sigmaP1/(sigmaP1+aa2)
      sigmaP1<- ((ss1+ss2+ss3)+aa3*aa3*(sum(swx_B1))+aa4*aa4*(sum(swx_B2)))/((n_samP1+n_samF1+n_samP2)+aa3*n_samB1+aa4*n_samB2)
      aa5<- abs(sigmaP1-aaa0)
      aaa0<- sigmaP1
      if (n_iter>20) break
    }
    sigmaF1<- sigmaP1;sigmaP2<- sigmaP1
    sigma1D<- matrix((sigmaP1+aa1),d21,1)
    sigma2D<- matrix((sigmaP1+aa2),d22,1)
    if(sum(sigma1D < 1e-30)>=1){break}
    if(sum(sigma2D < 1e-30)>=1){break}
    ####################### the stop criterion for iteration ############################################################
    L1<- logL(n_samP1,1,1,meanP1,sigmaP1,dataP1)+logL(n_samF1,1,1,meanF1,sigmaF1,dataF1)+logL(n_samP2,1,1,meanP2,sigmaP2,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1D,sigma1D,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2D,sigma2D,dataB2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*8

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1, 1,0,-1,1,0,0,-1, 0,1,0,0,1,1,0, 1,0,-1,0.5,0.5,-0.5,-0.5, 0,1,0,0.5,0.5,0.5,0.5),7,5)
  b_line1 <- matrix(c(meanP1,meanF1,meanP2,mean1D,mean2D))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigma1 - sigma1D[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigma1
  mm_1<- sigma1D[1]-sigmaP1
  if(mm_1 < 0) {mm_1 <- 0}
  nnn_1 <- mm_1/sigma1

  jj_2 <- sigma2 - sigma2D[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigma2
  mm_2<- sigma2D[1]-sigmaP1
  if(mm_2 < 0) {mm_2 <- 0}
  nnn_2 <- mm_2/sigma2

  ######### hypothesis testing #########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - meanP1)/sqrt(as.vector(sigmaP1))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[2])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),D_P1))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])

  ###############  F1  ###################
  dataF1<-sort(dataF1);bmw_F1 <- matrix(0,n_samF1,1); bmwsl_F1 <- matrix(0,n_samF1,1)

  gg_F1 <- (dataF1 - meanF1)/sqrt(as.vector(sigmaF1))
  bmw_F1[which(gg_F1>=0)] <- pnorm(gg_F1[gg_F1>=0])
  bmw_F1[which(gg_F1<0)] <- 1 - pnorm(abs(gg_F1[gg_F1<0]))
  bmwsl_F1[,1] <- bmw_F1

  P2_F1 <- rowSums(bmwsl_F1)
  nn<-dim(as.matrix(unique(P2_F1)))[1]
  if(nn<n_samF1){P2_F1<-P2_F1+runif(n_samF1)/1e4}

  dd_F1 <- as.matrix(c(sum(P2_F1),sum(P2_F1^2),sum((P2_F1-0.5)^2)))
  WW2_F1 <- 1/(12*n_samF1) + sum((P2_F1 - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  u_F1 <- as.matrix(c(12*n_samF1*((dd_F1[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((dd_F1[2]/n_samF1-1/3)^2),180*n_samF1*((dd_F1[3]/n_samF1-1/12)^2)))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[2])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),D_F1))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - meanP2)/sqrt(as.vector(sigmaP2))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[2])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),D_P2))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])

  ###############  B1  ###################
  dataB1<-sort(dataB1);bmw_B1 <- matrix(0,n_samB1,1); bmwsl_B1 <- matrix(0,n_samB1,d21)
  for(i in 1:d21){
    gg_B1 <- (dataB1 - mean1D[i])/sqrt(sigma1D[i])
    bmw_B1[which(gg_B1>=0)] <- pnorm(gg_B1[gg_B1>=0])
    bmw_B1[which(gg_B1<0)] <- 1 - pnorm(abs(gg_B1[gg_B1<0]))
    bmwsl_B1[,i] <- bmw_B1*mix_pi_1[i]
  }
  P2_B1 <- rowSums(bmwsl_B1)
  nn<-dim(as.matrix(unique(P2_B1)))[1]
  if(nn<n_samB1){P2_B1<-P2_B1+runif(n_samB1)/1e4}

  dd_B1 <- as.matrix(c(sum(P2_B1),sum(P2_B1^2),sum((P2_B1-0.5)^2)))
  WW2_B1 <- 1/(12*n_samB1) + sum((P2_B1 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  u_B1 <- as.matrix(c(12*n_samB1*((dd_B1[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((dd_B1[2]/n_samB1-1/3)^2),180*n_samB1*((dd_B1[3]/n_samB1-1/12)^2)))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[2])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),D_B1))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])

  ###############  B2  ###################
  dataB2<-sort(dataB2);bmw_B2 <- matrix(0,n_samB2,1); bmwsl_B2 <- matrix(0,n_samB2,d22)
  for(i in 1:d22){
    gg_B2 <- (dataB2 - mean2D[i])/sqrt(sigma2D[i])
    bmw_B2[which(gg_B2>=0)] <- pnorm(gg_B2[gg_B2>=0])
    bmw_B2[which(gg_B2<0)] <- 1 - pnorm(abs(gg_B2[gg_B2<0]))
    bmwsl_B2[,i] <- bmw_B2*mix_pi_2[i]
  }
  P2_B2 <- rowSums(bmwsl_B2)
  nn<-dim(as.matrix(unique(P2_B2)))[1]
  if(nn<n_samB2){P2_B2<-P2_B2+runif(n_samB2)/1e4}

  dd_B2 <- as.matrix(c(sum(P2_B2),sum(P2_B2^2),sum((P2_B2-0.5)^2)))
  WW2_B2 <- 1/(12*n_samB2) + sum((P2_B2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  u_B2 <- as.matrix(c(12*n_samB2*((dd_B2[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((dd_B2[2]/n_samB2-1/3)^2),180*n_samB2*((dd_B2[3]/n_samB2-1/12)^2)))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[2])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),D_B2))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX1-AD-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigmaP1,4),round(t(mean1D),4)," "," ",round(sigma1D[1],4),round(t(mix_pi_1),4)," "," ",
                       round(t(mean2D),4)," "," ",round(sigma2D[1],4),round(t(mix_pi_2),4)," "," ",
                       round(B1[1],4),round(B1[2],4)," ",round(B1[3],4)," ",round(B1[4],4),round(B1[5],4),round(jj_1,4),round(ll_1*100,4),round(mm_1,4),round(nnn_1*100,4) ,round(jj_2,4),round(ll_2*100,4),round(mm_2,4),round(nnn_2*100,4) ,
                       round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],
                       round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],
                       round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],
                       round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],
                       round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
}
########################### (D2)#####################################
G5BCModelFun[[11]] <- function(K1,logL,df11,df21,df31,df41,df51){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]))
  n_samP1<-dim(dataP1)[1];n_samP2<-dim(dataP2)[1];n_samF1<-dim(dataF1)[1];n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2)
  sumx4<-sum(dataB1);sumx5<-sum(dataB2)
  ss1<- (n_samP1-1)*as.numeric(var(dataP1));meanP1<- mean(dataP1)
  ss2<- (n_samF1-1)*as.numeric(var(dataF1));meanF1<- mean(dataF1)
  ss3<- (n_samP2-1)*as.numeric(var(dataP2));meanP2<- mean(dataP2)
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  sigmaP1<- sigma0;sigmaF1<- sigma0;sigmaP2<- sigma0
  mean1<- mean(dataB1); sigma1<- as.numeric(var(dataB1))
  mean2<- mean(dataB2); sigma2<- as.numeric(var(dataB2))

  m_esp<-0.0001
  ####################### MX1-A-AD Model ############ (D2)#####################
  d21<-2; d22<-2
  mi_1<- matrix(0.5,d21,1); a1<- sqrt(sigma1/(n_samB1-1))
  sigma1D<- matrix(sigma1/3,d21,1); mean1D<- as.matrix(c(mean1+a1,mean1-a1))

  mi_2<- matrix(0.5,d22,1); a2<- sqrt(sigma2/(n_samB2-1))
  sigma2D<- matrix(sigma2/3,d22,1); mean2D<- as.matrix(c(mean2+a2,mean2-a2))

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1D,sigma1D,dataB1)+logL(n_samB2,d22,mi_2,mean2D,sigma2D,dataB2)
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1D[i],sqrt(sigma1D[i]))/dmixnorm(dataB1,mean1D,sqrt(sigma1D),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2D[i],sqrt(sigma2D[i]))/dmixnorm(dataB2,mean2D,sqrt(sigma2D),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    ############ CM1-step for means ##############
    n01 <- n_samB1*mix_pi_1; n01[n01<0.000001] <- 0.000001
    n02 <- n_samB2*mix_pi_2; n02[n02<0.000001] <- 0.000001

    hh<- matrix(0,3,3)
    hh[1,1]<- sigmaP1/n_samP1+sigmaF1/n_samF1+sigma1D[1]/n01[1]+sigma1D[2]/n01[2]
    hh[1,2]<- sigmaF1/n_samF1
    hh[1,3]<- -sigma1D[1]/n01[1]+sigma1D[2]/n01[2]
    hh[2,2]<- sigmaF1/n_samF1+sigmaP2/n_samP2+sigma2D[1]/n02[1]+sigma2D[2]/n02[2]
    hh[2,3]<- sigma2D[1]/n02[1]-sigma2D[2]/n02[2]
    hh[3,3]<- sigma1D[1]/n01[1]+sigma1D[2]/n01[2]+sigma2D[1]/n02[1]+sigma2D[2]/n02[2]
    for(i in 2:3){
      for(j in 1:(i-1)){
        hh[i,j]<- hh[j,i]
      }
    }
    b_line<- matrix(0,3,1)
    b_line[1]<- sumx1/n_samP1+sumx2/n_samF1-sumwx_B1[1]/n01[1]-sumwx_B1[2]/n01[2]
    b_line[2]<- sumx2/n_samF1+sumx3/n_samP2-sumwx_B2[1]/n02[1]-sumwx_B2[2]/n02[2]
    b_line[3]<- sumwx_B1[1]/n01[1]-sumwx_B1[2]/n01[2]-sumwx_B2[1]/n02[1]+sumwx_B2[2]/n02[2]
    B <- solve(hh,b_line)

    meanP1<- (sumx1-sigmaP1*B[1])/n_samP1
    meanF1<- (sumx2-sigmaF1*(B[1]+B[2]))/n_samF1
    meanP2<- (sumx3-sigmaP2*B[2])/n_samP2
    mean1D[1]<- (sumwx_B1[1]+sigma1D[1]*(B[1]-B[3]))/n01[1]; mean1D[2]<- (sumwx_B1[2]+sigma1D[2]*(B[1]+B[3]))/n01[2]
    mean2D[1]<- (sumwx_B2[1]+sigma2D[1]*(B[2]+B[3]))/n02[1]; mean2D[2]<- (sumwx_B2[2]+sigma2D[2]*(B[2]-B[3]))/n02[2]

    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2); ss2<- sum((dataF1-meanF1)^2); ss3<- sum((dataP2-meanP2)^2)
    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1D[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2D[i])^2 }

    sigma1D[1]<- sum(swx_B1)/n_samB1; sigma1D<- matrix(sigma1D[1],d21,1)
    sigma2D[1]<- sum(swx_B2)/n_samB1; sigma2D<- matrix(sigma2D[1],d22,1)
    #################  CM3 to estimate the variance #################
    aa1<- sigma1D[1]-sigmaP1; aa2<- sigma2D[1]-sigmaP1
    aa1[aa1<0.0] <- 0.0; aa2[aa2<0.0] <- 0.0;
    n_iter<- 0.0; aaa0<- sigmaP1 ;aa5<- 1000
    while(aa5>0.0001){
      n_iter<- n_iter+1
      aa3<- sigmaP1/(sigmaP1+aa1); aa4<- sigmaP1/(sigmaP1+aa2)
      sigmaP1<- ((ss1+ss2+ss3)+aa3*aa3*(sum(swx_B1))+aa4*aa4*(sum(swx_B2)))/((n_samP1+n_samF1+n_samP2)+aa3*n_samB1+aa4*n_samB2)
      aa5<- abs(sigmaP1-aaa0)
      aaa0<- sigmaP1
      if (n_iter>20) break
    }
    sigmaF1<- sigmaP1;sigmaP2<- sigmaP1
    sigma1D<- matrix((sigmaP1+aa1),d21,1)
    sigma2D<- matrix((sigmaP1+aa2),d22,1)
    if(sum(sigma1D < 1e-30)>=1){break}
    if(sum(sigma2D < 1e-30)>=1){break}
    ####################### the stop criterion for iteration ############################################################
    L1<- logL(n_samP1,1,1,meanP1,sigmaP1,dataP1)+logL(n_samF1,1,1,meanF1,sigmaF1,dataF1)+logL(n_samP2,1,1,meanP2,sigmaP2,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1D,sigma1D,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2D,sigma2D,dataB2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*7

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1, 1,0,-1,1,0,0,-1, 1,0,-1,0.5,0.5,-0.5,-0.5, 0,1,0,0.5,0.5,0.5,0.5),7,4)
  b_line1 <- matrix(c(meanP1,meanF1,meanP2,mean1D,mean2D))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigma1 - sigma1D[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigma1
  mm_1<- sigma1D[1]-sigmaP1
  if(mm_1 < 0) {mm_1 <- 0}
  nnn_1 <- mm_1/sigma1

  jj_2 <- sigma2 - sigma2D[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigma2
  mm_2<- sigma2D[1]-sigmaP1
  if(mm_2 < 0) {mm_2 <- 0}
  nnn_2 <- mm_2/sigma2

  ######### hypothesis testing #########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - meanP1)/sqrt(as.vector(sigmaP1))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[2])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),D_P1))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])

  ###############  F1  ###################
  dataF1<-sort(dataF1);bmw_F1 <- matrix(0,n_samF1,1); bmwsl_F1 <- matrix(0,n_samF1,1)

  gg_F1 <- (dataF1 - meanF1)/sqrt(as.vector(sigmaF1))
  bmw_F1[which(gg_F1>=0)] <- pnorm(gg_F1[gg_F1>=0])
  bmw_F1[which(gg_F1<0)] <- 1 - pnorm(abs(gg_F1[gg_F1<0]))
  bmwsl_F1[,1] <- bmw_F1

  P2_F1 <- rowSums(bmwsl_F1)
  nn<-dim(as.matrix(unique(P2_F1)))[1]
  if(nn<n_samF1){P2_F1<-P2_F1+runif(n_samF1)/1e4}

  dd_F1 <- as.matrix(c(sum(P2_F1),sum(P2_F1^2),sum((P2_F1-0.5)^2)))
  WW2_F1 <- 1/(12*n_samF1) + sum((P2_F1 - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  u_F1 <- as.matrix(c(12*n_samF1*((dd_F1[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((dd_F1[2]/n_samF1-1/3)^2),180*n_samF1*((dd_F1[3]/n_samF1-1/12)^2)))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[2])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),D_F1))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - meanP2)/sqrt(as.vector(sigmaP2))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[2])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),D_P2))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])

  ###############  B1  ###################
  dataB1<-sort(dataB1);bmw_B1 <- matrix(0,n_samB1,1); bmwsl_B1 <- matrix(0,n_samB1,d21)
  for(i in 1:d21){
    gg_B1 <- (dataB1 - mean1D[i])/sqrt(sigma1D[i])
    bmw_B1[which(gg_B1>=0)] <- pnorm(gg_B1[gg_B1>=0])
    bmw_B1[which(gg_B1<0)] <- 1 - pnorm(abs(gg_B1[gg_B1<0]))
    bmwsl_B1[,i] <- bmw_B1*mix_pi_1[i]
  }
  P2_B1 <- rowSums(bmwsl_B1)
  nn<-dim(as.matrix(unique(P2_B1)))[1]
  if(nn<n_samB1){P2_B1<-P2_B1+runif(n_samB1)/1e4}

  dd_B1 <- as.matrix(c(sum(P2_B1),sum(P2_B1^2),sum((P2_B1-0.5)^2)))
  WW2_B1 <- 1/(12*n_samB1) + sum((P2_B1 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  u_B1 <- as.matrix(c(12*n_samB1*((dd_B1[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((dd_B1[2]/n_samB1-1/3)^2),180*n_samB1*((dd_B1[3]/n_samB1-1/12)^2)))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[2])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),D_B1))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])

  ###############  B2  ###################
  dataB2<-sort(dataB2);bmw_B2 <- matrix(0,n_samB2,1); bmwsl_B2 <- matrix(0,n_samB2,d22)
  for(i in 1:d22){
    gg_B2 <- (dataB2 - mean2D[i])/sqrt(sigma2D[i])
    bmw_B2[which(gg_B2>=0)] <- pnorm(gg_B2[gg_B2>=0])
    bmw_B2[which(gg_B2<0)] <- 1 - pnorm(abs(gg_B2[gg_B2<0]))
    bmwsl_B2[,i] <- bmw_B2*mix_pi_2[i]
  }
  P2_B2 <- rowSums(bmwsl_B2)
  nn<-dim(as.matrix(unique(P2_B2)))[1]
  if(nn<n_samB2){P2_B2<-P2_B2+runif(n_samB2)/1e4}

  dd_B2 <- as.matrix(c(sum(P2_B2),sum(P2_B2^2),sum((P2_B2-0.5)^2)))
  WW2_B2 <- 1/(12*n_samB2) + sum((P2_B2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  u_B2 <- as.matrix(c(12*n_samB2*((dd_B2[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((dd_B2[2]/n_samB2-1/3)^2),180*n_samB2*((dd_B2[3]/n_samB2-1/12)^2)))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[2])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),D_B2))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX1-A-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigmaP1,4),round(t(mean1D),4)," "," ",round(sigma1D[1],4),round(t(mix_pi_1),4)," "," ",
                       round(t(mean2D),4)," "," ",round(sigma2D[1],4),round(t(mix_pi_2),4)," "," ",
                       round(B1[1],4),round(B1[2],4)," "," "," ",round(B1[3],4),round(B1[4],4),round(jj_1,4),round(ll_1*100,4),round(mm_1,4),round(nnn_1*100,4) ,round(jj_2,4),round(ll_2*100,4),round(mm_2,4),round(nnn_2*100,4) ,
                       round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],
                       round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],
                       round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],
                       round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],
                       round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
}
########################### (D3)#####################################
G5BCModelFun[[12]] <- function(K1,logL,df11,df21,df31,df41,df51){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]))
  n_samP1<-dim(dataP1)[1];n_samP2<-dim(dataP2)[1];n_samF1<-dim(dataF1)[1];n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2)
  sumx4<-sum(dataB1);sumx5<-sum(dataB2)
  ss1<- (n_samP1-1)*as.numeric(var(dataP1));meanP1<- mean(dataP1)
  ss2<- (n_samF1-1)*as.numeric(var(dataF1));meanF1<- mean(dataF1)
  ss3<- (n_samP2-1)*as.numeric(var(dataP2));meanP2<- mean(dataP2)
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  sigmaP1<- sigma0;sigmaF1<- sigma0;sigmaP2<- sigma0
  mean1<- mean(dataB1); sigma1<- as.numeric(var(dataB1))
  mean2<- mean(dataB2); sigma2<- as.numeric(var(dataB2))

  m_esp<-0.0001
  ####################### MX1-EAD-AD Model ############ (D3)#####################
  d21<-1; d22<-2
  mi_1<- as.matrix(1)
  sigma1D<- as.matrix(sigma1/3) ; mean1D<- as.matrix(mean1)

  mi_2<- matrix(0.5,d22,1); a2<- sqrt(sigma2/(n_samB2-1))
  sigma2D<- matrix(sigma2/3,d22,1); mean2D<- as.matrix(c(mean2+a2,mean2-a2))

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1D,sigma1D,dataB1)+logL(n_samB2,d22,mi_2,mean2D,sigma2D,dataB2)
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1D[i],sqrt(sigma1D[i]))/dmixnorm(dataB1,mean1D,sqrt(sigma1D),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2D[i],sqrt(sigma2D[i]))/dmixnorm(dataB2,mean2D,sqrt(sigma2D),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    ############ CM1-step for means ##############
    n01 <- n_samB1*mix_pi_1; n01[n01<0.000001] <- 0.000001
    n02 <- n_samB2*mix_pi_2; n02[n02<0.000001] <- 0.000001

    aa1<- sigmaP1/n_samP1+sigmaP1/n_samF1+4.0*sigma1D[1]/n_samP1
    aa2<- sigmaF1/n_samF1
    aa3<- sigmaF1/n_samF1+sigmaP2/n_samP2+sigma2D[1]/n02[1]+sigma2D[2]/n02[2]
    aa4<- aa1*aa3-aa2*aa2
    s0<- matrix(0,2,1)
    s0[1]<- sumx1/n_samP1+sumx2/n_samF1-2.0*sumwx_B1[1]/n01[1]
    s0[2]<- sumx2/n_samF1+sumx3/n_samP2-sumwx_B2[1]/n02[1]-sumwx_B2[2]/n02[2]
    aa5<- s0[1]*aa3-s0[2]*aa2
    aa6<- aa1*s0[2]-s0[1]*aa2
    rr<- matrix(0,2,1)
    rr[1]<- aa5/aa4;rr[2]<- aa6/aa4

    meanP1<- (sumx1-sigmaP1*rr[1])/n_samP1
    meanF1<- (sumx2-sigmaF1*(rr[1]+rr[2]))/n_samF1
    meanP2<- (sumx3-sigmaP2*rr[2])/n_samP2
    mean1D[1]<- (sumwx_B1[1]+sigma1D[1]*rr[1]*2)/n01[1]
    mean2D[1]<- (sumwx_B2[1]+sigma2D[1]*rr[2])/n02[1]; mean2D[2]<- (sumwx_B2[2]+sigma2D[2]*rr[2])/n02[2]

    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2); ss2<- sum((dataF1-meanF1)^2); ss3<- sum((dataP2-meanP2)^2)
    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1D[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2D[i])^2 }

    sigma1D[1]<- sum(swx_B1)/n_samB1
    sigma2D[1]<- sum(swx_B2)/n_samB1; sigma2D<- matrix(sigma2D[1],d22,1)
    #################  CM3 to estimate the variance (sigma) #################
    aa1<- (sigma1D[1]-sigmaP1);aa2<- (sigma2D[1]-sigmaP1)
    aa1[aa1<0.0] <- 0.0; aa2[aa2<0.0] <- 0.0
    n_iter<- 0.0; aaa0<- sigmaP1 ;aa5<- 1000
    while(aa5>0.0001){
      n_iter<- n_iter+1
      aa3<- sigmaP1/(sigmaP1+aa1); aa4<- sigmaP1/(sigmaP1+aa2)
      sigmaP1<- ((ss1+ss2+ss3)+aa3*aa3*(sum(swx_B1))+aa4*aa4*(sum(swx_B2)))/((n_samP1+n_samF1+n_samP2)+aa3*n_samB1+aa4*n_samB2)
      aa5<- abs(sigmaP1-aaa0)
      aaa0<- sigmaP1
      if (n_iter>20) break
    }
    sigmaF1<- sigmaP1;sigmaP2<- sigmaP1
    sigma1D<- sigmaP1+aa1
    sigma2D<- matrix((sigmaP1+aa2),d22,1)
    if(sum(sigma1D < 1e-30)>=1){break}
    if(sum(sigma2D < 1e-30)>=1){break}
    ####################### the stop criterion for iteration ############################################################
    L1<- logL(n_samP1,1,1,meanP1,sigmaP1,dataP1)+logL(n_samF1,1,1,meanF1,sigmaF1,dataF1)+logL(n_samP2,1,1,meanP2,sigmaP2,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1D,sigma1D,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2D,sigma2D,dataB2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*7

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1, 1,1,-1,1,1,-1, 1,0,-1,0.5,-0.5,-0.5, 0,1,0,0.5,0.5,0.5),6,4)
  b_line1 <- matrix(c(meanP1,meanF1,meanP2,mean1D,mean2D))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigma1 - sigma1D[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigma1
  mm_1<- sigma1D[1]-sigmaP1
  if(mm_1 < 0) {mm_1 <- 0}
  nnn_1 <- mm_1/sigma1

  jj_2 <- sigma2 - sigma2D[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigma2
  mm_2<- sigma2D[1]-sigmaP1
  if(mm_2 < 0) {mm_2 <- 0}
  nnn_2 <- mm_2/sigma2

  ######### hypothesis testing #########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - meanP1)/sqrt(as.vector(sigmaP1))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[2])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),D_P1))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])

  ###############  F1  ###################
  dataF1<-sort(dataF1);bmw_F1 <- matrix(0,n_samF1,1); bmwsl_F1 <- matrix(0,n_samF1,1)

  gg_F1 <- (dataF1 - meanF1)/sqrt(as.vector(sigmaF1))
  bmw_F1[which(gg_F1>=0)] <- pnorm(gg_F1[gg_F1>=0])
  bmw_F1[which(gg_F1<0)] <- 1 - pnorm(abs(gg_F1[gg_F1<0]))
  bmwsl_F1[,1] <- bmw_F1

  P2_F1 <- rowSums(bmwsl_F1)
  nn<-dim(as.matrix(unique(P2_F1)))[1]
  if(nn<n_samF1){P2_F1<-P2_F1+runif(n_samF1)/1e4}

  dd_F1 <- as.matrix(c(sum(P2_F1),sum(P2_F1^2),sum((P2_F1-0.5)^2)))
  WW2_F1 <- 1/(12*n_samF1) + sum((P2_F1 - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  u_F1 <- as.matrix(c(12*n_samF1*((dd_F1[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((dd_F1[2]/n_samF1-1/3)^2),180*n_samF1*((dd_F1[3]/n_samF1-1/12)^2)))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[2])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),D_F1))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - meanP2)/sqrt(as.vector(sigmaP2))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[2])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),D_P2))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])

  ###############  B1  ###################
  dataB1<-sort(dataB1);bmw_B1 <- matrix(0,n_samB1,1); bmwsl_B1 <- matrix(0,n_samB1,d21)
  for(i in 1:d21){
    gg_B1 <- (dataB1 - mean1D[i])/sqrt(as.vector(sigma1D[i]))
    bmw_B1[which(gg_B1>=0)] <- pnorm(gg_B1[gg_B1>=0])
    bmw_B1[which(gg_B1<0)] <- 1 - pnorm(abs(gg_B1[gg_B1<0]))
    bmwsl_B1[,i] <- bmw_B1*mix_pi_1[i]
  }
  P2_B1 <- rowSums(bmwsl_B1)
  nn<-dim(as.matrix(unique(P2_B1)))[1]
  if(nn<n_samB1){P2_B1<-P2_B1+runif(n_samB1)/1e4}

  dd_B1 <- as.matrix(c(sum(P2_B1),sum(P2_B1^2),sum((P2_B1-0.5)^2)))
  WW2_B1 <- 1/(12*n_samB1) + sum((P2_B1 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  u_B1 <- as.matrix(c(12*n_samB1*((dd_B1[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((dd_B1[2]/n_samB1-1/3)^2),180*n_samB1*((dd_B1[3]/n_samB1-1/12)^2)))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[2])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),D_B1))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])

  ###############  B2  ###################
  dataB2<-sort(dataB2);bmw_B2 <- matrix(0,n_samB2,1); bmwsl_B2 <- matrix(0,n_samB2,d22)
  for(i in 1:d22){
    gg_B2 <- (dataB2 - mean2D[i])/sqrt(sigma2D[i])
    bmw_B2[which(gg_B2>=0)] <- pnorm(gg_B2[gg_B2>=0])
    bmw_B2[which(gg_B2<0)] <- 1 - pnorm(abs(gg_B2[gg_B2<0]))
    bmwsl_B2[,i] <- bmw_B2*mix_pi_2[i]
  }
  P2_B2 <- rowSums(bmwsl_B2)
  nn<-dim(as.matrix(unique(P2_B2)))[1]
  if(nn<n_samB2){P2_B2<-P2_B2+runif(n_samB2)/1e4}

  dd_B2 <- as.matrix(c(sum(P2_B2),sum(P2_B2^2),sum((P2_B2-0.5)^2)))
  WW2_B2 <- 1/(12*n_samB2) + sum((P2_B2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  u_B2 <- as.matrix(c(12*n_samB2*((dd_B2[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((dd_B2[2]/n_samB2-1/3)^2),180*n_samB2*((dd_B2[3]/n_samB2-1/12)^2)))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[2])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),D_B2))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX1-EAD-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigmaP1,4),round(t(mean1D),4)," "," "," ",round(sigma1D[1],4),round(t(mix_pi_1),4)," "," "," ",
                       round(t(mean2D),4)," "," ",round(sigma2D[1],4),round(t(mix_pi_2),4)," "," ",
                       round(B1[1],4),round(B1[2],4)," ",round(B1[2],4)," ",round(B1[3],4),round(B1[4],4),round(jj_1,4),round(ll_1*100,4),round(mm_1,4),round(nnn_1*100,4) ,round(jj_2,4),round(ll_2*100,4),round(mm_2,4),round(nnn_2*100,4) ,
                       round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],
                       round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],
                       round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],
                       round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],
                       round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
}
########################### (D4)#####################################
G5BCModelFun[[13]] <- function(K1,logL,df11,df21,df31,df41,df51){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]))
  n_samP1<-dim(dataP1)[1];n_samP2<-dim(dataP2)[1];n_samF1<-dim(dataF1)[1];n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2)
  sumx4<-sum(dataB1);sumx5<-sum(dataB2)
  ss1<- (n_samP1-1)*as.numeric(var(dataP1));meanP1<- mean(dataP1)
  ss2<- (n_samF1-1)*as.numeric(var(dataF1));meanF1<- mean(dataF1)
  ss3<- (n_samP2-1)*as.numeric(var(dataP2));meanP2<- mean(dataP2)
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  sigmaP1<- sigma0;sigmaF1<- sigma0;sigmaP2<- sigma0
  mean1<- mean(dataB1); sigma1<- as.numeric(var(dataB1))
  mean2<- mean(dataB2); sigma2<- as.numeric(var(dataB2))

  m_esp<-0.0001
  ####################### MX1-NCD-AD Model ############ (D4)#####################
  d21<-2; d22<-1
  mi_1<- matrix(0.5,d21,1); a1<- sqrt(sigma1/(n_samB1-1))
  sigma1D<- matrix(sigma1/3,d21,1); mean1D<- as.matrix(c(mean1+a1,mean1-a1))

  mi_2<- as.matrix(1)
  sigma2D<- as.matrix(sigma2/3) ; mean2D<- as.matrix(mean2)

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1D,sigma1D,dataB1)+logL(n_samB2,d22,mi_2,mean2D,sigma2D,dataB2)
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1D[i],sqrt(sigma1D[i]))/dmixnorm(dataB1,mean1D,sqrt(sigma1D),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2D[i],sqrt(sigma2D[i]))/dmixnorm(dataB2,mean2D,sqrt(sigma2D),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    ############ CM1-step for means ##############
    n01 <- n_samB1*mix_pi_1; n01[n01<0.000001] <- 0.000001
    n02 <- n_samB2*mix_pi_2; n02[n02<0.000001] <- 0.000001

    aa1<- sigmaP1/n_samP1+sigmaP1/n_samF1+sigma1D[1]/n01[1]+sigma1D[2]/n01[2]
    aa2<- sigmaF1/n_samF1
    aa3<- sigmaF1/n_samF1+sigmaP2/n_samP2+4.0*sigma2D[1]/n02[1]
    aa4<- aa1*aa3-aa2*aa2
    s0<- matrix(0,2,1)
    s0[1]<- sumx1/n_samP1+sumx2/n_samF1-sumwx_B1[1]/n01[1]-sumwx_B1[2]/n01[2]
    s0[2]<- sumx2/n_samF1+sumx3/n_samP2-2.0*sumwx_B2[1]/n02[1]
    aa5<- s0[1]*aa3-s0[2]*aa2
    aa6<- aa1*s0[2]-s0[1]*aa2
    rr<- matrix(0,2,1)
    rr[1]<- aa5/aa4; rr[2]<- aa6/aa4

    meanP1<- (sumx1-sigmaP1*rr[1])/n_samP1
    meanF1<- (sumx2-sigmaF1*(rr[1]+rr[2]))/n_samF1
    meanP2<- (sumx3-sigmaP2*rr[2])/n_samP2
    mean1D[1]<- (sumwx_B1[1]+sigma1D[1]*rr[1])/n01[1]; mean1D[2]<- (sumwx_B1[2]+sigma1D[2]*rr[1])/n01[2]
    mean2D[1]<- (sumwx_B2[1]+sigma2D[1]*rr[2]*2.0)/n02[1]

    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2); ss2<- sum((dataF1-meanF1)^2); ss3<- sum((dataP2-meanP2)^2)
    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1D[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2D[i])^2 }

    sigma1D[1]<- sum(swx_B1)/n_samB1; sigma1D<- matrix(sigma1D[1],d21,1)
    sigma2D<- sum(swx_B2)/n_samB1
    #################  CM3 to estimate the variance #################
    aa1<- sigma1D[1]-sigmaP1; aa2<- sigma2D[1]-sigmaP1
    aa1[aa1<0.0] <- 0.0; aa2[aa2<0.0] <- 0.0;
    n_iter<- 0.0; aaa0<- sigmaP1 ;aa5<- 1000
    while(aa5>0.0001){
      n_iter<- n_iter+1
      aa3<- sigmaP1/(sigmaP1+aa1); aa4<- sigmaP1/(sigmaP1+aa2)
      sigmaP1<- ((ss1+ss2+ss3)+aa3*aa3*(sum(swx_B1))+aa4*aa4*(sum(swx_B2)))/((n_samP1+n_samF1+n_samP2)+aa3*n_samB1+aa4*n_samB2)
      aa5<- abs(sigmaP1-aaa0)
      aaa0<- sigmaP1
      if (n_iter>20) break
    }
    sigmaF1<- sigmaP1;sigmaP2<- sigmaP1
    sigma1D<- matrix((sigmaP1+aa1),d21,1)
    sigma2D<- matrix((sigmaP1+aa2),d22,1)
    if(sum(sigma1D < 1e-30)>=1){break}
    if(sum(sigma2D < 1e-30)>=1){break}
    ####################### the stop criterion for iteration ############################################################
    L1<- logL(n_samP1,1,1,meanP1,sigmaP1,dataP1)+logL(n_samF1,1,1,meanF1,sigmaF1,dataF1)+logL(n_samP2,1,1,meanP2,sigmaP2,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1D,sigma1D,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2D,sigma2D,dataB2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*7

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1, 1,-1,-1,1,-1,-1, 1,0,-1,0.5,0.5,-0.5, 0,1,0,0.5,0.5,0.5),6,4)
  b_line1 <- matrix(c(meanP1,meanF1,meanP2,mean1D,mean2D))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigma1 - sigma1D[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigma1
  mm_1<- sigma1D[1]-sigmaP1
  if(mm_1 < 0) {mm_1 <- 0}
  nnn_1 <- mm_1/sigma1

  jj_2 <- sigma2 - sigma2D[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigma2
  mm_2<- sigma2D[1]-sigmaP1
  if(mm_2 < 0) {mm_2 <- 0}
  nnn_2 <- mm_2/sigma2

  ######### hypothesis testing #########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - meanP1)/sqrt(as.vector(sigmaP1))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[2])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),D_P1))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])

  ###############  F1  ###################
  dataF1<-sort(dataF1);bmw_F1 <- matrix(0,n_samF1,1); bmwsl_F1 <- matrix(0,n_samF1,1)

  gg_F1 <- (dataF1 - meanF1)/sqrt(as.vector(sigmaF1))
  bmw_F1[which(gg_F1>=0)] <- pnorm(gg_F1[gg_F1>=0])
  bmw_F1[which(gg_F1<0)] <- 1 - pnorm(abs(gg_F1[gg_F1<0]))
  bmwsl_F1[,1] <- bmw_F1

  P2_F1 <- rowSums(bmwsl_F1)
  nn<-dim(as.matrix(unique(P2_F1)))[1]
  if(nn<n_samF1){P2_F1<-P2_F1+runif(n_samF1)/1e4}

  dd_F1 <- as.matrix(c(sum(P2_F1),sum(P2_F1^2),sum((P2_F1-0.5)^2)))
  WW2_F1 <- 1/(12*n_samF1) + sum((P2_F1 - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  u_F1 <- as.matrix(c(12*n_samF1*((dd_F1[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((dd_F1[2]/n_samF1-1/3)^2),180*n_samF1*((dd_F1[3]/n_samF1-1/12)^2)))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[2])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),D_F1))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - meanP2)/sqrt(as.vector(sigmaP2))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[2])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),D_P2))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])

  ###############  B1  ###################
  dataB1<-sort(dataB1);bmw_B1 <- matrix(0,n_samB1,1); bmwsl_B1 <- matrix(0,n_samB1,d21)
  for(i in 1:d21){
    gg_B1 <- (dataB1 - mean1D[i])/sqrt(sigma1D[i])
    bmw_B1[which(gg_B1>=0)] <- pnorm(gg_B1[gg_B1>=0])
    bmw_B1[which(gg_B1<0)] <- 1 - pnorm(abs(gg_B1[gg_B1<0]))
    bmwsl_B1[,i] <- bmw_B1*mix_pi_1[i]
  }
  P2_B1 <- rowSums(bmwsl_B1)
  nn<-dim(as.matrix(unique(P2_B1)))[1]
  if(nn<n_samB1){P2_B1<-P2_B1+runif(n_samB1)/1e4}

  dd_B1 <- as.matrix(c(sum(P2_B1),sum(P2_B1^2),sum((P2_B1-0.5)^2)))
  WW2_B1 <- 1/(12*n_samB1) + sum((P2_B1 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  u_B1 <- as.matrix(c(12*n_samB1*((dd_B1[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((dd_B1[2]/n_samB1-1/3)^2),180*n_samB1*((dd_B1[3]/n_samB1-1/12)^2)))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[2])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),D_B1))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])

  ###############  B2  ###################
  dataB2<-sort(dataB2);bmw_B2 <- matrix(0,n_samB2,1); bmwsl_B2 <- matrix(0,n_samB2,d22)
  for(i in 1:d22){
    gg_B2 <- (dataB2 - mean2D[i])/sqrt(as.vector(sigma2D[i]))
    bmw_B2[which(gg_B2>=0)] <- pnorm(gg_B2[gg_B2>=0])
    bmw_B2[which(gg_B2<0)] <- 1 - pnorm(abs(gg_B2[gg_B2<0]))
    bmwsl_B2[,i] <- bmw_B2*mix_pi_2[i]
  }
  P2_B2 <- rowSums(bmwsl_B2)
  nn<-dim(as.matrix(unique(P2_B2)))[1]
  if(nn<n_samB2){P2_B2<-P2_B2+runif(n_samB2)/1e4}

  dd_B2 <- as.matrix(c(sum(P2_B2),sum(P2_B2^2),sum((P2_B2-0.5)^2)))
  WW2_B2 <- 1/(12*n_samB2) + sum((P2_B2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  u_B2 <- as.matrix(c(12*n_samB2*((dd_B2[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((dd_B2[2]/n_samB2-1/3)^2),180*n_samB2*((dd_B2[3]/n_samB2-1/12)^2)))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[2])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),D_B2))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX1-NCD-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigmaP1,4),round(t(mean1D),4)," "," ",round(sigma1D[1],4),round(t(mix_pi_1),4)," "," ",
                       round(t(mean2D),4)," "," "," ",round(sigma2D[1],4),round(t(mix_pi_2),4)," "," "," ",
                       round(B1[1],4),round(B1[2],4)," ",round(-B1[2],4)," ",round(B1[3],4),round(B1[4],4),round(jj_1,4),round(ll_1*100,4),round(mm_1,4),round(nnn_1*100,4) ,round(jj_2,4),round(ll_2*100,4),round(mm_2,4),round(nnn_2*100,4) ,
                       round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],
                       round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],
                       round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],
                       round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],
                       round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
}

########################### (E1)#####################################
G5BCModelFun[[14]] <- function(K1,logL,df11,df21,df31,df41,df51){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]))
  n_samP1<-dim(dataP1)[1];n_samP2<-dim(dataP2)[1];n_samF1<-dim(dataF1)[1];n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2)
  sumx4<-sum(dataB1);sumx5<-sum(dataB2)
  ss1<- (n_samP1-1)*as.numeric(var(dataP1));meanP1<- mean(dataP1)
  ss2<- (n_samF1-1)*as.numeric(var(dataF1));meanF1<- mean(dataF1)
  ss3<- (n_samP2-1)*as.numeric(var(dataP2));meanP2<- mean(dataP2)
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  sigmaP1<- sigma0;sigmaF1<- sigma0;sigmaP2<- sigma0
  mean1<- mean(dataB1); sigma1<- as.numeric(var(dataB1))
  mean2<- mean(dataB2); sigma2<- as.numeric(var(dataB2))

  m_esp<-0.0001
  ####################### MX2-AD-AD Model ############ (E1)#####################
  d21<-4; d22<-4
  mi_1<- matrix(0.25,d21,1); a1<- sqrt(sigma1/(n_samB1-1))
  sigma1E<- matrix(sigma1/3,d21,1); mean1E<- as.matrix(c(mean1+2*a1,mean1+0.8*a1,mean1-0.8*a1,mean1-2*a1))

  mi_2<- matrix(0.25,d22,1); a2<- sqrt(sigma2/(n_samB2-1))
  sigma2E<- matrix(sigma2/3,d22,1); mean2E<- as.matrix(c(mean2+2*a2,mean2+0.8*a2,mean2-0.8*a2,mean2-2*a2))

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1E,sigma1E,dataB1)+logL(n_samB2,d22,mi_2,mean2E,sigma2E,dataB2)
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1E[i],sqrt(sigma1E[i]))/dmixnorm(dataB1,mean1E,sqrt(sigma1E),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2E[i],sqrt(sigma2E[i]))/dmixnorm(dataB2,mean2E,sqrt(sigma2E),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    ############ CM1-step for means ##############
    n01 <- n_samB1*mix_pi_1; n01[n01<0.000001] <- 0.000001
    n02 <- n_samB2*mix_pi_2; n02[n02<0.000001] <- 0.000001

    hh<- matrix(0,4,4)
    hh[1,1]<- sigmaP1/n_samP1+sigmaP1/n_samF1+sigma1E[1]/n01[1]+sigma1E[4]/n01[4]
    hh[1,2]<- sigmaP1/n_samF1
    hh[1,3]<- -sigma1E[1]/n01[1]-sigma1E[4]/n01[4]
    hh[1,4]<- 0
    hh[2,2]<- sigmaP1/n_samF1+sigmaP1/n_samP2+sigma2E[1]/n02[1]+sigma2E[4]/n02[4]
    hh[2,3]<- 0
    hh[2,4]<- -sigma2E[1]/n02[1]-sigma2E[4]/n02[4]
    hh[3,3]<- sigma1E[1]/n01[1]+sigma1E[2]/n01[2]+sigma1E[3]/n01[3]+sigma1E[4]/n01[4]
    hh[3,4]<- 0
    hh[4,4]<- sigma2E[1]/n02[1]+sigma2E[2]/n02[2]+sigma2E[3]/n02[3]+sigma2E[4]/n02[4]
    for(i in 2:4){
      for(j in 1:(i-1)){
        hh[i,j]<- hh[j,i]
      }
    }
    b_line<- matrix(0,4,1)
    b_line[1]<- sumx1/n_samP1+sumx2/n_samF1-sumwx_B1[1]/n01[1]-sumwx_B1[4]/n01[4]
    b_line[2]<- sumx2/n_samF1+sumx3/n_samP2-sumwx_B2[1]/n02[1]-sumwx_B2[4]/n02[4]
    b_line[3]<- sumwx_B1[1]/n01[1]-sumwx_B1[2]/n01[2]-sumwx_B1[3]/n01[3]+sumwx_B1[4]/n01[4]
    b_line[4]<- sumwx_B2[1]/n02[1]-sumwx_B2[2]/n02[2]-sumwx_B2[3]/n02[3]+sumwx_B2[4]/n02[4]
    B <- solve(hh,b_line)

    meanP1<- (sumx1-sigmaP1*B[1])/n_samP1
    meanF1<- (sumx2-sigmaP1*(B[1]+B[2]))/n_samF1
    meanP2<- (sumx3-sigmaP1*B[2])/n_samP2
    mean1E[1]<- (sumwx_B1[1]+sigma1E[1]*(B[1]-B[3]))/n01[1]
    mean1E[2]<- (sumwx_B1[2]+sigma1E[2]*B[3])/n01[2]
    mean1E[3]<- (sumwx_B1[3]+sigma1E[3]*B[3])/n01[3]
    mean1E[4]<- (sumwx_B1[4]+sigma1E[4]*(B[1]-B[3]))/n01[4]
    mean2E[1]<- (sumwx_B2[1]+sigma2E[1]*(B[2]-B[4]))/n02[1]
    mean2E[2]<- (sumwx_B2[2]+sigma2E[2]*B[4])/n02[2]
    mean2E[3]<- (sumwx_B2[3]+sigma2E[3]*B[4])/n02[3]
    mean2E[4]<- (sumwx_B2[4]+sigma2E[4]*(B[2]-B[4]))/n02[4]

    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2); ss2<- sum((dataF1-meanF1)^2); ss3<- sum((dataP2-meanP2)^2)
    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1E[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2E[i])^2 }

    s0<- matrix(0,4,1)
    s0[1]<- sum(swx_B1); sigma1E<- as.matrix((s0[1]/n_samB1),d21,1)
    s0[2]<- sum(swx_B2); sigma2E<- as.matrix((s0[2]/n_samB2),d22,1)
    aa1<- sigma1E[1]-sigmaP1;aa2<- sigma2E[1]-sigmaP1
    aa1[aa1<0] <- 0; aa2[aa2<0] <- 0
    aaa0<- sigmaP1; n_iter<- 0;aa5<- 1000
    s0[3]<- ss1+ss2+ss3; s0[4]<- n_samP1+n_samF1+n_samP2
    while(aa5>0.0001){
      n_iter<- n_iter+1
      aa3<- sigmaP1/(aa1+sigmaP1); aa4<- sigmaP1/(aa2+sigmaP1)
      sigmaP1<- (s0[3]+aa3*aa3*s0[1]+aa4*aa4*s0[2])/(s0[4]+aa3*n_samB1+aa4*n_samB2)
      aa5<- abs(sigmaP1-aaa0); aaa0<- sigmaP1
      if (n_iter>20) break
    }
    sigmaF1<- sigmaP1;sigmaP2<- sigmaP1
    sigma1E<- matrix((sigmaP1+aa1),4,1)
    sigma2E<- matrix((sigmaP1+aa2),4,1)
    if(sum(sigma1E < 1e-30)>=1){break}
    if(sum(sigma2E < 1e-30)>=1){break}
    ####################### the stop criterion for iteration ############################################################
    L1<- logL(n_samP1,1,1,meanP1,sigmaP1,dataP1)+logL(n_samF1,1,1,meanF1,sigmaF1,dataF1)+logL(n_samP2,1,1,meanP2,sigmaP2,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1E,sigma1E,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2E,sigma2E,dataB2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*8

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1,1,1,1,1, 1,0,-1,1,1,0,0,0,0,-1,-1, 1,0,-1,1,0,1,0,0,-1,0,-1, 0,1,0,0,0,1,1,1,1,0,0,
                0,1,0,0,1,0,1,1,0,1,0, 1,0,-1,0.5,0.5,0.5,0.5,-0.5,-0.5,-0.5,-0.5, 0,1,0,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5),11,7)
  b_line1 <- matrix(c(meanP1,meanF1,meanP2,mean1E,mean2E))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigma1 - sigma1E[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigma1
  mm_1<- sigma1E[1]-sigmaP1
  if(mm_1 < 0) {mm_1 <- 0}
  nnn_1 <- mm_1/sigma1

  jj_2 <- sigma2 - sigma2E[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigma2
  mm_2<- sigma2E[1]-sigmaP1
  if(mm_2 < 0) {mm_2 <- 0}
  nnn_2 <- mm_2/sigma2

  ######### hypothesis testing #########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - meanP1)/sqrt(as.vector(sigmaP1))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[2])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),D_P1))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])

  ###############  F1  ###################
  dataF1<-sort(dataF1);bmw_F1 <- matrix(0,n_samF1,1); bmwsl_F1 <- matrix(0,n_samF1,1)

  gg_F1 <- (dataF1 - meanF1)/sqrt(as.vector(sigmaF1))
  bmw_F1[which(gg_F1>=0)] <- pnorm(gg_F1[gg_F1>=0])
  bmw_F1[which(gg_F1<0)] <- 1 - pnorm(abs(gg_F1[gg_F1<0]))
  bmwsl_F1[,1] <- bmw_F1

  P2_F1 <- rowSums(bmwsl_F1)
  nn<-dim(as.matrix(unique(P2_F1)))[1]
  if(nn<n_samF1){P2_F1<-P2_F1+runif(n_samF1)/1e4}

  dd_F1 <- as.matrix(c(sum(P2_F1),sum(P2_F1^2),sum((P2_F1-0.5)^2)))
  WW2_F1 <- 1/(12*n_samF1) + sum((P2_F1 - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  u_F1 <- as.matrix(c(12*n_samF1*((dd_F1[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((dd_F1[2]/n_samF1-1/3)^2),180*n_samF1*((dd_F1[3]/n_samF1-1/12)^2)))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[2])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),D_F1))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - meanP2)/sqrt(as.vector(sigmaP2))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[2])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),D_P2))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])

  ###############  B1  ###################
  dataB1<-sort(dataB1);bmw_B1 <- matrix(0,n_samB1,1); bmwsl_B1 <- matrix(0,n_samB1,d21)
  for(i in 1:d21){
    gg_B1 <- (dataB1 - mean1E[i])/sqrt(sigma1E[i])
    bmw_B1[which(gg_B1>=0)] <- pnorm(gg_B1[gg_B1>=0])
    bmw_B1[which(gg_B1<0)] <- 1 - pnorm(abs(gg_B1[gg_B1<0]))
    bmwsl_B1[,i] <- bmw_B1*mix_pi_1[i]
  }
  P2_B1 <- rowSums(bmwsl_B1)
  nn<-dim(as.matrix(unique(P2_B1)))[1]
  if(nn<n_samB1){P2_B1<-P2_B1+runif(n_samB1)/1e4}

  dd_B1 <- as.matrix(c(sum(P2_B1),sum(P2_B1^2),sum((P2_B1-0.5)^2)))
  WW2_B1 <- 1/(12*n_samB1) + sum((P2_B1 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  u_B1 <- as.matrix(c(12*n_samB1*((dd_B1[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((dd_B1[2]/n_samB1-1/3)^2),180*n_samB1*((dd_B1[3]/n_samB1-1/12)^2)))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[2])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),D_B1))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])

  ###############  B2  ###################
  dataB2<-sort(dataB2);bmw_B2 <- matrix(0,n_samB2,1); bmwsl_B2 <- matrix(0,n_samB2,d22)
  for(i in 1:d22){
    gg_B2 <- (dataB2 - mean2E[i])/sqrt(sigma2E[i])
    bmw_B2[which(gg_B2>=0)] <- pnorm(gg_B2[gg_B2>=0])
    bmw_B2[which(gg_B2<0)] <- 1 - pnorm(abs(gg_B2[gg_B2<0]))
    bmwsl_B2[,i] <- bmw_B2*mix_pi_2[i]
  }
  P2_B2 <- rowSums(bmwsl_B2)
  nn<-dim(as.matrix(unique(P2_B2)))[1]
  if(nn<n_samB2){P2_B2<-P2_B2+runif(n_samB2)/1e4}

  dd_B2 <- as.matrix(c(sum(P2_B2),sum(P2_B2^2),sum((P2_B2-0.5)^2)))
  WW2_B2 <- 1/(12*n_samB2) + sum((P2_B2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  u_B2 <- as.matrix(c(12*n_samB2*((dd_B2[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((dd_B2[2]/n_samB2-1/3)^2),180*n_samB2*((dd_B2[3]/n_samB2-1/12)^2)))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[2])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),D_B2))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX2-AD-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigmaP1,4),round(t(mean1E),4),round(sigma1E[1],4),round(t(mix_pi_1),4),
                       round(t(mean2E),4),round(sigma2E[1],4),round(t(mix_pi_2),4),
                       round(B1[1],4),round(B1[2],4),round(B1[3],4),round(B1[4],4),round(B1[5],4),round(B1[6],4),round(B1[7],4),round(jj_1,4),round(ll_1*100,4),round(mm_1,4),round(nnn_1*100,4) ,round(jj_2,4),round(ll_2*100,4),round(mm_2,4),round(nnn_2*100,4) ,
                       round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],
                       round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],
                       round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],
                       round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],
                       round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
}
########################### (E2)#####################################
G5BCModelFun[[15]] <- function(K1,logL,df11,df21,df31,df41,df51){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]))
  n_samP1<-dim(dataP1)[1];n_samP2<-dim(dataP2)[1];n_samF1<-dim(dataF1)[1];n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2)
  sumx4<-sum(dataB1);sumx5<-sum(dataB2)
  ss1<- (n_samP1-1)*as.numeric(var(dataP1));meanP1<- mean(dataP1)
  ss2<- (n_samF1-1)*as.numeric(var(dataF1));meanF1<- mean(dataF1)
  ss3<- (n_samP2-1)*as.numeric(var(dataP2));meanP2<- mean(dataP2)
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  sigmaP1<- sigma0;sigmaF1<- sigma0;sigmaP2<- sigma0
  mean1<- mean(dataB1); sigma1<- as.numeric(var(dataB1))
  mean2<- mean(dataB2); sigma2<- as.numeric(var(dataB2))

  m_esp<-0.0001
  ####################### MX2-A-AD Model ############ (E2)#####################
  d21<-4; d22<-4
  mi_1<- matrix(0.25,d21,1); a1<- sqrt(sigma1/(n_samB1-1))
  sigma1E<- matrix(sigma1/3,d21,1); mean1E<- as.matrix(c(mean1+2*a1,mean1+0.8*a1,mean1-0.8*a1,mean1-2*a1))

  mi_2<- matrix(0.25,d22,1); a2<- sqrt(sigma2/(n_samB2-1))
  sigma2E<- matrix(sigma2/3,d22,1); mean2E<- as.matrix(c(mean2+2*a2,mean2+0.8*a2,mean2-0.8*a2,mean2-2*a2))

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1E,sigma1E,dataB1)+logL(n_samB2,d22,mi_2,mean2E,sigma2E,dataB2)
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1E[i],sqrt(sigma1E[i]))/dmixnorm(dataB1,mean1E,sqrt(sigma1E),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2E[i],sqrt(sigma2E[i]))/dmixnorm(dataB2,mean2E,sqrt(sigma2E),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    ############ CM1-step for means ##############
    n01 <- n_samB1*mix_pi_1; n01[n01<0.000001] <- 0.000001
    n02 <- n_samB2*mix_pi_2; n02[n02<0.000001] <- 0.000001

    hh<- matrix(0,6,6)
    hh[1,1]<- sigmaP1/n_samP1+sigmaP1/n_samF1+sigma1E[1]/n01[1]+sigma1E[4]/n01[4]
    hh[1,2]<- sigmaP1/n_samF1
    hh[1,3]<- -sigma1E[1]/n01[1]-sigma1E[4]/n01[4]
    hh[1,4]<- 0
    hh[1,5]<- 0
    hh[1,6]<- -sigma1E[1]/n01[1]+sigma1E[4]/n01[4]
    hh[2,2]<- sigmaP1/n_samF1+sigmaP1/n_samP2+sigma2E[1]/n02[1]+sigma2E[4]/n02[4]
    hh[2,3]<- 0
    hh[2,4]<- -sigma2E[1]/n02[1]-sigma2E[4]/n02[4]
    hh[2,5]<- 0
    hh[2,6]<- sigma2E[1]/n02[1]-sigma2E[4]/n02[4]
    hh[3,3]<- sigma1E[1]/n01[1]+sigma1E[2]/n01[2]+sigma1E[3]/n01[3]+sigma1E[4]/n01[4]
    hh[3,4]<- 0
    hh[3,5]<- -sigma1E[2]/n01[2]+sigma1E[3]/n01[3]
    hh[3,6]<- sigma1E[1]/n01[1]-sigma1E[4]/n01[4]
    hh[4,4]<- sigma2E[1]/n02[1]+sigma2E[2]/n02[2]+sigma2E[3]/n02[3]+sigma2E[4]/n02[4]
    hh[4,5]<- sigma2E[2]/n02[2]-sigma2E[3]/n02[3]
    hh[4,6]<- -sigma2E[1]/n02[1]+sigma2E[4]/n02[4]
    hh[5,5]<- sigma1E[2]/n01[2]+sigma1E[3]/n01[3]+sigma2E[2]/n02[2]+sigma2E[3]/n02[3]
    hh[5,6]<- 0
    hh[6,6]<- sigma1E[1]/n01[1]+sigma1E[4]/n01[4]+sigma2E[1]/n02[1]+sigma2E[4]/n02[4]
    for(i in 2:6){
      for(j in 1:(i-1)){
        hh[i,j]<- hh[j,i]
      }
    }
    b_line<- matrix(0,6,1)
    b_line[1]<- sumx1/n_samP1+sumx2/n_samF1-sumwx_B1[1]/n01[1]-sumwx_B1[4]/n01[4]
    b_line[2]<- sumx2/n_samF1+sumx3/n_samP2-sumwx_B2[1]/n02[1]-sumwx_B2[4]/n02[4]
    b_line[3]<- sumwx_B1[1]/n01[1]-sumwx_B1[2]/n01[2]-sumwx_B1[3]/n01[3]+sumwx_B1[4]/n01[4]
    b_line[4]<- sumwx_B2[1]/n02[1]-sumwx_B2[2]/n02[2]-sumwx_B2[3]/n02[3]+sumwx_B2[4]/n02[4]
    b_line[5]<- sumwx_B1[2]/n01[2]-sumwx_B1[3]/n01[3]-sumwx_B2[2]/n02[2]+sumwx_B2[3]/n02[3]
    b_line[6]<- sumwx_B1[1]/n01[1]-sumwx_B1[4]/n01[4]-sumwx_B2[1]/n02[1]+sumwx_B2[4]/n02[4]
    B <- solve(hh,b_line)

    meanP1<- (sumx1-sigmaP1*B[1])/n_samP1
    meanF1<- (sumx2-sigmaP1*(B[1]+B[2]))/n_samF1
    meanP2<- (sumx3-sigmaP1*B[2])/n_samP2
    mean1E[1]<- (sumwx_B1[1]+sigma1E[1]*(B[1]-B[3]-B[6]))/n01[1]
    mean1E[2]<- (sumwx_B1[2]+sigma1E[1]*(B[3]-B[5]))/n01[2]
    mean1E[3]<- (sumwx_B1[3]+sigma1E[1]*(B[3]+B[5]))/n01[3]
    mean1E[4]<- (sumwx_B1[4]+sigma1E[1]*(B[1]-B[3]+B[6]))/n01[4]
    mean2E[1]<- (sumwx_B2[1]+sigma2E[1]*(B[2]-B[4]+B[6]))/n02[1]
    mean2E[2]<- (sumwx_B2[2]+sigma2E[1]*(B[4]+B[5]))/n02[2]
    mean2E[3]<- (sumwx_B2[3]+sigma2E[1]*(B[4]-B[5]))/n02[3]
    mean2E[4]<- (sumwx_B2[4]+sigma2E[1]*(B[2]-B[4]-B[6]))/n02[4]

    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2); ss2<- sum((dataF1-meanF1)^2); ss3<- sum((dataP2-meanP2)^2)
    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1E[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2E[i])^2 }

    s0<- matrix(0,6,1)
    s0[1]<- sum(swx_B1); sigma1E<- as.matrix((s0[1]/n_samB1),d21,1)
    s0[2]<- sum(swx_B2); sigma2E<- as.matrix((s0[2]/n_samB2),d22,1)
    aa1<- sigma1E[1]-sigmaP1;aa2<- sigma2E[1]-sigmaP1
    aa1[aa1<0] <- 0; aa2[aa2<0] <- 0
    aaa0<- sigmaP1; n_iter<- 0;aa5<- 1000
    s0[3]<- ss1+ss2+ss3; s0[4]<- n_samP1+n_samF1+n_samP2
    while(aa5>0.0001){
      n_iter<- n_iter+1
      aa3<- sigmaP1/(aa1+sigmaP1); aa4<- sigmaP1/(aa2+sigmaP1)
      sigmaP1<- (s0[3]+aa3*aa3*s0[1]+aa4*aa4*s0[2])/(s0[4]+aa3*n_samB1+aa4*n_samB2)
      aa5<- abs(sigmaP1-aaa0); aaa0<- sigmaP1
      if (n_iter>20) break
    }
    sigmaF1<- sigmaP1;sigmaP2<- sigmaP1
    sigma1E<- matrix((sigmaP1+aa1),d21,1)
    sigma2E<- matrix((sigmaP1+aa2),d22,1)
    if(sum(sigma1E < 1e-30)>=1){break}
    if(sum(sigma2E < 1e-30)>=1){break}
    ####################### the stop criterion for iteration ############################################################
    L1<- logL(n_samP1,1,1,meanP1,sigmaP1,dataP1)+logL(n_samF1,1,1,meanF1,sigmaF1,dataF1)+logL(n_samP2,1,1,meanP2,sigmaP2,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1E,sigma1E,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2E,sigma2E,dataB2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*7

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1,1,1,1,1, 1,0,-1,1,1,0,0,0,0,-1,-1, 1,0,-1,1,0,1,0,0,-1,0,-1, 1,0,-1,0.5,0.5,0.5,0.5,-0.5,-0.5,-0.5,-0.5, 0,1,0,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5),11,5)
  b_line1 <- matrix(c(meanP1,meanF1,meanP2,mean1E,mean2E))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigma1 - sigma1E[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigma1
  mm_1<- sigma1E[1]-sigmaP1
  if(mm_1 < 0) {mm_1 <- 0}
  nnn_1 <- mm_1/sigma1

  jj_2 <- sigma2 - sigma2E[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigma2
  mm_2<- sigma2E[1]-sigmaP1
  if(mm_2 < 0) {mm_2 <- 0}
  nnn_2 <- mm_2/sigma2

  ######### hypothesis testing #########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - meanP1)/sqrt(as.vector(sigmaP1))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[2])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),D_P1))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])

  ###############  F1  ###################
  dataF1<-sort(dataF1);bmw_F1 <- matrix(0,n_samF1,1); bmwsl_F1 <- matrix(0,n_samF1,1)

  gg_F1 <- (dataF1 - meanF1)/sqrt(as.vector(sigmaF1))
  bmw_F1[which(gg_F1>=0)] <- pnorm(gg_F1[gg_F1>=0])
  bmw_F1[which(gg_F1<0)] <- 1 - pnorm(abs(gg_F1[gg_F1<0]))
  bmwsl_F1[,1] <- bmw_F1

  P2_F1 <- rowSums(bmwsl_F1)
  nn<-dim(as.matrix(unique(P2_F1)))[1]
  if(nn<n_samF1){P2_F1<-P2_F1+runif(n_samF1)/1e4}

  dd_F1 <- as.matrix(c(sum(P2_F1),sum(P2_F1^2),sum((P2_F1-0.5)^2)))
  WW2_F1 <- 1/(12*n_samF1) + sum((P2_F1 - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  u_F1 <- as.matrix(c(12*n_samF1*((dd_F1[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((dd_F1[2]/n_samF1-1/3)^2),180*n_samF1*((dd_F1[3]/n_samF1-1/12)^2)))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[2])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),D_F1))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - meanP2)/sqrt(as.vector(sigmaP2))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[2])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),D_P2))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])

  ###############  B1  ###################
  dataB1<-sort(dataB1);bmw_B1 <- matrix(0,n_samB1,1); bmwsl_B1 <- matrix(0,n_samB1,d21)
  for(i in 1:d21){
    gg_B1 <- (dataB1 - mean1E[i])/sqrt(sigma1E[i])
    bmw_B1[which(gg_B1>=0)] <- pnorm(gg_B1[gg_B1>=0])
    bmw_B1[which(gg_B1<0)] <- 1 - pnorm(abs(gg_B1[gg_B1<0]))
    bmwsl_B1[,i] <- bmw_B1*mix_pi_1[i]
  }
  P2_B1 <- rowSums(bmwsl_B1)
  nn<-dim(as.matrix(unique(P2_B1)))[1]
  if(nn<n_samB1){P2_B1<-P2_B1+runif(n_samB1)/1e4}

  dd_B1 <- as.matrix(c(sum(P2_B1),sum(P2_B1^2),sum((P2_B1-0.5)^2)))
  WW2_B1 <- 1/(12*n_samB1) + sum((P2_B1 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  u_B1 <- as.matrix(c(12*n_samB1*((dd_B1[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((dd_B1[2]/n_samB1-1/3)^2),180*n_samB1*((dd_B1[3]/n_samB1-1/12)^2)))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[2])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),D_B1))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])

  ###############  B2  ###################
  dataB2<-sort(dataB2);bmw_B2 <- matrix(0,n_samB2,1); bmwsl_B2 <- matrix(0,n_samB2,d22)
  for(i in 1:d22){
    gg_B2 <- (dataB2 - mean2E[i])/sqrt(sigma2E[i])
    bmw_B2[which(gg_B2>=0)] <- pnorm(gg_B2[gg_B2>=0])
    bmw_B2[which(gg_B2<0)] <- 1 - pnorm(abs(gg_B2[gg_B2<0]))
    bmwsl_B2[,i] <- bmw_B2*mix_pi_2[i]
  }
  P2_B2 <- rowSums(bmwsl_B2)
  nn<-dim(as.matrix(unique(P2_B2)))[1]
  if(nn<n_samB2){P2_B2<-P2_B2+runif(n_samB2)/1e4}

  dd_B2 <- as.matrix(c(sum(P2_B2),sum(P2_B2^2),sum((P2_B2-0.5)^2)))
  WW2_B2 <- 1/(12*n_samB2) + sum((P2_B2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  u_B2 <- as.matrix(c(12*n_samB2*((dd_B2[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((dd_B2[2]/n_samB2-1/3)^2),180*n_samB2*((dd_B2[3]/n_samB2-1/12)^2)))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[2])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),D_B2))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX2-A-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigmaP1,4),round(t(mean1E),4),round(sigma1E[1],4),round(t(mix_pi_1),4),
                       round(t(mean2E),4),round(sigma2E[1],4),round(t(mix_pi_2),4),
                       round(B1[1],4),round(B1[2],4),round(B1[3],4)," "," ",round(B1[4],4),round(B1[5],4),round(jj_1,4),round(ll_1*100,4),round(mm_1,4),round(nnn_1*100,4) ,round(jj_2,4),round(ll_2*100,4),round(mm_2,4),round(nnn_2*100,4) ,
                       round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],
                       round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],
                       round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],
                       round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],
                       round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
}
########################### (E3)#####################################
G5BCModelFun[[16]] <- function(K1,logL,df11,df21,df31,df41,df51){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]))
  n_samP1<-dim(dataP1)[1];n_samP2<-dim(dataP2)[1];n_samF1<-dim(dataF1)[1];n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2)
  sumx4<-sum(dataB1);sumx5<-sum(dataB2)
  ss1<- (n_samP1-1)*as.numeric(var(dataP1));meanP1<- mean(dataP1)
  ss2<- (n_samF1-1)*as.numeric(var(dataF1));meanF1<- mean(dataF1)
  ss3<- (n_samP2-1)*as.numeric(var(dataP2));meanP2<- mean(dataP2)
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  sigmaP1<- sigma0;sigmaF1<- sigma0;sigmaP2<- sigma0
  mean1<- mean(dataB1); sigma1<- as.numeric(var(dataB1))
  mean2<- mean(dataB2); sigma2<- as.numeric(var(dataB2))

  m_esp<-0.0001
  ####################### MX2-EA-AD Model ############ (E3)#####################
  d21<-3; d22<-3
  mi_1<- as.matrix(c(0.25,0.5,0.25)); a1<- sqrt(sigma1/(n_samB1-1))
  sigma1E<- matrix(sigma1/3,d21,1); mean1E<- as.matrix(c(mean1+2*a1,mean1,mean1-2*a1))

  mi_2<- as.matrix(c(0.25,0.5,0.25)); a2<- sqrt(sigma2/(n_samB2-1))
  sigma2E<- matrix(sigma2/3,d22,1); mean2E<- as.matrix(c(mean2+2*a2,mean2,mean2-2*a2))

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1E,sigma1E,dataB1)+logL(n_samB2,d22,mi_2,mean2E,sigma2E,dataB2)
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1E[i],sqrt(sigma1E[i]))/dmixnorm(dataB1,mean1E,sqrt(sigma1E),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2E[i],sqrt(sigma2E[i]))/dmixnorm(dataB2,mean2E,sqrt(sigma2E),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    ############ CM1-step for means ##############
    n01 <- n_samB1*mix_pi_1; n01[n01<0.000001] <- 0.000001
    n02 <- n_samB2*mix_pi_2; n02[n02<0.000001] <- 0.000001

    hh<- matrix(0,5,5)
    hh[1,1]<- sigmaP1/n_samP1+sigmaP1/n_samF1+4.0*sigma1E[2]/n01[2]
    hh[1,2]<- sigmaP1/n_samF1
    hh[1,3]<- 4.0*sigma1E[2]/n01[2]
    hh[1,4]<- 0
    hh[1,5]<- 0
    hh[2,2]<- sigmaP1/n_samF1+sigmaP1/n_samP2+4.0*sigma2E[2]/n02[2]
    hh[2,3]<- 0
    hh[2,4]<- 4.0*sigma2E[2]/n02[2]
    hh[2,5]<- 0
    hh[3,3]<- sigma1E[1]/n01[1]+4.0*sigma1E[2]/n01[2]+sigma1E[3]/n01[3]
    hh[3,4]<- 0
    hh[3,5]<- sigma1E[1]/n01[1]-sigma1E[3]/n01[3]
    hh[4,4]<- sigma2E[1]/n02[1]+4.0*sigma2E[2]/n02[2]+sigma2E[3]/n02[3]
    hh[4,5]<- -sigma2E[1]/n02[1]+sigma2E[3]/n02[3]
    hh[5,5]<- sigma1E[1]/n01[1]+sigma1E[3]/n01[3]+sigma2E[1]/n02[1]+sigma2E[3]/n02[3]
    for(i in 2:5){
      for(j in 1:(i-1)){
        hh[i,j]<- hh[j,i]
      }
    }
    b_line<- matrix(0,5,1)
    b_line[1]<- sumx1/n_samP1+sumx2/n_samF1-2.0*sumwx_B1[2]/n01[2]
    b_line[2]<- sumx2/n_samF1+sumx3/n_samP2-2.0*sumwx_B2[2]/n02[2]
    b_line[3]<- sumwx_B1[1]/n01[1]-2.0*sumwx_B1[2]/n01[2]+sumwx_B1[3]/n01[3]
    b_line[4]<- sumwx_B2[1]/n02[1]-2.0*sumwx_B2[2]/n02[2]+sumwx_B2[3]/n02[3]
    b_line[5]<- sumwx_B1[1]/n01[1]-sumwx_B1[3]/n01[3]-sumwx_B2[1]/n02[1]+sumwx_B2[3]/n02[3]

    B <- solve(hh,b_line)

    meanP1<- (sumx1-sigmaP1*B[1])/n_samP1
    meanF1<- (sumx2-sigmaP1*(B[1]+B[2]))/n_samF1
    meanP2<- (sumx3-sigmaP1*B[2])/n_samP2
    mean1E[1]<- (sumwx_B1[1]-sigma1E[1]*(B[3]+B[5]))/n01[1]
    mean1E[2]<- (sumwx_B1[2]+sigma1E[1]*(2.0*B[1]+2.0*B[3]))/n01[2]
    mean1E[3]<- (sumwx_B1[3]-sigma1E[1]*(B[3]-B[5]))/n01[3]
    mean2E[1]<- (sumwx_B2[1]-sigma2E[1]*(B[4]-B[5]))/n02[1]
    mean2E[2]<- (sumwx_B2[2]+sigma2E[1]*(2.0*B[2]+2.0*B[4]))/n02[2]
    mean2E[3]<- (sumwx_B2[3]-sigma2E[1]*(B[4]+B[5]))/n02[3]

    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2); ss2<- sum((dataF1-meanF1)^2); ss3<- sum((dataP2-meanP2)^2)
    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1E[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2E[i])^2 }

    s0<- matrix(0,4,1)
    s0[1]<- sum(swx_B1); sigma1E<- as.matrix((s0[1]/n_samB1),d21,1)
    s0[2]<- sum(swx_B2); sigma2E<- as.matrix((s0[2]/n_samB2),d22,1)
    aa1<- sigma1E[1]-sigmaP1;aa2<- sigma2E[1]-sigmaP1
    aa1[aa1<0] <- 0; aa2[aa2<0] <- 0
    aaa0<- sigmaP1; n_iter<- 0;aa5<- 1000
    s0[3]<- ss1+ss2+ss3; s0[4]<- n_samP1+n_samF1+n_samP2
    while(aa5>0.0001){
      n_iter<- n_iter+1
      aa3<- sigmaP1/(aa1+sigmaP1); aa4<- sigmaP1/(aa2+sigmaP1)
      sigmaP1<- (s0[3]+aa3*aa3*s0[1]+aa4*aa4*s0[2])/(s0[4]+aa3*n_samB1+aa4*n_samB2)
      aa5<- abs(sigmaP1-aaa0); aaa0<- sigmaP1
      if (n_iter>20) break
    }
    sigmaF1<- sigmaP1;sigmaP2<- sigmaP1
    sigma1E<- matrix((sigmaP1+aa1),d21,1)
    sigma2E<- matrix((sigmaP1+aa2),d22,1)
    if(sum(sigma1E < 1e-30)>=1){break}
    if(sum(sigma2E < 1e-30)>=1){break}
    ####################### the stop criterion for iteration ############################################################
    L1<- logL(n_samP1,1,1,meanP1,sigmaP1,dataP1)+logL(n_samF1,1,1,meanF1,sigmaF1,dataF1)+logL(n_samP2,1,1,meanP2,sigmaP2,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1E,sigma1E,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2E,sigma2E,dataB2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*7

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1,1,1, 2,0,-2,2,1,0,0,-1,-2, 1,0,-1,0.5,0.5,0.5,-0.5,-0.5,-0.5, 0,1,0,0.5,0.5,0.5,0.5,0.5,0.5),9,4)
  b_line1 <- matrix(c(meanP1,meanF1,meanP2,mean1E,mean2E))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigma1 - sigma1E[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigma1
  mm_1<- sigma1E[1]-sigmaP1
  if(mm_1 < 0) {mm_1 <- 0}
  nnn_1 <- mm_1/sigma1

  jj_2 <- sigma2 - sigma2E[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigma2
  mm_2<- sigma2E[1]-sigmaP1
  if(mm_2 < 0) {mm_2 <- 0}
  nnn_2 <- mm_2/sigma2

  ######### hypothesis testing #########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - meanP1)/sqrt(as.vector(sigmaP1))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[2])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),D_P1))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])

  ###############  F1  ###################
  dataF1<-sort(dataF1);bmw_F1 <- matrix(0,n_samF1,1); bmwsl_F1 <- matrix(0,n_samF1,1)

  gg_F1 <- (dataF1 - meanF1)/sqrt(as.vector(sigmaF1))
  bmw_F1[which(gg_F1>=0)] <- pnorm(gg_F1[gg_F1>=0])
  bmw_F1[which(gg_F1<0)] <- 1 - pnorm(abs(gg_F1[gg_F1<0]))
  bmwsl_F1[,1] <- bmw_F1

  P2_F1 <- rowSums(bmwsl_F1)
  nn<-dim(as.matrix(unique(P2_F1)))[1]
  if(nn<n_samF1){P2_F1<-P2_F1+runif(n_samF1)/1e4}

  dd_F1 <- as.matrix(c(sum(P2_F1),sum(P2_F1^2),sum((P2_F1-0.5)^2)))
  WW2_F1 <- 1/(12*n_samF1) + sum((P2_F1 - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  u_F1 <- as.matrix(c(12*n_samF1*((dd_F1[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((dd_F1[2]/n_samF1-1/3)^2),180*n_samF1*((dd_F1[3]/n_samF1-1/12)^2)))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[2])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),D_F1))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - meanP2)/sqrt(as.vector(sigmaP2))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[2])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),D_P2))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])

  ###############  B1  ###################
  dataB1<-sort(dataB1);bmw_B1 <- matrix(0,n_samB1,1); bmwsl_B1 <- matrix(0,n_samB1,d21)
  for(i in 1:d21){
    gg_B1 <- (dataB1 - mean1E[i])/sqrt(sigma1E[i])
    bmw_B1[which(gg_B1>=0)] <- pnorm(gg_B1[gg_B1>=0])
    bmw_B1[which(gg_B1<0)] <- 1 - pnorm(abs(gg_B1[gg_B1<0]))
    bmwsl_B1[,i] <- bmw_B1*mix_pi_1[i]
  }
  P2_B1 <- rowSums(bmwsl_B1)
  nn<-dim(as.matrix(unique(P2_B1)))[1]
  if(nn<n_samB1){P2_B1<-P2_B1+runif(n_samB1)/1e4}

  dd_B1 <- as.matrix(c(sum(P2_B1),sum(P2_B1^2),sum((P2_B1-0.5)^2)))
  WW2_B1 <- 1/(12*n_samB1) + sum((P2_B1 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  u_B1 <- as.matrix(c(12*n_samB1*((dd_B1[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((dd_B1[2]/n_samB1-1/3)^2),180*n_samB1*((dd_B1[3]/n_samB1-1/12)^2)))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[2])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),D_B1))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])

  ###############  B2  ###################
  dataB2<-sort(dataB2);bmw_B2 <- matrix(0,n_samB2,1); bmwsl_B2 <- matrix(0,n_samB2,d22)
  for(i in 1:d22){
    gg_B2 <- (dataB2 - mean2E[i])/sqrt(sigma2E[i])
    bmw_B2[which(gg_B2>=0)] <- pnorm(gg_B2[gg_B2>=0])
    bmw_B2[which(gg_B2<0)] <- 1 - pnorm(abs(gg_B2[gg_B2<0]))
    bmwsl_B2[,i] <- bmw_B2*mix_pi_2[i]
  }
  P2_B2 <- rowSums(bmwsl_B2)
  nn<-dim(as.matrix(unique(P2_B2)))[1]
  if(nn<n_samB2){P2_B2<-P2_B2+runif(n_samB2)/1e4}

  dd_B2 <- as.matrix(c(sum(P2_B2),sum(P2_B2^2),sum((P2_B2-0.5)^2)))
  WW2_B2 <- 1/(12*n_samB2) + sum((P2_B2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  u_B2 <- as.matrix(c(12*n_samB2*((dd_B2[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((dd_B2[2]/n_samB2-1/3)^2),180*n_samB2*((dd_B2[3]/n_samB2-1/12)^2)))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[2])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),D_B2))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX2-EA-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigmaP1,4),round(t(mean1E),4)," ",round(sigma1E[1],4),round(t(mix_pi_1),4)," ",
                       round(t(mean2E),4)," ",round(sigma2E[1],4),round(t(mix_pi_2),4)," ",
                       round(B1[1],4),round(B1[2],4),round(B1[2],4)," "," ",round(B1[3],4),round(B1[4],4),round(jj_1,4),round(ll_1*100,4),round(mm_1,4),round(nnn_1*100,4) ,round(jj_2,4),round(ll_2*100,4),round(mm_2,4),round(nnn_2*100,4) ,
                       round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],
                       round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],
                       round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],
                       round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],
                       round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
}
########################### (E4)#####################################
G5BCModelFun[[17]] <- function(K1,logL,df11,df21,df31,df41,df51){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]))
  n_samP1<-dim(dataP1)[1];n_samP2<-dim(dataP2)[1];n_samF1<-dim(dataF1)[1];n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2)
  sumx4<-sum(dataB1);sumx5<-sum(dataB2)
  ss1<- (n_samP1-1)*as.numeric(var(dataP1));meanP1<- mean(dataP1)
  ss2<- (n_samF1-1)*as.numeric(var(dataF1));meanF1<- mean(dataF1)
  ss3<- (n_samP2-1)*as.numeric(var(dataP2));meanP2<- mean(dataP2)
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  sigmaP1<- sigma0;sigmaF1<- sigma0;sigmaP2<- sigma0
  mean1<- mean(dataB1); sigma1<- as.numeric(var(dataB1))
  mean2<- mean(dataB2); sigma2<- as.numeric(var(dataB2))

  m_esp<-0.0001
  ####################### MX2-CD-AD Model ############ (E4)#####################
  d21<-1; d22<-4
  mi_1<- as.matrix(1)
  sigma1E<- as.matrix(sigma1/3) ; mean1E<- as.matrix(mean1)

  mi_2<- matrix(0.25,d22,1); a2<- sqrt(sigma2/(n_samB2-1))
  sigma2E<- matrix(sigma2/3,d22,1); mean2E<- as.matrix(c(mean2+2*a2,mean2+0.5*a2,mean2-0.5*a2,mean2-2*a2))

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1E,sigma1E,dataB1)+logL(n_samB2,d22,mi_2,mean2E,sigma2E,dataB2)
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1E[i],sqrt(sigma1E[i]))/dmixnorm(dataB1,mean1E,sqrt(sigma1E),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2E[i],sqrt(sigma2E[i]))/dmixnorm(dataB2,mean2E,sqrt(sigma2E),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    ############ CM1-step for means ##############
    n01 <- n_samB1*mix_pi_1; n01[n01<0.000001] <- 0.000001
    n02 <- n_samB2*mix_pi_2; n02[n02<0.000001] <- 0.000001

    hh<- matrix(0,3,3)
    hh[1,1]<- sigmaP1/n_samP1+sigmaP1/n_samF1+4.0*sigma1E[1]/n_samB1
    hh[1,2]<- sigmaP1/n_samF1
    hh[1,3]<- 0
    hh[2,2]<- sigmaP1/n_samF1+sigmaP1/n_samP2+sigma2E[1]/n02[1]+sigma2E[4]/n02[4]
    hh[2,3]<- -sigma2E[1]/n02[1]-sigma2E[4]/n02[4]
    hh[3,3]<- sigma2E[1]/n02[1]+sigma2E[2]/n02[2]+sigma2E[3]/n02[3]+sigma2E[4]/n02[4]
    for(i in 2:3){
      for(j in 1:(i-1)){
        hh[i,j]<- hh[j,i]
      }
    }
    b_line<- matrix(0,3,1)
    b_line[1]<- sumx1/n_samP1+sumx2/n_samF1-2.0*sumwx_B1[1]/n01[1]
    b_line[2]<- sumx2/n_samF1+sumx3/n_samP2-sumwx_B2[1]/n02[1]-sumwx_B2[4]/n02[4]
    b_line[3]<- sumwx_B2[1]/n02[1]-sumwx_B2[2]/n02[2]-sumwx_B2[3]/n02[3]+sumwx_B2[4]/n02[4]

    B <- solve(hh,b_line)

    meanP1<- (sumx1-sigmaP1*B[1])/n_samP1
    meanF1<- (sumx2-sigmaP1*(B[1]+B[2]))/n_samF1
    meanP2<- (sumx3-sigmaP1*B[2])/n_samP2
    mean1E[1]<- (sumwx_B1[1]+sigma1E[1]*B[1]*2.0)/n01[1]
    mean2E[1]<- (sumwx_B2[1]+sigma2E[1]*(B[2]-B[3]))/n02[1]
    mean2E[2]<- (sumwx_B2[2]+sigma2E[1]*B[3])/n02[2]
    mean2E[3]<- (sumwx_B2[3]+sigma2E[1]*B[3])/n02[3]
    mean2E[4]<- (sumwx_B2[4]+sigma2E[1]*(B[2]-B[3]))/n02[4]
    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2); ss2<- sum((dataF1-meanF1)^2); ss3<- sum((dataP2-meanP2)^2)
    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1E[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2E[i])^2 }

    s0<- matrix(0,4,1)
    s0[1]<- sum(swx_B1); sigma1E<- s0[1]/n_samB1
    s0[2]<- sum(swx_B2); sigma2E<- as.matrix((s0[2]/n_samB2),d22,1)
    aa1<- sigma1E[1]-sigmaP1;aa2<- sigma2E[1]-sigmaP1
    aa1[aa1<0] <- 0; aa2[aa2<0] <- 0
    aaa0<- sigmaP1; n_iter<- 0;aa5<- 1000
    s0[3]<- ss1+ss2+ss3; s0[4]<- n_samP1+n_samF1+n_samP2
    while(aa5>0.0001){
      n_iter<- n_iter+1
      aa3<- sigmaP1/(aa1+sigmaP1); aa4<- sigmaP1/(aa2+sigmaP1)
      sigmaP1<- (s0[3]+aa3*aa3*s0[1]+aa4*aa4*s0[2])/(s0[4]+aa3*n_samB1+aa4*n_samB2)
      aa5<- abs(sigmaP1-aaa0); aaa0<- sigmaP1
      if (n_iter>20) break
    }
    sigmaF1<- sigmaP1;sigmaP2<- sigmaP1
    sigma1E<- matrix((sigmaP1+aa1),d21,1)
    sigma2E<- matrix((sigmaP1+aa2),d22,1)
    if(sum(sigma1E < 1e-30)>=1){break}
    if(sum(sigma2E < 1e-30)>=1){break}
    ####################### the stop criterion for iteration ############################################################
    L1<- logL(n_samP1,1,1,meanP1,sigmaP1,dataP1)+logL(n_samF1,1,1,meanF1,sigmaF1,dataF1)+logL(n_samP2,1,1,meanP2,sigmaP2,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1E,sigma1E,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2E,sigma2E,dataB2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*7

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1,1, 1,1,-1,1,1,1,-1,-1, 1,1,-1,1,1,-1,1,-1, 1,0,-1,0.5,-0.5,-0.5,-0.5,-0.5, 0,1,0,0.5,0.5,0.5,0.5,0.5),8,5)
  b_line1 <- matrix(c(meanP1,meanF1,meanP2,mean1E,mean2E))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigma1 - sigma1E[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigma1
  mm_1<- sigma1E[1]-sigmaP1
  if(mm_1 < 0) {mm_1 <- 0}
  nnn_1 <- mm_1/sigma1

  jj_2 <- sigma2 - sigma2E[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigma2
  mm_2<- sigma2E[1]-sigmaP1
  if(mm_2 < 0) {mm_2 <- 0}
  nnn_2 <- mm_2/sigma2

  ######### hypothesis testing #########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - meanP1)/sqrt(as.vector(sigmaP1))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[2])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),D_P1))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])

  ###############  F1  ###################
  dataF1<-sort(dataF1);bmw_F1 <- matrix(0,n_samF1,1); bmwsl_F1 <- matrix(0,n_samF1,1)

  gg_F1 <- (dataF1 - meanF1)/sqrt(as.vector(sigmaF1))
  bmw_F1[which(gg_F1>=0)] <- pnorm(gg_F1[gg_F1>=0])
  bmw_F1[which(gg_F1<0)] <- 1 - pnorm(abs(gg_F1[gg_F1<0]))
  bmwsl_F1[,1] <- bmw_F1

  P2_F1 <- rowSums(bmwsl_F1)
  nn<-dim(as.matrix(unique(P2_F1)))[1]
  if(nn<n_samF1){P2_F1<-P2_F1+runif(n_samF1)/1e4}

  dd_F1 <- as.matrix(c(sum(P2_F1),sum(P2_F1^2),sum((P2_F1-0.5)^2)))
  WW2_F1 <- 1/(12*n_samF1) + sum((P2_F1 - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  u_F1 <- as.matrix(c(12*n_samF1*((dd_F1[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((dd_F1[2]/n_samF1-1/3)^2),180*n_samF1*((dd_F1[3]/n_samF1-1/12)^2)))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[2])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),D_F1))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - meanP2)/sqrt(as.vector(sigmaP2))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[2])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),D_P2))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])

  ###############  B1  ###################
  dataB1<-sort(dataB1);bmw_B1 <- matrix(0,n_samB1,1); bmwsl_B1 <- matrix(0,n_samB1,d21)
  for(i in 1:d21){
    gg_B1 <- (dataB1 - mean1E[i])/sqrt(as.vector(sigma1E[i]))
    bmw_B1[which(gg_B1>=0)] <- pnorm(gg_B1[gg_B1>=0])
    bmw_B1[which(gg_B1<0)] <- 1 - pnorm(abs(gg_B1[gg_B1<0]))
    bmwsl_B1[,i] <- bmw_B1*mix_pi_1[i]
  }
  P2_B1 <- rowSums(bmwsl_B1)
  nn<-dim(as.matrix(unique(P2_B1)))[1]
  if(nn<n_samB1){P2_B1<-P2_B1+runif(n_samB1)/1e4}

  dd_B1 <- as.matrix(c(sum(P2_B1),sum(P2_B1^2),sum((P2_B1-0.5)^2)))
  WW2_B1 <- 1/(12*n_samB1) + sum((P2_B1 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  u_B1 <- as.matrix(c(12*n_samB1*((dd_B1[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((dd_B1[2]/n_samB1-1/3)^2),180*n_samB1*((dd_B1[3]/n_samB1-1/12)^2)))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[2])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),D_B1))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])

  ###############  B2  ###################
  dataB2<-sort(dataB2);bmw_B2 <- matrix(0,n_samB2,1); bmwsl_B2 <- matrix(0,n_samB2,d22)
  for(i in 1:d22){
    gg_B2 <- (dataB2 - mean2E[i])/sqrt(sigma2E[i])
    bmw_B2[which(gg_B2>=0)] <- pnorm(gg_B2[gg_B2>=0])
    bmw_B2[which(gg_B2<0)] <- 1 - pnorm(abs(gg_B2[gg_B2<0]))
    bmwsl_B2[,i] <- bmw_B2*mix_pi_2[i]
  }
  P2_B2 <- rowSums(bmwsl_B2)
  nn<-dim(as.matrix(unique(P2_B2)))[1]
  if(nn<n_samB2){P2_B2<-P2_B2+runif(n_samB2)/1e4}

  dd_B2 <- as.matrix(c(sum(P2_B2),sum(P2_B2^2),sum((P2_B2-0.5)^2)))
  WW2_B2 <- 1/(12*n_samB2) + sum((P2_B2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  u_B2 <- as.matrix(c(12*n_samB2*((dd_B2[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((dd_B2[2]/n_samB2-1/3)^2),180*n_samB2*((dd_B2[3]/n_samB2-1/12)^2)))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[2])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),D_B2))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX2-CD-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigmaP1,4),round(t(mean1E),4)," "," "," ",round(sigma1E[1],4),round(t(mix_pi_1),4)," "," "," ",
                       round(t(mean2E),4),round(sigma2E[1],4),round(t(mix_pi_2),4),
                       round(B1[1],4),round(B1[2],4),round(B1[3],4),round(B1[2],4),round(B1[3],4),round(B1[4],4),round(B1[5],4),round(jj_1,4),round(ll_1*100,4),round(mm_1,4),round(nnn_1*100,4) ,round(jj_2,4),round(ll_2*100,4),round(mm_2,4),round(nnn_2*100,4) ,
                       round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],
                       round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],
                       round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],
                       round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],
                       round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
}
########################### (E5)#####################################
G5BCModelFun[[18]] <- function(K1,logL,df11,df21,df31,df41,df51){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]))
  n_samP1<-dim(dataP1)[1];n_samP2<-dim(dataP2)[1];n_samF1<-dim(dataF1)[1];n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2)
  sumx4<-sum(dataB1);sumx5<-sum(dataB2)
  ss1<- (n_samP1-1)*as.numeric(var(dataP1));meanP1<- mean(dataP1)
  ss2<- (n_samF1-1)*as.numeric(var(dataF1));meanF1<- mean(dataF1)
  ss3<- (n_samP2-1)*as.numeric(var(dataP2));meanP2<- mean(dataP2)
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  sigmaP1<- sigma0;sigmaF1<- sigma0;sigmaP2<- sigma0
  mean1<- mean(dataB1); sigma1<- as.numeric(var(dataB1))
  mean2<- mean(dataB2); sigma2<- as.numeric(var(dataB2))

  m_esp<-0.0001
  ####################### MX2-EAD-AD Model ############ (E5)#####################
  d21<-1; d22<-3
  mi_1<- as.matrix(1)
  sigma1E<- as.matrix(sigma1/3) ; mean1E<- as.matrix(mean1)

  mi_2<- as.matrix(c(0.25,0.5,0.25)); a2<- sqrt(sigma2/(n_samB2-1))
  sigma2E<- matrix(sigma2/3,d22,1); mean2E<- as.matrix(c(mean2+2*a2,mean2,mean2-2*a2))

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1E,sigma1E,dataB1)+logL(n_samB2,d22,mi_2,mean2E,sigma2E,dataB2)
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1E[i],sqrt(sigma1E[i]))/dmixnorm(dataB1,mean1E,sqrt(sigma1E),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2E[i],sqrt(sigma2E[i]))/dmixnorm(dataB2,mean2E,sqrt(sigma2E),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    ############ CM1-step for means ##############
    n01 <- n_samB1*mix_pi_1; n01[n01<0.000001] <- 0.000001
    n02 <- n_samB2*mix_pi_2; n02[n02<0.000001] <- 0.000001

    hh<- matrix(0,3,3)
    hh[1,1]<- sigmaP1/n_samP1+sigmaP1/n_samF1+4.0*sigma1E[1]/n_samB1
    hh[1,2]<- sigmaP1/n_samF1
    hh[1,3]<- 0
    hh[2,2]<- sigmaP1/n_samF1+sigmaP1/n_samP2+4.0*sigma2E[2]/n02[2]
    hh[2,3]<- 4.0*sigma2E[2]/n02[2]
    hh[3,3]<- sigma2E[1]/n02[1]+4.0*sigma2E[2]/n02[2]+sigma2E[3]/n02[3]
    for(i in 2:3){
      for(j in 1:(i-1)){
        hh[i,j]<- hh[j,i]
      }
    }
    b_line<- matrix(0,3,1)
    b_line[1]<- sumx1/n_samP1+sumx2/n_samF1-2.0*sumwx_B1[1]/n_samB1
    b_line[2]<- sumx2/n_samF1+sumx3/n_samP2-2.0*sumwx_B2[2]/n02[2]
    b_line[3]<- sumwx_B2[1]/n02[1]-2.0*sumwx_B2[2]/n02[2]+sumwx_B2[3]/n02[3]

    B <- solve(hh,b_line)

    meanP1<- (sumx1-sigmaP1*B[1])/n_samP1
    meanF1<- (sumx2-sigmaF1*(B[1]+B[2]))/n_samF1
    meanP2<- (sumx3-sigmaP2*B[2])/n_samP2
    mean1E[1]<- (sumwx_B1[1]+sigma1E[1]*B[1]*2.0)/n01[1]
    mean2E[1]<- (sumwx_B2[1]-sigma2E[1]*B[3])/n02[1]
    mean2E[2]<- (sumwx_B2[2]+sigma2E[1]*(2.0*B[2]+2.0*B[3]))/n02[2]
    mean2E[3]<- (sumwx_B2[3]-sigma2E[1]*B[3])/n02[3]
    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2); ss2<- sum((dataF1-meanF1)^2); ss3<- sum((dataP2-meanP2)^2)
    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1E[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2E[i])^2 }

    s0<- matrix(0,4,1)
    s0[1]<- sum(swx_B1); sigma1E<- s0[1]/n_samB1
    s0[2]<- sum(swx_B2); sigma2E<- as.matrix((s0[2]/n_samB2),d22,1)
    aa1<- sigma1E[1]-sigmaP1;aa2<- sigma2E[1]-sigmaP1
    aa1[aa1<0] <- 0; aa2[aa2<0] <- 0
    aaa0<- sigmaP1; n_iter<- 0;aa5<- 1000
    s0[3]<- ss1+ss2+ss3; s0[4]<- n_samP1+n_samF1+n_samP2
    while(aa5>0.0001){
      n_iter<- n_iter+1
      aa3<- sigmaP1/(aa1+sigmaP1); aa4<- sigmaP1/(aa2+sigmaP1)
      sigmaP1<- (s0[3]+aa3*aa3*s0[1]+aa4*aa4*s0[2])/(s0[4]+aa3*n_samB1+aa4*n_samB2)
      aa5<- abs(sigmaP1-aaa0); aaa0<- sigmaP1
      if (n_iter>20) break
    }
    sigmaF1<- sigmaP1;sigmaP2<- sigmaP1
    sigma1E<- matrix((sigmaP1+aa1),d21,1)
    sigma2E<- matrix((sigmaP1+aa2),d22,1)
    if(sum(sigma1E < 1e-30)>=1){break}
    if(sum(sigma2E < 1e-30)>=1){break}
    ####################### the stop criterion for iteration ############################################################
    L1<- logL(n_samP1,1,1,meanP1,sigmaP1,dataP1)+logL(n_samF1,1,1,meanF1,sigmaF1,dataF1)+logL(n_samP2,1,1,meanP2,sigmaP2,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1E,sigma1E,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2E,sigma2E,dataB2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*7

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1, 2,2,-2,2,2,0,-2, 1,0,-1,0.5,-0.5,-0.5,-0.5, 0,1,0,0.5,0.5,0.5,0.5),7,4)
  b_line1 <- matrix(c(meanP1,meanF1,meanP2,mean1E,mean2E))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigma1 - sigma1E[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigma1
  mm_1<- sigma1E[1]-sigmaP1
  if(mm_1 < 0) {mm_1 <- 0}
  nnn_1 <- mm_1/sigma1

  jj_2 <- sigma2 - sigma2E[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigma2
  mm_2<- sigma2E[1]-sigmaP1
  if(mm_2 < 0) {mm_2 <- 0}
  nnn_2 <- mm_2/sigma2

  ######### hypothesis testing #########
  ###############  P1  ###################
  dataP1<-sort(dataP1);bmw_P1 <- matrix(0,n_samP1,1); bmwsl_P1 <- matrix(0,n_samP1,1)

  gg_P1 <- (dataP1 - meanP1)/sqrt(as.vector(sigmaP1))
  bmw_P1[which(gg_P1>=0)] <- pnorm(gg_P1[gg_P1>=0])
  bmw_P1[which(gg_P1<0)] <- 1 - pnorm(abs(gg_P1[gg_P1<0]))
  bmwsl_P1[,1] <- bmw_P1

  P2_P1 <- rowSums(bmwsl_P1)
  nn<-dim(as.matrix(unique(P2_P1)))[1]
  if(nn<n_samP1){P2_P1<-P2_P1+runif(n_samP1)/1e4}

  dd_P1 <- as.matrix(c(sum(P2_P1),sum(P2_P1^2),sum((P2_P1-0.5)^2)))
  WW2_P1 <- 1/(12*n_samP1) + sum((P2_P1 - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  u_P1 <- as.matrix(c(12*n_samP1*((dd_P1[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((dd_P1[2]/n_samP1-1/3)^2),180*n_samP1*((dd_P1[3]/n_samP1-1/12)^2)))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[2])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),D_P1))
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])

  ###############  F1  ###################
  dataF1<-sort(dataF1);bmw_F1 <- matrix(0,n_samF1,1); bmwsl_F1 <- matrix(0,n_samF1,1)

  gg_F1 <- (dataF1 - meanF1)/sqrt(as.vector(sigmaF1))
  bmw_F1[which(gg_F1>=0)] <- pnorm(gg_F1[gg_F1>=0])
  bmw_F1[which(gg_F1<0)] <- 1 - pnorm(abs(gg_F1[gg_F1<0]))
  bmwsl_F1[,1] <- bmw_F1

  P2_F1 <- rowSums(bmwsl_F1)
  nn<-dim(as.matrix(unique(P2_F1)))[1]
  if(nn<n_samF1){P2_F1<-P2_F1+runif(n_samF1)/1e4}

  dd_F1 <- as.matrix(c(sum(P2_F1),sum(P2_F1^2),sum((P2_F1-0.5)^2)))
  WW2_F1 <- 1/(12*n_samF1) + sum((P2_F1 - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  u_F1 <- as.matrix(c(12*n_samF1*((dd_F1[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((dd_F1[2]/n_samF1-1/3)^2),180*n_samF1*((dd_F1[3]/n_samF1-1/12)^2)))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[2])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),D_F1))
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])

  ###############  P2  ###################
  dataP2<-sort(dataP2);bmw_P2 <- matrix(0,n_samP2,1); bmwsl_P2 <- matrix(0,n_samP2,1)

  gg_P2 <- (dataP2 - meanP2)/sqrt(as.vector(sigmaP2))
  bmw_P2[which(gg_P2>=0)] <- pnorm(gg_P2[gg_P2>=0])
  bmw_P2[which(gg_P2<0)] <- 1 - pnorm(abs(gg_P2[gg_P2<0]))
  bmwsl_P2[,1] <- bmw_P2

  P2_P2 <- rowSums(bmwsl_P2)
  nn<-dim(as.matrix(unique(P2_P2)))[1]
  if(nn<n_samP2){P2_P2<-P2_P2+runif(n_samP2)/1e4}

  dd_P2 <- as.matrix(c(sum(P2_P2),sum(P2_P2^2),sum((P2_P2-0.5)^2)))
  WW2_P2 <- 1/(12*n_samP2) + sum((P2_P2 - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  u_P2 <- as.matrix(c(12*n_samP2*((dd_P2[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((dd_P2[2]/n_samP2-1/3)^2),180*n_samP2*((dd_P2[3]/n_samP2-1/12)^2)))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[2])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),D_P2))
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])

  ###############  B1  ###################
  dataB1<-sort(dataB1);bmw_B1 <- matrix(0,n_samB1,1); bmwsl_B1 <- matrix(0,n_samB1,d21)
  for(i in 1:d21){
    gg_B1 <- (dataB1 - mean1E[i])/sqrt(as.vector(sigma1E[i]))
    bmw_B1[which(gg_B1>=0)] <- pnorm(gg_B1[gg_B1>=0])
    bmw_B1[which(gg_B1<0)] <- 1 - pnorm(abs(gg_B1[gg_B1<0]))
    bmwsl_B1[,i] <- bmw_B1*mix_pi_1[i]
  }
  P2_B1 <- rowSums(bmwsl_B1)
  nn<-dim(as.matrix(unique(P2_B1)))[1]
  if(nn<n_samB1){P2_B1<-P2_B1+runif(n_samB1)/1e4}

  dd_B1 <- as.matrix(c(sum(P2_B1),sum(P2_B1^2),sum((P2_B1-0.5)^2)))
  WW2_B1 <- 1/(12*n_samB1) + sum((P2_B1 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  u_B1 <- as.matrix(c(12*n_samB1*((dd_B1[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((dd_B1[2]/n_samB1-1/3)^2),180*n_samB1*((dd_B1[3]/n_samB1-1/12)^2)))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[2])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),D_B1))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])

  ###############  B2  ###################
  dataB2<-sort(dataB2);bmw_B2 <- matrix(0,n_samB2,1); bmwsl_B2 <- matrix(0,n_samB2,d22)
  for(i in 1:d22){
    gg_B2 <- (dataB2 - mean2E[i])/sqrt(sigma2E[i])
    bmw_B2[which(gg_B2>=0)] <- pnorm(gg_B2[gg_B2>=0])
    bmw_B2[which(gg_B2<0)] <- 1 - pnorm(abs(gg_B2[gg_B2<0]))
    bmwsl_B2[,i] <- bmw_B2*mix_pi_2[i]
  }
  P2_B2 <- rowSums(bmwsl_B2)
  nn<-dim(as.matrix(unique(P2_B2)))[1]
  if(nn<n_samB2){P2_B2<-P2_B2+runif(n_samB2)/1e4}

  dd_B2 <- as.matrix(c(sum(P2_B2),sum(P2_B2^2),sum((P2_B2-0.5)^2)))
  WW2_B2 <- 1/(12*n_samB2) + sum((P2_B2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  u_B2 <- as.matrix(c(12*n_samB2*((dd_B2[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((dd_B2[2]/n_samB2-1/3)^2),180*n_samB2*((dd_B2[3]/n_samB2-1/12)^2)))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[2])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),D_B2))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX2-EAD-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigmaP1,4),round(t(mean1E),4)," "," "," ",round(sigma1E[1],4),round(t(mix_pi_1),4)," "," "," ",
                       round(t(mean2E),4)," ",round(sigma2E[1],4),round(t(mix_pi_2),4)," ",
                       round(B1[1],4),round(B1[2],4),round(B1[2],4),round(B1[2],4),round(B1[2],4),round(B1[3],4),round(B1[4],4),round(jj_1,4),round(ll_1*100,4),round(mm_1,4),round(nnn_1*100,4) ,round(jj_2,4),round(ll_2*100,4),round(mm_2,4),round(nnn_2*100,4) ,
                       round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],
                       round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],
                       round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],
                       round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],
                       round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
}



K1G5BC <- function(x){
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

logLG5BC <- function(nm,nng,mi,mn,s,d1) { sum2 <- sum(log(dmixnorm(d1,mn,sqrt(s),mi)));return (sum2) }


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
allresult=foreach(i=1:18,.combine = 'rbind')%dopar%{
  requireNamespace("KScorrect")
  G5BCModelFun[[i]](K1G5BC,logLG5BC,df11,df21,df31,df41,df51)[[1]]
}
stopCluster(cl)

mi_1<-NULL;mi_2<-NULL

}else{


allresultq=switch(model,"1MG-AD" = G5BCModelFun[[1]](K1G5BC,logLG5BC,df11,df21,df31,df41,df51),"1MG-A"=G5BCModelFun[[2]](K1G5BC,logLG5BC,df11,df21,df31,df41,df51),"1MG-EAD"=G5BCModelFun[[3]](K1G5BC,logLG5BC,df11,df21,df31,df41,df51),"1MG-NCD"=G5BCModelFun[[4]](K1G5BC,logLG5BC,df11,df21,df31,df41,df51),
         "2MG-AD"=G5BCModelFun[[5]](K1G5BC,logLG5BC,df11,df21,df31,df41,df51),"2MG-A"=G5BCModelFun[[6]](K1G5BC,logLG5BC,df11,df21,df31,df41,df51),"2MG-EA"=G5BCModelFun[[7]](K1G5BC,logLG5BC,df11,df21,df31,df41,df51),"2MG-CD"=G5BCModelFun[[8]](K1G5BC,logLG5BC,df11,df21,df31,df41,df51),
         "2MG-EAD"=G5BCModelFun[[9]](K1G5BC,logLG5BC,df11,df21,df31,df41,df51),"MX1-AD-AD"=G5BCModelFun[[10]](K1G5BC,logLG5BC,df11,df21,df31,df41,df51),"MX1-A-AD"=G5BCModelFun[[11]](K1G5BC,logLG5BC,df11,df21,df31,df41,df51),"MX1-EAD-AD"=G5BCModelFun[[12]](K1G5BC,logLG5BC,df11,df21,df31,df41,df51),
         "MX1-NCD-AD"=G5BCModelFun[[13]](K1G5BC,logLG5BC,df11,df21,df31,df41,df51),"MX2-AD-AD"=G5BCModelFun[[14]](K1G5BC,logLG5BC,df11,df21,df31,df41,df51),"MX2-A-AD"=G5BCModelFun[[15]](K1G5BC,logLG5BC,df11,df21,df31,df41,df51),"MX2-EA-AD"=G5BCModelFun[[16]](K1G5BC,logLG5BC,df11,df21,df31,df41,df51),
         "MX2-CD-AD"=G5BCModelFun[[17]](K1G5BC,logLG5BC,df11,df21,df31,df41,df51),"MX2-EAD-AD"=G5BCModelFun[[18]](K1G5BC,logLG5BC,df11,df21,df31,df41,df51))


allresult<-allresultq[[1]]

mi_1<-allresultq[[2]];mi_2<-allresultq[[3]]

}
colnames(allresult) <- G5BCcolname
out<-list(allresult,mi_1,mi_2)
return(out)
}





