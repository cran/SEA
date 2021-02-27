G6Fun<-function(df,model){

data<-sapply(df,as.character)
dP1<-data[-1,which(data[1,]=="P1")];P1<-as.numeric(dP1[which(is.na(as.numeric(dP1))==FALSE)]);df11<-as.data.frame(P1)
dF1<-data[-1,which(data[1,]=="F1")];F1<-as.numeric(dF1[which(is.na(as.numeric(dF1))==FALSE)]);df21<-as.data.frame(F1)
dP2<-data[-1,which(data[1,]=="P2")];P2<-as.numeric(dP2[which(is.na(as.numeric(dP2))==FALSE)]);df31<-as.data.frame(P2)
dB1<-data[-1,which(data[1,]=="B1")];B1<-as.numeric(dB1[which(is.na(as.numeric(dB1))==FALSE)]);df41<-as.data.frame(B1)
dB2<-data[-1,which(data[1,]=="B2")];B2<-as.numeric(dB2[which(is.na(as.numeric(dB2))==FALSE)]);df51<-as.data.frame(B2)
dF2<-data[-1,which(data[1,]=="F2")];F2<-as.numeric(dF2[which(is.na(as.numeric(dF2))==FALSE)]);df61<-as.data.frame(F2)

G6colname <- c("Model","Log_Max_likelihood_Value","AIC","mean[P1]","mean[F1]","mean[P2]","Var(P1 & P2 & F1)","B1-mean[1]","B1-mean[2]","B1-mean[3]","B1-mean[4]",
               "B1-Var(Residual+Polygene)","B1-Proportion[1]","B1-Proportion[2]","B1-Proportion[3]","B1-Proportion[4]","B2-mean[1]","B2-mean[2]","B2-mean[3]","B2-mean[4]",
               "B2-Var(Residual+Polygene)","B2-Proportion[1]","B2-Proportion[2]","B2-Proportion[3]","B2-Proportion[4]","F2-mean[1]","F2-mean[2]","F2-mean[3]","F2-mean[4]",
               "F2-mean[5]","F2-mean[6]","F2-mean[7]","F2-mean[8]","F2-mean[9]","F2-Var(Residual+Polygene)","F2-Proportion[1]","F2-Proportion[2]","F2-Proportion[3]","F2-Proportion[4]",
               "F2-Proportion[5]","F2-Proportion[6]","F2-Proportion[7]","F2-Proportion[8]","F2-Proportion[9]","m(m1)","m2","m3","m4","m5","m6","da","db","ha","hb","i","jab","jba","l",
               "[d]","[h]","B1-MajorGene Var","B1-Heritability(MajorGene)(%)","B1-Polygenes Var","B1-Heritability(Polygenes)(%)","B2-MajorGene Var","B2-Heritability(MajorGene)(%)",
               "B2-Polygenes Var","B2-Heritability(Polygenes)(%)",	"F2-MajorGene Var","F2-Heritability(MajorGene)(%)","F2-Polygenes Var","F2-Heritability(Polygenes)(%)",
               "U1 square-P1","P(U1 square-P1)","U2 square-P1","P(U2 square-P1)","U3 square-P1","P(U3 square-P1)","nW square-P1","P(nW square-P1)","Dn-P1","P(Dn-P1)","U1 square-F1","P(U1 square-F1)","U2 square-F1","P(U2 square-F1)","U3 square-F1","P(U3 square-F1)","nW square-F1","P(nW square-F1)","Dn-F1","P(Dn-F1)",
               "U1 square-P2","P(U1 square-P2)","U2 square-P2","P(U2 square-P2)","U3 square-P2","P(U3 square-P2)","nW square-P2","P(nW square-P2)","Dn-P2","P(Dn-P2)","U1 square-B1","P(U1 square-B1)","U2 square-B1","P(U2 square-B1)","U3 square-B1","P(U3 square-B1)","nW square-B1","P(nW square-B1)","Dn-B1","P(Dn-B1)",
               "U1 square-B2","P(U1 square-B2)","U2 square-B2","P(U2 square-B2)","U3 square-B2","P(U3 square-B2)","nW square-B2","P(nW square-B2)","Dn-B2","P(Dn-B2)","U1 square-F2","P(U1 square-F2)","U2 square-F2","P(U2 square-F2)","U3 square-F2","P(U3 square-F2)","nW square-F2","P(nW square-F2)","Dn-F2","P(Dn-F2)")

G6ModelFun<-list(NA)
###################define each model function##################
########################### (A1)#####################################
G6ModelFun[[1]] <- function(K1,logL,df11,df21,df31,df41,df51,df61){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1];n_samP2<-dim(dataP2)[1];n_samF1<-dim(dataF1)[1]
  n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1];n_samF2<-dim(dataF2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2)
  sumx4<-sum(dataB1);sumx5<-sum(dataB2);sumx6<-sum(dataF2)
  sigmaP1<- as.numeric(var(dataP1)); ss1<- (n_samP1-1)*sigmaP1
  sigmaF1<- as.numeric(var(dataF1)); ss2<- (n_samF1-1)*sigmaF1
  sigmaP2<- as.numeric(var(dataP2)); ss3<- (n_samP2-1)*sigmaP2
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  mean4<- mean(dataB1); sigmaB1<- as.numeric(var(dataB1)); sigma40<- sigmaB1
  mean5<- mean(dataB2); sigmaB2<- as.numeric(var(dataB2)); sigma50<- sigmaB2
  mean6<- mean(dataF2); sigmaF2<- as.numeric(var(dataF2)); sigma60<- sigmaF2

  m_esp<-0.0001
  ####################### 1MG-AD Model ############ (A1)#####################
  d21<-2; d22<-2; d23<-3
  sigma<- sigma0
  meanP1<- mean(dataP1);meanF1<- mean(dataF1);meanP2<- mean(dataP2)

  mi_1<- matrix(0.5,d21,1)
  sigma1A<- matrix(sigma,d21,1);mean1A<- as.matrix(c(meanP1,meanF1))

  mi_2<- matrix(0.5,d22,1)
  sigma2A<- matrix(sigma,d22,1);mean2A<- as.matrix(c(meanF1,meanP2))

  mi_3<- as.matrix(c(0.25,0.5,0.25))
  sigma3A<- matrix(sigma,d23,1);mean3A<- as.matrix(c(meanP1,meanF1,meanP2))

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1A,sigma1A,dataB1)+logL(n_samB2,d22,mi_2,mean2A,sigma2A,dataB2)+logL(n_samF2,d23,mi_3,mean3A,sigma3A,dataF2)
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)
  WW_F2 <- matrix(0,d23,n_samF2); swx_F2 <- matrix(0,d23,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1A[i],sqrt(sigma1A[i]))/dmixnorm(dataB1,mean1A,sqrt(sigma1A),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2A[i],sqrt(sigma2A[i]))/dmixnorm(dataB2,mean2A,sqrt(sigma2A),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    for(i in 1:d23) { WW_F2[i,] <- mi_3[i]*dnorm(dataF2,mean3A[i],sqrt(sigma3A[i]))/dmixnorm(dataF2,mean3A,sqrt(sigma3A),mi_3) }
    mix_pi_3 <- as.matrix(rowSums(WW_F2)/n_samF2)
    sumwx_F2 <- WW_F2%*%dataF2

    ############ CM1-step for means ##############
    s0<- matrix(0,6,1)
    s0[1]<- sumx1+sumwx_B1[1]+sumwx_F2[1]
    s0[2]<- sumx2+sumwx_B1[2]+sumwx_B2[1]+sumwx_F2[2]
    s0[3]<- sumx3+sumwx_B2[2]+sumwx_F2[3]
    n0<- matrix(0,6,1)
    n0[1]<- n_samP1+mix_pi_1[1]*n_samB1+mix_pi_3[1]*n_samF2
    n0[2]<- n_samF1+mix_pi_1[2]*n_samB1+mix_pi_2[1]*n_samB2+mix_pi_3[2]*n_samF2
    n0[3]<- n_samP2+mix_pi_2[2]*n_samB2+mix_pi_3[3]*n_samF2
    meanP1<- s0[1]/n0[1];meanF1<- s0[2]/n0[2];meanP2<- s0[3]/n0[3]
    mean1A<- as.matrix(c(meanP1,meanF1));mean2A<- as.matrix(c(meanF1,meanP2));mean3A<- as.matrix(c(meanP1,meanF1,meanP2))

    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2)
    ss2<- sum((dataF1-meanF1)^2)
    ss3<- sum((dataP2-meanP2)^2)

    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1A[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2A[i])^2 }
    for(i in 1:d23) {  swx_F2[i] <- WW_F2[i,]%*%(dataF2-mean3A[i])^2 }

    s0[5]<- ss1+ss2+ss3+sum(swx_B1)+sum(swx_B2)+sum(swx_F2)
    n0[5]<- n_samP1+n_samF1+n_samP2+n_samB1+n_samB2+n_samF2
    sigma<- s0[5]/n0[5]
    sigma1A<- matrix(sigma,d21,1);sigma2A<- matrix(sigma,d22,1);sigma3A<- matrix(sigma,d23,1)
    ########criteria for iterations to stop#######
    L1<- logL(n_samP1,1,1,meanP1,sigma,dataP1)+logL(n_samF1,1,1,meanF1,sigma,dataF1)+logL(n_samP2,1,1,meanP2,sigma,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1A,sigma1A,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2A,sigma2A,dataB2)+logL(n_samF2,d23,mix_pi_3,mean3A,sigma3A,dataF2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*4

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1, 1,0,-1, 0,1,0),3,3)
  b_line1 <- matrix(c(mean1A[1],mean2A))
  B1 <- solve(aa,b_line1)

  jj_1 <- sigmaB1 - sigma1A[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigmaB1

  jj_2 <- sigmaB2 - sigma2A[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigmaB2

  jj_3 <- sigmaF2 - sigma3A[1]
  if(jj_3 < 0) {jj_3 <- 0}
  ll_3 <- jj_3/sigmaF2

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
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

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
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),(1-pkolm(D_F1,n_samF1))))

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
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

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
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),(1-pkolm(D_B1,n_samB1))))

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
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),(1-pkolm(D_B2,n_samB2))))

  ###############  F2  ###################
  dataF2<-sort(dataF2);bmw_F2 <- matrix(0,n_samF2,1); bmwsl_F2 <- matrix(0,n_samF2,d23)
  for(i in 1:d23){
    gg_F2 <- (dataF2 - mean3A[i])/sqrt(sigma3A[i])
    bmw_F2[which(gg_F2>=0)] <- pnorm(gg_F2[gg_F2>=0])
    bmw_F2[which(gg_F2<0)] <- 1 - pnorm(abs(gg_F2[gg_F2<0]))
    bmwsl_F2[,i] <- bmw_F2*mix_pi_3[i]
  }
  P2_F2 <- rowSums(bmwsl_F2)
  nn<-dim(as.matrix(unique(P2_F2)))[1]
  if(nn<n_samF2){P2_F2<-P2_F2+runif(n_samF2)/1e4}

  dd_F2 <- as.matrix(c(sum(P2_F2),sum(P2_F2^2),sum((P2_F2-0.5)^2)))
  WW2_F2 <- 1/(12*n_samF2) + sum((P2_F2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  u_F2 <- as.matrix(c(12*n_samF2*((dd_F2[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((dd_F2[2]/n_samF2-1/3)^2),180*n_samF2*((dd_F2[3]/n_samF2-1/12)^2)))
  D_F2 <- as.numeric(ks.test(P2_F2,"punif")[[1]][1])
  tt_F2 <- as.matrix(c((1 - pchisq(u_F2[1],1)),(1 - pchisq(u_F2[2],1)),(1 - pchisq(u_F2[3],1)),K1(WW2_F2),(1-pkolm(D_F2,n_samF2))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)
  tt_F2[which(tt_F2>=10e-4)]<-round(tt_F2[which(tt_F2>=10e-4)],4);tt_F2[which(tt_F2<10e-4)]<-format(tt_F2[which(tt_F2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("1MG-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigma,4),round(t(mean1A),4)," "," ",round(sigma1A[1],4),round(t(mix_pi_1),4)," "," ",
                       round(t(mean2A),4)," "," ",round(sigma2A[1],4),round(t(mix_pi_2),4)," "," ",round(t(mean3A),4)," "," "," "," "," "," ",round(sigma3A[1],4),round(t(mix_pi_3),4)," "," "," "," "," "," ",
                       round(B1[1],4)," "," "," "," "," ",round(B1[2],4)," ",round(B1[3],4)," "," "," "," "," "," "," ",round(jj_1,4),round(ll_1*100,4)," "," ",
                       round(jj_2,4),round(ll_2*100,4)," "," ",round(jj_3,4),round(ll_3*100,4)," "," ",round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),
                       tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),
                       tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),
                       tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5],round(u_F2[1],4),tt_F2[1],round(u_F2[2],4),
                       tt_F2[2],round(u_F2[3],4),tt_F2[3],round(WW2_F2,4),tt_F2[4],round(D_F2,4),tt_F2[5])
  output<-as.matrix(output)

  OUTPUT<-list(output,mi_1,mi_2,mi_3)
  return(OUTPUT)
}
############################################ 1MG-A Model #########################################
G6ModelFun[[2]] <- function(K1,logL,df11,df21,df31,df41,df51,df61){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1]; n_samF1<-dim(dataF1)[1]; n_samP2<-dim(dataP2)[1]; n_samB1<-dim(dataB1)[1]; n_samB2<-dim(dataB2)[1]; n_samF2<-dim(dataF2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2)
  sumx4<-sum(dataB1);sumx5<-sum(dataB2);sumx6<-sum(dataF2)
  sigmaP1<- as.numeric(var(dataP1)); ss1<- (n_samP1-1)*sigmaP1
  sigmaF1<- as.numeric(var(dataF1)); ss2<- (n_samF1-1)*sigmaF1
  sigmaP2<- as.numeric(var(dataP2)); ss3<- (n_samP2-1)*sigmaP2
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  mean4<- mean(dataB1); sigmaB1<- as.numeric(var(dataB1)); sigma40<- sigmaB1
  mean5<- mean(dataB2); sigmaB2<- as.numeric(var(dataB2)); sigma50<- sigmaB2
  mean6<- mean(dataF2); sigmaF2<- as.numeric(var(dataF2)); sigma60<- sigmaF2

  m_esp<-0.0001
  ####################### 1MG-A Model ################  (A2) ###############
  d21<-2; d22<-2; d23<-3
  sigma<- sigma0
  meanP1<- mean(dataP1);meanF1<- mean(dataF1);meanP2<- mean(dataP2)

  mi_1<- matrix(0.5,d21,1)
  sigma1A<- matrix(sigma,d21,1);mean1A<- as.matrix(c(meanP1,meanF1))

  mi_2<- matrix(0.5,d22,1)
  sigma2A<- matrix(sigma,d22,1);mean2A<- as.matrix(c(meanF1,meanP2))

  mi_3<- as.matrix(c(0.25,0.5,0.25))
  sigma3A<- matrix(sigma,d23,1);mean3A<- as.matrix(c(meanP1,meanF1,meanP2))

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1A,sigma1A,dataB1)+logL(n_samB2,d22,mi_2,mean2A,sigma2A,dataB2)+logL(n_samF2,d23,mi_3,mean3A,sigma3A,dataF2)

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)
  WW_F2 <- matrix(0,d23,n_samF2); swx_F2 <- matrix(0,d23,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1A[i],sqrt(sigma1A[i]))/dmixnorm(dataB1,mean1A,sqrt(sigma1A),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2A[i],sqrt(sigma2A[i]))/dmixnorm(dataB2,mean2A,sqrt(sigma2A),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    for(i in 1:d23) { WW_F2[i,] <- mi_3[i]*dnorm(dataF2,mean3A[i],sqrt(sigma3A[i]))/dmixnorm(dataF2,mean3A,sqrt(sigma3A),mi_3) }
    mix_pi_3 <- as.matrix(rowSums(WW_F2)/n_samF2)
    sumwx_F2 <- WW_F2%*%dataF2

    ############ CM1-step for means ##############
    aaa0<- 0;aa1<- 1000
    s0<- matrix(0,6,1)
    s0[1]<- sumx1+sumwx_B1[1]+sumwx_F2[1]
    s0[2]<- sumx2+sumwx_B1[2]+sumwx_B2[1]+sumwx_F2[2]
    s0[3]<- sumx3+sumwx_B2[2]+sumwx_F2[3]
    n0<- matrix(0,6,1)
    n0[1]<- n_samP1+mix_pi_1[1]*n_samB1+mix_pi_3[1]*n_samF2
    n0[2]<- n_samF1+mix_pi_1[2]*n_samB1+mix_pi_2[1]*n_samB2+mix_pi_3[2]*n_samF2
    n0[3]<- n_samP2+mix_pi_2[2]*n_samB2+mix_pi_3[3]*n_samF2
    while (aa1>0.0001){
      aa3<- s0[1]/n0[1]-2.0*s0[2]/n0[2]+s0[3]/n0[3]
      aa4<- sigma*(1.0/n0[1]+4.0/n0[2]+1.0/n0[3])
      aaa1<- aa3/aa4
      meanP1<- (s0[1]-aaa1*sigma)/n0[1]
      meanF1<- (s0[2]+2.0*aaa1*sigma)/n0[2]
      meanP2<- (s0[3]-aaa1*sigma)/n0[3]
      aa1<- abs(aaa1-aaa0)
      aaa0<- aaa1
    }
    mean1A<- as.matrix(c(meanP1,meanF1));mean2A<- as.matrix(c(meanF1,meanP2));mean3A<- as.matrix(c(meanP1,meanF1,meanP2))

    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2)
    ss2<- sum((dataF1-meanF1)^2)
    ss3<- sum((dataP2-meanP2)^2)

    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1A[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2A[i])^2 }
    for(i in 1:d23) {  swx_F2[i] <- WW_F2[i,]%*%(dataF2-mean3A[i])^2 }

    s0[5]<- ss1+ss2+ss3+sum(swx_B1)+sum(swx_B2)+sum(swx_F2)
    n0[5]<- n_samP1+n_samF1+n_samP2+n_samB1+n_samB2+n_samF2
    sigma<- s0[5]/n0[5]
    sigma1A<- matrix(sigma,d21,1);sigma2A<- matrix(sigma,d22,1);sigma3A<- matrix(sigma,d23,1)
    ########criteria for iterations to stop#######
    L1<- logL(n_samP1,1,1,meanP1,sigma,dataP1)+logL(n_samF1,1,1,meanF1,sigma,dataF1)+logL(n_samP2,1,1,meanP2,sigma,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1A,sigma1A,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2A,sigma2A,dataB2)+logL(n_samF2,d23,mix_pi_3,mean3A,sigma3A,dataF2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*3

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1, 1,0,-1),3,2)
  b_line1 <- matrix(c(mean1A[1],mean2A))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigmaB1 - sigma1A[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigmaB1

  jj_2 <- sigmaB2 - sigma2A[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigmaB2

  jj_3 <- sigmaF2 - sigma3A[1]
  if(jj_3 < 0) {jj_3 <- 0}
  ll_3 <- jj_3/sigmaF2

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
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

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
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),(1-pkolm(D_F1,n_samF1))))

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
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

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
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),(1-pkolm(D_B1,n_samB1))))

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
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),(1-pkolm(D_B2,n_samB2))))

  ###############  F2  ###################
  dataF2<-sort(dataF2);bmw_F2 <- matrix(0,n_samF2,1); bmwsl_F2 <- matrix(0,n_samF2,d23)
  for(i in 1:d23){
    gg_F2 <- (dataF2 - mean3A[i])/sqrt(sigma3A[i])
    bmw_F2[which(gg_F2>=0)] <- pnorm(gg_F2[gg_F2>=0])
    bmw_F2[which(gg_F2<0)] <- 1 - pnorm(abs(gg_F2[gg_F2<0]))
    bmwsl_F2[,i] <- bmw_F2*mix_pi_3[i]
  }
  P2_F2 <- rowSums(bmwsl_F2)
  nn<-dim(as.matrix(unique(P2_F2)))[1]
  if(nn<n_samF2){P2_F2<-P2_F2+runif(n_samF2)/1e4}

  dd_F2 <- as.matrix(c(sum(P2_F2),sum(P2_F2^2),sum((P2_F2-0.5)^2)))
  WW2_F2 <- 1/(12*n_samF2) + sum((P2_F2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  u_F2 <- as.matrix(c(12*n_samF2*((dd_F2[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((dd_F2[2]/n_samF2-1/3)^2),180*n_samF2*((dd_F2[3]/n_samF2-1/12)^2)))
  D_F2 <- as.numeric(ks.test(P2_F2,"punif")[[1]][1])
  tt_F2 <- as.matrix(c((1 - pchisq(u_F2[1],1)),(1 - pchisq(u_F2[2],1)),(1 - pchisq(u_F2[3],1)),K1(WW2_F2),(1-pkolm(D_F2,n_samF2))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)
  tt_F2[which(tt_F2>=10e-4)]<-round(tt_F2[which(tt_F2>=10e-4)],4);tt_F2[which(tt_F2<10e-4)]<-format(tt_F2[which(tt_F2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("1MG-A",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigma,4),round(t(mean1A),4)," "," ",round(sigma1A[1],4),round(t(mix_pi_1),4)," "," ",
                       round(t(mean2A),4)," "," ",round(sigma2A[1],4),round(t(mix_pi_2),4)," "," ",round(t(mean3A),4)," "," "," "," "," "," ",round(sigma3A[1],4),round(t(mix_pi_3),4)," "," "," "," "," "," ",
                       round(B1[1],4)," "," "," "," "," ",round(B1[2],4)," "," "," "," "," "," "," "," "," ",round(jj_1,4),round(ll_1*100,4)," "," ",
                       round(jj_2,4),round(ll_2*100,4)," "," ",round(jj_3,4),round(ll_3*100,4)," "," ",round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),
                       tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),
                       tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),
                       tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5],round(u_F2[1],4),tt_F2[1],round(u_F2[2],4),
                       tt_F2[2],round(u_F2[3],4),tt_F2[3],round(WW2_F2,4),tt_F2[4],round(D_F2,4),tt_F2[5])

  output<-as.matrix(output)

  OUTPUT<-list(output,mi_1,mi_2,mi_3)
  return(OUTPUT)
}
############################################ 1MG-EAD Model #########################################
G6ModelFun[[3]] <- function(K1,logL,df11,df21,df31,df41,df51,df61){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1]; n_samF1<-dim(dataF1)[1]; n_samP2<-dim(dataP2)[1]; n_samB1<-dim(dataB1)[1]; n_samB2<-dim(dataB2)[1]; n_samF2<-dim(dataF2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2)
  sumx4<-sum(dataB1);sumx5<-sum(dataB2);sumx6<-sum(dataF2)
  sigmaP1<- as.numeric(var(dataP1)); ss1<- (n_samP1-1)*sigmaP1
  sigmaF1<- as.numeric(var(dataF1)); ss2<- (n_samF1-1)*sigmaF1
  sigmaP2<- as.numeric(var(dataP2)); ss3<- (n_samP2-1)*sigmaP2
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  mean4<- mean(dataB1); sigmaB1<- as.numeric(var(dataB1)); sigma40<- sigmaB1
  mean5<- mean(dataB2); sigmaB2<- as.numeric(var(dataB2)); sigma50<- sigmaB2
  mean6<- mean(dataF2); sigmaF2<- as.numeric(var(dataF2)); sigma60<- sigmaF2

  m_esp<-0.0001
  ####################### 1MG-EAD Model ######################  (A3) ########### A3
  d21<-1; d22<-2; d23<-2
  sigma<- sigma0
  meanP1<- mean(dataP1);meanF1<- mean(dataF1);meanP2<- mean(dataP2)
  meanF1<- (meanP1+meanF1)/2; meanP1<- meanF1

  mi_1<- as.matrix(1)
  sigma1A<- as.matrix(sigma) ;mean1A<- as.matrix(meanP1)

  mi_2<- matrix(0.5,d22,1)
  sigma2A<- matrix(sigma,d22,1);mean2A<- as.matrix(c(meanP1,meanP2))

  mi_3<- as.matrix(c(0.75,0.25))
  sigma3A<- matrix(sigma,d23,1);mean3A<- as.matrix(c(meanP1,meanP2))

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1A,sigma1A,dataB1)+logL(n_samB2,d22,mi_2,mean2A,sigma2A,dataB2)+logL(n_samF2,d23,mi_3,mean3A,sigma3A,dataF2)

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)
  WW_F2 <- matrix(0,d23,n_samF2); swx_F2 <- matrix(0,d23,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1A[i],sqrt(sigma1A[i]))/dmixnorm(dataB1,mean1A,sqrt(sigma1A),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2A[i],sqrt(sigma2A[i]))/dmixnorm(dataB2,mean2A,sqrt(sigma2A),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    for(i in 1:d23) { WW_F2[i,] <- mi_3[i]*dnorm(dataF2,mean3A[i],sqrt(sigma3A[i]))/dmixnorm(dataF2,mean3A,sqrt(sigma3A),mi_3) }
    mix_pi_3 <- as.matrix(rowSums(WW_F2)/n_samF2)
    sumwx_F2 <- WW_F2%*%dataF2

    ############ CM1-step for means ##############
    aaa0<- 0
    s0<- matrix(0,6,1)
    s0[1]<- sumx1+sumx2+sumx4+sumwx_B2[1]+sumwx_F2[1]
    s0[3]<- sumx3+sumwx_B2[2]+sumwx_F2[2]
    n0<- matrix(0,6,1)
    n0[1]<- n_samP1+n_samF1+n_samB1+mix_pi_2[1]*n_samB2+mix_pi_3[1]*n_samF2
    n0[3]<- n_samP2+mix_pi_2[2]*n_samB2+mix_pi_3[2]*n_samF2
    meanP1<- s0[1]/n0[1]
    meanP2<- s0[3]/n0[3]
    meanF1<- meanP1
    mean1A[1]<- meanP1; mean4<- meanP1
    mean2A<- as.matrix(c(meanP1,meanP2))
    mean3A<- mean2A

    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2)
    ss2<- sum((dataF1-meanF1)^2)
    ss3<- sum((dataP2-meanP2)^2)
    ss4<- sum((dataB1-mean4)^2)

    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2A[i])^2 }
    for(i in 1:d23) {  swx_F2[i] <- WW_F2[i,]%*%(dataF2-mean3A[i])^2 }

    s0[4]<- ss1+ss2+ss3+ss4+sum(swx_B2)+sum(swx_F2)
    n0[4]<- n_samP1+n_samF1+n_samP2+n_samB1+n_samB2+n_samF2
    sigma<- s0[4]/n0[4]
    sigma1A<- matrix(sigma,d21,1);sigma2A<- matrix(sigma,d22,1);sigma3A<- matrix(sigma,d23,1)
    ########criteria for iterations to stop#######
    L1<- logL(n_samP1,1,1,meanP1,sigma,dataP1)+logL(n_samF1,1,1,meanF1,sigma,dataF1)+logL(n_samP2,1,1,meanP2,sigma,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1A,sigma1A,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2A,sigma2A,dataB2)+logL(n_samF2,d23,mix_pi_3,mean3A,sigma3A,dataF2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*3

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1, 1,-1),2,2)
  b_line1 <- matrix(c(mean2A))
  B1 <- solve(aa,b_line1)

  jj_1 <- sigmaB1 - sigma1A[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigmaB1

  jj_2 <- sigmaB2 - sigma2A[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigmaB2

  jj_3 <- sigmaF2 - sigma3A[1]
  if(jj_3 < 0) {jj_3 <- 0}
  ll_3 <- jj_3/sigmaF2

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
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

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
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),(1-pkolm(D_F1,n_samF1))))

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
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

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
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),(1-pkolm(D_B1,n_samB1))))

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
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),(1-pkolm(D_B2,n_samB2))))

  ###############  F2  ###################
  dataF2<-sort(dataF2);bmw_F2 <- matrix(0,n_samF2,1); bmwsl_F2 <- matrix(0,n_samF2,d23)
  for(i in 1:d23){
    gg_F2 <- (dataF2 - mean3A[i])/sqrt(sigma3A[i])
    bmw_F2[which(gg_F2>=0)] <- pnorm(gg_F2[gg_F2>=0])
    bmw_F2[which(gg_F2<0)] <- 1 - pnorm(abs(gg_F2[gg_F2<0]))
    bmwsl_F2[,i] <- bmw_F2*mix_pi_3[i]
  }
  P2_F2 <- rowSums(bmwsl_F2)
  nn<-dim(as.matrix(unique(P2_F2)))[1]
  if(nn<n_samF2){P2_F2<-P2_F2+runif(n_samF2)/1e4}

  dd_F2 <- as.matrix(c(sum(P2_F2),sum(P2_F2^2),sum((P2_F2-0.5)^2)))
  WW2_F2 <- 1/(12*n_samF2) + sum((P2_F2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  u_F2 <- as.matrix(c(12*n_samF2*((dd_F2[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((dd_F2[2]/n_samF2-1/3)^2),180*n_samF2*((dd_F2[3]/n_samF2-1/12)^2)))
  D_F2 <- as.numeric(ks.test(P2_F2,"punif")[[1]][1])
  tt_F2 <- as.matrix(c((1 - pchisq(u_F2[1],1)),(1 - pchisq(u_F2[2],1)),(1 - pchisq(u_F2[3],1)),K1(WW2_F2),(1-pkolm(D_F2,n_samF2))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)
  tt_F2[which(tt_F2>=10e-4)]<-round(tt_F2[which(tt_F2>=10e-4)],4);tt_F2[which(tt_F2<10e-4)]<-format(tt_F2[which(tt_F2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("1MG-EAD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigma,4),round(t(mean1A),4)," "," "," ",round(sigma1A[1],4),round(t(mix_pi_1),4)," "," "," ",
                       round(t(mean2A),4)," "," ",round(sigma2A[1],4),round(t(mix_pi_2),4)," "," ",round(t(mean3A),4)," "," "," "," "," "," "," ",round(sigma3A[1],4),round(t(mix_pi_3),4)," "," "," "," "," "," "," ",
                       round(B1[1],4)," "," "," "," "," ",round(B1[2],4)," ",round(B1[2],4)," "," "," "," "," "," "," ",round(jj_1,4),round(ll_1*100,4)," "," ",
                       round(jj_2,4),round(ll_2*100,4)," "," ",round(jj_3,4),round(ll_3*100,4)," "," ",round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),
                       tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),
                       tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),
                       tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5],round(u_F2[1],4),tt_F2[1],round(u_F2[2],4),
                       tt_F2[2],round(u_F2[3],4),tt_F2[3],round(WW2_F2,4),tt_F2[4],round(D_F2,4),tt_F2[5])

  output<-as.matrix(output)


  OUTPUT<-list(output,mi_1,mi_2,mi_3)
  return(OUTPUT)
}
############################################ 1MG-NCD Model #########################################
G6ModelFun[[4]] <- function(K1,logL,df11,df21,df31,df41,df51,df61){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1]; n_samF1<-dim(dataF1)[1]; n_samP2<-dim(dataP2)[1]; n_samB1<-dim(dataB1)[1]; n_samB2<-dim(dataB2)[1]; n_samF2<-dim(dataF2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2)
  sumx4<-sum(dataB1);sumx5<-sum(dataB2);sumx6<-sum(dataF2)
  sigmaP1<- as.numeric(var(dataP1)); ss1<- (n_samP1-1)*sigmaP1
  sigmaF1<- as.numeric(var(dataF1)); ss2<- (n_samF1-1)*sigmaF1
  sigmaP2<- as.numeric(var(dataP2)); ss3<- (n_samP2-1)*sigmaP2
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  mean4<- mean(dataB1); sigmaB1<- as.numeric(var(dataB1)); sigma40<- sigmaB1
  mean5<- mean(dataB2); sigmaB2<- as.numeric(var(dataB2)); sigma50<- sigmaB2
  mean6<- mean(dataF2); sigmaF2<- as.numeric(var(dataF2)); sigma60<- sigmaF2

  m_esp<-0.0001
  ####################### 1MG-NCD Model ########################  (A4) #########
  d21<-2; d22<-1; d23<-2
  sigma<- sigma0
  meanP1<- mean(dataP1);meanF1<- mean(dataF1);meanP2<- mean(dataP2)
  meanP2<- (meanP2+meanF1)/2; meanF1<- meanP2

  mi_1<- matrix(0.5,d21,1)
  sigma1A<- matrix(sigma,d21,1);mean1A<- as.matrix(c(meanP1,meanP2))

  mi_2<- as.matrix(1)
  sigma2A<- as.matrix(sigma) ;mean2A<- as.matrix(meanP2)

  mi_3<- as.matrix(c(0.25,0.75))
  sigma3A<- matrix(sigma,d23,1);mean3A<- as.matrix(c(meanP1,meanP2))

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1A,sigma1A,dataB1)+logL(n_samB2,d22,mi_2,mean2A,sigma2A,dataB2)+logL(n_samF2,d23,mi_3,mean3A,sigma3A,dataF2)

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)
  WW_F2 <- matrix(0,d23,n_samF2); swx_F2 <- matrix(0,d23,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1A[i],sqrt(sigma1A[i]))/dmixnorm(dataB1,mean1A,sqrt(sigma1A),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2A[i],sqrt(sigma2A[i]))/dmixnorm(dataB2,mean2A,sqrt(sigma2A),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    for(i in 1:d23) { WW_F2[i,] <- mi_3[i]*dnorm(dataF2,mean3A[i],sqrt(sigma3A[i]))/dmixnorm(dataF2,mean3A,sqrt(sigma3A),mi_3) }
    mix_pi_3 <- as.matrix(rowSums(WW_F2)/n_samF2)
    sumwx_F2 <- WW_F2%*%dataF2

    ############ CM1-step for means ##############
    aaa0<- 0;aa1<- 1000
    s0<- matrix(0,6,1)
    s0[1]<- sumx1+sumwx_B1[1]+sumwx_F2[1]
    s0[3]<- sumx2+sumx3+sumwx_B1[2]+sumx5+sumwx_F2[2]
    n0<- matrix(0,6,1)
    n0[1]<- n_samP1+mix_pi_1[1]*n_samB1+mix_pi_3[1]*n_samF2
    n0[3]<- n_samF1+n_samP2+mix_pi_1[2]*n_samB1+n_samB2+mix_pi_3[2]*n_samF2
    meanP1<- s0[1]/n0[1]; meanP2<- s0[3]/n0[3]
    mean1A<- as.matrix(c(meanP1,meanP2));mean2A<- meanP2; mean5<- meanP2;mean3A<- as.matrix(c(meanP1,meanP2))

    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2)
    ss2<- sum((dataF1-meanF1)^2)
    ss3<- sum((dataP2-meanP2)^2)
    ss5<- sum((dataB2-mean5)^2)

    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1A[i])^2 }
    for(i in 1:d23) {  swx_F2[i] <- WW_F2[i,]%*%(dataF2-mean3A[i])^2 }

    s0[4]<- ss1+ss2+ss3+ss5+sum(swx_B1)+sum(swx_F2)
    n0[4]<- n_samP1+n_samF1+n_samP2+n_samB1+n_samB2+n_samF2
    sigma<- s0[4]/n0[4]
    sigma1A<- matrix(sigma,d21,1);sigma2A<- sigma ;sigma3A<- matrix(sigma,d23,1)
    ########criteria for iterations to stop#######
    L1<- logL(n_samP1,1,1,meanP1,sigma,dataP1)+logL(n_samF1,1,1,meanF1,sigma,dataF1)+logL(n_samP2,1,1,meanP2,sigma,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1A,sigma1A,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2A,sigma2A,dataB2)+logL(n_samF2,d23,mix_pi_3,mean3A,sigma3A,dataF2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*3

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1, 1,-1),2,2)
  b_line1 <- matrix(c(mean3A))
  B1 <- solve(aa,b_line1)

  jj_1 <- sigmaB1 - sigma1A[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigmaB1

  jj_2 <- sigmaB2 - sigma2A[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigmaB2

  jj_3 <- sigmaF2 - sigma3A[1]
  if(jj_3 < 0) {jj_3 <- 0}
  ll_3 <- jj_3/sigmaF2

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
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

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
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),(1-pkolm(D_F1,n_samF1))))

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
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

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
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),(1-pkolm(D_B1,n_samB1))))

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
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),(1-pkolm(D_B2,n_samB2))))

  ###############  F2  ###################
  dataF2<-sort(dataF2);bmw_F2 <- matrix(0,n_samF2,1); bmwsl_F2 <- matrix(0,n_samF2,d23)
  for(i in 1:d23){
    gg_F2 <- (dataF2 - mean3A[i])/sqrt(sigma3A[i])
    bmw_F2[which(gg_F2>=0)] <- pnorm(gg_F2[gg_F2>=0])
    bmw_F2[which(gg_F2<0)] <- 1 - pnorm(abs(gg_F2[gg_F2<0]))
    bmwsl_F2[,i] <- bmw_F2*mix_pi_3[i]
  }
  P2_F2 <- rowSums(bmwsl_F2)
  nn<-dim(as.matrix(unique(P2_F2)))[1]
  if(nn<n_samF2){P2_F2<-P2_F2+runif(n_samF2)/1e4}

  dd_F2 <- as.matrix(c(sum(P2_F2),sum(P2_F2^2),sum((P2_F2-0.5)^2)))
  WW2_F2 <- 1/(12*n_samF2) + sum((P2_F2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  u_F2 <- as.matrix(c(12*n_samF2*((dd_F2[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((dd_F2[2]/n_samF2-1/3)^2),180*n_samF2*((dd_F2[3]/n_samF2-1/12)^2)))
  D_F2 <- as.numeric(ks.test(P2_F2,"punif")[[1]][1])
  tt_F2 <- as.matrix(c((1 - pchisq(u_F2[1],1)),(1 - pchisq(u_F2[2],1)),(1 - pchisq(u_F2[3],1)),K1(WW2_F2),(1-pkolm(D_F2,n_samF2))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)
  tt_F2[which(tt_F2>=10e-4)]<-round(tt_F2[which(tt_F2>=10e-4)],4);tt_F2[which(tt_F2<10e-4)]<-format(tt_F2[which(tt_F2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("1MG-NCD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigma,4),round(t(mean1A),4)," "," ",round(sigma1A[1],4),round(t(mix_pi_1),4)," "," ",
                       round(t(mean2A),4)," "," "," ",round(sigma2A[1],4),round(t(mix_pi_2),4)," "," "," ",round(t(mean3A),4)," "," "," "," "," "," "," ",round(sigma3A[1],4),round(t(mix_pi_3),4)," "," "," "," "," "," "," ",
                       round(B1[1],4)," "," "," "," "," ",round(B1[2],4)," ",round(-B1[2],4)," "," "," "," "," "," "," ",round(jj_1,4),round(ll_1*100,4)," "," ",
                       round(jj_2,4),round(ll_2*100,4)," "," ",round(jj_3,4),round(ll_3*100,4)," "," ",round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),
                       tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),
                       tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),
                       tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5],round(u_F2[1],4),tt_F2[1],round(u_F2[2],4),
                       tt_F2[2],round(u_F2[3],4),tt_F2[3],round(WW2_F2,4),tt_F2[4],round(D_F2,4),tt_F2[5])

  output<-as.matrix(output)



  OUTPUT<-list(output,mi_1,mi_2,mi_3)
  return(OUTPUT)
}
############################################ 2MG-ADI Model #########################################
G6ModelFun[[5]] <- function(K1,logL,df11,df21,df31,df41,df51,df61){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1]; n_samF1<-dim(dataF1)[1]; n_samP2<-dim(dataP2)[1]; n_samB1<-dim(dataB1)[1]; n_samB2<-dim(dataB2)[1]; n_samF2<-dim(dataF2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2)
  sumx4<-sum(dataB1);sumx5<-sum(dataB2);sumx6<-sum(dataF2)
  sigmaP1<- as.numeric(var(dataP1)); ss1<- (n_samP1-1)*sigmaP1
  sigmaF1<- as.numeric(var(dataF1)); ss2<- (n_samF1-1)*sigmaF1
  sigmaP2<- as.numeric(var(dataP2)); ss3<- (n_samP2-1)*sigmaP2
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  mean4<- mean(dataB1); sigmaB1<- as.numeric(var(dataB1)); sigma40<- sigmaB1
  mean5<- mean(dataB2); sigmaB2<- as.numeric(var(dataB2)); sigma50<- sigmaB2
  mean6<- mean(dataF2); sigmaF2<- as.numeric(var(dataF2)); sigma60<- sigmaF2

  m_esp<-0.0001
  ####################### 2MG-ADI Model ######################  (B1) ##########
  d21<-4; d22<-4; d23<-9
  sigma<- sigma0
  meanP1<- mean(dataP1);meanF1<- mean(dataF1);meanP2<- mean(dataP2)

  mi_1<- matrix(0.25,d21,1)
  a1<- sqrt(sigmaB1/(n_samB1-1))
  if (meanP1<meanP2) {a1= -a1}
  mean1B<- as.matrix(c(meanP1,mean4+0.5*a1,mean4-0.5*a1,meanF1))

  mi_2<- matrix(0.25,d22,1)
  a2<- sqrt(sigmaB2/(n_samB2-1))
  if (meanP1<meanP2) { a2=-a2 }
  mean2B<- as.matrix(c(meanF1,mean5+0.5*a2,mean5-0.5*a2,meanP2))

  mi_3<- as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  mean3B<- as.matrix(c(meanP1,mean1B[2],0.5*(meanP1+mean1B[2]),mean1B[3],meanF1,mean2B[2],0.5*(mean2B[3]+meanP2),mean2B[3],meanP2))

  sigma1B<- matrix(sigma,d21,1)
  sigma2B<- matrix(sigma,d22,1)
  sigma3B<- matrix(sigma,d23,1)

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1B,sigma1B,dataB1)+logL(n_samB2,d22,mi_2,mean2B,sigma2B,dataB2)+logL(n_samF2,d23,mi_3,mean3B,sigma3B,dataF2)

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)
  WW_F2 <- matrix(0,d23,n_samF2); swx_F2 <- matrix(0,d23,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1B[i],sqrt(sigma1B[i]))/dmixnorm(dataB1,mean1B,sqrt(sigma1B),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2B[i],sqrt(sigma2B[i]))/dmixnorm(dataB2,mean2B,sqrt(sigma2B),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    for(i in 1:d23) { WW_F2[i,] <- mi_3[i]*dnorm(dataF2,mean3B[i],sqrt(sigma3B[i]))/dmixnorm(dataF2,mean3B,sqrt(sigma3B),mi_3) }
    mix_pi_3 <- as.matrix(rowSums(WW_F2)/n_samF2)
    sumwx_F2 <- WW_F2%*%dataF2

    ############ CM1-step for means ##############
    s0<- matrix(0,12,1)
    s0[1]<- sumx1+sumwx_B1[1]+sumwx_F2[1]
    s0[2]<- sumx2+sumwx_B1[4]+sumwx_B2[1]+sumwx_F2[5]
    s0[3]<- sumx3+sumwx_B2[4]+sumwx_F2[9]
    s0[4]<- sumwx_B1[2]+sumwx_F2[2]
    s0[5]<- sumwx_B1[3]+sumwx_F2[4]
    s0[6]<- sumwx_B2[2]+sumwx_F2[6]
    s0[7]<- sumwx_B2[3]+sumwx_F2[8]
    s0[8]<- sumwx_F2[3]
    s0[9]<- sumwx_F2[7]

    n0<- matrix(0,12,1)
    n0[1]<- n_samP1+mix_pi_1[1]*n_samB1+mix_pi_3[1]*n_samF2
    n0[2]<- n_samF1+mix_pi_1[4]*n_samB1+mix_pi_2[1]*n_samB2+mix_pi_3[5]*n_samF2
    n0[3]<- n_samP2+mix_pi_2[4]*n_samB2+mix_pi_3[9]*n_samF2
    n0[4]<- mix_pi_1[2]*n_samB1+mix_pi_3[2]*n_samF2
    n0[5]<- mix_pi_1[3]*n_samB1+mix_pi_3[4]*n_samF2
    n0[6]<- mix_pi_2[2]*n_samB2+mix_pi_3[6]*n_samF2
    n0[7]<- mix_pi_2[3]*n_samB2+mix_pi_3[8]*n_samF2
    n0[8]<- mix_pi_3[3]*n_samF2
    n0[9]<- mix_pi_3[7]*n_samF2

    meanP1<- s0[1]/n0[1]
    meanF1<- s0[2]/n0[2]
    meanP2<- s0[3]/n0[3]
    mean1B[2]<- s0[4]/n0[4];mean1B[3]<- s0[5]/n0[5]
    mean2B[2]<- s0[6]/n0[6];mean2B[3]<- s0[7]/n0[7]
    mean3B[3]<- s0[8]/n0[8];mean3B[7]<- s0[9]/n0[9]
    mean1B[1]<- meanP1;mean3B[1]<- meanP1
    mean1B[4]<- meanF1;mean2B[1]<- meanF1;mean3B[5]<- meanF1
    mean2B[4]<- meanP2;mean3B[9]<- meanP2
    mean3B[2]<- mean1B[2];mean3B[4]<- mean1B[3]
    mean3B[6]<- mean2B[2];mean3B[8]<- mean2B[3]

    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2)
    ss2<- sum((dataF1-meanF1)^2)
    ss3<- sum((dataP2-meanP2)^2)


    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1B[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2B[i])^2 }
    for(i in 1:d23) {  swx_F2[i] <- WW_F2[i,]%*%(dataF2-mean3B[i])^2 }

    s0[11]<- ss1+ss2+ss3+sum(swx_B1)+sum(swx_B2)+sum(swx_F2)
    n0[11]<- n_samP1+n_samF1+n_samP2+n_samB1+n_samB2+n_samF2
    sigma<- s0[11]/n0[11]
    sigma1B<- matrix(sigma,d21,1);sigma2B<- matrix(sigma,d22,1);sigma3B<- matrix(sigma,d23,1)
    ########criteria for iterations to stop#######
    L1<- logL(n_samP1,1,1,meanP1,sigma,dataP1)+logL(n_samF1,1,1,meanF1,sigma,dataF1)+logL(n_samP2,1,1,meanP2,sigma,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1B,sigma1B,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2B,sigma2B,dataB2)+logL(n_samF2,d23,mix_pi_3,mean3B,sigma3B,dataF2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*10

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1,1,1, 1,0,-1,1,0,0,-1,1,-1, 1,0,-1,0,1,-1,0,-1,1 ,0,1,0,0,1,1,0,0,0,
                0,1,0,1,0,0,1,0,0, 1,0,1,0,0,0,0,-1,-1, 0,0,0,1,0,0,-1,0,0, 0,0,0,0,1,-1,0,0,0, 0,1,0,0,0,0,0,0,0),9,9)
  b_line1 <- matrix(c(meanP1,meanF1,meanP2,mean1B[2],mean1B[3],mean2B[2],mean2B[3],mean3B[3],mean3B[7]))
  B1 <- solve(aa,b_line1)

  jj_1 <- sigmaB1 - sigma1B[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigmaB1

  jj_2 <- sigmaB2 - sigma2B[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigmaB2

  jj_3 <- sigmaF2 - sigma3B[1]
  if(jj_3 < 0) {jj_3 <- 0}
  ll_3 <- jj_3/sigmaF2

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
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

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
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),(1-pkolm(D_F1,n_samF1))))

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
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

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
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),(1-pkolm(D_B1,n_samB1))))

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
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),(1-pkolm(D_B2,n_samB2))))

  ###############  F2  ###################
  dataF2<-sort(dataF2);bmw_F2 <- matrix(0,n_samF2,1); bmwsl_F2 <- matrix(0,n_samF2,d23)
  for(i in 1:d23){
    gg_F2 <- (dataF2 - mean3B[i])/sqrt(sigma3B[i])
    bmw_F2[which(gg_F2>=0)] <- pnorm(gg_F2[gg_F2>=0])
    bmw_F2[which(gg_F2<0)] <- 1 - pnorm(abs(gg_F2[gg_F2<0]))
    bmwsl_F2[,i] <- bmw_F2*mix_pi_3[i]
  }
  P2_F2 <- rowSums(bmwsl_F2)
  nn<-dim(as.matrix(unique(P2_F2)))[1]
  if(nn<n_samF2){P2_F2<-P2_F2+runif(n_samF2)/1e4}

  dd_F2 <- as.matrix(c(sum(P2_F2),sum(P2_F2^2),sum((P2_F2-0.5)^2)))
  WW2_F2 <- 1/(12*n_samF2) + sum((P2_F2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  u_F2 <- as.matrix(c(12*n_samF2*((dd_F2[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((dd_F2[2]/n_samF2-1/3)^2),180*n_samF2*((dd_F2[3]/n_samF2-1/12)^2)))
  D_F2 <- as.numeric(ks.test(P2_F2,"punif")[[1]][1])
  tt_F2 <- as.matrix(c((1 - pchisq(u_F2[1],1)),(1 - pchisq(u_F2[2],1)),(1 - pchisq(u_F2[3],1)),K1(WW2_F2),(1-pkolm(D_F2,n_samF2))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)
  tt_F2[which(tt_F2>=10e-4)]<-round(tt_F2[which(tt_F2>=10e-4)],4);tt_F2[which(tt_F2<10e-4)]<-format(tt_F2[which(tt_F2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("2MG-ADI",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigma,4),round(t(mean1B),4),round(sigma1B[1],4),round(t(mix_pi_1),4),
                       round(t(mean2B),4),round(sigma2B[1],4),round(t(mix_pi_2),4),round(t(mean3B),4),round(sigma3B[1],4),round(t(mix_pi_3),4),
                       round(B1[1],4)," "," "," "," "," ",round(B1[2],4),round(B1[3],4),round(B1[4],4),round(B1[5],4),round(B1[6],4),round(B1[7],4),round(B1[8],4),round(B1[9],4)," "," ",
                       round(jj_1,4),round(ll_1*100,4)," "," ",round(jj_2,4),round(ll_2*100,4)," "," ",round(jj_3,4),round(ll_3*100,4)," "," ",round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),
                       tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),
                       tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),
                       tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5],round(u_F2[1],4),tt_F2[1],round(u_F2[2],4),
                       tt_F2[2],round(u_F2[3],4),tt_F2[3],round(WW2_F2,4),tt_F2[4],round(D_F2,4),tt_F2[5])

  output<-as.matrix(output)




  OUTPUT<-list(output,mi_1,mi_2,mi_3)
  return(OUTPUT)
}
############################################ 2MG-AD Model #########################################
G6ModelFun[[6]] <- function(K1,logL,df11,df21,df31,df41,df51,df61){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1]; n_samF1<-dim(dataF1)[1]; n_samP2<-dim(dataP2)[1]; n_samB1<-dim(dataB1)[1]; n_samB2<-dim(dataB2)[1]; n_samF2<-dim(dataF2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2)
  sumx4<-sum(dataB1);sumx5<-sum(dataB2);sumx6<-sum(dataF2)
  sigmaP1<- as.numeric(var(dataP1)); ss1<- (n_samP1-1)*sigmaP1
  sigmaF1<- as.numeric(var(dataF1)); ss2<- (n_samF1-1)*sigmaF1
  sigmaP2<- as.numeric(var(dataP2)); ss3<- (n_samP2-1)*sigmaP2
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  mean4<- mean(dataB1); sigmaB1<- as.numeric(var(dataB1)); sigma40<- sigmaB1
  mean5<- mean(dataB2); sigmaB2<- as.numeric(var(dataB2)); sigma50<- sigmaB2
  mean6<- mean(dataF2); sigmaF2<- as.numeric(var(dataF2)); sigma60<- sigmaF2

  m_esp<-0.0001
  ####################### 2MG-AD Model ###################### (B2) ############
  d21<-4; d22<-4; d23<-9
  sigma<- sigma0
  meanP1<- mean(dataP1);meanF1<- mean(dataF1);meanP2<- mean(dataP2)

  mi_1<- matrix(0.25,d21,1)
  a1<- sqrt(sigmaB1/(n_samB1-1))
  if (meanP1<meanP2) {a1= -a1}
  mean1B<- as.matrix(c(meanP1,mean4+0.5*a1,mean4-0.5*a1,meanF1))

  mi_2<- matrix(0.25,d22,1)
  a2<- sqrt(sigmaB2/(n_samB2-1))
  if (meanP1<meanP2) { a2= -a2 }
  mean2B<- as.matrix(c(meanF1,mean5+0.5*a2,mean5-0.5*a2,meanP2))

  mi_3<- as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  mean3B<- as.matrix(c(meanP1,mean1B[2],0.5*(meanP1+mean1B[2]),mean1B[3],meanF1,mean2B[2],0.5*(mean2B[3]+meanP2),mean2B[3],meanP2))

  sigma1B<- matrix(sigma,d21,1)
  sigma2B<- matrix(sigma,d22,1)
  sigma3B<- matrix(sigma,d23,1)

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1B,sigma1B,dataB1)+logL(n_samB2,d22,mi_2,mean2B,sigma2B,dataB2)+logL(n_samF2,d23,mi_3,mean3B,sigma3B,dataF2)

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)
  WW_F2 <- matrix(0,d23,n_samF2); swx_F2 <- matrix(0,d23,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1B[i],sqrt(sigma1B[i]))/dmixnorm(dataB1,mean1B,sqrt(sigma1B),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2B[i],sqrt(sigma2B[i]))/dmixnorm(dataB2,mean2B,sqrt(sigma2B),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    for(i in 1:d23) { WW_F2[i,] <- mi_3[i]*dnorm(dataF2,mean3B[i],sqrt(sigma3B[i]))/dmixnorm(dataF2,mean3B,sqrt(sigma3B),mi_3) }
    mix_pi_3 <- as.matrix(rowSums(WW_F2)/n_samF2)
    sumwx_F2 <- WW_F2%*%dataF2

    ############ CM1-step for means ##############
    s0<- matrix(0,12,1)
    s0[1]<- sumx1+sumwx_B1[1]+sumwx_F2[1]
    s0[2]<- sumx2+sumwx_B1[4]+sumwx_B2[1]+sumwx_F2[5]
    s0[3]<- sumx3+sumwx_B2[4]+sumwx_F2[9]
    s0[4]<- sumwx_B1[2]+sumwx_F2[2]
    s0[5]<- sumwx_B1[3]+sumwx_F2[4]
    s0[6]<- sumwx_B2[2]+sumwx_F2[6]
    s0[7]<- sumwx_B2[3]+sumwx_F2[8]
    s0[8]<- sumwx_F2[3]
    s0[9]<- sumwx_F2[7]

    n0<- matrix(0,12,1)
    n0[1]<- n_samP1+mix_pi_1[1]*n_samB1+mix_pi_3[1]*n_samF2
    n0[2]<- n_samF1+mix_pi_1[4]*n_samB1+mix_pi_2[1]*n_samB2+mix_pi_3[5]*n_samF2
    n0[3]<- n_samP2+mix_pi_2[4]*n_samB2+mix_pi_3[9]*n_samF2
    n0[4]<- mix_pi_1[2]*n_samB1+mix_pi_3[2]*n_samF2
    n0[5]<- mix_pi_1[3]*n_samB1+mix_pi_3[4]*n_samF2
    n0[6]<- mix_pi_2[2]*n_samB2+mix_pi_3[6]*n_samF2
    n0[7]<- mix_pi_2[3]*n_samB2+mix_pi_3[8]*n_samF2
    n0[8]<- mix_pi_3[3]*n_samF2
    n0[9]<- mix_pi_3[7]*n_samF2

    aaa1<- 1000; n_iter<- 0;AA<- matrix(0,4,1)
    while (aaa1>0.0001){
      n_iter<- n_iter+1
      hh<- matrix(0,4,4)
      hh[1,1]<- sigma*(1.0/n0[1]+1.0/n0[3]+1.0/n0[8]+1.0/n0[9])
      hh[1,2]<- sigma/n0[1]
      hh[1,3]<- sigma/n0[3]
      hh[1,4]<- sigma*(1.0/n0[8]-1.0/n0[9])
      hh[2,2]<- sigma*(1.0/n0[1]+1.0/n0[2]+1.0/n0[4]+1.0/n0[5])
      hh[2,3]<- sigma/n0[2]
      hh[2,4]<- sigma*(-1.0/n0[4]+1.0/n0[5])
      hh[3,3]<- sigma*(1.0/n0[2]+1.0/n0[3]+1.0/n0[6]+1.0/n0[7])
      hh[3,4]<- sigma*(-1.0/n0[6]+1.0/n0[7])
      hh[4,4]<- sigma*(1.0/n0[4]+1.0/n0[5]+1.0/n0[6]+1.0/n0[7]+1.0/n0[8]+1.0/n0[9])
      for(i in 2:4){
        for(j in 1:(i-1)){
          hh[i,j]<- hh[j,i]
        }
      }
      b_line<- matrix(0,4,1)
      b_line[1]<- s0[1]/n0[1]+s0[3]/n0[3]-s0[8]/n0[8]-s0[9]/n0[9]
      b_line[2]<- s0[1]/n0[1]+s0[2]/n0[2]-s0[4]/n0[4]-s0[5]/n0[5]
      b_line[3]<- s0[2]/n0[2]+s0[3]/n0[3]-s0[6]/n0[6]-s0[7]/n0[7]
      b_line[4]<- s0[4]/n0[4]-s0[5]/n0[5]+s0[6]/n0[6]-s0[7]/n0[7]-s0[8]/n0[8]+s0[9]/n0[9]
      B <- solve(hh,b_line)

      meanP1<- (s0[1]-sigma*(B[1]+B[2]))/n0[1]
      meanF1<- (s0[2]-sigma*(B[2]+B[3]))/n0[2]
      meanP2<- (s0[3]-sigma*(B[1]+B[3]))/n0[3]

      mean1B[2]<- (s0[4]+(B[2]-B[4])*sigma)/n0[4]; mean1B[3]<- (s0[5]+(B[2]+B[4])*sigma)/n0[5]
      mean2B[2]<- (s0[6]+sigma*(B[3]-B[4]))/n0[6]; mean2B[3]<- (s0[7]+sigma*(B[3]+B[4]))/n0[7]
      mean3B[3]<- (s0[8]+sigma*(B[1]+B[4]))/n0[8]; mean3B[7]<- (s0[9]+sigma*(B[1]-B[4]))/n0[9]

      mean1B[1]<- meanP1;mean3B[1]<- meanP1
      mean1B[4]<- meanF1;mean2B[1]<- meanF1;mean3B[5]<- meanF1
      mean2B[4]<- meanP2;mean3B[9]<- meanP2
      mean3B[2]<- mean1B[2];mean3B[4]<- mean1B[3]
      mean3B[6]<- mean3B[2];mean3B[8]<- mean2B[3]

      aaa1<- max(abs(AA-B))
      AA<- B
      #if (n_iter>20) break
    }

    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2)
    ss2<- sum((dataF1-meanF1)^2)
    ss3<- sum((dataP2-meanP2)^2)

    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1B[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2B[i])^2 }
    for(i in 1:d23) {  swx_F2[i] <- WW_F2[i,]%*%(dataF2-mean3B[i])^2 }

    s0[11]<- ss1+ss2+ss3+sum(swx_B1)+sum(swx_B2)+sum(swx_F2)
    n0[11]<- n_samP1+n_samF1+n_samP2+n_samB1+n_samB2+n_samF2
    sigma<- s0[11]/n0[11]
    sigma1B<- matrix(sigma,d21,1);sigma2B<- matrix(sigma,d22,1);sigma3B<- matrix(sigma,d23,1)
    ########criteria for iterations to stop#######
    L1<- logL(n_samP1,1,1,meanP1,sigma,dataP1)+logL(n_samF1,1,1,meanF1,sigma,dataF1)+logL(n_samP2,1,1,meanP2,sigma,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1B,sigma1B,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2B,sigma2B,dataB2)+logL(n_samF2,d23,mix_pi_3,mean3B,sigma3B,dataF2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*6

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1,1,1, 1,0,-1,1,0,0,-1,1,-1, 1,0,-1,0,1,-1,0,-1,1, 0,1,0,0,1,1,0,0,0, 0,1,0,1,0,0,1,0,0),9,5)
  b_line1 <- matrix(c(meanP1,meanF1,meanP2,mean1B[2],mean1B[3],mean2B[2],mean2B[3],mean3B[3],mean3B[7]))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigmaB1 - sigma1B[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigmaB1

  jj_2 <- sigmaB2 - sigma2B[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigmaB2

  jj_3 <- sigmaF2 - sigma3B[1]
  if(jj_3 < 0) {jj_3 <- 0}
  ll_3 <- jj_3/sigmaF2

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
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

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
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),(1-pkolm(D_F1,n_samF1))))

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
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

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
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),(1-pkolm(D_B1,n_samB1))))

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
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),(1-pkolm(D_B2,n_samB2))))

  ###############  F2  ###################
  dataF2<-sort(dataF2);bmw_F2 <- matrix(0,n_samF2,1); bmwsl_F2 <- matrix(0,n_samF2,d23)
  for(i in 1:d23){
    gg_F2 <- (dataF2 - mean3B[i])/sqrt(sigma3B[i])
    bmw_F2[which(gg_F2>=0)] <- pnorm(gg_F2[gg_F2>=0])
    bmw_F2[which(gg_F2<0)] <- 1 - pnorm(abs(gg_F2[gg_F2<0]))
    bmwsl_F2[,i] <- bmw_F2*mix_pi_3[i]
  }
  P2_F2 <- rowSums(bmwsl_F2)
  nn<-dim(as.matrix(unique(P2_F2)))[1]
  if(nn<n_samF2){P2_F2<-P2_F2+runif(n_samF2)/1e4}

  dd_F2 <- as.matrix(c(sum(P2_F2),sum(P2_F2^2),sum((P2_F2-0.5)^2)))
  WW2_F2 <- 1/(12*n_samF2) + sum((P2_F2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  u_F2 <- as.matrix(c(12*n_samF2*((dd_F2[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((dd_F2[2]/n_samF2-1/3)^2),180*n_samF2*((dd_F2[3]/n_samF2-1/12)^2)))
  D_F2 <- as.numeric(ks.test(P2_F2,"punif")[[1]][1])
  tt_F2 <- as.matrix(c((1 - pchisq(u_F2[1],1)),(1 - pchisq(u_F2[2],1)),(1 - pchisq(u_F2[3],1)),K1(WW2_F2),(1-pkolm(D_F2,n_samF2))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)
  tt_F2[which(tt_F2>=10e-4)]<-round(tt_F2[which(tt_F2>=10e-4)],4);tt_F2[which(tt_F2<10e-4)]<-format(tt_F2[which(tt_F2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("2MG-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigma,4),round(t(mean1B),4),round(sigma1B[1],4),round(t(mix_pi_1),4),
                       round(t(mean2B),4),round(sigma2B[1],4),round(t(mix_pi_2),4),round(t(mean3B),4),round(sigma3B[1],4),round(t(mix_pi_3),4),
                       round(B1[1],4)," "," "," "," "," ",round(B1[2],4),round(B1[3],4),round(B1[4],4),round(B1[5],4)," "," "," "," "," "," ",
                       round(jj_1,4),round(ll_1*100,4)," "," ",round(jj_2,4),round(ll_2*100,4)," "," ",round(jj_3,4),round(ll_3*100,4)," "," ",round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),
                       tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),
                       tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),
                       tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5],round(u_F2[1],4),tt_F2[1],round(u_F2[2],4),
                       tt_F2[2],round(u_F2[3],4),tt_F2[3],round(WW2_F2,4),tt_F2[4],round(D_F2,4),tt_F2[5])

  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2,mi_3)
  return(OUTPUT)
}
############################################ 2MG-A Model #########################################
G6ModelFun[[7]] <- function(K1,logL,df11,df21,df31,df41,df51,df61){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1]; n_samF1<-dim(dataF1)[1]; n_samP2<-dim(dataP2)[1]; n_samB1<-dim(dataB1)[1]; n_samB2<-dim(dataB2)[1]; n_samF2<-dim(dataF2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2)
  sumx4<-sum(dataB1);sumx5<-sum(dataB2);sumx6<-sum(dataF2)
  sigmaP1<- as.numeric(var(dataP1)); ss1<- (n_samP1-1)*sigmaP1
  sigmaF1<- as.numeric(var(dataF1)); ss2<- (n_samF1-1)*sigmaF1
  sigmaP2<- as.numeric(var(dataP2)); ss3<- (n_samP2-1)*sigmaP2
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  mean4<- mean(dataB1); sigmaB1<- as.numeric(var(dataB1)); sigma40<- sigmaB1
  mean5<- mean(dataB2); sigmaB2<- as.numeric(var(dataB2)); sigma50<- sigmaB2
  mean6<- mean(dataF2); sigmaF2<- as.numeric(var(dataF2)); sigma60<- sigmaF2

  m_esp<-0.0001
  ####################### 2MG-A Model ############# (B3) ######################
  d21<-4; d22<-4; d23<-9
  sigma<- sigma0
  meanP1<- mean(dataP1);meanF1<- mean(dataF1);meanP2<- mean(dataP2)

  mi_1<- matrix(0.25,d21,1)
  a1<- sqrt(sigmaB1/(n_samB1-1))
  if (meanP1<meanP2) {a1= -a1}
  mean1B<- as.matrix(c(meanP1,mean4+0.5*a1,mean4-0.5*a1,meanF1))

  mi_2<- matrix(0.25,d22,1)
  a2<- sqrt(sigmaB2/(n_samB2-1))
  if (meanP1<meanP2) { a2= -a2 }
  mean2B<- as.matrix(c(meanF1,mean5+0.5*a2,mean5-0.5*a2,meanP2))

  mi_3<- as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  mean3B<- as.matrix(c(meanP1,mean1B[2],0.5*(meanP1+mean1B[2]),mean1B[3],meanF1,mean2B[2],0.5*(mean2B[3]+meanP2),mean2B[3],meanP2))

  sigma1B<- matrix(sigma,d21,1)
  sigma2B<- matrix(sigma,d22,1)
  sigma3B<- matrix(sigma,d23,1)

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1B,sigma1B,dataB1)+logL(n_samB2,d22,mi_2,mean2B,sigma2B,dataB2)+logL(n_samF2,d23,mi_3,mean3B,sigma3B,dataF2)

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)
  WW_F2 <- matrix(0,d23,n_samF2); swx_F2 <- matrix(0,d23,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1B[i],sqrt(sigma1B[i]))/dmixnorm(dataB1,mean1B,sqrt(sigma1B),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2B[i],sqrt(sigma2B[i]))/dmixnorm(dataB2,mean2B,sqrt(sigma2B),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    for(i in 1:d23) { WW_F2[i,] <- mi_3[i]*dnorm(dataF2,mean3B[i],sqrt(sigma3B[i]))/dmixnorm(dataF2,mean3B,sqrt(sigma3B),mi_3) }
    mix_pi_3 <- as.matrix(rowSums(WW_F2)/n_samF2)
    sumwx_F2 <- WW_F2%*%dataF2

    ############ CM1-step for means ##############
    s0<- matrix(0,12,1)
    s0[1]<- sumx1+sumwx_B1[1]+sumwx_F2[1]
    s0[2]<- sumx2+sumwx_B1[4]+sumwx_B2[1]+sumwx_F2[5]
    s0[3]<- sumx3+sumwx_B2[4]+sumwx_F2[9]
    s0[4]<- sumwx_B1[2]+sumwx_F2[2]
    s0[5]<- sumwx_B1[3]+sumwx_F2[4]
    s0[6]<- sumwx_B2[2]+sumwx_F2[6]
    s0[7]<- sumwx_B2[3]+sumwx_F2[8]
    s0[8]<- sumwx_F2[3]
    s0[9]<- sumwx_F2[7]

    n0<- matrix(0,12,1)
    n0[1]<- n_samP1+mix_pi_1[1]*n_samB1+mix_pi_3[1]*n_samF2
    n0[2]<- n_samF1+mix_pi_1[4]*n_samB1+mix_pi_2[1]*n_samB2+mix_pi_3[5]*n_samF2
    n0[3]<- n_samP2+mix_pi_2[4]*n_samB2+mix_pi_3[9]*n_samF2
    n0[4]<- mix_pi_1[2]*n_samB1+mix_pi_3[2]*n_samF2
    n0[5]<- mix_pi_1[3]*n_samB1+mix_pi_3[4]*n_samF2
    n0[6]<- mix_pi_2[2]*n_samB2+mix_pi_3[6]*n_samF2
    n0[7]<- mix_pi_2[3]*n_samB2+mix_pi_3[8]*n_samF2
    n0[8]<- mix_pi_3[3]*n_samF2
    n0[9]<- mix_pi_3[7]*n_samF2

    aa1<- 1000;n_iter<- 0; AA<- matrix(0,6,1)
    while (aa1>0.0001){
      n_iter<- n_iter+1
      hh<- matrix(0,6,6)
      hh[1,1]<- sigma*(1.0/n0[1]+4.0/n0[2]+1.0/n0[3])
      hh[1,2]<- sigma*(1.0/n0[1]+1.0/n0[3])
      hh[1,3]<- sigma*(1.0/n0[1]-2.0/n0[2])
      hh[1,4]<- sigma*(-2.0/n0[2]+1.0/n0[3])
      hh[1,5]<- hh[1,6]<- 0
      hh[2,2]<- sigma*(1.0/n0[1]+1.0/n0[3]+1.0/n0[8]+1.0/n0[9])
      hh[2,3]<- sigma/n0[1]
      hh[2,4]<- sigma/n0[3]
      hh[2,5]<- 0
      hh[2,6]<- sigma*(1.0/n0[8]-1.0/n0[9])
      hh[3,3]<- sigma*(1.0/n0[1]+1.0/n0[2]+1.0/n0[4]+1.0/n0[5])
      hh[3,4]<- sigma/n0[2]
      hh[3,5]<- sigma*(-1.0/n0[4]+1.0/n0[5])
      hh[3,6]<- 0
      hh[4,4]<- sigma*(1.0/n0[2]+1.0/n0[3]+1.0/n0[6]+1.0/n0[7])
      hh[4,5]<- sigma*(1.0/n0[6]-1.0/n0[7])
      hh[4,6]<- 2.0*sigma*(-1.0/n0[6]+1.0/n0[7])
      hh[5,5]<- sigma*(1.0/n0[4]+1.0/n0[5]+1.0/n0[6]+1.0/n0[7])
      hh[5,6]<- -2.0*sigma*(1.0/n0[6]+1.0/n0[7])
      hh[6,6]<- sigma*(4.0/n0[6]+4.0/n0[7]+1.0/n0[8]+1.0/n0[9])
      for(i in 2:6){
        for(j in 1:(i-1)){
          hh[i,j]<- hh[j,i]
        }
      }
      b_line<- matrix(0,6,1)
      b_line[1]<- s0[1]/n0[1]-2.0*s0[2]/n0[2]+s0[3]/n0[3]
      b_line[2]<- s0[1]/n0[1]+s0[3]/n0[3]-s0[8]/n0[8]-s0[9]/n0[9]
      b_line[3]<- s0[1]/n0[1]+s0[2]/n0[2]-s0[4]/n0[4]-s0[5]/n0[5]
      b_line[4]<- s0[2]/n0[2]+s0[3]/n0[3]-s0[6]/n0[6]-s0[7]/n0[7]
      b_line[5]<- s0[4]/n0[4]-s0[5]/n0[5]-s0[6]/n0[6]+s0[7]/n0[7]
      b_line[6]<- 2.0*s0[6]/n0[6]-2.0*s0[7]/n0[7]-s0[8]/n0[8]+s0[9]/n0[9]
      B <- solve(hh,b_line)

      meanP1<- (s0[1]-sigma*(B[1]+B[2]+B[3]))/n0[1]
      meanF1<- (s0[2]-sigma*(2.0*B[1]-B[3]-B[4]))/n0[2]
      meanP2<- (s0[3]-sigma*(B[1]+B[2]+B[4]))/n0[3]

      mean1B[2]<- (s0[4]+(B[3]-B[5])*sigma)/n0[4]; mean1B[3]<- (s0[5]+(B[3]+B[5])*sigma)/n0[5]
      mean2B[2]<- (s0[6]+sigma*(B[4]+B[5]-2.0*B[6]))/n0[6]; mean2B[3]<- (s0[7]+sigma*(B[4]-B[5]+2.0*B[6]))/n0[7]
      mean3B[3]<- (s0[8]+sigma*(B[2]+B[6]))/n0[8]; mean3B[7]<- (s0[9]+sigma*(B[2]-B[6]))/n0[9]

      mean1B[1]<- meanP1;mean1B[4]<- meanF1
      mean2B[1]<- meanF1;mean2B[4]<- meanP2

      mean3B[1]<- meanP1;mean3B[2]<- mean1B[2];mean3B[4]<- mean1B[3]
      mean3B[5]<- meanF1;mean3B[6]<- mean2B[2];mean3B[8]<- mean2B[3]
      mean3B[9]<- meanP2

      aa1<- max(abs(AA-B))
      AA<- B
      if (n_iter>20) break

    }

    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2)
    ss2<- sum((dataF1-meanF1)^2)
    ss3<- sum((dataP2-meanP2)^2)

    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1B[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2B[i])^2 }
    for(i in 1:d23) {  swx_F2[i] <- WW_F2[i,]%*%(dataF2-mean3B[i])^2 }

    s0[11]<- ss1+ss2+ss3+sum(swx_B1)+sum(swx_B2)+sum(swx_F2)
    n0[11]<- n_samP1+n_samF1+n_samP2+n_samB1+n_samB2+n_samF2
    sigma<- s0[11]/n0[11]
    sigma1B<- matrix(sigma,d21,1);sigma2B<- matrix(sigma,d22,1);sigma3B<- matrix(sigma,d23,1)
    ########criteria for iterations to stop#######
    L1<- logL(n_samP1,1,1,meanP1,sigma,dataP1)+logL(n_samF1,1,1,meanF1,sigma,dataF1)+logL(n_samP2,1,1,meanP2,sigma,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1B,sigma1B,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2B,sigma2B,dataB2)+logL(n_samF2,d23,mix_pi_3,mean3B,sigma3B,dataF2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*4

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1,1,1, 1,0,-1,1,0,0,-1,1,-1, 1,0,-1,0,1,-1,0,-1,1 ),9,3)
  b_line1 <- matrix(c(meanP1,meanF1,meanP2,mean1B[2],mean1B[3],mean2B[2],mean2B[3],mean3B[3],mean3B[7]))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigmaB1 - sigma1B[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigmaB1

  jj_2 <- sigmaB2 - sigma2B[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigmaB2

  jj_3 <- sigmaF2 - sigma3B[1]
  if(jj_3 < 0) {jj_3 <- 0}
  ll_3 <- jj_3/sigmaF2

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
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

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
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),(1-pkolm(D_F1,n_samF1))))

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
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

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
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),(1-pkolm(D_B1,n_samB1))))

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
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),(1-pkolm(D_B2,n_samB2))))

  ###############  F2  ###################
  dataF2<-sort(dataF2);bmw_F2 <- matrix(0,n_samF2,1); bmwsl_F2 <- matrix(0,n_samF2,d23)
  for(i in 1:d23){
    gg_F2 <- (dataF2 - mean3B[i])/sqrt(sigma3B[i])
    bmw_F2[which(gg_F2>=0)] <- pnorm(gg_F2[gg_F2>=0])
    bmw_F2[which(gg_F2<0)] <- 1 - pnorm(abs(gg_F2[gg_F2<0]))
    bmwsl_F2[,i] <- bmw_F2*mix_pi_3[i]
  }
  P2_F2 <- rowSums(bmwsl_F2)
  nn<-dim(as.matrix(unique(P2_F2)))[1]
  if(nn<n_samF2){P2_F2<-P2_F2+runif(n_samF2)/1e4}

  dd_F2 <- as.matrix(c(sum(P2_F2),sum(P2_F2^2),sum((P2_F2-0.5)^2)))
  WW2_F2 <- 1/(12*n_samF2) + sum((P2_F2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  u_F2 <- as.matrix(c(12*n_samF2*((dd_F2[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((dd_F2[2]/n_samF2-1/3)^2),180*n_samF2*((dd_F2[3]/n_samF2-1/12)^2)))
  D_F2 <- as.numeric(ks.test(P2_F2,"punif")[[1]][1])
  tt_F2 <- as.matrix(c((1 - pchisq(u_F2[1],1)),(1 - pchisq(u_F2[2],1)),(1 - pchisq(u_F2[3],1)),K1(WW2_F2),(1-pkolm(D_F2,n_samF2))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)
  tt_F2[which(tt_F2>=10e-4)]<-round(tt_F2[which(tt_F2>=10e-4)],4);tt_F2[which(tt_F2<10e-4)]<-format(tt_F2[which(tt_F2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("2MG-A",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigma,4),round(t(mean1B),4),round(sigma1B[1],4),round(t(mix_pi_1),4),
                       round(t(mean2B),4),round(sigma2B[1],4),round(t(mix_pi_2),4),round(t(mean3B),4),round(sigma3B[1],4),round(t(mix_pi_3),4),
                       round(B1[1],4)," "," "," "," "," ",round(B1[2],4),round(B1[3],4)," "," "," "," "," "," "," "," ",
                       round(jj_1,4),round(ll_1*100,4)," "," ",round(jj_2,4),round(ll_2*100,4)," "," ",round(jj_3,4),round(ll_3*100,4)," "," ",round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),
                       tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),
                       tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),
                       tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5],round(u_F2[1],4),tt_F2[1],round(u_F2[2],4),
                       tt_F2[2],round(u_F2[3],4),tt_F2[3],round(WW2_F2,4),tt_F2[4],round(D_F2,4),tt_F2[5])

  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2,mi_3)
  return(OUTPUT)
}
############################################ 2MG-EA Model #########################################
G6ModelFun[[8]] <- function(K1,logL,df11,df21,df31,df41,df51,df61){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1]; n_samF1<-dim(dataF1)[1]; n_samP2<-dim(dataP2)[1]; n_samB1<-dim(dataB1)[1]; n_samB2<-dim(dataB2)[1]; n_samF2<-dim(dataF2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2)
  sumx4<-sum(dataB1);sumx5<-sum(dataB2);sumx6<-sum(dataF2)
  sigmaP1<- as.numeric(var(dataP1)); ss1<- (n_samP1-1)*sigmaP1
  sigmaF1<- as.numeric(var(dataF1)); ss2<- (n_samF1-1)*sigmaF1
  sigmaP2<- as.numeric(var(dataP2)); ss3<- (n_samP2-1)*sigmaP2
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  mean4<- mean(dataB1); sigmaB1<- as.numeric(var(dataB1)); sigma40<- sigmaB1
  mean5<- mean(dataB2); sigmaB2<- as.numeric(var(dataB2)); sigma50<- sigmaB2
  mean6<- mean(dataF2); sigmaF2<- as.numeric(var(dataF2)); sigma60<- sigmaF2

  m_esp<-0.0001
  ####################### 2MG-EA Model ############# (B4) ######################
  d21<-3; d22<-3; d23<-5
  sigma<- sigma0
  meanP1<- mean(dataP1);meanF1<- mean(dataF1);meanP2<- mean(dataP2)

  mi_1<- as.matrix(c(0.25,0.5,0.25))
  sigma1B<- matrix(sigma,d21,1)
  mean1B<- as.matrix(c(meanP1,0.5*(meanP1+meanF1),meanF1))

  mi_2<- as.matrix(c(0.25,0.5,0.25))
  sigma2B<- matrix(sigma,d22,1)
  mean2B<- as.matrix(c(meanF1,0.5*(meanF1+meanP2),meanP2))

  mi_3<- as.matrix(c(0.0625,0.25,0.375,0.25,0.0625))
  sigma3B<- matrix(sigma,d23,1)
  mean3B<- as.matrix(c(meanP1,mean1B[2],meanF1,mean2B[2],meanP2))

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1B,sigma1B,dataB1)+logL(n_samB2,d22,mi_2,mean2B,sigma2B,dataB2)+logL(n_samF2,d23,mi_3,mean3B,sigma3B,dataF2)

  ##########iteration process###########
  iteration <- 0; stopa <- 1000; AA<- matrix(0,3,1)
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)
  WW_F2 <- matrix(0,d23,n_samF2); swx_F2 <- matrix(0,d23,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1B[i],sqrt(sigma1B[i]))/dmixnorm(dataB1,mean1B,sqrt(sigma1B),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2B[i],sqrt(sigma2B[i]))/dmixnorm(dataB2,mean2B,sqrt(sigma2B),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    for(i in 1:d23) { WW_F2[i,] <- mi_3[i]*dnorm(dataF2,mean3B[i],sqrt(sigma3B[i]))/dmixnorm(dataF2,mean3B,sqrt(sigma3B),mi_3) }
    mix_pi_3 <- as.matrix(rowSums(WW_F2)/n_samF2)
    sumwx_F2 <- WW_F2%*%dataF2

    ############ CM1-step for means ##############
    s0<- matrix(0,12,1)
    s0[1]<- sumx1+sumwx_B1[1]+sumwx_F2[1]
    s0[2]<- sumx2+sumwx_B1[3]+sumwx_B2[1]+sumwx_F2[3]
    s0[3]<- sumx3+sumwx_B2[3]+sumwx_F2[5]
    s0[4]<- sumwx_B1[2]+sumwx_F2[2]
    s0[5]<- sumwx_B2[2]+sumwx_F2[4]

    n0<- matrix(0,12,1)
    n0[1]<- n_samP1+mix_pi_1[1]*n_samB1+mix_pi_3[1]*n_samF2
    n0[2]<- n_samF1+mix_pi_1[3]*n_samB1+mix_pi_2[1]*n_samB2+mix_pi_3[3]*n_samF2
    n0[3]<- n_samP2+mix_pi_2[3]*n_samB2+mix_pi_3[5]*n_samF2
    n0[4]<- mix_pi_1[2]*n_samB1+mix_pi_3[2]*n_samF2
    n0[5]<- mix_pi_2[2]*n_samB2+mix_pi_3[4]*n_samF2
    aaa1<- 1000;n_iter<- 0
    while (aaa1>0.0001){
      n_iter<- n_iter+1
      hh<- matrix(0,3,3)
      hh[1,1]<- sigma*(1.0/n0[1]+4.0/n0[2]+1.0/n0[3])
      hh[1,2]<- sigma*(1.0/n0[1]-1.0/n0[3])
      hh[1,3]<- 4.0*sigma/n0[2]
      hh[2,2]<- sigma*(1.0/n0[1]+1.0/n0[3]+4.0/n0[4]+4.0/n0[5])
      hh[2,3]<- sigma*(-2.0/n0[4]+2.0/n0[5])
      hh[3,3]<- sigma*(1.0/n0[4]+4.0/n0[2]+1.0/n0[5])
      for(i in 2:3){
        for(j in 1:(i-1)){
          hh[i,j]<- hh[j,i]
        }
      }
      b_line<- matrix(0,3,1)
      b_line[1]<- s0[1]/n0[1]-2.0*s0[2]/n0[2]+s0[3]/n0[3]
      b_line[2]<- s0[1]/n0[1]-s0[3]/n0[3]-2.0*s0[4]/n0[4]+2.0*s0[5]/n0[5]
      b_line[3]<- s0[4]/n0[4]-2.0*s0[2]/n0[2]+s0[5]/n0[5]
      B <- solve(hh,b_line)

      meanP1<- (s0[1]-sigma*(B[1]+B[2]))/n0[1]
      meanF1<- (s0[2]+sigma*(2.0*B[1]+2.0*B[3]))/n0[2]
      meanP2<- (s0[3]-sigma*(B[1]-B[2]))/n0[3]

      mean1B[2]<- (s0[4]+sigma*(2.0*B[2]-B[3]))/n0[4]
      mean2B[2]<- (s0[5]-(2.0*B[2]+B[3])*sigma)/n0[5]

      aa1<- max(abs(AA-B))
      AA<- B
      if (n_iter>20) break
    }
    mean1B[1]<- meanP1;mean1B[3]<- meanF1
    mean2B[1]<- meanF1;mean2B[3]<- meanP2
    mean3B[1]<- meanP1;mean3B[2]<- mean1B[2];mean3B[3]<- meanF1
    mean3B[4]<- mean2B[2];mean3B[5]<- meanP2
    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2)
    ss2<- sum((dataF1-meanF1)^2)
    ss3<- sum((dataP2-meanP2)^2)

    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1B[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2B[i])^2 }
    for(i in 1:d23) {  swx_F2[i] <- WW_F2[i,]%*%(dataF2-mean3B[i])^2 }

    s0[6]<- ss1+ss2+ss3+sum(swx_B1)+sum(swx_B2)+sum(swx_F2)
    n0[6]<- n_samP1+n_samF1+n_samP2+n_samB1+n_samB2+n_samF2
    sigma<- s0[6]/n0[6]
    sigma1B<- matrix(sigma,d21,1);sigma2B<- matrix(sigma,d22,1);sigma3B<- matrix(sigma,d23,1)
    ########criteria for iterations to stop#######
    L1<- logL(n_samP1,1,1,meanP1,sigma,dataP1)+logL(n_samF1,1,1,meanF1,sigma,dataF1)+logL(n_samP2,1,1,meanP2,sigma,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1B,sigma1B,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2B,sigma2B,dataB2)+logL(n_samF2,d23,mix_pi_3,mean3B,sigma3B,dataF2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*3

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1, 2,0,-2,1,-1 ),5,2)
  b_line1 <- matrix(c(meanP1,meanF1,meanP2,mean1B[2],mean2B[2]))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigmaB1 - sigma1B[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigmaB1

  jj_2 <- sigmaB2 - sigma2B[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigmaB2

  jj_3 <- sigmaF2 - sigma3B[1]
  if(jj_3 < 0) {jj_3 <- 0}
  ll_3 <- jj_3/sigmaF2

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
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

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
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),(1-pkolm(D_F1,n_samF1))))

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
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

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
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),(1-pkolm(D_B1,n_samB1))))

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
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),(1-pkolm(D_B2,n_samB2))))

  ###############  F2  ###################
  dataF2<-sort(dataF2);bmw_F2 <- matrix(0,n_samF2,1); bmwsl_F2 <- matrix(0,n_samF2,d23)
  for(i in 1:d23){
    gg_F2 <- (dataF2 - mean3B[i])/sqrt(sigma3B[i])
    bmw_F2[which(gg_F2>=0)] <- pnorm(gg_F2[gg_F2>=0])
    bmw_F2[which(gg_F2<0)] <- 1 - pnorm(abs(gg_F2[gg_F2<0]))
    bmwsl_F2[,i] <- bmw_F2*mix_pi_3[i]
  }
  P2_F2 <- rowSums(bmwsl_F2)
  nn<-dim(as.matrix(unique(P2_F2)))[1]
  if(nn<n_samF2){P2_F2<-P2_F2+runif(n_samF2)/1e4}

  dd_F2 <- as.matrix(c(sum(P2_F2),sum(P2_F2^2),sum((P2_F2-0.5)^2)))
  WW2_F2 <- 1/(12*n_samF2) + sum((P2_F2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  u_F2 <- as.matrix(c(12*n_samF2*((dd_F2[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((dd_F2[2]/n_samF2-1/3)^2),180*n_samF2*((dd_F2[3]/n_samF2-1/12)^2)))
  D_F2 <- as.numeric(ks.test(P2_F2,"punif")[[1]][1])
  tt_F2 <- as.matrix(c((1 - pchisq(u_F2[1],1)),(1 - pchisq(u_F2[2],1)),(1 - pchisq(u_F2[3],1)),K1(WW2_F2),(1-pkolm(D_F2,n_samF2))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)
  tt_F2[which(tt_F2>=10e-4)]<-round(tt_F2[which(tt_F2>=10e-4)],4);tt_F2[which(tt_F2<10e-4)]<-format(tt_F2[which(tt_F2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("2MG-EA",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigma,4),round(t(mean1B),4)," ",round(sigma1B[1],4),round(t(mix_pi_1),4)," ",
                       round(t(mean2B),4)," ",round(sigma2B[1],4),round(t(mix_pi_2),4)," ",round(t(mean3B),4)," "," "," "," ",round(sigma3B[1],4),round(t(mix_pi_3),4)," "," "," "," ",
                       round(B1[1],4)," "," "," "," "," ",round(B1[2],4),round(B1[2],4)," "," "," "," "," "," "," "," ",
                       round(jj_1,4),round(ll_1*100,4)," "," ",round(jj_2,4),round(ll_2*100,4)," "," ",round(jj_3,4),round(ll_3*100,4)," "," ",round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),
                       tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),
                       tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),
                       tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5],round(u_F2[1],4),tt_F2[1],round(u_F2[2],4),
                       tt_F2[2],round(u_F2[3],4),tt_F2[3],round(WW2_F2,4),tt_F2[4],round(D_F2,4),tt_F2[5])

  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2,mi_3)
  return(OUTPUT)
}
############################################ 2MG-CD Model #########################################
G6ModelFun[[9]] <- function(K1,logL,df11,df21,df31,df41,df51,df61){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1]; n_samF1<-dim(dataF1)[1]; n_samP2<-dim(dataP2)[1]; n_samB1<-dim(dataB1)[1]; n_samB2<-dim(dataB2)[1]; n_samF2<-dim(dataF2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2)
  sumx4<-sum(dataB1);sumx5<-sum(dataB2);sumx6<-sum(dataF2)
  sigmaP1<- as.numeric(var(dataP1)); ss1<- (n_samP1-1)*sigmaP1
  sigmaF1<- as.numeric(var(dataF1)); ss2<- (n_samF1-1)*sigmaF1
  sigmaP2<- as.numeric(var(dataP2)); ss3<- (n_samP2-1)*sigmaP2
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  mean4<- mean(dataB1); sigmaB1<- as.numeric(var(dataB1)); sigma40<- sigmaB1
  mean5<- mean(dataB2); sigmaB2<- as.numeric(var(dataB2)); sigma50<- sigmaB2
  mean6<- mean(dataF2); sigmaF2<- as.numeric(var(dataF2)); sigma60<- sigmaF2

  m_esp<-0.0001
  ####################### 2MG-CD Model ################### (B5) ##################
  d21<-1; d22<-4; d23<-4
  sigma<- sigma0
  meanP1<- mean(dataP1);meanF1<- mean(dataF1);meanP2<- mean(dataP2)

  mi_1<- as.matrix(1);sigma1B<- as.matrix(sigma)
  mean1B<- as.matrix(meanP1)

  mi_2<- as.matrix(c(0.25,0.25,0.25,0.25));sigma2B<- matrix(sigma,d22,1)
  a2<- sqrt(sigmaB2/(n_samB2-1))
  if (meanP1<meanP2) { a2=-a2 }
  mean2B<- as.matrix(c(meanP1,mean5-a2,mean5-2.0*a2,meanP2))

  mi_3<- as.matrix(c(0.5625,0.1875,0.1875,0.0625));sigma3B<- matrix(sigma,d23,1)
  mean3B<- as.matrix(c(meanP1,mean2B[2],mean2B[3],meanP2))

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1B,sigma1B,dataB1)+logL(n_samB2,d22,mi_2,mean2B,sigma2B,dataB2)+logL(n_samF2,d23,mi_3,mean3B,sigma3B,dataF2)

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)
  WW_F2 <- matrix(0,d23,n_samF2); swx_F2 <- matrix(0,d23,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1B[i],sqrt(sigma1B[i]))/dmixnorm(dataB1,mean1B,sqrt(sigma1B),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2B[i],sqrt(sigma2B[i]))/dmixnorm(dataB2,mean2B,sqrt(sigma2B),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    for(i in 1:d23) { WW_F2[i,] <- mi_3[i]*dnorm(dataF2,mean3B[i],sqrt(sigma3B[i]))/dmixnorm(dataF2,mean3B,sqrt(sigma3B),mi_3) }
    mix_pi_3 <- as.matrix(rowSums(WW_F2)/n_samF2)
    sumwx_F2 <- WW_F2%*%dataF2

    ############ CM1-step for means ##############
    s0<- matrix(0,12,1)
    s0[1]<- sumx1+sumx2+sumx4+sumwx_B2[1]+sumwx_F2[1]
    s0[2]<- sumx3+sumwx_B2[4]+sumwx_F2[4]
    s0[3]<- sumwx_B2[2]+sumwx_F2[2]
    s0[4]<- sumwx_B2[3]+sumwx_F2[3]

    n0<- matrix(0,12,1)
    n0[1]<- n_samP1+n_samF1+n_samB1+mix_pi_2[1]*n_samB2+mix_pi_3[1]*n_samF2
    n0[2]<- n_samP2+mix_pi_2[4]*n_samB2+mix_pi_3[4]*n_samF2
    n0[3]<- mix_pi_2[2]*n_samB2+mix_pi_3[2]*n_samF2
    n0[4]<- mix_pi_2[3]*n_samB2+mix_pi_3[3]*n_samF2

    aa1<- s0[1]/n0[1]+s0[2]/n0[2]-s0[3]/n0[3]-s0[4]/n0[4]
    aa2<- sigma*(1.0/n0[1]+1.0/n0[2]+1.0/n0[3]+1.0/n0[4])
    aaa1<- aa1/aa2

    meanP1<- (s0[1]-sigma*aaa1)/n0[1]
    meanP2<- (s0[2]-sigma*aaa1)/n0[2]

    mean2B[2]<- (s0[3]+aaa1*sigma)/n0[3]
    mean2B[3]<- (s0[4]+aaa1*sigma)/n0[4]

    mean1B[1]<- meanP1;mean4<- meanP1;meanF1<- meanP1
    mean2B[1]<- meanP1;mean2B[4]<- meanP2
    mean3B[1]<- meanP1;mean3B[2]<- mean2B[2];mean3B[3]<- mean2B[3];mean3B[4]<- meanP2

    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2)
    ss2<- sum((dataF1-meanF1)^2)
    ss3<- sum((dataP2-meanP2)^2)

    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1B[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2B[i])^2 }
    for(i in 1:d23) {  swx_F2[i] <- WW_F2[i,]%*%(dataF2-mean3B[i])^2 }

    s0[10]<- ss1+ss2+ss3+sum(swx_B1)+sum(swx_B2)+sum(swx_F2)
    n0[10]<- n_samP1+n_samF1+n_samP2+n_samB1+n_samB2+n_samF2
    sigma<- s0[10]/n0[10]
    sigma1B<- matrix(sigma,d21,1);sigma2B<- matrix(sigma,d22,1);sigma3B<- matrix(sigma,d23,1)
    ########criteria for iterations to stop#######
    L1<- logL(n_samP1,1,1,meanP1,sigma,dataP1)+logL(n_samF1,1,1,meanF1,sigma,dataF1)+logL(n_samP2,1,1,meanP2,sigma,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1B,sigma1B,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2B,sigma2B,dataB2)+logL(n_samF2,d23,mix_pi_3,mean3B,sigma3B,dataF2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*4

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1, 1,-1,1,-1, 1,-1,-1,1 ),4,3)
  b_line1 <- matrix(c(meanP1,meanP2,mean2B[2],mean2B[3]))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigmaB1 - sigma1B[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigmaB1

  jj_2 <- sigmaB2 - sigma2B[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigmaB2

  jj_3 <- sigmaF2 - sigma3B[1]
  if(jj_3 < 0) {jj_3 <- 0}
  ll_3 <- jj_3/sigmaF2

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
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

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
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),(1-pkolm(D_F1,n_samF1))))

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
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

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
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),(1-pkolm(D_B1,n_samB1))))

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
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),(1-pkolm(D_B2,n_samB2))))

  ###############  F2  ###################
  dataF2<-sort(dataF2);bmw_F2 <- matrix(0,n_samF2,1); bmwsl_F2 <- matrix(0,n_samF2,d23)
  for(i in 1:d23){
    gg_F2 <- (dataF2 - mean3B[i])/sqrt(sigma3B[i])
    bmw_F2[which(gg_F2>=0)] <- pnorm(gg_F2[gg_F2>=0])
    bmw_F2[which(gg_F2<0)] <- 1 - pnorm(abs(gg_F2[gg_F2<0]))
    bmwsl_F2[,i] <- bmw_F2*mix_pi_3[i]
  }
  P2_F2 <- rowSums(bmwsl_F2)
  nn<-dim(as.matrix(unique(P2_F2)))[1]
  if(nn<n_samF2){P2_F2<-P2_F2+runif(n_samF2)/1e4}

  dd_F2 <- as.matrix(c(sum(P2_F2),sum(P2_F2^2),sum((P2_F2-0.5)^2)))
  WW2_F2 <- 1/(12*n_samF2) + sum((P2_F2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  u_F2 <- as.matrix(c(12*n_samF2*((dd_F2[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((dd_F2[2]/n_samF2-1/3)^2),180*n_samF2*((dd_F2[3]/n_samF2-1/12)^2)))
  D_F2 <- as.numeric(ks.test(P2_F2,"punif")[[1]][1])
  tt_F2 <- as.matrix(c((1 - pchisq(u_F2[1],1)),(1 - pchisq(u_F2[2],1)),(1 - pchisq(u_F2[3],1)),K1(WW2_F2),(1-pkolm(D_F2,n_samF2))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)
  tt_F2[which(tt_F2>=10e-4)]<-round(tt_F2[which(tt_F2>=10e-4)],4);tt_F2[which(tt_F2<10e-4)]<-format(tt_F2[which(tt_F2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("2MG-CD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigma,4),round(t(mean1B),4)," "," "," ",round(sigma1B[1],4),round(t(mix_pi_1),4)," "," "," ",
                       round(t(mean2B),4),round(sigma2B[1],4),round(t(mix_pi_2),4),round(t(mean3B),4)," "," "," "," "," ",round(sigma3B[1],4),round(t(mix_pi_3),4)," "," "," "," "," ",
                       round(B1[1],4)," "," "," "," "," ",round(B1[2],4),round(B1[2],4),round(B1[3],4),round(B1[3],4)," "," "," "," "," "," ",
                       round(jj_1,4),round(ll_1*100,4)," "," ",round(jj_2,4),round(ll_2*100,4)," "," ",round(jj_3,4),round(ll_3*100,4)," "," ",round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),
                       tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),
                       tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),
                       tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5],round(u_F2[1],4),tt_F2[1],round(u_F2[2],4),
                       tt_F2[2],round(u_F2[3],4),tt_F2[3],round(WW2_F2,4),tt_F2[4],round(D_F2,4),tt_F2[5])

  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2,mi_3)
  return(OUTPUT)
}
############################################ 2MG-EAD Model #########################################
G6ModelFun[[10]] <- function(K1,logL,df11,df21,df31,df41,df51,df61){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1]; n_samF1<-dim(dataF1)[1]; n_samP2<-dim(dataP2)[1]; n_samB1<-dim(dataB1)[1]; n_samB2<-dim(dataB2)[1]; n_samF2<-dim(dataF2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2)
  sumx4<-sum(dataB1);sumx5<-sum(dataB2);sumx6<-sum(dataF2)
  sigmaP1<- as.numeric(var(dataP1)); ss1<- (n_samP1-1)*sigmaP1
  sigmaF1<- as.numeric(var(dataF1)); ss2<- (n_samF1-1)*sigmaF1
  sigmaP2<- as.numeric(var(dataP2)); ss3<- (n_samP2-1)*sigmaP2
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  mean4<- mean(dataB1); sigmaB1<- as.numeric(var(dataB1)); sigma40<- sigmaB1
  mean5<- mean(dataB2); sigmaB2<- as.numeric(var(dataB2)); sigma50<- sigmaB2
  mean6<- mean(dataF2); sigmaF2<- as.numeric(var(dataF2)); sigma60<- sigmaF2

  m_esp<-0.0001
  ####################### 2MG-EAD Model ############# (B6) ####################
  d21<-1; d22<-3; d23<-3
  sigma<- sigma0
  meanP1<- mean(dataP1);meanF1<- mean(dataF1);meanP2<- mean(dataP2)

  mi_1<- as.matrix(1); sigma1B<- as.matrix(sigma)
  mean1B<- as.matrix(meanP1)

  mi_2<- as.matrix(c(0.25,0.5,0.25));sigma2B<- matrix(sigma,d22,1)
  mean2B<- as.matrix(c(meanP1,0.5*(meanP1+meanP2),meanP2))

  mi_3<- as.matrix(c(0.5625,0.375,0.0625));sigma3B<- matrix(sigma,d23,1)
  mean3B<- as.matrix(c(meanP1,mean2B[2],meanP2))

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1B,sigma1B,dataB1)+logL(n_samB2,d22,mi_2,mean2B,sigma2B,dataB2)+logL(n_samF2,d23,mi_3,mean3B,sigma3B,dataF2)

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)
  WW_F2 <- matrix(0,d23,n_samF2); swx_F2 <- matrix(0,d23,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1B[i],sqrt(sigma1B[i]))/dmixnorm(dataB1,mean1B,sqrt(sigma1B),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2B[i],sqrt(sigma2B[i]))/dmixnorm(dataB2,mean2B,sqrt(sigma2B),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    for(i in 1:d23) { WW_F2[i,] <- mi_3[i]*dnorm(dataF2,mean3B[i],sqrt(sigma3B[i]))/dmixnorm(dataF2,mean3B,sqrt(sigma3B),mi_3) }
    mix_pi_3 <- as.matrix(rowSums(WW_F2)/n_samF2)
    sumwx_F2 <- WW_F2%*%dataF2

    ############ CM1-step for means ##############
    s0<- matrix(0,12,1)
    s0[1]<- sumx1+sumx2+sumx4+sumwx_B2[1]+sumwx_F2[1]
    s0[2]<- sumx3+sumwx_B2[3]+sumwx_F2[3]
    s0[3]<- sumwx_B2[2]+sumwx_F2[2]

    n0<- matrix(0,12,1)
    n0[1]<- n_samP1+n_samF1+n_samB1+mix_pi_2[1]*n_samB2+mix_pi_3[1]*n_samF2
    n0[2]<- n_samP2+mix_pi_2[3]*n_samB2+mix_pi_3[3]*n_samF2
    n0[3]<- mix_pi_2[2]*n_samB2+mix_pi_3[2]*n_samF2

    aa1<- s0[1]/n0[1]+s0[2]/n0[2]-2.0*s0[3]/n0[3]
    aa2<- sigma*(1.0/n0[1]+1.0/n0[2]+4.0/n0[3])
    aaa1<- aa1/aa2

    meanP1<- (s0[1]-sigma*aaa1)/n0[1]
    meanP2<- (s0[2]-sigma*aaa1)/n0[2]

    mean2B[2]<- (s0[3]+2.0*aaa1*sigma)/n0[3]
    meanF1<- meanP1;mean1B[1]<- meanP1;mean4<- meanP1
    mean2B[1]<- meanP1;mean2B[3]<- meanP2
    mean3B[1]<- meanP1;mean3B[2]<- mean2B[2];mean3B[3]<- meanP2

    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2)
    ss2<- sum((dataF1-meanF1)^2)
    ss3<- sum((dataP2-meanP2)^2)

    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1B[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2B[i])^2 }
    for(i in 1:d23) {  swx_F2[i] <- WW_F2[i,]%*%(dataF2-mean3B[i])^2 }

    s0[7]<- ss1+ss2+ss3+sum(swx_B1)+sum(swx_B2)+sum(swx_F2)
    n0[7]<- n_samP1+n_samF1+n_samP2+n_samB1+n_samB2+n_samF2
    sigma<- s0[7]/n0[7]
    sigma1B<- matrix(sigma,d21,1);sigma2B<- matrix(sigma,d22,1);sigma3B<- matrix(sigma,d23,1)
    ########criteria for iterations to stop#######
    L1<- logL(n_samP1,1,1,meanP1,sigma,dataP1)+logL(n_samF1,1,1,meanF1,sigma,dataF1)+logL(n_samP2,1,1,meanP2,sigma,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1B,sigma1B,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2B,sigma2B,dataB2)+logL(n_samF2,d23,mix_pi_3,mean3B,sigma3B,dataF2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*3

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1, 2,-2,0),3,2)
  b_line1 <- matrix(c(meanP1,meanP2,mean2B[2]))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigmaB1 - sigma1B[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigmaB1

  jj_2 <- sigmaB2 - sigma2B[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigmaB2

  jj_3 <- sigmaF2 - sigma3B[1]
  if(jj_3 < 0) {jj_3 <- 0}
  ll_3 <- jj_3/sigmaF2

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
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

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
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),(1-pkolm(D_F1,n_samF1))))

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
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

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
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),(1-pkolm(D_B1,n_samB1))))

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
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),(1-pkolm(D_B2,n_samB2))))

  ###############  F2  ###################
  dataF2<-sort(dataF2);bmw_F2 <- matrix(0,n_samF2,1); bmwsl_F2 <- matrix(0,n_samF2,d23)
  for(i in 1:d23){
    gg_F2 <- (dataF2 - mean3B[i])/sqrt(sigma3B[i])
    bmw_F2[which(gg_F2>=0)] <- pnorm(gg_F2[gg_F2>=0])
    bmw_F2[which(gg_F2<0)] <- 1 - pnorm(abs(gg_F2[gg_F2<0]))
    bmwsl_F2[,i] <- bmw_F2*mix_pi_3[i]
  }
  P2_F2 <- rowSums(bmwsl_F2)
  nn<-dim(as.matrix(unique(P2_F2)))[1]
  if(nn<n_samF2){P2_F2<-P2_F2+runif(n_samF2)/1e4}

  dd_F2 <- as.matrix(c(sum(P2_F2),sum(P2_F2^2),sum((P2_F2-0.5)^2)))
  WW2_F2 <- 1/(12*n_samF2) + sum((P2_F2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  u_F2 <- as.matrix(c(12*n_samF2*((dd_F2[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((dd_F2[2]/n_samF2-1/3)^2),180*n_samF2*((dd_F2[3]/n_samF2-1/12)^2)))
  D_F2 <- as.numeric(ks.test(P2_F2,"punif")[[1]][1])
  tt_F2 <- as.matrix(c((1 - pchisq(u_F2[1],1)),(1 - pchisq(u_F2[2],1)),(1 - pchisq(u_F2[3],1)),K1(WW2_F2),(1-pkolm(D_F2,n_samF2))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)
  tt_F2[which(tt_F2>=10e-4)]<-round(tt_F2[which(tt_F2>=10e-4)],4);tt_F2[which(tt_F2<10e-4)]<-format(tt_F2[which(tt_F2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("2MG-EAD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigma,4),round(t(mean1B),4)," "," "," ",round(sigma1B[1],4),round(t(mix_pi_1),4)," "," "," ",
                       round(t(mean2B),4)," ",round(sigma2B[1],4),round(t(mix_pi_2),4)," ",round(t(mean3B),4)," "," "," "," "," "," ",round(sigma3B[1],4),round(t(mix_pi_3),4)," "," "," "," "," "," ",
                       round(B1[1],4)," "," "," "," "," ",round(B1[2],4),round(B1[2],4),round(B1[2],4),round(B1[2],4)," "," "," "," "," "," ",
                       round(jj_1,4),round(ll_1*100,4)," "," ",round(jj_2,4),round(ll_2*100,4)," "," ",round(jj_3,4),round(ll_3*100,4)," "," ",round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),
                       tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),
                       tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),
                       tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5],round(u_F2[1],4),tt_F2[1],round(u_F2[2],4),
                       tt_F2[2],round(u_F2[3],4),tt_F2[3],round(WW2_F2,4),tt_F2[4],round(D_F2,4),tt_F2[5])

  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2,mi_3)
  return(OUTPUT)
}

############################################ PG-ADI Model #########################################
G6ModelFun[[11]] <- function(K1,logL,df11,df21,df31,df41,df51,df61){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1]; n_samF1<-dim(dataF1)[1]; n_samP2<-dim(dataP2)[1]; n_samB1<-dim(dataB1)[1]; n_samB2<-dim(dataB2)[1]; n_samF2<-dim(dataF2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2)
  sumx4<-sum(dataB1);sumx5<-sum(dataB2);sumx6<-sum(dataF2)
  sigmaP1<- as.numeric(var(dataP1)); ss1<- (n_samP1-1)*sigmaP1
  sigmaF1<- as.numeric(var(dataF1)); ss2<- (n_samF1-1)*sigmaF1
  sigmaP2<- as.numeric(var(dataP2)); ss3<- (n_samP2-1)*sigmaP2
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  mean4<- mean(dataB1); sigmaB1<- as.numeric(var(dataB1)); sigma40<- sigmaB1
  mean5<- mean(dataB2); sigmaB2<- as.numeric(var(dataB2)); sigma50<- sigmaB2
  mean6<- mean(dataF2); sigmaF2<- as.numeric(var(dataF2)); sigma60<- sigmaF2

  m_esp<-0.0001
  #######################  PG-ADI Model ################  (C0) ####################
  d21<-1; d22<-1; d23<-1
  sigma<- sigma0
  meanP1<- mean(dataP1);meanF1<- mean(dataF1);meanP2<- mean(dataP2)

  mix_pi_1<- 1; sigma1C<- sigmaB1; mean1C<- mean4
  mix_pi_2<- 1; sigma2C<- sigmaB2; mean2C<- mean5
  mix_pi_3<- 1; sigma3C<- sigmaF2; mean3C<- mean6

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1C,sigma40,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2C,sigma50,dataB2)+logL(n_samF2,d23,mix_pi_3,mean3C,sigma60,dataF2)

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)
  WW_F2 <- matrix(0,d23,n_samF2); swx_F2 <- matrix(0,d23,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2)
    ss2<- sum((dataF1-meanF1)^2)
    ss3<- sum((dataP2-meanP2)^2)
    ss4<- sum((dataB1-mean4)^2)
    ss5<- sum((dataB2-mean5)^2)
    ss6<- sum((dataF2-mean6)^2)
    abc1<- ss1+ss2+ss3
    abc2<- n_samP1+n_samF1+n_samP2
    aaa0<- sigma;aa3<- 1000
    while (aa3>0.0001){
      aa1<- n_samB1*sigma/ss4
      aa2<- n_samB2*sigma/ss5
      aa3<- n_samF2*sigma/ss6
      sigma<- (abc1+aa1*aa1*ss4+aa2*aa2*ss5+aa3*aa3*ss6)/(abc2+aa1*n_samB1+aa2*n_samB2+aa3*n_samF2)
      aa3<- abs(sigma-aaa0)
      aaa0<- sigma
    }
    sigma40<- ss4/n_samB1; sigma50<- ss5/n_samB2; sigma60<- ss6/n_samF2
    ########criteria for iterations to stop#######
    L1<- logL(n_samP1,1,1,meanP1,sigma,dataP1)+logL(n_samF1,1,1,meanF1,sigma,dataF1)+logL(n_samP2,1,1,meanP2,sigma,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1C,sigma40,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2C,sigma50,dataB2)+logL(n_samF2,d23,mix_pi_3,mean3C,sigma60,dataF2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*10

  #########first order genetic parameter process##########
  B1<- as.matrix(c(meanP1,meanF1,meanP2,mean1C,mean2C,mean3C))
  mm1<- sigma40-sigma
  if (mm1<0 | mm1>=sigma40) { mm1<- 0 }
  nn1<- mm1/sigma40

  mm2<- sigma50-sigma
  if (mm2<0 | mm2>=sigma50) { mm2<- 0 }
  nn2<- mm2/sigma50

  mm3<- sigma60-sigma
  if (mm3<0 | mm3>=sigma60) { mm3<- 0 }
  nn3<- mm3/sigma60

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
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

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
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),(1-pkolm(D_F1,n_samF1))))

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
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  B1  ###################
  dataB1<-sort(dataB1);bmw_B1 <- matrix(0,n_samB1,1); bmwsl_B1 <- matrix(0,n_samB1,d21)
  for(i in 1:d21){
    gg_B1 <- (dataB1 - mean1C[i])/sqrt(as.vector(sigma1C[i]))
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
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),(1-pkolm(D_B1,n_samB1))))

  ###############  B2  ###################
  dataB2<-sort(dataB2);bmw_B2 <- matrix(0,n_samB2,1); bmwsl_B2 <- matrix(0,n_samB2,d22)
  for(i in 1:d22){
    gg_B2 <- (dataB2 - mean2C[i])/sqrt(as.vector(sigma2C[i]))
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
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),(1-pkolm(D_B2,n_samB2))))

  ###############  F2  ###################
  dataF2<-sort(dataF2);bmw_F2 <- matrix(0,n_samF2,1); bmwsl_F2 <- matrix(0,n_samF2,d23)
  for(i in 1:d23){
    gg_F2 <- (dataF2 - mean3C[i])/sqrt(as.vector(sigma3C[i]))
    bmw_F2[which(gg_F2>=0)] <- pnorm(gg_F2[gg_F2>=0])
    bmw_F2[which(gg_F2<0)] <- 1 - pnorm(abs(gg_F2[gg_F2<0]))
    bmwsl_F2[,i] <- bmw_F2*mix_pi_3[i]
  }
  P2_F2 <- rowSums(bmwsl_F2)
  nn<-dim(as.matrix(unique(P2_F2)))[1]
  if(nn<n_samF2){P2_F2<-P2_F2+runif(n_samF2)/1e4}

  dd_F2 <- as.matrix(c(sum(P2_F2),sum(P2_F2^2),sum((P2_F2-0.5)^2)))
  WW2_F2 <- 1/(12*n_samF2) + sum((P2_F2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  u_F2 <- as.matrix(c(12*n_samF2*((dd_F2[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((dd_F2[2]/n_samF2-1/3)^2),180*n_samF2*((dd_F2[3]/n_samF2-1/12)^2)))
  D_F2 <- as.numeric(ks.test(P2_F2,"punif")[[1]][1])
  tt_F2 <- as.matrix(c((1 - pchisq(u_F2[1],1)),(1 - pchisq(u_F2[2],1)),(1 - pchisq(u_F2[3],1)),K1(WW2_F2),(1-pkolm(D_F2,n_samF2))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)
  tt_F2[which(tt_F2>=10e-4)]<-round(tt_F2[which(tt_F2>=10e-4)],4);tt_F2[which(tt_F2<10e-4)]<-format(tt_F2[which(tt_F2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("PG-ADI",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigma,4),round(t(mean1C),4)," "," "," ",round(sigma1C[1],4),round(t(mix_pi_1),4)," "," "," ",
                       round(t(mean2C),4)," "," "," ",round(sigma2C[1],4),round(t(mix_pi_2),4)," "," "," ",round(t(mean3C),4)," "," "," "," "," "," "," "," ",round(sigma3C[1],4),round(t(mix_pi_3),4),
                       " "," "," "," "," "," "," "," ",round(B1[1],4),round(B1[2],4),round(B1[3],4),round(B1[4],4),round(B1[5],4),round(B1[6],4)," "," "," "," "," "," "," "," "," "," ",
                       " "," ",round(mm1,4),round(nn1*100,4)," "," ",round(mm2,4),round(nn2*100,4)," "," ",round(mm3,4),round(nn3*100,4),round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),
                       tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),
                       tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),
                       tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5],round(u_F2[1],4),tt_F2[1],round(u_F2[2],4),
                       tt_F2[2],round(u_F2[3],4),tt_F2[3],round(WW2_F2,4),tt_F2[4],round(D_F2,4),tt_F2[5])

  output<-as.matrix(output)
  OUTPUT<-list(output)
  return(OUTPUT)
}
############################################ PG-AD Model #########################################
G6ModelFun[[12]] <- function(K1,logL,df11,df21,df31,df41,df51,df61){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1]; n_samF1<-dim(dataF1)[1]; n_samP2<-dim(dataP2)[1]; n_samB1<-dim(dataB1)[1]; n_samB2<-dim(dataB2)[1]; n_samF2<-dim(dataF2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2)
  sumx4<-sum(dataB1);sumx5<-sum(dataB2);sumx6<-sum(dataF2)
  sigmaP1<- as.numeric(var(dataP1)); ss1<- (n_samP1-1)*sigmaP1
  sigmaF1<- as.numeric(var(dataF1)); ss2<- (n_samF1-1)*sigmaF1
  sigmaP2<- as.numeric(var(dataP2)); ss3<- (n_samP2-1)*sigmaP2
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  mean4<- mean(dataB1); sigmaB1<- as.numeric(var(dataB1)); sigma40<- sigmaB1
  mean5<- mean(dataB2); sigmaB2<- as.numeric(var(dataB2)); sigma50<- sigmaB2
  mean6<- mean(dataF2); sigmaF2<- as.numeric(var(dataF2)); sigma60<- sigmaF2

  m_esp<-0.0001
  #######################  PG-AD Model ################### (C1) ##############
  d21<-1; d22<-1; d23<-1
  sigma<- sigma0
  meanP1<- mean(dataP1);meanF1<- mean(dataF1);meanP2<- mean(dataP2)

  mix_pi_1<- 1; mean1C<- mean4
  mix_pi_2<- 1; mean2C<- mean5
  mix_pi_3<- 1; mean3C<- mean6

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1C,sigma40,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2C,sigma50,dataB2)+logL(n_samF2,d23,mix_pi_3,mean3C,sigma60,dataF2)

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)
  WW_F2 <- matrix(0,d23,n_samF2); swx_F2 <- matrix(0,d23,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    ################# iteratively CM1-step for means ####################
    sigma40<- sum((dataB1-mean4)^2)
    sigma40<- sigma40/n_samB1
    sigma50<- sum((dataB2-mean5)^2)
    sigma50<- sigma50/n_samB2
    sigma60<- sum((dataF2-mean6)^2)
    sigma60<- sigma60/n_samF2

    hh<- matrix(0,3,3)
    hh[1,1]<- sigma*(1.0/n_samP1+1.0/n_samF1)+4.0*sigma40/n_samB1
    hh[1,2]<- sigma/n_samF1
    hh[1,3]<- -2.0*sigma40/n_samB1
    hh[2,2]<- sigma*(1.0/n_samF1+1.0/n_samP2)+4.0*sigma50/n_samB2
    hh[2,3]<- -2.0*sigma50/n_samB2
    hh[3,3]<- sigma40/n_samB1+sigma50/n_samB2+4.0*sigma60/n_samF2
    for(i in 2:3){
      for(j in 1:(i-1)){
        hh[i,j]<- hh[j,i]
      }
    }
    b_line<- matrix(0,3,1)
    b_line[1]<- sumx1/n_samP1+sumx2/n_samF1-2.0*sumx4/n_samB1
    b_line[2]<- sumx2/n_samF1+sumx3/n_samP2-2.0*sumx5/n_samB2
    b_line[3]<- sumx4/n_samB1+sumx5/n_samB2-2.0*sumx6/n_samF2
    B <- solve(hh,b_line)

    meanP1<- (sumx1-sigma*B[1])/n_samP1
    meanF1<- (sumx2-sigma*(B[1]+B[2]))/n_samF1
    meanP2<- (sumx3-sigma*B[2])/n_samP2
    mean1C[1]<- (sumx4+sigma40*(2.0*B[1]-B[3]))/n_samB1
    mean2C[1]<- (sumx5+(2.0*B[2]-B[3])*sigma50)/n_samB2
    mean3C[1]<- (sumx6+2.0*B[3]*sigma60)/n_samF2

    ################# iteratively CM2-step for variance ####################
    aaa0<- sigma
    ss1<- sum((dataP1-meanP1)^2)
    ss2<- sum((dataF1-meanF1)^2)
    ss3<- sum((dataP2-meanP2)^2)
    abc1<- ss1+ss2+ss3
    abc2<- n_samP1+n_samF1+n_samP2
    ss4<- sum((dataB1-mean4)^2)
    ss5<- sum((dataB2-mean5)^2)
    ss6<- sum((dataF2-mean6)^2)
    sigma40<- ss4/n_samB1;sigma_4<- sigma40-sigma
    sigma50<- ss5/n_samB2;sigma_5<- sigma50-sigma
    sigma60<- ss6/n_samF2;sigma_6<- sigma60-sigma
    ################ iteratively CM3-step for variance (sigma) ################
    sigma_4[(sigma_4<0)]<- 0.00001; sigma_5[(sigma_5<0)]<- 0.00001; sigma_6[(sigma_6<0)]<- 0.00001
    aa3<- 1000; n_iter<- 0
    while (aa3>0.0001){
      n_iter<- n_iter+1
      aa1<- sigma/(sigma_4+sigma)
      if (aa1>=1) { aa1<- 1 }
      aa2<- sigma/(sigma_5+sigma)
      if (aa2>=1) { aa2<- 1 }
      aa3<- sigma/(sigma_6+sigma)
      if (aa3>=1) { aa3<- 1 }
      aa4<- abc1+aa1*aa1*ss4+aa2*aa2*ss5+aa3*aa3*ss6
      aa5<- abc2+aa1*n_samB1+aa2*n_samB2+aa3*n_samF2
      sigma<- aa4/aa5
      aa3<- abs(sigma-aaa0)
      aaa0<- sigma
      if (n_iter>20) break
    }
    sigma40<- sigma_4+sigma;sigma50<- sigma_5+sigma;sigma60<- sigma_6+sigma
    ################### the stop criterion for iteration #########################
    L1<- logL(n_samP1,1,1,meanP1,sigma,dataP1)+logL(n_samF1,1,1,meanF1,sigma,dataF1)+logL(n_samP2,1,1,meanP2,sigma,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1C,sigma40,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2C,sigma50,dataB2)+logL(n_samF2,d23,mix_pi_3,mean3C,sigma60,dataF2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*7
  sigma1C<- sigma40; sigma2C<- sigma50; sigma3C<- sigma60;
  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1, 1,0,-1,0.5,-0.5,0, 0,1,0,0.5,0.5,0.5),6,3)
  b_line1 <- matrix(c(meanP1,meanF1,meanP2,mean1C[1],mean2C[1],mean3C[1]))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)
  mm1<- sigma40-sigma
  if (mm1<0 | mm1>=sigma40) { mm1<- 0 }
  nn1<- mm1/sigma40

  mm2<- sigma50-sigma
  if (mm2<0 | mm2>=sigma50) { mm2<- 0 }
  nn2<- mm2/sigma50

  mm3<- sigma60-sigma
  if (mm3<0 | mm3>=sigma60) { mm3<- 0 }
  nn3<- mm3/sigma60

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
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

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
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),(1-pkolm(D_F1,n_samF1))))

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
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

  ###############  B1  ###################
  dataB1<-sort(dataB1);bmw_B1 <- matrix(0,n_samB1,1); bmwsl_B1 <- matrix(0,n_samB1,d21)
  for(i in 1:d21){
    gg_B1 <- (dataB1 - mean1C[i])/sqrt(as.vector(sigma1C[i]))
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
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),(1-pkolm(D_B1,n_samB1))))

  ###############  B2  ###################
  dataB2<-sort(dataB2);bmw_B2 <- matrix(0,n_samB2,1); bmwsl_B2 <- matrix(0,n_samB2,d22)
  for(i in 1:d22){
    gg_B2 <- (dataB2 - mean2C[i])/sqrt(as.vector(sigma2C[i]))
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
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),(1-pkolm(D_B2,n_samB2))))

  ###############  F2  ###################
  dataF2<-sort(dataF2);bmw_F2 <- matrix(0,n_samF2,1); bmwsl_F2 <- matrix(0,n_samF2,d23)
  for(i in 1:d23){
    gg_F2 <- (dataF2 - mean3C[i])/sqrt(as.vector(sigma3C[i]))
    bmw_F2[which(gg_F2>=0)] <- pnorm(gg_F2[gg_F2>=0])
    bmw_F2[which(gg_F2<0)] <- 1 - pnorm(abs(gg_F2[gg_F2<0]))
    bmwsl_F2[,i] <- bmw_F2*mix_pi_3[i]
  }
  P2_F2 <- rowSums(bmwsl_F2)
  nn<-dim(as.matrix(unique(P2_F2)))[1]
  if(nn<n_samF2){P2_F2<-P2_F2+runif(n_samF2)/1e4}

  dd_F2 <- as.matrix(c(sum(P2_F2),sum(P2_F2^2),sum((P2_F2-0.5)^2)))
  WW2_F2 <- 1/(12*n_samF2) + sum((P2_F2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  u_F2 <- as.matrix(c(12*n_samF2*((dd_F2[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((dd_F2[2]/n_samF2-1/3)^2),180*n_samF2*((dd_F2[3]/n_samF2-1/12)^2)))
  D_F2 <- as.numeric(ks.test(P2_F2,"punif")[[1]][1])
  tt_F2 <- as.matrix(c((1 - pchisq(u_F2[1],1)),(1 - pchisq(u_F2[2],1)),(1 - pchisq(u_F2[3],1)),K1(WW2_F2),(1-pkolm(D_F2,n_samF2))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)
  tt_F2[which(tt_F2>=10e-4)]<-round(tt_F2[which(tt_F2>=10e-4)],4);tt_F2[which(tt_F2<10e-4)]<-format(tt_F2[which(tt_F2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("PG-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigma,4),round(t(mean1C),4)," "," "," ",round(sigma1C[1],4),round(t(mix_pi_1),4)," "," "," ",
                       round(t(mean2C),4)," "," "," ",round(sigma2C[1],4),round(t(mix_pi_2),4)," "," "," ",round(t(mean3C),4)," "," "," "," "," "," "," "," ",round(sigma3C[1],4),round(t(mix_pi_3),4),
                       " "," "," "," "," "," "," "," ",round(B1[1],4)," "," "," "," "," "," "," "," "," "," "," "," "," ",round(B1[2],4),round(B1[3],4),
                       " "," ",round(mm1,4),round(nn1*100,4)," "," ",round(mm2,4),round(nn2*100,4)," "," ",round(mm3,4),round(nn3*100,4),round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),
                       tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),
                       tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),
                       tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5],round(u_F2[1],4),tt_F2[1],round(u_F2[2],4),
                       tt_F2[2],round(u_F2[3],4),tt_F2[3],round(WW2_F2,4),tt_F2[4],round(D_F2,4),tt_F2[5])

  output<-as.matrix(output)
  OUTPUT<-list(output)
  return(OUTPUT)
}
############################################ MX1-AD-ADI Model #########################################
G6ModelFun[[13]] <- function(K1,logL,df11,df21,df31,df41,df51,df61){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1]; n_samF1<-dim(dataF1)[1]; n_samP2<-dim(dataP2)[1]; n_samB1<-dim(dataB1)[1]; n_samB2<-dim(dataB2)[1]; n_samF2<-dim(dataF2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2)
  sumx4<-sum(dataB1);sumx5<-sum(dataB2);sumx6<-sum(dataF2)
  sigmaP1<- as.numeric(var(dataP1)); ss1<- (n_samP1-1)*sigmaP1
  sigmaF1<- as.numeric(var(dataF1)); ss2<- (n_samF1-1)*sigmaF1
  sigmaP2<- as.numeric(var(dataP2)); ss3<- (n_samP2-1)*sigmaP2
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  mean4<- mean(dataB1); sigmaB1<- as.numeric(var(dataB1)); sigma40<- sigmaB1
  mean5<- mean(dataB2); sigmaB2<- as.numeric(var(dataB2)); sigma50<- sigmaB2
  mean6<- mean(dataF2); sigmaF2<- as.numeric(var(dataF2)); sigma60<- sigmaF2

  m_esp<-0.0001
  ####################### MX1-AD-ADI Model #################   (D0) ###################
  d21<-2; d22<-2; d23<-3
  sigma<- sigma0
  meanP1<- mean(dataP1);meanF1<- mean(dataF1);meanP2<- mean(dataP2)
  abb<- sigmaF2/(5.0*sigma0)

  mi_1<- matrix(0.5,d21,1);sigma1D<- matrix(sigmaB1/abb,d21,1)
  a1<- sqrt(sigmaB1/n_samB1)
  if (meanP1<meanP2) {a1= -a1}
  mean1D<- as.matrix(c(mean4+2*a1,mean4))

  mi_2<- matrix(0.5,d22,1);sigma2D<- matrix(sigmaB2/abb,d22,1)
  a2<- sqrt(sigmaB2/n_samB2)
  if (meanP1<meanP2) { a2= -a2 }
  mean2D<- as.matrix(c(mean5,mean5-2*a2))

  mi_3<- as.matrix(c(0.25,0.5,0.25));sigma3D<- matrix(sigmaF2/abb,d23,1)
  a3<- sqrt(sigmaF2/n_samF2)
  if (meanP1<meanP2) { a3= -a3 }
  mean3D<- as.matrix(c(mean6+2*a3,mean6,mean6-2*a3))

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1D,sigma1D,dataB1)+logL(n_samB2,d22,mi_2,mean2D,sigma2D,dataB2)+logL(n_samF2,d23,mi_3,mean3D,sigma3D,dataF2)

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)
  WW_F2 <- matrix(0,d23,n_samF2); swx_F2 <- matrix(0,d23,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1D[i],sqrt(sigma1D[i]))/dmixnorm(dataB1,mean1D,sqrt(sigma1D),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2D[i],sqrt(sigma2D[i]))/dmixnorm(dataB2,mean2D,sqrt(sigma2D),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    for(i in 1:d23) { WW_F2[i,] <- mi_3[i]*dnorm(dataF2,mean3D[i],sqrt(sigma3D[i]))/dmixnorm(dataF2,mean3D,sqrt(sigma3D),mi_3) }
    mix_pi_3 <- as.matrix(rowSums(WW_F2)/n_samF2)
    sumwx_F2 <- WW_F2%*%dataF2

    ############ CM1-step for means ##############
    n0<- as.matrix(c(mix_pi_1*n_samB1,mix_pi_2*n_samB2,mix_pi_3*n_samF2))
    n0[abs(n0)<0.0001] <- 0.0001

    aaa1<- 1000;n_iter<- 0; AA<- matrix(0,2,1)
    s0<- matrix(0,6,1)
    s0[1]<- sumwx_B1[1]/n0[1]-sumwx_B1[2]/n0[2]-sumwx_F2[1]/n0[5]+sumwx_F2[2]/n0[6]
    s0[2]<- sumwx_B2[1]/n0[3]-sumwx_B2[2]/n0[4]-sumwx_F2[2]/n0[6]+sumwx_F2[3]/n0[7]

    while (aaa1>0.0001){
      n_iter<- n_iter+1
      abc1<- sigma1D[1]/n0[1]+sigma1D[2]/n0[2]+sigma3D[1]/n0[5]+sigma3D[2]/n0[6]
      abc2<- -sigma3D[2]/n0[6]
      abc3<- sigma2D[1]/n0[3]+sigma2D[2]/n0[4]+sigma3D[2]/n0[6]+sigma3D[3]/n0[7]
      aa2<- abc1*abc3-abc2*abc2; aa3<- s0[1]*abc3-s0[2]*abc2; aa4<- s0[2]*abc1-s0[1]*abc2
      rr<- matrix(0,2,1)
      rr[1]<- aa3/aa2;rr[2]<- aa4/aa2
      mean1D[1]<- (sumwx_B1[1]-rr[1]*sigma1D[1])/n0[1]
      mean1D[2]<- (sumwx_B1[2]+rr[1]*sigma1D[2])/n0[2]
      mean2D[1]<- (sumwx_B2[1]-rr[2]*sigma2D[1])/n0[3]
      mean2D[2]<- (sumwx_B2[2]+rr[2]*sigma2D[2])/n0[4]
      mean3D[1]<- (sumwx_F2[1]+rr[1]*sigma3D[1])/n0[5]
      mean3D[2]<- (sumwx_F2[2]+sigma3D[2]*(-rr[1]+rr[2]))/n0[6]
      mean3D[3]<- (sumwx_F2[3]-rr[2]*sigma3D[3])/n0[7]

      aaa1<- max(abs(AA-rr))
      AA<- rr
      if (n_iter>20) break
    }
    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2)
    ss2<- sum((dataF1-meanF1)^2)
    ss3<- sum((dataP2-meanP2)^2)

    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1D[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2D[i])^2 }
    for(i in 1:d23) {  swx_F2[i] <- WW_F2[i,]%*%(dataF2-mean3D[i])^2 }

    sigma1D<- matrix(sum(swx_B1)/n_samB1,d21,1)
    sigma40<- sigma1D[1]-sigma
    sigma2D<- matrix(sum(swx_B2)/n_samB2,d22,1)
    sigma50<- sigma2D[1]-sigma
    sigma3D<- matrix(sum(swx_F2)/n_samF2,d23,1)
    sigma60<- sigma3D[1]-sigma
    ################ iteratively CM3-step for variance (sigma) ################
    sigma40[(sigma40<0)]<- 0.00001; sigma50[(sigma50<0)]<- 0.00001; sigma60[(sigma60<0)]<- 0.00001
    s0[3]<- ss1+ss2+ss3;s0[4]<- n_samP1+n_samF1+n_samP2
    aaa0<- 0;n_iter<- 0;aa3<- 1000
    while (aa3>0.0001){
      n_iter<- n_iter+1
      abc1<- sigma/(sigma+sigma40);abc2<- sigma/(sigma+sigma50);abc3<- sigma/(sigma+sigma60)
      aa4<- s0[3]+abc1*abc1*sum(swx_B1)+abc2*abc2*sum(swx_B2)+abc3*abc3*sum(swx_F2)
      aa5<- s0[4]+abc1*n_samB1+abc2*n_samB2+abc3*n_samF2
      sigma<- aa4/aa5;aa3<- abs(sigma-aaa0);aaa0<- sigma
      if (n_iter>20) break
    }
    sigma1D<- matrix(sigma+sigma40,d21,1);sigma2D<- matrix(sigma+sigma50,d22,1);sigma3D<- matrix(sigma+sigma60,d23,1)
    ####################### the stop criterion for iteration #########################################
    L1<- logL(n_samP1,1,1,meanP1,sigma,dataP1)+logL(n_samF1,1,1,meanF1,sigma,dataF1)+logL(n_samP2,1,1,meanP2,sigma,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1D,sigma1D,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2D,sigma2D,dataB2)+logL(n_samF2,d23,mix_pi_3,mean3D,sigma3D,dataF2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*12

  #########first order genetic parameter process##########
  aa<- matrix(c(1,0,0,0,0,0,0,0,0,0, 0,1,0,0,0,0,0,0,0,0, 0,0,1,0,0,0,0,0,0,0, 0,0,0,1,1,0,0,0,0,0,
                0,0,0,0,0,1,1,0,0,0, 0,0,0,0,0,0,0,1,1,1, 1,0,-1,1,0,0,-1,1,0,-1, 0,1,0,0,1,1,0,0,1,0),10,8)
  b_line1 <- matrix(c(meanP1,meanF1,meanP2,mean1D,mean2D,mean3D))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigmaB1 - sigma1D[1]
  if(jj_1 < 0 | jj_1>=sigmaB1) {jj_1 <- 0}
  ll_1 <- jj_1/sigmaB1
  mm1<- sigma1D[1]-sigma
  if (mm1<0 | mm1>=sigmaB1) { mm1<- 0 }
  nn1<- mm1/sigmaB1

  jj_2 <- sigmaB2 - sigma2D[1]
  if(jj_2 < 0 | jj_2>=sigmaB2) {jj_2 <- 0}
  ll_2 <- jj_2/sigmaB2
  mm2<- sigma2D[1]-sigma
  if (mm2<0 | mm2>=sigmaB2) { mm2<- 0 }
  nn2<- mm2/sigmaB2

  jj_3 <- sigmaF2 - sigma3D[1]
  if(jj_3 < 0 | jj_1>=sigmaF2) {jj_3 <- 0}
  ll_3 <- jj_3/sigmaF2
  mm3<- sigma3D[1]-sigma
  if (mm3<0 | mm3>=sigmaF2) { mm3<- 0 }
  nn3<- mm3/sigmaF2


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
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

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
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),(1-pkolm(D_F1,n_samF1))))

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
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

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
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),(1-pkolm(D_B1,n_samB1))))

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
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),(1-pkolm(D_B2,n_samB2))))

  ###############  F2  ###################
  dataF2<-sort(dataF2);bmw_F2 <- matrix(0,n_samF2,1); bmwsl_F2 <- matrix(0,n_samF2,d23)
  for(i in 1:d23){
    gg_F2 <- (dataF2 - mean3D[i])/sqrt(sigma3D[i])
    bmw_F2[which(gg_F2>=0)] <- pnorm(gg_F2[gg_F2>=0])
    bmw_F2[which(gg_F2<0)] <- 1 - pnorm(abs(gg_F2[gg_F2<0]))
    bmwsl_F2[,i] <- bmw_F2*mix_pi_3[i]
  }
  P2_F2 <- rowSums(bmwsl_F2)
  nn<-dim(as.matrix(unique(P2_F2)))[1]
  if(nn<n_samF2){P2_F2<-P2_F2+runif(n_samF2)/1e4}

  dd_F2 <- as.matrix(c(sum(P2_F2),sum(P2_F2^2),sum((P2_F2-0.5)^2)))
  WW2_F2 <- 1/(12*n_samF2) + sum((P2_F2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  u_F2 <- as.matrix(c(12*n_samF2*((dd_F2[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((dd_F2[2]/n_samF2-1/3)^2),180*n_samF2*((dd_F2[3]/n_samF2-1/12)^2)))
  D_F2 <- as.numeric(ks.test(P2_F2,"punif")[[1]][1])
  tt_F2 <- as.matrix(c((1 - pchisq(u_F2[1],1)),(1 - pchisq(u_F2[2],1)),(1 - pchisq(u_F2[3],1)),K1(WW2_F2),(1-pkolm(D_F2,n_samF2))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)
  tt_F2[which(tt_F2>=10e-4)]<-round(tt_F2[which(tt_F2>=10e-4)],4);tt_F2[which(tt_F2<10e-4)]<-format(tt_F2[which(tt_F2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX1-AD-ADI",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigma,4),round(t(mean1D),4)," "," ",round(sigma1D[1],4),round(t(mix_pi_1),4)," "," ",
                       round(t(mean2D),4)," "," ",round(sigma2D[1],4),round(t(mix_pi_2),4)," "," ",round(t(mean3D),4)," "," "," "," "," "," ",round(sigma3D[1],4),round(t(mix_pi_3),4),
                       " "," "," "," "," "," ",round(B1[1],4),round(B1[2],4),round(B1[3],4),round(B1[4],4),round(B1[5],4),round(B1[6],4),round(B1[7],4)," ",round(B1[8],4)," "," "," "," "," "," "," ",
                       round(jj_1,4),round(ll_1*100,4),round(mm1,4),round(nn1*100,4),round(jj_2,4),round(ll_2*100,4),round(mm2,4),round(nn2*100,4),round(jj_3,4),round(ll_3*100,4),round(mm3,4),round(nn3*100,4),round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),
                       tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),
                       tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),
                       tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5],round(u_F2[1],4),tt_F2[1],round(u_F2[2],4),
                       tt_F2[2],round(u_F2[3],4),tt_F2[3],round(WW2_F2,4),tt_F2[4],round(D_F2,4),tt_F2[5])

  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2,mi_3)
  return(OUTPUT)
}
############################################ MX1-AD-AD Model #########################################
G6ModelFun[[14]] <- function(K1,logL,df11,df21,df31,df41,df51,df61){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1]; n_samF1<-dim(dataF1)[1]; n_samP2<-dim(dataP2)[1]; n_samB1<-dim(dataB1)[1]; n_samB2<-dim(dataB2)[1]; n_samF2<-dim(dataF2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2)
  sumx4<-sum(dataB1);sumx5<-sum(dataB2);sumx6<-sum(dataF2)
  sigmaP1<- as.numeric(var(dataP1)); ss1<- (n_samP1-1)*sigmaP1
  sigmaF1<- as.numeric(var(dataF1)); ss2<- (n_samF1-1)*sigmaF1
  sigmaP2<- as.numeric(var(dataP2)); ss3<- (n_samP2-1)*sigmaP2
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  mean4<- mean(dataB1); sigmaB1<- as.numeric(var(dataB1)); sigma40<- sigmaB1
  mean5<- mean(dataB2); sigmaB2<- as.numeric(var(dataB2)); sigma50<- sigmaB2
  mean6<- mean(dataF2); sigmaF2<- as.numeric(var(dataF2)); sigma60<- sigmaF2

  m_esp<-0.0001
  ####################### MX1-AD-AD Model ###################### (D1) ##############
  d21<-2; d22<-2; d23<-3
  sigma<- sigma0
  meanP1<- mean(dataP1);meanF1<- mean(dataF1);meanP2<- mean(dataP2)
  abb<- sigmaF2/(5.0*sigma0)

  mi_1<- matrix(0.5,d21,1)
  sigma4<- sigmaB1/abb ;sigma1D<- matrix(sigma4,d21,1)
  a1<- sqrt(sigmaB1/n_samB1)
  if (meanP1<meanP2) {a1= -a1}
  mean1D<- as.matrix(c(mean4+2*a1,mean4))

  mi_2<- matrix(0.5,d22,1)
  sigma5<- sigmaB2/abb ;sigma2D<- matrix(sigma5,d22,1)
  a2<- sqrt(sigmaB2/n_samB2)
  if (meanP1<meanP2) { a2= -a2 }
  mean2D<- as.matrix(c(mean5,mean5-2*a2))

  mi_3<- as.matrix(c(0.25,0.5,0.25))
  sigma6<- sigmaF2/abb;sigma3D<- matrix(sigma6,d23,1)
  a3<- sqrt(sigmaF2/n_samF2)
  if (meanP1<meanP2) { a3= -a3 }
  mean3D<- as.matrix(c(mean6+2*a3,mean6,mean6-2*a3))

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1D,sigma1D,dataB1)+logL(n_samB2,d22,mi_2,mean2D,sigma2D,dataB2)+logL(n_samF2,d23,mi_3,mean3D,sigma3D,dataF2)

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)
  WW_F2 <- matrix(0,d23,n_samF2); swx_F2 <- matrix(0,d23,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1D[i],sqrt(sigma1D[i]))/dmixnorm(dataB1,mean1D,sqrt(sigma1D),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2D[i],sqrt(sigma2D[i]))/dmixnorm(dataB2,mean2D,sqrt(sigma2D),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    for(i in 1:d23) { WW_F2[i,] <- mi_3[i]*dnorm(dataF2,mean3D[i],sqrt(sigma3D[i]))/dmixnorm(dataF2,mean3D,sqrt(sigma3D),mi_3) }
    mix_pi_3 <- as.matrix(rowSums(WW_F2)/n_samF2)
    sumwx_F2 <- WW_F2%*%dataF2

    ############ CM1-step for means ##############
    n0<- as.matrix(c(mix_pi_1*n_samB1,mix_pi_2*n_samB2,mix_pi_3*n_samF2))
    n0[abs(n0)<0.0001] <- 0.0001
    aaa0<- 0; aaa1<- 1000; n_iter<- 0; AA<- matrix(0,5,1)
    while (aaa1>0.0001){
      n_iter<- n_iter+1
      hh<- matrix(0,5,5)
      hh[1,1]<- sigma*(1.0/n_samP1+1.0/n_samF1)+sigma1D[1]/n0[1]+sigma1D[2]/n0[2]
      hh[1,2]<- sigma/n_samF1
      hh[1,3]<- -sigma1D[1]/n0[1]+sigma1D[2]/n0[2]
      hh[1,4]<- 0
      hh[1,5]<- -sigma1D[2]/n0[2]
      hh[2,2]<- sigma*(1.0/n_samF1+1.0/n_samP2)+sigma2D[1]/n0[3]+sigma2D[2]/n0[4]
      hh[2,3]<- 0
      hh[2,4]<- -sigma2D[1]/n0[3]+sigma2D[2]/n0[4]
      hh[2,5]<- -sigma2D[1]/n0[3]
      hh[3,3]<- sigma1D[1]/n0[1]+sigma1D[2]/n0[2]+sigma3D[1]/n0[5]+sigma3D[2]/n0[6]
      hh[3,4]<- -sigma3D[2]/n0[6]
      hh[3,5]<- -sigma1D[2]/n0[2]-2.0*sigma3D[2]/n0[6]
      hh[4,4]<- sigma2D[1]/n0[3]+sigma2D[2]/n0[4]+sigma3D[2]/n0[6]+sigma3D[3]/n0[7]
      hh[4,5]<- sigma2D[1]/n0[5]+2.0*sigma3D[2]/n0[6]
      hh[5,5]<- sigma1D[2]/n0[2]+sigma2D[1]/n0[3]+4.0*sigma3D[2]/n0[6]
      for(i in 2:5){
        for(j in 1:(i-1)){
          hh[i,j]<- hh[j,i]
        }
      }
      b_line<- matrix(0,5,1)
      b_line[1]<- sumx1/n_samP1+sumx2/n_samF1-sumwx_B1[1]/n0[1]-sumwx_B1[2]/n0[2]
      b_line[2]<- sumx2/n_samF1+sumx3/n_samP2-sumwx_B2[1]/n0[3]-sumwx_B2[2]/n0[4]
      b_line[3]<- sumwx_B1[1]/n0[1]-sumwx_B1[2]/n0[2]-sumwx_F2[1]/n0[5]+sumwx_F2[2]/n0[6]
      b_line[4]<- sumwx_B2[1]/n0[3]-sumwx_B2[2]/n0[4]-sumwx_F2[2]/n0[6]+sumwx_F2[3]/n0[7]
      b_line[5]<- sumwx_B1[2]/n0[2]+sumwx_B2[1]/n0[3]-2.0*sumwx_F2[2]/n0[6]
      B <- solve(hh,b_line)

      meanP1<- (sumx1-sigma*B[1])/n_samP1
      meanF1<- (sumx2-sigma*(B[1]+B[2]))/n_samF1
      meanP2<- (sumx3-sigma*B[2])/n_samP2
      mean1D[1]<- (sumwx_B1[1]+sigma1D[1]*(B[1]-B[3]))/n0[1]
      mean1D[2]<- (sumwx_B1[2]+sigma1D[2]*(B[1]+B[3]-B[5]))/n0[2]
      mean2D[1]<- (sumwx_B2[1]+sigma2D[1]*(B[2]-B[4]-B[5]))/n0[3]
      mean2D[2]<- (sumwx_B2[2]+(B[2]+B[4])*sigma2D[2])/n0[4]
      mean3D[1]<- (sumwx_F2[1]+sigma3D[1]*B[3])/n0[5]
      mean3D[2]<- (sumwx_F2[2]+sigma3D[2]*(-B[3]+B[4]+2.0*B[5]))/n0[6]
      mean3D[3]<- (sumwx_F2[3]-sigma3D[3]*B[4])/n0[7]

      aaa1<- max(abs(AA-B))
      AA<- B
      if (n_iter>20) break
    }
    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2)
    ss2<- sum((dataF1-meanF1)^2)
    ss3<- sum((dataP2-meanP2)^2)

    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1D[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2D[i])^2 }
    for(i in 1:d23) {  swx_F2[i] <- WW_F2[i,]%*%(dataF2-mean3D[i])^2 }

    sigma1D<- matrix(sum(swx_B1)/n_samB1,d21,1)
    sigma40<- sigma1D[1]-sigma
    sigma2D<- matrix(sum(swx_B2)/n_samB2,d22,1)
    sigma50<- sigma2D[1]-sigma
    sigma3D<- matrix(sum(swx_F2)/n_samF2,d23,1)
    sigma60<- sigma3D[1]-sigma
    ################ iteratively CM3-step for variance (sigma) ################
    sigma40[(sigma40<0)]<- 0.00001; sigma50[(sigma50<0)]<- 0.00001; sigma60[(sigma60<0)]<- 0.00001
    s0<- matrix(0,2,1)
    s0[1]<- ss1+ss2+ss3;s0[2]<- n_samP1+n_samF1+n_samP2
    aaa0<- 0;n_iter<- 0;aa3<- 1000
    while (aa3>0.0001){
      n_iter<- n_iter+1
      abc1<- sigma/(sigma+sigma40);abc2<- sigma/(sigma+sigma50);abc3<- sigma/(sigma+sigma60)
      aa4<- s0[1]+abc1*abc1*sum(swx_B1)+abc2*abc2*sum(swx_B2)+abc3*abc3*sum(swx_F2)
      aa5<- s0[2]+abc1*n_samB1+abc2*n_samB2+abc3*n_samF2
      sigma<- aa4/aa5;aa3<- abs(sigma-aaa0);aaa0<- sigma
      if (n_iter>20) break
    }
    sigma1D<- matrix(sigma+sigma40,d21,1);sigma2D<- matrix(sigma+sigma50,d22,1);sigma3D<- matrix(sigma+sigma60,d23,1)
    ####################### the stop criterion for iteration #########################################
    L1<- logL(n_samP1,1,1,meanP1,sigma,dataP1)+logL(n_samF1,1,1,meanF1,sigma,dataF1)+logL(n_samP2,1,1,meanP2,sigma,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1D,sigma1D,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2D,sigma2D,dataB2)+logL(n_samF2,d23,mix_pi_3,mean3D,sigma3D,dataF2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*9

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1,1,1,1, 1,0,-1,1,0,0,-1,1,0,-1, 0,1,0,0,1,1,0,0,1,0, 1,0,-1,0.5,0.5,-0.5,-0.5,0,0,0,
                0,1,0,0.5,0.5,0.5,0.5,0.5,0.5,0.5),10,5)
  b_line1 <- matrix(c(meanP1,meanF1,meanP2,mean1D,mean2D,mean3D))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigmaB1 - sigma1D[1]
  if(jj_1 < 0 | jj_1>=sigmaB1) {jj_1 <- 0}
  ll_1 <- jj_1/sigmaB1
  mm1<- sigma1D[1]-sigma
  if (mm1<0 | mm1>=sigmaB1) { mm1<- 0 }
  nn1<- mm1/sigmaB1

  jj_2 <- sigmaB2 - sigma2D[1]
  if(jj_2 < 0 | jj_2>=sigmaB2) {jj_2 <- 0}
  ll_2 <- jj_2/sigmaB2
  mm2<- sigma2D[1]-sigma
  if (mm2<0 | mm2>=sigmaB2) { mm2<- 0 }
  nn2<- mm2/sigmaB2

  jj_3 <- sigmaF2 - sigma3D[1]
  if(jj_3 < 0 | jj_1>=sigmaF2) {jj_3 <- 0}
  ll_3 <- jj_3/sigmaF2
  mm3<- sigma3D[1]-sigma
  if (mm3<0 | mm3>=sigmaF2) { mm3<- 0 }
  nn3<- mm3/sigmaF2

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
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

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
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),(1-pkolm(D_F1,n_samF1))))

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
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

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
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),(1-pkolm(D_B1,n_samB1))))

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
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),(1-pkolm(D_B2,n_samB2))))

  ###############  F2  ###################
  dataF2<-sort(dataF2);bmw_F2 <- matrix(0,n_samF2,1); bmwsl_F2 <- matrix(0,n_samF2,d23)
  for(i in 1:d23){
    gg_F2 <- (dataF2 - mean3D[i])/sqrt(sigma3D[i])
    bmw_F2[which(gg_F2>=0)] <- pnorm(gg_F2[gg_F2>=0])
    bmw_F2[which(gg_F2<0)] <- 1 - pnorm(abs(gg_F2[gg_F2<0]))
    bmwsl_F2[,i] <- bmw_F2*mix_pi_3[i]
  }
  P2_F2 <- rowSums(bmwsl_F2)
  nn<-dim(as.matrix(unique(P2_F2)))[1]
  if(nn<n_samF2){P2_F2<-P2_F2+runif(n_samF2)/1e4}

  dd_F2 <- as.matrix(c(sum(P2_F2),sum(P2_F2^2),sum((P2_F2-0.5)^2)))
  WW2_F2 <- 1/(12*n_samF2) + sum((P2_F2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  u_F2 <- as.matrix(c(12*n_samF2*((dd_F2[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((dd_F2[2]/n_samF2-1/3)^2),180*n_samF2*((dd_F2[3]/n_samF2-1/12)^2)))
  D_F2 <- as.numeric(ks.test(P2_F2,"punif")[[1]][1])
  tt_F2 <- as.matrix(c((1 - pchisq(u_F2[1],1)),(1 - pchisq(u_F2[2],1)),(1 - pchisq(u_F2[3],1)),K1(WW2_F2),(1-pkolm(D_F2,n_samF2))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)
  tt_F2[which(tt_F2>=10e-4)]<-round(tt_F2[which(tt_F2>=10e-4)],4);tt_F2[which(tt_F2<10e-4)]<-format(tt_F2[which(tt_F2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX1-AD-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigma,4),round(t(mean1D),4)," "," ",round(sigma1D[1],4),round(t(mix_pi_1),4)," "," ",
                       round(t(mean2D),4)," "," ",round(sigma2D[1],4),round(t(mix_pi_2),4)," "," ",round(t(mean3D),4)," "," "," "," "," "," ",round(sigma3D[1],4),round(t(mix_pi_3),4),
                       " "," "," "," "," "," ",round(B1[1],4)," "," "," "," "," ",round(B1[2],4)," ",round(B1[3],4)," "," "," "," "," ",round(B1[4],4),round(B1[5],4),round(jj_1,4),round(ll_1*100,4),round(mm1,4),
                       round(nn1*100,4),round(jj_2,4),round(ll_2*100,4),round(mm2,4),round(nn2*100,4),round(jj_3,4),round(ll_3*100,4),round(mm3,4),round(nn3*100,4),round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),
                       tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),
                       tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),
                       tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5],round(u_F2[1],4),tt_F2[1],round(u_F2[2],4),
                       tt_F2[2],round(u_F2[3],4),tt_F2[3],round(WW2_F2,4),tt_F2[4],round(D_F2,4),tt_F2[5])

  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2,mi_3)
  return(OUTPUT)
}
############################################ MX1-A-AD Model #########################################
G6ModelFun[[15]] <- function(K1,logL,df11,df21,df31,df41,df51,df61){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1]; n_samF1<-dim(dataF1)[1]; n_samP2<-dim(dataP2)[1]; n_samB1<-dim(dataB1)[1]; n_samB2<-dim(dataB2)[1]; n_samF2<-dim(dataF2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2)
  sumx4<-sum(dataB1);sumx5<-sum(dataB2);sumx6<-sum(dataF2)
  sigmaP1<- as.numeric(var(dataP1)); ss1<- (n_samP1-1)*sigmaP1
  sigmaF1<- as.numeric(var(dataF1)); ss2<- (n_samF1-1)*sigmaF1
  sigmaP2<- as.numeric(var(dataP2)); ss3<- (n_samP2-1)*sigmaP2
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  mean4<- mean(dataB1); sigmaB1<- as.numeric(var(dataB1)); sigma40<- sigmaB1
  mean5<- mean(dataB2); sigmaB2<- as.numeric(var(dataB2)); sigma50<- sigmaB2
  mean6<- mean(dataF2); sigmaF2<- as.numeric(var(dataF2)); sigma60<- sigmaF2

  m_esp<-0.0001
  ####################### MX1-A-AD Model #################### (D2) ################
  d21<-2; d22<-2; d23<-3
  sigma<- sigma0
  meanP1<- mean(dataP1);meanF1<- mean(dataF1);meanP2<- mean(dataP2)
  abb<- sigmaF2/(5.0*sigma0)

  mi_1<- matrix(0.5,d21,1);sigma4<- sigmaB1/abb ;sigma1D<- matrix(sigma4,d21,1)
  a1<- sqrt(sigmaB1)
  if (meanP1<meanP2) {a1= -a1}
  mean1D<- as.matrix(c(mean4+2*a1,mean4))

  mi_2<- matrix(0.5,d22,1);sigma5<- sigmaB2/abb ;sigma2D<- matrix(sigma5,d22,1)
  a2<- sqrt(sigmaB2)
  if (meanP1<meanP2) { a2= -a2 }
  mean2D<- as.matrix(c(mean5,mean5-2*a2))

  mi_3<- as.matrix(c(0.25,0.5,0.25));sigma6<- sigmaF2/abb;sigma3D<- matrix(sigma6,d23,1)
  a3<- sqrt(sigmaF2)
  if (meanP1<meanP2) { a3= -a3 }
  mean3D<- as.matrix(c(mean6+2*a3,mean6,mean6-2*a3))

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1D,sigma1D,dataB1)+logL(n_samB2,d22,mi_2,mean2D,sigma2D,dataB2)+logL(n_samF2,d23,mi_3,mean3D,sigma3D,dataF2)

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)
  WW_F2 <- matrix(0,d23,n_samF2); swx_F2 <- matrix(0,d23,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1D[i],sqrt(sigma1D[i]))/dmixnorm(dataB1,mean1D,sqrt(sigma1D),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2D[i],sqrt(sigma2D[i]))/dmixnorm(dataB2,mean2D,sqrt(sigma2D),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    for(i in 1:d23) { WW_F2[i,] <- mi_3[i]*dnorm(dataF2,mean3D[i],sqrt(sigma3D[i]))/dmixnorm(dataF2,mean3D,sqrt(sigma3D),mi_3) }
    mix_pi_3 <- as.matrix(rowSums(WW_F2)/n_samF2)
    sumwx_F2 <- WW_F2%*%dataF2

    ############ CM1-step for means ##############
    n0<- as.matrix(c(mix_pi_1*n_samB1,mix_pi_2*n_samB2,mix_pi_3*n_samF2))
    n0[abs(n0)<0.0001] <- 0.0001
    aaa0<- 0; aaa1<- 1000; n_iter<- 0; AA<- matrix(0,6,1)
    while (aaa1>0.0001){
      n_iter<- n_iter+1
      hh<- matrix(0,6,6)
      hh[1,1]<- sigma*(1.0/n_samP1+1.0/n_samF1)+sigma1D[1]/n0[1]+sigma1D[2]/n0[2]
      hh[1,2]<- sigma/n_samF1
      hh[1,3]<- -sigma1D[1]/n0[1]+sigma1D[2]/n0[2]
      hh[1,4]<- 0
      hh[1,5]<- -sigma1D[2]/n0[2]
      hh[1,6]<- sigma/n_samP1+2.0*sigma1D[1]/n0[1]
      hh[2,2]<- sigma*(1.0/n_samF1+1.0/n_samP2)+sigma2D[1]/n0[3]+sigma2D[2]/n0[4]
      hh[2,3]<- 0
      hh[2,4]<- -sigma2D[1]/n0[3]+sigma2D[2]/n0[4]
      hh[2,5]<- -sigma2D[1]/n0[3]
      hh[2,6]<- -(sigma/n_samP2+2.0*sigma2D[1]/n0[3])
      hh[3,3]<- sigma1D[1]/n0[1]+sigma1D[2]/n0[2]+sigma3D[1]/n0[5]+sigma3D[2]/n0[6]
      hh[3,4]<- -sigma3D[2]/n0[6]
      hh[3,5]<- -(sigma1D[2]/n0[2]+2*sigma3D[2]/n0[6])
      hh[3,6]<- -2.0*sigma1D[1]/n0[1]
      hh[4,4]<- sigma2D[1]/n0[3]+sigma2D[2]/n0[4]+sigma3D[2]/n0[6]+sigma3D[3]/n0[7]
      hh[4,5]<- sigma2D[1]/n0[3]+2.0*sigma3D[2]/n0[6]
      hh[4,6]<- 2.0*sigma2D[1]/n0[3]
      hh[5,5]<- sigma1D[2]/n0[2]+sigma2D[1]/n0[3]+4.0*sigma3D[2]/n0[6]
      hh[5,6]<- 2.0*sigma2D[1]/n0[3]
      hh[6,6]<- sigma*(1.0/n_samP1+1.0/n_samP2)+4.0*sigma1D[1]/n0[1]+4.0*sigma2D[1]/n0[3]
      for(i in 2:6){
        for(j in 1:(i-1)){
          hh[i,j]<- hh[j,i]
        }
      }
      b_line<- matrix(0,6,1)
      b_line[1]<- sumx1/n_samP1+sumx2/n_samF1-sumwx_B1[1]/n0[1]-sumwx_B1[2]/n0[2]
      b_line[2]<- sumx2/n_samF1+sumx3/n_samP2-sumwx_B2[1]/n0[3]-sumwx_B2[2]/n0[4]
      b_line[3]<- sumwx_B1[1]/n0[1]-sumwx_B1[2]/n0[2]-sumwx_F2[1]/n0[5]+sumwx_F2[2]/n0[6]
      b_line[4]<- sumwx_B2[1]/n0[3]-sumwx_B2[2]/n0[4]-sumwx_F2[2]/n0[6]+sumwx_F2[3]/n0[7]
      b_line[5]<- sumwx_B1[2]/n0[2]+sumwx_B2[1]/n0[3]-2.0*sumwx_F2[2]/n0[6]
      b_line[6]<- sumx1/n_samP1-sumx3/n_samP2-2.0*sumwx_B1[1]/n0[1]+2.0*sumwx_B2[1]/n0[3]
      B <- solve(hh,b_line)

      meanP1<- (sumx1-sigma*(B[1]+B[6]))/n_samP1
      meanF1<- (sumx2-sigma*(B[1]+B[2]))/n_samF1
      meanP2<- (sumx3+sigma*(-B[2]+B[6]))/n_samP2
      mean1D[1]<- (sumwx_B1[1]+sigma1D[1]*(B[1]-B[3]+2.0*B[6]))/n0[1]
      mean1D[2]<- (sumwx_B1[2]+sigma1D[2]*(B[1]+B[3]-B[5]))/n0[2]
      mean2D[1]<- (sumwx_B2[1]+sigma2D[1]*(B[2]-B[4]-B[5]-2.0*B[6]))/n0[3]
      mean2D[2]<- (sumwx_B2[2]+(B[2]+B[4])*sigma2D[2])/n0[4]
      mean3D[1]<- (sumwx_F2[1]+sigma3D[1]*B[3])/n0[5]
      mean3D[2]<- (sumwx_F2[2]+sigma3D[2]*(-B[3]+B[4]+2.0*B[5]))/n0[6]
      mean3D[3]<- (sumwx_F2[3]-sigma3D[3]*B[4])/n0[7]

      aaa1<- max(abs(AA-B))
      AA<- B
      if (n_iter>20) break
    }
    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2)
    ss2<- sum((dataF1-meanF1)^2)
    ss3<- sum((dataP2-meanP2)^2)

    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1D[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2D[i])^2 }
    for(i in 1:d23) {  swx_F2[i] <- WW_F2[i,]%*%(dataF2-mean3D[i])^2 }

    sigma1D<- matrix(sum(swx_B1)/n_samB1,d21,1)
    sigma40<- sigma1D[1]-sigma
    sigma2D<- matrix(sum(swx_B2)/n_samB2,d22,1)
    sigma50<- sigma2D[1]-sigma
    sigma3D<- matrix(sum(swx_F2)/n_samF2,d23,1)
    sigma60<- sigma3D[1]-sigma
    ################ iteratively CM3-step for variance (sigma) ################
    sigma40[(sigma40<0)]<- 0.00001; sigma50[(sigma50<0)]<- 0.00001; sigma60[(sigma60<0)]<- 0.00001
    s0<- matrix(0,2,1)
    s0[1]<- ss1+ss2+ss3; s0[2]<- n_samP1+n_samF1+n_samP2
    aaa0<- 0;n_iter<- 0;aa3<- 1000
    while (aa3>0.0001){
      n_iter=n_iter+1
      abc1<- sigma/(sigma+sigma40);abc2<- sigma/(sigma+sigma50);abc3<- sigma/(sigma+sigma60)
      aa4<- s0[1]+abc1*abc1*sum(swx_B1)+abc2*abc2*sum(swx_B2)+abc3*abc3*sum(swx_F2)
      aa5<- s0[2]+abc1*n_samB1+abc2*n_samB2+abc3*n_samF2
      sigma<- aa4/aa5;aa3<- abs(sigma-aaa0);aaa0<- sigma
      if (n_iter>20) break
    }
    sigma1D<- matrix(sigma+sigma40,d21,1);sigma2D<- matrix(sigma+sigma50,d22,1);sigma3D<- matrix(sigma+sigma60,d23,1)
    ####################### the stop criterion for iteration #########################################
    L1<- logL(n_samP1,1,1,meanP1,sigma,dataP1)+logL(n_samF1,1,1,meanF1,sigma,dataF1)+logL(n_samP2,1,1,meanP2,sigma,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1D,sigma1D,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2D,sigma2D,dataB2)+logL(n_samF2,d23,mix_pi_3,mean3D,sigma3D,dataF2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*8

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1,1,1,1, 1,0,-1,1,0,0,-1,1,0,-1, 1,0,-1,0.5,0.5,-0.5,-0.5,0,0,0, 0,1,0,0.5,0.5,0.5,0.5,0.5,0.5,0.5),10,4)
  b_line1 <- matrix(c(meanP1,meanF1,meanP2,mean1D,mean2D,mean3D))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigmaB1 - sigma1D[1]
  if(jj_1 < 0 | jj_1>=sigmaB1) {jj_1 <- 0}
  ll_1 <- jj_1/sigmaB1
  mm1<- sigma1D[1]-sigma
  if (mm1<0 | mm1>=sigmaB1) { mm1<- 0 }
  nn1<- mm1/sigmaB1

  jj_2 <- sigmaB2 - sigma2D[1]
  if(jj_2 < 0 | jj_2>=sigmaB2) {jj_2 <- 0}
  ll_2 <- jj_2/sigmaB2
  mm2<- sigma2D[1]-sigma
  if (mm2<0 | mm2>=sigmaB2) { mm2<- 0 }
  nn2<- mm2/sigmaB2

  jj_3 <- sigmaF2 - sigma3D[1]
  if(jj_3 < 0 | jj_1>=sigmaF2) {jj_3 <- 0}
  ll_3 <- jj_3/sigmaF2
  mm3<- sigma3D[1]-sigma
  if (mm3<0 | mm3>=sigmaF2) { mm3<- 0 }
  nn3<- mm3/sigmaF2

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
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

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
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),(1-pkolm(D_F1,n_samF1))))

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
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

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
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),(1-pkolm(D_B1,n_samB1))))

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
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),(1-pkolm(D_B2,n_samB2))))

  ###############  F2  ###################
  dataF2<-sort(dataF2);bmw_F2 <- matrix(0,n_samF2,1); bmwsl_F2 <- matrix(0,n_samF2,d23)
  for(i in 1:d23){
    gg_F2 <- (dataF2 - mean3D[i])/sqrt(sigma3D[i])
    bmw_F2[which(gg_F2>=0)] <- pnorm(gg_F2[gg_F2>=0])
    bmw_F2[which(gg_F2<0)] <- 1 - pnorm(abs(gg_F2[gg_F2<0]))
    bmwsl_F2[,i] <- bmw_F2*mix_pi_3[i]
  }
  P2_F2 <- rowSums(bmwsl_F2)
  nn<-dim(as.matrix(unique(P2_F2)))[1]
  if(nn<n_samF2){P2_F2<-P2_F2+runif(n_samF2)/1e4}

  dd_F2 <- as.matrix(c(sum(P2_F2),sum(P2_F2^2),sum((P2_F2-0.5)^2)))
  WW2_F2 <- 1/(12*n_samF2) + sum((P2_F2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  u_F2 <- as.matrix(c(12*n_samF2*((dd_F2[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((dd_F2[2]/n_samF2-1/3)^2),180*n_samF2*((dd_F2[3]/n_samF2-1/12)^2)))
  D_F2 <- as.numeric(ks.test(P2_F2,"punif")[[1]][1])
  tt_F2 <- as.matrix(c((1 - pchisq(u_F2[1],1)),(1 - pchisq(u_F2[2],1)),(1 - pchisq(u_F2[3],1)),K1(WW2_F2),(1-pkolm(D_F2,n_samF2))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)
  tt_F2[which(tt_F2>=10e-4)]<-round(tt_F2[which(tt_F2>=10e-4)],4);tt_F2[which(tt_F2<10e-4)]<-format(tt_F2[which(tt_F2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX1-A-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigma,4),round(t(mean1D),4)," "," ",round(sigma1D[1],4),round(t(mix_pi_1),4)," "," ",
                       round(t(mean2D),4)," "," ",round(sigma2D[1],4),round(t(mix_pi_2),4)," "," ",round(t(mean3D),4)," "," "," "," "," "," ",round(sigma3D[1],4),round(t(mix_pi_3),4),
                       " "," "," "," "," "," ",round(B1[1],4)," "," "," "," "," ",round(B1[2],4)," "," "," "," "," "," "," ",round(B1[3],4),round(B1[4],4),round(jj_1,4),round(ll_1*100,4),round(mm1,4),
                       round(nn1*100,4),round(jj_2,4),round(ll_2*100,4),round(mm2,4),round(nn2*100,4),round(jj_3,4),round(ll_3*100,4),round(mm3,4),round(nn3*100,4),round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),
                       tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),
                       tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),
                       tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5],round(u_F2[1],4),tt_F2[1],round(u_F2[2],4),
                       tt_F2[2],round(u_F2[3],4),tt_F2[3],round(WW2_F2,4),tt_F2[4],round(D_F2,4),tt_F2[5])

  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2,mi_3)
  return(OUTPUT)
}
############################################ MX1-EAD-AD Model #########################################
G6ModelFun[[16]] <- function(K1,logL,df11,df21,df31,df41,df51,df61){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1]; n_samF1<-dim(dataF1)[1]; n_samP2<-dim(dataP2)[1]; n_samB1<-dim(dataB1)[1]; n_samB2<-dim(dataB2)[1]; n_samF2<-dim(dataF2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2)
  sumx4<-sum(dataB1);sumx5<-sum(dataB2);sumx6<-sum(dataF2)
  sigmaP1<- as.numeric(var(dataP1)); ss1<- (n_samP1-1)*sigmaP1
  sigmaF1<- as.numeric(var(dataF1)); ss2<- (n_samF1-1)*sigmaF1
  sigmaP2<- as.numeric(var(dataP2)); ss3<- (n_samP2-1)*sigmaP2
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  mean4<- mean(dataB1); sigmaB1<- as.numeric(var(dataB1)); sigma40<- sigmaB1
  mean5<- mean(dataB2); sigmaB2<- as.numeric(var(dataB2)); sigma50<- sigmaB2
  mean6<- mean(dataF2); sigmaF2<- as.numeric(var(dataF2)); sigma60<- sigmaF2

  m_esp<-0.0001
  ####################### MX1-EAD-AD Model ##################### (D3) ###################
  d21<-1; d22<-2; d23<-2
  sigma<- sigma0
  meanP1<- mean(dataP1);meanF1<- mean(dataF1);meanP2<- mean(dataP2)
  abb<- sigmaF2/(5.0*sigma0)

  mi_1<- as.matrix(1); sigma1D<- as.matrix(sigmaB1/abb)
  a1<- sqrt(sigmaB1/n_samB1)
  if (meanP1<meanP2) {a1= -a1}
  mean1D<- as.matrix(mean4+2*a1)

  mi_2<- matrix(0.5,d22,1);sigma2D<- matrix(sigmaB2/abb,d22,1)
  a2<- sqrt(sigmaB2/n_samB2)
  if (meanP1<meanP2) { a2= -a2 }
  mean2D<- as.matrix(c(mean5,mean5-2*a2))

  mi_3<- as.matrix(c(0.75,0.25)); sigma3D<- matrix(sigmaF2/abb,d23,1)
  mean3D<- as.matrix(c(mean1D,mean6))

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1D,sigma1D,dataB1)+logL(n_samB2,d22,mi_2,mean2D,sigma2D,dataB2)+logL(n_samF2,d23,mi_3,mean3D,sigma3D,dataF2)

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)
  WW_F2 <- matrix(0,d23,n_samF2); swx_F2 <- matrix(0,d23,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1D[i],sqrt(sigma1D[i]))/dmixnorm(dataB1,mean1D,sqrt(sigma1D),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2D[i],sqrt(sigma2D[i]))/dmixnorm(dataB2,mean2D,sqrt(sigma2D),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    for(i in 1:d23) { WW_F2[i,] <- mi_3[i]*dnorm(dataF2,mean3D[i],sqrt(sigma3D[i]))/dmixnorm(dataF2,mean3D,sqrt(sigma3D),mi_3) }
    mix_pi_3 <- as.matrix(rowSums(WW_F2)/n_samF2)
    sumwx_F2 <- WW_F2%*%dataF2

    ############ CM1-step for means ##############
    n0<- as.matrix(c(mix_pi_1*n_samB1,mix_pi_1*n_samB1,mix_pi_2*n_samB2,mix_pi_3*n_samF2))
    n0[abs(n0)<0.0001] <- 0.0001
    aaa0<- 0;aaa1<- 1000;n_iter<- 0;  AA<- matrix(0,4,1)
    while (aaa1>0.0001){
      n_iter<- n_iter+1
      hh<- matrix(0,4,4)
      hh[1,1]<- sigma*(1.0/n_samP1+1.0/n_samF1)+4.0*sigma1D[1]/n_samB1
      hh[1,2]<- sigma/n_samF1
      hh[1,3]<- 0
      hh[1,4]<- -2.0*sigma1D[1]/n_samB1
      hh[2,2]<- sigma*(1.0/n_samF1+1.0/n_samP2)+sigma2D[1]/n0[3]+sigma2D[2]/n0[4]
      hh[2,3]<- -sigma2D[1]/n0[3]+sigma2D[2]/n0[4]
      hh[2,4]<- -sigma2D[1]/n0[3]
      hh[3,3]<- sigma2D[1]/n0[3]+sigma2D[2]/n0[4]+sigma3D[1]/n0[5]+sigma3D[2]/n0[6]
      hh[3,4]<- sigma2D[1]/n0[3]+2.0*sigma3D[1]/n0[5]
      hh[4,4]<- sigma1D[1]/n0[1]+sigma2D[1]/n0[3]+4.0*sigma3D[1]/n0[5]
      for(i in 2:4){
        for(j in 1:(i-1)){
          hh[i,j]<- hh[j,i]
        }
      }
      b_line<- matrix(0,4,1)
      b_line[1]<- sumx1/n_samP1+sumx2/n_samF1-2.0*sumwx_B1[1]/n_samB1
      b_line[2]<- sumx2/n_samF1+sumx3/n_samP2-sumwx_B2[1]/n0[3]-sumwx_B2[2]/n0[4]
      b_line[3]<- sumwx_B2[1]/n0[3]-sumwx_B2[2]/n0[4]-sumwx_F2[1]/n0[5]+sumwx_F2[2]/n0[6]
      b_line[4]<- sumwx_B1[1]/n_samB1+sumwx_B2[1]/n0[3]-2.0*sumwx_F2[1]/n0[5]
      B <- solve(hh,b_line)

      meanP1<- (sumx1-sigma*B[1])/n_samP1
      meanF1<- (sumx2-sigma*(B[1]+B[2]))/n_samF1
      meanP2<- (sumx3-sigma*B[2])/n_samP2
      mean1D[1]<- (sumwx_B1[1]+sigma1D[1]*(2.0*B[1]-B[4]))/n0[1]
      mean2D[1]<- (sumwx_B2[1]+sigma2D[1]*(B[2]-B[3]-B[4]))/n0[3]
      mean2D[2]<- (sumwx_B2[2]+(B[2]+B[3])*sigma2D[2])/n0[4]
      mean3D[1]<- (sumwx_F2[1]+sigma3D[1]*(B[3]+2.0*B[4]))/n0[5]
      mean3D[2]<- (sumwx_F2[2]-sigma3D[2]*B[3])/n0[6]

      aaa1<- max(abs(AA-B))
      AA<- B
      if (n_iter>20) break
    }
    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2)
    ss2<- sum((dataF1-meanF1)^2)
    ss3<- sum((dataP2-meanP2)^2)

    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1D[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2D[i])^2 }
    for(i in 1:d23) {  swx_F2[i] <- WW_F2[i,]%*%(dataF2-mean3D[i])^2 }

    sigma1D<- matrix(sum(swx_B1)/n_samB1,d21,1)
    sigma40<- sigma1D[1]-sigma
    sigma2D<- matrix(sum(swx_B2)/n_samB2,d22,1)
    sigma50<- sigma2D[1]-sigma
    sigma3D<- matrix(sum(swx_F2)/n_samF2,d23,1)
    sigma60<- sigma3D[1]-sigma
    ################ iteratively CM3-step for variance (sigma) ################
    sigma40[(sigma40<0)]<- 0.00001; sigma50[(sigma50<0)]<- 0.00001; sigma60[(sigma60<0)]<- 0.00001
    s0<- matrix(0,2,1)
    s0[1]<- ss1+ss2+ss3; s0[2]<- n_samP1+n_samF1+n_samP2
    aaa0<- 0;n_iter<- 0;aa3<- 1000
    while (aa3>0.0001){
      n_iter<- n_iter+1
      abc1<- sigma/(sigma+sigma40);abc2<- sigma/(sigma+sigma50);abc3<- sigma/(sigma+sigma60)
      aa4<- s0[1]+abc1*abc1*sum(swx_B1)+abc2*abc2*sum(swx_B2)+abc3*abc3*sum(swx_F2)
      aa5<- s0[2]+abc1*n_samB1+abc2*n_samB2+abc3*n_samF2
      sigma<- aa4/aa5;aa3<- abs(sigma-aaa0);aaa0<- sigma
      if (n_iter>20) break
    }
    sigma1D<- matrix(sigma+sigma40,d21,1);sigma2D<- matrix(sigma+sigma50,d22,1);sigma3D<- matrix(sigma+sigma60,d23,1)
    ####################### the stop criterion for iteration #########################################
    L1<- logL(n_samP1,1,1,meanP1,sigma,dataP1)+logL(n_samF1,1,1,meanF1,sigma,dataF1)+logL(n_samP2,1,1,meanP2,sigma,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1D,sigma1D,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2D,sigma2D,dataB2)+logL(n_samF2,d23,mix_pi_3,mean3D,sigma3D,dataF2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*8

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1,1, 1,1,-1,1,1,-1,1,-1, 1,0,-1,0.5,-0.5,-0.5,0,0, 0,1,0,0.5,0.5,0.5,0.5,0.5),8,4)
  b_line1 <- matrix(c(meanP1,meanF1,meanP2,mean1D,mean2D,mean3D))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigmaB1 - sigma1D[1]
  if(jj_1 < 0 | jj_1>=sigmaB1) {jj_1 <- 0}
  ll_1 <- jj_1/sigmaB1
  mm1<- sigma1D[1]-sigma
  if (mm1<0 | mm1>=sigmaB1) { mm1<- 0 }
  nn1<- mm1/sigmaB1

  jj_2 <- sigmaB2 - sigma2D[1]
  if(jj_2 < 0 | jj_2>=sigmaB2) {jj_2 <- 0}
  ll_2 <- jj_2/sigmaB2
  mm2<- sigma2D[1]-sigma
  if (mm2<0 | mm2>=sigmaB2) { mm2<- 0 }
  nn2<- mm2/sigmaB2

  jj_3 <- sigmaF2 - sigma3D[1]
  if(jj_3 < 0 | jj_1>=sigmaF2) {jj_3 <- 0}
  ll_3 <- jj_3/sigmaF2
  mm3<- sigma3D[1]-sigma
  if (mm3<0 | mm3>=sigmaF2) { mm3<- 0 }
  nn3<- mm3/sigmaF2


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
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

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
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),(1-pkolm(D_F1,n_samF1))))

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
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

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
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),(1-pkolm(D_B1,n_samB1))))

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
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),(1-pkolm(D_B2,n_samB2))))

  ###############  F2  ###################
  dataF2<-sort(dataF2);bmw_F2 <- matrix(0,n_samF2,1); bmwsl_F2 <- matrix(0,n_samF2,d23)
  for(i in 1:d23){
    gg_F2 <- (dataF2 - mean3D[i])/sqrt(sigma3D[i])
    bmw_F2[which(gg_F2>=0)] <- pnorm(gg_F2[gg_F2>=0])
    bmw_F2[which(gg_F2<0)] <- 1 - pnorm(abs(gg_F2[gg_F2<0]))
    bmwsl_F2[,i] <- bmw_F2*mix_pi_3[i]
  }
  P2_F2 <- rowSums(bmwsl_F2)
  nn<-dim(as.matrix(unique(P2_F2)))[1]
  if(nn<n_samF2){P2_F2<-P2_F2+runif(n_samF2)/1e4}

  dd_F2 <- as.matrix(c(sum(P2_F2),sum(P2_F2^2),sum((P2_F2-0.5)^2)))
  WW2_F2 <- 1/(12*n_samF2) + sum((P2_F2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  u_F2 <- as.matrix(c(12*n_samF2*((dd_F2[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((dd_F2[2]/n_samF2-1/3)^2),180*n_samF2*((dd_F2[3]/n_samF2-1/12)^2)))
  D_F2 <- as.numeric(ks.test(P2_F2,"punif")[[1]][1])
  tt_F2 <- as.matrix(c((1 - pchisq(u_F2[1],1)),(1 - pchisq(u_F2[2],1)),(1 - pchisq(u_F2[3],1)),K1(WW2_F2),(1-pkolm(D_F2,n_samF2))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)
  tt_F2[which(tt_F2>=10e-4)]<-round(tt_F2[which(tt_F2>=10e-4)],4);tt_F2[which(tt_F2<10e-4)]<-format(tt_F2[which(tt_F2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX1-EAD-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigma,4),round(t(mean1D),4)," "," "," ",round(sigma1D[1],4),round(t(mix_pi_1),4)," "," "," ",
                       round(t(mean2D),4)," "," ",round(sigma2D[1],4),round(t(mix_pi_2),4)," "," ",round(t(mean3D),4)," "," "," "," "," "," "," ",round(sigma3D[1],4),round(t(mix_pi_3),4),
                       " "," "," "," "," "," "," ",round(B1[1],4)," "," "," "," "," ",round(B1[2],4)," ",round(B1[2],4)," "," "," "," "," ",round(B1[3],4),round(B1[4],4),round(jj_1,4),round(ll_1*100,4),round(mm1,4),
                       round(nn1*100,4),round(jj_2,4),round(ll_2*100,4),round(mm2,4),round(nn2*100,4),round(jj_3,4),round(ll_3*100,4),round(mm3,4),round(nn3*100,4),round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),
                       tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),
                       tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),
                       tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5],round(u_F2[1],4),tt_F2[1],round(u_F2[2],4),
                       tt_F2[2],round(u_F2[3],4),tt_F2[3],round(WW2_F2,4),tt_F2[4],round(D_F2,4),tt_F2[5])

  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2,mi_3)
  return(OUTPUT)
}
############################################MX1-NCD-AD Model #########################################
G6ModelFun[[17]] <- function(K1,logL,df11,df21,df31,df41,df51,df61){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1]; n_samF1<-dim(dataF1)[1]; n_samP2<-dim(dataP2)[1]; n_samB1<-dim(dataB1)[1]; n_samB2<-dim(dataB2)[1]; n_samF2<-dim(dataF2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2)
  sumx4<-sum(dataB1);sumx5<-sum(dataB2);sumx6<-sum(dataF2)
  sigmaP1<- as.numeric(var(dataP1)); ss1<- (n_samP1-1)*sigmaP1
  sigmaF1<- as.numeric(var(dataF1)); ss2<- (n_samF1-1)*sigmaF1
  sigmaP2<- as.numeric(var(dataP2)); ss3<- (n_samP2-1)*sigmaP2
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  mean4<- mean(dataB1); sigmaB1<- as.numeric(var(dataB1)); sigma40<- sigmaB1
  mean5<- mean(dataB2); sigmaB2<- as.numeric(var(dataB2)); sigma50<- sigmaB2
  mean6<- mean(dataF2); sigmaF2<- as.numeric(var(dataF2)); sigma60<- sigmaF2

  m_esp<-0.0001
  ####################### MX1-NCD-AD Model ################### (D4) #####################
  d21<-2; d22<-1; d23<-2
  sigma<- sigma0
  meanP1<- mean(dataP1);meanF1<- mean(dataF1);meanP2<- mean(dataP2)
  abb<- sigma60/(5.0*sigma0)

  mi_1<- matrix(0.5,d21,1)
  a1<- sqrt(sigma40/n_samB1)
  if (meanP1<meanP2) {a1= -a1}
  sigma1D<- matrix(sigmaB1/abb,d21,1)
  mean1D<- as.matrix(c(mean4+0.5*a1,0.25*(meanP1+meanP2)+0.5*meanF1))

  mi_2<- as.matrix(1)
  sigma2D<- as.matrix(sigmaB2/abb)
  mean2D<- as.matrix(mean5)

  mi_3<- as.matrix(c(0.25,0.75))
  sigma3D<- matrix(sigmaF2/abb,d23,1)
  mean3D<- as.matrix(c(mean1D[1],mean6))

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1D,sigma1D,dataB1)+logL(n_samB2,d22,mi_2,mean2D,sigma2D,dataB2)+logL(n_samF2,d23,mi_3,mean3D,sigma3D,dataF2)

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)
  WW_F2 <- matrix(0,d23,n_samF2); swx_F2 <- matrix(0,d23,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1D[i],sqrt(sigma1D[i]))/dmixnorm(dataB1,mean1D,sqrt(sigma1D),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2D[i],sqrt(sigma2D[i]))/dmixnorm(dataB2,mean2D,sqrt(sigma2D),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    for(i in 1:d23) { WW_F2[i,] <- mi_3[i]*dnorm(dataF2,mean3D[i],sqrt(sigma3D[i]))/dmixnorm(dataF2,mean3D,sqrt(sigma3D),mi_3) }
    mix_pi_3 <- as.matrix(rowSums(WW_F2)/n_samF2)
    sumwx_F2 <- WW_F2%*%dataF2

    ############ CM1-step for means ##############
    n0<- as.matrix(c(mix_pi_1*n_samB1,mix_pi_2*n_samB2,mix_pi_2*n_samB2,mix_pi_3*n_samF2))
    n0[abs(n0)<0.0001] <- 0.0001
    aaa1<- 1000; n_iter<- 0; AA<- matrix(0,4,1)
    while (aaa1>0.0001){
      n_iter<- n_iter+1
      hh<- matrix(0,4,4)
      hh[1,1]<- sigma*(1.0/n_samP1+1.0/n_samF1)+sigma1D[1]/n0[1]+sigma1D[2]/n0[2]
      hh[1,2]<- sigma/n_samF1
      hh[1,3]<- -sigma1D[1]/n0[1]+sigma1D[2]/n0[2]
      hh[1,4]<- -sigma1D[2]/n0[2]
      hh[2,2]<- sigma*(1.0/n_samF1+1.0/n_samP2)+4.0*sigma2D[1]/n_samB2
      hh[2,3]<- 0
      hh[2,4]<- -2.0*sigma2D[1]/n_samB2
      hh[3,3]<- sigma1D[1]/n0[1]+sigma1D[2]/n0[2]+sigma3D[1]/n0[5]+sigma3D[2]/n0[6]
      hh[3,4]<- -sigma1D[2]/n0[2]-2.0*sigma3D[2]/n0[6]
      hh[4,4]<- sigma1D[2]/n0[2]+sigma2D[1]/n_samB2+4.0*sigma3D[2]/n0[6]
      for(i in 2:4){
        for(j in 1:(i-1)){
          hh[i,j]<- hh[j,i]
        }
      }
      b_line<- matrix(0,4,1)
      b_line[1]<- sumx1/n_samP1+sumx2/n_samF1-sumwx_B1[1]/n0[1]-sumwx_B1[2]/n0[2]
      b_line[2]<- sumx2/n_samF1+sumx3/n_samP2-2.0*sumwx_B2[1]/n_samB2
      b_line[3]<- sumwx_B1[1]/n0[1]-sumwx_B1[2]/n0[2]-sumwx_F2[1]/n0[5]+sumwx_F2[2]/n0[6]
      b_line[4]<- sumwx_B1[2]/n0[2]+sumwx_B2[1]/n_samB2-2.0*sumwx_F2[2]/n0[6]
      B <- solve(hh,b_line)

      meanP1<- (sumx1-sigma*B[1])/n_samP1
      meanF1<- (sumx2-sigma*(B[1]+B[2]))/n_samF1
      meanP2<- (sumx3-sigma*B[2])/n_samP2
      mean1D[1]<- (sumwx_B1[1]+sigma1D[1]*(B[1]-B[3]))/n0[1]
      mean1D[2]<- (sumwx_B1[2]+sigma1D[2]*(B[1]+B[3]-B[4]))/n0[2]
      mean2D[1]<- (sumwx_B2[1]+sigma2D[1]*(2.0*B[2]-B[4]))/n0[3]
      mean3D[1]<- (sumwx_F2[1]+sigma3D[1]*B[3])/n0[5]
      mean3D[2]<- (sumwx_F2[2]+sigma3D[2]*(-B[3]+2*B[4]))/n0[6]

      aaa1<- max(abs(AA-B))
      AA<- B
      if (n_iter>20) break
    }
    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2)
    ss2<- sum((dataF1-meanF1)^2)
    ss3<- sum((dataP2-meanP2)^2)

    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1D[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2D[i])^2 }
    for(i in 1:d23) {  swx_F2[i] <- WW_F2[i,]%*%(dataF2-mean3D[i])^2 }

    sigma1D<- matrix(sum(swx_B1)/n_samB1,d21,1)
    sigma40<- sigma1D[1]-sigma
    sigma50<- sum((dataB2-mean2D[1])^2)
    sigma2D<- matrix(sigma50/n_samB2,d22,1)
    sigma50<- sigma2D[1]-sigma
    sigma3D<- matrix(sum(swx_F2)/n_samF2,d23,1)
    sigma60<- sigma3D[1]-sigma
    ################ iteratively CM3-step for variance (sigma) ################
    sigma40[(sigma40<0)]<- 0.00001; sigma50[(sigma50<0)]<- 0.00001; sigma60[(sigma60<0)]<- 0.00001
    s0<- matrix(0,2,1)
    s0[1]<- ss1+ss2+ss3; s0[2]<- n_samP1+n_samF1+n_samP2
    aaa0<- 0;aa3<- 1000
    while (aa3>0.0001){
      abc1<- sigma/(sigma+sigma40);abc2<- sigma/(sigma+sigma50);abc3<- sigma/(sigma+sigma60)
      aa4<- s0[1]+abc1*abc1*sum(swx_B1)+abc2*abc2*sum(swx_B2)+abc3*abc3*sum(swx_F2)
      aa5<- s0[2]+abc1*n_samB1+abc2*n_samB2+abc3*n_samF2
      sigma<- aa4/aa5;aa3<- abs(sigma-aaa0);aaa0<- sigma
    }
    sigma1D<- matrix(sigma+sigma40,d21,1);sigma2D<- matrix(sigma+sigma50,d22,1);sigma3D<- matrix(sigma+sigma60,d23,1)
    ####################### the stop criterion for iteration #########################################
    L1<- logL(n_samP1,1,1,meanP1,sigma,dataP1)+logL(n_samF1,1,1,meanF1,sigma,dataF1)+logL(n_samP2,1,1,meanP2,sigma,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1D,sigma1D,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2D,sigma2D,dataB2)+logL(n_samF2,d23,mix_pi_3,mean3D,sigma3D,dataF2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}

  }
  abc <- L1
  AIC <- -2*abc + 2*8

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1,1, 1,-1,-1,1,-1,-1,1,-1, 1,0,-1,0.5,0.5,-0.5,0,0, 0,1,0,0.5,0.5,0.5,0.5,0.5),8,4)
  b_line1 <- matrix(c(meanP1,meanF1,meanP2,mean1D,mean2D,mean3D))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigmaB1 - sigma1D[1]
  if(jj_1 < 0 | jj_1>=sigmaB1) {jj_1 <- 0}
  ll_1 <- jj_1/sigmaB1
  mm1<- sigma1D[1]-sigma
  if (mm1<0 | mm1>=sigmaB1) { mm1<- 0 }
  nn1<- mm1/sigmaB1

  jj_2 <- sigmaB2 - sigma2D[1]
  if(jj_2 < 0 | jj_2>=sigmaB2) {jj_2 <- 0}
  ll_2 <- jj_2/sigmaB2
  mm2<- sigma2D[1]-sigma
  if (mm2<0 | mm2>=sigmaB2) { mm2<- 0 }
  nn2<- mm2/sigmaB2

  jj_3 <- sigmaF2 - sigma3D[1]
  if(jj_3 < 0 | jj_1>=sigmaF2) {jj_3 <- 0}
  ll_3 <- jj_3/sigmaF2
  mm3<- sigma3D[1]-sigma
  if (mm3<0 | mm3>=sigmaF2) { mm3<- 0 }
  nn3<- mm3/sigmaF2

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
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

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
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),(1-pkolm(D_F1,n_samF1))))

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
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

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
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),(1-pkolm(D_B1,n_samB1))))

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
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),(1-pkolm(D_B2,n_samB2))))

  ###############  F2  ###################
  dataF2<-sort(dataF2);bmw_F2 <- matrix(0,n_samF2,1); bmwsl_F2 <- matrix(0,n_samF2,d23)
  for(i in 1:d23){
    gg_F2 <- (dataF2 - mean3D[i])/sqrt(sigma3D[i])
    bmw_F2[which(gg_F2>=0)] <- pnorm(gg_F2[gg_F2>=0])
    bmw_F2[which(gg_F2<0)] <- 1 - pnorm(abs(gg_F2[gg_F2<0]))
    bmwsl_F2[,i] <- bmw_F2*mix_pi_3[i]
  }
  P2_F2 <- rowSums(bmwsl_F2)
  nn<-dim(as.matrix(unique(P2_F2)))[1]
  if(nn<n_samF2){P2_F2<-P2_F2+runif(n_samF2)/1e4}

  dd_F2 <- as.matrix(c(sum(P2_F2),sum(P2_F2^2),sum((P2_F2-0.5)^2)))
  WW2_F2 <- 1/(12*n_samF2) + sum((P2_F2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  u_F2 <- as.matrix(c(12*n_samF2*((dd_F2[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((dd_F2[2]/n_samF2-1/3)^2),180*n_samF2*((dd_F2[3]/n_samF2-1/12)^2)))
  D_F2 <- as.numeric(ks.test(P2_F2,"punif")[[1]][1])
  tt_F2 <- as.matrix(c((1 - pchisq(u_F2[1],1)),(1 - pchisq(u_F2[2],1)),(1 - pchisq(u_F2[3],1)),K1(WW2_F2),(1-pkolm(D_F2,n_samF2))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)
  tt_F2[which(tt_F2>=10e-4)]<-round(tt_F2[which(tt_F2>=10e-4)],4);tt_F2[which(tt_F2<10e-4)]<-format(tt_F2[which(tt_F2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX1-NCD-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigma,4),round(t(mean1D),4)," "," ",round(sigma1D[1],4),round(t(mix_pi_1),4)," "," ",
                       round(t(mean2D),4)," "," "," ",round(sigma2D[1],4),round(t(mix_pi_2),4)," "," "," ",round(t(mean3D),4)," "," "," "," "," "," "," ",round(sigma3D[1],4),round(t(mix_pi_3),4),
                       " "," "," "," "," "," "," ",round(B1[1],4)," "," "," "," "," ",round(B1[2],4)," ",round(-B1[2],4)," "," "," "," "," ",round(B1[3],4),round(B1[4],4),round(jj_1,4),round(ll_1*100,4),round(mm1,4),
                       round(nn1*100,4),round(jj_2,4),round(ll_2*100,4),round(mm2,4),round(nn2*100,4),round(jj_3,4),round(ll_3*100,4),round(mm3,4),round(nn3*100,4),round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),
                       tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),
                       tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),
                       tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5],round(u_F2[1],4),tt_F2[1],round(u_F2[2],4),
                       tt_F2[2],round(u_F2[3],4),tt_F2[3],round(WW2_F2,4),tt_F2[4],round(D_F2,4),tt_F2[5])

  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2,mi_3)
  return(OUTPUT)
}
############################################ MX2-ADI-ADI Model #########################################
G6ModelFun[[18]] <- function(K1,logL,df11,df21,df31,df41,df51,df61){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1]; n_samF1<-dim(dataF1)[1]; n_samP2<-dim(dataP2)[1]; n_samB1<-dim(dataB1)[1]; n_samB2<-dim(dataB2)[1]; n_samF2<-dim(dataF2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2)
  sumx4<-sum(dataB1);sumx5<-sum(dataB2);sumx6<-sum(dataF2)
  sigmaP1<- as.numeric(var(dataP1)); ss1<- (n_samP1-1)*sigmaP1
  sigmaF1<- as.numeric(var(dataF1)); ss2<- (n_samF1-1)*sigmaF1
  sigmaP2<- as.numeric(var(dataP2)); ss3<- (n_samP2-1)*sigmaP2
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  mean4<- mean(dataB1); sigmaB1<- as.numeric(var(dataB1)); sigma40<- sigmaB1
  mean5<- mean(dataB2); sigmaB2<- as.numeric(var(dataB2)); sigma50<- sigmaB2
  mean6<- mean(dataF2); sigmaF2<- as.numeric(var(dataF2)); sigma60<- sigmaF2

  m_esp<-0.0001
  ####################### MX2-ADI-ADI Model #################### (E0) ################
  d21<-4; d22<-4; d23<-9
  sigma<- sigma0
  meanP1<- mean(dataP1);meanF1<- mean(dataF1);meanP2<- mean(dataP2)

  mi_1<- matrix(0.25,d21,1)
  a1<- sqrt(sigmaB1/n_samB1)
  if (meanP1<meanP2) {a1= -a1}
  sigma1E<- matrix(sigmaB1,d21,1)
  mean1E<- as.matrix(c(mean4+3*a1,mean4+a1,mean4-a1,mean4-3*a1))

  mi_2<- matrix(0.25,d22,1)
  a2<- sqrt(sigmaB2/n_samB2)
  if (meanP1<meanP2) { a2= -a2 }
  sigma2E<- matrix(sigmaB2,d22,1)
  mean2E<- as.matrix(c(mean5+3*a2,mean5+a2,mean5-a2,mean5-3*a2))

  mi_3<- as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  a3<- sqrt(sigmaF2/n_samF2)
  if (meanP1<meanP2) { a3= -a3 }
  sigma3E<- matrix(sigmaF2,d23,1)
  mean3E<- as.matrix(c(mean6+3*a3,mean6+1.8*a3,mean6,mean6-0.9*a3,mean6,mean6-0.9*a3,mean6,mean6-1.8*a3,mean6-3*a3))

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1E,sigma1E,dataB1)+logL(n_samB2,d22,mi_2,mean2E,sigma2E,dataB2)+logL(n_samF2,d23,mi_3,mean3E,sigma3E,dataF2)

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)
  WW_F2 <- matrix(0,d23,n_samF2); swx_F2 <- matrix(0,d23,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1E[i],sqrt(sigma1E[i]))/dmixnorm(dataB1,mean1E,sqrt(sigma1E),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2E[i],sqrt(sigma2E[i]))/dmixnorm(dataB2,mean2E,sqrt(sigma2E),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    for(i in 1:d23) { WW_F2[i,] <- mi_3[i]*dnorm(dataF2,mean3E[i],sqrt(sigma3E[i]))/dmixnorm(dataF2,mean3E,sqrt(sigma3E),mi_3) }
    mix_pi_3 <- as.matrix(rowSums(WW_F2)/n_samF2)
    sumwx_F2 <- WW_F2%*%dataF2

    ############ CM1-step for means ##############
    n0<- as.matrix(c(mix_pi_1*n_samB1,mix_pi_2*n_samB2))
    n0[abs(n0)< 0.00000001] <- 0.000001

    s0<- as.matrix(c(mix_pi_3*n_samF2))
    s0[abs(s0)< 0.00000001] <- 0.000001

    aaa0<- 0; aaa1<- 1000; n_iter<- 0; AA<- matrix(0,6,1)
    while (aaa1>0.0001){
      n_iter<- n_iter+1
      hh<- matrix(0,6,6)
      hh[1,1]<- sigma1E[1]/n0[1]+sigma1E[2]/n0[2]+sigma3E[1]/s0[1]+sigma3E[2]/s0[2]
      hh[1,2]<- 0
      hh[1,3]<- 0
      hh[1,4]<- 0
      hh[1,5]<- 0
      hh[1,6]<- sigma1E[1]/n0[1]+sigma3E[1]/s0[1]
      hh[2,2]<- sigma1E[3]/n0[3]+sigma1E[4]/n0[4]+sigma3E[4]/s0[4]+sigma3E[5]/s0[5]
      hh[2,3]<- -sigma3E[5]/s0[5]
      hh[2,4]<- 0
      hh[2,5]<- -sigma3E[5]/s0[5]
      hh[2,6]<- sigma1E[4]/n0[4]+sigma3E[5]/s0[5]
      hh[3,3]<- sigma2E[1]/n0[5]+sigma2E[2]/n0[6]+sigma3E[5]/s0[5]+sigma3E[6]/s0[6]
      hh[3,4]<- 0
      hh[3,5]<- sigma2E[1]/n0[5]+sigma3E[5]/s0[5]
      hh[3,6]<- -sigma3E[5]/s0[5]
      hh[4,4]<- sigma2E[3]/n0[7]+sigma2E[4]/n0[8]+sigma3E[8]/s0[8]+sigma3E[9]/s0[9]
      hh[4,5]<- sigma2E[4]/n0[8]+sigma3E[9]/s0[9]
      hh[4,6]<- 0
      hh[5,5]<- sigma2E[1]/n0[5]+sigma2E[4]/n0[8]+sigma3E[5]/s0[5]+sigma3E[9]/s0[9]
      hh[5,6]<- -sigma3E[5]/s0[5]
      hh[6,6]<- sigma1E[1]/n0[1]+sigma1E[4]/n0[4]+sigma3E[1]/s0[1]+sigma3E[5]/s0[5]
      for(i in 2:6){
        for(j in 1:(i-1)){
          hh[i,j]<- hh[j,i]
        }
      }
      b_line<- matrix(0,6,1)
      b_line[1]<- sumwx_B1[1]/n0[1]-sumwx_B1[2]/n0[2]-sumwx_F2[1]/s0[1]+sumwx_F2[2]/s0[2]
      b_line[2]<- sumwx_B1[3]/n0[3]-sumwx_B1[4]/n0[4]-sumwx_F2[4]/s0[4]+sumwx_F2[5]/s0[5]
      b_line[3]<- sumwx_B2[1]/n0[5]-sumwx_B2[2]/n0[6]-sumwx_F2[5]/s0[5]+sumwx_F2[6]/s0[6]
      b_line[4]<- sumwx_B2[3]/n0[7]-sumwx_B2[4]/n0[8]-sumwx_F2[8]/s0[8]+sumwx_F2[9]/s0[9]
      b_line[5]<- sumwx_B2[1]/n0[5]-sumwx_B2[4]/n0[8]-sumwx_F2[5]/s0[5]+sumwx_F2[9]/s0[9]
      b_line[6]<- sumwx_B1[1]/n0[1]-sumwx_B1[4]/n0[4]-sumwx_F2[1]/s0[1]+sumwx_F2[5]/s0[5]
      B <- solve(hh,b_line)

      mean1E[1]<- (sumwx_B1[1]-sigma1E[1]*(B[1]+B[6]))/n0[1]
      mean1E[2]<- (sumwx_B1[2]+sigma1E[2]*B[1])/n0[2]
      mean1E[3]<- (sumwx_B1[3]-sigma1E[3]*B[2])/n0[3]
      mean1E[4]<- (sumwx_B1[4]+sigma1E[4]*(B[2]+B[6]))/n0[4]

      mean2E[1]<- (sumwx_B2[1]-sigma2E[1]*(B[3]+B[5]))/n0[5]
      mean2E[2]<- (sumwx_B2[2]+sigma2E[2]*B[3])/n0[6]
      mean2E[3]<- (sumwx_B2[3]-sigma2E[3]*B[4])/n0[7]
      mean2E[4]<- (sumwx_B2[4]+sigma2E[4]*(B[4]+B[5]))/n0[8]

      mean3E[1]<- (sumwx_F2[1]+sigma3E[1]*(B[1]+B[6]))/s0[1]
      mean3E[2]<- (sumwx_F2[2]-sigma3E[2]*B[1])/s0[2]
      mean3E[3]<- sumwx_F2[3]/s0[3]
      mean3E[4]<- (sumwx_F2[4]+sigma3E[4]*B[2])/s0[4]
      mean3E[5]<- (sumwx_F2[5]+sigma3E[5]*(-B[2]+B[3]+B[5]-B[6]))/s0[5]
      mean3E[6]<- (sumwx_F2[6]-sigma3E[6]*B[3])/s0[6]
      mean3E[7]<- sumwx_F2[7]/s0[7]
      mean3E[8]<- (sumwx_F2[8]+sigma3E[8]*B[4])/s0[8]
      mean3E[9]<- (sumwx_F2[9]-sigma3E[9]*(B[4]+B[5]))/s0[9]

      aaa1<- max(abs(AA-B))
      AA<- B
      if (n_iter>20) break
    }
    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2)
    ss2<- sum((dataF1-meanF1)^2)
    ss3<- sum((dataP2-meanP2)^2)

    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1E[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2E[i])^2 }
    for(i in 1:d23) {  swx_F2[i] <- WW_F2[i,]%*%(dataF2-mean3E[i])^2 }
    aa1<-sum(swx_B1);aa2<-sum(swx_B2);aa3<-sum(swx_F2)
    sigma1E<- matrix(aa1/n_samB1,d21,1)
    sigma40<- sigma1E[1]-sigma
    sigma2E<- matrix(aa2/n_samB2,d22,1)
    sigma50<- sigma2E[1]-sigma
    sigma3E<- matrix(aa3/n_samF2,d23,1)
    sigma60<- sigma3E[1]-sigma
    ################ iteratively CM3-step for variance (sigma) ################
    sigma40[(sigma40<0)]<- 0.00001; sigma50[(sigma50<0)]<- 0.00001; sigma60[(sigma60<0)]<- 0.00001
    ab1<- ss1+ss2+ss3; ab2<- n_samP1+n_samF1+n_samP2
    aaa0<- sigma;n_iter<- 0; aaa1<- 1000
    while ( aaa1 >0.0001){
      n_iter<- n_iter+1
      aa4<- sigma/(sigma+sigma40); aa4[(aa4>1.0)]<- 1
      aa5<- sigma/(sigma+sigma50); aa5[(aa5>1.0)]<- 1
      aa6<- sigma/(sigma+sigma60); aa6[(aa6>1.0)]<- 1
      sigma<- (ab1+aa4*aa4*aa1+aa5*aa5*aa2+aa6*aa6*aa3)/(ab2+aa4*n_samB1+aa5*n_samB2+aa6*n_samF2)
      aaa1<- abs(sigma-aaa0)
      aaa0<- sigma
      if (n_iter>20) break
    }
    sigma1E<- matrix(sigma+sigma40,d21,1);sigma2E<- matrix(sigma+sigma50,d22,1);sigma3E<- matrix(sigma+sigma60,d23,1)
    ####################### the stop criterion for iteration ############################################################
    L1<- logL(n_samP1,1,1,meanP1,sigma,dataP1)+logL(n_samF1,1,1,meanF1,sigma,dataF1)+logL(n_samP2,1,1,meanP2,sigma,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1E,sigma1E,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2E,sigma2E,dataB2)+logL(n_samF2,d23,mix_pi_3,mean3E,sigma3E,dataF2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*18

  #########first order genetic parameter process##########
  aa<- matrix(c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1, 1,0,-1,1,1,0,0,0,0,-1,-1,1,1,1,0,0,0,-1,-1,-1, 1,0,-1,1,0,1,0,0,-1,0,-1,1,0,-1,1,0,-1,1,0,-1,
                0,1,0,0,0,1,1,1,1,0,0,0,0,0,1,1,1,0,0,0, 0,1,0,0,1,0,1,1,0,1,0,0,1,0,0,1,0,0,1,0, 1,0,1,1,0,0,0,0,0,0,1,1,0,-1,0,0,0,-1,0,1, 0,0,0,0,1,0,0,0,0,-1,0,0,1,0,0,0,0,0,-1,0,
                0,0,0,0,0,1,0,0,-1,0,0,0,0,0,1,0,-1,0,0,0, 0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0 ),20,14)
  b_line1 <- matrix(c(meanP1,meanF1,meanP2,mean1E,mean2E,mean3E))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigmaB1 - sigma1E[1]
  if(jj_1 < 0 | jj_1>=sigmaB1) {jj_1 <- 0}
  ll_1 <- jj_1/sigmaB1
  mm1<- sigma1E[1]-sigma
  if (mm1<0 | mm1>=sigmaB1) { mm1<- 0 }
  nn1<- mm1/sigmaB1

  jj_2 <- sigmaB2 - sigma2E[1]
  if(jj_2 < 0 | jj_2>=sigmaB2) {jj_2 <- 0}
  ll_2 <- jj_2/sigmaB2
  mm2<- sigma2E[1]-sigma
  if (mm2<0 | mm2>=sigmaB2) { mm2<- 0 }
  nn2<- mm2/sigmaB2

  jj_3 <- sigmaF2 - sigma3E[1]
  if(jj_3 < 0 | jj_3>=sigmaF2) {jj_3 <- 0}
  ll_3 <- jj_3/sigmaF2
  mm3<- sigma3E[1]-sigma
  if (mm3<0 | mm3>=sigmaF2) { mm3<- 0 }
  nn3<- mm3/sigmaF2

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
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

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
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),(1-pkolm(D_F1,n_samF1))))

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
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

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
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),(1-pkolm(D_B1,n_samB1))))

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
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),(1-pkolm(D_B2,n_samB2))))

  ###############  F2  ###################
  dataF2<-sort(dataF2);bmw_F2 <- matrix(0,n_samF2,1); bmwsl_F2 <- matrix(0,n_samF2,d23)
  for(i in 1:d23){
    gg_F2 <- (dataF2 - mean3E[i])/sqrt(sigma3E[i])
    bmw_F2[which(gg_F2>=0)] <- pnorm(gg_F2[gg_F2>=0])
    bmw_F2[which(gg_F2<0)] <- 1 - pnorm(abs(gg_F2[gg_F2<0]))
    bmwsl_F2[,i] <- bmw_F2*mix_pi_3[i]
  }
  P2_F2 <- rowSums(bmwsl_F2)
  nn<-dim(as.matrix(unique(P2_F2)))[1]
  if(nn<n_samF2){P2_F2<-P2_F2+runif(n_samF2)/1e4}

  dd_F2 <- as.matrix(c(sum(P2_F2),sum(P2_F2^2),sum((P2_F2-0.5)^2)))
  WW2_F2 <- 1/(12*n_samF2) + sum((P2_F2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  u_F2 <- as.matrix(c(12*n_samF2*((dd_F2[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((dd_F2[2]/n_samF2-1/3)^2),180*n_samF2*((dd_F2[3]/n_samF2-1/12)^2)))
  D_F2 <- as.numeric(ks.test(P2_F2,"punif")[[1]][1])
  tt_F2 <- as.matrix(c((1 - pchisq(u_F2[1],1)),(1 - pchisq(u_F2[2],1)),(1 - pchisq(u_F2[3],1)),K1(WW2_F2),(1-pkolm(D_F2,n_samF2))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)
  tt_F2[which(tt_F2>=10e-4)]<-round(tt_F2[which(tt_F2>=10e-4)],4);tt_F2[which(tt_F2<10e-4)]<-format(tt_F2[which(tt_F2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame(" MX2-ADI-ADI",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigma,4),round(t(mean1E),4),round(sigma1E[1],4),round(t(mix_pi_1),4),
                       round(t(mean2E),4),round(sigma2E[1],4),round(t(mix_pi_2),4),round(t(mean3E),4),round(sigma3E[1],4),round(t(mix_pi_3),4),
                       round(B1[1],4),round(B1[2],4),round(B1[3],4),round(B1[4],4),round(B1[5],4),round(B1[6],4),round(B1[7],4),round(B1[8],4),round(B1[9],4),round(B1[10],4),round(B1[11],4),
                       round(B1[12],4),round(B1[13],4),round(B1[14],4)," "," ",round(jj_1,4),round(ll_1*100,4),round(mm1,4),
                       round(nn1*100,4),round(jj_2,4),round(ll_2*100,4),round(mm2,4),round(nn2*100,4),round(jj_3,4),round(ll_3*100,4),round(mm3,4),round(nn3*100,4),round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),
                       tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),
                       tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),
                       tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5],round(u_F2[1],4),tt_F2[1],round(u_F2[2],4),
                       tt_F2[2],round(u_F2[3],4),tt_F2[3],round(WW2_F2,4),tt_F2[4],round(D_F2,4),tt_F2[5])

  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2,mi_3)
  return(OUTPUT)
}
############################################ MX2-ADI-AD Model #########################################
G6ModelFun[[19]] <- function(K1,logL,df11,df21,df31,df41,df51,df61){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1]; n_samF1<-dim(dataF1)[1]; n_samP2<-dim(dataP2)[1]; n_samB1<-dim(dataB1)[1]; n_samB2<-dim(dataB2)[1]; n_samF2<-dim(dataF2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2)
  sumx4<-sum(dataB1);sumx5<-sum(dataB2);sumx6<-sum(dataF2)
  sigmaP1<- as.numeric(var(dataP1)); ss1<- (n_samP1-1)*sigmaP1
  sigmaF1<- as.numeric(var(dataF1)); ss2<- (n_samF1-1)*sigmaF1
  sigmaP2<- as.numeric(var(dataP2)); ss3<- (n_samP2-1)*sigmaP2
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  mean4<- mean(dataB1); sigmaB1<- as.numeric(var(dataB1)); sigma40<- sigmaB1
  mean5<- mean(dataB2); sigmaB2<- as.numeric(var(dataB2)); sigma50<- sigmaB2
  mean6<- mean(dataF2); sigmaF2<- as.numeric(var(dataF2)); sigma60<- sigmaF2

  m_esp<-0.0001
  ####################### MX2-ADI-AD Model #################### (E1) ################
  d21<-4; d22<-4; d23<-9
  sigma<- sigma0
  meanP1<- mean(dataP1);meanF1<- mean(dataF1);meanP2<- mean(dataP2)

  mi_1<- matrix(0.25,d21,1)
  a1<- sqrt(sigmaB1/n_samB1)
  if (meanP1<meanP2) {a1= -a1}
  sigma1E<- matrix(sigmaB1,d21,1)
  mean1E<- as.matrix(c(mean4+2.5*a1,mean4+1.5*a1,mean4-1.5*a1,mean4-2.5*a1))

  mi_2<- matrix(0.25,d22,1)
  a2<- sqrt(sigmaB2/n_samB2)
  if (meanP1<meanP2) { a2= -a2 }
  sigma2E<- matrix(sigmaB2,d22,1)
  mean2E<- as.matrix(c(mean5+2.5*a2,mean5+1.5*a2,mean5-1.5*a2,mean5-2.5*a2))

  mi_3<- as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  a3<- sqrt(sigmaF2/n_samF2)
  if (meanP1<meanP2) { a3= -a3 }
  sigma3E<- matrix(sigmaF2,d23,1)
  mean3E<- as.matrix(c(mean6+3*a3,mean6+1.8*a3,mean6,mean6+1.2*a3,mean6+0.5*a3,mean6-0.5*a3,mean6,mean6-1.8*a3,mean6-3*a3))

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1E,sigma1E,dataB1)+logL(n_samB2,d22,mi_2,mean2E,sigma2E,dataB2)+logL(n_samF2,d23,mi_3,mean3E,sigma3E,dataF2)

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)
  WW_F2 <- matrix(0,d23,n_samF2); swx_F2 <- matrix(0,d23,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1E[i],sqrt(sigma1E[i]))/dmixnorm(dataB1,mean1E,sqrt(sigma1E),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2E[i],sqrt(sigma2E[i]))/dmixnorm(dataB2,mean2E,sqrt(sigma2E),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    for(i in 1:d23) { WW_F2[i,] <- mi_3[i]*dnorm(dataF2,mean3E[i],sqrt(sigma3E[i]))/dmixnorm(dataF2,mean3E,sqrt(sigma3E),mi_3) }
    mix_pi_3 <- as.matrix(rowSums(WW_F2)/n_samF2)
    sumwx_F2 <- WW_F2%*%dataF2

    ############ CM1-step for means ##############
    n0<- as.matrix(c(mix_pi_1*n_samB1,mix_pi_2*n_samB2))
    n0[abs(n0)<0.00000001] <- 0.000001

    s0<- as.matrix(c(mix_pi_3*n_samF2))
    s0[abs(s0)<0.00000001] <- 0.000001

    aaa0<- 0; aaa1<- 1000; n_iter<- 0;  AA<- matrix(0,9,1)
    while (aaa1>0.0001){
      n_iter<- n_iter+1
      hh<- matrix(0,9,9)
      hh[1,1]<- sigma/n_samP1+sigma/n_samP2+sigma1E[1]/n0[1]+sigma1E[4]/n0[4]+sigma2E[1]/n0[5]+sigma2E[4]/n0[8]
      hh[1,2]<- sigma/n_samP1-sigma/n_samP2+sigma1E[1]/n0[1]-sigma2E[4]/n0[8]
      hh[1,3]<- -sigma1E[1]/n0[1]
      hh[1,4]<- sigma1E[4]/n0[4]
      hh[1,5]<- sigma2E[1]/n0[5]
      hh[1,6]<- -sigma2E[4]/n0[8]
      hh[1,7]<- sigma2E[1]/n0[5]-sigma2E[4]/n0[8]
      hh[1,8]<- -sigma1E[1]/n0[1]+sigma1E[4]/n0[4]
      hh[1,9]<- hh[9,1]<- 0
      hh[2,2]<- sigma/n_samP1+4.0*sigma/n_samF1+sigma/n_samP2+sigma1E[1]/n0[1]+sigma2E[4]/n0[8]+4.0*sigma3E[5]/s0[5]
      hh[2,3]<- -sigma1E[1]/n0[1]
      hh[2,4]<- -2.0*sigma3E[5]/s0[5]
      hh[2,5]<- 2.0*sigma3E[5]/s0[5]
      hh[2,6]<- sigma2E[4]/n0[8]
      hh[2,7]<- sigma2E[4]/n0[8]+2.0*sigma3E[5]/s0[5]
      hh[2,8]<- -sigma1E[1]/n0[1]-2.0*sigma3E[5]/s0[5]
      hh[2,9]<- 0
      hh[3,3]<- sigma1E[1]/n0[1]+sigma1E[2]/n0[2]+sigma3E[1]/s0[1]+sigma3E[2]/s0[2]
      hh[3,4]<- 0
      hh[3,5]<- 0
      hh[3,6]<- 0
      hh[3,7]<- 0
      hh[3,8]<- sigma1E[1]/n0[1]+sigma3E[1]/s0[1]
      hh[3,9]<- -sigma1E[2]/n0[2]-sigma3E[2]/s0[2]
      hh[4,4]<- sigma1E[3]/n0[3]+sigma1E[4]/n0[4]+sigma3E[4]/s0[4]+sigma3E[5]/s0[5]
      hh[4,5]<- -sigma3E[5]/s0[5]
      hh[4,6]<- 0
      hh[4,7]<- -sigma3E[5]/s0[5]
      hh[4,8]<- sigma1E[4]/n0[4]+sigma3E[5]/s0[5]
      hh[4,9]<- 0
      hh[5,5]<- sigma2E[1]/n0[5]+sigma2E[2]/n0[6]+sigma3E[5]/s0[5]+sigma3E[6]/s0[6]
      hh[5,6]<- 0
      hh[5,7]<- sigma2E[1]/n0[5]+sigma3E[5]/s0[5]
      hh[5,8]<- -sigma3E[5]/s0[5]
      hh[5,9]<- 0
      hh[6,6]<- sigma2E[3]/n0[7]+sigma2E[4]/n0[8]+sigma3E[8]/s0[8]+sigma3E[9]/s0[9]
      hh[6,7]<- sigma2E[4]/n0[8]+sigma3E[9]/s0[9]
      hh[6,8]<- 0
      hh[6,9]<- sigma2E[3]/n0[7]+sigma3E[8]/s0[8]
      hh[7,7]<- sigma2E[1]/n0[5]+sigma2E[4]/n0[8]+sigma3E[5]/s0[5]+sigma3E[9]/s0[9]
      hh[7,8]<- -sigma3E[5]/s0[5]
      hh[7,9]<- 0
      hh[8,8]<- sigma1E[1]/n0[1]+sigma1E[4]/n0[4]+sigma3E[1]/s0[1]+sigma3E[5]/s0[5]
      hh[8,9]<- 0
      hh[9,9]<- sigma1E[2]/n0[2]+sigma2E[3]/n0[7]+sigma3E[2]/s0[2]+sigma3E[8]/s0[8]
      for(i in 2:9){
        for(j in 1:(i-1)){
          hh[i,j]<- hh[j,i]
        }
      }
      b_line<- matrix(0,9,1)
      b_line[1]<- sumx1/n_samP1-sumx3/n_samP2-sumwx_B1[1]/n0[1]-sumwx_B1[4]/n0[4]+sumwx_B2[1]/n0[5]+sumwx_B2[4]/n0[8]
      b_line[2]<- sumx1/n_samP1+2.0*sumx2/n_samF1+sumx3/n_samP2-sumwx_B1[1]/n0[1]-sumwx_B2[4]/n0[8]-2.0*sumwx_F2[5]/s0[5]
      b_line[3]<- sumwx_B1[1]/n0[1]-sumwx_B1[2]/n0[2]-sumwx_F2[1]/s0[1]+sumwx_F2[2]/s0[2]
      b_line[4]<- sumwx_B1[3]/n0[3]-sumwx_B1[4]/n0[4]-sumwx_F2[4]/s0[4]+sumwx_F2[5]/s0[5]
      b_line[5]<- sumwx_B2[1]/n0[5]-sumwx_B2[2]/n0[6]-sumwx_F2[5]/s0[5]+sumwx_F2[6]/s0[6]
      b_line[6]<- sumwx_B2[3]/n0[7]-sumwx_B2[4]/n0[8]-sumwx_F2[8]/s0[8]+sumwx_F2[9]/s0[9]
      b_line[7]<- sumwx_B2[1]/n0[5]-sumwx_B2[4]/n0[8]-sumwx_F2[5]/s0[5]+sumwx_F2[9]/s0[9]
      b_line[8]<- sumwx_B1[1]/n0[1]-sumwx_B1[4]/n0[4]-sumwx_F2[1]/s0[1]+sumwx_F2[5]/s0[5]
      b_line[9]<- sumwx_B1[2]/n0[2]+sumwx_B2[3]/n0[7]-sumwx_F2[2]/s0[2]-sumwx_F2[8]/s0[8]
      B <- solve(hh,b_line)

      meanP1<- (sumx1-sigma*(B[1]+B[2]))/n_samP1
      meanF1<- (sumx2-sigma*2.0*B[2])/n_samF1
      meanP2<- (sumx3+sigma*(B[1]-B[2]))/n_samP2

      mean1E[1]<- (sumwx_B1[1]+sigma1E[1]*(B[1]+B[2]-B[3]-B[8]))/n0[1]
      mean1E[2]<- (sumwx_B1[2]+sigma1E[2]*(B[3]-B[9]))/n0[2]
      mean1E[3]<- (sumwx_B1[3]-sigma1E[3]*B[4])/n0[3]
      mean1E[4]<- (sumwx_B1[4]+sigma1E[4]*(B[1]+B[4]+B[8]))/n0[4]

      mean2E[1]<- (sumwx_B2[1]-sigma2E[1]*(B[1]+B[5]+B[7]))/n0[5]
      mean2E[2]<- (sumwx_B2[2]+sigma2E[2]*B[5])/n0[6]
      mean2E[3]<- (sumwx_B2[3]-sigma2E[3]*(B[6]+B[9]))/n0[7]
      mean2E[4]<- (sumwx_B2[4]+sigma2E[4]*(-B[1]+B[2]+B[6]+B[7]))/n0[8]

      mean3E[1]<- (sumwx_F2[1]+sigma3E[1]*(B[3]+B[8]))/s0[1]
      mean3E[2]<- (sumwx_F2[2]+sigma3E[2]*(-B[3]+B[9]))/s0[2]
      mean3E[3]<- sumwx_F2[3]/s0[3]
      mean3E[4]<- (sumwx_F2[4]+sigma3E[4]*B[4])/s0[4]
      mean3E[5]<- (sumwx_F2[5]+sigma3E[5]*(2.0*B[2]-B[4]+B[5]+B[7]-B[8]))/s0[5]
      mean3E[6]<- (sumwx_F2[6]-sigma3E[6]*B[5])/s0[6]
      mean3E[7]<- sumwx_F2[7]/s0[7]
      mean3E[8]<- (sumwx_F2[8]+sigma3E[8]*(B[6]+B[9]))/s0[8]
      mean3E[9]<- (sumwx_F2[9]-sigma3E[9]*(B[6]+B[7]))/s0[9]

      aaa1<- max(abs(AA-B))
      AA<- B
      if (n_iter>20) break
    }
    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2)
    ss2<- sum((dataF1-meanF1)^2)
    ss3<- sum((dataP2-meanP2)^2)

    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1E[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2E[i])^2 }
    for(i in 1:d23) {  swx_F2[i] <- WW_F2[i,]%*%(dataF2-mean3E[i])^2 }

    aa1<-sum(swx_B1); aa2<-sum(swx_B2); aa3<-sum(swx_F2)
    sigma1E<- matrix(aa1/n_samB1,d21,1)
    sigma40<- sigma1E[1]-sigma
    sigma2E<- matrix(aa2/n_samB2,d22,1)
    sigma50<- sigma2E[1]-sigma
    sigma3E<- matrix(aa3/n_samF2,d23,1)
    sigma60<- sigma3E[1]-sigma
    ################ iteratively CM3-step for variance (sigma) ################
    sigma40[(sigma40<0)]<- 0.00001; sigma50[(sigma50<0)]<- 0.00001; sigma60[(sigma60<0)]<- 0.00001
    ab1<- ss1+ss2+ss3; ab2<- n_samP1+n_samF1+n_samP2
    aaa0<- sigma ;n_iter<- 0;aaa1<- 1000
    while (aaa1>0.0001){
      n_iter<- n_iter+1
      aa4<- sigma/(sigma+sigma40); aa4[(aa4>1.0)]<- 1
      aa5<- sigma/(sigma+sigma50); aa5[(aa5>1.0)]<- 1
      aa6<- sigma/(sigma+sigma60); aa6[(aa6>1.0)]<- 1
      sigma<- (ab1+aa4*aa4*aa1+aa5*aa5*aa2+aa6*aa6*aa3)/(ab2+aa4*n_samB1+aa5*n_samB2+aa6*n_samF2)
      aaa1<- abs(sigma-aaa0);aaa0<- sigma
      if (n_iter>20) break
    }
    sigma1E<- matrix(sigma+sigma40,d21,1);sigma2E<- matrix(sigma+sigma50,d22,1);sigma3E<- matrix(sigma+sigma60,d23,1)
    ####################### the stop criterion for iteration ############################################################
    L1<- logL(n_samP1,1,1,meanP1,sigma,dataP1)+logL(n_samF1,1,1,meanF1,sigma,dataF1)+logL(n_samP2,1,1,meanP2,sigma,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1E,sigma1E,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2E,sigma2E,dataB2)+logL(n_samF2,d23,mix_pi_3,mean3E,sigma3E,dataF2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*15

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, 1,0,-1,1,1,0,0,0,0,-1,-1,1,1,1,0,0,0,-1,-1,-1, 1,0,-1,1,0,1,0,0,-1,0,-1,1,0,-1,1,0,-1,1,0,-1, 0,1,0,0,0,1,1,1,1,0,0,0,0,0,1,1,1,0,0,0,
                0,1,0,0,1,0,1,1,0,1,0,0,1,0,0,1,0,0,1,0, 1,0,1,1,0,0,0,0,0,0,1,1,0,-1,0,0,0,-1,0,1, 0,0,0,0,1,0,0,0,0,-1,0,0,1,0,0,0,0,0,-1,0, 0,0,0,0,0,1,0,0,-1,0,0,0,0,0,1,0,-1,0,0,0,
                0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0, 1,0,-1,0.5,0.5,0.5,0.5,-0.5,-0.5,-0.5,-0.5,0,0,0,0,0,0,0,0,0, 0,1,0,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5 ),20,11)
  b_line1 <- matrix(c(meanP1,meanF1,meanP2,mean1E,mean2E,mean3E))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigmaB1 - sigma1E[1]
  if(jj_1 < 0 | jj_1>=sigmaB1) {jj_1 <- 0}
  ll_1 <- jj_1/sigmaB1
  mm1<- sigma1E[1]-sigma
  if (mm1<0 | mm1>=sigmaB1) { mm1<- 0 }
  nn1<- mm1/sigmaB1

  jj_2 <- sigmaB2 - sigma2E[1]
  if(jj_2 < 0 | jj_2>=sigmaB2) {jj_2 <- 0}
  ll_2 <- jj_2/sigmaB2
  mm2<- sigma2E[1]-sigma
  if (mm2<0 | mm2>=sigmaB2) { mm2<- 0 }
  nn2<- mm2/sigmaB2

  jj_3 <- sigmaF2 - sigma3E[1]
  if(jj_3 < 0 | jj_3>=sigmaF2) {jj_3 <- 0}
  ll_3 <- jj_3/sigmaF2
  mm3<- sigma3E[1]-sigma
  if (mm3<0 | mm3>=sigmaF2) { mm3<- 0 }
  nn3<- mm3/sigmaF2

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
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

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
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),(1-pkolm(D_F1,n_samF1))))

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
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

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
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),(1-pkolm(D_B1,n_samB1))))

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
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),(1-pkolm(D_B2,n_samB2))))

  ###############  F2  ###################
  dataF2<-sort(dataF2);bmw_F2 <- matrix(0,n_samF2,1); bmwsl_F2 <- matrix(0,n_samF2,d23)
  for(i in 1:d23){
    gg_F2 <- (dataF2 - mean3E[i])/sqrt(sigma3E[i])
    bmw_F2[which(gg_F2>=0)] <- pnorm(gg_F2[gg_F2>=0])
    bmw_F2[which(gg_F2<0)] <- 1 - pnorm(abs(gg_F2[gg_F2<0]))
    bmwsl_F2[,i] <- bmw_F2*mix_pi_3[i]
  }
  P2_F2 <- rowSums(bmwsl_F2)
  nn<-dim(as.matrix(unique(P2_F2)))[1]
  if(nn<n_samF2){P2_F2<-P2_F2+runif(n_samF2)/1e4}

  dd_F2 <- as.matrix(c(sum(P2_F2),sum(P2_F2^2),sum((P2_F2-0.5)^2)))
  WW2_F2 <- 1/(12*n_samF2) + sum((P2_F2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  u_F2 <- as.matrix(c(12*n_samF2*((dd_F2[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((dd_F2[2]/n_samF2-1/3)^2),180*n_samF2*((dd_F2[3]/n_samF2-1/12)^2)))
  D_F2 <- as.numeric(ks.test(P2_F2,"punif")[[1]][1])
  tt_F2 <- as.matrix(c((1 - pchisq(u_F2[1],1)),(1 - pchisq(u_F2[2],1)),(1 - pchisq(u_F2[3],1)),K1(WW2_F2),(1-pkolm(D_F2,n_samF2))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)
  tt_F2[which(tt_F2>=10e-4)]<-round(tt_F2[which(tt_F2>=10e-4)],4);tt_F2[which(tt_F2<10e-4)]<-format(tt_F2[which(tt_F2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame(" MX2-ADI-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigma,4),round(t(mean1E),4),round(sigma1E[1],4),round(t(mix_pi_1),4),
                       round(t(mean2E),4),round(sigma2E[1],4),round(t(mix_pi_2),4),round(t(mean3E),4),round(sigma3E[1],4),round(t(mix_pi_3),4),
                       round(B1[1],4)," "," "," "," "," ",round(B1[2],4),round(B1[3],4),round(B1[4],4),round(B1[5],4),round(B1[6],4),round(B1[7],4),round(B1[8],4),
                       round(B1[9],4),round(B1[10],4),round(B1[11],4),round(jj_1,4),round(ll_1*100,4),round(mm1,4),
                       round(nn1*100,4),round(jj_2,4),round(ll_2*100,4),round(mm2,4),round(nn2*100,4),round(jj_3,4),round(ll_3*100,4),round(mm3,4),round(nn3*100,4),round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),
                       tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),
                       tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),
                       tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5],round(u_F2[1],4),tt_F2[1],round(u_F2[2],4),
                       tt_F2[2],round(u_F2[3],4),tt_F2[3],round(WW2_F2,4),tt_F2[4],round(D_F2,4),tt_F2[5])

  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2,mi_3)
  return(OUTPUT)
}
############################################ MX2-AD-AD Model #########################################
G6ModelFun[[20]] <- function(K1,logL,df11,df21,df31,df41,df51,df61){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1]; n_samF1<-dim(dataF1)[1]; n_samP2<-dim(dataP2)[1]; n_samB1<-dim(dataB1)[1]; n_samB2<-dim(dataB2)[1]; n_samF2<-dim(dataF2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2)
  sumx4<-sum(dataB1);sumx5<-sum(dataB2);sumx6<-sum(dataF2)
  sigmaP1<- as.numeric(var(dataP1)); ss1<- (n_samP1-1)*sigmaP1
  sigmaF1<- as.numeric(var(dataF1)); ss2<- (n_samF1-1)*sigmaF1
  sigmaP2<- as.numeric(var(dataP2)); ss3<- (n_samP2-1)*sigmaP2
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  mean4<- mean(dataB1); sigmaB1<- as.numeric(var(dataB1)); sigma40<- sigmaB1
  mean5<- mean(dataB2); sigmaB2<- as.numeric(var(dataB2)); sigma50<- sigmaB2
  mean6<- mean(dataF2); sigmaF2<- as.numeric(var(dataF2)); sigma60<- sigmaF2

  m_esp<-0.0001
  ####################### MX2-AD-AD Model #################### (E2) ################
  d21<-4; d22<-4; d23<-9
  sigma<- sigma0
  meanP1<- mean(dataP1);meanF1<- mean(dataF1);meanP2<- mean(dataP2)
  abb<- sigmaF2/(10*sigma0)

  mi_1<- matrix(0.25,d21,1)
  a1<- sqrt(sigmaB1/n_samB1)
  if (meanP1<meanP2) {a1= -a1}
  sigma1E<- matrix((sigmaB1/abb),d21,1)
  mean1E<- as.matrix(c(mean4+3*a1,mean4+a1,mean4-a1,mean4-3*a1))

  mi_2<- matrix(0.25,d22,1)
  a2<- sqrt(sigmaB2/n_samB2)
  if (meanP1<meanP2) { a2= -a2 }
  sigma2E<- matrix((sigmaB2/abb),d22,1)
  mean2E<- as.matrix(c(mean5+3*a2,mean5+a2,mean5-a2,mean5-3*a2))

  mi_3<- as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  a3<- sqrt(sigmaF2/n_samF2)
  if (meanP1<meanP2) { a3= -a3 }
  sigma3E<- matrix((sigmaF2/abb),d23,1)
  mean3E<- as.matrix(c(mean6+3*a3,mean6+2.3*a3,mean6+1.8*a3,mean6+0.8*a3,mean6,mean6-0.8*a3,mean6-1.8*a3,mean6-2.3*a3,mean6-3*a3))

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1E,sigma1E,dataB1)+logL(n_samB2,d22,mi_2,mean2E,sigma2E,dataB2)+logL(n_samF2,d23,mi_3,mean3E,sigma3E,dataF2)

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)
  WW_F2 <- matrix(0,d23,n_samF2); swx_F2 <- matrix(0,d23,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1E[i],sqrt(sigma1E[i]))/dmixnorm(dataB1,mean1E,sqrt(sigma1E),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2E[i],sqrt(sigma2E[i]))/dmixnorm(dataB2,mean2E,sqrt(sigma2E),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    for(i in 1:d23) { WW_F2[i,] <- mi_3[i]*dnorm(dataF2,mean3E[i],sqrt(sigma3E[i]))/dmixnorm(dataF2,mean3E,sqrt(sigma3E),mi_3) }
    mix_pi_3 <- as.matrix(rowSums(WW_F2)/n_samF2)
    sumwx_F2 <- WW_F2%*%dataF2

    ############ CM1-step for means ##############
    n0<- as.matrix(c(mix_pi_1*n_samB1,mix_pi_2*n_samB2))
    n0[abs(n0)<0.00000001] <- 0.000001

    s0<- as.matrix(c(mix_pi_3*n_samF2))
    s0[abs(s0)<0.00000001] <- 0.000001

    aaa0<- 0;aaa1<- 1000; n_iter<- 0; AA<- matrix(0,13,1)
    while (aaa1>0.0001){
      n_iter<- n_iter+1
      hh<- matrix(0,13,13)
      hh[1,1]<- sigma/n_samP1+sigma/n_samP2+sigma1E[1]/n0[1]+sigma1E[4]/n0[4]+sigma2E[1]/n0[5]+sigma2E[4]/n0[8]
      hh[1,2]<- sigma/n_samP1-sigma/n_samP2+sigma1E[1]/n0[1]-sigma2E[4]/n0[8]
      hh[1,3]<- sigma/n_samP1+sigma/n_samP2
      hh[1,4]<- sigma/n_samP1+sigma/n_samP2
      hh[1,5]<- -sigma1E[1]/n0[1]
      hh[1,6]<- sigma1E[4]/n0[4]
      hh[1,7]<- sigma2E[1]/n0[5]
      hh[1,8]<- -sigma2E[4]/n0[8]
      hh[1,9]<- sigma2E[1]/n0[5]-sigma2E[4]/n0[8]
      hh[1,10]<- -sigma1E[1]/n0[1]+sigma1E[4]/n0[4]
      hh[1,11]<- hh[1,12]<- hh[1,13]<- 0
      hh[2,2]<- sigma*(1.0/n_samP1+4.0/n_samF1+1.0/n_samP2)+sigma1E[1]/n0[1]+sigma2E[4]/n0[8]+4.0*sigma3E[5]/s0[5]
      hh[2,3]<- sigma*(1.0/n_samP1-1.0/n_samP2)
      hh[2,4]<- sigma*(1.0/n_samP1-1.0/n_samP2)
      hh[2,5]<- -sigma1E[1]/n0[1]
      hh[2,6]<- -2.0*sigma3E[5]/s0[5]
      hh[2,7]<- 2.0*sigma3E[5]/s0[5]
      hh[2,8]<- sigma2E[4]/n0[8]
      hh[2,9]<- sigma2E[4]/n0[8]+2.0*sigma3E[5]/s0[5]
      hh[2,10]<- -sigma1E[1]/n0[1]-2.0*sigma3E[5]/s0[5]
      hh[2,11]<- hh[2,12]<- 0
      hh[2,13]<- 4.0*sigma3E[5]/s0[5]
      hh[3,3]<- sigma/n_samP1+sigma/n_samP2+4.0*sigma1E[3]/n0[3]+4.0*sigma2E[2]/n0[6]+sigma3E[3]/s0[3]+sigma3E[7]/s0[7]
      hh[3,4]<- sigma/n_samP1+sigma/n_samP2-sigma3E[3]/s0[3]-sigma3E[7]/s0[7]
      hh[3,5]<- 0
      hh[3,6]<- -2.0*sigma1E[3]/n0[3]
      hh[3,7]<- -2.0*sigma2E[2]/n0[6]
      hh[3,8]<- hh[3,9]<- hh[3,10]<- hh[3,11]<- 0
      hh[3,12]<- sigma3E[3]/s0[3]-sigma3E[7]/s0[7]
      hh[3,13]<- -sigma3E[7]/s0[7]
      hh[4,4]<- sigma/n_samP1+sigma/n_samP2+4.0*sigma1E[2]/n0[2]+4.0*sigma2E[3]/n0[7]+sigma3E[3]/s0[3]+sigma3E[7]/s0[7]
      hh[4,5]<- 2.0*sigma1E[2]/n0[2]
      hh[4,6]<- hh[4,7]<- 0
      hh[4,8]<- 2.0*sigma2E[3]/n0[7]
      hh[4,9]<- hh[4,10]<- 0
      hh[4,11]<- -2.0*sigma1E[2]/n0[2]+2.0*sigma2E[3]/n0[7]
      hh[4,12]<- -sigma3E[3]/s0[3]+sigma3E[7]/s0[7]
      hh[4,13]<- sigma3E[7]/s0[7]
      hh[5,5]<- sigma1E[1]/n0[1]+sigma1E[2]/n0[2]+sigma3E[1]/s0[1]+sigma3E[2]/s0[2]
      hh[5,6]<- hh[5,7]<- hh[5,8]<- hh[5,9]<- 0
      hh[5,10]<- sigma1E[1]/n0[1]+sigma3E[1]/s0[1]
      hh[5,11]<- -sigma1E[2]/n0[2]-sigma3E[2]/s0[2]
      hh[5,12]<- -sigma3E[1]/s0[1]
      hh[5,13]<- 0
      hh[6,6]<- sigma1E[3]/n0[3]+sigma1E[4]/n0[4]+sigma3E[4]/s0[4]+sigma3E[5]/s0[5]
      hh[6,7]<- -sigma3E[5]/s0[5]
      hh[6,8]<- 0
      hh[6,9]<- -sigma3E[5]/s0[5]
      hh[6,10]<- sigma1E[4]/n0[4]+sigma3E[5]/s0[5]
      hh[6,11]<- hh[6,12]<- 0
      hh[6,13]<- -sigma3E[4]/s0[4]-2.0*sigma3E[5]/s0[5]
      hh[7,7]<- sigma2E[1]/n0[5]+sigma2E[2]/n0[6]+sigma3E[5]/s0[5]+sigma3E[6]/s0[6]
      hh[7,8]<- 0
      hh[7,9]<- sigma2E[1]/n0[5]+sigma3E[5]/s0[5]
      hh[7,10]<- -sigma3E[5]/s0[5]
      hh[7,11]<- hh[7,12]<- 0
      hh[7,13]<- 2.0*sigma3E[5]/s0[5]+sigma3E[6]/s0[6]
      hh[8,8]<- sigma2E[3]/n0[7]+sigma2E[4]/n0[8]+sigma3E[8]/s0[8]+sigma3E[9]/s0[9]
      hh[8,9]<- sigma2E[4]/n0[8]+sigma3E[9]/s0[9]
      hh[8,10]<- 0
      hh[8,11]<- sigma2E[3]/n0[7]+sigma3E[8]/s0[8]
      hh[8,12]<- sigma3E[9]/s0[9]
      hh[8,13]<- -2.0*sigma3E[8]/s0[8]-sigma3E[9]/s0[9]
      hh[9,9]<- sigma2E[1]/n0[5]+sigma2E[4]/n0[8]+sigma3E[5]/s0[5]+sigma3E[9]/s0[9]
      hh[9,10]<- -sigma3E[5]/s0[5]
      hh[9,11]<- 0
      hh[9,12]<- sigma3E[9]/s0[9]
      hh[9,13]<- 2.0*sigma3E[5]/s0[5]-sigma3E[9]/s0[9]
      hh[10,10]<- sigma1E[1]/n0[1]+sigma1E[4]/n0[4]+sigma3E[1]/s0[1]+sigma3E[5]/s0[5]
      hh[10,11]<- 0
      hh[10,12]<- -sigma3E[1]/s0[1]
      hh[10,13]<- -2.0*sigma3E[5]/s0[5]
      hh[11,11]<- sigma1E[2]/n0[2]+sigma2E[3]/n0[7]+sigma3E[2]/s0[2]+sigma3E[8]/s0[8]
      hh[11,12]<- 0
      hh[11,13]<- -2.0*sigma3E[8]/s0[8]
      hh[12,12]<- sigma3E[1]/s0[1]+sigma3E[3]/s0[3]+sigma3E[7]/s0[7]+sigma3E[9]/s0[9]
      hh[12,13]<- sigma3E[7]/s0[7]-sigma3E[9]/s0[9]
      hh[13,13]<- sigma3E[4]/s0[4]+4.0*sigma3E[5]/s0[5]+sigma3E[6]/s0[6]+sigma3E[7]/s0[7]+4.0*sigma3E[8]/s0[8]+sigma3E[9]/s0[9]
      for(i in 2:13){
        for(j in 1:(i-1)){
          hh[i,j]<- hh[j,i]
        }
      }
      b_line<- matrix(0,13,1)
      b_line[1]<- sumx1/n_samP1-sumx3/n_samP2-sumwx_B1[1]/n0[1]-sumwx_B1[4]/n0[4]+sumwx_B2[1]/n0[5]+sumwx_B2[4]/n0[8]
      b_line[2]<- sumx1/n_samP1+2.0*sumx2/n_samF1+sumx3/n_samP2-sumwx_B1[1]/n0[1]-sumwx_B2[4]/n0[8]-2.0*sumwx_F2[5]/s0[5]
      b_line[3]<- sumx1/n_samP1-sumx3/n_samP2-2.0*sumwx_B1[3]/n0[3]+2.0*sumwx_B2[2]/n0[6]-sumwx_F2[3]/s0[3]+sumwx_F2[7]/s0[7]
      b_line[4]<- sumx1/n_samP1-sumx3/n_samP2-2.0*sumwx_B1[2]/n0[2]+2.0*sumwx_B2[3]/n0[7]+sumwx_F2[3]/s0[3]-sumwx_F2[7]/s0[7]
      b_line[5]<- sumwx_B1[1]/n0[1]-sumwx_B1[2]/n0[2]-sumwx_F2[1]/s0[1]+sumwx_F2[2]/s0[2]

      b_line[6]<- sumwx_B1[3]/n0[3]-sumwx_B1[4]/n0[4]-sumwx_F2[4]/s0[4]+sumwx_F2[5]/s0[5]
      b_line[7]<- sumwx_B2[1]/n0[5]-sumwx_B2[2]/n0[6]-sumwx_F2[5]/s0[5]+sumwx_F2[6]/s0[6]
      b_line[8]<- sumwx_B2[3]/n0[7]-sumwx_B2[4]/n0[8]-sumwx_F2[8]/s0[8]+sumwx_F2[9]/s0[9]
      b_line[9]<- sumwx_B2[1]/n0[5]-sumwx_B2[4]/n0[8]-sumwx_F2[5]/s0[5]+sumwx_F2[9]/s0[9]

      b_line[10]<- sumwx_B1[1]/n0[1]-sumwx_B1[4]/n0[4]-sumwx_F2[1]/s0[1]+sumwx_F2[5]/s0[5]
      b_line[11]<- sumwx_B1[2]/n0[2]+sumwx_B2[3]/n0[7]-sumwx_F2[2]/s0[2]-sumwx_F2[8]/s0[8]
      b_line[12]<- sumwx_F2[1]/s0[1]-sumwx_F2[3]/s0[3]-sumwx_F2[7]/s0[7]+sumwx_F2[9]/s0[9]
      b_line[13]<- sumwx_F2[4]/s0[4]-2.0*sumwx_F2[5]/s0[5]+sumwx_F2[6]/s0[6]-sumwx_F2[7]/s0[7]+2.0*sumwx_F2[8]/s0[8]-sumwx_F2[9]/s0[9]
      B <- solve(hh,b_line)

      meanP1<- (sumx1-sigma*(B[1]+B[2]+B[3]+B[4]))/n_samP1
      meanF1<- (sumx2-2.0*sigma*B[2])/n_samF1
      meanP2<- (sumx3+sigma*(B[1]-B[2]+B[3]+B[4]))/n_samP2
      mean1E[1]<- (sumwx_B1[1]+sigma1E[1]*(B[1]+B[2]-B[5]-B[10]))/n0[1]
      mean1E[2]<- (sumwx_B1[2]+sigma1E[2]*(2.0*B[4]+B[5]-B[11]))/n0[2]
      mean1E[3]<- (sumwx_B1[3]+sigma1E[3]*(2.0*B[3]-B[6]))/n0[3]
      mean1E[4]<- (sumwx_B1[4]+sigma1E[4]*(B[1]+B[6]+B[10]))/n0[4]

      mean2E[1]<- (sumwx_B2[1]-sigma2E[1]*(B[1]+B[7]+B[9]))/n0[5]
      mean2E[2]<- (sumwx_B2[2]-sigma2E[2]*(2.0*B[3]-B[7]))/n0[6]
      mean2E[3]<- (sumwx_B2[3]-sigma2E[3]*(2.0*B[4]+B[8]+B[11]))/n0[7]
      mean2E[4]<- (sumwx_B2[4]-sigma2E[4]*(B[1]-B[2]-B[8]-B[9]))/n0[8]

      mean3E[1]<- (sumwx_F2[1]+sigma3E[1]*(B[5]+B[10]-B[12]))/s0[1]
      mean3E[2]<- (sumwx_F2[2]-sigma3E[2]*(B[5]-B[11]))/s0[2]
      mean3E[3]<- (sumwx_F2[3]+sigma3E[3]*(B[3]-B[4]+B[12]))/s0[3]
      mean3E[4]<- (sumwx_F2[4]+sigma3E[4]*(B[6]-B[13]))/s0[4]
      mean3E[5]<- (sumwx_F2[5]+sigma3E[5]*(2.0*B[2]-B[6]+B[7]+B[9]-B[10]+2.0*B[13]))/s0[5]
      mean3E[6]<- (sumwx_F2[6]-sigma3E[6]*(B[7]+B[13]))/s0[6]
      mean3E[7]<- (sumwx_F2[7]+sigma3E[7]*(-B[3]+B[4]+B[12]+B[13]))/s0[7]
      mean3E[8]<- (sumwx_F2[8]+sigma3E[8]*(B[8]+B[11]-2.0*B[13]))/s0[8]
      mean3E[9]<- (sumwx_F2[9]-sigma3E[9]*(B[8]+B[9]+B[12]-B[13]))/s0[9]

      aaa1<- max(abs(AA-B))
      AA<- B
      if (n_iter>20) break
    }
    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2)
    ss2<- sum((dataF1-meanF1)^2)
    ss3<- sum((dataP2-meanP2)^2)

    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1E[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2E[i])^2 }
    for(i in 1:d23) {  swx_F2[i] <- WW_F2[i,]%*%(dataF2-mean3E[i])^2 }

    aa1<-sum(swx_B1);aa2<-sum(swx_B2);aa3<-sum(swx_F2)
    sigma1E<- matrix(aa1/n_samB1,d21,1)
    sigma40<- sigma1E[1]-sigma
    sigma2E<- matrix(aa2/n_samB2,d22,1)
    sigma50<- sigma2E[1]-sigma
    sigma3E<- matrix(aa3/n_samF2,d23,1)
    sigma60<- sigma3E[1]-sigma
    ################ iteratively CM3-step for variance (sigma) ################
    sigma40[(sigma40<0)]<- 0.00001; sigma50[(sigma50<0)]<- 0.00001; sigma60[(sigma60<0)]<- 0.00001
    ab1<- ss1+ss2+ss3; ab2<- n_samP1+n_samF1+n_samP2
    aaa0<- sigma ;n_iter<- 0;aaa1<- 1000
    while (aaa1>0.0001){
      n_iter<- n_iter+1
      aa4<- sigma/(sigma+sigma40); aa4[(aa4>1.0)]<- 1
      aa5<- sigma/(sigma+sigma50); aa5[(aa5>1.0)]<- 1
      aa6<- sigma/(sigma+sigma60); aa6[(aa6>1.0)]<- 1
      sigma<- (ab1+aa4*aa4*aa1+aa5*aa5*aa2+aa6*aa6*aa3)/(ab2+aa4*n_samB1+aa5*n_samB2+aa6*n_samF2)
      aaa1<- abs(sigma-aaa0)
      aaa0<- sigma
      if (n_iter>20) break
    }
    sigma1E<- matrix(sigma+sigma40,d21,1);sigma2E<- matrix(sigma+sigma50,d22,1);sigma3E<- matrix(sigma+sigma60,d23,1)
    ####################### the stop criterion for iteration ############################################################
    L1<- logL(n_samP1,1,1,meanP1,sigma,dataP1)+logL(n_samF1,1,1,meanF1,sigma,dataF1)+logL(n_samP2,1,1,meanP2,sigma,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1E,sigma1E,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2E,sigma2E,dataB2)+logL(n_samF2,d23,mix_pi_3,mean3E,sigma3E,dataF2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*11

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, 1,0,-1,1,1,0,0,0,0,-1,-1,1,1,1,0,0,0,-1,-1,-1, 1,0,-1,1,0,1,0,0,-1,0,-1,1,0,-1,1,0,-1,1,0,-1, 0,1,0,0,0,1,1,1,1,0,0,0,0,0,1,1,1,0,0,0,
                0,1,0,0,1,0,1,1,0,1,0,0,1,0,0,1,0,0,1,0, 1,0,-1,0.5,0.5,0.5,0.5,-0.5,-0.5,-0.5,-0.5,0,0,0,0,0,0,0,0,0, 0,1,0,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5 ),20,7)

  b_line1 <- matrix(c(meanP1,meanF1,meanP2,mean1E,mean2E,mean3E))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigmaB1 - sigma1E[1]
  if(jj_1 < 0 | jj_1>=sigmaB1) {jj_1 <- 0}
  ll_1 <- jj_1/sigmaB1
  mm1<- sigma1E[1]-sigma
  if (mm1<0 | mm1>=sigmaB1) { mm1<- 0 }
  nn1<- mm1/sigmaB1

  jj_2 <- sigmaB2 - sigma2E[1]
  if(jj_2 < 0 | jj_2>=sigmaB2) {jj_2 <- 0}
  ll_2 <- jj_2/sigmaB2
  mm2<- sigma2E[1]-sigma
  if (mm2<0 | mm2>=sigmaB2) { mm2<- 0 }
  nn2<- mm2/sigmaB2

  jj_3 <- sigmaF2 - sigma3E[1]
  if(jj_3 < 0 | jj_3>=sigmaF2) {jj_3 <- 0}
  ll_3 <- jj_3/sigmaF2
  mm3<- sigma3E[1]-sigma
  if (mm3<0 | mm3>=sigmaF2) { mm3<- 0 }
  nn3<- mm3/sigmaF2

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
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

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
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),(1-pkolm(D_F1,n_samF1))))

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
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

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
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),(1-pkolm(D_B1,n_samB1))))

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
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),(1-pkolm(D_B2,n_samB2))))

  ###############  F2  ###################
  dataF2<-sort(dataF2);bmw_F2 <- matrix(0,n_samF2,1); bmwsl_F2 <- matrix(0,n_samF2,d23)
  for(i in 1:d23){
    gg_F2 <- (dataF2 - mean3E[i])/sqrt(sigma3E[i])
    bmw_F2[which(gg_F2>=0)] <- pnorm(gg_F2[gg_F2>=0])
    bmw_F2[which(gg_F2<0)] <- 1 - pnorm(abs(gg_F2[gg_F2<0]))
    bmwsl_F2[,i] <- bmw_F2*mix_pi_3[i]
  }
  P2_F2 <- rowSums(bmwsl_F2)
  nn<-dim(as.matrix(unique(P2_F2)))[1]
  if(nn<n_samF2){P2_F2<-P2_F2+runif(n_samF2)/1e4}

  dd_F2 <- as.matrix(c(sum(P2_F2),sum(P2_F2^2),sum((P2_F2-0.5)^2)))
  WW2_F2 <- 1/(12*n_samF2) + sum((P2_F2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  u_F2 <- as.matrix(c(12*n_samF2*((dd_F2[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((dd_F2[2]/n_samF2-1/3)^2),180*n_samF2*((dd_F2[3]/n_samF2-1/12)^2)))
  D_F2 <- as.numeric(ks.test(P2_F2,"punif")[[1]][1])
  tt_F2 <- as.matrix(c((1 - pchisq(u_F2[1],1)),(1 - pchisq(u_F2[2],1)),(1 - pchisq(u_F2[3],1)),K1(WW2_F2),(1-pkolm(D_F2,n_samF2))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)
  tt_F2[which(tt_F2>=10e-4)]<-round(tt_F2[which(tt_F2>=10e-4)],4);tt_F2[which(tt_F2<10e-4)]<-format(tt_F2[which(tt_F2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame(" MX2-AD-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigma,4),round(t(mean1E),4),round(sigma1E[1],4),round(t(mix_pi_1),4),
                       round(t(mean2E),4),round(sigma2E[1],4),round(t(mix_pi_2),4),round(t(mean3E),4),round(sigma3E[1],4),round(t(mix_pi_3),4),
                       round(B1[1],4)," "," "," "," "," ",round(B1[2],4),round(B1[3],4),round(B1[4],4),round(B1[5],4)," "," "," "," ",
                       round(B1[6],4),round(B1[7],4),round(jj_1,4),round(ll_1*100,4),round(mm1,4),round(nn1*100,4),
                       round(jj_2,4),round(ll_2*100,4),round(mm2,4),round(nn2*100,4),round(jj_3,4),round(ll_3*100,4),round(mm3,4),round(nn3*100,4),round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),
                       tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),
                       tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),
                       tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5],round(u_F2[1],4),tt_F2[1],round(u_F2[2],4),
                       tt_F2[2],round(u_F2[3],4),tt_F2[3],round(WW2_F2,4),tt_F2[4],round(D_F2,4),tt_F2[5])

  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2,mi_3)
  return(OUTPUT)
}
############################################ MX2-A-AD Model #########################################
G6ModelFun[[21]] <- function(K1,logL,df11,df21,df31,df41,df51,df61){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1]; n_samF1<-dim(dataF1)[1]; n_samP2<-dim(dataP2)[1]; n_samB1<-dim(dataB1)[1]; n_samB2<-dim(dataB2)[1]; n_samF2<-dim(dataF2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2)
  sumx4<-sum(dataB1);sumx5<-sum(dataB2);sumx6<-sum(dataF2)
  sigmaP1<- as.numeric(var(dataP1)); ss1<- (n_samP1-1)*sigmaP1
  sigmaF1<- as.numeric(var(dataF1)); ss2<- (n_samF1-1)*sigmaF1
  sigmaP2<- as.numeric(var(dataP2)); ss3<- (n_samP2-1)*sigmaP2
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  mean4<- mean(dataB1); sigmaB1<- as.numeric(var(dataB1)); sigma40<- sigmaB1
  mean5<- mean(dataB2); sigmaB2<- as.numeric(var(dataB2)); sigma50<- sigmaB2
  mean6<- mean(dataF2); sigmaF2<- as.numeric(var(dataF2)); sigma60<- sigmaF2

  m_esp<-0.0001
  ####################### MX2-A-AD Model #################### (E3) ################
  d21<-4; d22<-4; d23<-9
  sigma<- sigma0
  meanP1<- mean(dataP1);meanF1<- mean(dataF1);meanP2<- mean(dataP2)
  abb<- sigmaF2/(10*sigma0)

  mi_1<- matrix(0.25,d21,1)
  a1<- sqrt(sigmaB1)
  if (meanP1<meanP2) {a1= -a1}
  sigma1E<- matrix(sigmaB1/abb,d21,1)
  mean1E<- as.matrix(c(mean4+3*a1,mean4+a1,mean4-a1,mean4-3*a1))

  mi_2<- matrix(0.25,d22,1)
  a2<- sqrt(sigmaB2)
  if (meanP1<meanP2) { a2= -a2 }
  sigma2E<- matrix(sigmaB2/abb,d22,1)
  mean2E<- as.matrix(c(mean5+3*a2,mean5+a2,mean5-a2,mean5+3*a2))

  mi_3<- as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  a3<- sqrt(sigmaF2)
  if (meanP1<meanP2) { a3= -a3 }
  sigma3E<- matrix(sigmaF2/abb,d23,1)
  mean3E<- as.matrix(c(mean6+3.2*a3,mean6+2.4*a3,mean6+1.6*a3,mean6+0.8*a3,mean6,mean6-0.8*a3,mean6-1.6*a3,mean6-2.4*a3,mean6-3*a3))

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1E,sigma1E,dataB1)+logL(n_samB2,d22,mi_2,mean2E,sigma2E,dataB2)+logL(n_samF2,d23,mi_3,mean3E,sigma3E,dataF2)

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)
  WW_F2 <- matrix(0,d23,n_samF2); swx_F2 <- matrix(0,d23,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1E[i],sqrt(sigma1E[i]))/dmixnorm(dataB1,mean1E,sqrt(sigma1E),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2E[i],sqrt(sigma2E[i]))/dmixnorm(dataB2,mean2E,sqrt(sigma2E),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    for(i in 1:d23) { WW_F2[i,] <- mi_3[i]*dnorm(dataF2,mean3E[i],sqrt(sigma3E[i]))/dmixnorm(dataF2,mean3E,sqrt(sigma3E),mi_3) }
    mix_pi_3 <- as.matrix(rowSums(WW_F2)/n_samF2)
    sumwx_F2 <- WW_F2%*%dataF2

    ############ CM1-step for means ##############
    n0<- as.matrix(c(mix_pi_1*n_samB1,mix_pi_2*n_samB2))
    n0[abs(n0)<0.00000001] <- 0.000001

    s0<- as.matrix(c(mix_pi_3*n_samF2))
    s0[abs(s0)<0.00000001] <- 0.000001

    aaa0<- 0; aaa1<- 1000; n_iter<- 0; AA<- matrix(0,15,1)
    while (aaa1>0.001){
      n_iter<- n_iter+1
      hh<- matrix(0,15,15)
      hh[1,1]<- sigma*(1.0/n_samP1+4.0/n_samF1+1.0/n_samP2)+16.0*sigma3E[5]/s0[5]
      hh[1,2]<- sigma*(1.0/n_samP1-1.0/n_samP2)
      hh[1,3]<- sigma*(1.0/n_samP1-1.0/n_samP2)
      hh[1,4]<- sigma*(1.0/n_samP1-1.0/n_samP2)
      hh[1,5]<- 0
      hh[1,6]<- -4.0*sigma3E[5]/s0[5]
      hh[1,7]<- -4.0*sigma3E[5]/s0[5]
      hh[1,8]<- 4.0*sigma3E[5]/s0[5]
      hh[1,9]<- 0
      hh[1,10]<- 4.0*sigma3E[5]/s0[5]
      hh[1,11]<- 8.0*sigma3E[5]/s0[5]
      hh[1,12]<- 8.0*sigma3E[5]/s0[5]
      hh[1,13]<- 8.0*sigma3E[5]/s0[5]
      hh[1,14]<- 8.0*sigma3E[5]/s0[5]
      hh[1,15]<- 0

      hh[2,2]<- sigma*(1.0/n_samP1+1.0/n_samP2)+sigma1E[1]/n0[1]+sigma1E[4]/n0[4]+sigma2E[1]/n0[5]+sigma2E[4]/n0[8]
      hh[2,3]<- sigma*(1.0/n_samP1+1.0/n_samP2)
      hh[2,4]<- sigma*(1.0/n_samP1+1.0/n_samP2)
      hh[2,5]<- -sigma1E[1]/n0[1]
      hh[2,6]<- sigma1E[4]/n0[4]
      hh[2,7]<- -sigma1E[1]/n0[1]+sigma1E[4]/n0[4]
      hh[2,8]<- sigma2E[1]/n0[5]
      hh[2,9]<- -sigma2E[4]/n0[8]
      hh[2,10]<- sigma2E[1]/n0[5]-sigma2E[4]/n0[8]
      hh[2,11]<- hh[2,12]<- hh[2,13]<- hh[2,14]<- hh[2,15]<- 0

      hh[3,3]<- sigma*(1.0/n_samP1+1.0/n_samP2)+4.0*sigma1E[3]/n0[3]+4.0*sigma2E[2]/n0[6]+sigma3E[3]/s0[3]+sigma3E[7]/s0[7]
      hh[3,4]<- sigma/n_samP1+sigma/n_samP2-sigma3E[3]/s0[3]-sigma3E[7]/s0[7]
      hh[3,5]<- 0
      hh[3,6]<- -2.0*sigma1E[3]/n0[3]
      hh[3,7]<- 0
      hh[3,8]<- -2.0*sigma2E[2]/n0[6]
      hh[3,9]<- hh[3,10]<- hh[3,11]<- hh[3,12]<- 0
      hh[3,13]<- -sigma3E[3]/s0[3]+sigma3E[4]/s0[4]
      hh[3,14]<- 0
      hh[3,15]<- sigma3E[7]/s0[7]

      hh[4,4]<- sigma*(1.0/n_samP1+1.0/n_samP2)+4.0*sigma1E[2]/n0[2]+4.0*sigma2E[3]/n0[7]+sigma3E[3]/s0[3]+sigma3E[7]/s0[7]
      hh[4,5]<- 2.0*sigma1E[2]/n0[2]
      hh[4,6]<- hh[4,7]<- hh[4,8]<- 0
      hh[4,9]<- 2.0*sigma2E[3]/n0[7]
      hh[4,10]<- hh[4,11]<- hh[4,12]<- 0
      hh[4,13]<- sigma3E[3]/s0[3]-sigma3E[7]/s0[7]
      hh[4,14]<- -2.0*sigma1E[2]/n0[2]+2.0*sigma2E[3]/n0[7]
      hh[4,15]<- -sigma3E[7]/s0[7]

      hh[5,5]<- sigma1E[1]/n0[1]+sigma1E[2]/n0[2]+sigma3E[1]/s0[1]+sigma3E[2]/s0[2]
      hh[5,6]<- 0
      hh[5,7]<- sigma1E[1]/n0[1]+sigma3E[1]/s0[1]
      hh[5,8]<- hh[5,9]<- hh[5,10]<- 0
      hh[5,11]<- -sigma3E[1]/s0[1]
      hh[5,12]<- sigma3E[2]/s0[2]
      hh[5,13]<- 0
      hh[5,14]<- -sigma1E[2]/n0[2]
      hh[5,15]<- 0

      hh[6,6]<- sigma1E[3]/n0[3]+sigma1E[4]/n0[4]+sigma3E[4]/s0[4]+sigma3E[5]/s0[5]
      hh[6,7]<- sigma1E[4]/n0[4]+sigma3E[5]/s0[5]
      hh[6,8]<- -sigma3E[5]/s0[5]
      hh[6,9]<- 0
      hh[6,10]<- -sigma3E[5]/s0[5]
      hh[6,11]<- -2.0*sigma3E[5]/s0[5]
      hh[6,12]<- -2.0*sigma3E[5]/s0[5]
      hh[6,13]<- -2.0*sigma3E[5]/s0[5]
      hh[6,14]<- -2.0*sigma3E[5]/s0[5]
      hh[6,15]<- 0

      hh[7,7]<- sigma1E[1]/n0[1]+sigma1E[4]/n0[4]+sigma3E[1]/s0[1]+sigma3E[5]/s0[5]
      hh[7,8]<- -sigma3E[5]/s0[5]
      hh[7,9]<- 0
      hh[7,10]<- -sigma3E[5]/s0[5]
      hh[7,11]<- -sigma3E[1]/s0[1]-2.0*sigma3E[5]/s0[5]
      hh[7,12]<- -2.0*sigma3E[5]/s0[5]
      hh[7,13]<- -2.0*sigma3E[5]/s0[5]
      hh[7,14]<- -2.0*sigma3E[5]/s0[5]
      hh[7,15]<- 0

      hh[8,8]<- sigma2E[1]/n0[5]+sigma2E[2]/n0[6]+sigma3E[5]/s0[5]+sigma3E[6]/s0[6]
      hh[8,9]<- 0
      hh[8,10]<- sigma2E[1]/n0[5]+sigma3E[5]/s0[5]
      hh[8,11]<- 2.0*sigma3E[5]/s0[5]
      hh[8,12]<- 2.0*sigma3E[5]/s0[5]
      hh[8,13]<- 2.0*sigma3E[5]/s0[5]
      hh[8,14]<- 2.0*sigma3E[5]/s0[5]
      hh[8,15]<- 0

      hh[9,9]<- sigma2E[3]/n0[7]+sigma2E[4]/n0[8]+sigma3E[8]/s0[8]+sigma3E[9]/s0[9]
      hh[9,10]<- sigma2E[4]/n0[8]+sigma3E[9]/s0[9]
      hh[9,11]<- sigma3E[9]/s0[9]
      hh[9,12]<- -sigma3E[8]/s0[8]
      hh[9,13]<- 0
      hh[9,14]<- sigma2E[3]/n0[7]
      hh[9,15]<- 2.0*sigma3E[8]/s0[8]+sigma3E[9]/s0[9]

      hh[10,10]<- sigma2E[1]/n0[5]+sigma2E[4]/n0[8]+sigma3E[5]/s0[5]+sigma3E[9]/s0[9]
      hh[10,11]<- 2.0*sigma3E[5]/s0[5]+sigma3E[9]/s0[9]
      hh[10,12]<- 2.0*sigma3E[5]/s0[5]
      hh[10,13]<- 2.0*sigma3E[5]/s0[5]
      hh[10,14]<- 2.0*sigma3E[5]/s0[5]
      hh[10,15]<- sigma3E[9]/s0[9]

      hh[11,11]<- sigma3E[1]/s0[1]+4.0*sigma3E[5]/s0[5]+sigma3E[9]/s0[9]
      hh[11,12]<- 4.0*sigma3E[5]/s0[5]
      hh[11,13]<- 4.0*sigma3E[5]/s0[5]
      hh[11,14]<- 4.0*sigma3E[5]/s0[5]
      hh[11,15]<- sigma3E[9]/s0[9]

      hh[12,12]<- sigma3E[2]/s0[2]+4.0*sigma3E[5]/s0[5]+sigma3E[8]/s0[8]
      hh[12,13]<- 4.0*sigma3E[5]/s0[5]
      hh[12,14]<- 4.0*sigma3E[5]/s0[5]
      hh[12,15]<- -2.0*sigma3E[8]/s0[8]

      hh[13,13]<- sigma3E[3]/s0[3]+4.0*sigma3E[5]/s0[5]+sigma3E[7]/s0[7]
      hh[13,14]<- 4.0*sigma3E[5]/s0[5]
      hh[13,15]<- sigma3E[7]/s0[7]

      hh[14,14]<- sigma1E[2]/n0[2]+sigma2E[3]/n0[7]+4.0*sigma3E[5]/s0[5]
      hh[14,15]<- 0

      hh[15,15]<- sigma3E[7]/s0[7]+4.0*sigma3E[8]/s0[8]+sigma3E[9]/s0[9]
      for(i in 2:15){
        for(j in 1:(i-1)){
          hh[i,j]<- hh[j,i]
        }
      }
      b_line<- matrix(0,15,1)
      b_line[1]<- sumx1/n_samP1+2.0*sumx2/n_samF1+sumx3/n_samP2-4.0*sumwx_F2[5]/s0[5]
      b_line[2]<- sumx1/n_samP1-sumx3/n_samP2-sumwx_B1[1]/n0[1]-sumwx_B1[4]/n0[4]+sumwx_B2[1]/n0[5]+sumwx_B2[4]/n0[8]
      b_line[3]<- sumx1/n_samP1-sumx3/n_samP2-2.0*sumwx_B1[3]/n0[3]+2.0*sumwx_B2[2]/n0[6]-sumwx_F2[3]/s0[3]+sumwx_F2[7]/s0[7]

      b_line[4]<- sumx1/n_samP1-sumx3/n_samP2-2.0*sumwx_B1[2]/n0[2]+2.0*sumwx_B2[3]/n0[7]+sumwx_F2[3]/s0[3]-sumwx_F2[7]/s0[7]
      b_line[5]<- sumwx_B1[1]/n0[1]-sumwx_B1[2]/n0[2]-sumwx_F2[1]/s0[1]+sumwx_F2[2]/s0[2]
      b_line[6]<- sumwx_B1[3]/n0[3]-sumwx_B1[4]/n0[4]-sumwx_F2[4]/s0[4]+sumwx_F2[5]/s0[5]

      b_line[7]<- sumwx_B1[1]/n0[1]-sumwx_B1[4]/n0[4]-sumwx_F2[1]/s0[1]+sumwx_F2[5]/s0[5]
      b_line[8]<- sumwx_B2[1]/n0[5]-sumwx_B2[2]/n0[6]-sumwx_F2[5]/s0[5]+sumwx_F2[6]/s0[6]
      b_line[9]<- sumwx_B2[3]/n0[7]-sumwx_B2[4]/n0[8]-sumwx_F2[8]/s0[8]+sumwx_F2[9]/s0[9]

      b_line[10]<- sumwx_B2[1]/n0[5]-sumwx_B2[4]/n0[8]-sumwx_F2[5]/s0[5]+sumwx_F2[9]/s0[9]
      b_line[11]<- sumwx_F2[1]/s0[1]-2.0*sumwx_F2[5]/s0[5]+sumwx_F2[9]/s0[9]
      b_line[12]<- sumwx_F2[2]/s0[2]-2.0*sumwx_F2[5]/s0[5]+sumwx_F2[8]/s0[8]

      b_line[13]<- sumwx_F2[3]/s0[3]-2.0*sumwx_F2[5]/s0[5]+sumwx_F2[7]/s0[7]
      b_line[14]<- sumwx_B1[2]/n0[2]+sumwx_B2[3]/n0[7]-2.0*sumwx_F2[5]/s0[5]
      b_line[15]<- sumwx_F2[7]/s0[7]-2.0*sumwx_F2[8]/s0[8]+sumwx_F2[9]/s0[9]

      B <- solve(hh,b_line)

      meanP1<- (sumx1-sigma*(B[1]+B[2]+B[3]+B[4]))/n_samP1
      meanF1<- (sumx2-sigma*2.0*B[1])/n_samF1
      meanP2<- (sumx3-sigma*(B[1]-B[2]-B[3]-B[4]))/n_samP2

      mean1E[1]<- (sumwx_B1[1]+sigma1E[1]*(B[2]-B[5]-B[7]))/n0[1]
      mean1E[2]<- (sumwx_B1[2]+sigma1E[2]*(2.0*B[4]+B[5]-B[14]))/n0[2]
      mean1E[3]<- (sumwx_B1[3]+sigma1E[3]*(2.0*B[3]-B[6]))/n0[3]
      mean1E[4]<- (sumwx_B1[4]+sigma1E[4]*(B[2]+B[6]+B[7]))/n0[4]

      mean2E[1]<- (sumwx_B2[1]-sigma2E[1]*(B[2]+B[8]+B[10]))/n0[5]
      mean2E[2]<- (sumwx_B2[2]-sigma2E[2]*(2.0*B[3]-B[8]))/n0[6]
      mean2E[3]<- (sumwx_B2[3]-sigma2E[3]*(2.0*B[4]+B[9]+B[14]))/n0[7]
      mean2E[4]<- (sumwx_B2[4]-sigma2E[4]*(B[2]-B[9]-B[10]))/n0[8]

      mean3E[1]<- (sumwx_F2[1]+sigma3E[1]*(B[5]+B[7]-B[11]))/s0[1]
      mean3E[2]<- (sumwx_F2[2]-sigma3E[2]*(B[5]+B[12]))/s0[2]
      mean3E[3]<- (sumwx_F2[3]+sigma3E[3]*(B[3]-B[4]-B[13]))/s0[3]
      mean3E[4]<- (sumwx_F2[4]+sigma3E[4]*B[6])/s0[4]
      mean3E[5]<- (sumwx_F2[5]+sigma3E[5]*(4.0*B[1]-B[6]-B[7]+B[8]+B[10]+2.0*B[11]+2.0*B[12]+2.0*B[13]+2.0*B[14]))/s0[5]
      mean3E[6]<- (sumwx_F2[6]-sigma3E[6]*B[8])/s0[6]
      mean3E[7]<- (sumwx_F2[7]-sigma3E[7]*(B[3]-B[4]+B[13]+B[15]))/s0[7]
      mean3E[8]<- (sumwx_F2[8]+sigma3E[8]*(B[9]-B[12]+2.0*B[15]))/s0[8]
      mean3E[9]<- (sumwx_F2[9]-sigma3E[9]*(B[9]+B[10]+B[11]+B[15]))/s0[9]

      aaa1<- max(abs(AA-B))
      AA<- B
      if (n_iter>20) break
    }
    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2)
    ss2<- sum((dataF1-meanF1)^2)
    ss3<- sum((dataP2-meanP2)^2)

    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1E[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2E[i])^2 }
    for(i in 1:d23) {  swx_F2[i] <- WW_F2[i,]%*%(dataF2-mean3E[i])^2 }
    aa1<-sum(swx_B1);aa2<-sum(swx_B2);aa3<-sum(swx_F2)
    sigma1E<- matrix(aa1/n_samB1,d21,1)
    sigma40<- sigma1E[1]-sigma
    sigma2E<- matrix(aa2/n_samB2,d22,1)
    sigma50<- sigma2E[1]-sigma
    sigma3E<- matrix(aa3/n_samF2,d23,1)
    sigma60<- sigma3E[1]-sigma
    ################ iteratively CM3-step for variance (sigma) ################
    sigma40[(sigma40<0)]<- 0.00001; sigma50[(sigma50<0)]<- 0.00001; sigma60[(sigma60<0)]<- 0.00001
    ab1<- ss1+ss2+ss3; ab2<- n_samP1+n_samF1+n_samP2
    aaa0<- sigma ;n_iter<- 0; aaa1<-1000
    while (aaa1>0.0001){
      n_iter<- n_iter+1
      aa4<- sigma/(sigma+sigma40); aa4[(aa4>1.0)]<- 1
      aa5<- sigma/(sigma+sigma50); aa5[(aa5>1.0)]<- 1
      aa6<- sigma/(sigma+sigma60); aa6[(aa6>1.0)]<- 1
      sigma<- (ab1+aa4*aa4*aa1+aa5*aa5*aa2+aa6*aa6*aa3)/(ab2+aa4*n_samB1+aa5*n_samB2+aa6*n_samF2);
      aaa1<- abs(sigma-aaa0)
      aaa0<- sigma
      if (n_iter>20) break
    }
    sigma1E<- matrix(sigma+sigma40,d21,1);sigma2E<- matrix(sigma+sigma50,d22,1);sigma3E<- matrix(sigma+sigma60,d23,1)
    ####################### the stop criterion for iteration ############################################################
    L1<- logL(n_samP1,1,1,meanP1,sigma,dataP1)+logL(n_samF1,1,1,meanF1,sigma,dataF1)+logL(n_samP2,1,1,meanP2,sigma,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1E,sigma1E,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2E,sigma2E,dataB2)+logL(n_samF2,d23,mix_pi_3,mean3E,sigma3E,dataF2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*9

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, 1,0,-1,1,1,0,0,0,0,-1,-1,1,1,1,0,0,0,-1,-1,-1, 1,0,-1,1,0,1,0,0,-1,0,-1,1,0,-1,1,0,-1,1,0,-1,
                1,0,-1,0.5,0.5,0.5,0.5,-0.5,-0.5,-0.5,-0.5,0,0,0,0,0,0,0,0,0, 0,1,0,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5 ),20,5)

  b_line1 <- matrix(c(meanP1,meanF1,meanP2,mean1E,mean2E,mean3E))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigmaB1 - sigma1E[1]
  if(jj_1 < 0 | jj_1>=sigmaB1) {jj_1 <- 0}
  ll_1 <- jj_1/sigmaB1
  mm1<- sigma1E[1]-sigma
  if (mm1<0 | mm1>=sigmaB1) { mm1<- 0 }
  nn1<- mm1/sigmaB1

  jj_2 <- sigmaB2 - sigma2E[1]
  if(jj_2 < 0 | jj_2>=sigmaB2) {jj_2 <- 0}
  ll_2 <- jj_2/sigmaB2
  mm2<- sigma2E[1]-sigma
  if (mm2<0 | mm2>=sigmaB2) { mm2<- 0 }
  nn2<- mm2/sigmaB2

  jj_3 <- sigmaF2 - sigma3E[1]
  if(jj_3 < 0 | jj_3>=sigmaF2) {jj_3 <- 0}
  ll_3 <- jj_3/sigmaF2
  mm3<- sigma3E[1]-sigma
  if (mm3<0 | mm3>=sigmaF2) { mm3<- 0 }
  nn3<- mm3/sigmaF2

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
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

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
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),(1-pkolm(D_F1,n_samF1))))

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
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

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
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),(1-pkolm(D_B1,n_samB1))))

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
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),(1-pkolm(D_B2,n_samB2))))

  ###############  F2  ###################
  dataF2<-sort(dataF2);bmw_F2 <- matrix(0,n_samF2,1); bmwsl_F2 <- matrix(0,n_samF2,d23)
  for(i in 1:d23){
    gg_F2 <- (dataF2 - mean3E[i])/sqrt(sigma3E[i])
    bmw_F2[which(gg_F2>=0)] <- pnorm(gg_F2[gg_F2>=0])
    bmw_F2[which(gg_F2<0)] <- 1 - pnorm(abs(gg_F2[gg_F2<0]))
    bmwsl_F2[,i] <- bmw_F2*mix_pi_3[i]
  }
  P2_F2 <- rowSums(bmwsl_F2)
  nn<-dim(as.matrix(unique(P2_F2)))[1]
  if(nn<n_samF2){P2_F2<-P2_F2+runif(n_samF2)/1e4}

  dd_F2 <- as.matrix(c(sum(P2_F2),sum(P2_F2^2),sum((P2_F2-0.5)^2)))
  WW2_F2 <- 1/(12*n_samF2) + sum((P2_F2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  u_F2 <- as.matrix(c(12*n_samF2*((dd_F2[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((dd_F2[2]/n_samF2-1/3)^2),180*n_samF2*((dd_F2[3]/n_samF2-1/12)^2)))
  D_F2 <- as.numeric(ks.test(P2_F2,"punif")[[1]][1])
  tt_F2 <- as.matrix(c((1 - pchisq(u_F2[1],1)),(1 - pchisq(u_F2[2],1)),(1 - pchisq(u_F2[3],1)),K1(WW2_F2),(1-pkolm(D_F2,n_samF2))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)
  tt_F2[which(tt_F2>=10e-4)]<-round(tt_F2[which(tt_F2>=10e-4)],4);tt_F2[which(tt_F2<10e-4)]<-format(tt_F2[which(tt_F2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame(" MX2-A-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigma,4),round(t(mean1E),4),round(sigma1E[1],4),round(t(mix_pi_1),4),
                       round(t(mean2E),4),round(sigma2E[1],4),round(t(mix_pi_2),4),round(t(mean3E),4),round(sigma3E[1],4),round(t(mix_pi_3),4),
                       round(B1[1],4)," "," "," "," "," ",round(B1[2],4),round(B1[3],4)," "," "," "," "," "," ",round(B1[4],4),round(B1[5],4),round(jj_1,4),round(ll_1*100,4),round(mm1,4),round(nn1*100,4),
                       round(jj_2,4),round(ll_2*100,4),round(mm2,4),round(nn2*100,4),round(jj_3,4),round(ll_3*100,4),round(mm3,4),round(nn3*100,4),round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),
                       tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),
                       tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),
                       tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5],round(u_F2[1],4),tt_F2[1],round(u_F2[2],4),
                       tt_F2[2],round(u_F2[3],4),tt_F2[3],round(WW2_F2,4),tt_F2[4],round(D_F2,4),tt_F2[5])

  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2,mi_3)
  return(OUTPUT)
}
############################################ MX2-EA-AD Model #########################################
G6ModelFun[[22]] <- function(K1,logL,df11,df21,df31,df41,df51,df61){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1]; n_samF1<-dim(dataF1)[1]; n_samP2<-dim(dataP2)[1]; n_samB1<-dim(dataB1)[1]; n_samB2<-dim(dataB2)[1]; n_samF2<-dim(dataF2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2)
  sumx4<-sum(dataB1);sumx5<-sum(dataB2);sumx6<-sum(dataF2)
  sigmaP1<- as.numeric(var(dataP1)); ss1<- (n_samP1-1)*sigmaP1
  sigmaF1<- as.numeric(var(dataF1)); ss2<- (n_samF1-1)*sigmaF1
  sigmaP2<- as.numeric(var(dataP2)); ss3<- (n_samP2-1)*sigmaP2
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  mean4<- mean(dataB1); sigmaB1<- as.numeric(var(dataB1)); sigma40<- sigmaB1
  mean5<- mean(dataB2); sigmaB2<- as.numeric(var(dataB2)); sigma50<- sigmaB2
  mean6<- mean(dataF2); sigmaF2<- as.numeric(var(dataF2)); sigma60<- sigmaF2

  m_esp<-0.0001
  ####################### MX2-EA-AD Model #################### (E4) ################
  d21<-3; d22<-3; d23<-5
  sigma<- sigma0
  meanP1<- mean(dataP1);meanF1<- mean(dataF1);meanP2<- mean(dataP2)
  abb<- sigmaB1
  abb[sigmaB2>abb]<- sigmaB2; abb[sigmaF2>abb]<- sigmaF2
  abb<- abb/sigmaB1

  mi_1<- as.matrix(c(0.25,0.5,0.25))
  a1<- sqrt(sigmaB1)
  if (meanP1<meanP2) {a1= -a1}
  sigma1E<- matrix(sigmaB1,d21,1)
  mean1E<- as.matrix(c(mean4+2*a1,mean4,mean4-2*a1))

  mi_2<- as.matrix(c(0.25,0.5,0.25))
  a2<- sqrt(sigmaB2)
  if (meanP1<meanP2) { a2= -a2 }
  sigma2E<- matrix(sigmaB2/abb,d22,1)
  mean2E<- as.matrix(c(mean5+2*a2,mean5,mean5-2*a2))

  mi_3<- as.matrix(c(0.0625,0.25,0.375,0.25,0.0625))
  a3<- sqrt(sigmaF2)
  if (meanP1<meanP2) { a3= -a3 }
  sigma3E<- matrix(sigmaF2/abb,d23,1)
  mean3E<- as.matrix(c(mean6+3*a3,mean6+1.5*a3,mean6,mean6-1.5*a3,mean6-3*a3))

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1E,sigma1E,dataB1)+logL(n_samB2,d22,mi_2,mean2E,sigma2E,dataB2)+logL(n_samF2,d23,mi_3,mean3E,sigma3E,dataF2)

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)
  WW_F2 <- matrix(0,d23,n_samF2); swx_F2 <- matrix(0,d23,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1E[i],sqrt(sigma1E[i]))/dmixnorm(dataB1,mean1E,sqrt(sigma1E),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2E[i],sqrt(sigma2E[i]))/dmixnorm(dataB2,mean2E,sqrt(sigma2E),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2

    for(i in 1:d23) { WW_F2[i,] <- mi_3[i]*dnorm(dataF2,mean3E[i],sqrt(sigma3E[i]))/dmixnorm(dataF2,mean3E,sqrt(sigma3E),mi_3) }
    mix_pi_3 <- as.matrix(rowSums(WW_F2)/n_samF2)
    s0<-as.matrix(rowSums(WW_F2));s0[abs(s0)<0.0001]<- 0.00001
    sumwx_F2 <- WW_F2%*%dataF2

    ############ CM1-step for means ##############
    n0<- as.matrix(c(mix_pi_1*n_samB1,mix_pi_2*n_samB2))
    n0[abs(n0)<0.00000001] <- 0.000001

    aaa0<- 0; aaa1<- 1000; n_iter<- 0; AA<- matrix(0,10,1)
    while (aaa1>0.0001){
      n_iter<- n_iter+1
      hh<- matrix(0,10,10)
      hh[1,1]<- sigma*(1.0/n_samP1+4.0/n_samF1+1.0/n_samP2)+16.0*sigma3E[3]/s0[3]
      hh[1,2]<- sigma*(1.0/n_samP1+4.0/n_samF1+1.0/n_samP2)
      hh[1,3]<- sigma*(1.0/n_samP1+4.0/n_samF1+1.0/n_samP2)+8.0*sigma3E[3]/s0[3]
      hh[1,4]<- sigma/n_samP1-sigma/n_samP2
      hh[1,5]<- sigma/n_samP1-sigma/n_samP2
      hh[1,6]<- 0
      hh[1,7]<- -4.0*sigma3E[3]/s0[3]
      hh[1,8]<- 4.0*sigma3E[3]/s0[3]
      hh[1,9]<- 0
      hh[1,10]<- 12.0*sigma3E[3]/s0[3]
      hh[2,2]<- sigma*(1.0/n_samP1+4.0/n_samF1+1.0/n_samP2)+4.0*sigma3E[1]/s0[1]+4.0*sigma3E[5]/s0[5]
      hh[2,3]<- sigma*(1.0/n_samP1+4.0/n_samF1+1.0/n_samP2)
      hh[2,4]<- sigma*(1.0/n_samP1-1.0/n_samP2)
      hh[2,5]<- sigma*(1.0/n_samP1-1.0/n_samP2)
      hh[2,6]<- 2.0*sigma3E[1]/s0[1]
      hh[2,7]<- 0
      hh[2,8]<- 0
      hh[2,9]<- -2.0*sigma3E[5]/s0[5]
      hh[2,10]<- 2.0*sigma3E[5]/s0[5]
      hh[3,3]<- sigma*(1.0/n_samP1+4.0/n_samF1+1.0/n_samP2)+sigma1E[1]/n0[1]+sigma2E[3]/n0[6]+4.0*sigma3E[3]/s0[3]
      hh[3,4]<- sigma*(1.0/n_samP1-1.0/n_samP2)+sigma1E[1]/n0[1]-sigma2E[3]/n0[6]
      hh[3,5]<- sigma*(1.0/n_samP1-1.0/n_samP2)
      hh[3,6]<- -sigma1E[1]/n0[1]
      hh[3,7]<- -2.0*sigma3E[3]/s0[3]
      hh[3,8]<- 2.0*sigma3E[3]/s0[3]
      hh[3,9]<- sigma2E[3]/n0[6]
      hh[3,10]<- 6.0*sigma3E[3]/s0[3]
      hh[4,4]<- sigma*(1.0/n_samP1+1.0/n_samP2)+sigma1E[1]/n0[1]+sigma1E[3]/n0[3]+sigma2E[1]/n0[4]+sigma2E[3]/n0[6]
      hh[4,5]<- sigma*(1.0/n_samP1+1.0/n_samP2)
      hh[4,6]<- -sigma1E[1]/n0[1]
      hh[4,7]<- sigma1E[3]/n0[3]
      hh[4,8]<- sigma2E[1]/n0[4]
      hh[4,9]<- -sigma2E[3]/n0[6]
      hh[4,10]<- 0
      hh[5,5]<- sigma*(1.0/n_samP1+1.0/n_samP2)+4.0*sigma1E[2]/n0[2]+4.0*sigma2E[2]/n0[5]
      hh[5,6]<- 2.0*sigma1E[2]/n0[2]
      hh[5,7]<- -2.0*sigma1E[2]/n0[2]
      hh[5,8]<- -2.0*sigma2E[2]/n0[5]
      hh[5,9]<- 2.0*sigma2E[2]/n0[5]
      hh[5,10]<- 0
      hh[6,6]<- sigma1E[1]/n0[1]+sigma1E[2]/n0[2]+sigma3E[1]/s0[1]+sigma3E[2]/s0[2]
      hh[6,7]<- -sigma1E[2]/n0[2]-sigma3E[2]/s0[2]
      hh[6,8]<- hh[6,9]<- 0
      hh[6,10]<- sigma3E[2]/s0[2]
      hh[7,7]<- sigma1E[2]/n0[2]+sigma1E[3]/n0[3]+sigma3E[2]/s0[2]+sigma3E[3]/s0[3]
      hh[7,8]<- -sigma3E[3]/s0[3]
      hh[7,9]<- 0
      hh[7,10]<- -sigma3E[2]/s0[2]-3.0*sigma3E[3]/s0[3]
      hh[8,8]<- sigma2E[1]/n0[4]+sigma2E[2]/n0[5]+sigma3E[3]/s0[3]+sigma3E[4]/s0[4]
      hh[8,9]<- -sigma2E[2]/n0[5]-sigma3E[4]/s0[4]
      hh[8,10]<- 3.0*sigma3E[3]/s0[3]+3.0*sigma3E[4]/s0[4]
      hh[9,9]<- sigma2E[2]/n0[5]+sigma2E[3]/n0[6]+sigma3E[4]/s0[4]+sigma3E[5]/s0[5]
      hh[9,10]<- -3.0*sigma3E[4]/s0[4]-sigma3E[5]/s0[5]
      hh[10,10]<- sigma3E[2]/s0[2]+9.0*sigma3E[3]/s0[3]+9.0*sigma3E[4]/s0[4]+sigma3E[5]/s0[5]
      for(i in 2:10){
        for(j in 1:(i-1)){
          hh[i,j]<- hh[j,i]
        }
      }
      b_line<- matrix(0,10,1)
      b_line[1]<- sumx1/n_samP1+2.0*sumx2/n_samF1+sumx3/n_samP2-4.0*sumwx_F2[3]/s0[3]
      b_line[2]<- sumx1/n_samP1+2.0*sumx2/n_samF1+sumx3/n_samP2-2.0*sumwx_F2[1]/s0[1]-2.0*sumwx_F2[5]/s0[5]
      b_line[3]<- sumx1/n_samP1+2.0*sumx2/n_samF1+sumx3/n_samP2-sumwx_B1[1]/n0[1]-sumwx_B2[3]/n0[6]-2.0*sumwx_F2[3]/s0[3]
      b_line[4]<- sumx1/n_samP1-sumx3/n_samP2-sumwx_B1[1]/n0[1]-sumwx_B1[3]/n0[3]+sumwx_B2[1]/n0[4]+sumwx_B2[3]/n0[6]
      b_line[5]<- sumx1/n_samP1-sumx3/n_samP2-2.0*sumwx_B1[2]/n0[2]+2.0*sumwx_B2[2]/n0[5]

      b_line[6]<- sumwx_B1[1]/n0[1]-sumwx_B1[2]/n0[2]-sumwx_F2[1]/s0[1]+sumwx_F2[2]/s0[2]
      b_line[7]<- sumwx_B1[2]/n0[2]-sumwx_B1[3]/n0[3]-sumwx_F2[2]/s0[2]+sumwx_F2[3]/s0[3]
      b_line[8]<- sumwx_B2[1]/n0[4]-sumwx_B2[2]/n0[5]-sumwx_F2[3]/s0[3]+sumwx_F2[4]/s0[4]
      b_line[9]<- sumwx_B2[2]/n0[5]-sumwx_B2[3]/n0[6]-sumwx_F2[4]/s0[4]+sumwx_F2[5]/s0[5]
      b_line[10]<- sumwx_F2[2]/s0[2]-3.0*sumwx_F2[3]/s0[3]+3.0*sumwx_F2[4]/s0[4]-sumwx_F2[5]/s0[5]

      B <- solve(hh,b_line)

      meanP1<- (sumx1-sigma*(B[1]+B[2]+B[3]+B[4]+B[5]))/n_samP1
      meanF1<- (sumx2-sigma*2.0*(B[1]+B[2]+B[3]))/n_samF1
      meanP2<- (sumx3-sigma*(B[1]+B[2]+B[3]-B[4]-B[5]))/n_samP2

      mean1E[1]<- (sumwx_B1[1]+sigma1E[1]*(B[3]+B[4]-B[6]))/n0[1]
      mean1E[2]<- (sumwx_B1[2]+sigma1E[2]*(2.0*B[5]+B[6]-B[7]))/n0[2]
      mean1E[3]<- (sumwx_B1[3]+sigma1E[3]*(B[4]+B[7]))/n0[3]

      mean2E[1]<- (sumwx_B2[1]-sigma2E[1]*(B[4]+B[8]))/n0[4]
      mean2E[2]<- (sumwx_B2[2]-sigma2E[2]*(2.0*B[5]-B[8]+B[9]))/n0[5]
      mean2E[3]<- (sumwx_B2[3]+sigma2E[3]*(B[3]-B[4]+B[9]))/n0[6]

      mean3E[1]<- (sumwx_F2[1]+sigma3E[1]*(2.0*B[2]+B[6]))/s0[1]
      mean3E[2]<- (sumwx_F2[2]-sigma3E[2]*(B[6]-B[7]+B[10]))/s0[2]
      mean3E[3]<- (sumwx_F2[3]+sigma3E[3]*(4.0*B[1]+2.0*B[3]-B[7]+B[8]+3.0*B[10]))/s0[3]
      mean3E[4]<- (sumwx_F2[4]-sigma3E[4]*(B[8]-B[9]+3.0*B[10]))/s0[4]
      mean3E[5]<- (sumwx_F2[5]+sigma3E[5]*(2.0*B[2]-B[9]+B[10]))/s0[5]

      aaa1<- max(abs(AA-B))
      AA<- B
      if (n_iter>20) break
    }
    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2)
    ss2<- sum((dataF1-meanF1)^2)
    ss3<- sum((dataP2-meanP2)^2)

    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1E[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2E[i])^2 }
    for(i in 1:d23) {  swx_F2[i] <- WW_F2[i,]%*%(dataF2-mean3E[i])^2 }

    aa1<-sum(swx_B1);aa2<-sum(swx_B2);aa3<-sum(swx_F2)
    sigma1E<- matrix(aa1/n_samB1,d21,1)
    sigma40<- sigma1E[1]-sigma
    sigma2E<- matrix(aa2/n_samB2,d22,1)
    sigma50<- sigma2E[1]-sigma
    sigma3E<- matrix(aa3/n_samF2,d23,1)
    sigma60<- sigma3E[1]-sigma
    ################ iteratively CM3-step for variance (sigma) ################
    sigma40[(sigma40<0)]<- 0.00001; sigma50[(sigma50<0)]<- 0.00001; sigma60[(sigma60<0)]<- 0.00001
    ab1<- ss1+ss2+ss3; ab2<- n_samP1+n_samF1+n_samP2
    aaa0<- sigma ;n_iter<- 0;aaa1<- 1000
    while (aaa1>0.0001){
      n_iter=n_iter+1
      aa4<- sigma/(sigma+sigma40); aa4[(aa4>1.0)]<- 1
      aa5<- sigma/(sigma+sigma50); aa5[(aa5>1.0)]<- 1
      aa6<- sigma/(sigma+sigma60); aa6[(aa6>1.0)]<- 1
      sigma<- (ab1+aa4*aa4*aa1+aa5*aa5*aa2+aa6*aa6*aa3)/(ab2+aa4*n_samB1+aa5*n_samB2+aa6*n_samF2)
      aaa1<- abs(sigma-aaa0)
      aaa0<- sigma
      if (n_iter>20) break
    }
    sigma1E<- matrix(sigma+sigma40,d21,1);sigma2E<- matrix(sigma+sigma50,d22,1);sigma3E<- matrix(sigma+sigma60,d23,1)
    ####################### the stop criterion for iteration ############################################################
    L1<- logL(n_samP1,1,1,meanP1,sigma,dataP1)+logL(n_samF1,1,1,meanF1,sigma,dataF1)+logL(n_samP2,1,1,meanP2,sigma,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1E,sigma1E,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2E,sigma2E,dataB2)+logL(n_samF2,d23,mix_pi_3,mean3E,sigma3E,dataF2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*8

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1, 2,0,-2,2,1,0,0,-1,-2,2,1,0,-1,-2,
                1,0,-1,0.5,0.5,0.5,-0.5,-0.5,-0.5,0,0,0,0,0, 0,1,0,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5),14,4)

  b_line1 <- matrix(c(meanP1,meanF1,meanP2,mean1E,mean2E,mean3E))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigmaB1 - sigma1E[1]
  if(jj_1 < 0 | jj_1>=sigmaB1) {jj_1 <- 0}
  ll_1 <- jj_1/sigmaB1
  mm1<- sigma1E[1]-sigma
  if (mm1<0 | mm1>=sigmaB1) { mm1<- 0 }
  nn1<- mm1/sigmaB1

  jj_2 <- sigmaB2 - sigma2E[1]
  if(jj_2 < 0 | jj_2>=sigmaB2) {jj_2 <- 0}
  ll_2 <- jj_2/sigmaB2
  mm2<- sigma2E[1]-sigma
  if (mm2<0 | mm2>=sigmaB2) { mm2<- 0 }
  nn2<- mm2/sigmaB2

  jj_3 <- sigmaF2 - sigma3E[1]
  if(jj_3 < 0 | jj_3>=sigmaF2) {jj_3 <- 0}
  ll_3 <- jj_3/sigmaF2
  mm3<- sigma3E[1]-sigma
  if (mm3<0 | mm3>=sigmaF2) { mm3<- 0 }
  nn3<- mm3/sigmaF2

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
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

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
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),(1-pkolm(D_F1,n_samF1))))

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
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

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
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),(1-pkolm(D_B1,n_samB1))))

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
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),(1-pkolm(D_B2,n_samB2))))

  ###############  F2  ###################
  dataF2<-sort(dataF2);bmw_F2 <- matrix(0,n_samF2,1); bmwsl_F2 <- matrix(0,n_samF2,d23)
  for(i in 1:d23){
    gg_F2 <- (dataF2 - mean3E[i])/sqrt(sigma3E[i])
    bmw_F2[which(gg_F2>=0)] <- pnorm(gg_F2[gg_F2>=0])
    bmw_F2[which(gg_F2<0)] <- 1 - pnorm(abs(gg_F2[gg_F2<0]))
    bmwsl_F2[,i] <- bmw_F2*mix_pi_3[i]
  }
  P2_F2 <- rowSums(bmwsl_F2)
  nn<-dim(as.matrix(unique(P2_F2)))[1]
  if(nn<n_samF2){P2_F2<-P2_F2+runif(n_samF2)/1e4}

  dd_F2 <- as.matrix(c(sum(P2_F2),sum(P2_F2^2),sum((P2_F2-0.5)^2)))
  WW2_F2 <- 1/(12*n_samF2) + sum((P2_F2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  u_F2 <- as.matrix(c(12*n_samF2*((dd_F2[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((dd_F2[2]/n_samF2-1/3)^2),180*n_samF2*((dd_F2[3]/n_samF2-1/12)^2)))
  D_F2 <- as.numeric(ks.test(P2_F2,"punif")[[1]][1])
  tt_F2 <- as.matrix(c((1 - pchisq(u_F2[1],1)),(1 - pchisq(u_F2[2],1)),(1 - pchisq(u_F2[3],1)),K1(WW2_F2),(1-pkolm(D_F2,n_samF2))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)
  tt_F2[which(tt_F2>=10e-4)]<-round(tt_F2[which(tt_F2>=10e-4)],4);tt_F2[which(tt_F2<10e-4)]<-format(tt_F2[which(tt_F2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame(" MX2-EA-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigma,4),round(t(mean1E),4)," ",round(sigma1E[1],4),round(t(mix_pi_1),4)," ",
                       round(t(mean2E),4)," ",round(sigma2E[1],4),round(t(mix_pi_2),4)," ",round(t(mean3E),4)," "," "," "," ",round(sigma3E[1],4),round(t(mix_pi_3),4)," "," "," "," ",
                       round(B1[1],4)," "," "," "," "," ",round(B1[2],4),round(B1[2],4)," "," "," "," "," "," ",round(B1[3],4),round(B1[4],4),round(jj_1,4),round(ll_1*100,4),round(mm1,4),round(nn1*100,4),
                       round(jj_2,4),round(ll_2*100,4),round(mm2,4),round(nn2*100,4),round(jj_3,4),round(ll_3*100,4),round(mm3,4),round(nn3*100,4),round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),
                       tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),
                       tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),
                       tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5],round(u_F2[1],4),tt_F2[1],round(u_F2[2],4),
                       tt_F2[2],round(u_F2[3],4),tt_F2[3],round(WW2_F2,4),tt_F2[4],round(D_F2,4),tt_F2[5])

  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2,mi_3)
  return(OUTPUT)
}

############################################ MX2-CD-AD Model #########################################
G6ModelFun[[23]] <- function(K1,logL,df11,df21,df31,df41,df51,df61){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1]; n_samF1<-dim(dataF1)[1]; n_samP2<-dim(dataP2)[1]; n_samB1<-dim(dataB1)[1]; n_samB2<-dim(dataB2)[1]; n_samF2<-dim(dataF2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2)
  sumx4<-sum(dataB1);sumx5<-sum(dataB2);sumx6<-sum(dataF2)
  sigmaP1<- as.numeric(var(dataP1)); ss1<- (n_samP1-1)*sigmaP1
  sigmaF1<- as.numeric(var(dataF1)); ss2<- (n_samF1-1)*sigmaF1
  sigmaP2<- as.numeric(var(dataP2)); ss3<- (n_samP2-1)*sigmaP2
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  mean4<- mean(dataB1); sigmaB1<- as.numeric(var(dataB1)); sigma40<- sigmaB1
  mean5<- mean(dataB2); sigmaB2<- as.numeric(var(dataB2)); sigma50<- sigmaB2
  mean6<- mean(dataF2); sigmaF2<- as.numeric(var(dataF2)); sigma60<- sigmaF2

  m_esp<-0.0001
  ####################### MX2-CD-AD Model #################### (E5) ###############
  d21<-1; d22<-4; d23<-4
  sigma<- sigma0
  meanP1<- mean(dataP1);meanF1<- mean(dataF1);meanP2<- mean(dataP2)
  abb<- sigmaF2/(5*sigma0)

  mi_1<- as.matrix(1)
  sigma1E<- as.matrix(sigmaB1/abb)
  mean1E<- as.matrix(mean4)

  mi_2<- matrix(0.25,4,1)
  a2<- sqrt(sigmaB2/n_samB2)
  if (meanP1<meanP2) { a2= -a2 }
  sigma2E<- matrix(sigmaB2/abb,d22,1)
  mean2E<- as.matrix(c(mean5+3*a2,mean5+a2,mean5-a2,mean5-3*a2))

  mi_3<- as.matrix(c(0.5625,0.1875,0.1875,0.0625))
  a3<- sqrt(sigmaF2/n_samF2)
  if (meanP1<meanP2) { a3= -a3 }
  sigma3E<- matrix(sigmaF2/abb,d23,1)
  mean3E<- as.matrix(c(mean6+3*a3,mean6+a3,mean6-a3,mean6-3*a3))

  sigmaP1<- sigma0;sigmaF1<- sigma0;sigmaP2<- sigma0

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1E,sigma1E,dataB1)+logL(n_samB2,d22,mi_2,mean2E,sigma2E,dataB2)+logL(n_samF2,d23,mi_3,mean3E,sigma3E,dataF2)

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)
  WW_F2 <- matrix(0,d23,n_samF2); swx_F2 <- matrix(0,d23,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1E[i],sqrt(sigma1E[i]))/dmixnorm(dataB1,mean1E,sqrt(sigma1E),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2E[i],sqrt(sigma2E[i]))/dmixnorm(dataB2,mean2E,sqrt(sigma2E),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    n0<-as.matrix(rowSums(WW_B2));n0[abs(n0)<0.0001]<- 0.00001
    sumwx_B2 <- WW_B2%*%dataB2

    for(i in 1:d23) { WW_F2[i,] <- mi_3[i]*dnorm(dataF2,mean3E[i],sqrt(sigma3E[i]))/dmixnorm(dataF2,mean3E,sqrt(sigma3E),mi_3) }
    mix_pi_3 <- as.matrix(rowSums(WW_F2)/n_samF2)
    s0<-as.matrix(rowSums(WW_F2));s0[abs(s0)<0.0001]<- 0.00001
    sumwx_F2 <- WW_F2%*%dataF2

    ###### CM1-step for means.
    aaa0<- 0; aaa1<- 1000; n_iter<- 0; AA<- matrix(0,7,1)
    while (aaa1>0.0001){
      n_iter<- n_iter+1
      hh<- matrix(0,7,7)
      hh[1,1]<- sigma*(1.0/n_samP1+1.0/n_samP2) + sigma1E[1]/n_samB1+sigma2E[4]/n0[4]+4.0*sigma3E[1]/s0[1]
      hh[1,2]<- sigma/n_samP1+2.0*sigma/n_samF1- sigma/n_samP2 +2*sigma1E[1]/n_samB1 - sigma2E[4]/n0[4]
      hh[1,3]<- hh[1,4]<- 0
      hh[1,5]<- sigma2E[4]/n0[4]+2.0*sigma3E[1]/s0[1]
      hh[1,6]<- -sigma1E[1]/n_samB1+2.0*sigma3E[1]/s0[1]
      hh[1,7]<- -2.0*sigma3E[1]/s0[1]
      hh[2,2]<- sigma*(1.0/n_samP1+1.0/n_samP2)+4.0*sigma1E[1]/n_samB1+sigma2E[1]/n0[1]+sigma2E[4]/n0[4]
      hh[2,3]<- sigma2E[1]/n0[1]
      hh[2,4]<- sigma2E[1]/n0[1]
      hh[2,5]<- sigma2E[1]/n0[1]-sigma2E[4]/n0[4]
      hh[2,6]<- -2.0*sigma1E[1]/n_samB1
      hh[2,7]<- 0
      hh[3,3]<- sigma2E[1]/n0[1]+sigma2E[2]/n0[2]+sigma3E[3]/s0[3]+sigma3E[4]/s0[4]
      hh[3,4]<- sigma2E[1]/n0[1]+sigma3E[4]/s0[4]
      hh[3,5]<- sigma2E[1]/n0[1]+sigma3E[4]/s0[4]
      hh[3,6]<- sigma3E[3]/s0[3]
      hh[3,7]<- sigma3E[3]/s0[3]+sigma3E[4]/s0[4]
      hh[4,4]<- sigma2E[1]/n0[1]+sigma2E[3]/n0[3]+sigma3E[2]/s0[2]+sigma3E[4]/s0[4]
      hh[4,5]<- sigma2E[1]/n0[1]+sigma3E[4]/s0[4]
      hh[4,6]<- -sigma2E[3]/n0[3]
      hh[4,7]<- sigma3E[2]/s0[2]+sigma3E[4]/s0[4]
      hh[5,5]<- sigma2E[1]/n0[1]+sigma2E[4]/n0[4]+sigma3E[1]/s0[1]+sigma3E[4]/s0[4]
      hh[5,6]<- sigma3E[1]/s0[1]
      hh[5,7]<- -sigma3E[1]/s0[1]+sigma3E[4]/s0[4]
      hh[6,6]<- sigma1E[1]/n_samB1+sigma2E[3]/n0[3]+sigma3E[1]/s0[1]+sigma3E[3]/s0[3]
      hh[6,7]<- -sigma3E[1]/s0[1]+sigma3E[3]/s0[3]
      hh[7,7]<- sigma3E[1]/s0[1]+sigma3E[2]/s0[2]+sigma3E[3]/s0[3]+sigma3E[4]/s0[4]
      for(i in 2:7){
        for(j in 1:(i-1)){
          hh[i,j]<- hh[j,i]
        }
      }
      b_line<- matrix(0,7,1)
      b_line[1]<- sumx1/n_samP1+2.0*sumx2/n_samF1+sumx3/n_samP2-sumwx_B1[1]/n_samB1-sumwx_B2[4]/n0[4]-2.0*sumwx_F2[1]/s0[1]
      b_line[2]<- sumx1/n_samP1-sumx3/n_samP2-2.0*sumwx_B1[1]/n_samB1+sumwx_B2[1]/n0[1]+sumwx_B2[4]/n0[4]
      b_line[3]<- sumwx_B2[1]/n0[1]-sumwx_B2[2]/n0[2]-sumwx_F2[3]/s0[3]+sumwx_F2[4]/s0[4]
      b_line[4]<- sumwx_B2[1]/n0[1]-sumwx_B2[3]/n0[3]-sumwx_F2[2]/s0[2]+sumwx_F2[4]/s0[4]
      b_line[5]<- sumwx_B2[1]/n0[1]-sumwx_B2[4]/n0[4]-sumwx_F2[1]/s0[1]+sumwx_F2[4]/s0[4]
      b_line[6]<- sumwx_B1[1]/n_samB1+sumwx_B2[3]/n0[3]-sumwx_F2[1]/s0[1]-sumwx_F2[3]/s0[3]
      b_line[7]<- sumwx_F2[1]/s0[1]-sumwx_F2[2]/s0[2]-sumwx_F2[3]/s0[3]+sumwx_F2[4]/s0[4]

      B <- solve(hh,b_line)

      meanP1<- (sumx1-sigma*(B[1]+B[2]))/n_samP1
      meanF1<- (sumx2-sigma*2.0*B[2])/n_samF1
      meanP2<- (sumx3-sigma*(B[1]-B[2]))/n_samP2

      mean1E[1]<- (sumwx_B1[1]+sigma1E[1]*(B[1]+2.0*B[2]-B[6]))/n_samB1

      mean2E[1]<- (sumwx_B2[1]-sigma2E[1]*(B[2]+B[3]+B[4]+B[5]))/n0[1]
      mean2E[2]<- (sumwx_B2[2]+sigma2E[2]*B[3])/n0[2]
      mean2E[3]<- (sumwx_B2[3]+sigma2E[3]*(B[4]-B[6]))/n0[3]
      mean2E[4]<- (sumwx_B2[4]+sigma2E[4]*(B[1]-B[2]+B[5]))/n0[4]

      mean3E[1]<- (sumwx_F2[1]+sigma3E[1]*(2.0*B[1]+B[5]+B[6]-B[7]))/s0[1]
      mean3E[2]<- (sumwx_F2[2]+sigma3E[2]*(B[4]+B[7]))/s0[2]
      mean3E[3]<- (sumwx_F2[3]+sigma3E[3]*(B[3]+B[6]+B[7]))/s0[3]
      mean3E[4]<- (sumwx_F2[4]-sigma3E[4]*(B[3]+B[4]+B[5]+B[7]))/s0[4]

      aaa1<- max(abs(AA-B))
      AA<- B
      if (n_iter>20) break
    }
    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2)
    ss2<- sum((dataF1-meanF1)^2)
    ss3<- sum((dataP2-meanP2)^2)

    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1E[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2E[i])^2 }
    for(i in 1:d23) {  swx_F2[i] <- WW_F2[i,]%*%(dataF2-mean3E[i])^2 }

    aa1<-sum(swx_B1);aa2<-sum(swx_B2);aa3<-sum(swx_F2)
    sigma1E<- matrix(aa1/n_samB1,d21,1)
    sigma40<- sigma1E[1]-sigma
    sigma2E<- matrix(aa2/n_samB2,d22,1)
    sigma50<- sigma2E[1]-sigma
    sigma3E<- matrix(aa3/n_samF2,d23,1)
    sigma60<- sigma3E[1]-sigma
    ################ iteratively CM3-step for variance (sigma) ################
    sigma40[(sigma40<0)]<- 0.00001; sigma50[(sigma50<0)]<- 0.00001; sigma60[(sigma60<0)]<- 0.00001
    ab1<- ss1+ss2+ss3; ab2<- n_samP1+n_samF1+n_samP2
    aaa0<- sigma ;n_iter<- 0;aaa1<- 1000
    while (aaa1>0.0001){
      n_iter<- n_iter+1
      aa4<- sigma/(sigma+sigma40); aa4[(aa4>1.0)]<- 1
      aa5<- sigma/(sigma+sigma50); aa5[(aa5>1.0)]<- 1
      aa6<- sigma/(sigma+sigma60); aa6[(aa6>1.0)]<- 1
      sigma<- (ab1+aa4*aa4*aa1+aa5*aa5*aa2+aa6*aa6*aa3)/(ab2+aa4*n_samB1+aa5*n_samB2+aa6*n_samF2);
      aaa1<- abs(sigma-aaa0)
      aaa0<- sigma
      if (n_iter>20) break
    }
    sigma1E<- matrix(sigma+sigma40,d21,1);sigma2E<- matrix(sigma+sigma50,d22,1);sigma3E<- matrix(sigma+sigma60,d23,1)
    ####################### the stop criterion for iteration ############################################################
    L1<- logL(n_samP1,1,1,meanP1,sigma,dataP1)+logL(n_samF1,1,1,meanF1,sigma,dataF1)+logL(n_samP2,1,1,meanP2,sigma,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1E,sigma1E,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2E,sigma2E,dataB2)+logL(n_samF2,d23,mix_pi_3,mean3E,sigma3E,dataF2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*9

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1,1,1,1,1,1, 1,1,-1,1,1,1,-1,-1,1,1,-1,-1, 1,1,-1,1,1,-1,1,-1,1,-1,1,-1,
                1,0,-1,0.5,-0.5,-0.5,-0.5,-0.5,0,0,0,0, 0,1,0,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5),12,5)

  b_line1 <- matrix(c(meanP1,meanF1,meanP2,mean1E,mean2E,mean3E))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigmaB1 - sigma1E[1]
  if(jj_1 < 0 | jj_1>=sigmaB1) {jj_1 <- 0}
  ll_1 <- jj_1/sigmaB1
  mm1<- sigma1E[1]-sigma
  if (mm1<0 | mm1>=sigmaB1) { mm1<- 0 }
  nn1<- mm1/sigmaB1

  jj_2 <- sigmaB2 - sigma2E[1]
  if(jj_2 < 0 | jj_2>=sigmaB2) {jj_2 <- 0}
  ll_2 <- jj_2/sigmaB2
  mm2<- sigma2E[1]-sigma
  if (mm2<0 | mm2>=sigmaB2) { mm2<- 0 }
  nn2<- mm2/sigmaB2

  jj_3 <- sigmaF2 - sigma3E[1]
  if(jj_3 < 0 | jj_3>=sigmaF2) {jj_3 <- 0}
  ll_3 <- jj_3/sigmaF2
  mm3<- sigma3E[1]-sigma
  if (mm3<0 | mm3>=sigmaF2) { mm3<- 0 }
  nn3<- mm3/sigmaF2

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
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

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
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),(1-pkolm(D_F1,n_samF1))))

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
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

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
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),(1-pkolm(D_B1,n_samB1))))

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
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),(1-pkolm(D_B2,n_samB2))))

  ###############  F2  ###################
  dataF2<-sort(dataF2);bmw_F2 <- matrix(0,n_samF2,1); bmwsl_F2 <- matrix(0,n_samF2,d23)
  for(i in 1:d23){
    gg_F2 <- (dataF2 - mean3E[i])/sqrt(sigma3E[i])
    bmw_F2[which(gg_F2>=0)] <- pnorm(gg_F2[gg_F2>=0])
    bmw_F2[which(gg_F2<0)] <- 1 - pnorm(abs(gg_F2[gg_F2<0]))
    bmwsl_F2[,i] <- bmw_F2*mix_pi_3[i]
  }
  P2_F2 <- rowSums(bmwsl_F2)
  nn<-dim(as.matrix(unique(P2_F2)))[1]
  if(nn<n_samF2){P2_F2<-P2_F2+runif(n_samF2)/1e4}

  dd_F2 <- as.matrix(c(sum(P2_F2),sum(P2_F2^2),sum((P2_F2-0.5)^2)))
  WW2_F2 <- 1/(12*n_samF2) + sum((P2_F2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  u_F2 <- as.matrix(c(12*n_samF2*((dd_F2[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((dd_F2[2]/n_samF2-1/3)^2),180*n_samF2*((dd_F2[3]/n_samF2-1/12)^2)))
  D_F2 <- as.numeric(ks.test(P2_F2,"punif")[[1]][1])
  tt_F2 <- as.matrix(c((1 - pchisq(u_F2[1],1)),(1 - pchisq(u_F2[2],1)),(1 - pchisq(u_F2[3],1)),K1(WW2_F2),(1-pkolm(D_F2,n_samF2))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)
  tt_F2[which(tt_F2>=10e-4)]<-round(tt_F2[which(tt_F2>=10e-4)],4);tt_F2[which(tt_F2<10e-4)]<-format(tt_F2[which(tt_F2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame(" MX2-CD-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigma,4),round(t(mean1E),4)," "," "," ",round(sigma1E[1],4),round(t(mix_pi_1),4)," "," "," ",
                       round(t(mean2E),4),round(sigma2E[1],4),round(t(mix_pi_2),4),round(t(mean3E),4)," "," "," "," "," ",round(sigma3E[1],4),round(t(mix_pi_3),4)," "," "," "," "," ",
                       round(B1[1],4)," "," "," "," "," ",round(B1[2],4),round(B1[2],4),round(B1[3],4),round(B1[3],4)," "," "," "," ",round(B1[4],4),round(B1[5],4),round(jj_1,4),round(ll_1*100,4),round(mm1,4),round(nn1*100,4),
                       round(jj_2,4),round(ll_2*100,4),round(mm2,4),round(nn2*100,4),round(jj_3,4),round(ll_3*100,4),round(mm3,4),round(nn3*100,4),round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),
                       tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),
                       tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),
                       tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5],round(u_F2[1],4),tt_F2[1],round(u_F2[2],4),
                       tt_F2[2],round(u_F2[3],4),tt_F2[3],round(WW2_F2,4),tt_F2[4],round(D_F2,4),tt_F2[5])

  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2,mi_3)
  return(OUTPUT)
}

############################################ MX2-EAD-AD Model #########################################
G6ModelFun[[24]] <- function(K1,logL,df11,df21,df31,df41,df51,df61){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]));dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]));dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1]; n_samF1<-dim(dataF1)[1]; n_samP2<-dim(dataP2)[1]; n_samB1<-dim(dataB1)[1]; n_samB2<-dim(dataB2)[1]; n_samF2<-dim(dataF2)[1]
  sumx1<-sum(dataP1);sumx2<-sum(dataF1);sumx3<-sum(dataP2)
  sumx4<-sum(dataB1);sumx5<-sum(dataB2);sumx6<-sum(dataF2)
  sigmaP1<- as.numeric(var(dataP1)); ss1<- (n_samP1-1)*sigmaP1
  sigmaF1<- as.numeric(var(dataF1)); ss2<- (n_samF1-1)*sigmaF1
  sigmaP2<- as.numeric(var(dataP2)); ss3<- (n_samP2-1)*sigmaP2
  sigma0<- (ss1+ss2+ss3)/(n_samP1+n_samF1+n_samP2-3)
  mean4<- mean(dataB1); sigmaB1<- as.numeric(var(dataB1)); sigma40<- sigmaB1
  mean5<- mean(dataB2); sigmaB2<- as.numeric(var(dataB2)); sigma50<- sigmaB2
  mean6<- mean(dataF2); sigmaF2<- as.numeric(var(dataF2)); sigma60<- sigmaF2

  m_esp<-0.0001
  ####################### MX2-EAD-AD Model #################### (E6) ##############
  d21<-1; d22<-3; d23<-3
  sigma<- sigma0
  meanP1<- mean(dataP1);meanF1<- mean(dataF1);meanP2<- mean(dataP2)
  abb<- sigmaF2/(5*sigma0)

  mi_1<- as.matrix(1)
  sigma1E<- as.matrix(sigmaB1/abb)
  mean1E<- as.matrix(mean4)

  mi_2<- as.matrix(c(0.25,0.5,0.25))
  a2<- sqrt(sigmaB2/n_samB2)
  if (meanP1<meanP2) { a2= -a2 }
  sigma2E<- matrix(sigmaB2/abb,d22,1)
  mean2E<- as.matrix(c(mean5+2.5*a2,mean5,mean5-2.5*a2))

  mi_3<- as.matrix(c(0.5625,0.375,0.0625))
  a3<- sqrt(sigmaF2/n_samF2)
  if (meanP1<meanP2) { a3= -a3 }
  sigma3E<- matrix(sigmaF2/abb,d23,1)
  mean3E<- as.matrix(c(mean6+2.5*a3,mean6,mean6-2.5*a3))

  L0 <- logL(n_samP1,1,1,meanP1,sigma0,dataP1)+logL(n_samF1,1,1,meanF1,sigma0,dataF1)+logL(n_samP2,1,1,meanP2,sigma0,dataP2)+logL(n_samB1,d21,mi_1,mean1E,sigma1E,dataB1)+logL(n_samB2,d22,mi_2,mean2E,sigma2E,dataB2)+logL(n_samF2,d23,mi_3,mean3E,sigma3E,dataF2)

  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)
  WW_F2 <- matrix(0,d23,n_samF2); swx_F2 <- matrix(0,d23,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean1E[i],sqrt(sigma1E[i]))/dmixnorm(dataB1,mean1E,sqrt(sigma1E),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean2E[i],sqrt(sigma2E[i]))/dmixnorm(dataB2,mean2E,sqrt(sigma2E),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    n0<-as.matrix(rowSums(WW_B2));n0[abs(n0)<0.0001]<- 0.00001
    sumwx_B2 <- WW_B2%*%dataB2

    for(i in 1:d23) { WW_F2[i,] <- mi_3[i]*dnorm(dataF2,mean3E[i],sqrt(sigma3E[i]))/dmixnorm(dataF2,mean3E,sqrt(sigma3E),mi_3) }
    mix_pi_3 <- as.matrix(rowSums(WW_F2)/n_samF2)
    s0<-as.matrix(rowSums(WW_F2));s0[abs(s0)<0.0001]<- 0.00001
    sumwx_F2 <- WW_F2%*%dataF2

    ############ CM1-step for means ##############
    aaa0<- 0;aaa1<- 1000;n_iter<- 0;AA<- matrix(0,6,1)
    while (aaa1>0.0001){
      n_iter<- n_iter+1
      hh<- matrix(0,6,6)
      hh[1,1]<- sigma*(1.0/n_samP1+4.0/n_samF1+1.0/n_samP2)+sigma1E[1]/n_samB1+sigma2E[3]/n0[3]+4.0*sigma3E[1]/s0[1]
      hh[1,2]<- sigma/n_samP1-sigma/n_samP2+2.0*sigma1E[1]/n_samB1-sigma2E[3]/n0[3]
      hh[1,3]<- 0
      hh[1,4]<- sigma2E[3]/n0[3]+2.0*sigma3E[1]/s0[1]
      hh[1,5]<- -sigma1E[1]/n_samB1+2.0*sigma3E[1]/s0[1]
      hh[1,6]<- -2.0*sigma3E[1]/s0[1]
      hh[2,2]<- sigma*(1.0/n_samP1+1.0/n_samP2)+4.0*sigma1E[1]/n_samB1+sigma2E[1]/n0[1]+sigma2E[3]/n0[3]
      hh[2,3]<- sigma2E[1]/n0[1]
      hh[2,4]<- sigma2E[1]/n0[1]-sigma2E[3]/n0[3]
      hh[2,5]<- -2.0*sigma1E[1]/n_samB1
      hh[2,6]<- 0
      hh[3,3]<- sigma2E[1]/n0[1]+sigma2E[2]/n0[2]+sigma3E[2]/s0[2]+sigma3E[3]/s0[3]
      hh[3,4]<- sigma2E[1]/n0[1]+sigma3E[3]/s0[3]
      hh[3,5]<- -sigma2E[2]/n0[2]+sigma3E[2]/s0[2]
      hh[3,6]<- 2.0*sigma3E[2]/s0[2]+sigma3E[3]/s0[3]
      hh[4,4]<- sigma2E[1]/n0[1]+sigma2E[3]/n0[3]+sigma3E[1]/s0[1]+sigma3E[3]/s0[3]
      hh[4,5]<- sigma3E[1]/s0[1]
      hh[4,6]<- -sigma3E[1]/s0[1]+sigma3E[3]/s0[3]
      hh[5,5]<- sigma1E[1]/n_samB1+sigma2E[2]/n0[2]+sigma3E[1]/s0[1]+sigma3E[2]/s0[2]
      hh[5,6]<- -sigma3E[1]/s0[1]+2.0*sigma3E[2]/s0[2]
      hh[6,6]<- sigma3E[1]/s0[1]+4.0*sigma3E[2]/s0[2]+sigma3E[3]/s0[3]
      for(i in 2:6){
        for(j in 1:(i-1)){
          hh[i,j]<- hh[j,i]
        }
      }
      b_line<- matrix(0,6,1)
      b_line[1]<- sumx1/n_samP1+2.0*sumx2/n_samF1+sumx3/n_samP2-sumwx_B1[1]/n_samB1-sumwx_B2[3]/n0[3]-2.0*sumwx_F2[1]/s0[1]
      b_line[2]<- sumx1/n_samP1-sumx3/n_samP2-2.0*sumwx_B1[1]/n_samB1+sumwx_B2[1]/n0[1]+sumwx_B2[3]/n0[3]
      b_line[3]<- sumwx_B2[1]/n0[1]-sumwx_B2[2]/n0[2]-sumwx_F2[2]/s0[2]+sumwx_F2[3]/s0[3]
      b_line[4]<- sumwx_B2[1]/n0[1]-sumwx_B2[3]/n0[3]-sumwx_F2[1]/s0[1]+sumwx_F2[3]/s0[3]
      b_line[5]<- sumwx_B1[1]/n_samB1+sumwx_B2[2]/n0[2]-sumwx_F2[1]/s0[1]-sumwx_F2[2]/s0[2]
      b_line[6]<- sumwx_F2[1]/s0[1]-2.0*sumwx_F2[2]/s0[2]+sumwx_F2[3]/s0[3]

      B <- solve(hh,b_line)

      meanP1<- (sumx1-sigma*(B[1]+B[2]))/n_samP1
      meanF1<- (sumx2-sigma*2.0*B[1])/n_samF1
      meanP2<- (sumx3-sigma*(B[1]-B[2]))/n_samP2

      mean1E[1]<- (sumwx_B1[1]+sigma1E[1]*(B[1]+2.0*B[2]-B[5]))/n_samB1

      mean2E[1]<- (sumwx_B2[1]-sigma2E[1]*(B[2]+B[3]+B[4]))/n0[1]
      mean2E[2]<- (sumwx_B2[2]+sigma2E[2]*(B[3]-B[5]))/n0[2]
      mean2E[3]<- (sumwx_B2[3]+sigma2E[3]*(B[1]-B[2]+B[4]))/n0[3]

      mean3E[1]<- (sumwx_F2[1]+sigma3E[1]*(2.0*B[1]+B[4]+B[5]-B[6]))/s0[1]
      mean3E[2]<- (sumwx_F2[2]+sigma3E[2]*(B[3]+B[5]+2.0*B[6]))/s0[2]
      mean3E[3]<- (sumwx_F2[3]-sigma3E[3]*(B[3]+B[4]+B[6]))/s0[3]

      aaa1<- max(abs(AA-B))
      AA<- B
      if (n_iter>20) break
    }
    ################# iteratively CM2-step for variance ####################
    ss1<- sum((dataP1-meanP1)^2)
    ss2<- sum((dataF1-meanF1)^2)
    ss3<- sum((dataP2-meanP2)^2)

    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean1E[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean2E[i])^2 }
    for(i in 1:d23) {  swx_F2[i] <- WW_F2[i,]%*%(dataF2-mean3E[i])^2 }

    aa1<- sum(swx_B1);aa2<- sum(swx_B2);aa3<- sum(swx_F2)
    sigma1E<- matrix(aa1/n_samB1,d21,1)
    sigma40<- sigma1E[1]-sigma
    sigma2E<- matrix(aa2/n_samB2,d22,1)
    sigma50<- sigma2E[1]-sigma
    sigma3E<- matrix(aa3/n_samF2,d23,1)
    sigma60<- sigma3E[1]-sigma
    ################ iteratively CM3-step for variance (sigma) ################
    sigma40[(sigma40<0)]<- 0.00001; sigma50[(sigma50<0)]<- 0.00001; sigma60[(sigma60<0)]<- 0.00001
    ab1<- ss1+ss2+ss3; ab2<- n_samP1+n_samF1+n_samP2
    aaa0<- sigma ;n_iter<- 0;aaa1<- 1000
    while (aaa1>0.0001){
      n_iter<- n_iter+1
      aa4<- sigma/(sigma+sigma40); aa4[(aa4>1.0)]<- 1
      aa5<- sigma/(sigma+sigma50); aa5[(aa5>1.0)]<- 1
      aa6<- sigma/(sigma+sigma60); aa6[(aa6>1.0)]<- 1
      sigma<- (ab1+aa4*aa4*aa1+aa5*aa5*aa2+aa6*aa6*aa3)/(ab2+aa4*n_samB1+aa5*n_samB2+aa6*n_samF2);
      aaa1<- abs(sigma-aaa0)
      aaa0<- sigma
      if (n_iter>20) break
    }
    sigma1E<- matrix(sigma+sigma40,d21,1);sigma2E<- matrix(sigma+sigma50,d22,1);sigma3E<- matrix(sigma+sigma60,d23,1)
    ####################### the stop criterion for iteration ############################################################
    L1<- logL(n_samP1,1,1,meanP1,sigma,dataP1)+logL(n_samF1,1,1,meanF1,sigma,dataF1)+logL(n_samP2,1,1,meanP2,sigma,dataP2)+logL(n_samB1,d21,mix_pi_1,mean1E,sigma1E,dataB1)+logL(n_samB2,d22,mix_pi_2,mean2E,sigma2E,dataB2)+logL(n_samF2,d23,mix_pi_3,mean3E,sigma3E,dataF2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc <- L1
  AIC <- -2*abc + 2*8

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,1,1,1,1,1,1, 2,2,-2,2,2,0,-2,2,0,-2, 1,0,-1,0.5,-0.5,-0.5,-0.5,0,0,0, 0,1,0,0.5,0.5,0.5,0.5,0.5,0.5,0.5),10,4)

  b_line1 <- matrix(c(meanP1,meanF1,meanP2,mean1E,mean2E,mean3E))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigmaB1 - sigma1E[1]
  if(jj_1 < 0 | jj_1>=sigmaB1) {jj_1 <- 0}
  ll_1 <- jj_1/sigmaB1
  mm1<- sigma1E[1]-sigma
  if (mm1<0 | mm1>=sigmaB1) { mm1<- 0 }
  nn1<- mm1/sigmaB1

  jj_2 <- sigmaB2 - sigma2E[1]
  if(jj_2 < 0 | jj_2>=sigmaB2) {jj_2 <- 0}
  ll_2 <- jj_2/sigmaB2
  mm2<- sigma2E[1]-sigma
  if (mm2<0 | mm2>=sigmaB2) { mm2<- 0 }
  nn2<- mm2/sigmaB2

  jj_3 <- sigmaF2 - sigma3E[1]
  if(jj_3 < 0 | jj_3>=sigmaF2) {jj_3 <- 0}
  ll_3 <- jj_3/sigmaF2
  mm3<- sigma3E[1]-sigma
  if (mm3<0 | mm3>=sigmaF2) { mm3<- 0 }
  nn3<- mm3/sigmaF2

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
  D_P1 <- as.numeric(ks.test(P2_P1,"punif")[[1]][1])
  tt_P1 <- as.matrix(c((1 - pchisq(u_P1[1],1)),(1 - pchisq(u_P1[2],1)),(1 - pchisq(u_P1[3],1)),K1(WW2_P1),(1-pkolm(D_P1,n_samP1))))

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
  D_F1 <- as.numeric(ks.test(P2_F1,"punif")[[1]][1])
  tt_F1 <- as.matrix(c((1 - pchisq(u_F1[1],1)),(1 - pchisq(u_F1[2],1)),(1 - pchisq(u_F1[3],1)),K1(WW2_F1),(1-pkolm(D_F1,n_samF1))))

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
  D_P2 <- as.numeric(ks.test(P2_P2,"punif")[[1]][1])
  tt_P2 <- as.matrix(c((1 - pchisq(u_P2[1],1)),(1 - pchisq(u_P2[2],1)),(1 - pchisq(u_P2[3],1)),K1(WW2_P2),(1-pkolm(D_P2,n_samP2))))

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
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),(1-pkolm(D_B1,n_samB1))))

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
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),(1-pkolm(D_B2,n_samB2))))

  ###############  F2  ###################
  dataF2<-sort(dataF2);bmw_F2 <- matrix(0,n_samF2,1); bmwsl_F2 <- matrix(0,n_samF2,d23)
  for(i in 1:d23){
    gg_F2 <- (dataF2 - mean3E[i])/sqrt(sigma3E[i])
    bmw_F2[which(gg_F2>=0)] <- pnorm(gg_F2[gg_F2>=0])
    bmw_F2[which(gg_F2<0)] <- 1 - pnorm(abs(gg_F2[gg_F2<0]))
    bmwsl_F2[,i] <- bmw_F2*mix_pi_3[i]
  }
  P2_F2 <- rowSums(bmwsl_F2)
  nn<-dim(as.matrix(unique(P2_F2)))[1]
  if(nn<n_samF2){P2_F2<-P2_F2+runif(n_samF2)/1e4}

  dd_F2 <- as.matrix(c(sum(P2_F2),sum(P2_F2^2),sum((P2_F2-0.5)^2)))
  WW2_F2 <- 1/(12*n_samF2) + sum((P2_F2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  u_F2 <- as.matrix(c(12*n_samF2*((dd_F2[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((dd_F2[2]/n_samF2-1/3)^2),180*n_samF2*((dd_F2[3]/n_samF2-1/12)^2)))
  D_F2 <- as.numeric(ks.test(P2_F2,"punif")[[1]][1])
  tt_F2 <- as.matrix(c((1 - pchisq(u_F2[1],1)),(1 - pchisq(u_F2[2],1)),(1 - pchisq(u_F2[3],1)),K1(WW2_F2),(1-pkolm(D_F2,n_samF2))))

  tt_P1[which(tt_P1>=10e-4)]<-round(tt_P1[which(tt_P1>=10e-4)],4);tt_P1[which(tt_P1<10e-4)]<-format(tt_P1[which(tt_P1<10e-4)],scientific=TRUE,digit=4)
  tt_F1[which(tt_F1>=10e-4)]<-round(tt_F1[which(tt_F1>=10e-4)],4);tt_F1[which(tt_F1<10e-4)]<-format(tt_F1[which(tt_F1<10e-4)],scientific=TRUE,digit=4)
  tt_P2[which(tt_P2>=10e-4)]<-round(tt_P2[which(tt_P2>=10e-4)],4);tt_P2[which(tt_P2<10e-4)]<-format(tt_P2[which(tt_P2<10e-4)],scientific=TRUE,digit=4)
  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)
  tt_F2[which(tt_F2>=10e-4)]<-round(tt_F2[which(tt_F2>=10e-4)],4);tt_F2[which(tt_F2<10e-4)]<-format(tt_F2[which(tt_F2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame(" MX2-EAD-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigma,4),round(t(mean1E),4)," "," "," ",round(sigma1E[1],4),round(t(mix_pi_1),4)," "," "," ",
                       round(t(mean2E),4)," ",round(sigma2E[1],4),round(t(mix_pi_2),4)," ",round(t(mean3E),4)," "," "," "," "," "," ",round(sigma3E[1],4),round(t(mix_pi_3),4)," "," "," "," "," "," ",
                       round(B1[1],4)," "," "," "," "," ",round(B1[2],4),round(B1[2],4),round(B1[2],4),round(B1[2],4)," "," "," "," ",round(B1[3],4),round(B1[4],4),round(jj_1,4),round(ll_1*100,4),round(mm1,4),round(nn1*100,4),
                       round(jj_2,4),round(ll_2*100,4),round(mm2,4),round(nn2*100,4),round(jj_3,4),round(ll_3*100,4),round(mm3,4),round(nn3*100,4),round(u_P1[1],4),tt_P1[1],round(u_P1[2],4),
                       tt_P1[2],round(u_P1[3],4),tt_P1[3],round(WW2_P1,4),tt_P1[4],round(D_P1,4),tt_P1[5],round(u_F1[1],4),tt_F1[1],round(u_F1[2],4),
                       tt_F1[2],round(u_F1[3],4),tt_F1[3],round(WW2_F1,4),tt_F1[4],round(D_F1,4),tt_F1[5],round(u_P2[1],4),tt_P2[1],round(u_P2[2],4),
                       tt_P2[2],round(u_P2[3],4),tt_P2[3],round(WW2_P2,4),tt_P2[4],round(D_P2,4),tt_P2[5],round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),
                       tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),
                       tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5],round(u_F2[1],4),tt_F2[1],round(u_F2[2],4),
                       tt_F2[2],round(u_F2[3],4),tt_F2[3],round(WW2_F2,4),tt_F2[4],round(D_F2,4),tt_F2[5])



  output<-as.matrix(output)

  OUTPUT<-list(output,mi_1,mi_2,mi_3)
  return(OUTPUT)
}



K1G6 <- function(x){
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

logLG6 <- function(nm,nng,mi,mn,s,d1) { sum2 <- sum(log(dmixnorm(d1,mn,sqrt(s),mi)));return (sum2) }


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
  allresult=foreach(i=1:24,.combine = 'rbind')%dopar%{
    requireNamespace("KScorrect")
    requireNamespace("kolmim")
    G6ModelFun[[i]](K1G6,logLG6,df11,df21,df31,df41,df51,df61)[[1]]
  }
  stopCluster(cl)
  mi_1<-NULL;mi_2<-NULL;mi_3<-NULL

}else{

 allresultq=switch(model,"1MG-AD" = G6ModelFun[[1]](K1G6,logLG6,df11,df21,df31,df41,df51,df61),"1MG-A"=G6ModelFun[[2]](K1G6,logLG6,df11,df21,df31,df41,df51,df61),"1MG-EAD"=G6ModelFun[[3]](K1G6,logLG6,df11,df21,df31,df41,df51,df61),"1MG-NCD"=G6ModelFun[[4]](K1G6,logLG6,df11,df21,df31,df41,df51,df61),"2MG-ADI"=G6ModelFun[[5]](K1G6,logLG6,df11,df21,df31,df41,df51,df61),
                   "2MG-AD"=G6ModelFun[[6]](K1G6,logLG6,df11,df21,df31,df41,df51,df61),"2MG-A"=G6ModelFun[[7]](K1G6,logLG6,df11,df21,df31,df41,df51,df61),"2MG-EA"=G6ModelFun[[8]](K1G6,logLG6,df11,df21,df31,df41,df51,df61),"2MG-CD"=G6ModelFun[[9]](K1G6,logLG6,df11,df21,df31,df41,df51,df61),"2MG-EAD"=G6ModelFun[[10]](K1G6,logLG6,df11,df21,df31,df41,df51,df61),
                   "PG-ADI"=G6ModelFun[[11]](K1G6,logLG6,df11,df21,df31,df41,df51,df61),"PG-AD"=G6ModelFun[[12]](K1G6,logLG6,df11,df21,df31,df41,df51,df61),"MX1-AD-ADI"=G6ModelFun[[13]](K1G6,logLG6,df11,df21,df31,df41,df51,df61),"MX1-AD-AD"=G6ModelFun[[14]](K1G6,logLG6,df11,df21,df31,df41,df51,df61),"MX1-A-AD"=G6ModelFun[[15]](K1G6,logLG6,df11,df21,df31,df41,df51,df61),
                   "MX1-EAD-AD"=G6ModelFun[[16]](K1G6,logLG6,df11,df21,df31,df41,df51,df61),"MX1-NCD-AD"=G6ModelFun[[17]](K1G6,logLG6,df11,df21,df31,df41,df51,df61),"MX2-ADI-ADI"=G6ModelFun[[18]](K1G6,logLG6,df11,df21,df31,df41,df51,df61),"MX2-ADI-AD"=G6ModelFun[[19]](K1G6,logLG6,df11,df21,df31,df41,df51,df61),"MX2-AD-AD"=G6ModelFun[[20]](K1G6,logLG6,df11,df21,df31,df41,df51,df61),
                   "MX2-A-AD"=G6ModelFun[[21]](K1G6,logLG6,df11,df21,df31,df41,df51,df61),"MX2-EA-AD"=G6ModelFun[[22]](K1G6,logLG6,df11,df21,df31,df41,df51,df61),"MX2-CD-AD"=G6ModelFun[[23]](K1G6,logLG6,df11,df21,df31,df41,df51,df61),"MX2-EAD-AD"=G6ModelFun[[24]](K1G6,logLG6,df11,df21,df31,df41,df51,df61))

 allresult<-allresultq[[1]]
 if(model=="PG-AD"||model=="PG-ADI"){
   mi_1<-NULL;mi_2<-NULL;mi_3<-NULL
 }else{
   mi_1<-allresultq[[2]];mi_2<-allresultq[[3]];mi_3<-allresultq[[4]]
 }
}
colnames(allresult) <- G6colname
out<-list(allresult,mi_1,mi_2,mi_3)
return(out)
}




