G6FFun<-function(df,model,G6Ftext2){

data<-sapply(df,as.character)
dP1<-data[-1,which(data[1,]=="P1")];P1<-as.numeric(dP1[which(is.na(as.numeric(dP1))==FALSE)]);df11<-as.data.frame(P1)
dF1<-data[-1,which(data[1,]=="F1")];F1<-as.numeric(dF1[which(is.na(as.numeric(dF1))==FALSE)]);df21<-as.data.frame(F1)
dP2<-data[-1,which(data[1,]=="P2")];P2<-as.numeric(dP2[which(is.na(as.numeric(dP2))==FALSE)]);df31<-as.data.frame(P2)
dB12<-data[-1,which(data[1,]=="B12")];B12<-as.numeric(dB12[which(is.na(as.numeric(dB12))==FALSE)]);df41<-as.data.frame(B12)
dB22<-data[-1,which(data[1,]=="B22")];B22<-as.numeric(dB22[which(is.na(as.numeric(dB22))==FALSE)]);df51<-as.data.frame(B22)
dF23<-data[-1,which(data[1,]=="F23")];F23<-as.numeric(dF23[which(is.na(as.numeric(dF23))==FALSE)]);df61<-as.data.frame(F23)

G6Fcolname<-c("Model","Log_Max_likelihood_value","AIC","mean[P1]","mean[F1]","mean[P2]","Var(P1 & P2 & F1)","B1:2-mean[1]","B1:2-mean[2]","B1:2-mean[3]","B1:2-mean[4]",
              "B1:2-Var[1]","B1:2-Var[2]","B1:2-Var[3]","B1:2-Var[4]","B1:2-Proportion[1]","B1:2-Proportion[2]","B1:2-Proportion[3]","B1:2-Proportion[4]",
              "B2:2-mean[1]","B2:2-mean[2]","B2:2-mean[3]","B2:2-mean[4]","B2:2-Var[1]","B2:2-Var[2]","B2:2-Var[3]","B2:2-Var[4]",
              "B2:2-Proportion[1]","B2:2-Proportion[2]","B2:2-Proportion[3]","B2:2-Proportion[4]",
              "F2:3-mean[1]","F2:3-mean[2]","F2:3-mean[3]","F2:3-mean[4]","F2:3-mean[5]","F2:3-mean[6]","F2:3-mean[7]","F2:3-mean[8]","F2:3-mean[9]",
              "F2:3-Var[1]","F2:3-Var[2]","F2:3-Var[3]","F2:3-Var[4]","F2:3-Var[5]","F2:3-Var[6]","F2:3-Var[7]","F2:3-Var[8]","F2:3-Var[9]",
              "F2:3-Proportion[1]","F2:3-Proportion[2]","F2:3-Proportion[3]","F2:3-Proportion[4]","F2:3-Proportion[5]","F2:3-Proportion[6]","F2:3-Proportion[7]","F2:3-Proportion[8]","F2:3-Proportion[9]",
              "m1(m)","m2","m3","m4","m5","m6","da","db","ha","hb","i","jab","jba","l","[d]","[h]",
              "B1:2-Major-Gene Var","B1:2-Heritability(Major-Gene)(%)","B1:2-Polygenes Var","B1:2-Heritability(Polygenes)(%)",
              "B2:2-Major-Gene Var","B2:2-Heritability(Major-Gene)(%)","B2:2-Polygenes Var","B2:2-Heritability(Polygenes)(%)",
              "F2:3-Major-Gene Var","F2:3-Heritability(Major-Gene)(%)","F2:3-Polygenes Var","F2:3-Heritability(Polygenes)(%)",
              "U1 square(P1)","P(U1 square(P1))","U2 square(P1)","P(U2 square(P1))","U3 square(P1)","P(U3 square(P1))","nW square(P1)","P(nW square(P1))","Dn(P1)","P(Dn(P1))",
              "U1 square(F1)","P(U1 square(F1))","U2 square(F1)","P(U2 square(F1))","U3 square(F1)","P(U3 square(F1))","nW square(F1)","P(nW square(F1))","Dn(F1)","P(Dn(F1))",
              "U1 square(P2)","P(U1 square(P2))","U2 square(P2)","P(U2 square(P2))","U3 square(P2)","P(U3 square(P2))","nW square(P2)","P(nW square(P2))","Dn(P2)","P(Dn(P2))",
              "U1 square(B1:2)","P(U1 square(B1:2))","U2 square(B1:2)","P(U2 square(B1:2))","U3 square(B1:2)","P(U3 square(B1:2))","nW square(B1:2)","P(nW square(B1:2))","Dn(B1:2)","P(Dn(B1:2))",
              "U1 square(B2:2)","P(U1 square(B2:2))","U2 square(B2:2)","P(U2 square(B2:2))","U3 square(B2:2)","P(U3 square(B2:2))","nW square(B2:2)","P(nW square(B2:2))","Dn(B2:2)","P(Dn(B2:2))",
              "U1 square(F2:3)","P(U1 square(F2:3))","U2 square(F2:3)","P(U2 square(F2:3))","U3 square(F2:3)","P(U3 square(F2:3))","nW square(F2:3)","P(nW square(F2:3))","Dn(F2:3)","P(Dn(F2:3))")
G6FModelFun<-list(NA)
###################define each model function############################
##########################1MG_AD(A_1)#################################
G6FModelFun[[1]] <- function(K1,logL,df11,df21,df31,df41,df51,df61,G6Ftext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1];n_samF1<-dim(dataF1)[1];n_samP2<-dim(dataP2)[1]
  n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1];n_samF2<-dim(dataF2)[1]
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2),sum(dataF2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2),sum(dataF2^2)))
  ss<-matrix(0,3,1);ss[1]<-s[1]-sumx[1]^2/n_samP1;ss[2]<-s[2]-sumx[2]^2/n_samF1;ss[3]<-s[3]-sumx[3]^2/n_samP2
  mean<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataB1),mean(dataB2),mean(dataF2)))
  sigma0<-sum(ss)/(n_samP1+n_samF1+n_samP2-3)
  sigmaB1<-var(dataB1);sigma40<-sigmaB1
  sigmaB2<-var(dataB2);sigma50<-sigmaB2
  sigmaF2<-var(dataF2);sigma60<-sigmaF2
  sigmaP1<-sigmaF1<-sigmaP2<-sigma0
  m_esp <- 0.0001;n_fam <- as.numeric(G6Ftext2)
  ###############procedure start###########################
  mi1<-matrix(0.5,2,1);sigma1<-matrix(0,2,1)
  mi2<-matrix(0.5,2,1);sigma2<-matrix(0,2,1)
  mi3<-as.matrix(c(0.25,0.5,0.25));sigma3<-matrix(0,3,1)
  sigma<-sigma0;sigma1[1]<-sigma
  a1<-sqrt(sigmaB1/n_samB1)
  if (mean[1]<mean[3]) a1<--a1
  mean1<-as.matrix(c(mean[4]+2*a1,mean[4]-2*a1))
  sigma2[2]<-sigma
  a2<-sqrt(sigmaB2/n_samB2)
  if (mean[1]<mean[3]) a2<--a2
  mean2<-as.matrix(c(mean[5]+2*a2,mean[5]-2*a2))
  sigma3[1]<-sigma3[3]<-sigma
  a3<-sqrt(sigmaF2/n_samF2)
  if (mean[1]<mean[3]) a3<--a3
  mean3<-as.matrix(c(mean[6]+2*a3,mean[6],mean[6]-2*a3))
  b1<-0.5*(mean[1]-mean[3])  #additive effect.
  b2<-(-6*mean[1]+10*mean[2]-6*mean[3]+2*mean1[2])/11   #dominance effect.
  sigma1[2]<-sigma+(0.5*b1^2+0.25*b2^2)/n_fam
  sigma2[1]<-sigma1[2];sigma3[2]<-sigma1[2]
  L0<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma0))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma0))))+
    sum(log(dnorm(dataP2,mean[3],sqrt(sigma0))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mi1)))+
    sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mi3)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,2,n_samB1);      swx1 <- matrix(0,2,1)
  W2 <- matrix(0,2,n_samB2);      swx2 <- matrix(0,2,1)
  W3 <- matrix(0,3,n_samF2);      swx3 <- matrix(0,3,1)
  s0<-matrix(0,6,1);              n0<-matrix(0,6,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:2) { W1[i,] <- mi1[i]*dnorm(dataB1,mean1[i],sqrt(sigma1[i]))/dmixnorm(dataB1,mean1,sqrt(sigma1),mi1) }
    mix_pi1 <- as.matrix(rowSums(W1)/n_samB1)
    sumwx1 <- W1%*%dataB1
    for(i in 1:2) { W2[i,] <- mi2[i]*dnorm(dataB2,mean2[i],sqrt(sigma2[i]))/dmixnorm(dataB2,mean2,sqrt(sigma2),mi2) }
    mix_pi2 <- as.matrix(rowSums(W2)/n_samB2)
    sumwx2 <- W2%*%dataB2
    for(i in 1:3) { W3[i,] <- mi3[i]*dnorm(dataF2,mean3[i],sqrt(sigma3[i]))/dmixnorm(dataF2,mean3,sqrt(sigma3),mi3) }
    mix_pi3 <- as.matrix(rowSums(W3)/n_samF2)
    sumwx3 <- W3%*%dataF2
    #########obtain means#################################
    aaa0<-0
    s0[1]<-sumx[1]+sumwx1[1]+sumwx3[1];s0[2]<-sumx[2]
    s0[3]<-sumx[3]+sumwx2[2]+sumwx3[3];s0[4]<-sumwx1[2]+sumwx2[1]+sumwx3[2]
    n0[1]<-n_samP1+mix_pi1[1]*n_samB1+mix_pi3[1]*n_samF2;n0[2]<-n_samF1
    n0[3]<-n_samP2+mix_pi2[2]*n_samB2+mix_pi3[3]*n_samF2;n0[4]<-mix_pi1[2]*n_samB1+mix_pi2[1]*n_samB2+mix_pi3[2]*n_samF2
    n0[c(1:4)][abs(n0[c(1:4)])<0.00000001]<-0.000001
    aa3<-s0[1]/n0[1]+s0[3]/n0[3]+2*s0[2]/n0[2]-4*s0[4]/n0[4];aa4<-sigma*(1/n0[1]+1/n0[3]+4/n0[2])
    aa1<-1000;n_iter<-0
    while (aa1>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-0.5*(mean[1]-mean[3])
      aa2<-(-6*mean[1]+10*mean[2]-6*mean[3]+2*mean1[2])/11
      sigma1[2]<-sigma+(0.5*aa1^2+0.25*aa2^2)/n_fam
      aa2<-aa4+16*sigma1[2]/n0[4]
      aaa1<-aa3/aa2        # coefficient in restricted condition.
      mean[1]<-(s0[1]-aaa1*sigma)/n0[1]
      mean[2]<-(s0[2]-2*aaa1*sigma)/n0[2]
      mean[3]<-(s0[3]-aaa1*sigma)/n0[3]
      mean1[2]<-(s0[4]+4*aaa1*sigma1[2])/n0[4]
      aa1<-abs(aaa1-aaa0)
      aaa0<-aaa1
      if(n_iter>20)break
    }
    mean1[1]<-mean3[1]<-mean[1];mean2[1]<-mean3[2]<-mean1[2];mean2[2]<-mean3[3]<-mean[3]
    ######################obtain variance#########################
    ss1<-sum((dataP1-mean[1])^2);ss3<-sum((dataP2-mean[3])^2);ss2<-sum((dataF1-mean[2])^2)
    for(i in 1:2) {swx1[i] <- W1[i,]%*%(dataB1-mean1[i])^2 };for(i in 1:2) {swx2[i] <- W2[i,]%*%(dataB2-mean2[i])^2 };for(i in 1:3) {swx3[i] <- W3[i,]%*%(dataF2-mean3[i])^2 }
    s0[5]<-ss1+ss2+ss3+swx1[1]+swx2[2]+swx3[1]+swx3[3]
    n0[5]<-n_samP1+n_samF1+n_samP2+mix_pi1[1]*n_samB1+mix_pi2[2]*n_samB2+(mix_pi3[1]+mix_pi3[3])*n_samF2
    s0[6]<-swx1[2]+swx2[1]+swx3[2];n0[6]<-mix_pi1[2]*n_samB1+mix_pi2[1]*n_samB2+mix_pi3[2]*n_samF2
    aaa0<-sigma
    a<-0.5*(mean[1]-mean[3])
    aa2<-(-6*mean[1]+10*mean[2]-6*mean[3]+2*mean1[2])/11
    a<-(0.5*a*a+0.25*aa2*aa2)/n_fam
    aa3<-1000;n_iter<-0
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      aa2<-sigma/(sigma+a)
      sigma<-(s0[5]+aa2*aa2*s0[6])/(n0[5]+aa2*n0[6])
      aa3<-abs(sigma-aaa0)
      aaa0<-sigma
      if(n_iter>20)break
    }
    sigma1[2]<-sigma+a;sigma2[1]<-sigma1[2];sigma3[2]<-sigma1[2]
    sigma1[1]<-sigma2[2]<-sigma3[1]<-sigma3[3]<-sigma
    if(sum(sigma < 1e-30)>=1){break}
    ##############################criteria for iterations to stop########################
    L1<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,mean[3],sqrt(sigma))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mix_pi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mix_pi3)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc<-L0
  AIC<--2*abc+2*4
  meanP1<-mean[1];meanF1<-mean[2];meanP2<-mean[3]
  sigma0<-sigma
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*n_samP1)
  P1bmw <- matrix(0,n_samP1,1)
  P1gg <- (dataP1 - mean[1])/sqrt(as.vector(sigma))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < n_samP1){P1bmw <- P1bmw+runif(n_samP1)/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  P1u<- as.matrix(c(12*n_samP1*((P1dd[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((P1dd[2]/n_samP1-1/3)^2),180*n_samP1*((P1dd[3]/n_samP1-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[2])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),P1D))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*n_samF1)
  F1bmw <- matrix(0,n_samF1,1)
  F1gg <- (dataF1 - mean[2])/sqrt(as.vector(sigma))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < n_samF1){F1bmw <- F1bmw+runif(n_samF1)/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  F1u<- as.matrix(c(12*n_samF1*((F1dd[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((F1dd[2]/n_samF1-1/3)^2),180*n_samF1*((F1dd[3]/n_samF1-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[2])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),F1D))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)
  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*n_samP2)
  P2bmw <- matrix(0,n_samP2,1)
  P2gg <- (dataP2 - mean[3])/sqrt(as.vector(sigma))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < n_samP2){P2bmw <- P2bmw+runif(n_samP2)/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  P2u<- as.matrix(c(12*n_samP2*((P2dd[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((P2dd[2]/n_samP2-1/3)^2),180*n_samP2*((P2dd[3]/n_samP2-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[2])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),P2D))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1);
  B1w1<-1/(12*n_samB1)
  B1bmw <- matrix(0,n_samB1,1); B1bmwsl <- matrix(0,n_samB1,2)
  for(i in 1:2){
    B1gg <- (dataB1 - mean1[i])/sqrt(sigma1[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi1[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B1P2)))[1]
  if(nn < n_samB1){B1P2 <- B1P2+runif(n_samB1)/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*n_samB1) + sum((B1P2 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  B1u <- as.matrix(c(12*n_samB1*((B1dd[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((B1dd[2]/n_samB1-1/3)^2),180*n_samB1*((B1dd[3]/n_samB1-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[2])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),B1D))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])

  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2);
  B2w1<-1/(12*n_samB2)
  B2bmw <- matrix(0,n_samB2,1); B2bmwsl <- matrix(0,n_samB2,2)
  for(i in 1:2){
    B2gg <- (dataB2 - mean2[i])/sqrt(sigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B2P2)))[1]
  if(nn < n_samB2){B2P2 <- B2P2+runif(n_samB2)/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*n_samB2) + sum((B2P2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  B2u <- as.matrix(c(12*n_samB2*((B2dd[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((B2dd[2]/n_samB2-1/3)^2),180*n_samB2*((B2dd[3]/n_samB2-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[2])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),B2D))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])

  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2);
  F2w1<-1/(12*n_samF2)
  F2bmw <- matrix(0,n_samF2,1); F2bmwsl <- matrix(0,n_samF2,3)
  for(i in 1:3){
    F2gg <- (dataF2 - mean3[i])/sqrt(sigma3[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi3[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < n_samF2){F2P2 <- F2P2+runif(n_samF2)/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*n_samF2) + sum((F2P2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  F2u <- as.matrix(c(12*n_samF2*((F2dd[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((F2dd[2]/n_samF2-1/3)^2),180*n_samF2*((F2dd[3]/n_samF2-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[2])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),F2D))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])

  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)

  ############first order genetic parameters#########################
  aa<-matrix(c(1,1,0,1,0,1,1,-1,0,1,0,0.5),4,3,byrow=T)
  mm<-as.matrix(c(mean[1],mean[2],mean[3],mean1[2]))
  B<-solve(crossprod(aa,aa))%*%crossprod(aa,mm)
  ############second order genetic parameters##########################
  jj1<-sigmaB1-sigma1[1]
  if (jj1<0 ) {jj1<-0}
  ll1<-jj1/sigmaB1
  jj2<-sigmaB2-sigma2[2]
  if (jj2<0) {jj2<-0}
  ll2<-jj2/sigmaB2
  jj3<-sigmaF2-sigma3[1]
  if (jj3<0) {jj3<-0}
  ll3<-jj3/sigmaF2

  output <- data.frame("1MG-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4), round(sigma0,4),round(t(mean1),4)," "," ",round(t(sigma1),4)," "," ",
                       round(t(mix_pi1),4)," "," ",round(t(mean2),4)," "," ",round(t(sigma2),4)," "," ",round(t(mix_pi2),4)," "," ",
                       round(t(mean3),4)," "," "," "," "," "," ",round(t(sigma3),4)," "," "," "," "," "," ",round(t(mix_pi3),4)," "," "," "," "," "," ",
                       round(B[1],4)," "," "," "," "," ",round(B[2],4)," ",round(B[3],4)," "," "," "," "," "," "," ",
                       round(jj1,4),round(ll1*100,4)," "," ",round(jj2,4),round(ll2*100,4)," "," ",round(jj3,4),round(ll3*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi1,mi2,mi3)
  return(OUTPUT)
}

###########1MG-A(A-2)###############################
G6FModelFun[[2]] <- function(K1,logL,df11,df21,df31,df41,df51,df61,G6Ftext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1];n_samF1<-dim(dataF1)[1];n_samP2<-dim(dataP2)[1]
  n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1];n_samF2<-dim(dataF2)[1]
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2),sum(dataF2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2),sum(dataF2^2)))
  ss<-matrix(0,3,1);ss[1]<-s[1]-sumx[1]^2/n_samP1;ss[2]<-s[2]-sumx[2]^2/n_samF1;ss[3]<-s[3]-sumx[3]^2/n_samP2
  mean<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataB1),mean(dataB2),mean(dataF2)))
  sigma0<-sum(ss)/(n_samP1+n_samF1+n_samP2-3)
  sigmaB1<-var(dataB1);sigma40<-sigmaB1
  sigmaB2<-var(dataB2);sigma50<-sigmaB2
  sigmaF2<-var(dataF2);sigma60<-sigmaF2
  sigmaP1<-sigmaF1<-sigmaP2<-sigma0
  m_esp <- 0.0001;n_fam <- as.numeric(G6Ftext2)
  ###############procedure start###########################
  mi1<-matrix(0.5,2,1);sigma1<-matrix(0,2,1)
  mi2<-matrix(0.5,2,1);sigma2<-matrix(0,2,1)
  mi3<-as.matrix(c(0.25,0.5,0.25));sigma3<-matrix(0,3,1)
  sigma<-sigma0;sigma1[1]<-sigma
  a1<-sqrt(sigmaB1/n_samB1)
  if (mean[1]<mean[3]) a1<--a1
  mean1<-as.matrix(c(mean[4]+2*a1,mean[4]-2*a1))
  sigma2[2]<-sigma
  a2<-sqrt(sigmaB2/n_samB2)
  if (mean[1]<mean[3]) a2<--a2
  mean2<-as.matrix(c(mean[5]+2*a2,mean[5]-2*a2))
  sigma3[1]<-sigma;sigma3[3]<-sigma
  a3<-sqrt(sigmaF2/n_samF2)
  if (mean[1]<mean[3]) a3<--a3
  mean3<-as.matrix(c(mean[6]+2*a3,mean[6],mean[6]-2*a3))
  b1<-0.5*(mean[1]-mean[3])  #additive effect.
  sigma1[2]<-sigma+(0.5*b1^2)/n_fam;sigma2[1]<-sigma1[2];sigma3[2]<-sigma1[2]
  L0<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma0))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma0))))+
    sum(log(dnorm(dataP2,mean[3],sqrt(sigma0))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mi1)))+
    sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mi3)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,2,n_samB1);      swx1 <- matrix(0,2,1)
  W2 <- matrix(0,2,n_samB2);      swx2 <- matrix(0,2,1)
  W3 <- matrix(0,3,n_samF2);      swx3 <- matrix(0,3,1)
  n0<-matrix(0,6,1);s0<-matrix(0,6,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:2) { W1[i,] <- mi1[i]*dnorm(dataB1,mean1[i],sqrt(sigma1[i]))/dmixnorm(dataB1,mean1,sqrt(sigma1),mi1) }
    mix_pi1 <- as.matrix(rowSums(W1)/n_samB1)
    sumwx1 <- W1%*%dataB1
    for(i in 1:2) { W2[i,] <- mi2[i]*dnorm(dataB2,mean2[i],sqrt(sigma2[i]))/dmixnorm(dataB2,mean2,sqrt(sigma2),mi2) }
    mix_pi2 <- as.matrix(rowSums(W2)/n_samB2)
    sumwx2 <- W2%*%dataB2
    for(i in 1:3) { W3[i,] <- mi3[i]*dnorm(dataF2,mean3[i],sqrt(sigma3[i]))/dmixnorm(dataF2,mean3,sqrt(sigma3),mi3) }
    mix_pi3 <- as.matrix(rowSums(W3)/n_samF2)
    sumwx3 <- W3%*%dataF2
    #########obtain means############################
    aaa0<-0
    s0[1]<-sumx[1]+sumwx1[1]+sumwx3[1];s0[2]<-sumwx1[2]+sumwx2[1]+sumwx3[2]
    s0[3]<-sumx[3]+sumwx2[2]+sumwx3[3];n0[1]<-n_samP1+mix_pi1[1]*n_samB1+mix_pi3[1]*n_samF2
    n0[2]<-mix_pi1[2]*n_samB1+mix_pi2[1]*n_samB2+mix_pi3[2]*n_samF2;n0[3]<-n_samP2+mix_pi2[2]*n_samB2+mix_pi3[3]*n_samF2
    n0[c(1:3)][abs(n0[c(1:3)])<0.00000001]<-0.000001
    aa2<-1000;n_iter<-0
    while (aa2>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-0.5*(mean[1]-mean[3])
      sigma1[2]<-sigma+0.5*aa1^2/n_fam
      s0[2]<-sumx[2]+s0[2]*sigma/sigma1[2]
      n0[2]<-n_samF1+n0[2]*sigma/sigma1[2]
      aa3<-s0[1]/n0[1]-2*s0[2]/n0[2]+s0[3]/n0[3]
      aa4<-sigma*(1/n0[1]+4/n0[2]+1/n0[3])
      aaa1<-aa3/aa4
      mean[1]<-(s0[1]-aaa1*sigma)/n0[1]
      mean[2]<-(s0[2]+2*aaa1*sigma)/n0[2]
      mean[3]<-(s0[3]-aaa1*sigma)/n0[3]
      aa2<-abs(aaa1-aaa0)
      aaa0<-aaa1
      if(n_iter>20)break
    }
    mean1[1]<-mean3[1]<-mean[1];mean1[2]<-mean2[1]<-mean3[2]<-mean[2];mean2[2]<-mean3[3]<-mean[3]
    ###########obtain variance##############################
    ss1<-sum((dataP1-mean[1])^2);ss3<-sum((dataP2-mean[3])^2);ss2<-sum((dataF1-mean[2])^2)
    for(i in 1:2) {swx1[i] <- W1[i,]%*%(dataB1-mean1[i])^2 };for(i in 1:2) {swx2[i] <- W2[i,]%*%(dataB2-mean2[i])^2 };for(i in 1:3) {swx3[i] <- W3[i,]%*%(dataF2-mean3[i])^2 }
    s0[5]<-ss1+ss2+ss3+swx1[1]+swx2[2]+swx3[1]+swx3[3]
    n0[5]<-n_samP1+n_samF1+n_samP2+mix_pi1[1]*n_samB1+mix_pi2[2]*n_samB2+(mix_pi3[1]+mix_pi3[3])*n_samF2
    s0[6]<-swx1[2]+swx2[1]+swx3[2]
    n0[6]<-mix_pi1[2]*n_samB1+mix_pi2[1]*n_samB2+mix_pi3[2]*n_samF2
    aaa0<-sigma
    aa1<-0.5*(mean[1]-mean[3])
    aa1<-0.5*aa1*aa1/n_fam
    aa3<-1000;n_iter<-0
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      aa2<-sigma/(sigma+aa1)
      sigma<-(s0[5]+aa2^2*s0[6])/(n0[5]+aa2*n0[6])
      aa3<-abs(sigma-aaa0)
      aaa0<-sigma
      if(n_iter>20)break
    }
    sigma1[2]<-sigma+aa1
    sigma2[1]<-sigma3[2]<-sigma1[2]
    sigma1[1]<-sigma2[2]<-sigma3[1]<-sigma3[3]<-sigma
    if(sum(sigma < 1e-30)>=1){break}
    #############################criteria for iterations to stop########################
    L1<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,mean[3],sqrt(sigma))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mix_pi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mix_pi3)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc<-L0
  AIC<--2*abc+2*3
  meanP1<-mean[1];meanF1<-mean[2];meanP2<-mean[3];sigma0<-sigma
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*n_samP1)
  P1bmw <- matrix(0,n_samP1,1)
  P1gg <- (dataP1 - mean[1])/sqrt(as.vector(sigma))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < n_samP1){P1bmw <- P1bmw+runif(n_samP1)/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  P1u<- as.matrix(c(12*n_samP1*((P1dd[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((P1dd[2]/n_samP1-1/3)^2),180*n_samP1*((P1dd[3]/n_samP1-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[2])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),P1D))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*n_samF1)
  F1bmw <- matrix(0,n_samF1,1)
  F1gg <- (dataF1 - mean[2])/sqrt(as.vector(sigma))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < n_samF1){F1bmw <- F1bmw+runif(n_samF1)/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  F1u<- as.matrix(c(12*n_samF1*((F1dd[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((F1dd[2]/n_samF1-1/3)^2),180*n_samF1*((F1dd[3]/n_samF1-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[2])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),F1D))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*n_samP2)
  P2bmw <- matrix(0,n_samP2,1)
  P2gg <- (dataP2 - mean[3])/sqrt(as.vector(sigma))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < n_samP2){P2bmw <- P2bmw+runif(n_samP2)/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  P2u<- as.matrix(c(12*n_samP2*((P2dd[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((P2dd[2]/n_samP2-1/3)^2),180*n_samP2*((P2dd[3]/n_samP2-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[2])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),P2D))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1);
  B1w1<-1/(12*n_samB1)
  B1bmw <- matrix(0,n_samB1,1); B1bmwsl <- matrix(0,n_samB1,2)
  for(i in 1:2){
    B1gg <- (dataB1 - mean1[i])/sqrt(sigma1[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi1[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B1P2)))[1]
  if(nn < n_samB1){B1P2 <- B1P2+runif(n_samB1)/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*n_samB1) + sum((B1P2 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  B1u <- as.matrix(c(12*n_samB1*((B1dd[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((B1dd[2]/n_samB1-1/3)^2),180*n_samB1*((B1dd[3]/n_samB1-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[2])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),B1D))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])

  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2);
  B2w1<-1/(12*n_samB2)
  B2bmw <- matrix(0,n_samB2,1); B2bmwsl <- matrix(0,n_samB2,2)
  for(i in 1:2){
    B2gg <- (dataB2 - mean2[i])/sqrt(sigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B2P2)))[1]
  if(nn < n_samB2){B2P2 <- B2P2+runif(n_samB2)/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*n_samB2) + sum((B2P2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  B2u <- as.matrix(c(12*n_samB2*((B2dd[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((B2dd[2]/n_samB2-1/3)^2),180*n_samB2*((B2dd[3]/n_samB2-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[2])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),B2D))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])

  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2);
  F2w1<-1/(12*n_samF2)
  F2bmw <- matrix(0,n_samF2,1); F2bmwsl <- matrix(0,n_samF2,3)
  for(i in 1:3){
    F2gg <- (dataF2 - mean3[i])/sqrt(sigma3[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi3[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < n_samF2){F2P2 <- F2P2+runif(n_samF2)/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*n_samF2) + sum((F2P2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  F2u <- as.matrix(c(12*n_samF2*((F2dd[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((F2dd[2]/n_samF2-1/3)^2),180*n_samF2*((F2dd[3]/n_samF2-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[2])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),F2D))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])

  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)

  ########first order parameters################
  aa<-matrix(c(1,1,1,0,1,-1),3,2,byrow=T)
  mm<-mean[c(1:3)]
  B<-solve(crossprod(aa,aa))%*%crossprod(aa,mm)
  #######second order parameters################
  jj1<-sigmaB1-sigma1[1]
  if (jj1<0 ) {jj1<-0}
  ll1<-jj1/sigmaB1
  jj2<-sigmaB2-sigma2[2]
  if (jj2<0) {jj2<-0}
  ll2<-jj2/sigmaB2
  jj3<-sigmaF2-sigma3[1]
  if (jj3<0) {jj3<-0}
  ll3<-jj3/sigmaF2

  output <- data.frame("1MG-A",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4), round(sigma0,4),round(t(mean1),4)," "," ",round(t(sigma1),4)," "," ",
                       round(t(mix_pi1),4)," "," ",round(t(mean2),4)," "," ",round(t(sigma2),4)," "," ",round(t(mix_pi2),4)," "," ",
                       round(t(mean3),4)," "," "," "," "," "," ",round(t(sigma3),4)," "," "," "," "," "," ",round(t(mix_pi3),4)," "," "," "," "," "," ",
                       round(B[1],4)," "," "," "," "," ",round(B[2],4)," "," "," "," "," "," "," "," "," ",
                       round(jj1,4),round(ll1*100,4)," "," ",round(jj2,4),round(ll2*100,4)," "," ",round(jj3,4),round(ll3*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi1,mi2,mi3)
  return(OUTPUT)
}

##################1MG-EAD(A-3)#########################
G6FModelFun[[3]] <- function(K1,logL,df11,df21,df31,df41,df51,df61,G6Ftext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1];n_samF1<-dim(dataF1)[1];n_samP2<-dim(dataP2)[1]
  n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1];n_samF2<-dim(dataF2)[1]
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2),sum(dataF2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2),sum(dataF2^2)))
  ss<-matrix(0,3,1);ss[1]<-s[1]-sumx[1]^2/n_samP1;ss[2]<-s[2]-sumx[2]^2/n_samF1;ss[3]<-s[3]-sumx[3]^2/n_samP2
  mean<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataB1),mean(dataB2),mean(dataF2)))
  sigma0<-sum(ss)/(n_samP1+n_samF1+n_samP2-3)
  sigmaB1<-var(dataB1);sigma40<-sigmaB1
  sigmaB2<-var(dataB2);sigma50<-sigmaB2
  sigmaF2<-var(dataF2);sigma60<-sigmaF2
  sigmaP1<-sigmaF1<-sigmaP2<-sigma0
  m_esp <- 0.0001;n_fam <- as.numeric(G6Ftext2)
  ###############procedure start###########################
  mi1<-matrix(0.5,2,1);sigma1<-matrix(0,2,1)
  mi2<-matrix(0.5,2,1);sigma2<-matrix(0,2,1)
  mi3<-as.matrix(c(0.25,0.5,0.25));sigma3<-matrix(0,3,1)
  sigma<-sigma0;sigma1[1]<-sigma
  a1<-sqrt(sigmaB1/n_samB1)
  if (mean[1]<mean[3]) a1<--a1
  mean1<-as.matrix(c(mean[4]+2*a1,mean[4]-2*a1))
  sigma2[2]<-sigma
  a2<-sqrt(sigmaB2/n_samB2)
  if (mean[1]<mean[3]) a2<--a2
  mean2<-as.matrix(c(mean[5]+2*a2,mean[5]-2*a2))
  sigma3[1]<-sigma3[3]<-sigma
  a3<-sqrt(sigmaF2/n_samF2)
  if (mean[1]<mean[3]) a3<--a3
  mean3<-as.matrix(c(mean[6]+2*a3,mean[6],mean[6]-2*a3))
  a1<-(5*mean[1]-7*mean[3]+2*mean1[2])/13  #additive effect.
  sigma1[2]<-sigma+(0.75*a1^2)/n_fam
  sigma2[1]<-sigma3[2]<-sigma1[2]
  L0<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma0))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma0))))+
    sum(log(dnorm(dataP2,mean[3],sqrt(sigma0))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mi1)))+
    sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mi3)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,2,n_samB1);      swx1 <- matrix(0,2,1)
  W2 <- matrix(0,2,n_samB2);      swx2 <- matrix(0,2,1)
  W3 <- matrix(0,3,n_samF2);      swx3 <- matrix(0,3,1)
  s0<-matrix(0,5,1);n0<-matrix(0,5,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:2) { W1[i,] <- mi1[i]*dnorm(dataB1,mean1[i],sqrt(sigma1[i]))/dmixnorm(dataB1,mean1,sqrt(sigma1),mi1) }
    mix_pi1 <- as.matrix(rowSums(W1)/n_samB1)
    sumwx1 <- W1%*%dataB1
    for(i in 1:2) { W2[i,] <- mi2[i]*dnorm(dataB2,mean2[i],sqrt(sigma2[i]))/dmixnorm(dataB2,mean2,sqrt(sigma2),mi2) }
    mix_pi2 <- as.matrix(rowSums(W2)/n_samB2)
    sumwx2 <- W2%*%dataB2
    for(i in 1:3) { W3[i,] <- mi3[i]*dnorm(dataF2,mean3[i],sqrt(sigma3[i]))/dmixnorm(dataF2,mean3,sqrt(sigma3),mi3) }
    mix_pi3 <- as.matrix(rowSums(W3)/n_samF2)
    sumwx3 <- W3%*%dataF2
    ######obtain means##################
    aaa0<-0
    s0[1]<-sumx[1]+sumx[2]+sumwx1[1]+sumwx3[1];s0[2]<-sumx[3]+sumwx2[2]+sumwx3[3]
    s0[3]<-sumwx1[2]+sumwx2[1]+sumwx3[2];n0[1]<-n_samP1+n_samF1+mix_pi1[1]*n_samB1+mix_pi3[1]*n_samF2
    n0[2]<-n_samP2+mix_pi2[2]*n_samB2+mix_pi3[3]*n_samF2;n0[3]<-mix_pi1[2]*n_samB1+mix_pi2[1]*n_samB2+mix_pi3[2]*n_samF2
    n0[c(1:3)][abs(n0[c(1:3)])<0.00000001]<-0.000001
    aa2<-3*s0[1]/n0[1]+s0[2]/n0[2]-4*s0[3]/n0[3];aa3<-9*sigma/n0[1]+sigma/n0[2]
    aa1<-1000;n_iter<-0
    while (aa1>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-(5*mean[1]-7*mean[3]+2*mean1[2])/13
      sigma1[2]<-sigma+0.75*aa1^2/n_fam
      aaa1<-aa2/(aa3+16*sigma1[2]/n0[3])   #coefficient in restricted condition.
      mean[1]<-(s0[1]-3*aaa1*sigma)/n0[1]
      mean[3]<-(s0[2]-aaa1*sigma)/n0[2]
      mean1[2]<-(s0[3]+4*aaa1*sigma1[2])/n0[3]
      aa1<-abs(aaa1-aaa0)
      aaa0<-aaa1
      if(n_iter>20)break
    }
    mean[2]<-mean1[1]<-mean3[1]<-mean[1];mean2[1]<-mean3[2]<-mean1[2];mean2[2]<-mean3[3]<-mean[3]
    ##########obtain variance##############################
    ss1<-sum((dataP1-mean[1])^2);ss3<-sum((dataP2-mean[3])^2);ss2<-sum((dataF1-mean[2])^2)
    for(i in 1:2) {swx1[i] <- W1[i,]%*%(dataB1-mean1[i])^2 };for(i in 1:2) {swx2[i] <- W2[i,]%*%(dataB2-mean2[i])^2 };for(i in 1:3) {swx3[i] <- W3[i,]%*%(dataF2-mean3[i])^2 }
    s0[4]<-ss1+ss2+ss3+swx1[1]+swx2[2]+swx3[1]+swx3[3]
    n0[4]<-n_samP1+n_samF1+n_samP2+mix_pi1[1]*n_samB1+mix_pi2[2]*n_samB2+(mix_pi3[1]+mix_pi3[3])*n_samF2
    aaa0<-sigma
    s0[5]<-swx1[2]+swx2[1]+swx3[2]
    n0[5]<-mix_pi1[2]*n_samB1+mix_pi2[1]*n_samB2+mix_pi3[2]*n_samF2
    aa1<-(5*mean[1]-7*mean[3]+2*mean1[2])/13
    aa1<-0.75*aa1^2/n_fam
    aa3<-1000;n_iter<-0
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      aa2<-sigma/(sigma+aa1)
      sigma<-(s0[4]+aa2^2*s0[5])/(n0[4]+aa2*n0[5])
      aa3<-abs(sigma-aaa0)
      aaa0<-sigma
      if(n_iter>20)break
    }
    sigma1[2]<-sigma+aa1
    sigma2[1]<-sigma3[2]<-sigma1[2]
    sigma1[1]<-sigma2[2]<-sigma3[1]<-sigma3[3]<-sigma
    if(sum(sigma < 1e-30)>=1){break}
    ##############################criteria for iterations to stop########################
    L1<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,mean[3],sqrt(sigma))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mix_pi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mix_pi3)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc<-L0
  AIC<--2*abc+2*3
  meanP1<-mean[1];meanF1<-mean[2];meanP2<-mean[3]
  sigma0<-sigma
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*n_samP1)
  P1bmw <- matrix(0,n_samP1,1)
  P1gg <- (dataP1 - mean[1])/sqrt(as.vector(sigma))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < n_samP1){P1bmw <- P1bmw+runif(n_samP1)/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  P1u<- as.matrix(c(12*n_samP1*((P1dd[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((P1dd[2]/n_samP1-1/3)^2),180*n_samP1*((P1dd[3]/n_samP1-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[2])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),P1D))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*n_samF1)
  F1bmw <- matrix(0,n_samF1,1)
  F1gg <- (dataF1 - mean[2])/sqrt(as.vector(sigma))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < n_samF1){F1bmw <- F1bmw+runif(n_samF1)/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  F1u<- as.matrix(c(12*n_samF1*((F1dd[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((F1dd[2]/n_samF1-1/3)^2),180*n_samF1*((F1dd[3]/n_samF1-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[2])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),F1D))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*n_samP2)
  P2bmw <- matrix(0,n_samP2,1)
  P2gg <- (dataP2 - mean[3])/sqrt(as.vector(sigma))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < n_samP2){P2bmw <- P2bmw+runif(n_samP2)/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  P2u<- as.matrix(c(12*n_samP2*((P2dd[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((P2dd[2]/n_samP2-1/3)^2),180*n_samP2*((P2dd[3]/n_samP2-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[2])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),P2D))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1);
  B1w1<-1/(12*n_samB1)
  B1bmw <- matrix(0,n_samB1,1); B1bmwsl <- matrix(0,n_samB1,2)
  for(i in 1:2){
    B1gg <- (dataB1 - mean1[i])/sqrt(sigma1[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi1[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B1P2)))[1]
  if(nn < n_samB1){B1P2 <- B1P2+runif(n_samB1)/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*n_samB1) + sum((B1P2 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  B1u <- as.matrix(c(12*n_samB1*((B1dd[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((B1dd[2]/n_samB1-1/3)^2),180*n_samB1*((B1dd[3]/n_samB1-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[2])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),B1D))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])

  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2);
  B2w1<-1/(12*n_samB2)
  B2bmw <- matrix(0,n_samB2,1); B2bmwsl <- matrix(0,n_samB2,2)
  for(i in 1:2){
    B2gg <- (dataB2 - mean2[i])/sqrt(sigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B2P2)))[1]
  if(nn < n_samB2){B2P2 <- B2P2+runif(n_samB2)/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*n_samB2) + sum((B2P2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  B2u <- as.matrix(c(12*n_samB2*((B2dd[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((B2dd[2]/n_samB2-1/3)^2),180*n_samB2*((B2dd[3]/n_samB2-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[2])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),B2D))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])

  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2);
  F2w1<-1/(12*n_samF2)
  F2bmw <- matrix(0,n_samF2,1); F2bmwsl <- matrix(0,n_samF2,3)
  for(i in 1:3){
    F2gg <- (dataF2 - mean3[i])/sqrt(sigma3[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi3[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < n_samF2){F2P2 <- F2P2+runif(n_samF2)/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*n_samF2) + sum((F2P2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  F2u <- as.matrix(c(12*n_samF2*((F2dd[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((F2dd[2]/n_samF2-1/3)^2),180*n_samF2*((F2dd[3]/n_samF2-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[2])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),F2D))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])

  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)
  ########first order parameters#################
  aa<-matrix(c(1,1,1,-1,1,0.5),3,2,byrow=T)
  mm<-as.matrix(c(mean[1],mean[3],mean1[2]))
  B<-solve(crossprod(aa,aa))%*%crossprod(aa,mm)
  ##########second order parameters########
  jj1<-sigmaB1-sigma1[1]
  if (jj1<0 ) {jj1<-0}
  ll1<-jj1/sigmaB1
  jj2<-sigmaB2-sigma2[2]
  if (jj2<0) {jj2<-0}
  ll2<-jj2/sigmaB2
  jj3<-sigmaF2-sigma3[1]
  if (jj3<0) {jj3<-0}
  ll3<-jj3/sigmaF2

  output <- data.frame("1MG-EAD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4), round(sigma0,4),round(t(mean1),4)," "," ",round(t(sigma1),4)," "," ",
                       round(t(mix_pi1),4)," "," ",round(t(mean2),4)," "," ",round(t(sigma2),4)," "," ",round(t(mix_pi2),4)," "," ",
                       round(t(mean3),4)," "," "," "," "," "," ",round(t(sigma3),4)," "," "," "," "," "," ",round(t(mix_pi3),4)," "," "," "," "," "," ",
                       round(B[1],4)," "," "," "," "," ",round(B[2],4)," ",round(B[2],4)," "," "," "," "," "," "," ",
                       round(jj1,4),round(ll1*100,4)," "," ",round(jj2,4),round(ll2*100,4)," "," ",round(jj3,4),round(ll3*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi1,mi2,mi3)
  return(OUTPUT)
}

############1MG-NCD(A-4)####################
G6FModelFun[[4]] <- function(K1,logL,df11,df21,df31,df41,df51,df61,G6Ftext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1];n_samF1<-dim(dataF1)[1];n_samP2<-dim(dataP2)[1]
  n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1];n_samF2<-dim(dataF2)[1]
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2),sum(dataF2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2),sum(dataF2^2)))
  ss<-matrix(0,3,1);ss[1]<-s[1]-sumx[1]^2/n_samP1;ss[2]<-s[2]-sumx[2]^2/n_samF1;ss[3]<-s[3]-sumx[3]^2/n_samP2
  mean<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataB1),mean(dataB2),mean(dataF2)))
  sigma0<-sum(ss)/(n_samP1+n_samF1+n_samP2-3)
  sigmaB1<-var(dataB1);sigma40<-sigmaB1
  sigmaB2<-var(dataB2);sigma50<-sigmaB2
  sigmaF2<-var(dataF2);sigma60<-sigmaF2
  sigmaP1<-sigmaF1<-sigmaP2<-sigma0
  m_esp <- 0.0001;n_fam <- as.numeric(G6Ftext2)
  ###############procedure start###########################
  mi1<-matrix(0.5,2,1);sigma1<-matrix(0,2,1)
  mi2<-matrix(0.5,2,1);sigma2<-matrix(0,2,1)
  mi3<-as.matrix(c(0.25,0.5,0.25));sigma3<-matrix(0,3,1)
  sigma<-sigma0;sigma1[1]<-sigma
  a1<-sqrt(sigmaB1/n_samB1)
  if (mean[1]<mean[3]) a1<--a1
  mean1<-as.matrix(c(mean[4]+2*a1,mean[4]-2*a1))
  sigma2[2]<-sigma
  a2<-sqrt(sigmaB2/n_samB2)
  if (mean[1]<mean[3]) a2<--a2
  mean2<-as.matrix(c(mean[5]+2*a2,mean[5]-2*a2))
  sigma3[1]<-sigma3[3]<-sigma
  a3<-sqrt(sigmaF2/n_samF2)
  if (mean[1]<mean[3]) a3<--a3
  mean3<-as.matrix(c(mean[6]+2*a3,mean[6],mean[6]-2*a3))
  b1<-(7*mean[1]-5*mean[3]-2*mean1[2])/13  # additive effect.
  sigma1[2]<-sigma+0.75*b1^2/n_fam;sigma2[1]<-sigma3[2]<-sigma1[2]
  L0<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma0))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma0))))+
    sum(log(dnorm(dataP2,mean[3],sqrt(sigma0))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mi1)))+
    sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mi3)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,2,n_samB1);      swx1 <- matrix(0,2,1)
  W2 <- matrix(0,2,n_samB2);      swx2 <- matrix(0,2,1)
  W3 <- matrix(0,3,n_samF2);      swx3 <- matrix(0,3,1)
  s0<-matrix(0,5,1);n0<-matrix(0,5,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:2) { W1[i,] <- mi1[i]*dnorm(dataB1,mean1[i],sqrt(sigma1[i]))/dmixnorm(dataB1,mean1,sqrt(sigma1),mi1) }
    mix_pi1 <- as.matrix(rowSums(W1)/n_samB1)
    sumwx1 <- W1%*%dataB1
    for(i in 1:2) { W2[i,] <- mi2[i]*dnorm(dataB2,mean2[i],sqrt(sigma2[i]))/dmixnorm(dataB2,mean2,sqrt(sigma2),mi2) }
    mix_pi2 <- as.matrix(rowSums(W2)/n_samB2)
    sumwx2 <- W2%*%dataB2
    for(i in 1:3) { W3[i,] <- mi3[i]*dnorm(dataF2,mean3[i],sqrt(sigma3[i]))/dmixnorm(dataF2,mean3,sqrt(sigma3),mi3) }
    mix_pi3 <- as.matrix(rowSums(W3)/n_samF2)
    sumwx3 <- W3%*%dataF2
    #######obtain means#################
    aaa0<-0
    s0[1]<-sumx[1]+sumwx1[1]+sumwx3[1];s0[2]<-sumx[2]+sumx[3]+sumwx2[2]+sumwx3[3]
    s0[3]<-sumwx1[2]+sumwx2[1]+sumwx3[2];n0[1]<-n_samP1+mix_pi1[1]*n_samB1+mix_pi3[1]*n_samF2
    n0[2]<-n_samF1+n_samP2+mix_pi2[2]*n_samB2+mix_pi3[3]*n_samF2;n0[3]<-mix_pi1[2]*n_samB1+mix_pi2[1]*n_samB2+mix_pi3[2]*n_samF2
    n0[c(1:3)][abs(n0[c(1:3)])<0.00000001]<-0.000001
    aa2<-s0[1]/n0[1]+3*s0[2]/n0[2]-4*s0[3]/n0[3]
    aa3<-sigma*(1/n0[1]+9/n0[2])
    aa1<-1000;n_iter<-0
    while (aa1>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-(7*mean[1]-5*mean[3]-2*mean1[2])/13
      sigma1[2]<-sigma+0.75*aa1^2/n_fam
      aaa1<-aa2/(aa3+16*sigma1[2]/n0[3])  # coefficient in restricted condition.
      mean[1]<-(s0[1]-aaa1*sigma)/n0[1]
      mean[3]<-(s0[2]-3*aaa1*sigma)/n0[2]
      mean1[2]<-(s0[3]+4*aaa1*sigma1[2])/n0[3]
      aa1<-abs(aaa1-aaa0)
      aaa0<-aaa1
      if(n_iter>20)break
    }
    mean[2]<-mean2[2]<-mean3[3]<-mean[3];mean1[1]<-mean3[1]<-mean[1];mean2[1]<-mean3[2]<-mean1[2]
    ##############################################
    ss1<-sum((dataP1-mean[1])^2);ss3<-sum((dataP2-mean[3])^2);ss2<-sum((dataF1-mean[2])^2)
    for(i in 1:2) {swx1[i] <- W1[i,]%*%(dataB1-mean1[i])^2 };for(i in 1:2) {swx2[i] <- W2[i,]%*%(dataB2-mean2[i])^2 };for(i in 1:3) {swx3[i] <- W3[i,]%*%(dataF2-mean3[i])^2 }
    s0[4]<-ss1+ss2+ss3+swx1[1]+swx2[2]+swx3[1]+swx3[3]
    n0[4]<-n_samP1+n_samF1+n_samP2+mix_pi1[1]*n_samB1+mix_pi2[2]*n_samB2+(mix_pi3[1]+mix_pi3[3])*n_samF2
    s0[5]<-swx1[2]+swx2[1]+swx3[2]
    n0[5]<-mix_pi1[2]*n_samB1+mix_pi2[1]*n_samB2+mix_pi3[2]*n_samF2
    aaa0<-sigma
    aa1<-(7*mean[1]-5*mean[3]-2*mean1[2])/13
    aa1<-0.75*aa1^2/n_fam
    aa3<-1000;n_iter<-0
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      aa2<-sigma/(sigma+aa1)
      sigma<-(s0[4]+aa2^2*s0[5])/(n0[4]+aa2*n0[5])
      aa3<-abs(sigma-aaa0)
      aaa0<-sigma
      if(n_iter>20)break
    }
    sigma1[2]<-sigma+aa1;sigma2[1]<-sigma3[2]<-sigma1[2]
    sigma1[1]<-sigma2[2]<-sigma3[1]<-sigma3[3]<-sigma
    if(sum(sigma < 1e-30)>=1){break}
    ##############################criteria for iterations to stop########################
    L1<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,mean[3],sqrt(sigma))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mix_pi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mix_pi3)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*3
  meanP1<-mean[1];meanF1<-mean[2];meanP2<-mean[3]
  sigma0<-sigma
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*n_samP1)
  P1bmw <- matrix(0,n_samP1,1)
  P1gg <- (dataP1 - mean[1])/sqrt(as.vector(sigma))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < n_samP1){P1bmw <- P1bmw+runif(n_samP1)/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  P1u<- as.matrix(c(12*n_samP1*((P1dd[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((P1dd[2]/n_samP1-1/3)^2),180*n_samP1*((P1dd[3]/n_samP1-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[2])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),P1D))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*n_samF1)
  F1bmw <- matrix(0,n_samF1,1)
  F1gg <- (dataF1 - mean[2])/sqrt(as.vector(sigma))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < n_samF1){F1bmw <- F1bmw+runif(n_samF1)/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  F1u<- as.matrix(c(12*n_samF1*((F1dd[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((F1dd[2]/n_samF1-1/3)^2),180*n_samF1*((F1dd[3]/n_samF1-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[2])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),F1D))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*n_samP2)
  P2bmw <- matrix(0,n_samP2,1)
  P2gg <- (dataP2 - mean[3])/sqrt(as.vector(sigma))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < n_samP2){P2bmw <- P2bmw+runif(n_samP2)/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  P2u<- as.matrix(c(12*n_samP2*((P2dd[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((P2dd[2]/n_samP2-1/3)^2),180*n_samP2*((P2dd[3]/n_samP2-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[2])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),P2D))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1);
  B1w1<-1/(12*n_samB1)
  B1bmw <- matrix(0,n_samB1,1); B1bmwsl <- matrix(0,n_samB1,2)
  for(i in 1:2){
    B1gg <- (dataB1 - mean1[i])/sqrt(sigma1[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi1[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B1P2)))[1]
  if(nn < n_samB1){B1P2 <- B1P2+runif(n_samB1)/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*n_samB1) + sum((B1P2 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  B1u <- as.matrix(c(12*n_samB1*((B1dd[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((B1dd[2]/n_samB1-1/3)^2),180*n_samB1*((B1dd[3]/n_samB1-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[2])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),B1D))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])

  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2);
  B2w1<-1/(12*n_samB2)
  B2bmw <- matrix(0,n_samB2,1); B2bmwsl <- matrix(0,n_samB2,2)
  for(i in 1:2){
    B2gg <- (dataB2 - mean2[i])/sqrt(sigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B2P2)))[1]
  if(nn < n_samB2){B2P2 <- B2P2+runif(n_samB2)/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*n_samB2) + sum((B2P2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  B2u <- as.matrix(c(12*n_samB2*((B2dd[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((B2dd[2]/n_samB2-1/3)^2),180*n_samB2*((B2dd[3]/n_samB2-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[2])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),B2D))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])

  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2);
  F2w1<-1/(12*n_samF2)
  F2bmw <- matrix(0,n_samF2,1); F2bmwsl <- matrix(0,n_samF2,3)
  for(i in 1:3){
    F2gg <- (dataF2 - mean3[i])/sqrt(sigma3[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi3[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < n_samF2){F2P2 <- F2P2+runif(n_samF2)/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*n_samF2) + sum((F2P2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  F2u <- as.matrix(c(12*n_samF2*((F2dd[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((F2dd[2]/n_samF2-1/3)^2),180*n_samF2*((F2dd[3]/n_samF2-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[2])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),F2D))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])

  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)

  ###########first order parameters#######################
  aa<-matrix(c(1,1,1,-1,1,-0.5),3,2,byrow=T)
  mm<-as.matrix(c(mean[1],mean[3],mean1[2]))
  B<-solve(crossprod(aa,aa))%*%crossprod(aa,mm)
  ###########second order parameters#######################
  jj1<-sigmaB1-sigma1[1]
  if (jj1<0 ) {jj1<-0}
  ll1<-jj1/sigmaB1
  jj2<-sigmaB2-sigma2[2]
  if (jj2<0) {jj2<-0}
  ll2<-jj2/sigmaB2
  jj3<-sigmaF2-sigma3[1]
  if (jj3<0) {jj3<-0}
  ll3<-jj3/sigmaF2

  output <- data.frame("1MG-NCD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4), round(sigma0,4),round(t(mean1),4)," "," ",round(t(sigma1),4)," "," ",
                       round(t(mix_pi1),4)," "," ",round(t(mean2),4)," "," ",round(t(sigma2),4)," "," ",round(t(mix_pi2),4)," "," ",
                       round(t(mean3),4)," "," "," "," "," "," ",round(t(sigma3),4)," "," "," "," "," "," ",round(t(mix_pi3),4)," "," "," "," "," "," ",
                       round(B[1],4)," "," "," "," "," ",round(B[2],4)," ",round(-B[2],4)," "," "," "," "," "," "," ",
                       round(jj1,4),round(ll1*100,4)," "," ",round(jj2,4),round(ll2*100,4)," "," ",round(jj3,4),round(ll3*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi1,mi2,mi3)
  return(OUTPUT)
}

##############2MG-ADI(B-1)######################
G6FModelFun[[5]] <- function(K1,logL,df11,df21,df31,df41,df51,df61,G6Ftext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1];n_samF1<-dim(dataF1)[1];n_samP2<-dim(dataP2)[1]
  n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1];n_samF2<-dim(dataF2)[1]
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2),sum(dataF2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2),sum(dataF2^2)))
  ss<-matrix(0,3,1);ss[1]<-s[1]-sumx[1]^2/n_samP1;ss[2]<-s[2]-sumx[2]^2/n_samF1;ss[3]<-s[3]-sumx[3]^2/n_samP2
  mean<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataB1),mean(dataB2),mean(dataF2)))
  sigma0<-sum(ss)/(n_samP1+n_samF1+n_samP2-3)
  sigmaB1<-var(dataB1);sigma40<-sigmaB1
  sigmaB2<-var(dataB2);sigma50<-sigmaB2
  sigmaF2<-var(dataF2);sigma60<-sigmaF2
  sigmaP1<-sigmaF1<-sigmaP2<-sigma0
  m_esp <- 0.0001;n_fam <- as.numeric(G6Ftext2)
  ###############procedure start###########################
  mi1<-matrix(0.25,4,1);sigma1<-matrix(0,4,1)
  mi2<-matrix(0.25,4,1);sigma2<-matrix(0,4,1)
  mi3<-as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  sigma3<-matrix(0,9,1)
  sigma<-sigma0
  a1<-sqrt(sigma/n_samB1)
  if (mean[1]<mean[3]) a1<--a1
  mean1<-as.matrix(c(mean[4]+2.4*a1,mean[4]+0.8*a1,mean[4]-0.8*a1,mean[4]-2.4*a1))
  a2<-sqrt(sigmaB2/n_samB2)
  if (mean[1]<mean[3]) a2<--a2
  mean2<-as.matrix(c(mean[5]+2.4*a2,mean[5]+0.8*a2,mean[5]-0.8*a2,mean[5]-2.4*a2))
  a3<-sqrt(sigmaF2/n_samF2)
  if (mean[1]<mean[3]) a3<--a3
  mean3<-as.matrix(c(mean[6]+2.4*a3,mean[6]+1.8*a3,mean[6]+1.2*a3,mean[6]+0.6*a3,mean[6],mean[6]-0.6*a3,mean[6]-1.2*a3,mean[6]-1.8*a3,mean[6]-2.4*a3))
  gs<-matrix(0,8,1)
  gs[1]<-0.25*(mean[1]-mean[3]+mean3[3]-mean3[7])       	 #da.
  gs[2]<-0.25*(mean[1]-mean[3]-mean3[3]+mean3[7])  	 #db.
  gs[3]<-(-8*mean[1]-2*mean[2]-8*mean[3]-2*mean1[2]+15*mean1[3]+8*mean1[4]+15*mean2[2]-2*mean2[3]-8*mean3[3]-8*mean3[7])/17
  # ha.
  gs[4]<-(-8*mean[1]-2*mean[2]-8*mean[3]+15*mean1[2]-2*mean1[3]+8*mean1[4]-2*mean2[2]+15*mean2[3]-8*mean3[3]-8*mean3[7])/17
  # hb.
  gs[5]<-0.25*(mean[1]+mean[3]-mean3[3]-mean3[7])  #  i.
  gs[6]<--0.5*mean[1]+0.5*mean[3]+mean1[2]-mean2[3]-0.5*mean3[3]+0.5*mean3[7]
  #  jab.
  gs[7]<--0.5*mean[1]+0.5*mean[3]+mean1[3]-mean2[2]+0.5*mean3[3]-0.5*mean3[7]
  #  jba.
  gs[8]<-(12*mean[1]+20*mean[2]+12*mean[3]-14*mean1[2]-14*mean1[3]-12*mean1[4]-14*mean2[2]-14*mean2[3]+12*mean3[3]+12*mean3[7])/17
  #  l.
  g_aa1<-(0.5*(gs[2]+gs[5])^2+0.25*(gs[4]+gs[6])^2)/n_fam
  #   0.5(db+i)**2+0.25(hb+jab)**2.
  g_aa2<-(0.5*(gs[1]+gs[5])^2+0.25*(gs[3]+gs[7])^2)/n_fam
  #   0.5(da+i)**2+0.25(ha+jba)**2.
  g_aa3<-(0.5*(gs[1]-gs[5])^2+0.25*(gs[3]-gs[7])^2)/n_fam
  #   0.5(da-i)**2+0.25(ha-jba)**2.
  g_aa4<-(0.5*(gs[2]-gs[5])^2+0.25*(gs[4]-gs[6])^2)/n_fam
  #   0.5(db-i)**2+0.25(hb-jab)**2.
  g_aa5<-0.25*(gs[1]^2+gs[2]^2+gs[5]^2+(gs[1]+gs[6])^2+(gs[2]+gs[7])^2+(gs[3]+gs[8]/2)^2+(gs[4]+gs[8]/2)^2+gs[8]^2/4)/n_fam

  sigma1[1]<-sigma;sigma1[2]<-sigma+g_aa1;sigma1[3]<-sigma+g_aa2;sigma1[4]<-sigma+g_aa5
  sigma2[1]<-sigma1[4];sigma2[2]<-sigma+g_aa3;sigma2[3]<-sigma+g_aa4;sigma2[4]<-sigma
  sigma3[1]<-sigma3[3]<-sigma3[7]<-sigma3[9]<-sigma
  sigma3[2]<-sigma1[2];sigma3[4]<-sigma1[3];sigma3[5]<-sigma1[4];sigma3[6]<-sigma2[2];sigma3[8]<-sigma2[3]
  L0<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma0))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma0))))+
    sum(log(dnorm(dataP2,mean[3],sqrt(sigma0))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mi1)))+
    sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mi3)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,4,n_samB1);      swx1 <- matrix(0,4,1)
  W2 <- matrix(0,4,n_samB2);      swx2 <- matrix(0,4,1)
  W3 <- matrix(0,9,n_samF2);      swx3 <- matrix(0,9,1)
  s0<-matrix(0,16,1);n0<-matrix(0,16,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:4) { W1[i,] <- mi1[i]*dnorm(dataB1,mean1[i],sqrt(sigma1[i]))/dmixnorm(dataB1,mean1,sqrt(sigma1),mi1) }
    mix_pi1 <- as.matrix(rowSums(W1)/n_samB1)
    sumwx1 <- W1%*%dataB1
    for(i in 1:4) { W2[i,] <- mi2[i]*dnorm(dataB2,mean2[i],sqrt(sigma2[i]))/dmixnorm(dataB2,mean2,sqrt(sigma2),mi2) }
    mix_pi2 <- as.matrix(rowSums(W2)/n_samB2)
    sumwx2 <- W2%*%dataB2
    for(i in 1:9) { W3[i,] <- mi3[i]*dnorm(dataF2,mean3[i],sqrt(sigma3[i]))/dmixnorm(dataF2,mean3,sqrt(sigma3),mi3) }
    mix_pi3 <- as.matrix(rowSums(W3)/n_samF2)
    sumwx3 <- W3%*%dataF2
    ######obtain means#####################
    aaa0<-0
    s0[1]<-sumx[1]+sumwx1[1]+sumwx3[1];s0[2]<-sumx[2]
    s0[3]<-sumx[3]+sumwx2[4]+sumwx3[9];s0[4]<-sumwx1[2]+sumwx3[2]
    s0[5]<-sumwx1[3]+sumwx3[4];s0[6]<-sumwx1[4]+sumwx2[1]+sumwx3[5]
    s0[7]<-sumwx2[2]+sumwx3[6];s0[8]<-sumwx2[3]+sumwx3[8]
    s0[9]<-sumwx3[3];s0[10]<-sumwx3[7]
    n0[1]<-n_samP1+mix_pi1[1]*n_samB1+mix_pi3[1]*n_samF2;n0[2]<-n_samF1
    n0[3]<-n_samP2+mix_pi2[4]*n_samB2+mix_pi3[9]*n_samF2;n0[4]<-mix_pi1[2]*n_samB1+mix_pi3[2]*n_samF2
    n0[5]<-mix_pi1[3]*n_samB1+mix_pi3[4]*n_samF2;n0[6]<-mix_pi1[4]*n_samB1+mix_pi2[1]*n_samB2+mix_pi3[5]*n_samF2
    n0[7]<-mix_pi2[2]*n_samB2+mix_pi3[6]*n_samF2;n0[8]<-mix_pi2[3]*n_samB2+mix_pi3[8]*n_samF2
    n0[9]<-mix_pi3[3]*n_samF2;n0[10]<-mix_pi3[7]*n_samF2
    s0[c(1:10)][abs(s0[c(1:10)])<0.000001]<-0.000001
    n0[c(1:10)][abs(n0[c(1:10)])<0.000001]<-0.000001
    aa3<-s0[1]/n0[1]-4*s0[2]/n0[2]+s0[3]/n0[3]-4*s0[4]/n0[4]-4*s0[5]/n0[5]+16*s0[6]/n0[6]-4*s0[7]/n0[7]-4*s0[8]/n0[8]+s0[9]/n0[9]+s0[10]/n0[10]
    aa4<-sigma*(1/n0[1]+16/n0[2]+1/n0[3]+1/n0[9]+1/n0[10])
    aa1<-1000;n_iter<-0
    while (aa1>0.0001)
    {
      n_iter<-n_iter+1
      gs[1]<-0.25*(mean[1]-mean[3]+mean3[3]-mean3[7])       	 #da.
      gs[2]<-0.25*(mean[1]-mean[3]-mean3[3]+mean3[7])  	       #db.
      gs[3]<-(-8*mean[1]-2*mean[2]-8*mean[3]-2*mean1[2]+15*mean1[3]+
                8*mean1[4]+15*mean2[2]-2*mean2[3]-8*mean3[3]-8*mean3[7])/17
      # ha.
      gs[4]<-(-8*mean[1]-2*mean[2]-8*mean[3]+15*mean1[2]-2*mean1[3]+
                8*mean1[4]-2*mean2[2]+15*mean2[3]-8*mean3[3]-8*mean3[7])/17
      # hb.
      gs[5]<-0.25*(mean[1]+mean[3]-mean3[3]-mean3[7])  #  i.
      gs[6]<--0.5*mean[1]+0.5*mean[3]+mean1[2]-mean2[3]-0.5*mean3[3]+0.5*mean3[7]
      #  jab.
      gs[7]<--0.5*mean[1]+0.5*mean[3]+mean1[3]-mean2[2]+0.5*mean3[3]-0.5*mean3[7]
      #  jba.
      gs[8]<-(12*mean[1]+20*mean[2]+12*mean[3]-14*mean1[2]-14*mean1[3]-
                12*mean1[4]-14*mean2[2]-14*mean2[3]+12*mean3[3]+12*mean3[7])/17
      #  l.
      g_aa1<-(0.5*(gs[2]+gs[5])^2+0.25*(gs[4]+gs[6])^2)/n_fam
      #   0.5(db+i)**2+0.25(hb+jab)**2.
      g_aa2<-(0.5*(gs[1]+gs[5])^2+0.25*(gs[3]+gs[7])^2)/n_fam
      #   0.5(da+i)**2+0.25(ha+jba)**2.
      g_aa3<-(0.5*(gs[1]-gs[5])^2+0.25*(gs[3]-gs[7])^2)/n_fam
      #   0.5(da-i)**2+0.25(ha-jba)**2.
      g_aa4<-(0.5*(gs[2]-gs[5])^2+0.25*(gs[4]-gs[6])^2)/n_fam
      #   0.5(db-i)**2+0.25(hb-jab)**2.
      g_aa5<-0.25*(gs[1]^2+gs[2]^2+gs[5]^2+(gs[1]+gs[6])^2+(gs[2]+gs[7])^2+(gs[3]+gs[8]/2)^2+(gs[4]+gs[8]/2)^2+gs[8]^2/4)/n_fam
      sigma1[2]<-sigma+g_aa1;sigma1[3]<-sigma+g_aa2;sigma1[4]<-sigma+g_aa5
      sigma2[1]<-sigma1[4];sigma2[2]<-sigma+g_aa3;sigma2[3]<-sigma+g_aa4;sigma2[4]<-sigma
      sigma3[1]<-sigma3[3]<-sigma3[7]<-sigma3[9]<-sigma
      sigma3[2]<-sigma1[2];sigma3[4]<-sigma1[3];sigma3[5]<-sigma1[4];sigma3[6]<-sigma2[2];sigma3[8]<-sigma2[3]
      aa2<-aa4+16*sigma1[2]/n0[4]+16*sigma1[3]/n0[5]+256*sigma1[4]/n0[6]+16*sigma2[2]/n0[7]+16*sigma2[3]/n0[8]
      aaa1<-aa3/aa2         # coefficient in restricted condition.
      mean[1]<-(s0[1]-aaa1*sigma)/n0[1]
      mean[2]<-(s0[2]+4*aaa1*sigma)/n0[2]
      mean[3]<-(s0[3]-aaa1*sigma)/n0[3]
      mean1[2]<-(s0[4]+4*aaa1*sigma1[2])/n0[4]
      mean1[3]<-(s0[5]+4*aaa1*sigma1[3])/n0[5]
      mean1[4]<-(s0[6]-16*aaa1*sigma1[4])/n0[6]
      mean2[2]<-(s0[7]+4*aaa1*sigma2[2])/n0[7]
      mean2[3]<-(s0[8]+4*aaa1*sigma2[3])/n0[8]
      mean3[3]<-(s0[9]-aaa1*sigma)/n0[9]
      mean3[7]<-(s0[10]-aaa1*sigma)/n0[10]
      mean1[1]<-mean[1]
      mean2[1]<-mean1[4];mean2[4]<-mean[3]
      mean3[1]<-mean[1];mean3[2]<-mean1[2];mean3[4]<-mean1[3]
      mean3[5]<-mean1[4];mean3[6]<-mean2[2];mean3[8]<-mean2[3];mean3[9]<-mean[3]
      aa1<-abs(aaa1-aaa0)
      aaa0<-aaa1
      if(n_iter>20)break
    }
    #######obtain variance###############################
    ss1<-sum((dataP1-mean[1])^2);ss3<-sum((dataP2-mean[3])^2);ss2<-sum((dataF1-mean[2])^2)
    for(i in 1:4) {swx1[i] <- W1[i,]%*%(dataB1-mean1[i])^2} ;for(i in 1:4) {swx2[i] <- W2[i,]%*%(dataB2-mean2[i])^2 };for(i in 1:9) {swx3[i] <- W3[i,]%*%(dataF2-mean3[i])^2 }
    s0[11]<-ss1+ss2+ss3+swx1[1]+swx2[4]+swx3[1]+swx3[3]+swx3[7]+swx3[9]
    n0[11]<-n_samP1+n_samF1+n_samP2+mix_pi1[1]*n_samB1+mix_pi2[4]*n_samB2+(mix_pi3[1]+mix_pi3[3]+mix_pi3[7]+mix_pi3[9])*n_samF2
    s0[12]<-swx1[2]+swx3[2]
    s0[13]<-swx1[3]+swx3[4]
    s0[14]<-swx1[4]+swx2[1]+swx3[5]
    s0[15]<-swx2[2]+swx3[6]
    s0[16]<-swx2[3]+swx3[8]
    n0[12]<-mix_pi1[2]*n_samB1+mix_pi3[2]*n_samF2
    n0[13]<-mix_pi1[3]*n_samB1+mix_pi3[4]*n_samF2
    n0[14]<-mix_pi1[4]*n_samB1+mix_pi2[1]*n_samB2+mix_pi3[5]*n_samF2
    n0[15]<-mix_pi2[2]*n_samB2+mix_pi3[6]*n_samF2
    n0[16]<-mix_pi2[3]*n_samB2+mix_pi3[8]*n_samF2

    gs[1]<-0.25*(mean[1]-mean[3]+mean3[3]-mean3[7])       	 #da.
    gs[2]<-0.25*(mean[1]-mean[3]-mean3[3]+mean3[7])  	       #db.
    gs[3]<-(-8*mean[1]-2*mean[2]-8*mean[3]-2*mean1[2]+15*mean1[3]+
              8*mean1[4]+15*mean2[2]-2*mean2[3]-8*mean3[3]-8*mean3[7])/17
    # ha.
    gs[4]<-(-8*mean[1]-2*mean[2]-8*mean[3]+15*mean1[2]-2*mean1[3]+
              8*mean1[4]-2*mean2[2]+15*mean2[3]-8*mean3[3]-8*mean3[7])/17
    # hb.
    gs[5]<-0.25*(mean[1]+mean[3]-mean3[3]-mean3[7])  #  i.
    gs[6]<--0.5*mean[1]+0.5*mean[3]+mean1[2]-mean2[3]-0.5*mean3[3]+0.5*mean3[7]
    #  jab.
    gs[7]<--0.5*mean[1]+0.5*mean[3]+mean1[3]-mean2[2]+0.5*mean3[3]-0.5*mean3[7]
    #  jba.
    gs[8]<-(12*mean[1]+20*mean[2]+12*mean[3]-14*mean1[2]-14*mean1[3]-
              12*mean1[4]-14*mean2[2]-14*mean2[3]+12*mean3[3]+12*mean3[7])/17
    #  l.
    g_aa1<-(0.5*(gs[2]+gs[5])^2+0.25*(gs[4]+gs[6])^2)/n_fam
    #   0.5(db+i)**2+0.25(hb+jab)**2.
    g_aa2<-(0.5*(gs[1]+gs[5])^2+0.25*(gs[3]+gs[7])^2)/n_fam
    #   0.5(da+i)**2+0.25(ha+jba)**2.
    g_aa3<-(0.5*(gs[1]-gs[5])^2+0.25*(gs[3]-gs[7])^2)/n_fam
    #   0.5(da-i)**2+0.25(ha-jba)**2.
    g_aa4<-(0.5*(gs[2]-gs[5])^2+0.25*(gs[4]-gs[6])^2)/n_fam
    #   0.5(db-i)**2+0.25(hb-jab)**2.
    g_aa5<-0.25*(gs[1]^2+gs[2]^2+gs[5]^2+(gs[1]+gs[6])^2+(gs[2]+gs[7])^2+(gs[3]+gs[8]/2)^2+(gs[4]+gs[8]/2)^2+gs[8]^2/4)/n_fam
    aaa0<-sigma;aa3<-1000;n_iter<-0
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma/(sigma+g_aa1);aa2<-sigma/(sigma+g_aa2)
      aa3<-sigma/(sigma+g_aa5);aa4<-sigma/(sigma+g_aa3)
      aa5<-sigma/(sigma+g_aa4)
      sigma<-(s0[11]+aa1^2*s0[12]+aa2^2*s0[13]+aa3^2*s0[14]+aa4^2*s0[15]+aa5^2*s0[16])/(n0[11]+aa1*n0[12]+aa2*n0[13]+aa3*n0[14]+aa4*n0[15]+aa5*n0[16])
      aa3<-abs(sigma-aaa0)
      aaa0<-sigma
      if(n_iter>20)break
    }
    sigma1[1]<-sigma;sigma1[2]<-sigma+g_aa1;sigma1[3]<-sigma+g_aa2;sigma1[4]<-sigma+g_aa5
    sigma2[1]<-sigma1[4];sigma2[2]<-sigma+g_aa3;sigma2[3]<-sigma+g_aa4;sigma2[4]<-sigma
    sigma3[1]<-sigma3[3]<-sigma3[7]<-sigma3[9]<-sigma
    sigma3[2]<-sigma1[2];sigma3[4]<-sigma1[3];sigma3[5]<-sigma1[4]
    sigma3[6]<-sigma2[2];sigma3[8]<-sigma2[3]
    if(sum(sigma < 1e-30)>=1){break}
    ##############################criteria for iterations to stop########################
    L1<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,mean[3],sqrt(sigma))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mix_pi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mix_pi3)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*10
  meanP1<-mean[1];meanF1<-mean[2];meanP2<-mean[3]
  sigma0<-sigma
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*n_samP1)
  P1bmw <- matrix(0,n_samP1,1)
  P1gg <- (dataP1 - mean[1])/sqrt(as.vector(sigma))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < n_samP1){P1bmw <- P1bmw+runif(n_samP1)/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  P1u<- as.matrix(c(12*n_samP1*((P1dd[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((P1dd[2]/n_samP1-1/3)^2),180*n_samP1*((P1dd[3]/n_samP1-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[2])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),P1D))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*n_samF1)
  F1bmw <- matrix(0,n_samF1,1)
  F1gg <- (dataF1 - mean[2])/sqrt(as.vector(sigma))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < n_samF1){F1bmw <- F1bmw+runif(n_samF1)/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  F1u<- as.matrix(c(12*n_samF1*((F1dd[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((F1dd[2]/n_samF1-1/3)^2),180*n_samF1*((F1dd[3]/n_samF1-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[2])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),F1D))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*n_samP2)
  P2bmw <- matrix(0,n_samP2,1)
  P2gg <- (dataP2 - mean[3])/sqrt(as.vector(sigma))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < n_samP2){P2bmw <- P2bmw+runif(n_samP2)/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  P2u<- as.matrix(c(12*n_samP2*((P2dd[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((P2dd[2]/n_samP2-1/3)^2),180*n_samP2*((P2dd[3]/n_samP2-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[2])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),P2D))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1);
  B1w1<-1/(12*n_samB1)
  B1bmw <- matrix(0,n_samB1,1); B1bmwsl <- matrix(0,n_samB1,4)
  for(i in 1:4){
    B1gg <- (dataB1 - mean1[i])/sqrt(sigma1[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi1[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B1P2)))[1]
  if(nn < n_samB1){B1P2 <- B1P2+runif(n_samB1)/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*n_samB1) + sum((B1P2 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  B1u <- as.matrix(c(12*n_samB1*((B1dd[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((B1dd[2]/n_samB1-1/3)^2),180*n_samB1*((B1dd[3]/n_samB1-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[2])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),B1D))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])

  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2);
  B2w1<-1/(12*n_samB2)
  B2bmw <- matrix(0,n_samB2,1); B2bmwsl <- matrix(0,n_samB2,4)
  for(i in 1:4){
    B2gg <- (dataB2 - mean2[i])/sqrt(sigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B2P2)))[1]
  if(nn < n_samB2){B2P2 <- B2P2+runif(n_samB2)/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*n_samB2) + sum((B2P2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  B2u <- as.matrix(c(12*n_samB2*((B2dd[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((B2dd[2]/n_samB2-1/3)^2),180*n_samB2*((B2dd[3]/n_samB2-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[2])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),B2D))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])

  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2);
  F2w1<-1/(12*n_samF2)
  F2bmw <- matrix(0,n_samF2,1); F2bmwsl <- matrix(0,n_samF2,9)
  for(i in 1:9){
    F2gg <- (dataF2 - mean3[i])/sqrt(sigma3[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi3[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < n_samF2){F2P2 <- F2P2+runif(n_samF2)/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*n_samF2) + sum((F2P2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  F2u <- as.matrix(c(12*n_samF2*((F2dd[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((F2dd[2]/n_samF2-1/3)^2),180*n_samF2*((F2dd[3]/n_samF2-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[2])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),F2D))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])

  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)

  #########first order parameters###################################
  aa<-matrix(c(1,1,1,0,0,1,0,0,0,1,0,0,1,1,0,0,0,1,1,-1,-1,0,0,1,0,0,0,
               1,1,0,0,0.5,0,0.5,0,0,1,0,1,0.5,0,0,0,0.5,0,1,0,0,0.5,0.5,0,0,0,0.25,
               1,0,-1,0.5,0,0,0,-0.5,0,1,-1,0,0,0.5,0,-0.5,0,0,1,1,-1,0,0,-1,0,0,0,
               1,-1,1,0,0,-1,0,0,0),10,9,byrow=T)
  mm<-as.matrix(c(mean[1],mean[2],mean[3],mean1[2],mean1[3],mean1[4],mean2[2],mean2[3],mean3[3],mean3[7]))
  B<-solve(crossprod(aa,aa))%*%crossprod(aa,mm)
  #######second ordedr parameters##############################
  jj1<-sigmaB1-sigma1[1]
  if (jj1<0 ) {jj1<-0}
  ll1<-jj1/sigmaB1
  jj2<-sigmaB2-sigma2[4]
  if (jj2<0) {jj2<-0}
  ll2<-jj2/sigmaB2
  jj3<-sigmaF2-sigma3[1]
  if (jj3<0) {jj3<-0}
  ll3<-jj3/sigmaF2

  output <- data.frame("2MG-ADI",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4), round(sigma0,4),round(t(mean1),4),round(t(sigma1),4),
                       round(t(mix_pi1),4),round(t(mean2),4),round(t(sigma2),4),round(t(mix_pi2),4),
                       round(t(mean3),4),round(t(sigma3),4),round(t(mix_pi3),4),
                       round(B[1],4)," "," "," "," "," ",round(B[2],4),round(B[3],4),round(B[4],4),round(B[5],4),round(B[6],4),round(B[7],4),round(B[8],4),round(B[9],4)," "," ",
                       round(jj1,4),round(ll1*100,4)," "," ",round(jj2,4),round(ll2*100,4)," "," ",round(jj3,4),round(ll3*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi1,mi2,mi3)
  return(OUTPUT)
}

#################2MG-AD(B-2)######################################
G6FModelFun[[6]] <- function(K1,logL,df11,df21,df31,df41,df51,df61,G6Ftext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1];n_samF1<-dim(dataF1)[1];n_samP2<-dim(dataP2)[1]
  n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1];n_samF2<-dim(dataF2)[1]
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2),sum(dataF2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2),sum(dataF2^2)))
  ss<-matrix(0,3,1);ss[1]<-s[1]-sumx[1]^2/n_samP1;ss[2]<-s[2]-sumx[2]^2/n_samF1;ss[3]<-s[3]-sumx[3]^2/n_samP2
  mean<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataB1),mean(dataB2),mean(dataF2)))
  sigma0<-sum(ss)/(n_samP1+n_samF1+n_samP2-3)
  sigmaB1<-var(dataB1);sigma40<-sigmaB1
  sigmaB2<-var(dataB2);sigma50<-sigmaB2
  sigmaF2<-var(dataF2);sigma60<-sigmaF2
  sigmaP1<-sigmaF1<-sigmaP2<-sigma0
  m_esp <- 0.0001;n_fam <- as.numeric(G6Ftext2)
  ###############procedure start###########################
  mi1<-matrix(0.25,4,1);sigma1<-matrix(0,4,1)
  mi2<-matrix(0.25,4,1);sigma2<-matrix(0,4,1)
  mi3<-as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625));sigma3<-matrix(0,9,1)
  sigma<-sigma0
  a1<-sqrt(sigmaB1/n_samB1)
  if (mean[1]<mean[3]) a1<--a1
  mean1<-as.matrix(c(mean[4]+2.4*a1,mean[4]+0.8*a1,mean[4]-0.8*a1,mean[4]-2.4*a1))
  a2<-sqrt(sigmaB2/n_samB2)
  if (mean[1]<mean[3]) a2<--a2
  mean2<-as.matrix(c(mean[5]+2.4*a2,mean[5]+0.8*a2,mean[5]-0.8*a2,mean[5]-2.4*a2))
  a3<-sqrt(sigmaF2/n_samF2)
  if (mean[1]<mean[3]) a3<--a3
  mean3<-as.matrix(c(mean[6]+2.4*a3,mean[6]+1.8*a3,mean[6]+1.2*a3,mean[6]+0.6*a3,mean[6],mean[6]-0.6*a3,mean[6]-1.2*a3,mean[6]-1.8*a3,mean[6]-2.4*a3))
  gs<-matrix(0,4,1)
  gs[1]<-(mean[1]-mean[3]+mean1[2]-mean2[3]+mean3[3]-mean3[7])/6
  #da.
  gs[2]<-(mean[1]-mean[3]+mean1[3]-mean2[2]-mean3[3]+mean3[7])/6
  #db.
  gs[3]<--mean[1]/7+3*mean[2]/7-mean[3]/7-0.5*mean1[2]+0.5*mean1[3]+
    mean1[4]/7+0.5*mean2[2]-0.5*mean2[3]-mean3[3]/7-mean3[7]/7
  # ha.
  gs[4]<--mean[1]/7+3*mean[2]/7-mean[3]/7+0.5*mean1[2]-0.5*mean1[3]+
    mean1[4]/7-0.5*mean2[2]+0.5*mean2[3]-mean3[3]/7-mean3[7]/7
  # hb.

  g_aa1<-(0.5*gs[2]^2+0.25*gs[4]^2)/n_fam
  #   0.5*db**2+0.25*hb**2.
  g_aa2<-(0.5*gs[1]^2+0.25*gs[3]^2)/n_fam
  #   0.5*da**2+0.25*ha**2.
  g_aa3<-g_aa1+g_aa2
  #   0.5(da**2+db**2)+0.25(ha**2+hb**2).
  sigma1[1]<-sigma;sigma1[2]<-sigma+g_aa1;sigma1[3]<-sigma+g_aa2;sigma1[4]<-sigma+g_aa3
  sigma2[1]<-sigma1[4];sigma2[2]<-sigma1[3];sigma2[3]<-sigma1[2];sigma2[4]<-sigma
  sigma3[1]<-sigma3[3]<-sigma3[7]<-sigma3[9]<-sigma
  sigma3[2]<-sigma3[8]<-sigma1[2]
  sigma3[4]<-sigma3[6]<-sigma1[3]
  sigma3[5]<-sigma1[4]
  L0<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma0))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma0))))+
    sum(log(dnorm(dataP2,mean[3],sqrt(sigma0))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mi1)))+
    sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mi3)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,4,n_samB1);      swx1 <- matrix(0,4,1)
  W2 <- matrix(0,4,n_samB2);      swx2 <- matrix(0,4,1)
  W3 <- matrix(0,9,n_samF2);      swx3 <- matrix(0,9,1)
  s0<-matrix(0,14,1);n0<-matrix(0,14,1)
  hh<-matrix(0,5,5);b_line<-matrix(0,5,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:4) { W1[i,] <- mi1[i]*dnorm(dataB1,mean1[i],sqrt(sigma1[i]))/dmixnorm(dataB1,mean1,sqrt(sigma1),mi1) }
    mix_pi1 <- as.matrix(rowSums(W1)/n_samB1)
    sumwx1 <- W1%*%dataB1
    for(i in 1:4) { W2[i,] <- mi2[i]*dnorm(dataB2,mean2[i],sqrt(sigma2[i]))/dmixnorm(dataB2,mean2,sqrt(sigma2),mi2) }
    mix_pi2 <- as.matrix(rowSums(W2)/n_samB2)
    sumwx2 <- W2%*%dataB2
    for(i in 1:9) { W3[i,] <- mi3[i]*dnorm(dataF2,mean3[i],sqrt(sigma3[i]))/dmixnorm(dataF2,mean3,sqrt(sigma3),mi3) }
    mix_pi3 <- as.matrix(rowSums(W3)/n_samF2)
    sumwx3 <- W3%*%dataF2
    ######obtain means########################
    aaa0<-0
    s0[1]<-sumx[1]+sumwx1[1]+sumwx3[1];s0[2]<-sumx[2]
    s0[3]<-sumx[3]+sumwx2[4]+sumwx3[9];s0[4]<-sumwx1[2]+sumwx3[2]
    s0[5]<-sumwx1[3]+sumwx3[4];s0[6]<-sumwx1[4]+sumwx2[1]+sumwx3[5]
    s0[7]<-sumwx2[2]+sumwx3[6];s0[8]<-sumwx2[3]+sumwx3[8]
    s0[9]<-sumwx3[3];s0[10]<-sumwx3[7]
    n0[1]<-n_samP1+mix_pi1[1]*n_samB1+mix_pi3[1]*n_samF2;n0[2]<-n_samF1
    n0[3]<-n_samP2+mix_pi2[4]*n_samB2+mix_pi3[9]*n_samF2;n0[4]<-mix_pi1[2]*n_samB1+mix_pi3[2]*n_samF2
    n0[5]<-mix_pi1[3]*n_samB1+mix_pi3[4]*n_samF2;n0[6]<-mix_pi1[4]*n_samB1+mix_pi2[1]*n_samB2+mix_pi3[5]*n_samF2
    n0[7]<-mix_pi2[2]*n_samB2+mix_pi3[6]*n_samF2;n0[8]<-mix_pi2[3]*n_samB2+mix_pi3[8]*n_samF2
    n0[9]<-mix_pi3[3]*n_samF2;n0[10]<-mix_pi3[7]*n_samF2
    s0[c(1:10)][abs(s0[c(1:10)])<0.000001]<-0.000001;n0[c(1:10)][abs(n0[c(1:10)])<0.000001]<-0.000001
    n_iter<-0;aaa1<-1000;AA<-matrix(0,5,1)
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      gs[1]<-(mean[1]-mean[3]+mean1[2]-mean2[3]+mean3[3]-mean3[7])/6
      #da.
      gs[2]<-(mean[1]-mean[3]+mean1[3]-mean2[2]-mean3[3]+mean3[7])/6
      #db.
      gs[3]<--mean[1]/7+3*mean[2]/7-mean[3]/7-0.5*mean1[2]+0.5*mean1[3]+
        mean1[4]/7+0.5*mean2[2]-0.5*mean2[3]-mean3[3]/7-mean3[7]/7
      # ha.
      gs[4]<--mean[1]/7+3*mean[2]/7-mean[3]/7+0.5*mean1[2]-0.5*mean1[3]+
        mean1[4]/7-0.5*mean2[2]+0.5*mean2[3]-mean3[3]/7-mean3[7]/7
      # hb.
      g_aa1<-(0.5*gs[2]^2+0.25*gs[4]^2)/n_fam
      #   0.5*db**2+0.25*hb**2.
      g_aa2<-(0.5*gs[1]^2+0.25*gs[3]^2)/n_fam
      #   0.5*da**2+0.25*ha**2.
      g_aa3<-g_aa1+g_aa2
      #   0.5(da**2+db**2)+0.25(ha**2+hb**2).
      sigma1[2]<-sigma+g_aa1;sigma1[3]<-sigma+g_aa2;sigma1[4]<-sigma+g_aa3
      sigma2[1]<-sigma1[4];sigma2[2]<-sigma1[3];sigma2[3]<-sigma1[2];sigma2[4]<-sigma
      sigma3[1]<-sigma3[3]<-sigma3[7]<-sigma3[9]<-sigma
      sigma3[2]<-sigma3[8]<-sigma1[2]
      sigma3[4]<-sigma1[3];sigma3[5]<-sigma1[4];sigma3[6]<-sigma1[3]
      #############################################################
      hh[1,1]<-sigma*(1/n0[1]+1/n0[3]+1/n0[9]+1/n0[10])
      hh[1,2]<-sigma*(1/n0[1]-1/n0[3]-1/n0[9]+1/n0[10])
      hh[1,3]<-sigma*(1/n0[1]-1/n0[3]+1/n0[9]-1/n0[10])
      hh[1,4]<-sigma*(5/n0[1]+19/n0[3])
      hh[1,5]<-sigma*(1/n0[1]+1/n0[3])
      hh[2,2]<-sigma*(1/n0[1]+1/n0[3]+1/n0[9]+1/n0[10])+sigma1[2]*(4/n0[4]+4/n0[8])
      hh[2,3]<-sigma*(1/n0[1]+1/n0[3]-1/n0[9]-1/n0[10])
      hh[2,4]<-sigma*(5/n0[1]-19/n0[3])-28*sigma1[2]/n0[8]
      hh[2,5]<-sigma*(1/n0[1]-1/n0[3])
      hh[3,3]<-sigma*(1/n0[1]+1/n0[3]+1/n0[9]+1/n0[10])+sigma1[3]*(4/n0[5]+4/n0[7])
      hh[3,4]<-sigma*(5/n0[1]-19/n0[3])-28*sigma1[3]/n0[7]
      hh[3,5]<-sigma*(1/n0[1]-1/n0[3])
      hh[4,4]<-sigma*(25/n0[1]+100/n0[2]+361/n0[3])+196*sigma1[3]/n0[7]+196*sigma1[2]/n0[8]+36*sigma1[4]/n0[6]
      hh[4,5]<-sigma*(5/n0[1]+20/n0[2]+19/n0[3])+24*sigma1[4]/n0[6]
      hh[5,5]<-sigma*(1/n0[1]+4/n0[2]+1/n0[3])+16*sigma1[4]/n0[6]
      for(i in 2:5)
      {
        for(j in 1:(i-1))
        {
          hh[i,j]<-hh[j,i]
        }
      }
      ################################################################
      b_line[1]<-s0[1]/n0[1]+s0[3]/n0[3]-s0[9]/n0[9]-s0[10]/n0[10]
      b_line[2]<-s0[1]/n0[1]-s0[3]/n0[3]-2*s0[4]/n0[4]+2*s0[8]/n0[8]+s0[9]/n0[9]-s0[10]/n0[10]
      b_line[3]<-s0[1]/n0[1]-s0[3]/n0[3]-2*s0[5]/n0[5]+2*s0[7]/n0[7]-s0[9]/n0[9]+s0[10]/n0[10]
      b_line[4]<-5*s0[1]/n0[1]+10*s0[2]/n0[2]+19*s0[3]/n0[3]-14*s0[7]/n0[7]-14*s0[8]/n0[8]-6*s0[6]/n0[6]
      b_line[5]<-s0[1]/n0[1]+2*s0[2]/n0[2]+s0[3]/n0[3]-4*s0[6]/n0[6]
      B<-solve(hh,b_line)
      mean[1]<-(s0[1]-sigma*(B[1]+B[2]+B[3]+5*B[4]+B[5]))/n0[1]
      mean[2]<-(s0[2]-sigma*(10*B[4]+2*B[5]))/n0[2]
      mean[3]<-(s0[3]-sigma*(B[1]-B[2]-B[3]+19*B[4]+B[5]))/n0[3]
      mean1[2]<-(s0[4]+2*B[2]*sigma1[2])/n0[4]
      mean1[3]<-(s0[5]+2*B[3]*sigma1[3])/n0[5]
      mean1[4]<-(s0[6]+sigma1[4]*(6*B[4]+4*B[5]))/n0[6]
      mean2[2]<-(s0[7]-sigma1[3]*(2*B[3]-14*B[4]))/n0[7]
      mean2[3]<-(s0[8]-sigma1[2]*(2*B[2]-14*B[4]))/n0[8]
      mean3[3]<-(s0[9]+sigma*(B[1]-B[2]+B[3]))/n0[9]
      mean3[7]<-(s0[10]+sigma*(B[1]+B[2]-B[3]))/n0[10]
      mean1[1]<-mean[1];mean2[1]<-mean1[4];mean2[4]<-mean[3]
      mean3[1]<-mean[1];mean3[2]<-mean1[2];mean3[4]<-mean1[3]
      mean3[5]<-mean1[4];mean3[6]<-mean2[2];mean3[8]<-mean2[3];mean3[9]<-mean[3]
      aaa1<-max(abs(B-AA))
      AA<-B
      if(n_iter>20)break

    }
    ############obtain variance#######################
    ss1<-sum((dataP1-mean[1])^2);ss3<-sum((dataP2-mean[3])^2);ss2<-sum((dataF1-mean[2])^2)
    for(i in 1:4) {swx1[i] <- W1[i,]%*%(dataB1-mean1[i])^2 };for(i in 1:4) {swx2[i] <- W2[i,]%*%(dataB2-mean2[i])^2 };for(i in 1:9) {swx3[i] <- W3[i,]%*%(dataF2-mean3[i])^2 }
    s0[11]<-ss1+ss2+ss3+swx1[1]+swx2[4]+swx3[1]+swx3[3]+swx3[7]+swx3[9]
    n0[11]<-n_samP1+n_samF1+n_samP2+mix_pi1[1]*n_samB1+mix_pi2[4]*n_samB2+(mix_pi3[1]+mix_pi3[3]+mix_pi3[7]+mix_pi3[9])*n_samF2
    s0[12]<-swx1[2]+swx2[3]+swx3[2]+swx3[8];s0[13]<-swx1[3]+swx2[2]+swx3[4]+swx3[6];s0[14]<-swx1[4]+swx2[1]+swx3[5]
    n0[12]<-mix_pi1[2]*n_samB1+mix_pi2[3]*n_samB2+mix_pi3[2]*n_samF2+mix_pi3[8]*n_samF2
    n0[13]<-mix_pi1[3]*n_samB1+mix_pi2[2]*n_samB2+mix_pi3[4]*n_samF2+mix_pi3[6]*n_samF2
    n0[14]<-mix_pi1[4]*n_samB1+mix_pi2[1]*n_samB2+mix_pi3[5]*n_samF2

    gs[1]<-(mean[1]-mean[3]+mean1[2]-mean2[3]+mean3[3]-mean3[7])/6
    #da.
    gs[2]<-(mean[1]-mean[3]+mean1[3]-mean2[2]-mean3[3]+mean3[7])/6
    #db.
    gs[3]<--mean[1]/7+3*mean[2]/7-mean[3]/7-0.5*mean1[2]+0.5*mean1[3]+
      mean1[4]/7+0.5*mean2[2]-0.5*mean2[3]-mean3[3]/7-mean3[7]/7
    # ha.
    gs[4]<--mean[1]/7+3*mean[2]/7-mean[3]/7+0.5*mean1[2]-0.5*mean1[3]+
      mean1[4]/7-0.5*mean2[2]+0.5*mean2[3]-mean3[3]/7-mean3[7]/7
    # hb.
    g_aa1<-(0.5*gs[2]^2+0.25*gs[4]^2)/n_fam
    #   0.5*db**2+0.25*hb**2.
    g_aa2<-(0.5*gs[1]^2+0.25*gs[3]^2)/n_fam
    #   0.5*da**2+0.25*ha**2.
    g_aa3<-g_aa1+g_aa2
    #   0.5(da**2+db**2)+0.25(ha**2+hb**2).
    aaa0<-sigma;aa3<-1000;n_iter<-0
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma/(sigma+g_aa1)
      aa2<-sigma/(sigma+g_aa2)
      aa3<-sigma/(sigma+g_aa3)
      sigma<-(s0[11]+aa1^2*s0[12]+aa2^2*s0[13]+aa3^2*s0[14])/(n0[11]+aa1*n0[12]+aa2*n0[13]+aa3*n0[14])
      aa3<-abs(sigma-aaa0)
      aaa0<-sigma
      if(n_iter>20)break
    }
    sigma1[1]<-sigma;sigma1[2]<-sigma+g_aa1;sigma1[3]<-sigma+g_aa2;sigma1[4]<-sigma+g_aa3
    sigma2[1]<-sigma1[4];sigma2[2]<-sigma1[3];sigma2[3]<-sigma1[2];sigma2[4]<-sigma
    sigma3[1]<-sigma3[3]<-sigma3[7]<-sigma3[9]<-sigma
    sigma3[2]<-sigma3[8]<-sigma1[2]
    sigma3[4]<-sigma3[6]<-sigma1[3]
    sigma3[5]<-sigma1[4]
    if(sum(sigma < 1e-30)>=1){break}
    ##############################criteria for iterations to stop########################
    L1<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,mean[3],sqrt(sigma))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mix_pi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mix_pi3)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*6
  meanP1<-mean[1];meanF1<-mean[2];meanP2<-mean[3]
  sigma0<-sigma
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*n_samP1)
  P1bmw <- matrix(0,n_samP1,1)
  P1gg <- (dataP1 - mean[1])/sqrt(as.vector(sigma))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < n_samP1){P1bmw <- P1bmw+runif(n_samP1)/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  P1u<- as.matrix(c(12*n_samP1*((P1dd[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((P1dd[2]/n_samP1-1/3)^2),180*n_samP1*((P1dd[3]/n_samP1-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[2])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),P1D))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*n_samF1)
  F1bmw <- matrix(0,n_samF1,1)
  F1gg <- (dataF1 - mean[2])/sqrt(as.vector(sigma))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < n_samF1){F1bmw <- F1bmw+runif(n_samF1)/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  F1u<- as.matrix(c(12*n_samF1*((F1dd[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((F1dd[2]/n_samF1-1/3)^2),180*n_samF1*((F1dd[3]/n_samF1-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[2])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),F1D))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*n_samP2)
  P2bmw <- matrix(0,n_samP2,1)
  P2gg <- (dataP2 - mean[3])/sqrt(as.vector(sigma))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < n_samP2){P2bmw <- P2bmw+runif(n_samP2)/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  P2u<- as.matrix(c(12*n_samP2*((P2dd[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((P2dd[2]/n_samP2-1/3)^2),180*n_samP2*((P2dd[3]/n_samP2-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[2])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),P2D))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1);
  B1w1<-1/(12*n_samB1)
  B1bmw <- matrix(0,n_samB1,1); B1bmwsl <- matrix(0,n_samB1,4)
  for(i in 1:4){
    B1gg <- (dataB1 - mean1[i])/sqrt(sigma1[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi1[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B1P2)))[1]
  if(nn < n_samB1){B1P2 <- B1P2+runif(n_samB1)/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*n_samB1) + sum((B1P2 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  B1u <- as.matrix(c(12*n_samB1*((B1dd[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((B1dd[2]/n_samB1-1/3)^2),180*n_samB1*((B1dd[3]/n_samB1-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[2])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),B1D))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])

  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2);
  B2w1<-1/(12*n_samB2)
  B2bmw <- matrix(0,n_samB2,1); B2bmwsl <- matrix(0,n_samB2,4)
  for(i in 1:4){
    B2gg <- (dataB2 - mean2[i])/sqrt(sigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B2P2)))[1]
  if(nn < n_samB2){B2P2 <- B2P2+runif(n_samB2)/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*n_samB2) + sum((B2P2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  B2u <- as.matrix(c(12*n_samB2*((B2dd[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((B2dd[2]/n_samB2-1/3)^2),180*n_samB2*((B2dd[3]/n_samB2-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[2])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),B2D))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])

  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2);
  F2w1<-1/(12*n_samF2)
  F2bmw <- matrix(0,n_samF2,1); F2bmwsl <- matrix(0,n_samF2,9)
  for(i in 1:9){
    F2gg <- (dataF2 - mean3[i])/sqrt(sigma3[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi3[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < n_samF2){F2P2 <- F2P2+runif(n_samF2)/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*n_samF2) + sum((F2P2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  F2u <- as.matrix(c(12*n_samF2*((F2dd[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((F2dd[2]/n_samF2-1/3)^2),180*n_samF2*((F2dd[3]/n_samF2-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[2])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),F2D))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])

  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)

  #########first order parameters###################
  aa<-matrix(c(1,1,1,0,0,1,0,0,1,1,1,-1,-1,0,0,1,1,0,0,0.5,1,0,1,0.5,0,1,0,0,0.5,0.5,
               1,0,-1,0.5,0,1,-1,0,0,0.5,1,1,-1,0,0,1,-1,1,0,0),10,5,byrow=T)
  mm<-as.matrix(c(mean[1],mean[2],mean[3],mean1[2],mean1[3],mean1[4],mean2[2],mean2[3],mean3[3],mean3[7]))
  B1<-solve(crossprod(aa,aa))%*%crossprod(aa,mm)
  #########second order parameters####################
  jj1<-sigmaB1-sigma1[1]
  if (jj1<0 ) {jj1<-0}
  ll1<-jj1/sigmaB1
  jj2<-sigmaB2-sigma2[4]
  if (jj2<0) {jj2<-0}
  ll2<-jj2/sigmaB2
  jj3<-sigmaF2-sigma3[1]
  if (jj3<0) {jj3<-0}
  ll3<-jj3/sigmaF2


  output <- data.frame("2MG-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4), round(sigma0,4),round(t(mean1),4),round(t(sigma1),4),
                       round(t(mix_pi1),4),round(t(mean2),4),round(t(sigma2),4),round(t(mix_pi2),4),
                       round(t(mean3),4),round(t(sigma3),4),round(t(mix_pi3),4),
                       round(B1[1],4)," "," "," "," "," ",round(B1[2],4),round(B1[3],4),round(B1[4],4),round(B1[5],4)," "," "," "," "," "," ",
                       round(jj1,4),round(ll1*100,4)," "," ",round(jj2,4),round(ll2*100,4)," "," ",round(jj3,4),round(ll3*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi1,mi2,mi3)
  return(OUTPUT)
}

###################2MG-A(B-3)#############################
G6FModelFun[[7]] <- function(K1,logL,df11,df21,df31,df41,df51,df61,G6Ftext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1];n_samF1<-dim(dataF1)[1];n_samP2<-dim(dataP2)[1]
  n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1];n_samF2<-dim(dataF2)[1]
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2),sum(dataF2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2),sum(dataF2^2)))
  ss<-matrix(0,3,1);ss[1]<-s[1]-sumx[1]^2/n_samP1;ss[2]<-s[2]-sumx[2]^2/n_samF1;ss[3]<-s[3]-sumx[3]^2/n_samP2
  mean<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataB1),mean(dataB2),mean(dataF2)))
  sigma0<-sum(ss)/(n_samP1+n_samF1+n_samP2-3)
  sigmaB1<-var(dataB1);sigma40<-sigmaB1
  sigmaB2<-var(dataB2);sigma50<-sigmaB2
  sigmaF2<-var(dataF2);sigma60<-sigmaF2
  sigmaP1<-sigmaF1<-sigmaP2<-sigma0
  m_esp <- 0.0001;n_fam <- as.numeric(G6Ftext2)
  ###############procedure start###########################
  mi1<-matrix(0.25,4,1);sigma1<-matrix(0,4,1)
  mi2<-matrix(0.25,4,1);sigma2<-matrix(0,4,1)
  mi3<-as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  sigma3<-matrix(0,9,1)
  sigma<-sigma0
  a1<-sqrt(sigmaB1/n_samB1)
  if (mean[1]<mean[3]) a1<--a1
  mean1<-as.matrix(c(mean[4]+2.4*a1,mean[4]+0.8*a1,mean[4]-0.8*a1,mean[4]-2.4*a1))
  a2<-sqrt(sigmaB2/n_samB2)
  if (mean[1]<mean[3]) a2<--a2
  mean2<-as.matrix(c(mean[5]+2.4*a2,mean[5]+0.8*a2,mean[5]-0.8*a2,mean[5]-2.4*a2))
  a3<-sqrt(sigmaF2/n_samF2)
  if (mean[1]<mean[3]) a3<--a3
  mean3<-as.matrix(c(mean[6]+2.4*a3,mean[6]+1.8*a3,mean[6]+1.2*a3,mean[6]+0.6*a3,mean[6],mean[6]-0.6*a3,mean[6]-1.2*a3,mean[6]-1.8*a3,mean[6]-2.4*a3))
  gs<-matrix(0,2,1)
  gs[1]<-(mean[1]-mean[3]+mean1[2]-mean2[3]+mean3[3]-mean3[7])/6
  #da.
  gs[2]<-(mean[1]-mean[3]+mean1[3]-mean2[2]-mean3[3]+mean3[7])/6
  #db.
  sigma1[1]<-sigma
  g_aa1<-0.5*gs[2]^2/n_fam #   0.5*db**2.
  g_aa2<-0.5*gs[1]^2/n_fam #   0.5*da**2.
  g_aa3<-g_aa1+g_aa2#   0.5(da**2+db**2).
  sigma1[2]<-sigma+g_aa1;sigma1[3]<-sigma+g_aa2;sigma1[4]<-sigma+g_aa3
  sigma2[1]<-sigma1[4];sigma2[2]<-sigma1[3];sigma2[3]<-sigma1[2];sigma2[4]<-sigma
  sigma3[1]<-sigma3[3]<-sigma3[7]<-sigma3[9]<-sigma
  sigma3[2]<-sigma3[8]<-sigma1[2]
  sigma3[4]<-sigma3[6]<-sigma1[3]
  sigma3[5]<-sigma1[4]
  L0<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma0))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma0))))+
    sum(log(dnorm(dataP2,mean[3],sqrt(sigma0))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mi1)))+
    sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mi3)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,4,n_samB1);      swx1 <- matrix(0,4,1)
  W2 <- matrix(0,4,n_samB2);      swx2 <- matrix(0,4,1)
  W3 <- matrix(0,9,n_samF2);      swx3 <- matrix(0,9,1)
  s0<-matrix(0,14,1);n0<-matrix(0,14,1)
  hh<-matrix(0,6,6);b_line<-matrix(0,6,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:4) { W1[i,] <- mi1[i]*dnorm(dataB1,mean1[i],sqrt(sigma1[i]))/dmixnorm(dataB1,mean1,sqrt(sigma1),mi1) }
    mix_pi1 <- as.matrix(rowSums(W1)/n_samB1)
    sumwx1 <- W1%*%dataB1
    for(i in 1:4) { W2[i,] <- mi2[i]*dnorm(dataB2,mean2[i],sqrt(sigma2[i]))/dmixnorm(dataB2,mean2,sqrt(sigma2),mi2) }
    mix_pi2 <- as.matrix(rowSums(W2)/n_samB2)
    sumwx2 <- W2%*%dataB2
    for(i in 1:9) { W3[i,] <- mi3[i]*dnorm(dataF2,mean3[i],sqrt(sigma3[i]))/dmixnorm(dataF2,mean3,sqrt(sigma3),mi3) }
    mix_pi3 <- as.matrix(rowSums(W3)/n_samF2)
    sumwx3 <- W3%*%dataF2
    ##########obtain means###########################
    aaa0<-0
    s0[1]<-sumx[1]+sumwx1[1]+sumwx3[1];s0[2]<-sigma1[4]*sumx[2]+sigma*(sumwx1[4]+sumwx2[1]+sumwx3[5])
    s0[3]<-sumx[3]+sumwx2[4]+sumwx3[9];s0[4]<-sumwx1[2]+sumwx3[2]
    s0[5]<-sumwx1[3]+sumwx3[4];s0[7]<-sumwx2[2]+sumwx3[6]
    s0[8]<-sumwx2[3]+sumwx3[8];s0[9]<-sumwx3[3];s0[10]<-sumwx3[7]
    n0[1]<-n_samP1+mix_pi1[1]*n_samB1+mix_pi3[1]*n_samF2
    n0[2]<-n_samF1*sigma1[4]+sigma*(mix_pi1[4]*n_samB1+mix_pi2[1]*n_samB2+mix_pi3[5]*n_samF2)
    n0[3]<-n_samP2+mix_pi2[4]*n_samB2+mix_pi3[9]*n_samF2
    n0[4]<-mix_pi1[2]*n_samB1+mix_pi3[2]*n_samF2;n0[5]<-mix_pi1[3]*n_samB1+mix_pi3[4]*n_samF2
    n0[7]<-mix_pi2[2]*n_samB2+mix_pi3[6]*n_samF2;n0[8]<-mix_pi2[3]*n_samB2+mix_pi3[8]*n_samF2
    n0[9]<-mix_pi3[3]*n_samF2;n0[10]<-mix_pi3[7]*n_samF2
    s0[c(1:10)][abs(s0[c(1:10)])<0.000001]<-0.000001;n0[c(1:10)][abs(n0[c(1:10)])<0.000001]<-0.000001
    AA<-matrix(0,6,1);aaa1<-1000;n_iter<-0
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      gs[1]<-(mean[1]-mean[3]+mean1[2]-mean2[3]+mean3[3]-mean3[7])/6
      #da.
      gs[2]<-(mean[1]-mean[3]+mean1[3]-mean2[2]-mean3[3]+mean3[7])/6
      #db.
      sigma1[1]<-sigma
      g_aa1<-0.5*gs[2]^2/n_fam #   0.5*db**2.
      g_aa2<-0.5*gs[1]^2/n_fam #   0.5*da**2.
      g_aa3<-g_aa1+g_aa2#   0.5(da**2+db**2).
      sigma1[2]<-sigma+g_aa1;sigma1[3]<-sigma+g_aa2;sigma1[4]<-sigma+g_aa3
      sigma2[1]<-sigma1[4];sigma2[2]<-sigma1[3];sigma2[3]<-sigma1[2];sigma2[4]<-sigma
      sigma3[1]<-sigma3[3]<-sigma3[7]<-sigma3[9]<-sigma
      sigma3[2]<-sigma3[8]<-sigma1[2]
      sigma3[4]<-sigma3[6]<-sigma1[3]
      sigma3[5]<-sigma1[4]
      #############################################################
      hh[1,1]<-sigma*(1/n0[1]+1/n0[3]+1/n0[9]+1/n0[10])
      hh[1,2]<-sigma*(1/n0[1]-1/n0[3]-1/n0[9]+1/n0[10])
      hh[1,3]<-sigma*(1/n0[1]-1/n0[3]+1/n0[9]-1/n0[10])
      hh[1,4]<-0
      hh[1,5]<-0
      hh[1,6]<-sigma*(1/n0[9]+1/n0[10])
      hh[2,2]<-sigma*(1/n0[1]+1/n0[3]+1/n0[9]+1/n0[10])+sigma1[2]*(4/n0[4]+4/n0[8])
      hh[2,3]<-sigma*(1/n0[1]+1/n0[3]-1/n0[9]-1/n0[10])
      hh[2,4]<-2*sigma1[2]*(1/n0[4]-1/n0[8])
      hh[2,5]<-sigma1[2]*(-2/n0[4]+2/n0[8])
      hh[2,6]<-sigma*(-1/n0[9]+1/n0[10])
      hh[3,3]<-sigma*(1/n0[1]+1/n0[3]+1/n0[9]+1/n0[10])+sigma1[3]*(4/n0[5]+4/n0[7])
      hh[3,4]<-2*sigma1[3]*(1/n0[5]-1/n0[7])
      hh[3,5]<-2*sigma1[3]*(1/n0[5]-1/n0[7])
      hh[3,6]<-sigma*(1/n0[9]-1/n0[10])
      hh[4,4]<-16*sigma*sigma1[4]/n0[2]+sigma1[2]/n0[4]+sigma1[3]/n0[5]+sigma1[3]/n0[7]+sigma1[2]/n0[8]
      hh[4,5]<-sigma1[3]*(1/n0[5]+1/n0[7])-sigma1[2]*(1/n0[4]+1/n0[8])
      hh[4,6]<-8*sigma*sigma1[4]/n0[2]
      hh[5,5]<-sigma1[2]*(1/n0[4]+1/n0[8])+sigma1[3]*(1/n0[5]+1/n0[7])
      hh[5,6]<-0
      hh[6,6]<-sigma*(4*sigma1[4]/n0[2]+1/n0[9]+1/n0[10])
      for(i in 2:6)
      {
        for(j in 1:(i-1))
        {
          hh[i,j]<-hh[j,i]
        }
      }
      ########################################################
      b_line[1]<-s0[1]/n0[1]+s0[3]/n0[3]-s0[9]/n0[9]-s0[10]/n0[10]
      b_line[2]<-s0[1]/n0[1]-s0[3]/n0[3]-2*s0[4]/n0[4]+2*s0[8]/n0[8]+s0[9]/n0[9]-s0[10]/n0[10]
      b_line[3]<-s0[1]/n0[1]-s0[3]/n0[3]-2*s0[5]/n0[5]+2*s0[7]/n0[7]-s0[9]/n0[9]+s0[10]/n0[10]
      b_line[4]<-4*s0[2]/n0[2]-s0[4]/n0[4]-s0[5]/n0[5]-s0[7]/n0[7]-s0[8]/n0[8]
      b_line[5]<-s0[4]/n0[4]-s0[5]/n0[5]-s0[7]/n0[7]+s0[8]/n0[8]
      b_line[6]<-2*s0[2]/n0[2]-s0[9]/n0[9]-s0[10]/n0[10]
      B<-solve(hh,b_line)
      mean[1]<-(s0[1]-sigma*(B[1]+B[2]+B[3]))/n0[1]
      mean[2]<-(s0[2]-sigma*sigma1[4]*(4*B[4]+2*B[6]))/n0[2]
      mean[3]<-(s0[3]-sigma*(B[1]-B[2]-B[3]))/n0[3]
      mean1[2]<-(s0[4]+(2*B[2]+B[4]-B[5])*sigma1[2])/n0[4]
      mean1[3]<-(s0[5]+(2*B[3]+B[4]+B[5])*sigma1[3])/n0[5]
      mean2[2]<-(s0[7]-sigma1[3]*(2*B[3]-B[4]-B[5]))/n0[7]
      mean2[3]<-(s0[8]-sigma1[2]*(2*B[2]-B[4]+B[5]))/n0[8]
      mean3[3]<-(s0[9]+sigma*(B[1]-B[2]+B[3]+B[6]))/n0[9]
      mean3[7]<-(s0[10]+sigma*(B[1]+B[2]-B[3]+B[6]))/n0[10]
      mean1[1]<-mean[1];mean1[4]<-mean2[1]<-mean[2]
      mean2[4]<-mean[3]
      mean3[1]<-mean[1];mean3[2]<-mean1[2];mean3[4]<-mean1[3];mean3[5]<-mean[2]
      mean3[6]<-mean2[2];mean3[8]<-mean2[3];mean3[9]<-mean[3]
      aaa1<-max(abs(B-AA))
      AA<-B
      if(n_iter>20)break
    }
    #########obtain variance###############################
    ss1<-sum((dataP1-mean[1])^2);ss3<-sum((dataP2-mean[3])^2);ss2<-sum((dataF1-mean[2])^2)
    for(i in 1:4) {swx1[i] <- W1[i,]%*%(dataB1-mean1[i])^2 };for(i in 1:4) {swx2[i] <- W2[i,]%*%(dataB2-mean2[i])^2 };for(i in 1:9) {swx3[i] <- W3[i,]%*%(dataF2-mean3[i])^2 }
    s0[11]<-ss1+ss2+ss3+swx1[1]+swx2[4]+swx3[1]+swx3[3]+swx3[7]+swx3[9]
    n0[11]<-n_samP1+n_samF1+n_samP2+mix_pi1[1]*n_samB1+mix_pi2[4]*n_samB2+(mix_pi3[1]+mix_pi3[3]+mix_pi3[7]+mix_pi3[9])*n_samF2
    s0[12]<-swx1[2]+swx2[3]+swx3[2]+swx3[8]
    s0[13]<-swx1[3]+swx2[2]+swx3[4]+swx3[6]
    s0[14]<-swx1[4]+swx2[1]+swx3[5]
    n0[12]<-mix_pi1[2]*n_samB1+mix_pi2[3]*n_samB2+mix_pi3[2]*n_samF2+mix_pi3[8]*n_samF2
    n0[13]<-mix_pi1[3]*n_samB1+mix_pi2[2]*n_samB2+mix_pi3[4]*n_samF2+mix_pi3[6]*n_samF2
    n0[14]<-mix_pi1[4]*n_samB1+mix_pi2[1]*n_samB2+mix_pi3[5]*n_samF2
    gs[1]<-(mean[1]-mean[3]+mean1[2]-mean2[3]+mean3[3]-mean3[7])/6
    #da.
    gs[2]<-(mean[1]-mean[3]+mean1[3]-mean2[2]-mean3[3]+mean3[7])/6
    #db.
    g_aa1<-0.5*gs[2]^2/n_fam #   0.5*db**2.
    g_aa2<-0.5*gs[1]^2/n_fam #   0.5*da**2.
    g_aa3<-g_aa1+g_aa2#   0.5(da**2+db**2).
    aaa0<-sigma;aa4<-1000;n_iter<-0
    while (aa4>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma/(sigma+g_aa1)
      aa2<-sigma/(sigma+g_aa2)
      aa3<-sigma/(sigma+g_aa3)
      sigma<-(s0[11]+aa1*aa1*s0[12]+aa2*aa2*s0[13]+aa3*aa3*s0[14])/(n0[11]+aa1*n0[12]+aa2*n0[13]+aa3*n0[14])
      aa4<-abs(sigma-aaa0)
      aaa0<-sigma
      if(n_iter>20)break
    }
    sigma1[1]<-sigma;sigma1[2]<-sigma+g_aa1
    sigma1[3]<-sigma+g_aa2;sigma1[4]<-sigma+g_aa3
    sigma2[1]<-sigma1[4];sigma2[2]<-sigma1[3]
    sigma2[3]<-sigma1[2];sigma2[4]<-sigma
    sigma3[1]<-sigma3[3]<-sigma3[7]<-sigma3[9]<-sigma
    sigma3[2]<-sigma3[8]<-sigma1[2]
    sigma3[4]<-sigma3[6]<-sigma1[3]
    sigma3[5]<-sigma1[4]
    if(sum(sigma < 1e-30)>=1){break}
    ##############################criteria for iterations to stop########################
    L1<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,mean[3],sqrt(sigma))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mix_pi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mix_pi3)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*4
  meanP1<-mean[1];meanF1<-mean[2];meanP2<-mean[3]
  sigma0<-sigma
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*n_samP1)
  P1bmw <- matrix(0,n_samP1,1)
  P1gg <- (dataP1 - mean[1])/sqrt(as.vector(sigma))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < n_samP1){P1bmw <- P1bmw+runif(n_samP1)/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  P1u<- as.matrix(c(12*n_samP1*((P1dd[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((P1dd[2]/n_samP1-1/3)^2),180*n_samP1*((P1dd[3]/n_samP1-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[2])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),P1D))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*n_samF1)
  F1bmw <- matrix(0,n_samF1,1)
  F1gg <- (dataF1 - mean[2])/sqrt(as.vector(sigma))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < n_samF1){F1bmw <- F1bmw+runif(n_samF1)/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  F1u<- as.matrix(c(12*n_samF1*((F1dd[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((F1dd[2]/n_samF1-1/3)^2),180*n_samF1*((F1dd[3]/n_samF1-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[2])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),F1D))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*n_samP2)
  P2bmw <- matrix(0,n_samP2,1)
  P2gg <- (dataP2 - mean[3])/sqrt(as.vector(sigma))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < n_samP2){P2bmw <- P2bmw+runif(n_samP2)/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  P2u<- as.matrix(c(12*n_samP2*((P2dd[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((P2dd[2]/n_samP2-1/3)^2),180*n_samP2*((P2dd[3]/n_samP2-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[2])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),P2D))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1);
  B1w1<-1/(12*n_samB1)
  B1bmw <- matrix(0,n_samB1,1); B1bmwsl <- matrix(0,n_samB1,4)
  for(i in 1:4){
    B1gg <- (dataB1 - mean1[i])/sqrt(sigma1[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi1[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B1P2)))[1]
  if(nn < n_samB1){B1P2 <- B1P2+runif(n_samB1)/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*n_samB1) + sum((B1P2 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  B1u <- as.matrix(c(12*n_samB1*((B1dd[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((B1dd[2]/n_samB1-1/3)^2),180*n_samB1*((B1dd[3]/n_samB1-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[2])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),B1D))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])

  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2);
  B2w1<-1/(12*n_samB2)
  B2bmw <- matrix(0,n_samB2,1); B2bmwsl <- matrix(0,n_samB2,4)
  for(i in 1:4){
    B2gg <- (dataB2 - mean2[i])/sqrt(sigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B2P2)))[1]
  if(nn < n_samB2){B2P2 <- B2P2+runif(n_samB2)/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*n_samB2) + sum((B2P2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  B2u <- as.matrix(c(12*n_samB2*((B2dd[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((B2dd[2]/n_samB2-1/3)^2),180*n_samB2*((B2dd[3]/n_samB2-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[2])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),B2D))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])

  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2);
  F2w1<-1/(12*n_samF2)
  F2bmw <- matrix(0,n_samF2,1); F2bmwsl <- matrix(0,n_samF2,9)
  for(i in 1:9){
    F2gg <- (dataF2 - mean3[i])/sqrt(sigma3[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi3[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < n_samF2){F2P2 <- F2P2+runif(n_samF2)/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*n_samF2) + sum((F2P2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  F2u <- as.matrix(c(12*n_samF2*((F2dd[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((F2dd[2]/n_samF2-1/3)^2),180*n_samF2*((F2dd[3]/n_samF2-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[2])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),F2D))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])

  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)

  #########first order parameters########################
  aa<-matrix(c(1,1,1,1,0,0,1,-1,-1,1,1,0,1,0,1,
               1,0,-1,1,-1,0,1,1,-1,1,-1,1),9,3,byrow=T)
  mm<-as.matrix(c(mean[1],mean[2],mean[3],mean1[2],mean1[3],mean2[2],mean2[3],mean3[3],mean3[7]))
  B1<-solve(crossprod(aa,aa))%*%crossprod(aa,mm)
  ########second order parameters########################
  jj1<-sigmaB1-sigma1[1]
  if (jj1<0 ) {jj1<-0}
  ll1<-jj1/sigmaB1
  jj2<-sigmaB2-sigma2[4]
  if (jj2<0) {jj2<-0}
  ll2<-jj2/sigmaB2
  jj3<-sigmaF2-sigma3[1]
  if (jj3<0) {jj3<-0}
  ll3<-jj3/sigmaF2

  output <- data.frame("2MG-A",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4), round(sigma0,4),round(t(mean1),4),round(t(sigma1),4),
                       round(t(mix_pi1),4),round(t(mean2),4),round(t(sigma2),4),round(t(mix_pi2),4),
                       round(t(mean3),4),round(t(sigma3),4),round(t(mix_pi3),4),
                       round(B1[1],4)," "," "," "," "," ",round(B1[2],4),round(B1[3],4)," "," "," "," "," "," "," "," ",
                       round(jj1,4),round(ll1*100,4)," "," ",round(jj2,4),round(ll2*100,4)," "," ",round(jj3,4),round(ll3*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi1,mi2,mi3)
  return(OUTPUT)
}

#####################2MG-EA(B-4)#########################
G6FModelFun[[8]] <- function(K1,logL,df11,df21,df31,df41,df51,df61,G6Ftext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1];n_samF1<-dim(dataF1)[1];n_samP2<-dim(dataP2)[1]
  n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1];n_samF2<-dim(dataF2)[1]
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2),sum(dataF2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2),sum(dataF2^2)))
  ss<-matrix(0,3,1);ss[1]<-s[1]-sumx[1]^2/n_samP1;ss[2]<-s[2]-sumx[2]^2/n_samF1;ss[3]<-s[3]-sumx[3]^2/n_samP2
  mean<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataB1),mean(dataB2),mean(dataF2)))
  sigma0<-sum(ss)/(n_samP1+n_samF1+n_samP2-3)
  sigmaB1<-var(dataB1);sigma40<-sigmaB1
  sigmaB2<-var(dataB2);sigma50<-sigmaB2
  sigmaF2<-var(dataF2);sigma60<-sigmaF2
  sigmaP1<-sigmaF1<-sigmaP2<-sigma0
  m_esp <- 0.0001;n_fam <- as.numeric(G6Ftext2)
  ###############procedure start###########################
  mi1<-as.matrix(c(0.25,0.5,0.25))
  mean1<-as.matrix(c(mean[1],mean[4],mean[2]))
  sigma1<-matrix(0,3,1)
  mi2<-as.matrix(c(0.25,0.5,0.25))
  mean2<-as.matrix(c(mean[2],mean[5],mean[3]))
  sigma2<-matrix(0,3,1)
  mi3<-as.matrix(c(0.0625,0.125,0.25,0.25,0.25,0.0625))
  mean3<-as.matrix(c(mean[1],mean[2],mean1[2],mean[2],mean2[2],mean[3]))
  sigma3<-matrix(0,6,1)
  sigma<-sigma0
  b1<-a1<-0.2*mean[1]-0.2*mean[3]+0.1*mean1[2]-0.1*mean2[2]#da
  a1<-b1^2/n_fam
  sigma1[1]<-sigma;sigma1[2]<-sigma+0.5*a1;sigma1[3]<-sigma+a1
  sigma2[1]<-sigma1[3];sigma2[2]<-sigma1[2];sigma2[3]<-sigma
  sigma3[1]<-sigma3[2]<-sigma3[6]<-sigma
  sigma3[3]<-sigma3[5]<-sigma1[2]
  sigma3[4]<-sigma1[3]
  L0<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma0))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma0))))+
    sum(log(dnorm(dataP2,mean[3],sqrt(sigma0))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mi1)))+
    sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mi3)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,3,n_samB1);      swx1 <- matrix(0,3,1)
  W2 <- matrix(0,3,n_samB2);      swx2 <- matrix(0,3,1)
  W3 <- matrix(0,6,n_samF2);      swx3 <- matrix(0,6,1)
  s0<-matrix(0,8,1);n0<-matrix(0,8,1)
  hh<-matrix(0,3,3);b_line<-matrix(0,3,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:3) { W1[i,] <- mi1[i]*dnorm(dataB1,mean1[i],sqrt(sigma1[i]))/dmixnorm(dataB1,mean1,sqrt(sigma1),mi1) }
    mix_pi1 <- as.matrix(rowSums(W1)/n_samB1)
    sumwx1 <- W1%*%dataB1
    for(i in 1:3) { W2[i,] <- mi2[i]*dnorm(dataB2,mean2[i],sqrt(sigma2[i]))/dmixnorm(dataB2,mean2,sqrt(sigma2),mi2) }
    mix_pi2 <- as.matrix(rowSums(W2)/n_samB2)
    sumwx2 <- W2%*%dataB2
    for(i in 1:6) { W3[i,] <- mi3[i]*dnorm(dataF2,mean3[i],sqrt(sigma3[i]))/dmixnorm(dataF2,mean3,sqrt(sigma3),mi3) }
    mix_pi3 <- as.matrix(rowSums(W3)/n_samF2)
    sumwx3 <- W3%*%dataF2
    ########obtain means########################
    aaa0<-0
    s0[1]<-sumx[1]+sumwx1[1]+sumwx3[1];s0[2]<-(sumx[2]+sumwx3[2])*sigma1[3]+sigma*(sumwx1[3]+sumwx2[1]+sumwx3[4])
    s0[3]<-sumx[3]+sumwx2[3]+sumwx3[6];s0[4]<-sumwx1[2]+sumwx3[3]
    s0[5]<-sumwx2[2]+sumwx3[5]
    n0[1]<-n_samP1+mix_pi1[1]*n_samB1+mix_pi3[1]*n_samF2
    n0[2]<-(n_samF1+mix_pi3[2]*n_samF2)*sigma1[3]+sigma*(mix_pi1[3]*n_samB1+mix_pi2[1]*n_samB2+mix_pi3[4]*n_samF2)
    n0[3]<-n_samP2+mix_pi2[3]*n_samB2+mix_pi3[6]*n_samF2
    n0[4]<-mix_pi1[2]*n_samB1+mix_pi3[3]*n_samF2
    n0[5]<-mix_pi2[2]*n_samB2+mix_pi3[5]*n_samF2
    n0[c(1:5)][abs(n0[c(1:5)])<0.00000001]<-0.000001
    AA<-matrix(0,3,1);aa3<-0;aa4<-0;aaa1<-1000
    while(aaa1>0.0001)
    {
      aa4<-aa4+1
      aa6<-0.2*mean[1]-0.2*mean[3]+0.1*mean1[2]-0.1*mean2[2]
      aa6<-aa6^2/n_fam
      sigma1[2]<-sigma+0.5*aa6
      sigma1[3]<-sigma+aa6
      hh[1,1]<-sigma/n0[1]+4*sigma*sigma1[3]/n0[2]+sigma/n0[3]
      hh[1,2]<-4*sigma*sigma1[3]/n0[2]
      hh[1,3]<-sigma/n0[1]+6*sigma*sigma1[3]/n0[2]
      hh[2,2]<-sigma1[2]/n0[4]+4*sigma*sigma1[3]/n0[2]+sigma1[2]/n0[5]
      hh[2,3]<-6*sigma*sigma1[3]/n0[2]+2*sigma1[2]/n0[5]
      hh[3,3]<-sigma/n0[1]+9*sigma*sigma1[3]/n0[2]+4*sigma1[2]/n0[5]
      for(i in 2:3)
      {
        for(j in 1:(i-1))
        {
          hh[i,j]<-hh[j,i]
        }
      }
      b_line[1]<-s0[1]/n0[1]-2*s0[2]/n0[2]+s0[3]/n0[3]
      b_line[2]<-s0[4]/n0[4]-2*s0[2]/n0[2]+s0[5]/n0[5]
      b_line[3]<-s0[1]/n0[1]-3*s0[2]/n0[2]+2*s0[5]/n0[5]
      B<-solve(hh,b_line)
      mean[1]<-(s0[1]-sigma*(B[1]+B[3]))/n0[1]
      mean[2]<-(s0[2]+sigma*sigma1[3]*(2*B[1]+2*B[2]+3*B[3]))/n0[2]
      mean[3]<-(s0[3]-B[1]*sigma)/n0[3]
      mean1[2]<-(s0[4]-B[2]*sigma1[2])/n0[4]
      mean2[2]<-(s0[5]-(B[2]+2*B[3])*sigma1[2])/n0[5]
      aaa1<-max(abs(B-AA))
      AA<-B
      if (aa4>20) break
    }
    mean1[1]<-mean[1];mean1[3]<-mean[2]
    mean2[1]<-mean[2];mean2[3]<-mean[3]
    mean3[1]<-mean[1]
    mean3[2]<-mean3[4]<-mean[2]
    mean3[3]<-mean1[2];mean3[5]<-mean2[2];mean3[6]<-mean[3]
    #######obtain variance#################################
    ss1<-sum((dataP1-mean[1])^2);ss3<-sum((dataP2-mean[3])^2);ss2<-sum((dataF1-mean[2])^2)
    for(i in 1:3) {swx1[i] <- W1[i,]%*%(dataB1-mean1[i])^2 };for(i in 1:3) {swx2[i] <- W2[i,]%*%(dataB2-mean2[i])^2 };for(i in 1:6) {swx3[i] <- W3[i,]%*%(dataF2-mean3[i])^2 }
    s0[6]<-ss1+ss2+ss3+swx1[1]+swx2[3]+swx3[1]+swx3[2]+swx3[6]
    n0[6]<-n_samP1+n_samF1+n_samP2+mix_pi1[1]*n_samB1+mix_pi2[3]*n_samB2+(mix_pi3[1]+mix_pi3[2]+mix_pi3[6])*n_samF2
    s0[7]<-swx1[2]+swx2[2]+swx3[3]+swx3[5]
    s0[8]<-swx1[3]+swx2[1]+swx3[4]
    n0[7]<-mix_pi1[2]*n_samB1+mix_pi2[2]*n_samB2+(mix_pi3[3]+mix_pi3[5])*n_samF2
    n0[8]<-mix_pi1[3]*n_samB1+mix_pi2[1]*n_samB2+mix_pi3[4]*n_samF2
    aa6<-0.2*mean[1]-0.2*mean[3]+0.1*mean1[2]-0.1*mean2[2]
    aa6<-aa6*aa6/n_fam
    aaa0<-sigma;aa4<-0;aa3<-1000
    while (aa3>0.0001)
    {
      aa4<-aa4+1
      aa1<-sigma/(sigma+0.5*aa6)
      aa2<-sigma/(sigma+aa6)
      sigma<-(s0[6]+aa1^2*s0[7]+aa2^2*s0[8])/(n0[6]+aa1*n0[7]+aa2*n0[8])
      aa3<-abs(sigma-aaa0)
      aaa0<-sigma
      if (aa4>20) break
    }
    sigma1[1]<-sigma;sigma1[2]<-sigma+0.5*aa6;sigma1[3]<-sigma+aa6
    sigma2[1]<-sigma1[3];sigma2[2]<-sigma1[2];sigma2[3]<-sigma
    sigma3[1]<-sigma3[2]<-sigma3[6]<-sigma
    sigma3[3]<-sigma3[5]<-sigma1[2]
    sigma3[4]<-sigma1[3]
    if(sum(sigma < 1e-30)>=1){break}
    ##############################criteria for iterations to stop########################
    L1<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,mean[3],sqrt(sigma))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mix_pi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mix_pi3)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*3
  meanP1<-mean[1];meanF1<-mean[2];meanP2<-mean[3]
  sigma0<-sigma
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*n_samP1)
  P1bmw <- matrix(0,n_samP1,1)
  P1gg <- (dataP1 - mean[1])/sqrt(as.vector(sigma))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < n_samP1){P1bmw <- P1bmw+runif(n_samP1)/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  P1u<- as.matrix(c(12*n_samP1*((P1dd[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((P1dd[2]/n_samP1-1/3)^2),180*n_samP1*((P1dd[3]/n_samP1-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[2])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),P1D))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*n_samF1)
  F1bmw <- matrix(0,n_samF1,1)
  F1gg <- (dataF1 - mean[2])/sqrt(as.vector(sigma))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < n_samF1){F1bmw <- F1bmw+runif(n_samF1)/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  F1u<- as.matrix(c(12*n_samF1*((F1dd[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((F1dd[2]/n_samF1-1/3)^2),180*n_samF1*((F1dd[3]/n_samF1-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[2])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),F1D))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*n_samP2)
  P2bmw <- matrix(0,n_samP2,1)
  P2gg <- (dataP2 - mean[3])/sqrt(as.vector(sigma))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < n_samP2){P2bmw <- P2bmw+runif(n_samP2)/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  P2u<- as.matrix(c(12*n_samP2*((P2dd[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((P2dd[2]/n_samP2-1/3)^2),180*n_samP2*((P2dd[3]/n_samP2-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[2])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),P2D))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1);
  B1w1<-1/(12*n_samB1)
  B1bmw <- matrix(0,n_samB1,1); B1bmwsl <- matrix(0,n_samB1,3)
  for(i in 1:3){
    B1gg <- (dataB1 - mean1[i])/sqrt(sigma1[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi1[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B1P2)))[1]
  if(nn < n_samB1){B1P2 <- B1P2+runif(n_samB1)/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*n_samB1) + sum((B1P2 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  B1u <- as.matrix(c(12*n_samB1*((B1dd[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((B1dd[2]/n_samB1-1/3)^2),180*n_samB1*((B1dd[3]/n_samB1-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[2])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),B1D))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])

  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2);
  B2w1<-1/(12*n_samB2)
  B2bmw <- matrix(0,n_samB2,1); B2bmwsl <- matrix(0,n_samB2,3)
  for(i in 1:3){
    B2gg <- (dataB2 - mean2[i])/sqrt(sigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B2P2)))[1]
  if(nn < n_samB2){B2P2 <- B2P2+runif(n_samB2)/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*n_samB2) + sum((B2P2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  B2u <- as.matrix(c(12*n_samB2*((B2dd[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((B2dd[2]/n_samB2-1/3)^2),180*n_samB2*((B2dd[3]/n_samB2-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[2])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),B2D))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])

  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2);
  F2w1<-1/(12*n_samF2)
  F2bmw <- matrix(0,n_samF2,1); F2bmwsl <- matrix(0,n_samF2,6)
  for(i in 1:6){
    F2gg <- (dataF2 - mean3[i])/sqrt(sigma3[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi3[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < n_samF2){F2P2 <- F2P2+runif(n_samF2)/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*n_samF2) + sum((F2P2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  F2u <- as.matrix(c(12*n_samF2*((F2dd[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((F2dd[2]/n_samF2-1/3)^2),180*n_samF2*((F2dd[3]/n_samF2-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[2])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),F2D))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])

  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)

  ########first order parameters##################
  aa<-matrix(c(1,2,1,0,1,-2,1,1,1,-1),5,2,byrow=T)
  mm<-as.matrix(c(mean[1],mean[2],mean[3],mean1[2],mean2[2]))
  #######second order parameters##############
  B1<-solve(crossprod(aa,aa))%*%crossprod(aa,mm)
  jj1<-sigmaB1-sigma1[1]
  if (jj1<0 ) {jj1<-0}
  ll1<-jj1/sigmaB1
  jj2<-sigmaB2-sigma2[3]
  if (jj2<0)  {jj2<-0}
  ll2<-jj2/sigmaB2
  jj3<-sigmaF2-sigma3[1]
  if (jj3<0) {jj3<-0}
  ll3<-jj3/sigmaF2

  output <- data.frame("2MG-EA",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4), round(sigma0,4),round(t(mean1),4)," ",round(t(sigma1),4)," ",
                       round(t(mix_pi1),4)," ",round(t(mean2),4)," ",round(t(sigma2),4)," ",round(t(mix_pi2),4)," ",
                       round(t(mean3),4)," "," "," ",round(t(sigma3),4)," "," "," ",round(t(mix_pi3),4)," "," "," ",
                       round(B1[1],4)," "," "," "," "," ",round(B1[2],4),round(B1[2],4)," "," "," "," "," "," "," "," ",
                       round(jj1,4),round(ll1*100,4)," "," ",round(jj2,4),round(ll2*100,4)," "," ",round(jj3,4),round(ll3*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi1,mi2,mi3)
  return(OUTPUT)
}

#########################2MG-CD(B-5)##################################
G6FModelFun[[9]] <- function(K1,logL,df11,df21,df31,df41,df51,df61,G6Ftext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1];n_samF1<-dim(dataF1)[1];n_samP2<-dim(dataP2)[1]
  n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1];n_samF2<-dim(dataF2)[1]
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2),sum(dataF2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2),sum(dataF2^2)))
  ss<-matrix(0,3,1);ss[1]<-s[1]-sumx[1]^2/n_samP1;ss[2]<-s[2]-sumx[2]^2/n_samF1;ss[3]<-s[3]-sumx[3]^2/n_samP2
  mean<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataB1),mean(dataB2),mean(dataF2)))
  sigma0<-sum(ss)/(n_samP1+n_samF1+n_samP2-3)
  sigmaB1<-var(dataB1);sigma40<-sigmaB1
  sigmaB2<-var(dataB2);sigma50<-sigmaB2
  sigmaF2<-var(dataF2);sigma60<-sigmaF2
  sigmaP1<-sigmaF1<-sigmaP2<-sigma0
  m_esp <- 0.0001;n_fam <- as.numeric(G6Ftext2)
  ###############procedure start###########################
  mi1<-matrix(0.25,4,1);sigma1<-matrix(0,4,1)
  mi2<-matrix(0.25,4,1);sigma2<-matrix(0,4,1)
  mi3<-as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  sigma3<-matrix(0,9,1);sigma<-sigma0
  a1<-sqrt(sigmaB1/n_samB1)
  if (mean[1]<mean[3]) a1<--a1
  mean1<-as.matrix(c(mean[1],mean[4]-0.5*a1,mean[4]-1.5*a1,mean[4]-2.5*a1))
  a2<-sqrt(sigmaB2/n_samB2)
  if (mean[1]<mean[3]) a2<--a2
  mean2<-as.matrix(c(mean1[4],mean[5]-a2,mean[5]-2*a2,mean[3]))
  a3<-sqrt(sigmaF2/n_samF2)
  if (mean[1]<mean[3]) a3<--a3
  mean3<-as.matrix(c(mean[1],mean1[2],mean[6],mean1[3],mean1[4],mean2[2],mean[6],mean2[3],mean[3]))
  gs<-matrix(0,2,1)
  gs[1]<-(5*mean[1]-7*mean[3]+5*mean1[2]+2*mean1[3]+2*mean1[4]+2*mean2[2]-7*mean2[3]+5*mean3[3]-7*mean3[7])/39
  #da.
  gs[2]<-(5*mean[1]-7*mean[3]+2*mean1[2]+5*mean1[3]+2*mean1[4]-7*mean2[2]+2*mean2[3]-7*mean3[3]+5*mean3[7])/39
  #db.
  g_aa1<-0.75*gs[2]^2/n_fam #   0.75*db**2.
  g_aa2<-0.75*gs[1]^2/n_fam #   0.75*da**2.
  g_aa3<-g_aa1+g_aa2#   0.75*(da**2+db**2)/num_l.
  sigma1[1]<-sigma;sigma1[2]<-sigma+g_aa1;sigma1[3]<-sigma+g_aa2;sigma1[4]<-sigma+g_aa3
  sigma2[1]<-sigma1[4];sigma2[2]<-sigma1[3];sigma2[3]<-sigma1[2];sigma2[4]<-sigma
  sigma3[1]<-sigma3[3]<-sigma3[7]<-sigma3[9]<-sigma
  sigma3[2]<-sigma3[8]<-sigma1[2]
  sigma3[4]<-sigma3[6]<-sigma1[3]
  sigma3[5]<-sigma1[4]
  L0<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma0))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma0))))+
    sum(log(dnorm(dataP2,mean[3],sqrt(sigma0))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mi1)))+
    sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mi3)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,4,n_samB1);      swx1 <- matrix(0,4,1)
  W2 <- matrix(0,4,n_samB2);      swx2 <- matrix(0,4,1)
  W3 <- matrix(0,9,n_samF2);      swx3 <- matrix(0,9,1)
  s0<-matrix(0,13,1);n0<-matrix(0,13,1)
  hh<-matrix(0,6,6);b_line<-matrix(0,6,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:4) { W1[i,] <- mi1[i]*dnorm(dataB1,mean1[i],sqrt(sigma1[i]))/dmixnorm(dataB1,mean1,sqrt(sigma1),mi1) }
    mix_pi1 <- as.matrix(rowSums(W1)/n_samB1)
    sumwx1 <- W1%*%dataB1
    for(i in 1:4) { W2[i,] <- mi2[i]*dnorm(dataB2,mean2[i],sqrt(sigma2[i]))/dmixnorm(dataB2,mean2,sqrt(sigma2),mi2) }
    mix_pi2 <- as.matrix(rowSums(W2)/n_samB2)
    sumwx2 <- W2%*%dataB2
    for(i in 1:9) { W3[i,] <- mi3[i]*dnorm(dataF2,mean3[i],sqrt(sigma3[i]))/dmixnorm(dataF2,mean3,sqrt(sigma3),mi3) }
    mix_pi3 <- as.matrix(rowSums(W3)/n_samF2)
    sumwx3 <- W3%*%dataF2
    ######obtain variance##################
    aaa0<-0
    s0[1]<-sumx[1]+sumx[2]+sumwx1[1]+sumwx3[1];s0[2]<-sumx[3]+sumwx2[4]+sumwx3[9]
    s0[3]<-sumwx1[2]+sumwx3[2];s0[4]<-sumwx1[3]+sumwx3[4]
    s0[5]<-sumwx1[4]+sumwx2[1]+sumwx3[5];s0[6]<-sumwx2[2]+sumwx3[6]
    s0[7]<-sumwx2[3]+sumwx3[8];s0[8]<-sumwx3[3];s0[9]<-sumwx3[7]
    n0[1]<-n_samP1+n_samF1+mix_pi1[1]*n_samB1+mix_pi3[1]*n_samF2
    n0[2]<-n_samP2+mix_pi2[4]*n_samB2+mix_pi3[9]*n_samF2
    n0[3]<-mix_pi1[2]*n_samB1+mix_pi3[2]*n_samF2
    n0[4]<-mix_pi1[3]*n_samB1+mix_pi3[4]*n_samF2
    n0[5]<-mix_pi1[4]*n_samB1+mix_pi2[1]*n_samB2+mix_pi3[5]*n_samF2
    n0[6]<-mix_pi2[2]*n_samB2+mix_pi3[6]*n_samF2
    n0[7]<-mix_pi2[3]*n_samB2+mix_pi3[8]*n_samF2
    n0[8]<-mix_pi3[3]*n_samF2
    n0[9]<-mix_pi3[7]*n_samF2
    n0[c(1:9)][abs(n0[c(1:9)])<0.00000001]<-0.000001
    AA<-matrix(0,6,1);aa7<-0;aaa1<-1000
    while (aaa1>0.0001)
    {
      aa7<-aa7+1
      gs[1]<-(5*mean[1]-7*mean[3]+5*mean1[2]+2*mean1[3]+2*mean1[4]+2*mean2[2]-7*mean2[3]+5*mean3[3]-7*mean3[7])/39 #da.
      gs[2]<-(5*mean[1]-7*mean[3]+2*mean1[2]+5*mean1[3]+2*mean1[4]-7*mean2[2]+2*mean2[3]-7*mean3[3]+5*mean3[7])/39 #db.
      g_aa1<-0.75*gs[2]^2/n_fam #   0.75*db**2.
      g_aa2<-0.75*gs[1]^2/n_fam #   0.75*da**2.
      g_aa3<-g_aa1+g_aa2 #   0.75*(da**2+db**2)/n.
      sigma1[2]<-sigma+g_aa1;sigma1[3]<-sigma+g_aa2;sigma1[4]<-sigma+g_aa3
      sigma2[1]<-sigma1[4];sigma2[2]<-sigma1[3];sigma2[3]<-sigma1[2];sigma2[4]<-sigma
      sigma3[1]<-sigma3[3]<-sigma3[7]<-sigma3[9]<-sigma
      sigma3[2]<-sigma3[8]<-sigma1[2]
      sigma3[4]<-sigma3[6]<-sigma1[3]
      sigma3[5]<-sigma1[4]
      ##################################################
      hh[1,1]<-sigma*(1/n0[1]+1/n0[2]+1/n0[8]+1/n0[9])
      hh[1,2]<-sigma*(1/n0[1]-1/n0[2]-1/n0[8]+1/n0[9])
      hh[1,3]<-sigma*(1/n0[1]-1/n0[2]+1/n0[8]-1/n0[9])
      hh[1,4]<-sigma*(3/n0[1]+1/n0[2])
      hh[1,5]<-sigma/n0[2]
      hh[1,6]<-0
      hh[2,2]<-sigma*(1/n0[1]+1/n0[2]+1/n0[8]+1/n0[9])+sigma1[2]*(4/n0[3]+4/n0[7])
      hh[2,3]<-sigma*(1/n0[1]+1/n0[2]-1/n0[8]-1/n0[9])
      hh[2,4]<-sigma*(3/n0[1]-1/n0[2])
      hh[2,5]<--sigma/n0[2]-2*sigma1[2]/n0[7]
      hh[2,6]<--6*sigma1[2]/n0[3]+2*sigma1[2]/n0[7]
      hh[3,3]<-sigma*(1/n0[1]+1/n0[2]+1/n0[8]+1/n0[9])+sigma1[3]*(4/n0[4]+4/n0[6])
      hh[3,4]<-sigma*(3/n0[1]-1/n0[2])
      hh[3,5]<--sigma/n0[2]-2*sigma1[3]/n0[6]
      hh[3,6]<-6*sigma1[3]/n0[4]-2*sigma1[3]/n0[6]
      hh[4,4]<-sigma*(9/n0[1]+1/n0[2])+16*sigma1[4]/n0[5]
      hh[4,5]<-sigma/n0[2]-4*sigma1[4]/n0[5]
      hh[4,6]<-0
      hh[5,5]<-sigma/n0[2]+sigma1[4]/n0[5]+sigma1[3]/n0[6]+sigma1[2]/n0[7]
      hh[5,6]<-sigma1[3]/n0[6]-sigma1[2]/n0[7]
      hh[6,6]<-9*sigma1[2]/n0[3]+9*sigma1[3]/n0[4]+sigma1[3]/n0[6]+sigma1[2]/n0[7]
      for(i in 2:6)
      {
        for(j in 1:(i-1))
        {
          hh[i,j]<-hh[j,i]
        }
      }
      #################################################
      b_line[1]<-s0[1]/n0[1]+s0[2]/n0[2]-s0[8]/n0[8]-s0[9]/n0[9]
      b_line[2]<-s0[1]/n0[1]-s0[2]/n0[2]-2*s0[3]/n0[3]+2*s0[7]/n0[7]+s0[8]/n0[8]-s0[9]/n0[9]
      b_line[3]<-s0[1]/n0[1]-s0[2]/n0[2]-2*s0[4]/n0[4]+2*s0[6]/n0[6]-s0[8]/n0[8]+s0[9]/n0[9]
      b_line[4]<-3*s0[1]/n0[1]+s0[2]/n0[2]-4*s0[5]/n0[5]
      b_line[5]<-s0[2]/n0[2]+s0[5]/n0[5]-s0[6]/n0[6]-s0[7]/n0[7]
      b_line[6]<-3*s0[3]/n0[3]-3*s0[4]/n0[4]-s0[6]/n0[6]+s0[7]/n0[7]
      B<-solve(hh,b_line)
      mean[1]<-(s0[1]-sigma*(B[1]+B[2]+B[3]+3*B[4]))/n0[1]
      mean[3]<-(s0[2]-sigma*(B[1]-B[2]-B[3]+B[4]+B[5]))/n0[2]
      mean1[2]<-(s0[3]+(2*B[2]-3*B[6])*sigma1[2])/n0[3]
      mean1[3]<-(s0[4]+(2*B[3]+3*B[6])*sigma1[3])/n0[4]
      mean1[4]<-(s0[5]+sigma1[4]*(4*B[4]-B[5]))/n0[5]
      mean2[2]<-(s0[6]-sigma1[3]*(2*B[3]-B[5]-B[6]))/n0[6]
      mean2[3]<-(s0[7]-sigma1[2]*(2*B[2]-B[5]+B[6]))/n0[7]
      mean3[3]<-(s0[8]+sigma*(B[1]-B[2]+B[3]))/n0[8]
      mean3[7]<-(s0[9]+sigma*(B[1]+B[2]-B[3]))/n0[9]
      mean[2]<-mean1[1]<-mean[1]
      mean2[1]<-mean1[4];mean2[4]<-mean[3]
      mean3[1]<-mean[1];mean3[2]<-mean1[2];mean3[4]<-mean1[3]
      mean3[5]<-mean1[4];mean3[6]<-mean2[2];mean3[8]<-mean2[3];mean3[9]<-mean[3]
      aaa1<-max(abs(B-AA))
      AA<-B
      if (aa7>20) break
    }
    ##########obtain variance###########################
    ss1<-sum((dataP1-mean[1])^2);ss3<-sum((dataP2-mean[3])^2);ss2<-sum((dataF1-mean[2])^2)
    for(i in 1:4) {swx1[i] <- W1[i,]%*%(dataB1-mean1[i])^2 };for(i in 1:4) {swx2[i] <- W2[i,]%*%(dataB2-mean2[i])^2 };for(i in 1:9) {swx3[i] <- W3[i,]%*%(dataF2-mean3[i])^2 }
    s0[10]<-ss1+ss2+ss3+swx1[1]+swx2[4]+swx3[1]+swx3[3]+swx3[7]+swx3[9]
    n0[10]<-n_samP1+n_samF1+n_samP2+mix_pi1[1]*n_samB1+mix_pi2[4]*n_samB2+(mix_pi3[1]+mix_pi3[3]+mix_pi3[7]+mix_pi3[9])*n_samF2
    s0[11]<-swx1[2]+swx2[3]+swx3[2]+swx3[8]
    s0[12]<-swx1[3]+swx2[2]+swx3[4]+swx3[6]
    s0[13]<-swx1[4]+swx2[1]+swx3[5]
    n0[11]<-mix_pi1[2]*n_samB1+mix_pi2[3]*n_samB2+mix_pi3[2]*n_samF2+mix_pi3[8]*n_samF2
    n0[12]<-mix_pi1[3]*n_samB1+mix_pi2[2]*n_samB2+mix_pi3[4]*n_samF2+mix_pi3[6]*n_samF2
    n0[13]<-mix_pi1[4]*n_samB1+mix_pi2[1]*n_samB2+mix_pi3[5]*n_samF2
    gs[1]<-(5*mean[1]-7*mean[3]+5*mean1[2]+2*mean1[3]+2*mean1[4]+2*mean2[2]-7*mean2[3]+5*mean3[3]-7*mean3[7])/39 #da.
    gs[2]<-(5*mean[1]-7*mean[3]+2*mean1[2]+5*mean1[3]+2*mean1[4]-7*mean2[2]+2*mean2[3]-7*mean3[3]+5*mean3[7])/39 #db.
    g_aa1<-0.75*gs[2]^2/n_fam #   0.75*db**2.
    g_aa2<-0.75*gs[1]^2/n_fam #   0.75*da**2.
    g_aa3<-g_aa1+g_aa2 #   0.75*(da**2+db**2)/n.
    aaa0<-sigma;aa7<-0;aa3<-1000
    while (aa3>0.0001)
    {
      aa7<-aa7+1
      aa1<-sigma/(sigma+g_aa1)
      aa2<-sigma/(sigma+g_aa2)
      aa3<-sigma/(sigma+g_aa3)
      sigma<-(s0[10]+aa1^2*s0[11]+aa2^2*s0[12]+aa3^2*s0[13])/(n0[10]+aa1*n0[11]+aa2*n0[12]+aa3*n0[13])
      aa3<-abs(sigma-aaa0)
      aaa0<-sigma
      if (aa7>20) break
    }
    sigma1[1]<-sigma;sigma1[2]<-sigma+g_aa1;sigma1[3]<-sigma+g_aa2;sigma1[4]<-sigma+g_aa3
    sigma2[1]<-sigma1[4];sigma2[2]<-sigma1[3];sigma2[3]<-sigma1[2];sigma2[4]<-sigma
    sigma3[1]<-sigma3[3]<-sigma3[7]<-sigma3[9]<-sigma
    sigma3[2]<-sigma3[8]<-sigma1[2]
    sigma3[4]<-sigma3[6]<-sigma1[3]
    sigma3[5]<-sigma1[4]
    if(sum(sigma < 1e-30)>=1){break}
    ##############################criteria for iterations to stop########################
    L1<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,mean[3],sqrt(sigma))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mix_pi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mix_pi3)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*4
  meanP1<-mean[1];meanF1<-mean[2];meanP2<-mean[3]
  sigma0<-sigma
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*n_samP1)
  P1bmw <- matrix(0,n_samP1,1)
  P1gg <- (dataP1 - mean[1])/sqrt(as.vector(sigma))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < n_samP1){P1bmw <- P1bmw+runif(n_samP1)/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  P1u<- as.matrix(c(12*n_samP1*((P1dd[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((P1dd[2]/n_samP1-1/3)^2),180*n_samP1*((P1dd[3]/n_samP1-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[2])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),P1D))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*n_samF1)
  F1bmw <- matrix(0,n_samF1,1)
  F1gg <- (dataF1 - mean[2])/sqrt(as.vector(sigma))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < n_samF1){F1bmw <- F1bmw+runif(n_samF1)/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  F1u<- as.matrix(c(12*n_samF1*((F1dd[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((F1dd[2]/n_samF1-1/3)^2),180*n_samF1*((F1dd[3]/n_samF1-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[2])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),F1D))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*n_samP2)
  P2bmw <- matrix(0,n_samP2,1)
  P2gg <- (dataP2 - mean[3])/sqrt(as.vector(sigma))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < n_samP2){P2bmw <- P2bmw+runif(n_samP2)/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  P2u<- as.matrix(c(12*n_samP2*((P2dd[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((P2dd[2]/n_samP2-1/3)^2),180*n_samP2*((P2dd[3]/n_samP2-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[2])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),P2D))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1);
  B1w1<-1/(12*n_samB1)
  B1bmw <- matrix(0,n_samB1,1); B1bmwsl <- matrix(0,n_samB1,4)
  for(i in 1:4){
    B1gg <- (dataB1 - mean1[i])/sqrt(sigma1[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi1[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B1P2)))[1]
  if(nn < n_samB1){B1P2 <- B1P2+runif(n_samB1)/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*n_samB1) + sum((B1P2 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  B1u <- as.matrix(c(12*n_samB1*((B1dd[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((B1dd[2]/n_samB1-1/3)^2),180*n_samB1*((B1dd[3]/n_samB1-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[2])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),B1D))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])

  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2);
  B2w1<-1/(12*n_samB2)
  B2bmw <- matrix(0,n_samB2,1); B2bmwsl <- matrix(0,n_samB2,4)
  for(i in 1:4){
    B2gg <- (dataB2 - mean2[i])/sqrt(sigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B2P2)))[1]
  if(nn < n_samB2){B2P2 <- B2P2+runif(n_samB2)/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*n_samB2) + sum((B2P2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  B2u <- as.matrix(c(12*n_samB2*((B2dd[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((B2dd[2]/n_samB2-1/3)^2),180*n_samB2*((B2dd[3]/n_samB2-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[2])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),B2D))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])

  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2);
  F2w1<-1/(12*n_samF2)
  F2bmw <- matrix(0,n_samF2,1); F2bmwsl <- matrix(0,n_samF2,9)
  for(i in 1:9){
    F2gg <- (dataF2 - mean3[i])/sqrt(sigma3[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi3[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < n_samF2){F2P2 <- F2P2+runif(n_samF2)/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*n_samF2) + sum((F2P2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  F2u <- as.matrix(c(12*n_samF2*((F2dd[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((F2dd[2]/n_samF2-1/3)^2),180*n_samF2*((F2dd[3]/n_samF2-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[2])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),F2D))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])

  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)

  ################first order parameters########################
  aa<-matrix(c(1,1,1,1,-1,-1,1,1,0.5,1,0.5,1,1,0.5,0.5,
               1,0.5,-1,1,-1,0.5,1,1,-1,1,-1,1),9,3,byrow=T)
  mm<-as.matrix(c(mean[1],mean[3],mean1[2],mean1[3],mean1[4],mean2[2],mean2[3],mean3[3],mean3[7]))
  ############second order parameters############################
  B1<-solve(crossprod(aa,aa))%*%crossprod(aa,mm)
  jj1<-sigmaB1-sigma1[1]
  if (jj1<0 ) {jj1<-0}
  ll1<-jj1/sigmaB1
  jj2<-sigmaB2-sigma2[4]
  if (jj2<0) {jj2<-0}
  ll2<-jj2/sigmaB2
  jj3<-sigmaF2-sigma3[1]
  if (jj3<0) {jj3<-0}
  ll3<-jj3/sigmaF2

  output <- data.frame("2MG-CD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4), round(sigma0,4),round(t(mean1),4),round(t(sigma1),4),
                       round(t(mix_pi1),4),round(t(mean2),4),round(t(sigma2),4),round(t(mix_pi2),4),
                       round(t(mean3),4),round(t(sigma3),4),round(t(mix_pi3),4),
                       round(B1[1],4)," "," "," "," "," ",round(B1[2],4),round(B1[3],4),round(B1[2],4),round(B1[3],4)," "," "," "," "," "," ",
                       round(jj1,4),round(ll1*100,4)," "," ",round(jj2,4),round(ll2*100,4)," "," ",round(jj3,4),round(ll3*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi1,mi2,mi3)
  return(OUTPUT)
}

###############2MG-EAD(B-6)#########################
G6FModelFun[[10]] <- function(K1,logL,df11,df21,df31,df41,df51,df61,G6Ftext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1];n_samF1<-dim(dataF1)[1];n_samP2<-dim(dataP2)[1]
  n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1];n_samF2<-dim(dataF2)[1]
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2),sum(dataF2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2),sum(dataF2^2)))
  ss<-matrix(0,3,1);ss[1]<-s[1]-sumx[1]^2/n_samP1;ss[2]<-s[2]-sumx[2]^2/n_samF1;ss[3]<-s[3]-sumx[3]^2/n_samP2
  mean<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataB1),mean(dataB2),mean(dataF2)))
  sigma0<-sum(ss)/(n_samP1+n_samF1+n_samP2-3)
  sigmaB1<-var(dataB1);sigma40<-sigmaB1
  sigmaB2<-var(dataB2);sigma50<-sigmaB2
  sigmaF2<-var(dataF2);sigma60<-sigmaF2
  sigmaP1<-sigmaF1<-sigmaP2<-sigma0
  m_esp <- 0.0001;n_fam <- as.numeric(G6Ftext2)
  ###############procedure start###########################
  mi1<-as.matrix(c(0.25,0.5,0.25));sigma1<-matrix(0,3,1)
  mi2<-as.matrix(c(0.25,0.5,0.25));sigma2<-matrix(0,3,1)
  mi3<-as.matrix(c(0.0625,0.125,0.25,0.25,0.25,0.0625))
  sigma3<-matrix(0,6,1);sigma<-sigma0
  a1<-sqrt(sigma40/(n_samB1-1))
  if (mean[1]<mean[3]) {a1<--a1}
  mean1<-as.matrix(c(mean[1],mean[4]-a1,mean[4]))
  mean2<-as.matrix(c(mean1[3],mean[5],mean[3]))
  mean3<-as.matrix(c(mean[1],mean[6],mean1[2],mean1[3],mean2[2],mean[3]))
  a1<-(10*mean[1]-14*mean[3]+7*mean1[2]+4*mean1[3]-5*mean2[2]-2*mean3[2])/65
  a1<-a1^2/n_fam
  sigma1[1]<-sigma;sigma1[2]<-sigma+0.75*a1;sigma1[3]<-sigma+1.5*a1
  sigma2[1]<-sigma1[3];sigma2[2]<-sigma1[2];sigma2[3]<-sigma
  sigma3[1]<-sigma3[2]<-sigma3[6]<-sigma
  sigma3[3]<-sigma3[5]<-sigma1[2]
  sigma3[4]<-sigma1[3]
  L0<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma0))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma0))))+
    sum(log(dnorm(dataP2,mean[3],sqrt(sigma0))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mi1)))+
    sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mi3)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,3,n_samB1);      swx1 <- matrix(0,3,1)
  W2 <- matrix(0,3,n_samB2);      swx2 <- matrix(0,3,1)
  W3 <- matrix(0,6,n_samF2);      swx3 <- matrix(0,6,1)
  s0<-matrix(0,9,1);n0<-matrix(0,9,1)
  hh<-matrix(0,4,4);b_line<-matrix(0,4,1)
  gs<-matrix(0,1,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:3) { W1[i,] <- mi1[i]*dnorm(dataB1,mean1[i],sqrt(sigma1[i]))/dmixnorm(dataB1,mean1,sqrt(sigma1),mi1) }
    mix_pi1 <- as.matrix(rowSums(W1)/n_samB1)
    sumwx1 <- W1%*%dataB1
    for(i in 1:3) { W2[i,] <- mi2[i]*dnorm(dataB2,mean2[i],sqrt(sigma2[i]))/dmixnorm(dataB2,mean2,sqrt(sigma2),mi2) }
    mix_pi2 <- as.matrix(rowSums(W2)/n_samB2)
    sumwx2 <- W2%*%dataB2
    for(i in 1:6) { W3[i,] <- mi3[i]*dnorm(dataF2,mean3[i],sqrt(sigma3[i]))/dmixnorm(dataF2,mean3,sqrt(sigma3),mi3) }
    mix_pi3 <- as.matrix(rowSums(W3)/n_samF2)
    sumwx3 <- W3%*%dataF2
    ################obtain means####################
    aaa0<-0
    s0[1]<-sumx[1]+sumx[2]+sumwx1[1]+sumwx3[1];s0[2]<-sumx[3]+sumwx2[3]+sumwx3[6]
    s0[3]<-sumwx1[2]+sumwx3[3];s0[4]<-sumwx1[3]+sumwx2[1]+sumwx3[4]
    s0[5]<-sumwx2[2]+sumwx3[5];s0[6]<-sumwx3[2]
    n0[1]<-n_samP1+n_samF1+mix_pi1[1]*n_samB1+mix_pi3[1]*n_samF2
    n0[2]<-n_samP2+mix_pi2[3]*n_samB2+mix_pi3[6]*n_samF2
    n0[3]<-mix_pi1[2]*n_samB1+mix_pi3[3]*n_samF2
    n0[4]<-mix_pi1[3]*n_samB1+mix_pi2[1]*n_samB2+mix_pi3[4]*n_samF2
    n0[5]<-mix_pi2[2]*n_samB2+mix_pi3[5]*n_samF2
    n0[6]<-mix_pi3[2]*n_samF2
    n0[c(1:6)][abs(n0[c(1:6)])<0.00000001]<-0.000001
    AA<-matrix(0,4,1);ab5<-0;aaa1<-1000
    while(aaa1>0.0001)
    {
      ab5<-ab5+1
      gs[1]<-(10*mean[1]-14*mean[3]+7*mean1[2]+4*mean1[3]-5*mean2[2]-2*mean3[2])/65 #da.
      g_aa1<-0.75*gs[1]*gs[1]/n_fam #   0.75*d*d.
      g_aa2<-1.5*gs[1]*gs[1]/n_fam  #   1.5*d*d.
      sigma1[1]<-sigma;sigma1[2]<-sigma+g_aa1;sigma1[3]<-sigma+g_aa2
      sigma2[1]<-sigma1[3];sigma2[2]<-sigma1[2];sigma2[3]<-sigma
      sigma3[1]<-sigma3[2]<-sigma3[6]<-sigma
      sigma3[3]<-sigma3[5]<-sigma1[2]
      sigma3[4]<-sigma1[3]
      hh[1,1]<-sigma*(1/n0[1]+1/n0[2]+4/n0[6])
      hh[1,2]<-sigma*(1/n0[1]-1/n0[2])
      hh[1,3]<-sigma*(3/n0[1]+1/n0[2])
      hh[1,4]<-sigma/n0[2]
      hh[2,2]<-sigma*(1/n0[1]+1/n0[2])+sigma1[2]*(4/n0[3]+4/n0[5])
      hh[2,3]<-sigma*(3/n0[1]-1/n0[2])
      hh[2,4]<--sigma/n0[2]-4*sigma1[2]/n0[5]
      hh[3,3]<-sigma*(9/n0[1]+1/n0[2])+16*sigma1[3]/n0[4]
      hh[3,4]<-sigma/n0[2]-4*sigma1[3]/n0[4]
      hh[4,4]<-sigma/n0[2]+sigma1[3]/n0[4]+4*sigma1[2]/n0[5]
      for(i in 2:4)
      {
        for(j in 1:(i-1))
        {
          hh[i,j]<-hh[j,i]
        }
      }
      b_line[1]<-s0[1]/n0[1]+s0[2]/n0[2]-2*s0[6]/n0[6]
      b_line[2]<-s0[1]/n0[1]-s0[2]/n0[2]-2*s0[3]/n0[3]+2*s0[5]/n0[5]
      b_line[3]<-3*s0[1]/n0[1]+s0[2]/n0[2]-4*s0[4]/n0[4]
      b_line[4]<-s0[2]/n0[2]+s0[4]/n0[4]-2*s0[5]/n0[5]
      B<-solve(hh,b_line)
      mean[1]<-(s0[1]-sigma*(B[1]+B[2]+3*B[3]))/n0[1]
      mean[3]<-(s0[2]-sigma*(B[1]-B[2]+B[3]+B[4]))/n0[2]
      mean1[2]<-(s0[3]+2*B[2]*sigma1[2])/n0[3]
      mean1[3]<-(s0[4]+(4*B[3]-B[4])*sigma1[3])/n0[4]
      mean2[2]<-(s0[5]-sigma1[2]*(2*B[2]-2*B[4]))/n0[5]
      mean3[2]<-(s0[6]+2*sigma*B[1])/n0[6]
      mean[2]<-mean1[1]<-mean[1]
      mean2[1]<-mean1[3];mean2[3]<-mean[3]
      mean3[1]<-mean[1];mean3[3]<-mean1[2]
      mean3[4]<-mean1[3];mean3[5]<-mean2[2];mean3[6]<-mean[3]
      aaa1<-max(abs(B-AA))
      AA<-B
      if (ab5>20) break
    }
    ############obtain variance##################################
    ss1<-sum((dataP1-mean[1])^2);ss3<-sum((dataP2-mean[3])^2);ss2<-sum((dataF1-mean[2])^2)
    for(i in 1:3) {swx1[i] <- W1[i,]%*%(dataB1-mean1[i])^2 };for(i in 1:3) {swx2[i] <- W2[i,]%*%(dataB2-mean2[i])^2 };for(i in 1:6) {swx3[i] <- W3[i,]%*%(dataF2-mean3[i])^2 }
    s0[7]<-ss1+ss2+ss3+swx1[1]+swx2[3]+swx3[1]+swx3[2]+swx3[6]
    n0[7]<-n_samP1+n_samF1+n_samP2+mix_pi1[1]*n_samB1+mix_pi2[3]*n_samB2+(mix_pi3[1]+mix_pi3[2]+mix_pi3[6])*n_samF2
    s0[8]<-swx1[2]+swx2[2]+swx3[3]+swx3[5]
    s0[9]<-swx1[3]+swx2[1]+swx3[4]
    n0[8]<-mix_pi1[2]*n_samB1+mix_pi2[2]*n_samB2+mix_pi3[3]*n_samF2+mix_pi3[5]*n_samF2
    n0[9]<-mix_pi1[3]*n_samB1+mix_pi2[1]*n_samB2+mix_pi3[4]*n_samF2
    aaa0<-sigma
    gs[1]<-(10*mean[1]-14*mean[3]+7*mean1[2]+4*mean1[3]-5*mean2[2]-2*mean3[2])/65 #da.
    g_aa1<-0.75*gs[1]^2/n_fam
    g_aa2<-1.5*gs[1]^2/n_fam
    aa4<-1000;ab5<-0
    while (aa4>0.0001)
    {
      ab5<-ab5+1
      aa1<-sigma/(sigma+g_aa1)
      aa2<-sigma/(sigma+g_aa2)
      sigma<-(s0[7]+aa1^2*s0[8]+aa2^2*s0[9])/(n0[7]+aa1*n0[8]+aa2*n0[9])
      aa4<-abs(sigma-aaa0)
      aaa0<-sigma
      if (ab5>20) break
    }
    sigma1[1]<-sigma;sigma1[2]<-sigma+g_aa1;sigma1[3]<-sigma+g_aa2
    sigma2[1]<-sigma1[3];sigma2[2]<-sigma1[2];sigma2[3]<-sigma
    sigma3[1]<-sigma3[2]<-sigma3[6]<-sigma
    sigma3[3]<-sigma3[5]<-sigma1[2]
    sigma3[4]<-sigma1[3]
    if(sum(sigma < 1e-30)>=1){break}
    ##############################criteria for iterations to stop########################
    L1<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,mean[3],sqrt(sigma))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mix_pi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mix_pi3)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*3
  meanP1<-mean[1];meanF1<-mean[2];meanP2<-mean[3]
  sigma0<-sigma
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*n_samP1)
  P1bmw <- matrix(0,n_samP1,1)
  P1gg <- (dataP1 - mean[1])/sqrt(as.vector(sigma))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < n_samP1){P1bmw <- P1bmw+runif(n_samP1)/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  P1u<- as.matrix(c(12*n_samP1*((P1dd[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((P1dd[2]/n_samP1-1/3)^2),180*n_samP1*((P1dd[3]/n_samP1-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[2])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),P1D))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*n_samF1)
  F1bmw <- matrix(0,n_samF1,1)
  F1gg <- (dataF1 - mean[2])/sqrt(as.vector(sigma))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < n_samF1){F1bmw <- F1bmw+runif(n_samF1)/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  F1u<- as.matrix(c(12*n_samF1*((F1dd[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((F1dd[2]/n_samF1-1/3)^2),180*n_samF1*((F1dd[3]/n_samF1-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[2])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),F1D))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*n_samP2)
  P2bmw <- matrix(0,n_samP2,1)
  P2gg <- (dataP2 - mean[3])/sqrt(as.vector(sigma))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < n_samP2){P2bmw <- P2bmw+runif(n_samP2)/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  P2u<- as.matrix(c(12*n_samP2*((P2dd[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((P2dd[2]/n_samP2-1/3)^2),180*n_samP2*((P2dd[3]/n_samP2-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[2])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),P2D))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1);
  B1w1<-1/(12*n_samB1)
  B1bmw <- matrix(0,n_samB1,1); B1bmwsl <- matrix(0,n_samB1,3)
  for(i in 1:3){
    B1gg <- (dataB1 - mean1[i])/sqrt(sigma1[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi1[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B1P2)))[1]
  if(nn < n_samB1){B1P2 <- B1P2+runif(n_samB1)/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*n_samB1) + sum((B1P2 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  B1u <- as.matrix(c(12*n_samB1*((B1dd[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((B1dd[2]/n_samB1-1/3)^2),180*n_samB1*((B1dd[3]/n_samB1-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[2])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),B1D))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])

  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2);
  B2w1<-1/(12*n_samB2)
  B2bmw <- matrix(0,n_samB2,1); B2bmwsl <- matrix(0,n_samB2,3)
  for(i in 1:3){
    B2gg <- (dataB2 - mean2[i])/sqrt(sigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B2P2)))[1]
  if(nn < n_samB2){B2P2 <- B2P2+runif(n_samB2)/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*n_samB2) + sum((B2P2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  B2u <- as.matrix(c(12*n_samB2*((B2dd[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((B2dd[2]/n_samB2-1/3)^2),180*n_samB2*((B2dd[3]/n_samB2-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[2])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),B2D))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])

  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2);
  F2w1<-1/(12*n_samF2)
  F2bmw <- matrix(0,n_samF2,1); F2bmwsl <- matrix(0,n_samF2,6)
  for(i in 1:6){
    F2gg <- (dataF2 - mean3[i])/sqrt(sigma3[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi3[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < n_samF2){F2P2 <- F2P2+runif(n_samF2)/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*n_samF2) + sum((F2P2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  F2u <- as.matrix(c(12*n_samF2*((F2dd[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((F2dd[2]/n_samF2-1/3)^2),180*n_samF2*((F2dd[3]/n_samF2-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[2])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),F2D))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])

  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)

  ########first order parameters#############################
  aa<-matrix(c(1,2,1,-2,1,0.5,1,1,1,-0.5,1,0),6,2,byrow=T)
  mm<-as.matrix(c(mean[1],mean[3],mean1[2],mean1[3],mean2[2],mean3[2]))
  B1<-solve(crossprod(aa,aa))%*%crossprod(aa,mm)
  #######second order parameters################################
  jj1<-sigmaB1-sigma1[1]
  if (jj1<0 ) {jj1<-0}
  ll1<-jj1/sigmaB1
  jj2<-sigmaB2-sigma2[3]
  if (jj2<0) {jj2<-0}
  ll2<-jj2/sigmaB2
  jj3<-sigmaF2-sigma3[1]
  if (jj3<0) {jj3<-0}
  ll3<-jj3/sigmaF2

  output <- data.frame("2MG-EAD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4), round(sigma0,4),round(t(mean1),4)," ",round(t(sigma1),4)," ",
                       round(t(mix_pi1),4)," ",round(t(mean2),4)," ",round(t(sigma2),4)," ",round(t(mix_pi2),4)," ",
                       round(t(mean3),4)," "," "," ",round(t(sigma3),4)," "," "," ",round(t(mix_pi3),4)," "," "," ",
                       round(B1[1],4)," "," "," "," "," ",round(B1[2],4),round(B1[2],4),round(B1[2],4),round(B1[2],4)," "," "," "," "," "," ",
                       round(jj1,4),round(ll1*100,4)," "," ",round(jj2,4),round(ll2*100,4)," "," ",round(jj3,4),round(ll3*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi1,mi2,mi3)
  return(OUTPUT)
}


##################PG-ADI(C-0)#########################################
G6FModelFun[[11]] <- function(K1,logL,df11,df21,df31,df41,df51,df61,G6Ftext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1];n_samF1<-dim(dataF1)[1];n_samP2<-dim(dataP2)[1]
  n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1];n_samF2<-dim(dataF2)[1]
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2),sum(dataF2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2),sum(dataF2^2)))
  ss<-matrix(0,3,1);ss[1]<-s[1]-sumx[1]^2/n_samP1;ss[2]<-s[2]-sumx[2]^2/n_samF1;ss[3]<-s[3]-sumx[3]^2/n_samP2
  mean<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataB1),mean(dataB2),mean(dataF2)))
  sigma0<-sum(ss)/(n_samP1+n_samF1+n_samP2-3)
  sigmaB1<-var(dataB1);sigma40<-sigmaB1
  sigmaB2<-var(dataB2);sigma50<-sigmaB2
  sigmaF2<-var(dataF2);sigma60<-sigmaF2
  sigmaP1<-sigmaF1<-sigmaP2<-sigma0
  m_esp <- 0.0001;n_fam <- as.numeric(G6Ftext2)
  ###############procedure start###########################
  sigma<-sigma0
  mix_pi1<-matrix(1,1,1);mean1<-matrix(mean[4],1,1);sigma1<-matrix(0,1,1)
  mix_pi2<-matrix(1,1,1);mean2<-matrix(mean[5],1,1);sigma2<-matrix(0,1,1)
  mix_pi3<-matrix(1,1,1);mean3<-matrix(mean[6],1,1);sigma3<-matrix(0,1,1)
  L0<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma0))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma0))))+
    sum(log(dnorm(dataP2,mean[3],sqrt(sigma0))))+sum(log(dnorm(dataB1,mean1[1],sqrt(sigma40))))+
    sum(log(dnorm(dataB2,mean2[1],sqrt(sigma50))))+sum(log(dnorm(dataF2,mean3[1],sqrt(sigma60))))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ss1<-sum((dataP1-mean[1])^2);ss3<-sum((dataP2-mean[3])^2)
    ss2<-sum((dataF1-mean[2])^2);ss4<-sigma40*(n_samB1-1)
    ss5<-sigma50*(n_samB2-1);ss6<-sigma60*(n_samF2-1)
    abc1<-ss1+ss2+ss3;abc2<-n_samP1+n_samF1+n_samP2
    aaa0<-sigma;aa3<-1000;n_iter<-0
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma/sigma40
      aa2<-sigma/sigma50
      aa3<-sigma/sigma60
      sigma<-(abc1+aa1^2*ss4+aa2^2*ss5+aa3^2*ss6)/(abc2+aa1*n_samB1+aa2*n_samB2+aa3*n_samF2)
      aa3<-abs(sigma-aaa0)
      aaa0<-sigma
      if(n_iter>20)break
    }
    sigma_4<-sigma40-sigma;
    if (sigma_4<0) {sigma_4<-0;sigma40<-sigma}
    sigma40<-sigma_4+sigma
    sigma_5<-sigma50-sigma;
    if (sigma_5<0) {sigma_5<-0;sigma50<-sigma}
    sigma50<-sigma_5+sigma
    sigma_6<-sigma60-sigma;
    if (sigma_6<0) {sigma_6<-0;sigma60<-sigma}
    sigma60<-sigma_6+sigma
    if(sum(sigma < 1e-30)>=1){break}
    ##############################criteria for iterations to stop########################
    L1<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,mean[3],sqrt(sigma))))+sum(log(dnorm(dataB1,mean1[1],sqrt(sigma40))))+
      sum(log(dnorm(dataB2,mean2[1],sqrt(sigma50))))+sum(log(dnorm(dataF2,mean3[1],sqrt(sigma60))))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc<-L0
  AIC<--2*abc+2*10
  meanP1<-mean[1];meanF1<-mean[2];meanP2<-mean[3]
  sigma0<-sigma
  sigma1[1]<-sigma40;sigma2[1]<-sigma50;sigma3[1]<-sigma60
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*n_samP1)
  P1bmw <- matrix(0,n_samP1,1)
  P1gg <- (dataP1 - mean[1])/sqrt(as.vector(sigma))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < n_samP1){P1bmw <- P1bmw+runif(n_samP1)/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  P1u<- as.matrix(c(12*n_samP1*((P1dd[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((P1dd[2]/n_samP1-1/3)^2),180*n_samP1*((P1dd[3]/n_samP1-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[2])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),P1D))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*n_samF1)
  F1bmw <- matrix(0,n_samF1,1)
  F1gg <- (dataF1 - mean[2])/sqrt(as.vector(sigma))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < n_samF1){F1bmw <- F1bmw+runif(n_samF1)/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  F1u<- as.matrix(c(12*n_samF1*((F1dd[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((F1dd[2]/n_samF1-1/3)^2),180*n_samF1*((F1dd[3]/n_samF1-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[2])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),F1D))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*n_samP2)
  P2bmw <- matrix(0,n_samP2,1)
  P2gg <- (dataP2 - mean[3])/sqrt(as.vector(sigma))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < n_samP2){P2bmw <- P2bmw+runif(n_samP2)/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  P2u<- as.matrix(c(12*n_samP2*((P2dd[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((P2dd[2]/n_samP2-1/3)^2),180*n_samP2*((P2dd[3]/n_samP2-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[2])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),P2D))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1);
  B1w1<-1/(12*n_samB1)
  B1bmw <- matrix(0,n_samB1,1); B1bmwsl <- matrix(0,n_samB1,1)
  for(i in 1:1){
    B1gg <- (dataB1 - mean1[i])/sqrt(as.vector(sigma1[i]))
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi1[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B1P2)))[1]
  if(nn < n_samB1){B1P2 <- B1P2+runif(n_samB1)/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*n_samB1) + sum((B1P2 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  B1u <- as.matrix(c(12*n_samB1*((B1dd[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((B1dd[2]/n_samB1-1/3)^2),180*n_samB1*((B1dd[3]/n_samB1-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[2])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),B1D))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])

  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2);
  B2w1<-1/(12*n_samB2)
  B2bmw <- matrix(0,n_samB2,1); B2bmwsl <- matrix(0,n_samB2,1)
  for(i in 1:1){
    B2gg <- (dataB2 - mean2[i])/sqrt(as.vector(sigma2[i]))
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B2P2)))[1]
  if(nn < n_samB2){B2P2 <- B2P2+runif(n_samB2)/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*n_samB2) + sum((B2P2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  B2u <- as.matrix(c(12*n_samB2*((B2dd[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((B2dd[2]/n_samB2-1/3)^2),180*n_samB2*((B2dd[3]/n_samB2-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[2])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),B2D))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])

  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2);
  F2w1<-1/(12*n_samF2)
  F2bmw <- matrix(0,n_samF2,1); F2bmwsl <- matrix(0,n_samF2,1)
  for(i in 1:1){
    F2gg <- (dataF2 - mean3[i])/sqrt(as.vector(sigma3[i]))
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi3[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < n_samF2){F2P2 <- F2P2+runif(n_samF2)/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*n_samF2) + sum((F2P2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  F2u <- as.matrix(c(12*n_samF2*((F2dd[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((F2dd[2]/n_samF2-1/3)^2),180*n_samF2*((F2dd[3]/n_samF2-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[2])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),F2D))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])

  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)

  #########first order parameters#####################
  m1<-meanP1;m2<-meanF1;m3<-meanP2
  m4<-mean1[1];m5<-mean2[1];m6<-mean3[1]
  #########second order parameters##################
  mm1<-sigma40-sigma
  if (mm1<0 || mm1>=sigma40) {mm1<-0}
  nn1<-mm1/sigma40
  mm2<-sigma50-sigma
  if (mm2<0 || mm2>=sigma50) {mm2<-0}
  nn2<-mm2/sigma50
  mm3<-sigma60-sigma
  if (mm3<0 || mm3>=sigmaF2) {mm3<-0}
  nn3<-mm3/sigma60

  output <- data.frame("PG-ADI",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4),round(sigma0,4),round(t(mean1),4)," "," "," ",round(t(sigma1),4)," "," "," ",
                       round(t(mix_pi1),4)," "," "," ",round(t(mean2),4)," "," "," ",round(t(sigma2),4)," "," "," ",round(t(mix_pi2),4)," "," "," ",
                       round(t(mean3),4)," "," "," "," "," "," "," "," ",round(t(sigma3),4)," "," "," "," "," "," "," "," ",round(t(mix_pi3),4)," "," "," "," "," "," "," "," ",
                       round(m1,4),round(m2,4),round(m3,4),round(m4,4),round(m5,4),round(m6,4)," "," "," "," "," "," "," "," "," "," "," "," ",
                       round(mm1,4),round(nn1*100,4)," "," ",round(mm2,4),round(nn2*100,4)," "," ",round(mm3,4),round(nn3*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output)
  return(OUTPUT)
}


#######################PG-AD(C-1)#############################
G6FModelFun[[12]] <- function(K1,logL,df11,df21,df31,df41,df51,df61,G6Ftext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1];n_samF1<-dim(dataF1)[1];n_samP2<-dim(dataP2)[1]
  n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1];n_samF2<-dim(dataF2)[1]
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2),sum(dataF2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2),sum(dataF2^2)))
  ss<-matrix(0,3,1);ss[1]<-s[1]-sumx[1]^2/n_samP1;ss[2]<-s[2]-sumx[2]^2/n_samF1;ss[3]<-s[3]-sumx[3]^2/n_samP2
  mean<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataB1),mean(dataB2),mean(dataF2)))
  sigma0<-sum(ss)/(n_samP1+n_samF1+n_samP2-3)
  sigmaB1<-var(dataB1);sigma40<-sigmaB1
  sigmaB2<-var(dataB2);sigma50<-sigmaB2
  sigmaF2<-var(dataF2);sigma60<-sigmaF2
  sigmaP1<-sigmaF1<-sigmaP2<-sigma0
  m_esp <- 0.0001;n_fam <- as.numeric(G6Ftext2)
  ###############procedure start###########################
  sigma<-sigma0
  mix_pi1<-matrix(1,1,1);mean1<-matrix(mean[4],1,1);sigma1<-matrix(0,1,1)
  mix_pi2<-matrix(1,1,1);mean2<-matrix(mean[5],1,1);sigma2<-matrix(0,1,1)
  mix_pi3<-matrix(1,1,1);mean3<-matrix(mean[6],1,1);sigma3<-matrix(0,1,1)
  L0<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma0))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma0))))+
    sum(log(dnorm(dataP2,mean[3],sqrt(sigma0))))+sum(log(dnorm(dataB1,mean1[1],sqrt(sigma40))))+
    sum(log(dnorm(dataB2,mean2[1],sqrt(sigma50))))+sum(log(dnorm(dataF2,mean3[1],sqrt(sigma60))))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  hh<-matrix(0,3,3);b_line<-matrix(0,3,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    aaa0<-0
    sigma40<-sum((dataB1-mean[4])^2);sigma40<-sigma40/(n_samB1-1)
    sigma50<-sum((dataB2-mean[5])^2);sigma50<-sigma50/(n_samB2-1)
    sigma60<-sum((dataF2-mean[6])^2);sigma60<-sigma60/(n_samF2-1)
    hh[1,1]<-sigma*(1/n_samP1+1/n_samP2)+4*sigma40/n_samB1+4*sigma50/n_samB2
    hh[1,2]<-3*sigma*(1/n_samP1-1/n_samP2)
    hh[1,3]<--2*sigma40/n_samB1+2*sigma50/n_samB2
    hh[2,2]<-sigma*(9/n_samP1+4/n_samF1+9/n_samP2)+64*sigma60/n_samF2
    hh[2,3]<-16*sigma60/n_samF2
    hh[3,3]<-sigma40/n_samB1+sigma50/n_samB2+4*sigma60/n_samF2
    for(i in 2:3)
    {
      for(j in 1:(i-1))
      {
        hh[i,j]<-hh[j,i]
      }
    }
    ###################################################
    b_line[1]<-sumx[1]/n_samP1-sumx[3]/n_samP2-2*sumx[4]/n_samB1+2*sumx[5]/n_samB2
    b_line[2]<-3*sumx[1]/n_samP1+2*sumx[2]/n_samF1+3*sumx[3]/n_samP2-8*sumx[6]/n_samF2
    b_line[3]<-sumx[4]/n_samB1+sumx[5]/n_samB2-2*sumx[6]/n_samF2
    B<-solve(hh,b_line)
    mean[1]<-(sumx[1]-sigma*(B[1]+3*B[2]))/n_samP1
    mean[2]<-(sumx[2]-2*sigma*B[2])/n_samF1
    mean[3]<-(sumx[3]+sigma*(B[1]-3*B[2]))/n_samP2
    mean1[1]<-(sumx[4]+sigma40*(2*B[1]-B[3]))/n_samB1
    mean2[1]<-(sumx[5]-(2*B[1]+B[3])*sigma50)/n_samB2
    mean3[1]<-(sumx[6]+(8*B[2]+2*B[3])*sigma60)/n_samF2
    aaa0<-sigma
    ss1<-sum((dataP1-mean[1])^2);ss3<-sum((dataP2-mean[3])^2);ss2<-sum((dataF1-mean[2])^2)
    abc1<-ss1+ss2+ss3;abc2<-n_samP1+n_samF1+n_samP2
    ss4<-sum((dataB1-mean[4])^2);ss5<-sum((dataB2-mean[5])^2);ss6<-sum((dataF2-mean[6])^2)
    sigma40<-ss4/(n_samB1-1)
    sigma_4<-sigma40-sigma;if (sigma_4<0) {sigma_4<-0}
    sigma40<-sigma_4+sigma;sigma50<-ss5/(n_samB2-1)
    sigma_5<-sigma50-sigma;if (sigma_5<0) {sigma_5<-0}
    sigma50<-sigma_5+sigma;sigma60<-ss6/(n_samF2-1)
    sigma_6<-sigma60-sigma;if (sigma_6<0) {sigma_6<-0}
    sigma60<-sigma_6+sigma
    aa3<-1000;n_iter<-0
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma/sigma40
      if (aa1>=1) {aa1<-1}
      aa2<-sigma/sigma50
      if (aa2>=1) {aa2<-1}
      aa3<-sigma/sigma60
      if (aa3>=1) {aa3<-1}
      aa4<-abc1+aa1^2*ss4+aa2^2*ss5+aa3^2*ss6
      aa5<-abc2+aa1*n_samB1+aa2*n_samB2+aa3*n_samF2
      sigma<-aa4/aa5
      aa3<-abs(sigma-aaa0)
      aaa0<-sigma
      if(n_iter>20)break
    }
    sigma40<-sigma_4+sigma;sigma50<-sigma_5+sigma;sigma60<-sigma_6+sigma
    if(sum(sigma < 1e-30)>=1){break}
    ##############################criteria for iterations to stop########################
    L1<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,mean[3],sqrt(sigma))))+sum(log(dnorm(dataB1,mean1[1],sqrt(sigma40))))+
      sum(log(dnorm(dataB2,mean2[1],sqrt(sigma50))))+sum(log(dnorm(dataF2,mean3[1],sqrt(sigma60))))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*7
  meanP1<-mean[1];meanF1<-mean[2];meanP2<-mean[3]
  sigma0<-sigma
  sigma1[1]<-sigma40;sigma2[1]<-sigma50;sigma3[1]<-sigma60
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*n_samP1)
  P1bmw <- matrix(0,n_samP1,1)
  P1gg <- (dataP1 - mean[1])/sqrt(as.vector(sigma))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < n_samP1){P1bmw <- P1bmw+runif(n_samP1)/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  P1u<- as.matrix(c(12*n_samP1*((P1dd[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((P1dd[2]/n_samP1-1/3)^2),180*n_samP1*((P1dd[3]/n_samP1-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[2])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),P1D))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*n_samF1)
  F1bmw <- matrix(0,n_samF1,1)
  F1gg <- (dataF1 - mean[2])/sqrt(as.vector(sigma))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < n_samF1){F1bmw <- F1bmw+runif(n_samF1)/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  F1u<- as.matrix(c(12*n_samF1*((F1dd[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((F1dd[2]/n_samF1-1/3)^2),180*n_samF1*((F1dd[3]/n_samF1-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[2])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),F1D))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*n_samP2)
  P2bmw <- matrix(0,n_samP2,1)
  P2gg <- (dataP2 - mean[3])/sqrt(as.vector(sigma))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < n_samP2){P2bmw <- P2bmw+runif(n_samP2)/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  P2u<- as.matrix(c(12*n_samP2*((P2dd[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((P2dd[2]/n_samP2-1/3)^2),180*n_samP2*((P2dd[3]/n_samP2-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[2])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),P2D))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1);
  B1w1<-1/(12*n_samB1)
  B1bmw <- matrix(0,n_samB1,1); B1bmwsl <- matrix(0,n_samB1,1)
  for(i in 1:1){
    B1gg <- (dataB1 - mean1[i])/sqrt(as.vector(sigma1[i]))
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi1[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B1P2)))[1]
  if(nn < n_samB1){B1P2 <- B1P2+runif(n_samB1)/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*n_samB1) + sum((B1P2 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  B1u <- as.matrix(c(12*n_samB1*((B1dd[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((B1dd[2]/n_samB1-1/3)^2),180*n_samB1*((B1dd[3]/n_samB1-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[2])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),B1D))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])

  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2);
  B2w1<-1/(12*n_samB2)
  B2bmw <- matrix(0,n_samB2,1); B2bmwsl <- matrix(0,n_samB2,1)
  for(i in 1:1){
    B2gg <- (dataB2 - mean2[i])/sqrt(as.vector(sigma2[i]))
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B2P2)))[1]
  if(nn < n_samB2){B2P2 <- B2P2+runif(n_samB2)/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*n_samB2) + sum((B2P2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  B2u <- as.matrix(c(12*n_samB2*((B2dd[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((B2dd[2]/n_samB2-1/3)^2),180*n_samB2*((B2dd[3]/n_samB2-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[2])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),B2D))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])

  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2);
  F2w1<-1/(12*n_samF2)
  F2bmw <- matrix(0,n_samF2,1); F2bmwsl <- matrix(0,n_samF2,1)
  for(i in 1:1){
    F2gg <- (dataF2 - mean3[i])/sqrt(as.vector(sigma3[i]))
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi3[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < n_samF2){F2P2 <- F2P2+runif(n_samF2)/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*n_samF2) + sum((F2P2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  F2u <- as.matrix(c(12*n_samF2*((F2dd[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((F2dd[2]/n_samF2-1/3)^2),180*n_samF2*((F2dd[3]/n_samF2-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[2])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),F2D))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])

  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)

  ############first order parameters###########
  aa<-matrix(c(1,1,0,1,0,1,1,-1,0,1,0.5,0.25,1,-0.5,0.25,1,0,0.25),6,3,byrow=T)
  mm<-as.matrix(c(mean[1],mean[2],mean[3],mean1[1],mean2[1],mean3[1]))
  B1<-solve(crossprod(aa,aa))%*%crossprod(aa,mm)
  ########second order parameters##########
  mm1<-sigma40-sigma
  if (mm1<0 || mm1>=sigma40) {mm1<-0}
  nn1<-mm1/sigma40
  mm2<-sigma50-sigma
  if (mm2<0 || mm2>=sigma50) {mm2<-0}
  nn2<-mm2/sigma50
  mm3<-sigma60-sigma
  if (mm3<0 || mm3>=sigma60) {mm3<-0}
  nn3<-mm3/sigma60

  output <- data.frame("PG-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4), round(sigma0,4),round(t(mean1),4)," "," "," ",round(t(sigma1),4)," "," "," ",
                       round(t(mix_pi1),4)," "," "," ",round(t(mean2),4)," "," "," ",round(t(sigma2),4)," "," "," ",round(t(mix_pi2),4)," "," "," ",
                       round(t(mean3),4)," "," "," "," "," "," "," "," ",round(t(sigma3),4)," "," "," "," "," "," "," "," ",round(t(mix_pi3),4)," "," "," "," "," "," "," "," ",
                       round(B1[1],4)," "," "," "," "," "," "," "," "," "," "," "," "," ",round(B1[2],4),round(B1[3],4)," "," ",
                       round(mm1,4),round(nn1*100,4)," "," ",round(mm2,4),round(nn2*100,4)," "," ",round(mm3,4),round(nn3*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output)
  return(OUTPUT)
}

############MX1-AD-ADI(D-0)################################
G6FModelFun[[13]] <- function(K1,logL,df11,df21,df31,df41,df51,df61,G6Ftext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1];n_samF1<-dim(dataF1)[1];n_samP2<-dim(dataP2)[1]
  n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1];n_samF2<-dim(dataF2)[1]
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2),sum(dataF2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2),sum(dataF2^2)))
  ss<-matrix(0,3,1);ss[1]<-s[1]-sumx[1]^2/n_samP1;ss[2]<-s[2]-sumx[2]^2/n_samF1;ss[3]<-s[3]-sumx[3]^2/n_samP2
  mean<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataB1),mean(dataB2),mean(dataF2)))
  sigma0<-sum(ss)/(n_samP1+n_samF1+n_samP2-3)
  sigmaB1<-var(dataB1);sigma40<-sigmaB1
  sigmaB2<-var(dataB2);sigma50<-sigmaB2
  sigmaF2<-var(dataF2);sigma60<-sigmaF2
  sigmaP1<-sigmaF1<-sigmaP2<-sigma0
  m_esp <- 0.0001;n_fam <- as.numeric(G6Ftext2)
  ###############procedure start###########################
  mi1<-matrix(0.5,2,1);sigma1<-matrix(0,2,1)
  mi2<-matrix(0.5,2,1);sigma2<-matrix(0,2,1)
  mi3<-as.matrix(c(0.25,0.5,0.25));sigma3<-matrix(0,3,1)
  sigma<-sigma0
  a1<-sqrt(sigmaB1/n_samB1);if (mean[1]<mean[3]) a1<--a1
  mean1<-as.matrix(c(mean[4]+a1,mean[4]))
  a2<-sqrt(sigmaB2/n_samB2);if (mean[1]<mean[3]) a2<--a2
  mean2<-as.matrix(c(mean[5],mean[5]-a2))
  a3<-sqrt(sigmaF2/n_samF2);if (mean[1]<mean[3]) a3<--a3
  mean3<-as.matrix(c(mean[6]+a3,mean[6],mean[6]-a3))
  b1<-(mean1[1]-mean1[2]+mean2[1]-mean2[2]+2*mean3[1]-2*mean3[3])/6# additive effect.
  b2<--0.6*mean1[1]+0.6*mean1[2]+0.6*mean2[1]-0.6*mean2[2]-0.4*mean3[1]+0.8*mean3[2]-0.4*mean3[3]
  b3<-(0.5*b1^2+0.25*b2^2)/n_fam
  sigma1[1]<-sigmaB1/2;sigma1[2]<-sigma1[1]+b3
  sigma2[2]<-sigmaB2/2;sigma2[1]<-sigma2[2]+b3
  sigma3[1]<-sigmaF2/2;sigma3[3]<-sigma3[1];sigma3[2]<-sigma3[1]+b3
  L0<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma0))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma0))))+
    sum(log(dnorm(dataP2,mean[3],sqrt(sigma0))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mi1)))+
    sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mi3)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,2,n_samB1);      swx1 <- matrix(0,2,1)
  W2 <- matrix(0,2,n_samB2);      swx2 <- matrix(0,2,1)
  W3 <- matrix(0,3,n_samF2);      swx3 <- matrix(0,3,1)
  n0<-matrix(0,9,1);s0<-matrix(0,9,1); rr<-matrix(0,2,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:2) { W1[i,] <- mi1[i]*dnorm(dataB1,mean1[i],sqrt(sigma1[i]))/dmixnorm(dataB1,mean1,sqrt(sigma1),mi1) }
    mix_pi1 <- as.matrix(rowSums(W1)/n_samB1)
    sumwx1 <- W1%*%dataB1
    for(i in 1:2) { W2[i,] <- mi2[i]*dnorm(dataB2,mean2[i],sqrt(sigma2[i]))/dmixnorm(dataB2,mean2,sqrt(sigma2),mi2) }
    mix_pi2 <- as.matrix(rowSums(W2)/n_samB2)
    sumwx2 <- W2%*%dataB2
    for(i in 1:3) { W3[i,] <- mi3[i]*dnorm(dataF2,mean3[i],sqrt(sigma3[i]))/dmixnorm(dataF2,mean3,sqrt(sigma3),mi3) }
    mix_pi3 <- as.matrix(rowSums(W3)/n_samF2)
    sumwx3 <- W3%*%dataF2
    n0[1]<-mix_pi1[1]*n_samB1;n0[2]<-mix_pi1[2]*n_samB1
    n0[3]<-mix_pi2[1]*n_samB2;n0[4]<-mix_pi2[2]*n_samB2
    n0[5]<-mix_pi3[1]*n_samF2;n0[6]<-mix_pi3[2]*n_samF2
    n0[7]<-mix_pi3[3]*n_samF2
    n0[c(1:7)][abs(n0[c(1:7)])<0.0001]<-0.0001
    ###########obtain means##########################
    aaa0<-0    # CM1-step for means.
    s0[1]<-sumwx1[1]/n0[1]-sumwx1[2]/n0[2]-sumwx3[1]/n0[5]+sumwx3[2]/n0[6]
    s0[2]<-sumwx2[1]/n0[3]-sumwx2[2]/n0[4]-sumwx3[2]/n0[6]+sumwx3[3]/n0[7]
    abc5<-0;abc6<-0;n_iter<-0;aaa1<-1000;AA<-matrix(0,2,1)
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-(mean1[1]-mean1[2]+mean2[1]-mean2[2]+2*mean3[1]-2*mean3[3])/6
      aa2<--0.6*mean1[1]+0.6*mean1[2]+0.6*mean2[1]-0.6*mean2[2]-0.4*mean3[1]+0.8*mean3[2]-0.4*mean3[3]
      aa1<-(0.5*aa1*aa1+0.25*aa2*aa2)/n_fam
      sigma1[2]<-sigma1[1]+aa1
      sigma2[1]<-sigma2[2]+aa1
      sigma3[2]<-sigma3[1]+aa1
      abc1<-sigma1[1]/n0[1]+sigma1[2]/n0[2]+sigma3[1]/n0[5]+sigma3[2]/n0[6]
      abc2<--sigma3[2]/n0[6]
      abc3<-sigma2[1]/n0[3]+sigma2[2]/n0[4]+sigma3[2]/n0[6]+sigma3[3]/n0[7]
      aa2<-abc1*abc3-abc2^2
      aa3<-s0[1]*abc3-s0[2]*abc2
      aa4<-s0[2]*abc1-s0[1]*abc2
      rr[1]<-aa3/aa2;rr[2]<-aa4/aa2
      mean1[1]<-(sumwx1[1]-rr[1]*sigma1[1])/n0[1];mean1[2]<-(sumwx1[2]+rr[1]*sigma1[2])/n0[2]
      mean2[1]<-(sumwx2[1]-rr[2]*sigma2[1])/n0[3];mean2[2]<-(sumwx2[2]+rr[2]*sigma2[2])/n0[4]
      mean3[1]<-(sumwx3[1]+rr[1]*sigma3[1])/n0[5];mean3[2]<-(sumwx3[2]+sigma3[2]*(-rr[1]+rr[2]))/n0[6]
      mean3[3]<-(sumwx3[3]-rr[2]*sigma3[3])/n0[7]
      aaa1<-max(abs(rr-AA))
      AA<-rr
      if (n_iter>20) break
    }
    ##########obtain variance#############################
    ss1<-sum((dataP1-mean[1])^2);ss3<-sum((dataP2-mean[3])^2);ss2<-sum((dataF1-mean[2])^2)
    for(i in 1:2) {swx1[i] <- W1[i,]%*%(dataB1-mean1[i])^2 };for(i in 1:2) {swx2[i] <- W2[i,]%*%(dataB2-mean2[i])^2 };for(i in 1:3) {swx3[i] <- W3[i,]%*%(dataF2-mean3[i])^2 }
    n0[8]<-mix_pi1[1]*n_samB1;n0[9]<-mix_pi1[2]*n_samB1
    aaa0<-sigma1[1];n_iter<-0;aa3<-1000
    aa1<-(mean1[1]-mean1[2]+mean2[1]-mean2[2]+2*mean3[1]-2*mean3[3])/6
    aa2<--0.6*mean1[1]+0.6*mean1[2]+0.6*mean2[1]-0.6*mean2[2]-0.4*mean3[1]+0.8*mean3[2]-0.4*mean3[3]
    aa1<-(0.5*aa1^2+0.25*aa2^2)/n_fam
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      ab2<-sigma1[1]/(sigma1[1]+aa1)
      sigma1[1]<-(swx1[1]+ab2^2*swx1[2])/(n0[8]+ab2*n0[9])
      aa3<-abs(sigma1[1]-aaa0)
      aaa0<-sigma1[1]
      if (n_iter>20) break
    }
    sigma40<-sigma1[1]-sigma;
    if (sigma40<0) {sigma40<-0;sigma1[1]<-sigma}
    sigma1[1]<-sigma40+sigma;sigma1[2]<-sigma1[1]+aa1
    n0[8]<-mix_pi2[1]*n_samB2;n0[9]<-mix_pi2[2]*n_samB2
    aaa0<-sigma2[2];n_iter<-0;aa3<-1000
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      ab3<-sigma2[2]/(sigma2[2]+aa1)
      sigma2[2]<-(ab3^2*swx2[1]+swx2[2])/(ab3*n0[8]+n0[9])
      aa3<-abs(sigma2[2]-aaa0)
      aaa0<-sigma2[2]
      if (n_iter>20) break
    }
    sigma50<-sigma2[2]-sigma;
    if (sigma50<0) {sigma50<-0;sigma2[2]<-sigma}
    sigma2[2]<-sigma50+sigma;sigma2[1]<-sigma2[2]+aa1
    n0[8]<-(mix_pi3[3]+mix_pi3[1])*n_samF2;n0[9]<-mix_pi3[2]*n_samF2
    aaa0<-sigma3[1];n_iter<-0;aa3<-1000
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      ab4<-sigma3[1]/(sigma3[1]+aa1)
      sigma3[1]<-(swx3[1]+ab4^2*swx3[2]+swx3[3])/(n0[8]+ab4*n0[9])
      aa3<-abs(sigma3[1]-aaa0)
      aaa0<-sigma3[1]
      if (n_iter>20) break
    }
    sigma60<-sigma3[1]-sigma;
    if (sigma60<0) {sigma60<-0;sigma3[1]<-sigma}
    sigma3[1]<-sigma60+sigma;sigma3[3]<-sigma3[1];sigma3[2]<-sigma3[1]+aa1
    # iteratively CM3-step for variance (sigma).
    s0[1]<-ss1+ss2+ss3;s0[2]<-n_samP1+n_samF1+n_samP2
    aaa0<-0;n_iter<-0;aa3<-1000
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      abc1<-sigma/(sigma+sigma40)
      abc2<-sigma/(sigma+sigma40+aa1)
      abc3<-sigma/(sigma+sigma50+aa1)
      abc4<-sigma/(sigma+sigma50)
      abc5<-sigma/(sigma+sigma60)
      abc6<-sigma/(sigma+sigma60+aa1)
      aa4<-s0[1]+abc1^2*swx1[1]+abc2^2*swx1[2]+abc3^2*swx2[1]+abc4^2*swx2[2]+abc5^2*(swx3[1]+swx3[3])+abc6^2*swx3[2]
      aa5<-s0[2]+abc1*n0[1]+abc2*n0[2]+abc3*n0[3]+abc4*n0[4]+abc5*(n0[5]+n0[7])+abc6*n0[6]
      sigma<-aa4/aa5
      aa3<-abs(sigma-aaa0)
      aaa0<-sigma
      if (n_iter>20) break
    }
    sigma1[1]<-sigma+sigma40;sigma1[2]<-sigma1[1]+aa1
    sigma2[1]<-sigma+sigma50+aa1;sigma2[2]<-sigma+sigma50
    sigma3[1]<-sigma+sigma60;sigma3[3]<-sigma3[1];sigma3[2]<-sigma3[1]+aa1
    if(sum(sigma < 1e-30)>=1){break}
    ##############################criteria for iterations to stop########################
    L1<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,mean[3],sqrt(sigma))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mix_pi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mix_pi3)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*12
  meanP1<-mean[1];meanF1<-mean[2];meanP2<-mean[3]
  sigma0<-sigma
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*n_samP1)
  P1bmw <- matrix(0,n_samP1,1)
  P1gg <- (dataP1 - mean[1])/sqrt(as.vector(sigma))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < n_samP1){P1bmw <- P1bmw+runif(n_samP1)/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  P1u<- as.matrix(c(12*n_samP1*((P1dd[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((P1dd[2]/n_samP1-1/3)^2),180*n_samP1*((P1dd[3]/n_samP1-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[2])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),P1D))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*n_samF1)
  F1bmw <- matrix(0,n_samF1,1)
  F1gg <- (dataF1 - mean[2])/sqrt(as.vector(sigma))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < n_samF1){F1bmw <- F1bmw+runif(n_samF1)/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  F1u<- as.matrix(c(12*n_samF1*((F1dd[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((F1dd[2]/n_samF1-1/3)^2),180*n_samF1*((F1dd[3]/n_samF1-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[2])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),F1D))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*n_samP2)
  P2bmw <- matrix(0,n_samP2,1)
  P2gg <- (dataP2 - mean[3])/sqrt(as.vector(sigma))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < n_samP2){P2bmw <- P2bmw+runif(n_samP2)/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  P2u<- as.matrix(c(12*n_samP2*((P2dd[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((P2dd[2]/n_samP2-1/3)^2),180*n_samP2*((P2dd[3]/n_samP2-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[2])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),P2D))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1);
  B1w1<-1/(12*n_samB1)
  B1bmw <- matrix(0,n_samB1,1); B1bmwsl <- matrix(0,n_samB1,2)
  for(i in 1:2){
    B1gg <- (dataB1 - mean1[i])/sqrt(sigma1[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi1[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B1P2)))[1]
  if(nn < n_samB1){B1P2 <- B1P2+runif(n_samB1)/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*n_samB1) + sum((B1P2 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  B1u <- as.matrix(c(12*n_samB1*((B1dd[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((B1dd[2]/n_samB1-1/3)^2),180*n_samB1*((B1dd[3]/n_samB1-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[2])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),B1D))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])

  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2);
  B2w1<-1/(12*n_samB2)
  B2bmw <- matrix(0,n_samB2,1); B2bmwsl <- matrix(0,n_samB2,2)
  for(i in 1:2){
    B2gg <- (dataB2 - mean2[i])/sqrt(sigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B2P2)))[1]
  if(nn < n_samB2){B2P2 <- B2P2+runif(n_samB2)/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*n_samB2) + sum((B2P2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  B2u <- as.matrix(c(12*n_samB2*((B2dd[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((B2dd[2]/n_samB2-1/3)^2),180*n_samB2*((B2dd[3]/n_samB2-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[2])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),B2D))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])

  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2);
  F2w1<-1/(12*n_samF2)
  F2bmw <- matrix(0,n_samF2,1); F2bmwsl <- matrix(0,n_samF2,3)
  for(i in 1:3){
    F2gg <- (dataF2 - mean3[i])/sqrt(sigma3[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi3[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < n_samF2){F2P2 <- F2P2+runif(n_samF2)/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*n_samF2) + sum((F2P2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  F2u <- as.matrix(c(12*n_samF2*((F2dd[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((F2dd[2]/n_samF2-1/3)^2),180*n_samF2*((F2dd[3]/n_samF2-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[2])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),F2D))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])

  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)

  ########first order parameters##################
  aa<-matrix(c(1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,-1,0,0,0,0,1,0,0,1,0,
               0,0,0,1,0,0,0,0.5,0,0,0,0,1,0,0,0.5,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,1,0,
               0,0,0,0,0,1,0,0.5,0,0,0,0,0,1,-1,0),10,8,byrow=T)
  mm<-as.matrix(c(mean[1],mean[2],mean[3],mean1[1],mean1[2],mean2[1],mean2[2],mean3[1],mean3[2],mean3[3]))
  B<-solve(crossprod(aa,aa))%*%crossprod(aa,mm)
  #####second order parameters################
  jj1<-sigmaB1-sigma1[1]
  if (jj1<0 || jj1>=sigmaB1) {jj1<-0}
  ll1<-jj1/sigmaB1
  mm1<-sigma1[1]-sigma
  if (mm1<0 || mm1>=sigmaB1) {mm1<-0}
  nn1<-mm1/sigmaB1
  jj2<-sigmaB2-sigma2[2]
  if (jj2<0 || jj2>=sigmaB2) {jj2<-0}
  ll2<-jj2/sigmaB2
  mm2<-sigma2[2]-sigma
  if (mm2<0 || mm2>=sigmaB2) {mm2<-0}
  nn2<-mm2/sigmaB2
  jj3<-sigmaF2-sigma3[1]
  if (jj3<0 || jj3>=sigmaF2) {jj3<-0}
  ll3<-jj3/sigmaF2
  mm3<-sigma3[1]-sigma
  if (mm3<0 || mm3>=sigmaF2) {mm3<-0}
  nn3<-mm3/sigmaF2

  output <- data.frame("MX1-AD-ADI",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4), round(sigma0,4),round(t(mean1),4)," "," ",round(t(sigma1),4)," "," ",
                       round(t(mix_pi1),4)," "," ",round(t(mean2),4)," "," ",round(t(sigma2),4)," "," ",round(t(mix_pi2),4)," "," ",
                       round(t(mean3),4)," "," "," "," "," "," ",round(t(sigma3),4)," "," "," "," "," "," ",round(t(mix_pi3),4)," "," "," "," "," "," ",
                       round(B[1],4),round(B[2],4),round(B[3],4),round(B[4],4),round(B[5],4),round(B[6],4),round(B[7],4)," ",round(B[8],4)," "," "," "," "," "," "," ",
                       round(jj1,4),round(ll1*100,4),round(mm1,4),round(nn1*100,4),round(jj2,4),round(ll2*100,4),round(mm2,4),round(nn2*100,4),round(jj3,4),round(ll3*100,4),round(mm3,4),round(nn3*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi1,mi2,mi3)
  return(OUTPUT)
}

##############MX1-AD-AD(D-1)############################
G6FModelFun[[14]] <- function(K1,logL,df11,df21,df31,df41,df51,df61,G6Ftext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1];n_samF1<-dim(dataF1)[1];n_samP2<-dim(dataP2)[1]
  n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1];n_samF2<-dim(dataF2)[1]
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2),sum(dataF2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2),sum(dataF2^2)))
  ss<-matrix(0,3,1);ss[1]<-s[1]-sumx[1]^2/n_samP1;ss[2]<-s[2]-sumx[2]^2/n_samF1;ss[3]<-s[3]-sumx[3]^2/n_samP2
  mean<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataB1),mean(dataB2),mean(dataF2)))
  sigma0<-sum(ss)/(n_samP1+n_samF1+n_samP2-3)
  sigmaB1<-var(dataB1);sigma40<-sigmaB1
  sigmaB2<-var(dataB2);sigma50<-sigmaB2
  sigmaF2<-var(dataF2);sigma60<-sigmaF2
  sigmaP1<-sigmaF1<-sigmaP2<-sigma0
  m_esp <- 0.0001;n_fam <- as.numeric(G6Ftext2)
  ###############procedure start###########################
  mi1<-matrix(0.5,2,1);sigma1<-matrix(0,2,1)
  mi2<-matrix(0.5,2,1);sigma2<-matrix(0,2,1)
  mi3<-as.matrix(c(0.25,0.5,0.25));sigma3<-matrix(0,3,1)
  sigma<-sigma0
  a1<-sqrt(sigmaB1/n_samB1);if (mean[1]<mean[3]) a1<--a1
  mean1<-as.matrix(c(mean[4]+a1,mean[4]))
  a2<-sqrt(sigmaB2/n_samB2);if (mean[1]<mean[3]) a2<--a2
  mean2<-as.matrix(c(mean[5],mean[5]-a2))
  a3<-sqrt(sigmaF2/n_samF2);if (mean[1]<mean[3]) a3<--a3
  mean3<-as.matrix(c(mean[6]+a3,mean[6],mean[6]-a3))
  b1<-(mean1[1]-mean1[2]+mean2[1]-mean2[2]+2*mean3[1]-2*mean3[3])/6
  # additive effect.
  b2<-(3*mean[1]+2*mean[2]+3*mean[3]-24.5*mean1[1]+30*mean1[2]+30*
         mean2[1]-24.5*mean2[2]-24.5*mean3[1]+30*mean3[2]-24.5*mean3[3])/47
  b3<-(0.5*b1^2+0.25*b2^2)/n_fam
  sigma1[1]<-sigmaB1/2;sigma1[2]<-sigma1[1]+b3
  sigma2[2]<-sigmaB2/2;sigma2[1]<-sigma2[2]+b3
  sigma3[1]<-sigmaF2/2;sigma3[3]<-sigma3[1];sigma3[2]<-sigma3[1]+b3
  L0<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma0))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma0))))+
    sum(log(dnorm(dataP2,mean[3],sqrt(sigma0))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mi1)))+
    sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mi3)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,2,n_samB1);      swx1 <- matrix(0,2,1)
  W2 <- matrix(0,2,n_samB2);      swx2 <- matrix(0,2,1)
  W3 <- matrix(0,3,n_samF2);      swx3 <- matrix(0,3,1)
  n0<-matrix(0,9,1);s0<-matrix(0,9,1)
  hh<-matrix(0,5,5);b_line<-matrix(0,5,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:2) { W1[i,] <- mi1[i]*dnorm(dataB1,mean1[i],sqrt(sigma1[i]))/dmixnorm(dataB1,mean1,sqrt(sigma1),mi1) }
    mix_pi1 <- as.matrix(rowSums(W1)/n_samB1)
    sumwx1 <- W1%*%dataB1
    for(i in 1:2) { W2[i,] <- mi2[i]*dnorm(dataB2,mean2[i],sqrt(sigma2[i]))/dmixnorm(dataB2,mean2,sqrt(sigma2),mi2) }
    mix_pi2 <- as.matrix(rowSums(W2)/n_samB2)
    sumwx2 <- W2%*%dataB2
    for(i in 1:3) { W3[i,] <- mi3[i]*dnorm(dataF2,mean3[i],sqrt(sigma3[i]))/dmixnorm(dataF2,mean3,sqrt(sigma3),mi3) }
    mix_pi3 <- as.matrix(rowSums(W3)/n_samF2)
    sumwx3 <- W3%*%dataF2
    #######obtain means########################
    n0[1]<-mix_pi1[1]*n_samB1;n0[2]<-mix_pi1[2]*n_samB1
    n0[3]<-mix_pi2[1]*n_samB2;n0[4]<-mix_pi2[2]*n_samB2
    n0[5]<-mix_pi3[1]*n_samF2;n0[6]<-mix_pi3[2]*n_samF2
    n0[7]<-mix_pi3[3]*n_samF2
    n0[c(1:7)][abs(n0[c(1:7)])<0.000001]<-0.000001
    AA<-matrix(0,5,1);aaa0<-0;n_iter<-0;aaa1<-1000
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-(mean1[1]-mean1[2]+mean2[1]-mean2[2]+2*mean3[1]-2*mean3[3])/6
      aa2<-(3.0*mean[1]+2*mean[2]+3*mean[3]-24.5*mean1[1]+30*mean1[2]
            +30*mean2[1]-24.5*mean2[2]-24.5*mean3[1]+30*mean3[2]-24.5*mean3[3])/47
      aa1<-(0.5*aa1*aa1+0.25*aa2*aa2)/n_fam
      sigma1[2]<-sigma1[1]+aa1
      sigma2[1]<-sigma2[2]+aa1
      sigma3[2]<-sigma3[1]+aa1
      hh[1,1]<-sigma*(9/n_samP1+4/n_samF1+9/n_samP2)+4*sigma3[1]/n0[5]+16*sigma3[2]/n0[6]+4*sigma3[3]/n0[7]
      hh[1,2]<-3*sigma*(1/n_samP1-1/n_samP2)
      hh[1,3]<-sigma3[1]*(2/n0[5]-2/n0[7])
      hh[1,4]<-sigma3[1]*(2/n0[5]+2/n0[7])
      hh[1,5]<-8*sigma3[2]/n0[6]
      hh[2,2]<-sigma*(1.0/n_samP1+1.0/n_samP2)+sigma1[1]/n0[1]+sigma1[2]/n0[2]+sigma2[1]/n0[3]+sigma2[2]/n0[4]
      hh[2,3]<--sigma1[1]/n0[1]+sigma1[2]/n0[2]+sigma2[1]/n0[3]-sigma2[2]/n0[4]
      hh[2,4]<--sigma1[1]/n0[1]+sigma2[2]/n0[4]
      hh[2,5]<--sigma1[2]/n0[2]+sigma2[1]/n0[3]
      hh[3,3]<-sigma1[1]/n0[1]+sigma1[2]/n0[2]+sigma2[1]/n0[3]+sigma2[2]/n0[4]+sigma3[1]/n0[5]+sigma3[3]/n0[7]
      hh[3,4]<-sigma1[1]/n0[1]-sigma2[2]/n0[4]+sigma3[1]/n0[4]-sigma3[3]/n0[7]
      hh[3,5]<--sigma1[2]/n0[2]+sigma2[1]/n0[3]
      hh[4,4]<-sigma1[1]/n0[1]+sigma2[2]/n0[4]+sigma3[1]/n0[5]+sigma3[3]/n0[7]
      hh[4,5]<-0
      hh[5,5]<-sigma1[2]/n0[2]+sigma2[1]/n0[3]+4*sigma3[2]/n0[6]
      for(i in 2:5)
      {
        for(j in 1:(i-1))
        {
          hh[i,j]<-hh[j,i]
        }
      }
      ##############################################################
      b_line[1]<-3*sumx[1]/n_samP1+2*sumx[2]/n_samF1+3*sumx[3]/n_samP2-2*sumwx3[1]/n0[5]-4*sumwx3[2]/n0[6]-2*sumwx3[3]/n0[7]
      b_line[2]<-sumx[1]/n_samP1-sumx[3]/n_samP2-sumwx1[1]/n0[1]-sumwx1[2]/n0[2]+sumwx2[1]/n0[3]+sumwx2[2]/n0[4]
      b_line[3]<-sumwx1[1]/n0[1]-sumwx1[2]/n0[2]+sumwx2[1]/n0[3]-sumwx2[2]/n0[4]-sumwx3[1]/n0[5]+sumwx3[3]/n0[7]
      b_line[4]<-sumwx1[1]/n0[1]+sumwx2[2]/n0[4]-sumwx3[1]/n0[5]-sumwx3[3]/n0[7]
      b_line[5]<-sumwx1[2]/n0[2]+sumwx2[1]/n0[3]-2*sumwx3[2]/n0[6]
      B<-solve(hh,b_line)
      mean[1]<-(sumx[1]-sigma*(3*B[1]+B[2]))/n_samP1
      mean[2]<-(sumx[2]-2*sigma*B[1])/n_samF1
      mean[3]<-(sumx[3]+sigma*(-3*B[1]+B[2]))/n_samP2
      mean1[1]<-(sumwx1[1]+sigma1[1]*(B[2]-B[3]-B[4]))/n0[1]
      mean1[2]<-(sumwx1[2]+sigma1[2]*(B[2]+B[3]-B[5]))/n0[2]
      mean2[1]<-(sumwx2[1]-sigma2[1]*(B[2]+B[3]+B[5]))/n0[3]
      mean2[2]<-(sumwx2[2]+(-B[2]+B[3]-B[4])*sigma2[2])/n0[4]
      mean3[1]<-(sumwx3[1]+sigma3[1]*(2*B[1]+B[3]+B[4]))/n0[5]
      mean3[2]<-(sumwx3[2]+sigma3[2]*(4*B[1]+2*B[5]))/n0[6]
      mean3[3]<-(sumwx3[3]+sigma3[3]*(2*B[1]-B[3]+B[4]))/n0[7]
      aaa1<-max(abs(B-AA))
      AA<-B
      if (n_iter>20) break
    }
    ########obtain variance#############################
    ss1<-sum((dataP1-mean[1])^2);ss3<-sum((dataP2-mean[3])^2);ss2<-sum((dataF1-mean[2])^2)
    for(i in 1:2) {swx1[i] <- W1[i,]%*%(dataB1-mean1[i])^2 };for(i in 1:2) {swx2[i] <- W2[i,]%*%(dataB2-mean2[i])^2 };for(i in 1:3) {swx3[i] <- W3[i,]%*%(dataF2-mean3[i])^2 }
    aa1<-(mean1[1]-mean1[2]+mean2[1]-mean2[2]+2*mean3[1]-2*mean3[3])/6
    aa2<-(3.0*mean[1]+2*mean[2]+3*mean[3]-24.5*mean1[1]+30*mean1[2]+30*mean2[1]-24.5*mean2[2]-24.5*mean3[1]+30*mean3[2]-24.5*mean3[3])/47
    aa1<-(0.5*aa1^2+0.25*aa2^2)/n_fam
    aaa0<-sigma1[1];n_iter<-0;aa3<-1000
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      ab2<-sigma1[1]/(sigma1[1]+aa1)
      sigma1[1]<-(swx1[1]+ab2^2*swx1[2])/(n0[1]+ab2*n0[2])
      aa3<-abs(sigma1[1]-aaa0)
      aaa0<-sigma1[1]
      if (n_iter>20) break
    }
    sigma40<-sigma1[1]-sigma;
    if (sigma40<0) {sigma40<-0;sigma1[1]<-sigma}
    sigma1[1]<-sigma40+sigma;sigma1[2]<-sigma1[1]+aa1
    n0[8]<-mix_pi2[1]*n_samB2;n0[9]<-mix_pi2[2]*n_samB2
    aaa0<-sigma2[2];n_iter<-0;aa3<-1000
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      ab3<-sigma2[2]/(sigma2[2]+aa1)
      sigma2[2]<-(ab3^2*swx2[1]+swx2[2])/(ab3*n0[8]+n0[9])
      aa3<-abs(sigma2[2]-aaa0)
      aaa0<-sigma2[2]
      if (n_iter>20) break
    }
    sigma50<-sigma2[2]-sigma;
    if (sigma50<0) {sigma50<-0;sigma2[2]<-sigma}
    sigma2[2]<-sigma50+sigma;sigma2[1]<-sigma2[2]+aa1
    n0[8]<-(mix_pi3[3]+mix_pi3[1])*n_samF2;n0[9]<-mix_pi3[2]*n_samF2
    aaa0<-sigma3[1];n_iter<-0;aa3<-1000
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      ab4<-sigma3[1]/(sigma3[1]+aa1)
      sigma3[1]<-(swx3[1]+ab4^2*swx3[2]+swx3[3])/(n0[8]+ab4*n0[9])
      aa3<-abs(sigma3[1]-aaa0)
      aaa0<-sigma3[1]
      if (n_iter>20) break
    }
    sigma60<-sigma3[1]-sigma;
    if (sigma60<0) {sigma60<-0;sigma3[1]<-sigma}
    sigma3[1]<-sigma60+sigma;sigma3[2]<-sigma3[1]+aa1
    sigma3[3]<-sigma3[1]
    # iteratively CM3-step for variance (sigma).
    s0[1]<-ss1+ss2+ss3;s0[2]<-n_samP1+n_samF1+n_samP2
    aaa0<-0;n_iter<-0;aa3<-1000
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      abc1<-sigma/(sigma+sigma40)
      abc2<-sigma/(sigma+sigma40+aa1)
      abc3<-sigma/(sigma+sigma50+aa1)
      abc4<-sigma/(sigma+sigma50)
      abc5<-sigma/(sigma+sigma60)
      abc6<-sigma/(sigma+sigma60+aa1)
      aa4<-s0[1]+abc1^2*swx1[1]+abc2^2*swx1[2]+abc3^2*swx2[1]+abc4^2*swx2[2]+abc5^2*(swx3[1]+swx3[3])+abc6^2*swx3[2]
      aa5<-s0[2]+abc1*n0[1]+abc2*n0[2]+abc3*n0[3]+abc4*n0[4]+abc5*(n0[5]+n0[7])+abc6*n0[6]
      sigma<-aa4/aa5
      aa3<-abs(sigma-aaa0)
      aaa0<-sigma
      if (n_iter>20) break
    }
    sigma1[1]<-sigma+sigma40;sigma1[2]<-sigma1[1]+aa1
    sigma2[1]<-sigma+sigma50+aa1;sigma2[2]<-sigma+sigma50
    sigma3[1]<-sigma+sigma60;sigma3[3]<-sigma3[1];sigma3[2]<-sigma3[1]+aa1
    if(sum(sigma < 1e-30)>=1){break}
    ##############################criteria for iterations to stop########################
    L1<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,mean[3],sqrt(sigma))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mix_pi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mix_pi3)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*9
  meanP1<-mean[1];meanF1<-mean[2];meanP2<-mean[3]
  sigma0<-sigma
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*n_samP1)
  P1bmw <- matrix(0,n_samP1,1)
  P1gg <- (dataP1 - mean[1])/sqrt(as.vector(sigma))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < n_samP1){P1bmw <- P1bmw+runif(n_samP1)/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  P1u<- as.matrix(c(12*n_samP1*((P1dd[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((P1dd[2]/n_samP1-1/3)^2),180*n_samP1*((P1dd[3]/n_samP1-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[2])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),P1D))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*n_samF1)
  F1bmw <- matrix(0,n_samF1,1)
  F1gg <- (dataF1 - mean[2])/sqrt(as.vector(sigma))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < n_samF1){F1bmw <- F1bmw+runif(n_samF1)/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  F1u<- as.matrix(c(12*n_samF1*((F1dd[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((F1dd[2]/n_samF1-1/3)^2),180*n_samF1*((F1dd[3]/n_samF1-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[2])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),F1D))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*n_samP2)
  P2bmw <- matrix(0,n_samP2,1)
  P2gg <- (dataP2 - mean[3])/sqrt(as.vector(sigma))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < n_samP2){P2bmw <- P2bmw+runif(n_samP2)/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  P2u<- as.matrix(c(12*n_samP2*((P2dd[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((P2dd[2]/n_samP2-1/3)^2),180*n_samP2*((P2dd[3]/n_samP2-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[2])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),P2D))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1);
  B1w1<-1/(12*n_samB1)
  B1bmw <- matrix(0,n_samB1,1); B1bmwsl <- matrix(0,n_samB1,2)
  for(i in 1:2){
    B1gg <- (dataB1 - mean1[i])/sqrt(sigma1[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi1[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B1P2)))[1]
  if(nn < n_samB1){B1P2 <- B1P2+runif(n_samB1)/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*n_samB1) + sum((B1P2 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  B1u <- as.matrix(c(12*n_samB1*((B1dd[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((B1dd[2]/n_samB1-1/3)^2),180*n_samB1*((B1dd[3]/n_samB1-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[2])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),B1D))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])

  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2);
  B2w1<-1/(12*n_samB2)
  B2bmw <- matrix(0,n_samB2,1); B2bmwsl <- matrix(0,n_samB2,2)
  for(i in 1:2){
    B2gg <- (dataB2 - mean2[i])/sqrt(sigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B2P2)))[1]
  if(nn < n_samB2){B2P2 <- B2P2+runif(n_samB2)/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*n_samB2) + sum((B2P2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  B2u <- as.matrix(c(12*n_samB2*((B2dd[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((B2dd[2]/n_samB2-1/3)^2),180*n_samB2*((B2dd[3]/n_samB2-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[2])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),B2D))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])

  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2);
  F2w1<-1/(12*n_samF2)
  F2bmw <- matrix(0,n_samF2,1); F2bmwsl <- matrix(0,n_samF2,3)
  for(i in 1:3){
    F2gg <- (dataF2 - mean3[i])/sqrt(sigma3[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi3[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < n_samF2){F2P2 <- F2P2+runif(n_samF2)/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*n_samF2) + sum((F2P2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  F2u <- as.matrix(c(12*n_samF2*((F2dd[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((F2dd[2]/n_samF2-1/3)^2),180*n_samF2*((F2dd[3]/n_samF2-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[2])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),F2D))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])

  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)

  ########first order parameters##########################
  aa<-matrix(c(1,1,0,1,0,1,0,1,0,1,1,-1,0,-1,0,1,1,0,0.5,0.25,
               1,0,0.5,0.5,0.25,1,0,0.5,-0.5,0.25,1,-1,0,-0.5,0.25,
               1,1,0,0,0.25,1,0,0.5,0,0.25,1,-1,0,0,0.25),10,5,byrow=T)
  mm<-as.matrix(c(mean[1],mean[2],mean[3],mean1[1],mean1[2],mean2[1],mean2[2],mean3[1],mean3[2],mean3[3]))
  B1<-solve(crossprod(aa,aa))%*%crossprod(aa,mm)
  #########second order parameters#######################
  jj1<-sigmaB1-sigma1[1]
  if (jj1<0 || jj1>=sigmaB1) {jj1<-0}
  ll1<-jj1/sigmaB1
  mm1<-sigma1[1]-sigma
  if (mm1<0 || mm1>=sigmaB1) {mm1<-0}
  nn1<-mm1/sigmaB1
  jj2<-sigmaB2-sigma2[2]
  if (jj2<0 || jj2>=sigmaB2) {jj2<-0}
  ll2<-jj2/sigmaB2
  mm2<-sigma2[2]-sigma
  if (mm2<0 || mm2>=sigmaB2) {mm2<-0}
  nn2<-mm2/sigmaB2
  jj3<-sigmaF2-sigma3[1]
  if (jj3<0 || jj3>=sigmaF2) {jj3<-0}
  ll3<-jj3/sigmaF2
  mm3<-sigma3[1]-sigma
  if (mm3<0 || mm3>=sigmaF2) {mm3<-0}
  nn3<-mm3/sigmaF2

  output <- data.frame("MX1-AD-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4), round(sigma0,4),round(t(mean1),4)," "," ",round(t(sigma1),4)," "," ",
                       round(t(mix_pi1),4)," "," ",round(t(mean2),4)," "," ",round(t(sigma2),4)," "," ",round(t(mix_pi2),4)," "," ",
                       round(t(mean3),4)," "," "," "," "," "," ",round(t(sigma3),4)," "," "," "," "," "," ",round(t(mix_pi3),4)," "," "," "," "," "," ",
                       round(B1[1],4)," "," "," "," "," ",round(B1[2],4)," ",round(B1[3],4)," "," "," "," "," ",round(B1[4],4),round(B1[5],4),
                       round(jj1,4),round(ll1*100,4),round(mm1,4),round(nn1*100,4),round(jj2,4),round(ll2*100,4),round(mm2,4),round(nn2*100,4),round(jj3,4),round(ll3*100,4),round(mm3,4),round(nn3*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi1,mi2,mi3)
  return(OUTPUT)
}

##################MX1-A-AD(D-2)################################
G6FModelFun[[15]] <- function(K1,logL,df11,df21,df31,df41,df51,df61,G6Ftext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1];n_samF1<-dim(dataF1)[1];n_samP2<-dim(dataP2)[1]
  n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1];n_samF2<-dim(dataF2)[1]
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2),sum(dataF2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2),sum(dataF2^2)))
  ss<-matrix(0,3,1);ss[1]<-s[1]-sumx[1]^2/n_samP1;ss[2]<-s[2]-sumx[2]^2/n_samF1;ss[3]<-s[3]-sumx[3]^2/n_samP2
  mean<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataB1),mean(dataB2),mean(dataF2)))
  sigma0<-sum(ss)/(n_samP1+n_samF1+n_samP2-3)
  sigmaB1<-var(dataB1);sigma40<-sigmaB1
  sigmaB2<-var(dataB2);sigma50<-sigmaB2
  sigmaF2<-var(dataF2);sigma60<-sigmaF2
  sigmaP1<-sigmaF1<-sigmaP2<-sigma0
  m_esp <- 0.0001;n_fam <- as.numeric(G6Ftext2)
  ###############procedure start###########################
  mi1<-matrix(0.5,2,1);sigma1<-matrix(0,2,1)
  mi2<-matrix(0.5,2,1);sigma2<-matrix(0,2,1)
  mi3<-as.matrix(c(0.25,0.5,0.25));sigma3<-matrix(0,3,1)
  sigma<-sigma0
  a1<-sqrt(sigmaB1);if (mean[1]<mean[3]) a1<--a1
  mean1<-as.matrix(c(mean[4]+a1,mean[4]))
  a2<-sqrt(sigmaB2);if (mean[1]<mean[3]) a2<--a2
  mean2<-as.matrix(c(mean[5],mean[5]-a2))
  a3<-sqrt(sigmaF2);if (mean[1]<mean[3]) a3<--a3
  mean3<-as.matrix(c(mean[6]+a3,mean[6],mean[6]-a3))
  b1<-(mean1[1]-mean1[2]+mean2[1]-mean2[2]+2*mean3[1]-2*mean3[3])/6
  # additive effect.
  b3<-(0.5*b1^2)/n_fam
  sigma1[1]<-sigmaB1/2;sigma1[2]<-sigma1[1]+b3
  sigma2[2]<-sigmaB2/2;sigma2[1]<-sigma2[2]+b3
  sigma3[1]<-sigmaF2/2;sigma3[3]<-sigma3[1];sigma3[2]<-sigma3[1]+b3
  L0<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma0))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma0))))+
    sum(log(dnorm(dataP2,mean[3],sqrt(sigma0))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mi1)))+
    sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mi3)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,2,n_samB1);      swx1 <- matrix(0,2,1)
  W2 <- matrix(0,2,n_samB2);      swx2 <- matrix(0,2,1)
  W3 <- matrix(0,3,n_samF2);      swx3 <- matrix(0,3,1)
  hh<-matrix(0,6,6);b_line<-matrix(0,6,1)
  n0<-matrix(0,9,1);s0<-matrix(0,2,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:2) { W1[i,] <- mi1[i]*dnorm(dataB1,mean1[i],sqrt(sigma1[i]))/dmixnorm(dataB1,mean1,sqrt(sigma1),mi1) }
    mix_pi1 <- as.matrix(rowSums(W1)/n_samB1)
    sumwx1 <- W1%*%dataB1
    for(i in 1:2) { W2[i,] <- mi2[i]*dnorm(dataB2,mean2[i],sqrt(sigma2[i]))/dmixnorm(dataB2,mean2,sqrt(sigma2),mi2) }
    mix_pi2 <- as.matrix(rowSums(W2)/n_samB2)
    sumwx2 <- W2%*%dataB2
    for(i in 1:3) { W3[i,] <- mi3[i]*dnorm(dataF2,mean3[i],sqrt(sigma3[i]))/dmixnorm(dataF2,mean3,sqrt(sigma3),mi3) }
    mix_pi3 <- as.matrix(rowSums(W3)/n_samF2)
    sumwx3 <- W3%*%dataF2
    ########obtain means#####################
    n0[1]<-mix_pi1[1]*n_samB1;n0[2]<-mix_pi1[2]*n_samB1
    n0[3]<-mix_pi2[1]*n_samB2;n0[4]<-mix_pi2[2]*n_samB2
    n0[5]<-mix_pi3[1]*n_samF2;n0[6]<-mix_pi3[2]*n_samF2
    n0[7]<-mix_pi3[3]*n_samF2
    n0[c(1:7)][abs(n0[c(1:7)])<0.0001]<-0.0001
    aaa0<-0    # CM1-step for means.
    AA<-matrix(0,6,1);n_iter<-0;aaa1<-1000
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-(mean1[1]-mean1[2]+mean2[1]-mean2[2]+2*mean3[1]-2*mean3[3])/6
      aa1<-0.5*aa1^2/n_fam
      sigma1[2]<-sigma1[1]+aa1
      sigma2[1]<-sigma2[2]+aa1
      sigma3[2]<-sigma3[1]+aa1
      hh[1,1]<-sigma*(9/n_samP1+4/n_samF1+9/n_samP2)+64*sigma3[2]/n0[6]
      hh[1,2]<-3*sigma*(1/n_samP1-1/n_samP2)
      hh[1,3]<-16*sigma3[2]/n0[6]
      hh[1,4]<-0
      hh[1,5]<-16*sigma3[2]/n0[6]
      hh[1,6]<-16*sigma3[2]/n0[6]
      hh[2,2]<-sigma*(1/n_samP1+1/n_samP2)+sigma1[1]/n0[1]+sigma1[2]/n0[2]+sigma2[1]/n0[3]+sigma2[2]/n0[4]
      hh[2,3]<--sigma1[1]/n0[1]+sigma2[2]/n0[4]
      hh[2,4]<--sigma1[1]/n0[1]+sigma1[2]/n0[2]+sigma2[1]/n0[3]-sigma2[2]/n0[4]
      hh[2,5]<--sigma1[2]/n0[2]+sigma2[1]/n0[3]
      hh[2,6]<-0
      hh[3,3]<-sigma1[1]/n0[1]+sigma2[2]/n0[4]+4*sigma3[2]/n0[6]
      hh[3,4]<-sigma1[1]/n0[1]-sigma2[2]/n0[4]
      hh[3,5]<-4*sigma3[2]/n0[6]
      hh[3,6]<-4*sigma3[2]/n0[6]
      hh[4,4]<-sigma1[1]/n0[1]+sigma1[2]/n0[2]+sigma2[1]/n0[3]+sigma2[2]/n0[4]+sigma3[1]/n0[5]+sigma3[3]/n0[7]
      hh[4,5]<--sigma1[2]/n0[2]+sigma2[1]/n0[3]
      hh[4,6]<--sigma3[1]/n0[5]+sigma3[3]/n0[7]
      hh[5,5]<-sigma1[2]/n0[2]+sigma2[1]/n0[3]+4*sigma3[2]/n0[6]
      hh[5,6]<-4*sigma3[2]/n0[6]
      hh[6,6]<-sigma3[1]/n0[5]+sigma3[3]/n0[7]+4*sigma3[2]/n0[6]
      for(i in 2:6)
      {
        for(j in 1:(i-1))
        {
          hh[i,j]<-hh[j,i]
        }
      }
      ##############################################################
      b_line[1]<-3*sumx[1]/n_samP1+2*sumx[2]/n_samF1+3*sumx[3]/n_samP2-8*sumwx3[2]/n0[6]
      b_line[2]<-sumx[1]/n_samP1-sumx[3]/n_samP2-sumwx1[1]/n0[1]-sumwx1[2]/n0[2]+sumwx2[1]/n0[3]+sumwx2[2]/n0[4]
      b_line[3]<-sumwx1[1]/n0[1]+sumwx2[2]/n0[4]-2*sumwx3[2]/n0[6]
      b_line[4]<-sumwx1[1]/n0[1]-sumwx1[2]/n0[2]+sumwx2[1]/n0[3]-sumwx2[2]/n0[4]-sumwx3[1]/n0[5]+sumwx3[3]/n0[7]
      b_line[5]<-sumwx1[2]/n0[2]+sumwx2[1]/n0[3]-2*sumwx3[2]/n0[6]
      b_line[6]<-sumwx3[1]/n0[5]+sumwx3[3]/n0[7]-2*sumwx3[2]/n0[6]
      B<-solve(hh,b_line)
      # to solve the restricted condition functions.
      # to estimate the means.
      mean[1]<-(sumx[1]-sigma*(3*B[1]+B[2]))/n_samP1
      mean[2]<-(sumx[2]-2*sigma*B[1])/n_samF1
      mean[3]<-(sumx[3]+sigma*(-3*B[1]+B[2]))/n_samP2
      mean1[1]<-(sumwx1[1]+sigma1[1]*(B[2]-B[3]-B[4]))/n0[1]
      mean1[2]<-(sumwx1[2]+sigma1[2]*(B[2]+B[4]-B[5]))/n0[2]
      mean2[1]<-(sumwx2[1]-sigma2[1]*(B[2]+B[4]+B[5]))/n0[3]
      mean2[2]<-(sumwx2[2]+(-B[2]-B[3]+B[4])*sigma2[2])/n0[4]
      mean3[1]<-(sumwx3[1]+sigma3[1]*(B[4]-B[6]))/n0[5]
      mean3[2]<-(sumwx3[2]+sigma3[2]*(8*B[1]+2*B[3]+2*B[5]+2*B[6]))/n0[6]
      mean3[3]<-(sumwx3[3]-sigma3[3]*(B[4]+B[6]))/n0[7]
      aaa1<-max(abs(B-AA))
      AA<-B
      if (n_iter>20) break
    }
    ##########obtain variance###########################
    ss1<-sum((dataP1-mean[1])^2);ss3<-sum((dataP2-mean[3])^2);ss2<-sum((dataF1-mean[2])^2)
    for(i in 1:2) {swx1[i] <- W1[i,]%*%(dataB1-mean1[i])^2 };for(i in 1:2) {swx2[i] <- W2[i,]%*%(dataB2-mean2[i])^2 };for(i in 1:3) {swx3[i] <- W3[i,]%*%(dataF2-mean3[i])^2 }
    aa1<-(mean1[1]-mean1[2]+mean2[1]-mean2[2]+2*mean3[1]-2*mean3[3])/6
    aa1<-0.5*aa1*aa1/n_fam
    aaa0<-sigma1[1];n_iter<-0;aa3<-1000
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      ab2<-sigma1[1]/(sigma1[1]+aa1)
      sigma1[1]<-(swx1[1]+ab2^2*swx1[2])/(n0[1]+ab2*n0[2])
      aa3<-abs(sigma1[1]-aaa0)
      aaa0<-sigma1[1]
      if (n_iter>20) break
    }
    sigma40<-sigma1[1]-sigma;
    if (sigma40<0) {sigma40<-0;sigma1[1]<-sigma}
    sigma1[1]<-sigma40+sigma;sigma1[2]<-sigma1[1]+aa1
    n0[8]<-mix_pi2[1]*n_samB2;n0[9]<-mix_pi2[2]*n_samB2
    aa3<-1000;aaa0<-sigma2[2];n_iter<-0
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      ab3<-sigma2[2]/(sigma2[2]+aa1)
      sigma2[2]<-(ab3^2*swx2[1]+swx2[2])/(ab3*n0[8]+n0[9])
      aa3<-abs(sigma2[2]-aaa0)
      aaa0<-sigma2[2]
      if (n_iter>20) break
    }
    sigma50<-sigma2[2]-sigma;
    if (sigma50<0) {sigma50<-0;sigma2[2]<-sigma}
    sigma2[2]<-sigma50+sigma;sigma2[1]<-sigma2[2]+aa1
    n0[8]<-(mix_pi3[3]+mix_pi3[1])*n_samF2;n0[9]<-mix_pi3[2]*n_samF2
    aaa0<-sigma3[1];n_iter<-0;aa3<-1000
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      ab4<-sigma3[1]/(sigma3[1]+aa1)
      sigma3[1]<-(swx3[1]+ab4^2*swx3[2]+swx3[3])/(n0[8]+ab4*n0[9])
      aa3<-abs(sigma3[1]-aaa0)
      aaa0<-sigma3[1]
      if (n_iter>20) break
    }
    sigma60<-sigma3[1]-sigma;
    if (sigma60<0) {sigma60<-0;sigma3[1]<-sigma}
    sigma3[1]<-sigma60+sigma;sigma3[2]<-sigma3[1]+aa1
    sigma3[3]<-sigma3[1]
    # iteratively CM3-step for variance (sigma).
    s0[1]<-ss1+ss2+ss3;s0[2]<-n_samP1+n_samF1+n_samP2
    aaa0<-0;n_iter<-0;aa3<-1000
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      abc1<-sigma/(sigma+sigma40)
      abc2<-sigma/(sigma+sigma40+aa1)
      abc3<-sigma/(sigma+sigma50+aa1)
      abc4<-sigma/(sigma+sigma50)
      abc5<-sigma/(sigma+sigma60)
      abc6<-sigma/(sigma+sigma60+aa1)
      aa4<-s0[1]+abc1^2*swx1[1]+abc2^2*swx1[2]+abc3^2*
        swx2[1]+abc4^2*swx2[2]+abc5^2*(swx3[1]+swx3[3])+abc6^2*swx3[2]
      aa5<-s0[2]+abc1*n0[1]+abc2*n0[2]+abc3*n0[3]+abc4*n0[4]+abc5*(n0[5]+n0[7])+abc6*n0[6]
      sigma<-aa4/aa5
      aa3<-abs(sigma-aaa0)
      aaa0<-sigma
      if (n_iter>20) break
    }
    sigma1[1]<-sigma+sigma40;sigma1[2]<-sigma1[1]+aa1
    sigma2[1]<-sigma+sigma50+aa1;sigma2[2]<-sigma+sigma50
    sigma3[1]<-sigma+sigma60;sigma3[3]<-sigma3[1];sigma3[2]<-sigma3[1]+aa1
    if(sum(sigma < 1e-30)>=1){break}
    ##############################criteria for iterations to stop########################
    L1<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,mean[3],sqrt(sigma))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mix_pi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mix_pi3)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*8
  meanP1<-mean[1];meanF1<-mean[2];meanP2<-mean[3]
  sigma0<-sigma
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*n_samP1)
  P1bmw <- matrix(0,n_samP1,1)
  P1gg <- (dataP1 - mean[1])/sqrt(as.vector(sigma))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < n_samP1){P1bmw <- P1bmw+runif(n_samP1)/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  P1u<- as.matrix(c(12*n_samP1*((P1dd[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((P1dd[2]/n_samP1-1/3)^2),180*n_samP1*((P1dd[3]/n_samP1-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[2])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),P1D))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*n_samF1)
  F1bmw <- matrix(0,n_samF1,1)
  F1gg <- (dataF1 - mean[2])/sqrt(as.vector(sigma))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < n_samF1){F1bmw <- F1bmw+runif(n_samF1)/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  F1u<- as.matrix(c(12*n_samF1*((F1dd[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((F1dd[2]/n_samF1-1/3)^2),180*n_samF1*((F1dd[3]/n_samF1-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[2])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),F1D))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*n_samP2)
  P2bmw <- matrix(0,n_samP2,1)
  P2gg <- (dataP2 - mean[3])/sqrt(as.vector(sigma))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < n_samP2){P2bmw <- P2bmw+runif(n_samP2)/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  P2u<- as.matrix(c(12*n_samP2*((P2dd[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((P2dd[2]/n_samP2-1/3)^2),180*n_samP2*((P2dd[3]/n_samP2-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[2])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),P2D))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1);
  B1w1<-1/(12*n_samB1)
  B1bmw <- matrix(0,n_samB1,1); B1bmwsl <- matrix(0,n_samB1,2)
  for(i in 1:2){
    B1gg <- (dataB1 - mean1[i])/sqrt(sigma1[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi1[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B1P2)))[1]
  if(nn < n_samB1){B1P2 <- B1P2+runif(n_samB1)/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*n_samB1) + sum((B1P2 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  B1u <- as.matrix(c(12*n_samB1*((B1dd[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((B1dd[2]/n_samB1-1/3)^2),180*n_samB1*((B1dd[3]/n_samB1-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[2])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),B1D))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])

  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2);
  B2w1<-1/(12*n_samB2)
  B2bmw <- matrix(0,n_samB2,1); B2bmwsl <- matrix(0,n_samB2,2)
  for(i in 1:2){
    B2gg <- (dataB2 - mean2[i])/sqrt(sigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B2P2)))[1]
  if(nn < n_samB2){B2P2 <- B2P2+runif(n_samB2)/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*n_samB2) + sum((B2P2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  B2u <- as.matrix(c(12*n_samB2*((B2dd[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((B2dd[2]/n_samB2-1/3)^2),180*n_samB2*((B2dd[3]/n_samB2-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[2])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),B2D))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])

  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2);
  F2w1<-1/(12*n_samF2)
  F2bmw <- matrix(0,n_samF2,1); F2bmwsl <- matrix(0,n_samF2,3)
  for(i in 1:3){
    F2gg <- (dataF2 - mean3[i])/sqrt(sigma3[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi3[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < n_samF2){F2P2 <- F2P2+runif(n_samF2)/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*n_samF2) + sum((F2P2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  F2u <- as.matrix(c(12*n_samF2*((F2dd[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((F2dd[2]/n_samF2-1/3)^2),180*n_samF2*((F2dd[3]/n_samF2-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[2])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),F2D))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])

  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)

  #######first order parameters#############################
  aa<-matrix(c(1,1,1,0,1,0,0,1,1,-1,-1,0,1,1,0.5,0.25,1,0,0.5,0.25,
               1,0,-0.5,0.25,1,-1,-0.5,0.25,1,1,0,0.25,1,0,0,0.25,
               1,-1,0,0.25),10,4,byrow=T)
  mm<-as.matrix(c(mean[1],mean[2],mean[3],mean1[1],mean1[2],mean2[1],mean2[2],mean3[1],mean3[2],mean3[3]))
  B1<-solve(crossprod(aa,aa))%*%crossprod(aa,mm)
  ########second order parameters#############################
  jj1<-sigmaB1-sigma1[1]
  if (jj1<0 || jj1>=sigmaB1) {jj1<-0}
  ll1<-jj1/sigmaB1
  mm1<-sigma1[1]-sigma
  if (mm1<0 || mm1>=sigmaB1) {mm1<-0}
  nn1<-mm1/sigmaB1
  jj2<-sigmaB2-sigma2[2]
  if (jj2<0 || jj2>=sigmaB2) {jj2<-0}
  ll2<-jj2/sigmaB2
  mm2<-sigma2[2]-sigma
  if (mm2<0 || mm2>=sigmaB2) {mm2<-0}
  nn2<-mm2/sigmaB2
  jj3<-sigmaF2-sigma3[1]
  if (jj3<0 || jj3>=sigmaF2) {jj3<-0}
  ll3<-jj3/sigmaF2
  mm3<-sigma3[1]-sigma
  if (mm3<0 || mm3>=sigmaF2) {mm3<-0}
  nn3<-mm3/sigmaF2

  output <- data.frame("MX1-A-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4), round(sigma0,4),round(t(mean1),4)," "," ",round(t(sigma1),4)," "," ",
                       round(t(mix_pi1),4)," "," ",round(t(mean2),4)," "," ",round(t(sigma2),4)," "," ",round(t(mix_pi2),4)," "," ",
                       round(t(mean3),4)," "," "," "," "," "," ",round(t(sigma3),4)," "," "," "," "," "," ",round(t(mix_pi3),4)," "," "," "," "," "," ",
                       round(B1[1],4)," "," "," "," "," ",round(B1[2],4)," "," "," "," "," "," "," ",round(B1[3],4),round(B1[4],4),
                       round(jj1,4),round(ll1*100,4),round(mm1,4),round(nn1*100,4),round(jj2,4),round(ll2*100,4),round(mm2,4),round(nn2*100,4),round(jj3,4),round(ll3*100,4),round(mm3,4),round(nn3*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi1,mi2,mi3)
  return(OUTPUT)
}

#################MX1-EAD-AD(D-3)############################
G6FModelFun[[16]] <- function(K1,logL,df11,df21,df31,df41,df51,df61,G6Ftext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1];n_samF1<-dim(dataF1)[1];n_samP2<-dim(dataP2)[1]
  n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1];n_samF2<-dim(dataF2)[1]
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2),sum(dataF2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2),sum(dataF2^2)))
  ss<-matrix(0,3,1);ss[1]<-s[1]-sumx[1]^2/n_samP1;ss[2]<-s[2]-sumx[2]^2/n_samF1;ss[3]<-s[3]-sumx[3]^2/n_samP2
  mean<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataB1),mean(dataB2),mean(dataF2)))
  sigma0<-sum(ss)/(n_samP1+n_samF1+n_samP2-3)
  sigmaB1<-var(dataB1);sigma40<-sigmaB1
  sigmaB2<-var(dataB2);sigma50<-sigmaB2
  sigmaF2<-var(dataF2);sigma60<-sigmaF2
  sigmaP1<-sigmaF1<-sigmaP2<-sigma0
  m_esp <- 0.0001;n_fam <- as.numeric(G6Ftext2)
  ###############procedure start###########################
  mi1<-matrix(0.5,2,1);sigma1<-matrix(0,2,1)
  mi2<-matrix(0.5,2,1);sigma2<-matrix(0,2,1)
  mi3<-as.matrix(c(0.25,0.5,0.25));sigma3<-matrix(0,3,1)
  sigma<-sigma0
  a1<-sqrt(sigmaB1/n_samB1);if (mean[1]<mean[3]) a1<--a1
  mean1<-as.matrix(c(mean[4]+0.5*a1,0.25*(mean[1]+mean[3])+0.5*mean[2]))
  a2<-sqrt(sigmaB2/n_samB2);if (mean[1]<mean[3]) a2<--a2
  mean2<-as.matrix(c(mean[5],mean[5]-2*a2))
  a3<-sqrt(sigmaF2/n_samF2);if (mean[1]<mean[3]) a3<--a3
  mean3<-as.matrix(c(mean1[1],mean[6],mean2[2]))
  b1<-(12*mean[1]+8*mean[2]+12*mean[3]+120*mean1[1]-98*mean1[2]+338*mean2[1]-316*mean2[2]+338*mean3[1]+120*mean3[2]-534*mean3[3])/1496
  # additive effect.
  b3<-0.75*b1^2/n_fam
  sigma1[1]<-sigmaB1/2;sigma1[2]<-sigma1[1]+2*b3
  sigma2[2]<-sigmaB2/2;sigma2[1]<-sigma2[2]+2*b3
  sigma3[1]<-sigmaF2/2;sigma3[3]<-sigma3[1];sigma3[2]<-sigma3[1]+2*b3
  L0<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma0))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma0))))+
    sum(log(dnorm(dataP2,mean[3],sqrt(sigma0))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mi1)))+
    sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mi3)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,2,n_samB1);      swx1 <- matrix(0,2,1)
  W2 <- matrix(0,2,n_samB2);      swx2 <- matrix(0,2,1)
  W3 <- matrix(0,3,n_samF2);      swx3 <- matrix(0,3,1)
  hh<-matrix(0,6,6);b_line<-matrix(0,6,1)
  n0<-matrix(0,9,1);s0<-matrix(0,2,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:2) { W1[i,] <- mi1[i]*dnorm(dataB1,mean1[i],sqrt(sigma1[i]))/dmixnorm(dataB1,mean1,sqrt(sigma1),mi1) }
    mix_pi1 <- as.matrix(rowSums(W1)/n_samB1)
    sumwx1 <- W1%*%dataB1
    for(i in 1:2) { W2[i,] <- mi2[i]*dnorm(dataB2,mean2[i],sqrt(sigma2[i]))/dmixnorm(dataB2,mean2,sqrt(sigma2),mi2) }
    mix_pi2 <- as.matrix(rowSums(W2)/n_samB2)
    sumwx2 <- W2%*%dataB2
    for(i in 1:3) { W3[i,] <- mi3[i]*dnorm(dataF2,mean3[i],sqrt(sigma3[i]))/dmixnorm(dataF2,mean3,sqrt(sigma3),mi3) }
    mix_pi3 <- as.matrix(rowSums(W3)/n_samF2)
    sumwx3 <- W3%*%dataF2
    #######obtain means############################
    n0[1]<-mix_pi1[1]*n_samB1;n0[2]<-mix_pi1[2]*n_samB1
    n0[3]<-mix_pi2[1]*n_samB2;n0[4]<-mix_pi2[2]*n_samB2
    n0[5]<-mix_pi3[1]*n_samF2;n0[6]<-mix_pi3[2]*n_samF2
    n0[7]<-mix_pi3[3]*n_samF2
    n0[c(1:7)][abs(n0[c(1:7)])<0.00000001]<-0.000001
    aaa0<-0 ;AA<-matrix(0,6,1);aaa1<-1000;n_iter<-0
    while(aaa1>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-(12*mean[1]+8*mean[2]+12*mean[3]+120*mean1[1]-98*mean1[2]+
              338*mean2[1]-316*mean2[2]+338*mean3[1]+120*mean3[2]-534*mean3[3])/1496
      aa1<-0.75*aa1^2/n_fam
      sigma1[2]<-sigma1[1]+aa1
      sigma2[1]<-sigma2[2]+aa1
      sigma3[2]<-sigma3[1]+aa1
      hh[1,1]<-sigma*(9/n_samP1+4/n_samF1+9/n_samP2)+4*sigma3[1]/n0[5]+16*sigma3[2]/n0[6]+4*sigma3[3]/n0[7]
      hh[1,2]<-sigma*(3/n_samP1-3/n_samP2)
      hh[1,3]<-sigma3[1]*(2/n0[5]-2/n0[7])
      hh[1,4]<-sigma3[1]*(2/n0[5]+2/n0[7])
      hh[1,5]<-8*sigma3[2]/n0[6]
      hh[1,6]<-6*sigma3[1]/n0[5]-16*sigma3[2]/n0[6]+2*sigma3[3]/n0[7]
      hh[2,2]<-sigma*(1/n_samP1+1/n_samP2)+sigma1[1]/n0[1]+sigma1[2]/n0[2]+sigma2[1]/n0[3]+sigma2[2]/n0[4]
      hh[2,3]<--sigma1[1]/n0[1]+sigma1[2]/n0[2]+sigma2[1]/n0[3]-sigma2[2]/n0[4]
      hh[2,4]<--sigma1[1]/n0[1]+sigma2[2]/n0[4]
      hh[2,5]<--sigma1[2]/n0[2]+sigma2[1]/n0[3]
      hh[2,6]<-0
      hh[3,3]<-sigma1[1]/n0[1]+sigma1[2]/n0[2]+sigma2[1]/n0[3]+sigma2[2]/n0[4]+sigma3[1]/n0[5]+sigma3[3]/n0[7]
      hh[3,4]<-sigma1[1]/n0[1]-sigma2[2]/n0[4]+sigma3[1]/n0[5]-sigma3[3]/n0[7]
      hh[3,5]<--sigma1[2]/n0[2]+sigma2[1]/n0[3]
      hh[3,6]<-sigma3[1]*(3/n0[5]-1/n0[7])
      hh[4,4]<-sigma1[1]/n0[1]+sigma2[2]/n0[4]+sigma3[1]/n0[5]+sigma3[3]/n0[7]
      hh[4,5]<-0
      hh[4,6]<-3*sigma3[1]/n0[5]+sigma3[3]/n0[7]
      hh[5,5]<-sigma1[2]/n0[2]+sigma2[1]/n0[3]+4*sigma3[2]/n0[6]
      hh[5,6]<--8*sigma3[2]/n0[6]
      hh[6,6]<-9*sigma3[1]/n0[5]+sigma3[3]/n0[7]+16*sigma3[2]/n0[6]
      for(i in 2:6)
      {
        for(j in 1:(i-1))
        {
          hh[i,j]<-hh[j,i]
        }
      }
      ##############################################################
      b_line[1]<-3*sumx[1]/n_samP1+2*sumx[2]/n_samF1+3*sumx[3]/n_samP2-2*sumwx3[1]/n0[5]-4*sumwx3[2]/n0[6]-2*sumwx3[3]/n0[7]
      b_line[2]<-sumx[1]/n_samP1-sumx[3]/n_samP2-sumwx1[1]/n0[1]-sumwx1[2]/n0[2]+sumwx2[1]/n0[3]+sumwx2[2]/n0[4]
      b_line[3]<-sumwx1[1]/n0[1]-sumwx1[2]/n0[2]+sumwx2[1]/n0[3]-sumwx2[2]/n0[4]-sumwx3[1]/n0[5]+sumwx3[3]/n0[7]
      b_line[4]<-sumwx1[1]/n0[1]+sumwx2[2]/n0[4]-sumwx3[1]/n0[5]-sumwx3[3]/n0[7]
      b_line[5]<-sumwx1[2]/n0[2]+sumwx2[1]/n0[3]-2*sumwx3[2]/n0[6]
      b_line[6]<-4*sumwx3[2]/n0[6]-3*sumwx3[1]/n0[5]-sumwx3[3]/n0[7]
      ##############################################################
      B<-solve(hh,b_line)
      mean[1]<-(sumx[1]-sigma*(3*B[1]+B[2]))/n_samP1
      mean[2]<-(sumx[2]-2*sigma*B[1])/n_samF1
      mean[3]<-(sumx[3]+sigma*(-3*B[1]+B[2]))/n_samP2
      mean1[1]<-(sumwx1[1]+sigma1[1]*(B[2]-B[3]-B[4]))/n0[1]
      mean1[2]<-(sumwx1[2]+sigma1[2]*(B[2]+B[3]-B[5]))/n0[2]
      mean2[1]<-(sumwx2[1]-sigma2[1]*(B[2]+B[3]+B[5]))/n0[3]
      mean2[2]<-(sumwx2[2]+(-B[2]+B[3]-B[4])*sigma2[2])/n0[4]
      mean3[1]<-(sumwx3[1]+sigma3[1]*(2*B[1]+B[3]+B[4]+3*B[6]))/n0[5]
      mean3[2]<-(sumwx3[2]+sigma3[2]*(4*B[1]+2*B[5]-4*B[6]))/n0[6]
      mean3[3]<-(sumwx3[3]+sigma3[1]*(2*B[1]-B[3]+B[4]+B[6]))/n0[7]
      aaa1<-max(abs(B-AA))
      AA<-B
      if(n_iter>20)break
    }
    ##############obtain variance###########################
    ss1<-sum((dataP1-mean[1])^2);ss3<-sum((dataP2-mean[3])^2);ss2<-sum((dataF1-mean[2])^2)
    for(i in 1:2) {swx1[i] <- W1[i,]%*%(dataB1-mean1[i])^2 };for(i in 1:2) {swx2[i] <- W2[i,]%*%(dataB2-mean2[i])^2 };for(i in 1:3) {swx3[i] <- W3[i,]%*%(dataF2-mean3[i])^2 }
    aa1<-(12*mean[1]+8*mean[2]+12*mean[3]+120*mean1[1]-98*mean1[2]+
            338*mean2[1]-316*mean2[2]+338*mean3[1]+120*mean3[2]-534*mean3[3])/1496
    aa1<-0.75*aa1^2/n_fam
    aaa0<-sigma1[1];aa3<-1000;n_iter<-0
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      ab2<-sigma1[1]/(sigma1[1]+aa1)
      sigma1[1]<-(swx1[1]+ab2^2*swx1[2])/(n0[1]+ab2*n0[2])
      aa3<-abs(sigma1[1]-aaa0)
      aaa0<-sigma1[1]
      if(n_iter>20)break
    }
    sigma40<-sigma1[1]-sigma;
    if (sigma40<0) {sigma40<-0;sigma1[1]<-sigma}
    sigma1[1]<-sigma40+sigma;sigma1[2]<-sigma1[1]+aa1
    n0[8]<-mix_pi2[1]*n_samB2;n0[9]<-mix_pi2[2]*n_samB2
    aaa0<-sigma2[2]
    aa3<-1000;n_iter<-0
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      ab3<-sigma2[2]/(sigma2[2]+aa1)
      sigma2[2]<-(ab3*ab3*swx2[1]+swx2[2])/(ab3*n0[8]+n0[9])
      aa3<-abs(sigma2[2]-aaa0)
      aaa0<-sigma2[2]
      if(n_iter>20)break
    }
    sigma50<-sigma2[2]-sigma;
    if (sigma50<0) {sigma50<-0;sigma2[2]<-sigma}
    sigma2[2]<-sigma50+sigma;sigma2[1]<-sigma2[2]+aa1
    n0[8]<-(mix_pi3[3]+mix_pi3[1])*n_samF2;n0[9]<-mix_pi3[2]*n_samF2
    aaa0<-sigma3[1];aa3<-1000;n_iter<-0
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      ab4<-sigma3[1]/(sigma3[1]+aa1)
      sigma3[1]<-(swx3[1]+ab4^2*swx3[2]+swx3[3])/(n0[8]+ab4*n0[9])
      aa3<-abs(sigma3[1]-aaa0)
      aaa0<-sigma3[1]
      if(n_iter>20)break
    }
    sigma60<-sigma3[1]-sigma;
    if (sigma60<0) {sigma60<-0;sigma3[1]<-sigma}
    sigma3[1]<-sigma60+sigma;sigma3[2]<-sigma3[1]+aa1
    sigma3[3]<-sigma3[1]
    # iteratively CM3-step for variance (sigma).
    s0[1]<-ss1+ss2+ss3;s0[2]<-n_samP1+n_samF1+n_samP2
    aaa0<-0;aa3<-1000;n_iter<-0
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      abc1<-sigma/(sigma+sigma40)
      abc2<-sigma/(sigma+sigma40+aa1)
      abc3<-sigma/(sigma+sigma50+aa1)
      abc4<-sigma/(sigma+sigma50)
      abc5<-sigma/(sigma+sigma60)
      abc6<-sigma/(sigma+sigma60+aa1)
      aa4<-s0[1]+abc1^2*swx1[1]+abc2^2*swx1[2]+abc3^2*swx2[1]+abc4^2*swx2[2]+abc5^2*(swx3[1]+swx3[3])+abc6^2*swx3[2]
      aa5<-s0[2]+abc1*n0[1]+abc2*n0[2]+abc3*n0[3]+abc4*n0[4]+abc5*(n0[5]+n0[7])+abc6*n0[6]
      sigma<-aa4/aa5
      aa3<-abs(sigma-aaa0)
      aaa0<-sigma
      if(n_iter>20)break
    }
    sigma1[1]<-sigma+sigma40;sigma1[2]<-sigma1[1]+aa1
    sigma2[1]<-sigma+sigma50+aa1;sigma2[2]<-sigma+sigma50
    sigma3[1]<-sigma+sigma60;sigma3[3]<-sigma3[1];sigma3[2]<-sigma3[1]+aa1
    if(sum(sigma < 1e-30)>=1){break}
    ##############################criteria for iterations to stop########################
    L1<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,mean[3],sqrt(sigma))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mix_pi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mix_pi3)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc<-L0
  AIC<--2*abc+2*8
  meanP1<-mean[1];meanF1<-mean[2];meanP2<-mean[3]
  sigma0<-sigma
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*n_samP1)
  P1bmw <- matrix(0,n_samP1,1)
  P1gg <- (dataP1 - mean[1])/sqrt(as.vector(sigma))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < n_samP1){P1bmw <- P1bmw+runif(n_samP1)/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  P1u<- as.matrix(c(12*n_samP1*((P1dd[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((P1dd[2]/n_samP1-1/3)^2),180*n_samP1*((P1dd[3]/n_samP1-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[2])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),P1D))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*n_samF1)
  F1bmw <- matrix(0,n_samF1,1)
  F1gg <- (dataF1 - mean[2])/sqrt(as.vector(sigma))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < n_samF1){F1bmw <- F1bmw+runif(n_samF1)/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  F1u<- as.matrix(c(12*n_samF1*((F1dd[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((F1dd[2]/n_samF1-1/3)^2),180*n_samF1*((F1dd[3]/n_samF1-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[2])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),F1D))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*n_samP2)
  P2bmw <- matrix(0,n_samP2,1)
  P2gg <- (dataP2 - mean[3])/sqrt(as.vector(sigma))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < n_samP2){P2bmw <- P2bmw+runif(n_samP2)/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  P2u<- as.matrix(c(12*n_samP2*((P2dd[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((P2dd[2]/n_samP2-1/3)^2),180*n_samP2*((P2dd[3]/n_samP2-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[2])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),P2D))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1);
  B1w1<-1/(12*n_samB1)
  B1bmw <- matrix(0,n_samB1,1); B1bmwsl <- matrix(0,n_samB1,2)
  for(i in 1:2){
    B1gg <- (dataB1 - mean1[i])/sqrt(sigma1[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi1[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B1P2)))[1]
  if(nn < n_samB1){B1P2 <- B1P2+runif(n_samB1)/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*n_samB1) + sum((B1P2 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  B1u <- as.matrix(c(12*n_samB1*((B1dd[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((B1dd[2]/n_samB1-1/3)^2),180*n_samB1*((B1dd[3]/n_samB1-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[2])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),B1D))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])

  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2);
  B2w1<-1/(12*n_samB2)
  B2bmw <- matrix(0,n_samB2,1); B2bmwsl <- matrix(0,n_samB2,2)
  for(i in 1:2){
    B2gg <- (dataB2 - mean2[i])/sqrt(sigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B2P2)))[1]
  if(nn < n_samB2){B2P2 <- B2P2+runif(n_samB2)/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*n_samB2) + sum((B2P2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  B2u <- as.matrix(c(12*n_samB2*((B2dd[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((B2dd[2]/n_samB2-1/3)^2),180*n_samB2*((B2dd[3]/n_samB2-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[2])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),B2D))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])

  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2);
  F2w1<-1/(12*n_samF2)
  F2bmw <- matrix(0,n_samF2,1); F2bmwsl <- matrix(0,n_samF2,3)
  for(i in 1:3){
    F2gg <- (dataF2 - mean3[i])/sqrt(sigma3[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi3[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < n_samF2){F2P2 <- F2P2+runif(n_samF2)/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*n_samF2) + sum((F2P2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  F2u <- as.matrix(c(12*n_samF2*((F2dd[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((F2dd[2]/n_samF2-1/3)^2),180*n_samF2*((F2dd[3]/n_samF2-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[2])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),F2D))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])

  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)

  #############first order parameters############################
  aa<-matrix(c(1,1,1,0,1,1,0,1,1,-1,-1,0,1,1,0.5,0.25,1,0.5,0.5,0.25,
               1,0.5,-0.5,0.25,1,-1,-0.5,0.25,1,1,0,0.25,1,0.5,0,0.25,
               1,-1,0,0.25),10,4,byrow=T)
  mm<-as.matrix(c(mean[1],mean[2],mean[3],mean1[1],mean1[2],mean2[1],mean2[2],mean3[1],mean3[2],mean3[3]))
  B1<-solve(crossprod(aa,aa))%*%crossprod(aa,mm)
  ##########second order parameters################################
  jj1<-sigmaB1-sigma1[1]
  if (jj1<0 || jj1>=sigmaB1) {jj1<-0}
  ll1<-jj1/sigmaB1
  mm1<-sigma1[1]-sigma
  if (mm1<0 || mm1>=sigmaB1) {mm1<-0}
  nn1<-mm1/sigmaB1
  jj2<-sigmaB2-sigma2[2]
  if (jj2<0 || jj2>=sigmaB2) {jj2<-0}
  ll2<-jj2/sigmaB2
  mm2<-sigma2[2]-sigma
  if (mm2<0 || mm2>=sigmaB2) {mm2<-0}
  nn2<-mm2/sigmaB2
  jj3<-sigmaF2-sigma3[1]
  if (jj3<0 || jj3>=sigmaF2) {jj3<-0}
  ll3<-jj3/sigmaF2
  mm3<-sigma3[1]-sigma
  if (mm3<0 || mm3>=sigmaF2) {mm3<-0}
  nn3<-mm3/sigmaF2

  output <- data.frame("MX1-EAD-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4), round(sigma0,4),round(t(mean1),4)," "," ",round(t(sigma1),4)," "," ",
                       round(t(mix_pi1),4)," "," ",round(t(mean2),4)," "," ",round(t(sigma2),4)," "," ",round(t(mix_pi2),4)," "," ",
                       round(t(mean3),4)," "," "," "," "," "," ",round(t(sigma3),4)," "," "," "," "," "," ",round(t(mix_pi3),4)," "," "," "," "," "," ",
                       round(B1[1],4)," "," "," "," "," ",round(B1[2],4)," ",round(B1[2],4)," "," "," "," "," ",round(B1[3],4),round(B1[4],4),
                       round(jj1,4),round(ll1*100,4),round(mm1,4),round(nn1*100,4),round(jj2,4),round(ll2*100,4),round(mm2,4),round(nn2*100,4),round(jj3,4),round(ll3*100,4),round(mm3,4),round(nn3*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi1,mi2,mi3)
  return(OUTPUT)
}

################MX1-NCD-AD(D-4)###############################
G6FModelFun[[17]] <- function(K1,logL,df11,df21,df31,df41,df51,df61,G6Ftext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1];n_samF1<-dim(dataF1)[1];n_samP2<-dim(dataP2)[1]
  n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1];n_samF2<-dim(dataF2)[1]
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2),sum(dataF2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2),sum(dataF2^2)))
  ss<-matrix(0,3,1);ss[1]<-s[1]-sumx[1]^2/n_samP1;ss[2]<-s[2]-sumx[2]^2/n_samF1;ss[3]<-s[3]-sumx[3]^2/n_samP2
  mean<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataB1),mean(dataB2),mean(dataF2)))
  sigma0<-sum(ss)/(n_samP1+n_samF1+n_samP2-3)
  sigmaB1<-var(dataB1);sigma40<-sigmaB1
  sigmaB2<-var(dataB2);sigma50<-sigmaB2
  sigmaF2<-var(dataF2);sigma60<-sigmaF2
  sigmaP1<-sigmaF1<-sigmaP2<-sigma0
  m_esp <- 0.0001;n_fam <- as.numeric(G6Ftext2)
  ###############procedure start###########################
  mi1<-matrix(0.5,2,1);sigma1<-matrix(0,2,1)
  mi2<-matrix(0.5,2,1);sigma2<-matrix(0,2,1)
  mi3<-as.matrix(c(0.25,0.5,0.25));sigma3<-matrix(0,3,1)
  sigma<-sigma0
  a1<-sqrt(sigmaB1/n_samB1);if (mean[1]<mean[3]) a1<--a1
  mean1<-as.matrix(c(mean[4]+0.5*a1,0.25*(mean[1]+mean[3])+0.5*mean[2]))
  a2<-sqrt(sigmaB2/n_samB2);if (mean[1]<mean[3]) a2<--a2
  mean2<-as.matrix(c(mean[5],mean[5]-2*a2))
  a3<-sqrt(sigmaF2/n_samF2);if (mean[1]<mean[3]) a3<--a3
  mean3<-as.matrix(c(mean1[1],mean[6],mean2[2]))
  b1<-(-12*mean[1]-8*mean[2]-12*mean[3]+316*mean1[1]-338*mean1[2]+98*mean2[1]-120*mean2[2]+534*mean3[1]-120*mean3[2]-338*mean3[3])/1496
  # additive effect.
  b3<-0.75*b1^2/n_fam
  sigma1[1]<-sigmaB1/2;sigma1[2]<-sigma1[1]+2*b3
  sigma2[2]<-sigmaB2/2;sigma2[1]<-sigma2[2]+2*b3
  sigma3[1]<-sigmaF2/2;sigma3[3]<-sigma3[1];sigma3[2]<-sigma3[1]+2*b3
  L0<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma0))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma0))))+
    sum(log(dnorm(dataP2,mean[3],sqrt(sigma0))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mi1)))+
    sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mi3)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,2,n_samB1);      swx1 <- matrix(0,2,1)
  W2 <- matrix(0,2,n_samB2);      swx2 <- matrix(0,2,1)
  W3 <- matrix(0,3,n_samF2);      swx3 <- matrix(0,3,1)
  hh<-matrix(0,6,6);b_line<-matrix(0,6,1)
  n0<-matrix(0,9,1);s0<-matrix(0,2,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:2) { W1[i,] <- mi1[i]*dnorm(dataB1,mean1[i],sqrt(sigma1[i]))/dmixnorm(dataB1,mean1,sqrt(sigma1),mi1) }
    mix_pi1 <- as.matrix(rowSums(W1)/n_samB1)
    sumwx1 <- W1%*%dataB1
    for(i in 1:2) { W2[i,] <- mi2[i]*dnorm(dataB2,mean2[i],sqrt(sigma2[i]))/dmixnorm(dataB2,mean2,sqrt(sigma2),mi2) }
    mix_pi2 <- as.matrix(rowSums(W2)/n_samB2)
    sumwx2 <- W2%*%dataB2
    for(i in 1:3) { W3[i,] <- mi3[i]*dnorm(dataF2,mean3[i],sqrt(sigma3[i]))/dmixnorm(dataF2,mean3,sqrt(sigma3),mi3) }
    mix_pi3 <- as.matrix(rowSums(W3)/n_samF2)
    sumwx3 <- W3%*%dataF2
    ########obtain means###################
    n0[1]<-mix_pi1[1]*n_samB1;n0[2]<-mix_pi1[2]*n_samB1
    n0[3]<-mix_pi2[1]*n_samB2;n0[4]<-mix_pi2[2]*n_samB2
    n0[5]<-mix_pi3[1]*n_samF2;n0[6]<-mix_pi3[2]*n_samF2
    n0[7]<-mix_pi3[3]*n_samF2
    n0[c(1:7)][abs(n0[c(1:7)])<0.00000001]<-0.000001
    aaa0<-0 ;AA<-matrix(0,6,1);aaa1<-1000;n_iter<-0
    while(aaa1>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-(-12*mean[1]-8*mean[2]-12*mean[3]+316*mean1[1]-338*mean1[2]+98*mean2[1]-120*mean2[2]+534*mean3[1]-120*mean3[2]-338*mean3[3])/1496
      aa1<-0.75*aa1^2/n_fam
      sigma1[2]<-sigma1[1]+aa1
      sigma2[1]<-sigma2[2]+aa1
      sigma3[2]<-sigma3[1]+aa1
      hh[1,1]<-sigma*(9/n_samP1+4/n_samF1+9/n_samP2)+4*sigma3[1]/n0[5]+16*sigma3[2]/n0[6]+4*sigma3[3]/n0[7]
      hh[1,2]<-sigma*(3/n_samP1-3/n_samP2)
      hh[1,3]<-sigma3[1]*(2/n0[5]-2/n0[7])
      hh[1,4]<-sigma3[1]*(2/n0[5]+2/n0[7])
      hh[1,5]<-8*sigma3[2]/n0[6]
      hh[1,6]<-2*sigma3[1]/n0[5]-16*sigma3[2]/n0[6]+6*sigma3[3]/n0[7]
      hh[2,2]<-sigma*(1/n_samP1+1/n_samP2)+sigma1[1]/n0[1]+sigma1[2]/n0[2]+sigma2[1]/n0[3]+sigma2[2]/n0[4]
      hh[2,3]<--sigma1[1]/n0[1]+sigma1[2]/n0[2]+sigma2[1]/n0[3]-sigma2[2]/n0[4]
      hh[2,4]<--sigma1[1]/n0[1]+sigma2[2]/n0[4]
      hh[2,5]<--sigma1[2]/n0[2]+sigma2[1]/n0[3]
      hh[2,6]<-0
      hh[3,3]<-sigma1[1]/n0[1]+sigma1[2]/n0[2]+sigma2[1]/n0[3]+sigma2[2]/n0[4]+sigma3[1]/n0[5]+sigma3[3]/n0[7]
      hh[3,4]<-sigma1[1]/n0[1]-sigma2[2]/n0[4]+sigma3[1]/n0[5]-sigma3[3]/n0[7]
      hh[3,5]<--sigma1[2]/n0[2]+sigma2[1]/n0[3]
      hh[3,6]<-sigma3[1]*(1/n0[5]-3/n0[7])
      hh[4,4]<-sigma1[1]/n0[1]+sigma2[2]/n0[4]+sigma3[1]/n0[5]+sigma3[3]/n0[7]
      hh[4,5]<-0
      hh[4,6]<-sigma3[1]/n0[5]+3*sigma3[3]/n0[7]
      hh[5,5]<-sigma1[2]/n0[2]+sigma2[1]/n0[3]+4*sigma3[2]/n0[6]
      hh[5,6]<--8*sigma3[2]/n0[6]
      hh[6,6]<-sigma3[1]/n0[5]+16*sigma3[2]/n0[6]+9*sigma3[3]/n0[7]
      for(i in 2:6)
      {
        for(j in 1:(i-1))
        {
          hh[i,j]<-hh[j,i]
        }
      }
      ########################################################
      b_line[1]<-3*sumx[1]/n_samP1+2*sumx[2]/n_samF1+3*sumx[3]/n_samP2-2*sumwx3[1]/n0[5]-4*sumwx3[2]/n0[6]-2*sumwx3[3]/n0[7]
      b_line[2]<-sumx[1]/n_samP1-sumx[3]/n_samP2-sumwx1[1]/n0[1]-sumwx1[2]/n0[2]+sumwx2[1]/n0[3]+sumwx2[2]/n0[4]
      b_line[3]<-sumwx1[1]/n0[1]-sumwx1[2]/n0[2]+sumwx2[1]/n0[3]-sumwx2[2]/n0[4]-sumwx3[1]/n0[5]+sumwx3[3]/n0[7]
      b_line[4]<-sumwx1[1]/n0[1]+sumwx2[2]/n0[4]-sumwx3[1]/n0[5]-sumwx3[3]/n0[7]
      b_line[5]<-sumwx1[2]/n0[2]+sumwx2[1]/n0[3]-2*sumwx3[2]/n0[6]
      b_line[6]<-4*sumwx3[2]/n0[6]-sumwx3[1]/n0[5]-3*sumwx3[3]/n0[7]
      B<-solve(hh,b_line)
      # to solve the restricted condition functions.
      # to estimate the means.
      mean[1]<-(sumx[1]-sigma*(3*B[1]+B[2]))/n_samP1
      mean[2]<-(sumx[2]-2*sigma*B[1])/n_samF1
      mean[3]<-(sumx[3]+sigma*(-3*B[1]+B[2]))/n_samP2
      mean1[1]<-(sumwx1[1]+sigma1[1]*(B[2]-B[3]-B[4]))/n0[1]
      mean1[2]<-(sumwx1[2]+sigma1[2]*(B[2]+B[3]-B[5]))/n0[2]
      mean2[1]<-(sumwx2[1]-sigma2[1]*(B[2]+B[3]+B[5]))/n0[3]
      mean2[2]<-(sumwx2[2]+(-B[2]+B[3]-B[4])*sigma2[2])/n0[4]
      mean3[1]<-(sumwx3[1]+sigma3[1]*(2*B[1]+B[3]+B[4]+B[6]))/n0[5]
      mean3[2]<-(sumwx3[2]+sigma3[2]*(4*B[1]+2*B[5]-4*B[6]))/n0[6]
      mean3[3]<-(sumwx3[3]+sigma3[3]*(2*B[1]-B[3]+B[4]+3*B[6]))/n0[7]
      aaa1<-max(abs(B-AA))
      AA<-B
      if(n_iter>20)break
    }
    ##########################################################################
    ss1<-sum((dataP1-mean[1])^2);ss3<-sum((dataP2-mean[3])^2);ss2<-sum((dataF1-mean[2])^2)
    for(i in 1:2) {swx1[i] <- W1[i,]%*%(dataB1-mean1[i])^2 };for(i in 1:2) {swx2[i] <- W2[i,]%*%(dataB2-mean2[i])^2 };for(i in 1:3) {swx3[i] <- W3[i,]%*%(dataF2-mean3[i])^2 }
    aa1<-(-12*mean[1]-8*mean[2]-12*mean[3]+316*mean1[1]-338*mean1[2]+98*mean2[1]-120*mean2[2]+534*mean3[1]-120*mean3[2]-338*mean3[3])/1496
    aa1<-0.75*aa1^2/n_fam
    aaa0<-sigma1[1];aa3<-1000;n_iter<-0
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      ab2<-sigma1[1]/(sigma1[1]+aa1)
      sigma1[1]<-(swx1[1]+ab2^2*swx1[2])/(n0[1]+ab2*n0[2])
      aa3<-abs(sigma1[1]-aaa0)
      aaa0<-sigma1[1]
      if(n_iter>20)break
    }
    sigma40<-sigma1[1]-sigma;
    if (sigma40<0) {sigma40<-0;sigma1[1]<-sigma}
    sigma1[1]<-sigma40+sigma;sigma1[2]<-sigma1[1]+aa1
    n0[8]<-mix_pi2[1]*n_samB2;n0[9]<-mix_pi2[2]*n_samB2
    aaa0<-sigma2[2];aa3<-1000;n_iter<-0
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      ab3<-sigma2[2]/(sigma2[2]+aa1)
      sigma2[2]<-(ab3^2*swx2[1]+swx2[2])/(ab3*n0[8]+n0[9])
      aa3<-abs(sigma2[2]-aaa0)
      aaa0<-sigma2[2]
      if(n_iter>20)break
    }
    sigma50<-sigma2[2]-sigma;
    if (sigma50<0) {sigma50<-0;sigma2[2]<-sigma}
    sigma2[2]<-sigma50+sigma;sigma2[1]<-sigma2[2]+aa1
    n0[8]<-(mix_pi3[3]+mix_pi3[1])*n_samF2;n0[9]<-mix_pi3[2]*n_samF2
    aaa0<-sigma3[1];aa3<-1000;n_iter<-0
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      ab4<-sigma3[1]/(sigma3[1]+aa1)
      sigma3[1]<-(swx3[1]+ab4^2*swx3[2]+swx3[3])/(n0[8]+ab4*n0[9])
      aa3<-abs(sigma3[1]-aaa0)
      aaa0<-sigma3[1]
      if(n_iter>20)break
    }
    sigma60<-sigma3[1]-sigma;
    if (sigma60<0) {sigma60<-0;sigma3[1]<-sigma}
    sigma3[1]<-sigma60+sigma;sigma3[2]<-sigma3[1]+aa1;sigma3[3]<-sigma3[1]
    # iteratively CM3-step for variance (sigma).
    s0[1]<-ss1+ss2+ss3;s0[2]<-n_samP1+n_samF1+n_samP2
    aaa0<-0;aa3<-1000;n_iter<-0
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      abc1<-sigma/(sigma+sigma40)
      abc2<-sigma/(sigma+sigma40+aa1)
      abc3<-sigma/(sigma+sigma50+aa1)
      abc4<-sigma/(sigma+sigma50)
      abc5<-sigma/(sigma+sigma60)
      abc6<-sigma/(sigma+sigma60+aa1)
      aa4<-s0[1]+abc1^2*swx1[1]+abc2^2*swx1[2]+abc3^2*swx2[1]+abc4^2*swx2[2]+abc5^2*(swx3[1]+swx3[3])+abc6^2*swx3[2]
      aa5<-s0[2]+abc1*n0[1]+abc2*n0[2]+abc3*n0[3]+abc4*n0[4]+abc5*(n0[5]+n0[7])+abc6*n0[6]
      sigma<-aa4/aa5
      aa3<-abs(sigma-aaa0)
      aaa0<-sigma
      if(n_iter>20)break
    }
    sigma1[1]<-sigma+sigma40;sigma1[2]<-sigma1[1]+aa1
    sigma2[1]<-sigma+sigma50+aa1;sigma2[2]<-sigma+sigma50
    sigma3[1]<-sigma+sigma60;sigma3[3]<-sigma3[1];sigma3[2]<-sigma3[1]+aa1
    if(sum(sigma < 1e-30)>=1){break}
    ##############################criteria for iterations to stop########################
    L1<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,mean[3],sqrt(sigma))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mix_pi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mix_pi3)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*8
  meanP1<-mean[1];meanF1<-mean[2];meanP2<-mean[3]
  sigma0<-sigma
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*n_samP1)
  P1bmw <- matrix(0,n_samP1,1)
  P1gg <- (dataP1 - mean[1])/sqrt(as.vector(sigma))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < n_samP1){P1bmw <- P1bmw+runif(n_samP1)/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  P1u<- as.matrix(c(12*n_samP1*((P1dd[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((P1dd[2]/n_samP1-1/3)^2),180*n_samP1*((P1dd[3]/n_samP1-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[2])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),P1D))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*n_samF1)
  F1bmw <- matrix(0,n_samF1,1)
  F1gg <- (dataF1 - mean[2])/sqrt(as.vector(sigma))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < n_samF1){F1bmw <- F1bmw+runif(n_samF1)/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  F1u<- as.matrix(c(12*n_samF1*((F1dd[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((F1dd[2]/n_samF1-1/3)^2),180*n_samF1*((F1dd[3]/n_samF1-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[2])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),F1D))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*n_samP2)
  P2bmw <- matrix(0,n_samP2,1)
  P2gg <- (dataP2 - mean[3])/sqrt(as.vector(sigma))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < n_samP2){P2bmw <- P2bmw+runif(n_samP2)/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  P2u<- as.matrix(c(12*n_samP2*((P2dd[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((P2dd[2]/n_samP2-1/3)^2),180*n_samP2*((P2dd[3]/n_samP2-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[2])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),P2D))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1);
  B1w1<-1/(12*n_samB1)
  B1bmw <- matrix(0,n_samB1,1); B1bmwsl <- matrix(0,n_samB1,2)
  for(i in 1:2){
    B1gg <- (dataB1 - mean1[i])/sqrt(sigma1[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi1[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B1P2)))[1]
  if(nn < n_samB1){B1P2 <- B1P2+runif(n_samB1)/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*n_samB1) + sum((B1P2 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  B1u <- as.matrix(c(12*n_samB1*((B1dd[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((B1dd[2]/n_samB1-1/3)^2),180*n_samB1*((B1dd[3]/n_samB1-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[2])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),B1D))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])

  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2);
  B2w1<-1/(12*n_samB2)
  B2bmw <- matrix(0,n_samB2,1); B2bmwsl <- matrix(0,n_samB2,2)
  for(i in 1:2){
    B2gg <- (dataB2 - mean2[i])/sqrt(sigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B2P2)))[1]
  if(nn < n_samB2){B2P2 <- B2P2+runif(n_samB2)/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*n_samB2) + sum((B2P2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  B2u <- as.matrix(c(12*n_samB2*((B2dd[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((B2dd[2]/n_samB2-1/3)^2),180*n_samB2*((B2dd[3]/n_samB2-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[2])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),B2D))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])

  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2);
  F2w1<-1/(12*n_samF2)
  F2bmw <- matrix(0,n_samF2,1); F2bmwsl <- matrix(0,n_samF2,3)
  for(i in 1:3){
    F2gg <- (dataF2 - mean3[i])/sqrt(sigma3[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi3[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < n_samF2){F2P2 <- F2P2+runif(n_samF2)/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*n_samF2) + sum((F2P2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  F2u <- as.matrix(c(12*n_samF2*((F2dd[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((F2dd[2]/n_samF2-1/3)^2),180*n_samF2*((F2dd[3]/n_samF2-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[2])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),F2D))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])

  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)

  ########first order parameters##############################
  aa<-matrix(c(1,1,1,0,1,-1,0,1,1,-1,-1,0,1,1,0.5,0.25,1,-0.5,0.5,0.25,
               1,-0.5,-0.5,0.25,1,-1,-0.5,0.25,1,1,0,0.25,1,-0.5,0,0.25,
               1,-1,0,0.25),10,4,byrow=T)
  mm<-as.matrix(c(mean[1],mean[2],mean[3],mean1[1],mean1[2],mean2[1],mean2[2],mean3[1],mean3[2],mean3[3]))
  B1<-solve(crossprod(aa,aa))%*%crossprod(aa,mm)
  ########second order parameters##############################
  jj1<-sigmaB1-sigma1[1]
  if (jj1<0 || jj1>=sigmaB1) {jj1<-0}
  ll1<-jj1/sigmaB1
  mm1<-sigma1[1]-sigma
  if (mm1<0 || mm1>=sigmaB1) {mm1<-0}
  nn1<-mm1/sigmaB1
  jj2<-sigmaB2-sigma2[2]
  if (jj2<0 || jj2>=sigmaB2) {jj2<-0}
  ll2<-jj2/sigmaB2
  mm2<-sigma2[2]-sigma
  if (mm2<0 || mm2>=sigmaB2) {mm2<-0}
  nn2<-mm2/sigmaB2
  jj3<-sigmaF2-sigma3[1]
  if (jj3<0 || jj3>=sigmaF2) {jj3<-0}
  ll3<-jj3/sigmaF2
  mm3<-sigma3[1]-sigma
  if (mm3<0 || mm3>=sigmaF2) {mm3<-0}
  nn3<-mm3/sigmaF2

  output <- data.frame("MX1-NCD-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4), round(sigma0,4),round(t(mean1),4)," "," ",round(t(sigma1),4)," "," ",
                       round(t(mix_pi1),4)," "," ",round(t(mean2),4)," "," ",round(t(sigma2),4)," "," ",round(t(mix_pi2),4)," "," ",
                       round(t(mean3),4)," "," "," "," "," "," ",round(t(sigma3),4)," "," "," "," "," "," ",round(t(mix_pi3),4)," "," "," "," "," "," ",
                       round(B1[1],4)," "," "," "," "," ",round(B1[2],4)," ",round(-B1[2],4)," "," "," "," "," ",round(B1[3],4),round(B1[4],4),
                       round(jj1,4),round(ll1*100,4),round(mm1,4),round(nn1*100,4),round(jj2,4),round(ll2*100,4),round(mm2,4),round(nn2*100,4),round(jj3,4),round(ll3*100,4),round(mm3,4),round(nn3*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi1,mi2,mi3)
  return(OUTPUT)
}

##############MX2-ADI-ADI(E-0)#######################
G6FModelFun[[18]] <- function(K1,logL,df11,df21,df31,df41,df51,df61,G6Ftext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1];n_samF1<-dim(dataF1)[1];n_samP2<-dim(dataP2)[1]
  n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1];n_samF2<-dim(dataF2)[1]
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2),sum(dataF2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2),sum(dataF2^2)))
  ss<-matrix(0,3,1);ss[1]<-s[1]-sumx[1]^2/n_samP1;ss[2]<-s[2]-sumx[2]^2/n_samF1;ss[3]<-s[3]-sumx[3]^2/n_samP2
  mean<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataB1),mean(dataB2),mean(dataF2)))
  sigma0<-sum(ss)/(n_samP1+n_samF1+n_samP2-3)
  sigmaB1<-var(dataB1);sigma40<-sigmaB1
  sigmaB2<-var(dataB2);sigma50<-sigmaB2
  sigmaF2<-var(dataF2);sigma60<-sigmaF2
  sigmaP1<-sigmaF1<-sigmaP2<-sigma0
  m_esp <- 0.0001;n_fam <- as.numeric(G6Ftext2)
  ###############procedure start###########################
  mi1<-matrix(0.25,4,1);sigma1<-matrix(0,4,1)
  mi2<-matrix(0.25,4,1);sigma2<-matrix(0,4,1)
  mi3<-as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  sigma3<-matrix(0,9,1)
  sigma<-sigma0
  a1<-sqrt(sigma40/n_samB1);if (mean[1]<mean[3]) a1<--a1
  mean1<-as.matrix(c(mean[4]+2*a1,mean[4]+a1,mean[4]-a1,mean[4]-2*a1))
  a2<-sqrt(sigma50/n_samB2);if (mean[1]<mean[3]) a2<--a2
  mean2<-as.matrix(c(mean[5]+2*a2,mean[5]+a2,mean[5]-a2,mean[5]-2*a2))
  a3<-sqrt(sigma60/n_samF2);if (mean[1]<mean[3]) a3<--a3
  mean3<-as.matrix(c(mean[6]+2.8*a3,mean[6]+2.1*a3,mean[6]+1.4*a3,mean[6]+0.7*a3,mean[6],mean[6]-0.7*a3,mean[6]-1.4*a3,mean[6]-2.1*a3,mean[6]-2.8*a3))
  sigma1[1]<-sigmaB1/2;sigma2[4]<-sigmaB2/2;sigma3[1]<-sigmaF2/2
  aa<-matrix(c(1,0,0,0,0,0,1,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,1,0,0,1,0,0,0,-1,-1,0,0,
               1,0,0,0,0,0,0,1,0,0,1,1,0,0,1,0,0,0,0,0,0,1,0,0,1,0,0,0.5,0,0.5,0,0,0,0,0,1,0,0,
               0,1,0.5,0,0,0,0.5,0,0,0,0,1,0,0,0,0,0.5,0.5,0,0,0,0.25,0,0,0,0,1,0,0,0,0.5,0.5,0,
               0,0,0.25,0,0,0,0,1,0,0,-1,0.5,0,0,0,-0.5,0,0,0,0,0,1,0,-1,0,0,0.5,0,-0.5,0,0,0,0,0,
               0,1,0,-1,-1,0,0,1,0,0,0,0,0,0,0,0,1,1,1,0,0,1,0,0,0,0,0,0,0,0,1,1,0,0,0.5,0,0.5,0,0,
               0,0,0,0,0,1,1,-1,0,0,-1,0,0,0,0,0,0,0,0,1,0,1,0.5,0,0,0,0.5,0,0,0,0,0,0,1,0,0,0.5,0.5,
               0,0,0,0.25,0,0,0,0,0,1,0,-1,0.5,0,0,0,-0.5,0,0,0,0,0,0,1,-1,1,0,0,-1,0,0,0,0,0,0,0,0,1,
               -1,0,0,0.5,0,-0.5,0,0,0,0,0,0,0,1,-1,-1,0,0,1,0,0,0),20,14,byrow=T)
  mm<-as.matrix(c(mean[1],mean[2],mean[3],mean1[1],mean1[2],mean1[3],mean1[4],mean2[1],mean2[2],mean2[3],mean2[4],
                  mean3[1],mean3[2],mean3[3],mean3[4],mean3[5],mean3[6],mean3[7],mean3[8],mean3[9]))
  B<-solve(crossprod(aa,aa))%*%crossprod(aa,mm)
  gs<-matrix(0,8,1)
  gs[1]<-B[7];gs[2]<-B[8];gs[3]<-B[9];gs[4]<-B[10];gs[5]<-B[11];gs[6]<-B[12];gs[7]<-B[13];gs[8]<-B[14]
  g_aa1<-(0.5*(gs[2]+gs[5])^2+0.25*(gs[4]+gs[6])^2)/n_fam
  #   0.5(db+i)**2+0.25(hb+jab)**2.
  g_aa2<-(0.5*(gs[1]+gs[5])^2+0.25*(gs[3]+gs[7])^2)/n_fam
  #   0.5(da+i)**2+0.25(ha+jba)**2.
  g_aa3<-(0.5*(gs[1]-gs[5])^2+0.25*(gs[3]-gs[7])^2)/n_fam
  #   0.5(da-i)**2+0.25(ha-jba)**2.
  g_aa4<-(0.5*(gs[2]-gs[5])^2+0.25*(gs[4]-gs[6])^2)/n_fam
  #   0.5(db-i)**2+0.25(hb-jab)**2.
  g_aa5<-0.25*(gs[1]^2+gs[2]^2+gs[5]^2+(gs[1]+gs[6])^2+(gs[2]+gs[7])^2+(gs[3]+gs[8]/2)^2+(gs[4]+gs[8]/2)^2+gs[8]^2/4)/n_fam
  sigma1[2]<-sigma1[1]+g_aa1;sigma1[3]<-sigma1[1]+g_aa2;sigma1[4]<-sigma1[1]+g_aa5
  sigma2[1]<-sigma2[4]+g_aa5;sigma2[2]<-sigma2[4]+g_aa3;sigma2[3]<-sigma2[4]+g_aa4
  sigma3[2]<-sigma3[1]+g_aa1;sigma3[3]<-sigma3[1];sigma3[4]<-sigma3[1]+g_aa2
  sigma3[5]<-sigma3[1]+g_aa5;sigma3[6]<-sigma3[1]+g_aa3;sigma3[7]<-sigma3[1]
  sigma3[8]<-sigma3[1]+g_aa4;sigma3[9]<-sigma3[1]
  L0<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma0))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma0))))+
    sum(log(dnorm(dataP2,mean[3],sqrt(sigma0))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mi1)))+
    sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mi3)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,4,n_samB1);      swx1 <- matrix(0,4,1)
  W2 <- matrix(0,4,n_samB2);      swx2 <- matrix(0,4,1)
  W3 <- matrix(0,9,n_samF2);      swx3 <- matrix(0,9,1)
  hh<-matrix(0,6,6); b_line1<-matrix(0,20,1);b_line2<-matrix(0,6,1)
  n0<-matrix(0,18,1);s0<-matrix(0,18,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:4) { W1[i,] <- mi1[i]*dnorm(dataB1,mean1[i],sqrt(sigma1[i]))/dmixnorm(dataB1,mean1,sqrt(sigma1),mi1) }
    mix_pi1 <- as.matrix(rowSums(W1)/n_samB1)
    sumwx1 <- W1%*%dataB1
    for(i in 1:4) { W2[i,] <- mi2[i]*dnorm(dataB2,mean2[i],sqrt(sigma2[i]))/dmixnorm(dataB2,mean2,sqrt(sigma2),mi2) }
    mix_pi2 <- as.matrix(rowSums(W2)/n_samB2)
    sumwx2 <- W2%*%dataB2
    for(i in 1:9) { W3[i,] <- mi3[i]*dnorm(dataF2,mean3[i],sqrt(sigma3[i]))/dmixnorm(dataF2,mean3,sqrt(sigma3),mi3) }
    mix_pi3 <- as.matrix(rowSums(W3)/n_samF2)
    sumwx3 <- W3%*%dataF2
    ###########obtain means################
    n0[1]<-mix_pi1[1]*n_samB1;n0[2]<-mix_pi1[2]*n_samB1;n0[3]<-mix_pi1[3]*n_samB1
    n0[4]<-mix_pi1[4]*n_samB1;n0[5]<-mix_pi2[1]*n_samB2;n0[6]<-mix_pi2[2]*n_samB2
    n0[7]<-mix_pi2[3]*n_samB2;n0[8]<-mix_pi2[4]*n_samB2
    s0[1]<-mix_pi3[1]*n_samF2;s0[2]<-mix_pi3[2]*n_samF2;s0[3]<-mix_pi3[3]*n_samF2
    s0[4]<-mix_pi3[4]*n_samF2;s0[5]<-mix_pi3[5]*n_samF2;s0[6]<-mix_pi3[6]*n_samF2
    s0[7]<-mix_pi3[7]*n_samF2;s0[8]<-mix_pi3[8]*n_samF2;s0[9]<-mix_pi3[9]*n_samF2
    s0[c(1:9)][abs(s0[c(1:9)])<0.00000001]<-0.000001
    n0[c(1:8)][abs(n0[c(1:8)])<0.00000001]<-0.000001
    n_iter<-0;aaa1<-1000;AA<-matrix(0,6,1)
    while(aaa1>0.0001)
    {
      n_iter<-n_iter+1
      aa11<-aa
      mm11<-as.matrix(c(mean[1],mean[2],mean[3],mean1[1],mean1[2],mean1[3],mean1[4],mean2[1],mean2[2],mean2[3],mean2[4],
                        mean3[1],mean3[2],mean3[3],mean3[4],mean3[5],mean3[6],mean3[7],mean3[8],mean3[9]))
      B1<-solve(t(aa11)%*%aa11)%*%(t(aa11)%*%mm11)
      gs[1]<-B1[7];gs[2]<-B1[8];gs[3]<-B1[9];gs[4]<-B1[10];gs[5]<-B1[11];gs[6]<-B1[12];gs[7]<-B1[13];gs[8]<-B1[14];
      g_aa1<-(0.5*(gs[2]+gs[5])^2+0.25*(gs[4]+gs[6])^2)/n_fam
      #   0.5(db+i)**2+0.25(hb+jab)**2.
      g_aa2<-(0.5*(gs[1]+gs[5])^2+0.25*(gs[3]+gs[7])^2)/n_fam
      #   0.5(da+i)**2+0.25(ha+jba)**2.
      g_aa3<-(0.5*(gs[1]-gs[5])^2+0.25*(gs[3]-gs[7])^2)/n_fam
      #   0.5(da-i)**2+0.25(ha-jba)**2.
      g_aa4<-(0.5*(gs[2]-gs[5])^2+0.25*(gs[4]-gs[6])^2)/n_fam
      #   0.5(db-i)**2+0.25(hb-jab)**2.
      g_aa5<-0.25*(gs[1]^2+gs[2]^2+gs[5]^2+(gs[1]+gs[6])^2+(gs[2]+gs[7])^2+(gs[3]+gs[8]/2)^2+(gs[4]+gs[8]/2)^2+gs[8]^2/4)/n_fam
      sigma1[2]<-sigma1[1]+g_aa1;sigma1[3]<-sigma1[1]+g_aa2;sigma1[4]<-sigma1[1]+g_aa5
      sigma2[1]<-sigma2[4]+g_aa5;sigma2[2]<-sigma2[4]+g_aa3;sigma2[3]<-sigma2[4]+g_aa4
      sigma3[2]<-sigma3[1]+g_aa1;sigma3[3]<-sigma3[7]<-sigma3[9]<-sigma3[1]
      sigma3[4]<-sigma3[1]+g_aa2;sigma3[5]<-sigma3[1]+g_aa5;sigma3[6]<-sigma3[1]+g_aa3
      sigma3[8]<-sigma3[1]+g_aa4
      ##################################################################
      hh[1,1]<-sigma1[1]/n0[1]+sigma1[2]/n0[2]+sigma3[1]/s0[1]+sigma3[2]/s0[2]
      hh[1,2]<-0
      hh[1,3]<-sigma1[1]/n0[1]+sigma3[1]/s0[1]
      hh[1,4]<-0
      hh[1,5]<-0
      hh[1,6]<-0
      hh[2,2]<-sigma1[3]/n0[3]+sigma1[4]/n0[4]+sigma3[4]/s0[4]+sigma3[5]/s0[5]
      hh[2,3]<-sigma1[4]/n0[4]+sigma3[5]/s0[5]
      hh[2,4]<--sigma3[5]/s0[5]
      hh[2,5]<-0
      hh[2,6]<--sigma3[5]/s0[5]
      hh[3,3]<-sigma1[1]/n0[1]+sigma1[4]/n0[4]+sigma3[1]/s0[1]+sigma3[5]/s0[5]
      hh[3,4]<--sigma3[5]/s0[5]
      hh[3,5]<-0
      hh[3,6]<--sigma3[5]/s0[5]
      hh[4,4]<-sigma2[1]/n0[5]+sigma2[2]/n0[6]+sigma3[5]/s0[5]+sigma3[6]/s0[6]
      hh[4,5]<-0
      hh[4,6]<-sigma2[1]/n0[5]+sigma3[5]/s0[5]
      hh[5,5]<-sigma2[3]/n0[7]+sigma2[4]/n0[8]+sigma3[8]/s0[8]+sigma3[9]/s0[9]
      hh[5,6]<-sigma2[4]/n0[8]+sigma3[9]/s0[9]
      hh[6,6]<-sigma2[1]/n0[5]+sigma2[4]/n0[8]+sigma3[1]/s0[1]+sigma3[9]/s0[9]
      for(i in 2:6)
      {
        for(j in 1:(i-1))
        {
          hh[i,j]<-hh[j,i]
        }
      }
      ##################################################
      b_line2[1]<-sumwx1[1]/n0[1]-sumwx1[2]/n0[2]-sumwx3[1]/s0[1]+sumwx3[2]/s0[2]
      b_line2[2]<-sumwx1[3]/n0[3]-sumwx1[4]/n0[4]-sumwx3[4]/s0[4]+sumwx3[5]/s0[5]
      b_line2[3]<-sumwx1[1]/n0[1]-sumwx1[4]/n0[4]-sumwx3[1]/s0[1]+sumwx3[5]/s0[5]
      b_line2[4]<-sumwx2[1]/n0[5]-sumwx2[2]/n0[6]-sumwx3[5]/s0[5]+sumwx3[6]/s0[6]
      b_line2[5]<-sumwx2[3]/n0[7]-sumwx2[4]/n0[8]-sumwx3[8]/s0[8]+sumwx3[9]/s0[9]
      b_line2[6]<-sumwx2[1]/n0[5]-sumwx2[4]/n0[8]-sumwx3[5]/s0[5]+sumwx3[9]/s0[9]
      B2<-solve(hh,b_line2)
      mean[1]<-sumx[1]/n_samP1;mean[2]<-sumx[2]/n_samF1;mean[3]<-sumx[3]/n_samP2
      mean1[1]<-(sumwx1[1]-sigma1[1]*(B2[1]+B2[3]))/n0[1]
      mean1[2]<-(sumwx1[2]+sigma1[2]*B2[1])/n0[2]
      mean1[3]<-(sumwx1[3]-sigma1[3]*B2[2])/n0[3]
      mean1[4]<-(sumwx1[4]+sigma1[4]*(B2[2]+B2[3]))/n0[4]
      mean2[1]<-(sumwx2[1]-sigma2[1]*(B2[4]+B2[6]))/n0[5]
      mean2[2]<-(sumwx2[2]+sigma2[2]*B2[4])/n0[6]
      mean2[3]<-(sumwx2[3]-sigma2[3]*B2[5])/n0[7]
      mean2[4]<-(sumwx2[4]+sigma2[4]*(B2[5]+B2[6]))/n0[8]
      mean3[1]<-(sumwx3[1]+sigma3[1]*(B2[1]+B2[3]))/s0[1]
      mean3[2]<-(sumwx3[2]-sigma3[2]*B2[1])/s0[2]
      mean3[3]<-sumwx3[3]/s0[3]
      mean3[7]<-sumwx3[7]/s0[7]
      mean3[4]<-(sumwx3[4]+sigma3[4]*B2[2])/s0[4]
      mean3[5]<-(sumwx3[5]+sigma3[5]*(-B2[2]-B2[3]+B2[4]+B2[6]))/s0[5]
      mean3[6]<-(sumwx3[6]-sigma3[6]*B2[4])/s0[6]
      mean3[8]<-(sumwx3[8]+sigma3[8]*B2[5])/s0[8]
      mean3[9]<-(sumwx3[9]-sigma3[9]*(B2[5]+B2[6]))/s0[9]
      aaa1<-max(abs(B2-AA))
      AA<-B2
      if (n_iter>20) break
    }
    ############obtain variance######################
    ss1<-sum((dataP1-mean[1])^2);ss3<-sum((dataP2-mean[3])^2);ss2<-sum((dataF1-mean[2])^2)
    for(i in 1:4) {swx1[i] <- W1[i,]%*%(dataB1-mean1[i])^2 };for(i in 1:4) {swx2[i] <- W2[i,]%*%(dataB2-mean2[i])^2 };for(i in 1:9) {swx3[i] <- W3[i,]%*%(dataF2-mean3[i])^2 }
    aaa0<-sigma1[1];n_iter<-0;aaa1<-1000
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      aa2<-sigma1[1]/(sigma1[1]+g_aa1)
      aa3<-sigma1[1]/(sigma1[1]+g_aa2)
      aa4<-sigma1[1]/(sigma1[1]+g_aa5)
      as1<-swx1[1]+swx1[2]*aa2^2+swx1[3]*aa3^2+swx1[4]*aa4^2
      as2<-n0[1]+aa2*n0[2]+aa3*n0[3]+aa4*n0[4]
      sigma1[1]<-as1/as2
      aaa1<-abs(sigma1[1]-aaa0)
      aaa0<-sigma1[1]
      if (n_iter>20) break
    }
    sigma40<-sigma1[1]-sigma;
    if (sigma40<0) {sigma40<-0;sigma1[1]<-sigma}
    sigma1[1]<-sigma40+sigma;sigma1[2]<-sigma1[1]+g_aa1
    sigma1[3]<-sigma1[1]+g_aa2;sigma1[4]<-sigma1[1]+g_aa5
    # to estimate sigma50.
    aaa0<-sigma2[4];n_iter<-0;aaa1<-1000
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma2[4]/(sigma2[4]+g_aa5)
      aa2<-sigma2[4]/(sigma2[4]+g_aa3)
      aa3<-sigma2[4]/(sigma2[4]+g_aa4)
      as3<-swx2[1]*aa1^2+swx2[2]*aa2^2+swx2[3]*aa3^2+swx2[4]
      as4<-aa1*n0[5]+aa2*n0[6]+aa3*n0[7]+n0[8]
      sigma2[4]<-as3/as4
      aaa1<-abs(sigma2[4]-aaa0)
      aaa0<-sigma2[4]
      if (n_iter>20) break
    }
    sigma50<-sigma2[4]-sigma;
    if (sigma50<0) {sigma50<-0;sigma2[4]<-sigma}
    sigma2[4]<-sigma50+sigma;sigma2[1]<-sigma2[4]+g_aa5
    sigma2[2]<-sigma2[4]+g_aa3;sigma2[3]<-sigma2[4]+g_aa4
    # to estimate sigma60.
    aaa0<-sigma3[1];aa6<-swx3[1]+swx3[3]+swx3[7]+swx3[9]
    aa7<-s0[1]+s0[3]+s0[7]+s0[9]
    n_iter<-0;aaa1<-1000
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma3[1]/(sigma3[1]+g_aa1)
      aa2<-sigma3[1]/(sigma3[1]+g_aa2)
      aa3<-sigma3[1]/(sigma3[1]+g_aa3)
      aa4<-sigma3[1]/(sigma3[1]+g_aa4)
      aa5<-sigma3[1]/(sigma3[1]+g_aa5)
      as5<-aa6+swx3[2]*aa1^2+swx3[4]*aa2^2+swx3[5]*aa5^2+swx3[6]*aa3^2+swx3[8]*aa4^2
      as6<-aa7+aa1*s0[2]+aa2*s0[4]+aa5*s0[5]+aa3*s0[6]+aa4*s0[8]
      sigma3[1]<-as5/as6
      aaa1<-abs(sigma3[1]-aaa0)
      aaa0<-sigma3[1]
      if (n_iter>20) break
    }
    sigma60<-sigma3[1]-sigma;
    if (sigma60<0) {sigma60<-0;sigma3[1]<-sigma}
    sigma3[1]<-sigma60+sigma;sigma3[2]<-sigma3[1]+g_aa1
    sigma3[3]<-sigma3[7]<-sigma3[9]<-sigma3[1]
    sigma3[4]<-sigma3[1]+g_aa2;sigma3[5]<-sigma3[1]+g_aa5
    sigma3[6]<-sigma3[1]+g_aa3;sigma3[8]<-sigma3[1]+g_aa4
    # CM3 to estimate the variance (sigma).
    ab1<-ss1+ss2+ss3;ab2<-n_samP1+n_samF1+n_samP2
    n_iter<-0;aaa0<-sigma;aaa1<-1000
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      n0[11]<-sigma/(sigma+sigma40)
      n0[12]<-sigma/(sigma+sigma40+g_aa1)
      n0[13]<-sigma/(sigma+sigma40+g_aa2)
      n0[14]<-sigma/(sigma+sigma40+g_aa5)
      s0[11]<-sigma/(sigma+sigma50+g_aa5)
      s0[12]<-sigma/(sigma+sigma50+g_aa3)
      s0[13]<-sigma/(sigma+sigma50+g_aa4)
      s0[14]<-sigma/(sigma+sigma50)
      ab3<-sum(swx1[c(1:4)]*n0[c(11:14)]^2+swx2[c(1:4)]*s0[c(11:14)]^2)
      ab4<-sum(n0[c(1:4)]*n0[c(11:14)]+n0[c(5:8)]*s0[c(11:14)])
      n0[11]<-sigma/(sigma+sigma60);n0[13]<-n0[17]<-n0[19]<-n0[11]
      n0[12]<-sigma/(sigma+sigma60+g_aa1);n0[14]<-sigma/(sigma+sigma60+g_aa2)
      n0[15]<-sigma/(sigma+sigma60+g_aa5);n0[16]<-sigma/(sigma+sigma60+g_aa3)
      n0[18]<-sigma/(sigma+sigma60+g_aa4)
      ab3<-ab3+sum(swx3[c(1:9)]*n0[c(11:19)]^2)
      ab4<-ab4+sum(s0[c(1:9)]*n0[11:19])
      sigma<-(ab1+ab3)/(ab2+ab4)
      aaa1<-abs(sigma-aaa0)
      aaa0<-sigma
      if (n_iter>20) break
    }
    sigma1[1]<-sigma+sigma40;sigma1[2]<-sigma1[1]+g_aa1
    sigma1[3]<-sigma1[1]+g_aa2;sigma1[4]<-sigma1[1]+g_aa5
    sigma2[4]<-sigma+sigma50;sigma2[1]<-sigma2[4]+g_aa5
    sigma2[2]<-sigma2[4]+g_aa3;sigma2[3]<-sigma2[4]+g_aa4
    sigma3[1]<-sigma+sigma60
    sigma3[3]<-sigma3[7]<-sigma3[9]<-sigma3[1]
    sigma3[2]<-sigma3[1]+g_aa1;sigma3[4]<-sigma3[1]+g_aa2
    sigma3[5]<-sigma3[1]+g_aa5;sigma3[6]<-sigma3[1]+g_aa3
    sigma3[8]<-sigma3[1]+g_aa4
    if(sum(sigma < 1e-30)>=1){break}
    ##############################criteria for iterations to stop########################
    L1<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,mean[3],sqrt(sigma))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mix_pi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mix_pi3)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*18
  meanP1<-mean[1];meanF1<-mean[2];meanP2<-mean[3]
  sigma0<-sigma
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*n_samP1)
  P1bmw <- matrix(0,n_samP1,1)
  P1gg <- (dataP1 - mean[1])/sqrt(as.vector(sigma))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < n_samP1){P1bmw <- P1bmw+runif(n_samP1)/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  P1u<- as.matrix(c(12*n_samP1*((P1dd[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((P1dd[2]/n_samP1-1/3)^2),180*n_samP1*((P1dd[3]/n_samP1-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[2])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),P1D))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*n_samF1)
  F1bmw <- matrix(0,n_samF1,1)
  F1gg <- (dataF1 - mean[2])/sqrt(as.vector(sigma))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < n_samF1){F1bmw <- F1bmw+runif(n_samF1)/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  F1u<- as.matrix(c(12*n_samF1*((F1dd[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((F1dd[2]/n_samF1-1/3)^2),180*n_samF1*((F1dd[3]/n_samF1-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[2])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),F1D))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*n_samP2)
  P2bmw <- matrix(0,n_samP2,1)
  P2gg <- (dataP2 - mean[3])/sqrt(as.vector(sigma))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < n_samP2){P2bmw <- P2bmw+runif(n_samP2)/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  P2u<- as.matrix(c(12*n_samP2*((P2dd[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((P2dd[2]/n_samP2-1/3)^2),180*n_samP2*((P2dd[3]/n_samP2-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[2])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),P2D))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1);
  B1w1<-1/(12*n_samB1)
  B1bmw <- matrix(0,n_samB1,1); B1bmwsl <- matrix(0,n_samB1,4)
  for(i in 1:4){
    B1gg <- (dataB1 - mean1[i])/sqrt(sigma1[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi1[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B1P2)))[1]
  if(nn < n_samB1){B1P2 <- B1P2+runif(n_samB1)/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*n_samB1) + sum((B1P2 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  B1u <- as.matrix(c(12*n_samB1*((B1dd[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((B1dd[2]/n_samB1-1/3)^2),180*n_samB1*((B1dd[3]/n_samB1-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[2])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),B1D))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])

  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2);
  B2w1<-1/(12*n_samB2)
  B2bmw <- matrix(0,n_samB2,1); B2bmwsl <- matrix(0,n_samB2,4)
  for(i in 1:4){
    B2gg <- (dataB2 - mean2[i])/sqrt(sigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B2P2)))[1]
  if(nn < n_samB2){B2P2 <- B2P2+runif(n_samB2)/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*n_samB2) + sum((B2P2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  B2u <- as.matrix(c(12*n_samB2*((B2dd[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((B2dd[2]/n_samB2-1/3)^2),180*n_samB2*((B2dd[3]/n_samB2-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[2])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),B2D))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])

  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2);
  F2w1<-1/(12*n_samF2)
  F2bmw <- matrix(0,n_samF2,1); F2bmwsl <- matrix(0,n_samF2,9)
  for(i in 1:9){
    F2gg <- (dataF2 - mean3[i])/sqrt(sigma3[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi3[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < n_samF2){F2P2 <- F2P2+runif(n_samF2)/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*n_samF2) + sum((F2P2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  F2u <- as.matrix(c(12*n_samF2*((F2dd[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((F2dd[2]/n_samF2-1/3)^2),180*n_samF2*((F2dd[3]/n_samF2-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[2])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),F2D))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])

  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)

  ############first order parameters#############################
  aa3<-aa
  b_line3<-as.matrix(c(mean[1],mean[2],mean[3],mean1[1],mean1[2],mean1[3],mean1[4],mean2[1],mean2[2],mean2[3],mean2[4],
                       mean3[1],mean3[2],mean3[3],mean3[4],mean3[5],mean3[6],mean3[7],mean3[8],mean3[9]))
  B3<-solve(crossprod(aa3,aa3))%*%crossprod(aa3,b_line3)
  ########second order parameters########################
  jj1<-sigmaB1-sigma1[1]
  if (jj1<0 || jj1>=sigmaB1) {jj1<-0}
  ll1<-jj1/sigmaB1
  mm1<-sigma1[1]-sigma
  if (mm1<0 || mm1>=sigmaB1) {mm1<-0}
  nn1<-mm1/sigmaB1
  jj2<-sigmaB2-sigma2[4]
  if (jj2<0 || jj2>=sigmaB2) {jj2<-0}
  ll2<-jj2/sigmaB2
  mm2<-sigma2[4]-sigma
  if (mm2<0 || mm2>=sigmaB2) {mm2<-0}
  nn2<-mm2/sigmaB2
  jj3<-sigmaF2-sigma3[1]
  if (jj3<0 || jj3>=sigmaF2) {jj3<-0}
  ll3<-jj3/sigmaF2
  mm3<-sigma3[1]-sigma
  if (mm3<0 || mm3>=sigmaF2) {mm3<-0}
  nn3<-mm3/sigmaF2

  output <- data.frame("MX2-ADI-ADI",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4), round(sigma0,4),round(t(mean1),4),round(t(sigma1),4),
                       round(t(mix_pi1),4),round(t(mean2),4),round(t(sigma2),4),round(t(mix_pi2),4),
                       round(t(mean3),4),round(t(sigma3),4),round(t(mix_pi3),4),
                       round(B3[1],4),round(B3[2],4),round(B3[3],4),round(B3[4],4),round(B3[5],4),round(B3[6],4),round(B3[7],4),round(B3[8],4),round(B3[9],4),round(B3[10],4),round(B3[11],4),round(B3[12],4),round(B3[13],4),round(B3[14],4)," "," ",
                       round(jj1,4),round(ll1*100,4),round(mm1,4),round(nn1*100,4),round(jj2,4),round(ll2*100,4),round(mm2,4),round(nn2*100,4),round(jj3,4),round(ll3*100,4),round(mm3,4),round(nn3*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi1,mi2,mi3)
  return(OUTPUT)
}

######################MX2-ADI-AD(E-1)#####################################
G6FModelFun[[19]] <- function(K1,logL,df11,df21,df31,df41,df51,df61,G6Ftext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1];n_samF1<-dim(dataF1)[1];n_samP2<-dim(dataP2)[1]
  n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1];n_samF2<-dim(dataF2)[1]
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2),sum(dataF2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2),sum(dataF2^2)))
  ss<-matrix(0,3,1);ss[1]<-s[1]-sumx[1]^2/n_samP1;ss[2]<-s[2]-sumx[2]^2/n_samF1;ss[3]<-s[3]-sumx[3]^2/n_samP2
  mean<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataB1),mean(dataB2),mean(dataF2)))
  sigma0<-sum(ss)/(n_samP1+n_samF1+n_samP2-3)
  sigmaB1<-var(dataB1);sigma40<-sigmaB1
  sigmaB2<-var(dataB2);sigma50<-sigmaB2
  sigmaF2<-var(dataF2);sigma60<-sigmaF2
  sigmaP1<-sigmaF1<-sigmaP2<-sigma0
  m_esp <- 0.0001;n_fam <- as.numeric(G6Ftext2)
  ###############procedure start###########################
  mi1<-matrix(0.25,4,1);sigma1<-matrix(0,4,1)
  mi2<-matrix(0.25,4,1);sigma2<-matrix(0,4,1)
  mi3<-as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  sigma3<-matrix(0,9,1)
  sigma<-sigma0
  a1<-sqrt(sigma40/n_samB1);if (mean[1]<mean[3]) a1<--a1
  mean1<-as.matrix(c(mean[4]+2.5*a1,mean[4]+0.5*a1,mean[4]-0.5*a1,mean[4]-2.5*a1))
  a2<-sqrt(sigma50/n_samB2);if (mean[1]<mean[3]) a2<--a2
  mean2<-as.matrix(c(mean[5]+2.5*a2,mean[5]+0.5*a2,mean[5]-0.5*a2,mean[5]-2.5*a2))
  a3<-sqrt(sigma60/n_samF2);if (mean[1]<mean[3]) a3<--a3
  mean3<-as.matrix(c(mean[6]+3*a3,mean[6]+2.1*a3,mean[6]+1.4*a3,mean[6]+0.7*a3,mean[6],mean[6]-0.7*a3,mean[6]-1.4*a3,mean[6]-2.1*a3,mean[6]-3*a3))
  sigma1[1]<-sigmaB1/2;sigma2[4]<-sigmaB2/2;sigma3[1]<-sigmaF2/2
  gs<-matrix(0,8,1)
  gs[1]<-0.01844*mean[1]-0.0013*mean[2]-0.02236*mean[3]+0.08161*mean1[1]-
    0.02293*mean1[2]-0.03212*mean1[3]-0.06119*mean1[4]+0.06467*mean2[1]+
    0.03081*mean2[2]+0.03655*mean2[3]-0.08506*mean2[4]+0.14454*mean3[1]+
    0.03309*mean3[2]+0.24459*mean3[3]+0.03081*mean3[4]+0.00174*mean3[5]-
    0.03212*mean3[6]-0.25541*mean3[7]-0.02639*mean3[8]-0.14799*mean3[9]
  #da.
  gs[2]<-0.01844*mean[1]-0.0013*mean[2]-0.02236*mean[3]+0.08161*mean1[1]-
    0.02293*mean1[2]-0.03212*mean1[3]-0.06119*mean1[4]+0.06467*mean2[1]+
    0.03081*mean2[2]+0.03655*mean2[3]-0.08506*mean2[4]+0.14454*mean3[1]+
    0.03309*mean3[2]-0.25541*mean3[3]+0.03081*mean3[4]+0.00174*mean3[5]-
    0.03212*mean3[6]+0.24459*mean3[7]-0.02639*mean3[8]-0.14799*mean3[9]
  #db.
  gs[3]<--0.19992*mean[1]-0.12552*mean[2]-0.17662*mean[3]-0.13771*mean1[1]+
    0.0197*mean1[2]+0.43142*mean1[3]+0.1557*mean1[4]+0.179*mean2[1]+
    0.44307*mean2[2]-0.06186*mean2[3]-0.13771*mean2[4]-0.12606*mean3[1]-
    0.15507*mean3[2]-0.4637*mean3[3]+0.44307*mean3[4]+0.16735*mean3[5]+
    0.43142*mean3[6]-0.4637*mean3[7]-0.07351*mean3[8]-0.14936*mean3[9]
  #ha.
  gs[4]<--0.15292*mean[1]-0.08351*mean[2]-0.09759*mean[3]-0.11036*mean1[1]+
    0.53258*mean1[2]-0.05558*mean1[3]+0.08368*mean1[4]+0.139*mean2[1]-
    0.02792*mean2[2]+0.33893*mean2[3]-0.11036*mean2[4]-0.0827*mean3[1]+
    0.11763*mean3[2]-0.34598*mean3[3]-0.02792*mean3[4]+0.11134*mean3[5]-
    0.05558*mean3[6]-0.34598*mean3[7]+0.31127*mean3[8]-0.13802*mean3[9]
  #hb.
  gs[5]<-0.03146*mean[1]+0.02343*mean[2]+0.03883*mean[3]+0.10843*mean1[1]+
    0.03146*mean1[2]+0.00987*mean1[3]-0.03492*mean1[4]-0.02756*mean2[1]+
    0.01356*mean2[2]+0.00569*mean2[3]+0.10843*mean2[4]+0.11211*mean3[1]-
    0.02376*mean3[2]-0.24799*mean3[3]+0.01356*mean3[4]-0.03124*mean3[5]+
    0.00987*mean3[6]-0.24799*mean3[7]+0.00201*mean3[8]+0.10475*mean3[9]
  #i.
  gs[6]<--0.12844*mean[1]+0.02088*mean[2]+0.19107*mean[3]-0.13908*mean1[1]+
    0.36686*mean1[2]+0.0139*mean1[3]-0.02092*mean1[4]-0.03475*mean2[1]+
    0.00698*mean2[2]-0.58473*mean2[3]+0.19426*mean2[4]-0.14599*mean3[1]+
    0.47059*mean3[2]-0.41351*mean3[3]+0.00698*mean3[4]-0.02784*mean3[5]+
    0.0139*mean3[6]+0.58649*mean3[7]-0.57782*mean3[8]+0.20117*mean3[9]
  #jab.
  gs[7]<-(-mean[1]+mean[3]-mean1[1]+3*mean1[3]-3*mean2[2]+mean2[4]-
            mean3[1]+3*mean3[3]+3*mean3[4]-3*mean3[6]-3*mean3[7]+
            mean3[9])/6.0	 #jba.
  gs[8]<-1.01848*mean[1]+0.65486*mean[2]+0.9461*mean[3]-0.1674*mean1[1]-
    0.7552*mean1[2]-0.65447*mean1[3]+0.49638*mean1[4]+0.424*mean2[1]-
    0.69066*mean2[2]-0.50187*mean2[3]-0.1674*mean2[4]-0.20359*mean3[1]-
    0.21234*mean3[2]+0.64749*mean3[3]-0.69066*mean3[4]+0.46019*mean3[5]-
    0.65447*mean3[6]+0.64749*mean3[7]-0.46567*mean3[8]-0.13121*mean3[9]
  #l.
  g_aa1<-(0.5*(gs[2]+gs[5])^2+0.25*(gs[4]+gs[6])^2)/n_fam
  #   0.5(db+i)**2+0.25(hb+jab)**2.
  g_aa2<-(0.5*(gs[1]+gs[5])^2+0.25*(gs[3]+gs[7])^2)/n_fam
  #   0.5(da+i)**2+0.25(ha+jba)**2.
  g_aa3<-(0.5*(gs[1]-gs[5])^2+0.25*(gs[3]-gs[7])^2)/n_fam
  #   0.5(da-i)**2+0.25(ha-jba)**2.
  g_aa4<-(0.5*(gs[2]-gs[5])^2+0.25*(gs[4]-gs[6])^2)/n_fam
  #   0.5(db-i)**2+0.25(hb-jab)**2.
  g_aa5<-0.25*(gs[1]^2+gs[2]^2+gs[5]^2+(gs[1]+gs[6])^2+(gs[2]+gs[7])^2+(gs[3]+gs[8]/2)^2+(gs[4]+gs[8]/2)^2+gs[8]^2/4)/n_fam
  sigma1[2]<-sigma1[1]+g_aa1;sigma1[3]<-sigma1[1]+g_aa2;sigma1[4]<-sigma1[1]+g_aa5
  sigma2[1]<-sigma2[4]+g_aa5;sigma2[2]<-sigma2[4]+g_aa3;sigma2[3]<-sigma2[4]+g_aa4
  sigma3[2]<-sigma3[1]+g_aa1
  sigma3[3]<-sigma3[7]<-sigma3[9]<-sigma3[1]
  sigma3[4]<-sigma3[1]+g_aa2;sigma3[5]<-sigma3[1]+g_aa5;sigma3[6]<-sigma3[1]+g_aa3;sigma3[8]<-sigma3[1]+g_aa4
  L0<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma0))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma0))))+
    sum(log(dnorm(dataP2,mean[3],sqrt(sigma0))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mi1)))+
    sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mi3)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,4,n_samB1);      swx1 <- matrix(0,4,1)
  W2 <- matrix(0,4,n_samB2);      swx2 <- matrix(0,4,1)
  W3 <- matrix(0,9,n_samF2);      swx3 <- matrix(0,9,1)
  hh<-matrix(0,9,9);b_line<-matrix(0,9,1)
  n0<-matrix(0,18,1);s0<-matrix(0,18,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:4) { W1[i,] <- mi1[i]*dnorm(dataB1,mean1[i],sqrt(sigma1[i]))/dmixnorm(dataB1,mean1,sqrt(sigma1),mi1) }
    mix_pi1 <- as.matrix(rowSums(W1)/n_samB1)
    sumwx1 <- W1%*%dataB1
    for(i in 1:4) { W2[i,] <- mi2[i]*dnorm(dataB2,mean2[i],sqrt(sigma2[i]))/dmixnorm(dataB2,mean2,sqrt(sigma2),mi2) }
    mix_pi2 <- as.matrix(rowSums(W2)/n_samB2)
    sumwx2 <- W2%*%dataB2
    for(i in 1:9) { W3[i,] <- mi3[i]*dnorm(dataF2,mean3[i],sqrt(sigma3[i]))/dmixnorm(dataF2,mean3,sqrt(sigma3),mi3) }
    mix_pi3 <- as.matrix(rowSums(W3)/n_samF2)
    sumwx3 <- W3%*%dataF2
    # CM1-step for means.
    n0[c(1:4)]<-mix_pi1[c(1:4)]*n_samB1;n0[c(5:8)]<-mix_pi2[c(1:4)]*n_samB2
    s0[c(1:9)]<-mix_pi3[c(1:9)]*n_samF2
    s0[c(1:9)][abs(s0[c(1:9)])<0.00000001]<-0.000001
    n0[c(1:8)][abs(n0[c(1:8)])<0.00000001]<-0.000001
    aaa0<-0;AA<-matrix(0,9,1);n_iter<-0;aaa1<-1000
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      gs[1]<-0.01844*mean[1]-0.0013*mean[2]-0.02236*mean[3]+0.08161*mean1[1]-
        0.02293*mean1[2]-0.03212*mean1[3]-0.06119*mean1[4]+0.06467*mean2[1]+
        0.03081*mean2[2]+0.03655*mean2[3]-0.08506*mean2[4]+0.14454*mean3[1]+
        0.03309*mean3[2]+0.24459*mean3[3]+0.03081*mean3[4]+0.00174*mean3[5]-
        0.03212*mean3[6]-0.25541*mean3[7]-0.02639*mean3[8]-0.14799*mean3[9]
      #da.
      gs[2]<-0.01844*mean[1]-0.0013*mean[2]-0.02236*mean[3]+0.08161*mean1[1]-
        0.02293*mean1[2]-0.03212*mean1[3]-0.06119*mean1[4]+0.06467*mean2[1]+
        0.03081*mean2[2]+0.03655*mean2[3]-0.08506*mean2[4]+0.14454*mean3[1]+
        0.03309*mean3[2]-0.25541*mean3[3]+0.03081*mean3[4]+0.00174*mean3[5]-
        0.03212*mean3[6]+0.24459*mean3[7]-0.02639*mean3[8]-0.14799*mean3[9]
      #db.
      gs[3]<--0.19992*mean[1]-0.12552*mean[2]-0.17662*mean[3]-0.13771*mean1[1]+
        0.0197*mean1[2]+0.43142*mean1[3]+0.1557*mean1[4]+0.179*mean2[1]+
        0.44307*mean2[2]-0.06186*mean2[3]-0.13771*mean2[4]-0.12606*mean3[1]-
        0.15507*mean3[2]-0.4637*mean3[3]+0.44307*mean3[4]+0.16735*mean3[5]+
        0.43142*mean3[6]-0.4637*mean3[7]-0.07351*mean3[8]-0.14936*mean3[9]
      #ha.
      gs[4]<--0.15292*mean[1]-0.08351*mean[2]-0.09759*mean[3]-0.11036*mean1[1]+
        0.53258*mean1[2]-0.05558*mean1[3]+0.08368*mean1[4]+0.139*mean2[1]-
        0.02792*mean2[2]+0.33893*mean2[3]-0.11036*mean2[4]-0.0827*mean3[1]+
        0.11763*mean3[2]-0.34598*mean3[3]-0.02792*mean3[4]+0.11134*mean3[5]-
        0.05558*mean3[6]-0.34598*mean3[7]+0.31127*mean3[8]-0.13802*mean3[9]
      #hb.
      gs[5]<-0.03146*mean[1]+0.02343*mean[2]+0.03883*mean[3]+0.10843*mean1[1]+
        0.03146*mean1[2]+0.00987*mean1[3]-0.03492*mean1[4]-0.02756*mean2[1]+
        0.01356*mean2[2]+0.00569*mean2[3]+0.10843*mean2[4]+0.11211*mean3[1]-
        0.02376*mean3[2]-0.24799*mean3[3]+0.01356*mean3[4]-0.03124*mean3[5]+
        0.00987*mean3[6]-0.24799*mean3[7]+0.00201*mean3[8]+0.10475*mean3[9]
      #i.
      gs[6]<--0.12844*mean[1]+0.02088*mean[2]+0.19107*mean[3]-0.13908*mean1[1]+
        0.36686*mean1[2]+0.0139*mean1[3]-0.02092*mean1[4]-0.03475*mean2[1]+
        0.00698*mean2[2]-0.58473*mean2[3]+0.19426*mean2[4]-0.14599*mean3[1]+
        0.47059*mean3[2]-0.41351*mean3[3]+0.00698*mean3[4]-0.02784*mean3[5]+
        0.0139*mean3[6]+0.58649*mean3[7]-0.57782*mean3[8]+0.20117*mean3[9]
      #jab.
      gs[7]<-(-mean[1]+mean[3]-mean1[1]+3*mean1[3]-3*mean2[2]+mean2[4]-
                mean3[1]+3*mean3[3]+3*mean3[4]-3*mean3[6]-3*mean3[7]+
                mean3[9])/6	 #jba.
      gs[8]<-1.01848*mean[1]+0.65486*mean[2]+0.9461*mean[3]-0.1674*mean1[1]-
        0.7552*mean1[2]-0.65447*mean1[3]+0.49638*mean1[4]+0.424*mean2[1]-
        0.69066*mean2[2]-0.50187*mean2[3]-0.1674*mean2[4]-0.20359*mean3[1]-
        0.21234*mean3[2]+0.64749*mean3[3]-0.69066*mean3[4]+0.46019*mean3[5]-
        0.65447*mean3[6]+0.64749*mean3[7]-0.46567*mean3[8]-0.13121*mean3[9]
      #l.
      g_aa1<-(0.5*(gs[2]+gs[5])^2+0.25*(gs[4]+gs[6])^2)/n_fam
      #   0.5(db+i)**2+0.25(hb+jab)**2.
      g_aa2<-(0.5*(gs[1]+gs[5])^2+0.25*(gs[3]+gs[7])^2)/n_fam
      #   0.5(da+i)**2+0.25(ha+jba)**2.
      g_aa3<-(0.5*(gs[1]-gs[5])^2+0.25*(gs[3]-gs[7])^2)/n_fam
      #   0.5(da-i)**2+0.25(ha-jba)**2.
      g_aa4<-(0.5*(gs[2]-gs[5])^2+0.25*(gs[4]-gs[6])^2)/n_fam
      #   0.5(db-i)**2+0.25(hb-jab)**2.
      g_aa5<-0.25*(gs[1]^2+gs[2]^2+gs[5]^2+(gs[1]+gs[6])^2+(gs[2]+gs[7])^2+(gs[3]+gs[8]/2)^2+(gs[4]+gs[8]/2)^2+gs[8]^2/4)/n_fam
      hh[1,1]<-sigma1[1]/n0[1]+sigma1[2]/n0[2]+sigma3[1]/s0[1]+sigma3[2]/s0[2]
      hh[1,2]<-0
      hh[1,3]<-sigma1[1]/n0[1]+sigma3[1]/s0[1]
      hh[1,4]<-0
      hh[1,5]<-0
      hh[1,6]<-0
      hh[1,7]<--sigma1[1]/n0[1]
      hh[1,8]<-0
      hh[1,9]<-7*sigma3[1]/s0[1]+4*sigma3[2]/s0[2]
      hh[2,2]<-sigma1[3]/n0[3]+sigma1[4]/n0[4]+sigma3[4]/s0[4]+sigma3[5]/s0[5]
      hh[2,3]<-sigma1[4]/n0[4]+sigma3[5]/s0[5]
      hh[2,4]<--sigma3[5]/s0[5]
      hh[2,5]<-0
      hh[2,6]<--sigma3[5]/s0[5]
      hh[2,7]<-sigma1[4]/n0[4]
      hh[2,8]<--sigma1[4]/n0[4]-2*sigma3[5]/s0[5]
      hh[2,9]<--4*sigma3[4]/s0[4]-16*sigma3[5]/s0[5]
      hh[3,3]<-sigma1[1]/n0[1]+sigma1[4]/n0[4]+sigma3[1]/s0[1]+sigma3[5]/s0[5]
      hh[3,4]<--sigma3[5]/s0[5]
      hh[3,5]<-0
      hh[3,6]<--sigma3[5]/s0[5]
      hh[3,7]<--sigma1[1]/n0[1]+sigma1[4]/n0[4]
      hh[3,8]<--sigma1[4]/n0[4]-2*sigma3[5]/s0[5]
      hh[3,9]<-7*sigma3[1]/s0[1]-16*sigma3[5]/s0[5]
      hh[4,4]<-sigma2[1]/n0[5]+sigma2[2]/n0[6]+sigma3[5]/s0[5]+sigma3[6]/s0[6]
      hh[4,5]<-0
      hh[4,6]<-sigma2[1]/n0[5]+sigma3[5]/s0[5]
      hh[4,7]<-sigma2[1]/n0[5]
      hh[4,8]<-sigma2[1]/n0[5]+2*sigma3[5]/s0[5]
      hh[4,9]<-16*sigma3[5]/s0[5]+4*sigma3[6]/s0[6]
      hh[5,5]<-sigma2[3]/n0[7]+sigma2[4]/n0[8]+sigma3[8]/s0[8]+sigma3[9]/s0[9]
      hh[5,6]<-sigma2[4]/n0[8]+sigma3[9]/s0[9]
      hh[5,7]<--sigma2[4]/n0[8]
      hh[5,8]<-0
      hh[5,9]<--4*sigma3[8]/s0[8]-7*sigma3[9]/s0[9]
      hh[6,6]<-sigma2[1]/n0[5]+sigma2[4]/n0[8]+sigma3[5]/s0[5]+sigma3[9]/s0[9]
      hh[6,7]<-sigma2[1]/n0[5]-sigma2[4]/n0[8]
      hh[6,8]<-sigma2[1]/n0[5]+2*sigma3[5]/s0[5]
      hh[6,9]<-16*sigma3[5]/s0[5]-7*sigma3[9]/s0[9]
      hh[7,7]<-sigma*(1/n_samP1+1/n_samP2)+sigma1[1]/n0[1]+sigma1[4]/n0[4]+sigma2[1]/n0[5]+sigma2[4]/n0[8]
      hh[7,8]<--sigma1[4]/n0[4]+sigma2[1]/n0[5]
      hh[7,9]<-6*(sigma/n_samP1-sigma/n_samP2)
      hh[8,8]<-sigma1[4]/n0[4]+sigma2[1]/n0[5]+4*sigma3[5]/s0[5]
      hh[8,9]<-32*sigma3[5]/s0[5]
      hh[9,9]<-sigma*(36/n_samP1+16/n_samF1+36/n_samP2)+49*sigma3[1]/s0[1]+16*
        sigma3[2]/s0[2]+sigma3[3]/s0[3]+16*sigma3[4]/s0[4]+256*sigma3[5]/s0[5]+16*sigma3[6]/s0[6]+sigma3[7]/s0[7]+16*sigma3[8]/s0[8]+49*sigma3[9]/s0[9]
      for(i in 2:9)
      {
        for(j in 1:(i-1))
        {
          hh[i,j]<-hh[j,i]
        }
      }
      ##########################################################################
      b_line[1]<-sumwx1[1]/n0[1]-sumwx1[2]/n0[2]-sumwx3[1]/s0[1]+sumwx3[2]/s0[2]
      b_line[2]<-sumwx1[3]/n0[3]-sumwx1[4]/n0[4]-sumwx3[4]/s0[4]+sumwx3[5]/s0[5]
      b_line[3]<-sumwx1[1]/n0[1]-sumwx1[4]/n0[4]-sumwx3[1]/s0[1]+sumwx3[5]/s0[5]
      b_line[4]<-sumwx2[1]/n0[5]-sumwx2[2]/n0[6]-sumwx3[5]/s0[5]+sumwx3[6]/s0[6]
      b_line[5]<-sumwx2[3]/n0[7]-sumwx2[4]/n0[8]-sumwx3[8]/s0[8]+sumwx3[9]/s0[9]
      b_line[6]<-sumwx2[1]/n0[5]-sumwx2[4]/n0[8]-sumwx3[5]/s0[5]+sumwx3[9]/s0[9]
      b_line[7]<-sumx[1]/n_samP1-sumx[3]/n_samP2-sumwx1[1]/n0[1]-sumwx1[4]/n0[4]+sumwx2[1]/n0[5]+sumwx2[4]/n0[8]
      b_line[8]<-sumwx1[4]/n0[4]+sumwx2[1]/n0[5]-2*sumwx3[5]/s0[5]
      b_line[9]<-6*sumx[1]/n_samP1+4*sumx[2]/n_samF1+6*sumx[3]/n_samP2-7*sumwx3[1]/s0[1]+4*sumwx3[2]/s0[2]-sumwx3[3]/s0[3]+4*sumwx3[4]/s0[4]-16*
        sumwx3[5]/s0[5]+4*sumwx3[6]/s0[6]-sumwx3[7]/s0[7]+4*sumwx3[8]/s0[8]-7*sumwx3[9]/s0[9]
      B<-solve(hh,b_line)
      ####################################################
      mean[1]<-(sumx[1]-sigma*(B[7]+6*B[9]))/n_samP1
      mean[2]<-(sumx[2]-sigma*4*B[9])/n_samF1
      mean[3]<-(sumx[3]+sigma*(B[7]-6*B[9]))/n_samP2
      mean1[1]<-(sumwx1[1]+sigma1[1]*(-B[1]-B[3]+B[7]))/n0[1]
      mean1[2]<-(sumwx1[2]+sigma1[2]*B[1])/n0[2]
      mean1[3]<-(sumwx1[3]-sigma1[3]*B[2])/n0[3]
      mean1[4]<-(sumwx1[4]+sigma1[4]*(B[2]+B[3]+B[7]-B[8]))/n0[4]
      mean2[1]<-(sumwx2[1]-sigma2[1]*(B[4]+B[6]+B[7]+B[8]))/n0[5]
      mean2[2]<-(sumwx2[2]+sigma2[2]*B[4])/n0[6]
      mean2[3]<-(sumwx2[3]-sigma2[3]*B[5])/n0[7]
      mean2[4]<-(sumwx2[4]+sigma2[4]*(B[5]+B[6]-B[7]))/n0[8]
      mean3[1]<-(sumwx3[1]+sigma3[1]*(B[1]+B[3]+7*B[9]))/s0[1]
      mean3[2]<-(sumwx3[2]-sigma3[2]*(B[1]+4*B[9]))/s0[2]
      mean3[3]<-(sumwx3[3]+sigma3[3]*B[9])/s0[3]
      mean3[7]<-(sumwx3[7]+sigma3[7]*B[9])/s0[7]
      mean3[4]<-(sumwx3[4]+sigma3[4]*(B[2]-4*B[9]))/s0[4]
      mean3[5]<-(sumwx3[5]+sigma3[5]*(-B[2]-B[3]+B[4]+B[6]+2*B[8]+16*B[9]))/s0[5]
      mean3[6]<-(sumwx3[6]-sigma3[6]*(B[4]+4*B[9]))/s0[6]
      mean3[8]<-(sumwx3[8]+sigma3[8]*(B[5]-4*B[9]))/s0[8]
      mean3[9]<-(sumwx3[9]+sigma3[9]*(-B[5]-B[6]+7*B[9]))/s0[9]
      aaa1<-max(abs(B-AA))
      AA<-B
      if (n_iter>20) break
    }
    ############obtain variance###########################
    ss1<-sum((dataP1-mean[1])^2);ss3<-sum((dataP2-mean[3])^2);ss2<-sum((dataF1-mean[2])^2)
    for(i in 1:4) {swx1[i] <- W1[i,]%*%(dataB1-mean1[i])^2 };for(i in 1:4) {swx2[i] <- W2[i,]%*%(dataB2-mean2[i])^2 };for(i in 1:9) {swx3[i] <- W3[i,]%*%(dataF2-mean3[i])^2 }
    # to estimate sigma40.
    aaa0<-sigma1[1];n_iter<-0;aaa1<-1000
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      aa2<-sigma1[1]/(sigma1[1]+g_aa1)
      aa3<-sigma1[1]/(sigma1[1]+g_aa2)
      aa4<-sigma1[1]/(sigma1[1]+g_aa5)
      as1<-swx1[1]+swx1[2]*aa2^2+swx1[3]*aa3^2+swx1[4]*aa4^2
      as2<-n0[1]+aa2*n0[2]+aa3*n0[3]+aa4*n0[4]
      sigma1[1]<-as1/as2
      aaa1<-abs(sigma1[1]-aaa0)
      aaa0<-sigma1[1]
      if (n_iter>20) break
    }
    sigma40<-sigma1[1]-sigma;
    if (sigma40<0) {sigma40<-0;sigma1[1]<-sigma}
    sigma1[1]<-sigma40+sigma;sigma1[2]<-sigma1[1]+g_aa1
    sigma1[3]<-sigma1[1]+g_aa2;sigma1[4]<-sigma1[1]+g_aa5
    # to estimate sigma50.
    aaa0<-sigma2[4];n_iter<-0;aaa1<-1000
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma2[4]/(sigma2[4]+g_aa5)
      aa2<-sigma2[4]/(sigma2[4]+g_aa3)
      aa3<-sigma2[4]/(sigma2[4]+g_aa4)
      as3<-swx2[1]*aa1^2+swx2[2]*aa2^2+swx2[3]*aa3^2+swx2[4]
      as4<-aa1*n0[5]+aa2*n0[6]+aa3*n0[7]+n0[8]
      sigma2[4]<-as3/as4
      aaa1<-abs(sigma2[4]-aaa0)
      aaa0<-sigma2[4]
      if (n_iter>20) break
    }
    sigma50<-sigma2[4]-sigma;
    if (sigma50<0) {sigma50<-0;sigma2[4]<-sigma}
    sigma2[4]<-sigma50+sigma;sigma2[1]<-sigma2[4]+g_aa5
    sigma2[2]<-sigma2[4]+g_aa3;sigma2[3]<-sigma2[4]+g_aa4
    # to estimate sigma60.
    aaa0<-sigma3[1];aa6<-swx3[1]+swx3[3]+swx3[7]+swx3[9]
    aa7<-s0[1]+s0[3]+s0[7]+s0[9]
    n_iter<-0;aaa1<-1000
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma3[1]/(sigma3[1]+g_aa1)
      aa2<-sigma3[1]/(sigma3[1]+g_aa2)
      aa3<-sigma3[1]/(sigma3[1]+g_aa3)
      aa4<-sigma3[1]/(sigma3[1]+g_aa4)
      aa5<-sigma3[1]/(sigma3[1]+g_aa5)
      as5<-aa6+swx3[2]*aa1^2+swx3[4]*aa2^2+swx3[5]*aa5^2+swx3[6]*aa3^2+swx3[8]*aa4^2
      as6<-aa7+aa1*s0[2]+aa2*s0[4]+aa5*s0[5]+aa3*s0[6]+aa4*s0[8]
      sigma3[1]<-as5/as6
      aaa1<-abs(sigma3[1]-aaa0)
      aaa0<-sigma3[1]
      if (n_iter>20) break
    }
    sigma60<-sigma3[1]-sigma;
    if (sigma60<0) {sigma60<-0;sigma3[1]<-sigma}
    sigma3[1]<-sigma60+sigma;sigma3[2]<-sigma3[1]+g_aa1
    sigma3[3]<-sigma3[7]<-sigma3[9]<-sigma3[1]
    sigma3[4]<-sigma3[1]+g_aa2;sigma3[5]<-sigma3[1]+g_aa5
    sigma3[6]<-sigma3[1]+g_aa3;sigma3[8]<-sigma3[1]+g_aa4
    # CM3 to estimate the variance (sigma).
    ab1<-ss1+ss2+ss3;ab2<-n_samP1+n_samF1+n_samP2
    n_iter<-0;aaa0<-sigma;aaa1<-1000
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      n0[11]<-sigma/(sigma+sigma40)
      n0[12]<-sigma/(sigma+sigma40+g_aa1)
      n0[13]<-sigma/(sigma+sigma40+g_aa2)
      n0[14]<-sigma/(sigma+sigma40+g_aa5)
      s0[11]<-sigma/(sigma+sigma50+g_aa5)
      s0[12]<-sigma/(sigma+sigma50+g_aa3)
      s0[13]<-sigma/(sigma+sigma50+g_aa4)
      s0[14]<-sigma/(sigma+sigma50)
      ab3<-sum(swx1[c(1:4)]*n0[c(11:14)]^2+swx2[c(1:4)]*s0[c(11:14)]^2)
      ab4<-sum(n0[c(1:4)]*n0[c(11:14)]+n0[c(5:8)]*s0[c(11:14)])
      n0[11]<-sigma/(sigma+sigma60);n0[13]<-n0[17]<-n0[19]<-n0[11]
      n0[12]<-sigma/(sigma+sigma60+g_aa1);n0[14]<-sigma/(sigma+sigma60+g_aa2)
      n0[15]<-sigma/(sigma+sigma60+g_aa5);n0[16]<-sigma/(sigma+sigma60+g_aa3)
      n0[18]<-sigma/(sigma+sigma60+g_aa4)
      ab3<-ab3+sum(swx3[c(1:9)]*n0[c(11:19)]^2)
      ab4<-ab4+sum(s0[c(1:9)]*n0[11:19])
      sigma<-(ab1+ab3)/(ab2+ab4)
      aaa1<-abs(sigma-aaa0)
      aaa0<-sigma
      if (n_iter>20) break
    }
    sigma1[1]<-sigma+sigma40;sigma1[2]<-sigma1[1]+g_aa1
    sigma1[3]<-sigma1[1]+g_aa2;sigma1[4]<-sigma1[1]+g_aa5
    sigma2[4]<-sigma+sigma50;sigma2[1]<-sigma2[4]+g_aa5
    sigma2[2]<-sigma2[4]+g_aa3;sigma2[3]<-sigma2[4]+g_aa4
    sigma3[1]<-sigma+sigma60;sigma3[3]<-sigma3[7]<-sigma3[9]<-sigma3[1]
    sigma3[2]<-sigma3[1]+g_aa1;sigma3[4]<-sigma3[1]+g_aa2
    sigma3[5]<-sigma3[1]+g_aa5;sigma3[6]<-sigma3[1]+g_aa3
    sigma3[8]<-sigma3[1]+g_aa4
    if(sum(sigma < 1e-30)>=1){break}
    ##############################criteria for iterations to stop########################
    L1<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,mean[3],sqrt(sigma))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mix_pi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mix_pi3)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*15
  meanP1<-mean[1];meanF1<-mean[2];meanP2<-mean[3]
  sigma0<-sigma

  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*n_samP1)
  P1bmw <- matrix(0,n_samP1,1)
  P1gg <- (dataP1 - mean[1])/sqrt(as.vector(sigma))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < n_samP1){P1bmw <- P1bmw+runif(n_samP1)/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  P1u<- as.matrix(c(12*n_samP1*((P1dd[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((P1dd[2]/n_samP1-1/3)^2),180*n_samP1*((P1dd[3]/n_samP1-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[2])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),P1D))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*n_samF1)
  F1bmw <- matrix(0,n_samF1,1)
  F1gg <- (dataF1 - mean[2])/sqrt(as.vector(sigma))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < n_samF1){F1bmw <- F1bmw+runif(n_samF1)/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  F1u<- as.matrix(c(12*n_samF1*((F1dd[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((F1dd[2]/n_samF1-1/3)^2),180*n_samF1*((F1dd[3]/n_samF1-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[2])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),F1D))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*n_samP2)
  P2bmw <- matrix(0,n_samP2,1)
  P2gg <- (dataP2 - mean[3])/sqrt(as.vector(sigma))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < n_samP2){P2bmw <- P2bmw+runif(n_samP2)/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  P2u<- as.matrix(c(12*n_samP2*((P2dd[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((P2dd[2]/n_samP2-1/3)^2),180*n_samP2*((P2dd[3]/n_samP2-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[2])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),P2D))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1);
  B1w1<-1/(12*n_samB1)
  B1bmw <- matrix(0,n_samB1,1); B1bmwsl <- matrix(0,n_samB1,4)
  for(i in 1:4){
    B1gg <- (dataB1 - mean1[i])/sqrt(sigma1[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi1[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B1P2)))[1]
  if(nn < n_samB1){B1P2 <- B1P2+runif(n_samB1)/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*n_samB1) + sum((B1P2 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  B1u <- as.matrix(c(12*n_samB1*((B1dd[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((B1dd[2]/n_samB1-1/3)^2),180*n_samB1*((B1dd[3]/n_samB1-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[2])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),B1D))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])

  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2);
  B2w1<-1/(12*n_samB2)
  B2bmw <- matrix(0,n_samB2,1); B2bmwsl <- matrix(0,n_samB2,4)
  for(i in 1:4){
    B2gg <- (dataB2 - mean2[i])/sqrt(sigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B2P2)))[1]
  if(nn < n_samB2){B2P2 <- B2P2+runif(n_samB2)/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*n_samB2) + sum((B2P2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  B2u <- as.matrix(c(12*n_samB2*((B2dd[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((B2dd[2]/n_samB2-1/3)^2),180*n_samB2*((B2dd[3]/n_samB2-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[2])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),B2D))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])

  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2);
  F2w1<-1/(12*n_samF2)
  F2bmw <- matrix(0,n_samF2,1); F2bmwsl <- matrix(0,n_samF2,9)
  for(i in 1:9){
    F2gg <- (dataF2 - mean3[i])/sqrt(sigma3[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi3[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < n_samF2){F2P2 <- F2P2+runif(n_samF2)/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*n_samF2) + sum((F2P2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  F2u <- as.matrix(c(12*n_samF2*((F2dd[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((F2dd[2]/n_samF2-1/3)^2),180*n_samF2*((F2dd[3]/n_samF2-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[2])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),F2D))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])

  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)


  ###################first order parameters###########################
  aa<-matrix(c(1,1,1,0,0,1,0,0,0,1,0,1,0,0,1,1,0,0,0,1,0,1,1,-1,-1,0,0,1,0,0,
               0,-1,0,1,1,1,0,0,1,0,0,0,0.5,0.25,1,1,0,0,0.5,0,0.5,0,0,0.5,0.25,
               1,0,1,0.5,0,0,0,0.5,0,0.5,0.25,1,0,0,0.5,0.5,0,0,0,0.25,0.5,0.25,
               1,0,0,0.5,0.5,0,0,0,0.25,-0.5,0.25,1,0,-1,0.5,0,0,0,-0.5,0,-0.5,0.25,
               1,-1,0,0,0.5,0,-0.5,0,0,-0.5,0.25,1,-1,-1,0,0,1,0,0,0,-0.5,0.25,1,1,1,
               0,0,1,0,0,0,0,0.25,1,1,0,0,0.5,0,0.5,0,0,0,0.25,1,1,-1,0,0,-1,0,0,0,0,
               0.25,1,0,1,0.5,0,0,0,0.5,0,0,0.25,1,0,0,0.5,0.5,0,0,0,0.25,0,0.25,1,0,-1,
               0.5,0,0,0,-0.5,0,0,0.25,1,-1,1,0,0,-1,0,0,0,0,0.25,1,-1,0,0,0.5,0,-0.5,0,
               0,0,0.25,1,-1,-1,0,0,1,0,0,0,0,0.25),20,11,byrow=T)
  b_line1<-as.matrix(c(mean[1],mean[2],mean[3],mean1[1],mean1[2],mean1[3],mean1[4],mean2[1],mean2[2],mean2[3],mean2[4],
                       mean3[1],mean3[2],mean3[3],mean3[4],mean3[5],mean3[6],mean3[7],mean3[8],mean3[9]))

  B1<-solve(crossprod(aa,aa))%*%crossprod(aa,b_line1)
  #########second order parameters###################
  jj1<-sigmaB1-sigma1[1]
  if (jj1<0 || jj1>=sigmaB1) {jj1<-0}
  ll1<-jj1/sigmaB1
  mm1<-sigma1[1]-sigma
  if (mm1<0 || mm1>=sigmaB1) {mm1<-0}
  nn1<-mm1/sigmaB1
  jj2<-sigmaB2-sigma2[4]
  if (jj2<0 || jj2>=sigmaB2) {jj2<-0}
  ll2<-jj2/sigmaB2
  mm2<-sigma2[4]-sigma
  if (mm2<0 || mm2>=sigmaB2) {mm2<-0}
  nn2<-mm2/sigmaB2
  jj3<-sigmaF2-sigma3[1]
  if (jj3<0 || jj3>=sigmaF2) {jj3<-0}
  ll3<-jj3/sigmaF2
  mm3<-sigma3[1]-sigma
  if (mm3<0 || mm3>=sigmaF2) {mm3<-0}
  nn3<-mm3/sigmaF2

  output <- data.frame("MX2-ADI-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4), round(sigma0,4),round(t(mean1),4),round(t(sigma1),4),
                       round(t(mix_pi1),4),round(t(mean2),4),round(t(sigma2),4),round(t(mix_pi2),4),
                       round(t(mean3),4),round(t(sigma3),4),round(t(mix_pi3),4),
                       round(B1[1],4)," "," "," "," "," ",round(B1[2],4),round(B1[3],4),round(B1[4],4),round(B1[5],4),round(B1[6],4),round(B1[7],4),round(B1[8],4),round(B1[9],4),round(B1[10],4),round(B1[11],4),
                       round(jj1,4),round(ll1*100,4),round(mm1,4),round(nn1*100,4),round(jj2,4),round(ll2*100,4),round(mm2,4),round(nn2*100,4),round(jj3,4),round(ll3*100,4),round(mm3,4),round(nn3*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi1,mi2,mi3)
  return(OUTPUT)
}


#################MX2-AD-AD(E-2)#################################
G6FModelFun[[20]] <- function(K1,logL,df11,df21,df31,df41,df51,df61,G6Ftext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1];n_samF1<-dim(dataF1)[1];n_samP2<-dim(dataP2)[1]
  n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1];n_samF2<-dim(dataF2)[1]
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2),sum(dataF2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2),sum(dataF2^2)))
  ss<-matrix(0,3,1);ss[1]<-s[1]-sumx[1]^2/n_samP1;ss[2]<-s[2]-sumx[2]^2/n_samF1;ss[3]<-s[3]-sumx[3]^2/n_samP2
  mean<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataB1),mean(dataB2),mean(dataF2)))
  sigma0<-sum(ss)/(n_samP1+n_samF1+n_samP2-3)
  sigmaB1<-var(dataB1);sigma40<-sigmaB1
  sigmaB2<-var(dataB2);sigma50<-sigmaB2
  sigmaF2<-var(dataF2);sigma60<-sigmaF2
  sigmaP1<-sigmaF1<-sigmaP2<-sigma0
  m_esp <- 0.0001;n_fam <- as.numeric(G6Ftext2)
  ###############procedure start###########################
  mi1<-matrix(0.25,4,1);sigma1<-matrix(0,4,1)
  mi2<-matrix(0.25,4,1);sigma2<-matrix(0,4,1)
  mi3<-as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  sigma3<-matrix(0,9,1)
  sigma<-sigma0
  a1<-sqrt(sigmaB1/n_samB1);if (mean[1]<mean[3]) a1<--a1
  mean1<-as.matrix(c(mean[1],mean[4]-0.5*a1,mean[4]-a1,mean[4]-1.5*a1))
  a2<-sqrt(sigmaB2/n_samB2);if (mean[1]<mean[3]) a2<--a2
  mean2<-as.matrix(c(mean[5]+2*a2,mean[5],mean[5]-a2,mean[5]-2*a2))
  a3<-sqrt(sigmaF2/n_samF2);if (mean[1]<mean[3]) a3<--a3
  mean3<-matrix(0,9,1)
  mean3[1]<-mean[1];mean3[2]<-mean[1]-0.5*a3;mean3[3]<-(mean[1]+mean[3])/2
  mean3[4]<-mean[6]+0.6*a3;mean3[5]<-mean[6];mean3[6]<-mean3[4]
  mean3[7]<-mean3[3];mean3[8]<-mean[6]-a3;mean3[9]<-mean[3]
  sigma1[1]<-sigmaB1/2;sigma2[4]<-sigmaB2/2;sigma3[1]<-sigmaF2/2
  gs<-matrix(0,4,1)
  gs[1]<-(mean1[1]+mean1[2]-mean1[3]-mean1[4]+mean2[1]+mean2[2]-
            mean2[3]-mean2[4]+2*mean3[1]+2*mean3[2]+2*mean3[3]-2*mean3[7]
          -2*mean3[8]-2*mean3[9])/16
  #da.
  gs[2]<-(mean1[1]-mean1[2]+mean1[3]-mean1[4]+mean2[1]-mean2[2]+
            mean2[3]-mean2[4]+2*mean3[1]-2*mean3[3]+2*mean3[4]-2*
            mean3[6]+2*mean3[7]-2*mean3[9])/16
  #db.
  gs[3]<-0.03846*mean[1]+0.02564*mean[2]+0.03846*mean[3]-0.19872*
    mean1[1]-0.21474*mean1[2]+0.28526*mean1[3]+0.26923*mean1[4]+
    0.26923*mean2[1]+0.28526*mean2[2]-0.21474*mean2[3]-0.19872*
    mean2[4]-0.19872*mean3[1]-0.21474*mean3[2]-0.19872*mean3[3]+
    0.28526*mean3[4]+0.26923*mean3[5]+0.28526*mean3[6]-0.19872*
    mean3[7]-0.21474*mean3[8]-0.19872*mean3[9]
  #ha.
  gs[4]<-0.03846*mean[1]+0.02564*mean[2]+0.03846*mean[3]-0.19872*
    mean1[1]+0.28526*mean1[2]-0.21474*mean1[3]+0.26923*mean1[4]+
    0.26923*mean2[1]-0.21474*mean2[2]+0.28526*mean2[3]-0.19872*
    mean2[4]-0.19872*mean3[1]+0.28526*mean3[2]-0.19872*mean3[3]-
    0.21474*mean3[4]+0.26923*mean3[5]-0.21474*mean3[6]-0.19872*
    mean3[7]+0.28526*mean3[8]-0.19872*mean3[9]
  #hb.
  g_aa1<-(0.5*gs[2]^2+0.25*gs[4]^2)/n_fam #   0.5db**2+0.25hb**2.
  g_aa2<-(0.5*gs[1]^2+0.25*gs[3]^2)/n_fam #   0.5da**2+0.25ha**2.
  g_aa3<-g_aa1+g_aa2
  sigma1[2]<-sigma1[1]+g_aa1;sigma1[3]<-sigma1[1]+g_aa2
  sigma1[4]<-sigma1[1]+g_aa3;sigma2[1]<-sigma2[4]+g_aa3
  sigma2[2]<-sigma2[4]+g_aa2;sigma2[3]<-sigma2[4]+g_aa1
  sigma3[2]<-sigma3[1]+g_aa1;sigma3[3]<-sigma3[7]<-sigma3[9]<-sigma3[1]
  sigma3[4]<-sigma3[1]+g_aa2;sigma3[5]<-sigma3[1]+g_aa3
  sigma3[6]<-sigma3[1]+g_aa2;sigma3[8]<-sigma3[1]+g_aa1
  L0<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma0))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma0))))+
    sum(log(dnorm(dataP2,mean[3],sqrt(sigma0))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mi1)))+
    sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mi3)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,4,n_samB1);      swx1 <- matrix(0,4,1)
  W2 <- matrix(0,4,n_samB2);      swx2 <- matrix(0,4,1)
  W3 <- matrix(0,9,n_samF2);      swx3 <- matrix(0,9,1)
  hh<-matrix(0,13,13);b_line<-matrix(0,13,1)
  n0<-matrix(0,18,1);s0<-matrix(0,18,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:4) { W1[i,] <- mi1[i]*dnorm(dataB1,mean1[i],sqrt(sigma1[i]))/dmixnorm(dataB1,mean1,sqrt(sigma1),mi1) }
    mix_pi1 <- as.matrix(rowSums(W1)/n_samB1)
    sumwx1 <- W1%*%dataB1
    for(i in 1:4) { W2[i,] <- mi2[i]*dnorm(dataB2,mean2[i],sqrt(sigma2[i]))/dmixnorm(dataB2,mean2,sqrt(sigma2),mi2) }
    mix_pi2 <- as.matrix(rowSums(W2)/n_samB2)
    sumwx2 <- W2%*%dataB2
    for(i in 1:9) { W3[i,] <- mi3[i]*dnorm(dataF2,mean3[i],sqrt(sigma3[i]))/dmixnorm(dataF2,mean3,sqrt(sigma3),mi3) }
    mix_pi3 <- as.matrix(rowSums(W3)/n_samF2)
    sumwx3 <- W3%*%dataF2
    ########obtain means#############################
    n0[c(1:4)]<-mix_pi1[c(1:4)]*n_samB1;n0[c(5:8)]<-mix_pi2[c(1:4)]*n_samB2
    s0[c(1:9)]<-mix_pi3[c(1:9)]*n_samF2
    s0[c(1:9)][abs(s0[c(1:9)])<0.00000001]<-0.000001
    n0[c(1:8)][abs(n0[c(1:8)])<0.00000001]<-0.000001
    aaa0<-0;AA<-matrix(0,13,1);n_iter<-0;aaa1<-1000
    while (aaa1>0.001)
    {
      n_iter<-n_iter+1
      gs[1]<-(mean1[1]+mean1[2]-mean1[3]-mean1[4]+mean2[1]+mean2[2]-
                mean2[3]-mean2[4]+2*mean3[1]+2*mean3[2]+2*mean3[3]-2*mean3[7]-
                2*mean3[8]-2*mean3[9])/16
      #da.
      gs[2]<-(mean1[1]-mean1[2]+mean1[3]-mean1[4]+mean2[1]-mean2[2]+
                mean2[3]-mean2[4]+2*mean3[1]-2*mean3[3]+2*mean3[4]-2*
                mean3[6]+2*mean3[7]-2*mean3[9])/16
      #db.
      gs[3]<-0.03846*mean[1]+0.02564*mean[2]+0.03846*mean[3]-0.19872*
        mean1[1]-0.21474*mean1[2]+0.28526*mean1[3]+0.26923*mean1[4]+
        0.26923*mean2[1]+0.28526*mean2[2]-0.21474*mean2[3]-0.19872*
        mean2[4]-0.19872*mean3[1]-0.21474*mean3[2]-0.19872*mean3[3]+
        0.28526*mean3[4]+0.26923*mean3[5]+0.28526*mean3[6]-0.19872*
        mean3[7]-0.21474*mean3[8]-0.19872*mean3[9]
      #ha.
      gs[4]<-0.03846*mean[1]+0.02564*mean[2]+0.03846*mean[3]-0.19872*
        mean1[1]+0.28526*mean1[2]-0.21474*mean1[3]+0.26923*mean1[4]+
        0.26923*mean2[1]-0.21474*mean2[2]+0.28526*mean2[3]-0.19872*
        mean2[4]-0.19872*mean3[1]+0.28526*mean3[2]-0.19872*mean3[3]-
        0.21474*mean3[4]+0.26923*mean3[5]-0.21474*mean3[6]-0.19872*
        mean3[7]+0.28526*mean3[8]-0.19872*mean3[9]
      #hb.
      g_aa1<-(0.5*gs[2]*gs[1]+0.25*gs[4]*gs[4])/n_fam #   0.5db**2+0.25hb**2.
      g_aa2<-(0.5*gs[1]*gs[1]+0.25*gs[3]*gs[3])/n_fam #   0.5da**2+0.25ha**2.
      g_aa3<-g_aa1+g_aa2
      sigma1[2]<-sigma1[1]+g_aa1;sigma1[3]<-sigma1[1]+g_aa2
      sigma1[4]<-sigma1[1]+g_aa3;sigma2[1]<-sigma2[4]+g_aa3
      sigma2[2]<-sigma2[4]+g_aa2;sigma2[3]<-sigma2[4]+g_aa1
      sigma3[2]<-sigma3[1]+g_aa1;sigma3[3]<-sigma3[7]<-sigma3[9]<-sigma3[1]
      sigma3[4]<-sigma3[1]+g_aa2;sigma3[5]<-sigma3[1]+g_aa3
      sigma3[6]<-sigma3[1]+g_aa2;sigma3[8]<-sigma3[1]+g_aa1
      ########################################################
      hh[1,1]<-sigma1[1]/n0[1]+sigma1[2]/n0[2]+sigma3[1]/s0[1]+sigma3[2]/s0[2]
      hh[1,2]<-0
      hh[1,3]<-sigma1[1]/n0[1]+sigma3[1]/s0[1]
      hh[1,4]<-0
      hh[1,5]<-0
      hh[1,6]<-0
      hh[1,7]<--sigma1[1]/n0[1]
      hh[1,8]<-0
      hh[1,9]<-2*sigma3[2]/s0[2]
      hh[1,10]<--sigma3[1]/s0[1]
      hh[1,11]<-0
      hh[1,12]<--sigma1[1]/n0[1]+sigma1[2]/n0[2]
      hh[1,13]<--sigma3[1]/s0[1]-sigma3[2]/s0[2]
      hh[2,2]<-sigma1[3]/n0[3]+sigma1[4]/n0[4]+sigma3[4]/s0[4]+sigma3[5]/s0[5]
      hh[2,3]<-sigma1[4]/n0[4]+sigma3[5]/s0[5]
      hh[2,4]<--sigma3[5]/s0[5]
      hh[2,5]<-0
      hh[2,6]<-hh[2,4]
      hh[2,7]<-sigma1[4]/n0[4]
      hh[2,8]<--sigma1[4]/n0[4]-2*sigma3[5]/s0[5]
      hh[2,9]<--2*sigma3[4]/s0[4]-8*sigma3[5]/s0[5]
      hh[2,10]<-0
      hh[2,11]<--2*sigma1[3]/n0[3]+4*sigma3[4]/s0[4]
      hh[2,12]<--sigma1[3]/n0[3]+sigma1[4]/n0[4]
      hh[2,13]<-sigma3[4]/s0[4]
      hh[3,3]<-sigma1[1]/n0[1]+sigma1[4]/n0[4]+sigma3[1]/s0[1]+sigma3[5]/s0[5]
      hh[3,4]<--sigma3[5]/s0[5]
      hh[3,5]<-0
      hh[3,6]<-hh[3,4]
      hh[3,7]<--sigma1[1]/n0[1]+sigma1[4]/n0[4]
      hh[3,8]<--sigma1[4]/n0[4]-2*sigma3[5]/s0[5]
      hh[3,9]<--8*sigma3[5]/s0[5]
      hh[3,10]<--sigma3[1]/s0[1]
      hh[3,11]<-0
      hh[3,12]<--sigma1[1]/n0[1]+sigma1[4]/n0[4]
      hh[3,13]<--sigma3[1]/s0[1]
      hh[4,4]<-sigma2[1]/n0[5]+sigma2[2]/n0[6]+sigma3[5]/s0[5]+sigma3[6]/s0[6]
      hh[4,5]<-0
      hh[4,6]<-sigma2[1]/n0[5]+sigma3[5]/s0[5]
      hh[4,7]<-sigma2[1]/n0[5]
      hh[4,8]<-sigma2[1]/n0[5]+2*sigma3[5]/s0[5]
      hh[4,9]<-8*sigma3[5]/s0[5]+2*sigma3[6]/s0[6]
      hh[4,10]<-0
      hh[4,11]<--2*sigma2[2]/n0[6]+4*sigma3[6]/s0[6]
      hh[4,12]<--sigma2[1]/n0[5]+sigma2[2]/n0[6]
      hh[4,13]<-sigma3[6]/s0[6]
      hh[5,5]<-sigma2[3]/n0[7]+sigma2[4]/n0[8]+sigma3[8]/s0[8]+sigma3[9]/s0[9]
      hh[5,6]<-sigma2[4]/n0[8]+sigma3[9]/s0[9]
      hh[5,7]<--sigma2[4]/n0[8]
      hh[5,8]<-0
      hh[5,9]<--2*sigma3[8]/s0[8]
      hh[5,10]<-sigma3[9]/s0[9]
      hh[5,11]<--4*sigma3[9]/s0[9]
      hh[5,12]<--sigma2[3]/n0[7]+sigma2[4]/n0[8]
      hh[5,13]<--sigma3[8]/s0[8]-sigma3[9]/s0[9]
      hh[6,6]<-sigma2[1]/n0[5]+sigma2[4]/n0[8]+sigma3[5]/s0[5]+sigma3[9]/s0[9]
      hh[6,7]<-sigma2[1]/n0[5]-sigma2[4]/n0[8]
      hh[6,8]<-sigma2[1]/n0[5]+2*sigma3[5]/s0[5]
      hh[6,9]<-8*sigma3[5]/s0[5]
      hh[6,10]<-sigma3[9]/s0[9]
      hh[6,11]<--4*sigma3[9]/s0[9]
      hh[6,12]<--sigma2[1]/n0[5]+sigma2[4]/n0[8]
      hh[6,13]<--sigma3[9]/s0[9]
      hh[7,7]<-sigma*(1/n_samP1+1/n_samP2)+sigma1[1]/n0[1]+sigma1[4]/n0[4]+sigma2[1]/n0[5]+sigma2[4]/n0[8]
      hh[7,8]<--sigma1[4]/n0[4]+sigma2[1]/n0[5]
      hh[7,9]<-sigma*(3/n_samP1-3/n_samP2)
      hh[7,10]<-0
      hh[7,11]<-sigma*(1/n_samP1+1/n_samP2)
      hh[7,12]<-sigma*(3/n_samP1-3/n_samP2)+sigma1[1]/n0[1]+sigma1[4]/n0[4]-sigma2[1]/n0[5]-sigma2[4]/n0[8]
      hh[7,13]<-0
      hh[8,8]<-sigma1[4]/n0[4]+sigma2[1]/n0[5]+4*sigma3[5]/s0[5]
      hh[8,9]<-16*sigma3[5]/s0[5]
      hh[8,10]<-0
      hh[8,11]<-0
      hh[8,12]<--sigma1[4]/n0[4]-sigma2[1]/n0[5]
      hh[8,13]<-0
      hh[9,9]<-sigma*(9/n_samP1+4/n_samF1+9/n_samP2)+4*sigma3[2]/s0[2]+16*sigma3[3]/s0[3]+4*sigma3[4]/s0[4]+64*sigma3[5]/s0[5]+4*sigma3[6]/s0[6]+16*sigma3[7]/s0[7]+4*sigma3[8]/s0[8]
      hh[9,10]<-4*sigma3[3]/s0[3]+4*sigma3[7]/s0[7]
      hh[9,11]<-sigma*(3/n_samP1-3/n_samP2)+4*sigma3[3]/s0[3]-8*sigma3[4]/s0[4]+8*sigma3[6]/s0[6]-20*sigma3[7]/s0[7]
      hh[9,12]<-sigma*(9/n_samP1+4/n_samF1+9/n_samP2)
      hh[9,13]<--2*sigma3[2]/s0[2]-2*sigma3[4]/s0[4]+2*sigma3[6]/s0[6]+2*sigma3[8]/s0[8]
      hh[10,10]<-sigma3[1]/s0[1]+sigma3[3]/s0[3]+sigma3[7]/s0[7]+sigma3[9]/s0[9]
      hh[10,11]<-sigma3[3]/s0[3]-5*sigma3[7]/s0[7]-4*sigma3[9]/s0[9]
      hh[10,12]<-0
      hh[10,13]<-sigma3[1]/s0[1]-sigma3[9]/s0[9]
      hh[11,11]<-sigma*(1/n_samP1+1/n_samP2)+4*sigma1[3]/n0[3]+4*sigma2[2]/n0[6]+sigma3[3]/s0[3]+16*sigma3[4]/s0[4]+16*sigma3[6]/s0[6]+25*sigma3[7]/s0[7]+16*sigma3[9]/s0[9]
      hh[11,12]<-sigma*(3/n_samP1-3/n_samP2)+2*sigma1[3]/n0[3]-2*sigma2[2]/n0[6]
      hh[11,13]<-4*(sigma3[4]/s0[4]+sigma3[6]/s0[6]+sigma3[9]/s0[9])
      hh[12,12]<-sigma*(9/n_samP1+4/n_samF1+9/n_samP2)+sigma1[1]/n0[1]+sigma1[2]/n0[2]+sigma1[3]/n0[3]+sigma1[4]/n0[4]+sigma2[1]/n0[5]+sigma2[2]/n0[6]+sigma2[3]/n0[7]+sigma2[4]/n0[8]
      hh[12,13]<-0
      hh[13,13]<-sigma3[1]/s0[1]+sigma3[2]/s0[2]+sigma3[4]/s0[4]+sigma3[6]/s0[6]+sigma3[8]/s0[8]+sigma3[9]/s0[9]
      for(i in 2:13)
      {
        for(j in 1:(i-1))
        {
          hh[i,j]<-hh[j,i]
        }
      }
      b_line[1]<-sumwx1[1]/n0[1]-sumwx1[2]/n0[2]-sumwx3[1]/s0[1]+sumwx3[2]/s0[2]
      b_line[2]<-sumwx1[3]/n0[3]-sumwx1[4]/n0[4]-sumwx3[4]/s0[4]+sumwx3[5]/s0[5]
      b_line[3]<-sumwx1[1]/n0[1]-sumwx1[4]/n0[4]-sumwx3[1]/s0[1]+sumwx3[5]/s0[5]
      b_line[4]<-sumwx2[1]/n0[5]-sumwx2[2]/n0[6]-sumwx3[5]/s0[5]+sumwx3[6]/s0[6]
      b_line[5]<-sumwx2[3]/n0[7]-sumwx2[4]/n0[8]-sumwx3[8]/s0[8]+sumwx3[9]/s0[9]
      b_line[6]<-sumwx2[1]/n0[5]-sumwx2[4]/n0[8]-sumwx3[5]/s0[5]+sumwx3[9]/s0[9]
      b_line[7]<-sumx[1]/n_samP1-sumx[3]/n_samP2-sumwx1[1]/n0[1]-sumwx1[4]/n0[4]+sumwx2[1]/n0[5]+sumwx2[4]/n0[8]
      b_line[8]<-sumwx1[4]/n0[4]+sumwx2[1]/n0[5]-2*sumwx3[5]/s0[5]
      b_line[9]<-3*sumx[1]/n_samP1+2*sumx[2]/n_samF1+3*sumx[3]/n_samP2+2*sumwx3[2]/s0[2]-4*sumwx3[3]/s0[3]+2*sumwx3[4]/s0[4]+2*sumwx3[6]/s0[6]-4*sumwx3[7]/s0[7]+2*sumwx3[8]/s0[8]-8*sumwx3[5]/s0[5]
      b_line[10]<-sumwx3[1]/s0[1]-sumwx3[3]/s0[3]-sumwx3[7]/s0[7]+sumwx3[9]/s0[9]
      b_line[11]<-sumx[1]/n_samP1-sumx[3]/n_samP2-2*sumwx1[3]/n0[3]+2*sumwx2[2]/n0[6]-sumwx3[3]/s0[3]-4*sumwx3[4]/s0[4]+4*sumwx3[6]/s0[6]+5*sumwx3[7]/s0[7]-4*sumwx3[9]/s0[9]
      b_line[12]<-3*sumx[1]/n_samP1+2*sumx[2]/n_samF1+3*sumx[3]/n_samP2-sumwx1[1]/n0[1]-sumwx1[2]/n0[2]-sumwx1[3]/n0[3]-sumwx1[4]/n0[4]-sumwx2[1]/n0[5]-sumwx2[2]/n0[6]-sumwx2[3]/n0[7]-sumwx2[4]/n0[8]
      b_line[13]<-sumwx3[1]/s0[1]-sumwx3[2]/s0[2]-sumwx3[4]/s0[4]+sumwx3[6]/s0[6]+sumwx3[8]/s0[8]-sumwx3[9]/s0[9]
      B<-solve(hh,b_line)
      mean[1]<-(sumx[1]+sigma*(-B[7]-3*B[9]-B[11]-3*B[12]))/n_samP1
      mean[2]<-(sumx[2]+sigma*(-2*B[9]-2*B[12]))/n_samF1
      mean[3]<-(sumx[3]+sigma*(B[7]-3*B[9]+B[11]-3*B[12]))/n_samP2
      mean1[1]<-(sumwx1[1]+sigma1[1]*(-B[1]-B[3]+B[7]+B[12]))/n0[1]
      mean1[2]<-(sumwx1[2]+sigma1[2]*(B[1]+B[12]))/n0[2]
      mean1[3]<-(sumwx1[3]+sigma1[3]*(-B[2]+2*B[11]+B[12]))/n0[3]
      mean1[4]<-(sumwx1[4]+sigma1[4]*(B[2]+B[3]+B[7]-B[8]+B[12]))/n0[4]
      mean2[1]<-(sumwx2[1]+sigma2[1]*(-B[4]-B[6]-B[7]-B[8]+B[12]))/n0[5]
      mean2[2]<-(sumwx2[2]+sigma2[2]*(B[4]-2*B[11]+B[12]))/n0[6]
      mean2[3]<-(sumwx2[3]+sigma2[3]*(-B[5]+B[12]))/n0[7]
      mean2[4]<-(sumwx2[4]+sigma2[4]*(B[5]+B[6]-B[7]+B[12]))/n0[8]
      mean3[1]<-(sumwx3[1]+sigma3[1]*(B[1]+B[3]-B[10]-B[13]))/s0[1]
      mean3[2]<-(sumwx3[2]+sigma3[2]*(-B[1]-2*B[9]+B[13]))/s0[2]
      mean3[3]<-(sumwx3[3]+sigma3[3]*(4*B[9]+B[10]+B[11]))/s0[3]
      mean3[7]<-(sumwx3[7]+sigma3[7]*(4*B[9]+B[10]-5*B[11]))/s0[7]
      mean3[4]<-(sumwx3[4]+sigma3[4]*(B[2]-2*B[9]+4*B[11]+B[13]))/s0[4]
      mean3[5]<-(sumwx3[5]+sigma3[5]*(-B[2]-B[3]+B[4]+B[6]+2*B[8]+8*B[9]))/s0[5]
      mean3[6]<-(sumwx3[6]+sigma3[6]*(-B[4]-2*B[9]-4*B[11]-B[13]))/s0[6]
      mean3[8]<-(sumwx3[8]+sigma3[8]*(B[5]-2*B[9]-B[13]))/s0[8]
      mean3[9]<-(sumwx3[9]+sigma3[9]*(-B[5]-B[6]-B[10]+4*B[11]+B[13]))/s0[9]
      aaa1<-max(abs(B-AA))
      AA<-B
      if (n_iter>20) break
    }
    ##############obtain variance#################
    ss1<-sum((dataP1-mean[1])^2);ss3<-sum((dataP2-mean[3])^2);ss2<-sum((dataF1-mean[2])^2)
    for(i in 1:4) {swx1[i] <- W1[i,]%*%(dataB1-mean1[i])^2 };for(i in 1:4) {swx2[i] <- W2[i,]%*%(dataB2-mean2[i])^2 };for(i in 1:9) {swx3[i] <- W3[i,]%*%(dataF2-mean3[i])^2 }
    # to estimate sigma40.
    aaa0<-sigma1[1];n_iter<-0;aaa1<-1000
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      aa2<-sigma1[1]/(sigma1[1]+g_aa1)
      aa3<-sigma1[1]/(sigma1[1]+g_aa2)
      aa4<-sigma1[1]/(sigma1[1]+g_aa3)
      as1<-swx1[1]+swx1[2]*aa2^2+swx1[3]*aa3^2+swx1[4]*aa4^2
      as2<-n0[1]+aa2*n0[2]+aa3*n0[3]+aa4*n0[4]
      sigma1[1]<-as1/as2
      aaa1<-abs(sigma1[1]-aaa0)
      aaa0<-sigma1[1]
      if (n_iter>20) break
    }
    sigma40<-sigma1[1]-sigma;
    if (sigma40<0) {sigma40<-0;sigma1[1]<-sigma}
    sigma1[1]<-sigma40+sigma;sigma1[2]<-sigma1[1]+g_aa1
    sigma1[3]<-sigma1[1]+g_aa2;sigma1[4]<-sigma1[1]+g_aa3
    # to estimate sigma50.
    aaa0<-sigma2[4];n_iter<-0;aaa1<-1000
    while (aaa1>0.001)
    {
      n_iter<-n_iter+1
      aa1<-sigma2[4]/(sigma2[4]+g_aa3)
      aa2<-sigma2[4]/(sigma2[4]+g_aa2)
      aa3<-sigma2[4]/(sigma2[4]+g_aa1)
      as3<-swx2[1]*aa1^2+swx2[2]*aa2^2+swx2[3]*aa3^2+swx2[4]
      as4<-aa1*n0[5]+aa2*n0[6]+aa3*n0[7]+n0[8]
      sigma2[4]<-as3/as4
      aaa1<-abs(sigma2[4]-aaa0)
      aaa0<-sigma2[4]
      if (n_iter>20) break
    }
    sigma50<-sigma2[4]-sigma;
    if (sigma50<0) {sigma50<-0;sigma2[4]<-sigma}
    sigma2[4]<-sigma50+sigma;sigma2[1]<-sigma2[4]+g_aa3
    sigma2[2]<-sigma2[4]+g_aa2;sigma2[3]<-sigma2[4]+g_aa1
    # to estimate sigma60.
    aaa0<-sigma3[1];aa6<-swx3[1]+swx3[3]+swx3[7]+swx3[9]
    aa7<-s0[1]+s0[3]+s0[7]+s0[9]
    n_iter<-0;aaa1<-1000
    while (aaa1>0.001)
    {
      n_iter<-n_iter+1
      aa1<-sigma3[1]/(sigma3[1]+g_aa1)
      aa2<-sigma3[1]/(sigma3[1]+g_aa2)
      aa3<-sigma3[1]/(sigma3[1]+g_aa3)
      as5<-aa6+(swx3[2]+swx3[8])*aa1^2+(swx3[4]+swx3[6])*aa2^2+swx3[5]*aa3^2
      as6<-aa7+aa1*(s0[2]+s0[8])+aa2*(s0[4]+s0[6])+aa3*s0[5]
      sigma3[1]<-as5/as6
      aaa1<-abs(sigma3[1]-aaa0)
      aaa0<-sigma3[1]
      if (n_iter>20) break
    }
    sigma60<-sigma3[1]-sigma;
    if (sigma60<0) {sigma60<-0;sigma3[1]<-sigma}
    sigma3[1]<-sigma60+sigma;sigma3[2]<-sigma3[1]+g_aa1
    sigma3[3]<-sigma3[7]<-sigma3[9]<-sigma3[1];sigma3[4]<-sigma3[1]+g_aa2
    sigma3[5]<-sigma3[1]+g_aa3;sigma3[6]<-sigma3[1]+g_aa2;sigma3[8]<-sigma3[1]+g_aa1
    # CM3 to estimate the variance (sigma).
    ab1<-ss1+ss2+ss3;ab2<-n_samP1+n_samF1+n_samP2
    n_iter<-0;aaa0<-sigma;aaa1<-1000
    while (aaa1>0.001)
    {
      n_iter<-n_iter+1
      n0[11]<-sigma/(sigma+sigma40)
      n0[12]<-sigma/(sigma+sigma40+g_aa1)
      n0[13]<-sigma/(sigma+sigma40+g_aa2)
      n0[14]<-sigma/(sigma+sigma40+g_aa3)
      s0[11]<-sigma/(sigma+sigma50+g_aa3)
      s0[12]<-sigma/(sigma+sigma50+g_aa2)
      s0[13]<-sigma/(sigma+sigma50+g_aa1)
      s0[14]<-sigma/(sigma+sigma50)
      ab3<-sum(swx1[c(1:4)]*n0[c(11:14)]^2+swx2[c(1:4)]*s0[c(11:14)]^2)
      ab4<-sum(n0[c(1:4)]*n0[c(11:14)]+n0[c(5:8)]*s0[c(11:14)])
      n0[11]<-sigma/(sigma+sigma60);n0[13]<-n0[17]<-n0[19]<-n0[11]
      n0[12]<-sigma/(sigma+sigma60+g_aa1);n0[14]<-sigma/(sigma+sigma60+g_aa2)
      n0[15]<-sigma/(sigma+sigma60+g_aa3);n0[16]<-sigma/(sigma+sigma60+g_aa2)
      n0[18]<-sigma/(sigma+sigma60+g_aa1)
      ab3<-ab3+sum(swx3[c(1:9)]*n0[c(11:19)]^2)
      ab4<-ab4+sum(s0[c(1:9)]*n0[11:19])
      sigma<-(ab1+ab3)/(ab2+ab4)
      aaa1<-abs(sigma-aaa0)
      aaa0<-sigma
      if (n_iter>20) break
    }
    sigma1[1]<-sigma+sigma40;sigma1[2]<-sigma1[1]+g_aa1
    sigma1[3]<-sigma1[1]+g_aa2;sigma1[4]<-sigma1[1]+g_aa3
    sigma2[4]<-sigma+sigma50;sigma2[1]<-sigma2[4]+g_aa3
    sigma2[2]<-sigma2[4]+g_aa2;sigma2[3]<-sigma2[4]+g_aa1
    sigma3[1]<-sigma+sigma60;sigma3[3]<-sigma3[7]<-sigma3[9]<-sigma3[1]
    sigma3[2]<-sigma3[1]+g_aa1;sigma3[4]<-sigma3[1]+g_aa2
    sigma3[5]<-sigma3[1]+g_aa3;sigma3[6]<-sigma3[1]+g_aa2
    sigma3[8]<-sigma3[1]+g_aa1
    if(sum(sigma < 1e-30)>=1){break}
    ##############################criteria for iterations to stop########################
    L1<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,mean[3],sqrt(sigma))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mix_pi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mix_pi3)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*11
  meanP1<-mean[1];meanF1<-mean[2];meanP2<-mean[3]
  sigma0<-sigma

  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*n_samP1)
  P1bmw <- matrix(0,n_samP1,1)
  P1gg <- (dataP1 - mean[1])/sqrt(as.vector(sigma))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < n_samP1){P1bmw <- P1bmw+runif(n_samP1)/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  P1u<- as.matrix(c(12*n_samP1*((P1dd[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((P1dd[2]/n_samP1-1/3)^2),180*n_samP1*((P1dd[3]/n_samP1-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[2])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),P1D))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*n_samF1)
  F1bmw <- matrix(0,n_samF1,1)
  F1gg <- (dataF1 - mean[2])/sqrt(as.vector(sigma))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < n_samF1){F1bmw <- F1bmw+runif(n_samF1)/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  F1u<- as.matrix(c(12*n_samF1*((F1dd[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((F1dd[2]/n_samF1-1/3)^2),180*n_samF1*((F1dd[3]/n_samF1-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[2])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),F1D))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*n_samP2)
  P2bmw <- matrix(0,n_samP2,1)
  P2gg <- (dataP2 - mean[3])/sqrt(as.vector(sigma))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < n_samP2){P2bmw <- P2bmw+runif(n_samP2)/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  P2u<- as.matrix(c(12*n_samP2*((P2dd[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((P2dd[2]/n_samP2-1/3)^2),180*n_samP2*((P2dd[3]/n_samP2-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[2])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),P2D))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1);
  B1w1<-1/(12*n_samB1)
  B1bmw <- matrix(0,n_samB1,1); B1bmwsl <- matrix(0,n_samB1,4)
  for(i in 1:4){
    B1gg <- (dataB1 - mean1[i])/sqrt(sigma1[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi1[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B1P2)))[1]
  if(nn < n_samB1){B1P2 <- B1P2+runif(n_samB1)/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*n_samB1) + sum((B1P2 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  B1u <- as.matrix(c(12*n_samB1*((B1dd[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((B1dd[2]/n_samB1-1/3)^2),180*n_samB1*((B1dd[3]/n_samB1-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[2])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),B1D))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])

  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2);
  B2w1<-1/(12*n_samB2)
  B2bmw <- matrix(0,n_samB2,1); B2bmwsl <- matrix(0,n_samB2,4)
  for(i in 1:4){
    B2gg <- (dataB2 - mean2[i])/sqrt(sigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B2P2)))[1]
  if(nn < n_samB2){B2P2 <- B2P2+runif(n_samB2)/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*n_samB2) + sum((B2P2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  B2u <- as.matrix(c(12*n_samB2*((B2dd[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((B2dd[2]/n_samB2-1/3)^2),180*n_samB2*((B2dd[3]/n_samB2-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[2])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),B2D))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])

  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2);
  F2w1<-1/(12*n_samF2)
  F2bmw <- matrix(0,n_samF2,1); F2bmwsl <- matrix(0,n_samF2,9)
  for(i in 1:9){
    F2gg <- (dataF2 - mean3[i])/sqrt(sigma3[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi3[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < n_samF2){F2P2 <- F2P2+runif(n_samF2)/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*n_samF2) + sum((F2P2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  F2u <- as.matrix(c(12*n_samF2*((F2dd[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((F2dd[2]/n_samF2-1/3)^2),180*n_samF2*((F2dd[3]/n_samF2-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[2])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),F2D))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])

  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)

  ############first order parameters#########################################
  aa<-matrix(c(1,1,1,0,0,1,0,1,0,0,1,1,0,1,1,-1,-1,0,0,-1,0,1,1,1,0,0,0.5,0.25,1,1,0,0,0.5,0.5,0.25,1,0,1,0.5,0,0.5,0.25,
               1,0,0,0.5,0.5,0.5,0.25,1,0,0,0.5,0.5,-0.5,0.25,1,0,-1,0.5,0,-0.5,0.25,1,-1,0,0,0.5,-0.5,0.25,1,-1,-1,0,0,-0.5,
               0.25,1,1,1,0,0,0,0.25,1,1,0,0,0.5,0,0.25,1,1,-1,0,0,0,0.25,1,0,1,0.5,0,0,0.25,1,0,0,0.5,0.5,0,0.25,1,0,-1,0.5,
               0,0,0.25,1,-1,1,0,0,0,0.25,1,-1,0,0,0.5,0,0.25,1,-1,-1,0,0,0,0.25),20,7,byrow=T)
  b_line1<-as.matrix(c(mean[1],mean[2],mean[3],mean1[1],mean1[2],mean1[3],mean1[4],mean2[1],mean2[2],mean2[3],mean2[4],
                       mean3[1],mean3[2],mean3[3],mean3[4],mean3[5],mean3[6],mean3[7],mean3[8],mean3[9]))
  B1<-solve(crossprod(aa,aa))%*%crossprod(aa,b_line1)
  #########second order parameters####################
  jj1<-sigmaB1-sigma1[1]
  if (jj1<0 || jj1>=sigmaB1) {jj1<-0}
  ll1<-jj1/sigmaB1
  mm1<-sigma1[1]-sigma
  if (mm1<0 || mm1>=sigmaB1) {mm1<-0}
  nn1<-mm1/sigmaB1
  jj2<-sigmaB2-sigma2[4]
  if (jj2<0 || jj2>=sigmaB2) {jj2<-0}
  ll2<-jj2/sigmaB2
  mm2<-sigma2[4]-sigma
  if (mm2<0 || mm2>=sigmaB2) {mm2<-0}
  nn2<-mm2/sigmaB2
  jj3<-sigmaF2-sigma3[1]
  if (jj3<0 || jj3>=sigmaF2) {jj3<-0}
  ll3<-jj3/sigmaF2
  mm3<-sigma3[1]-sigma
  if (mm3<0 || mm3>=sigmaF2) {mm3<-0}
  nn3<-mm3/sigmaF2

  output <- data.frame("MX2-AD-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4), round(sigma0,4),round(t(mean1),4),round(t(sigma1),4),
                       round(t(mix_pi1),4),round(t(mean2),4),round(t(sigma2),4),round(t(mix_pi2),4),
                       round(t(mean3),4),round(t(sigma3),4),round(t(mix_pi3),4),
                       round(B1[1],4)," "," "," "," "," ",round(B1[2],4),round(B1[3],4),round(B1[4],4),round(B1[5],4)," "," "," "," ",round(B1[6],4),round(B1[7],4),
                       round(jj1,4),round(ll1*100,4),round(mm1,4),round(nn1*100,4),round(jj2,4),round(ll2*100,4),round(mm2,4),round(nn2*100,4),round(jj3,4),round(ll3*100,4),round(mm3,4),round(nn3*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi1,mi2,mi3)
  return(OUTPUT)
}

#####################MX2-A-AD(E-3)####################################
G6FModelFun[[21]] <- function(K1,logL,df11,df21,df31,df41,df51,df61,G6Ftext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1];n_samF1<-dim(dataF1)[1];n_samP2<-dim(dataP2)[1]
  n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1];n_samF2<-dim(dataF2)[1]
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2),sum(dataF2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2),sum(dataF2^2)))
  ss<-matrix(0,3,1);ss[1]<-s[1]-sumx[1]^2/n_samP1;ss[2]<-s[2]-sumx[2]^2/n_samF1;ss[3]<-s[3]-sumx[3]^2/n_samP2
  mean<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataB1),mean(dataB2),mean(dataF2)))
  sigma0<-sum(ss)/(n_samP1+n_samF1+n_samP2-3)
  sigmaB1<-var(dataB1);sigma40<-sigmaB1
  sigmaB2<-var(dataB2);sigma50<-sigmaB2
  sigmaF2<-var(dataF2);sigma60<-sigmaF2
  sigmaP1<-sigmaF1<-sigmaP2<-sigma0
  m_esp <- 0.0001;n_fam <- as.numeric(G6Ftext2)
  ###############procedure start###########################
  mi1<-matrix(0.25,4,1);sigma1<-matrix(0,4,1)
  mi2<-matrix(0.25,4,1);sigma2<-matrix(0,4,1)
  mi3<-as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  sigma3<-matrix(0,9,1)
  sigma<-sigma0
  a1<-sqrt(sigma40/n_samB1);if (mean[1]<mean[3]) a1<--a1
  mean1<-as.matrix(c(mean[1],mean[4]-0.5*a1,mean[4]-0.8*a1,mean[2]))
  a2<-sqrt(sigma50/n_samB2);if (mean[1]<mean[3]) a2<--a2
  mean2<-as.matrix(c(mean[2],mean[5]-0.5*a2,mean[5]-a2,mean[3]))
  a3<-sqrt(sigma60/n_samF2);if (mean[1]<mean[3]) a3<--a3
  mean3<-matrix(0,9,1)
  mean3[1]<-mean3[2]<-mean[1];mean3[3]<-(mean[1]+mean[3])/2
  mean3[4]<-mean[6]+0.6*a3;mean3[5]<-mean[6]
  mean3[6]<-mean3[4];mean3[7]<-mean3[3]
  mean3[8]<-mean[6]-a3;mean3[9]<-mean[3]
  sigma1[1]<-sigmaB1/2;sigma2[4]<-sigmaB2/2;sigma3[1]<-sigmaF2/2
  gs<-matrix(0,2,1)
  gs[1]<-(mean1[1]+mean1[2]-mean1[3]-mean1[4]+mean2[1]+mean2[2]-
            mean2[3]-mean2[4]+2*mean3[1]+2*mean3[2]+2*mean3[3]-2*mean3[7]-
            2*mean3[8]-2*mean3[9])/16
  #da.
  gs[2]<-(mean1[1]-mean1[2]+mean1[3]-mean1[4]+mean2[1]-mean2[2]+
            mean2[3]-mean2[4]+2*mean3[1]-2*mean3[3]+2*mean3[4]-2*
            mean3[6]+2*mean3[7]-2*mean3[9])/16
  #db.
  g_aa1<-0.5*gs[2]*gs[2]/n_fam #   0.5db**2.
  g_aa2<-0.5*gs[1]*gs[1]/n_fam #   0.5da**2.
  g_aa3<-g_aa1+g_aa2
  sigma1[2]<-sigma1[1]+g_aa1;sigma1[3]<-sigma1[1]+g_aa2
  sigma1[4]<-sigma1[1]+g_aa3;sigma2[1]<-sigma2[4]+g_aa3
  sigma2[2]<-sigma2[4]+g_aa2;sigma2[3]<-sigma2[4]+g_aa1
  sigma3[2]<-sigma3[1]+g_aa1;sigma3[3]<-sigma3[7]<-sigma3[9]<-sigma3[1]
  sigma3[4]<-sigma3[1]+g_aa2;sigma3[5]<-sigma3[1]+g_aa3
  sigma3[6]<-sigma3[1]+g_aa2;sigma3[8]<-sigma3[1]+g_aa1
  L0<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma0))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma0))))+
    sum(log(dnorm(dataP2,mean[3],sqrt(sigma0))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mi1)))+
    sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mi3)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,4,n_samB1);      swx1 <- matrix(0,4,1)
  W2 <- matrix(0,4,n_samB2);      swx2 <- matrix(0,4,1)
  W3 <- matrix(0,9,n_samF2);      swx3 <- matrix(0,9,1)
  hh<-matrix(0,15,15);b_line<-matrix(0,15,1)
  n0<-matrix(0,18,1);s0<-matrix(0,18,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:4) { W1[i,] <- mi1[i]*dnorm(dataB1,mean1[i],sqrt(sigma1[i]))/dmixnorm(dataB1,mean1,sqrt(sigma1),mi1) }
    mix_pi1 <- as.matrix(rowSums(W1)/n_samB1)
    sumwx1 <- W1%*%dataB1
    for(i in 1:4) { W2[i,] <- mi2[i]*dnorm(dataB2,mean2[i],sqrt(sigma2[i]))/dmixnorm(dataB2,mean2,sqrt(sigma2),mi2) }
    mix_pi2 <- as.matrix(rowSums(W2)/n_samB2)
    sumwx2 <- W2%*%dataB2
    for(i in 1:9) { W3[i,] <- mi3[i]*dnorm(dataF2,mean3[i],sqrt(sigma3[i]))/dmixnorm(dataF2,mean3,sqrt(sigma3),mi3) }
    mix_pi3 <- as.matrix(rowSums(W3)/n_samF2)
    sumwx3 <- W3%*%dataF2
    ############obtain means####################
    n0[c(1:4)]<-mix_pi1[c(1:4)]*n_samB1;n0[c(5:8)]<-mix_pi2[c(1:4)]*n_samB2
    s0[c(1:9)]<-mix_pi3[c(1:9)]*n_samF2
    s0[c(1:9)][abs(s0[c(1:9)])<0.00000001]<-0.000001
    n0[c(1:8)][abs(n0[c(1:8)])<0.00000001]<-0.000001
    aaa0<-0 ;AA<-matrix(0,15,1);n_iter<-0;aaa1<-1000
    while(aaa1>0.001)
    {
      n_iter<-n_iter+1
      gs[1]<-(mean1[1]+mean1[2]-mean1[3]-mean1[4]+mean2[1]+mean2[2]-
                mean2[3]-mean2[4]+2*mean3[1]+2*mean3[2]+2*mean3[3]-2*mean3[7]-
                2*mean3[8]-2*mean3[9])/16
      #da.
      gs[2]<-(mean1[1]-mean1[2]+mean1[3]-mean1[4]+mean2[1]-mean2[2]+
                mean2[3]-mean2[4]+2*mean3[1]-2*mean3[3]+2*mean3[4]-2*
                mean3[6]+2*mean3[7]-2*mean3[9])/16
      #db.
      g_aa1<-0.5*gs[2]*gs[2]/n_fam #   0.5db**2.
      g_aa2<-0.5*gs[1]*gs[1]/n_fam #   0.5da**2.
      g_aa3<-g_aa1+g_aa2
      sigma1[2]<-sigma1[1]+g_aa1;sigma1[3]<-sigma1[1]+g_aa2
      sigma1[4]<-sigma1[1]+g_aa3;sigma2[1]<-sigma2[4]+g_aa3
      sigma2[2]<-sigma2[4]+g_aa2;sigma2[3]<-sigma2[4]+g_aa1
      sigma3[2]<-sigma3[1]+g_aa1;sigma3[3]<-sigma3[7]<-sigma3[9]<-sigma3[1]
      sigma3[4]<-sigma3[1]+g_aa2;sigma3[5]<-sigma3[1]+g_aa3
      sigma3[6]<-sigma3[1]+g_aa2;sigma3[8]<-sigma3[1]+g_aa1
      hh[1,1]<-sigma*(9/n_samP1+4/n_samF1+9/n_samP2)+64*sigma3[5]/s0[5]
      hh[1,2]<-0
      hh[1,3]<--8*sigma3[5]/s0[5]
      hh[1,4]<-hh[1,3]
      hh[1,5]<--hh[1,3]
      hh[1,6]<-0
      hh[1,7]<-hh[1,5]
      hh[1,8]<-sigma*(3/n_samP1-3/n_samP2)
      hh[1,9]<-16*sigma3[5]/s0[5]
      hh[1,10]<-0
      hh[1,11]<-0
      hh[1,12]<-0
      hh[1,13]<-0
      hh[1,14]<-48*sigma3[5]/s0[5]
      hh[1,15]<-0
      hh[2,2]<-sigma1[1]/n0[1]+sigma1[2]/n0[2]+sigma3[1]/s0[1]+sigma3[2]/s0[2]
      hh[2,3]<-0
      hh[2,4]<-sigma1[1]/n0[1]+sigma3[1]/s0[1]
      hh[2,5]<-0
      hh[2,6]<-0
      hh[2,7]<-hh[2,6]
      hh[2,8]<--sigma1[1]/n0[1]
      hh[2,9]<-0
      hh[2,10]<-sigma3[2]/s0[2]
      hh[2,11]<--sigma3[1]/s0[1]
      hh[2,12]<-0
      hh[2,13]<-hh[2,11]
      hh[2,14]<-sigma1[1]/n0[1]-2*sigma1[2]/n0[2]
      hh[2,15]<--sigma1[2]/n0[2]
      hh[3,3]<-sigma1[3]/n0[3]+sigma1[4]/n0[4]+sigma3[4]/s0[4]+sigma3[5]/s0[5]
      hh[3,4]<-sigma1[4]/n0[4]+sigma3[5]/s0[5]
      hh[3,5]<--sigma3[5]/s0[5]
      hh[3,6]<-0
      hh[3,7]<--sigma3[5]/s0[5]
      hh[3,8]<-sigma1[4]/n0[4]
      hh[3,9]<--sigma1[4]/n0[4]-2*sigma3[5]/s0[5]
      hh[3,10]<-0
      hh[3,11]<-0
      hh[3,12]<--sigma3[4]/s0[4]
      hh[3,13]<-2*sigma3[4]/s0[4]
      hh[3,14]<--6*sigma3[5]/s0[5]
      hh[3,15]<--sigma1[3]/n0[3]
      hh[4,4]<-sigma1[1]/n0[1]+sigma1[4]/n0[4]+sigma3[1]/s0[1]+sigma3[5]/s0[5]
      hh[4,5]<--sigma3[5]/s0[5]
      hh[4,6]<-0
      hh[4,7]<-hh[4,5]
      hh[4,8]<--sigma1[1]/n0[1]+sigma1[4]/n0[4]
      hh[4,9]<--sigma1[4]/n0[4]-2*sigma3[5]/s0[5]
      hh[4,10]<-0
      hh[4,11]<--sigma3[1]/s0[1]
      hh[4,12]<-0
      hh[4,13]<--sigma3[1]/s0[1]
      hh[4,14]<-sigma1[1]/n0[1]-6*sigma3[5]/s0[5]
      hh[4,15]<-0
      hh[5,5]<-sigma2[1]/n0[5]+sigma2[2]/n0[6]+sigma3[5]/s0[5]+sigma3[6]/s0[6]
      hh[5,6]<-0
      hh[5,7]<-sigma2[1]/n0[5]+sigma3[5]/s0[5]
      hh[5,8]<-sigma2[1]/n0[5]
      hh[5,9]<-sigma2[1]/n0[5]+2*sigma3[5]/s0[5]
      hh[5,10]<-0
      hh[5,11]<-0
      hh[5,12]<--sigma3[6]/s0[6]
      hh[5,13]<-0
      hh[5,14]<-6*sigma3[5]/s0[5]
      hh[5,15]<-sigma2[2]/n0[6]
      hh[6,6]<-sigma2[3]/n0[7]+sigma2[4]/n0[8]+sigma3[8]/s0[8]+sigma3[9]/s0[9]
      hh[6,7]<-sigma2[4]/n0[8]+sigma3[9]/s0[9]
      hh[6,8]<--sigma2[4]/n0[8]
      hh[6,9]<-0
      hh[6,10]<--sigma3[8]/s0[8]
      hh[6,11]<-sigma3[9]/s0[9]
      hh[6,12]<-hh[6,11]
      hh[6,13]<--2*sigma3[8]/s0[8]-sigma3[9]/s0[9]
      hh[6,14]<-2*sigma2[3]/n0[7]-sigma2[4]/n0[8]
      hh[6,15]<-sigma2[3]/n0[7]
      hh[7,7]<-sigma2[1]/n0[5]+sigma2[4]/n0[8]+sigma3[5]/s0[5]+sigma3[9]/s0[9]
      hh[7,8]<-sigma2[1]/n0[5]-sigma2[4]/n0[8]
      hh[7,9]<-sigma2[1]/n0[5]+2*sigma3[5]/s0[5]
      hh[7,10]<-0
      hh[7,11]<-sigma3[9]/s0[9]
      hh[7,12]<-hh[7,11]
      hh[7,13]<--hh[7,11]
      hh[7,14]<--sigma2[4]/n0[8]+6*sigma3[5]/s0[5]
      hh[7,15]<-0
      hh[8,8]<-sigma*(1/n_samP1+1/n_samP2)+sigma1[1]/n0[1]+sigma1[4]/n0[4]+sigma2[1]/n0[5]+sigma2[4]/n0[8]
      hh[8,9]<--sigma1[4]/n0[4]+sigma2[1]/n0[5]
      hh[8,10]<-0
      hh[8,11]<-hh[8,10]
      hh[8,12]<-hh[8,10]
      hh[8,13]<-hh[8,10]
      hh[8,14]<--sigma1[1]/n0[1]+sigma2[4]/n0[8]
      hh[8,15]<-0
      hh[9,9]<-sigma1[4]/n0[4]+sigma2[1]/n0[5]+4*sigma3[5]/s0[5]
      hh[9,10]<-0
      hh[9,11]<-0
      hh[9,12]<-0
      hh[9,13]<-0
      hh[9,14]<-12*sigma3[5]/s0[5]
      hh[9,15]<-0
      hh[10,10]<-sigma3[2]/s0[2]+sigma3[3]/s0[3]+sigma3[7]/s0[7]+sigma3[8]/s0[8]
      hh[10,11]<-sigma3[3]/s0[3]+sigma3[7]/s0[7]
      hh[10,12]<-sigma3[7]/s0[7]
      hh[10,13]<-2*sigma3[8]/s0[8]
      hh[10,14]<-0
      hh[10,15]<-0
      hh[11,11]<-sigma3[1]/s0[1]+sigma3[3]/s0[3]+sigma3[7]/s0[7]+sigma3[9]/s0[9]
      hh[11,12]<-sigma3[7]/s0[7]+sigma3[9]/s0[9]
      hh[11,13]<-sigma3[1]/s0[1]-sigma3[9]/s0[9]
      hh[11,14]<-0
      hh[11,15]<-0
      hh[12,12]<-sigma3[4]/s0[4]+sigma3[6]/s0[6]+sigma3[7]/s0[7]+sigma3[9]/s0[9]
      hh[12,13]<--2*sigma3[4]/s0[4]-sigma3[9]/s0[9]
      hh[12,14]<-0
      hh[12,15]<-0
      hh[13,13]<-sigma3[1]/s0[1]+4*sigma3[4]/s0[4]+4*sigma3[8]/s0[8]+sigma3[9]/s0[9]
      hh[13,14]<-0
      hh[13,15]<-0
      hh[14,14]<-sigma1[1]/n0[1]+4*sigma1[2]/n0[2]+4*sigma2[3]/n0[7]+sigma2[4]/n0[8]+36*sigma3[5]/s0[5]
      hh[14,15]<-2*sigma1[2]/n0[2]+2*sigma2[3]/n0[7]
      hh[15,15]<-sigma1[2]/n0[2]+sigma1[3]/n0[3]+sigma2[2]/n0[6]+sigma2[3]/n0[7]
      for(i in 2:15)
      {
        for(j in 1:(i-1))
        {
          hh[i,j]<-hh[j,i]
        }
      }
      b_line[1]<-3*sumx[1]/n_samP1+2*sumx[2]/n_samF1+3*sumx[3]/n_samP2-8*sumwx3[5]/s0[5]
      b_line[2]<-sumwx1[1]/n0[1]-sumwx1[2]/n0[2]-sumwx3[1]/s0[1]+sumwx3[2]/s0[2]
      b_line[3]<-sumwx1[3]/n0[3]-sumwx1[4]/n0[4]-sumwx3[4]/s0[4]+sumwx3[5]/s0[5]
      b_line[4]<-sumwx1[1]/n0[1]-sumwx1[4]/n0[4]-sumwx3[1]/s0[1]+sumwx3[5]/s0[5]
      b_line[5]<-sumwx2[1]/n0[5]-sumwx2[2]/n0[6]-sumwx3[5]/s0[5]+sumwx3[6]/s0[6]
      b_line[6]<-sumwx2[3]/n0[7]-sumwx2[4]/n0[8]-sumwx3[8]/s0[8]+sumwx3[9]/s0[9]
      b_line[7]<-sumwx2[1]/n0[5]-sumwx2[4]/n0[8]-sumwx3[5]/s0[5]+sumwx3[9]/s0[9]
      b_line[8]<-sumx[1]/n_samP1-sumx[3]/n_samP2-sumwx1[1]/n0[1]-sumwx1[4]/n0[4]+sumwx2[1]/n0[5]+sumwx2[4]/n0[8]
      b_line[9]<-sumwx1[4]/n0[4]+sumwx2[1]/n0[5]-2*sumwx3[5]/s0[5]
      b_line[10]<-sumwx3[2]/s0[2]-sumwx3[3]/s0[3]-sumwx3[7]/s0[7]+sumwx3[8]/s0[8]
      b_line[11]<-sumwx3[1]/s0[1]-sumwx3[3]/s0[3]-sumwx3[7]/s0[7]+sumwx3[9]/s0[9]
      b_line[12]<-sumwx3[4]/s0[4]-sumwx3[6]/s0[6]-sumwx3[7]/s0[7]+sumwx3[9]/s0[9]
      b_line[13]<-sumwx3[1]/s0[1]-2*sumwx3[4]/s0[4]+2*sumwx3[8]/s0[8]-sumwx3[9]/s0[9]
      b_line[14]<-sumwx1[1]/n0[1]+2*sumwx1[2]/n0[2]+2*sumwx2[3]/n0[7]+sumwx2[4]/n0[8]-6*sumwx3[5]/s0[5]
      b_line[15]<-sumwx1[2]/n0[2]-sumwx1[3]/n0[3]-sumwx2[2]/n0[6]+sumwx2[3]/n0[7]
      B<-solve(hh,b_line)
      mean[1]<-(sumx[1]+sigma*(-3*B[1]-B[8]))/n_samP1
      mean[2]<-(sumx[2]+sigma*(-2*B[1]))/n_samF1
      mean[3]<-(sumx[3]+sigma*(-3*B[1]+B[8]))/n_samP2
      mean1[1]<-(sumwx1[1]+sigma1[1]*(-B[2]-B[4]+B[8]-B[14]))/n0[1]
      mean1[2]<-(sumwx1[2]+sigma1[2]*(B[2]-2*B[14]-B[15]))/n0[2]
      mean1[3]<-(sumwx1[3]+sigma1[3]*(-B[3]+B[15]))/n0[3]
      mean1[4]<-(sumwx1[4]+sigma1[4]*(B[3]+B[4]+B[8]-B[9]))/n0[4]
      mean2[1]<-(sumwx2[1]+sigma2[1]*(-B[5]-B[7]-B[8]-B[9]))/n0[5]
      mean2[2]<-(sumwx2[2]+sigma2[2]*(B[5]+B[15]))/n0[6]
      mean2[3]<-(sumwx2[3]+sigma2[3]*(-B[6]-2*B[14]-B[15]))/n0[7]
      mean2[4]<-(sumwx2[4]+sigma2[4]*(B[6]+B[7]-B[8]-B[14]))/n0[8]
      mean3[1]<-(sumwx3[1]+sigma3[1]*(B[2]+B[4]-B[11]-B[13]))/s0[1]
      mean3[2]<-(sumwx3[2]+sigma3[2]*(-B[2]-B[10]))/s0[2]
      mean3[3]<-(sumwx3[3]+sigma3[3]*(B[10]+B[11]))/s0[3]
      mean3[7]<-(sumwx3[7]+sigma3[7]*(B[10]+B[11]+B[12]))/s0[7]
      mean3[4]<-(sumwx3[4]+sigma3[4]*(B[3]-B[12]+2*B[13]))/s0[4]
      mean3[5]<-(sumwx3[5]+sigma3[5]*(8*B[1]-B[3]-B[4]+B[5]+B[7]+2*B[9]+6*B[14]))/s0[5]
      mean3[6]<-(sumwx3[6]+sigma3[6]*(-B[5]+B[12]))/s0[6]
      mean3[8]<-(sumwx3[8]+sigma3[8]*(B[6]-B[10]-2*B[13]))/s0[8]
      mean3[9]<-(sumwx3[9]+sigma3[9]*(-B[6]-B[7]-B[11]-B[12]+B[13]))/s0[9]
      aaa1<-max(abs(B-AA))
      AA<-B
      if (n_iter>20) break
    }
    ##########obtain variance###############################
    ss1<-sum((dataP1-mean[1])^2);ss3<-sum((dataP2-mean[3])^2);ss2<-sum((dataF1-mean[2])^2)
    for(i in 1:4) {swx1[i] <- W1[i,]%*%(dataB1-mean1[i])^2 };for(i in 1:4) {swx2[i] <- W2[i,]%*%(dataB2-mean2[i])^2 };for(i in 1:9) {swx3[i] <- W3[i,]%*%(dataF2-mean3[i])^2 }
    # to estimate sigma40.
    aaa0<-sigma1[1];n_iter<-0;aaa1<-1000
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      aa2<-sigma1[1]/(sigma1[1]+g_aa1)
      aa3<-sigma1[1]/(sigma1[1]+g_aa2)
      aa4<-sigma1[1]/(sigma1[1]+g_aa3)
      as1<-swx1[1]+swx1[2]*aa2^2+swx1[3]*aa3^2+swx1[4]*aa4^2
      as2<-n0[1]+aa2*n0[2]+aa3*n0[3]+aa4*n0[4]
      sigma1[1]<-as1/as2
      aaa1<-abs(sigma1[1]-aaa0)
      aaa0<-sigma1[1]
      if (n_iter>20) break
    }
    sigma40<-sigma1[1]-sigma;
    if (sigma40<0) {sigma40<-0;sigma1[1]<-sigma}
    sigma1[1]<-sigma40+sigma;sigma1[2]<-sigma1[1]+g_aa1
    sigma1[3]<-sigma1[1]+g_aa2;sigma1[4]<-sigma1[1]+g_aa3
    # to estimate sigma50.
    aaa0<-sigma2[4];n_iter<-0;aaa1<-1000
    while (aaa1>0.001)
    {
      n_iter<-n_iter+1
      aa1<-sigma2[4]/(sigma2[4]+g_aa3)
      aa2<-sigma2[4]/(sigma2[4]+g_aa2)
      aa3<-sigma2[4]/(sigma2[4]+g_aa1)
      as3<-swx2[1]*aa1^2+swx2[2]*aa2^2+swx2[3]*aa3^2+swx2[4]
      as4<-aa1*n0[5]+aa2*n0[6]+aa3*n0[7]+n0[8]
      sigma2[4]<-as3/as4
      aaa1<-abs(sigma2[4]-aaa0)
      aaa0<-sigma2[4]
      if (n_iter>20) break
    }
    sigma50<-sigma2[4]-sigma;
    if (sigma50<0) {sigma50<-0;sigma2[4]<-sigma}
    sigma2[4]<-sigma50+sigma;sigma2[1]<-sigma2[4]+g_aa3
    sigma2[2]<-sigma2[4]+g_aa2;sigma2[3]<-sigma2[4]+g_aa1
    # to estimate sigma60.
    aaa0<-sigma3[1];aa6<-swx3[1]+swx3[3]+swx3[7]+swx3[9]
    aa7<-s0[1]+s0[3]+s0[7]+s0[9]
    n_iter<-0;aaa1<-1000
    while (aaa1>0.001)
    {
      n_iter<-n_iter+1
      aa1<-sigma3[1]/(sigma3[1]+g_aa1)
      aa2<-sigma3[1]/(sigma3[1]+g_aa2)
      aa3<-sigma3[1]/(sigma3[1]+g_aa3)
      as5<-aa6+(swx3[2]+swx3[8])*aa1^2+(swx3[4]+swx3[6])*aa2^2+swx3[5]*aa3^2
      as6<-aa7+aa1*(s0[2]+s0[8])+aa2*(s0[4]+s0[6])+aa3*s0[5]
      sigma3[1]<-as5/as6
      aaa1<-abs(sigma3[1]-aaa0)
      aaa0<-sigma3[1]
      if (n_iter>20) break
    }
    sigma60<-sigma3[1]-sigma;
    if (sigma60<0) {sigma60<-0;sigma3[1]<-sigma}
    sigma3[1]<-sigma60+sigma;sigma3[2]<-sigma3[1]+g_aa1
    sigma3[3]<-sigma3[7]<-sigma3[9]<-sigma3[1];sigma3[4]<-sigma3[1]+g_aa2
    sigma3[5]<-sigma3[1]+g_aa3;sigma3[6]<-sigma3[1]+g_aa2;sigma3[8]<-sigma3[1]+g_aa1
    # CM3 to estimate the variance (sigma).
    ab1<-ss1+ss2+ss3;ab2<-n_samP1+n_samF1+n_samP2
    n_iter<-0;aaa0<-sigma;aaa1<-1000
    while (aaa1>0.001)
    {
      n_iter<-n_iter+1
      n0[11]<-sigma/(sigma+sigma40)
      n0[12]<-sigma/(sigma+sigma40+g_aa1)
      n0[13]<-sigma/(sigma+sigma40+g_aa2)
      n0[14]<-sigma/(sigma+sigma40+g_aa3)
      s0[11]<-sigma/(sigma+sigma50+g_aa3)
      s0[12]<-sigma/(sigma+sigma50+g_aa2)
      s0[13]<-sigma/(sigma+sigma50+g_aa1)
      s0[14]<-sigma/(sigma+sigma50)
      ab3<-sum(swx1[c(1:4)]*n0[c(11:14)]^2+swx2[c(1:4)]*s0[c(11:14)]^2)
      ab4<-sum(n0[c(1:4)]*n0[c(11:14)]+n0[c(5:8)]*s0[c(11:14)])
      n0[11]<-sigma/(sigma+sigma60);n0[13]<-n0[17]<-n0[19]<-n0[11]
      n0[12]<-sigma/(sigma+sigma60+g_aa1);n0[14]<-sigma/(sigma+sigma60+g_aa2)
      n0[15]<-sigma/(sigma+sigma60+g_aa3);n0[16]<-sigma/(sigma+sigma60+g_aa2)
      n0[18]<-sigma/(sigma+sigma60+g_aa1)
      ab3<-ab3+sum(swx3[c(1:9)]*n0[c(11:19)]^2)
      ab4<-ab4+sum(s0[c(1:9)]*n0[11:19])
      sigma<-(ab1+ab3)/(ab2+ab4)
      aaa1<-abs(sigma-aaa0)
      aaa0<-sigma
      if (n_iter>20) break
    }
    sigma1[1]<-sigma+sigma40;sigma1[2]<-sigma1[1]+g_aa1
    sigma1[3]<-sigma1[1]+g_aa2;sigma1[4]<-sigma1[1]+g_aa3
    sigma2[4]<-sigma+sigma50;sigma2[1]<-sigma2[4]+g_aa3
    sigma2[2]<-sigma2[4]+g_aa2;sigma2[3]<-sigma2[4]+g_aa1
    sigma3[1]<-sigma+sigma60;sigma3[3]<-sigma3[7]<-sigma3[9]<-sigma3[1]
    sigma3[2]<-sigma3[1]+g_aa1;sigma3[4]<-sigma3[1]+g_aa2
    sigma3[5]<-sigma3[1]+g_aa3;sigma3[6]<-sigma3[1]+g_aa2;sigma3[8]<-sigma3[1]+g_aa1
    if(sum(sigma < 1e-30)>=1){break}
    ##############################criteria for iterations to stop########################
    L1<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,mean[3],sqrt(sigma))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mix_pi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mix_pi3)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*9
  meanP1<-mean[1];meanF1<-mean[2];meanP2<-mean[3]
  sigma0<-sigma
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*n_samP1)
  P1bmw <- matrix(0,n_samP1,1)
  P1gg <- (dataP1 - mean[1])/sqrt(as.vector(sigma))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < n_samP1){P1bmw <- P1bmw+runif(n_samP1)/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  P1u<- as.matrix(c(12*n_samP1*((P1dd[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((P1dd[2]/n_samP1-1/3)^2),180*n_samP1*((P1dd[3]/n_samP1-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[2])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),P1D))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*n_samF1)
  F1bmw <- matrix(0,n_samF1,1)
  F1gg <- (dataF1 - mean[2])/sqrt(as.vector(sigma))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < n_samF1){F1bmw <- F1bmw+runif(n_samF1)/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  F1u<- as.matrix(c(12*n_samF1*((F1dd[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((F1dd[2]/n_samF1-1/3)^2),180*n_samF1*((F1dd[3]/n_samF1-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[2])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),F1D))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*n_samP2)
  P2bmw <- matrix(0,n_samP2,1)
  P2gg <- (dataP2 - mean[3])/sqrt(as.vector(sigma))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < n_samP2){P2bmw <- P2bmw+runif(n_samP2)/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  P2u<- as.matrix(c(12*n_samP2*((P2dd[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((P2dd[2]/n_samP2-1/3)^2),180*n_samP2*((P2dd[3]/n_samP2-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[2])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),P2D))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1);
  B1w1<-1/(12*n_samB1)
  B1bmw <- matrix(0,n_samB1,1); B1bmwsl <- matrix(0,n_samB1,4)
  for(i in 1:4){
    B1gg <- (dataB1 - mean1[i])/sqrt(sigma1[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi1[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B1P2)))[1]
  if(nn < n_samB1){B1P2 <- B1P2+runif(n_samB1)/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*n_samB1) + sum((B1P2 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  B1u <- as.matrix(c(12*n_samB1*((B1dd[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((B1dd[2]/n_samB1-1/3)^2),180*n_samB1*((B1dd[3]/n_samB1-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[2])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),B1D))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])

  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2);
  B2w1<-1/(12*n_samB2)
  B2bmw <- matrix(0,n_samB2,1); B2bmwsl <- matrix(0,n_samB2,4)
  for(i in 1:4){
    B2gg <- (dataB2 - mean2[i])/sqrt(sigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B2P2)))[1]
  if(nn < n_samB2){B2P2 <- B2P2+runif(n_samB2)/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*n_samB2) + sum((B2P2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  B2u <- as.matrix(c(12*n_samB2*((B2dd[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((B2dd[2]/n_samB2-1/3)^2),180*n_samB2*((B2dd[3]/n_samB2-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[2])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),B2D))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])

  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2);
  F2w1<-1/(12*n_samF2)
  F2bmw <- matrix(0,n_samF2,1); F2bmwsl <- matrix(0,n_samF2,9)
  for(i in 1:9){
    F2gg <- (dataF2 - mean3[i])/sqrt(sigma3[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi3[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < n_samF2){F2P2 <- F2P2+runif(n_samF2)/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*n_samF2) + sum((F2P2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  F2u <- as.matrix(c(12*n_samF2*((F2dd[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((F2dd[2]/n_samF2-1/3)^2),180*n_samF2*((F2dd[3]/n_samF2-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[2])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),F2D))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])

  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)

  ##################first order parameters######################
  aa<-matrix(c(1,1,1,1,0,1,0,0,0,1,1,-1,-1,-1,0,1,1,1,0.5,0.25,1,1,0,0.5,0.25,1,0,1,0.5,0.25,
               1,0,0,0.5,0.25,1,0,0,-0.5,0.25,1,0,-1,-0.5,0.25,1,-1,0,-0.5,0.25,1,-1,-1,-0.5,
               0.25,1,1,1,0,0.25,1,1,0,0,0.25,1,1,-1,0,0.25,1,0,1,0,0.25,1,0,0,0,0.25,1,0,-1,0,
               0.25,1,-1,1,0,0.25,1,-1,0,0,0.25,1,-1,-1,0,0.25),20,5,byrow=T)
  b_line1<-as.matrix(c(mean[1],mean[2],mean[3],mean1[1],mean1[2],mean1[3],mean1[4],mean2[1],mean2[2],mean2[3],mean2[4],
                       mean3[1],mean3[2],mean3[3],mean3[4],mean3[5],mean3[6],mean3[7],mean3[8],mean3[9]))
  B1<-solve(crossprod(aa,aa))%*%crossprod(aa,b_line1)
  ###############second order parameters########################
  jj1<-sigmaB1-sigma1[1]
  if (jj1<0 || jj1>=sigmaB1) {jj1<-0}
  ll1<-jj1/sigmaB1
  mm1<-sigma1[1]-sigma
  if (mm1<0 || mm1>=sigmaB1) {mm1<-0}
  nn1<-mm1/sigmaB1
  jj2<-sigmaB2-sigma2[4]
  if (jj2<0 || jj2>=sigmaB2) {jj2<-0}
  ll2<-jj2/sigmaB2
  mm2<-sigma2[4]-sigma
  if (mm2<0 || mm2>=sigmaB2) {mm2<-0}
  nn2<-mm2/sigmaB2
  jj3<-sigmaF2-sigma3[1]
  if (jj3<0 || jj3>=sigmaF2) {jj3<-0}
  ll3<-jj3/sigmaF2
  mm3<-sigma3[1]-sigma
  if (mm3<0 || mm3>=sigmaF2) {mm3<-0}
  nn3<-mm3/sigmaF2

  output <- data.frame("MX2-A-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4), round(sigma0,4),round(t(mean1),4),round(t(sigma1),4),
                       round(t(mix_pi1),4),round(t(mean2),4),round(t(sigma2),4),round(t(mix_pi2),4),
                       round(t(mean3),4),round(t(sigma3),4),round(t(mix_pi3),4),
                       round(B1[1],4)," "," "," "," "," ",round(B1[2],4),round(B1[3],4)," "," "," "," "," "," ",round(B1[4],4),round(B1[5],4),
                       round(jj1,4),round(ll1*100,4),round(mm1,4),round(nn1*100,4),round(jj2,4),round(ll2*100,4),round(mm2,4),round(nn2*100,4),round(jj3,4),round(ll3*100,4),round(mm3,4),round(nn3*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi1,mi2,mi3)
  return(OUTPUT)
}

###############MX2-EA-AD(E-4)##################################
G6FModelFun[[22]] <- function(K1,logL,df11,df21,df31,df41,df51,df61,G6Ftext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1];n_samF1<-dim(dataF1)[1];n_samP2<-dim(dataP2)[1]
  n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1];n_samF2<-dim(dataF2)[1]
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2),sum(dataF2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2),sum(dataF2^2)))
  ss<-matrix(0,3,1);ss[1]<-s[1]-sumx[1]^2/n_samP1;ss[2]<-s[2]-sumx[2]^2/n_samF1;ss[3]<-s[3]-sumx[3]^2/n_samP2
  mean<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataB1),mean(dataB2),mean(dataF2)))
  sigma0<-sum(ss)/(n_samP1+n_samF1+n_samP2-3)
  sigmaB1<-var(dataB1);sigma40<-sigmaB1
  sigmaB2<-var(dataB2);sigma50<-sigmaB2
  sigmaF2<-var(dataF2);sigma60<-sigmaF2
  sigmaP1<-sigmaF1<-sigmaP2<-sigma0
  m_esp <- 0.0001;n_fam <- as.numeric(G6Ftext2)
  ###############procedure start###########################
  mi1<-as.matrix(c(0.25,0.5,0.25));sigma1<-matrix(0,3,1)
  mi2<-as.matrix(c(0.25,0.5,0.25));sigma2<-matrix(0,3,1)
  mi3<-as.matrix(c(0.0625,0.25,0.125,0.25,0.25,0.0625))
  sigma3<-matrix(0,6,1)
  sigma<-sigma0
  a1<-sqrt(sigmaB1/n_samB1);if (mean[1]<mean[3]) {a1<--a1}
  mean1<-as.matrix(c(mean[1],mean[4],mean[4]-1.2*a1))
  a2<-sqrt(sigmaB2/n_samB2);if (mean[1]<mean[3]) {a2<--a2}
  mean2<-as.matrix(c(mean[2]+0.5*a1,mean[5],mean[5]-1.2*a2))
  a3<-sqrt(sigmaF2/n_samF2);if (mean[1]<mean[3]){a3<--a3}
  mean3<-as.matrix(c(mean[1],mean[4],mean[2],mean[2],mean[5],mean[2]))
  sigma1[1]<-sigmaB1/2;sigma2[3]<-sigmaB2/2;sigma3[1]<-sigmaF2/2
  gs<-matrix(0,1,1)
  gs[1]<-(mean1[1]-mean1[3]+mean2[1]-mean2[3]+2*mean3[1]+mean3[2]-mean3[5]-2*mean3[6])/14     #d.
  g_aa1<-0.5*gs[1]^2/n_fam  #   0.5d**2.
  g_aa2<-gs[1]^2/n_fam      #   d**2.
  sigma1[2]<-sigma1[1]+g_aa1;sigma1[3]<-sigma1[1]+g_aa2
  sigma2[1]<-sigma2[3]+g_aa2;sigma2[2]<-sigma2[3]+g_aa1
  sigma3[2]<-sigma3[1]+g_aa1;sigma3[3]<-sigma3[6]<-sigma3[1]
  sigma3[4]<-sigma3[1]+g_aa2;sigma3[5]<-sigma3[1]+g_aa1
  L0<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma0))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma0))))+
    sum(log(dnorm(dataP2,mean[3],sqrt(sigma0))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mi1)))+
    sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mi3)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,3,n_samB1);      swx1 <- matrix(0,3,1)
  W2 <- matrix(0,3,n_samB2);      swx2 <- matrix(0,3,1)
  W3 <- matrix(0,6,n_samF2);      swx3 <- matrix(0,6,1)
  hh<-matrix(0,10,10);b_line<-matrix(0,10,1)
  n0<-matrix(0,15,1);s0<-matrix(0,15,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:3) { W1[i,] <- mi1[i]*dnorm(dataB1,mean1[i],sqrt(sigma1[i]))/dmixnorm(dataB1,mean1,sqrt(sigma1),mi1) }
    mix_pi1 <- as.matrix(rowSums(W1)/n_samB1)
    sumwx1 <- W1%*%dataB1
    for(i in 1:3) { W2[i,] <- mi2[i]*dnorm(dataB2,mean2[i],sqrt(sigma2[i]))/dmixnorm(dataB2,mean2,sqrt(sigma2),mi2) }
    mix_pi2 <- as.matrix(rowSums(W2)/n_samB2)
    sumwx2 <- W2%*%dataB2
    for(i in 1:6) { W3[i,] <- mi3[i]*dnorm(dataF2,mean3[i],sqrt(sigma3[i]))/dmixnorm(dataF2,mean3,sqrt(sigma3),mi3) }
    mix_pi3 <- as.matrix(rowSums(W3)/n_samF2)
    sumwx3 <- W3%*%dataF2
    n0[c(1:3)]<-mix_pi1[c(1:3)]*n_samB1;n0[c(4:6)]<-mix_pi2[c(1:3)]*n_samB2
    s0[c(1:6)]<-mix_pi3[c(1:6)]*n_samF2
    s0[c(1:6)][abs(s0[c(1:6)])<0.00000001]<-0.000001
    n0[c(1:6)][abs(n0[c(1:6)])<0.00000001]<-0.000001
    ################# CM1-step for means####################
    aaa0<-0;AA<-matrix(0,10,1);n_iter<-0
    aaa1<-1000
    while(aaa1>0.001)
    {
      n_iter<-n_iter+1;
      gs[1]<-(mean1[1]-mean1[3]+mean2[1]-mean2[3]+2*mean3[1]+mean3[2]-mean3[5]-2*mean3[6])/14     #d.
      g_aa1<-0.5*gs[1]^2/n_fam
      g_aa2<-gs[1]^2/n_fam
      sigma1[2]<-sigma1[1]+g_aa1;sigma1[3]<-sigma1[1]+g_aa2
      sigma2[1]<-sigma2[3]+g_aa2;sigma2[2]<-sigma2[3]+g_aa1
      sigma3[2]<-sigma3[1]+g_aa1;sigma3[3]<-sigma3[6]<-sigma3[1]
      sigma3[4]<-sigma3[1]+g_aa2;sigma3[5]<-sigma3[1]+g_aa1
      aa11<-sigma3[4]*sumwx3[3]+sigma3[3]*sumwx3[4]
      aa12<-sigma3[4]*s0[3]+sigma3[3]*s0[4]
      aa13<-sigma3[3]*sigma3[4]
      hh[1,1]<-sigma*(9/n_samP1+4/n_samF1+9/n_samP2)+64*aa13/aa12
      hh[1,2]<-sigma*(3/n_samP1-3/n_samP2)
      hh[1,3]<-hh[1,2]
      hh[1,4]<-hh[1,2]
      hh[1,5]<--8*aa13/aa12
      hh[1,6]<-0
      hh[1,7]<-16*aa13/aa12
      hh[1,8]<-hh[1,7]
      hh[1,9]<-0.5*hh[1,7]
      hh[1,10]<-hh[1,7]
      hh[2,2]<-sigma*(1/n_samP1+1/n_samP2)+4*sigma1[1]/n0[1]+36*sigma1[2]/n0[2]+36*sigma3[5]/s0[5]+4*sigma3[6]/s0[6]
      hh[2,3]<-sigma*(1/n_samP1+1/n_samP2)-2*sigma1[1]/n0[1]
      hh[2,4]<-sigma*(1/n_samP1+1/n_samP2)+12*sigma1[2]/n0[2]
      hh[2,5]<-2*sigma1[1]/n0[1]
      hh[2,6]<-2*sigma1[1]/n0[1]+6*sigma1[2]/n0[2]
      hh[2,7]<--6*sigma1[2]/n0[2]
      hh[2,8]<-0
      hh[2,9]<-6*sigma3[5]/s0[5]
      hh[2,10]<--2*sigma3[6]/s0[6]
      hh[3,3]<-sigma*(1/n_samP1+1/n_samP2)+sigma1[1]/n0[1]+sigma1[3]/n0[3]+sigma2[1]/n0[4]+sigma2[3]/n0[6]
      hh[3,4]<-sigma*(1/n_samP1+1/n_samP2)
      hh[3,5]<--sigma1[1]/n0[1]+sigma1[3]/n0[3]
      hh[3,6]<--sigma1[1]/n0[1]
      hh[3,7]<-0
      hh[3,8]<--sigma1[3]/n0[3]+sigma2[1]/n0[4]
      hh[3,9]<-sigma2[1]/n0[4]
      hh[3,10]<-0
      hh[4,4]<-sigma*(1/n_samP1+1/n_samP2)+4*sigma1[2]/n0[2]+4*sigma2[2]/n0[5]
      hh[4,5]<-0
      hh[4,6]<-2*sigma1[2]/n0[2]
      hh[4,7]<--2*sigma1[2]/n0[2]+2*sigma2[2]/n0[5]
      hh[4,8]<-0
      hh[4,9]<--2*sigma2[2]/n0[5]
      hh[4,10]<-0
      hh[5,5]<-sigma1[1]/n0[1]+sigma1[3]/n0[3]+sigma3[1]/s0[1]+aa13/aa12
      hh[5,6]<-sigma1[1]/n0[1]+sigma3[1]/s0[1]
      hh[5,7]<--2*aa13/aa12
      hh[5,8]<--sigma1[3]/n0[3]-2*aa13/aa12
      hh[5,9]<--aa13/aa12
      hh[5,10]<--sigma3[1]/s0[1]-2*aa13/aa12
      hh[6,6]<-sigma1[1]/n0[1]+sigma1[2]/n0[2]+sigma3[1]/s0[1]+sigma3[2]/s0[2]
      hh[6,7]<--sigma1[2]/n0[2]
      hh[6,8]<-0
      hh[6,9]<-0
      hh[6,10]<--sigma3[1]/s0[1]
      hh[7,7]<-sigma1[2]/n0[2]+sigma2[2]/n0[5]+4*aa13/aa12
      hh[7,8]<-4*aa13/aa12
      hh[7,9]<--sigma2[2]/n0[5]+2*aa13/aa12
      hh[7,10]<-4*aa13/aa12
      hh[8,8]<-sigma1[3]/n0[3]+sigma2[1]/n0[4]+4*aa13/aa12
      hh[8,9]<-sigma2[1]/n0[4]+2*aa13/aa12
      hh[8,10]<-4*aa13/aa12
      hh[9,9]<-sigma2[1]/n0[4]+sigma2[2]/n0[5]+aa13/aa12+sigma3[5]/s0[5]
      hh[9,10]<-2*aa13/aa12
      hh[10,10]<-sigma3[1]/s0[1]+sigma3[6]/s0[6]+4*aa13/aa12
      for(i in 2:10)
      {
        for(j in 1:(i-1))
        {
          hh[i,j]<-hh[j,i]
        }
      }
      b_line[1]<-3*sumx[1]/n_samP1+2*sumx[2]/n_samF1+3*sumx[3]/n_samP2-8*aa11/aa12;
      b_line[2]<-sumx[1]/n_samP1-sumx[3]/n_samP2+2*sumwx1[1]/n0[1]-6*sumwx1[2]/n0[2]+6*sumwx3[5]/s0[5]-2*sumwx3[6]/s0[6];
      b_line[3]<-sumx[1]/n_samP1-sumx[3]/n_samP2-sumwx1[1]/n0[1]-sumwx1[3]/n0[3]+sumwx2[1]/n0[4]+sumwx2[3]/n0[6];
      b_line[4]<-sumx[1]/n_samP1-sumx[3]/n_samP2-2*sumwx1[2]/n0[2]+2*sumwx2[2]/n0[5];
      b_line[5]<-sumwx1[1]/n0[1]-sumwx1[3]/n0[3]-sumwx3[1]/s0[1]+aa11/aa12;
      b_line[6]<-sumwx1[1]/n0[1]-sumwx1[2]/n0[2]-sumwx3[1]/s0[1]+sumwx3[2]/s0[2];
      b_line[7]<-sumwx1[2]/n0[2]+sumwx2[2]/n0[5]-2*aa11/aa12;
      b_line[8]<-sumwx1[3]/n0[3]+sumwx2[1]/n0[4]-2*aa11/aa12;
      b_line[9]<-sumwx2[1]/n0[4]-sumwx2[2]/n0[5]-aa11/aa12+sumwx3[5]/s0[5];
      b_line[10]<-sumwx3[1]/s0[1]-2*aa11/aa12+sumwx3[6]/s0[6];
      B<-solve(hh,b_line)
      mean[1]<-(sumx[1]+sigma*(-3*B[1]-B[2]-B[3]-B[4]))/n_samP1
      mean[2]<-(sumx[2]-sigma*2*B[1])/n_samF1
      mean[3]<-(sumx[3]+sigma*(-3*B[1]+B[2]+B[3]+B[4]))/n_samP2
      mean1[1]<-(sumwx1[1]+sigma1[1]*(-2*B[2]+B[3]-B[5]-B[6]))/n0[1]
      mean1[2]<-(sumwx1[2]+sigma1[2]*(6*B[2]+2*B[4]+B[6]-B[7]))/n0[2]
      mean1[3]<-(sumwx1[3]+sigma1[3]*(B[3]+B[5]-B[8]))/n0[3]
      mean2[1]<-(sumwx2[1]+sigma2[1]*(-B[3]-B[8]-B[9]))/n0[4]
      mean2[2]<-(sumwx2[2]+sigma2[2]*(-2*B[4]-B[7]+B[9]))/n0[5]
      mean2[3]<-(sumwx2[3]-sigma2[3]*B[3])/n0[6]
      mean3[1]<-(sumwx3[1]+sigma3[1]*(B[5]+B[6]-B[10]))/s0[1]
      mean3[2]<-(sumwx3[2]+sigma3[2]*(-B[6]))/s0[2]
      mean3[3]<-(sigma3[4]*sumwx3[3]+sigma3[3]*sumwx3[4]+sigma3[3]*sigma3[4]*(8*B[1]-B[5]+2*B[7]+2*B[8]+B[9]+2*B[10]))/(sigma3[4]*s0[3]+sigma3[3]*s0[4])
      mean3[4]<-mean3[3]
      mean3[5]<-(sumwx3[5]+sigma3[5]*(-6*B[2]-B[9]))/s0[5]
      mean3[6]<-(sumwx3[6]+sigma3[6]*(2*B[2]-B[10]))/s0[6]
      aaa1<-max(abs(B-AA))
      AA<-B
      if (n_iter>20) break
    }
    ########obtain variance#############################
    ss1<-sum((dataP1-mean[1])^2);ss3<-sum((dataP2-mean[3])^2);ss2<-sum((dataF1-mean[2])^2)
    for(i in 1:3) {swx1[i] <- W1[i,]%*%(dataB1-mean1[i])^2 };for(i in 1:3) {swx2[i] <- W2[i,]%*%(dataB2-mean2[i])^2 };for(i in 1:6) {swx3[i] <- W3[i,]%*%(dataF2-mean3[i])^2 }
    aaa0<-sigma1[1];n_iter<-0;aaa1<-1000
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      aa2<-sigma1[1]/(sigma1[1]+g_aa1)
      aa3<-sigma1[1]/(sigma1[1]+g_aa2)
      as1<-swx1[1]+swx1[2]*aa2^2+swx1[3]*aa3^2
      as2<-n0[1]+aa2*n0[2]+aa3*n0[3]
      sigma1[1]<-as1/as2
      aaa1<-abs(sigma1[1]-aaa0)
      aaa0<-sigma1[1]
      if (n_iter>20) break
    }
    sigma40<-sigma1[1]-sigma;
    if (sigma40<0) {sigma40<-0;sigma1[1]<-sigma}
    sigma1[1]<-sigma+sigma40;sigma1[2]<-sigma1[1]+g_aa1;sigma1[3]<-sigma1[1]+g_aa2
    # to estimate sigma50.
    aaa0<-sigma2[3];n_iter<-0;aaa1<-1000
    while (aaa1>0.001)
    {
      n_iter<-n_iter+1
      aa1<-sigma2[3]/(sigma2[3]+g_aa2)
      aa2<-sigma2[3]/(sigma2[3]+g_aa1)
      as3<-swx2[1]*aa1^2+swx2[2]*aa2^2+swx2[3]
      as4<-aa1*n0[4]+aa2*n0[5]+n0[6]
      sigma2[3]<-as3/as4
      aaa1<-abs(sigma2[3]-aaa0)
      aaa0<-sigma2[3]
      if (n_iter>20) break
    }
    sigma50<-sigma2[3]-sigma;
    if (sigma50<0) {sigma50<-0;sigma2[3]<-sigma}
    sigma2[3]<-sigma+sigma50;sigma2[1]<-sigma2[3]+g_aa2;sigma2[2]<-sigma2[3]+g_aa1
    # to estimate sigma60.
    aaa0<-sigma3[1];aa6<-swx3[1]+swx3[3]+swx3[6];aa7<-s0[1]+s0[3]+s0[6]
    n_iter<-0;aaa1<-1000
    while (aaa1>0.001)
    {
      n_iter<-n_iter+1
      aa1<-sigma3[1]/(sigma3[1]+g_aa1)
      aa2<-sigma3[1]/(sigma3[1]+g_aa2)
      as5<-aa6+(swx3[2]+swx3[5])*aa1^2+swx3[4]*aa2^2
      as6<-aa7+aa1*(s0[2]+s0[5])+aa2*s0[4]
      sigma3[1]<-as5/as6
      aaa1<-abs(sigma3[1]-aaa0)
      aaa0<-sigma3[1]
      if (n_iter>20) break
    }
    sigma60<-sigma3[1]-sigma;
    if (sigma60<0) {sigma60<-0;sigma3[1]<-sigma}
    sigma3[1]<-sigma+sigma60;sigma3[2]<-sigma3[1]+g_aa1
    sigma3[3]<-sigma3[6]<-sigma3[1];sigma3[4]<-sigma3[1]+g_aa2;sigma3[5]<-sigma3[1]+g_aa1
    # CM3 to estimate the variance (sigma).
    ab1<-ss1+ss2+ss3;ab2<-n_samP1+n_samF1+n_samP2
    n_iter<-0;aaa0<-sigma;aaa1<-1000
    while (aaa1>0.001)
    {
      n_iter<-n_iter+1
      n0[11]<-sigma/(sigma+sigma40)
      n0[12]<-sigma/(sigma+sigma40+g_aa1)
      n0[13]<-sigma/(sigma+sigma40+g_aa2)
      s0[11]<-sigma/(sigma+sigma50+g_aa2)
      s0[12]<-sigma/(sigma+sigma50+g_aa1)
      s0[13]<-sigma/(sigma+sigma50)
      ab3<-sum(swx1[c(1:3)]*n0[c(11:13)]^2+swx2[c(1:3)]*s0[c(11:13)]^2)
      ab4<-sum(n0[c(1:3)]*n0[c(11:13)]+n0[c(4:6)]*s0[c(11:13)])
      n0[11]<-sigma/(sigma+sigma60);n0[13]<-n0[16]<-n0[11]
      n0[12]<-sigma/(sigma+sigma60+g_aa1)
      n0[14]<-sigma/(sigma+sigma60+g_aa2)
      n0[15]<-sigma/(sigma+sigma60+g_aa1)
      ab3<-ab3+sum(swx3[c(1:6)]*n0[c(11:16)]^2)
      ab4<-ab4+sum(s0[c(1:6)]*n0[11:16])
      sigma<-(ab1+ab3)/(ab2+ab4)
      aaa1<-abs(sigma-aaa0)
      aaa0<-sigma
      if (n_iter>20) break
    }
    sigma1[1]<-sigma+sigma40;sigma1[2]<-sigma1[1]+g_aa1
    sigma1[3]<-sigma1[1]+g_aa2;sigma2[3]<-sigma+sigma50
    sigma2[1]<-sigma2[3]+g_aa2;sigma2[2]<-sigma2[3]+g_aa1
    sigma3[1]<-sigma+sigma60;sigma3[2]<-sigma3[1]+g_aa1
    sigma3[3]<-sigma3[6]<-sigma3[1];sigma3[4]<-sigma3[1]+g_aa2
    sigma3[5]<-sigma3[1]+g_aa1
    if(sum(sigma < 1e-30)>=1){break}
    ##############################criteria for iterations to stop########################
    L1<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,mean[3],sqrt(sigma))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mix_pi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mix_pi3)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*8
  meanP1<-mean[1];meanF1<-mean[2];meanP2<-mean[3]
  sigma0<-sigma
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*n_samP1)
  P1bmw <- matrix(0,n_samP1,1)
  P1gg <- (dataP1 - mean[1])/sqrt(as.vector(sigma))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < n_samP1){P1bmw <- P1bmw+runif(n_samP1)/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  P1u<- as.matrix(c(12*n_samP1*((P1dd[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((P1dd[2]/n_samP1-1/3)^2),180*n_samP1*((P1dd[3]/n_samP1-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[2])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),P1D))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*n_samF1)
  F1bmw <- matrix(0,n_samF1,1)
  F1gg <- (dataF1 - mean[2])/sqrt(as.vector(sigma))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < n_samF1){F1bmw <- F1bmw+runif(n_samF1)/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  F1u<- as.matrix(c(12*n_samF1*((F1dd[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((F1dd[2]/n_samF1-1/3)^2),180*n_samF1*((F1dd[3]/n_samF1-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[2])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),F1D))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*n_samP2)
  P2bmw <- matrix(0,n_samP2,1)
  P2gg <- (dataP2 - mean[3])/sqrt(as.vector(sigma))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < n_samP2){P2bmw <- P2bmw+runif(n_samP2)/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  P2u<- as.matrix(c(12*n_samP2*((P2dd[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((P2dd[2]/n_samP2-1/3)^2),180*n_samP2*((P2dd[3]/n_samP2-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[2])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),P2D))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1);
  B1w1<-1/(12*n_samB1)
  B1bmw <- matrix(0,n_samB1,1); B1bmwsl <- matrix(0,n_samB1,3)
  for(i in 1:3){
    B1gg <- (dataB1 - mean1[i])/sqrt(sigma1[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi1[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B1P2)))[1]
  if(nn < n_samB1){B1P2 <- B1P2+runif(n_samB1)/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*n_samB1) + sum((B1P2 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  B1u <- as.matrix(c(12*n_samB1*((B1dd[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((B1dd[2]/n_samB1-1/3)^2),180*n_samB1*((B1dd[3]/n_samB1-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[2])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),B1D))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])

  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2);
  B2w1<-1/(12*n_samB2)
  B2bmw <- matrix(0,n_samB2,1); B2bmwsl <- matrix(0,n_samB2,3)
  for(i in 1:3){
    B2gg <- (dataB2 - mean2[i])/sqrt(sigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B2P2)))[1]
  if(nn < n_samB2){B2P2 <- B2P2+runif(n_samB2)/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*n_samB2) + sum((B2P2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  B2u <- as.matrix(c(12*n_samB2*((B2dd[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((B2dd[2]/n_samB2-1/3)^2),180*n_samB2*((B2dd[3]/n_samB2-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[2])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),B2D))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])

  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2);
  F2w1<-1/(12*n_samF2)
  F2bmw <- matrix(0,n_samF2,1); F2bmwsl <- matrix(0,n_samF2,6)
  for(i in 1:6){
    F2gg <- (dataF2 - mean3[i])/sqrt(sigma3[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi3[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < n_samF2){F2P2 <- F2P2+runif(n_samF2)/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*n_samF2) + sum((F2P2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  F2u <- as.matrix(c(12*n_samF2*((F2dd[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((F2dd[2]/n_samF2-1/3)^2),180*n_samF2*((F2dd[3]/n_samF2-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[2])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),F2D))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])

  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)

  ###########first order parameters##################
  aa<-matrix(c(1,2,1,0,1,0,0,1,1,-2,-1,0,1,2,0.5,0.25,1,1,0.5,0.25,
               1,0,0.5,0.25,1,0,-0.5,0.25,1,-1,-0.5,0.25,1,-2,-0.5,0.25,
               1,2,0,0.25,1,1,0,0.25,1,0,0,0.25,1,-1,0,0.25,1,-2,0,0.25),14,4,byrow=T)
  b_line1<-as.matrix(c(mean[1],mean[2],mean[3],mean1[1],mean1[2],mean1[3],mean2[1],mean2[2],mean2[3],mean3[1],mean3[2],mean3[3],mean3[5],mean3[6]))
  B1<-solve(crossprod(aa,aa))%*%crossprod(aa,b_line1)
  ##########second order parameters#####################
  jj1<-sigmaB1-sigma1[1]
  if (jj1<0 || jj1>=sigmaB1) {jj1<-0}
  ll1<-jj1/sigmaB1
  mm1<-sigma1[1]-sigma
  if (mm1<0 || mm1>=sigmaB1) {mm1<-0}
  nn1<-mm1/sigmaB1
  jj2<-sigmaB2-sigma2[3]
  if (jj2<0 || jj2>=sigmaB2) {jj2<-0}
  ll2<-jj2/sigmaB2
  mm2<-sigma2[3]-sigma
  if (mm2<0 || mm2>=sigmaB2) {mm2<-0}
  nn2<-mm2/sigmaB2
  jj3<-sigmaF2-sigma3[1]
  if (jj3<0 || jj3>=sigmaF2) {jj3<-0}
  ll3<-jj3/sigmaF2
  mm3<-sigma3[1]-sigma
  if (mm3<0 || mm3>=sigmaF2) {mm3<-0}
  nn3<-mm3/sigmaF2

  output <- data.frame("MX2-EA-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4), round(sigma0,4),round(t(mean1),4)," ",round(t(sigma1),4)," ",
                       round(t(mix_pi1),4)," ",round(t(mean2),4)," ",round(t(sigma2),4)," ",round(t(mix_pi2),4)," ",
                       round(t(mean3),4)," "," "," ",round(t(sigma3),4)," "," "," ",round(t(mix_pi3),4)," "," "," ",
                       round(B1[1],4)," "," "," "," "," ",round(B1[2],4),round(B1[2],4)," "," "," "," "," "," ",round(B1[3],4),round(B1[4],4),
                       round(jj1,4),round(ll1*100,4),round(mm1,4),round(nn1*100,4),round(jj2,4),round(ll2*100,4),round(mm2,4),round(nn2*100,4),round(jj3,4),round(ll3*100,4),round(mm3,4),round(nn3*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi1,mi2,mi3)
  return(OUTPUT)
}

#####################MX2-CD-AD(E-5)##########################################
G6FModelFun[[23]] <- function(K1,logL,df11,df21,df31,df41,df51,df61,G6Ftext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1];n_samF1<-dim(dataF1)[1];n_samP2<-dim(dataP2)[1]
  n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1];n_samF2<-dim(dataF2)[1]
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2),sum(dataF2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2),sum(dataF2^2)))
  ss<-matrix(0,3,1);ss[1]<-s[1]-sumx[1]^2/n_samP1;ss[2]<-s[2]-sumx[2]^2/n_samF1;ss[3]<-s[3]-sumx[3]^2/n_samP2
  mean<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataB1),mean(dataB2),mean(dataF2)))
  sigma0<-sum(ss)/(n_samP1+n_samF1+n_samP2-3)
  sigmaB1<-var(dataB1);sigma40<-sigmaB1
  sigmaB2<-var(dataB2);sigma50<-sigmaB2
  sigmaF2<-var(dataF2);sigma60<-sigmaF2
  sigmaP1<-sigmaF1<-sigmaP2<-sigma0
  m_esp <- 0.0001;n_fam <- as.numeric(G6Ftext2)
  ###############procedure start###########################
  mi1<-matrix(0.25,4,1);sigma1<-matrix(0,4,1)
  mi2<-matrix(0.25,4,1);sigma2<-matrix(0,4,1)
  mi3<-as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  sigma3<-matrix(0,9,1);sigma<-sigma0
  a1<-sqrt(sigmaB1/n_samB1);if (mean[1]<mean[3]) a1<--a1
  mean1<-as.matrix(c(mean[1],mean[4]-0.5*a1,mean[4]-a1,mean[4]-1.5*a1))
  a2<-sqrt(sigma50/n_samB2);if (mean[1]<mean[3]) a2<--a2
  mean2<-as.matrix(c(mean[5]+a2,mean[5],mean[5]-a2,mean[5]-2*a2))
  a3<-sqrt(sigma60/n_samF2);if (mean[1]<mean[3]) a3<--a3
  mean3<-matrix(0,9,1)
  mean3[1]<-mean[1];mean3[2]<-mean[1]-0.5*a3
  mean3[3]<-(mean[1]+mean[3])/2;mean3[4]<-mean[6]+0.6*a3
  mean3[5]<-mean[6];mean3[6]<-mean3[4]-0.6*a3
  mean3[7]<-mean3[3];mean3[8]<-mean[6]-a3;mean3[9]<-mean[3]
  sigma1[1]<-sigmaB1/2;sigma2[4]<-sigmaB2/2;sigma3[1]<-sigmaF2/2
  gs<-matrix(0,2,1)
  gs[1]<-0.00453*mean[1]+0.00302*mean[2]+0.00453*mean[3]+0.03172*
    mean1[1]+0.03193*mean1[2]-0.02362*mean1[3]-0.02341*mean1[4]+
    0.08686*mean2[1]+0.08749*mean2[2]-0.07918*mean2[3]-0.07855*
    mean2[4]+0.08686*mean3[1]+0.08707*mean3[2]+0.0877*mean3[3]+
    0.03151*mean3[4]+0.03172*mean3[5]+0.03235*mean3[6]-0.13453*
    mean3[7]-0.13432*mean3[8]-0.13369*mean3[9];
  #da.
  gs[2]<-0.00453*mean[1]+0.00302*mean[2]+0.00453*mean[3]+0.03172*
    mean1[1]-0.02362*mean1[2]+0.03193*mean1[3]-0.02341*mean1[4]+
    0.08686*mean2[1]-0.07918*mean2[2]+0.08749*mean2[3]-0.07855*
    mean2[4]+0.08686*mean3[1]+0.03151*mean3[2]-0.13453*mean3[3]+
    0.08707*mean3[4]+0.03172*mean3[5]-0.13432*mean3[6]+0.0877*
    mean3[7]+0.03235*mean3[8]-0.13369*mean3[9];
  #db.
  g_aa1<-0.75*gs[2]^2/n_fam #   0.75*db*db.
  g_aa2<-0.75*gs[1]^2/n_fam #   0.75*da*da.
  g_aa3<-0.75*(gs[1]^2+gs[2]^2)/n_fam
  sigma1[2]<-sigma1[1]+g_aa1;sigma1[3]<-sigma1[1]+g_aa2
  sigma1[4]<-sigma1[1]+g_aa3;sigma2[1]<-sigma2[4]+g_aa3
  sigma2[2]<-sigma2[4]+g_aa2;sigma2[3]<-sigma2[4]+g_aa1
  sigma3[2]<-sigma3[1]+g_aa1;sigma3[3]<-sigma3[7]<-sigma3[9]<-sigma3[1]
  sigma3[4]<-sigma3[1]+g_aa2;sigma3[5]<-sigma3[1]+g_aa3
  sigma3[6]<-sigma3[4];sigma3[8]<-sigma3[2]
  L0<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma0))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma0))))+
    sum(log(dnorm(dataP2,mean[3],sqrt(sigma0))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mi1)))+
    sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mi3)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,4,n_samB1);      swx1 <- matrix(0,4,1)
  W2 <- matrix(0,4,n_samB2);      swx2 <- matrix(0,4,1)
  W3 <- matrix(0,9,n_samF2);      swx3 <- matrix(0,9,1)
  hh<-matrix(0,15,15);b_line<-matrix(0,15,1)
  n0<-matrix(0,18,1);s0<-matrix(0,18,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:4) { W1[i,] <- mi1[i]*dnorm(dataB1,mean1[i],sqrt(sigma1[i]))/dmixnorm(dataB1,mean1,sqrt(sigma1),mi1) }
    mix_pi1 <- as.matrix(rowSums(W1)/n_samB1)
    sumwx1 <- W1%*%dataB1
    for(i in 1:4) { W2[i,] <- mi2[i]*dnorm(dataB2,mean2[i],sqrt(sigma2[i]))/dmixnorm(dataB2,mean2,sqrt(sigma2),mi2) }
    mix_pi2 <- as.matrix(rowSums(W2)/n_samB2)
    sumwx2 <- W2%*%dataB2
    for(i in 1:9) { W3[i,] <- mi3[i]*dnorm(dataF2,mean3[i],sqrt(sigma3[i]))/dmixnorm(dataF2,mean3,sqrt(sigma3),mi3) }
    mix_pi3 <- as.matrix(rowSums(W3)/n_samF2)
    sumwx3 <- W3%*%dataF2
    #####################CM1-step for means##################
    n0[c(1:4)]<-mix_pi1[c(1:4)]*n_samB1;n0[c(5:8)]<-mix_pi2[c(1:4)]*n_samB2
    s0[c(1:9)]<-mix_pi3[c(1:9)]*n_samF2
    s0[c(1:9)][abs(s0[c(1:9)])<0.00000001]<-0.000001
    n0[c(1:8)][abs(n0[c(1:8)])<0.00000001]<-0.000001
    aaa0<-0;AA<-matrix(0,15,1);n_iter<-0;aaa1<-1000
    while(aaa1>0.001)
    {
      n_iter<-n_iter+1
      gs[1]<-0.00453*mean[1]+0.00302*mean[2]+0.00453*mean[3]+0.03172*
        mean1[1]+0.03193*mean1[2]-0.02362*mean1[3]-0.02341*mean1[4]+
        0.08686*mean2[1]+0.08749*mean2[2]-0.07918*mean2[3]-0.07855*
        mean2[4]+0.08686*mean3[1]+0.08707*mean3[2]+0.0877*mean3[3]+
        0.03151*mean3[4]+0.03172*mean3[5]+0.03235*mean3[6]-0.13453*
        mean3[7]-0.13432*mean3[8]-0.13369*mean3[9];
      #da.
      gs[2]<-0.00453*mean[1]+0.00302*mean[2]+0.00453*mean[3]+0.03172*
        mean1[1]-0.02362*mean1[2]+0.03193*mean1[3]-0.02341*mean1[4]+
        0.08686*mean2[1]-0.07918*mean2[2]+0.08749*mean2[3]-0.07855*
        mean2[4]+0.08686*mean3[1]+0.03151*mean3[2]-0.13453*mean3[3]+
        0.08707*mean3[4]+0.03172*mean3[5]-0.13432*mean3[6]+0.0877*
        mean3[7]+0.03235*mean3[8]-0.13369*mean3[9];
      #db.
      g_aa1<-0.75*gs[2]^2/n_fam #   0.75*db*db.
      g_aa2<-0.75*gs[1]^2/n_fam #   0.75*da*da.
      g_aa3<-0.75*(gs[1]^2+gs[2]^2)/n_fam
      sigma1[2]<-sigma1[1]+g_aa1;sigma1[3]<-sigma1[1]+g_aa2
      sigma1[4]<-sigma1[1]+g_aa3;sigma2[1]<-sigma2[4]+g_aa3
      sigma2[2]<-sigma2[4]+g_aa2;sigma2[3]<-sigma2[4]+g_aa1
      sigma3[2]<-sigma3[1]+g_aa1;sigma3[3]<-sigma3[7]<-sigma3[9]<-sigma3[1]
      sigma3[4]<-sigma3[1]+g_aa2;sigma3[5]<-sigma3[1]+g_aa3
      sigma3[6]<-sigma3[4];sigma3[8]<-sigma3[2]
      hh[1,1]<-sigma1[1]/n0[1]+sigma1[2]/n0[2]+sigma3[1]/s0[1]+sigma3[2]/s0[2]
      hh[1,2]<-0
      hh[1,3]<-sigma1[1]/n0[1]+sigma3[1]/s0[1]
      hh[1,4]<-0
      hh[1,5]<-0
      hh[1,6]<-0
      hh[1,7]<--sigma1[1]/n0[1]
      hh[1,8]<-0
      hh[1,9]<-0
      hh[1,10]<--sigma3[1]/s0[1]
      hh[1,11]<-0
      hh[1,12]<--sigma1[1]/n0[1]+sigma1[2]/n0[2]
      hh[1,13]<--sigma3[1]/s0[1]-sigma3[2]/s0[2]
      hh[1,14]<--5*sigma1[2]/n0[2]
      hh[1,15]<-4*sigma3[1]/s0[1]
      hh[2,2]<-sigma1[3]/n0[3]+sigma1[4]/n0[4]+sigma3[4]/s0[4]+sigma3[5]/s0[5]
      hh[2,3]<-sigma1[4]/n0[4]+sigma3[5]/s0[5]
      hh[2,4]<--sigma3[5]/s0[5]
      hh[2,5]<-0
      hh[2,6]<-hh[2,4]
      hh[2,7]<-sigma1[4]/n0[4]
      hh[2,8]<--sigma1[4]/n0[4]-2*sigma3[5]/s0[5]
      hh[2,9]<--4*sigma3[5]/s0[5]
      hh[2,10]<-0
      hh[2,11]<--2*sigma1[3]/n0[3]
      hh[2,12]<--sigma1[3]/n0[3]+sigma1[4]/n0[4]
      hh[2,13]<-sigma3[4]/s0[4]
      hh[2,14]<--5*sigma1[3]/n0[3]
      hh[2,15]<-2*sigma1[4]/n0[4]+4*sigma3[5]/s0[5]
      hh[3,3]<-sigma1[1]/n0[1]+sigma1[4]/n0[4]+sigma3[1]/s0[1]+sigma3[5]/s0[5]
      hh[3,4]<--sigma3[5]/s0[5]
      hh[3,5]<-0
      hh[3,6]<-hh[3,4]
      hh[3,7]<--sigma1[1]/n0[1]+sigma1[4]/n0[4]
      hh[3,8]<--sigma1[4]/n0[4]-2*sigma3[5]/s0[5]
      hh[3,9]<--4*sigma3[5]/s0[5]
      hh[3,10]<--sigma3[1]/s0[1]
      hh[3,11]<-0
      hh[3,12]<--sigma1[1]/n0[1]+sigma1[4]/n0[4]
      hh[3,13]<--sigma3[1]/s0[1]
      hh[3,14]<-0
      hh[3,15]<-2*sigma1[4]/n0[4]+4*sigma3[1]/s0[1]+4*sigma3[5]/s0[5]
      hh[4,4]<-sigma2[1]/n0[5]+sigma2[2]/n0[6]+sigma3[5]/s0[5]+sigma3[6]/s0[6]
      hh[4,5]<-0
      hh[4,6]<-sigma2[1]/n0[5]+sigma3[5]/s0[5]
      hh[4,7]<-sigma2[1]/n0[5]
      hh[4,8]<-sigma2[1]/n0[5]+2*sigma3[5]/s0[5]
      hh[4,9]<-4*sigma3[5]/s0[5]
      hh[4,10]<-0
      hh[4,11]<--2*sigma2[2]/n0[6]
      hh[4,12]<--sigma2[1]/n0[5]+sigma2[2]/n0[6]
      hh[4,13]<-sigma3[6]/s0[6]
      hh[4,14]<-3*sigma2[2]/n0[6]
      hh[4,15]<-2*sigma2[1]/n0[5]-4*sigma3[5]/s0[5]
      hh[5,5]<-sigma2[3]/n0[7]+sigma2[4]/n0[8]+sigma3[8]/s0[8]+sigma3[9]/s0[9]
      hh[5,6]<-sigma2[4]/n0[8]+sigma3[9]/s0[9]
      hh[5,7]<--sigma2[4]/n0[8]
      hh[5,8]<-0
      hh[5,9]<-0
      hh[5,10]<-sigma3[9]/s0[9]
      hh[5,11]<-0
      hh[5,12]<--sigma2[3]/n0[7]+sigma2[4]/n0[8]
      hh[5,13]<--sigma3[8]/s0[8]-sigma3[9]/s0[9]
      hh[5,14]<-3*sigma2[3]/n0[7]
      hh[5,15]<-0
      hh[6,6]<-sigma2[1]/n0[5]+sigma2[4]/n0[8]+sigma3[5]/s0[5]+sigma3[9]/s0[9]
      hh[6,7]<-sigma2[1]/n0[5]-sigma2[4]/n0[8]
      hh[6,8]<-sigma2[1]/n0[5]+2*sigma3[5]/s0[5]
      hh[6,9]<-4*sigma3[5]/s0[5]
      hh[6,10]<-sigma3[9]/s0[9]
      hh[6,11]<-0
      hh[6,12]<--sigma2[1]/n0[5]+sigma2[4]/n0[8]
      hh[6,13]<--sigma3[9]/s0[9]
      hh[6,14]<-0
      hh[6,15]<-2*sigma2[1]/n0[5]-4*sigma3[5]/s0[5]
      hh[7,7]<-sigma*(1/n_samP1+1/n_samP2)+sigma1[1]/n0[1]+sigma1[4]/n0[4]+sigma2[1]/n0[5]+sigma2[4]/n0[8]
      hh[7,8]<--sigma1[4]/n0[4]+sigma2[1]/n0[5]
      hh[7,9]<-sigma*(3/n_samP1-3/n_samP2)
      hh[7,10]<-0
      hh[7,11]<-sigma*(1/n_samP1+1/n_samP2)
      hh[7,12]<-sigma*(3/n_samP1-3/n_samP2)+sigma1[1]/n0[1]+sigma1[4]/n0[4]-sigma2[1]/n0[5]-sigma2[4]/n0[8]
      hh[7,13]<-0
      hh[7,14]<-0
      hh[7,15]<-sigma*(1/n_samP1+1/n_samP2)+2*sigma1[4]/n0[4]+2*sigma2[1]/n0[5]
      hh[8,8]<-sigma1[4]/n0[4]+sigma2[1]/n0[5]+4*sigma3[5]/s0[5]
      hh[8,9]<-8*sigma3[5]/s0[5]
      hh[8,10]<-0
      hh[8,11]<-0
      hh[8,12]<--sigma1[4]/n0[4]-sigma2[1]/n0[5]
      hh[8,13]<-0
      hh[8,14]<-0
      hh[8,15]<--2*sigma1[4]/n0[4]+2*sigma2[1]/n0[5]-8*sigma3[5]/s0[5]
      hh[9,9]<-sigma*(9/n_samP1+4/n_samF1+9/n_samP2)+4*sigma3[3]/s0[3]+16*sigma3[5]/s0[5]+4*sigma3[7]/s0[7]
      hh[9,10]<-2*sigma3[3]/s0[3]+2*sigma3[7]/s0[7]
      hh[9,11]<-sigma*(3/n_samP1-3/n_samP2)+2*sigma3[3]/s0[3]-2*sigma3[7]/s0[7]
      hh[9,12]<-sigma*(9/n_samP1+4/n_samF1+9/n_samP2)
      hh[9,13]<-0
      hh[9,14]<--2*sigma3[3]/s0[3]+2*sigma3[7]/s0[7]
      hh[9,15]<-sigma*(3/n_samP1-3/n_samP2)-16*sigma3[5]/s0[5]
      hh[10,10]<-sigma3[1]/s0[1]+sigma3[3]/s0[3]+sigma3[7]/s0[7]+sigma3[9]/s0[9]
      hh[10,11]<-sigma3[3]/s0[3]-sigma3[7]/s0[7]
      hh[10,12]<-0
      hh[10,13]<-sigma3[1]/s0[1]-sigma3[9]/s0[9]
      hh[10,14]<--sigma3[3]/s0[3]+sigma3[7]/s0[7]
      hh[10,15]<--4*sigma3[1]/s0[1]
      hh[11,11]<-sigma*(1/n_samP1+1/n_samP2)+4*sigma1[3]/n0[3]+4*sigma2[2]/n0[6]+sigma3[3]/s0[3]+sigma3[7]/s0[7]
      hh[11,12]<-sigma*(3/n_samP1-3/n_samP2)+2*sigma1[3]/n0[3]-2*sigma2[2]/n0[6]
      hh[11,13]<-0
      hh[11,14]<-10*sigma1[3]/n0[3]-6*sigma2[2]/n0[6]-sigma3[3]/s0[3]-sigma3[7]/s0[7]
      hh[11,15]<-sigma*(1/n_samP1+1/n_samP2)
      hh[12,12]<-sigma*(9/n_samP1+4/n_samF1+9/n_samP2)+sigma1[1]/n0[1]+sigma1[2]/n0[2]+sigma1[3]/n0[3]+sigma1[4]/n0[4]+sigma2[1]/n0[5]+sigma2[2]/n0[6]+sigma2[3]/n0[7]+sigma2[4]/n0[8]
      hh[12,13]<-0
      hh[12,14]<--5*sigma1[2]/n0[2]+5*sigma1[3]/n0[3]+3*sigma2[2]/n0[6]-3*sigma2[3]/n0[7]
      hh[12,15]<-sigma*(3/n_samP1-3/n_samP2)+2*sigma1[4]/n0[4]-2*sigma2[1]/n0[5]
      hh[13,13]<-sigma3[1]/s0[1]+sigma3[2]/s0[2]+sigma3[4]/s0[4]+sigma3[6]/s0[6]+sigma3[8]/s0[8]+sigma3[9]/s0[9]
      hh[13,14]<-0
      hh[13,15]<--4*sigma3[1]/s0[1]
      hh[14,14]<-25*sigma1[2]/n0[2]+25*sigma1[3]/n0[3]+9*sigma2[2]/n0[6]+9*sigma2[3]/n0[7]+sigma3[3]/s0[3]+sigma3[7]/s0[7]
      hh[14,15]<-0
      hh[15,15]<-sigma*(1/n_samP1+1/n_samP2)+4*sigma1[4]/n0[4]+4*sigma2[1]/n0[5]+ 16*sigma3[1]/s0[1]+16*sigma3[5]/s0[5]
      for(i in 2:15)
      {
        for(j in 1:(i-1))
        {
          hh[i,j]<-hh[j,i]
        }
      }
      ###########################################################################
      b_line[1]<-sumwx1[1]/n0[1]-sumwx1[2]/n0[2]-sumwx3[1]/s0[1]+sumwx3[2]/s0[2]
      b_line[2]<-sumwx1[3]/n0[3]-sumwx1[4]/n0[4]-sumwx3[4]/s0[4]+sumwx3[5]/s0[5]
      b_line[3]<-sumwx1[1]/n0[1]-sumwx1[4]/n0[4]-sumwx3[1]/s0[1]+sumwx3[5]/s0[5]
      b_line[4]<-sumwx2[1]/n0[5]-sumwx2[2]/n0[6]-sumwx3[5]/s0[5]+sumwx3[6]/s0[6]
      b_line[5]<-sumwx2[3]/n0[7]-sumwx2[4]/n0[8]-sumwx3[8]/s0[8]+sumwx3[9]/s0[9]
      b_line[6]<-sumwx2[1]/n0[5]-sumwx2[4]/n0[8]-sumwx3[5]/s0[5]+sumwx3[9]/s0[9]
      b_line[7]<-sumx[1]/n_samP1-sumx[3]/n_samP2-sumwx1[1]/n0[1]-sumwx1[4]/n0[4]+sumwx2[1]/n0[5]+sumwx2[4]/n0[8]
      b_line[8]<-sumwx1[4]/n0[4]+sumwx2[1]/n0[5]-2*sumwx3[5]/s0[5]
      b_line[9]<-3*sumx[1]/n_samP1+2*sumx[2]/n_samF1+3*sumx[3]/n_samP2-2*sumwx3[3]/s0[3]-4*sumwx3[5]/s0[5]-2*sumwx3[7]/s0[7]
      b_line[10]<-sumwx3[1]/s0[1]-sumwx3[3]/s0[3]-sumwx3[7]/s0[7]+sumwx3[9]/s0[9]
      b_line[11]<-sumx[1]/n_samP1-sumx[3]/n_samP2-2*sumwx1[3]/n0[3]+2*sumwx2[2]/n0[6]-sumwx3[3]/s0[3]+sumwx3[7]/s0[7]
      b_line[12]<-3*sumx[1]/n_samP1+2*sumx[2]/n_samF1+3*sumx[3]/n_samP2-sumwx1[1]/n0[1]-sumwx1[2]/n0[2]-sumwx1[3]/n0[3]-sumwx1[4]/n0[4]-sumwx2[1]/n0[5]-sumwx2[2]/n0[6]-sumwx2[3]/n0[7]-sumwx2[4]/n0[8]
      b_line[13]<-sumwx3[1]/s0[1]-sumwx3[2]/s0[2]-sumwx3[4]/s0[4]+sumwx3[6]/s0[6]+sumwx3[8]/s0[8]-sumwx3[9]/s0[9]
      b_line[14]<-5*sumwx1[2]/n0[2]-5*sumwx1[3]/n0[3]-3*sumwx2[2]/n0[6]+3*sumwx2[3]/n0[7]+sumwx3[3]/s0[3]-sumwx3[7]/s0[7]
      b_line[15]<-sumx[1]/n_samP1-sumx[3]/n_samP2-2*sumwx1[4]/n0[4]+2*sumwx2[1]/n0[5]-4*sumwx3[1]/s0[1]+4*sumwx3[5]/s0[5]
      B<-solve(hh,b_line)
      mean[1]<-(sumx[1]+sigma*(-B[7]-3*B[9]-B[11]-3*B[12]-B[15]))/n_samP1
      mean[2]<-(sumx[2]+sigma*(-2*B[9]-2*B[12]))/n_samF1
      mean[3]<-(sumx[3]+sigma*(B[7]-3*B[9]+B[11]-3*B[12]+B[15]))/n_samP2
      mean1[1]<-(sumwx1[1]+sigma1[1]*(-B[1]-B[3]+B[7]+B[12]))/n0[1]
      mean1[2]<-(sumwx1[2]+sigma1[2]*(B[1]+B[12]-5*B[14]))/n0[2]
      mean1[3]<-(sumwx1[3]+sigma1[3]*(-B[2]+2*B[11]+B[12]+5*B[14]))/n0[3]
      mean1[4]<-(sumwx1[4]+sigma1[4]*(B[2]+B[3]+B[7]-B[8]+B[12]+2*B[15]))/n0[4]
      mean2[1]<-(sumwx2[1]+sigma2[1]*(-B[4]-B[6]-B[7]-B[8]+B[12]-2*B[15]))/n0[5]
      mean2[2]<-(sumwx2[2]+sigma2[2]*(B[4]-2*B[11]+B[12]+3*B[14]))/n0[6]
      mean2[3]<-(sumwx2[3]+sigma2[3]*(-B[5]+B[12]-3*B[14]))/n0[7]
      mean2[4]<-(sumwx2[4]+sigma2[4]*(B[5]+B[6]-B[7]+B[12]))/n0[8]
      mean3[1]<-(sumwx3[1]+sigma3[1]*(B[1]+B[3]-B[10]-B[13]+4*B[15]))/s0[1]
      mean3[2]<-(sumwx3[2]+sigma3[2]*(-B[1]+B[13]))/s0[2]
      mean3[3]<-(sumwx3[3]+sigma3[3]*(2*B[9]+B[10]+B[11]-B[14]))/s0[3]
      mean3[7]<-(sumwx3[7]+sigma3[7]*(2*B[9]+B[10]-B[11]+B[14]))/s0[7]
      mean3[4]<-(sumwx3[4]+sigma3[4]*(B[2]+B[13]))/s0[4]
      mean3[5]<-(sumwx3[5]+sigma3[5]*(-B[2]-B[3]+B[4]+B[6]+2*B[8]+4*B[9]-4*B[15]))/s0[5]
      mean3[6]<-(sumwx3[6]+sigma3[6]*(-B[4]-B[13]))/s0[6]
      mean3[8]<-(sumwx3[8]+sigma3[8]*(B[5]-B[13]))/s0[8]
      mean3[9]<-(sumwx3[9]+sigma3[9]*(-B[5]-B[6]-B[10]+B[13]))/s0[9]
      aaa1<-max(abs(B-AA))
      AA<-B
      if (n_iter>20) break
    }

    ################obtain variance############################
    ss1<-sum((dataP1-mean[1])^2);ss3<-sum((dataP2-mean[3])^2);ss2<-sum((dataF1-mean[2])^2)
    for(i in 1:4) {swx1[i] <- W1[i,]%*%(dataB1-mean1[i])^2 };for(i in 1:4) {swx2[i] <- W2[i,]%*%(dataB2-mean2[i])^2 }; for(i in 1:9) {swx3[i] <- W3[i,]%*%(dataF2-mean3[i])^2 }
    aaa0<-sigma1[1];n_iter<-0;aaa1<-1000
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      aa2<-sigma1[1]/(sigma1[1]+g_aa1)
      aa3<-sigma1[1]/(sigma1[1]+g_aa2)
      aa4<-sigma1[1]/(sigma1[1]+g_aa3)
      as1<-swx1[1]+swx1[2]*aa2^2+swx1[3]*aa3^2+swx1[4]*aa4^2
      as2<-n0[1]+aa2*n0[2]+aa3*n0[3]+aa4*n0[4]
      sigma1[1]<-as1/as2
      aaa1<-abs(sigma1[1]-aaa0)
      aaa0<-sigma1[1]
      if (n_iter>20) break
    }
    sigma40<-sigma1[1]-sigma;
    if (sigma40<0) {sigma40<-0;sigma1[1]<-sigma}
    sigma1[1]<-sigma40+sigma;sigma1[2]<-sigma1[1]+g_aa1
    sigma1[3]<-sigma1[1]+g_aa2;sigma1[4]<-sigma1[1]+g_aa3
    # to estimate sigma50.
    aaa0<-sigma2[4];n_iter<-0;aaa1<-1000
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma2[4]/(sigma2[4]+g_aa3)
      aa2<-sigma2[4]/(sigma2[4]+g_aa2)
      aa3<-sigma2[4]/(sigma2[4]+g_aa1)
      as3<-swx2[1]*aa1^2+swx2[2]*aa2^2+swx2[3]*aa3^2+swx2[4]
      as4<-aa1*n0[5]+aa2*n0[6]+aa3*n0[7]+n0[8]
      sigma2[4]<-as3/as4
      aaa1<-abs(sigma2[4]-aaa0)
      aaa0<-sigma2[4]
      if (n_iter>20) break
    }
    sigma50<-sigma2[4]-sigma;
    if (sigma50<0) {sigma50<-0;sigma2[4]<-sigma}
    sigma2[4]<-sigma50+sigma;sigma2[1]<-sigma2[4]+g_aa3
    sigma2[2]<-sigma2[4]+g_aa2;sigma2[3]<-sigma2[4]+g_aa1
    # to estimate sigma60.
    aaa0<-sigma3[1]
    aa6<-swx3[1]+swx3[3]+swx3[7]+swx3[9];aa7<-s0[1]+s0[3]+s0[7]+s0[9]
    n_iter<-0
    aaa1<-1000
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma3[1]/(sigma3[1]+g_aa1)
      aa2<-sigma3[1]/(sigma3[1]+g_aa2)
      aa3<-sigma3[1]/(sigma3[1]+g_aa3)
      as5<-aa6+(swx3[2]+swx3[8])*aa1^2+(swx3[4]+swx3[6])*aa2^2+swx3[5]*aa3^2
      as6<-aa7+aa1*(s0[2]+s0[8])+aa2*(s0[4]+s0[6])+aa3*s0[5]
      sigma3[1]<-as5/as6
      aaa1<-abs(sigma3[1]-aaa0)
      aaa0<-sigma3[1]
      if (n_iter>20) break
    }
    sigma60<-sigma3[1]-sigma;
    if (sigma60<0) {sigma60<-0;sigma3[1]<-sigma}
    sigma3[1]<-sigma60+sigma;sigma3[2]<-sigma3[1]+g_aa1
    sigma3[3]<-sigma3[7]<-sigma3[9]<-sigma3[1]
    sigma3[4]<-sigma3[1]+g_aa2;sigma3[5]<-sigma3[1]+g_aa3
    sigma3[6]<-sigma3[4];sigma3[8]<-sigma3[2]
    # CM3 to estimate the variance (sigma).
    ab1<-ss1+ss2+ss3;ab2<-n_samP1+n_samF1+n_samP2
    n_iter<-0;aaa0<-sigma;aaa1<-1000
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      n0[11]<-sigma/(sigma+sigma40)
      n0[12]<-sigma/(sigma+sigma40+g_aa1)
      n0[13]<-sigma/(sigma+sigma40+g_aa2)
      n0[14]<-sigma/(sigma+sigma40+g_aa3)
      s0[11]<-sigma/(sigma+sigma50+g_aa3)
      s0[12]<-sigma/(sigma+sigma50+g_aa2)
      s0[13]<-sigma/(sigma+sigma50+g_aa1)
      s0[14]<-sigma/(sigma+sigma50)
      ab3<-sum(swx1[c(1:4)]*n0[c(11:14)]^2+swx2[c(1:4)]*s0[c(11:14)]^2)
      ab4<-sum(n0[c(1:4)]*n0[c(11:14)]+n0[c(5:8)]*s0[c(11:14)])
      n0[11]<-sigma/(sigma+sigma60);n0[13]<-n0[17]<-n0[19]<-n0[11]
      n0[12]<-sigma/(sigma+sigma60+g_aa1);n0[18]<-n0[12]
      n0[14]<-sigma/(sigma+sigma60+g_aa2);n0[16]<-n0[14]
      n0[15]<-sigma/(sigma+sigma60+g_aa3)
      ab3<-ab3+sum(swx3[c(1:9)]*n0[c(11:19)]^2)
      ab4<-ab4+sum(s0[c(1:9)]*n0[11:19])
      sigma<-(ab1+ab3)/(ab2+ab4);aaa1<-abs(sigma-aaa0);aaa0<-sigma
      if (n_iter>20) break
    }

    sigma1[1]<-sigma+sigma40;sigma1[2]<-sigma1[1]+g_aa1
    sigma1[3]<-sigma1[1]+g_aa2;sigma1[4]<-sigma1[1]+g_aa3
    sigma2[4]<-sigma+sigma50;sigma2[1]<-sigma2[4]+g_aa3
    sigma2[2]<-sigma2[4]+g_aa2;sigma2[3]<-sigma2[4]+g_aa1
    sigma3[1]<-sigma+sigma60;sigma3[3]<-sigma3[7]<-sigma3[9]<-sigma3[1]
    sigma3[2]<-sigma3[1]+g_aa1;sigma3[8]<-sigma3[2]
    sigma3[4]<-sigma3[1]+g_aa2;sigma3[6]<-sigma3[4]
    sigma3[5]<-sigma3[1]+g_aa3
    if(sum(sigma < 1e-30)>=1){break}
    ##############################criteria for iterations to stop########################
    L1<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,mean[3],sqrt(sigma))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mix_pi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mix_pi3)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc<-L0
  AIC<--2*abc+2*9
  meanP1<-mean[1];meanF1<-mean[2];meanP2<-mean[3]
  sigma0<-sigma

  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*n_samP1)
  P1bmw <- matrix(0,n_samP1,1)
  P1gg <- (dataP1 - mean[1])/sqrt(as.vector(sigma))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < n_samP1){P1bmw <- P1bmw+runif(n_samP1)/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  P1u<- as.matrix(c(12*n_samP1*((P1dd[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((P1dd[2]/n_samP1-1/3)^2),180*n_samP1*((P1dd[3]/n_samP1-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[2])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),P1D))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*n_samF1)
  F1bmw <- matrix(0,n_samF1,1)
  F1gg <- (dataF1 - mean[2])/sqrt(as.vector(sigma))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < n_samF1){F1bmw <- F1bmw+runif(n_samF1)/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  F1u<- as.matrix(c(12*n_samF1*((F1dd[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((F1dd[2]/n_samF1-1/3)^2),180*n_samF1*((F1dd[3]/n_samF1-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[2])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),F1D))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*n_samP2)
  P2bmw <- matrix(0,n_samP2,1)
  P2gg <- (dataP2 - mean[3])/sqrt(as.vector(sigma))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < n_samP2){P2bmw <- P2bmw+runif(n_samP2)/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  P2u<- as.matrix(c(12*n_samP2*((P2dd[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((P2dd[2]/n_samP2-1/3)^2),180*n_samP2*((P2dd[3]/n_samP2-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[2])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),P2D))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1);
  B1w1<-1/(12*n_samB1)
  B1bmw <- matrix(0,n_samB1,1); B1bmwsl <- matrix(0,n_samB1,4)
  for(i in 1:4){
    B1gg <- (dataB1 - mean1[i])/sqrt(sigma1[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi1[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B1P2)))[1]
  if(nn < n_samB1){B1P2 <- B1P2+runif(n_samB1)/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*n_samB1) + sum((B1P2 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  B1u <- as.matrix(c(12*n_samB1*((B1dd[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((B1dd[2]/n_samB1-1/3)^2),180*n_samB1*((B1dd[3]/n_samB1-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[2])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),B1D))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])

  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2);
  B2w1<-1/(12*n_samB2)
  B2bmw <- matrix(0,n_samB2,1); B2bmwsl <- matrix(0,n_samB2,4)
  for(i in 1:4){
    B2gg <- (dataB2 - mean2[i])/sqrt(sigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B2P2)))[1]
  if(nn < n_samB2){B2P2 <- B2P2+runif(n_samB2)/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*n_samB2) + sum((B2P2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  B2u <- as.matrix(c(12*n_samB2*((B2dd[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((B2dd[2]/n_samB2-1/3)^2),180*n_samB2*((B2dd[3]/n_samB2-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[2])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),B2D))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])

  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2);
  F2w1<-1/(12*n_samF2)
  F2bmw <- matrix(0,n_samF2,1); F2bmwsl <- matrix(0,n_samF2,9)
  for(i in 1:9){
    F2gg <- (dataF2 - mean3[i])/sqrt(sigma3[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi3[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < n_samF2){F2P2 <- F2P2+runif(n_samF2)/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*n_samF2) + sum((F2P2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  F2u <- as.matrix(c(12*n_samF2*((F2dd[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((F2dd[2]/n_samF2-1/3)^2),180*n_samF2*((F2dd[3]/n_samF2-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[2])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),F2D))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])

  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)


  ################first order parameters###############
  aa<-matrix(c(1,1,1,1,0,1,1,1,0,1,1,-1,-1,-1,0,1,1,1,0.5,0.25,1,1,0.5,0.5,0.25,1,0.5,1,0.5,
               0.25,1,0.5,0.5,0.5,0.25,1,0.5,0.5,-0.5,0.25,1,0.5,-1,-0.5,0.25,1,-1,0.5,-0.5,
               0.25,1,-1,-1,-0.5,0.25,1,1,1,0,0.25,1,1,0.5,0,0.25,1,1,-1,0,0.25,1,0.5,1,0,0.25,
               1,0.5,0.5,0,0.25,1,0.5,-1,0,0.25,1,-1,1,0,0.25,1,-1,0.5,0,0.25,1,-1,-1,0,0.25),20,5,byrow=T)
  b_line1<-as.matrix(c(mean[1],mean[2],mean[3],mean1[1],mean1[2],mean1[3],mean1[4],mean2[1],mean2[2],mean2[3],mean2[4],
                       mean3[1],mean3[2],mean3[3],mean3[4],mean3[5],mean3[6],mean3[7],mean3[8],mean3[9]))
  B1<-solve(crossprod(aa,aa))%*%crossprod(aa,b_line1)
  #########second order parameters##################
  jj1<-sigmaB1-sigma1[1]
  if (jj1<0 || jj1>=sigmaB1) {jj1<-0}
  ll1<-jj1/sigmaB1
  mm1<-sigma1[1]-sigma
  if (mm1<0 || mm1>=sigmaB1) {mm1<-0}
  nn1<-mm1/sigmaB1
  jj2<-sigmaB2-sigma2[4]
  if (jj2<0 || jj2>=sigmaB2) {jj2<-0}
  ll2<-jj2/sigmaB2
  mm2<-sigma2[4]-sigma
  if (mm2<0 || mm2>=sigmaB2) {mm2<-0}
  nn2<-mm2/sigmaB2
  jj3<-sigmaF2-sigma3[1]
  if (jj3<0 || jj3>=sigmaF2) {jj3<-0}
  ll3<-jj3/sigmaF2
  mm3<-sigma3[1]-sigma
  if (mm3<0 || mm3>=sigmaF2) {mm3<-0}
  nn3<-mm3/sigmaF2

  output <- data.frame("MX2-CD-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4), round(sigma0,4),round(t(mean1),4),round(t(sigma1),4),
                       round(t(mix_pi1),4),round(t(mean2),4),round(t(sigma2),4),round(t(mix_pi2),4),
                       round(t(mean3),4),round(t(sigma3),4),round(t(mix_pi3),4),
                       round(B1[1],4)," "," "," "," "," ",round(B1[2],4),round(B1[3],4),round(B1[2],4),round(B1[3],4)," "," "," "," ",round(B1[4],4),round(B1[5],4),
                       round(jj1,4),round(ll1*100,4),round(mm1,4),round(nn1*100,4),round(jj2,4),round(ll2*100,4),round(mm2,4),round(nn2*100,4),round(jj3,4),round(ll3*100,4),round(mm3,4),round(nn3*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi1,mi2,mi3)
  return(OUTPUT)
}

#################MX2-EAD-AD(E-6)#########################
G6FModelFun[[24]] <- function(K1,logL,df11,df21,df31,df41,df51,df61,G6Ftext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]));dataF2 <- as.matrix(as.numeric(df61[,1]))
  n_samP1<-dim(dataP1)[1];n_samF1<-dim(dataF1)[1];n_samP2<-dim(dataP2)[1]
  n_samB1<-dim(dataB1)[1];n_samB2<-dim(dataB2)[1];n_samF2<-dim(dataF2)[1]
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2),sum(dataF2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2),sum(dataF2^2)))
  ss<-matrix(0,3,1);ss[1]<-s[1]-sumx[1]^2/n_samP1;ss[2]<-s[2]-sumx[2]^2/n_samF1;ss[3]<-s[3]-sumx[3]^2/n_samP2
  mean<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataB1),mean(dataB2),mean(dataF2)))
  sigma0<-sum(ss)/(n_samP1+n_samF1+n_samP2-3)
  sigmaB1<-var(dataB1);sigma40<-sigmaB1
  sigmaB2<-var(dataB2);sigma50<-sigmaB2
  sigmaF2<-var(dataF2);sigma60<-sigmaF2
  sigmaP1<-sigmaF1<-sigmaP2<-sigma0
  m_esp <- 0.0001;n_fam <- as.numeric(G6Ftext2)
  ###############procedure start###########################
  mi1<-as.matrix(c(0.25,0.5,0.25));sigma1<-matrix(0,3,1)
  mi2<-as.matrix(c(0.25,0.5,0.25));sigma2<-matrix(0,3,1)
  mi3<-as.matrix(c(0.0625,0.25,0.125,0.25,0.25,0.0625))
  sigma3<-matrix(0,6,1);sigma<-sigma0
  a1<-sqrt(sigmaB1/n_samB1);if (mean[1]<mean[3]) {a1<--a1}
  mean1<-as.matrix(c(mean[1],mean[4],mean[4]-a1))
  a2<-sqrt(sigmaB2/n_samB2);if (mean[1]<mean[3]) {a2<--a2}
  mean2<-as.matrix(c(mean1[3],mean[5],mean[5]-a2))
  a3<-sqrt(sigmaF2/n_samF2);if (mean[1]<mean[3]){a3<--a3}
  mean3<-as.matrix(c(mean[6]+2.5*a3,mean[6]+1.5*a3,mean[6]+0.2*a3,mean[6]+a3,mean[6]-0.5*a3,mean[6]-2.5*a3))
  sigma1[1]<-sigmaB1/2;sigma2[3]<-sigmaB2/2;sigma3[1]<-sigmaF2/2
  gs<-matrix(0,1,1)
  gs[1]<-0.00459*mean[1]+0.00306*mean[2]+0.00459*mean[3]+0.03559*mean1[1]+
    0.00421*mean1[2]-0.02717*mean1[3]+0.09835*mean2[1]+0.00421*mean2[2]-
    0.08993*mean2[3]+0.09835*mean3[1]+0.06697*mean3[2]-0.02717*mean3[3]+
    0.03559*mean3[4]-0.05855*mean3[5]-0.15270*mean3[6]
  #d.
  g_aa1<-0.75*gs[1]^2/n_fam  #   0.75*d*d.
  g_aa2<-1.5*gs[1]^2/n_fam   #   1.5*d*d.
  sigma1[2]<-sigma1[1]+g_aa1;sigma1[3]<-sigma1[1]+g_aa2
  sigma2[1]<-sigma2[3]+g_aa2;sigma2[2]<-sigma2[3]+g_aa1
  sigma3[2]<-sigma3[1]+g_aa1;sigma3[3]<-sigma3[6]<-sigma3[1]
  sigma3[4]<-sigma3[1]+g_aa2;sigma3[5]<-sigma3[1]+g_aa1
  L0<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma0))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma0))))+
    sum(log(dnorm(dataP2,mean[3],sqrt(sigma0))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mi1)))+
    sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mi3)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,3,n_samB1);      swx1 <- matrix(0,3,1)
  W2 <- matrix(0,3,n_samB2);      swx2 <- matrix(0,3,1)
  W3 <- matrix(0,6,n_samF2);      swx3 <- matrix(0,6,1)
  hh<-matrix(0,11,11);b_line<-matrix(0,11,1)
  n0<-matrix(0,15,1);s0<-matrix(0,15,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:3) { W1[i,] <- mi1[i]*dnorm(dataB1,mean1[i],sqrt(sigma1[i]))/dmixnorm(dataB1,mean1,sqrt(sigma1),mi1) }
    mix_pi1 <- as.matrix(rowSums(W1)/n_samB1)
    sumwx1 <- W1%*%dataB1
    for(i in 1:3) { W2[i,] <- mi2[i]*dnorm(dataB2,mean2[i],sqrt(sigma2[i]))/dmixnorm(dataB2,mean2,sqrt(sigma2),mi2) }
    mix_pi2 <- as.matrix(rowSums(W2)/n_samB2)
    sumwx2 <- W2%*%dataB2
    for(i in 1:6) { W3[i,] <- mi3[i]*dnorm(dataF2,mean3[i],sqrt(sigma3[i]))/dmixnorm(dataF2,mean3,sqrt(sigma3),mi3) }
    mix_pi3 <- as.matrix(rowSums(W3)/n_samF2)
    sumwx3 <- W3%*%dataF2
    ######obtain means#################################
    n0[c(1:3)]<-mix_pi1[c(1:3)]*n_samB1;n0[c(4:6)]<-mix_pi2[c(1:3)]*n_samB2
    s0[c(1:6)]<-mix_pi3[c(1:6)]*n_samF2
    s0[c(1:6)][abs(s0[c(1:6)])<0.00000001]<-0.000001
    n0[c(1:6)][abs(n0[c(1:6)])<0.00000001]<-0.000001
    aaa0<-0;AA<-matrix(0,11,1);n_iter<-0;aaa1<-1000
    while(aaa1>0.001)
    {
      n_iter<-n_iter+1;
      gs[1]<-0.00459*mean[1]+0.00306*mean[2]+0.00459*mean[3]+0.03559*mean1[1]+
        0.00421*mean1[2]-0.02717*mean1[3]+0.09835*mean2[1]+0.00421*mean2[2]-
        0.08993*mean2[3]+0.09835*mean3[1]+0.06697*mean3[2]-0.02717*mean3[3]+
        0.03559*mean3[4]-0.05855*mean3[5]-0.15270*mean3[6]
      #d.
      g_aa1<-0.75*gs[1]*gs[1]/n_fam  #   0.75*d*d.
      g_aa2<-1.5*gs[1]*gs[1]/n_fam   #   1.5*d*d.
      sigma1[2]<-sigma1[1]+g_aa1;sigma1[3]<-sigma1[1]+g_aa2
      sigma2[1]<-sigma2[3]+g_aa2;sigma2[2]<-sigma2[3]+g_aa1
      sigma3[2]<-sigma3[1]+g_aa1;sigma3[3]<-sigma3[6]<-sigma3[1]
      sigma3[4]<-sigma3[1]+g_aa2;sigma3[5]<-sigma3[1]+g_aa1
      hh[1,1]<-sigma*(1/n_samP1+1/n_samF1+4/n_samP2)+sigma1[3]/n0[3]+sigma2[1]/n0[4]+4*sigma3[1]/s0[1]+4*sigma3[6]/s0[6]
      hh[1,2]<-sigma*(3/n_samP1+2/n_samF1+6/n_samP2)
      hh[1,3]<-sigma*(1/n_samP1-2/n_samP2)+4*sigma3[6]/s0[6]
      hh[1,4]<-sigma*(1/n_samP1-2/n_samP2)-sigma1[3]/n0[3]-sigma2[1]/n0[4]
      hh[1,5]<-sigma1[3]/n0[3]
      hh[1,6]<--6*sigma3[1]/s0[1]
      hh[1,7]<-2*sigma3[1]/s0[1]+2*sigma3[6]/s0[6]
      hh[1,8]<-0
      hh[1,9]<-sigma1[3]/n0[3]-sigma2[1]/n0[4]
      hh[1,10]<-0
      hh[1,11]<--2*sigma3[6]/s0[6]
      hh[2,2]<-sigma*(9/n_samP1+4/n_samF1+9/n_samP2)+4*sigma1[1]/n0[1]+4*sigma2[3]/n0[6]+16*sigma3[4]/s0[4]
      hh[2,3]<-sigma*(3/n_samP1-3/n_samP2)-4*sigma1[1]/n0[1]
      hh[2,4]<-sigma*(3/n_samP1-3/n_samP2)+2*sigma1[1]/n0[1]-2*sigma2[3]/n0[6]
      hh[2,5]<--2*sigma1[1]/n0[1]
      hh[2,6]<-4*sigma3[4]/s0[4]
      hh[2,7]<--2*sigma1[1]/n0[1]-2*sigma2[3]/n0[6]
      hh[2,8]<-0
      hh[2,9]<-8*sigma3[4]/s0[4]
      hh[2,10]<-4*sigma3[4]/s0[4]
      hh[2,11]<--4*sigma3[4]/s0[4]
      hh[3,3]<-sigma*(1/n_samP1+1/n_samP2)+4*sigma1[1]/n0[1]+4*sigma3[6]/s0[6]+36*sigma3[5]/s0[5]+36*sigma1[2]/n0[2]
      hh[3,4]<-sigma*(1/n_samP1-1/n_samP2)-2*sigma1[1]/n0[1]
      hh[3,5]<-2*sigma1[1]/n0[1]+12*sigma1[2]/n0[2]
      hh[3,6]<-6*sigma3[5]/s0[5]
      hh[3,7]<-2*sigma1[1]/n0[1]+2*sigma3[6]/s0[6]
      hh[3,8]<--6*sigma1[2]/n0[2]-6*sigma3[5]/s0[5]
      hh[3,9]<-0
      hh[3,10]<-6*sigma3[5]/s0[5]
      hh[3,11]<--2*sigma3[6]/s0[6]-12*sigma3[5]/s0[5]
      hh[4,4]<-sigma*(1/n_samP1+1/n_samP2)+sigma1[1]/n0[1]+sigma1[3]/n0[3]+sigma2[1]/n0[4]+sigma2[3]/n0[6]
      hh[4,5]<--sigma1[1]/n0[1]-sigma1[3]/n0[3]
      hh[4,6]<-0
      hh[4,7]<--sigma1[1]/n0[1]+sigma2[3]/n0[6]
      hh[4,8]<-0
      hh[4,9]<--sigma1[3]/n0[3]+sigma2[1]/n0[4]
      hh[4,10]<-0
      hh[4,11]<-0
      hh[5,5]<-sigma1[1]/n0[1]+4*sigma1[2]/n0[2]+sigma1[3]/n0[3]
      hh[5,6]<-0
      hh[5,7]<-sigma1[1]/n0[1]
      hh[5,8]<--2*sigma1[2]/n0[2]
      hh[5,9]<-sigma1[3]/n0[3]
      hh[5,10]<-0
      hh[5,11]<-0
      hh[6,6]<-9*sigma3[1]/s0[1]+9*sigma3[2]/s0[2]+sigma3[4]/s0[4]+sigma3[5]/s0[5]
      hh[6,7]<--3*sigma3[1]/s0[1]
      hh[6,8]<-3*sigma3[2]/s0[2]-sigma3[5]/s0[5]
      hh[6,9]<-2*sigma3[4]/s0[4]
      hh[6,10]<--3*sigma3[2]/s0[2]+sigma3[4]/s0[4]+sigma3[5]/s0[5]
      hh[6,11]<--sigma3[4]/s0[4]-2*sigma3[5]/s0[5]
      hh[7,7]<-sigma1[1]/n0[1]+sigma2[3]/n0[6]+sigma3[1]/s0[1]+sigma3[6]/s0[6]
      hh[7,8]<-0
      hh[7,9]<-0
      hh[7,10]<-0
      hh[7,11]<--sigma3[6]/s0[6]
      hh[8,8]<-sigma1[2]/n0[2]+sigma2[2]/n0[5]+sigma3[2]/s0[2]+sigma3[5]/s0[5]
      hh[8,9]<-0
      hh[8,10]<--sigma3[2]/s0[2]-sigma3[5]/s0[5]
      hh[8,11]<-2*sigma3[5]/s0[5]
      hh[9,9]<-sigma1[3]/n0[3]+sigma2[1]/n0[4]+4*sigma3[4]/s0[4]
      hh[9,10]<-2*sigma3[4]/s0[4]
      hh[9,11]<--2*sigma3[4]/s0[4]
      hh[10,10]<-sigma3[2]/s0[2]+sigma3[3]/s0[3]+sigma3[4]/s0[4]+sigma3[5]/s0[5]
      hh[10,11]<--sigma3[4]/s0[4]-2*sigma3[5]/s0[5]
      hh[11,11]<-sigma3[4]/s0[4]+sigma3[6]/s0[6]+4*sigma3[5]/s0[5]
      for(i in 2:11)
      {
        for(j in 1:(i-1))
        {
          hh[i,j]<-hh[j,i]
        }
      }
      b_line[1]<-sumx[1]/n_samP1+sumx[2]/n_samF1+2*sumx[3]/n_samP2+sumwx1[3]/n0[3]-sumwx2[1]/n0[4]-2*sumwx3[1]/s0[1]-2*sumwx3[6]/s0[6]
      b_line[2]<-3*sumx[1]/n_samP1+2*sumx[2]/n_samF1+3*sumx[3]/n_samP2-2*sumwx1[1]/n0[1]-2*sumwx2[3]/n0[6]-4*sumwx3[4]/s0[4]
      b_line[3]<-sumx[1]/n_samP1-sumx[3]/n_samP2+2*sumwx1[1]/n0[1]-2*sumwx3[6]/s0[6]+6*sumwx3[5]/s0[5]-6*sumwx1[2]/n0[2]
      b_line[4]<-sumx[1]/n_samP1-sumx[3]/n_samP2-sumwx1[1]/n0[1]-sumwx1[3]/n0[3]+sumwx2[1]/n0[4]+sumwx2[3]/n0[6]
      b_line[5]<-sumwx1[1]/n0[1]-2*sumwx1[2]/n0[2]+sumwx1[3]/n0[3]
      b_line[6]<-3*sumwx3[1]/s0[1]-3*sumwx3[2]/s0[2]-sumwx3[4]/s0[4]+sumwx3[5]/s0[5]
      b_line[7]<-sumwx1[1]/n0[1]+sumwx2[3]/n0[6]-sumwx3[1]/s0[1]-sumwx3[6]/s0[6]
      b_line[8]<-sumwx1[2]/n0[2]+sumwx2[2]/n0[5]-sumwx3[2]/s0[2]-sumwx3[5]/s0[5]
      b_line[9]<-sumwx1[3]/n0[3]+sumwx2[1]/n0[4]-2*sumwx3[4]/s0[4]
      b_line[10]<-sumwx3[2]/s0[2]-sumwx3[3]/s0[3]-sumwx3[4]/s0[4]+sumwx3[5]/s0[5]
      b_line[11]<-sumwx3[4]/s0[4]-2*sumwx3[5]/s0[5]+sumwx3[6]/s0[6]
      B<-solve(hh,b_line)
      mean[1]<-(sumx[1]-sigma*(B[1]+3*B[2]+B[3]+B[4]))/n_samP1
      mean[2]<-(sumx[2]-sigma*(B[1]+2*B[2]))/n_samF1
      mean[3]<-(sumx[3]-sigma*(2*B[1]+3*B[2]-B[3]-B[4]))/n_samP2
      mean1[1]<-(sumwx1[1]+sigma1[1]*(2*B[2]-2*B[3]+B[4]-B[5]-B[7]))/n0[1]
      mean1[2]<-(sumwx1[2]+sigma1[2]*(6*B[3]+2*B[5]-B[8]))/n0[2]
      mean1[3]<-(sumwx1[3]+sigma1[3]*(-B[1]+B[4]-B[5]-B[9]))/n0[3]
      mean2[1]<-(sumwx2[1]+sigma2[1]*(B[1]-B[4]-B[9]))/n0[4]
      mean2[2]<-(sumwx2[2]+sigma2[2]*(-B[8]))/n0[5]
      mean2[3]<-(sumwx2[3]+sigma2[3]*(2*B[2]-B[4]-B[7]))/n0[6]
      mean3[1]<-(sumwx3[1]+sigma3[1]*(2*B[1]-3*B[6]+B[7]))/s0[1]
      mean3[2]<-(sumwx3[2]+sigma3[2]*(3*B[6]+B[8]-B[10]))/s0[2]
      mean3[3]<-(sumwx3[3]+sigma3[3]*B[10])/s0[3]
      mean3[4]<-(sumwx3[4]+sigma3[4]*(4*B[2]+B[6]+2*B[9]+B[10]-B[11]))/s0[4]
      mean3[5]<-(sumwx3[5]+sigma3[5]*(-6*B[3]-B[6]+B[8]-B[10]+2*B[11]))/s0[5]
      mean3[6]<-(sumwx3[6]+sigma3[6]*(2*B[1]+2*B[3]+B[7]-B[11]))/s0[6]
      aaa1<-max(abs(B-AA))
      AA<-B
      if (n_iter>20) break
    }
    ##########obtain variance############################
    ss1<-sum((dataP1-mean[1])^2);ss3<-sum((dataP2-mean[3])^2);ss2<-sum((dataF1-mean[2])^2)
    for(i in 1:3) {swx1[i] <- W1[i,]%*%(dataB1-mean1[i])^2 };for(i in 1:3) {swx2[i] <- W2[i,]%*%(dataB2-mean2[i])^2 };for(i in 1:6) {swx3[i] <- W3[i,]%*%(dataF2-mean3[i])^2 }
    aaa0<-sigma1[1];n_iter<-0;aaa1<-1000
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      aa2<-sigma1[1]/(sigma1[1]+g_aa1)
      aa3<-sigma1[1]/(sigma1[1]+g_aa2)
      as1<-swx1[1]+swx1[2]*aa2^2+swx1[3]*aa3^2
      as2<-n0[1]+aa2*n0[2]+aa3*n0[3]
      sigma1[1]<-as1/as2
      aaa1<-abs(sigma1[1]-aaa0)
      aaa0<-sigma1[1]
      if (n_iter>20) break
    }
    sigma40<-sigma1[1]-sigma;
    if (sigma40<0) {sigma40<-0;sigma1[1]<-sigma}
    sigma1[1]<-sigma40+sigma;sigma1[2]<-sigma1[1]+g_aa1;sigma1[3]<-sigma1[1]+g_aa2
    # to estimate sigma50.
    aaa0<-sigma2[3];n_iter<-0;aaa1<-1000
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma2[3]/(sigma2[3]+g_aa2)
      aa2<-sigma2[3]/(sigma2[3]+g_aa1)
      as3<-swx2[1]*aa1^2+swx2[2]*aa2^2+swx2[3]
      as4<-aa1*n0[4]+aa2*n0[5]+n0[6]
      sigma2[3]<-as3/as4
      aaa1<-abs(sigma2[3]-aaa0)
      aaa0<-sigma2[3]
      if (n_iter>20) break
    }
    sigma50<-sigma2[3]-sigma;
    if (sigma50<0) {sigma50<-0;sigma2[3]<-sigma}
    sigma2[3]<-sigma+sigma50;sigma2[1]<-sigma2[3]+g_aa2;sigma2[2]<-sigma2[3]+g_aa1
    # to estimate sigma60.
    aaa0<-sigma3[1];aa6<-swx3[1]+swx3[3]+swx3[6];aa7<-s0[1]+s0[3]+s0[6]
    n_iter<-0;aaa1<-1000
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma3[1]/(sigma3[1]+g_aa1)
      aa2<-sigma3[1]/(sigma3[1]+g_aa2)
      as5<-aa6+(swx3[2]+swx3[5])*aa1^2+swx3[4]*aa2^2
      as6<-aa7+aa1*(s0[2]+s0[5])+aa2*s0[4]
      sigma3[1]<-as5/as6
      aaa1<-abs(sigma3[1]-aaa0)
      aaa0<-sigma3[1]
      if (n_iter>20) break
    }
    sigma60<-sigma3[1]-sigma;
    if (sigma60<0) {sigma60<-0;sigma3[1]<-sigma}
    sigma3[1]<-sigma+sigma60;sigma3[2]<-sigma3[1]+g_aa1
    sigma3[3]<-sigma3[6]<-sigma3[1]
    sigma3[4]<-sigma3[1]+g_aa2
    sigma3[5]<-sigma3[2]
    # CM3 to estimate the variance (sigma).
    ab1<-ss1+ss2+ss3;ab2<-n_samP1+n_samF1+n_samP2
    n_iter<-0;aaa0<-sigma;aaa1<-1000
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      n0[11]<-sigma/(sigma+sigma40)
      n0[12]<-sigma/(sigma+sigma40+g_aa1)
      n0[13]<-sigma/(sigma+sigma40+g_aa2)
      s0[11]<-sigma/(sigma+sigma50+g_aa2)
      s0[12]<-sigma/(sigma+sigma50+g_aa1)
      s0[13]<-sigma/(sigma+sigma50)
      ab3<-sum(swx1[c(1:3)]*n0[c(11:13)]^2+swx2[c(1:3)]*s0[c(11:13)]^2)
      ab4<-sum(n0[c(1:3)]*n0[c(11:13)]+n0[c(4:6)]*s0[c(11:13)])
      n0[11]<-sigma/(sigma+sigma60)
      n0[13]<-n0[16]<-n0[11]
      n0[12]<-sigma/(sigma+sigma60+g_aa1)
      n0[14]<-sigma/(sigma+sigma60+g_aa2)
      n0[15]<-n0[12]
      ab3<-ab3+sum(swx3[c(1:6)]*n0[c(11:16)]^2)
      ab4<-ab4+sum(s0[c(1:6)]*n0[11:16])
      sigma<-(ab1+ab3)/(ab2+ab4)
      aaa1<-abs(sigma-aaa0)
      aaa0<-sigma
      if (n_iter>20) break
    }
    sigma1[1]<-sigma+sigma40;sigma1[2]<-sigma1[1]+g_aa1;sigma1[3]<-sigma1[1]+g_aa2
    sigma2[3]<-sigma+sigma50;sigma2[1]<-sigma2[3]+g_aa2;sigma2[2]<-sigma2[3]+g_aa1
    sigma3[1]<-sigma+sigma60;sigma3[3]<-sigma3[6]<-sigma3[1];sigma3[2]<-sigma3[1]+g_aa1
    sigma3[5]<-sigma3[2];sigma3[4]<-sigma3[1]+g_aa2
    if(sum(sigma < 1e-30)>=1){break}
    #############################criteria for iterations to stop########################
    L1<-sum(log(dnorm(dataP1,mean[1],sqrt(sigma))))+sum(log(dnorm(dataF1,mean[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,mean[3],sqrt(sigma))))+sum(log(dmixnorm(dataB1,mean1,sqrt(sigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mean2,sqrt(sigma2),mix_pi2)))+sum(log(dmixnorm(dataF2,mean3,sqrt(sigma3),mix_pi3)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*8
  meanP1<-mean[1];meanF1<-mean[2];meanP2<-mean[3]
  sigma0<-sigma

  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*n_samP1)
  P1bmw <- matrix(0,n_samP1,1)
  P1gg <- (dataP1 - mean[1])/sqrt(as.vector(sigma))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < n_samP1){P1bmw <- P1bmw+runif(n_samP1)/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:n_samP1)) - 0.5)/n_samP1)^2)
  P1u<- as.matrix(c(12*n_samP1*((P1dd[1]/n_samP1-0.5)^2),((45*n_samP1)/4)*((P1dd[2]/n_samP1-1/3)^2),180*n_samP1*((P1dd[3]/n_samP1-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[2])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),P1D))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*n_samF1)
  F1bmw <- matrix(0,n_samF1,1)
  F1gg <- (dataF1 - mean[2])/sqrt(as.vector(sigma))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < n_samF1){F1bmw <- F1bmw+runif(n_samF1)/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:n_samF1)) - 0.5)/n_samF1)^2)
  F1u<- as.matrix(c(12*n_samF1*((F1dd[1]/n_samF1-0.5)^2),((45*n_samF1)/4)*((F1dd[2]/n_samF1-1/3)^2),180*n_samF1*((F1dd[3]/n_samF1-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[2])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),F1D))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*n_samP2)
  P2bmw <- matrix(0,n_samP2,1)
  P2gg <- (dataP2 - mean[3])/sqrt(as.vector(sigma))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < n_samP2){P2bmw <- P2bmw+runif(n_samP2)/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:n_samP2)) - 0.5)/n_samP2)^2)
  P2u<- as.matrix(c(12*n_samP2*((P2dd[1]/n_samP2-0.5)^2),((45*n_samP2)/4)*((P2dd[2]/n_samP2-1/3)^2),180*n_samP2*((P2dd[3]/n_samP2-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[2])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),P2D))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1);
  B1w1<-1/(12*n_samB1)
  B1bmw <- matrix(0,n_samB1,1); B1bmwsl <- matrix(0,n_samB1,3)
  for(i in 1:3){
    B1gg <- (dataB1 - mean1[i])/sqrt(sigma1[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi1[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B1P2)))[1]
  if(nn < n_samB1){B1P2 <- B1P2+runif(n_samB1)/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*n_samB1) + sum((B1P2 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  B1u <- as.matrix(c(12*n_samB1*((B1dd[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((B1dd[2]/n_samB1-1/3)^2),180*n_samB1*((B1dd[3]/n_samB1-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[2])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),B1D))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])

  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2);
  B2w1<-1/(12*n_samB2)
  B2bmw <- matrix(0,n_samB2,1); B2bmwsl <- matrix(0,n_samB2,3)
  for(i in 1:3){
    B2gg <- (dataB2 - mean2[i])/sqrt(sigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(B2P2)))[1]
  if(nn < n_samB2){B2P2 <- B2P2+runif(n_samB2)/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*n_samB2) + sum((B2P2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  B2u <- as.matrix(c(12*n_samB2*((B2dd[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((B2dd[2]/n_samB2-1/3)^2),180*n_samB2*((B2dd[3]/n_samB2-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[2])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),B2D))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])

  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2);
  F2w1<-1/(12*n_samF2)
  F2bmw <- matrix(0,n_samF2,1); F2bmwsl <- matrix(0,n_samF2,6)
  for(i in 1:6){
    F2gg <- (dataF2 - mean3[i])/sqrt(sigma3[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi3[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < n_samF2){F2P2 <- F2P2+runif(n_samF2)/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*n_samF2) + sum((F2P2 - (as.matrix(c(1:n_samF2)) - 0.5)/n_samF2)^2)
  F2u <- as.matrix(c(12*n_samF2*((F2dd[1]/n_samF2-0.5)^2),((45*n_samF2)/4)*((F2dd[2]/n_samF2-1/3)^2),180*n_samF2*((F2dd[3]/n_samF2-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[2])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),F2D))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])

  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)

  #########first order parameters######################
  aa<-matrix(c(1,2,1,0,1,0,0,1,1,-2,-1,0,1,2,0.5,0.25,1,1.5,0.5,0.25,1,1,0.5,0.25,1,1,
               -0.5,0.25,1,-0.5,-0.5,0.25,1,-2,-0.5,0.25,1,2,0,0.25,1,1.5,0,0.25,1,0,0,
               0.25,1,1,0,0.25,1,-0.5,0,0.25,1,-2,0,0.25),15,4,byrow=T)
  b_line1<-as.matrix(c(mean[1],mean[2],mean[3],mean1[1],mean1[2],mean1[3],mean2[1],mean2[2],mean2[3],mean3[1],mean3[2],mean3[3],mean3[4],mean3[5],mean3[6]))
  B1<-solve(crossprod(aa,aa))%*%crossprod(aa,b_line1)
  ##########second order parameters######################
  jj1<-sigmaB1-sigma1[1]
  if (jj1<0 || jj1>=sigmaB1) {jj1<-0}
  ll1<-jj1/sigmaB1
  mm1<-sigma1[1]-sigma
  if (mm1<0 || mm1>=sigmaB1) {mm1<-0}
  nn1<-mm1/sigmaB1
  jj2<-sigmaB2-sigma2[3]
  if (jj2<0 || jj2>=sigmaB2) {jj2<-0}
  ll2<-jj2/sigmaB2
  mm2<-sigma2[3]-sigma
  if (mm2<0 || mm2>=sigmaB2) {mm2<-0}
  nn2<-mm2/sigmaB2
  jj3<-sigmaF2-sigma3[1]
  if (jj3<0 || jj3>=sigmaF2) {jj3<-0}
  ll3<-jj3/sigmaF2
  mm3<-sigma3[1]-sigma
  if (mm3<0 || mm3>=sigmaF2) {mm3<-0}
  nn3<-mm3/sigmaF2

  output <- data.frame("MX2-EAD-AD",round(abc,4),round(AIC,4),round(meanP1,4),round(meanF1,4),round(meanP2,4), round(sigma0,4),round(t(mean1),4)," ",round(t(sigma1),4)," ",
                       round(t(mix_pi1),4)," ",round(t(mean2),4)," ",round(t(sigma2),4)," ",round(t(mix_pi2),4)," ",
                       round(t(mean3),4)," "," "," ",round(t(sigma3),4)," "," "," ",round(t(mix_pi3),4)," "," "," ",
                       round(B1[1],4)," "," "," "," "," ",round(B1[2],4),round(B1[2],4),round(B1[2],4),round(B1[2],4)," "," "," "," ",round(B1[3],4),round(B1[4],4),
                       round(jj1,4),round(ll1*100,4),round(mm1,4),round(nn1*100,4),round(jj2,4),round(ll2*100,4),round(mm2,4),round(nn2*100,4),round(jj3,4),round(ll3*100,4),round(mm3,4),round(nn3*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi1,mi2,mi3)
  return(OUTPUT)
}



K1G6F <- function(x){
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

logLG6F <- function(nm,nng,mi,mn,s,d1) { sum2 <- sum(log(dmixnorm(d1,mn,sqrt(s),mi)));return (sum2) }


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
  G6FModelFun[[i]](K1G6F,logLG6F,df11,df21,df31,df41,df51,df61,G6Ftext2)[[1]]
}
stopCluster(cl)

mi1<-NULL;mi2<-NULL;mi3<-NULL

}else{

allresultq=switch(model,"1MG-AD" = G6FModelFun[[1]](K1G6F,logLG6F,df11,df21,df31,df41,df51,df61,G6Ftext2),"1MG-A"=G6FModelFun[[2]](K1G6F,logLG6F,df11,df21,df31,df41,df51,df61,G6Ftext2),"1MG-EAD"=G6FModelFun[[3]](K1G6F,logLG6F,df11,df21,df31,df41,df51,df61,G6Ftext2),"1MG-NCD"=G6FModelFun[[4]](K1G6F,logLG6F,df11,df21,df31,df41,df51,df61,G6Ftext2),"2MG-ADI"=G6FModelFun[[5]](K1G6F,logLG6F,df11,df21,df31,df41,df51,df61,G6Ftext2),
                        "2MG-AD"=G6FModelFun[[6]](K1G6F,logLG6F,df11,df21,df31,df41,df51,df61,G6Ftext2),"2MG-A"=G6FModelFun[[7]](K1G6F,logLG6F,df11,df21,df31,df41,df51,df61,G6Ftext2),"2MG-EA"=G6FModelFun[[8]](K1G6F,logLG6F,df11,df21,df31,df41,df51,df61,G6Ftext2),"2MG-CD"=G6FModelFun[[9]](K1G6F,logLG6F,df11,df21,df31,df41,df51,df61,G6Ftext2),"2MG-EAD"=G6FModelFun[[10]](K1G6F,logLG6F,df11,df21,df31,df41,df51,df61,G6Ftext2),
                        "PG-ADI"=G6FModelFun[[11]](K1G6F,logLG6F,df11,df21,df31,df41,df51,df61,G6Ftext2),"PG-AD"=G6FModelFun[[12]](K1G6F,logLG6F,df11,df21,df31,df41,df51,df61,G6Ftext2),"MX1-AD-ADI"=G6FModelFun[[13]](K1G6F,logLG6F,df11,df21,df31,df41,df51,df61,G6Ftext2),"MX1-AD-AD"=G6FModelFun[[14]](K1G6F,logLG6F,df11,df21,df31,df41,df51,df61,G6Ftext2),"MX1-A-AD"=G6FModelFun[[15]](K1G6F,logLG6F,df11,df21,df31,df41,df51,df61,G6Ftext2),
                        "MX1-EAD-AD"=G6FModelFun[[16]](K1G6F,logLG6F,df11,df21,df31,df41,df51,df61,G6Ftext2),"MX1-NCD-AD"=G6FModelFun[[17]](K1G6F,logLG6F,df11,df21,df31,df41,df51,df61,G6Ftext2),"MX2-ADI-ADI"=G6FModelFun[[18]](K1G6F,logLG6F,df11,df21,df31,df41,df51,df61,G6Ftext2),"MX2-ADI-AD"=G6FModelFun[[19]](K1G6F,logLG6F,df11,df21,df31,df41,df51,df61,G6Ftext2),"MX2-AD-AD"=G6FModelFun[[20]](K1G6F,logLG6F,df11,df21,df31,df41,df51,df61,G6Ftext2),
                        "MX2-A-AD"=G6FModelFun[[21]](K1G6F,logLG6F,df11,df21,df31,df41,df51,df61,G6Ftext2),"MX2-EA-AD"=G6FModelFun[[22]](K1G6F,logLG6F,df11,df21,df31,df41,df51,df61,G6Ftext2),"MX2-CD-AD"=G6FModelFun[[23]](K1G6F,logLG6F,df11,df21,df31,df41,df51,df61,G6Ftext2),"MX2-EAD-AD"=G6FModelFun[[24]](K1G6F,logLG6F,df11,df21,df31,df41,df51,df61,G6Ftext2))



allresult<-allresultq[[1]]
if(model=="PG-AD"||model=="PG-ADI"){
  mi1<-NULL;mi2<-NULL;mi3<-NULL
}else{
  mi1<-allresultq[[2]];mi2<-allresultq[[3]];mi3<-allresultq[[4]]
}
}
colnames(allresult) <- G6Fcolname
out<-list(allresult,mi1,mi2,mi3)
return(out)
}




