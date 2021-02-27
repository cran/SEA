BCFun<-function(df,model){

data<-sapply(df,as.character)
dB1<-data[-1,which(data[1,]=="B1")];B1<-as.numeric(dB1[which(is.na(as.numeric(dB1))==FALSE)]);df11<-as.data.frame(B1)
dB2<-data[-1,which(data[1,]=="B2")];B2<-as.numeric(dB2[which(is.na(as.numeric(dB2))==FALSE)]);df21<-as.data.frame(B2)


#########################################Column names##################################################################
BCcolname <- c("Model","Log_Max_likelihood_value","AIC","B1-mean[1]","B1-mean[2]","B1-mean[3]","B1-mean[4]","B1-Var(Residual+Polygene)",
               "B1-Proportion[1]","B1-Proportion[2]","B1-Proportion[3]","B1-Proportion[4]","B2-mean[1]","B2-mean[2]","B2-mean[3]","B2-mean[4]",
               "B2-Var(Residual+Polygene)","B2-Proportion[1]","B2-Proportion[2]","B2-Proportion[3]","B2-Proportion[4]","m1","m2","da","db","ha",
               "hb","B1-Major-Gene Var ","B1-Heritability(Major-Gene)(%)","B2-Major-Gene Var ","B2-Heritability(Major-Gene)(%)","U1 square-B1","p(U1 square-B1)",
               "U2 square-B1","p(U2 square-B1)","U3 square-B1","p(U3 square-B1)","nW square-B1","p(nW square-B1)","Dn-B1","p(Dn-B1)","U1 square-B2","p(U1 square-B2)",
               "U2 square-B2","p(U2 square-B2)","U3 square-B2","p(U3 square-B2)","nW square-B2","p(nW square-B2)","Dn-B2","p(Dn-B2)")

BCModelFun<-list(NA)
###################define each model function##################
####################(A-0)Model##########################################
BCModelFun[[1]] <- function(K1,logL,df11,df21){
  dataB1 <- as.matrix(as.numeric(df11[,1])); dataB2 <- as.matrix(as.numeric(df21[,1]))
  n_samB1 <- dim(dataB1)[1];n_samB2 <- dim(dataB2)[1]
  mean1<-mean(dataB1);mean2<-mean(dataB2)
  sigma1<- as.numeric(var(dataB1)); sigma2<- as.numeric(var(dataB2))

  m_esp<-0.0001
  ####################(A-0)Model##########################################
  d2<- 1
  mix_pi_1<-1.0;mix_pi_2<-1.0
  mean11<-mean1;mean22<-mean2
  sigma11<-sigma1;sigma22<-sigma2

  abc <-logL(n_samB1,d2,mix_pi_1,mean11,sigma11,dataB1)+logL(n_samB2,d2,mix_pi_2,mean22,sigma22,dataB2)
  AIC <- -2.0*abc+2.0*2.0
  ######### hypothesis testing #########
  ###############  B1  ###################
  dataB1<-sort(dataB1)
  bmw_B1 <- matrix(0,n_samB1,1); bmwsl_B1 <- matrix(0,n_samB1,1)
  gg_B1 <- (dataB1 - mean11)/sqrt(as.vector(sigma11))
  bmw_B1[which(gg_B1>=0)] <- pnorm(gg_B1[gg_B1>=0])
  bmw_B1[which(gg_B1<0)] <- 1 - pnorm(abs(gg_B1[gg_B1<0]))
  bmwsl_B1[,1] <- bmw_B1

  P2_B1 <- rowSums(bmwsl_B1)
  nn<-dim(as.matrix(unique(P2_B1)))[1]
  if(nn<n_samB1){P2_B1<-P2_B1+runif(n_samB1)/1e4}

  dd_B1 <- as.matrix(c(sum(P2_B1),sum(P2_B1^2),sum((P2_B1-0.5)^2)))
  WW2_B1 <- 1/(12*n_samB1) + sum((P2_B1 - (as.matrix(c(1:n_samB1)) - 0.5)/n_samB1)^2)
  u_B1 <- as.matrix(c(12*n_samB1*((dd_B1[1]/n_samB1-0.5)^2),((45*n_samB1)/4)*((dd_B1[2]/n_samB1-1/3)^2),180*n_samB1*((dd_B1[3]/n_samB1-1/12)^2)))
  D_B1 <- as.numeric(ks.test(P2_B1,"punif")[[1]][1])
  tt_B1 <- as.matrix(c((1 - pchisq(u_B1[1],1)),(1 - pchisq(u_B1[2],1)),(1 - pchisq(u_B1[3],1)),K1(WW2_B1),(1-pkolm(D_B1,n_samB1))))

  ###############  B2  ###################
  dataB2<-sort(dataB2);bmw_B2 <- matrix(0,n_samB2,1); bmwsl_B2 <- matrix(0,n_samB2,1)
  gg_B2 <- (dataB2 - mean22)/sqrt(as.vector(sigma22))
  bmw_B2[which(gg_B2>=0)] <- pnorm(gg_B2[gg_B2>=0])
  bmw_B2[which(gg_B2<0)] <- 1 - pnorm(abs(gg_B2[gg_B2<0]))
  bmwsl_B2[,1] <- bmw_B2

  P2_B2 <- rowSums(bmwsl_B2)
  nn<-dim(as.matrix(unique(P2_B2)))[1]
  if(nn<n_samB2){P2_B2<-P2_B2+runif(n_samB2)/1e4}

  dd_B2 <- as.matrix(c(sum(P2_B2),sum(P2_B2^2),sum((P2_B2-0.5)^2)))
  WW2_B2 <- 1/(12*n_samB2) + sum((P2_B2 - (as.matrix(c(1:n_samB2)) - 0.5)/n_samB2)^2)
  u_B2 <- as.matrix(c(12*n_samB2*((dd_B2[1]/n_samB2-0.5)^2),((45*n_samB2)/4)*((dd_B2[2]/n_samB2-1/3)^2),180*n_samB2*((dd_B2[3]/n_samB2-1/12)^2)))
  D_B2 <- as.numeric(ks.test(P2_B2,"punif")[[1]][1])
  tt_B2 <- as.matrix(c((1 - pchisq(u_B2[1],1)),(1 - pchisq(u_B2[2],1)),(1 - pchisq(u_B2[3],1)),K1(WW2_B2),(1-pkolm(D_B2,n_samB2))))

  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)


  output <- data.frame("0MG",round(abc,4),round(AIC,4),round(t(mean11),4)," "," "," ",round(sigma11[1],4),round(t(mix_pi_1),4)," "," "," ",round(t(mean22),4)," "," "," ",round(sigma22[1],4),round(t(mix_pi_2),4),
                       " "," "," "," "," "," "," "," "," "," "," "," "," ",round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],
                       round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],
                       round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output)
  return(OUTPUT)
}
####################(A-1)Model##########################################
BCModelFun[[2]] <- function(K1,logL,df11,df21){
  dataB1 <- as.matrix(as.numeric(df11[,1])); dataB2 <- as.matrix(as.numeric(df21[,1]))
  n_samB1 <- dim(dataB1)[1];n_samB2 <- dim(dataB2)[1]
  mean1<-mean(dataB1);mean2<-mean(dataB2)
  sigma1<- as.numeric(var(dataB1)); sigma2<- as.numeric(var(dataB2))
  m_esp<-0.0001
  ####################(A-1)Model##########################################
  d2<-2
  mi_1 <- as.matrix(c(0.5,0.5))
  sigma11 <- matrix((sigma1/5),d2,1)
  a1 <- sqrt(sigma1/n_samB1)
  mean11 <- as.matrix(c((mean1+2*a1),(mean1-2*a1)))

  mi_2 <- as.matrix(c(0.5,0.5))
  sigma22 <- matrix((sigma2/5),d2,1)
  a2 <- sqrt(sigma2/n_samB2)
  mean22 <- as.matrix(c((mean2+2*a2),(mean2-2*a2)))

  L0 <- logL(n_samB1,d2,mi_1,mean11,sigma11,dataB1)+logL(n_samB2,d2,mi_2,mean22,sigma22,dataB2)
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d2,n_samB1); swx_B1 <- matrix(0,d2,1)
  WW_B2 <- matrix(0,d2,n_samB2); swx_B2 <- matrix(0,d2,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d2) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean11[i],sqrt(sigma11[i]))/dmixnorm(dataB1,mean11,sqrt(sigma11),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d2) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean22[i],sqrt(sigma22[i]))/dmixnorm(dataB2,mean22,sqrt(sigma22),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2
    ############ CM1-step for means ##############
    n0_1 <- n_samB1*mix_pi_1
    n0_1[n0_1<0.000001] <- 0.000001

    n0_2 <- n_samB2*mix_pi_2
    n0_2[n0_2<0.000001] <- 0.000001

    mean11<- sumwx_B1/n0_1
    mean22<- sumwx_B2/n0_2

    ##########obtain variance##########
    for(i in 1:d2) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean11[i])^2 }
    for(i in 1:d2) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean22[i])^2 }

    sigma11<- matrix((sum(swx_B1)/n_samB1),d2,1)
    sigma22<- matrix((sum(swx_B2)/n_samB2),d2,1)

    ########criteria for iterations to stop#######

    L1 <- logL(n_samB1,d2,mix_pi_1,mean11,sigma11,dataB1)+logL(n_samB2,d2,mix_pi_2,mean22,sigma22,dataB2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*6

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,0,0, 0,0,1,1, 1,0,0,-1, 0,1,1,0),4,4)
  b_line1 <- matrix(c(mean11,mean22))
  B1 <- solve(aa,b_line1)

  jj_1 <- sigma1 - sigma11[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigma1

  jj_2 <- sigma2 - sigma22[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigma2

  ######### hypothesis testing #########
  ###############  B1  ###################
  dataB1<-sort(dataB1);bmw_B1 <- matrix(0,n_samB1,1); bmwsl_B1 <- matrix(0,n_samB1,d2)
  for(i in 1:d2){
    gg_B1 <- (dataB1 - mean11[i])/sqrt(sigma11[i])
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
  dataB2<-sort(dataB2);bmw_B2 <- matrix(0,n_samB2,1); bmwsl_B2 <- matrix(0,n_samB2,d2)
  for(i in 1:d2){
    gg_B2 <- (dataB2 - mean22[i])/sqrt(sigma22[i])
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

  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("1MG-AD",round(abc,4),round(AIC,4),round(t(mean11),4)," "," ",round(sigma11[1],4),round(t(mix_pi_1),4)," "," ",round(t(mean22),4)," "," ",round(sigma22[1],4),round(t(mix_pi_2),4),
                       " "," ",round(B1[1],4),round(B1[2],4),round(B1[3],4)," ",round(B1[4],4)," ",round(jj_1,4),round(ll_1*100,4),round(jj_2,4),round(ll_2*100,4),round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],
                       round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],
                       round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
}
####################(A-2)Model##########################################
BCModelFun[[3]] <- function(K1,logL,df11,df21){
  dataB1 <- as.matrix(as.numeric(df11[,1])); dataB2 <- as.matrix(as.numeric(df21[,1]))
  n_samB1 <- dim(dataB1)[1];n_samB2 <- dim(dataB2)[1]
  mean1<-mean(dataB1);mean2<-mean(dataB2)
  sigma1<- as.numeric(var(dataB1)); sigma2<- as.numeric(var(dataB2))
  m_esp<-0.0001
  ####################(A-2)Model##########################################
  d2<-2
  mi_1 <- as.matrix(c(0.5,0.5))
  sigma11 <- matrix((sigma1/5),d2,1)
  a1 <- sqrt(sigma1/n_samB1)
  mean11 <- as.matrix(c((mean1+2*a1),(mean1-2*a1)))

  mi_2 <- as.matrix(c(0.5,0.5))
  sigma22 <- matrix((sigma2/5),d2,1)
  a2 <- sqrt(sigma2/n_samB2)
  mean22 <- as.matrix(c((mean2+2*a2),(mean2-2*a2)))

  L0 <- logL(n_samB1,d2,mi_1,mean11,sigma11,dataB1)+logL(n_samB2,d2,mi_2,mean22,sigma22,dataB2)
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d2,n_samB1); swx_B1 <- matrix(0,d2,1)
  WW_B2 <- matrix(0,d2,n_samB2); swx_B2 <- matrix(0,d2,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d2) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean11[i],sqrt(sigma11[i]))/dmixnorm(dataB1,mean11,sqrt(sigma11),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d2) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean22[i],sqrt(sigma22[i]))/dmixnorm(dataB2,mean22,sqrt(sigma22),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2
    ############ CM1-step for means ##############
    n0_1 <- n_samB1*mix_pi_1
    n0_1[n0_1<0.000001] <- 0.000001

    n0_2 <- n_samB2*mix_pi_2
    n0_2[n0_2<0.000001] <- 0.000001

    aa1<- sum(sigma11/n0_1)+ sum(sigma22/n0_2)
    aa2<- sumwx_B1[1]/n0_1[1]-sumwx_B1[2]/n0_1[2]-sumwx_B2[1]/n0_2[1]+sumwx_B2[2]/n0_2[2]
    aa3<- aa2/aa1

    mean11[1]<- (sumwx_B1[1]-sigma11[1]*aa3)/n0_1[1]
    mean11[2]<- (sumwx_B1[2]+sigma11[2]*aa3)/n0_1[2]
    mean22[1]<- (sumwx_B2[1]+sigma22[1]*aa3)/n0_2[1]
    mean22[2]<- (sumwx_B2[2]-sigma22[2]*aa3)/n0_2[2]
    ##########obtain variance##########
    for(i in 1:d2) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean11[i])^2 }
    for(i in 1:d2) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean22[i])^2 }

    sigma11<- matrix((sum(swx_B1)/n_samB1),d2,1)
    sigma22<- matrix((sum(swx_B2)/n_samB2),d2,1)

    ########criteria for iterations to stop#######

    L1 <- logL(n_samB1,d2,mix_pi_1,mean11,sigma11,dataB1)+logL(n_samB2,d2,mix_pi_2,mean22,sigma22,dataB2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*5

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,0,0, 0,0,1,1, 1,0,0,-1),4,3)
  b_line1 <- matrix(c(mean11,mean22))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigma1 - sigma11[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigma1

  jj_2 <- sigma2 - sigma22[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigma2

  ######### hypothesis testing #########
  ###############  B1  ###################
  dataB1<-sort(dataB1);bmw_B1 <- matrix(0,n_samB1,1); bmwsl_B1 <- matrix(0,n_samB1,d2)
  for(i in 1:d2){
    gg_B1 <- (dataB1 - mean11[i])/sqrt(sigma11[i])
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
  dataB2<-sort(dataB2);bmw_B2 <- matrix(0,n_samB2,1); bmwsl_B2 <- matrix(0,n_samB2,d2)
  for(i in 1:d2){
    gg_B2 <- (dataB2 - mean22[i])/sqrt(sigma22[i])
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

  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("1MG-A",round(abc,4),round(AIC,4),round(t(mean11),4)," "," ",round(sigma11[1],4),round(t(mix_pi_1),4)," "," ",round(t(mean22),4)," "," ",round(sigma22[1],4),round(t(mix_pi_2),4),
                       " "," ",round(B1[1],4),round(B1[2],4),round(B1[3],4)," "," "," ",round(jj_1,4),round(ll_1*100,4),round(jj_2,4),round(ll_2*100,4),round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],
                       round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],
                       round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
}
####################(A-3)Model##########################################
BCModelFun[[4]] <- function(K1,logL,df11,df21){
  dataB1 <- as.matrix(as.numeric(df11[,1])); dataB2 <- as.matrix(as.numeric(df21[,1]))
  n_samB1 <- dim(dataB1)[1];n_samB2 <- dim(dataB2)[1]
  mean1<-mean(dataB1);mean2<-mean(dataB2)
  sigma1<- as.numeric(var(dataB1)); sigma2<- as.numeric(var(dataB2))
  m_esp<-0.0001
  ####################(A-3)Model##########################################
  d21<- 1
  mi_1 <- as.matrix(1)
  sigma11 <- as.matrix(sigma1/5)
  mean11 <- as.matrix(mean1)

  d22<- 2
  mi_2 <- as.matrix(c(0.5,0.5))
  sigma22 <- matrix((sigma2/5),d22,1)
  a2 <- sqrt(sigma2/n_samB2)
  mean22 <- as.matrix(c((mean2+2*a2),(mean2-2*a2)))

  L0 <- logL(n_samB1,d21,mi_1,mean11,sigma11,dataB1)+logL(n_samB2,d22,mi_2,mean22,sigma22,dataB2)
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean11[i],sqrt(sigma11[i]))/dmixnorm(dataB1,mean11,sqrt(sigma11),mi_1)}
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean22[i],sqrt(sigma22[i]))/dmixnorm(dataB2,mean22,sqrt(sigma22),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2
    ############ CM1-step for means ##############
    n0_1 <- n_samB1

    n0_2 <- n_samB2*mix_pi_2
    n0_2[n0_2<0.000001] <- 0.000001

    mean11<- sumwx_B1/n0_1
    mean22<- sumwx_B2/n0_2

    ##########obtain variance##########
    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean11[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean22[i])^2 }

    sigma11<- sum(swx_B1)/n_samB1
    sigma22<- matrix((sum(swx_B2)/n_samB2),d22,1)

    ########criteria for iterations to stop#######

    L1 <- logL(n_samB1,d21,mix_pi_1,mean11,sigma11,dataB1)+logL(n_samB2,d22,mix_pi_2,mean22,sigma22,dataB2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa< 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*6

  #########first order genetic parameter process##########
  aa<- matrix(c(1,0,0, 0,1,1, 1,1,-1),3,3)
  b_line1 <- matrix(c(mean11,mean22))
  B1 <- solve(aa,b_line1)

  jj_1 <- 0
  ll_1 <- jj_1/sigma1

  jj_2 <- sigma2 - sigma22[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigma2

  ######### hypothesis testing #########
  ###############  B1  ###################
  dataB1<-sort(dataB1);bmw_B1 <- matrix(0,n_samB1,1); bmwsl_B1 <- matrix(0,n_samB1,d21)
  for(i in 1:d21){
    gg_B1 <- (dataB1 - mean11[i])/sqrt(as.vector(sigma11[i]))
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
    gg_B2 <- (dataB2 - mean22[i])/sqrt(sigma22[i])
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

  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("1MG-EAD",round(abc,4),round(AIC,4),round(t(mean11),4)," "," "," ",round(sigma11[1],4),round(t(mix_pi_1),4)," "," "," ",round(t(mean22),4)," "," ",round(sigma22[1],4),round(t(mix_pi_2),4),
                       " "," ",round(B1[1],4),round(B1[2],4),round(B1[3],4)," ",round(B1[3],4)," ",round(jj_1,4),round(ll_1*100,4),round(jj_2,4),round(ll_2*100,4),round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],
                       round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],
                       round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
}
####################(A-4)Model##########################################
BCModelFun[[5]] <- function(K1,logL,df11,df21){
  dataB1 <- as.matrix(as.numeric(df11[,1])); dataB2 <- as.matrix(as.numeric(df21[,1]))
  n_samB1 <- dim(dataB1)[1];n_samB2 <- dim(dataB2)[1]
  mean1<-mean(dataB1);mean2<-mean(dataB2)
  sigma1<- as.numeric(var(dataB1)); sigma2<- as.numeric(var(dataB2))
  m_esp<-0.0001
  ####################(A-4)Model##########################################
  d21<-2
  mi_1 <- as.matrix(c(0.5,0.5))
  sigma11 <- matrix((sigma1/5),d21,1)
  a1 <- sqrt(sigma1/n_samB1)
  mean11 <- as.matrix(c((mean1+2*a1),(mean1-2*a1)))

  d22<-1
  mi_2 <- as.matrix(1)
  sigma22 <- as.matrix(sigma2/5)
  mean22 <- as.matrix(mean2)

  L0 <- logL(n_samB1,d21,mi_1,mean11,sigma11,dataB1)+logL(n_samB2,d22,mi_2,mean22,sigma22,dataB2)
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean11[i],sqrt(sigma11[i]))/dmixnorm(dataB1,mean11,sqrt(sigma11),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean22[i],sqrt(sigma22[i]))/dmixnorm(dataB2,mean22,sqrt(sigma22),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2
    ############ CM1-step for means ##############
    n0_1 <- n_samB1*mix_pi_1
    n0_1[n0_1<0.000001] <- 0.000001

    n0_2 <- n_samB2

    mean11<- sumwx_B1/n0_1
    mean22<- sumwx_B2/n_samB2
    ##########obtain variance##########
    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean11[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean22[1])^2 }

    sigma11<- matrix((sum(swx_B1)/n_samB1),d21,1)
    sigma22<- sum(swx_B2)/n_samB2
    ########criteria for iterations to stop#######

    L1 <- logL(n_samB1,d21,mix_pi_1,mean11,sigma11,dataB1)+logL(n_samB2,d22,mix_pi_2,mean22,sigma22,dataB2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*6
  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,0, 0,0,1, 1,-1,-1),3,3)
  b_line1 <- matrix(c(mean11,mean22))
  B1 <- solve(aa,b_line1)

  jj_1 <- sigma1 - sigma11[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigma1

  jj_2 <- 0
  ll_2 <- jj_2/sigma2

  ######### hypothesis testing #########
  ###############  B1  ###################
  dataB1<-sort(dataB1);bmw_B1 <- matrix(0,n_samB1,1); bmwsl_B1 <- matrix(0,n_samB1,d21)
  for(i in 1:d21){
    gg_B1 <- (dataB1 - mean11[i])/sqrt(sigma11[i])
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
    gg_B2 <- (dataB2 - mean22[i])/sqrt(as.vector(sigma22[i]))
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

  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("1MG-NCD",round(abc,4),round(AIC,4),round(t(mean11),4)," "," ",round(sigma11[1],4),round(t(mix_pi_1),4)," "," ",round(t(mean22),4)," "," "," ",round(sigma22[1],4),round(t(mix_pi_2),4),
                       " "," "," ",round(B1[1],4),round(B1[2],4),round(B1[3],4)," ",round(-B1[3],4)," ",round(jj_1,4),round(ll_1*100,4),round(jj_2,4),round(ll_2*100,4),round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],
                       round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],
                       round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
}
####################(B-1)Model##########################################
BCModelFun[[6]] <- function(K1,logL,df11,df21){
  dataB1 <- as.matrix(as.numeric(df11[,1])); dataB2 <- as.matrix(as.numeric(df21[,1]))
  n_samB1 <- dim(dataB1)[1];n_samB2 <- dim(dataB2)[1]
  mean1<-mean(dataB1);mean2<-mean(dataB2)
  sigma1<- as.numeric(var(dataB1)); sigma2<- as.numeric(var(dataB2))
  m_esp<-0.0001
  ####################(B-1)Model##########################################
  d2<-4
  mi_1 <- as.matrix(c(0.25,0.25,0.25,0.25))
  sigma11 <- matrix((sigma1/2),d2,1)
  a1 <- sqrt(sigma1/n_samB1)
  mean11 <- as.matrix(c((mean1+2.4*a1),(mean1+0.8*a1),(mean1-0.8*a1),(mean1-2.4*a1)))

  mi_2 <- as.matrix(c(0.25,0.25,0.25,0.25))
  sigma22 <- matrix((sigma2/5),d2,1)
  a2 <- sqrt(sigma2/n_samB2)
  mean22 <- as.matrix(c((mean2+2.4*a2),(mean2+0.8*a2),(mean2-0.8*a2),(mean2-2.4*a2)))

  L0 <- logL(n_samB1,d2,mi_1,mean11,sigma11,dataB1)+logL(n_samB2,d2,mi_2,mean22,sigma22,dataB2)
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d2,n_samB1); swx_B1 <- matrix(0,d2,1)
  WW_B2 <- matrix(0,d2,n_samB2); swx_B2 <- matrix(0,d2,1)
  rr <- matrix(0,2,1)
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d2) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean11[i],sqrt(sigma11[i]))/dmixnorm(dataB1,mean11,sqrt(sigma11),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d2) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean22[i],sqrt(sigma22[i]))/dmixnorm(dataB2,mean22,sqrt(sigma22),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2
    ############ CM1-step for means ##############
    n0_1 <- n_samB1*mix_pi_1
    n0_1[n0_1<0.000001] <- 0.000001

    n0_2 <- n_samB2*mix_pi_2
    n0_2[n0_2<0.000001] <- 0.000001

    aa1<- sumwx_B1[1]/n0_1[1]-sumwx_B1[2]/n0_1[2]-sumwx_B1[3]/n0_1[3]+sumwx_B1[4]/n0_1[4]
    aa2<- sumwx_B2[1]/n0_2[1]-sumwx_B2[2]/n0_2[2]-sumwx_B2[3]/n0_2[3]+sumwx_B2[4]/n0_2[4]
    aa3<- sum(sigma11/n0_1);aa4<- sum(sigma22/n0_2)

    rr[1]<- aa1/aa3;rr[2]<- aa2/aa4

    mean11[1]<- (sumwx_B1[1]-sigma11[1]*rr[1])/n0_1[1]
    mean11[2]<- (sumwx_B1[2]+sigma11[2]*rr[1])/n0_1[2]
    mean11[3]<- (sumwx_B1[3]+sigma11[3]*rr[1])/n0_1[3]
    mean11[4]<- (sumwx_B1[4]-sigma11[4]*rr[1])/n0_1[4]
    mean22[1]<- (sumwx_B2[1]-sigma22[1]*rr[2])/n0_2[1]
    mean22[2]<- (sumwx_B2[2]+sigma22[2]*rr[2])/n0_2[2]
    mean22[3]<- (sumwx_B2[3]+sigma22[3]*rr[2])/n0_2[3]
    mean22[4]<- (sumwx_B2[4]-sigma22[4]*rr[2])/n0_2[4]

    ##########obtain variance##########
    for(i in 1:d2) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean11[i])^2 }
    for(i in 1:d2) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean22[i])^2 }

    sigma11<- matrix((sum(swx_B1)/n_samB1),d2,1)
    sigma22<- matrix((sum(swx_B2)/n_samB2),d2,1)

    ########criteria for iterations to stop#######

    L1 <- logL(n_samB1,d2,mix_pi_1,mean11,sigma11,dataB1)+logL(n_samB2,d2,mix_pi_2,mean22,sigma22,dataB2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*8

  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,0,0,0,0, 0,0,0,0,1,1,1,1, 1,1,0,0,0,0,-1,-1, 1,0,1,0,0,-1,0,-1, 0,0,1,1,1,1,0,0 ,0,1,0,1,1,0,1,0),8,6)
  b_line1 <- matrix(c(mean11,mean22))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigma1 - sigma11[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigma1

  jj_2 <- sigma2 - sigma22[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigma2

  ######### hypothesis testing #########
  ###############  B1  ###################
  dataB1<-sort(dataB1);bmw_B1 <- matrix(0,n_samB1,1); bmwsl_B1 <- matrix(0,n_samB1,d2)
  for(i in 1:d2){
    gg_B1 <- (dataB1 - mean11[i])/sqrt(sigma11[i])
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
  dataB2<-sort(dataB2);bmw_B2 <- matrix(0,n_samB2,1); bmwsl_B2 <- matrix(0,n_samB2,d2)
  for(i in 1:d2){
    gg_B2 <- (dataB2 - mean22[i])/sqrt(sigma22[i])
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

  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("2MG-AD",round(abc,4),round(AIC,4),round(t(mean11),4),round(sigma11[1],4),round(t(mix_pi_1),4),round(t(mean22),4),round(sigma22[1],4),round(t(mix_pi_2),4),
                       round(B1[1],4),round(B1[2],4),round(B1[3],4),round(B1[4],4),round(B1[5],4),round(B1[6],4),round(jj_1,4),round(ll_1*100,4),round(jj_2,4),round(ll_2*100,4),round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],
                       round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],
                       round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
}
####################(B-2)Model##########################################
BCModelFun[[7]] <- function(K1,logL,df11,df21){
  dataB1 <- as.matrix(as.numeric(df11[,1])); dataB2 <- as.matrix(as.numeric(df21[,1]))
  n_samB1 <- dim(dataB1)[1];n_samB2 <- dim(dataB2)[1]
  mean1<-mean(dataB1);mean2<-mean(dataB2)
  sigma1<- as.numeric(var(dataB1)); sigma2<- as.numeric(var(dataB2))
  m_esp<-0.0001
  ####################(B-2)Model##########################################
  d2<-4
  mi_1 <- as.matrix(c(0.25,0.25,0.25,0.25))
  sigma11 <- matrix((sigma1/2),d2,1)
  a1 <- sqrt(sigma1/n_samB1)
  mean11 <- as.matrix(c((mean1+2.4*a1),(mean1+1.1*a1),(mean1+0.5*a1),(mean1-2*a1)))

  mi_2 <- as.matrix(c(0.25,0.25,0.25,0.25))
  sigma22 <- matrix((sigma2/5),d2,1)
  a2 <- sqrt(sigma2/n_samB2)
  mean22 <- as.matrix(c((mean2+2.4*a2),(mean2+1.1*a2),(mean2+0.5*a2),(mean2-2*a2)))

  L0 <- logL(n_samB1,d2,mi_1,mean11,sigma11,dataB1)+logL(n_samB2,d2,mi_2,mean22,sigma22,dataB2)
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d2,n_samB1); swx_B1 <- matrix(0,d2,1)
  WW_B2 <- matrix(0,d2,n_samB2); swx_B2 <- matrix(0,d2,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d2) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean11[i],sqrt(sigma11[i]))/dmixnorm(dataB1,mean11,sqrt(sigma11),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d2) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean22[i],sqrt(sigma22[i]))/dmixnorm(dataB2,mean22,sqrt(sigma22),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2
    ############ CM1-step for means ##############
    n0_1 <- n_samB1*mix_pi_1
    n0_1[n0_1<0.000001] <- 0.000001

    n0_2 <- n_samB2*mix_pi_2
    n0_2[n0_2<0.000001] <- 0.000001

    ############solve the linear equation############
    hh<-matrix(0,4,4)
    hh[1,1]<- sigma11[1]/n0_1[1]+sigma11[2]/n0_1[2]+sigma11[3]/n0_1[3]+sigma11[4]/n0_1[4]
    hh[1,2]<- 0
    hh[1,3]<- sigma11[1]/n0_1[1]-sigma11[4]/n0_1[4]
    hh[1,4]<- -sigma11[2]/n0_1[2]+sigma11[3]/n0_1[3]

    hh[2,2]<- sigma22[1]/n0_2[1]+sigma22[2]/n0_2[2]+sigma22[3]/n0_2[3]+sigma22[4]/n0_2[4]
    hh[2,3]<- -sigma22[1]/n0_2[1]+sigma22[4]/n0_2[4]
    hh[2,4]<- sigma22[2]/n0_2[2]-sigma22[3]/n0_2[3]

    hh[3,3]<- sigma11[1]/n0_1[1]+sigma11[4]/n0_1[4]+sigma22[1]/n0_2[1]+sigma22[4]/n0_2[4]
    hh[3,4]<- 0

    hh[4,4]<- sigma11[2]/n0_1[2]+sigma11[3]/n0_1[3]+sigma22[2]/n0_2[2]+sigma22[3]/n0_2[3]
    for(i in 2:4){
      for(j in 1:(i-1)){
        hh[i,j]<- hh[j,i]
      }
    }
    ##################################################
    b_line<-matrix(0,4,1)
    b_line[1]<- sumwx_B1[1]/n0_1[1]-sumwx_B1[2]/n0_1[2]-sumwx_B1[3]/n0_1[3]+sumwx_B1[4]/n0_1[4]
    b_line[2]<- sumwx_B2[1]/n0_2[1]-sumwx_B2[2]/n0_2[2]-sumwx_B2[3]/n0_2[3]+sumwx_B2[4]/n0_2[4]
    b_line[3]<- sumwx_B1[1]/n0_1[1]-sumwx_B1[4]/n0_1[4]-sumwx_B2[1]/n0_2[1]+sumwx_B2[4]/n0_2[4]
    b_line[4]<- sumwx_B1[2]/n0_1[2]-sumwx_B1[3]/n0_1[3]-sumwx_B2[2]/n0_2[2]+sumwx_B2[3]/n0_2[3]

    B <- solve(hh,b_line)
    ##################################################
    mean11[1]<- (sumwx_B1[1]-sigma11[1]*(B[1]+B[3]))/n0_1[1]
    mean11[2]<- (sumwx_B1[2]+sigma11[2]*(B[1]-B[4]))/n0_1[2]
    mean11[3]<- (sumwx_B1[3]+sigma11[3]*(B[1]+B[4]))/n0_1[3]
    mean11[4]<- (sumwx_B1[4]-sigma11[4]*(B[1]-B[3]))/n0_1[4]
    mean22[1]<- (sumwx_B2[1]-sigma22[1]*(B[2]-B[3]))/n0_2[1]
    mean22[2]<- (sumwx_B2[2]+sigma22[2]*(B[2]+B[4]))/n0_2[2]
    mean22[3]<- (sumwx_B2[3]+sigma22[3]*(B[2]-B[4]))/n0_2[3]
    mean22[4]<- (sumwx_B2[4]-sigma22[4]*(B[2]+B[3]))/n0_2[4]

    ##########obtain variance##########
    for(i in 1:d2) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean11[i])^2 }
    for(i in 1:d2) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean22[i])^2 }

    sigma11<- matrix((sum(swx_B1)/n_samB1),d2,1)
    sigma22<- matrix((sum(swx_B2)/n_samB2),d2,1)

    ########criteria for iterations to stop#######
    L1 <- logL(n_samB1,d2,mix_pi_1,mean11,sigma11,dataB1)+logL(n_samB2,d2,mix_pi_2,mean22,sigma22,dataB2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*6
  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,1,0,0,0,0, 0,0,0,0,1,1,1,1, 1,1,0,0,0,0,-1,-1, 1,0,1,0,0,-1,0,-1),8,4)
  b_line1 <- matrix(c(mean11,mean22))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigma1 - sigma11[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigma1

  jj_2 <- sigma2 - sigma22[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigma2

  ######### hypothesis testing #########
  ###############  B1  ###################
  dataB1<-sort(dataB1);bmw_B1 <- matrix(0,n_samB1,1); bmwsl_B1 <- matrix(0,n_samB1,d2)
  for(i in 1:d2){
    gg_B1 <- (dataB1 - mean11[i])/sqrt(sigma11[i])
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
  dataB2<-sort(dataB2);bmw_B2 <- matrix(0,n_samB2,1); bmwsl_B2 <- matrix(0,n_samB2,d2)
  for(i in 1:d2){
    gg_B2 <- (dataB2 - mean22[i])/sqrt(sigma22[i])
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

  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("2MG-A",round(abc,4),round(AIC,4),round(t(mean11),4),round(sigma11[1],4),round(t(mix_pi_1),4),round(t(mean22),4),round(sigma22[1],4),round(t(mix_pi_2),4),
                       round(B1[1],4),round(B1[2],4),round(B1[3],4),round(B1[4],4)," "," ",round(jj_1,4),round(ll_1*100,4),round(jj_2,4),round(ll_2*100,4),round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],
                       round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],
                       round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
}
####################(B-3)Model##########################################
BCModelFun[[8]] <- function(K1,logL,df11,df21){
  dataB1 <- as.matrix(as.numeric(df11[,1])); dataB2 <- as.matrix(as.numeric(df21[,1]))
  n_samB1 <- dim(dataB1)[1];n_samB2 <- dim(dataB2)[1]
  mean1<-mean(dataB1);mean2<-mean(dataB2)
  sigma1<- as.numeric(var(dataB1)); sigma2<- as.numeric(var(dataB2))
  m_esp<-0.0001
  ####################(B-3)Model##########################################
  d2<-3
  mi_1 <- as.matrix(c(0.25,0.5,0.25))
  sigma11 <- matrix((sigma1/2),d2,1)
  a1 <- sqrt(sigma1/n_samB1)
  mean11 <- as.matrix(c((mean1+2.4*a1),mean1,(mean1-2.4*a1)))

  mi_2 <- as.matrix(c(0.25,0.5,0.25))
  sigma22 <- matrix((sigma2/5),d2,1)
  a2 <- sqrt(sigma2/n_samB2)
  mean22 <- as.matrix(c((mean2+2.4*a2),mean2,(mean2-2.4*a2)))

  L0 <- logL(n_samB1,d2,mi_1,mean11,sigma11,dataB1)+logL(n_samB2,d2,mi_2,mean22,sigma22,dataB2)
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d2,n_samB1); swx_B1 <- matrix(0,d2,1)
  WW_B2 <- matrix(0,d2,n_samB2); swx_B2 <- matrix(0,d2,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d2) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean11[i],sqrt(sigma11[i]))/dmixnorm(dataB1,mean11,sqrt(sigma11),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d2) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean22[i],sqrt(sigma22[i]))/dmixnorm(dataB2,mean22,sqrt(sigma22),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2
    ############ CM1-step for means ##############
    n0_1 <- n_samB1*mix_pi_1
    n0_1[n0_1<0.000001] <- 0.000001

    n0_2 <- n_samB2*mix_pi_2
    n0_2[n0_2<0.000001] <- 0.000001

    ############solve the linear equation############
    hh<-matrix(0,3,3)
    hh[1,1]<- sigma11[1]/n0_1[1]+4*sigma11[2]/n0_1[2]+sigma11[3]/n0_1[3]
    hh[1,2]<- 0
    hh[1,3]<- sigma11[1]/n0_1[1]-sigma11[3]/n0_1[3]

    hh[2,2]<- sigma22[1]/n0_2[1]+4*sigma22[2]/n0_2[2]+sigma22[3]/n0_2[3]
    hh[2,3]<- -sigma22[1]/n0_2[1]+sigma22[3]/n0_2[3]

    hh[3,3]<- sigma11[1]/n0_1[1]+sigma11[3]/n0_1[3]+sigma22[1]/n0_2[1]+sigma22[3]/n0_2[3]
    for(i in 2:3){
      for(j in 1:(i-1)){
        hh[i,j]<- hh[j,i]
      }
    }
    ##################################################
    b_line<-matrix(0,3,1)
    b_line[1]<- sumwx_B1[1]/n0_1[1]-2*sumwx_B1[2]/n0_1[2]+sumwx_B1[3]/n0_1[3]
    b_line[2]<- sumwx_B2[1]/n0_2[1]-2*sumwx_B2[2]/n0_2[2]+sumwx_B2[3]/n0_2[3]
    b_line[3]<- sumwx_B1[1]/n0_1[1]-sumwx_B1[3]/n0_1[3]-sumwx_B2[1]/n0_2[1]+sumwx_B2[3]/n0_2[3]

    B <- solve(hh,b_line)
    ##################################################
    mean11[1]<- (sumwx_B1[1]-sigma11[1]*(B[1]+B[3]))/n0_1[1]
    mean11[2]<- (sumwx_B1[2]+sigma11[2]*2.0*B[1])/n0_1[2]
    mean11[3]<- (sumwx_B1[3]-sigma11[3]*(B[1]-B[3]))/n0_1[3]
    mean22[1]<- (sumwx_B2[1]-sigma22[1]*(B[2]-B[3]))/n0_2[1]
    mean22[2]<- (sumwx_B2[2]+sigma22[2]*2.0*B[2])/n0_2[2]
    mean22[3]<- (sumwx_B2[3]-sigma22[3]*(B[2]+B[3]))/n0_2[3]
    ##########obtain variance##########
    for(i in 1:d2) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean11[i])^2 }
    for(i in 1:d2) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean22[i])^2 }

    sigma11<- matrix((sum(swx_B1)/n_samB1),d2,1)
    sigma22<- matrix((sum(swx_B2)/n_samB2),d2,1)
    ########criteria for iterations to stop#######
    L1 <- logL(n_samB1,d2,mix_pi_1,mean11,sigma11,dataB1)+logL(n_samB2,d2,mix_pi_2,mean22,sigma22,dataB2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*5
  #########first order genetic parameter process##########
  aa<- matrix(c(1,1,1,0,0,0, 0,0,0,1,1,1, 2,1,0,0,-1,-2),6,3)
  b_line1 <- matrix(c(mean11,mean22))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- sigma1 - sigma11[1]
  if(jj_1 < 0) {jj_1 <- 0}
  ll_1 <- jj_1/sigma1

  jj_2 <- sigma2 - sigma22[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigma2
  ######### hypothesis testing #########
  ###############  B1  ###################
  dataB1<-sort(dataB1);bmw_B1 <- matrix(0,n_samB1,1); bmwsl_B1 <- matrix(0,n_samB1,d2)
  for(i in 1:d2){
    gg_B1 <- (dataB1 - mean11[i])/sqrt(sigma11[i])
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
  dataB2<-sort(dataB2);bmw_B2 <- matrix(0,n_samB2,1); bmwsl_B2 <- matrix(0,n_samB2,d2)
  for(i in 1:d2){
    gg_B2 <- (dataB2 - mean22[i])/sqrt(sigma22[i])
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

  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("2MG-EA",round(abc,4),round(AIC,4),round(t(mean11),4)," ",round(sigma11[1],4),round(t(mix_pi_1),4)," ",round(t(mean22),4)," ",round(sigma22[1],4),round(t(mix_pi_2),4)," ",
                       round(B1[1],4),round(B1[2],4),round(B1[3],4),round(B1[3],4)," "," ",round(jj_1,4),round(ll_1*100,4),round(jj_2,4),round(ll_2*100,4),round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],
                       round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],
                       round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
}
####################(B-4)Model##########################################
BCModelFun[[9]] <- function(K1,logL,df11,df21){
  dataB1 <- as.matrix(as.numeric(df11[,1])); dataB2 <- as.matrix(as.numeric(df21[,1]))
  n_samB1 <- dim(dataB1)[1];n_samB2 <- dim(dataB2)[1]
  mean1<-mean(dataB1);mean2<-mean(dataB2)
  sigma1<- as.numeric(var(dataB1)); sigma2<- as.numeric(var(dataB2))
  m_esp<-0.0001
  ####################(B-4)Model##########################################
  d21<-1
  mi_1 <- as.matrix(1)
  sigma11 <- as.matrix(sigma1/2)
  mean11 <- as.matrix(mean1)

  d22<-4
  mi_2 <- as.matrix(c(0.25,0.25,0.25,0.25))
  sigma22 <- matrix((sigma2/5),d22,1)
  a2 <- sqrt(sigma2/n_samB2)
  mean22 <- as.matrix(c((mean2+2.4*a2),(mean2+1.1*a2),(mean2+0.5*a2),(mean2-2.4*a2)))

  L0 <- logL(n_samB1,d21,mi_1,mean11,sigma11,dataB1)+logL(n_samB2,d22,mi_2,mean22,sigma22,dataB2)
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)
  rr <- matrix(0,2,1)
  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean11[i],sqrt(sigma11[i]))/dmixnorm(dataB1,mean11,sqrt(sigma11),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean22[i],sqrt(sigma22[i]))/dmixnorm(dataB2,mean22,sqrt(sigma22),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2
    ############ CM1-step for means ##############
    n0_1 <- n_samB1*mix_pi_1
    n0_1[n0_1<0.000001] <- 0.000001

    n0_2 <- n_samB2*mix_pi_2
    n0_2[n0_2<0.000001] <- 0.000001

    aa1<- sumwx_B2[1]/n0_2[1]-sumwx_B2[2]/n0_2[2]-sumwx_B2[3]/n0_2[3]+sumwx_B2[4]/n0_2[4]
    aa2<- sigma22[1]/n0_2[1]+sigma22[2]/n0_2[2]+sigma22[3]/n0_2[3]+sigma22[4]/n0_2[4]
    aa3<- aa1/aa2

    ##################################################
    mean11[1]<- sumwx_B1[1]/n_samB1

    mean22[1]<- (sumwx_B2[1]-sigma22[1]*aa3)/n0_2[1]
    mean22[2]<- (sumwx_B2[2]+sigma22[2]*aa3)/n0_2[2]
    mean22[3]<- (sumwx_B2[3]+sigma22[3]*aa3)/n0_2[3]
    mean22[4]<- (sumwx_B2[4]-sigma22[4]*aa3)/n0_2[4]
    ##########obtain variance##########
    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean11[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean22[i])^2 }

    sigma11<- sum(swx_B1)/n_samB1
    sigma22<- matrix((sum(swx_B2)/n_samB2),d22,1)

    ########criteria for iterations to stop#######

    L1 <- logL(n_samB1,d21,mix_pi_1,mean11,sigma11,dataB1)+logL(n_samB2,d22,mix_pi_2,mean22,sigma22,dataB2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*6
  #########first order genetic parameter process##########
  aa<- matrix(c(1,0,0,0,0, 0,1,1,1,1, 1,1,1,-1,-1, 1,1,-1,1,-1),5,4)
  b_line1 <- matrix(c(mean11,mean22))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- 0
  ll_1 <- jj_1/sigma1

  jj_2 <- sigma2 - sigma22[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigma2

  ######### hypothesis testing #########
  ###############  B1  ###################
  dataB1<-sort(dataB1);bmw_B1 <- matrix(0,n_samB1,1); bmwsl_B1 <- matrix(0,n_samB1,d21)
  for(i in 1:d21){
    gg_B1 <- (dataB1 - mean11[i])/sqrt(as.vector(sigma11[i]))
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
    gg_B2 <- (dataB2 - mean22[i])/sqrt(sigma22[i])
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

  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("2MG-CD",round(abc,4),round(AIC,4),round(t(mean11),4)," "," "," ",round(sigma11[1],4),round(t(mix_pi_1),4)," "," "," ",round(t(mean22),4),round(sigma22[1],4),round(t(mix_pi_2),4),
                       round(B1[1],4),round(B1[2],4),round(B1[3],4),round(B1[4],4),round(B1[3],4),round(B1[4],4),round(jj_1,4),round(ll_1*100,4),round(jj_2,4),round(ll_2*100,4),round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],
                       round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],
                       round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
}
####################(B-5)Model##########################################
BCModelFun[[10]] <- function(K1,logL,df11,df21){
  dataB1 <- as.matrix(as.numeric(df11[,1])); dataB2 <- as.matrix(as.numeric(df21[,1]))
  n_samB1 <- dim(dataB1)[1];n_samB2 <- dim(dataB2)[1]
  mean1<-mean(dataB1);mean2<-mean(dataB2)
  sigma1<- as.numeric(var(dataB1)); sigma2<- as.numeric(var(dataB2))
  m_esp<-0.0001
  ####################(B-5)Model##########################################
  d21<-1
  mi_1 <- as.matrix(1)
  sigma11 <- as.matrix(sigma1/2)
  mean11 <- as.matrix(mean1)

  d22<-3
  mi_2 <- as.matrix(c(0.25,0.5,0.25))
  sigma22 <- matrix((sigma2/5),d22,1)
  a2 <- sqrt(sigma2/n_samB2)
  mean22 <- as.matrix(c((mean2+2.4*a2),mean2,(mean2-2.4*a2)))

  L0 <- logL(n_samB1,d21,mi_1,mean11,sigma11,dataB1)+logL(n_samB2,d22,mi_2,mean22,sigma22,dataB2)
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW_B1 <- matrix(0,d21,n_samB1); swx_B1 <- matrix(0,d21,1)
  WW_B2 <- matrix(0,d22,n_samB2); swx_B2 <- matrix(0,d22,1)

  while(stopa > m_esp && iteration<=1000){
    iteration <- iteration + 1
    ############ E-step #############
    for(i in 1:d21) { WW_B1[i,] <- mi_1[i]*dnorm(dataB1,mean11[i],sqrt(sigma11[i]))/dmixnorm(dataB1,mean11,sqrt(sigma11),mi_1) }
    mix_pi_1 <- as.matrix(rowSums(WW_B1)/n_samB1)
    sumwx_B1 <- WW_B1%*%dataB1

    for(i in 1:d22) { WW_B2[i,] <- mi_2[i]*dnorm(dataB2,mean22[i],sqrt(sigma22[i]))/dmixnorm(dataB2,mean22,sqrt(sigma22),mi_2) }
    mix_pi_2 <- as.matrix(rowSums(WW_B2)/n_samB2)
    sumwx_B2 <- WW_B2%*%dataB2
    ############ CM1-step for means ##############
    n0_1 <- n_samB1*mix_pi_1
    n0_1[n0_1<0.000001] <- 0.000001

    n0_2 <- n_samB2*mix_pi_2
    n0_2[n0_2<0.000001] <- 0.000001

    aa1<- sumwx_B2[1]/n0_2[1]-2.0*sumwx_B2[2]/n0_2[2]+sumwx_B2[3]/n0_2[3]
    aa2<- sigma22[1]/n0_2[1]+4.0*sigma22[2]/n0_2[2]+sigma22[3]/n0_2[3]
    aa3<- aa1/aa2
    mean11[1]<- sumwx_B1[1]/n_samB1

    mean22[1]<- (sumwx_B2[1]-sigma22[1]*aa3)/n0_2[1]
    mean22[2]<- (sumwx_B2[2]+sigma22[2]*2.0*aa3)/n0_2[2]
    mean22[3]<- (sumwx_B2[3]-sigma22[3]*aa3)/n0_2[3]
    ##########obtain variance##########
    for(i in 1:d21) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean11[i])^2 }
    for(i in 1:d22) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean22[i])^2 }

    sigma11<- sum(swx_B1)/n_samB1
    sigma22<- matrix((sum(swx_B2)/n_samB2),d22,1)

    ########criteria for iterations to stop#######
    L1 <- logL(n_samB1,d21,mix_pi_1,mean11,sigma11,dataB1)+logL(n_samB2,d22,mix_pi_2,mean22,sigma22,dataB2)

    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }

  abc <- L1
  AIC <- -2*abc + 2*5
  #########first order genetic parameter process##########
  aa<- matrix(c(1,0,0,0, 0,1,1,1, 2,2,0,-2),4,3)
  b_line1 <- matrix(c(mean11,mean22))
  B1 <- solve(t(aa)%*%aa)%*%(t(aa)%*%b_line1)

  jj_1 <- 0
  ll_1 <- jj_1/sigma1

  jj_2 <- sigma2 - sigma22[1]
  if(jj_2 < 0) {jj_2 <- 0}
  ll_2 <- jj_2/sigma2

  ######### hypothesis testing #########
  ###############  B1  ###################
  dataB1<-sort(dataB1);bmw_B1 <- matrix(0,n_samB1,1); bmwsl_B1 <- matrix(0,n_samB1,d21)
  for(i in 1:d21){
    gg_B1 <- (dataB1 - mean11[i])/sqrt(as.vector(sigma11[i]))
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
    gg_B2 <- (dataB2 - mean22[i])/sqrt(sigma22[i])
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

  tt_B1[which(tt_B1>=10e-4)]<-round(tt_B1[which(tt_B1>=10e-4)],4);tt_B1[which(tt_B1<10e-4)]<-format(tt_B1[which(tt_B1<10e-4)],scientific=TRUE,digit=4)
  tt_B2[which(tt_B2>=10e-4)]<-round(tt_B2[which(tt_B2>=10e-4)],4);tt_B2[which(tt_B2<10e-4)]<-format(tt_B2[which(tt_B2<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("2MG-EAD",round(abc,4),round(AIC,4),round(t(mean11),4)," "," "," ",round(sigma11[1],4),round(t(mix_pi_1),4)," "," "," ",round(t(mean22),4)," ",round(sigma22[1],4),round(t(mix_pi_2),4)," ",
                       round(B1[1],4),round(B1[2],4),round(B1[3],4),round(B1[3],4),round(B1[3],4),round(B1[3],4),round(jj_1,4),round(ll_1*100,4),round(jj_2,4),round(ll_2*100,4),round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],
                       round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],
                       round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
}

K1BC <- function(x){
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

logLBC <- function(nm,nng,mi,mn,s,d1) { sum2 <- sum(log(dmixnorm(d1,mn,sqrt(s),mi)));return (sum2) }


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
  allresult=foreach(i=1:10,.combine = 'rbind')%dopar%{
    requireNamespace("KScorrect")
    requireNamespace("kolmim")
    BCModelFun[[i]](K1BC,logLBC,df11,df21)[[1]]
  }
  stopCluster(cl)
  mi_1<-NULL;mi_2<-NULL
}else{

allresultq=switch(model,"0MG"=BCModelFun[[1]](K1BC,logLBC,df11,df21),"1MG-AD" = BCModelFun[[2]](K1BC,logLBC,df11,df21),"1MG-A"=BCModelFun[[3]](K1BC,logLBC,df11,df21),"1MG-EAD"=BCModelFun[[4]](K1BC,logLBC,df11,df21),"1MG-NCD"=BCModelFun[[5]](K1BC,logLBC,df11,df21),
                        "2MG-AD"=BCModelFun[[6]](K1BC,logLBC,df11,df21),"2MG-A"=BCModelFun[[7]](K1BC,logLBC,df11,df21),"2MG-EA"=BCModelFun[[8]](K1BC,logLBC,df11,df21),"2MG-CD"=BCModelFun[[9]](K1BC,logLBC,df11,df21),"2MG-EAD"=BCModelFun[[10]](K1BC,logLBC,df11,df21))

allresult<-allresultq[[1]]
if(model!="0MG"){
  mi_1<-allresultq[[2]];mi_2<-allresultq[[3]]
}else{
  mi_1<-NULL;mi_2<-NULL
}
}
colnames(allresult) <- BCcolname

out<-list(allresult,mi_1,mi_2)
return(out)
}


