
BCFFun<-function(df,model,BCFtext2){
  
data<-sapply(df,as.character)
dB12<-data[-1,which(data[1,]=="B12")];B12<-as.numeric(dB12[which(is.na(as.numeric(dB12))==FALSE)]);df11<-as.data.frame(B12)
dB22<-data[-1,which(data[1,]=="B22")];B22<-as.numeric(dB22[which(is.na(as.numeric(dB22))==FALSE)]);df21<-as.data.frame(B22)
##################################################
BCFcolname <- c("Model","Log_Max_likelihood_value","AIC","mean(B1:2)[1]","mean(B1:2)[2]","mean(B1:2)[3]","mean(B1:2)[4]","Var(B1:2)[1]",
                "Var(B1:2)[2]","Var(B1:2)[3]","Var(B1:2)[4]","Proportion(B1:2)[1]","Proportion(B1:2)[2]","Proportion(B1:2)[3]","Proportion(B1:2)[4]","mean(B2:2)[1]",
                "mean(B2:2)[2]","mean(B2:2)[3]","mean(B2:2)[4]","Var(B2:2)[1]","Var(B2:2)[2]","Var(B2:2)[3]","Var(B2:2)[4]","Proportion(B2:2)[1]",
                "Proportion(B2:2)[2]","Proportion(B2:2)[3]","Proportion(B2:2)[4]","m1","m2","da(d)","db","ha(h)","hb",
                "Major-Gene Var(B1:2)","Heritability(Major-Gene(B1:2))(%)","Major-Gene Var(B2:2)","Heritability(Major-Gene(B2:2))(%)","U1 square(B1:2)","p(U1 square(B1:2))",
                "U2 square(B1:2)","p(U2 square(B1:2))","U3 square(B1:2)","p(U3 square(B1:2))","nW square(B1:2)","p(nW square(B1:2))","Dn(B1:2)","p(Dn(B1:2))","U1 square(B2:2)","p(U1 square(B2:2))","U2 square(B2:2)","p(U2 square(B2:2))",
                "U3 square(B2:2)","p(U3 square(B2:2))","nW square(B2:2)","p(nW square(B2:2))","Dn(B2:2)","p(Dn(B2:2))")

BCFModelFun<-list(NA)
###################define each model function##################
####################(A-0)Model##########################################
BCFModelFun[[1]] <- function(K1,logL,df11,df21,BCFtext2){
  dataB1 <- as.matrix(as.numeric(df11[,1])); dataB2 <- as.matrix(as.numeric(df21[,1]))
  n_samB1 <- dim(dataB1)[1];n_samB2 <- dim(dataB2)[1]
  sumx1<- sum(dataB1);mean1<- mean(dataB1);sigma1<- as.numeric(var(dataB1))
  sumx2<- sum(dataB2);mean2<- mean(dataB2);sigma2<- as.numeric(var(dataB2)) 
  m_esp<-0.0001 ;num_l<- as.numeric(BCFtext2)
  ####################0MG Model##########################################
  d2<- 1
  mean11<- mean1;mean22<- mean2
  sigma11<- sigma1;sigma22<- sigma2
  mix_pi_1<- 1;mix_pi_2<- 1
  
  L0<- logL(n_samB1,d2,mix_pi_1,mean11,sigma11,dataB1)+logL(n_samB2,d2,mix_pi_2,mean22,sigma22,dataB2)  
  aBCF<- L0
  AIC<- -2.0*aBCF+2.0*2.0
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
  
  output <- data.frame("0MG",round(aBCF,4),round(AIC,4),round(t(mean11),4)," "," "," ",round(t(sigma11),4)," "," "," ",round(t(mix_pi_1),4)," "," "," ",
                       round(t(mean22),4)," "," "," ",round(t(sigma22),4)," "," "," ",round(t(mix_pi_2),4)," "," "," ",
                       " "," "," "," "," "," "," "," "," "," ",
                       round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],
                       round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output)
  return(OUTPUT)
  
}

####################(A-1)Model##########################################
BCFModelFun[[2]] <- function(K1,logL,df11,df21,BCFtext2){
  dataB1 <- as.matrix(as.numeric(df11[,1])); dataB2 <- as.matrix(as.numeric(df21[,1]))
  n_samB1 <- dim(dataB1)[1];n_samB2 <- dim(dataB2)[1]
  sumx1<- sum(dataB1);mean1<-mean(dataB1);sigma1<- as.numeric(var(dataB1))
  sumx2<- sum(dataB2);mean2<-mean(dataB2);sigma2<- as.numeric(var(dataB2))
  
  m_esp<-0.0001 ;num_l<- as.numeric(BCFtext2)
  ####################1MG-AD Model########  (A1) ##############################
  d2<- 2
  mi_1<- matrix(0.5,d2,1);mi_2<- matrix(0.5,d2,1)
  sigma11<- matrix(sigma1/5,d2,1) ;sigma22<- matrix(sigma2/5,d2,1)
  a1<-sqrt(sigma1/(n_samB1-1))
  mean11<- matrix(c((mean1+2*a1),(mean1-2*a1)))
  a2<-sqrt(sigma2/(n_samB2-1))
  mean22<- matrix(c((mean2+2*a2),(mean2-2*a2)))
  ############  first order genetic parameter  ############
  hh<- matrix(c(1,1,0,0, 0,0,1,1, 1,0,0,-1, 0,0.5,0.5,0),4,4)
  b_line <- matrix(c(mean11,mean22))
  B1 <- solve(hh,b_line) ;a1<- B1[3];a2<- B1[4]
  gg1<- (0.5*a1*a1+0.25*a2*a2)/num_l
  sigma11[2]<- sigma11[1]+gg1; sigma22[1]<- sigma22[2]+gg1
  ############   likelihood values of the initial value of calculation ##############
  L0<- logL(n_samB1,d2,mi_1,mean11,sigma11,dataB1)+logL(n_samB2,d2,mi_2,mean22,sigma22,dataB2)  
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
    #########first order genetic parameter process##########
    b_line <- matrix(c(mean11,mean22))
    B11 <- solve(hh,b_line) 
    a1<- B11[3];a2<- B11[4]
    gg1<- (0.5*a1*a1+0.25*a2*a2)/num_l
    sigma11[2]<- sigma11[1]+gg1; sigma22[1]<- sigma22[2]+gg1
    ###########obtain variance##########
    for(i in 1:d2) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean11[i])^2 }
    for(i in 1:d2) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean22[i])^2 }
    
    aaa0<- sigma11[1]; n_iter<- 0;aa2<-1000
    while (aa2>0.0001){
      n_iter<- n_iter+1
      aa1<- sigma11[1]/(sigma11[1]+gg1)
      sigma11[1]<- (swx_B1[1]+aa1*aa1*swx_B1[2])/(n0_1[1]+aa1*n0_1[2])
      aa2<- abs(sigma11[1]-aaa0)
      aaa0<- sigma11[1]
      if (n_iter>20) break
    }
    sigma11[2]<- sigma11[1]+gg1
    aaa0<- sigma22[2]; n_iter<- 0; aa2<- 1000
    
    while (aa2>0.0001){
      n_iter<- n_iter+1
      aa1<- sigma22[2]/(sigma22[2]+gg1)
      sigma22[2]<- (swx_B2[2]+aa1*aa1*swx_B2[1])/(n0_2[2]+aa1*n0_2[1])
      aa2<- abs(sigma22[2]-aaa0)
      aaa0<- sigma22[2]
      if (n_iter>20) break
    }
    sigma22[1]<- sigma22[2]+gg1
    ########criteria for iterations to stop#######
    L1 <- logL(n_samB1,d2,mix_pi_1,mean11,sigma11,dataB1)+logL(n_samB2,d2,mix_pi_2,mean22,sigma22,dataB2) 
    
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  
  aBCF <- L1
  AIC <- -2*aBCF + 2*6
  
  #########first order genetic parameter process##########
  b_line <- matrix(c(mean11,mean22))
  B111 <- solve(hh,b_line) 
  #########second order genetic parameter process##########
  jj_1<- sigma1-sigma11[1]
  if(jj_1<0){jj_1<- 0}
  ll_1<- jj_1/sigma1
  
  jj_2<-sigma2-sigma22[2]
  if(jj_2<0){jj_2<- 0}
  ll_2<- jj_2/sigma2
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
  
  output <- data.frame("1MG-AD",round(aBCF,4),round(AIC,4),round(t(mean11),4)," "," ",round(t(sigma11),4)," "," ",round(t(mix_pi_1),4)," "," ",
                       round(t(mean22),4)," "," ",round(t(sigma22),4)," "," ",round(t(mix_pi_2),4)," "," ",
                       round(B111[1],4),round(B111[2],4),round(B111[3],4)," ",round(B111[4],4)," ",round(jj_1,4),round(ll_1*100,4),round(jj_2,4),round(ll_2*100,4),
                       round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],
                       round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
} 

####################(A-2)Model##########################################
BCFModelFun[[3]] <- function(K1,logL,df11,df21,BCFtext2){
  dataB1 <- as.matrix(as.numeric(df11[,1])); dataB2 <- as.matrix(as.numeric(df21[,1]))
  n_samB1 <- dim(dataB1)[1];n_samB2 <- dim(dataB2)[1]
  sumx1<- sum(dataB1);mean1<-mean(dataB1);sigma1<- as.numeric(var(dataB1))
  sumx2<- sum(dataB2);mean2<-mean(dataB2);sigma2<- as.numeric(var(dataB2))
  
  m_esp<-0.0001 ;num_l<- as.numeric(BCFtext2)
  ####################1MG-A Model########  (A2) ##############################
  d2<- 2
  mi_1<- matrix(0.5,d2,1);mi_2<- matrix(0.5,d2,1)
  sigma11<- matrix(sigma1/5,d2,1) ;sigma22<- matrix(sigma2/5,d2,1)
  a1<-sqrt(sigma1/(n_samB1-1))
  mean11<- matrix(c((mean1+2*a1),(mean1-2*a1)))
  a2<-sqrt(sigma2/(n_samB2-1))
  mean22<- matrix(c((mean2+2*a2),(mean2-2*a2)))
  
  ############  first order genetic parameter  ############
  hh<- matrix(c(1,1,0,0, 0,0,1,1, 1,0,0,-1),4,3)
  b_line <- matrix(c(mean11,mean22))
  B1 <- solve(t(hh)%*%hh)%*%(t(hh)%*%b_line) ;a1<- B1[3]
  gg1<- 0.5*a1*a1/num_l
  sigma11[2]<- sigma11[1]+gg1; sigma22[1]<- sigma22[2]+gg1
  ############  likelihood values of the initial value of calculation ###################
  L0<- logL(n_samB1,d2,mi_1,mean11,sigma11,dataB1)+logL(n_samB2,d2,mi_2,mean22,sigma22,dataB2)  
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
    
    aa1<- sigma11[1]/n0_1[1]+sigma11[2]/n0_1[2]+sigma22[1]/n0_2[1]+sigma22[2]/n0_2[2]
    aa2<- sumwx_B1[1]/n0_1[1]-sumwx_B1[2]/n0_1[2]-sumwx_B2[1]/n0_2[1]+sumwx_B2[2]/n0_2[2]
    aa3<- aa2/aa1
    mean11[1]<- (sumwx_B1[1]-sigma11[1]*aa3)/n0_1[1]
    mean11[2]<- (sumwx_B1[2]+sigma11[2]*aa3)/n0_1[2]
    mean22[1]<- (sumwx_B2[1]+sigma22[1]*aa3)/n0_2[1]
    mean22[2]<- (sumwx_B2[2]-sigma22[2]*aa3)/n0_2[2]
    #########first order genetic parameter process##########
    b_line <- matrix(c(mean11,mean22))
    B11 <- solve(t(hh)%*%hh)%*%(t(hh)%*%b_line)
    a1<- B11[3]
    gg1<- 0.5*a1*a1/num_l
    sigma11[2]<- sigma11[1]+gg1; sigma22[1]<- sigma22[2]+gg1
    ###########obtain variance##########
    for(i in 1:d2) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean11[i])^2 }
    for(i in 1:d2) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean22[i])^2 }
    
    aaa0<- sigma11[1]; n_iter<- 0;aa2<-1000
    while (aa2>0.0001){
      n_iter<- n_iter+1
      aa1<- sigma11[1]/(sigma11[1]+gg1)
      sigma11[1]<- (swx_B1[1]+aa1*aa1*swx_B1[2])/(n0_1[1]+aa1*n0_1[2])
      aa2<- abs(sigma11[1]-aaa0)
      aaa0<- sigma11[1]
      if (n_iter>20) break
    }
    sigma11[2]<- sigma11[1]+gg1
    aaa0<- sigma22[2]; n_iter<- 0; aa2<- 1000
    while (aa2>0.0001){
      n_iter<- n_iter+1
      aa1<- sigma22[2]/(sigma22[2]+gg1)
      sigma22[2]<- (swx_B2[2]+aa1*aa1*swx_B2[1])/(n0_2[2]+aa1*n0_2[1])
      aa2<- abs(sigma22[2]-aaa0)
      aaa0<- sigma22[2]
      if (n_iter>20) break
    }
    sigma22[1]<- sigma22[2]+gg1
    
    ########criteria for iterations to stop#######
    L1 <- logL(n_samB1,d2,mix_pi_1,mean11,sigma11,dataB1)+logL(n_samB2,d2,mix_pi_2,mean22,sigma22,dataB2) 
    
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  
  aBCF <- L1
  AIC <- -2*aBCF + 2*5
  
  #########first order genetic parameter process##########
  b_line <- matrix(c(mean11,mean22))
  B111 <- solve(t(hh)%*%hh)%*%(t(hh)%*%b_line)
  #########second order genetic parameter process##########
  jj_1<- sigma1-sigma11[1]
  if(jj_1<0){jj_1<- 0}
  ll_1<- jj_1/sigma1
  
  jj_2<-sigma2-sigma22[2]
  if(jj_2<0){jj_2<- 0}
  ll_2<- jj_2/sigma2
  
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
  
  output <- data.frame("1MG-A",round(aBCF,4),round(AIC,4),round(t(mean11),4)," "," ",round(t(sigma11),4)," "," ",round(t(mix_pi_1),4)," "," ",
                       round(t(mean22),4)," "," ",round(t(sigma22),4)," "," ",round(t(mix_pi_2),4)," "," ",
                       round(B111[1],4),round(B111[2],4),round(B111[3],4)," "," "," ",round(jj_1,4),round(ll_1*100,4),round(jj_2,4),round(ll_2*100,4),
                       round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],
                       round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
} 

####################(A-3)Model##########################################
BCFModelFun[[4]] <- function(K1,logL,df11,df21,BCFtext2){
  dataB1 <- as.matrix(as.numeric(df11[,1])); dataB2 <- as.matrix(as.numeric(df21[,1]))
  n_samB1 <- dim(dataB1)[1];n_samB2 <- dim(dataB2)[1]
  sumx1<- sum(dataB1);mean1<-mean(dataB1);sigma1<- as.numeric(var(dataB1))
  sumx2<- sum(dataB2);mean2<-mean(dataB2);sigma2<- as.numeric(var(dataB2))
  
  m_esp<-0.0001 ;num_l<- as.numeric(BCFtext2)
  ####################1MG-EAD Model########  (A3) ##############################
  d2<- 2
  mi_1<- matrix(0.5,d2,1);mi_2<- matrix(0.5,d2,1)
  sigma11<- matrix(sigma1/5,d2,1) ;sigma22<- matrix(sigma2/5,d2,1)
  a1<-sqrt(sigma1/(n_samB1-1))
  mean11<- matrix(c((mean1+2*a1),(mean1-2*a1)))
  a2<-sqrt(sigma2/(n_samB2-1))
  mean22<- matrix(c((mean2+2*a2),(mean2-2*a2)))
  
  ############  first order genetic parameter  ############
  hh<- matrix(c(1,1,0,0, 0,0,1,1, 1,0.5,0.5,-1),4,3)
  b_line <- matrix(c(mean11,mean22))
  B1 <- solve(t(hh)%*%hh)%*%(t(hh)%*%b_line) ;a1<- B1[3]
  gg1<- 0.75*a1*a1/num_l
  sigma11[2]<- sigma11[1]+gg1; sigma22[1]<- sigma22[2]+gg1
  ############  likelihood values of the initial value of calculation ##########
  L0<- logL(n_samB1,d2,mi_1,mean11,sigma11,dataB1)+logL(n_samB2,d2,mi_2,mean22,sigma22,dataB2)  
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
    
    aa1<- 9.0*sigma11[1]/n0_1[1]+9.0*sigma11[2]/n0_1[2]+sigma22[1]/n0_2[1]+sigma22[2]/n0_2[2]
    aa2<- 3.0*sumwx_B1[1]/n0_1[1]-3.0*sumwx_B1[2]/n0_1[2]-sumwx_B2[1]/n0_2[1]+sumwx_B2[2]/n0_2[2]
    aa3<- aa2/aa1
    mean11[1]<- (sumwx_B1[1]-sigma11[1]*aa3*3)/n0_1[1]
    mean11[2]<- (sumwx_B1[2]+sigma11[2]*aa3*3)/n0_1[2]
    mean22[1]<- (sumwx_B2[1]+sigma22[1]*aa3)/n0_2[1]
    mean22[2]<- (sumwx_B2[2]-sigma22[2]*aa3)/n0_2[2]
    #########first order genetic parameter process##########
    b_line <- matrix(c(mean11,mean22))
    B11 <- solve(t(hh)%*%hh)%*%(t(hh)%*%b_line)
    a1<- B11[3]
    gg1<- 0.5*a1*a1/num_l
    sigma11[2]<- sigma11[1]+gg1; sigma22[1]<- sigma22[2]+gg1
    ###########obtain variance##########
    for(i in 1:d2) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean11[i])^2 }
    for(i in 1:d2) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean22[i])^2 }
    
    aaa0<- sigma11[1]; n_iter<- 0;aa2<-1000
    while (aa2>0.0001){
      n_iter<- n_iter+1
      aa1<- sigma11[1]/(sigma11[1]+gg1)
      sigma11[1]<- (swx_B1[1]+aa1*aa1*swx_B1[2])/(n0_1[1]+aa1*n0_1[2])
      aa2<- abs(sigma11[1]-aaa0)
      aaa0<- sigma11[1]
      if (n_iter>20) break
    }
    sigma11[2]<- sigma11[1]+gg1
    aaa0<- sigma22[2]; n_iter<- 0; aa2<- 1000
    while (aa2>0.0001){
      n_iter<- n_iter+1
      aa1<- sigma22[2]/(sigma22[2]+gg1)
      sigma22[2]<- (swx_B2[2]+aa1*aa1*swx_B2[1])/(n0_2[2]+aa1*n0_2[1])
      aa2<- abs(sigma22[2]-aaa0)
      aaa0<- sigma22[2]
      if (n_iter>20) break
    }
    sigma22[1]<- sigma22[2]+gg1
    
    ########criteria for iterations to stop#######
    
    L1 <- logL(n_samB1,d2,mix_pi_1,mean11,sigma11,dataB1)+logL(n_samB2,d2,mix_pi_2,mean22,sigma22,dataB2) 
    
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  
  aBCF <- L1
  AIC <- -2*aBCF + 2*5
  
  #########first order genetic parameter process##########
  b_line <- matrix(c(mean11,mean22))
  B111 <- solve(t(hh)%*%hh)%*%(t(hh)%*%b_line)
  
  #########second order genetic parameter process##########
  jj_1<- sigma1-sigma11[1]
  if(jj_1<0){jj_1<- 0}
  ll_1<- jj_1/sigma1
  
  jj_2<-sigma2-sigma22[2]
  if(jj_2<0){jj_2<- 0}
  ll_2<- jj_2/sigma2
  
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
  
  output <- data.frame("1MG-EAD",round(aBCF,4),round(AIC,4),round(t(mean11),4)," "," ",round(t(sigma11),4)," "," ",round(t(mix_pi_1),4)," "," ",
                       round(t(mean22),4)," "," ",round(t(sigma22),4)," "," ",round(t(mix_pi_2),4)," "," ",
                       round(B111[1],4),round(B111[2],4),round(B111[3],4)," "," "," ",round(jj_1,4),round(ll_1*100,4),round(jj_2,4),round(ll_2*100,4),
                       round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],
                       round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
} 

####################(A-4)Model##########################################
BCFModelFun[[5]] <- function(K1,logL,df11,df21,BCFtext2){
  dataB1 <- as.matrix(as.numeric(df11[,1])); dataB2 <- as.matrix(as.numeric(df21[,1]))
  n_samB1 <- dim(dataB1)[1];n_samB2 <- dim(dataB2)[1]
  sumx1<- sum(dataB1);mean1<-mean(dataB1);sigma1<- as.numeric(var(dataB1))
  sumx2<- sum(dataB2);mean2<-mean(dataB2);sigma2<- as.numeric(var(dataB2))
  
  m_esp<-0.0001 ;num_l<- as.numeric(BCFtext2)
  #################### 1MG-NCD Model########  (A4) ##############################
  d2<- 2
  mi_1<- matrix(0.5,d2,1);mi_2<- matrix(0.5,d2,1)
  sigma11<- matrix(sigma1/5,d2,1) ;sigma22<- matrix(sigma2/5,d2,1)
  a1<-sqrt(sigma1/(n_samB1-1))
  mean11<- matrix(c((mean1+2*a1),(mean1-2*a1)))
  a2<-sqrt(sigma2/(n_samB2-1))
  mean22<- matrix(c((mean2+a2),(mean2-a2)))
  
  ############  first order genetic parameter  ############
  hh<- matrix(c(1,1,0,0, 0,0,1,1, 1,-0.5,-0.5,-1),4,3)
  b_line <- matrix(c(mean11,mean22))
  B1 <- solve(t(hh)%*%hh)%*%(t(hh)%*%b_line) ;a1<- B1[3]
  gg1<- 0.75*a1*a1/num_l
  sigma11[2]<- sigma11[1]+gg1; sigma22[1]<- sigma22[2]+gg1
  ############  likelihood values of the initial value of calculation #########
  L0<- logL(n_samB1,d2,mi_1,mean11,sigma11,dataB1)+logL(n_samB2,d2,mi_2,mean22,sigma22,dataB2)  
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
    
    aa1<- sigma11[1]/n0_1[1]+sigma11[2]/n0_1[2]+9.0*sigma22[1]/n0_2[1]+9.0*sigma22[2]/n0_2[2]
    aa2<- sumwx_B1[1]/n0_1[1]-sumwx_B1[2]/n0_1[2]-3.0*sumwx_B2[1]/n0_2[1]+3.0*sumwx_B2[2]/n0_2[2]
    aa3<- aa2/aa1
    mean11[1]<- (sumwx_B1[1]-sigma11[1]*aa3)/n0_1[1]
    mean11[2]<- (sumwx_B1[2]+sigma11[2]*aa3)/n0_1[2]
    mean22[1]<- (sumwx_B2[1]+sigma22[1]*aa3*3)/n0_2[1]
    mean22[2]<- (sumwx_B2[2]-sigma22[2]*aa3*3)/n0_2[2]
    #########first order genetic parameter process##########
    b_line <- matrix(c(mean11,mean22))
    B11 <- solve(t(hh)%*%hh)%*%(t(hh)%*%b_line)
    a1<- B11[3]
    gg1<- 0.5*a1*a1/num_l
    sigma11[2]<- sigma11[1]+gg1; sigma22[1]<- sigma22[2]+gg1
    ###########obtain variance##########
    for(i in 1:d2) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean11[i])^2 }
    for(i in 1:d2) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean22[i])^2 }
    
    aaa0<- sigma11[1]; n_iter<- 0;aa2<-1000
    while (aa2>0.0001){
      n_iter<- n_iter+1
      aa1<- sigma11[1]/(sigma11[1]+gg1)
      sigma11[1]<- (swx_B1[1]+aa1*aa1*swx_B1[2])/(n0_1[1]+aa1*n0_1[2])
      aa2<- abs(sigma11[1]-aaa0)
      aaa0<- sigma11[1]
      if (n_iter>20) break
    }
    sigma11[2]<- sigma11[1]+gg1
    aaa0<- sigma22[2]; n_iter<- 0; aa2<- 1000
    while (aa2>0.0001){
      n_iter<- n_iter+1
      aa1<- sigma22[2]/(sigma22[2]+gg1)
      sigma22[2]<- (swx_B2[2]+aa1*aa1*swx_B2[1])/(n0_2[2]+aa1*n0_2[1])
      aa2<- abs(sigma22[2]-aaa0)
      aaa0<- sigma22[2]
      if (n_iter>20) break
    }
    sigma22[1]<- sigma22[2]+gg1
    
    ########criteria for iterations to stop#######
    L1 <- logL(n_samB1,d2,mix_pi_1,mean11,sigma11,dataB1)+logL(n_samB2,d2,mix_pi_2,mean22,sigma22,dataB2) 
    
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  
  aBCF <- L1
  AIC <- -2*aBCF + 2*5
  
  #########first order genetic parameter process##########
  b_line <- matrix(c(mean11,mean22))
  B111 <- solve(t(hh)%*%hh)%*%(t(hh)%*%b_line)
  
  #########second order genetic parameter process##########
  jj_1<- sigma1-sigma11[1]
  if(jj_1<0){jj_1<- 0}
  ll_1<- jj_1/sigma1
  
  jj_2<-sigma2-sigma22[2]
  if(jj_2<0){jj_2<- 0}
  ll_2<- jj_2/sigma2
  
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
  
  output <- data.frame("1MG-NCD",round(aBCF,4),round(AIC,4),round(t(mean11),4)," "," ",round(t(sigma11),4)," "," ",round(t(mix_pi_1),4)," "," ",
                       round(t(mean22),4)," "," ",round(t(sigma22),4)," "," ",round(t(mix_pi_2),4)," "," ",
                       round(B111[1],4),round(B111[2],4),round(B111[3],4)," "," "," ",round(jj_1,4),round(ll_1*100,4),round(jj_2,4),round(ll_2*100,4),
                       round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],
                       round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
} 

####################(B2)Model##########################################
BCFModelFun[[6]] <- function(K1,logL,df11,df21,BCFtext2){
  dataB1 <- as.matrix(as.numeric(df11[,1])); dataB2 <- as.matrix(as.numeric(df21[,1]))
  n_samB1 <- dim(dataB1)[1];n_samB2 <- dim(dataB2)[1]
  sumx1<- sum(dataB1);mean1<-mean(dataB1);sigma1<- as.numeric(var(dataB1))
  sumx2<- sum(dataB2);mean2<-mean(dataB2);sigma2<- as.numeric(var(dataB2))
  
  m_esp<-0.0001 ;num_l<- as.numeric(BCFtext2)
  #################### 2MG-AD Model########  (B2) ##############################
  d2<- 4
  mi_1<- matrix(0.25,d2,1);mi_2<- matrix(0.25,d2,1)
  sigma11<- matrix(sigma1/5,d2,1) ;sigma22<- matrix(sigma2/5,d2,1)
  a1<-sqrt(sigma1/(n_samB1-1))
  mean11<- matrix(c((mean1+2.4*a1),(mean1+0.8*a1),(mean1-0.8*a1),(mean1-2.4*a1)))
  a2<-sqrt(sigma2/(n_samB2-1))
  mean22<- matrix(c((mean2+2.4*a2),(mean2+0.8*a2),(mean2-0.8*a2),(mean2-2.4*a2)))
  
  ############  first order genetic parameter  ############
  hh<- matrix(c(1,1,1,1,0,0,0,0, 0,0,0,0,1,1,1,1, 1,1,0,0,0,0,-1,-1, 1,0,1,0,0,-1,0,-1, 0,0,0.5,0.5,0.5,0.5,0,0, 0,0.5,0,0.5,0.5,0,0.5,0),8,6)
  b_line <- matrix(c(mean11,mean22))
  B1 <- solve(t(hh)%*%hh)%*%(t(hh)%*%b_line) ;a1<- B1[3];a2<- B1[4];a3<- B1[5];a4<- B1[6] ## da; db; ha; hb
  gg2<- (0.5*a2*a2+0.25*a4*a4)/num_l; gg3<- (0.5*a1*a1+0.25*a3*a3)/num_l; gg4<- gg2+gg3;
  sigma11[2]<- sigma11[1]+gg2; sigma11[3]<- sigma11[1]+gg3; sigma11[4]<- sigma11[1]+gg4
  sigma22[1]<- sigma22[4]+gg4; sigma22[2]<- sigma22[4]+gg3; sigma22[3]<- sigma22[4]+gg2
  
  ############  likelihood values of the initial value of calculation #####################
  L0<- logL(n_samB1,d2,mi_1,mean11,sigma11,dataB1)+logL(n_samB2,d2,mi_2,mean22,sigma22,dataB2)  
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
    
    aa1<- sumwx_B1[1]/n0_1[1]-sumwx_B1[2]/n0_1[2]-sumwx_B1[3]/n0_1[3]+sumwx_B1[4]/n0_1[4]
    aa2<- sumwx_B2[1]/n0_2[1]-sumwx_B2[2]/n0_2[2]-sumwx_B2[3]/n0_2[3]+sumwx_B2[4]/n0_2[4]
    aa3<- sigma11[1]/n0_1[1]+sigma11[2]/n0_1[2]+sigma11[3]/n0_1[3]+sigma11[4]/n0_1[4]
    aa4<- sigma22[1]/n0_2[1]+sigma22[2]/n0_2[3]+sigma22[3]/n0_2[3]+sigma22[4]/n0_2[4]
    rr<- matrix(0,2,1); rr[1]<- aa1/aa3; rr[2]<- aa2/aa4 
    mean11[1]<- (sumwx_B1[1]-sigma11[1]*rr[1])/n0_1[1]
    mean11[2]<- (sumwx_B1[2]+sigma11[2]*rr[1])/n0_1[2]
    mean11[3]<- (sumwx_B1[3]+sigma11[3]*rr[1])/n0_1[3]
    mean11[4]<- (sumwx_B1[4]-sigma11[4]*rr[1])/n0_1[4]
    
    mean22[1]<- (sumwx_B2[1]-sigma22[1]*rr[2])/n0_2[1]
    mean22[2]<- (sumwx_B2[2]+sigma22[2]*rr[2])/n0_2[2]
    mean22[3]<- (sumwx_B2[3]+sigma22[3]*rr[2])/n0_2[3]
    mean22[4]<- (sumwx_B2[4]-sigma22[4]*rr[2])/n0_2[4]
    #########first order genetic parameter process##########
    b_line <- matrix(c(mean11,mean22))
    B11 <- solve(t(hh)%*%hh)%*%(t(hh)%*%b_line)
    a1<- B11[3]; a2<- B11[4];a3<- B11[5];a4<- B11[6]  ## da; db; ha; hb
    gg2<- (0.5*a2*a2+0.25*a4*a4)/num_l; gg3<- (0.5*a1*a1+0.25*a3*a3)/num_l; gg4<- gg2+gg3
    
    sigma11[2]<- sigma11[1]+gg2; sigma11[3]<- sigma11[1]+gg3; sigma11[4]<- sigma11[1]+gg4
    sigma22[1]<- sigma22[4]+gg4; sigma22[2]<- sigma22[4]+gg3; sigma22[3]<- sigma22[4]+gg2
    
    ###########obtain variance##########
    for(i in 1:d2) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean11[i])^2 }
    for(i in 1:d2) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean22[i])^2 }
    
    aaa0<- sigma11[1]; n_iter<- 0;aa5<-1000
    while (aa5>0.0001){
      n_iter<- n_iter+1
      aa2<- sigma11[1]/(sigma11[1]+gg2)
      aa3<- sigma11[1]/(sigma11[1]+gg3)
      aa4<- sigma11[1]/(sigma11[1]+gg4)
      sigma11[1]<- (swx_B1[1]+aa2*aa2*swx_B1[2]+aa3*aa3*swx_B1[3]+aa4*aa4*swx_B1[4])/(n0_1[1]+aa2*n0_1[2]+aa3*n0_1[3]+aa4*n0_1[4])
      aa5<- abs(sigma11[1]-aaa0)
      aaa0<- sigma11[1]
      if (n_iter>20) break
    }
    sigma11[2]<- sigma11[1]+gg2; sigma11[3]<- sigma11[1]+gg3; sigma11[4]<- sigma11[1]+gg4
    aaa0<- sigma22[4]; n_iter<- 0; aa5<- 1000
    while (aa2>0.0001){
      n_iter<- n_iter+1
      aa2<- sigma22[4]/(sigma22[4]+gg4)
      aa3<- sigma22[4]/(sigma22[4]+gg3)
      aa4<- sigma22[4]/(sigma22[4]+gg2)
      sigma22[4]<- (aa2*aa2*swx_B2[1]+aa3*aa3*swx_B2[2]+aa4*aa4*swx_B2[3]+swx_B2[4])/(aa2*n0_2[1]+aa3*n0_2[2]+aa4*n0_2[3]+n0_2[4])
      aa5<- abs(sigma22[4]-aaa0)
      aaa0<- sigma22[4]
      if (n_iter>20) break
    }
    sigma22[1]<- sigma22[4]+gg4; sigma22[2]<- sigma22[4]+gg3; sigma22[3]<- sigma22[4]+gg2
    
    ########criteria for iterations to stop#######
    
    L1 <- logL(n_samB1,d2,mix_pi_1,mean11,sigma11,dataB1)+logL(n_samB2,d2,mix_pi_2,mean22,sigma22,dataB2) 
    
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  
  aBCF <- L1
  AIC <- -2*aBCF + 2*8
  
  #########first order genetic parameter process##########
  b_line <- matrix(c(mean11,mean22))
  B111 <- solve(t(hh)%*%hh)%*%(t(hh)%*%b_line)
  
  #########second order genetic parameter process##########
  jj_1<- sigma1-sigma11[1]
  if(jj_1<0){jj_1<- 0}
  ll_1<- jj_1/sigma1
  
  jj_2<-sigma2-sigma22[4]
  if(jj_2<0){jj_2<- 0}
  ll_2<- jj_2/sigma2
  
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
  
  output <- data.frame("2MG-AD",round(aBCF,4),round(AIC,4),round(t(mean11),4),round(t(sigma11),4),round(t(mix_pi_1),4),
                       round(t(mean22),4),round(t(sigma22),4),round(t(mix_pi_2),4),
                       round(B111[1],4),round(B111[2],4),round(B111[3],4),round(B111[4],4),round(B111[5],4),round(B111[6],4),round(jj_1,4),round(ll_1*100,4),round(jj_2,4),round(ll_2*100,4),
                       round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],
                       round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
} 

####################(B3)Model##########################################
BCFModelFun[[7]] <- function(K1,logL,df11,df21,BCFtext2){
  dataB1 <- as.matrix(as.numeric(df11[,1])); dataB2 <- as.matrix(as.numeric(df21[,1]))
  n_samB1 <- dim(dataB1)[1];n_samB2 <- dim(dataB2)[1]
  sumx1<- sum(dataB1);mean1<-mean(dataB1);sigma1<- as.numeric(var(dataB1))
  sumx2<- sum(dataB2);mean2<-mean(dataB2);sigma2<- as.numeric(var(dataB2))
  
  m_esp<-0.0001 ;num_l<- as.numeric(BCFtext2)
  #################### 2MG-A Model########  (B3) ##############################
  d2<- 4
  mi_1<- matrix(0.25,d2,1);mi_2<- matrix(0.25,d2,1)
  sigma11<- matrix(sigma1/5,d2,1) ;sigma22<- matrix(sigma2/5,d2,1)
  a1<-sqrt(sigma1/(n_samB1-1))
  mean11<- matrix(c((mean1+2.4*a1),(mean1+1.1*a1),(mean1+0.5*a1),(mean1-2*a1)))
  a2<-sqrt(sigma2/(n_samB2-1))
  mean22<- matrix(c((mean2+2.4*a2),(mean2+1.1*a2),(mean2+0.5*a2),(mean2-2*a2)))
  
  ############  first order genetic parameter  ############
  hh<- matrix(c(1,1,1,1,0,0,0,0, 0,0,0,0,1,1,1,1, 1,1,0,0,0,0,-1,-1, 1,0,1,0,0,-1,0,-1),8,4)
  b_line <- matrix(c(mean11,mean22))
  B1 <- solve(t(hh)%*%hh)%*%(t(hh)%*%b_line) ;a1<- B1[3]; a2<- B1[4] ## da; db
  gg2<- 0.5*a2*a2/num_l; gg3<- 0.5*a1*a1/num_l; gg4<- gg2+gg3;
  sigma11[2]<- sigma11[1]+gg2; sigma11[3]<- sigma11[1]+gg3; sigma11[4]<- sigma11[1]+gg4
  sigma22[1]<- sigma22[4]+gg4; sigma22[2]<- sigma22[4]+gg3; sigma22[3]<- sigma22[4]+gg2
  
  ############ likelihood values of the initial value of calculation
  L0<- logL(n_samB1,d2,mi_1,mean11,sigma11,dataB1)+logL(n_samB2,d2,mi_2,mean22,sigma22,dataB2)  
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
    
    ############ Solution of linear equations ############ 
    aa<- matrix(0,4,4)
    aa[1,1]<- sigma11[1]/n0_1[1]+sigma11[2]/n0_1[2]+sigma11[3]/n0_1[3]+sigma11[4]/n0_1[4]
    aa[1,2]<- 0
    aa[1,3]<- sigma11[1]/n0_1[1]-sigma11[4]/n0_1[4]
    aa[1,4]<- -sigma11[2]/n0_1[2]+sigma11[3]/n0_1[3]
    aa[2,2]<- sigma22[1]/n0_2[1]+sigma22[2]/n0_2[2]+sigma22[3]/n0_2[3]+sigma22[4]/n0_2[4]
    aa[2,3]<- -sigma22[1]/n0_2[1]+sigma22[4]/n0_2[4]
    aa[2,4]<- sigma22[2]/n0_2[2]-sigma22[3]/n0_2[3]
    aa[3,3]<- sigma11[1]/n0_1[1]+sigma11[4]/n0_1[4]+sigma22[1]/n0_2[1]+sigma22[4]/n0_2[4]
    aa[3,4]<- 0
    aa[4,4]<- sigma11[2]/n0_1[2]+sigma11[3]/n0_1[3]+sigma22[2]/n0_2[2]+sigma22[3]/n0_2[3]
    for(i in 2:4){
      for(j in 1:(i-1)){
        aa[i,j]<- aa[j,i]
      }
    }
    b_line1<- matrix(0,4,1)
    b_line1[1]<- sumwx_B1[1]/n0_1[1]-sumwx_B1[2]/n0_1[2]-sumwx_B1[3]/n0_1[3]+sumwx_B1[4]/n0_1[4];                            
    b_line1[2]<- sumwx_B2[1]/n0_2[1]-sumwx_B2[2]/n0_2[2]-sumwx_B2[3]/n0_2[3]+sumwx_B2[4]/n0_2[4];
    b_line1[3]<- sumwx_B1[1]/n0_1[1]-sumwx_B1[4]/n0_1[4]-sumwx_B2[1]/n0_2[1]+sumwx_B2[4]/n0_2[4];
    b_line1[4]<- sumwx_B1[2]/n0_1[2]-sumwx_B1[3]/n0_1[3]-sumwx_B2[2]/n0_2[2]+sumwx_B2[3]/n0_2[3];
    B <- solve(aa,b_line1)
    
    mean11[1]<- (sumwx_B1[1]-sigma11[1]*(B[1]+B[3]))/n0_1[1]
    mean11[2]<- (sumwx_B1[2]+sigma11[2]*(B[1]-B[4]))/n0_1[2]
    mean11[3]<- (sumwx_B1[3]+sigma11[3]*(B[1]+B[4]))/n0_1[3]
    mean11[4]<- (sumwx_B1[4]-sigma11[4]*(B[1]-B[3]))/n0_1[4]
    
    mean22[1]<- (sumwx_B2[1]-sigma22[1]*(B[2]-B[3]))/n0_2[1]
    mean22[2]<- (sumwx_B2[2]+sigma22[2]*(B[2]+B[4]))/n0_2[2]
    mean22[3]<- (sumwx_B2[3]+sigma22[3]*(B[2]-B[4]))/n0_2[3]
    mean22[4]<- (sumwx_B2[4]-sigma22[4]*(B[2]+B[3]))/n0_2[4]
    #########first order genetic parameter process##########
    b_line <- matrix(c(mean11,mean22))
    B11 <- solve(t(hh)%*%hh)%*%(t(hh)%*%b_line)
    a1<- B11[3]; a2<- B11[4]  ## da; db
    gg2<- 0.5*a2*a2/num_l; gg3<- 0.5*a1*a1/num_l; gg4<- gg2+gg3
    
    sigma11[2]<- sigma11[1]+gg2; sigma11[3]<- sigma11[1]+gg3; sigma11[4]<- sigma11[1]+gg4
    sigma22[1]<- sigma22[4]+gg4; sigma22[2]<- sigma22[4]+gg3; sigma22[3]<- sigma22[4]+gg2
    
    ###########obtain variance##########
    for(i in 1:d2) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean11[i])^2 }
    for(i in 1:d2) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean22[i])^2 }
    
    aaa0<- sigma11[1]; n_iter<- 0; aa5<-1000
    while (aa5>0.0001){
      n_iter<- n_iter+1
      aa2<- sigma11[1]/(sigma11[1]+gg2)
      aa3<- sigma11[1]/(sigma11[1]+gg3)
      aa4<- sigma11[1]/(sigma11[1]+gg4)
      sigma11[1]<- (swx_B1[1]+aa2*aa2*swx_B1[2]+aa3*aa3*swx_B1[3]+aa4*aa4*swx_B1[4])/(n0_1[1]+aa2*n0_1[2]+aa3*n0_1[3]+aa4*n0_1[4])
      aa5<- abs(sigma11[1]-aaa0)
      aaa0<- sigma11[1]
      if (n_iter>20) break
    }
    sigma11[2]<- sigma11[1]+gg2; sigma11[3]<- sigma11[1]+gg3; sigma11[4]<- sigma11[1]+gg4
    aaa0<- sigma22[4]; n_iter<- 0; aa5<- 1000
    while (aa2>0.0001){
      n_iter<- n_iter+1
      aa2<- sigma22[4]/(sigma22[4]+gg4)
      aa3<- sigma22[4]/(sigma22[4]+gg3)
      aa4<- sigma22[4]/(sigma22[4]+gg2)
      sigma22[4]<- (aa2*aa2*swx_B2[1]+aa3*aa3*swx_B2[2]+aa4*aa4*swx_B2[3]+swx_B2[4])/(aa2*n0_2[1]+aa3*n0_2[2]+aa4*n0_2[3]+n0_2[4])
      aa5<- abs(sigma22[4]-aaa0)
      aaa0<- sigma22[4]
      if (n_iter>20) break
    }
    sigma22[1]<- sigma22[4]+gg4; sigma22[2]<- sigma22[4]+gg3; sigma22[3]<- sigma22[4]+gg2
    
    ########criteria for iterations to stop#######
    
    L1 <- logL(n_samB1,d2,mix_pi_1,mean11,sigma11,dataB1)+logL(n_samB2,d2,mix_pi_2,mean22,sigma22,dataB2) 
    
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  
  aBCF <- L1
  AIC <- -2*aBCF + 2*6
  
  #########first order genetic parameter process##########
  b_line <- matrix(c(mean11,mean22))
  B111 <- solve(t(hh)%*%hh)%*%(t(hh)%*%b_line)
  
  #########second order genetic parameter process##########
  jj_1<- sigma1-sigma11[1]
  if(jj_1<0){jj_1<- 0}
  ll_1<- jj_1/sigma1
  
  jj_2<-sigma2-sigma22[4]
  if(jj_2<0){jj_2<- 0}
  ll_2<- jj_2/sigma2
  
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
  
  output <- data.frame("2MG-A",round(aBCF,4),round(AIC,4),round(t(mean11),4),round(t(sigma11),4),round(t(mix_pi_1),4),
                       round(t(mean22),4),round(t(sigma22),4),round(t(mix_pi_2),4),
                       round(B111[1],4),round(B111[2],4),round(B111[3],4),round(B111[4],4)," "," ",round(jj_1,4),round(ll_1*100,4),round(jj_2,4),round(ll_2*100,4),
                       round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],
                       round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
} 

####################(B4)Model##########################################
BCFModelFun[[8]] <- function(K1,logL,df11,df21,BCFtext2){
  dataB1 <- as.matrix(as.numeric(df11[,1])); dataB2 <- as.matrix(as.numeric(df21[,1]))
  n_samB1 <- dim(dataB1)[1];n_samB2 <- dim(dataB2)[1]
  sumx1<- sum(dataB1);mean1<-mean(dataB1);sigma1<- as.numeric(var(dataB1))
  sumx2<- sum(dataB2);mean2<-mean(dataB2);sigma2<- as.numeric(var(dataB2))
  
  m_esp<-0.0001 ;num_l<- as.numeric(BCFtext2)
  #################### 2MG-EA Model########  (B4) ##############################
  d2<- 3
  mi_1<- as.matrix(c(0.25,0.5,0.25));mi_2<- as.matrix(c(0.25,0.5,0.25))
  sigma11<- matrix(sigma1/5,d2,1) ;sigma22<- matrix(sigma2/5,d2,1)
  a1<-sqrt(sigma1/(n_samB1-1))
  mean11<- matrix(c((mean1+2.4*a1),mean1,(mean1-2.4*a1)))
  a2<-sqrt(sigma2/(n_samB2-1))
  mean22<- matrix(c((mean2+2.4*a2),mean2,(mean2-2.4*a2)))
  
  ############  first order genetic parameter  ############
  hh<- matrix(c(1,1,1,0,0,0, 0,0,0,1,1,1, 2,1,0,0,-1,-2),6,3)
  b_line <- matrix(c(mean11,mean22))
  B1 <- solve(t(hh)%*%hh)%*%(t(hh)%*%b_line) ;a1<- B1[3]  ## d
  gg2<- 0.5*a1*a1/num_l; gg3<- a1*a1/num_l
  sigma11[2]<- sigma11[1]+gg2; sigma11[3]<- sigma11[1]+gg3
  sigma22[1]<- sigma22[3]+gg3; sigma22[2]<- sigma22[3]+gg2
  
  ############ likelihood values of the initial value of calculation
  L0<- logL(n_samB1,d2,mi_1,mean11,sigma11,dataB1)+logL(n_samB2,d2,mi_2,mean22,sigma22,dataB2)  
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
    
    ############ Solution of linear equations ############ 
    aa<- matrix(0,3,3)
    aa[1,1]<- sigma11[1]/n0_1[1]+4.0*sigma11[2]/n0_1[2]+sigma11[3]/n0_1[3]
    aa[1,2]<- 0
    aa[1,3]<- sigma11[1]/n0_1[1]-sigma11[3]/n0_1[3]
    aa[2,2]<- sigma22[1]/n0_2[1]+4.0*sigma22[2]/n0_2[2]+sigma22[3]/n0_2[3]
    aa[2,3]<- -sigma22[1]/n0_2[1]+sigma22[3]/n0_2[3]
    aa[3,3]<- sigma11[1]/n0_1[1]+sigma11[3]/n0_1[3]+sigma22[1]/n0_2[1]+sigma22[3]/n0_2[3]
    for(i in 2:3){
      for(j in 1:(i-1)){
        aa[i,j]<- aa[j,i]
      }
    }
    b_line1<- matrix(0,3,1)
    b_line1[1]<- sumwx_B1[1]/n0_1[1]-2.0*sumwx_B1[2]/n0_1[2]+sumwx_B1[3]/n0_1[3]                           
    b_line1[2]<- sumwx_B2[1]/n0_2[1]-2.0*sumwx_B2[2]/n0_2[2]+sumwx_B2[3]/n0_2[3]
    b_line1[3]<- sumwx_B1[1]/n0_1[1]-sumwx_B1[3]/n0_1[3]-sumwx_B2[1]/n0_2[1]+sumwx_B2[3]/n0_2[3];
    
    B <- solve(aa,b_line1)
    
    mean11[1]<- (sumwx_B1[1]-sigma11[1]*(B[1]+B[3]))/n0_1[1]
    mean11[2]<- (sumwx_B1[2]+sigma11[2]*2.0*B[1])/n0_1[2]
    mean11[3]<- (sumwx_B1[3]-sigma11[3]*(B[1]-B[3]))/n0_1[3]
    
    mean22[1]<- (sumwx_B2[1]-sigma22[1]*(B[2]-B[3]))/n0_2[1]
    mean22[2]<- (sumwx_B2[2]+sigma22[2]*2.0*B[2])/n0_2[2]
    mean22[3]<- (sumwx_B2[3]-sigma22[3]*(B[2]+B[3]))/n0_2[3]
    #########first order genetic parameter process##########
    b_line <- matrix(c(mean11,mean22))
    B11 <- solve(t(hh)%*%hh)%*%(t(hh)%*%b_line)
    a1<- B11[3]  ## d
    gg2<- 0.5*a1*a1/num_l; gg3<- a1*a1/num_l
    
    sigma11[2]<- sigma11[1]+gg2; sigma11[3]<- sigma11[1]+gg3; 
    sigma22[1]<- sigma22[3]+gg3; sigma22[2]<- sigma22[3]+gg2; 
    
    ###########obtain variance##########
    for(i in 1:d2) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean11[i])^2 }
    for(i in 1:d2) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean22[i])^2 }
    
    aaa0<- sigma11[1]; n_iter<- 0; aa5<-1000
    while (aa5>0.0001){
      n_iter<- n_iter+1
      aa2<- sigma11[1]/(sigma11[1]+gg2)
      aa3<- sigma11[1]/(sigma11[1]+gg3)
      sigma11[1]<- (swx_B1[1]+aa2*aa2*swx_B1[2]+aa3*aa3*swx_B1[3])/(n0_1[1]+aa2*n0_1[2]+aa3*n0_1[3])
      aa5<- abs(sigma11[1]-aaa0)
      aaa0<- sigma11[1]
      if (n_iter>20) break
    }
    sigma11[2]<- sigma11[1]+gg2; sigma11[3]<- sigma11[1]+gg3
    aaa0<- sigma22[3]; n_iter<- 0; aa5<- 1000
    while (aa2>0.0001){
      n_iter<- n_iter+1
      aa2<- sigma22[3]/(sigma22[3]+gg3)
      aa3<- sigma22[3]/(sigma22[3]+gg2)
      
      sigma22[3]<- (aa2*aa2*swx_B2[1]+aa3*aa3*swx_B2[2]+swx_B2[3])/(aa2*n0_2[1]+aa3*n0_2[2]+n0_2[3])
      aa5<- abs(sigma22[3]-aaa0)
      aaa0<- sigma22[3]
      if (n_iter>20) break
    }
    sigma22[1]<- sigma22[3]+gg3; sigma22[2]<- sigma22[3]+gg2
    ########criteria for iterations to stop#######
    
    L1 <- logL(n_samB1,d2,mix_pi_1,mean11,sigma11,dataB1)+logL(n_samB2,d2,mix_pi_2,mean22,sigma22,dataB2) 
    
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  
  aBCF <- L1
  AIC <- -2*aBCF + 2*5
  
  #########first order genetic parameter process##########
  b_line <- matrix(c(mean11,mean22))
  B111 <- solve(t(hh)%*%hh)%*%(t(hh)%*%b_line)
  
  #########second order genetic parameter process##########
  jj_1<- sigma1-sigma11[1]
  if(jj_1<0){jj_1<- 0}
  ll_1<- jj_1/sigma1
  
  jj_2<-sigma2-sigma22[3]
  if(jj_2<0){jj_2<- 0}
  ll_2<- jj_2/sigma2
  
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
  
  output <- data.frame("2MG-EA",round(aBCF,4),round(AIC,4),round(t(mean11),4)," ",round(t(sigma11),4)," ",round(t(mix_pi_1),4)," ",
                       round(t(mean22),4)," ",round(t(sigma22),4)," ",round(t(mix_pi_2),4)," ",
                       round(B111[1],4),round(B111[2],4),round(B111[3],4)," "," "," ",round(jj_1,4),round(ll_1*100,4),round(jj_2,4),round(ll_2*100,4),
                       round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],
                       round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
} 

####################(B5)Model##########################################
BCFModelFun[[9]] <- function(K1,logL,df11,df21,BCFtext2){
  dataB1 <- as.matrix(as.numeric(df11[,1])); dataB2 <- as.matrix(as.numeric(df21[,1]))
  n_samB1 <- dim(dataB1)[1];n_samB2 <- dim(dataB2)[1]
  sumx1<- sum(dataB1);mean1<-mean(dataB1);sigma1<- as.numeric(var(dataB1))
  sumx2<- sum(dataB2);mean2<-mean(dataB2);sigma2<- as.numeric(var(dataB2))
  
  m_esp<-0.0001 ;num_l<- as.numeric(BCFtext2)
  #################### 2MG-CD Model########  (B5) ##############################
  d2<- 4
  mi_1<- matrix(0.25,d2,1);mi_2<- matrix(0.25,d2,1)
  sigma11<- matrix(sigma1/5,d2,1) ;sigma22<- matrix(sigma2/5,d2,1)
  a1<-sqrt(sigma1/(n_samB1-1))
  mean11<- matrix(c((mean1+2.4*a1),(mean1+1.1*a1),(mean1+0.5*a1),(mean1-2.4*a1)))
  a2<-sqrt(sigma2/(n_samB2-1))
  mean22<- matrix(c((mean2+2.4*a2),(mean2+1.1*a2),(mean2+0.5*a2),(mean2-2.4*a2)))
  
  ############  first order genetic parameter  ############
  hh<- matrix(c(1,1,1,1,0,0,0,0, 0,0,0,0,1,1,1,1, 1,1,0.5,0.5,0.5,0.5,-1,-1, 1,0.5,1,0.5,0.5,-1,0.5,-1),8,4)
  b_line <- matrix(c(mean11,mean22))
  B1 <- solve(t(hh)%*%hh)%*%(t(hh)%*%b_line) ;a1<- B1[3]; a2<- B1[4] ## da; db
  gg2<- 0.75*a2*a2/num_l; gg3<- 0.75*a1*a1/num_l; gg4<- gg2+gg3;
  sigma11[2]<- sigma11[1]+gg2; sigma11[3]<- sigma11[1]+gg3; sigma11[4]<- sigma11[1]+gg4
  sigma22[1]<- sigma22[4]+gg4; sigma22[2]<- sigma22[4]+gg3; sigma22[3]<- sigma22[4]+gg2
  
  ############ likelihood values of the initial value of calculation
  L0<- logL(n_samB1,d2,mi_1,mean11,sigma11,dataB1)+logL(n_samB2,d2,mi_2,mean22,sigma22,dataB2)  
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
    
    ############ Solution of linear equations ############ 
    aa<- matrix(0,4,4)
    aa[1,1]<- sigma11[1]/n0_1[1]+sigma11[2]/n0_1[2]+sigma11[3]/n0_1[3]+sigma11[4]/n0_1[4]
    aa[1,2]<- 0
    aa[1,3]<- -3.0*sigma11[2]/n0_1[2]+3.0*sigma11[3]/n0_1[3]
    aa[1,4]<- 3.0*sigma11[1]/n0_1[1]-3.0*sigma11[4]/n0_1[4]
    aa[2,2]<- sigma22[1]/n0_2[1]+sigma22[2]/n0_2[2]+sigma22[3]/n0_2[3]+sigma22[4]/n0_2[4]
    aa[2,3]<- sigma22[2]/n0_2[2]-sigma22[3]/n0_2[3]
    aa[2,4]<- -sigma22[1]/n0_2[1]+sigma22[4]/n0_2[4] 
    aa[3,3]<- 9.0*sigma11[2]/n0_1[2]+9.0*sigma11[3]/n0_1[3]+sigma22[2]/n0_2[2]+sigma22[3]/n0_2[3]
    aa[3,4]<- 0
    aa[4,4]<- 9.0*sigma11[1]/n0_1[1]+9.0*sigma11[4]/n0_1[4]+sigma22[1]/n0_2[1]+sigma22[4]/n0_2[4]
    for(i in 2:4){
      for(j in 1:(i-1)){
        aa[i,j]<- aa[j,i]
      }
    }
    b_line1<- matrix(0,4,1)
    b_line1[1]<- sumwx_B1[1]/n0_1[1]-sumwx_B1[2]/n0_1[2]-sumwx_B1[3]/n0_1[3]+sumwx_B1[4]/n0_1[4];                            
    b_line1[2]<- sumwx_B2[1]/n0_2[1]-sumwx_B2[2]/n0_2[2]-sumwx_B2[3]/n0_2[3]+sumwx_B2[4]/n0_2[4];
    b_line1[3]<- 3.0*sumwx_B1[2]/n0_1[2]-3.0*sumwx_B1[3]/n0_1[3]-sumwx_B2[2]/n0_2[2]+sumwx_B2[3]/n0_2[3];
    b_line1[4]<- 3.0*sumwx_B1[1]/n0_1[1]-3.0*sumwx_B1[4]/n0_1[4]-sumwx_B2[1]/n0_2[1]+sumwx_B2[4]/n0_2[4];
    B <- solve(aa,b_line1)
    
    mean11[1]<- (sumwx_B1[1]-sigma11[1]*(B[1]+3.0*B[4]))/n0_1[1]
    mean11[2]<- (sumwx_B1[2]+sigma11[2]*(B[1]-3.0*B[3]))/n0_1[2]
    mean11[3]<- (sumwx_B1[3]+sigma11[3]*(B[1]+3.0*B[3]))/n0_1[3]
    mean11[4]<- (sumwx_B1[4]-sigma11[4]*(B[1]-3.0*B[4]))/n0_1[4]
    
    mean22[1]<- (sumwx_B2[1]-sigma22[1]*(B[2]-B[4]))/n0_2[1]
    mean22[2]<- (sumwx_B2[2]+sigma22[2]*(B[2]+B[3]))/n0_2[2]
    mean22[3]<- (sumwx_B2[3]+sigma22[3]*(B[2]-B[3]))/n0_2[3]
    mean22[4]<- (sumwx_B2[4]-sigma22[4]*(B[2]+B[4]))/n0_2[4]
    #########first order genetic parameter process##########
    b_line <- matrix(c(mean11,mean22))
    B11 <- solve(t(hh)%*%hh)%*%(t(hh)%*%b_line)
    a1<- B11[3]; a2<- B11[4]  ## da; db
    gg2<- 0.75*a2*a2/num_l; gg3<- 0.75*a1*a1/num_l; gg4<- gg2+gg3
    
    sigma11[2]<- sigma11[1]+gg2; sigma11[3]<- sigma11[1]+gg3; sigma11[4]<- sigma11[1]+gg4
    sigma22[1]<- sigma22[4]+gg4; sigma22[2]<- sigma22[4]+gg3; sigma22[3]<- sigma22[4]+gg2
    
    ###########obtain variance##########
    for(i in 1:d2) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean11[i])^2 }
    for(i in 1:d2) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean22[i])^2 }
    
    aaa0<- sigma11[1]; n_iter<- 0; aa5<-1000
    while (aa5>0.0001){
      n_iter<- n_iter+1
      aa2<- sigma11[1]/(sigma11[1]+gg2)
      aa3<- sigma11[1]/(sigma11[1]+gg3)
      aa4<- sigma11[1]/(sigma11[1]+gg4)
      sigma11[1]<- (swx_B1[1]+aa2*aa2*swx_B1[2]+aa3*aa3*swx_B1[3]+aa4*aa4*swx_B1[4])/(n0_1[1]+aa2*n0_1[2]+aa3*n0_1[3]+aa4*n0_1[4])
      aa5<- abs(sigma11[1]-aaa0)
      aaa0<- sigma11[1]
      if (n_iter>20) break
    }
    sigma11[2]<- sigma11[1]+gg2; sigma11[3]<- sigma11[1]+gg3; sigma11[4]<- sigma11[1]+gg4
    aaa0<- sigma22[4]; n_iter<- 0; aa5<- 1000
    while (aa2>0.0001){
      n_iter<- n_iter+1
      aa2<- sigma22[4]/(sigma22[4]+gg4)
      aa3<- sigma22[4]/(sigma22[4]+gg3)
      aa4<- sigma22[4]/(sigma22[4]+gg2)
      sigma22[4]<- (aa2*aa2*swx_B2[1]+aa3*aa3*swx_B2[2]+aa4*aa4*swx_B2[3]+swx_B2[4])/(aa2*n0_2[1]+aa3*n0_2[2]+aa4*n0_2[3]+n0_2[4])
      aa5<- abs(sigma22[4]-aaa0)
      aaa0<- sigma22[4]
      if (n_iter>20) break
    }
    sigma22[1]<- sigma22[4]+gg4; sigma22[2]<- sigma22[4]+gg3; sigma22[3]<- sigma22[4]+gg2
    
    ########criteria for iterations to stop#######
    
    L1 <- logL(n_samB1,d2,mix_pi_1,mean11,sigma11,dataB1)+logL(n_samB2,d2,mix_pi_2,mean22,sigma22,dataB2) 
    
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  
  aBCF <- L1
  AIC <- -2*aBCF + 2*6
  
  #########first order genetic parameter process##########
  b_line <- matrix(c(mean11,mean22))
  B111 <- solve(t(hh)%*%hh)%*%(t(hh)%*%b_line)
  
  #########second order genetic parameter process##########
  jj_1<- sigma1-sigma11[1]
  if(jj_1<0){jj_1<- 0}
  ll_1<- jj_1/sigma1
  
  jj_2<-sigma2-sigma22[4]
  if(jj_2<0){jj_2<- 0}
  ll_2<- jj_2/sigma2
  
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
  
  output <- data.frame("2MG-CD",round(aBCF,4),round(AIC,4),round(t(mean11),4),round(t(sigma11),4),round(t(mix_pi_1),4),
                       round(t(mean22),4),round(t(sigma22),4),round(t(mix_pi_2),4),
                       round(B111[1],4),round(B111[2],4),round(B111[3],4),round(B111[4],4)," "," ",round(jj_1,4),round(ll_1*100,4),round(jj_2,4),round(ll_2*100,4),
                       round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],
                       round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
} 

####################(B6)Model##########################################
BCFModelFun[[10]] <- function(K1,logL,df11,df21,BCFtext2){
  dataB1 <- as.matrix(as.numeric(df11[,1])); dataB2 <- as.matrix(as.numeric(df21[,1]))
  n_samB1 <- dim(dataB1)[1];n_samB2 <- dim(dataB2)[1]
  sumx1<- sum(dataB1);mean1<-mean(dataB1);sigma1<- as.numeric(var(dataB1))
  sumx2<- sum(dataB2);mean2<-mean(dataB2);sigma2<- as.numeric(var(dataB2))
  
  m_esp<-0.0001 ;num_l<- as.numeric(BCFtext2)
  #################### 2MG-EAD Model########  (B6) ##############################
  d2<- 3
  mi_1<- as.matrix(c(0.25,0.5,0.25));mi_2<- as.matrix(c(0.25,0.5,0.25))
  sigma11<- matrix(sigma1/5,d2,1) ;sigma22<- matrix(sigma2/5,d2,1)
  a1<-sqrt(sigma1/(n_samB1-1))
  mean11<- matrix(c((mean1+2.4*a1),mean1,(mean1-2.4*a1)))
  a2<-sqrt(sigma2/(n_samB2-1))
  mean22<- matrix(c((mean2+2.4*a2),mean2,(mean2-2.4*a2)))
  
  ############  first order genetic parameter  ############
  hh<- matrix(c(1,1,1,0,0,0, 0,0,0,1,1,1, 2,1.5,1,1,-0.5,-2),6,3)
  b_line <- matrix(c(mean11,mean22))
  B1 <- solve(t(hh)%*%hh)%*%(t(hh)%*%b_line) ;a1<- B1[3]  ## d
  gg2<- 0.75*a1*a1/num_l; gg3<- 1.5*a1*a1/num_l
  sigma11[2]<- sigma11[1]+gg2; sigma11[3]<- sigma11[1]+gg3
  sigma22[1]<- sigma22[3]+gg3; sigma22[2]<- sigma22[3]+gg2
  
  ############ likelihood values of the initial value of calculation
  L0<- logL(n_samB1,d2,mi_1,mean11,sigma11,dataB1)+logL(n_samB2,d2,mi_2,mean22,sigma22,dataB2)  
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
    
    ############ Solution of linear equations ############ 
    aa<- matrix(0,3,3)
    aa[1,1]<- sigma11[1]/n0_1[1]+4.0*sigma11[2]/n0_1[2]+sigma11[3]/n0_1[3]
    aa[1,2]<- 0
    aa[1,3]<- -12.0*sigma11[2]/n0_1[2]-6.0*sigma11[3]/n0_1[3] 
    aa[2,2]<- sigma22[1]/n0_2[1]+4.0*sigma22[2]/n0_2[2]+sigma22[3]/n0_2[3]
    aa[2,3]<- -sigma22[1]/n0_2[1]+sigma22[3]/n0_2[3] 
    aa[3,3]<- 36.0*sigma11[2]/n0_1[2]+36.0*sigma11[3]/n0_1[3]+sigma22[1]/n0_2[1]+sigma22[3]/n0_2[3]
    for(i in 2:3){
      for(j in 1:(i-1)){
        aa[i,j]<- aa[j,i]
      }
    }
    b_line1<- matrix(0,3,1)
    b_line1[1]<- sumwx_B1[1]/n0_1[1]-2.0*sumwx_B1[2]/n0_1[2]+sumwx_B1[3]/n0_1[3]                           
    b_line1[2]<- sumwx_B2[1]/n0_2[1]-2.0*sumwx_B2[2]/n0_2[2]+sumwx_B2[3]/n0_2[3]
    b_line1[3]<- 6.0*sumwx_B1[2]/n0_1[2]-6.0*sumwx_B1[3]/n0_1[3]-sumwx_B2[1]/n0_2[1]+sumwx_B2[3]/n0_2[3];
    
    B <- solve(aa,b_line1)
    
    mean11[1]<- (sumwx_B1[1]-sigma11[1]*B[1])/n0_1[1]
    mean11[2]<- (sumwx_B1[2]+sigma11[2]*(B[1]*2.0-6.0*B[3]))/n0_1[2]
    mean11[3]<- (sumwx_B1[3]-sigma11[3]*(B[1]-6.0*B[3]))/n0_1[3]
    
    mean22[1]<- (sumwx_B2[1]-sigma22[1]*(B[2]-B[3]))/n0_2[1]
    mean22[2]<- (sumwx_B2[2]+sigma22[2]*2.0*B[2])/n0_2[2]
    mean22[3]<- (sumwx_B2[3]-sigma22[3]*(B[2]+B[3]))/n0_2[3]
    #########first order genetic parameter process##########
    b_line <- matrix(c(mean11,mean22))
    B11 <- solve(t(hh)%*%hh)%*%(t(hh)%*%b_line)
    a1<- B11[3]  ## d
    gg2<- 0.75*a1*a1/num_l; gg3<- 1.5*a1*a1/num_l
    
    sigma11[2]<- sigma11[1]+gg2; sigma11[3]<- sigma11[1]+gg3; 
    sigma22[1]<- sigma22[3]+gg3; sigma22[2]<- sigma22[3]+gg2; 
    
    ###########obtain variance##########
    for(i in 1:d2) {  swx_B1[i] <- WW_B1[i,]%*%(dataB1-mean11[i])^2 }
    for(i in 1:d2) {  swx_B2[i] <- WW_B2[i,]%*%(dataB2-mean22[i])^2 }
    
    aaa0<- sigma11[1]; n_iter<- 0; aa5<-1000
    while (aa5>0.0001){
      n_iter<- n_iter+1
      aa2<- sigma11[1]/(sigma11[1]+gg2)
      aa3<- sigma11[1]/(sigma11[1]+gg3)
      sigma11[1]<- (swx_B1[1]+aa2*aa2*swx_B1[2]+aa3*aa3*swx_B1[3])/(n0_1[1]+aa2*n0_1[2]+aa3*n0_1[3])
      aa5<- abs(sigma11[1]-aaa0)
      aaa0<- sigma11[1]
      if (n_iter>20) break
    }
    sigma11[2]<- sigma11[1]+gg2; sigma11[3]<- sigma11[1]+gg3
    aaa0<- sigma22[3]; n_iter<- 0; aa5<- 1000
    while (aa2>0.0001){
      n_iter<- n_iter+1
      aa2<- sigma22[3]/(sigma22[3]+gg3)
      aa3<- sigma22[3]/(sigma22[3]+gg2)
      
      sigma22[3]<- (aa2*aa2*swx_B2[1]+aa3*aa3*swx_B2[2]+swx_B2[3])/(aa2*n0_2[1]+aa3*n0_2[2]+n0_2[3])
      aa5<- abs(sigma22[3]-aaa0)
      aaa0<- sigma22[3]
      if (n_iter>20) break
    }
    sigma22[1]<- sigma22[3]+gg3; sigma22[2]<- sigma22[3]+gg2
    ########criteria for iterations to stop#######
    
    L1 <- logL(n_samB1,d2,mix_pi_1,mean11,sigma11,dataB1)+logL(n_samB2,d2,mix_pi_2,mean22,sigma22,dataB2) 
    
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  
  aBCF <- L1
  AIC <- -2*aBCF + 2*5
  
  #########first order genetic parameter process##########
  b_line <- matrix(c(mean11,mean22))
  B111 <- solve(t(hh)%*%hh)%*%(t(hh)%*%b_line)
  
  #########second order genetic parameter process##########
  jj_1<- sigma1-sigma11[1]
  if(jj_1<0){jj_1<- 0}
  ll_1<- jj_1/sigma1
  
  jj_2<- sigma2-sigma22[3]
  if(jj_2<0){jj_2<- 0}
  ll_2<- jj_2/sigma2
  
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
  
  output <- data.frame("2MG-EAD",round(aBCF,4),round(AIC,4),round(t(mean11),4)," ",round(t(sigma11),4)," ",round(t(mix_pi_1),4)," ",
                       round(t(mean22),4)," ",round(t(sigma22),4)," ",round(t(mix_pi_2),4)," ",
                       round(B111[1],4),round(B111[2],4),round(B111[3],4)," "," "," ",round(jj_1,4),round(ll_1*100,4),round(jj_2,4),round(ll_2*100,4),
                       round(u_B1[1],4),tt_B1[1],round(u_B1[2],4),tt_B1[2],round(u_B1[3],4),tt_B1[3],round(WW2_B1,4),tt_B1[4],round(D_B1,4),tt_B1[5],
                       round(u_B2[1],4),tt_B2[1],round(u_B2[2],4),tt_B2[2],round(u_B2[3],4),tt_B2[3],round(WW2_B2,4),tt_B2[4],round(D_B2,4),tt_B2[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi_1,mi_2)
  return(OUTPUT)
} 



K1BCF <- function(x){
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

logLBCF <- function(nm,nng,mi,mn,s,d1) { sum2 <- sum(log(dmixnorm(d1,mn,sqrt(s),mi)));return (sum2) } 



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
    BCFModelFun[[i]](K1BCF,logLBCF,df11,df21,BCFtext2)[[1]]
  }
  stopCluster(cl)
  mi_1<-NULL;mi_2<-NULL
}else{
   
  allresultq=switch(model,"0MG"=BCFModelFun[[1]](K1BCF,logLBCF,df11,df21,BCFtext2),"1MG-AD" = BCFModelFun[[2]](K1BCF,logLBCF,df11,df21,BCFtext2),"1MG-A"=BCFModelFun[[3]](K1BCF,logLBCF,df11,df21,BCFtext2),"1MG-EAD"=BCFModelFun[[4]](K1BCF,logLBCF,df11,df21,BCFtext2),"1MG-NCD"=BCFModelFun[[5]](K1BCF,logLBCF,df11,df21,BCFtext2),
               "2MG-AD"=BCFModelFun[[6]](K1BCF,logLBCF,df11,df21,BCFtext2),"2MG-A"=BCFModelFun[[7]](K1BCF,logLBCF,df11,df21,BCFtext2),"2MG-EA"=BCFModelFun[[8]](K1BCF,logLBCF,df11,df21,BCFtext2),"2MG-CD"=BCFModelFun[[9]](K1BCF,logLBCF,df11,df21,BCFtext2),"2MG-EAD"=BCFModelFun[[10]](K1BCF,logLBCF,df11,df21,BCFtext2))
  
  allresult<-allresultq[[1]]
  if(model!="0MG"){
    mi_1<-allresultq[[2]];mi_2<-allresultq[[3]]  
  }else{
    mi_1<-NULL;mi_2<-NULL
  }
}
  
colnames(allresult) <- BCFcolname

out<-list(allresult,mi_1,mi_2)
return(out)
}
















