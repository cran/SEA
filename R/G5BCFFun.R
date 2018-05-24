G5BCFFun<-function(df,model,G5BCFtext2){

data<-sapply(df,as.character)

dP1<-data[-1,which(data[1,]=="P1")];P1<-as.numeric(dP1[which(is.na(as.numeric(dP1))==FALSE)]);df11<-as.data.frame(P1)
dF1<-data[-1,which(data[1,]=="F1")];F1<-as.numeric(dF1[which(is.na(as.numeric(dF1))==FALSE)]);df21<-as.data.frame(F1)
dP2<-data[-1,which(data[1,]=="P2")];P2<-as.numeric(dP2[which(is.na(as.numeric(dP2))==FALSE)]);df31<-as.data.frame(P2)
dB12<-data[-1,which(data[1,]=="B12")];B12<-as.numeric(dB12[which(is.na(as.numeric(dB12))==FALSE)]);df41<-as.data.frame(B12)
dB22<-data[-1,which(data[1,]=="B22")];B22<-as.numeric(dB22[which(is.na(as.numeric(dB22))==FALSE)]);df51<-as.data.frame(B22)

G5BCFcolname<-c("Model","Log_Max_likelihood_Value","AIC","mean[P1]","mean[F1]","mean[P2]","mean(B1:2)[1]","mean(B1:2)[2]","mean(B1:2)[3]","mean(B1:2)[4]",
                "var(B1:2)[1]","var(B1:2)[2]","var(B1:2)[3]","var(B1:2)[4]","Proportion(B1:2)[1]","Proportion(B1:2)[2]","Proportion(B1:2)[3]","Proportion(B1:2)[4]",
                "mean(B2:2)[1]","mean(B2:2)[2]","mean(B2:2)[3]","mean(B2:2)[4]","var(B2:2)[1]","var(B2:2)[2]","var(B2:2)[3]","var(B1:2)[4]",
                "Proportion(B2:2)[1]","Proportion(B2:2)[2]","Proportion(B2:2)[3]","Proportion(B2:2)[4]","Var(Residual)",
                "m","da(d)","db","ha(h)","hb","[d]","[h]","Major-Gene Var(B1:2)","Heritability(Major-Gene(B1:2))(%)","Poly-Gene Var(B1:2)","Heritability(Poly-Gene(B1:2))(%)",
                "Major-Gene Var(B2:2)","Heritability(Major-Gene(B2:2))(%)","Poly-Gene Var(B2:2)","Heritability(Poly-Gene(B2:2))(%)",
                "U1 square(P1)","P(U1 square(P1))","U2 square(P1)","P(U2 square(P1))","U3 square(P1)","P(U3 square(P1))","nW square(P1)","P(nW square(P1))","Dn(P1)","P(Dn(P1))",
                "U1 square(F1)","P(U1 square(F1))","U2 square(F1)","P(U2 square(F1))","U3 square(F1)","P(U3 square(F1))","nW square(F1)","P(nW square(F1))","Dn(F1)","P(Dn(F1))",
                "U1 square(P2)","P(U1 square(P2))","U2 square(P2)","P(U2 square(P2))","U3 square(P2)","P(U3 square(P2))","nW square(P2)","P(nW square(P2))","Dn(P2)","P(Dn(P2))",
                "U1 square(B1:2)","P(U1 square(B1:2))","U2 square(B1:2)","P(U2 square(B1:2))","U3 square(B1:2)","P(U3 square(B1:2))","nW square(B1:2)","P(nW square(B1:2))","Dn(B1:2)","P(Dn(B1:2))",
                "U1 square(B2:2)","P(U1 square(B2:2))","U2 square(B2:2)","P(U2 square(B2:2))","U3 square(B2:2)","P(U3 square(B2:2))","nW square(B2:2)","P(nW square(B2:2))","Dn(B2:2)","P(Dn(B2:2))")

G5BCFModelFun<-list(NA)
###################define each model function############################
################A-1(1MG-AD)##############################################
G5BCFModelFun[[1]] <- function(K1,logL,df11,df21,df31,df41,df51,G5BCFtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]))
  nn<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataB1)[1],dim(dataB2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2)))
  ss<-s[c(1:3)]-sumx[c(1:3)]^2/nn[c(1:3)]
  m<-matrix(0,13,1)
  m[11]<-mean(dataP1);m[12]<-mean(dataF1);m[13]<-mean(dataP2);man0<-mean(dataB1);man1<-mean(dataB2)
  sigma0<-sum(ss)/(nn[1]+nn[2]+nn[3]-3)
  sigma1<-var(dataB1);sigma2<-var(dataB2)
  d1<-1;d2<-2;d3<-3;d4<-4
  m_esp <- 0.0001;num_l <- as.numeric(G5BCFtext2)
  ###############procedure start###########################
  mi<-matrix(0,6,1);mi[c(1,2,5,6)]<-0.5
  sigma<-matrix(0,11,1)
  sigma[1]<-sigma1/3;sigma[2]<-sigma[1];sigma[5]<-sigma2/3;sigma[6]<-sigma[5];sigma[11]<-sigma0
  a1<-sqrt(sigma1/(nn[4]-1))
  m[1]<-man0+a1; m[2]<-man0-a1
  a1<-sqrt(sigma2/(nn[5]-1))
  m[5]<-man1+a1;m[6]<-man1-a1
  mm1<-m[c(1,2)];mmi1<-mi[c(1,2)];ssigma1<-sigma[c(1,2)]
  mm2<-m[c(5,6)];mmi2<-mi[c(5,6)];ssigma2<-sigma[c(5,6)]
  L0<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
    sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)))+
    sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,2,nn[4]); swx1 <- matrix(0,2,1)
  W2 <- matrix(0,2,nn[5]); swx2 <- matrix(0,2,1)
  sumwx<-matrix(0,6,1);mix_pi<-matrix(0,6,1)
  n0<-matrix(0,6,1);s0<-matrix(0,4,1);swx<-matrix(0,6,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:2) { W1[i,] <- mi[i]*dnorm(dataB1,m[i],sqrt(sigma[i]))/dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)}
    mix_pi[c(1,2)] <- as.matrix(rowSums(W1)/nn[4])
    sumwx[c(1,2)] <- W1%*%dataB1
    for(i in 1:2) { W2[i,] <- mmi2[i]*dnorm(dataB2,mm2[i],sqrt(ssigma2[i]))/dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)}
    mix_pi[c(5,6)] <- as.matrix(rowSums(W2)/nn[5])
    sumwx[c(5,6)] <- W2%*%dataB2
    n0[c(1,2)]<-nn[4]*mix_pi[c(1,2)];n0[c(5,6)]<-nn[5]*mix_pi[c(5,6)]
    n0[c(1,2,5,6)][abs(n0[c(1,2,5,6)])<0.000001]<-0.000001
    s0[1]<-sumx[1]+sumwx[1];s0[2]<-nn[1]+n0[1];s0[3]<-sumx[3]+sumwx[6];s0[4]<-nn[3]+n0[6]
    aa3<-(s0[1]/s0[2]+2*sumx[2]/nn[2]+s0[3]/s0[4]-4*sumwx[2]/n0[2])/(sigma[11]/s0[2]+4*sigma[11]/nn[2]+sigma[11]/s0[4]+16*sigma[2]/n0[2])
    m[11]<-(s0[1]-sigma[11]*aa3)/s0[2]                # mean1.
    m[12]<-(sumx[2]-sigma[11]*aa3*2)/nn[2]            # mean2.
    m[13]<-(s0[3]-sigma[11]*aa3)/s0[4]                # mean3.
    m[2]<-(sumwx[2]+sigma[2]*aa3*4)/n0[2]
    m[1]<-m[11];m[5]<-m[2];m[6]<-m[13]
    ########first order genetic parameters###############   
    hh1<-matrix(c(1,1,1,1,1,0,-1,0,0,1,0,0.5),4,3)
    mm<-m[c(11,12,13,2)]
    B<-solve(crossprod(hh1,hh1))%*%crossprod(hh1,mm)
    a1<-B[2];a2<-B[3]
    a1<-(0.5*a1^2+0.25*a2^2)/num_l
    ########obtain variance#######################
    ss1<-sum((dataP1-m[11])^2);ss3<-sum((dataP2-m[13])^2);ss2<-sum((dataF1-m[12])^2)
    mm2<-m[c(5,6)]
    for(i in 1:d2) {swx1[i] <- W1[i,]%*%(dataB1-m[i])^2 } ;for(i in 1:d2) {swx2[i] <- W2[i,]%*%(dataB2-mm2[i])^2 }  
    swx[c(1,2)]<-swx1;swx[c(5,6)]<-swx2
    s0[1]<-ss1+ss2+ss3;s0[2]<-nn[1]+nn[2]+nn[3]
    aaa0<-sigma[11];n_iter<-0;aa2<-1000
    while (aa2>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma[11]/(sigma[11]+a1)
      sigma[11]<-(s0[1]+(swx[1]+swx[6]+aa1^2*(swx[2]+swx[5])))/(s0[2]+n0[1]+n0[6]+aa1*(n0[2]+n0[5]))
      aa2<-abs(sigma[11]-aaa0)
      aaa0<-sigma[11]
      if (n_iter>20) break
    }
    sigma[1]<-sigma[6]<-sigma[11];sigma[2]<-sigma[1]+a1;sigma[5]<-sigma[2]
    mm1<-m[c(1,2)];ssigma1<-sigma[c(1,2)];mm2<-m[c(5,6)];ssigma2<-sigma[c(5,6)]
    mix_pi1<-mix_pi[c(1,2)];mix_pi2<-mix_pi[c(5,6)]
    ########criteria for iterations to stop########
    L1<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
      sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mix_pi2)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
    if(iteration>300)break
  }
  abc<-L0
  AIC<--2*abc+2*8
  ########second order genetic parameters###############   
  B1jj<-sigma1-sigma[11]
  if(B1jj<0){B1jj<-0}
  B1ll<-B1jj/sigma1
  B2jj<-sigma2-sigma[11]
  if(B2jj<0) {B2jj<-0}
  B2ll<-B2jj/sigma2   
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*nn[1])
  P1bmw <- matrix(0,nn[1],1)
  P1gg <- (dataP1 - m[11])/sqrt(as.vector(sigma[11]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn1 < nn[1]){P1bmw <- P1bmw+runif(nn[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:nn[1])) - 0.5)/nn[1])^2)
  P1u<- as.matrix(c(12*nn[1]*((P1dd[1]/nn[1]-0.5)^2),((45*nn[1])/4)*((P1dd[2]/nn[1]-1/3)^2),180*nn[1]*((P1dd[3]/nn[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,nn[1]))))
  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*nn[2])
  F1bmw <- matrix(0,nn[2],1)
  F1gg <- (dataF1 - m[12])/sqrt(as.vector(sigma[11]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn1 < nn[2]){F1bmw <- F1bmw+runif(nn[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:nn[2])) - 0.5)/nn[2])^2)
  F1u<- as.matrix(c(12*nn[2]*((F1dd[1]/nn[2]-0.5)^2),((45*nn[2])/4)*((F1dd[2]/nn[2]-1/3)^2),180*nn[2]*((F1dd[3]/nn[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,nn[2]))))
  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)
  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*nn[3])
  P2bmw <- matrix(0,nn[3],1)
  P2gg <- (dataP2 - m[13])/sqrt(as.vector(sigma[11]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn1 < nn[3]){P2bmw <- P2bmw+runif(nn[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:nn[3])) - 0.5)/nn[3])^2)
  P2u<- as.matrix(c(12*nn[3]*((P2dd[1]/nn[3]-0.5)^2),((45*nn[3])/4)*((P2dd[2]/nn[3]-1/3)^2),180*nn[3]*((P2dd[3]/nn[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,nn[3]))))
  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1); 
  B1w1<-1/(12*nn[4])
  B1bmw <- matrix(0,nn[4],1); B1bmwsl <- matrix(0,nn[4],d2)
  for(i in 1:d2){
    B1gg <- (dataB1 - m[i])/sqrt(sigma[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B1P2)))[1]
  if(nn1 < nn[4]){B1P2 <- B1P2+runif(nn[4])/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*nn[4]) + sum((B1P2 - (as.matrix(c(1:nn[4])) - 0.5)/nn[4])^2)
  B1u <- as.matrix(c(12*nn[4]*((B1dd[1]/nn[4]-0.5)^2),((45*nn[4])/4)*((B1dd[2]/nn[4]-1/3)^2),180*nn[4]*((B1dd[3]/nn[4]-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),(1-pkolm(B1D,nn[4]))))
  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2); 
  B2w1<-1/(12*nn[5])
  B2bmw <- matrix(0,nn[5],1); B2bmwsl <- matrix(0,nn[5],d2)
  for(i in 1:d2){
    B2gg <- (dataB2 - mm2[i])/sqrt(ssigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B2P2)))[1]
  if(nn1 < nn[5]){B2P2 <- B2P2+runif(nn[5])/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*nn[5]) + sum((B2P2 - (as.matrix(c(1:nn[5])) - 0.5)/nn[5])^2)
  B2u <- as.matrix(c(12*nn[5]*((B2dd[1]/nn[5]-0.5)^2),((45*nn[5])/4)*((B2dd[2]/nn[5]-1/3)^2),180*nn[5]*((B2dd[3]/nn[5]-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),(1-pkolm(B2D,nn[5]))))
  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("1MG-AD",round(abc,4),round(AIC,4),round(m[11],4),round(m[12],4),round(m[13],4), round(m[1],4),round(m[2],4)," "," ",round(sigma[1],4),round(sigma[2],4),
                       " "," ",round(mix_pi[1],4),round(mix_pi[2],4)," "," ",round(m[5],4),round(m[6],4)," "," ",round(sigma[5],4),round(sigma[6],4)," "," ",        
                       round(mix_pi[5],4),round(mix_pi[6],4)," "," ",round(sigma[11],4),round(B[1],4),round(B[2],4)," ",round(B[3],4)," "," "," ",round(B1jj,4),round(B1ll*100,4)," "," ",round(B2jj,4),round(B2ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,as.matrix(mmi1),as.matrix(mmi2))
  return(OUTPUT)
}
######################1MG-A(A-2)##################################
G5BCFModelFun[[2]] <- function(K1,logL,df11,df21,df31,df41,df51,G5BCFtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]))
  nn<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataB1)[1],dim(dataB2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2)))
  ss<-s[c(1:3)]-sumx[c(1:3)]^2/nn[c(1:3)]
  m<-matrix(0,13,1);m[11]<-mean(dataP1);m[12]<-mean(dataF1);m[13]<-mean(dataP2);man0<-mean(dataB1);man1<-mean(dataB2)
  sigma0<-sum(ss)/(nn[1]+nn[2]+nn[3]-3)
  sigma1<-var(dataB1);sigma2<-var(dataB2)
  d1<-1;d2<-2;d3<-3;d4<-4
  m_esp <- 0.0001;num_l <- as.numeric(G5BCFtext2)
  ###############procedure start###########################  
  mi<-matrix(0,6,1);mi[c(1,2,5,6)]<-0.5
  sigma<-matrix(0,11,1)
  sigma[1]<-sigma1/3;sigma[2]<-sigma[1];sigma[5]<-sigma2/3;sigma[6]<-sigma[5];sigma[11]<-sigma0
  a1<-sqrt(sigma1/(nn[4]-1))
  m[1]<-man0+a1;m[2]<-man0-a1
  a1<-sqrt(sigma2/(nn[5]-1))
  m[5]<-man1+a1;m[6]<-man1-a1
  mm1<-m[c(1,2)];mmi1<-mi[c(1,2)];ssigma1<-sigma[c(1,2)]
  mm2<-m[c(5,6)];mmi2<-mi[c(5,6)];ssigma2<-sigma[c(5,6)]
  L0<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
    sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)))+
    sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,2,nn[4]); swx1 <- matrix(0,2,1)
  W2 <- matrix(0,2,nn[5]); swx2 <- matrix(0,2,1)
  sumwx<-matrix(0,6,1);mix_pi<-matrix(0,6,1)
  n0<-matrix(0,6,1);s0<-matrix(0,6,1);swx<-matrix(0,6,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:2) { W1[i,] <- mi[i]*dnorm(dataB1,m[i],sqrt(sigma[i]))/dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)}
    mix_pi[c(1,2)] <- as.matrix(rowSums(W1)/nn[4])
    sumwx[c(1,2)] <- W1%*%dataB1
    for(i in 1:2) { W2[i,] <- mmi2[i]*dnorm(dataB2,mm2[i],sqrt(ssigma2[i]))/dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)}
    mix_pi[c(5,6)] <- as.matrix(rowSums(W2)/nn[5])
    sumwx[c(5,6)] <- W2%*%dataB2
    n0[c(1,2)]<-nn[4]*mix_pi[c(1,2)];n0[c(5,6)]<-nn[5]*mix_pi[c(5,6)]
    n0[c(1,2,5,6)][abs(n0[c(1,2,5,6)])<0.000001]<-0.000001
    s0[1]<-sumx[1]+sumwx[1];s0[2]<-nn[1]+n0[1];s0[3]<-sumx[2]*sigma[2]+(sigma[11]/sigma[2])*(sumwx[2]+sumwx[5])
    s0[4]<-nn[2]*sigma[2]+(sigma[11]/sigma[2])*(n0[2]+n0[5]);s0[5]<-sumx[3]+sumwx[6];s0[6]<-nn[3]+n0[6]
    aa3<-(s0[1]/s0[2]-2*s0[3]/s0[4]+s0[5]/s0[6])/(sigma[11]/s0[2]+4*sigma[11]*sigma[2]/s0[4]+sigma[11]/s0[6])
    m[11]<-(s0[1]-sigma[11]*aa3)/s0[2]                # mean1.
    m[12]<-(s0[3]+sigma[11]*sigma[2]*aa3*2)/s0[4]   # mean2.
    m[13]<-(s0[5]-sigma[11]*aa3)/s0[6]                # mean3.
    m[1]<-m[11];m[2]<-m[12];m[5]<-m[12];m[6]<-m[13]
    ############first order genetic parameters#######
    hh<-matrix(c(1,1,1,1,0,-1),3,2)
    mm<-m[c(11,12,13)]
    B<-solve(crossprod(hh,hh))%*%crossprod(hh,mm)
    a1<-B[2];a1<-0.5*a1^2/num_l
    ########obtain variance#######################
    ss1<-sum((dataP1-m[11])^2);ss3<-sum((dataP2-m[13])^2);ss2<-sum((dataF1-m[12])^2)
    mm2<-m[c(5,6)]
    for(i in 1:d2) {swx1[i] <- W1[i,]%*%(dataB1-m[i])^2 } ;for(i in 1:d2) {swx2[i] <- W2[i,]%*%(dataB2-mm2[i])^2 }  
    swx[c(1,2)]<-swx1;swx[c(5,6)]<-swx2
    s0[1]<-ss1+ss2+ss3;s0[2]<-nn[1]+nn[2]+nn[3]
    aaa0<-sigma[11];n_iter<-0;aa2<-1000
    while (aa2>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma[11]/(sigma[11]+a1)
      sigma[11]<-(s0[1]+(swx[1]+swx[6]+aa1^2*(swx[2]+swx[5])))/(s0[2]+n0[1]+n0[6]+aa1*(n0[2]+n0[5]))
      aa2<-abs(sigma[11]-aaa0)
      aaa0<-sigma[11]
      if (n_iter>20) break
    }
    sigma[1]<-sigma[6]<-sigma[11];sigma[2]<-sigma[1]+a1;sigma[5]<-sigma[2]
    mm1<-m[c(1,2)];ssigma1<-sigma[c(1,2)];mm2<-m[c(5,6)];ssigma2<-sigma[c(5,6)]
    mix_pi1<-mix_pi[c(1,2)];mix_pi2<-mix_pi[c(5,6)]
    ########criteria for iterations to stop########
    L1<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
      sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mix_pi2)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
    if(iteration>300)break
  }
  abc<-L0
  AIC<--2*abc+2*7
  ########second order genetic parameters###############   
  B1jj<-sigma1-sigma[11]
  if(B1jj<0){B1jj<-0}
  B1ll<-B1jj/sigma1
  B2jj<-sigma2-sigma[11]
  if(B2jj<0) {B2jj<-0}
  B2ll<-B2jj/sigma2   
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*nn[1])
  P1bmw <- matrix(0,nn[1],1)
  P1gg <- (dataP1 - m[11])/sqrt(as.vector(sigma[11]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn1 < nn[1]){P1bmw <- P1bmw+runif(nn[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:nn[1])) - 0.5)/nn[1])^2)
  P1u<- as.matrix(c(12*nn[1]*((P1dd[1]/nn[1]-0.5)^2),((45*nn[1])/4)*((P1dd[2]/nn[1]-1/3)^2),180*nn[1]*((P1dd[3]/nn[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,nn[1]))))
  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*nn[2])
  F1bmw <- matrix(0,nn[2],1)
  F1gg <- (dataF1 - m[12])/sqrt(as.vector(sigma[11]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn1 < nn[2]){F1bmw <- F1bmw+runif(nn[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:nn[2])) - 0.5)/nn[2])^2)
  F1u<- as.matrix(c(12*nn[2]*((F1dd[1]/nn[2]-0.5)^2),((45*nn[2])/4)*((F1dd[2]/nn[2]-1/3)^2),180*nn[2]*((F1dd[3]/nn[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,nn[2]))))
  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)
  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*nn[3])
  P2bmw <- matrix(0,nn[3],1)
  P2gg <- (dataP2 - m[13])/sqrt(as.vector(sigma[11]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn1 < nn[3]){P2bmw <- P2bmw+runif(nn[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:nn[3])) - 0.5)/nn[3])^2)
  P2u<- as.matrix(c(12*nn[3]*((P2dd[1]/nn[3]-0.5)^2),((45*nn[3])/4)*((P2dd[2]/nn[3]-1/3)^2),180*nn[3]*((P2dd[3]/nn[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,nn[3]))))
  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1); 
  B1w1<-1/(12*nn[4])
  B1bmw <- matrix(0,nn[4],1); B1bmwsl <- matrix(0,nn[4],d2)
  for(i in 1:d2){
    B1gg <- (dataB1 - m[i])/sqrt(sigma[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B1P2)))[1]
  if(nn1 < nn[4]){B1P2 <- B1P2+runif(nn[4])/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*nn[4]) + sum((B1P2 - (as.matrix(c(1:nn[4])) - 0.5)/nn[4])^2)
  B1u <- as.matrix(c(12*nn[4]*((B1dd[1]/nn[4]-0.5)^2),((45*nn[4])/4)*((B1dd[2]/nn[4]-1/3)^2),180*nn[4]*((B1dd[3]/nn[4]-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),(1-pkolm(B1D,nn[4]))))
  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2); 
  B2w1<-1/(12*nn[5])
  B2bmw <- matrix(0,nn[5],1); B2bmwsl <- matrix(0,nn[5],d2)
  for(i in 1:d2){
    B2gg <- (dataB2 - mm2[i])/sqrt(ssigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B2P2)))[1]
  if(nn1 < nn[5]){B2P2 <- B2P2+runif(nn[5])/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*nn[5]) + sum((B2P2 - (as.matrix(c(1:nn[5])) - 0.5)/nn[5])^2)
  B2u <- as.matrix(c(12*nn[5]*((B2dd[1]/nn[5]-0.5)^2),((45*nn[5])/4)*((B2dd[2]/nn[5]-1/3)^2),180*nn[5]*((B2dd[3]/nn[5]-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),(1-pkolm(B2D,nn[5]))))
  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("1MG-A",round(abc,4),round(AIC,4),round(m[11],4),round(m[12],4),round(m[13],4), round(m[1],4),round(m[2],4)," "," ",round(sigma[1],4),round(sigma[2],4),
                       " "," ",round(mix_pi[1],4),round(mix_pi[2],4)," "," ",round(m[5],4),round(m[6],4)," "," ",round(sigma[5],4),round(sigma[6],4)," "," ",        
                       round(mix_pi[5],4),round(mix_pi[6],4)," "," ",round(sigma[11],4),round(B[1],4),round(B[2],4)," "," "," "," "," ",round(B1jj,4),round(B1ll*100,4)," "," ",round(B2jj,4),round(B2ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,as.matrix(mmi1),as.matrix(mmi2))
  return(OUTPUT)
}

###################1MG-EAD(A-3)#########################
G5BCFModelFun[[3]] <- function(K1,logL,df11,df21,df31,df41,df51,G5BCFtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]))
  nn<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataB1)[1],dim(dataB2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2)))
  ss<-s[c(1:3)]-sumx[c(1:3)]^2/nn[c(1:3)]
  m<-matrix(0,13,1);m[11]<-mean(dataP1);m[12]<-mean(dataF1);m[13]<-mean(dataP2);man0<-mean(dataB1);man1<-mean(dataB2)
  sigma0<-sum(ss)/(nn[1]+nn[2]+nn[3]-3)
  sigma1<-var(dataB1);sigma2<-var(dataB2)
  d1<-1;d2<-2;d3<-3;d4<-4
  m_esp <- 0.0001;num_l <- as.numeric(G5BCFtext2)
  ###############procedure start###########################  
  mi<-matrix(0,6,1);mi[c(1,2,5,6)]<-0.5
  sigma<-matrix(0,11,1);sigma[1]<-sigma1/3;sigma[2]<-sigma[1];sigma[5]<-sigma2/3;sigma[6]<-sigma[5];sigma[11]<-sigma0
  a1<-sqrt(sigma1/(nn[4]-1))
  m[1]<-man0+a1;m[2]<-man0-a1
  a1<-sqrt(sigma2/(nn[5]-1))
  m[5]<-man1+a1;m[6]<-man1-a1
  mm1<-m[c(1,2)];mmi1<-mi[c(1,2)];ssigma1<-sigma[c(1,2)]
  mm2<-m[c(5,6)];mmi2<-mi[c(5,6)];ssigma2<-sigma[c(5,6)]
  L0<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
    sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)))+
    sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,2,nn[4]); swx1 <- matrix(0,2,1)
  W2 <- matrix(0,2,nn[5]); swx2 <- matrix(0,2,1)
  sumwx<-matrix(0,6,1);mix_pi<-matrix(0,6,1)
  n0<-matrix(0,6,1);s0<-matrix(0,6,1);swx<-matrix(0,6,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:2) { W1[i,] <- mi[i]*dnorm(dataB1,m[i],sqrt(sigma[i]))/dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)}
    mix_pi[c(1,2)] <- as.matrix(rowSums(W1)/nn[4])
    sumwx[c(1,2)] <- W1%*%dataB1
    for(i in 1:2) { W2[i,] <- mmi2[i]*dnorm(dataB2,mm2[i],sqrt(ssigma2[i]))/dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)}
    mix_pi[c(5,6)] <- as.matrix(rowSums(W2)/nn[5])
    sumwx[c(5,6)] <- W2%*%dataB2
    n0[c(1,2)]<-nn[4]*mix_pi[c(1,2)];n0[c(5,6)]<-nn[5]*mix_pi[c(5,6)]
    n0[c(1,2,5,6)][abs(n0[c(1,2,5,6)])<0.000001]<-0.000001
    s0[1]<-sumx[1]+sumx[2]+sumwx[1];s0[2]<-nn[1]+nn[2]+n0[1];s0[3]<-sumx[3]+sumwx[6]
    s0[4]<-nn[3]+n0[6];s0[5]<-sumwx[2]+sumwx[5];s0[6]<-n0[2]+n0[5]
    aa3<-(3*s0[1]/s0[2]-4*s0[5]/s0[6]+s0[3]/s0[4])/(9*sigma[11]/s0[2]+16*sigma[2]/s0[6]+sigma[11]/s0[4])
    m[11]<-(s0[1]-sigma[11]*aa3*3)/s0[2]            #mean1.
    m[12]<-m[11]                                  #mean2.
    m[13]<-(s0[3]-sigma[11]*aa3)/s0[4]           #mean3.
    m[2]<-(s0[5]+sigma[2]*aa3*4)/s0[6]
    m[1]<-m[11]; m[5]<-m[2];m[6]<-m[13]
    ############first order genetic parameters######
    hh<-matrix(c(1,1,1,1,-1,0.5),3,2)
    mm<-m[c(11,13,2)]
    B<-solve(crossprod(hh,hh))%*%crossprod(hh,mm)
    a1<-B[2];a1<-0.75*a1^2/num_l
    ########obtain variance#######################
    ss1<-sum((dataP1-m[11])^2);ss3<-sum((dataP2-m[13])^2);ss2<-sum((dataF1-m[12])^2)
    mm2<-m[c(5,6)]
    for(i in 1:d2) {swx1[i] <- W1[i,]%*%(dataB1-m[i])^2 } ;for(i in 1:d2) {swx2[i] <- W2[i,]%*%(dataB2-mm2[i])^2 }  
    swx[c(1,2)]<-swx1;swx[c(5,6)]<-swx2
    s0[1]<-ss1+ss2+ss3;s0[2]<-nn[1]+nn[2]+nn[3]
    aaa0<-sigma[11]; n_iter<-0;aa2<-1000
    while (aa2>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma[11]/(sigma[11]+a1)
      sigma[11]<-(s0[1]+(swx[1]+swx[6]+aa1^2*(swx[2]+swx[5])))/(s0[2]+n0[1]+n0[6]+aa1*(n0[2]+n0[5]))
      aa2<-abs(sigma[11]-aaa0)
      aaa0<-sigma[11]
      if (n_iter>20) break
    }
    sigma[1]<-sigma[11];sigma[6]<-sigma[11];sigma[2]<-sigma[1]+a1;sigma[5]<-sigma[2]
    mm1<-m[c(1,2)];ssigma1<-sigma[c(1,2)]
    mm2<-m[c(5,6)];ssigma2<-sigma[c(5,6)]
    mix_pi1<-mix_pi[c(1,2)];mix_pi2<-mix_pi[c(5,6)]
    ########criteria for iterations to stop########
    L1<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
      sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mix_pi2)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
    if(iteration>300)break
  }
  abc<-L0
  AIC<--2*abc+2*7
  ########second order genetic parameters###############   
  B1jj<-sigma1-sigma[11]
  if(B1jj<0){B1jj<-0}
  B1ll<-B1jj/sigma1
  B2jj<-sigma2-sigma[11]
  if(B2jj<0) {B2jj<-0}
  B2ll<-B2jj/sigma2   
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*nn[1])
  P1bmw <- matrix(0,nn[1],1)
  P1gg <- (dataP1 - m[11])/sqrt(as.vector(sigma[11]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn1 < nn[1]){P1bmw <- P1bmw+runif(nn[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:nn[1])) - 0.5)/nn[1])^2)
  P1u<- as.matrix(c(12*nn[1]*((P1dd[1]/nn[1]-0.5)^2),((45*nn[1])/4)*((P1dd[2]/nn[1]-1/3)^2),180*nn[1]*((P1dd[3]/nn[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,nn[1]))))
  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*nn[2])
  F1bmw <- matrix(0,nn[2],1)
  F1gg <- (dataF1 - m[12])/sqrt(as.vector(sigma[11]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn1 < nn[2]){F1bmw <- F1bmw+runif(nn[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:nn[2])) - 0.5)/nn[2])^2)
  F1u<- as.matrix(c(12*nn[2]*((F1dd[1]/nn[2]-0.5)^2),((45*nn[2])/4)*((F1dd[2]/nn[2]-1/3)^2),180*nn[2]*((F1dd[3]/nn[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,nn[2]))))
  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)
  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*nn[3])
  P2bmw <- matrix(0,nn[3],1)
  P2gg <- (dataP2 - m[13])/sqrt(as.vector(sigma[11]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn1 < nn[3]){P2bmw <- P2bmw+runif(nn[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:nn[3])) - 0.5)/nn[3])^2)
  P2u<- as.matrix(c(12*nn[3]*((P2dd[1]/nn[3]-0.5)^2),((45*nn[3])/4)*((P2dd[2]/nn[3]-1/3)^2),180*nn[3]*((P2dd[3]/nn[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,nn[3]))))
  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1); 
  B1w1<-1/(12*nn[4])
  B1bmw <- matrix(0,nn[4],1); B1bmwsl <- matrix(0,nn[4],d2)
  for(i in 1:d2){
    B1gg <- (dataB1 - m[i])/sqrt(sigma[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B1P2)))[1]
  if(nn1 < nn[4]){B1P2 <- B1P2+runif(nn[4])/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*nn[4]) + sum((B1P2 - (as.matrix(c(1:nn[4])) - 0.5)/nn[4])^2)
  B1u <- as.matrix(c(12*nn[4]*((B1dd[1]/nn[4]-0.5)^2),((45*nn[4])/4)*((B1dd[2]/nn[4]-1/3)^2),180*nn[4]*((B1dd[3]/nn[4]-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),(1-pkolm(B1D,nn[4]))))
  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2); 
  B2w1<-1/(12*nn[5])
  B2bmw <- matrix(0,nn[5],1); B2bmwsl <- matrix(0,nn[5],d2)
  for(i in 1:d2){
    B2gg <- (dataB2 - mm2[i])/sqrt(ssigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B2P2)))[1]
  if(nn1 < nn[5]){B2P2 <- B2P2+runif(nn[5])/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*nn[5]) + sum((B2P2 - (as.matrix(c(1:nn[5])) - 0.5)/nn[5])^2)
  B2u <- as.matrix(c(12*nn[5]*((B2dd[1]/nn[5]-0.5)^2),((45*nn[5])/4)*((B2dd[2]/nn[5]-1/3)^2),180*nn[5]*((B2dd[3]/nn[5]-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),(1-pkolm(B2D,nn[5]))))
  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("1MG-EAD",round(abc,4),round(AIC,4),round(m[11],4),round(m[12],4),round(m[13],4), round(m[1],4),round(m[2],4)," "," ",round(sigma[1],4),round(sigma[2],4),
                       " "," ",round(mix_pi[1],4),round(mix_pi[2],4)," "," ",round(m[5],4),round(m[6],4)," "," ",round(sigma[5],4),round(sigma[6],4)," "," ",        
                       round(mix_pi[5],4),round(mix_pi[6],4)," "," ",round(sigma[11],4),round(B[1],4),round(B[2],4)," "," "," "," "," ",round(B1jj,4),round(B1ll*100,4)," "," ",round(B2jj,4),round(B2ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,as.matrix(mmi1),as.matrix(mmi2))
  return(OUTPUT)
}

####################1MG-NCD(A-4)##################################
G5BCFModelFun[[4]] <- function(K1,logL,df11,df21,df31,df41,df51,G5BCFtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]))
  nn<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataB1)[1],dim(dataB2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2)))
  ss<-s[c(1:3)]-sumx[c(1:3)]^2/nn[c(1:3)]
  m<-matrix(0,13,1);m[11]<-mean(dataP1);m[12]<-mean(dataF1);m[13]<-mean(dataP2);man0<-mean(dataB1);man1<-mean(dataB2)
  sigma0<-sum(ss)/(nn[1]+nn[2]+nn[3]-3)
  sigma1<-var(dataB1);sigma2<-var(dataB2)
  d1<-1;d2<-2;d3<-3;d4<-4
  m_esp <- 0.0001;num_l <- as.numeric(G5BCFtext2)
  ###############procedure start###########################  
  mi<-matrix(0,6,1);mi[c(1,2,5,6)]<-0.5
  sigma<-matrix(0,11,1)
  sigma[1]<-sigma1/3;sigma[2]<-sigma[1];sigma[5]<-sigma2/3;sigma[6]<-sigma[5];sigma[11]<-sigma0
  a1<-sqrt(sigma1/(nn[4]-1))
  m[1]<-man0+a1;m[2]<-man0-a1
  a1<-sqrt(sigma2/(nn[5]-1))
  m[5]<-man1+a1;m[6]<-man1-a1
  mm1<-m[c(1,2)];mmi1<-mi[c(1,2)];ssigma1<-sigma[c(1,2)]
  mm2<-m[c(5,6)];mmi2<-mi[c(5,6)];ssigma2<-sigma[c(5,6)]
  L0<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
    sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)))+
    sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,2,nn[4]); swx1 <- matrix(0,2,1)
  W2 <- matrix(0,2,nn[5]); swx2 <- matrix(0,2,1)
  sumwx<-matrix(0,6,1);mix_pi<-matrix(0,6,1)
  n0<-matrix(0,6,1);s0<-matrix(0,6,1);swx<-matrix(0,6,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:2) { W1[i,] <- mi[i]*dnorm(dataB1,m[i],sqrt(sigma[i]))/dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)}
    mix_pi[c(1,2)] <- as.matrix(rowSums(W1)/nn[4])
    sumwx[c(1,2)] <- W1%*%dataB1
    for(i in 1:2) { W2[i,] <- mmi2[i]*dnorm(dataB2,mm2[i],sqrt(ssigma2[i]))/dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)}
    mix_pi[c(5,6)] <- as.matrix(rowSums(W2)/nn[5])
    sumwx[c(5,6)] <- W2%*%dataB2
    n0[c(1,2)]<-nn[4]*mix_pi[c(1,2)];n0[c(5,6)]<-nn[5]*mix_pi[c(5,6)]
    n0[c(1,2,5,6)][abs(n0[c(1,2,5,6)])<0.000001]<-0.000001
    s0[1]<-sumx[1]+sumwx[1];s0[2]<-nn[1]+n0[1];s0[3]<-sumx[2]+sumx[3]+sumwx[6]
    s0[4]<-nn[2]+nn[3]+n0[6];s0[5]<-sumwx[2]+sumwx[5];s0[6]<-n0[2]+n0[5]
    aa3<-(s0[1]/s0[2]-4*s0[5]/s0[6]+3*s0[3]/s0[4])/(sigma[11]/s0[2]+16*sigma[2]/s0[6]+9*sigma[11]/s0[4])
    m[11]<-(s0[1]-sigma[11]*aa3)/s0[2]                # mean1.
    m[12]<-(s0[3]-sigma[11]*aa3*3)/s0[4]            # mean2.
    m[13]<-m[12]                                   # mean3.
    m[2]<-(s0[5]+sigma[2]*aa3*4)/s0[6]
    m[1]<-m[11];m[5]<-m[2];m[6]<-m[13]
    ############first order genetic parameters######
    hh<-matrix(c(1,1,1,1,-1,-0.5),3,2)
    mm<-m[c(1,3,2)]
    B<-solve(crossprod(hh,hh))%*%crossprod(hh,mm)
    a1<-B[2];a1<-0.75*a1^2/num_l
    ########obtain variance#######################
    ss1<-sum((dataP1-m[11])^2);ss3<-sum((dataP2-m[13])^2);ss2<-sum((dataF1-m[12])^2)
    for(i in 1:d2) {swx1[i] <- W1[i,]%*%(dataB1-m[i])^2 } ;for(i in 1:d2) {swx2[i] <- W2[i,]%*%(dataB2-mm2[i])^2 }  
    swx[c(1,2)]<-swx1;swx[c(5,6)]<-swx2
    s0[1]<-ss1+ss2+ss3;s0[2]<-nn[1]+nn[2]+nn[3]
    aaa0<-sigma[11];n_iter<-0;aa2<-1000
    while (aa2>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma[11]/(sigma[11]+a1)
      sigma[11]<-(s0[1]+(swx[1]+swx[6]+aa1^2*(swx[2]+swx[5])))/(s0[2]+n0[1]+n0[6]+aa1*(n0[2]+n0[5]))
      aa2<-abs(sigma[11]-aaa0)
      aaa0<-sigma[11]
      if (n_iter>20) break
    }
    sigma[1]<-sigma[6]<-sigma[11];sigma[2]<-sigma[1]+a1;sigma[5]<-sigma[2]
    mm1<-m[c(1,2)];ssigma1<-sigma[c(1,2)]
    mm2<-m[c(5,6)];ssigma2<-sigma[c(5,6)]
    mix_pi1<-mix_pi[c(1,2)];mix_pi2<-mix_pi[c(5,6)]
    ########criteria for iterations to stop########
    L1<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
      sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mix_pi2)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
    if(iteration>300)break
  }
  
  abc<-L0
  AIC<--2*abc+2*7
  ########second order genetic parameters###############   
  B1jj<-sigma1-sigma[11]
  if(B1jj<0){B1jj<-0}
  B1ll<-B1jj/sigma1
  B2jj<-sigma2-sigma[11]
  if(B2jj<0) {B2jj<-0}
  B2ll<-B2jj/sigma2   
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*nn[1])
  P1bmw <- matrix(0,nn[1],1)
  P1gg <- (dataP1 - m[11])/sqrt(as.vector(sigma[11]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn1 < nn[1]){P1bmw <- P1bmw+runif(nn[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:nn[1])) - 0.5)/nn[1])^2)
  P1u<- as.matrix(c(12*nn[1]*((P1dd[1]/nn[1]-0.5)^2),((45*nn[1])/4)*((P1dd[2]/nn[1]-1/3)^2),180*nn[1]*((P1dd[3]/nn[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,nn[1]))))
  
  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*nn[2])
  F1bmw <- matrix(0,nn[2],1)
  F1gg <- (dataF1 - m[12])/sqrt(as.vector(sigma[11]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn1 < nn[2]){F1bmw <- F1bmw+runif(nn[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:nn[2])) - 0.5)/nn[2])^2)
  F1u<- as.matrix(c(12*nn[2]*((F1dd[1]/nn[2]-0.5)^2),((45*nn[2])/4)*((F1dd[2]/nn[2]-1/3)^2),180*nn[2]*((F1dd[3]/nn[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,nn[2]))))
  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)
  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*nn[3])
  P2bmw <- matrix(0,nn[3],1)
  P2gg <- (dataP2 - m[13])/sqrt(as.vector(sigma[11]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn1 < nn[3]){P2bmw <- P2bmw+runif(nn[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:nn[3])) - 0.5)/nn[3])^2)
  P2u<- as.matrix(c(12*nn[3]*((P2dd[1]/nn[3]-0.5)^2),((45*nn[3])/4)*((P2dd[2]/nn[3]-1/3)^2),180*nn[3]*((P2dd[3]/nn[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,nn[3]))))
  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1); 
  B1w1<-1/(12*nn[4])
  B1bmw <- matrix(0,nn[4],1); B1bmwsl <- matrix(0,nn[4],d2)
  for(i in 1:d2){
    B1gg <- (dataB1 - m[i])/sqrt(sigma[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B1P2)))[1]
  if(nn1 < nn[4]){B1P2 <- B1P2+runif(nn[4])/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*nn[4]) + sum((B1P2 - (as.matrix(c(1:nn[4])) - 0.5)/nn[4])^2)
  B1u <- as.matrix(c(12*nn[4]*((B1dd[1]/nn[4]-0.5)^2),((45*nn[4])/4)*((B1dd[2]/nn[4]-1/3)^2),180*nn[4]*((B1dd[3]/nn[4]-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),(1-pkolm(B1D,nn[4]))))
  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2); 
  B2w1<-1/(12*nn[5])
  B2bmw <- matrix(0,nn[5],1); B2bmwsl <- matrix(0,nn[5],d2)
  for(i in 1:d2){
    B2gg <- (dataB2 - mm2[i])/sqrt(ssigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B2P2)))[1]
  if(nn1 < nn[5]){B2P2 <- B2P2+runif(nn[5])/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*nn[5]) + sum((B2P2 - (as.matrix(c(1:nn[5])) - 0.5)/nn[5])^2)
  B2u <- as.matrix(c(12*nn[5]*((B2dd[1]/nn[5]-0.5)^2),((45*nn[5])/4)*((B2dd[2]/nn[5]-1/3)^2),180*nn[5]*((B2dd[3]/nn[5]-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),(1-pkolm(B2D,nn[5]))))
  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("1MG-NCD",round(abc,4),round(AIC,4),round(m[11],4),round(m[12],4),round(m[13],4), round(m[1],4),round(m[2],4)," "," ",round(sigma[1],4),round(sigma[2],4),
                       " "," ",round(mix_pi[1],4),round(mix_pi[2],4)," "," ",round(m[5],4),round(m[6],4)," "," ",round(sigma[5],4),round(sigma[6],4)," "," ",        
                       round(mix_pi[5],4),round(mix_pi[6],4)," "," ",round(sigma[11],4),round(B[1],4),round(B[2],4)," "," "," "," "," ",round(B1jj,4),round(B1ll*100,4)," "," ",round(B2jj,4),round(B2ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,as.matrix(mmi1),as.matrix(mmi2))
  return(OUTPUT)
}

################2MG-AD(B-2)##########################
G5BCFModelFun[[5]] <- function(K1,logL,df11,df21,df31,df41,df51,G5BCFtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]))
  nn<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataB1)[1],dim(dataB2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2)))
  ss<-s[c(1:3)]-sumx[c(1:3)]^2/nn[c(1:3)]
  m<-matrix(0,13,1);m[11]<-mean(dataP1);m[12]<-mean(dataF1);m[13]<-mean(dataP2);man0<-mean(dataB1);man1<-mean(dataB2)
  sigma0<-sum(ss)/(nn[1]+nn[2]+nn[3]-3)
  sigma1<-var(dataB1);sigma2<-var(dataB2)
  d1<-1;d2<-2;d3<-3;d4<-4
  m_esp <- 0.0001;num_l <- as.numeric(G5BCFtext2)
  ###############procedure start###########################  
  mi<-matrix(0.25,8,1)
  a1<-sqrt(sigma1/(nn[4]-1))
  m[1]<-man0+2*a1;m[2]<-man0+0.8*a1;m[3]<-man0-0.8*a1;m[4]<-man0-2*a1
  a1<-sqrt(sigma2/(nn[5]-1))
  m[5]<-man1+2*a1;m[6]<-man1+0.8*a1;m[7]<-man1-0.8*a1;m[8]<-man1-2*a1
  sigma<-matrix(0,11,1);sigma[11]<-sigma0;sigma[1]<-sigma1/3;sigma[8]<-sigma2/3
  hh1<-matrix(c(1,1,1,1,1,1,1,1,1,0,-1,1,0,0,0,-1,1,0,-1,0,1,0,-1,0,
                0,1,0,0,0.5,0.5,0.5,0,0,1,0,0.5,0,0.5,0,0.5),8,5)
  mm<-m[c(11,12,13,2,3,4,6,7)]
  B1<-solve(crossprod(hh1,hh1))%*%crossprod(hh1,mm)
  a1<-B1[2]    # da.
  a2<-B1[3]    # db.
  a3<-B1[4]    # ha.
  a4<-B1[5]    # hb.
  g<-matrix(0,3,1)
  g[1]<-(0.5*a2^2+0.25*a4^2)/num_l;g[2]<-(0.5*a1^2+0.25*a3^2)/num_l;g[3]<-g[1]+g[2]
  sigma[2]<-sigma[1]+g[1];sigma[3]<-sigma[1]+g[2];sigma[4]<-sigma[1]+g[3]
  sigma[5]<-sigma[8]+g[3];sigma[6]<-sigma[8]+g[2];sigma[7]<-sigma[8]+g[1]
  mm1<-m[c(1:4)];ssigma1<-sigma[c(1:4)];mmi1<-mi[c(1:4)]
  mm2<-m[c(5:8)];ssigma2<-sigma[c(5:8)];mmi2<-mi[c(5:8)]
  L0<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
    sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)))+
    sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,4,nn[4]); swx1 <- matrix(0,4,1)
  W2 <- matrix(0,4,nn[5]); swx2 <- matrix(0,4,1)
  sumwx<-matrix(0,8,1);mix_pi<-matrix(0,8,1)
  n0<-matrix(0,8,1);s0<-matrix(0,6,1);swx<-matrix(0,8,1);
  hh<-matrix(0,3,3);b_line<-matrix(0,3,1)
  aa3<-matrix(0,4,1); aa1<-matrix(0,8,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:4) { W1[i,] <- mi[i]*dnorm(dataB1,m[i],sqrt(sigma[i]))/dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)}
    mix_pi[c(1:4)] <- as.matrix(rowSums(W1)/nn[4])
    sumwx[c(1:4)] <- W1%*%dataB1
    for(i in 1:4) { W2[i,] <- mmi2[i]*dnorm(dataB2,mm2[i],sqrt(ssigma2[i]))/dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)}
    mix_pi[c(5:8)] <- as.matrix(rowSums(W2)/nn[5])
    sumwx[c(5:8)] <- W2%*%dataB2
    n0[c(1:4)]<-nn[4]*mix_pi[c(1:4)]; n0[c(5:8)]<-nn[5]*mix_pi[c(5:8)]
    n0[c(1:8)][abs(n0[c(1:8)])<0.000001]<-0.000001
    s0[1]<-sumx[1]+sumwx[1];s0[2]<-nn[1]+n0[1];s0[3]<-sumx[3]+sumwx[8]
    s0[4]<-nn[3]+n0[8];s0[5]<-sumwx[4]+sumwx[5];s0[6]<-n0[4]+n0[5]
    ##############################################
    hh[1,1]<-sigma[11]/s0[2]+4*sigma[11]/nn[2]+sigma[11]/s0[4]+16*sigma[4]/s0[6]
    hh[1,2]<-sigma[11]/s0[2]+4*sigma[11]/nn[2]+5*sigma[11]/s0[4]
    hh[1,3]<-sigma[11]/s0[2]-4*sigma[4]/s0[6]
    hh[2,2]<-sigma[11]/s0[2]+4*sigma[11]/nn[2]+25*sigma[11]/s0[4]+16*sigma[6]/n0[6]+16*sigma[7]/n0[7]
    hh[2,3]<-sigma[11]/s0[2]
    hh[3,3]<-sigma[11]/s0[2]+sigma[2]/n0[2]+sigma[3]/n0[3]+sigma[4]/s0[6]
    for(i in 2:3)
    {
      for(j in 1:(i-1))
      {
        hh[i,j]<-hh[j,i]
      }
    }
    ################################################
    b_line[1]<-s0[1]/s0[2]+2*sumx[2]/nn[2]+s0[3]/s0[4]-4*s0[5]/s0[6]
    b_line[2]<-s0[1]/s0[2]+2*sumx[2]/nn[2]+5*s0[3]/s0[4]-4*sumwx[6]/n0[6]-4*sumwx[7]/n0[7]
    b_line[3]<-s0[1]/s0[2]-sumwx[2]/n0[2]-sumwx[3]/n0[3]+s0[5]/s0[6]
    B2<-solve(hh,b_line) 
    m[11]<-(s0[1]-sigma[11]*(B2[1]+B2[2]+B2[3]))/s0[2]
    m[12]<-(sumx[2]-sigma[11]*(2*B2[1]+2*B2[2]))/nn[2]
    m[13]<-(s0[3]-sigma[11]*(B2[1]+5*B2[2]))/s0[4]
    m[2]<-(sumwx[2]+sigma[2]*B2[3])/n0[2]
    m[3]<-(sumwx[3]+sigma[3]*B2[3])/n0[3]
    m[4]<-(s0[5]+sigma[4]*(4*B2[1]-B2[3]))/s0[6]
    m[6]<-(sumwx[6]+sigma[6]*B2[2]*4)/n0[6]
    m[7]<-(sumwx[7]+sigma[7]*B2[2]*4)/n0[7]
    m[1]<-m[11]; m[5]<-m[4];m[8]<-m[13]
    #################################################
    hh3<-matrix(c(1,1,1,1,1,1,1,1,1,0,-1,1,0,0,0,-1,1,0,-1,0,1,0,-1,0,
                  0,1,0,0,0.5,0.5,0.5,0,0,1,0,0.5,0,0.5,0,0.5),8,5)
    mm<-m[c(11,12,13,2,3,4,6,7)]
    B3<-solve(crossprod(hh3,hh3))%*%crossprod(hh3,mm)
    a1<-B3[2]    # da.
    a2<-B3[3]    # db.
    a3<-B3[4]    # ha.
    a4<-B3[5]    # hb.
    g[1]<-(0.5*a2^2+0.25*a4^2)/num_l;g[2]<-(0.5*a1^2+0.25*a3^2)/num_l;g[3]<-g[1]+g[2]
    sigma[2]<-sigma[1]+g[1];sigma[3]<-sigma[1]+g[2];sigma[4]<-sigma[1]+g[3]
    sigma[5]<-sigma[8]+g[3];sigma[6]<-sigma[8]+g[2];sigma[7]<-sigma[8]+g[1]
    ########obtain variance#######################
    ss1<-sum((dataP1-m[11])^2);ss3<-sum((dataP2-m[13])^2);ss2<-sum((dataF1-m[12])^2)
    mm2<-m[c(5:8)]
    for(i in 1:d4) {swx1[i] <- W1[i,]%*%(dataB1-m[i])^2 } ;for(i in 1:d4) {swx2[i] <- W2[i,]%*%(dataB2-mm2[i])^2 }  
    swx[c(1:4)]<-swx1;swx[c(5:8)]<-swx2
    s0[1]<-ss1+ss2+ss3;s0[2]<-nn[1]+nn[2]+nn[3]
    aaa0<-sigma[11];n_iter<-0;abc2<-1000
    while (abc2>0.0001)
    {
      n_iter<-n_iter+1
      aa3[1]<-sigma[11]/sigma[11]
      aa3[c(2:4)]<-sigma[11]/(sigma[11]+g[c(1:3)])
      s0[3]<-sum(aa3[c(1:4)]^2*swx[c(1:4)])
      s0[4]<-sum(aa3[c(1:4)]*n0[c(1:4)])
      aa1[8]<-sigma[11]/sigma[11]
      aa1[c(5:7)]<-sigma[11]/(sigma[11]+g[c(3:1)])
      s0[5]<-sum(aa1[c(5:8)]^2*swx[c(5:8)])
      s0[6]<-sum(aa1[c(5:8)]*n0[c(5:8)])
      sigma[11]<-(s0[1]+(s0[3]+s0[5]))/(s0[2]+s0[4]+s0[6])
      abc2<-abs(sigma[11]-aaa0)
      aaa0<-sigma[11]
      if (n_iter>20) break
    }
    sigma[1]<-sigma[11];sigma[c(2:4)]<-sigma[11]+g[c(1:3)];sigma[c(5:7)]<-sigma[11]+g[c(3:1)];sigma[8]<-sigma[11]
    mm1<-m[c(1:4)];ssigma1<-sigma[c(1:4)]
    mm2<-m[c(5:8)];ssigma2<-sigma[c(5:8)]
    mix_pi1<-mix_pi[c(1:4)];mix_pi2<-mix_pi[c(5:8)]
    ########criteria for iterations to stop####################
    L1<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
      sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mix_pi2)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
    if(iteration>300)break
  }
  abc<-L0
  AIC<--2*abc+2*6
  ########first order genetic parameters###############
  hh4<-hh3
  mm<-m[c(11,12,13,2,3,4,6,7)]
  B4<-solve(t(hh4)%*%hh4)%*%(t(hh4)%*%mm)
  ########second order genetic parameters###############   
  B1jj<-sigma1-sigma[11]
  if(B1jj<0){B1jj<-0}
  B1ll<-B1jj/sigma1
  B2jj<-sigma2-sigma[11]
  if(B2jj<0) {B2jj<-0}
  B2ll<-B2jj/sigma2   
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*nn[1])
  P1bmw <- matrix(0,nn[1],1)
  P1gg <- (dataP1 - m[11])/sqrt(as.vector(sigma[11]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn1 < nn[1]){P1bmw <- P1bmw+runif(nn[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:nn[1])) - 0.5)/nn[1])^2)
  P1u<- as.matrix(c(12*nn[1]*((P1dd[1]/nn[1]-0.5)^2),((45*nn[1])/4)*((P1dd[2]/nn[1]-1/3)^2),180*nn[1]*((P1dd[3]/nn[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,nn[1]))))
  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*nn[2])
  F1bmw <- matrix(0,nn[2],1)
  F1gg <- (dataF1 - m[12])/sqrt(as.vector(sigma[11]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn1 < nn[2]){F1bmw <- F1bmw+runif(nn[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:nn[2])) - 0.5)/nn[2])^2)
  F1u<- as.matrix(c(12*nn[2]*((F1dd[1]/nn[2]-0.5)^2),((45*nn[2])/4)*((F1dd[2]/nn[2]-1/3)^2),180*nn[2]*((F1dd[3]/nn[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,nn[2]))))
  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)
  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*nn[3])
  P2bmw <- matrix(0,nn[3],1)
  P2gg <- (dataP2 - m[13])/sqrt(as.vector(sigma[11]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn1 < nn[3]){P2bmw <- P2bmw+runif(nn[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:nn[3])) - 0.5)/nn[3])^2)
  P2u<- as.matrix(c(12*nn[3]*((P2dd[1]/nn[3]-0.5)^2),((45*nn[3])/4)*((P2dd[2]/nn[3]-1/3)^2),180*nn[3]*((P2dd[3]/nn[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,nn[3]))))
  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1); 
  B1w1<-1/(12*nn[4])
  B1bmw <- matrix(0,nn[4],1); B1bmwsl <- matrix(0,nn[4],d4)
  for(i in 1:d4){
    B1gg <- (dataB1 - m[i])/sqrt(sigma[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B1P2)))[1]
  if(nn1 < nn[4]){B1P2 <- B1P2+runif(nn[4])/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*nn[4]) + sum((B1P2 - (as.matrix(c(1:nn[4])) - 0.5)/nn[4])^2)
  B1u <- as.matrix(c(12*nn[4]*((B1dd[1]/nn[4]-0.5)^2),((45*nn[4])/4)*((B1dd[2]/nn[4]-1/3)^2),180*nn[4]*((B1dd[3]/nn[4]-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),(1-pkolm(B1D,nn[4]))))
  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2); 
  B2w1<-1/(12*nn[5])
  B2bmw <- matrix(0,nn[5],1); B2bmwsl <- matrix(0,nn[5],d4)
  for(i in 1:d4){
    B2gg <- (dataB2 - mm2[i])/sqrt(ssigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B2P2)))[1]
  if(nn1 < nn[5]){B2P2 <- B2P2+runif(nn[5])/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*nn[5]) + sum((B2P2 - (as.matrix(c(1:nn[5])) - 0.5)/nn[5])^2)
  B2u <- as.matrix(c(12*nn[5]*((B2dd[1]/nn[5]-0.5)^2),((45*nn[5])/4)*((B2dd[2]/nn[5]-1/3)^2),180*nn[5]*((B2dd[3]/nn[5]-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),(1-pkolm(B2D,nn[5]))))
  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("2MG-AD",round(abc,4),round(AIC,4),round(m[11],4),round(m[12],4),round(m[13],4),round(m[1],4),round(m[2],4),round(m[3],4),round(m[4],4),round(sigma[1],4),round(sigma[2],4),
                       round(sigma[3],4),round(sigma[4],4),round(mix_pi[1],4),round(mix_pi[2],4),round(mix_pi[3],4),round(mix_pi[4],4),round(m[5],4),round(m[6],4),round(m[7],4),round(m[8],4),round(sigma[5],4),round(sigma[6],4), round(sigma[7],4),round(sigma[8],4),      
                       round(mix_pi[5],4),round(mix_pi[6],4),round(mix_pi[7],4),round(mix_pi[8],4),round(sigma[11],4),round(B4[1],4),round(B4[2],4),round(B4[3],4),round(B4[4],4),round(B4[5],4)," "," ",round(B1jj,4),round(B1ll*100,4)," "," ",round(B2jj,4),round(B2ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,as.matrix(mmi1),as.matrix(mmi2))
  return(OUTPUT)
}

#############2MG-A(B-3)###########################
G5BCFModelFun[[6]] <- function(K1,logL,df11,df21,df31,df41,df51,G5BCFtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]))
  nn<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataB1)[1],dim(dataB2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2)))
  ss<-s[c(1:3)]-sumx[c(1:3)]^2/nn[c(1:3)]
  m<-matrix(0,13,1);m[11]<-mean(dataP1);m[12]<-mean(dataF1);m[13]<-mean(dataP2);man0<-mean(dataB1);man1<-mean(dataB2)
  sigma0<-sum(ss)/(nn[1]+nn[2]+nn[3]-3)
  sigma1<-var(dataB1);sigma2<-var(dataB2)
  d1<-1;d2<-2;d3<-3;d4<-4
  m_esp <- 0.0001;num_l <- as.numeric(G5BCFtext2)
  ###############procedure start###########################  
  mi<-matrix(0.25,8,1)
  a1<-sqrt(sigma1/(nn[4]-1))
  m[1]<-man0+2*a1;m[2]<-man0+0.8*a1;m[3]<-man0-0.8*a1;m[4]<-man0-2*a1
  a1<-sqrt(sigma2/(nn[5]-1))
  m[5]<-man1+2*a1;m[6]<-man1+0.8*a1;m[7]<-man1-0.8*a1;m[8]<-man1-2*a1
  sigma<-matrix(0,11,1)
  sigma[11]<-sigma0;sigma[1]<-sigma1/3;sigma[8]<-sigma2/3
  hh1<-matrix(c(1,1,1,1,1,1,1,1,1,0,-1,1,0,0,0,-1,1,0,-1,0,1,0,-1,0),8,3)
  mm<-m[c(11,12,13,2,3,4,6,7)]
  B1<-solve(crossprod(hh1,hh1))%*%crossprod(hh1,mm)
  a1<-B1[2]    # da.
  a2<-B1[3]    # db.
  g<-matrix(0,3,1)
  g[1]<-0.5*a2^2/num_l;g[2]<-0.5*a1^2/num_l;g[3]<-g[1]+g[2]
  sigma[2]<-sigma[1]+g[1];sigma[3]<-sigma[1]+g[2];sigma[4]<-sigma[1]+g[3]
  sigma[5]<-sigma[8]+g[3];sigma[6]<-sigma[8]+g[2];sigma[7]<-sigma[8]+g[1]
  mm1<-m[c(1:4)];ssigma1<-sigma[c(1:4)];mmi1<-mi[c(1:4)]
  mm2<-m[c(5:8)];ssigma2<-sigma[c(5:8)];mmi2<-mi[c(5:8)]
  L0<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
    sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)))+
    sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,4,nn[4]); swx1 <- matrix(0,4,1)
  W2 <- matrix(0,4,nn[5]); swx2 <- matrix(0,4,1)
  sumwx<-matrix(0,8,1);mix_pi<-matrix(0,8,1)
  n0<-matrix(0,8,1);s0<-matrix(0,6,1)
  hh<-matrix(0,4,4);b_line<-matrix(0,4,1);swx<-matrix(0,8,1)
  aa3<-matrix(0,4,1);aa1<-matrix(0,8,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:4) { W1[i,] <- mi[i]*dnorm(dataB1,m[i],sqrt(sigma[i]))/dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)}
    mix_pi[c(1:4)] <- as.matrix(rowSums(W1)/nn[4])
    sumwx[c(1:4)] <- W1%*%dataB1
    for(i in 1:4) { W2[i,] <- mmi2[i]*dnorm(dataB2,mm2[i],sqrt(ssigma2[i]))/dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)}
    mix_pi[c(5:8)] <- as.matrix(rowSums(W2)/nn[5])
    sumwx[c(5:8)] <- W2%*%dataB2
    n0[c(1:4)]<-nn[4]*mix_pi[c(1:4)];n0[c(5:8)]<-nn[5]*mix_pi[c(5:8)]
    n0[c(1:8)][abs(n0[c(1:8)])<0.000001]<-0.000001
    s0[1]<-sumx[1]+sumwx[1];s0[2]<-nn[1]+n0[1];s0[3]<-sumx[2]+(sumwx[4]+sumwx[5])*(sigma[11]/sigma[4])
    s0[4]<-nn[2]+(n0[4]+n0[5])*(sigma[11]/sigma[4]);s0[5]<-sumx[3]+sumwx[8];s0[6]<-nn[3]+n0[8]
    ################################################
    hh[1,1]<-9*sigma[11]/s0[2]+sigma[11]/s0[6]+4*sigma[2]/n0[2]+4*sigma[3]/n0[3]
    hh[1,2]<-3*sigma[11]/s0[2]+3*sigma[11]/s0[6]
    hh[1,3]<-3*sigma[11]/s0[2]+sigma[11]/s0[6]
    hh[1,4]<--2*sigma[2]/n0[2]+2*sigma[3]/n0[3]
    hh[2,2]<-sigma[11]/s0[2]+9*sigma[11]/s0[6]+4*sigma[6]/n0[6]+4*sigma[7]/n0[7]
    hh[2,3]<-sigma[11]/s0[2]+3*sigma[11]/s0[6]
    hh[2,4]<-2*sigma[6]/n0[6]-2*sigma[7]/n0[7]
    hh[3,3]<-sigma[11]/s0[2]+4*sigma[11]/s0[4]+sigma[11]/s0[6]
    hh[3,4]<-0
    hh[4,4]<-sigma[2]/n0[2]+sigma[3]/n0[3]+sigma[6]/n0[6]+sigma[7]/n0[7]
    for(i in 2:4)
    {
      for(j in 1:(i-1))
      {
        hh[i,j]<-hh[j,i]
      }
    }
    ############################################
    b_line[1]<-3*s0[1]/s0[2]+s0[5]/s0[6]-2*sumwx[2]/n0[2]-2*sumwx[3]/n0[3]
    b_line[2]<-s0[1]/s0[2]+3*s0[5]/s0[6]-2*sumwx[6]/n0[6]-2*sumwx[7]/n0[7]
    b_line[3]<-s0[1]/s0[2]-2*s0[3]/s0[4]+s0[5]/s0[6]
    b_line[4]<-sumwx[2]/n0[2]-sumwx[3]/n0[3]-sumwx[6]/n0[6]+sumwx[7]/n0[7]
    B2<-solve(hh,b_line)
    #############################################
    m[11]<-(s0[1]-sigma[11]*(3*B2[1]+B2[2]+B2[3]))/s0[2]
    m[12]<-(s0[3]+sigma[11]*2*B2[3])/s0[4]
    m[13]<-(s0[5]-sigma[11]*(B2[1]+3*B2[2]+B2[3]))/s0[6]
    m[2]<-(sumwx[2]+sigma[2]*(2*B2[1]-B2[4]))/n0[2]
    m[3]<-(sumwx[3]+sigma[3]*(2*B2[1]+B2[4]))/n0[3]
    m[6]<-(sumwx[6]+sigma[6]*(2*B2[2]+B2[4]))/n0[6]
    m[7]<-(sumwx[7]+sigma[7]*(2*B2[2]-B2[4]))/n0[7]
    m[1]<-m[11];  m[4]<-m[12];m[5]<-m[12];m[8]<-m[13]
    ##############################################
    hh3<-hh1
    mm<-m[c(11,12,13,2,3,4,6,7)]
    B3<-solve(crossprod(hh3,hh3))%*%crossprod(hh3,mm)
    a1<-B3[2]    # da.
    a2<-B3[3]    # db.
    g[1]<-0.5*a2^2/num_l;g[2]<-0.5*a1^2/num_l;g[3]<-g[1]+g[2]
    sigma[2]<-sigma[1]+g[1];sigma[3]<-sigma[1]+g[2];sigma[4]<-sigma[1]+g[3]
    sigma[5]<-sigma[8]+g[3];sigma[6]<-sigma[8]+g[2];sigma[7]<-sigma[8]+g[1]
    ########obtain variance#######################
    ss1<-sum((dataP1-m[11])^2);ss3<-sum((dataP2-m[13])^2);ss2<-sum((dataF1-m[12])^2)
    mm2<-m[c(5:8)]
    for(i in 1:d4) {swx1[i] <- W1[i,]%*%(dataB1-m[i])^2 } ;for(i in 1:d4) {swx2[i] <- W2[i,]%*%(dataB2-mm2[i])^2 }  
    swx[c(1:4)]<-swx1;swx[c(5:8)]<-swx2
    s0[1]<-ss1+ss2+ss3;s0[2]<-nn[1]+nn[2]+nn[3]
    aaa0<-sigma[11];n_iter<-0;abc2<-1000
    while (abc2>0.0001)
    {
      n_iter<-n_iter+1
      aa3[1]<-sigma[11]/sigma[11]
      aa3[c(2:4)]<-sigma[11]/(sigma[11]+g[c(1:3)])
      s0[3]<-sum(aa3[c(1:4)]^2*swx[c(1:4)])
      s0[4]<-sum(aa3[c(1:4)]*n0[c(1:4)])
      aa1[8]<-sigma[11]/sigma[11]
      aa1[c(5:7)]<-sigma[11]/(sigma[11]+g[c(3:1)])
      s0[5]<-sum(aa1[c(5:8)]^2*swx[c(5:8)])
      s0[6]<-sum(aa1[c(5:8)]*n0[c(5:8)])
      sigma[11]<-(s0[1]+(s0[3]+s0[5]))/(s0[2]+s0[4]+s0[6])
      abc2<-abs(sigma[11]-aaa0)
      aaa0<-sigma[11]
      if (n_iter>20) break
    }
    sigma[1]<-sigma[11];sigma[c(2:4)]<-sigma[11]+g[c(1:3)];sigma[c(5:7)]<-sigma[11]+g[c(3:1)];sigma[8]<-sigma[11]
    mm1<-m[c(1:4)];ssigma1<-sigma[c(1:4)]
    mm2<-m[c(5:8)];ssigma2<-sigma[c(5:8)]
    mix_pi1<-mix_pi[c(1:4)];mix_pi2<-mix_pi[c(5:8)]
    ########criteria for iterations to stop########
    L1<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
      sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mix_pi2)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
    if(iteration>300)break
  }
  abc<-L0
  AIC<--2*abc+2*4
  ########first order genetic parameters###############
  hh4<-hh1
  mm<-m[c(11,12,13,2,3,4,6,7)]
  B4<-solve(t(hh4)%*%hh4)%*%(t(hh4)%*%mm)
  ########second order genetic parameters###############   
  B1jj<-sigma1-sigma[11]
  if(B1jj<0){B1jj<-0}
  B1ll<-B1jj/sigma1
  B2jj<-sigma2-sigma[11]
  if(B2jj<0) {B2jj<-0}
  B2ll<-B2jj/sigma2   
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*nn[1])
  P1bmw <- matrix(0,nn[1],1)
  P1gg <- (dataP1 - m[11])/sqrt(as.vector(sigma[11]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn1 < nn[1]){P1bmw <- P1bmw+runif(nn[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:nn[1])) - 0.5)/nn[1])^2)
  P1u<- as.matrix(c(12*nn[1]*((P1dd[1]/nn[1]-0.5)^2),((45*nn[1])/4)*((P1dd[2]/nn[1]-1/3)^2),180*nn[1]*((P1dd[3]/nn[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,nn[1]))))
  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*nn[2])
  F1bmw <- matrix(0,nn[2],1)
  F1gg <- (dataF1 - m[12])/sqrt(as.vector(sigma[11]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn1 < nn[2]){F1bmw <- F1bmw+runif(nn[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:nn[2])) - 0.5)/nn[2])^2)
  F1u<- as.matrix(c(12*nn[2]*((F1dd[1]/nn[2]-0.5)^2),((45*nn[2])/4)*((F1dd[2]/nn[2]-1/3)^2),180*nn[2]*((F1dd[3]/nn[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,nn[2]))))
  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)
  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*nn[3])
  P2bmw <- matrix(0,nn[3],1)
  P2gg <- (dataP2 - m[13])/sqrt(as.vector(sigma[11]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn1 < nn[3]){P2bmw <- P2bmw+runif(nn[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:nn[3])) - 0.5)/nn[3])^2)
  P2u<- as.matrix(c(12*nn[3]*((P2dd[1]/nn[3]-0.5)^2),((45*nn[3])/4)*((P2dd[2]/nn[3]-1/3)^2),180*nn[3]*((P2dd[3]/nn[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,nn[3]))))
  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1); 
  B1w1<-1/(12*nn[4])
  B1bmw <- matrix(0,nn[4],1); B1bmwsl <- matrix(0,nn[4],d4)
  for(i in 1:d4){
    B1gg <- (dataB1 - m[i])/sqrt(sigma[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B1P2)))[1]
  if(nn1 < nn[4]){B1P2 <- B1P2+runif(nn[4])/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*nn[4]) + sum((B1P2 - (as.matrix(c(1:nn[4])) - 0.5)/nn[4])^2)
  B1u <- as.matrix(c(12*nn[4]*((B1dd[1]/nn[4]-0.5)^2),((45*nn[4])/4)*((B1dd[2]/nn[4]-1/3)^2),180*nn[4]*((B1dd[3]/nn[4]-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),(1-pkolm(B1D,nn[4]))))
  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2); 
  B2w1<-1/(12*nn[5])
  B2bmw <- matrix(0,nn[5],1); B2bmwsl <- matrix(0,nn[5],d4)
  for(i in 1:d4){
    B2gg <- (dataB2 - mm2[i])/sqrt(ssigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B2P2)))[1]
  if(nn1 < nn[5]){B2P2 <- B2P2+runif(nn[5])/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*nn[5]) + sum((B2P2 - (as.matrix(c(1:nn[5])) - 0.5)/nn[5])^2)
  B2u <- as.matrix(c(12*nn[5]*((B2dd[1]/nn[5]-0.5)^2),((45*nn[5])/4)*((B2dd[2]/nn[5]-1/3)^2),180*nn[5]*((B2dd[3]/nn[5]-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),(1-pkolm(B2D,nn[5]))))
  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("2MG-A",round(abc,4),round(AIC,4),round(m[11],4),round(m[12],4),round(m[13],4),round(m[1],4),round(m[2],4),round(m[3],4),round(m[4],4),round(sigma[1],4),round(sigma[2],4),
                       round(sigma[3],4),round(sigma[4],4),round(mix_pi[1],4),round(mix_pi[2],4),round(mix_pi[3],4),round(mix_pi[4],4),round(m[5],4),round(m[6],4),round(m[7],4),round(m[8],4),round(sigma[5],4),round(sigma[6],4), round(sigma[7],4),round(sigma[8],4),      
                       round(mix_pi[5],4),round(mix_pi[6],4),round(mix_pi[7],4),round(mix_pi[8],4),round(sigma[11],4),round(B4[1],4),round(B4[2],4),round(B4[3],4)," "," "," "," ",round(B1jj,4),round(B1ll*100,4)," "," ",round(B2jj,4),round(B2ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,as.matrix(mmi1),as.matrix(mmi2))
  return(OUTPUT)
}

###############2MG-EA(B-4)##############################################
G5BCFModelFun[[7]] <- function(K1,logL,df11,df21,df31,df41,df51,G5BCFtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]))
  nn<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataB1)[1],dim(dataB2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2)))
  ss<-s[c(1:3)]-sumx[c(1:3)]^2/nn[c(1:3)]
  m<-matrix(0,13,1);m[11]<-mean(dataP1);m[12]<-mean(dataF1);m[13]<-mean(dataP2);man0<-mean(dataB1);man1<-mean(dataB2)
  sigma0<-sum(ss)/(nn[1]+nn[2]+nn[3]-3)
  sigma1<-var(dataB1);sigma2<-var(dataB2)
  d1<-1;d2<-2;d3<-3;d4<-4
  m_esp <- 0.0001;num_l <- as.numeric(G5BCFtext2)
  ###############procedure start###########################  
  mi<-matrix(0,7,1);mi[c(1,3)]<-0.25;mi[c(2,6)]<-0.5;mi[c(5,7)]<-0.25
  a1<-sqrt(sigma1/(nn[4]-1))
  m[1]<-man0+2*a1;m[2]<-man0;m[3]<-man0-2*a1
  a1<-sqrt(sigma2/(nn[5]-1))
  m[5]<-man1+2*a1;m[6]<-man1;m[7]<-man1-2*a1
  sigma<-matrix(0,11,1)
  sigma[11]<-sigma0;sigma[1]<-sigma1/3;sigma[c(7,8)]<-sigma2/3
  hh1<-matrix(c(1,1,1,1,1,2,0,-2,1,-1),5,2)
  mm<-m[c(11,12,13,2,6)]
  B1<-solve(crossprod(hh1,hh1))%*%crossprod(hh1,mm)
  a1<-B1[2]    # d.
  g<-matrix(0,3,1)
  g[1]<-0.5*a1^2/num_l;g[2]<-a1^2/num_l
  sigma[2]<-sigma[1]+g[1];sigma[3]<-sigma[1]+g[2]
  sigma[5]<-sigma[7]+g[2];sigma[6]<-sigma[7]+g[1]
  mm1<-m[c(1:3)];ssigma1<-sigma[c(1:3)];mmi1<-mi[c(1:3)]
  mm2<-m[c(5:7)];ssigma2<-sigma[c(5:7)];mmi2<-mi[c(5:7)]
  L0<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
    sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)))+
    sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,3,nn[4]); swx1 <- matrix(0,3,1)
  W2 <- matrix(0,3,nn[5]); swx2 <- matrix(0,3,1)
  sumwx<-matrix(0,8,1);mix_pi<-matrix(0,8,1);swx<-matrix(0,8,1);
  n0<-matrix(0,8,1); s0<-matrix(0,6,1)
  hh<-matrix(0,3,3);b_line<-matrix(0,3,1)
  aa3<-matrix(0,3,1);aa1<-matrix(0,7,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:3) { W1[i,] <- mi[i]*dnorm(dataB1,m[i],sqrt(sigma[i]))/dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)}
    mix_pi[c(1:3)] <- as.matrix(rowSums(W1)/nn[4])
    sumwx[c(1:3)] <- W1%*%dataB1
    for(i in 1:3) { W2[i,] <- mmi2[i]*dnorm(dataB2,mm2[i],sqrt(ssigma2[i]))/dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)}
    mix_pi[c(5:7)] <- as.matrix(rowSums(W2)/nn[5])
    sumwx[c(5:7)] <- W2%*%dataB2
    n0[c(1:3)]<-nn[4]*mix_pi[c(1:3)];n0[c(5:7)]<-nn[5]*mix_pi[c(5:7)]
    n0[c(1,2,3,5,6,7)][abs(n0[c(1,2,3,5,6,7)])<0.000001]<-0.000001
    s0[1]<-sumx[1]+sumwx[1];s0[2]<-nn[1]+n0[1];s0[3]<-sumx[2]*sigma[3]+(sumwx[3]+sumwx[5])*sigma[11]
    s0[4]<-nn[2]*sigma[3]+(n0[3]+n0[5])*sigma[11];s0[5]<-sumx[3]+sumwx[7];s0[6]<-nn[3]+n0[7]
    sigma[9]<-sigma[3]*sigma[11]
    #########################################
    hh[1,1]<-sigma[11]/s0[2]+4*sigma[9]/s0[4]+sigma[11]/s0[6]
    hh[1,2]<--4*sigma[9]/s0[4]
    hh[1,3]<-sigma[11]/s0[2]-sigma[11]/s0[6]
    hh[2,2]<-4*sigma[9]/s0[4]+sigma[2]/n0[2]+sigma[6]/n0[6]
    hh[2,3]<-2*sigma[2]/n0[2]-2*sigma[6]/n0[6]
    hh[3,3]<-sigma[11]/s0[2]+sigma[11]/s0[6]+4*sigma[2]/n0[2]+4*sigma[6]/n0[6]
    for(i in 2:3)
    {
      for(j in 1:(i-1))
      {
        hh[i,j]<-hh[j,i]
      }
    }
    ###########################################
    b_line[1]<-s0[1]/s0[2]-2*s0[3]/s0[4]+s0[5]/s0[6]
    b_line[2]<-2*s0[3]/s0[4]-sumwx[2]/n0[2]-sumwx[6]/n0[6]
    b_line[3]<-s0[1]/s0[2]-s0[5]/s0[6]-2*sumwx[2]/n0[2]+2*sumwx[6]/n0[6]
    B2<-solve(hh,b_line)
    #############################################
    m[11]<-(s0[1]-sigma[11]*(B2[1]+B2[3]))/s0[2]
    m[12]<-(s0[3]+sigma[9]*(2*B2[1]-2*B2[2]))/s0[4]
    m[13]<-(s0[5]-sigma[11]*(B2[1]-B2[3]))/s0[6]
    m[2]<-(sumwx[2]+sigma[2]*(B2[2]+2*B2[3]))/n0[2]
    m[6]<-(sumwx[6]+sigma[6]*(B2[2]-2*B2[3]))/n0[6]
    m[1]<-m[11];m[3]<-m[12]; m[5]<-m[12];m[7]<-m[13]
    ##############################################
    hh3<-hh1
    mm<-m[c(11,12,13,2,6)]
    B3<-solve(crossprod(hh3,hh3))%*%crossprod(hh3,mm)
    a1<-B3[2]    # d.
    g[1]<-0.5*a1^2/num_l;g[2]<-a1^2/num_l
    sigma[2]<-sigma[1]+g[1];sigma[3]<-sigma[1]+g[2];sigma[5]<-sigma[7]+g[2];sigma[6]<-sigma[7]+g[1]
    ########obtain variance#######################
    ss1<-sum((dataP1-m[11])^2);ss3<-sum((dataP2-m[13])^2);ss2<-sum((dataF1-m[12])^2)
    mm2<-m[c(5:7)]
    for(i in 1:d3) {swx1[i] <- W1[i,]%*%(dataB1-m[i])^2 } ;for(i in 1:d3) {swx2[i] <- W2[i,]%*%(dataB2-mm2[i])^2 }  
    swx[c(1:3)]<-swx1;swx[c(5:7)]<-swx2
    s0[1]<-ss1+ss2+ss3;s0[2]<-nn[1]+nn[2]+nn[3]
    aaa0<-sigma[11];n_iter<-0;abc2<-1000
    while (abc2>0.0001)
    {
      n_iter<-n_iter+1
      aa3[1]<-sigma[11]/sigma[11]
      aa3[c(2:3)]<-sigma[11]/(sigma[11]+g[c(1:2)])
      s0[3]<-sum(aa3[c(1:3)]^2*swx[c(1:3)])
      s0[4]<-sum(aa3[c(1:3)]*n0[c(1:3)])
      aa1[7]<-sigma[11]/sigma[11]
      aa1[c(5:6)]<-sigma[11]/(sigma[11]+g[c(2:1)])
      s0[5]<-sum(aa1[c(5:7)]^2*swx[c(5:7)])
      s0[6]<-sum(aa1[c(5:7)]*n0[c(5:7)])
      sigma[11]<-(s0[1]+(s0[3]+s0[5]))/(s0[2]+s0[4]+s0[6])
      abc2<-abs(sigma[11]-aaa0)
      aaa0<-sigma[11]
      if (n_iter>20) break
    }
    sigma[1]<-sigma[11];sigma[c(2:3)]<-sigma[11]+g[c(1:2)];sigma[c(5:6)]<-sigma[11]+g[c(2:1)];sigma[7]<-sigma[11]
    mm1<-m[c(1:3)];ssigma1<-sigma[c(1:3)]
    mm2<-m[c(5:7)];ssigma2<-sigma[c(5:7)]
    mix_pi1<-mix_pi[c(1:3)];mix_pi2<-mix_pi[c(5:7)]
    ########criteria for iterations to stop############
    L1<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
      sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mix_pi2)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
    if(iteration>300)break
  }
  abc<-L0
  AIC<--2*abc+2*3
  ########first order genetic parameters###############
  hh4<-hh1
  mm<-m[c(11,12,13,2,6)]
  B4<-solve(t(hh4)%*%hh4)%*%(t(hh4)%*%mm)
  ########second order genetic parameters###############   
  B1jj<-sigma1-sigma[11]
  if(B1jj<0){B1jj<-0}
  B1ll<-B1jj/sigma1
  B2jj<-sigma2-sigma[11]
  if(B2jj<0) {B2jj<-0}
  B2ll<-B2jj/sigma2   
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*nn[1])
  P1bmw <- matrix(0,nn[1],1)
  P1gg <- (dataP1 - m[11])/sqrt(as.vector(sigma[11]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn1 < nn[1]){P1bmw <- P1bmw+runif(nn[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:nn[1])) - 0.5)/nn[1])^2)
  P1u<- as.matrix(c(12*nn[1]*((P1dd[1]/nn[1]-0.5)^2),((45*nn[1])/4)*((P1dd[2]/nn[1]-1/3)^2),180*nn[1]*((P1dd[3]/nn[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,nn[1]))))
  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*nn[2])
  F1bmw <- matrix(0,nn[2],1)
  F1gg <- (dataF1 - m[12])/sqrt(as.vector(sigma[11]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn1 < nn[2]){F1bmw <- F1bmw+runif(nn[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:nn[2])) - 0.5)/nn[2])^2)
  F1u<- as.matrix(c(12*nn[2]*((F1dd[1]/nn[2]-0.5)^2),((45*nn[2])/4)*((F1dd[2]/nn[2]-1/3)^2),180*nn[2]*((F1dd[3]/nn[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,nn[2]))))
  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)
  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*nn[3])
  P2bmw <- matrix(0,nn[3],1)
  P2gg <- (dataP2 - m[13])/sqrt(as.vector(sigma[11]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn1 < nn[3]){P2bmw <- P2bmw+runif(nn[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:nn[3])) - 0.5)/nn[3])^2)
  P2u<- as.matrix(c(12*nn[3]*((P2dd[1]/nn[3]-0.5)^2),((45*nn[3])/4)*((P2dd[2]/nn[3]-1/3)^2),180*nn[3]*((P2dd[3]/nn[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,nn[3]))))
  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1); 
  B1w1<-1/(12*nn[4])
  B1bmw <- matrix(0,nn[4],1); B1bmwsl <- matrix(0,nn[4],d3)
  for(i in 1:d3){
    B1gg <- (dataB1 - m[i])/sqrt(sigma[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B1P2)))[1]
  if(nn1 < nn[4]){B1P2 <- B1P2+runif(nn[4])/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*nn[4]) + sum((B1P2 - (as.matrix(c(1:nn[4])) - 0.5)/nn[4])^2)
  B1u <- as.matrix(c(12*nn[4]*((B1dd[1]/nn[4]-0.5)^2),((45*nn[4])/4)*((B1dd[2]/nn[4]-1/3)^2),180*nn[4]*((B1dd[3]/nn[4]-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),(1-pkolm(B1D,nn[4]))))
  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2); 
  B2w1<-1/(12*nn[5])
  B2bmw <- matrix(0,nn[5],1); B2bmwsl <- matrix(0,nn[5],d3)
  for(i in 1:d3){
    B2gg <- (dataB2 - mm2[i])/sqrt(ssigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B2P2)))[1]
  if(nn1 < nn[5]){B2P2 <- B2P2+runif(nn[5])/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*nn[5]) + sum((B2P2 - (as.matrix(c(1:nn[5])) - 0.5)/nn[5])^2)
  B2u <- as.matrix(c(12*nn[5]*((B2dd[1]/nn[5]-0.5)^2),((45*nn[5])/4)*((B2dd[2]/nn[5]-1/3)^2),180*nn[5]*((B2dd[3]/nn[5]-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),(1-pkolm(B2D,nn[5]))))
  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("2MG-EA",round(abc,4),round(AIC,4),round(m[11],4),round(m[12],4),round(m[13],4),round(m[1],4),round(m[2],4),round(m[3],4)," ",round(sigma[1],4),round(sigma[2],4),
                       round(sigma[3],4)," ",round(mix_pi[1],4),round(mix_pi[2],4),round(mix_pi[3],4)," ",round(m[5],4),round(m[6],4),round(m[7],4)," ",round(sigma[5],4),round(sigma[6],4), round(sigma[7],4)," ",      
                       round(mix_pi[5],4),round(mix_pi[6],4),round(mix_pi[7],4)," ",round(sigma[11],4),round(B4[1],4),round(B4[2],4)," "," "," "," "," ",round(B1jj,4),round(B1ll*100,4)," "," ",round(B2jj,4),round(B2ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,as.matrix(mmi1),as.matrix(mmi2))
  return(OUTPUT)
}

################2MG-CD(B-5)##########################
G5BCFModelFun[[8]] <- function(K1,logL,df11,df21,df31,df41,df51,G5BCFtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]))
  nn<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataB1)[1],dim(dataB2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2)))
  ss<-s[c(1:3)]-sumx[c(1:3)]^2/nn[c(1:3)]
  m<-matrix(0,13,1);m[11]<-mean(dataP1);m[12]<-mean(dataF1);m[13]<-mean(dataP2);man0<-mean(dataB1);man1<-mean(dataB2)
  sigma0<-sum(ss)/(nn[1]+nn[2]+nn[3]-3)
  sigma1<-var(dataB1);sigma2<-var(dataB2)
  d1<-1;d2<-2;d3<-3;d4<-4
  m_esp <- 0.0001;num_l <- as.numeric(G5BCFtext2)
  ###############procedure start###########################  
  mi<-matrix(0.25,8,1);a1<-sqrt(sigma1/(nn[4]-1))
  m[1]<-man0+2*a1;m[2]<-man0+0.8*a1;m[3]<-man0-0.8*a1;m[4]<-man0-2*a1
  a1<-sqrt(sigma2/(nn[5]-1))
  m[5]<-man1+2*a1;m[6]<-man1+0.8*a1;m[7]<-man1-0.8*a1;m[8]<-man1-2*a1
  sigma<-matrix(0,11,1);sigma[11]<-sigma0;sigma[1]<-sigma1/3;sigma[8]<-sigma2/3
  hh1<-matrix(c(1,1,1,1,1,1,1,1,-1,1,0.5,0.5,0.5,-1,1,-1,0.5,1,0.5,-1,0.5),7,3)
  mm<-m[c(11,13,2,3,4,6,7)]
  B1<-solve(crossprod(hh1,hh1))%*%crossprod(hh1,mm)
  a1<-B1[2]    # da.
  a2<-B1[3]    # db.
  g<-matrix(0,3,1)
  g[1]<-0.75*a2^2/num_l;g[2]<-0.75*a1^2/num_l;g[3]<-g[1]+g[2]
  sigma[2]<-sigma[1]+g[1];sigma[3]<-sigma[1]+g[2];sigma[4]<-sigma[1]+g[3]
  sigma[5]<-sigma[8]+g[3];sigma[6]<-sigma[8]+g[2];sigma[7]<-sigma[8]+g[1]
  mm1<-m[c(1:4)];ssigma1<-sigma[c(1:4)];mmi1<-mi[c(1:4)]
  mm2<-m[c(5:8)];ssigma2<-sigma[c(5:8)];mmi2<-mi[c(5:8)]
  L0<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
    sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)))+
    sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,4,nn[4]); swx1 <- matrix(0,4,1)
  W2 <- matrix(0,4,nn[5]); swx2 <- matrix(0,4,1)
  sumwx<-matrix(0,8,1);mix_pi<-matrix(0,8,1)
  n0<-matrix(0,8,1); s0<-matrix(0,6,1)
  hh<-matrix(0,4,4);b_line<-matrix(0,4,1)
  swx<-matrix(0,8,1);aa3<-matrix(0,4,1);aa1<-matrix(0,8,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:4) { W1[i,] <- mi[i]*dnorm(dataB1,m[i],sqrt(sigma[i]))/dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)}
    mix_pi[c(1:4)] <- as.matrix(rowSums(W1)/nn[4])
    sumwx[c(1:4)] <- W1%*%dataB1
    for(i in 1:4) { W2[i,] <- mmi2[i]*dnorm(dataB2,mm2[i],sqrt(ssigma2[i]))/dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)}
    mix_pi[c(5:8)] <- as.matrix(rowSums(W2)/nn[5])
    sumwx[c(5:8)] <- W2%*%dataB2
    n0[c(1:4)]<-nn[4]*mix_pi[c(1:4)];n0[c(5:8)]<-nn[5]*mix_pi[c(5:8)]
    n0[c(1:8)][abs(n0[c(1:8)])<0.000001]<-0.000001
    s0[1]<-sumx[1]+sumx[2]+sumwx[1];s0[2]<-nn[1]+nn[2]+n0[1]
    s0[3]<-sumx[3]+sumwx[8];s0[4]<-nn[3]+n0[8]
    s0[5]<-sumwx[4]+sumwx[5];s0[6]<-n0[4]+n0[5]
    ########################################
    hh[1,1]<-9*sigma[11]/s0[2]+sigma[11]/s0[4]+16*sigma[4]/s0[6]
    hh[1,2]<-9*sigma[11]/s0[2]+5*sigma[11]/s0[4]
    hh[1,3]<-3*sigma[11]/s0[2]-4*sigma[4]/s0[6]
    hh[1,4]<-0
    hh[2,2]<-9*sigma[11]/s0[2]+25*sigma[11]/s0[4]+16*sigma[6]/n0[6]+16*sigma[7]/n0[7]
    hh[2,3]<-3*sigma[11]/s0[2]
    hh[2,4]<-4*sigma[6]/n0[6]-4*sigma[7]/n0[7]
    hh[3,3]<-sigma[11]/s0[2]+sigma[2]/n0[2]+sigma[3]/n0[3]+sigma[4]/s0[6]
    hh[3,4]<--3*sigma[2]/n0[2]+3*sigma[3]/n0[3]
    hh[4,4]<-9*sigma[2]/n0[2]+9*sigma[3]/n0[3]+sigma[6]/n0[6]+sigma[7]/n0[7]
    for(i in 2:4)
    {
      for(j in 1:(i-1))
      {
        hh[i,j]<-hh[j,i]
      }
    }
    ###############################################
    b_line[1]<-3*s0[1]/s0[2]+s0[3]/s0[4]-4*s0[5]/s0[6]
    b_line[2]<-3*s0[1]/s0[2]+5*s0[3]/s0[4]-4*sumwx[6]/n0[6]-4*sumwx[7]/n0[7]
    b_line[3]<-s0[1]/s0[2]-sumwx[2]/n0[2]-sumwx[3]/n0[3]+s0[5]/s0[6]
    b_line[4]<-3*sumwx[2]/n0[2]-3*sumwx[3]/n0[3]-sumwx[6]/n0[6]+sumwx[7]/n0[7] 
    B2<-solve(hh,b_line)
    ###############################################
    m[11]<-(s0[1]-sigma[11]*(3*B2[1]+3*B2[2]+B2[3]))/s0[2]
    m[13]<-(s0[3]-sigma[11]*(B2[1]+5*B2[2]))/s0[4]
    m[2]<-(sumwx[2]+sigma[2]*(B2[3]-3*B2[4]))/n0[2]
    m[3]<-(sumwx[3]+sigma[3]*(B2[3]+3*B2[4]))/n0[3]
    m[4]<-(s0[5]+sigma[4]*(4*B2[1]-B2[3]))/s0[6]
    m[6]<-(sumwx[6]+sigma[6]*(B2[2]*4+B2[4]))/n0[6]
    m[7]<-(sumwx[7]+sigma[7]*(B2[2]*4-B2[4]))/n0[7]
    m[12]<-m[11];m[1]<-m[11];m[5]<-m[4];m[8]<-m[13]
    ##############################################
    hh3<-hh1
    mm<-m[c(11,13,2,3,4,6,7)]
    B3<-solve(crossprod(hh3,hh3))%*%crossprod(hh3,mm)
    a1<-B3[2]    # da.
    a2<-B3[3]    # db.
    g[1]<-0.75*a2^2/num_l;g[2]<-0.75*a1^2/num_l;g[3]<-g[1]+g[2]
    sigma[2]<-sigma[1]+g[1];sigma[3]<-sigma[1]+g[2];sigma[4]<-sigma[1]+g[3]
    sigma[5]<-sigma[8]+g[3];sigma[6]<-sigma[8]+g[2];sigma[7]<-sigma[8]+g[1]
    ########obtain variance#######################
    ss1<-sum((dataP1-m[11])^2);ss3<-sum((dataP2-m[13])^2);ss2<-sum((dataF1-m[12])^2)
    mm2<-m[c(5:8)]
    for(i in 1:d4) {swx1[i] <- W1[i,]%*%(dataB1-m[i])^2 } ;for(i in 1:d4) {swx2[i] <- W2[i,]%*%(dataB2-mm2[i])^2 }  
    swx[c(1:4)]<-swx1;swx[c(5:8)]<-swx2
    s0[1]<-ss1+ss2+ss3;s0[2]<-nn[1]+nn[2]+nn[3]
    aaa0<-sigma[11];n_iter<-0;abc2<-1000
    while (abc2>0.0001)
    {
      n_iter<-n_iter+1
      aa3[1]<-sigma[11]/sigma[11]
      aa3[c(2:4)]<-sigma[11]/(sigma[11]+g[c(1:3)])
      s0[3]<-sum(aa3[c(1:4)]^2*swx[c(1:4)])
      s0[4]<-sum(aa3[c(1:4)]*n0[c(1:4)])
      aa1[8]<-sigma[11]/sigma[11]
      aa1[c(5:7)]<-sigma[11]/(sigma[11]+g[c(3:1)])
      s0[5]<-sum(aa1[c(5:8)]^2*swx[c(5:8)])
      s0[6]<-sum(aa1[c(5:8)]*n0[c(5:8)])
      sigma[11]<-(s0[1]+(s0[3]+s0[5]))/(s0[2]+s0[4]+s0[6])
      abc2<-abs(sigma[11]-aaa0)
      aaa0<-sigma[11]
      if (n_iter>20) break
    }
    sigma[1]<-sigma[11];sigma[c(2:4)]<-sigma[11]+g[c(1:3)];sigma[c(5:7)]<-sigma[11]+g[c(3:1)];sigma[8]<-sigma[11]
    mm1<-m[c(1:4)];ssigma1<-sigma[c(1:4)]
    mm2<-m[c(5:8)];ssigma2<-sigma[c(5:8)]
    mix_pi1<-mix_pi[c(1:4)];mix_pi2<-mix_pi[c(5:8)]
    ########criteria for iterations to stop########
    L1<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
      sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mix_pi2)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
    if(iteration>300)break
  }
  abc<-L0
  AIC<--2*abc+2*4
  ########first order genetic parameters###############
  hh4<-hh1
  mm<-m[c(11,13,2,3,4,6,7)]
  B4<-solve(t(hh4)%*%hh4)%*%(t(hh4)%*%mm)
  ########second order genetic parameters###############   
  B1jj<-sigma1-sigma[11]
  if(B1jj<0){B1jj<-0}
  B1ll<-B1jj/sigma1
  B2jj<-sigma2-sigma[11]
  if(B2jj<0) {B2jj<-0}
  B2ll<-B2jj/sigma2   
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*nn[1])
  P1bmw <- matrix(0,nn[1],1)
  P1gg <- (dataP1 - m[11])/sqrt(as.vector(sigma[11]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn1 < nn[1]){P1bmw <- P1bmw+runif(nn[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:nn[1])) - 0.5)/nn[1])^2)
  P1u<- as.matrix(c(12*nn[1]*((P1dd[1]/nn[1]-0.5)^2),((45*nn[1])/4)*((P1dd[2]/nn[1]-1/3)^2),180*nn[1]*((P1dd[3]/nn[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,nn[1]))))
  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*nn[2])
  F1bmw <- matrix(0,nn[2],1)
  F1gg <- (dataF1 - m[12])/sqrt(as.vector(sigma[11]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn1 < nn[2]){F1bmw <- F1bmw+runif(nn[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:nn[2])) - 0.5)/nn[2])^2)
  F1u<- as.matrix(c(12*nn[2]*((F1dd[1]/nn[2]-0.5)^2),((45*nn[2])/4)*((F1dd[2]/nn[2]-1/3)^2),180*nn[2]*((F1dd[3]/nn[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,nn[2]))))
  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)
  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*nn[3])
  P2bmw <- matrix(0,nn[3],1)
  P2gg <- (dataP2 - m[13])/sqrt(as.vector(sigma[11]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn1 < nn[3]){P2bmw <- P2bmw+runif(nn[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:nn[3])) - 0.5)/nn[3])^2)
  P2u<- as.matrix(c(12*nn[3]*((P2dd[1]/nn[3]-0.5)^2),((45*nn[3])/4)*((P2dd[2]/nn[3]-1/3)^2),180*nn[3]*((P2dd[3]/nn[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,nn[3]))))
  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1); 
  B1w1<-1/(12*nn[4])
  B1bmw <- matrix(0,nn[4],1); B1bmwsl <- matrix(0,nn[4],d4)
  for(i in 1:d4){
    B1gg <- (dataB1 - m[i])/sqrt(sigma[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B1P2)))[1]
  if(nn1 < nn[4]){B1P2 <- B1P2+runif(nn[4])/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*nn[4]) + sum((B1P2 - (as.matrix(c(1:nn[4])) - 0.5)/nn[4])^2)
  B1u <- as.matrix(c(12*nn[4]*((B1dd[1]/nn[4]-0.5)^2),((45*nn[4])/4)*((B1dd[2]/nn[4]-1/3)^2),180*nn[4]*((B1dd[3]/nn[4]-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),(1-pkolm(B1D,nn[4]))))
  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2); 
  B2w1<-1/(12*nn[5])
  B2bmw <- matrix(0,nn[5],1); B2bmwsl <- matrix(0,nn[5],d4)
  for(i in 1:d4){
    B2gg <- (dataB2 - mm2[i])/sqrt(ssigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B2P2)))[1]
  if(nn1 < nn[5]){B2P2 <- B2P2+runif(nn[5])/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*nn[5]) + sum((B2P2 - (as.matrix(c(1:nn[5])) - 0.5)/nn[5])^2)
  B2u <- as.matrix(c(12*nn[5]*((B2dd[1]/nn[5]-0.5)^2),((45*nn[5])/4)*((B2dd[2]/nn[5]-1/3)^2),180*nn[5]*((B2dd[3]/nn[5]-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),(1-pkolm(B2D,nn[5]))))
  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("2MG-CD",round(abc,4),round(AIC,4),round(m[11],4),round(m[12],4),round(m[13],4),round(m[1],4),round(m[2],4),round(m[3],4),round(m[4],4),round(sigma[1],4),round(sigma[2],4),
                       round(sigma[3],4),round(sigma[4],4),round(mix_pi[1],4),round(mix_pi[2],4),round(mix_pi[3],4),round(mix_pi[4],4),round(m[5],4),round(m[6],4),round(m[7],4),round(m[8],4),round(sigma[5],4),round(sigma[6],4), round(sigma[7],4),round(sigma[8],4),      
                       round(mix_pi[5],4),round(mix_pi[6],4),round(mix_pi[7],4),round(mix_pi[8],4),round(sigma[11],4),round(B4[1],4),round(B4[2],4),round(B4[3],4)," "," "," "," ",round(B1jj,4),round(B1ll*100,4)," "," ",round(B2jj,4),round(B2ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,as.matrix(mmi1),as.matrix(mmi2))
  return(OUTPUT)
}

############2MG-EAD(B-6)######################
G5BCFModelFun[[9]] <- function(K1,logL,df11,df21,df31,df41,df51,G5BCFtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]))
  nn<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataB1)[1],dim(dataB2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2)))
  ss<-s[c(1:3)]-sumx[c(1:3)]^2/nn[c(1:3)]
  m<-matrix(0,13,1);m[11]<-mean(dataP1);m[12]<-mean(dataF1);m[13]<-mean(dataP2);man0<-mean(dataB1);man1<-mean(dataB2)
  sigma0<-sum(ss)/(nn[1]+nn[2]+nn[3]-3)
  sigma1<-var(dataB1);sigma2<-var(dataB2)
  d1<-1;d2<-2;d3<-3;d4<-4
  m_esp <- 0.0001;num_l <- as.numeric(G5BCFtext2)
  ###############procedure start###########################  
  mi<-matrix(0,7,1)
  mi[c(1,3)]<-0.25;mi[c(2,6)]<-0.5;mi[c(5,7)]<-0.25
  a1<-sqrt(sigma1/(nn[4]-1))
  m[1]<-man0+2*a1;m[2]<-man0;m[3]<-man0-2*a1
  a1<-sqrt(sigma2/(nn[5]-1))
  m[5]<-man1+2*a1;m[6]<-man1;m[7]<-man1-2*a1
  sigma<-matrix(0,11,1)
  sigma[11]<-sigma0;sigma[1]<-sigma1/3;sigma[c(7,8)]<-sigma2/3
  hh1<-matrix(c(1,1,1,1,1,2,-2,1.5,1,-0.5),5,2)
  mm<-m[c(11,13,2,3,6)]
  B1<-solve(crossprod(hh1,hh1))%*%crossprod(hh1,mm)
  a1<-B1[2]    # d.
  g<-matrix(0,2,1)
  g[1]<-0.75*a1^2/num_l;g[2]<-1.5*a1^2/num_l
  sigma[2]<-sigma[1]+g[1];sigma[3]<-sigma[1]+g[2]
  sigma[5]<-sigma[7]+g[2];sigma[6]<-sigma[7]+g[1]
  mm1<-m[c(1:3)];ssigma1<-sigma[c(1:3)];mmi1<-mi[c(1:3)]
  mm2<-m[c(5:7)];ssigma2<-sigma[c(5:7)];mmi2<-mi[c(5:7)]
  L0<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
    sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)))+
    sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,3,nn[4]); swx1 <- matrix(0,3,1)
  W2 <- matrix(0,3,nn[5]); swx2 <- matrix(0,3,1)
  sumwx<-matrix(0,8,1);mix_pi<-matrix(0,8,1)
  n0<-matrix(0,8,1);s0<-matrix(0,6,1)
  hh<-matrix(0,3,3);b_line<-matrix(0,3,1)
  swx<-matrix(0,8,1);aa3<-matrix(0,3,1);aa1<-matrix(0,7,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:3) { W1[i,] <- mi[i]*dnorm(dataB1,m[i],sqrt(sigma[i]))/dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)}
    mix_pi[c(1:3)] <- as.matrix(rowSums(W1)/nn[4])
    sumwx[c(1:3)] <- W1%*%dataB1
    for(i in 1:3) { W2[i,] <- mmi2[i]*dnorm(dataB2,mm2[i],sqrt(ssigma2[i]))/dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)}
    mix_pi[c(5:7)] <- as.matrix(rowSums(W2)/nn[5])
    sumwx[c(5:7)] <- W2%*%dataB2
    n0[c(1:3)]<-nn[4]*mix_pi[c(1:3)];n0[c(5:7)]<-nn[5]*mix_pi[c(5:7)]
    n0[c(1,2,3,5,6,7)][abs(n0[c(1,2,3,5,6,7)])<0.000001]<-0.000001
    s0[1]<-sumx[1]+sumx[2]+sumwx[1];s0[2]<-nn[1]+nn[2]+n0[1];s0[3]<-sumx[3]+sumwx[7]
    s0[4]<-nn[3]+n0[7];s0[5]<-sumwx[3]+sumwx[5];s0[6]<-n0[3]+n0[5]
    ###########################################################
    hh[1,1]<-49*sigma[11]/s0[2]+sigma[11]/s0[4]+64*sigma[2]/n0[2]
    hh[1,2]<-21*sigma[11]/s0[2]+5*sigma[11]/s0[4]
    hh[1,3]<--24*sigma[2]/n0[2]
    hh[2,1]<-hh[1,2]
    hh[2,2]<-9*sigma[11]/s0[2]+25*sigma[11]/s0[4]+64*sigma[6]/n0[6]
    hh[2,3]<--8*sigma[6]/n0[6]
    hh[3,1]<-hh[1,3]
    hh[3,2]<-hh[2,3]
    hh[3,3]<-9*sigma[2]/n0[2]+16*sigma[3]/s0[6]+sigma[6]/n0[6]
    for(i in 2:3)
    {
      for(j in 1:(i-1))
      {
        hh[i,j]<-hh[j,i]
      }
    }
    #########################################
    b_line[1]<-7*s0[1]/s0[2]+s0[3]/s0[4]-8*sumwx[2]/n0[2]
    b_line[2]<-3*s0[1]/s0[2]+5*s0[3]/s0[4]-8*sumwx[6]/n0[6]
    b_line[3]<-3*sumwx[2]/n0[2]-4*s0[5]/s0[6]+sumwx[6]/n0[6]
    B2<-solve(hh,b_line)
    #############################################
    m[11]<-(s0[1]-sigma[11]*(7*B2[1]+3*B2[2]))/s0[2]
    m[13]<-(s0[3]-sigma[11]*(B2[1]+5*B2[2]))/s0[4]
    m[2]<-(sumwx[2]+sigma[2]*(8*B2[1]-3*B2[3]))/n0[2]
    m[3]<-(s0[5]+sigma[3]*4*B2[3])/s0[6]
    m[6]<-(sumwx[6]+sigma[6]*(8*B2[2]-B2[3]))/n0[6]
    m[12]<-m[11];m[1]<-m[11]; m[5]<-m[3];m[7]<-m[13]
    ##############################################
    hh3<-hh1
    mm<-m[c(11,13,2,3,6)]
    B3<-solve(crossprod(hh3,hh3))%*%crossprod(hh3,mm)
    a1<-B3[2]    # d.
    g[1]<-0.75*a1^2/num_l;g[2]<-1.5*a1^2/num_l
    sigma[2]<-sigma[1]+g[1];sigma[3]<-sigma[1]+g[2]
    sigma[5]<-sigma[7]+g[2];sigma[6]<-sigma[7]+g[1]
    ########obtain variance#######################
    ss1<-sum((dataP1-m[11])^2);ss3<-sum((dataP2-m[13])^2);ss2<-sum((dataF1-m[12])^2)
    mm2<-m[c(5:7)]
    for(i in 1:d3) {swx1[i] <- W1[i,]%*%(dataB1-m[i])^2 } ;for(i in 1:d3) {swx2[i] <- W2[i,]%*%(dataB2-mm2[i])^2 }  
    swx[c(1:3)]<-swx1;swx[c(5:7)]<-swx2
    s0[1]<-ss1+ss2+ss3;s0[2]<-nn[1]+nn[2]+nn[3]
    aaa0<-sigma[11];n_iter<-0;abc2<-1000
    while (abc2>0.0001)
    {
      n_iter<-n_iter+1
      aa3[1]<-sigma[11]/sigma[11]
      aa3[c(2:3)]<-sigma[11]/(sigma[11]+g[c(1:2)])
      s0[3]<-sum(aa3[c(1:3)]^2*swx[c(1:3)])
      s0[4]<-sum(aa3[c(1:3)]*n0[c(1:3)])
      aa1[7]<-sigma[11]/sigma[11]
      aa1[c(5:6)]<-sigma[11]/(sigma[11]+g[c(2:1)])
      s0[5]<-sum(aa1[c(5:7)]^2*swx[c(5:7)])
      s0[6]<-sum(aa1[c(5:7)]*n0[c(5:7)])
      sigma[11]<-(s0[1]+(s0[3]+s0[5]))/(s0[2]+s0[4]+s0[6])
      abc2<-abs(sigma[11]-aaa0)
      aaa0<-sigma[11]
      if (n_iter>20) break
    }
    sigma[1]<-sigma[11];sigma[c(2:3)]<-sigma[11]+g[c(1:2)];sigma[c(5:6)]<-sigma[11]+g[c(2:1)];sigma[7]<-sigma[11]
    mm1<-m[c(1:3)];ssigma1<-sigma[c(1:3)]
    mm2<-m[c(5:7)];ssigma2<-sigma[c(5:7)]
    mix_pi1<-mix_pi[c(1:3)];mix_pi2<-mix_pi[c(5:7)]
    ########criteria for iterations to stop########
    L1<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
      sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mix_pi2)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
    if(iteration>300)break
  }
  abc<-L0
  AIC<--2*abc+2*3
  ########first order genetic parameters###############
  hh4<-hh1
  mm<-m[c(11,13,2,3,6)]
  B4<-solve(t(hh4)%*%hh4)%*%(t(hh4)%*%mm)
  ########second order genetic parameters###############   
  B1jj<-sigma1-sigma[11]
  if(B1jj<0){B1jj<-0}
  B1ll<-B1jj/sigma1
  B2jj<-sigma2-sigma[11]
  if(B2jj<0) {B2jj<-0}
  B2ll<-B2jj/sigma2   
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*nn[1])
  P1bmw <- matrix(0,nn[1],1)
  P1gg <- (dataP1 - m[11])/sqrt(as.vector(sigma[11]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn1 < nn[1]){P1bmw <- P1bmw+runif(nn[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:nn[1])) - 0.5)/nn[1])^2)
  P1u<- as.matrix(c(12*nn[1]*((P1dd[1]/nn[1]-0.5)^2),((45*nn[1])/4)*((P1dd[2]/nn[1]-1/3)^2),180*nn[1]*((P1dd[3]/nn[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,nn[1]))))
  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*nn[2])
  F1bmw <- matrix(0,nn[2],1)
  F1gg <- (dataF1 - m[12])/sqrt(as.vector(sigma[11]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn1 < nn[2]){F1bmw <- F1bmw+runif(nn[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:nn[2])) - 0.5)/nn[2])^2)
  F1u<- as.matrix(c(12*nn[2]*((F1dd[1]/nn[2]-0.5)^2),((45*nn[2])/4)*((F1dd[2]/nn[2]-1/3)^2),180*nn[2]*((F1dd[3]/nn[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,nn[2]))))
  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)
  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*nn[3])
  P2bmw <- matrix(0,nn[3],1)
  P2gg <- (dataP2 - m[13])/sqrt(as.vector(sigma[11]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn1 < nn[3]){P2bmw <- P2bmw+runif(nn[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:nn[3])) - 0.5)/nn[3])^2)
  P2u<- as.matrix(c(12*nn[3]*((P2dd[1]/nn[3]-0.5)^2),((45*nn[3])/4)*((P2dd[2]/nn[3]-1/3)^2),180*nn[3]*((P2dd[3]/nn[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,nn[3]))))
  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1); 
  B1w1<-1/(12*nn[4])
  B1bmw <- matrix(0,nn[4],1); B1bmwsl <- matrix(0,nn[4],d3)
  for(i in 1:d3){
    B1gg <- (dataB1 - m[i])/sqrt(sigma[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B1P2)))[1]
  if(nn1 < nn[4]){B1P2 <- B1P2+runif(nn[4])/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*nn[4]) + sum((B1P2 - (as.matrix(c(1:nn[4])) - 0.5)/nn[4])^2)
  B1u <- as.matrix(c(12*nn[4]*((B1dd[1]/nn[4]-0.5)^2),((45*nn[4])/4)*((B1dd[2]/nn[4]-1/3)^2),180*nn[4]*((B1dd[3]/nn[4]-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),(1-pkolm(B1D,nn[4]))))
  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2); 
  B2w1<-1/(12*nn[5])
  B2bmw <- matrix(0,nn[5],1); B2bmwsl <- matrix(0,nn[5],d3)
  for(i in 1:d3){
    B2gg <- (dataB2 - mm2[i])/sqrt(ssigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B2P2)))[1]
  if(nn1 < nn[5]){B2P2 <- B2P2+runif(nn[5])/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*nn[5]) + sum((B2P2 - (as.matrix(c(1:nn[5])) - 0.5)/nn[5])^2)
  B2u <- as.matrix(c(12*nn[5]*((B2dd[1]/nn[5]-0.5)^2),((45*nn[5])/4)*((B2dd[2]/nn[5]-1/3)^2),180*nn[5]*((B2dd[3]/nn[5]-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),(1-pkolm(B2D,nn[5]))))
  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("2MG-EAD",round(abc,4),round(AIC,4),round(m[11],4),round(m[12],4),round(m[13],4),round(m[1],4),round(m[2],4),round(m[3],4)," ",round(sigma[1],4),round(sigma[2],4),
                       round(sigma[3],4)," ",round(mix_pi[1],4),round(mix_pi[2],4),round(mix_pi[3],4)," ",round(m[5],4),round(m[6],4),round(m[7],4)," ",round(sigma[5],4),round(sigma[6],4), round(sigma[7],4)," ",      
                       round(mix_pi[5],4),round(mix_pi[6],4),round(mix_pi[7],4)," ",round(sigma[11],4),round(B4[1],4),round(B4[2],4)," "," "," "," "," ",round(B1jj,4),round(B1ll*100,4)," "," ",round(B2jj,4),round(B2ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,as.matrix(mmi1),as.matrix(mmi2))
  return(OUTPUT)
}

####################MX1-AD-AD(D-1)############################################
G5BCFModelFun[[10]] <- function(K1,logL,df11,df21,df31,df41,df51,G5BCFtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]))
  nn<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataB1)[1],dim(dataB2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2)))
  ss<-s[c(1:3)]-sumx[c(1:3)]^2/nn[c(1:3)]
  m<-matrix(0,13,1);m[11]<-mean(dataP1);m[12]<-mean(dataF1);m[13]<-mean(dataP2);man0<-mean(dataB1);man1<-mean(dataB2)
  sigma0<-sum(ss)/(nn[1]+nn[2]+nn[3]-3)
  sigma1<-var(dataB1);sigma2<-var(dataB2)
  d1<-1;d2<-2;d3<-3;d4<-4
  m_esp <- 0.0001;num_l <- as.numeric(G5BCFtext2)
  ###############procedure start###########################  
  mi<-matrix(0,6,1);mi[c(1,2,5,6)]<-0.5
  sigma<-matrix(0,11,1)
  sigma[1]<-sigma1/3;sigma[2]<-sigma[1];sigma[5]<-sigma2/3;sigma[6]<-sigma[5];sigma[11]<-sigma0
  a1<-sqrt(sigma1/(nn[4]-1))
  m[1]<-man0+a1;m[2]<-man0-a1
  a1<-sqrt(sigma2/(nn[5]-1))
  m[5]<-man1+a1;m[6]<-man1-a1
  mm1<-m[c(1,2)];mmi1<-mi[c(1,2)];ssigma1<-sigma[c(1,2)]
  mm2<-m[c(5,6)];mmi2<-mi[c(5,6)];ssigma2<-sigma[c(5,6)]
  L0<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
    sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)))+
    sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,2,nn[4]); swx1 <- matrix(0,2,1)
  W2 <- matrix(0,2,nn[5]); swx2 <- matrix(0,2,1)
  sumwx<-matrix(0,6,1);mix_pi<-matrix(0,6,1)
  n0<-matrix(0,6,1);s0<-matrix(0,2,1);  rr<-matrix(0,2,1)
  swx<-matrix(0,6,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:2) { W1[i,] <- mi[i]*dnorm(dataB1,m[i],sqrt(sigma[i]))/dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)}
    mix_pi[c(1,2)] <- as.matrix(rowSums(W1)/nn[4])
    sumwx[c(1,2)] <- W1%*%dataB1
    for(i in 1:2) { W2[i,] <- mmi2[i]*dnorm(dataB2,mm2[i],sqrt(ssigma2[i]))/dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)}
    mix_pi[c(5,6)] <- as.matrix(rowSums(W2)/nn[5])
    sumwx[c(5,6)] <- W2%*%dataB2
    n0[c(1,2)]<-nn[4]*mix_pi[c(1,2)];n0[c(5,6)]<-nn[5]*mix_pi[c(5,6)]
    n0[c(1,2,5,6)][abs(n0[c(1,2,5,6)])<0.000001]<-0.000001
    aa1<-25*sigma[11]/nn[1]+4*sigma[11]/nn[2]+sigma[11]/nn[3]+16*sigma[1]/n0[1]+16*sigma[2]/n0[2]
    aa2<-5*sigma[11]/nn[1]+4*sigma[11]/nn[2]+5*sigma[11]/nn[3]
    aa3<-sigma[11]/nn[1]+4*sigma[11]/nn[2]+25*sigma[11]/nn[3]+16*sigma[5]/n0[5]+16*sigma[6]/n0[6]
    aa4<-aa1*aa3-aa2^2
    s0[1]<-5*sumx[1]/nn[1]+2*sumx[2]/nn[2]+sumx[3]/nn[3]-4*sumwx[1]/n0[1]-4*sumwx[2]/n0[2]
    s0[2]<-sumx[1]/nn[1]+2*sumx[2]/nn[2]+5*sumx[3]/nn[3]-4*sumwx[5]/n0[5]-4*sumwx[6]/n0[6]
    aa5<-s0[1]*aa3-s0[2]*aa2
    aa6<-aa1*s0[2]-s0[1]*aa2
    rr[1]<-aa5/aa4;rr[2]<-aa6/aa4
    m[11]<-(sumx[1]-sigma[11]*(5*rr[1]+rr[2]))/nn[1]
    m[12]<-(sumx[2]-sigma[11]*(2*rr[1]+rr[2]*2))/nn[2]
    m[13]<-(sumx[3]-sigma[11]*(rr[1]+5*rr[2]))/nn[3]
    m[1]<-(sumwx[1]+sigma[1]*rr[1]*4)/n0[1]
    m[2]<-(sumwx[2]+sigma[2]*rr[1]*4)/n0[2]
    m[5]<-(sumwx[5]+sigma[5]*rr[2]*4)/n0[5]
    m[6]<-(sumwx[6]+sigma[6]*rr[2]*4)/n0[6]
    hh1<-matrix(c(1,1,1,1,1,1,1,1,0,-1,1,0,0,-1,0,1,0,0,0.5,0.5,0,1,0,-1,0.5,0.5,-0.5,-0.5,
                  0,1,0,0.5,0.5,0.25,0.25),7,5)
    mm<-m[c(11,12,13,1,2,5,6)]
    B<-solve(crossprod(hh1,hh1))%*%crossprod(hh1,mm)
    a1<-B[2];a2<-B[3]
    a1<-(0.5*a1^2+0.25*a2^2)/num_l
    ########obtain variance#######################
    ss1<-sum((dataP1-m[11])^2);ss3<-sum((dataP2-m[13])^2);ss2<-sum((dataF1-m[12])^2)
    mm2<-m[c(5,6)]
    for(i in 1:d2) {swx1[i] <- W1[i,]%*%(dataB1-m[i])^2 } ;for(i in 1:d2) {swx2[i] <- W2[i,]%*%(dataB2-mm2[i])^2 }  
    swx[c(1,2)]<-swx1;swx[c(5,6)]<-swx2
    n_iter<-0;aaa0<-sigma[1];aa2<-1000
    while (aa2>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma[1]/(sigma[1]+a1)
      sigma[1]<-(swx[1]+aa1^2*swx[2])/(n0[1]+aa1*n0[2])
      aa2<-abs(sigma[1]-aaa0)
      aaa0<-sigma[1]
      if (n_iter>20) break
    }
    sigma40<-sigma[1]-sigma[11]
    if (sigma40<0) {sigma40<-0;sigma[1]<-sigma[11]}
    sigma[2]<-sigma[1]+a1
    n_iter<-0;aaa0<-sigma[6];aa2<-1000
    while(aa2>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma[6]/(sigma[6]+a1)
      sigma[6]<-(swx[6]+aa1^2*swx[5])/(n0[6]+aa1*n0[5])
      aa2<-abs(sigma[6]-aaa0)
      aaa0<-sigma[6]
      if(n_iter>20) break
    }
    sigma50<-sigma[6]-sigma[11]
    if (sigma50<0) {sigma[6]<-sigma[11]}
    sigma[5]<-sigma[6]+a1
    aa1<-sigma40;aa2<-sigma50
    if (aa1<0) {aa1<-0};if (aa2<0) {aa2<-0}
    n_iter<-0;aaa0<-sigma[11]
    s0[1]<-ss1+ss2+ss3;s0[2]<-nn[1]+nn[2]+nn[3]
    a3<-1000
    while (a3>0.0001)
    {
      n_iter<-n_iter+1
      aa3<-sigma[11]/(sigma[11]+aa1)
      aa4<-sigma[11]/(sigma[11]+aa1+a1)
      aa5<-sigma[11]/(sigma[11]+aa2+a1)
      aa6<-sigma[11]/(sigma[11]+aa2)
      sigma[11]<-(s0[1]+(aa3^2*swx[1]+aa4^2*swx[2]+aa5^2*swx[5]+aa6^2*swx[6]))/(s0[2]+aa3*n0[1]+aa4*n0[2]+aa5*n0[5]+aa6*n0[6])
      a3<-abs(sigma[11]-aaa0)
      aaa0<-sigma[11]
      if (n_iter>20) break
    } 
    sigma[1]<-sigma[11]+sigma40;sigma[2]<-sigma[11]+sigma40+a1;sigma[5]<-sigma[11]+sigma50+a1;sigma[6]<-sigma[11]+sigma50
    mm1<-m[c(1,2)];ssigma1<-sigma[c(1,2)]
    mm2<-m[c(5,6)];ssigma2<-sigma[c(5,6)]
    mix_pi1<-mix_pi[c(1,2)];mix_pi2<-mix_pi[c(5,6)]
    ########criteria for iterations to stop########
    L1<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
      sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mix_pi2)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
    if(iteration>300)break
  }
  abc<-L0
  AIC<--2*abc+2*8
  ########second order genetic parameters###############   
  B1jj <- sigma1 - sigma[1]
  B1_gg<- sigma[1]-sigma[11]
  if(B1jj<0) {B1jj<-0}
  if(B1_gg<0 || B1_gg>sigma1){B1_gg<-0}
  B1ll <- B1jj/sigma1
  B1rr <- B1_gg/sigma1
  B2jj <- sigma2 - sigma[6]
  B2_gg <- sigma[6]-sigma[11]
  if(B2jj<0){B2jj<-0}
  if(B2_gg<0 || B2_gg>sigma2){B2_gg<-0}
  B2ll <- B2jj/sigma2
  B2rr <- B2_gg/sigma2
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*nn[1])
  P1bmw <- matrix(0,nn[1],1)
  P1gg <- (dataP1 - m[11])/sqrt(as.vector(sigma[11]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn1 < nn[1]){P1bmw <- P1bmw+runif(nn[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:nn[1])) - 0.5)/nn[1])^2)
  P1u<- as.matrix(c(12*nn[1]*((P1dd[1]/nn[1]-0.5)^2),((45*nn[1])/4)*((P1dd[2]/nn[1]-1/3)^2),180*nn[1]*((P1dd[3]/nn[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,nn[1]))))
  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*nn[2])
  F1bmw <- matrix(0,nn[2],1)
  F1gg <- (dataF1 - m[12])/sqrt(as.vector(sigma[11]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn1 < nn[2]){F1bmw <- F1bmw+runif(nn[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:nn[2])) - 0.5)/nn[2])^2)
  F1u<- as.matrix(c(12*nn[2]*((F1dd[1]/nn[2]-0.5)^2),((45*nn[2])/4)*((F1dd[2]/nn[2]-1/3)^2),180*nn[2]*((F1dd[3]/nn[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,nn[2]))))
  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)
  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*nn[3])
  P2bmw <- matrix(0,nn[3],1)
  P2gg <- (dataP2 - m[13])/sqrt(as.vector(sigma[11]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn1 < nn[3]){P2bmw <- P2bmw+runif(nn[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:nn[3])) - 0.5)/nn[3])^2)
  P2u<- as.matrix(c(12*nn[3]*((P2dd[1]/nn[3]-0.5)^2),((45*nn[3])/4)*((P2dd[2]/nn[3]-1/3)^2),180*nn[3]*((P2dd[3]/nn[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,nn[3]))))
  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1); 
  B1w1<-1/(12*nn[4])
  B1bmw <- matrix(0,nn[4],1); B1bmwsl <- matrix(0,nn[4],d2)
  for(i in 1:d2){
    B1gg <- (dataB1 - m[i])/sqrt(sigma[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B1P2)))[1]
  if(nn1 < nn[4]){B1P2 <- B1P2+runif(nn[4])/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*nn[4]) + sum((B1P2 - (as.matrix(c(1:nn[4])) - 0.5)/nn[4])^2)
  B1u <- as.matrix(c(12*nn[4]*((B1dd[1]/nn[4]-0.5)^2),((45*nn[4])/4)*((B1dd[2]/nn[4]-1/3)^2),180*nn[4]*((B1dd[3]/nn[4]-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),(1-pkolm(B1D,nn[4]))))
  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2); 
  B2w1<-1/(12*nn[5])
  B2bmw <- matrix(0,nn[5],1); B2bmwsl <- matrix(0,nn[5],d2)
  for(i in 1:d2){
    B2gg <- (dataB2 - mm2[i])/sqrt(ssigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B2P2)))[1]
  if(nn1 < nn[5]){B2P2 <- B2P2+runif(nn[5])/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*nn[5]) + sum((B2P2 - (as.matrix(c(1:nn[5])) - 0.5)/nn[5])^2)
  B2u <- as.matrix(c(12*nn[5]*((B2dd[1]/nn[5]-0.5)^2),((45*nn[5])/4)*((B2dd[2]/nn[5]-1/3)^2),180*nn[5]*((B2dd[3]/nn[5]-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),(1-pkolm(B2D,nn[5]))))
  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("MX1-AD-AD",round(abc,4),round(AIC,4),round(m[11],4),round(m[12],4),round(m[13],4), round(m[1],4),round(m[2],4)," "," ",round(sigma[1],4),round(sigma[2],4),
                       " "," ",round(mix_pi[1],4),round(mix_pi[2],4)," "," ",round(m[5],4),round(m[6],4)," "," ",round(sigma[5],4),round(sigma[6],4)," "," ",        
                       round(mix_pi[5],4),round(mix_pi[6],4)," "," ",round(sigma[11],4),round(B[1],4),round(B[2],4)," ",round(B[3],4)," ",round(B[4],4),round(B[5],4),round(B1jj,4),round(B1ll*100,4),round(B1_gg,4),round(B1rr*100,4),round(B2jj,4),round(B2ll*100,4),round(B2_gg,4),round(B2rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,as.matrix(mmi1),as.matrix(mmi2))
  return(OUTPUT)
}

##################MX1-A-AD(D-2)############################
G5BCFModelFun[[11]] <- function(K1,logL,df11,df21,df31,df41,df51,G5BCFtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]))
  nn<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataB1)[1],dim(dataB2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2)))
  ss<-s[c(1:3)]-sumx[c(1:3)]^2/nn[c(1:3)]
  m<-matrix(0,13,1);m[11]<-mean(dataP1);m[12]<-mean(dataF1);m[13]<-mean(dataP2);man0<-mean(dataB1);man1<-mean(dataB2)
  sigma0<-sum(ss)/(nn[1]+nn[2]+nn[3]-3)
  sigma1<-var(dataB1);sigma2<-var(dataB2)
  d1<-1;d2<-2;d3<-3;d4<-4
  m_esp <- 0.0001;num_l <- as.numeric(G5BCFtext2)
  ###############procedure start###########################  
  mi<-matrix(0,6,1);mi[c(1,2,5,6)]<-0.5
  sigma<-matrix(0,11,1);sigma[1]<-sigma1/3;sigma[2]<-sigma[1];sigma[5]<-sigma2/3;sigma[6]<-sigma[5];sigma[11]<-sigma0
  a1<-sqrt(sigma1/(nn[4]-1))
  m[1]<-man0+a1;m[2]<-man0-a1
  a1<-sqrt(sigma2/(nn[5]-1))
  m[5]<-man1+a1;m[6]<-man1-a1
  mm1<-m[c(1,2)];mmi1<-mi[c(1,2)];ssigma1<-sigma[c(1,2)]
  mm2<-m[c(5,6)];mmi2<-mi[c(5,6)];ssigma2<-sigma[c(5,6)]
  L0<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
    sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)))+
    sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,2,nn[4]); swx1 <- matrix(0,2,1)
  W2 <- matrix(0,2,nn[5]); swx2 <- matrix(0,2,1)
  sumwx<-matrix(0,6,1);mix_pi<-matrix(0,6,1)
  n0<-matrix(0,6,1); swx<-matrix(0,6,1)
  hh<-matrix(0,3,3);b_line<-matrix(0,3,1)
  s0<-matrix(0,2,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:2) { W1[i,] <- mi[i]*dnorm(dataB1,m[i],sqrt(sigma[i]))/dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)}
    mix_pi[c(1,2)] <- as.matrix(rowSums(W1)/nn[4])
    sumwx[c(1,2)] <- W1%*%dataB1
    for(i in 1:2) { W2[i,] <- mmi2[i]*dnorm(dataB2,mm2[i],sqrt(ssigma2[i]))/dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)}
    mix_pi[c(5,6)] <- as.matrix(rowSums(W2)/nn[5])
    sumwx[c(5,6)] <- W2%*%dataB2
    n0[c(1,2)]<-nn[4]*mix_pi[c(1,2)];n0[c(5,6)]<-nn[5]*mix_pi[c(5,6)]
    n0[c(1,2,5,6)][abs(n0[c(1,2,5,6)])<0.000001]<-0.000001
    ##############################################
    hh[1,1]<-25*sigma[11]/nn[1]+4*sigma[11]/nn[2]+sigma[11]/nn[3]+16*sigma[1]/n0[1]+16*sigma[2]/n0[2]
    hh[1,2]<-5*sigma[11]/nn[1]+4*sigma[11]/nn[2]+5*sigma[11]/nn[3]
    hh[1,3]<--4*sigma[1]/n0[1]+4*sigma[2]/n0[2]
    hh[2,2]<-sigma[11]/nn[1]+4*sigma[11]/nn[2]+25*sigma[11]/nn[3]+16*sigma[5]/n0[5]+16*sigma[6]/n0[6]
    hh[2,3]<-4*sigma[5]/n0[5]-4*sigma[6]/n0[6]
    hh[3,3]<-sigma[1]/n0[1]+sigma[2]/n0[2]+sigma[5]/n0[5]+sigma[6]/n0[6]
    for(i in 2:3)
    {
      for(j in 1:(i-1))
      {
        hh[i,j]<-hh[j,i]
      }
    }
    #########################################
    b_line[1]<-5*sumx[1]/nn[1]+2*sumx[2]/nn[2]+sumx[3]/nn[3]-4*sumwx[1]/n0[1]-4*sumwx[2]/n0[2]
    b_line[2]<-sumx[1]/nn[1]+2*sumx[2]/nn[2]+5*sumx[3]/nn[3]-4*sumwx[5]/n0[5]-4*sumwx[6]/n0[6]
    b_line[3]<-sumwx[1]/n0[1]-sumwx[2]/n0[2]-sumwx[5]/n0[5]+sumwx[6]/n0[6]
    B1<-solve(hh,b_line)
    #########################################
    m[11]<-(sumx[1]-sigma[11]*(5*B1[1]+B1[2]))/nn[1]
    m[12]<-(sumx[2]-sigma[11]*(2*B1[1]+B1[2]*2))/nn[2]
    m[13]<-(sumx[3]-sigma[11]*(B1[1]+5*B1[2]))/nn[3]
    m[1]<-(sumwx[1]+sigma[1]*(4*B1[1]-B1[3]))/n0[1]
    m[2]<-(sumwx[2]+sigma[2]*(4*B1[1]+B1[3]))/n0[2]
    m[5]<-(sumwx[5]+sigma[5]*(4*B1[2]+B1[3]))/n0[5]
    m[6]<-(sumwx[6]+sigma[6]*(4*B1[2]-B1[3]))/n0[6]
    hh2<-matrix(c(1,1,1,1,1,1,1,1,0,-1,1,0,0,-1,1,0,-1,0.5,0.5,-0.5,-0.5,0,1,0,0.5,0.5,0.25,0.25),7,4)
    mm<-m[c(11,12,13,1,2,5,6)]
    B2<-solve(crossprod(hh2,hh2))%*%crossprod(hh2,mm)
    a1<-B2[2];a1<-(0.5*a1^2)/num_l
    ########obtain variance#######################
    ss1<-sum((dataP1-m[11])^2);ss3<-sum((dataP2-m[13])^2);ss2<-sum((dataF1-m[12])^2)
    mm2<-m[c(5,6)]
    for(i in 1:d2) {swx1[i] <- W1[i,]%*%(dataB1-m[i])^2 } ;for(i in 1:d2) {swx2[i] <- W2[i,]%*%(dataB2-mm2[i])^2 }  
    swx[c(1,2)]<-swx1;swx[c(5,6)]<-swx2
    n_iter<-0;aaa0<-sigma[1];aa2<-1000
    while (aa2>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma[1]/(sigma[1]+a1)
      sigma[1]<-(swx[1]+aa1^2*swx[2])/(n0[1]+aa1*n0[2])
      aa2<-abs(sigma[1]-aaa0)
      aaa0<-sigma[1]
      if (n_iter>20) break
    }
    sigma40<-sigma[1]-sigma[11]
    if (sigma40<0) {sigma40<-0;sigma[1]<-sigma[11]}
    sigma[2]<-sigma[1]+a1
    n_iter<-0;aaa0<-sigma[6];aa2<-1000
    while(aa2>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma[6]/(sigma[6]+a1)
      sigma[6]<-(swx[6]+aa1^2*swx[5])/(n0[6]+aa1*n0[5])
      aa2<-abs(sigma[6]-aaa0)
      aaa0<-sigma[6]
      if(n_iter>20) break
    }
    sigma50<-sigma[6]-sigma[11]
    if (sigma50<0) {sigma50<-0;sigma[6]<-sigma[11]}
    sigma[5]<-sigma[6]+a1
    aa1<-sigma40;aa2<-sigma50
    if (aa1<0) {aa1<-0};if (aa2<0) {aa2<-0}
    n_iter<-0;aaa0<-sigma[11]
    s0[1]<-ss1+ss2+ss3;s0[2]<-nn[1]+nn[2]+nn[3];a3<-1000
    while (a3>0.0001)
    {
      n_iter<-n_iter+1
      aa3<-sigma[11]/(sigma[11]+aa1)
      aa4<-sigma[11]/(sigma[11]+aa1+a1)
      aa5<-sigma[11]/(sigma[11]+aa2+a1)
      aa6<-sigma[11]/(sigma[11]+aa2)
      sigma[11]<-(s0[1]+(aa3^2*swx[1]+aa4^2*swx[2]+aa5^2*swx[5]+aa6^2*swx[6]))/(s0[2]+aa3*n0[1]+aa4*n0[2]+aa5*n0[5]+aa6*n0[6])
      a3<-abs(sigma[11]-aaa0)
      aaa0<-sigma[11]
      if (n_iter>20) break
    } 
    sigma[1]<-sigma[11]+sigma40;sigma[2]<-sigma[11]+sigma40+a1;sigma[5]<-sigma[11]+sigma50+a1;sigma[6]<-sigma[11]+sigma50
    mm1<-m[c(1,2)];ssigma1<-sigma[c(1,2)]
    mm2<-m[c(5,6)];ssigma2<-sigma[c(5,6)]
    mix_pi1<-mix_pi[c(1,2)];mix_pi2<-mix_pi[c(5,6)]
    ########criteria for iterations to stop########
    L1<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
      sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mix_pi2)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
    if(iteration>300)break
  }
  abc<-L0
  AIC<--2*abc+2*7
  ########second order genetic parameters###############   
  B1jj <- sigma1 - sigma[1]
  B1_gg <- sigma[1]-sigma[11]
  if(B1jj<0) {B1jj<-0}
  if(B1_gg<0 || B1_gg>sigma1){B1_gg<-0}
  B1ll <- B1jj/sigma1
  B1rr <- B1_gg/sigma1
  B2jj <- sigma2 - sigma[6]
  B2_gg <- sigma[6]-sigma[11]
  if(B2jj<0){B2jj<-0}
  if(B2_gg<0 || B2_gg>sigma2){B2_gg<-0}
  B2ll <- B2jj/sigma2
  B2rr <- B2_gg/sigma2
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*nn[1])
  P1bmw <- matrix(0,nn[1],1)
  P1gg <- (dataP1 - m[11])/sqrt(as.vector(sigma[11]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn1 < nn[1]){P1bmw <- P1bmw+runif(nn[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:nn[1])) - 0.5)/nn[1])^2)
  P1u<- as.matrix(c(12*nn[1]*((P1dd[1]/nn[1]-0.5)^2),((45*nn[1])/4)*((P1dd[2]/nn[1]-1/3)^2),180*nn[1]*((P1dd[3]/nn[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,nn[1]))))
  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*nn[2])
  F1bmw <- matrix(0,nn[2],1)
  F1gg <- (dataF1 - m[12])/sqrt(as.vector(sigma[11]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn1 < nn[2]){F1bmw <- F1bmw+runif(nn[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:nn[2])) - 0.5)/nn[2])^2)
  F1u<- as.matrix(c(12*nn[2]*((F1dd[1]/nn[2]-0.5)^2),((45*nn[2])/4)*((F1dd[2]/nn[2]-1/3)^2),180*nn[2]*((F1dd[3]/nn[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,nn[2]))))
  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)
  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*nn[3])
  P2bmw <- matrix(0,nn[3],1)
  P2gg <- (dataP2 - m[13])/sqrt(as.vector(sigma[11]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn1 < nn[3]){P2bmw <- P2bmw+runif(nn[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:nn[3])) - 0.5)/nn[3])^2)
  P2u<- as.matrix(c(12*nn[3]*((P2dd[1]/nn[3]-0.5)^2),((45*nn[3])/4)*((P2dd[2]/nn[3]-1/3)^2),180*nn[3]*((P2dd[3]/nn[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,nn[3]))))
  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1); 
  B1w1<-1/(12*nn[4])
  B1bmw <- matrix(0,nn[4],1); B1bmwsl <- matrix(0,nn[4],d2)
  for(i in 1:d2){
    B1gg <- (dataB1 - m[i])/sqrt(sigma[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B1P2)))[1]
  if(nn1 < nn[4]){B1P2 <- B1P2+runif(nn[4])/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*nn[4]) + sum((B1P2 - (as.matrix(c(1:nn[4])) - 0.5)/nn[4])^2)
  B1u <- as.matrix(c(12*nn[4]*((B1dd[1]/nn[4]-0.5)^2),((45*nn[4])/4)*((B1dd[2]/nn[4]-1/3)^2),180*nn[4]*((B1dd[3]/nn[4]-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),(1-pkolm(B1D,nn[4]))))
  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2); 
  B2w1<-1/(12*nn[5])
  B2bmw <- matrix(0,nn[5],1); B2bmwsl <- matrix(0,nn[5],d2)
  for(i in 1:d2){
    B2gg <- (dataB2 - mm2[i])/sqrt(ssigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B2P2)))[1]
  if(nn1 < nn[5]){B2P2 <- B2P2+runif(nn[5])/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*nn[5]) + sum((B2P2 - (as.matrix(c(1:nn[5])) - 0.5)/nn[5])^2)
  B2u <- as.matrix(c(12*nn[5]*((B2dd[1]/nn[5]-0.5)^2),((45*nn[5])/4)*((B2dd[2]/nn[5]-1/3)^2),180*nn[5]*((B2dd[3]/nn[5]-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),(1-pkolm(B2D,nn[5]))))
  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("MX1-A-AD",round(abc,4),round(AIC,4),round(m[11],4),round(m[12],4),round(m[13],4), round(m[1],4),round(m[2],4)," "," ",round(sigma[1],4),round(sigma[2],4),
                       " "," ",round(mix_pi[1],4),round(mix_pi[2],4)," "," ",round(m[5],4),round(m[6],4)," "," ",round(sigma[5],4),round(sigma[6],4)," "," ",        
                       round(mix_pi[5],4),round(mix_pi[6],4)," "," ",round(sigma[11],4),round(B2[1],4),round(B2[2],4)," "," "," ",round(B2[3],4),round(B2[4],4),round(B1jj,4),round(B1ll*100,4),round(B1_gg,4),round(B1rr*100,4),round(B2jj,4),round(B2ll*100,4),round(B2_gg,4),round(B2rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,as.matrix(mmi1),as.matrix(mmi2))
  return(OUTPUT)
}

#################MX1-EAD-AD(D-3)#################################
G5BCFModelFun[[12]] <- function(K1,logL,df11,df21,df31,df41,df51,G5BCFtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]))
  nn<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataB1)[1],dim(dataB2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2)))
  ss<-s[c(1:3)]-sumx[c(1:3)]^2/nn[c(1:3)]
  m<-matrix(0,13,1);m[11]<-mean(dataP1);m[12]<-mean(dataF1);m[13]<-mean(dataP2);man0<-mean(dataB1);man1<-mean(dataB2)
  sigma0<-sum(ss)/(nn[1]+nn[2]+nn[3]-3)
  sigma1<-var(dataB1);sigma2<-var(dataB2)
  d1<-1;d2<-2;d3<-3;d4<-4
  m_esp <- 0.0001;num_l <- as.numeric(G5BCFtext2)
  ###############procedure start###########################  
  mi<-matrix(0,6,1);mi[c(1,2,5,6)]<-0.5
  sigma<-matrix(0,11,1);sigma[1]<-sigma1/3;sigma[2]<-sigma[1];sigma[5]<-sigma2/3;sigma[6]<-sigma[5];sigma[11]<-sigma0
  a1<-sqrt(sigma1/(nn[4]-1))
  m[1]<-man0+a1;m[2]<-man0-a1
  a1<-sqrt(sigma2/(nn[5]-1))
  m[5]<-man1+a1;m[6]<-man1-a1
  mm1<-m[c(1,2)];mmi1<-mi[c(1,2)];ssigma1<-sigma[c(1,2)]
  mm2<-m[c(5,6)];mmi2<-mi[c(5,6)];ssigma2<-sigma[c(5,6)]
  L0<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
    sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)))+
    sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,2,nn[4]); swx1 <- matrix(0,2,1)
  W2 <- matrix(0,2,nn[5]); swx2 <- matrix(0,2,1)
  sumwx<-matrix(0,6,1);mix_pi<-matrix(0,6,1)
  n0<-matrix(0,6,1);hh<-matrix(0,3,3);b_line<-matrix(0,3,1)
  swx<-matrix(0,6,1); s0<-matrix(0,2,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:2) { W1[i,] <- mi[i]*dnorm(dataB1,m[i],sqrt(sigma[i]))/dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)}
    mix_pi[c(1,2)] <- as.matrix(rowSums(W1)/nn[4])
    sumwx[c(1,2)] <- W1%*%dataB1
    for(i in 1:2) { W2[i,] <- mmi2[i]*dnorm(dataB2,mm2[i],sqrt(ssigma2[i]))/dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)}
    mix_pi[c(5,6)] <- as.matrix(rowSums(W2)/nn[5])
    sumwx[c(5,6)] <- W2%*%dataB2
    n0[c(1,2)]<-nn[4]*mix_pi[c(1,2)];n0[c(5,6)]<-nn[5]*mix_pi[c(5,6)]
    n0[c(1,2,5,6)][abs(n0[c(1,2,5,6)])<0.000001]<-0.000001
    ######################################
    hh[1,1]<-25*sigma[11]/nn[1]+4*sigma[11]/nn[2]+sigma[11]/nn[3]+16*sigma[1]/n0[1]+16*sigma[2]/n0[2]
    hh[1,2]<-5*sigma[11]/nn[1]+4*sigma[11]/nn[2]+5*sigma[11]/nn[3]
    hh[1,3]<--12*sigma[1]/n0[1]+12*sigma[2]/n0[2]
    hh[2,2]<-sigma[11]/nn[1]+4*sigma[11]/nn[2]+25*sigma[11]/nn[3]+16*sigma[5]/n0[5]+16*sigma[6]/n0[6]
    hh[2,3]<-4*sigma[5]/n0[5]-4*sigma[6]/n0[6]
    hh[3,3]<-9*sigma[1]/n0[1]+9*sigma[2]/n0[2]+sigma[5]/n0[5]+sigma[6]/n0[6]
    for(i in 2:3)
    {
      for(j in 1:(i-1))
      {
        hh[i,j]<-hh[j,i]
      }
    }
    ##########################################
    b_line[1]<-5*sumx[1]/nn[1]+2*sumx[2]/nn[2]+sumx[3]/nn[3]-4*sumwx[1]/n0[1]-4*sumwx[2]/n0[2]
    b_line[2]<-sumx[1]/nn[1]+2*sumx[2]/nn[2]+5*sumx[3]/nn[3]-4*sumwx[5]/n0[5]-4*sumwx[6]/n0[6]
    b_line[3]<-3*sumwx[1]/n0[1]-3*sumwx[2]/n0[2]-sumwx[5]/n0[5]+sumwx[6]/n0[6]
    B1<-solve(hh,b_line)
    #########################################
    m[11]<-(sumx[1]-sigma[11]*(5*B1[1]+B1[2]))/nn[1]
    m[12]<-(sumx[2]-sigma[11]*(2*B1[1]+B1[2]*2))/nn[2]
    m[13]<-(sumx[3]-sigma[11]*(B1[1]+5*B1[2]))/nn[3]
    m[1]<-(sumwx[1]+sigma[1]*(B1[1]*4-3*B1[3]))/n0[1]
    m[2]<-(sumwx[2]+sigma[2]*(B1[1]*4+3*B1[2]))/n0[2]
    m[5]<-(sumwx[5]+sigma[5]*(4*B1[2]+B1[3]))/n0[5]
    m[6]<-(sumwx[6]+sigma[6]*(4*B1[2]-B1[3]))/n0[6]
    hh2<-matrix(c(1,1,1,1,1,1,1,1,1,-1,1,0.5,0.5,-1,1,0,-1,0.5,0.5,-0.5,-0.5,0,1,0,0.5,0.5,0.25,0.25),7,4)
    mm<-m[c(11,12,13,1,2,5,6)]
    B2<-solve(crossprod(hh2,hh2))%*%crossprod(hh2,mm)
    a1<-B2[2];a1<-(0.75*a1^2)/num_l
    ########obtain variance#######################
    ss1<-sum((dataP1-m[11])^2);ss3<-sum((dataP2-m[13])^2);ss2<-sum((dataF1-m[12])^2)
    mm2<-m[c(5,6)]
    for(i in 1:d2) {swx1[i] <- W1[i,]%*%(dataB1-m[i])^2 } ;for(i in 1:d2) {swx2[i] <- W2[i,]%*%(dataB2-mm2[i])^2 }  
    swx[c(1,2)]<-swx1;swx[c(5,6)]<-swx2
    n_iter<-0;aaa0<-sigma[1];aa2<-1000
    while (aa2>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma[1]/(sigma[1]+a1)
      sigma[1]<-(swx[1]+aa1^2*swx[2])/(n0[1]+aa1*n0[2])
      aa2<-abs(sigma[1]-aaa0)
      aaa0<-sigma[1]
      if (n_iter>20) break
    }
    sigma40<-sigma[1]-sigma[11]
    if (sigma40<0) {sigma40<-0;sigma[1]<-sigma[11]}
    sigma[2]<-sigma[1]+a1
    n_iter<-0;aaa0<-sigma[6];aa2<-1000
    while(aa2>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma[6]/(sigma[6]+a1)
      sigma[6]<-(swx[6]+aa1^2*swx[5])/(n0[6]+aa1*n0[5])
      aa2<-abs(sigma[6]-aaa0)
      aaa0<-sigma[6]
      if(n_iter>20) break
    }
    sigma50<-sigma[6]-sigma[11]
    if (sigma50<0) {sigma50<-0;sigma[6]<-sigma[11]}
    sigma[5]<-sigma[6]+a1
    aa1<-sigma40;aa2<-sigma50
    if (aa1<0) {aa1<-0};if (aa2<0) {aa2<-0}
    n_iter<-0;aaa0<-sigma[11]
    s0[1]<-ss1+ss2+ss3;s0[2]<-nn[1]+nn[2]+nn[3];a3<-1000
    while (a3>0.0001)
    {
      n_iter<-n_iter+1
      aa3<-sigma[11]/(sigma[11]+aa1)
      aa4<-sigma[11]/(sigma[11]+aa1+a1)
      aa5<-sigma[11]/(sigma[11]+aa2+a1)
      aa6<-sigma[11]/(sigma[11]+aa2)
      sigma[11]<-(s0[1]+(aa3^2*swx[1]+aa4^2*swx[2]+aa5^2*swx[5]+aa6^2*swx[6]))/(s0[2]+aa3*n0[1]+aa4*n0[2]+aa5*n0[5]+aa6*n0[6])
      a3<-abs(sigma[11]-aaa0)
      aaa0<-sigma[11]
      if (n_iter>20) break
    } 
    sigma[1]<-sigma[11]+sigma40;sigma[2]<-sigma[11]+sigma40+a1
    sigma[5]<-sigma[11]+sigma50+a1;sigma[6]<-sigma[11]+sigma50
    mm1<-m[c(1,2)];ssigma1<-sigma[c(1,2)]
    mm2<-m[c(5,6)];ssigma2<-sigma[c(5,6)]
    mix_pi1<-mix_pi[c(1,2)];mix_pi2<-mix_pi[c(5,6)]
    ########criteria for iterations to stop########
    L1<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
      sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mix_pi2)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
    if(iteration>300)break
  }
  abc<-L0
  AIC<--2*abc+2*7
  ########second order genetic parameters###############   
  B1jj <- sigma1 - sigma[1]
  B1_gg <- sigma[1]-sigma[11]
  if(B1jj<0) {B1jj<-0}
  if(B1_gg<0 || B1_gg>sigma1){B1_gg<-0}
  B1ll <- B1jj/sigma1
  B1rr <- B1_gg/sigma1
  B2jj <- sigma2 - sigma[6]
  B2_gg <- sigma[6]-sigma[11]
  if(B2jj<0){B2jj<-0}
  if(B2_gg<0 || B2_gg>sigma2){B2_gg<-0}
  B2ll <- B2jj/sigma2
  B2rr <- B2_gg/sigma2
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*nn[1])
  P1bmw <- matrix(0,nn[1],1)
  P1gg <- (dataP1 - m[11])/sqrt(as.vector(sigma[11]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn1 < nn[1]){P1bmw <- P1bmw+runif(nn[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:nn[1])) - 0.5)/nn[1])^2)
  P1u<- as.matrix(c(12*nn[1]*((P1dd[1]/nn[1]-0.5)^2),((45*nn[1])/4)*((P1dd[2]/nn[1]-1/3)^2),180*nn[1]*((P1dd[3]/nn[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,nn[1]))))
  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*nn[2])
  F1bmw <- matrix(0,nn[2],1)
  F1gg <- (dataF1 - m[12])/sqrt(as.vector(sigma[11]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn1 < nn[2]){F1bmw <- F1bmw+runif(nn[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:nn[2])) - 0.5)/nn[2])^2)
  F1u<- as.matrix(c(12*nn[2]*((F1dd[1]/nn[2]-0.5)^2),((45*nn[2])/4)*((F1dd[2]/nn[2]-1/3)^2),180*nn[2]*((F1dd[3]/nn[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,nn[2]))))
  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)
  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*nn[3])
  P2bmw <- matrix(0,nn[3],1)
  P2gg <- (dataP2 - m[13])/sqrt(as.vector(sigma[11]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn1 < nn[3]){P2bmw <- P2bmw+runif(nn[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:nn[3])) - 0.5)/nn[3])^2)
  P2u<- as.matrix(c(12*nn[3]*((P2dd[1]/nn[3]-0.5)^2),((45*nn[3])/4)*((P2dd[2]/nn[3]-1/3)^2),180*nn[3]*((P2dd[3]/nn[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,nn[3]))))
  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1); 
  B1w1<-1/(12*nn[4])
  B1bmw <- matrix(0,nn[4],1); B1bmwsl <- matrix(0,nn[4],d2)
  for(i in 1:d2){
    B1gg <- (dataB1 - m[i])/sqrt(sigma[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B1P2)))[1]
  if(nn1 < nn[4]){B1P2 <- B1P2+runif(nn[4])/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*nn[4]) + sum((B1P2 - (as.matrix(c(1:nn[4])) - 0.5)/nn[4])^2)
  B1u <- as.matrix(c(12*nn[4]*((B1dd[1]/nn[4]-0.5)^2),((45*nn[4])/4)*((B1dd[2]/nn[4]-1/3)^2),180*nn[4]*((B1dd[3]/nn[4]-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),(1-pkolm(B1D,nn[4]))))
  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2); 
  B2w1<-1/(12*nn[5])
  B2bmw <- matrix(0,nn[5],1); B2bmwsl <- matrix(0,nn[5],d2)
  for(i in 1:d2){
    B2gg <- (dataB2 - mm2[i])/sqrt(ssigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B2P2)))[1]
  if(nn1 < nn[5]){B2P2 <- B2P2+runif(nn[5])/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*nn[5]) + sum((B2P2 - (as.matrix(c(1:nn[5])) - 0.5)/nn[5])^2)
  B2u <- as.matrix(c(12*nn[5]*((B2dd[1]/nn[5]-0.5)^2),((45*nn[5])/4)*((B2dd[2]/nn[5]-1/3)^2),180*nn[5]*((B2dd[3]/nn[5]-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),(1-pkolm(B2D,nn[5]))))
  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("MX1-EAD-AD",round(abc,4),round(AIC,4),round(m[11],4),round(m[12],4),round(m[13],4), round(m[1],4),round(m[2],4)," "," ",round(sigma[1],4),round(sigma[2],4),
                       " "," ",round(mix_pi[1],4),round(mix_pi[2],4)," "," ",round(m[5],4),round(m[6],4)," "," ",round(sigma[5],4),round(sigma[6],4)," "," ",        
                       round(mix_pi[5],4),round(mix_pi[6],4)," "," ",round(sigma[11],4),round(B2[1],4),round(B2[2],4)," "," "," ",round(B2[3],4),round(B2[4],4),round(B1jj,4),round(B1ll*100,4),round(B1_gg,4),round(B1rr*100,4),round(B2jj,4),round(B2ll*100,4),round(B2_gg,4),round(B2rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,as.matrix(mmi1),as.matrix(mmi2))
  return(OUTPUT)
}

################MX1-NCD-AD(D-4)########################
G5BCFModelFun[[13]] <- function(K1,logL,df11,df21,df31,df41,df51,G5BCFtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]))
  nn<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataB1)[1],dim(dataB2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2)))
  ss<-s[c(1:3)]-sumx[c(1:3)]^2/nn[c(1:3)]
  m<-matrix(0,13,1);m[11]<-mean(dataP1);m[12]<-mean(dataF1);m[13]<-mean(dataP2);man0<-mean(dataB1);man1<-mean(dataB2)
  sigma0<-sum(ss)/(nn[1]+nn[2]+nn[3]-3)
  sigma1<-var(dataB1);sigma2<-var(dataB2)
  d1<-1;d2<-2;d3<-3;d4<-4
  m_esp <- 0.0001;num_l <- as.numeric(G5BCFtext2)
  ###############procedure start###########################  
  mi<-matrix(0,6,1);mi[c(1,2,5,6)]<-0.5
  sigma<-matrix(0,11,1)
  sigma[1]<-sigma1/3;sigma[2]<-sigma[1];sigma[5]<-sigma2/3;sigma[6]<-sigma[5];sigma[11]<-sigma0
  a1<-sqrt(sigma1/(nn[4]-1))
  m[1]<-man0+a1;m[2]<-man0-a1
  a1<-sqrt(sigma2/(nn[5]-1))
  m[5]<-man1+a1;m[6]<-man1-a1
  mm1<-m[c(1,2)];mmi1<-mi[c(1,2)];ssigma1<-sigma[c(1,2)]
  mm2<-m[c(5,6)];mmi2<-mi[c(5,6)];ssigma2<-sigma[c(5,6)]
  L0<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
    sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)))+
    sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,2,nn[4]); swx1 <- matrix(0,2,1)
  W2 <- matrix(0,2,nn[5]); swx2 <- matrix(0,2,1)
  sumwx<-matrix(0,6,1);mix_pi<-matrix(0,6,1)
  n0<-matrix(0,6,1);swx<-matrix(0,6,1);s0<-matrix(0,2,1)
  hh<-matrix(0,3,3);b_line<-matrix(0,3,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:2) { W1[i,] <- mi[i]*dnorm(dataB1,m[i],sqrt(sigma[i]))/dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)}
    mix_pi[c(1,2)] <- as.matrix(rowSums(W1)/nn[4])
    sumwx[c(1,2)] <- W1%*%dataB1
    for(i in 1:2) { W2[i,] <- mmi2[i]*dnorm(dataB2,mm2[i],sqrt(ssigma2[i]))/dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)}
    mix_pi[c(5,6)] <- as.matrix(rowSums(W2)/nn[5])
    sumwx[c(5,6)] <- W2%*%dataB2
    n0[c(1,2)]<-nn[4]*mix_pi[c(1,2)]; n0[c(5,6)]<-nn[5]*mix_pi[c(5,6)]
    n0[c(1,2,5,6)][abs(n0[c(1,2,5,6)])<0.000001]<-0.000001
    ######################################
    hh[1,1]<-25*sigma[11]/nn[1]+4*sigma[11]/nn[2]+sigma[11]/nn[3]+16*sigma[1]/n0[1]+16*sigma[2]/n0[2]
    hh[1,2]<-5*sigma[11]/nn[1]+4*sigma[11]/nn[2]+5*sigma[11]/nn[3]
    hh[1,3]<--4*sigma[1]/n0[1]+4*sigma[2]/n0[2]
    hh[2,1]<-hh[1,2]
    hh[2,2]<-sigma[11]/nn[1]+4*sigma[11]/nn[2]+25*sigma[11]/nn[3]+16*sigma[5]/n0[5]+16*sigma[6]/n0[6]
    hh[2,3]<-12*sigma[5]/n0[5]-12*sigma[6]/n0[6]
    hh[3,1]<-hh[1,3]
    hh[3,2]<-hh[2,3]
    hh[3,3]<-sigma[1]/n0[1]+sigma[2]/n0[2]+9*sigma[5]/n0[5]+9*sigma[6]/n0[6]
    for(i in 2:3)
    {
      for(j in 1:(i-1))
      {
        hh[i,j]<-hh[j,i]
      }
    }
    #######################################
    b_line[1]<-5*sumx[1]/nn[1]+2*sumx[2]/nn[2]+sumx[3]/nn[3]-4*sumwx[1]/n0[1]-4*sumwx[2]/n0[2]
    b_line[2]<-sumx[1]/nn[1]+2*sumx[2]/nn[2]+5*sumx[3]/nn[3]-4*sumwx[5]/n0[5]-4*sumwx[6]/n0[6]
    b_line[3]<-sumwx[1]/n0[1]-sumwx[2]/n0[2]-3*sumwx[5]/n0[5]+3*sumwx[6]/n0[6]
    B1<-solve(hh,b_line)
    #########################################
    m[11]<-(sumx[1]-sigma[11]*(5*B1[1]+B1[2]))/nn[1]
    m[12]<-(sumx[2]-sigma[11]*(2*B1[1]+B1[2]*2))/nn[2]
    m[13]<-(sumx[3]-sigma[11]*(B1[1]+5*B1[2]))/nn[3]
    m[1]<-(sumwx[1]+sigma[1]*(4*B1[1]-B1[3]))/n0[1]
    m[2]<-(sumwx[2]+sigma[2]*(4*B1[1]+B1[3]))/n0[2]
    m[5]<-(sumwx[5]+sigma[5]*(4*B1[2]+3*B1[3]))/n0[5]
    m[6]<-(sumwx[6]+sigma[6]*(4*B1[2]-3*B1[3]))/n0[6]
    hh2<-matrix(c(1,1,1,1,1,1,1,1,-1,-1,1,-0.5,-0.5,-1,1,0,-1,0.5,0.5,-0.5,-0.5,
                  0,1,0,0.5,0.5,0.25,0.25),7,4)
    mm<-m[c(11,12,13,1,2,5,6)]
    B2<-solve(crossprod(hh2,hh2))%*%crossprod(hh2,mm)
    a1<-B2[2];a1<-(0.75*a1^2)/num_l
    ########obtain variance#######################
    ss1<-sum((dataP1-m[11])^2);ss3<-sum((dataP2-m[13])^2);ss2<-sum((dataF1-m[12])^2)
    mm2<-m[c(5,6)]
    for(i in 1:d2) {swx1[i] <- W1[i,]%*%(dataB1-m[i])^2 } ;for(i in 1:d2) {swx2[i] <- W2[i,]%*%(dataB2-mm2[i])^2 }  
    swx[c(1,2)]<-swx1;swx[c(5,6)]<-swx2
    n_iter<-0;aaa0<-sigma[1];aa2<-1000
    while (aa2>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma[1]/(sigma[1]+a1)
      sigma[1]<-(swx[1]+aa1^2*swx[2])/(n0[1]+aa1*n0[2])
      aa2<-abs(sigma[1]-aaa0)
      aaa0<-sigma[1]
      if (n_iter>20) break
    }
    sigma40<-sigma[1]-sigma[11]
    if (sigma40<0) {sigma40<-0;sigma[1]<-sigma[11]}
    sigma[2]<-sigma[1]+a1
    n_iter<-0;aaa0<-sigma[6];aa2<-1000
    while(aa2>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma[6]/(sigma[6]+a1)
      sigma[6]<-(swx[6]+aa1^2*swx[5])/(n0[6]+aa1*n0[5])
      aa2<-abs(sigma[6]-aaa0)
      aaa0<-sigma[6]
      if(n_iter>20) break
    }
    sigma50<-sigma[6]-sigma[11]
    if (sigma50<0) {sigma50<-0;sigma[6]<-sigma[11]}
    sigma[5]<-sigma[6]+a1
    aa1<-sigma40;aa2<-sigma50
    if (aa1<0) {aa1<-0};if (aa2<0) {aa2<-0}
    n_iter<-0;aaa0<-sigma[11]
    s0[1]<-ss1+ss2+ss3;s0[2]<-nn[1]+nn[2]+nn[3]
    a3<-1000
    while (a3>0.0001)
    {
      n_iter<-n_iter+1
      aa3<-sigma[11]/(sigma[11]+aa1)
      aa4<-sigma[11]/(sigma[11]+aa1+a1)
      aa5<-sigma[11]/(sigma[11]+aa2+a1)
      aa6<-sigma[11]/(sigma[11]+aa2)
      sigma[11]<-(s0[1]+(aa3^2*swx[1]+aa4^2*swx[2]+aa5^2*swx[5]+aa6^2*swx[6]))/(s0[2]+aa3*n0[1]+aa4*n0[2]+aa5*n0[5]+aa6*n0[6])
      a3<-abs(sigma[11]-aaa0)
      aaa0<-sigma[11]
      if (n_iter>20) break
    } 
    sigma[1]<-sigma[11]+sigma40;sigma[2]<-sigma[11]+sigma40+a1
    sigma[5]<-sigma[11]+sigma50+a1;sigma[6]<-sigma[11]+sigma50
    mm1<-m[c(1,2)];ssigma1<-sigma[c(1,2)]
    mm2<-m[c(5,6)];ssigma2<-sigma[c(5,6)]
    mix_pi1<-mix_pi[c(1,2)];mix_pi2<-mix_pi[c(5,6)]
    ########criteria for iterations to stop############
    L1<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
      sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mix_pi2)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
    if(iteration>300)break
  }
  abc<-L0
  AIC<--2*abc+2*7
  ########second order genetic parameters###############   
  B1jj <- sigma1 - sigma[1]
  B1_gg <- sigma[1]-sigma[11]
  if(B1jj<0) {B1jj<-0}
  if(B1_gg<0 || B1_gg>sigma1){B1_gg<-0}
  B1ll <- B1jj/sigma1
  B1rr <- B1_gg/sigma1
  B2jj <- sigma2 - sigma[6]
  B2_gg <- sigma[6]-sigma[11]
  if(B2jj<0){B2jj<-0}
  if(B2_gg<0 || B2_gg>sigma2){B2_gg<-0}
  B2ll <- B2jj/sigma2
  B2rr <- B2_gg/sigma2
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*nn[1])
  P1bmw <- matrix(0,nn[1],1)
  P1gg <- (dataP1 - m[11])/sqrt(as.vector(sigma[11]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn1 < nn[1]){P1bmw <- P1bmw+runif(nn[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:nn[1])) - 0.5)/nn[1])^2)
  P1u<- as.matrix(c(12*nn[1]*((P1dd[1]/nn[1]-0.5)^2),((45*nn[1])/4)*((P1dd[2]/nn[1]-1/3)^2),180*nn[1]*((P1dd[3]/nn[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,nn[1]))))
  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*nn[2])
  F1bmw <- matrix(0,nn[2],1)
  F1gg <- (dataF1 - m[12])/sqrt(as.vector(sigma[11]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn1 < nn[2]){F1bmw <- F1bmw+runif(nn[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:nn[2])) - 0.5)/nn[2])^2)
  F1u<- as.matrix(c(12*nn[2]*((F1dd[1]/nn[2]-0.5)^2),((45*nn[2])/4)*((F1dd[2]/nn[2]-1/3)^2),180*nn[2]*((F1dd[3]/nn[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,nn[2]))))
  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)
  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*nn[3])
  P2bmw <- matrix(0,nn[3],1)
  P2gg <- (dataP2 - m[13])/sqrt(as.vector(sigma[11]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn1 < nn[3]){P2bmw <- P2bmw+runif(nn[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:nn[3])) - 0.5)/nn[3])^2)
  P2u<- as.matrix(c(12*nn[3]*((P2dd[1]/nn[3]-0.5)^2),((45*nn[3])/4)*((P2dd[2]/nn[3]-1/3)^2),180*nn[3]*((P2dd[3]/nn[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,nn[3]))))
  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1); 
  B1w1<-1/(12*nn[4])
  B1bmw <- matrix(0,nn[4],1); B1bmwsl <- matrix(0,nn[4],d2)
  for(i in 1:d2){
    B1gg <- (dataB1 - m[i])/sqrt(sigma[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B1P2)))[1]
  if(nn1 < nn[4]){B1P2 <- B1P2+runif(nn[4])/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*nn[4]) + sum((B1P2 - (as.matrix(c(1:nn[4])) - 0.5)/nn[4])^2)
  B1u <- as.matrix(c(12*nn[4]*((B1dd[1]/nn[4]-0.5)^2),((45*nn[4])/4)*((B1dd[2]/nn[4]-1/3)^2),180*nn[4]*((B1dd[3]/nn[4]-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),(1-pkolm(B1D,nn[4]))))
  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2); 
  B2w1<-1/(12*nn[5])
  B2bmw <- matrix(0,nn[5],1); B2bmwsl <- matrix(0,nn[5],d2)
  for(i in 1:d2){
    B2gg <- (dataB2 - mm2[i])/sqrt(ssigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B2P2)))[1]
  if(nn1 < nn[5]){B2P2 <- B2P2+runif(nn[5])/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*nn[5]) + sum((B2P2 - (as.matrix(c(1:nn[5])) - 0.5)/nn[5])^2)
  B2u <- as.matrix(c(12*nn[5]*((B2dd[1]/nn[5]-0.5)^2),((45*nn[5])/4)*((B2dd[2]/nn[5]-1/3)^2),180*nn[5]*((B2dd[3]/nn[5]-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),(1-pkolm(B2D,nn[5]))))
  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("MX1-NCD-AD",round(abc,4),round(AIC,4),round(m[11],4),round(m[12],4),round(m[13],4), round(m[1],4),round(m[2],4)," "," ",round(sigma[1],4),round(sigma[2],4),
                       " "," ",round(mix_pi[1],4),round(mix_pi[2],4)," "," ",round(m[5],4),round(m[6],4)," "," ",round(sigma[5],4),round(sigma[6],4)," "," ",        
                       round(mix_pi[5],4),round(mix_pi[6],4)," "," ",round(sigma[11],4),round(B2[1],4),round(B2[2],4)," "," "," ",round(B2[3],4),round(B2[4],4),round(B1jj,4),round(B1ll*100,4),round(B1_gg,4),round(B1rr*100,4),round(B2jj,4),round(B2ll*100,4),round(B2_gg,4),round(B2rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,as.matrix(mmi1),as.matrix(mmi2))
  return(OUTPUT)
}

###########MX2-AD-AD(E-2)#####################################
G5BCFModelFun[[14]] <- function(K1,logL,df11,df21,df31,df41,df51,G5BCFtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]))
  nn<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataB1)[1],dim(dataB2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2)))
  ss<-s[c(1:3)]-sumx[c(1:3)]^2/nn[c(1:3)]
  m<-matrix(0,13,1);m[11]<-mean(dataP1);m[12]<-mean(dataF1);m[13]<-mean(dataP2);man0<-mean(dataB1);man1<-mean(dataB2)
  sigma0<-sum(ss)/(nn[1]+nn[2]+nn[3]-3)
  sigma1<-var(dataB1);sigma2<-var(dataB2)
  d1<-1;d2<-2;d3<-3;d4<-4
  m_esp <- 0.0001;num_l <- as.numeric(G5BCFtext2)
  ###############procedure start###########################  
  mi<-matrix(0.25,8,1)
  a1<-sqrt(sigma1/(nn[4]-1))
  m[1]<-man0+2*a1;m[2]<-man0+0.8*a1;m[3]<-man0-0.8*a1;m[4]<-man0-2*a1
  a1<-sqrt(sigma2/(nn[5]-1))
  m[5]<-man1+2*a1;m[6]<-man1+0.8*a1;m[7]<-man1-0.8*a1;m[8]<-man1-2*a1
  sigma<-matrix(0,11,1);sigma[11]<-sigma0;sigma[1]<-sigma1/3;sigma[8]<-sigma2/3
  hh1<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,0,-1,1,1,0,0,0,0,-1,-1,
                1,0,-1,1,0,1,0,0,-1,0,-1,0,1,0,0,0,0.5,0.5,0.5,0.5,0,0,
                0,1,0,0,0.5,0,0.5,0.5,0,0.5,0,1,0,-1,0.5,0.5,0.5,0.5,-0.5,-0.5,-0.5,-0.5,
                0,1,0,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25),11,7)
  mm<-m[c(11,12,13,1,2,3,4,5,6,7,8)]
  B1<-solve(crossprod(hh1,hh1))%*%crossprod(hh1,mm)
  a1<-B1[2]    # da.
  a2<-B1[3]    # db.
  a3<-B1[4]    # ha.
  a4<-B1[5]    # hb.
  g<-matrix(0,3,1)
  g[1]<-(0.5*a2^2+0.25*a4^2)/num_l;g[2]<-(0.5*a1^2+0.25*a3^2)/num_l;g[3]<-g[1]+g[2]
  sigma[2]<-sigma[1]+g[1];sigma[3]<-sigma[1]+g[2];sigma[4]<-sigma[1]+g[3]
  sigma[5]<-sigma[8]+g[3];sigma[6]<-sigma[8]+g[2];sigma[7]<-sigma[8]+g[1]
  mm1<-m[c(1:4)];ssigma1<-sigma[c(1:4)];mmi1<-mi[c(1:4)]
  mm2<-m[c(5:8)];ssigma2<-sigma[c(5:8)];mmi2<-mi[c(5:8)]
  L0<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
    sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)))+
    sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,4,nn[4]); swx1 <- matrix(0,4,1)
  W2 <- matrix(0,4,nn[5]); swx2 <- matrix(0,4,1)
  sumwx<-matrix(0,8,1);mix_pi<-matrix(0,8,1)
  n0<-matrix(0,8,1);g<-matrix(0,3,1);s0<-matrix(0,6,1)
  hh<-matrix(0,4,4);b_line<-matrix(0,4,1)
  swx<-matrix(0,8,1);aa3<-matrix(0,8,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:4) { W1[i,] <- mi[i]*dnorm(dataB1,m[i],sqrt(sigma[i]))/dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)}
    mix_pi[c(1:4)] <- as.matrix(rowSums(W1)/nn[4])
    sumwx[c(1:4)] <- W1%*%dataB1
    for(i in 1:4) { W2[i,] <- mmi2[i]*dnorm(dataB2,mm2[i],sqrt(ssigma2[i]))/dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)}
    mix_pi[c(5:8)] <- as.matrix(rowSums(W2)/nn[5])
    sumwx[c(5:8)] <- W2%*%dataB2
    n0[c(1:4)]<-nn[4]*mix_pi[c(1:4)];n0[c(5:8)]<-nn[5]*mix_pi[c(5:8)]
    n0[c(1:8)][abs(n0[c(1:8)])<0.000001]<-0.000001
    ##########################################
    hh[1,1]<-25*sigma[11]/nn[1]+4*sigma[11]/nn[2]+sigma[11]/nn[3]+16*sigma[1]/n0[1]+16*sigma[4]/n0[4]
    hh[1,2]<-5*sigma[11]/nn[1]+4*sigma[11]/nn[2]+5*sigma[11]/nn[3]
    hh[1,3]<--4*sigma[1]/n0[1]-4*sigma[4]/n0[4]
    hh[1,4]<-0
    hh[2,2]<-sigma[11]/nn[1]+4*sigma[11]/nn[2]+25*sigma[11]/nn[3]+16*sigma[6]/n0[6]+16*sigma[7]/n0[7]
    hh[2,3]<-0
    hh[2,4]<-4*sigma[6]/n0[6]+4*sigma[7]/n0[7]
    hh[3,3]<-sigma[1]/n0[1]+sigma[2]/n0[2]+sigma[3]/n0[3]+sigma[4]/n0[4]
    hh[3,4]<-0
    hh[4,4]<-sigma[5]/n0[5]+sigma[6]/n0[6]+sigma[7]/n0[7]+sigma[8]/n0[8]
    for(i in 2:4)
    {
      for(j in 1:(i-1))
      {
        hh[i,j]<-hh[j,i]
      }
    }
    #############################################
    b_line[1]<-5*sumx[1]/nn[1]+2*sumx[2]/nn[2]+sumx[3]/nn[3]-4*sumwx[1]/n0[1]-4*sumwx[4]/n0[4]                          
    b_line[2]<-sumx[1]/nn[1]+2*sumx[2]/nn[2]+5*sumx[3]/nn[3]-4*sumwx[6]/n0[6]-4*sumwx[7]/n0[7]
    b_line[3]<-sumwx[1]/n0[1]-sumwx[2]/n0[2]-sumwx[3]/n0[3]+sumwx[4]/n0[4]
    b_line[4]<-sumwx[5]/n0[5]-sumwx[6]/n0[6]-sumwx[7]/n0[7]+sumwx[8]/n0[8] 
    B2<-solve(hh,b_line) 
    m[11]<-(sumx[1]-sigma[11]*(5*B2[1]+B2[2]))/nn[1]
    m[12]<-(sumx[2]-sigma[11]*(2*B2[1]+B2[2]*2))/nn[2]
    m[13]<-(sumx[3]-sigma[11]*(B2[1]+5*B2[2]))/nn[3]
    m[1]<-(sumwx[1]+sigma[1]*(4*B2[1]-B2[3]))/n0[1]
    m[2]<-(sumwx[2]+sigma[2]*B2[3])/n0[2]
    m[3]<-(sumwx[3]+sigma[3]*B2[3])/n0[3]
    m[4]<-(sumwx[4]+sigma[4]*(4*B2[1]-B2[3]))/n0[4]
    m[5]<-(sumwx[5]-sigma[5]*B2[4])/n0[5]
    m[6]<-(sumwx[6]+sigma[6]*(4*B2[2]+B2[4]))/n0[6]
    m[7]<-(sumwx[7]+sigma[7]*(4*B2[2]+B2[4]))/n0[7]
    m[8]<-(sumwx[8]-sigma[8]*B2[4])/n0[8]
    hh3<-hh1
    mm<-m[c(11,12,13,1,2,3,4,5,6,7,8)] 
    B3<-solve(crossprod(hh3,hh3))%*%crossprod(hh3,mm)
    a1<-B3[2]    # da.
    a2<-B3[3]    # db.
    a3<-B3[4]    # ha.
    a4<-B3[5]    # hb.
    g[1]<-(0.5*a2^2+0.25*a4^2)/num_l;g[2]<-(0.5*a1^2+0.25*a3^2)/num_l;g[3]<-g[1]+g[2]
    sigma[2]<-sigma[1]+g[1];sigma[3]<-sigma[1]+g[2];sigma[4]<-sigma[1]+g[3]
    sigma[5]<-sigma[8]+g[3];sigma[6]<-sigma[8]+g[2];sigma[7]<-sigma[8]+g[1]
    ########obtain variance#######################
    ss1<-sum((dataP1-m[11])^2);ss3<-sum((dataP2-m[13])^2); ss2<-sum((dataF1-m[12])^2)
    mm2<-m[c(5:8)]
    for(i in 1:d4) {swx1[i] <- W1[i,]%*%(dataB1-m[i])^2 } ;for(i in 1:d4) {swx2[i] <- W2[i,]%*%(dataB2-mm2[i])^2 }  
    swx[c(1:4)]<-swx1;swx[c(5:8)]<-swx2
    aaa0<-sigma[1];n_iter<-0;abc1<-1000
    while (abc1>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma[1]/(sigma[1]+g[1])
      aa2<-sigma[1]/(sigma[1]+g[2])
      aa3<-sigma[1]/(sigma[1]+g[3])
      sigma[1]<-(swx[1]+aa1^2*swx[2]+aa2^2*swx[3]+aa3^2*swx[4])/(n0[1]+aa1*n0[2]+aa2*n0[3]+aa3*n0[4])
      abc1<-abs(sigma[1]-aaa0)
      aaa0<-sigma[1]
      if (n_iter>20) break
    } 
    sigma40<-sigma[1]-sigma[11]
    if (sigma40<0) {sigma40<-0;sigma[1]<-sigma[11]}
    sigma[2]<-sigma[1]+g[1];sigma[3]<-sigma[1]+g[2];sigma[4]<-sigma[1]+g[3]
    aaa0<-sigma[8];aa4<-1000;n_iter<-0
    while (aa4>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma[8]/(sigma[8]+g[1])
      aa2<-sigma[8]/(sigma[8]+g[2])
      aa3<-sigma[8]/(sigma[8]+g[3])
      sigma[8]<-(aa3^2*swx[5]+aa2^2*swx[6]+aa1^2*swx[7]+swx[8])/(aa3*n0[5]+aa2*n0[6]+aa1*n0[7]+n0[8])
      aa4<-abs(sigma[8]-aaa0)
      aaa0<-sigma[8]
      if (n_iter>20) break
    } 
    sigma50<-sigma[8]-sigma[11]
    if (sigma50<0) {sigma50<-0;sigma[8]<-sigma[11]}
    sigma[5]<-sigma[8]+g[3];sigma[6]<-sigma[8]+g[2];sigma[7]<-sigma[8]+g[1];
    aa1<-sigma40;aa2<-sigma50
    if (aa1<0) {aa1<-0};if (aa2<0) {aa2<-0}
    s0[1]<-ss1+ss2+ss3;s0[2]<-nn[1]+nn[2]+nn[3]
    aaa0<-sigma[11];n_iter<-0;abc2<-1000
    while (abc2>0.0001)
    {
      n_iter<-n_iter+1
      aa3[1]<-sigma[11]/(sigma[11]+aa1)
      aa3[c(2:4)]<-sigma[11]/(sigma[11]+aa1+g[c(1:3)])
      s0[3]<-sum(aa3[c(1:4)]^2*swx[c(1:4)])
      s0[4]<-sum(aa3[c(1:4)]*n0[c(1:4)])
      aa3[8]<-sigma[11]/(sigma[11]+aa2)
      aa3[c(5:7)]<-sigma[11]/(sigma[11]+aa2+g[c(3:1)])
      s0[5]<-sum(aa3[c(5:8)]^2*swx[c(5:8)])
      s0[6]<-sum(aa3[c(5:8)]*n0[c(5:8)])
      sigma[11]<-(s0[1]+(s0[3]+s0[5]))/(s0[2]+s0[4]+s0[6])
      abc2<-abs(sigma[11]-aaa0)
      aaa0<-sigma[11]
      if (n_iter>20) break
    }
    sigma[1]<-sigma[11]+sigma40;sigma[c(2:4)]<-sigma[11]+sigma40+g[c(1:3)];sigma[c(5:7)]<-sigma[11]+sigma50+g[c(3:1)];sigma[8]<-sigma[11]+sigma50
    mm1<-m[c(1:4)];ssigma1<-sigma[c(1:4)]
    mm2<-m[c(5:8)];ssigma2<-sigma[c(5:8)]
    mix_pi1<-mix_pi[c(1:4)];mix_pi2<-mix_pi[c(5:8)]
    ########criteria for iterations to stop###########
    L1<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
      sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mix_pi2)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
    if(iteration>300)break
  }
  abc<-L0
  AIC<--2*abc+2*10
  ########first order genetic parameters###############
  hh4<-hh1
  mm<-m[c(11,12,13,1,2,3,4,5,6,7,8)] 
  B4<-solve(t(hh4)%*%hh4)%*%(t(hh4)%*%mm)
  ########second order genetic parameters###############   
  B1jj <- sigma1 - sigma[1]	   
  if(B1jj<0) {B1jj<-0}	   
  B1ll <- B1jj/sigma1
  B1_gg <- sigma[1] - sigma[11]
  if(B1_gg<0|| B1_gg>sigma1){B1_gg<-0}	   
  B1rr <- B1_gg/sigma1	   
  B2jj <- sigma2 - sigma[8]
  if(B2jj<0) {B2jj<-0}
  B2ll <- B2jj/sigma2
  B2_gg <- sigma[8] - sigma[11]
  if(B2_gg<0 || B2_gg>sigma2) {B2_gg<-0}	   
  B2rr <- B2_gg/sigma2
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*nn[1])
  P1bmw <- matrix(0,nn[1],1)
  P1gg <- (dataP1 - m[11])/sqrt(as.vector(sigma[11]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn1 < nn[1]){P1bmw <- P1bmw+runif(nn[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:nn[1])) - 0.5)/nn[1])^2)
  P1u<- as.matrix(c(12*nn[1]*((P1dd[1]/nn[1]-0.5)^2),((45*nn[1])/4)*((P1dd[2]/nn[1]-1/3)^2),180*nn[1]*((P1dd[3]/nn[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,nn[1]))))
  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*nn[2])
  F1bmw <- matrix(0,nn[2],1)
  F1gg <- (dataF1 - m[12])/sqrt(as.vector(sigma[11]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn1 < nn[2]){F1bmw <- F1bmw+runif(nn[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:nn[2])) - 0.5)/nn[2])^2)
  F1u<- as.matrix(c(12*nn[2]*((F1dd[1]/nn[2]-0.5)^2),((45*nn[2])/4)*((F1dd[2]/nn[2]-1/3)^2),180*nn[2]*((F1dd[3]/nn[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,nn[2]))))
  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)
  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*nn[3])
  P2bmw <- matrix(0,nn[3],1)
  P2gg <- (dataP2 - m[13])/sqrt(as.vector(sigma[11]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn1 < nn[3]){P2bmw <- P2bmw+runif(nn[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:nn[3])) - 0.5)/nn[3])^2)
  P2u<- as.matrix(c(12*nn[3]*((P2dd[1]/nn[3]-0.5)^2),((45*nn[3])/4)*((P2dd[2]/nn[3]-1/3)^2),180*nn[3]*((P2dd[3]/nn[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,nn[3]))))
  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1); 
  B1w1<-1/(12*nn[4])
  B1bmw <- matrix(0,nn[4],1); B1bmwsl <- matrix(0,nn[4],d4)
  for(i in 1:d4){
    B1gg <- (dataB1 - m[i])/sqrt(sigma[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B1P2)))[1]
  if(nn1 < nn[4]){B1P2 <- B1P2+runif(nn[4])/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*nn[4]) + sum((B1P2 - (as.matrix(c(1:nn[4])) - 0.5)/nn[4])^2)
  B1u <- as.matrix(c(12*nn[4]*((B1dd[1]/nn[4]-0.5)^2),((45*nn[4])/4)*((B1dd[2]/nn[4]-1/3)^2),180*nn[4]*((B1dd[3]/nn[4]-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),(1-pkolm(B1D,nn[4]))))
  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2); 
  B2w1<-1/(12*nn[5])
  B2bmw <- matrix(0,nn[5],1); B2bmwsl <- matrix(0,nn[5],d4)
  for(i in 1:d4){
    B2gg <- (dataB2 - mm2[i])/sqrt(ssigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B2P2)))[1]
  if(nn1 < nn[5]){B2P2 <- B2P2+runif(nn[5])/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*nn[5]) + sum((B2P2 - (as.matrix(c(1:nn[5])) - 0.5)/nn[5])^2)
  B2u <- as.matrix(c(12*nn[5]*((B2dd[1]/nn[5]-0.5)^2),((45*nn[5])/4)*((B2dd[2]/nn[5]-1/3)^2),180*nn[5]*((B2dd[3]/nn[5]-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),(1-pkolm(B2D,nn[5]))))
  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("MX2-AD-AD",round(abc,4),round(AIC,4),round(m[11],4),round(m[12],4),round(m[13],4),round(m[1],4),round(m[2],4),round(m[3],4),round(m[4],4),round(sigma[1],4),round(sigma[2],4),
                       round(sigma[3],4),round(sigma[4],4),round(mix_pi[1],4),round(mix_pi[2],4),round(mix_pi[3],4),round(mix_pi[4],4),round(m[5],4),round(m[6],4),round(m[7],4),round(m[8],4),round(sigma[5],4),round(sigma[6],4), round(sigma[7],4),round(sigma[8],4),      
                       round(mix_pi[5],4),round(mix_pi[6],4),round(mix_pi[7],4),round(mix_pi[8],4),round(sigma[11],4),round(B4[1],4),round(B4[2],4),round(B4[3],4),round(B4[4],4),round(B4[5],4),round(B4[6],4),round(B4[7],4),round(B1jj,4),round(B1ll*100,4),round(B1_gg,4),round(B1rr*100,4),round(B2jj,4),round(B2ll*100,4),round(B2_gg,4),round(B2rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,as.matrix(mmi1),as.matrix(mmi2))
  return(OUTPUT)
}

#############MX2-A-AD(E-3)#########################################
G5BCFModelFun[[15]] <- function(K1,logL,df11,df21,df31,df41,df51,G5BCFtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]))
  nn<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataB1)[1],dim(dataB2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2)))
  ss<-s[c(1:3)]-sumx[c(1:3)]^2/nn[c(1:3)]
  m<-matrix(0,13,1);m[11]<-mean(dataP1);m[12]<-mean(dataF1);m[13]<-mean(dataP2);man0<-mean(dataB1);man1<-mean(dataB2)
  sigma0<-sum(ss)/(nn[1]+nn[2]+nn[3]-3)
  sigma1<-var(dataB1);sigma2<-var(dataB2)
  d1<-1;d2<-2;d3<-3;d4<-4
  m_esp <- 0.0001;num_l <- as.numeric(G5BCFtext2)
  ###############procedure start###########################  
  mi<-matrix(0.25,8,1)
  a1<-sqrt(sigma1/(nn[4]-1))
  m[1]<-man0+2*a1;m[2]<-man0+0.8*a1;m[3]<-man0-0.8*a1;m[4]<-man0-2*a1
  a1<-sqrt(sigma2/(nn[5]-1))
  m[5]<-man1+2*a1;m[6]<-man1+0.8*a1;m[7]<-man1-0.8*a1;m[8]<-man1-2*a1
  sigma<-matrix(0,11,1)
  sigma[11]<-sigma0;sigma[1]<-sigma1/3;sigma[8]<-sigma2/3
  hh1<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,0,-1,1,1,0,0,0,0,-1,-1,
                1,0,-1,1,0,1,0,0,-1,0,-1,1,0,-1,0.5,0.5,0.5,0.5,-0.5,-0.5,-0.5,-0.5,
                0,1,0,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25),11,5)
  mm<-m[c(11,12,13,1,2,3,4,5,6,7,8)]
  B1<-solve(crossprod(hh1,hh1))%*%crossprod(hh1,mm)
  a1<-B1[2]    # da.
  a2<-B1[3]    # db.
  g<-matrix(0,3,1)
  g[1]<-(0.5*a2^2)/num_l;g[2]<-(0.5*a1^2)/num_l;g[3]<-g[1]+g[2]
  sigma[2]<-sigma[1]+g[1];sigma[3]<-sigma[1]+g[2];sigma[4]<-sigma[1]+g[3]
  sigma[5]<-sigma[8]+g[3];sigma[6]<-sigma[8]+g[2];sigma[7]<-sigma[8]+g[1]
  mm1<-m[c(1:4)];ssigma1<-sigma[c(1:4)];mmi1<-mi[c(1:4)]
  mm2<-m[c(5:8)];ssigma2<-sigma[c(5:8)];mmi2<-mi[c(5:8)]
  L0<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
    sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)))+
    sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,4,nn[4]); swx1 <- matrix(0,4,1)
  W2 <- matrix(0,4,nn[5]); swx2 <- matrix(0,4,1)
  sumwx<-matrix(0,8,1);mix_pi<-matrix(0,8,1)
  hh<-matrix(0,6,6);b_line<-matrix(0,6,1)
  n0<-matrix(0,8,1);swx<-matrix(0,8,1);
  s0<-matrix(0,6,1);aa3<-matrix(0,8,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:4) { W1[i,] <- mi[i]*dnorm(dataB1,m[i],sqrt(sigma[i]))/dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)}
    mix_pi[c(1:4)] <- as.matrix(rowSums(W1)/nn[4])
    sumwx[c(1:4)] <- W1%*%dataB1
    for(i in 1:4) { W2[i,] <- mmi2[i]*dnorm(dataB2,mm2[i],sqrt(ssigma2[i]))/dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)}
    mix_pi[c(5:8)] <- as.matrix(rowSums(W2)/nn[5])
    sumwx[c(5:8)] <- W2%*%dataB2
    n0[c(1:4)]<-nn[4]*mix_pi[c(1:4)];n0[c(5:8)]<-nn[5]*mix_pi[c(5:8)]
    n0[c(1:8)][abs(n0[c(1:8)])<0.000001]<-0.000001
    ##########################################
    hh[1,1]<-25*sigma[11]/nn[1]+4*sigma[11]/nn[2]+sigma[11]/nn[3]+16*sigma[1]/n0[1]+16*sigma[4]/n0[4]
    hh[1,2]<-5*sigma[11]/nn[1]+4*sigma[11]/nn[2]+5*sigma[11]/nn[3]
    hh[1,3]<--4*sigma[1]/n0[1]-4*sigma[4]/n0[4]
    hh[1,4]<-0
    hh[1,5]<-0
    hh[1,6]<--4*sigma[1]/n0[1]+4*sigma[4]/n0[4]
    hh[2,2]<-sigma[11]/nn[1]+4*sigma[11]/nn[2]+25*sigma[11]/nn[3]+16*sigma[6]/n0[6]+16*sigma[7]/n0[7]
    hh[2,3]<-0
    hh[2,4]<-4*sigma[6]/n0[6]+4*sigma[7]/n0[7]
    hh[2,5]<-4*sigma[6]/n0[6]-4*sigma[7]/n0[7]
    hh[2,6]<-0
    hh[3,3]<-sigma[1]/n0[1]+sigma[2]/n0[2]+sigma[3]/n0[3]+sigma[4]/n0[4]
    hh[3,4]<-0
    hh[3,5]<--sigma[2]/n0[2]+sigma[3]/n0[3]
    hh[3,6]<-sigma[1]/n0[1]-sigma[4]/n0[4]
    hh[4,4]<-sigma[5]/n0[5]+sigma[6]/n0[6]+sigma[7]/n0[7]+sigma[8]/n0[8]
    hh[4,5]<-sigma[6]/n0[6]-sigma[7]/n0[7]
    hh[4,6]<--sigma[5]/n0[5]+sigma[8]/n0[8]
    hh[5,5]<-sigma[2]/n0[2]+sigma[3]/n0[3]+sigma[6]/n0[6]+sigma[7]/n0[7]
    hh[5,6]<-0
    hh[6,6]<-sigma[1]/n0[1]+sigma[4]/n0[4]+sigma[5]/n0[5]+sigma[8]/n0[8]
    for(i in 2:6)
    {
      for(j in 1:(i-1))
      {
        hh[i,j]<-hh[j,i]
      }
    }
    #########################################################
    b_line[1]<-5*sumx[1]/nn[1]+2*sumx[2]/nn[2]+sumx[3]/nn[3]-4*sumwx[1]/n0[1]-4*sumwx[4]/n0[4]
    b_line[2]<-sumx[1]/nn[1]+2*sumx[2]/nn[2]+5*sumx[3]/nn[3]-4*sumwx[6]/n0[6]-4*sumwx[7]/n0[7]
    b_line[3]<-sumwx[1]/n0[1]-sumwx[2]/n0[2]-sumwx[3]/n0[3]+sumwx[4]/n0[4]
    b_line[4]<-sumwx[5]/n0[5]-sumwx[6]/n0[6]-sumwx[7]/n0[7]+sumwx[8]/n0[8]
    b_line[5]<-sumwx[2]/n0[2]-sumwx[3]/n0[3]-sumwx[6]/n0[6]+sumwx[7]/n0[7]
    b_line[6]<-sumwx[1]/n0[1]-sumwx[4]/n0[4]-sumwx[5]/n0[5]+sumwx[8]/n0[8]
    B2<-solve(hh,b_line) 
    #########################################################
    m[11]<-(sumx[1]-sigma[11]*(5*B2[1]+B2[2]))/nn[1]
    m[12]<-(sumx[2]-sigma[11]*(2*B2[1]+B2[2]*2))/nn[2]
    m[13]<-(sumx[3]-sigma[11]*(B2[1]+5*B2[2]))/nn[3]
    m[1]<-(sumwx[1]+sigma[1]*(4*B2[1]-B2[3]-B2[6]))/n0[1]
    m[2]<-(sumwx[2]+sigma[2]*(B2[3]-B2[5]))/n0[2]
    m[3]<-(sumwx[3]+sigma[3]*(B2[3]+B2[5]))/n0[3]
    m[4]<-(sumwx[4]+sigma[4]*(4*B2[1]-B2[3]+B2[6]))/n0[4]
    m[5]<-(sumwx[5]+sigma[5]*(-B2[4]+B2[6]))/n0[5]
    m[6]<-(sumwx[6]+sigma[6]*(4*B2[2]+B2[4]+B2[5]))/n0[6]
    m[7]<-(sumwx[7]+sigma[7]*(4*B2[2]+B2[4]-B2[5]))/n0[7]
    m[8]<-(sumwx[8]+sigma[8]*(-B2[4]-B2[6]))/n0[8]
    hh3<-hh1
    mm<-m[c(11,12,13,1,2,3,4,5,6,7,8)]
    B3<-solve(crossprod(hh3,hh3))%*%crossprod(hh3,mm)
    a1<-B3[2]    # da.
    a2<-B3[3]    # db.
    g[1]<-(0.5*a2^2)/num_l;g[2]<-(0.5*a1^2)/num_l;g[3]<-g[1]+g[2]
    sigma[2]<-sigma[1]+g[1];sigma[3]<-sigma[1]+g[2];sigma[4]<-sigma[1]+g[3]
    sigma[5]<-sigma[8]+g[3];sigma[6]<-sigma[8]+g[2];sigma[7]<-sigma[8]+g[1]
    ########obtain variance#######################
    ss1<-sum((dataP1-m[11])^2);ss3<-sum((dataP2-m[13])^2);ss2<-sum((dataF1-m[12])^2)
    mm2<-m[c(5:8)]
    for(i in 1:d4) {swx1[i] <- W1[i,]%*%(dataB1-m[i])^2 }   ;for(i in 1:d4) {swx2[i] <- W2[i,]%*%(dataB2-mm2[i])^2 }  
    swx[c(1:4)]<-swx1;swx[c(5:8)]<-swx2
    aaa0<-sigma[1];n_iter<-0;abc1<-1000
    while (abc1>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma[1]/(sigma[1]+g[1])
      aa2<-sigma[1]/(sigma[1]+g[2])
      aa3<-sigma[1]/(sigma[1]+g[3])
      sigma[1]<-(swx[1]+aa1^2*swx[2]+aa2^2*swx[3]+aa3^2*swx[4])/(n0[1]+aa1*n0[2]+aa2*n0[3]+aa3*n0[4])
      abc1<-abs(sigma[1]-aaa0)
      aaa0<-sigma[1]
      if (n_iter>20) break
    } 
    sigma40<-sigma[1]-sigma[11]
    if (sigma40<0) {sigma40<-0;sigma[1]<-sigma[11]}
    sigma[2]<-sigma[1]+g[1];sigma[3]<-sigma[1]+g[2];sigma[4]<-sigma[1]+g[3]
    aaa0<-sigma[8]; aa4<-1000;n_iter<-0
    while (aa4>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma[8]/(sigma[8]+g[1])
      aa2<-sigma[8]/(sigma[8]+g[2])
      aa3<-sigma[8]/(sigma[8]+g[3])
      sigma[8]<-(aa3*aa3*swx[5]+aa2*aa2*swx[6]+aa1*aa1*swx[7]+swx[8])/(aa3*n0[5]+aa2*n0[6]+aa1*n0[7]+n0[8])
      aa4<-abs(sigma[8]-aaa0)
      aaa0<-sigma[8]
      if (n_iter>20) break
    } 
    sigma50<-sigma[8]-sigma[11]
    if (sigma50<0) {sigma50<-0;sigma[8]<-sigma[11]}
    sigma[5]<-sigma[8]+g[3];sigma[6]<-sigma[8]+g[2];sigma[7]<-sigma[8]+g[1]
    aa1<-sigma40;aa2<-sigma50
    if (aa1<0) {aa1<-0};if (aa2<0) {aa2<-0}
    s0[1]<-ss1+ss2+ss3;s0[2]<-nn[1]+nn[2]+nn[3]
    aaa0<-sigma[11];n_iter<-0;abc2<-1000
    while (abc2>0.0001)
    {
      n_iter<-n_iter+1
      aa3[1]<-sigma[11]/(sigma[11]+aa1)
      aa3[c(2:4)]<-sigma[11]/(sigma[11]+aa1+g[c(1:3)])
      s0[3]<-sum(aa3[c(1:4)]^2*swx[c(1:4)])
      s0[4]<-sum(aa3[c(1:4)]*n0[c(1:4)])
      aa3[8]<-sigma[11]/(sigma[11]+aa2)
      aa3[c(5:7)]<-sigma[11]/(sigma[11]+aa2+g[c(3:1)])
      s0[5]<-sum(aa3[c(5:8)]^2*swx[c(5:8)])
      s0[6]<-sum(aa3[c(5:8)]*n0[c(5:8)])
      sigma[11]<-(s0[1]+(s0[3]+s0[5]))/(s0[2]+s0[4]+s0[6])
      abc2<-abs(sigma[11]-aaa0)
      aaa0<-sigma[11]
      if (n_iter>20) break
    }
    sigma[1]<-sigma[11]+sigma40;sigma[c(2:4)]<-sigma[11]+sigma40+g[c(1:3)]
    sigma[c(5:7)]<-sigma[11]+sigma50+g[c(3:1)];sigma[8]<-sigma[11]+sigma50
    mm1<-m[c(1:4)];ssigma1<-sigma[c(1:4)]
    mm2<-m[c(5:8)];ssigma2<-sigma[c(5:8)]
    mix_pi1<-mix_pi[c(1:4)];mix_pi2<-mix_pi[c(5:8)]
    ########criteria for iterations to stop########
    L1<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
      sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mix_pi2)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
    if(iteration>300)break
  }
  abc<-L0
  AIC<--2*abc+2*8
  ########first order genetic parameters###############
  hh4<-hh1
  mm<-m[c(11,12,13,1,2,3,4,5,6,7,8)] 
  B4<-solve(t(hh4)%*%hh4)%*%(t(hh4)%*%mm)
  ########second order genetic parameters###############   
  B1jj <- sigma1 - sigma[1]	   
  if(B1jj<0) {B1jj<-0}	   
  B1ll <- B1jj/sigma1
  B1_gg <- sigma[1] - sigma[11]
  if(B1_gg<0|| B1_gg>sigma1){B1_gg<-0}	   
  B1rr <- B1_gg/sigma1	   
  B2jj <- sigma2 - sigma[8]
  if(B2jj<0) {B2jj<-0}
  B2ll <- B2jj/sigma2
  B2_gg <- sigma[8] - sigma[11]
  if(B2_gg<0 || B2_gg>sigma2) {B2_gg<-0}	   
  B2rr <- B2_gg/sigma2
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*nn[1])
  P1bmw <- matrix(0,nn[1],1)
  P1gg <- (dataP1 - m[11])/sqrt(as.vector(sigma[11]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn1 < nn[1]){P1bmw <- P1bmw+runif(nn[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:nn[1])) - 0.5)/nn[1])^2)
  P1u<- as.matrix(c(12*nn[1]*((P1dd[1]/nn[1]-0.5)^2),((45*nn[1])/4)*((P1dd[2]/nn[1]-1/3)^2),180*nn[1]*((P1dd[3]/nn[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,nn[1]))))
  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*nn[2])
  F1bmw <- matrix(0,nn[2],1)
  F1gg <- (dataF1 - m[12])/sqrt(as.vector(sigma[11]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn1 < nn[2]){F1bmw <- F1bmw+runif(nn[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:nn[2])) - 0.5)/nn[2])^2)
  F1u<- as.matrix(c(12*nn[2]*((F1dd[1]/nn[2]-0.5)^2),((45*nn[2])/4)*((F1dd[2]/nn[2]-1/3)^2),180*nn[2]*((F1dd[3]/nn[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,nn[2]))))
  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)
  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*nn[3])
  P2bmw <- matrix(0,nn[3],1)
  P2gg <- (dataP2 - m[13])/sqrt(as.vector(sigma[11]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn1 < nn[3]){P2bmw <- P2bmw+runif(nn[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:nn[3])) - 0.5)/nn[3])^2)
  P2u<- as.matrix(c(12*nn[3]*((P2dd[1]/nn[3]-0.5)^2),((45*nn[3])/4)*((P2dd[2]/nn[3]-1/3)^2),180*nn[3]*((P2dd[3]/nn[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,nn[3]))))
  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1); 
  B1w1<-1/(12*nn[4])
  B1bmw <- matrix(0,nn[4],1); B1bmwsl <- matrix(0,nn[4],d4)
  for(i in 1:d4){
    B1gg <- (dataB1 - m[i])/sqrt(sigma[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B1P2)))[1]
  if(nn1 < nn[4]){B1P2 <- B1P2+runif(nn[4])/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*nn[4]) + sum((B1P2 - (as.matrix(c(1:nn[4])) - 0.5)/nn[4])^2)
  B1u <- as.matrix(c(12*nn[4]*((B1dd[1]/nn[4]-0.5)^2),((45*nn[4])/4)*((B1dd[2]/nn[4]-1/3)^2),180*nn[4]*((B1dd[3]/nn[4]-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),(1-pkolm(B1D,nn[4]))))
  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2); 
  B2w1<-1/(12*nn[5])
  B2bmw <- matrix(0,nn[5],1); B2bmwsl <- matrix(0,nn[5],d4)
  for(i in 1:d4){
    B2gg <- (dataB2 - mm2[i])/sqrt(ssigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B2P2)))[1]
  if(nn1 < nn[5]){B2P2 <- B2P2+runif(nn[5])/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*nn[5]) + sum((B2P2 - (as.matrix(c(1:nn[5])) - 0.5)/nn[5])^2)
  B2u <- as.matrix(c(12*nn[5]*((B2dd[1]/nn[5]-0.5)^2),((45*nn[5])/4)*((B2dd[2]/nn[5]-1/3)^2),180*nn[5]*((B2dd[3]/nn[5]-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),(1-pkolm(B2D,nn[5]))))
  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("MX2-A-AD",round(abc,4),round(AIC,4),round(m[11],4),round(m[12],4),round(m[13],4),round(m[1],4),round(m[2],4),round(m[3],4),round(m[4],4),round(sigma[1],4),round(sigma[2],4),
                       round(sigma[3],4),round(sigma[4],4),round(mix_pi[1],4),round(mix_pi[2],4),round(mix_pi[3],4),round(mix_pi[4],4),round(m[5],4),round(m[6],4),round(m[7],4),round(m[8],4),round(sigma[5],4),round(sigma[6],4), round(sigma[7],4),round(sigma[8],4),      
                       round(mix_pi[5],4),round(mix_pi[6],4),round(mix_pi[7],4),round(mix_pi[8],4),round(sigma[11],4),round(B4[1],4),round(B4[2],4),round(B4[3],4)," "," ",round(B4[4],4),round(B4[5],4),round(B1jj,4),round(B1ll*100,4),round(B1_gg,4),round(B1rr*100,4),round(B2jj,4),round(B2ll*100,4),round(B2_gg,4),round(B2rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,as.matrix(mmi1),as.matrix(mmi2))
  return(OUTPUT)
}

###############MX2-EA-AD(E-4)###############################
G5BCFModelFun[[16]] <- function(K1,logL,df11,df21,df31,df41,df51,G5BCFtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]))
  nn<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataB1)[1],dim(dataB2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2)))
  ss<-s[c(1:3)]-sumx[c(1:3)]^2/nn[c(1:3)]
  m<-matrix(0,13,1);m[11]<-mean(dataP1);m[12]<-mean(dataF1);m[13]<-mean(dataP2);man0<-mean(dataB1);man1<-mean(dataB2)
  sigma0<-sum(ss)/(nn[1]+nn[2]+nn[3]-3)
  sigma1<-var(dataB1);sigma2<-var(dataB2)
  d1<-1;d2<-2;d3<-3;d4<-4
  m_esp <- 0.0001;num_l <- as.numeric(G5BCFtext2)
  ###############procedure start########################### 
  mi<-matrix(0,7,1);mi[c(1,3)]<-0.25;mi[c(2,6)]<-0.5;mi[c(5,7)]<-0.25
  a1<-sqrt(sigma1/(nn[4]-1))
  m[1]<-man0+2*a1;m[2]<-man0;m[3]<-man0-2*a1
  a1<-sqrt(sigma2/(nn[5]-1))
  m[5]<-man1+2*a1;m[6]<-man1;m[7]<-man1-2*a1
  sigma<-matrix(0,11,1);sigma[11]<-sigma0;sigma[1]<-sigma1/3;sigma[c(7,8)]<-sigma2/3
  hh1<-matrix(c(1,1,1,1,1,1,1,1,1,2,0,-2,2,1,0,0,-1,-2,1,0,-1,0.5,0.5,0.5,-0.5,-0.5,-0.5,
                0,1,0,0.25,0.25,0.25,0.25,0.25,0.25),9,4)
  mm<-m[c(11,12,13,1,2,3,5,6,7)]
  B1<-solve(crossprod(hh1,hh1))%*%crossprod(hh1,mm)
  a1<-B1[2]    # d.
  g<-matrix(0,2,1)
  g[1]<-0.5*a1^2/num_l;g[2]<-a1^2/num_l
  sigma[2]<-sigma[1]+g[1];sigma[3]<-sigma[1]+g[2];sigma[5]<-sigma[7]+g[2];sigma[6]<-sigma[7]+g[1]
  mm1<-m[c(1:3)];ssigma1<-sigma[c(1:3)];mmi1<-mi[c(1:3)]
  mm2<-m[c(5:7)];ssigma2<-sigma[c(5:7)];mmi2<-mi[c(5:7)]
  L0<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
    sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)))+
    sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,3,nn[4]); swx1 <- matrix(0,3,1)
  W2 <- matrix(0,3,nn[5]); swx2 <- matrix(0,3,1)
  sumwx<-matrix(0,8,1);mix_pi<-matrix(0,8,1)
  hh<-matrix(0,5,5);b_line<-matrix(0,5,1);s0<-matrix(0,6,1)
  n0<-matrix(0,8,1);swx<-matrix(0,8,1); aa3<-matrix(0,8,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:3) { W1[i,] <- mi[i]*dnorm(dataB1,m[i],sqrt(sigma[i]))/dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)}
    mix_pi[c(1:3)] <- as.matrix(rowSums(W1)/nn[4])
    sumwx[c(1:3)] <- W1%*%dataB1
    for(i in 1:3) { W2[i,] <- mmi2[i]*dnorm(dataB2,mm2[i],sqrt(ssigma2[i]))/dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)}
    mix_pi[c(5:7)] <- as.matrix(rowSums(W2)/nn[5])
    sumwx[c(5:7)] <- W2%*%dataB2
    n0[c(1:3)]<-nn[4]*mix_pi[c(1:3)];n0[c(5:7)]<-nn[5]*mix_pi[c(5:7)]
    n0[c(1,2,3,5,6,7)][abs(n0[c(1,2,3,5,6,7)])<0.000001]<-0.000001
    ########################################################
    hh[1,1]<-25*sigma[11]/nn[1]+4*sigma[11]/nn[2]+sigma[11]/nn[3]+64*sigma[2]/n0[2]
    hh[1,2]<-5*sigma[11]/nn[1]+4*sigma[11]/nn[2]+5*sigma[11]/nn[3]
    hh[1,3]<-16*sigma[2]/n0[2]
    hh[1,4]<-0
    hh[1,5]<-0
    hh[2,2]<-sigma[11]/nn[1]+4*sigma[11]/nn[2]+25*sigma[11]/nn[3]+64*sigma[6]/n0[6]
    hh[2,3]<-0
    hh[2,4]<-16*sigma[6]/n0[6]
    hh[2,5]<-0
    hh[3,3]<-sigma[1]/n0[1]+4*sigma[2]/n0[2]+sigma[3]/n0[3]
    hh[3,4]<-0
    hh[3,5]<-sigma[1]/n0[1]-sigma[3]/n0[3]
    hh[4,4]<-sigma[5]/n0[5]+4*sigma[6]/n0[6]+sigma[7]/n0[7]
    hh[4,5]<--sigma[5]/n0[5]+sigma[7]/n0[7]
    hh[5,5]<-sigma[1]/n0[1]+sigma[3]/n0[3]+sigma[5]/n0[5]+sigma[7]/n0[7]
    for(i in 2:5)
    {
      for(j in 1:(i-1))
      {
        hh[i,j]<-hh[j,i]
      }
    }
    #########################################
    b_line[1]<-5*sumx[1]/nn[1]+2*sumx[2]/nn[2]+sumx[3]/nn[3]-8*sumwx[2]/n0[2];
    b_line[2]<-sumx[1]/nn[1]+2*sumx[2]/nn[2]+5*sumx[3]/nn[3]-8*sumwx[6]/n0[6];
    b_line[3]<-sumwx[1]/n0[1]-2*sumwx[2]/n0[2]+sumwx[3]/n0[3];
    b_line[4]<-sumwx[5]/n0[5]-2*sumwx[6]/n0[6]+sumwx[7]/n0[7];
    b_line[5]<-sumwx[1]/n0[1]-sumwx[3]/n0[3]-sumwx[5]/n0[5]+sumwx[7]/n0[7];
    B2<-solve(hh,b_line)
    #############################################
    m[11]<-(sumx[1]-sigma[11]*(5*B2[1]+B2[2]))/nn[1]
    m[12]<-(sumx[2]-sigma[11]*(2*B2[1]+B2[2]*2))/nn[2]
    m[13]<-(sumx[3]-sigma[11]*(B2[1]+5*B2[2]))/nn[3]
    m[1]<-(sumwx[1]-sigma[1]*(B2[3]+B2[5]))/n0[1]
    m[2]<-(sumwx[2]+sigma[2]*(8*B2[1]+2*B2[3]))/n0[2]
    m[3]<-(sumwx[3]-sigma[3]*(B2[3]-B2[5]))/n0[3]
    m[5]<-(sumwx[5]-sigma[5]*(B2[4]-B2[5]))/n0[5]
    m[6]<-(sumwx[6]+sigma[6]*(8*B2[2]+2*B2[4]))/n0[6]
    m[7]<-(sumwx[7]-sigma[7]*(B2[4]+B2[5]))/n0[7]
    ##############################################
    hh3<-hh1
    mm<-m[c(11,12,13,1,2,3,5,6,7)]
    B3<-solve(crossprod(hh3,hh3))%*%crossprod(hh3,mm)
    a1<-B3[2]    # d.
    g[1]<-0.5*a1^2/num_l;g[2]<-a1^2/num_l
    sigma[2]<-sigma[1]+g[1];sigma[3]<-sigma[1]+g[2]
    ########obtain variance#######################
    ss1<-sum((dataP1-m[11])^2);ss3<-sum((dataP2-m[13])^2);ss2<-sum((dataF1-m[12])^2)
    mm2<-m[c(5:7)]
    for(i in 1:d3) {swx1[i] <- W1[i,]%*%(dataB1-m[i])^2 } ;for(i in 1:d3) {swx2[i] <- W2[i,]%*%(dataB2-mm2[i])^2 }  
    swx[c(1:3)]<-swx1;swx[c(5:7)]<-swx2
    aaa0<-sigma[1];n_iter<-0;abc1<-1000
    while (abc1>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma[1]/(sigma[1]+g[1])
      aa2<-sigma[1]/(sigma[1]+g[2])
      sigma[1]<-(swx[1]+aa1^2*swx[2]+aa2^2*swx[3])/(n0[1]+aa1*n0[2]+aa2*n0[3])
      abc1<-abs(sigma[1]-aaa0)
      aaa0<-sigma[1]
      if (n_iter>20) break
    } 
    sigma40<-sigma[1]-sigma[11]
    if (sigma40<0) {sigma40<-0;sigma[1]<-sigma[11]}
    sigma[2]<-sigma[1]+g[1];sigma[3]<-sigma[1]+g[2]
    aaa0<-sigma[7];aa4<-1000;n_iter<-0
    while (aa4>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma[7]/(sigma[7]+g[1])
      aa2<-sigma[7]/(sigma[7]+g[2])
      sigma[7]<-(aa2^2*swx[5]+aa1^2*swx[6]+swx[7])/(aa2*n0[5]+aa1*n0[6]+n0[7])
      aa4<-abs(sigma[7]-aaa0)
      aaa0<-sigma[7]
      if (n_iter>20) break
    } 
    sigma50<-sigma[7]-sigma[11]
    if (sigma50<0) {sigma50<-0;sigma[7]<-sigma[11]}
    sigma[5]<-sigma[7]+g[2];sigma[6]<-sigma[7]+g[1]
    aa1<-sigma40;aa2<-sigma50
    if (aa1<0) {aa1<-0};if (aa2<0) {aa2<-0}
    s0[1]<-ss1+ss2+ss3;s0[2]<-nn[1]+nn[2]+nn[3]
    aaa0<-sigma[11];n_iter<-0;abc2<-1000
    while (abc2>0.0001)
    {
      n_iter<-n_iter+1
      aa3[1]<-sigma[11]/(sigma[11]+aa1)
      aa3[c(2:3)]<-sigma[11]/(sigma[11]+aa1+g[c(1:2)])
      s0[3]<-sum(aa3[c(1:4)]^2*swx[c(1:4)])
      s0[4]<-sum(aa3[c(1:4)]*n0[c(1:4)])
      aa3[7]<-sigma[11]/(sigma[11]+aa2)
      aa3[c(5:6)]<-sigma[11]/(sigma[11]+aa2+g[c(2:1)])
      s0[5]<-sum(aa3[c(5:7)]^2*swx[c(5:7)])
      s0[6]<-sum(aa3[c(5:7)]*n0[c(5:7)])
      sigma[11]<-(s0[1]+(s0[3]+s0[5]))/(s0[2]+s0[4]+s0[6])
      abc2<-abs(sigma[11]-aaa0)
      aaa0<-sigma[11]
      if (n_iter>20) break
    }
    sigma[1]<-sigma[11]+sigma40;sigma[c(2:3)]<-sigma[11]+sigma40+g[c(1:2)]
    sigma[c(5:6)]<-sigma[11]+sigma50+g[c(2:1)];sigma[7]<-sigma[11]+sigma50
    mm1<-m[c(1:3)];ssigma1<-sigma[c(1:3)]
    mm2<-m[c(5:7)];ssigma2<-sigma[c(5:7)]
    mix_pi1<-mix_pi[c(1:3)];mix_pi2<-mix_pi[c(5:7)]
    ########criteria for iterations to stop########
    L1<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
      sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mix_pi2)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
    if(iteration>300)break
  }
  abc<-L0
  AIC<--2*abc+2*7
  ########first order genetic parameters###############
  hh4<-hh1
  mm<-m[c(11,12,13,1,2,3,5,6,7)]
  B4<-solve(t(hh4)%*%hh4)%*%(t(hh4)%*%mm)
  ########second order genetic parameters###############   
  B1jj <- sigma1 - sigma[1]	   
  if(B1jj<0) {B1jj<-0}	   
  B1ll <- B1jj/sigma1
  B1_gg <- sigma[1] - sigma[11]
  if(B1_gg<0|| B1_gg>sigma1){B1_gg<-0}	   
  B1rr <- B1_gg/sigma1	   
  B2jj <- sigma2 - sigma[7]
  if(B2jj<0) {B2jj<-0}
  B2ll <- B2jj/sigma2
  B2_gg <- sigma[7] - sigma[11]
  if(B2_gg<0 || B2_gg>sigma2) {B2_gg<-0}	   
  B2rr <- B2_gg/sigma2
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*nn[1])
  P1bmw <- matrix(0,nn[1],1)
  P1gg <- (dataP1 - m[11])/sqrt(as.vector(sigma[11]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn1 < nn[1]){P1bmw <- P1bmw+runif(nn[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:nn[1])) - 0.5)/nn[1])^2)
  P1u<- as.matrix(c(12*nn[1]*((P1dd[1]/nn[1]-0.5)^2),((45*nn[1])/4)*((P1dd[2]/nn[1]-1/3)^2),180*nn[1]*((P1dd[3]/nn[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,nn[1]))))
  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*nn[2])
  F1bmw <- matrix(0,nn[2],1)
  F1gg <- (dataF1 - m[12])/sqrt(as.vector(sigma[11]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn1 < nn[2]){F1bmw <- F1bmw+runif(nn[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:nn[2])) - 0.5)/nn[2])^2)
  F1u<- as.matrix(c(12*nn[2]*((F1dd[1]/nn[2]-0.5)^2),((45*nn[2])/4)*((F1dd[2]/nn[2]-1/3)^2),180*nn[2]*((F1dd[3]/nn[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,nn[2]))))
  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)
  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*nn[3])
  P2bmw <- matrix(0,nn[3],1)
  P2gg <- (dataP2 - m[13])/sqrt(as.vector(sigma[11]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn1 < nn[3]){P2bmw <- P2bmw+runif(nn[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:nn[3])) - 0.5)/nn[3])^2)
  P2u<- as.matrix(c(12*nn[3]*((P2dd[1]/nn[3]-0.5)^2),((45*nn[3])/4)*((P2dd[2]/nn[3]-1/3)^2),180*nn[3]*((P2dd[3]/nn[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,nn[3]))))
  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1); 
  B1w1<-1/(12*nn[4])
  B1bmw <- matrix(0,nn[4],1); B1bmwsl <- matrix(0,nn[4],d3)
  for(i in 1:d3){
    B1gg <- (dataB1 - m[i])/sqrt(sigma[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B1P2)))[1]
  if(nn1 < nn[4]){B1P2 <- B1P2+runif(nn[4])/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*nn[4]) + sum((B1P2 - (as.matrix(c(1:nn[4])) - 0.5)/nn[4])^2)
  B1u <- as.matrix(c(12*nn[4]*((B1dd[1]/nn[4]-0.5)^2),((45*nn[4])/4)*((B1dd[2]/nn[4]-1/3)^2),180*nn[4]*((B1dd[3]/nn[4]-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),(1-pkolm(B1D,nn[4]))))
  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2); 
  B2w1<-1/(12*nn[5])
  B2bmw <- matrix(0,nn[5],1); B2bmwsl <- matrix(0,nn[5],d3)
  for(i in 1:d3){
    B2gg <- (dataB2 - mm2[i])/sqrt(ssigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B2P2)))[1]
  if(nn1 < nn[5]){B2P2 <- B2P2+runif(nn[5])/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*nn[5]) + sum((B2P2 - (as.matrix(c(1:nn[5])) - 0.5)/nn[5])^2)
  B2u <- as.matrix(c(12*nn[5]*((B2dd[1]/nn[5]-0.5)^2),((45*nn[5])/4)*((B2dd[2]/nn[5]-1/3)^2),180*nn[5]*((B2dd[3]/nn[5]-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),(1-pkolm(B2D,nn[5]))))
  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("MX2-EA-AD",round(abc,4),round(AIC,4),round(m[11],4),round(m[12],4),round(m[13],4),round(m[1],4),round(m[2],4),round(m[3],4)," ",round(sigma[1],4),round(sigma[2],4),
                       round(sigma[3],4), " ",round(mix_pi[1],4),round(mix_pi[2],4),round(mix_pi[3],4)," ",round(m[5],4),round(m[6],4),round(m[7],4)," ",round(sigma[5],4),round(sigma[6],4), round(sigma[7],4)," ",   
                       round(mix_pi[5],4),round(mix_pi[6],4),round(mix_pi[7],4)," ",round(sigma[11],4),round(B4[1],4),round(B4[2],4)," "," "," ",round(B4[3],4),round(B4[4],4),round(B1jj,4),round(B1ll*100,4),round(B1_gg,4),round(B1rr*100,4),round(B2jj,4),round(B2ll*100,4),round(B2_gg,4),round(B2rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,as.matrix(mmi1),as.matrix(mmi2))
  return(OUTPUT)
}

#################MX2-CD-AD(E-5)#################################################
G5BCFModelFun[[17]] <- function(K1,logL,df11,df21,df31,df41,df51,G5BCFtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]))
  nn<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataB1)[1],dim(dataB2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2)))
  ss<-s[c(1:3)]-sumx[c(1:3)]^2/nn[c(1:3)]
  m<-matrix(0,13,1);m[11]<-mean(dataP1);m[12]<-mean(dataF1);m[13]<-mean(dataP2);man0<-mean(dataB1);man1<-mean(dataB2)
  sigma0<-sum(ss)/(nn[1]+nn[2]+nn[3]-3)
  sigma1<-var(dataB1);sigma2<-var(dataB2)
  d1<-1;d2<-2;d3<-3;d4<-4
  m_esp <- 0.0001;num_l <- as.numeric(G5BCFtext2)
  ###############procedure start###########################
  mi<-matrix(0.25,8,1)
  a1<-sqrt(sigma1/(nn[4]-1))
  m[1]<-man0+2*a1;m[2]<-man0+0.8*a1
  m[3]<-man0-0.8*a1;m[4]<-man0-2*a1
  a1<-sqrt(sigma2/(nn[5]-1))
  m[5]<-man1+2*a1;m[6]<-man1+0.8*a1
  m[7]<-man1-0.8*a1;m[8]<-man1-2*a1
  sigma<-matrix(0,11,1)
  sigma[11]<-sigma0;sigma[1]<-sigma1/3;sigma[8]<-sigma2/3
  hh1<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,-1,1,1,0.5,0.5,0.5,0.5,-1,-1,
                1,1,-1,1,0.5,1,0.5,0.5,-1,0.5,-1,1,0,-1,0.5,0.5,0.5,0.5,-0.5,-0.5,-0.5,-0.5,
                0,1,0,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25),11,5)
  mm<-m[c(11,12,13,1,2,3,4,5,6,7,8)]
  B1<-solve(crossprod(hh1,hh1))%*%crossprod(hh1,mm)
  a1<-B1[2]    # da.
  a2<-B1[3]    # db.
  g<-matrix(0,3,1)
  g[1]<-(0.75*a2^2)/num_l;g[2]<-(0.75*a1^2)/num_l;g[3]<-g[1]+g[2]
  sigma[2]<-sigma[1]+g[1];sigma[3]<-sigma[1]+g[2];sigma[4]<-sigma[1]+g[3]
  sigma[5]<-sigma[8]+g[3];sigma[6]<-sigma[8]+g[2];sigma[7]<-sigma[8]+g[1]
  mm1<-m[c(1:4)];ssigma1<-sigma[c(1:4)];mmi1<-mi[c(1:4)]
  mm2<-m[c(5:8)];ssigma2<-sigma[c(5:8)];mmi2<-mi[c(5:8)]
  L0<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
    sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)))+
    sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,4,nn[4]); swx1 <- matrix(0,4,1)
  W2 <- matrix(0,4,nn[5]); swx2 <- matrix(0,4,1)
  sumwx<-matrix(0,8,1);mix_pi<-matrix(0,8,1)
  n0<-matrix(0,8,1);hh<-matrix(0,6,6);b_line<-matrix(0,6,1)
  swx<-matrix(0,8,1);s0<-matrix(0,6,1); aa3<-matrix(0,8,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:4) { W1[i,] <- mi[i]*dnorm(dataB1,m[i],sqrt(sigma[i]))/dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)}
    mix_pi[c(1:4)] <- as.matrix(rowSums(W1)/nn[4])
    sumwx[c(1:4)] <- W1%*%dataB1
    for(i in 1:4) { W2[i,] <- mmi2[i]*dnorm(dataB2,mm2[i],sqrt(ssigma2[i]))/dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)}
    mix_pi[c(5:8)] <- as.matrix(rowSums(W2)/nn[5])
    sumwx[c(5:8)] <- W2%*%dataB2
    n0[c(1:4)]<-nn[4]*mix_pi[c(1:4)];n0[c(5:8)]<-nn[5]*mix_pi[c(5:8)]
    n0[c(1:8)][abs(n0[c(1:8)])<0.000001]<-0.000001
    ##########################################
    hh[1,1]<-25*sigma[11]/nn[1]+4*sigma[11]/nn[2]+sigma[11]/nn[3]+16*sigma[1]/n0[1]+16*sigma[4]/n0[4]
    hh[1,2]<-5*sigma[11]/nn[1]+4*sigma[11]/nn[2]+5*sigma[11]/nn[3]
    hh[1,3]<--4*sigma[1]/n0[1]-4*sigma[4]/n0[4]
    hh[1,4]<-0
    hh[1,5]<-0
    hh[1,6]<--12*sigma[1]/n0[1]+12*sigma[4]/n0[4]
    hh[2,2]<-sigma[11]/nn[1]+4*sigma[11]/nn[2]+25*sigma[11]/nn[3]+16*sigma[6]/n0[6]+16*sigma[7]/n0[7]
    hh[2,3]<-0
    hh[2,4]<-4*sigma[6]/n0[6]+4*sigma[7]/n0[7]
    hh[2,5]<-4*sigma[6]/n0[6]-4*sigma[7]/n0[7]
    hh[2,6]<-0
    hh[3,3]<-sigma[1]/n0[1]+sigma[2]/n0[2]+sigma[3]/n0[3]+sigma[4]/n0[4]
    hh[3,4]<-0
    hh[3,5]<--3*sigma[2]/n0[2]+3*sigma[3]/n0[3]
    hh[3,6]<-3*sigma[1]/n0[1]-3*sigma[4]/n0[4]
    hh[4,4]<-sigma[5]/n0[5]+sigma[6]/n0[6]+sigma[7]/n0[7]+sigma[8]/n0[8]
    hh[4,5]<-sigma[6]/n0[6]-sigma[7]/n0[7]
    hh[4,6]<--sigma[5]/n0[5]+sigma[8]/n0[8]
    hh[5,5]<-9*sigma[2]/n0[2]+9*sigma[3]/n0[3]+sigma[6]/n0[6]+sigma[7]/n0[7]
    hh[5,6]<-0
    hh[6,6]<-9*sigma[1]/n0[1]+9*sigma[4]/n0[4]+sigma[5]/n0[5]+sigma[8]/n0[8]
    for(i in 2:6)
    {
      for(j in 1:(i-1))
      {
        hh[i,j]<-hh[j,i]
      }
    }
    ############################################################
    b_line[1]<-5*sumx[1]/nn[1]+2*sumx[2]/nn[2]+sumx[3]/nn[3]-4*sumwx[1]/n0[1]-4*sumwx[4]/n0[4]                            
    b_line[2]<-sumx[1]/nn[1]+2*sumx[2]/nn[2]+5*sumx[3]/nn[3]-4*sumwx[6]/n0[6]-4*sumwx[7]/n0[7]
    b_line[3]<-sumwx[1]/n0[1]-sumwx[2]/n0[2]-sumwx[3]/n0[3]+sumwx[4]/n0[4]
    b_line[4]<-sumwx[5]/n0[5]-sumwx[6]/n0[6]-sumwx[7]/n0[7]+sumwx[8]/n0[8]   
    b_line[5]<-3*sumwx[2]/n0[2]-3*sumwx[3]/n0[3]-sumwx[6]/n0[6]+sumwx[7]/n0[7]
    b_line[6]<-3*sumwx[1]/n0[1]-3*sumwx[4]/n0[4]-sumwx[5]/n0[5]+sumwx[8]/n0[8] 
    B2<-solve(hh,b_line) 
    #########################################################
    m[11]<-(sumx[1]-sigma[11]*(5*B2[1]+B2[2]))/nn[1]
    m[12]<-(sumx[2]-sigma[11]*(2*B2[1]+B2[2]*2))/nn[2]
    m[13]<-(sumx[3]-sigma[11]*(B2[1]+5*B2[2]))/nn[3]
    m[1]<-(sumwx[1]+sigma[1]*(4*B2[1]-B2[3]-3*B2[6]))/n0[1]
    m[2]<-(sumwx[2]+sigma[2]*(B2[3]-3*B2[5]))/n0[2]
    m[3]<-(sumwx[3]+sigma[3]*(B2[3]+3*B2[5]))/n0[3]
    m[4]<-(sumwx[4]+sigma[4]*(4*B2[1]-B2[3]+3*B2[6]))/n0[4]
    m[5]<-(sumwx[5]-sigma[5]*(B2[4]-B2[6]))/n0[5]
    m[6]<-(sumwx[6]+sigma[6]*(4*B2[2]+B2[4]+B2[5]))/n0[6]
    m[7]<-(sumwx[7]+sigma[7]*(4*B2[2]+B2[4]-B2[5]))/n0[7]
    m[8]<-(sumwx[8]-sigma[8]*(B2[4]+B2[6]))/n0[8]
    hh3<-hh1
    mm<-m[c(11,12,13,1,2,3,4,5,6,7,8)]
    B3<-solve(crossprod(hh3,hh3))%*%crossprod(hh3,mm)
    a1<-B3[2]    # da.
    a2<-B3[3]    # db.
    g[1]<-(0.75*a2^2)/num_l;g[2]<-(0.75*a1^2)/num_l;g[3]<-g[1]+g[2]
    sigma[2]<-sigma[1]+g[1];sigma[3]<-sigma[1]+g[2];sigma[4]<-sigma[1]+g[3]
    sigma[5]<-sigma[8]+g[3];sigma[6]<-sigma[8]+g[2];sigma[7]<-sigma[8]+g[1]
    ########obtain variance#######################
    ss1<-sum((dataP1-m[11])^2);ss3<-sum((dataP2-m[13])^2);ss2<-sum((dataF1-m[12])^2)
    mm2<-m[c(5:8)]
    for(i in 1:d4) {swx1[i] <- W1[i,]%*%(dataB1-m[i])^2 } ;for(i in 1:d4) {swx2[i] <- W2[i,]%*%(dataB2-mm2[i])^2 }  
    swx[c(1:4)]<-swx1;swx[c(5:8)]<-swx2
    aaa0<-sigma[1];n_iter<-0;abc1<-1000
    while (abc1>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma[1]/(sigma[1]+g[1])
      aa2<-sigma[1]/(sigma[1]+g[2])
      aa3<-sigma[1]/(sigma[1]+g[3])
      sigma[1]<-(swx[1]+aa1^2*swx[2]+aa2^2*swx[3]+aa3^2*swx[4])/(n0[1]+aa1*n0[2]+aa2*n0[3]+aa3*n0[4])
      abc1<-abs(sigma[1]-aaa0)
      aaa0<-sigma[1]
      if (n_iter>20) break
    } 
    sigma40<-sigma[1]-sigma[11]
    if (sigma40<0) {sigma40<-0;sigma[1]<-sigma[11]}
    sigma[2]<-sigma[1]+g[1];sigma[3]<-sigma[1]+g[2];sigma[4]<-sigma[1]+g[3]
    aaa0<-sigma[8];aa4<-1000;n_iter<-0
    while (aa4>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma[8]/(sigma[8]+g[1])
      aa2<-sigma[8]/(sigma[8]+g[2])
      aa3<-sigma[8]/(sigma[8]+g[3])
      sigma[8]<-(aa3^2*swx[5]+aa2^2*swx[6]+aa1^2*swx[7]+swx[8])/(aa3*n0[5]+aa2*n0[6]+aa1*n0[7]+n0[8])
      aa4<-abs(sigma[8]-aaa0)
      aaa0<-sigma[8]
      if (n_iter>20) break
    } 
    sigma50<-sigma[8]-sigma[11]
    if (sigma50<0) {sigma50<-0;sigma[8]<-sigma[11]}
    sigma[5]<-sigma[8]+g[3];sigma[6]<-sigma[8]+g[2];sigma[7]<-sigma[8]+g[1]
    aa1<-sigma40;aa2<-sigma50
    if (aa1<0) {aa1<-0};if (aa2<0) {aa2<-0}
    s0[1]<-ss1+ss2+ss3;s0[2]<-nn[1]+nn[2]+nn[3]
    aaa0<-sigma[11];n_iter<-0;abc2<-1000
    while (abc2>0.0001)
    {
      n_iter<-n_iter+1
      aa3[1]<-sigma[11]/(sigma[11]+aa1)
      aa3[c(2:4)]<-sigma[11]/(sigma[11]+aa1+g[c(1:3)])
      s0[3]<-sum(aa3[c(1:4)]^2*swx[c(1:4)])
      s0[4]<-sum(aa3[c(1:4)]*n0[c(1:4)])
      aa3[8]<-sigma[11]/(sigma[11]+aa2)
      aa3[c(5:7)]<-sigma[11]/(sigma[11]+aa2+g[c(3:1)])
      s0[5]<-sum(aa3[c(5:8)]^2*swx[c(5:8)])
      s0[6]<-sum(aa3[c(5:8)]*n0[c(5:8)])
      sigma[11]<-(s0[1]+(s0[3]+s0[5]))/(s0[2]+s0[4]+s0[6])
      abc2<-abs(sigma[11]-aaa0)
      aaa0<-sigma[11]
      if (n_iter>20) break
    }
    sigma[1]<-sigma[11]+sigma40;sigma[c(2:4)]<-sigma[11]+sigma40+g[c(1:3)];sigma[c(5:7)]<-sigma[11]+sigma50+g[c(3:1)];sigma[8]<-sigma[11]+sigma50
    mm1<-m[c(1:4)];ssigma1<-sigma[c(1:4)]
    mm2<-m[c(5:8)];ssigma2<-sigma[c(5:8)]
    mix_pi1<-mix_pi[c(1:4)];mix_pi2<-mix_pi[c(5:8)]
    L1<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
      sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mix_pi2)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
    if(iteration>300)break
  }
  abc<-L0
  AIC<--2*abc+2*8
  ########first order genetic parameters###############
  hh4<-hh1
  mm<-m[c(11,12,13,1,2,3,4,5,6,7,8)]
  B4<-solve(t(hh4)%*%hh4)%*%(t(hh4)%*%mm)
  ########second order genetic parameters###############   
  B1jj <- sigma1 - sigma[1]	   
  if(B1jj<0) {B1jj<-0}	   
  B1ll <- B1jj/sigma1
  B1_gg <- sigma[1] - sigma[11]
  if(B1_gg<0|| B1_gg>sigma1){B1_gg<-0}	   
  B1rr <- B1_gg/sigma1	   
  B2jj <- sigma2 - sigma[8]
  if(B2jj<0) {B2jj<-0}
  B2ll <- B2jj/sigma2
  B2_gg <- sigma[8] - sigma[11]
  if(B2_gg<0 || B2_gg>sigma2) {B2_gg<-0}	   
  B2rr <- B2_gg/sigma2
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*nn[1])
  P1bmw <- matrix(0,nn[1],1)
  P1gg <- (dataP1 - m[11])/sqrt(as.vector(sigma[11]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn1 < nn[1]){P1bmw <- P1bmw+runif(nn[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:nn[1])) - 0.5)/nn[1])^2)
  P1u<- as.matrix(c(12*nn[1]*((P1dd[1]/nn[1]-0.5)^2),((45*nn[1])/4)*((P1dd[2]/nn[1]-1/3)^2),180*nn[1]*((P1dd[3]/nn[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,nn[1]))))
  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*nn[2])
  F1bmw <- matrix(0,nn[2],1)
  F1gg <- (dataF1 - m[12])/sqrt(as.vector(sigma[11]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn1 < nn[2]){F1bmw <- F1bmw+runif(nn[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:nn[2])) - 0.5)/nn[2])^2)
  F1u<- as.matrix(c(12*nn[2]*((F1dd[1]/nn[2]-0.5)^2),((45*nn[2])/4)*((F1dd[2]/nn[2]-1/3)^2),180*nn[2]*((F1dd[3]/nn[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,nn[2]))))
  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)
  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*nn[3])
  P2bmw <- matrix(0,nn[3],1)
  P2gg <- (dataP2 - m[13])/sqrt(as.vector(sigma[11]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn1 < nn[3]){P2bmw <- P2bmw+runif(nn[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:nn[3])) - 0.5)/nn[3])^2)
  P2u<- as.matrix(c(12*nn[3]*((P2dd[1]/nn[3]-0.5)^2),((45*nn[3])/4)*((P2dd[2]/nn[3]-1/3)^2),180*nn[3]*((P2dd[3]/nn[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,nn[3]))))
  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1); 
  B1w1<-1/(12*nn[4])
  B1bmw <- matrix(0,nn[4],1); B1bmwsl <- matrix(0,nn[4],d4)
  for(i in 1:d4){
    B1gg <- (dataB1 - m[i])/sqrt(sigma[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B1P2)))[1]
  if(nn1 < nn[4]){B1P2 <- B1P2+runif(nn[4])/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*nn[4]) + sum((B1P2 - (as.matrix(c(1:nn[4])) - 0.5)/nn[4])^2)
  B1u <- as.matrix(c(12*nn[4]*((B1dd[1]/nn[4]-0.5)^2),((45*nn[4])/4)*((B1dd[2]/nn[4]-1/3)^2),180*nn[4]*((B1dd[3]/nn[4]-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),(1-pkolm(B1D,nn[4]))))
  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2); 
  B2w1<-1/(12*nn[5])
  B2bmw <- matrix(0,nn[5],1); B2bmwsl <- matrix(0,nn[5],d4)
  for(i in 1:d4){
    B2gg <- (dataB2 - mm2[i])/sqrt(ssigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B2P2)))[1]
  if(nn1 < nn[5]){B2P2 <- B2P2+runif(nn[5])/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*nn[5]) + sum((B2P2 - (as.matrix(c(1:nn[5])) - 0.5)/nn[5])^2)
  B2u <- as.matrix(c(12*nn[5]*((B2dd[1]/nn[5]-0.5)^2),((45*nn[5])/4)*((B2dd[2]/nn[5]-1/3)^2),180*nn[5]*((B2dd[3]/nn[5]-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),(1-pkolm(B2D,nn[5]))))
  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("MX2-CD-AD",round(abc,4),round(AIC,4),round(m[11],4),round(m[12],4),round(m[13],4),round(m[1],4),round(m[2],4),round(m[3],4),round(m[4],4),round(sigma[1],4),round(sigma[2],4),
                       round(sigma[3],4),round(sigma[4],4),round(mix_pi[1],4),round(mix_pi[2],4),round(mix_pi[3],4),round(mix_pi[4],4),round(m[5],4),round(m[6],4),round(m[7],4),round(m[8],4),round(sigma[5],4),round(sigma[6],4), round(sigma[7],4),round(sigma[8],4),      
                       round(mix_pi[5],4),round(mix_pi[6],4),round(mix_pi[7],4),round(mix_pi[8],4),round(sigma[11],4),round(B4[1],4),round(B4[2],4),round(B4[3],4)," "," ",round(B4[4],4),round(B4[5],4),round(B1jj,4),round(B1ll*100,4),round(B1_gg,4),round(B1rr*100,4),round(B2jj,4),round(B2ll*100,4),round(B2_gg,4),round(B2rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,as.matrix(mmi1),as.matrix(mmi2))
  return(OUTPUT)
}

#################MX2-EAD-AD(E-6)############################
G5BCFModelFun[[18]] <- function(K1,logL,df11,df21,df31,df41,df51,G5BCFtext2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataB1 <- as.matrix(as.numeric(df41[,1]))
  dataB2 <- as.matrix(as.numeric(df51[,1]))
  nn<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataB1)[1],dim(dataB2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataB1),sum(dataB2)))
  s<-as.matrix(c(sum(dataP1^2),sum(dataF1^2),sum(dataP2^2),sum(dataB1^2),sum(dataB2^2)))
  ss<-s[c(1:3)]-sumx[c(1:3)]^2/nn[c(1:3)]
  m<-matrix(0,13,1);m[11]<-mean(dataP1);m[12]<-mean(dataF1);m[13]<-mean(dataP2);man0<-mean(dataB1);man1<-mean(dataB2)
  sigma0<-sum(ss)/(nn[1]+nn[2]+nn[3]-3)
  sigma1<-var(dataB1);sigma2<-var(dataB2)
  d1<-1;d2<-2;d3<-3;d4<-4
  m_esp <- 0.0001;num_l <- as.numeric(G5BCFtext2)
  ###############procedure start###########################
  mi<-matrix(0,7,1);mi[c(1,3)]<-0.25;mi[c(2,6)]<-0.5;mi[c(5,7)]<-0.25
  a1<-sqrt(sigma1/(nn[4]-1))
  m[1]<-man0+2*a1;m[2]<-man0;m[3]<-man0-2*a1
  a1<-sqrt(sigma2/(nn[5]-1))
  m[5]<-man1+2*a1;m[6]<-man1;m[7]<-man1-2*a1
  sigma<-matrix(0,11,1)
  sigma[11]<-sigma0;sigma[1]<-sigma1/3;sigma[c(7,8)]<-sigma2/3
  hh1<-matrix(c(1,1,1,1,1,1,1,1,1,2,0,-2,2,1.5,1,1,-0.5,-0.5,
                1,0,-1,0.5,0.5,0.5,-0.5,-0.5,-0.5,0,1,0,0.25,0.25,0.25,0.25,0.25,0.25),9,4)
  mm<-m[c(11,12,13,1,2,3,5,6,7)]
  B1<-solve(crossprod(hh1,hh1))%*%crossprod(hh1,mm)
  a1<-B1[2]    # d.
  g<-matrix(0,2,1);g[1]<-0.75*a1^2/num_l;g[2]<-1.5*a1^2/num_l
  sigma[2]<-sigma[1]+g[1];sigma[3]<-sigma[1]+g[2];sigma[5]<-sigma[7]+g[2];sigma[6]<-sigma[7]+g[1]
  mm1<-m[c(1:3)];ssigma1<-sigma[c(1:3)];mmi1<-mi[c(1:3)]
  mm2<-m[c(5:7)];ssigma2<-sigma[c(5:7)];mmi2<-mi[c(5:7)]
  L0<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
    sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)))+
    sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)))
  ##########iteration process###########
  ############E-step###################
  iteration <- 0; stopa <- 1000
  W1 <- matrix(0,3,nn[4]); swx1 <- matrix(0,3,1)
  W2 <- matrix(0,3,nn[5]); swx2 <- matrix(0,3,1)
  sumwx<-matrix(0,8,1);mix_pi<-matrix(0,8,1)
  n0<-matrix(0,8,1);hh<-matrix(0,5,5);b_line<-matrix(0,5,1)
  swx<-matrix(0,8,1);s0<-matrix(0,6,1);aa3<-matrix(0,8,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:3) { W1[i,] <- mi[i]*dnorm(dataB1,m[i],sqrt(sigma[i]))/dmixnorm(dataB1,mm1,sqrt(ssigma1),mmi1)}
    mix_pi[c(1:3)] <- as.matrix(rowSums(W1)/nn[4])
    sumwx[c(1:3)] <- W1%*%dataB1
    for(i in 1:3) { W2[i,] <- mmi2[i]*dnorm(dataB2,mm2[i],sqrt(ssigma2[i]))/dmixnorm(dataB2,mm2,sqrt(ssigma2),mmi2)}
    mix_pi[c(5:7)] <- as.matrix(rowSums(W2)/nn[5])
    sumwx[c(5:7)] <- W2%*%dataB2
    n0[c(1:3)]<-nn[4]*mix_pi[c(1:3)];n0[c(5:7)]<-nn[5]*mix_pi[c(5:7)]
    n0[c(1,2,3,5,6,7)][abs(n0[c(1,2,3,5,6,7)])<0.000001]<-0.000001
    #########################################################
    hh[1,1]<-25*sigma[11]/nn[1]+4*sigma[11]/nn[2]+sigma[11]/nn[3]+64*sigma[2]/n0[2]
    hh[1,2]<-5*sigma[11]/nn[1]+4*sigma[11]/nn[2]+5*sigma[11]/nn[3]
    hh[1,3]<-16*sigma[2]/n0[2]
    hh[1,4]<-0
    hh[1,5]<-0
    hh[2,2]<-sigma[11]/nn[1]+4*sigma[11]/nn[2]+25*sigma[11]/nn[3]+64*sigma[6]/n0[6]
    hh[2,3]<-0
    hh[2,4]<-16*sigma[6]/n0[6]
    hh[2,5]<-0
    hh[3,3]<-sigma[1]/n0[1]+4*sigma[2]/n0[2]+sigma[3]/n0[3]
    hh[3,4]<-0
    hh[3,5]<-3*sigma[1]/n0[1]-3*sigma[3]/n0[3]
    hh[4,4]<-sigma[5]/n0[5]+4*sigma[6]/n0[6]+sigma[7]/n0[7]
    hh[4,5]<--sigma[5]/n0[5]+sigma[7]/n0[7]
    hh[5,5]<-9*sigma[1]/n0[1]+9*sigma[3]/n0[3]+sigma[5]/n0[5]+sigma[7]/n0[7]
    for(i in 2:5)
    {
      for(j in 1:(i-1))
      {
        hh[i,j]<-hh[j,i]
      }
    }
    ###########################################
    b_line[1]<-5*sumx[1]/nn[1]+2*sumx[2]/nn[2]+sumx[3]/nn[3]-8*sumwx[2]/n0[2]
    b_line[2]<-sumx[1]/nn[1]+2*sumx[2]/nn[2]+5*sumx[3]/nn[3]-8*sumwx[6]/n0[6]
    b_line[3]<-sumwx[1]/n0[1]-2*sumwx[2]/n0[2]+sumwx[3]/n0[3]
    b_line[4]<-sumwx[5]/n0[5]-2*sumwx[6]/n0[6]+sumwx[7]/n0[7]
    b_line[5]<-3*sumwx[1]/n0[1]-3*sumwx[3]/n0[3]-sumwx[5]/n0[5]+sumwx[7]/n0[7]
    B2<-solve(hh,b_line)
    #############################################
    m[11]<-(sumx[1]-sigma[11]*(5*B2[1]+B2[2]))/nn[1]
    m[12]<-(sumx[2]-sigma[11]*(2*B2[1]+B2[2]*2))/nn[2]
    m[13]<-(sumx[3]-sigma[11]*(B2[1]+5*B2[2]))/nn[3]
    m[1]<-(sumwx[1]-sigma[1]*(B2[3]+3*B2[5]))/n0[1]
    m[2]<-(sumwx[2]+sigma[2]*(8*B2[1]+2*B2[3]))/n0[2]
    m[3]<-(sumwx[3]+sigma[3]*(-B2[3]+3*B2[5]))/n0[3]
    m[5]<-(sumwx[5]+sigma[5]*(-B2[4]+B2[5]))/n0[5]
    m[6]<-(sumwx[6]+sigma[6]*(8*B2[2]+2*B2[4]))/n0[6]
    m[7]<-(sumwx[7]-sigma[7]*(B2[4]+B2[5]))/n0[7]
    ##############################################
    hh3<-hh1
    mm<-m[c(11,12,13,1,2,3,5,6,7)]
    B3<-solve(crossprod(hh3,hh3))%*%crossprod(hh3,mm)
    a1<-B3[2]    # d.
    g[1]<-0.75*a1^2/num_l;g[2]<-1.5*a1^2/num_l
    sigma[2]<-sigma[1]+g[1];sigma[3]<-sigma[1]+g[2]
    ########obtain variance#######################
    ss1<-sum((dataP1-m[11])^2);ss3<-sum((dataP2-m[13])^2);ss2<-sum((dataF1-m[12])^2)
    mm2<-m[c(5:7)]
    for(i in 1:d3) {swx1[i] <- W1[i,]%*%(dataB1-m[i])^2 } ;for(i in 1:d3) {swx2[i] <- W2[i,]%*%(dataB2-mm2[i])^2 }  
    swx[c(1:3)]<-swx1;swx[c(5:7)]<-swx2
    aaa0<-sigma[1];n_iter<-0;abc1<-1000
    while (abc1>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma[1]/(sigma[1]+g[1])
      aa2<-sigma[1]/(sigma[1]+g[2])
      sigma[1]<-(swx[1]+aa1^2*swx[2]+aa2^2*swx[3])/(n0[1]+aa1*n0[2]+aa2*n0[3])
      abc1<-abs(sigma[1]-aaa0)
      aaa0<-sigma[1]
      if (n_iter>20) break
    } 
    sigma40<-sigma[1]-sigma[11]
    if (sigma40<0) {sigma40<-0;sigma[1]<-sigma[11]}
    sigma[2]<-sigma[1]+g[1];sigma[3]<-sigma[1]+g[2]
    aaa0<-sigma[7];aa4<-1000;n_iter<-0
    while (aa4>0.0001)
    {
      n_iter<-n_iter+1
      aa1<-sigma[7]/(sigma[7]+g[1])
      aa2<-sigma[7]/(sigma[7]+g[2])
      sigma[7]<-(aa2^2*swx[5]+aa1^2*swx[6]+swx[7])/(aa2*n0[5]+aa1*n0[6]+n0[7])
      aa4<-abs(sigma[7]-aaa0)
      aaa0<-sigma[7]
      if (n_iter>20) break
    } 
    sigma50<-sigma[7]-sigma[11]
    if (sigma50<0) {sigma50<-0;sigma[7]<-sigma[11]}
    sigma[5]<-sigma[7]+g[2];sigma[6]<-sigma[7]+g[1]
    aa1<-sigma40;aa2<-sigma50
    if (aa1<0) {aa1<-0};if (aa2<0) {aa2<-0}
    s0[1]<-ss1+ss2+ss3;s0[2]<-nn[1]+nn[2]+nn[3]
    aaa0<-sigma[11];n_iter<-0;abc2<-1000
    while (abc2>0.0001)
    {
      n_iter<-n_iter+1
      aa3[1]<-sigma[11]/(sigma[11]+aa1)
      aa3[c(2:3)]<-sigma[11]/(sigma[11]+aa1+g[c(1:2)])
      s0[3]<-sum(aa3[c(1:3)]^2*swx[c(1:3)])
      s0[4]<-sum(aa3[c(1:3)]*n0[c(1:3)])
      aa3[7]<-sigma[11]/(sigma[11]+aa2)
      aa3[c(5:6)]<-sigma[11]/(sigma[11]+aa2+g[c(2:1)])
      s0[5]<-sum(aa3[c(5:7)]^2*swx[c(5:7)])
      s0[6]<-sum(aa3[c(5:7)]*n0[c(5:7)])
      sigma[11]<-(s0[1]+(s0[3]+s0[5]))/(s0[2]+s0[4]+s0[6])
      abc2<-abs(sigma[11]-aaa0)
      aaa0<-sigma[11]
      if (n_iter>20) break
    }
    sigma[1]<-sigma[11]+sigma40;sigma[c(2:3)]<-sigma[11]+sigma40+g[c(1:2)];sigma[c(5:6)]<-sigma[11]+sigma50+g[c(2:1)];sigma[7]<-sigma[11]+sigma50
    mm1<-m[c(1:3)];ssigma1<-sigma[c(1:3)]
    mm2<-m[c(5:7)];ssigma2<-sigma[c(5:7)]
    mix_pi1<-mix_pi[c(1:3)];mix_pi2<-mix_pi[c(5:7)]
    ########criteria for iterations to stop########
    L1<-sum(log(dnorm(dataP1,m[11],sqrt(sigma[11]))))+sum(log(dnorm(dataF1,m[12],sqrt(sigma[11]))))+
      sum(log(dnorm(dataP2,m[13],sqrt(sigma[11]))))+sum(log(dmixnorm(dataB1,mm1,sqrt(ssigma1),mix_pi1)))+
      sum(log(dmixnorm(dataB2,mm2,sqrt(ssigma2),mix_pi2)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
    if(iteration>300)break
  }
  abc<-L0
  AIC<--2*abc+2*7
  ########first order genetic parameters###############
  hh4<-hh1
  mm<-m[c(11,12,13,1,2,3,5,6,7)]
  B4<-solve(t(hh4)%*%hh4)%*%(t(hh4)%*%mm)
  ########second order genetic parameters###############   
  B1jj <- sigma1 - sigma[1]	   
  if(B1jj<0) {B1jj<-0}	   
  B1ll <- B1jj/sigma1
  B1_gg <- sigma[1] - sigma[11]
  if(B1_gg<0|| B1_gg>sigma1){B1_gg<-0}	   
  B1rr <- B1_gg/sigma1	   
  B2jj <- sigma2 - sigma[7]
  if(B2jj<0) {B2jj<-0}
  B2ll <- B2jj/sigma2
  B2_gg <- sigma[7] - sigma[11]
  if(B2_gg<0 || B2_gg>sigma2) {B2_gg<-0}	   
  B2rr <- B2_gg/sigma2
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*nn[1])
  P1bmw <- matrix(0,nn[1],1)
  P1gg <- (dataP1 - m[11])/sqrt(as.vector(sigma[11]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn1 < nn[1]){P1bmw <- P1bmw+runif(nn[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:nn[1])) - 0.5)/nn[1])^2)
  P1u<- as.matrix(c(12*nn[1]*((P1dd[1]/nn[1]-0.5)^2),((45*nn[1])/4)*((P1dd[2]/nn[1]-1/3)^2),180*nn[1]*((P1dd[3]/nn[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,nn[1]))))
  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*nn[2])
  F1bmw <- matrix(0,nn[2],1)
  F1gg <- (dataF1 - m[12])/sqrt(as.vector(sigma[11]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn1 < nn[2]){F1bmw <- F1bmw+runif(nn[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:nn[2])) - 0.5)/nn[2])^2)
  F1u<- as.matrix(c(12*nn[2]*((F1dd[1]/nn[2]-0.5)^2),((45*nn[2])/4)*((F1dd[2]/nn[2]-1/3)^2),180*nn[2]*((F1dd[3]/nn[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,nn[2]))))
  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)
  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*nn[3])
  P2bmw <- matrix(0,nn[3],1)
  P2gg <- (dataP2 - m[13])/sqrt(as.vector(sigma[11]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn1 < nn[3]){P2bmw <- P2bmw+runif(nn[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:nn[3])) - 0.5)/nn[3])^2)
  P2u<- as.matrix(c(12*nn[3]*((P2dd[1]/nn[3]-0.5)^2),((45*nn[3])/4)*((P2dd[2]/nn[3]-1/3)^2),180*nn[3]*((P2dd[3]/nn[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,nn[3]))))
  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B1#########################################
  dataB1 <- sort(dataB1); 
  B1w1<-1/(12*nn[4])
  B1bmw <- matrix(0,nn[4],1); B1bmwsl <- matrix(0,nn[4],d3)
  for(i in 1:d3){
    B1gg <- (dataB1 - m[i])/sqrt(sigma[i])
    B1bmw[which(B1gg>=0)] <- pnorm(B1gg[B1gg>=0])
    B1bmw[which(B1gg<0)] <- 1 - pnorm(abs(B1gg[B1gg<0]))
    B1bmwsl[,i] <- B1bmw*mix_pi[i]
  }
  B1P2 <- rowSums(B1bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B1P2)))[1]
  if(nn1 < nn[4]){B1P2 <- B1P2+runif(nn[4])/1e4}
  ##########################################################
  B1dd <- as.matrix(c(sum(B1P2),sum(B1P2^2),sum((B1P2-0.5)^2)))
  B1WW2 <- 1/(12*nn[4]) + sum((B1P2 - (as.matrix(c(1:nn[4])) - 0.5)/nn[4])^2)
  B1u <- as.matrix(c(12*nn[4]*((B1dd[1]/nn[4]-0.5)^2),((45*nn[4])/4)*((B1dd[2]/nn[4]-1/3)^2),180*nn[4]*((B1dd[3]/nn[4]-1/12)^2)))
  B1D <- as.numeric(ks.test(B1P2,"punif")[[1]][1])
  B1tt <- as.matrix(c((1 - pchisq(B1u[1],1)),(1 - pchisq(B1u[2],1)),(1 - pchisq(B1u[3],1)),K1(B1WW2),(1-pkolm(B1D,nn[4]))))
  B1tt[which( B1tt>=10e-4)]<-round(B1tt[which(B1tt>=10e-4)],4);B1tt[which(B1tt<10e-4)]<-format(B1tt[which(B1tt<10e-4)],scientific=TRUE,digit=4)
  ####################################hypothesis testing for B2#########################################
  dataB2 <- sort(dataB2); 
  B2w1<-1/(12*nn[5])
  B2bmw <- matrix(0,nn[5],1); B2bmwsl <- matrix(0,nn[5],d3)
  for(i in 1:d3){
    B2gg <- (dataB2 - mm2[i])/sqrt(ssigma2[i])
    B2bmw[which(B2gg>=0)] <- pnorm(B2gg[B2gg>=0])
    B2bmw[which(B2gg<0)] <- 1 - pnorm(abs(B2gg[B2gg<0]))
    B2bmwsl[,i] <- B2bmw*mix_pi2[i]
  }
  B2P2 <- rowSums(B2bmwsl)
  #############deal with ties ####################
  nn1 <- dim(as.matrix(unique(B2P2)))[1]
  if(nn1 < nn[5]){B2P2 <- B2P2+runif(nn[5])/1e4}
  ##########################################################
  B2dd <- as.matrix(c(sum(B2P2),sum(B2P2^2),sum((B2P2-0.5)^2)))
  B2WW2 <- 1/(12*nn[5]) + sum((B2P2 - (as.matrix(c(1:nn[5])) - 0.5)/nn[5])^2)
  B2u <- as.matrix(c(12*nn[5]*((B2dd[1]/nn[5]-0.5)^2),((45*nn[5])/4)*((B2dd[2]/nn[5]-1/3)^2),180*nn[5]*((B2dd[3]/nn[5]-1/12)^2)))
  B2D <- as.numeric(ks.test(B2P2,"punif")[[1]][1])
  B2tt <- as.matrix(c((1 - pchisq(B2u[1],1)),(1 - pchisq(B2u[2],1)),(1 - pchisq(B2u[3],1)),K1(B2WW2),(1-pkolm(B2D,nn[5]))))
  B2tt[which( B2tt>=10e-4)]<-round(B2tt[which(B2tt>=10e-4)],4);B2tt[which(B2tt<10e-4)]<-format(B2tt[which(B2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("MX2-EAD-AD",round(abc,4),round(AIC,4),round(m[11],4),round(m[12],4),round(m[13],4),round(m[1],4),round(m[2],4),round(m[3],4)," ",round(sigma[1],4),round(sigma[2],4),
                       round(sigma[3],4), " ",round(mix_pi[1],4),round(mix_pi[2],4),round(mix_pi[3],4)," ",round(m[5],4),round(m[6],4),round(m[7],4)," ",round(sigma[5],4),round(sigma[6],4), round(sigma[7],4)," ",   
                       round(mix_pi[5],4),round(mix_pi[6],4),round(mix_pi[7],4)," ",round(sigma[11],4),round(B4[1],4),round(B4[2],4)," "," "," ",round(B4[3],4),round(B4[4],4),round(B1jj,4),round(B1ll*100,4),round(B1_gg,4),round(B1rr*100,4),round(B2jj,4),round(B2ll*100,4),round(B2_gg,4),round(B2rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(B1u[1],4),B1tt[1],round(B1u[2],4),B1tt[2],round(B1u[3],4),B1tt[3],round(B1WW2,4),B1tt[4],round(B1D,4),B1tt[5],
                       round(B2u[1],4),B2tt[1],round(B2u[2],4),B2tt[2],round(B2u[3],4),B2tt[3],round(B2WW2,4),B2tt[4],round(B2D,4),B2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,as.matrix(mmi1),as.matrix(mmi2))
  return(OUTPUT)
}



K1G5BCF <- function(x){
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

logLG5BCF <- function(nm,nng,mi,mn,s,d1) { sum2 <- sum(log(dmixnorm(d1,mn,sqrt(s),mi)));return (sum2) } 


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
  requireNamespace("kolmim")
  G5BCFModelFun[[i]](K1G5BCF,logLG5BCF,df11,df21,df31,df41,df51,G5BCFtext2)[[1]]
}
stopCluster(cl)

mmi1<-NULL;mmi2<-NULL

}else{
  
allresultq=switch(model,"1MG-AD" = G5BCFModelFun[[1]](K1G5BCF,logLG5BCF,df11,df21,df31,df41,df51,G5BCFtext2),"1MG-A"=G5BCFModelFun[[2]](K1G5BCF,logLG5BCF,df11,df21,df31,df41,df51,G5BCFtext2),"1MG-EAD"=G5BCFModelFun[[3]](K1G5BCF,logLG5BCF,df11,df21,df31,df41,df51,G5BCFtext2),"1MG-NCD"=G5BCFModelFun[[4]](K1G5BCF,logLG5BCF,df11,df21,df31,df41,df51,G5BCFtext2),
         "2MG-AD"=G5BCFModelFun[[5]](K1G5BCF,logLG5BCF,df11,df21,df31,df41,df51,G5BCFtext2),"2MG-A"=G5BCFModelFun[[6]](K1G5BCF,logLG5BCF,df11,df21,df31,df41,df51,G5BCFtext2),"2MG-EA"=G5BCFModelFun[[7]](K1G5BCF,logLG5BCF,df11,df21,df31,df41,df51,G5BCFtext2),"2MG-CD"=G5BCFModelFun[[8]](K1G5BCF,logLG5BCF,df11,df21,df31,df41,df51,G5BCFtext2),
         "2MG-EAD"=G5BCFModelFun[[9]](K1G5BCF,logLG5BCF,df11,df21,df31,df41,df51,G5BCFtext2),"MX1-AD-AD"=G5BCFModelFun[[10]](K1G5BCF,logLG5BCF,df11,df21,df31,df41,df51,G5BCFtext2),"MX1-A-AD"=G5BCFModelFun[[11]](K1G5BCF,logLG5BCF,df11,df21,df31,df41,df51,G5BCFtext2),"MX1-EAD-AD"=G5BCFModelFun[[12]](K1G5BCF,logLG5BCF,df11,df21,df31,df41,df51,G5BCFtext2),
          "MX1-NCD-AD"=G5BCFModelFun[[13]](K1G5BCF,logLG5BCF,df11,df21,df31,df41,df51,G5BCFtext2),"MX2-AD-AD"=G5BCFModelFun[[14]](K1G5BCF,logLG5BCF,df11,df21,df31,df41,df51,G5BCFtext2),"MX2-A-AD"=G5BCFModelFun[[15]](K1G5BCF,logLG5BCF,df11,df21,df31,df41,df51,G5BCFtext2),"MX2-EA-AD"=G5BCFModelFun[[16]](K1G5BCF,logLG5BCF,df11,df21,df31,df41,df51,G5BCFtext2),
         "MX2-CD-AD"=G5BCFModelFun[[17]](K1G5BCF,logLG5BCF,df11,df21,df31,df41,df51,G5BCFtext2),"MX2-EAD-AD"=G5BCFModelFun[[18]](K1G5BCF,logLG5BCF,df11,df21,df31,df41,df51,G5BCFtext2)) 
  
  
allresult<-allresultq[[1]]
if(model!="0MG"){
  mmi1<-allresultq[[2]];mmi2<-allresultq[[3]] 
}else{
  mmi1<-NULL;mmi2<-NULL
}
}
colnames(allresult) <- G5BCFcolname
out<-list(allresult,mmi1,mmi2)
return(out)
} 




















