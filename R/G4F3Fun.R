G4F3Fun<-function(df,model,G4F3text2){

data<-sapply(df,as.character)

dP1<-data[-1,which(data[1,]=="P1")];P1<-as.numeric(dP1[which(is.na(as.numeric(dP1))==FALSE)]);df11<-as.data.frame(P1)
dF1<-data[-1,which(data[1,]=="F1")];F1<-as.numeric(dF1[which(is.na(as.numeric(dF1))==FALSE)]);df21<-as.data.frame(F1)
dP2<-data[-1,which(data[1,]=="P2")];P2<-as.numeric(dP2[which(is.na(as.numeric(dP2))==FALSE)]);df31<-as.data.frame(P2)
dF23<-data[-1,which(data[1,]=="F23")];F23<-as.numeric(dF23[which(is.na(as.numeric(dF23))==FALSE)]);df41<-as.data.frame(F23)

##################################################
G4F3colname <- c("Model","Log_Max_likelihood_Value","AIC","meanP1","meanF1","meanP2","mean[1]","mean[2]","mean[3]","mean[4]","mean[5]","mean[6]","mean[7]","mean[8]","mean[9]",
                 "Var(Residual)","Var[1]","Var[2]","Var[3]","Var[4]","Var[5]","Var[6]","Var[7]","Var[8]","Var[9]","Proportion[1]","Proportion[2]","Proportion[3]","Proportion[4]","Proportion[5]","Proportion[6]","Proportion[7]","Proportion[8]","Proportion[9]",
                 "m1(m)","m2","m3","m4","da(d)","db","ha(h)","hb","i","jab","jba","l","[d]","[h]","Major-Gene Var","Heritability(Major-Gene)(%)","Poly-Gene Var","Heritability(Poly-Gene)(%)","U1 square(P1)","P(U1 square(P1))","U2 square(P1)","P(U2 square(P1))","U3 square(P1)","P(U3 square(P1))","nW square(P1)","P(nW square(P1))","Dn(P1)","P(Dn(P1))",
                 "U1 square(F1)","P(U1 square(F1))","U2 square(F1)","P(U2 square(F1))","U3 square(F1)","P(U3 square(F1))","nW square(F1)","P(nW square(F1))","Dn(F1)","P(Dn(F1))","U1 square(P2)","P(U1 square(P2))","U2 square(P2)","P(U2 square(P2))","U3 square(P2)","P(U3 square(P2))","nW square(P2)","P(nW square(P2))","Dn(P2)","P(Dn(P2))",
                 "U1 square(F3)","P(U1 square(F3))","U2 square(F3)","P(U2 square(F3))","U3 square(F3)","P(U3 square(F3))","nW square(F3)","P(nW square(F3))","Dn(F3)","P(Dn(F3))")

G4F3ModelFun<-list(NA)
###################define each model function############################
##########################1MG_AD(A_1)#################################
G4F3ModelFun[[1]] <- function(K1,logL,df11,df21,df31,df41,G4F3text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF3 <- as.matrix(as.numeric(df41[,1]))

  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF3)))
  sigmaf3<-var(dataF3);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a8<-sqrt(sigmaf3/(m_sam[4]))
  m_nf <- as.numeric(G4F3text2)
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF3)))
  ###############procedure start##############################
  d1<-3
  mi<-as.matrix(c(0.25,0.5,0.25))
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaf3/2,sigmaf3/2,sigmaf3/2))
  m[c(1:4)]<-m
  m[5]<-m[4] ;m[4]<-m[5]+2.0*a8;m[6]<-m[5]-2.0*a8
  pi<-m[c(4:6)];gh<-sigma[c(4:6)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mi)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,3,m_sam[4]); swx <- matrix(0,3,1);g_a<-matrix(0,1,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ########################E-step####################
    mm<-m[c(4,5,6)];ssigma<-sigma[c(4,5,6)]
    for(i in 1:d1) { WW[i,] <- mi[i]*dnorm(dataF3,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF3,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF3
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    #####################restrictions################################
    aa1<-(sumx[1]+m_nf*sumwx[1])/(m_sam[1]+m_nf*n0[1])+2*sumx[2]/m_sam[2]+(sumx[3]+ m_nf*sumwx[3])/(m_sam[3]+m_nf*n0[3])-4*sumwx[2]/n0[2]
    aa2<-sigma[1]/(m_sam[1]+m_nf*n0[1])+4*sigma[2]/m_sam[2]+sigma[3]/(m_sam[3]+m_nf*n0[3])+16*sigma[5]/n0[2]
    aa3<-aa1/aa2
    m[1]<-(sumx[1]+m_nf*sumwx[1]-sigma[1]*aa3)/(m_sam[1]+m_nf*n0[1])
    m[2]<-(sumx[2]-sigma[2]*aa3*2)/m_sam[2]
    m[3]<-(sumx[3]+m_nf*sumwx[3]-sigma[3]*aa3)/(m_sam[3]+m_nf*n0[3])
    m[5]<-(sumwx[2]+sigma[5]*aa3*4)/n0[2]
    m[4]<-m[1]; m[6]<-m[3]
    #########first order genetic parameter process##########
    hh1 <- matrix(c(1,1,1,1,1,0,-1,0,0,1,0,0.5),4,3)
    B1<-solve(crossprod(hh1,hh1))%*%crossprod(hh1,m[c(1,2,3,5)])
    g_a[1]<-(0.5*B1[2]^2+0.25*B1[3]^2)/m_nf
    ########################obtain variance##################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d1) {swx[i] <- WW[i,]%*%(dataF3-m[i+3])^2 }
    a1<-sigma[1];a2<-ss1+ss2+ss3;a3<-m_sam[1]+m_sam[2]+m_sam[3]
    n_iter<-0
    a6<-1000
    while(a6>0.0001){
      n_iter<-n_iter+1
      a4<-(sigma[1]/m_nf)/(sigma[1]/m_nf+g_a[1])
      a4[a4>1]<-1
      sigma[1]<-(a2+m_nf*(swx[1]+swx[3])+a4^2*m_nf*swx[2])/(a3+n0[1]+n0[3]+a4*n0[2])
      a6<-abs(sigma[1]-a1)
      a1<-sigma[1]
      if(n_iter>20)break
    }
    sigma[2]<-sigma[3]<-sigma[1];sigma[4]<-sigma[1]/m_nf
    sigma[6]<-sigma[4];sigma[5]<-sigma[4]+g_a[1]
    ########criteria for iterations to stop####################
    pi<-m[c(4:6)];gh<-sigma[c(4:6)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*4
  #################second oder genetic parameter process#########################################
  jj<-sigmaf3 - sigma[1]/m_nf
  if(jj<0){jj<-0}
  ll<-jj/sigmaf3
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma[1]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < m_sam[1]){P1bmw <- P1bmw+runif(m_sam[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:m_sam[1])) - 0.5)/m_sam[1])^2)
  P1u<- as.matrix(c(12*m_sam[1]*((P1dd[1]/m_sam[1]-0.5)^2),((45*m_sam[1])/4)*((P1dd[2]/m_sam[1]-1/3)^2),180*m_sam[1]*((P1dd[3]/m_sam[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,m_sam[1]))))

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*m_sam[2])
  F1bmw <- matrix(0,m_sam[2],1)
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma[2]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < m_sam[2]){F1bmw <- F1bmw+runif(m_sam[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:m_sam[2])) - 0.5)/m_sam[2])^2)
  F1u<- as.matrix(c(12*m_sam[2]*((F1dd[1]/m_sam[2]-0.5)^2),((45*m_sam[2])/4)*((F1dd[2]/m_sam[2]-1/3)^2),180*m_sam[2]*((F1dd[3]/m_sam[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,m_sam[2]))))

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*m_sam[3])
  P2bmw <- matrix(0,m_sam[3],1)
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma[3]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < m_sam[3]){P2bmw <- P2bmw+runif(m_sam[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:m_sam[3])) - 0.5)/m_sam[3])^2)
  P2u<- as.matrix(c(12*m_sam[3]*((P2dd[1]/m_sam[3]-0.5)^2),((45*m_sam[3])/4)*((P2dd[2]/m_sam[3]-1/3)^2),180*m_sam[3]*((P2dd[3]/m_sam[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,m_sam[3]))))

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3);
  F3w1<-1/(12*m_sam[4])
  F3bmw <- matrix(0,m_sam[4],1); F3bmwsl <- matrix(0,m_sam[4],d1)
  for(i in 1:d1){
    F3gg <- (dataF3 - pi[i])/sqrt(gh[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mix_pi[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[4]){F3P2 <- F3P2+runif(m_sam[4])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[4]) + sum((F3P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F3u <- as.matrix(c(12*m_sam[4]*((F3dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F3dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F3dd[3]/m_sam[4]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[4]))))

  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("1MG-AD",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," "," ",round(sigma[1],4),round(sigma[4],4),round(sigma[5],4),round(sigma[6],4)," "," "," "," "," "," ",
                       round(t(mix_pi),4)," "," "," "," "," "," ",round(B1[1],4)," "," "," ",round(B1[2],4)," ",round(B1[3],4)," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
####################################1MG_A(A_2)###############################################
G4F3ModelFun[[2]] <- function(K1,logL,df11,df21,df31,df41,G4F3text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF3 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF3)))
  sigmaf3<-var(dataF3);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a8<-sqrt(sigmaf3/(m_sam[4]))
  m_nf <- as.numeric(G4F3text2)
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF3)))
  ############procedure start#############
  mi<-as.matrix(c(0.25,0.5,0.25))
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaf3/2,1.5*sigmaf3/2,sigmaf3/2))
  m[c(1:4)]<-m
  m[5]<-m[4];m[4]<-m[5]+2.0*a8;m[6]<-m[5]-2.0*a8
  pi<-m[c(4:6)];gh<-sigma[c(4:6)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mi)))
  ##############iteration process###############
  d1<-3
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,3,m_sam[4]); swx <- matrix(0,3,1);g_a<-matrix(0,1,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4,5,6)];ssigma<-sigma[c(4,5,6)]
    for(i in 1:d1) { WW[i,] <- mi[i]*dnorm(dataF3,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF3,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF3
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    ##############restriction###################
    aa1<-(sumx[1]+m_nf*sumwx[1])/(m_sam[1]+m_nf*n0[1])-2*(sumx[2]+m_nf*sumwx[2])/(m_sam[2]+m_nf*n0[2])+(sumx[3]+m_nf*sumwx[3])/(m_sam[3]+m_nf*n0[3])
    aa2<-sigma[1]/(m_sam[1]+m_nf*n0[1])+4*sigma[2]/(m_sam[2]+m_nf*n0[2])+sigma[3]/(m_sam[3]+m_nf*n0[3])
    rr1<-aa1/aa2
    m[1]<-(sumx[1]+m_nf*sumwx[1]-sigma[1]*rr1)/(m_sam[1]+m_nf*n0[1])
    m[2]<-(sumx[2]+m_nf*sumwx[2]+sigma[1]*rr1*2)/(m_sam[2]+m_nf*n0[2])
    m[3]<-(sumx[3]+m_nf*sumwx[3]-sigma[1]*rr1)/(m_sam[3]+m_nf*n0[3])
    m[c(4:6)]<-m[c(1:3)]
    #########first order genetic parameter process#######
    hh1 <- matrix(c(1,1,1,1,0,-1),3,2)
    B2<-solve(crossprod(hh1,hh1))%*%crossprod(hh1,m[c(1,2,3)])
    g_a[1]<-(0.5*B2[2]^2)/m_nf
    #####################obtain variance##################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d1) {swx[i] <- WW[i,]%*%(dataF3-m[i+3])^2 }
    a1<-sigma[1];a2<-ss1+ss2+ss3;a3<-m_sam[1]+m_sam[2]+m_sam[3]
    n_iter<-0
    a6<-1000
    while(a6>0.0001){
      n_iter<-n_iter+1
      a4<-(sigma[1]/m_nf)/(sigma[1]/m_nf+g_a[1]);
      a4[a4>1]<-1
      sigma[1]<-(a2+m_nf*(swx[1]+swx[3])+a4^2*m_nf*swx[2])/(a3+n0[1]+n0[3]+a4*n0[2])
      a6<-abs(sigma[1]-a1)
      a1<-sigma[1]
      if(n_iter>20)break
    }
    sigma[2]<-sigma[3]<-sigma[1];sigma[4]<-sigma[1]/m_nf
    sigma[6]<-sigma[4];sigma[5]<-sigma[4]+g_a[1]
    ########criteria for iterations to stop####################
    pi<-m[c(4:6)];gh<-sigma[c(4:6)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*2
  ##############second oder genetic parameter process############################################
  jj<-sigmaf3-sigma[1]/m_nf
  if(jj<0){jj<-0}
  ll<-jj/sigmaf3
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma[1]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < m_sam[1]){P1bmw <- P1bmw+runif(m_sam[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:m_sam[1])) - 0.5)/m_sam[1])^2)
  P1u<- as.matrix(c(12*m_sam[1]*((P1dd[1]/m_sam[1]-0.5)^2),((45*m_sam[1])/4)*((P1dd[2]/m_sam[1]-1/3)^2),180*m_sam[1]*((P1dd[3]/m_sam[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,m_sam[1]))))

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*m_sam[2])
  F1bmw <- matrix(0,m_sam[2],1)
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma[2]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < m_sam[2]){F1bmw <- F1bmw+runif(m_sam[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:m_sam[2])) - 0.5)/m_sam[2])^2)
  F1u<- as.matrix(c(12*m_sam[2]*((F1dd[1]/m_sam[2]-0.5)^2),((45*m_sam[2])/4)*((F1dd[2]/m_sam[2]-1/3)^2),180*m_sam[2]*((F1dd[3]/m_sam[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,m_sam[2]))))

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*m_sam[3])
  P2bmw <- matrix(0,m_sam[3],1)
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma[3]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < m_sam[3]){P2bmw <- P2bmw+runif(m_sam[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:m_sam[3])) - 0.5)/m_sam[3])^2)
  P2u<- as.matrix(c(12*m_sam[3]*((P2dd[1]/m_sam[3]-0.5)^2),((45*m_sam[3])/4)*((P2dd[2]/m_sam[3]-1/3)^2),180*m_sam[3]*((P2dd[3]/m_sam[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,m_sam[3]))))

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3);
  F3w1<-1/(12*m_sam[4])
  F3bmw <- matrix(0,m_sam[4],1); F3bmwsl <- matrix(0,m_sam[4],d1)
  for(i in 1:d1){
    F3gg <- (dataF3 - pi[i])/sqrt(gh[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mix_pi[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[4]){F3P2 <- F3P2+runif(m_sam[4])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[4]) + sum((F3P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F3u <- as.matrix(c(12*m_sam[4]*((F3dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F3dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F3dd[3]/m_sam[4]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[4]))))

  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("1MG-A",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," "," ",round(sigma[1],4),round(sigma[4],4),round(sigma[5],4),round(sigma[6],4)," "," "," "," "," "," ",
                       round(t(mix_pi),4)," "," "," "," "," "," ",round(B2[1],4)," "," "," ",round(B2[2],4)," "," "," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
####################################1MG_EAD(A_3)##################################################
G4F3ModelFun[[3]] <- function(K1,logL,df11,df21,df31,df41,G4F3text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF3 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF3)))
  sigmaf3<-var(dataF3);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a8<-sqrt(sigmaf3/(m_sam[4]))
  m_nf <- as.numeric(G4F3text2)
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF3)))
  #################procedure start######################################
  mi<-as.matrix(c(0.25,0.5,0.25))
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaf3/2,1.5*sigmaf3/2,sigmaf3/2))
  m[c(1:4)]<-m
  m[5]<-m[4]+1.5*a8; m[6]<-m[4]-2*a8;m[4]<-m[4]+2*a8
  pi<-m[c(4:6)]; gh<-sigma[c(4:6)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mi)))
  ##############iteration process#####################
  d1<-3
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,3,m_sam[4]); swx <- matrix(0,3,1);g_a<-matrix(0,1,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4,5,6)];ssigma<-sigma[c(4,5,6)]
    for(i in 1:d1) { WW[i,] <- mi[i]*dnorm(dataF3,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF3,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF3
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    #####################restrictions########################
    aa1<-3*(sumx[1]+sumx[2]+m_nf*sumwx[1])/(m_sam[1]+m_sam[2]+m_nf*n0[1])+(sumx[3]+m_nf*sumwx[3])/(m_sam[3]+m_nf*n0[3])-4*sumwx[2]/n0[2]
    aa2<-9*sigma[1]/(m_sam[1]+m_sam[2]+m_nf*n0[1])+sigma[3]/(m_sam[3]+m_nf*n0[3])+16*sigma[5]/n0[2]
    rr1<-aa1/aa2
    m[1]<-(sumx[1]+sumx[2]+m_nf*sumwx[1]-sigma[1]*rr1*3)/(m_sam[1]+m_sam[2]+m_nf*n0[1])
    m[3]<-(sumx[3]+m_nf*sumwx[3]-sigma[3]*rr1)/(m_sam[3]+m_nf*n0[3])
    m[5]<-(sumwx[2]+sigma[5]*4*rr1)/n0[2]
    m[2]<-m[1];m[4]<-m[1];m[6]<-m[3]
    #########first order genetic parameter process##########
    hh1 <- matrix(c(1,1,1,1,-1,0.5),3,2)
    B3<-solve(crossprod(hh1,hh1))%*%crossprod(hh1,m[c(1,3,5)])
    g_a[1]<-(0.75*B3[2]^2)/m_nf
    ########################obtain variance##################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d1) {swx[i] <- WW[i,]%*%(dataF3-m[i+3])^2 }
    a1<-sigma[1];a2<-ss1+ss2+ss3;a3<-m_sam[1]+m_sam[2]+m_sam[3]
    n_iter<-0
    a6<-1000
    while(a6>0.0001){
      n_iter<-n_iter+1
      a4<-(sigma[1]/m_nf)/(sigma[1]/m_nf+g_a[1])
      a4[a4>1]<-1
      sigma[1]<-(a2+m_nf*(swx[1]+swx[3])+a4^2*m_nf*swx[2])/(a3+n0[1]+n0[3]+a4*n0[2])
      a6<-abs(sigma[1]-a1)
      a1<-sigma[1]
      if(n_iter>20)break
    }
    sigma[2]<-sigma[3]<-sigma[1];sigma[4]<-sigma[1]/m_nf
    sigma[6]<-sigma[4];sigma[5]<-sigma[4]+g_a[1]
    ################criteria for iterations to stop###################
    pi<-m[c(4:6)];gh<-sigma[c(4:6)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*2
  ##################second order genetic parameter pcocess########################################
  jj<-sigmaf3-sigma[1]/m_nf
  if(jj<0){jj<-0}
  ll<-jj/sigmaf3
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma[1]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < m_sam[1]){P1bmw <- P1bmw+runif(m_sam[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:m_sam[1])) - 0.5)/m_sam[1])^2)
  P1u<- as.matrix(c(12*m_sam[1]*((P1dd[1]/m_sam[1]-0.5)^2),((45*m_sam[1])/4)*((P1dd[2]/m_sam[1]-1/3)^2),180*m_sam[1]*((P1dd[3]/m_sam[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,m_sam[1]))))

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*m_sam[2])
  F1bmw <- matrix(0,m_sam[2],1)
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma[2]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < m_sam[2]){F1bmw <- F1bmw+runif(m_sam[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:m_sam[2])) - 0.5)/m_sam[2])^2)
  F1u<- as.matrix(c(12*m_sam[2]*((F1dd[1]/m_sam[2]-0.5)^2),((45*m_sam[2])/4)*((F1dd[2]/m_sam[2]-1/3)^2),180*m_sam[2]*((F1dd[3]/m_sam[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,m_sam[2]))))

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*m_sam[3])
  P2bmw <- matrix(0,m_sam[3],1)
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma[3]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < m_sam[3]){P2bmw <- P2bmw+runif(m_sam[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:m_sam[3])) - 0.5)/m_sam[3])^2)
  P2u<- as.matrix(c(12*m_sam[3]*((P2dd[1]/m_sam[3]-0.5)^2),((45*m_sam[3])/4)*((P2dd[2]/m_sam[3]-1/3)^2),180*m_sam[3]*((P2dd[3]/m_sam[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,m_sam[3]))))

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3);
  F3w1<-1/(12*m_sam[4])
  F3bmw <- matrix(0,m_sam[4],1); F3bmwsl <- matrix(0,m_sam[4],d1)
  for(i in 1:d1){
    F3gg <- (dataF3 - pi[i])/sqrt(gh[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mix_pi[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[4]){F3P2 <- F3P2+runif(m_sam[4])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[4]) + sum((F3P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F3u <- as.matrix(c(12*m_sam[4]*((F3dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F3dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F3dd[3]/m_sam[4]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[4]))))

  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("1MG-EAD",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," "," ",round(sigma[1],4),round(sigma[4],4),round(sigma[5],4),round(sigma[6],4)," "," "," "," "," "," ",
                       round(t(mix_pi),4)," "," "," "," "," "," ",round(B3[1],4)," "," "," ",round(B3[2],4)," "," "," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
######################################1MG_NCD(A_4)#####################################################
G4F3ModelFun[[4]] <- function(K1,logL,df11,df21,df31,df41,G4F3text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF3 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF3)))
  sigmaf3<-var(dataF3);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a8<-sqrt(sigmaf3/(m_sam[4]))
  m_nf <- as.numeric(G4F3text2)
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF3)))
  ####################procedure start###############################################
  mi<-as.matrix(c(0.25,0.5,0.25))
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaf3/2,1.5*sigmaf3/2,sigmaf3/2))
  m[c(1:4)]<-m
  m[5]<-m[4]-1.5*a8;m[6]<-m[4]-2*a8;m[4]<-m[4]+2*a8;
  pi<-m[c(4:6)];gh<-sigma[c(4:6)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mi)))
  ##########iteration process###########
  d1<-3
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,3,m_sam[4]); swx <- matrix(0,3,1);g_a<-matrix(0,1,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4,5,6)];ssigma<-sigma[c(4,5,6)]
    for(i in 1:d1) { WW[i,] <- mi[i]*dnorm(dataF3,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF3,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF3
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    #################restrictions###############################
    aa1<-(sumx[1]+m_nf*sumwx[1])/(m_sam[1]+m_nf*n0[1])+3*(sumx[2]+sumx[3]+m_nf*sumwx[3])/(m_sam[2]+m_sam[3]+m_nf*n0[3])-4*sumwx[2]/n0[2]
    aa2<-sigma[1]/(m_sam[1]+m_nf*n0[1])+9*sigma[3]/(m_sam[2]+m_sam[3]+m_nf*n0[3])+16*sigma[5]/n0[2];
    rr1<-aa1/aa2
    m[1]<-(sumx[1]+m_nf*sumwx[1]-sigma[1]*rr1)/(m_sam[1]+m_nf*n0[1])
    m[2]<-(sumx[2]+sumx[3]+m_nf*sumwx[3]-sigma[3]*rr1*3)/(m_sam[2]+m_sam[3]+m_nf*n0[3])
    m[5]<-(sumwx[2]+sigma[5]*4*rr1)/n0[2]
    m[3]<-m[2]; m[4]<-m[1]; m[6]<-m[3]
    #########first order genetic parameter process##########
    hh1 <- matrix(c(1,1,1,1,-1,-0.5),3,2)
    B4<-solve(crossprod(hh1,hh1))%*%crossprod(hh1,m[c(1,3,5)])
    g_a[1]<-(0.75*B4[2]^2)/m_nf
    ########################obtain variance##################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d1) {swx[i] <- WW[i,]%*%(dataF3-m[i+3])^2 }
    a1<-sigma[1];a2<-ss1+ss2+ss3;a3<-m_sam[1]+m_sam[2]+m_sam[3]
    n_iter<-0
    a6<-1000
    while(a6>0.0001){
      n_iter<-n_iter+1
      a4<-(sigma[1]/m_nf)/(sigma[1]/m_nf+g_a[1])
      a4[a4>1]<-1
      sigma[1]<-(a2+m_nf*(swx[1]+swx[3])+a4^2*m_nf*swx[2])/(a3+n0[1]+n0[3]+a4*n0[2])
      a6<-abs(sigma[1]-a1)
      a1<-sigma[1]
      if(n_iter>20)break
    }
    sigma[2]<-sigma[3]<-sigma[1];sigma[4]<-sigma[1]/m_nf
    sigma[6]<-sigma[4];sigma[5]<-sigma[4]+g_a[1]
    ##############criteria for iterations to stop##############################
    pi<-m[c(4:6)];gh<-sigma[c(4:6)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*2
  ##############second order genetic parameter process#################################
  jj<-sigmaf3-sigma[1]/m_nf
  if(jj<0){jj<-0}
  ll<-jj/sigmaf3
  ##################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma[1]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < m_sam[1]){P1bmw <- P1bmw+runif(m_sam[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:m_sam[1])) - 0.5)/m_sam[1])^2)
  P1u<- as.matrix(c(12*m_sam[1]*((P1dd[1]/m_sam[1]-0.5)^2),((45*m_sam[1])/4)*((P1dd[2]/m_sam[1]-1/3)^2),180*m_sam[1]*((P1dd[3]/m_sam[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,m_sam[1]))))

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*m_sam[2])
  F1bmw <- matrix(0,m_sam[2],1)
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma[2]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < m_sam[2]){F1bmw <- F1bmw+runif(m_sam[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:m_sam[2])) - 0.5)/m_sam[2])^2)
  F1u<- as.matrix(c(12*m_sam[2]*((F1dd[1]/m_sam[2]-0.5)^2),((45*m_sam[2])/4)*((F1dd[2]/m_sam[2]-1/3)^2),180*m_sam[2]*((F1dd[3]/m_sam[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,m_sam[2]))))

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*m_sam[3])
  P2bmw <- matrix(0,m_sam[3],1)
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma[3]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < m_sam[3]){P2bmw <- P2bmw+runif(m_sam[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:m_sam[3])) - 0.5)/m_sam[3])^2)
  P2u<- as.matrix(c(12*m_sam[3]*((P2dd[1]/m_sam[3]-0.5)^2),((45*m_sam[3])/4)*((P2dd[2]/m_sam[3]-1/3)^2),180*m_sam[3]*((P2dd[3]/m_sam[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,m_sam[3]))))

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3);
  F3w1<-1/(12*m_sam[4])
  F3bmw <- matrix(0,m_sam[4],1); F3bmwsl <- matrix(0,m_sam[4],d1)
  for(i in 1:d1){
    F3gg <- (dataF3 - pi[i])/sqrt(gh[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mix_pi[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[4]){F3P2 <- F3P2+runif(m_sam[4])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[4]) + sum((F3P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F3u <- as.matrix(c(12*m_sam[4]*((F3dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F3dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F3dd[3]/m_sam[4]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[4]))))

  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("1MG-NCD",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," "," ",round(sigma[1],4),round(sigma[4],4),round(sigma[5],4),round(sigma[6],4)," "," "," "," "," "," ",
                       round(t(mix_pi),4)," "," "," "," "," "," ",round(B4[1],4)," "," "," ",round(B4[2],4)," "," "," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
###################################2MG-ADI(B-1)######################################################
G4F3ModelFun[[5]] <- function(K1,logL,df11,df21,df31,df41,G4F3text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF3 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF3)))
  sigmaf3<-var(dataF3);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a9<-sqrt(sigmaf3)
  m_nf <- as.numeric(G4F3text2)
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF3)))
  ###################procedure start######################################
  d3<-9
  mi<-as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2))
  m[c(1:4)]<-m
  if (m[1]<m[3]) {a9<--a9}
  m[5]<-m[4]+2.4*a9;m[6]<-m[4]+1.6*a9;m[7]<-m[4]+0.8*a9
  m[8]<-m[4];m[9]<-m[4]-0.8*a9;m[10]<-m[4]-1.6*a9
  m[11]<-m[4]-2.4*a9;m[12]<-m[4]-3.2*a9;m[4]<-m[4]+3.2*a9
  pi<-m[c(4:12)];gh<-sigma[c(4:12)]
  m[4:12] <- c(128,134,80,126,102,117,114,91,90,90)
  m[1:3] <- m[4:6]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mi)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,9,m_sam[4]); swx <- matrix(0,9,1);s0<-matrix(0,4,1);g_a<-matrix(0,9,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:12)];ssigma<-sigma[c(4:12)]
    for(i in 1:d3) { WW[i,] <- mi[i]*dnorm(dataF3,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF3,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF3
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    ###############restrictions####################################
    s0[1]<-sumx[1]+m_nf*sumwx[1];s0[2]<-m_sam[1]+m_nf*n0[1]
    s0[3]<-sumx[3]+m_nf*sumwx[9];s0[4]<-m_sam[3]+m_nf*n0[9]
    aa1<-s0[1]/s0[2]-4*sumx[2]/m_sam[2]+s0[3]/s0[4]+(sumwx[3]/n0[3]+sumwx[7]/n0[7])+16*sumwx[5]/n0[5]-4*(sumwx[2]/n0[2]+sumwx[4]/n0[4]+sumwx[6]/n0[6]+sumwx[8]/n0[8])
    aa2<-(sigma[1]/s0[2]+sigma[1]/s0[4])+16*sigma[1]/m_sam[2]+(sigma[6]/n0[3]+sigma[10]/n0[7])+256*sigma[8]/n0[5]+16*(sigma[5]/n0[2]+sigma[7]/n0[4]+sigma[9]/n0[6]+sigma[11]/n0[8])
    aa3<-aa1/aa2
    m[1]<-(s0[1]-sigma[1]*aa3)/s0[2];m[2]<-(sumx[2]+sigma[1]*aa3*4)/m_sam[2]
    m[3]<-(s0[3]-sigma[1]*aa3)/s0[4];m[5]<-(sumwx[2]+sigma[5]*aa3*4)/n0[2]
    m[6]<-(sumwx[3]-sigma[6]*aa3)/n0[3];m[7]<-(sumwx[4]+sigma[7]*aa3*4)/n0[4]
    m[8]<-(sumwx[5]-sigma[8]*aa3*16)/n0[5];m[9]<-(sumwx[6]+sigma[9]*aa3*4)/n0[6]
    m[10]<-(sumwx[7]-sigma[10]*aa3)/n0[7];m[11]<-(sumwx[8]+sigma[11]*aa3*4)/n0[8]
    m[4]<-m[1];m[12]<-m[3]
    #########first order genetic parameter process##########
    hh1 <- matrix(c(1,1,1,1,1,1,1,1,1,1,1,0,-1,1,1,0,0,0,-1,-1,1,0,-1,0,-1,1,0,-1,1,0,0,1,0,0,0,0.5,0.5,0.5,0,0,
                    0,1,0,0.5,0,0,0.5,0,0,0.5,1,0,1,0,-1,0,0,0,-1,0,0,0,0,0.5,0,0,0,0,0,-0.5,0,0,0,0,0,0.5,0,-0.5,0,0,
                    0,1,0,0,0,0,0.25,0,0,0),10,9)
    B5<-solve(crossprod(hh1,hh1))%*%crossprod(hh1,m[c(1,2,3,5,6,7,8,9,10,11)])
    g_a[c(1,3,7,9)]<-0
    g_a[2]<-(0.5*(B5[3]+B5[6])^2+0.25*(B5[5]+B5[7])^2)/m_nf
    g_a[4]<-(0.5*(B5[2]+B5[6])^2+0.25*(B5[4]+B5[8])^2)/m_nf
    g_a[5]<-0.25*((B5[2]+B5[7])^2+(B5[3]+B5[8])^2+(B5[4]+0.5*B5[9])^2+(B5[5]+0.5*B5[9])^2+0.25*B5[9]^2)/m_nf
    g_a[6]<-(0.5*(B5[2]-B5[6])^2+0.25*(B5[4]-B5[8])^2)/m_nf
    g_a[8]<-(0.5*(B5[3]-B5[6])^2+0.25*(B5[5]-B5[7])^2)/m_nf
    ##############################obtain variance###################################################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d3) {swx[i] <- WW[i,]%*%(dataF3-m[i+3])^2 }
    ################error process##########################################
    a1<-sigma[1]
    n_iter<-0
    a4<-1000
    while(a4>0.0001){
      h1<-(sigma[1]/m_nf)/(sigma[1]/m_nf+g_a)
      s0[1]<-sum(swx*h1^2*m_nf)
      s0[2]<-sum(h1*n0)
      s0[1]<-s0[1]+ss1+ss2+ss3
      s0[2]<-s0[2]+m_sam[1]+m_sam[2]+m_sam[3]
      sigma[1]<-s0[1]/s0[2]
      a4<-abs(sigma[1]-a1)
      a1<-sigma[1]
      if(n_iter>20)break
    }
    sigma[2]<-sigma[3]<-sigma[1]
    sigma[c(4:12)]<-sigma[1]/m_nf+g_a
    #######################criteria for iterations to stop #############################
    pi<-m[c(4:12)];gh<-sigma[c(4:12)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*10
  ########################second order genetic parameter process##########################
  jj<-sigmaf3-sigma[1]/m_nf
  if(jj<0){jj<-0}
  ll<-jj/sigmaf3
  ###########################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma[1]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < m_sam[1]){P1bmw <- P1bmw+runif(m_sam[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:m_sam[1])) - 0.5)/m_sam[1])^2)
  P1u<- as.matrix(c(12*m_sam[1]*((P1dd[1]/m_sam[1]-0.5)^2),((45*m_sam[1])/4)*((P1dd[2]/m_sam[1]-1/3)^2),180*m_sam[1]*((P1dd[3]/m_sam[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,m_sam[1]))))

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*m_sam[2])
  F1bmw <- matrix(0,m_sam[2],1)
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma[2]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < m_sam[2]){F1bmw <- F1bmw+runif(m_sam[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:m_sam[2])) - 0.5)/m_sam[2])^2)
  F1u<- as.matrix(c(12*m_sam[2]*((F1dd[1]/m_sam[2]-0.5)^2),((45*m_sam[2])/4)*((F1dd[2]/m_sam[2]-1/3)^2),180*m_sam[2]*((F1dd[3]/m_sam[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,m_sam[2]))))

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*m_sam[3])
  P2bmw <- matrix(0,m_sam[3],1)
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma[3]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < m_sam[3]){P2bmw <- P2bmw+runif(m_sam[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:m_sam[3])) - 0.5)/m_sam[3])^2)
  P2u<- as.matrix(c(12*m_sam[3]*((P2dd[1]/m_sam[3]-0.5)^2),((45*m_sam[3])/4)*((P2dd[2]/m_sam[3]-1/3)^2),180*m_sam[3]*((P2dd[3]/m_sam[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,m_sam[3]))))

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3);
  F3w1<-1/(12*m_sam[4])
  F3bmw <- matrix(0,m_sam[4],1); F3bmwsl <- matrix(0,m_sam[4],d3)
  for(i in 1:d3){
    F3gg <- (dataF3 - pi[i])/sqrt(gh[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mix_pi[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[4]){F3P2 <- F3P2+runif(m_sam[4])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[4]) + sum((F3P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F3u <- as.matrix(c(12*m_sam[4]*((F3dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F3dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F3dd[3]/m_sam[4]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[4]))))

  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("2MG-ADI",round(abc,4),round(AIC,4),round(t(m),4),round(sigma[1],4),round(sigma[4],4),round(sigma[5],4),round(sigma[6],4),round(sigma[7],4),round(sigma[8],4),round(sigma[9],4),round(sigma[10],4),round(sigma[11],4),round(sigma[12],4),
                       round(t(mix_pi),4),round(B5[1],4)," "," "," ",round(B5[2],4),round(B5[3],4),round(B5[4],4),round(B5[5],4),round(B5[6],4),round(B5[7],4),round(B5[8],4),round(B5[9],4)," "," ",round(jj,4),round(ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
###################################2MG-AD(B-2)######################################################
G4F3ModelFun[[6]] <- function(K1,logL,df11,df21,df31,df41,G4F3text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF3 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF3)))
  sigmaf3<-var(dataF3);sigmaP1<-var(dataP1)
  sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a9<-sqrt(sigmaf3/(m_sam[4]*(m_sam[4]-1)))
  m_nf <- as.numeric(G4F3text2)
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF3)))
  ################procedure start#############################################
  d3<-9
  mi<-as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2))
  m[c(1:4)]<-m
  if (m[1]<m[3]) {a9<--a9}
  m[5]<-m[4]+2.4*a9;m[6]<-m[4]+1.6*a9;m[7]<-m[4]+0.8*a9
  m[8]<-m[4];m[9]<-m[4]-0.8*a9;m[10]<-m[4]-1.6*a9
  m[11]<-m[4]-2.4*a9;m[12]<-m[4]-3.2*a9;m[4]<-m[4]+3.2*a9
  pi<-m[c(4:12)];gh<-sigma[c(4:12)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mi)))
  ##########iteration process###########
  iteration <- 0; stopa <- 1000
  hh<-matrix(0,5,5);b_line<-matrix(0,5,1)
  WW <- matrix(0,9,m_sam[4]); swx <- matrix(0,9,1)
  s0<-matrix(0,4,1); g_a<-matrix(0,9,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:12)];ssigma<-sigma[c(4:12)]
    for(i in 1:d3) { WW[i,] <- mi[i]*dnorm(dataF3,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF3,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF3
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    ###################restrictions#############################
    s0[1]<-sumx[1]+m_nf*sumwx[1];s0[2]<-m_sam[1]+m_nf*n0[1]
    s0[3]<-sumx[3]+m_nf*sumwx[9];s0[4]<-m_sam[3]+m_nf*n0[9]
    #############################################################
    hh[1,1]<-sigma[1]/s0[2]+4.0*sigma[1]/m_sam[2]+sigma[1]/s0[4]+16.0*sigma[8]/n0[5]
    hh[1,2]<-sigma[1]/s0[2]-4.0*sigma[8]/n0[5]
    hh[1,3]<-sigma[1]/s0[2]
    hh[1,4]<-sigma[1]/s0[2]+sigma[1]/s0[4]
    hh[1,5]<-sigma[1]/s0[4]-4.0*sigma[8]/n0[5]
    hh[2,2]<-sigma[1]/s0[2]+sigma[5]/n0[2]+sigma[7]/n0[4]+sigma[8]/n0[5]
    hh[2,3]<-sigma[1]/s0[2]+sigma[5]/n0[2]
    hh[2,4]<-sigma[1]/s0[2]
    hh[2,5]<-sigma[8]/n0[5]
    hh[3,3]<-sigma[1]/s0[2]+sigma[5]/n0[2]+sigma[10]/n0[7]+sigma[11]/n0[8]
    hh[3,4]<-sigma[1]/s0[2]+sigma[10]/n0[7]
    hh[3,5]<--sigma[11]/n0[8]
    hh[4,4]<-sigma[1]/s0[2]+sigma[6]/n0[3]+sigma[10]/n0[3]+sigma[1]/s0[4]
    hh[4,5]<-sigma[1]/s0[4]
    hh[5,5]<-sigma[8]/n0[5]+sigma[9]/n0[6]+sigma[11]/n0[8]+sigma[1]/s0[4]
    for(i in 2:5)
    {
      for(j in 1:(i-1))
      {
        hh[i,j]<-hh[j,i]
      }
    }
    #################################################################
    b_line[1]<-s0[1]/s0[2]+2*sumx[2]/m_sam[2]+s0[3]/s0[4]-4*sumwx[5]/n0[5]
    b_line[2]<-s0[1]/s0[2]-sumwx[2]/n0[2]-sumwx[4]/n0[4]+sumwx[5]/n0[5]
    b_line[3]<-s0[1]/s0[2]-sumwx[2]/n0[2]-sumwx[7]/n0[7]+sumwx[8]/n0[8]
    b_line[4]<-s0[1]/s0[2]-sumwx[3]/n0[3]-sumwx[7]/n0[7]+s0[3]/s0[4]
    b_line[5]<-sumwx[5]/n0[5]-sumwx[6]/n0[6]-sumwx[8]/n0[8]+s0[3]/s0[4]
    B6 <- solve(hh,b_line)
    ################################################################
    m[1]<-(s0[1]-sigma[1]*(B6[1]+B6[2]+B6[3]+B6[4]))/s0[2]
    m[2]<-(sumx[2]-sigma[1]*B6[1]*2)/m_sam[2]
    m[3]<-(s0[3]-sigma[1]*(B6[1]+B6[4]+B6[5]))/s0[4]
    m[5]<-(sumwx[2]+sigma[5]*(B6[2]+B6[3]))/n0[2]
    m[6]<-(sumwx[3]+sigma[6]*B6[4])/n0[3]
    m[7]<-(sumwx[4]+sigma[7]*B6[2])/n0[4]
    m[8]<-(sumwx[5]+sigma[8]*(4*B6[1]-B6[2]-B6[5]))/n0[5]
    m[9]<-(sumwx[6]+sigma[9]*B6[5])/n0[6]
    m[10]<-(sumwx[7]+sigma[10]*(B6[3]+B6[4]))/n0[7]
    m[11]<-(sumwx[8]+sigma[11]*(-B6[3]+B6[5]))/n0[8]
    m[4]<-m[1];m[12]<-m[3]
    #########first order genetic parameter process####################
    hh61 <- matrix(c(1,1,1,1,1,1,1,1,1,1,1,0,-1,1,1,0,0,0,-1,-1,1,0,-1,0,-1,1,0,-1,1,0,0,1,0,0,0,0.5,0.5,0.5,0,0,
                     0,1,0,0.5,0,0,0.5,0,0,0.5),10,5)
    B61<-solve(crossprod(hh61,hh61))%*%crossprod(hh61,m[c(1,2,3,5,6,7,8,9,10,11)])
    g_a[c(1,3,7,9)]<-0
    g_a[2]<-(0.5*B61[3]^2+0.25*B61[5]^2)/m_nf
    g_a[4]<-(0.5*B61[2]^2+0.25*B61[4]^2)/m_nf
    g_a[5]<-(0.5*B61[2]^2+0.5*B61[3]^2+0.25*B61[4]^2+0.25*B61[5]^2)/m_nf
    g_a[6]<-(0.5*B61[2]^2+0.25*B61[4]^2)/m_nf
    g_a[8]<-(0.5*B61[3]^2+0.25*B61[5]^2)/m_nf
    ########################obtain variance###########################################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d3) {swx[i] <- WW[i,]%*%(dataF3-m[i+3])^2 }
    ################error process##########################################
    a1<-sigma[1];n_iter<-0;a4<-1000
    while(a4>0.0001){
      n_iter<-n_iter+1
      h1<-(sigma[1]/m_nf)/(sigma[1]/m_nf+g_a)
      s0[1]<-sum(swx*h1^2*m_nf)
      s0[2]<-sum(h1*n0)
      s0[1]<-s0[1]+ss1+ss2+ss3
      s0[2]<-s0[2]+m_sam[1]+m_sam[2]+m_sam[3]
      sigma[1]<-s0[1]/s0[2]
      a4<-abs(sigma[1]-a1)
      a1<-sigma[1]
      if(n_iter>20) break
    }
    sigma[2]<-sigma[3]<-sigma[1]
    sigma[c(4:12)]<-sigma[1]/m_nf+g_a
    ##############################criteria for iterations to stop########################
    pi<-m[c(4:12)];gh<-sigma[c(4:12)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*6
  ################second order genetic parameter process###############################
  jj<-sigmaf3-sigma[1]/m_nf
  if(jj<0){jj<-0}
  ll<-jj/sigmaf3
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma[1]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < m_sam[1]){P1bmw <- P1bmw+runif(m_sam[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:m_sam[1])) - 0.5)/m_sam[1])^2)
  P1u<- as.matrix(c(12*m_sam[1]*((P1dd[1]/m_sam[1]-0.5)^2),((45*m_sam[1])/4)*((P1dd[2]/m_sam[1]-1/3)^2),180*m_sam[1]*((P1dd[3]/m_sam[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,m_sam[1]))))

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*m_sam[2])
  F1bmw <- matrix(0,m_sam[2],1)
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma[2]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < m_sam[2]){F1bmw <- F1bmw+runif(m_sam[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:m_sam[2])) - 0.5)/m_sam[2])^2)
  F1u<- as.matrix(c(12*m_sam[2]*((F1dd[1]/m_sam[2]-0.5)^2),((45*m_sam[2])/4)*((F1dd[2]/m_sam[2]-1/3)^2),180*m_sam[2]*((F1dd[3]/m_sam[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,m_sam[2]))))

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*m_sam[3])
  P2bmw <- matrix(0,m_sam[3],1)
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma[3]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < m_sam[3]){P2bmw <- P2bmw+runif(m_sam[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:m_sam[3])) - 0.5)/m_sam[3])^2)
  P2u<- as.matrix(c(12*m_sam[3]*((P2dd[1]/m_sam[3]-0.5)^2),((45*m_sam[3])/4)*((P2dd[2]/m_sam[3]-1/3)^2),180*m_sam[3]*((P2dd[3]/m_sam[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,m_sam[3]))))

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3);
  F3w1<-1/(12*m_sam[4])
  F3bmw <- matrix(0,m_sam[4],1); F3bmwsl <- matrix(0,m_sam[4],d3)
  for(i in 1:d3){
    F3gg <- (dataF3 - pi[i])/sqrt(gh[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mix_pi[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[4]){F3P2 <- F3P2+runif(m_sam[4])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[4]) + sum((F3P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F3u <- as.matrix(c(12*m_sam[4]*((F3dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F3dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F3dd[3]/m_sam[4]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[4]))))

  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("2MG-AD",round(abc,4),round(AIC,4),round(t(m),4),round(sigma[1],4),round(sigma[4],4),round(sigma[5],4),round(sigma[6],4),round(sigma[7],4),round(sigma[8],4),round(sigma[9],4),round(sigma[10],4),round(sigma[11],4),round(sigma[12],4),
                       round(t(mix_pi),4),round(B61[1],4)," "," "," ",round(B61[2],4),round(B61[3],4),round(B61[4],4),round(B61[5],4)," "," "," "," "," "," ",round(jj,4),round(ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
######################################2MG-A(B-3)########################################################
G4F3ModelFun[[7]] <- function(K1,logL,df11,df21,df31,df41,G4F3text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF3 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF3)))
  sigmaf3<-var(dataF3);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a9<-sqrt(sigmaf3/(m_sam[4]*(m_sam[4]-1)))


  m_nf <- as.numeric(G4F3text2)
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF3)))
  ####################procedure start###############################
  d3<-9
  a9<-sqrt(sigmaf3/(m_sam[4]*(m_sam[4]-1)))
  mi<-as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2))
  m[c(1:4)]<-m
  if (m[1]<m[3]) {a9<--a9}
  m[5]<-m[4]+2.4*a9;m[6]<-m[4]+1.6*a9;m[7]<-m[4]+0.8*a9
  m[8]<-m[4];m[9]<-m[4]-0.8*a9;m[10]<-m[4]-1.6*a9
  m[11]<-m[4]-2.4*a9;m[12]<-m[4]-3.2*a9;m[4]<-m[4]+3.2*a9
  pi<-m[c(4:12)];gh<-sigma[c(4:12)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mi)))
  ##################################iteration process########################################################
  iteration <- 0; stopa <- 1000
  hh<-matrix(0,6,6);b_line<-matrix(0,6,1)
  WW <- matrix(0,9,m_sam[4]); swx <- matrix(0,9,1)
  s0<-matrix(0,6,1); g_a<-matrix(0,9,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:12)];ssigma<-sigma[c(4:12)]
    for(i in 1:d3) { WW[i,] <- mi[i]*dnorm(dataF3,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF3,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF3
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    #############restrictions###################################
    s0[1]<-sumx[1]+m_nf*sumwx[1];s0[2]<-m_sam[1]+m_nf*n0[1];s0[3]<-sumx[2]+m_nf*sumwx[5]
    s0[4]<-m_sam[2]+m_nf*n0[5];s0[5]<-sumx[3]+m_nf*sumwx[9];s0[6]<-m_sam[3]+m_nf*n0[9]
    ############################################################################
    hh[1,1]<-sigma[1]/s0[2]+sigma[1]/s0[6]+sigma[6]/n0[3]+sigma[10]/n0[7]
    hh[1,2]<-sigma[1]/s0[2]+sigma[10]/n0[7]
    hh[1,3]<--sigma[6]/n0[3]-sigma[10]/n0[7]
    hh[1,4]<-2.0*sigma[6]/n0[3]
    hh[1,5]<-sigma[1]/s0[6]-sigma[10]/n0[7]
    hh[1,6]<-sigma[1]/s0[2]-sigma[10]/n0[7]
    hh[2,2]<-sigma[1]/s0[2]+sigma[5]/n0[2]+sigma[10]/n0[7]+sigma[11]/n0[8]
    hh[2,3]<--sigma[10]/n0[7]
    hh[2,4]<--2.0*sigma[5]/n0[2]
    hh[2,5]<--sigma[10]/n0[7]-2.0*sigma[11]/n0[8]
    hh[2,6]<-sigma[1]/s0[2]-sigma[10]/n0[7]
    hh[3,3]<-4*sigma[1]/s0[4]+sigma[6]/n0[3]+sigma[10]/n0[7]
    hh[3,4]<--2*sigma[6]/n0[3]
    hh[3,5]<-sigma[10]/n0[7]
    hh[3,6]<-sigma[10]/n0[7]
    hh[4,4]<-4*sigma[5]/n0[2]+4*sigma[6]/n0[3]+sigma[7]/n0[4]+sigma[9]/n0[6]
    hh[4,5]<-0
    hh[4,6]<-2*sigma[7]/n0[4]
    hh[5,5]<-sigma[1]/s0[6]+sigma[10]/n0[7]+4.0*sigma[11]/n0[8]
    hh[5,6]<-sigma[10]/n0[7]
    hh[6,6]<-sigma[1]/s0[2]+4.0*sigma[7]/n0[4]+sigma[10]/n0[7]
    for(i in 2:6)
    {
      for(j in 1:(i-1))
      {
        hh[i,j]<-hh[j,i]
      }
    }
    ####################################################################################################
    b_line[1]<-s0[1]/s0[2]-sumwx[3]/n0[3]-sumwx[7]/n0[7]+s0[5]/s0[6]
    b_line[2]<-s0[1]/s0[2]-sumwx[2]/n0[2]-sumwx[7]/n0[7]+sumwx[8]/n0[8]
    b_line[3]<-sumwx[3]/n0[3]-2*s0[3]/s0[4]+sumwx[7]/n0[7]
    b_line[4]<-2*sumwx[2]/n0[2]-2*sumwx[3]/n0[3]-sumwx[4]/n0[4]+sumwx[6]/n0[6]
    b_line[5]<-s0[5]/s0[6]-2*sumwx[8]/n0[8]+sumwx[7]/n0[7]
    b_line[6]<-s0[1]/s0[2]-2*sumwx[4]/n0[4]+sumwx[7]/n0[7]
    B7 <- solve(hh,b_line)
    ###########################################################################################3
    m[1]<-(s0[1]-sigma[1]*(B7[1]+B7[2]+B7[6]))/s0[2]
    m[2]<-(s0[3]+sigma[1]*B7[3]*2.0)/s0[4]
    m[3]<-(s0[5]-sigma[1]*(B7[1]+B7[5]))/s0[6]
    m[5]<-(sumwx[2]+sigma[5]*(B7[2]-2.0*B7[4]))/n0[2]
    m[6]<-(sumwx[3]+sigma[6]*(B7[1]-B7[3]+2.0*B7[4]))/n0[3]
    m[7]<-(sumwx[4]+sigma[7]*(B7[4]+2.0*B7[6]))/n0[4]
    m[9]<-(sumwx[6]-sigma[9]*B7[4])/n0[6]
    m[10]<-(sumwx[7]+sigma[10]*(B7[1]+B7[2]-B7[3]-B7[5]-B7[6]))/n0[7]
    m[11]<-(sumwx[8]-sigma[11]*(B7[2]-2.0*B7[5]))/n0[8]
    m[4]<-m[1]; m[8]<-m[2];m[12]<-m[3]
    #########first order genetic parameter process##########
    hh71 <- matrix(c(1,1,1,1,1,1,1,1,1,1,1,0,-1,1,1,0,0,0,-1,-1,1,0,-1,0,-1,1,0,-1,1,0),10,3)
    B71<-solve(crossprod(hh71,hh71))%*%crossprod(hh71,m[c(1,2,3,5,6,7,8,9,10,11)])
    g_a[c(1,3,7,9)]<-0
    g_a[2]<-0.5*B71[3]^2/m_nf
    g_a[4]<-0.5*B71[2]^2/m_nf
    g_a[5]<-(0.5*B71[2]^2+0.5*B71[3]^2)/m_nf
    g_a[6]<-0.5*B71[2]^2/m_nf
    g_a[8]<-0.5*B71[3]^2/m_nf
    #################obtain variance#############################################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d3) {swx[i] <- WW[i,]%*%(dataF3-m[i+3])^2 }
    ################error process##########################################
    a1<-sigma[1];n_iter<-0;a4<-1000
    while(a4>0.0001){
      n_iter<-n_iter+1
      h1<-(sigma[1]/m_nf)/(sigma[1]/m_nf+g_a)
      s0[1]<-sum(swx*h1^2*m_nf)
      s0[2]<-sum(h1*n0)
      s0[1]<-s0[1]+ss1+ss2+ss3
      s0[2]<-s0[2]+m_sam[1]+m_sam[2]+m_sam[3]
      sigma[1]<-s0[1]/s0[2]
      a4<-abs(sigma[1]-a1)
      a1<-sigma[1]
      if(n_iter>20) break
    }
    sigma[2]<-sigma[3]<-sigma[1]
    sigma[c(4:12)]<-sigma[1]/m_nf+g_a
    ####################criteria for iterations to stop############################
    pi<-m[c(4:12)];gh<-sigma[c(4:12)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*4
  ##############second order genetic parameter process############################
  jj<-sigmaf3-sigma[1]/m_nf
  if(jj<0){jj<-0}
  ll<-jj/sigmaf3
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma[1]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < m_sam[1]){P1bmw <- P1bmw+runif(m_sam[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:m_sam[1])) - 0.5)/m_sam[1])^2)
  P1u<- as.matrix(c(12*m_sam[1]*((P1dd[1]/m_sam[1]-0.5)^2),((45*m_sam[1])/4)*((P1dd[2]/m_sam[1]-1/3)^2),180*m_sam[1]*((P1dd[3]/m_sam[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,m_sam[1]))))

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*m_sam[2])
  F1bmw <- matrix(0,m_sam[2],1)
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma[2]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < m_sam[2]){F1bmw <- F1bmw+runif(m_sam[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:m_sam[2])) - 0.5)/m_sam[2])^2)
  F1u<- as.matrix(c(12*m_sam[2]*((F1dd[1]/m_sam[2]-0.5)^2),((45*m_sam[2])/4)*((F1dd[2]/m_sam[2]-1/3)^2),180*m_sam[2]*((F1dd[3]/m_sam[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,m_sam[2]))))

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*m_sam[3])
  P2bmw <- matrix(0,m_sam[3],1)
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma[3]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < m_sam[3]){P2bmw <- P2bmw+runif(m_sam[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:m_sam[3])) - 0.5)/m_sam[3])^2)
  P2u<- as.matrix(c(12*m_sam[3]*((P2dd[1]/m_sam[3]-0.5)^2),((45*m_sam[3])/4)*((P2dd[2]/m_sam[3]-1/3)^2),180*m_sam[3]*((P2dd[3]/m_sam[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,m_sam[3]))))

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3);
  F3w1<-1/(12*m_sam[4])
  F3bmw <- matrix(0,m_sam[4],1); F3bmwsl <- matrix(0,m_sam[4],d3)
  for(i in 1:d3){
    F3gg <- (dataF3 - pi[i])/sqrt(gh[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mix_pi[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[4]){F3P2 <- F3P2+runif(m_sam[4])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[4]) + sum((F3P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F3u <- as.matrix(c(12*m_sam[4]*((F3dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F3dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F3dd[3]/m_sam[4]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[4]))))

  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("2MG-A",round(abc,4),round(AIC,4),round(t(m),4),round(sigma[1],4),round(sigma[4],4),round(sigma[5],4),round(sigma[6],4),round(sigma[7],4),round(sigma[8],4),round(sigma[9],4),round(sigma[10],4),round(sigma[11],4),round(sigma[12],4),
                       round(t(mix_pi),4),round(B71[1],4)," "," "," ",round(B71[2],4),round(B71[3],4)," "," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
##################################2MG-EA(B-4)#####################################################################
G4F3ModelFun[[8]] <- function(K1,logL,df11,df21,df31,df41,G4F3text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF3 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF3)))
  sigmaf3<-var(dataF3);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a9<-sqrt(sigmaf3/(m_sam[4]*(m_sam[4]-1)))
  m_nf <- as.numeric(G4F3text2)
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF3)))
  ################procedure start#################################
  d2<-6
  mi<-as.matrix(c(0.0625,0.25,0.125,0.25,0.25,0.0625))
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2))
  m[c(1:4)]<-m
  if (m[1]<m[3]) {a9<--a9}
  m[5]<-m[4]+1.5*a9;m[6]<-m[4];m[7]<-sumx[4]/m_sam[4]
  m[8]<-m[4]-1.5*a9;m[9]<-m[4]-3*a9;m[4]<-m[4]+3*a9
  pi<-m[c(4:9)];gh<-sigma[c(4:9)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+ sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mi)))
  ##################################iteration process########################################################
  iteration <- 0; stopa <- 1000
  hh<-matrix(0,3,3);b_line<-matrix(0,3,1)
  WW <- matrix(0,6,m_sam[4]); swx <- matrix(0,6,1)
  s0<-matrix(0,6,1);g_a<-matrix(0,6,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:9)];ssigma<-sigma[c(4:9)]
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataF3,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF3,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF3
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    ################restrictions#################################
    s0[1]<-sumx[1]+m_nf*sumwx[1];s0[2]<-m_sam[1]+m_nf*n0[1];s0[3]<-sumx[2]+m_nf*sumwx[3]+sigma[1]*sumwx[4]/sigma[7];
    s0[4]<-m_sam[2]+m_nf*n0[3]+sigma[1]*n0[4]/sigma[7];s0[5]<-sumx[3]+m_nf*sumwx[6];s0[6]<-m_sam[3]+m_nf*n0[6]
    ###########################################################################################
    hh[1,1]<-sigma[1]/s0[2]+4*sigma[1]/s0[4]+sigma[1]/s0[6]
    hh[1,2]<-4*sigma[1]/s0[4]
    hh[1,3]<-sigma[1]/s0[2]-2*sigma[1]/s0[4]
    hh[2,2]<-sigma[5]/n0[2]+4*sigma[1]/s0[4]+sigma[8]/n0[5]
    hh[2,3]<--2*sigma[1]/s0[4]-2*sigma[5]/n0[2]
    hh[3,3]<-sigma[1]/s0[2]+4*sigma[5]/n0[2]+sigma[1]/s0[4]
    for(i in 2:3)
    {
      for(j in 1:(i-1))
      {
        hh[i,j]<-hh[j,i]
      }
    }
    #########################################################################################
    b_line[1]<-s0[1]/s0[2]-2*s0[3]/s0[4]+s0[5]/s0[6]
    b_line[2]<-sumwx[2]/n0[2]-2*s0[3]/s0[4]+sumwx[5]/n0[5]
    b_line[3]<--2*sumwx[2]/n0[2]+s0[1]/s0[2]+s0[3]/s0[4]
    B8 <- solve(hh,b_line)
    #########################################################################################
    m[1]<-(s0[1]-sigma[1]*(B8[1]+B8[3]))/s0[2]
    m[2]<-(s0[3]+sigma[1]*(B8[1]*2.0+2.0*B8[2]-B8[3]))/s0[4]
    m[3]<-(s0[5]-sigma[1]*B8[1])/s0[6]
    m[5]<-(sumwx[2]-sigma[5]*(B8[2]-2.0*B8[3]))/n0[2]
    m[8]<-(sumwx[5]-sigma[8]*B8[2])/n0[5]
    m[4]<-m[1];m[6]<-m[2];m[7]<-m[2];m[9]<-m[3]
    ##########################first order genetic parameter process#######################
    hh81 <- matrix(c(1,1,1,1,1,2,0,-2,0,-1),5,2)
    B81<-solve(crossprod(hh81,hh81))%*%crossprod(hh81,m[c(1,2,3,6,8)])
    g_a[c(1,3,6)]<-0
    g_a[2]<-0.5*B81[2]^2/m_nf
    g_a[4]<-B81[2]^2/m_nf
    g_a[5]<-0.5*B81[2]^2/m_nf
    ###################obtain variance############################################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d2) {swx[i] <- WW[i,]%*%(dataF3-m[i+3])^2 }
    ################error process##########################################
    a1<-sigma[1];n_iter<-0;a4<-1000
    while(a4>0.0001){
      n_iter<-n_iter+1
      h1<-(sigma[1]/m_nf)/(sigma[1]/m_nf+g_a)
      s0[1]<-sum(swx*h1^2*m_nf)
      s0[2]<-sum(h1*n0)
      s0[1]<-s0[1]+ss1+ss2+ss3
      s0[2]<-s0[2]+m_sam[1]+m_sam[2]+m_sam[3]
      sigma[1]<-s0[1]/s0[2]
      a4<-abs(sigma[1]-a1)
      a1<-sigma[1]
      if(n_iter>20) break
    }
    sigma[2]<-sigma[3]<-sigma[1]
    sigma[c(4:9)]<-sigma[1]/m_nf+g_a
    ####################criteria for iterations to stop######################
    pi<-m[c(4:9)];gh<-sigma[c(4:9)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*4
  ###################second order genetic parameter process#######################################
  jj<-sigmaf3-sigma[1]/m_nf
  if(jj<0){jj<-0}
  ll<-jj/sigmaf3
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma[1]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < m_sam[1]){P1bmw <- P1bmw+runif(m_sam[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:m_sam[1])) - 0.5)/m_sam[1])^2)
  P1u<- as.matrix(c(12*m_sam[1]*((P1dd[1]/m_sam[1]-0.5)^2),((45*m_sam[1])/4)*((P1dd[2]/m_sam[1]-1/3)^2),180*m_sam[1]*((P1dd[3]/m_sam[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,m_sam[1]))))

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*m_sam[2])
  F1bmw <- matrix(0,m_sam[2],1)
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma[2]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < m_sam[2]){F1bmw <- F1bmw+runif(m_sam[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:m_sam[2])) - 0.5)/m_sam[2])^2)
  F1u<- as.matrix(c(12*m_sam[2]*((F1dd[1]/m_sam[2]-0.5)^2),((45*m_sam[2])/4)*((F1dd[2]/m_sam[2]-1/3)^2),180*m_sam[2]*((F1dd[3]/m_sam[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,m_sam[2]))))

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*m_sam[3])
  P2bmw <- matrix(0,m_sam[3],1)
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma[3]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < m_sam[3]){P2bmw <- P2bmw+runif(m_sam[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:m_sam[3])) - 0.5)/m_sam[3])^2)
  P2u<- as.matrix(c(12*m_sam[3]*((P2dd[1]/m_sam[3]-0.5)^2),((45*m_sam[3])/4)*((P2dd[2]/m_sam[3]-1/3)^2),180*m_sam[3]*((P2dd[3]/m_sam[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,m_sam[3]))))

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3);
  F3w1<-1/(12*m_sam[4])
  F3bmw <- matrix(0,m_sam[4],1); F3bmwsl <- matrix(0,m_sam[4],d2)
  for(i in 1:d2){
    F3gg <- (dataF3 - pi[i])/sqrt(gh[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mix_pi[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[4]){F3P2 <- F3P2+runif(m_sam[4])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[4]) + sum((F3P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F3u <- as.matrix(c(12*m_sam[4]*((F3dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F3dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F3dd[3]/m_sam[4]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[4]))))

  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("2MG-EA",round(abc,4),round(AIC,4),round(t(m),4)," "," ", " ",round(sigma[1],4),round(sigma[4],4),round(sigma[5],4),round(sigma[6],4),round(sigma[7],4),round(sigma[8],4),round(sigma[9],4)," "," "," ",
                       round(t(mix_pi),4)," "," "," ",round(B81[1],4)," "," "," ",round(B81[2],4)," "," "," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
#################################2MG-CD(B-5)#######################################
G4F3ModelFun[[9]] <- function(K1,logL,df11,df21,df31,df41,G4F3text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF3 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF3)))
  sigmaf3<-var(dataF3);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a9<-sqrt(sigmaf3/(m_sam[4]*(m_sam[4]-1)))
  m_nf <- as.numeric(G4F3text2)
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF3)))
  ############procedure start################################################
  d3<-9
  a9<-sqrt(sigmaf3/(m_sam[4]*(m_sam[4]-1)))
  mi<-as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2))
  m[c(1:4)]<-m
  if (m[1]<m[3]) {a9<--a9}
  m[5]<-m[4]+2.4*a9;m[6]<-m[4]+1.6*a9;m[7]<-m[4]+0.8*a9
  m[8]<-m[4];m[9]<-m[4]-0.8*a9;m[10]<-m[4]-1.6*a9
  m[11]<-m[4]-2.4*a9;m[12]<-m[4]-3.2*a9;m[4]<-m[4]+3.2*a9
  pi<-m[c(4:12)];gh<-sigma[c(4:12)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mi)))
  #####################iteration process####################################
  iteration <- 0; stopa <- 1000
  hh<-matrix(0,6,6);b_line<-matrix(0,6,1)
  WW <- matrix(0,9,m_sam[4]); swx <- matrix(0,9,1)
  s0<-matrix(0,4,1);g_a<-matrix(0,9,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:12)];ssigma<-sigma[c(4:12)]
    for(i in 1:d3) { WW[i,] <- mi[i]*dnorm(dataF3,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF3,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF3
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    ###############restrictions##############################################
    s0[1]<-sumx[1]+sumx[2]+m_nf*sumwx[1];s0[2]<-m_sam[1]+m_sam[2]+m_nf*n0[1]
    s0[3]<-sumx[3]+m_nf*sumwx[9];s0[4]<-m_sam[3]+m_nf*n0[9]
    ###################################################################
    hh[1,1]<-9.0*sigma[1]/s0[2]+16.0*sigma[5]/n0[2]+sigma[6]/n0[3]
    hh[1,2]<-0.0
    hh[1,3]<--12.0*sigma[5]/n0[2]
    hh[1,4]<-3.0*sigma[1]/s0[2]+4.0*sigma[5]/n0[2]
    hh[1,5]<-3.0*sigma[1]/s0[2]-sigma[6]/n0[3]
    hh[1,6]<-0.0
    hh[2,2]<-9.0*sigma[7]/n0[4]+16.0*sigma[8]/n0[5]+sigma[9]/n0[6]
    hh[2,3]<-16.0*sigma[8]/n0[5]
    hh[2,4]<-0.0
    hh[2,5]<--3.0*sigma[7]/n0[4]+sigma[9]/n0[6]
    hh[2,6]<--4.0*sigma[8]/n0[5]-sigma[9]/n0[6]
    hh[3,3]<-9.0*sigma[5]/n0[2]+16.0*sigma[8]/n0[5]+sigma[11]/n0[8]
    hh[3,4]<--3.0*sigma[5]/n0[2]+sigma[11]/n0[8]
    hh[3,5]<-0.0
    hh[3,6]<--4.0*sigma[8]/n0[5]-sigma[11]/n0[8]
    hh[4,4]<-sigma[1]/s0[2]+sigma[5]/n0[2]+sigma[10]/n0[7]+sigma[11]/n0[8]
    hh[4,5]<-sigma[1]/s0[2]
    hh[4,6]<--sigma[11]/n0[8]
    hh[5,5]<-sigma[1]/s0[2]+sigma[6]/n0[3]+sigma[7]/n0[4]+sigma[9]/n0[6]
    hh[5,6]<--sigma[9]/n0[6]
    hh[6,6]<-sigma[8]/n0[5]+sigma[9]/n0[6]+sigma[11]/n0[8]+sigma[1]/s0[4]
    for(i in 2:6)
    {
      for(j in 1:(i-1))
      {
        hh[i,j]<-hh[j,i]
      }
    }
    #############################################################################################
    b_line[1]<-3.0*s0[1]/s0[2]-4.0*sumwx[2]/n0[2]+sumwx[3]/n0[3]
    b_line[2]<-3.0*sumwx[4]/n0[4]-4.0*sumwx[5]/n0[5]+sumwx[6]/n0[6]
    b_line[3]<-3.0*sumwx[2]/n0[2]-4.0*sumwx[5]/n0[5]+sumwx[8]/n0[8]
    b_line[4]<-s0[1]/s0[2]-sumwx[2]/n0[2]-sumwx[7]/n0[7]+sumwx[8]/n0[8]
    b_line[5]<-s0[1]/s0[2]-sumwx[3]/n0[3]-sumwx[4]/n0[4]+sumwx[6]/n0[6]
    b_line[6]<-sumwx[5]/n0[5]-sumwx[6]/n0[6]-sumwx[8]/n0[8]+s0[3]/s0[4]
    B9<-solve(hh,b_line)
    #############################################################################################
    m[1]<-(s0[1]-sigma[1]*(B9[1]*3+B9[4]+B9[5]))/s0[2]
    m[3]<-(s0[3]-sigma[1]*B9[6])/s0[4]
    m[5]<-(sumwx[2]+sigma[5]*(4*B9[1]-3*B9[3]+B9[4]))/n0[2]
    m[6]<-(sumwx[3]+sigma[6]*(-B9[1]+B9[5]))/n0[3]
    m[7]<-(sumwx[4]+sigma[7]*(-3*B9[2]+B9[5]))/n0[4]
    m[8]<-(sumwx[5]+sigma[8]*(4*B9[2]+4*B9[3]-B9[6]))/n0[5]
    m[9]<-(sumwx[6]+sigma[9]*(-B9[2]-B9[5]+B9[6]))/n0[6]
    m[10]<-(sumwx[7]+sigma[10]*B9[4])/n0[7]
    m[11]<-(sumwx[8]+sigma[11]*(-B9[3]-B9[4]+B9[6]))/n0[8]
    m[2]<-m[1]; m[4]<-m[1];m[12]<-m[3]
    ##########################################################################################
    hh91<-matrix(c(1,1,1,1,1,1,1,1,1,1,-1,1,1,0.5,0.5,0.5,-1,-1,1,-1,0.5,-1,1,0.5,-1,1,0.5),9,3)
    B91<-solve(crossprod(hh91,hh91))%*%crossprod(hh91,m[c(1,3,5,6,7,8,9,10,11)])
    g_a[c(1,3,7,9)]<-0
    g_a[2]<-0.75*B91[3]^2/m_nf
    g_a[4]<-0.75*B91[2]^2/m_nf
    g_a[5]<-(0.75*B91[2]^2+0.75*B91[3]^2)/m_nf
    g_a[6]<-0.75*B91[2]^2/m_nf
    g_a[8]<-0.75*B91[3]^2/m_nf
    ############obtain variance#############################################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d3) {swx[i] <- WW[i,]%*%(dataF3-m[i+3])^2 }
    #################error process##########################################
    a1<-sigma[1];n_iter<-0;a4<-1000
    while(a4>0.0001){
      n_iter<-n_iter+1
      h1<-(sigma[1]/m_nf)/(sigma[1]/m_nf+g_a)
      s0[1]<-sum(swx*h1^2*m_nf)
      s0[2]<-sum(h1*n0)
      s0[1]<-s0[1]+ss1+ss2+ss3
      s0[2]<-s0[2]+m_sam[1]+m_sam[2]+m_sam[3]
      sigma[1]<-s0[1]/s0[2]
      a4<-abs(sigma[1]-a1)
      a1<-sigma[1]
      if(n_iter>20) break
    }
    sigma[2]<-sigma[3]<-sigma[1]
    sigma[c(4:12)]<-sigma[1]/m_nf+g_a
    ####################criteria for iterations to stop#######################
    pi<-m[c(4:12)];gh<-sigma[c(4:12)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*5
  #####################second order genetic parameter process#####################
  jj<-sigmaf3-sigma[1]/m_nf
  if(jj<0){jj<-0}
  ll<-jj/sigmaf3
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma[1]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < m_sam[1]){P1bmw <- P1bmw+runif(m_sam[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:m_sam[1])) - 0.5)/m_sam[1])^2)
  P1u<- as.matrix(c(12*m_sam[1]*((P1dd[1]/m_sam[1]-0.5)^2),((45*m_sam[1])/4)*((P1dd[2]/m_sam[1]-1/3)^2),180*m_sam[1]*((P1dd[3]/m_sam[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,m_sam[1]))))

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*m_sam[2])
  F1bmw <- matrix(0,m_sam[2],1)
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma[2]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < m_sam[2]){F1bmw <- F1bmw+runif(m_sam[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:m_sam[2])) - 0.5)/m_sam[2])^2)
  F1u<- as.matrix(c(12*m_sam[2]*((F1dd[1]/m_sam[2]-0.5)^2),((45*m_sam[2])/4)*((F1dd[2]/m_sam[2]-1/3)^2),180*m_sam[2]*((F1dd[3]/m_sam[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,m_sam[2]))))

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*m_sam[3])
  P2bmw <- matrix(0,m_sam[3],1)
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma[3]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < m_sam[3]){P2bmw <- P2bmw+runif(m_sam[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:m_sam[3])) - 0.5)/m_sam[3])^2)
  P2u<- as.matrix(c(12*m_sam[3]*((P2dd[1]/m_sam[3]-0.5)^2),((45*m_sam[3])/4)*((P2dd[2]/m_sam[3]-1/3)^2),180*m_sam[3]*((P2dd[3]/m_sam[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,m_sam[3]))))

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3);
  F3w1<-1/(12*m_sam[4])
  F3bmw <- matrix(0,m_sam[4],1); F3bmwsl <- matrix(0,m_sam[4],d3)
  for(i in 1:d3){
    F3gg <- (dataF3 - pi[i])/sqrt(gh[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mix_pi[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[4]){F3P2 <- F3P2+runif(m_sam[4])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[4]) + sum((F3P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F3u <- as.matrix(c(12*m_sam[4]*((F3dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F3dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F3dd[3]/m_sam[4]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[4]))))

  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("2MG-CD",round(abc,4),round(AIC,4),round(t(m),4),round(sigma[1],4),round(sigma[4],4),round(sigma[5],4),round(sigma[6],4),round(sigma[7],4),round(sigma[8],4),round(sigma[9],4),round(sigma[10],4),round(sigma[11],4),round(sigma[12],4),
                       round(t(mix_pi),4),round(B91[1],4)," "," "," ",round(B91[2],4),round(B91[3],4)," "," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
#############################2MG-EAD(B_6)######################################################
G4F3ModelFun[[10]] <- function(K1,logL,df11,df21,df31,df41,G4F3text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF3 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF3)))
  sigmaf3<-var(dataF3);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a9<-sqrt(sigmaf3/(m_sam[4]*(m_sam[4]-1)))
  m_nf <- as.numeric(G4F3text2)
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF3)))
  #####################procedure start########################################
  d2<-6
  mi<-as.matrix(c(0.0625,0.25,0.125,0.25,0.25,0.0625))
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2))
  m[c(1:4)]<-m
  if (m[1]<m[3]) {a9<--a9}
  m[5]<-m[4]+2.1*a9;m[6]<-m[4]+1.2*a9;m[7]<-m[4]
  m[8]<-m[4]-0.7*a9;m[9]<-m[4]-3*a9;m[4]<-m[4]+3*a9
  pi<-m[c(4:9)];gh<-sigma[c(4:9)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mi)))
  #########################iteration process#################################
  iteration <- 0; stopa <- 1000
  hh<-matrix(0,4,4);b_line<-matrix(0,4,1)
  WW <- matrix(0,6,m_sam[4]); swx <- matrix(0,6,1)
  s0<-matrix(0,6,1); g_a<-matrix(0,6,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:9)];ssigma<-sigma[c(4:9)]
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataF3,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF3,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF3
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    ######################restrictions###################
    s0[1]<-sumx[1]+m_nf*sumwx[1];s0[2]<-m_sam[1]+m_nf*n0[1];s0[3]<-sumx[2]+m_nf*sumwx[4];
    s0[4]<-m_sam[2]+m_nf*n0[4];s0[5]<-sumx[3]+m_nf*sumwx[6];s0[6]<-m_sam[3]+m_nf*n0[6];
    ######################################################################
    hh[1,1]<-sigma[1]/s0[2]+4*sigma[5]/n0[2]+sigma[6]/n0[3]
    hh[1,2]<-sigma[1]/s0[2]
    hh[1,3]<--2*sigma[5]/n0[2]-sigma[6]/n0[3]
    hh[1,4]<--2*sigma[6]/n0[3]
    hh[2,2]<-sigma[1]/s0[2]+4*sigma[1]/s0[4]+sigma[1]/s0[6]
    hh[2,3]<-2*sigma[1]/s0[4]
    hh[2,4]<--6*sigma[1]/s0[4]-sigma[1]/s0[6]
    hh[3,3]<-sigma[1]/s0[4]+sigma[5]/n0[2]+sigma[6]/n0[3]+sigma[8]/n0[5]
    hh[3,4]<-2*sigma[6]/n0[3]-3*sigma[1]/s0[4]
    hh[4,4]<-9.0*sigma[1]/s0[4]+sigma[1]/s0[6]+4.0*sigma[6]/n0[3]
    for(i in 2:4)
    {
      for(j in 1:(i-1))
      {
        hh[i,j]<-hh[j,i]
      }
    }
    ##########################################################################
    b_line[1]<-s0[1]/s0[2]-2*sumwx[2]/n0[2]+sumwx[3]/n0[3]
    b_line[2]<-s0[1]/s0[2]-2*s0[3]/s0[4]+s0[5]/s0[6]
    b_line[3]<-sumwx[2]/n0[2]-sumwx[3]/n0[3]-s0[3]/s0[4]+sumwx[5]/n0[5]
    b_line[4]<-3*s0[3]/s0[4]-s0[5]/s0[6]-2*sumwx[3]/n0[3]
    B10<-solve(hh,b_line)
    ###########################################################################
    m[1]<-(s0[1]-sigma[1]*(B10[1]+B10[2]))/s0[2]
    m[2]<-(s0[3]+sigma[1]*(B10[2]*2.0+B10[3]-3.0*B10[4]))/s0[4]
    m[3]<-(s0[5]-sigma[1]*(-B10[2]+B10[4]))/s0[6]
    m[5]<-(sumwx[2]+sigma[5]*(2.0*B10[1]-B10[3]))/n0[2]
    m[6]<-(sumwx[3]+sigma[6]*(-B10[1]+B10[3]+2.0*B10[4]))/n0[3]
    m[8]<-(sumwx[5]-sigma[8]*B10[3])/n0[5]
    m[4]<-m[1]; m[7]<-m[2];m[9]<-m[3]
    #########################################################################
    hh101 <- matrix(c(1,1,1,1,1,1,2,-2,1.5,1,0,-0.5),6,2)
    B101<-solve(crossprod(hh101,hh101))%*%crossprod(hh101,m[c(1,3,5,6,7,8)])
    g_a[c(1,4,6)]<-0
    g_a[2]<-0.75*B101[2]^2/m_nf
    g_a[3]<-1.50*B101[2]^2/m_nf
    g_a[5]<-0.75*B101[2]^2/m_nf
    #################obtain variance#############################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d2) {swx[i] <- WW[i,]%*%(dataF3-m[i+3])^2 }
    ################restrictions#######################################
    a1<-sigma[1];n_iter<-0; a4<-1000
    while(a4>0.0001){
      n_iter<-n_iter+1
      h1<-(sigma[1]/m_nf)/(sigma[1]/m_nf+g_a)
      s0[1]<-sum(swx*h1^2*m_nf)
      s0[2]<-sum(h1*n0)
      s0[1]<-s0[1]+ss1+ss2+ss3
      s0[2]<-s0[2]+m_sam[1]+m_sam[2]+m_sam[3]
      sigma[1]<-s0[1]/s0[2]
      a4<-abs(sigma[1]-a1)
      a1<-sigma[1]
      if(n_iter>20) break
    }
    sigma[2]<-sigma[3]<-sigma[1]
    sigma[c(4:9)]<-sigma[1]/m_nf+g_a
    ###############criteria for iterations to stop####################
    pi<-m[c(4:9)];gh<-sigma[c(4:9)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*4
  ####################second order genetic parameter process################################
  jj<-sigmaf3-sigma[1]/m_nf
  if(jj<0){jj<-0}
  ll<-jj/sigmaf3
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma[1]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < m_sam[1]){P1bmw <- P1bmw+runif(m_sam[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:m_sam[1])) - 0.5)/m_sam[1])^2)
  P1u<- as.matrix(c(12*m_sam[1]*((P1dd[1]/m_sam[1]-0.5)^2),((45*m_sam[1])/4)*((P1dd[2]/m_sam[1]-1/3)^2),180*m_sam[1]*((P1dd[3]/m_sam[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,m_sam[1]))))

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*m_sam[2])
  F1bmw <- matrix(0,m_sam[2],1)
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma[2]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < m_sam[2]){F1bmw <- F1bmw+runif(m_sam[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:m_sam[2])) - 0.5)/m_sam[2])^2)
  F1u<- as.matrix(c(12*m_sam[2]*((F1dd[1]/m_sam[2]-0.5)^2),((45*m_sam[2])/4)*((F1dd[2]/m_sam[2]-1/3)^2),180*m_sam[2]*((F1dd[3]/m_sam[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,m_sam[2]))))

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*m_sam[3])
  P2bmw <- matrix(0,m_sam[3],1)
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma[3]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < m_sam[3]){P2bmw <- P2bmw+runif(m_sam[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:m_sam[3])) - 0.5)/m_sam[3])^2)
  P2u<- as.matrix(c(12*m_sam[3]*((P2dd[1]/m_sam[3]-0.5)^2),((45*m_sam[3])/4)*((P2dd[2]/m_sam[3]-1/3)^2),180*m_sam[3]*((P2dd[3]/m_sam[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,m_sam[3]))))

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3);
  F3w1<-1/(12*m_sam[4])
  F3bmw <- matrix(0,m_sam[4],1); F3bmwsl <- matrix(0,m_sam[4],d2)
  for(i in 1:d2){
    F3gg <- (dataF3 - pi[i])/sqrt(gh[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mix_pi[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[4]){F3P2 <- F3P2+runif(m_sam[4])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[4]) + sum((F3P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F3u <- as.matrix(c(12*m_sam[4]*((F3dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F3dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F3dd[3]/m_sam[4]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[4]))))

  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("2MG-EAD",round(abc,4),round(AIC,4),round(t(m),4)," "," "," ",round(sigma[1],4),round(sigma[4],4),round(sigma[5],4),round(sigma[6],4),round(sigma[7],4),round(sigma[8],4),round(sigma[9],4)," "," "," ",
                       round(t(mix_pi),4)," "," "," ",round(B101[1],4)," "," "," ",round(B101[2],4)," "," "," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
#################################PG-ADI(C-0)#####################################################
G4F3ModelFun[[11]] <- function(K1,logL,df11,df21,df31,df41,G4F3text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF3 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF3)))
  sigmaf3<-var(dataF3);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a9<-sqrt(sigmaf3/(m_sam[4]*(m_sam[4]-1)))
  m_nf <- as.numeric(G4F3text2)
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF3)))
  #########procedure start#############################################
  d4<-1;mix_pi<-1;
  m[c(1:4)]<-m
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaf3))
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dnorm(dataF3,m[4],sqrt(sigma[4]))))
  ###########################iteration process#########################################
  iteration <- 0; stopa <- 1000
  swx <- matrix(0,1,1); s0<-matrix(0,2,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    m[c(1:4)]<-m
    ########################obtain variance##################################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    swx[1]<-sum((dataF3-m[4])^2);s0[1]<-ss1+ss2+ss3;s0[2]<-m_sam[1]+m_sam[2]+m_sam[3]
    sigma[4]<-swx[1]/m_sam[4];

    sigma40<-sigma[4]-sigma[1]/m_nf

    if (sigma40<0){sigma40<-0;sigma[4]<-sigma[1]/m_nf}
    a2<-sigma[1];n_iter<-0;a1<-1000
    while(a1>0.0001){
      n_iter<-n_iter+1
      a3<-(sigma[1]/m_nf)/(sigma[1]/m_nf+sigma40)
      sigma[1]<-(s0[1]+m_nf*a3^2*swx[1])/(s0[2]+a3*m_sam[4])
      a1<-abs(sigma[1]-a2)
      a2<-sigma[1]
      if(n_iter>20) break
    }
    sigma[2]<-sigma[3]<-sigma[1]
    sigma[4]<-sigma[1]/m_nf+sigma40
    ###################criteria for iterations to stop#######################################
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+ sum(log(dnorm(dataF3,m[4],sqrt(sigma[4]))))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*6
  ##############first order genetic parameter process##################################
  gg<-sigmaf3 - sigma[1]/m_nf
  if(gg<0){gg<-0}
  rr<-gg/sigmaf3
  ####################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma[1]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < m_sam[1]){P1bmw <- P1bmw+runif(m_sam[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:m_sam[1])) - 0.5)/m_sam[1])^2)
  P1u<- as.matrix(c(12*m_sam[1]*((P1dd[1]/m_sam[1]-0.5)^2),((45*m_sam[1])/4)*((P1dd[2]/m_sam[1]-1/3)^2),180*m_sam[1]*((P1dd[3]/m_sam[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,m_sam[1]))))

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*m_sam[2])
  F1bmw <- matrix(0,m_sam[2],1)
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma[2]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < m_sam[2]){F1bmw <- F1bmw+runif(m_sam[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:m_sam[2])) - 0.5)/m_sam[2])^2)
  F1u<- as.matrix(c(12*m_sam[2]*((F1dd[1]/m_sam[2]-0.5)^2),((45*m_sam[2])/4)*((F1dd[2]/m_sam[2]-1/3)^2),180*m_sam[2]*((F1dd[3]/m_sam[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,m_sam[2]))))

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*m_sam[3])
  P2bmw <- matrix(0,m_sam[3],1)
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma[3]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < m_sam[3]){P2bmw <- P2bmw+runif(m_sam[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:m_sam[3])) - 0.5)/m_sam[3])^2)
  P2u<- as.matrix(c(12*m_sam[3]*((P2dd[1]/m_sam[3]-0.5)^2),((45*m_sam[3])/4)*((P2dd[2]/m_sam[3]-1/3)^2),180*m_sam[3]*((P2dd[3]/m_sam[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,m_sam[3]))))

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F3#########################################
  dataF3<-sort(dataF3)
  F3w1<-1/(12*m_sam[4])
  F3bmw <- matrix(0,m_sam[4],1)
  F3gg <- (dataF3 - m[4])/sqrt(as.vector(sigma[4]))
  F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
  F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3bmw)))[1]
  if(nn < m_sam[4]){F3bmw <- F3bmw+runif(m_sam[4])/1e4}
  ##########################################################
  F3dd<-c((sum(F3bmw)),(sum(F3bmw^2)),sum((F3bmw-0.5)^2))
  F3w<-F3w1+sum((F3bmw - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F3u<- as.matrix(c(12*m_sam[4]*((F3dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F3dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F3dd[3]/m_sam[4]-1/12)^2)))
  F3D<-as.numeric(ks.test(F3bmw,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3w),(1-pkolm(F3D,m_sam[4]))))

  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("PG-ADI",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," "," "," "," ",round(sigma[1],4),round(sigma[4],4)," "," "," "," "," "," "," "," ",
                       round(t(mix_pi),4)," "," "," "," "," "," "," "," ",round(m[1],4),round(m[2],4),round(m[3],4),round(m[4],4)," "," "," "," "," "," "," "," "," "," "," "," ",round(gg,4),round(rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3w,4),F3tt[4],round(F3D,4),F3tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output)
  return(OUTPUT)
}
####################################PG-AD(C-1)####################################
G4F3ModelFun[[12]] <- function(K1,logL,df11,df21,df31,df41,G4F3text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF3 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF3)))
  sigmaf3<-var(dataF3);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a9<-sqrt(sigmaf3/(m_sam[4]*(m_sam[4]-1)))
  m_nf <- as.numeric(G4F3text2)
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF3)))
  ##############procedure start#########################
  d4<-1;mix_pi<-1;m[c(1:4)]<-m
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaf3))
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dnorm(dataF3,m[4],sqrt(sigma[4]))))
  ####################iteration process######################################################
  iteration <- 0; stopa <- 1000
  swx <- matrix(0,1,1); s0<-matrix(0,2,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ####################restrictions############################################
    aa1<-3*sumx[1]/m_sam[1]+2*sumx[2]/m_sam[2]+3*sumx[3]/m_sam[3]-8*sumx[4]/m_sam[4];
    aa2<-9*sigma[1]/m_sam[1]+4*sigma[1]/m_sam[2]+9*sigma[1]/m_sam[3]+64*sigma[4]/m_sam[4];
    rr1<-aa1/aa2;
    m[1]<-(sumx[1]-3*sigma[1]*rr1)/m_sam[1];m[2]<-(sumx[2]-2*sigma[1]*rr1)/m_sam[2];
    m[3]<-(sumx[3]-3*sigma[1]*rr1)/m_sam[3];m[4]<-(sumx[4]+8*sigma[4]*rr1)/m_sam[4];
    ###################obtain variance####################################################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    swx[1]<-sum((dataF3-m[4])^2)
    s0[1]<-ss1+ss2+ss3;s0[2]<-m_sam[1]+m_sam[2]+m_sam[3]
    sigma[4]<-swx[1]/m_sam[4];

    sigma40<-sigma[4]-sigma[1]/m_nf

    if (sigma40<0){sigma40<-0;sigma[4]<-sigma[1]/m_nf}
    a2<-sigma[1];n_iter<-0;a1<-1000
    while(a1>0.0001){
      n_iter<-n_iter+1
      a3<-(sigma[1]/m_nf)/(sigma[1]/m_nf+sigma40)
      sigma[1]<-(s0[1]+m_nf*a3^2*swx[1])/(s0[2]+a3*m_sam[4])
      a1<-abs(sigma[1]-a2)
      a2<-sigma[1]
      if(n_iter>20) break
    }
    sigma[2]<-sigma[3]<-sigma[1]
    sigma[4]<-sigma[1]/m_nf+sigma40
    #####################criteria for iterations to stop##########################
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+ sum(log(dnorm(dataF3,m[4],sqrt(sigma[4]))))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*5
  ###################first order genetic parameter process##################################
  hh11 <- matrix(c(1,1,1,1,1,0,-1,0,0,1,0,0.25),4,3)
  B11<-solve(crossprod(hh11,hh11))%*%crossprod(hh11,m)
  gg<-sigmaf3 - sigma[1]/m_nf
  if(gg<0){gg<-0}
  rr<-gg/sigmaf3
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma[1]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < m_sam[1]){P1bmw <- P1bmw+runif(m_sam[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:m_sam[1])) - 0.5)/m_sam[1])^2)
  P1u<- as.matrix(c(12*m_sam[1]*((P1dd[1]/m_sam[1]-0.5)^2),((45*m_sam[1])/4)*((P1dd[2]/m_sam[1]-1/3)^2),180*m_sam[1]*((P1dd[3]/m_sam[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,m_sam[1]))))

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*m_sam[2])
  F1bmw <- matrix(0,m_sam[2],1)
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma[2]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < m_sam[2]){F1bmw <- F1bmw+runif(m_sam[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:m_sam[2])) - 0.5)/m_sam[2])^2)
  F1u<- as.matrix(c(12*m_sam[2]*((F1dd[1]/m_sam[2]-0.5)^2),((45*m_sam[2])/4)*((F1dd[2]/m_sam[2]-1/3)^2),180*m_sam[2]*((F1dd[3]/m_sam[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,m_sam[2]))))

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*m_sam[3])
  P2bmw <- matrix(0,m_sam[3],1)
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma[3]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < m_sam[3]){P2bmw <- P2bmw+runif(m_sam[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:m_sam[3])) - 0.5)/m_sam[3])^2)
  P2u<- as.matrix(c(12*m_sam[3]*((P2dd[1]/m_sam[3]-0.5)^2),((45*m_sam[3])/4)*((P2dd[2]/m_sam[3]-1/3)^2),180*m_sam[3]*((P2dd[3]/m_sam[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,m_sam[3]))))

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F3#########################################
  dataF3<-sort(dataF3)
  F3w1<-1/(12*m_sam[4])
  F3bmw <- matrix(0,m_sam[4],1)
  F3gg <- (dataF3 - m[4])/sqrt(as.vector(sigma[4]))
  F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
  F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3bmw)))[1]
  if(nn < m_sam[4]){F3bmw <- F3bmw+runif(m_sam[4])/1e4}
  ##########################################################
  F3dd<-c((sum(F3bmw)),(sum(F3bmw^2)),sum((F3bmw-0.5)^2))
  F3w<-F3w1+sum((F3bmw - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F3u<- as.matrix(c(12*m_sam[4]*((F3dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F3dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F3dd[3]/m_sam[4]-1/12)^2)))
  F3D<-as.numeric(ks.test(F3bmw,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3w),(1-pkolm(F3D,m_sam[4]))))

  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("PG-AD",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," "," "," "," ",round(sigma[1],4),round(sigma[4],4)," "," "," "," "," "," "," "," ",
                       round(t(mix_pi),4)," "," "," "," "," "," "," "," ",round(B11[1],4)," "," "," "," "," "," "," "," "," "," "," ",round(B11[2],4),round(B11[3],4)," "," ",round(gg,4),round(rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3w,4),F3tt[4],round(F3D,4),F3tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output)
  return(OUTPUT)
}
###########################MX1-AD-ADI(D-0)###############################################
G4F3ModelFun[[13]] <- function(K1,logL,df11,df21,df31,df41,G4F3text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF3 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF3)))
  sigmaf3<-var(dataF3); sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a8<-sqrt(sigmaf3)
  m_nf <- as.numeric(G4F3text2)
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF3)))
  ##########procedure start################################
  d1<-3
  mi<-as.matrix(c(0.25,0.5,0.25))
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaf3/2,1.5*sigmaf3/2,sigmaf3/2))
  m[c(1:4)]<-m
  m[5]<-m[4];m[4]<-m[5]+2*a8;m[6]<-m[5]-2*a8
  pi<-m[c(4,5,6)] ;gh<-sigma[c(4,5,6)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mi)))
  ##################################iteration process########################################################
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,3,m_sam[4]); swx <- matrix(0,3,1)
  g_a<-matrix(0,1,1);s0<-matrix(0,2,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:6)];ssigma<-sigma[c(4:6)]
    for(i in 1:d1) { WW[i,] <- mi[i]*dnorm(dataF3,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF3,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF3
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    m[c(1:3)]<-sumx[c(1:3)]/m_sam[c(1:3)];m[c(4:6)]<-sumwx[c(1:3)]/n0[c(1:3)]
    ###################first order genetic parameter process#########################
    hh12 <- matrix(c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,1,1,1,0,-1,1,0,-1,0,1,0,0,0.5,0),6,6)
    B12 <- solve(hh12,m)
    g_a[1]<-(0.5*B12[5]^2+0.25*B12[6]^2)/m_nf
    #####################obtain variance############################################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d1) {swx[i] <- WW[i,]%*%(dataF3-m[i+3])^2 }
    s0[1]<-swx[1]+swx[3];s0[2]<-n0[1]+n0[3];
    a2<-sigma[4];n_iter<-0;a1<-1000
    while(a1>0.0001){
      n_iter<-n_iter+1
      a3<-sigma[4]/(sigma[4]+g_a[1])
      sigma[4]<-(s0[1]+a3^2*swx[2])/(s0[2]+a3*n0[2])
      a1<-abs(sigma[4]-a2)
      a2<-sigma[4]
      if(n_iter>20) break
    }

    sigma40<-sigma[4]-sigma[1]/m_nf
    if (sigma40<0){sigma40<-0;sigma[4]<-sigma[1]/m_nf}
    sigma[6]<-sigma[4];sigma[5]<-sigma[4]+g_a[1]
    a1<-sigma[1];a2<-ss1+ss2+ss3;a3<-m_sam[1]+m_sam[2]+m_sam[3]
    n_iter<-0;a6<-1000
    while(a6>0.0001){
      n_iter<-n_iter+1
      a4<-(sigma[1]/m_nf)/(sigma[1]/m_nf+sigma40)
      a5<-(sigma[1]/m_nf)/(sigma[1]/m_nf+sigma40+g_a[1])
      a4[a4>1]<-1;a5[a5>1]<-1
      sigma[1]<-(a2+a4^2*m_nf*(swx[1]+swx[3])+a5^2*m_nf*swx[2])/(a3+a4*(n0[1]+n0[3])+a5*n0[2])
      a6<-abs(sigma[1]-a1)
      a1<-sigma[1]
      if(n_iter>20) break
    }
    sigma[2]<-sigma[3]<-sigma[1];sigma[4]<-sigma[1]/m_nf+sigma40
    sigma[6]<-sigma[4];sigma[5]<-sigma[4]+g_a[1]
    #######################criteria for iteration to stop##########################
    pi<-m[c(4,5,6)];gh<-sigma[c(4,5,6)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*7
  ##########second order genetic parameter process###############################################
  jj<-sigmaf3-sigma[1]/m_nf
  if(jj<0){jj<-0}
  ll<-jj/sigmaf3
  jj<-sigmaf3 - sigma[4];
  gg<-sigma[4]-sigma[1]/m_nf;
  jj[jj<0]<-0
  gg[gg<0]<-0
  ll<-jj/sigmaf3
  rr<-gg/sigmaf3
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma[1]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < m_sam[1]){P1bmw <- P1bmw+runif(m_sam[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:m_sam[1])) - 0.5)/m_sam[1])^2)
  P1u<- as.matrix(c(12*m_sam[1]*((P1dd[1]/m_sam[1]-0.5)^2),((45*m_sam[1])/4)*((P1dd[2]/m_sam[1]-1/3)^2),180*m_sam[1]*((P1dd[3]/m_sam[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,m_sam[1]))))

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*m_sam[2])
  F1bmw <- matrix(0,m_sam[2],1)
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma[2]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < m_sam[2]){F1bmw <- F1bmw+runif(m_sam[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:m_sam[2])) - 0.5)/m_sam[2])^2)
  F1u<- as.matrix(c(12*m_sam[2]*((F1dd[1]/m_sam[2]-0.5)^2),((45*m_sam[2])/4)*((F1dd[2]/m_sam[2]-1/3)^2),180*m_sam[2]*((F1dd[3]/m_sam[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,m_sam[2]))))

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*m_sam[3])
  P2bmw <- matrix(0,m_sam[3],1)
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma[3]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < m_sam[3]){P2bmw <- P2bmw+runif(m_sam[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:m_sam[3])) - 0.5)/m_sam[3])^2)
  P2u<- as.matrix(c(12*m_sam[3]*((P2dd[1]/m_sam[3]-0.5)^2),((45*m_sam[3])/4)*((P2dd[2]/m_sam[3]-1/3)^2),180*m_sam[3]*((P2dd[3]/m_sam[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,m_sam[3]))))

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3);
  F3w1<-1/(12*m_sam[4])
  F3bmw <- matrix(0,m_sam[4],1); F3bmwsl <- matrix(0,m_sam[4],d1)
  for(i in 1:d1){
    F3gg <- (dataF3 - pi[i])/sqrt(gh[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mix_pi[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[4]){F3P2 <- F3P2+runif(m_sam[4])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[4]) + sum((F3P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F3u <- as.matrix(c(12*m_sam[4]*((F3dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F3dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F3dd[3]/m_sam[4]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[4]))))

  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX1-AD-ADI",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," "," ",round(sigma[1],4),round(sigma[4],4),round(sigma[5],4),round(sigma[6],4)," "," "," "," "," "," ",
                       round(t(mix_pi),4)," "," "," "," "," "," ",round(B12[1],4),round(B12[2],4),round(B12[3],4),round(B12[4],4),round(B12[5],4)," ",round(B12[6],4)," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4),round(gg,4),round(rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
#####################################MX1-AD-AD(D-1)##############################################################
G4F3ModelFun[[14]] <- function(K1,logL,df11,df21,df31,df41,G4F3text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF3 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF3)))
  sigmaf3<-var(dataF3);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a8<-sqrt(sigmaf3)
  m_nf <- as.numeric(G4F3text2)
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF3)))
  ############procedure start##########################################
  d1<-3
  mi<-as.matrix(c(0.25,0.5,0.25))
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaf3/2,1.5*sigmaf3/2,sigmaf3/2))
  m[c(1:4)]<-m
  m[5]<-m[4];m[4]<-m[5]+2*a8;m[6]<-m[5]-2*a8
  pi<-m[c(4,5,6)] ;gh<-sigma[c(4,5,6)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mi)))
  #######################iteration process########################
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,3,m_sam[4]); swx <- matrix(0,3,1)
  g_a<-matrix(0,1,1);s0<-matrix(0,2,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:6)];ssigma<-sigma[c(4:6)]
    for(i in 1:d1) { WW[i,] <- mi[i]*dnorm(dataF3,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF3,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF3
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    ##############restrictions#################################
    aa1<-3*sumx[1]/m_sam[1]+2*sumx[2]/m_sam[2]+3*sumx[3]/m_sam[3]-2*sumwx[1]/n0[1]-4*sumwx[2]/n0[2]-2*sumwx[3]/n0[3]
    aa2<-9*sigma[1]/m_sam[1]+4*sigma[2]/m_sam[2]+9*sigma[3]/m_sam[3]+4*sigma[4]/n0[1]+16*sigma[5]/n0[2]+4*sigma[6]/n0[3]
    aa3<-aa1/aa2
    m[1]<-(sumx[1]-sigma[1]*aa3*3)/m_sam[1];m[2]<-(sumx[2]-sigma[2]*aa3*2)/m_sam[2];m[3]<-(sumx[3]-sigma[3]*aa3*3)/m_sam[3];
    m[4]<-(sumwx[1]+sigma[4]*aa3*2)/n0[1];m[5]<-(sumwx[2]+sigma[5]*aa3*4)/n0[2];m[6]<-(sumwx[3]+sigma[6]*aa3*2)/n0[3];
    ########first order genetic parameter process############
    hh13<-matrix(c(1,1,1,1,1,1,1,0,-1,1,0,-1,0,1,0,0,0.5,0,1,0,-1,0,0,0,0,1,0,0.25,0.25,0.25),6,5)
    B13<-solve(crossprod(hh13,hh13))%*%crossprod(hh13,m)
    g_a[1]<-(0.5*B13[2]^2+0.25*B13[3]^2)/m_nf
    ##########obtain variance######################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d1) {swx[i] <- WW[i,]%*%(dataF3-m[i+3])^2 }
    s0[1]<-swx[1]+swx[3];s0[2]<-n0[1]+n0[3]
    a2<-sigma[4];n_iter<-0;a1<-1000
    while(a1>0.0001){
      n_iter<-n_iter+1
      a3<-sigma[4]/(sigma[4]+g_a[1])
      sigma[4]<-(s0[1]+a3^2*swx[2])/(s0[2]+a3*n0[2])
      a1<-abs(sigma[4]-a2)
      a2<-sigma[4]
      if(n_iter>20) break
    }

    sigma40<-sigma[4]-sigma[1]/m_nf;

    if (sigma40<0){sigma40<-0;sigma[4]<-sigma[1]/m_nf}
    sigma[6]<-sigma[4];sigma[5]<-sigma[4]+g_a[1]
    a1<-sigma[1];a2<-ss1+ss2+ss3;a3<-m_sam[1]+m_sam[2]+m_sam[3];
    n_iter<-0;a6<-1000
    while(a6>0.0001){
      n_iter<-n_iter+1
      a4<-(sigma[1]/m_nf)/(sigma[1]/m_nf+sigma40);
      a5<-(sigma[1]/m_nf)/(sigma[1]/m_nf+sigma40+g_a[1]);
      a4[a4>1]<-1;a5[a5>1]<-1
      sigma[1]<-(a2+a4^2*m_nf*(swx[1]+swx[3])+a5^2*m_nf*swx[2])/(a3+a4*(n0[1]+n0[3])+a5*n0[2])
      a6<-abs(sigma[1]-a1)
      a1<-sigma[1]
      if(n_iter>20) break
    }
    sigma[2]<-sigma[3]<-sigma[1];sigma[4]<-sigma[1]/m_nf+sigma40
    sigma[6]<-sigma[4];sigma[5]<-sigma[4]+g_a[1]
    ###################criteria for iterations to stop#####################
    pi<-m[c(4,5,6)] ;gh<-sigma[c(4,5,6)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*7
  ###################second order genetic parameter process#################
  jj<-sigmaf3-sigma[1]/m_nf
  if(jj<0){jj<-0}
  ll<-jj/sigmaf3
  jj<-sigmaf3 - sigma[4];
  gg<-sigma[4]-sigma[1]/m_nf;
  jj[jj<0]<-0
  gg[gg<0]<-0
  ll<-jj/sigmaf3
  rr<-gg/sigmaf3
  ######################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma[1]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < m_sam[1]){P1bmw <- P1bmw+runif(m_sam[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:m_sam[1])) - 0.5)/m_sam[1])^2)
  P1u<- as.matrix(c(12*m_sam[1]*((P1dd[1]/m_sam[1]-0.5)^2),((45*m_sam[1])/4)*((P1dd[2]/m_sam[1]-1/3)^2),180*m_sam[1]*((P1dd[3]/m_sam[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,m_sam[1]))))

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*m_sam[2])
  F1bmw <- matrix(0,m_sam[2],1)
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma[2]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < m_sam[2]){F1bmw <- F1bmw+runif(m_sam[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:m_sam[2])) - 0.5)/m_sam[2])^2)
  F1u<- as.matrix(c(12*m_sam[2]*((F1dd[1]/m_sam[2]-0.5)^2),((45*m_sam[2])/4)*((F1dd[2]/m_sam[2]-1/3)^2),180*m_sam[2]*((F1dd[3]/m_sam[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,m_sam[2]))))

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*m_sam[3])
  P2bmw <- matrix(0,m_sam[3],1)
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma[3]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < m_sam[3]){P2bmw <- P2bmw+runif(m_sam[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:m_sam[3])) - 0.5)/m_sam[3])^2)
  P2u<- as.matrix(c(12*m_sam[3]*((P2dd[1]/m_sam[3]-0.5)^2),((45*m_sam[3])/4)*((P2dd[2]/m_sam[3]-1/3)^2),180*m_sam[3]*((P2dd[3]/m_sam[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,m_sam[3]))))

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3);
  F3w1<-1/(12*m_sam[4])
  F3bmw <- matrix(0,m_sam[4],1); F3bmwsl <- matrix(0,m_sam[4],d1)
  for(i in 1:d1){
    F3gg <- (dataF3 - pi[i])/sqrt(gh[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mix_pi[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[4]){F3P2 <- F3P2+runif(m_sam[4])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[4]) + sum((F3P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F3u <- as.matrix(c(12*m_sam[4]*((F3dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F3dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F3dd[3]/m_sam[4]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[4]))))

  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX1-AD-AD",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," "," ",round(sigma[1],4),round(sigma[4],4),round(sigma[5],4),round(sigma[6],4)," "," "," "," "," "," ",
                       round(t(mix_pi),4)," "," "," "," "," "," ",round(B13[1],4)," "," "," ",round(B13[2],4)," ",round(B13[3],4)," "," "," "," "," ",round(B13[4],4),round(B13[5],4),round(jj,4),round(ll*100,4),round(gg,4),round(rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
#################################MX1-A-AD(D-2)#######################################################
G4F3ModelFun[[15]] <- function(K1,logL,df11,df21,df31,df41,G4F3text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF3 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF3)))
  sigmaf3<-var(dataF3);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a8<-sqrt(sigmaf3)
  m_nf <- as.numeric(G4F3text2)
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF3)))
  ############procedure start##########################################
  d1<-3
  mi<-as.matrix(c(0.25,0.5,0.25))
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaf3/2,1.5*sigmaf3/2,sigmaf3/2))
  m[c(1:4)]<-m
  m[5]<-m[4];m[4]<-m[5]+2*a8;m[6]<-m[5]-2*a8
  pi<-m[c(4,5,6)];gh<-sigma[c(4,5,6)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mi)))
  ####################################################################
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,3,m_sam[4]); swx <- matrix(0,3,1)
  rr<-matrix(0,2,1); g_a<-matrix(0,1,1); s0<-matrix(0,2,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:6)];ssigma<-sigma[c(4:6)]
    for(i in 1:d1) { WW[i,] <- mi[i]*dnorm(dataF3,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF3,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF3
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    ##################restrictions############################################
    aa1<-3*sumx[1]/m_sam[1]+2*sumx[2]/m_sam[2]+3*sumx[3]/m_sam[3]-8*sumwx[2]/n0[2]
    aa2<-sumwx[1]/n0[1]-2*sumwx[2]/n0[2]+sumwx[3]/n0[3]
    aa3<-9*sigma[1]/m_sam[1]+4*sigma[2]/m_sam[2]+9*sigma[3]/m_sam[3]+64*sigma[5]/n0[2]
    aa4<-16*sigma[5]/n0[2]
    aa5<-sigma[4]/n0[1]+4*sigma[5]/n0[2]+sigma[6]/n0[3]
    aa6<-aa3*aa5-aa4^2;aa7<-aa1*aa5-aa2*aa4;aa8<-aa2*aa3-aa1*aa4
    rr[1]<-aa7/aa6;rr[2]<-aa8/aa6
    m[1]<-(sumx[1]-sigma[1]*rr[1]*3)/m_sam[1];
    m[2]<-(sumx[2]-sigma[2]*rr[1]*2)/m_sam[2];
    m[3]<-(sumx[3]-sigma[3]*rr[1]*3)/m_sam[3];
    m[4]<-(sumwx[1]-sigma[4]*rr[2])/n0[1];
    m[5]<-(sumwx[2]+sigma[5]*(8*rr[1]+2*rr[2]))/n0[2];
    m[6]<-(sumwx[3]-sigma[6]*rr[2])/n0[3];
    #################first order genetic parameter process#####################################
    hh14<-matrix(c(1,1,1,1,1,1,1,0,-1,1,0,-1,1,0,-1,0,0,0,0,1,0,0.25,0.25,0.25),6,4)
    B14<-solve(crossprod(hh14,hh14))%*%crossprod(hh14,m)
    g_a[1]<-0.5*B14[2]^2/m_nf;
    ##############obtain variance###############################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d1) {swx[i] <- WW[i,]%*%(dataF3-m[i+3])^2 }
    s0[1]<-swx[1]+swx[3];s0[2]<-n0[1]+n0[3];a2<-sigma[4]
    n_iter<-0;a1<-1000
    while(a1>0.0001){
      n_iter<-n_iter+1
      a3<-sigma[4]/(sigma[4]+g_a[1])
      sigma[4]<-(s0[1]+a3^2*swx[2])/(s0[2]+a3*n0[2])
      a1<-abs(sigma[4]-a2)
      a2<-sigma[4]
      if(n_iter>20) break
    }

    sigma40<-sigma[4]-sigma[1]/m_nf

    if (sigma40<0){sigma40<-0;sigma[4]<-sigma[1]/m_nf}
    sigma[6]<-sigma[4];sigma[5]<-sigma[4]+g_a[1]
    a1<-sigma[1];a2<-ss1+ss2+ss3;a3<-m_sam[1]+m_sam[2]+m_sam[3]
    n_iter<-0;a6<-1000
    while(a6>0.0001){
      n_iter<-n_iter+1
      a4<-(sigma[1]/m_nf)/(sigma[1]/m_nf+sigma40)
      a5<-(sigma[1]/m_nf)/(sigma[1]/m_nf+sigma40+g_a[1])
      a4[a4>1]<-1;a5[a5>1]<-1
      sigma[1]<-(a2+a4^2*m_nf*(swx[1]+swx[3])+a5^2*m_nf*swx[2])/(a3+a4*(n0[1]+n0[3])+a5*n0[2])
      a6<-abs(sigma[1]-a1)
      a1<-sigma[1]
      if(n_iter>20) break
    }
    sigma[2]<-sigma[3]<-sigma[1]
    sigma[4]<-sigma[1]/m_nf+sigma40;sigma[6]<-sigma[4];sigma[5]<-sigma[4]+g_a[1]
    #################criteria for iterations to stop#####################
    pi<-m[c(4,5,6)] ;gh<-sigma[c(4,5,6)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*5
  ###########second order genetic parameter process#######################
  jj<-sigmaf3-sigma[1]/m_nf
  if(jj<0){jj<-0}
  ll<-jj/sigmaf3
  jj<-sigmaf3 - sigma[4];
  gg<-sigma[4]-sigma[1]/m_nf;
  jj[jj<0]<-0
  gg[gg<0]<-0
  ll<-jj/sigmaf3
  rr<-gg/sigmaf3
  ####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma[1]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < m_sam[1]){P1bmw <- P1bmw+runif(m_sam[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:m_sam[1])) - 0.5)/m_sam[1])^2)
  P1u<- as.matrix(c(12*m_sam[1]*((P1dd[1]/m_sam[1]-0.5)^2),((45*m_sam[1])/4)*((P1dd[2]/m_sam[1]-1/3)^2),180*m_sam[1]*((P1dd[3]/m_sam[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,m_sam[1]))))

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*m_sam[2])
  F1bmw <- matrix(0,m_sam[2],1)
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma[2]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < m_sam[2]){F1bmw <- F1bmw+runif(m_sam[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:m_sam[2])) - 0.5)/m_sam[2])^2)
  F1u<- as.matrix(c(12*m_sam[2]*((F1dd[1]/m_sam[2]-0.5)^2),((45*m_sam[2])/4)*((F1dd[2]/m_sam[2]-1/3)^2),180*m_sam[2]*((F1dd[3]/m_sam[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,m_sam[2]))))

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*m_sam[3])
  P2bmw <- matrix(0,m_sam[3],1)
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma[3]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < m_sam[3]){P2bmw <- P2bmw+runif(m_sam[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:m_sam[3])) - 0.5)/m_sam[3])^2)
  P2u<- as.matrix(c(12*m_sam[3]*((P2dd[1]/m_sam[3]-0.5)^2),((45*m_sam[3])/4)*((P2dd[2]/m_sam[3]-1/3)^2),180*m_sam[3]*((P2dd[3]/m_sam[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,m_sam[3]))))

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3);
  F3w1<-1/(12*m_sam[4])
  F3bmw <- matrix(0,m_sam[4],1); F3bmwsl <- matrix(0,m_sam[4],d1)
  for(i in 1:d1){
    F3gg <- (dataF3 - pi[i])/sqrt(gh[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mix_pi[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[4]){F3P2 <- F3P2+runif(m_sam[4])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[4]) + sum((F3P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F3u <- as.matrix(c(12*m_sam[4]*((F3dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F3dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F3dd[3]/m_sam[4]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[4]))))

  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX1-A-AD",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," "," ",round(sigma[1],4),round(sigma[4],4),round(sigma[5],4),round(sigma[6],4)," "," "," "," "," "," ",
                       round(t(mix_pi),4)," "," "," "," "," "," ",round(B14[1],4)," "," "," ",round(B14[2],4)," "," "," "," "," "," "," ",round(B14[3],4),round(B14[4],4),round(jj,4),round(ll*100,4),round(gg,4),round(rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
##############################MX1-EAD-AD(D-3)#############################################################
G4F3ModelFun[[16]] <- function(K1,logL,df11,df21,df31,df41,G4F3text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF3 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF3)))
  sigmaf3<-var(dataF3);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a8<-sqrt(sigmaf3/(m_sam[4]))
  m_nf <- as.numeric(G4F3text2)
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF3)))
  ##################procedure start#############################################
  d1<-3
  mi<-as.matrix(c(0.25,0.5,0.25))
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaf3/2,1.5*sigmaf3/2,sigmaf3/2))
  m[c(1:4)]<-m
  m[5]<-m[4]+1.5*a8;m[6]<-m[4]-2*a8;m[4]<-m[4]+2*a8
  pi<-m[c(4,5,6)] ;gh<-sigma[c(4,5,6)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+ sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mi)))
  ################iteration process########################
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,3,m_sam[4]); swx <- matrix(0,3,1)
  rr<-matrix(0,2,1);g_a<-matrix(0,1,1); s0<-matrix(0,2,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:6)];ssigma<-sigma[c(4:6)]
    for(i in 1:d1) { WW[i,] <- mi[i]*dnorm(dataF3,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF3,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF3
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    ##################restrictions###################
    aa1<-3*sumx[1]/m_sam[1]+2*sumx[2]/m_sam[2]+3*sumx[3]/m_sam[3]-2*sumwx[1]/n0[1]-4*sumwx[2]/n0[2]-2*sumwx[3]/n0[3]
    aa2<-3*sumwx[1]/n0[1]-4*sumwx[2]/n0[2]+sumwx[3]/n0[3]
    aa3<-9*sigma[1]/m_sam[1]+4*sigma[2]/m_sam[2]+9*sigma[3]/m_sam[3]+4*sigma[4]/n0[1]+16*sigma[5]/n0[2]+4*sigma[6]/n0[3]
    aa4<--6*sigma[4]/n0[1]+16*sigma[5]/n0[2]-2*sigma[6]/n0[3]
    aa5<-9*sigma[4]/n0[1]+16*sigma[5]/n0[2]+sigma[6]/n0[3]
    aa6<-aa3*aa5-aa4*aa4; aa7<-aa1*aa5-aa2*aa4; aa8<-aa2*aa3-aa1*aa4
    rr[1]<-aa7/aa6;rr[2]<-aa8/aa6
    m[1]<-(sumx[1]-sigma[1]*rr[1]*3)/m_sam[1]
    m[2]<-(sumx[2]-sigma[2]*rr[1]*2)/m_sam[2]
    m[3]<-(sumx[3]-sigma[3]*rr[1]*3)/m_sam[3]
    m[4]<-(sumwx[1]+sigma[4]*(2*rr[1]-3*rr[2]))/n0[1]
    m[5]<-(sumwx[2]+sigma[5]*(4*rr[1]+4*rr[2]))/n0[2]
    m[6]<-(sumwx[3]+sigma[6]*(2*rr[1]-rr[2]))/n0[3]
    ###########first order genetic parameter process##############
    hh15<-matrix(c(1,1,1,1,1,1,1,1,-1,1,0.5,-1,1,0,-1,0,0,0,0,1,0,0.25,0.25,0.25),6,4)
    B15<-solve(crossprod(hh15,hh15))%*%crossprod(hh15,m)
    g_a[1]<-0.75*B15[2]^2/m_nf
    ############obtain variance##########################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d1) {swx[i] <- WW[i,]%*%(dataF3-m[i+3])^2 }
    s0[1]<-swx[1]+swx[3];s0[2]<-n0[1]+n0[3];a2<-sigma[4];n_iter<-0;a1<-1000
    while(a1>0.0001){
      n_iter<-n_iter+1
      a3<-sigma[4]/(sigma[4]+g_a[1])
      sigma[4]<-(s0[1]+a3^2*swx[2])/(s0[2]+a3*n0[2])
      a1<-abs(sigma[4]-a2)
      a2<-sigma[4]
      if(n_iter>20) break
    }

    sigma40<-sigma[4]-sigma[1]/m_nf

    if (sigma40<0){sigma40<-0;sigma[4]<-sigma[1]/m_nf}
    sigma[6]<-sigma[4]
    sigma[5]<-sigma[4]+g_a[1]
    a1<-sigma[1]; a2<-ss1+ss2+ss3;a3<-m_sam[1]+m_sam[2]+m_sam[3]
    n_iter<-0;a6<-1000
    while(a6>0.0001){
      n_iter<-n_iter+1
      a4<-(sigma[1]/m_nf)/(sigma[1]/m_nf+sigma40)
      a5<-(sigma[1]/m_nf)/(sigma[1]/m_nf+sigma40+g_a[1])
      a4[a4>1]<-1;a5[a5>1]<-1
      sigma[1]<-(a2+a4^2*m_nf*(swx[1]+swx[3])+a5^2*m_nf*swx[2])/(a3+a4*(n0[1]+n0[3])+a5*n0[2])
      a6<-abs(sigma[1]-a1)
      a1<-sigma[1]
      if(n_iter>20) break
    }
    sigma[2]<-sigma[3]<-sigma[1]
    sigma[4]<-sigma[1]/m_nf+sigma40;sigma[6]<-sigma[4];sigma[5]<-sigma[4]+g_a[1]
    ##############criteria for iteration to stop###############################
    pi<-m[c(4,5,6)] ;gh<-sigma[c(4,5,6)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*5
  ################second order genetic parameter process##########################
  jj<-sigmaf3 - sigma[4]
  gg<-sigma[4]-sigma[1]/m_nf
  jj[jj<0]<-0
  gg[gg<0]<-0
  ll<-jj/sigmaf3
  rr<-gg/sigmaf3
  #############hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma[1]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < m_sam[1]){P1bmw <- P1bmw+runif(m_sam[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:m_sam[1])) - 0.5)/m_sam[1])^2)
  P1u<- as.matrix(c(12*m_sam[1]*((P1dd[1]/m_sam[1]-0.5)^2),((45*m_sam[1])/4)*((P1dd[2]/m_sam[1]-1/3)^2),180*m_sam[1]*((P1dd[3]/m_sam[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,m_sam[1]))))

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*m_sam[2])
  F1bmw <- matrix(0,m_sam[2],1)
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma[2]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < m_sam[2]){F1bmw <- F1bmw+runif(m_sam[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:m_sam[2])) - 0.5)/m_sam[2])^2)
  F1u<- as.matrix(c(12*m_sam[2]*((F1dd[1]/m_sam[2]-0.5)^2),((45*m_sam[2])/4)*((F1dd[2]/m_sam[2]-1/3)^2),180*m_sam[2]*((F1dd[3]/m_sam[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,m_sam[2]))))

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*m_sam[3])
  P2bmw <- matrix(0,m_sam[3],1)
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma[3]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < m_sam[3]){P2bmw <- P2bmw+runif(m_sam[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:m_sam[3])) - 0.5)/m_sam[3])^2)
  P2u<- as.matrix(c(12*m_sam[3]*((P2dd[1]/m_sam[3]-0.5)^2),((45*m_sam[3])/4)*((P2dd[2]/m_sam[3]-1/3)^2),180*m_sam[3]*((P2dd[3]/m_sam[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,m_sam[3]))))

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3);
  F3w1<-1/(12*m_sam[4])
  F3bmw <- matrix(0,m_sam[4],1); F3bmwsl <- matrix(0,m_sam[4],d1)
  for(i in 1:d1){
    F3gg <- (dataF3 - pi[i])/sqrt(gh[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mix_pi[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[4]){F3P2 <- F3P2+runif(m_sam[4])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[4]) + sum((F3P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F3u <- as.matrix(c(12*m_sam[4]*((F3dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F3dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F3dd[3]/m_sam[4]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[4]))))

  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX1-EAD-AD",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," "," ",round(sigma[1],4),round(sigma[4],4),round(sigma[5],4),round(sigma[6],4)," "," "," "," "," "," ",
                       round(t(mix_pi),4)," "," "," "," "," "," ",round(B15[1],4)," "," "," ",round(B15[2],4)," "," "," "," "," "," "," ",round(B15[3],4),round(B15[4],4),round(jj,4),round(ll*100,4),round(gg,4),round(rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
#############################MX1-NCD-AD(D-4)###################################################
G4F3ModelFun[[17]] <- function(K1,logL,df11,df21,df31,df41,G4F3text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF3 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF3)))
  sigmaf3<-var(dataF3);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a8<-sqrt(sigmaf3/(m_sam[4]))
  m_nf <- as.numeric(G4F3text2)
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF3)))
  ##################procedure start#############################################
  d1<-3
  mi<-as.matrix(c(0.25,0.5,0.25))
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaf3/2,1.5*sigmaf3/2,sigmaf3/2))
  m[c(1:4)]<-m
  m[5]<-m[4]-1.5*a8;m[6]<-m[4]-2*a8;m[4]<-m[4]+2*a8
  pi<-m[c(4,5,6)] ;gh<-sigma[c(4,5,6)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mi)))
  ################iteration process######################
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,3,m_sam[4]); swx <- matrix(0,3,1)
  rr<-matrix(0,2,1); g_a<-matrix(0,1,1);s0<-matrix(0,2,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:6)];ssigma<-sigma[c(4:6)]
    for(i in 1:d1) { WW[i,] <- mi[i]*dnorm(dataF3,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF3,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF3
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    ##############restrictions###################################
    aa1<-3*sumx[1]/m_sam[1]+2*sumx[2]/m_sam[2]+3*sumx[3]/m_sam[3]-2*sumwx[1]/n0[1]-4*sumwx[2]/n0[2]-2*sumwx[3]/n0[3]
    aa2<-sumwx[1]/n0[1]-4*sumwx[2]/n0[2]+3*sumwx[3]/n0[3]
    aa3<-9*sigma[1]/m_sam[1]+4*sigma[2]/m_sam[2]+9*sigma[3]/m_sam[3]+4*sigma[4]/n0[1]+16*sigma[5]/n0[2]+4*sigma[6]/n0[3]
    aa4<--2*sigma[4]/n0[1]+16*sigma[5]/n0[2]-6*sigma[6]/n0[3]
    aa5<-sigma[4]/n0[1]+16*sigma[5]/n0[2]+9*sigma[6]/n0[3]
    aa6<-aa3*aa5-aa4*aa4;aa7<-aa1*aa5-aa2*aa4;aa8<-aa2*aa3-aa1*aa4
    rr[1]<-aa7/aa6;rr[2]<-aa8/aa6;
    m[1]<-(sumx[1]-sigma[1]*rr[1]*3)/m_sam[1]
    m[2]<-(sumx[2]-sigma[2]*rr[1]*2)/m_sam[2]
    m[3]<-(sumx[3]-sigma[3]*rr[1]*3)/m_sam[3]
    m[4]<-(sumwx[1]+sigma[4]*(2*rr[1]-rr[2]))/n0[1]
    m[5]<-(sumwx[2]+sigma[5]*(4*rr[1]+4*rr[2]))/n0[2]
    m[6]<-(sumwx[3]+sigma[6]*(2*rr[1]-3*rr[2]))/n0[3]
    ########first order genetic parameter process ################
    hh16<-matrix(c(1,1,1,1,1,1,1,-1,-1,1,-0.5,-1,1,0,-1,0,0,0,0,1,0,0.25,0.25,0.25),6,4)
    B16<-solve(crossprod(hh16,hh16))%*%crossprod(hh16,m)
    g_a[1]<-0.75*B16[2]^2/m_nf
    ##########obtain variance##########################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d1) {swx[i] <- WW[i,]%*%(dataF3-m[i+3])^2 }
    s0[1]<-swx[1]+swx[3];s0[2]<-n0[1]+n0[3]
    a2<-sigma[4]; n_iter<-0;a1<-1000
    while(a1>0.0001){
      n_iter<-n_iter+1
      a3<-sigma[4]/(sigma[4]+g_a[1])
      sigma[4]<-(s0[1]+a3^2*swx[2])/(s0[2]+a3*n0[2])
      a1<-abs(sigma[4]-a2)
      a2<-sigma[4]
      if(n_iter>20) break
    }

    sigma40<-sigma[4]-sigma[1]/m_nf

    if (sigma40<0){sigma40<-0;sigma[4]<-sigma[1]/m_nf}
    sigma[6]<-sigma[4];sigma[5]<-sigma[4]+g_a[1]
    a1<-sigma[1];a2<-ss1+ss2+ss3;a3<-m_sam[1]+m_sam[2]+m_sam[3]
    n_iter<-0;a6<-1000
    while(a6>0.0001){
      n_iter<-n_iter+1
      a4<-(sigma[1]/m_nf)/(sigma[1]/m_nf+sigma40)
      a5<-(sigma[1]/m_nf)/(sigma[1]/m_nf+sigma40+g_a[1])
      a4[a4>1]<-1;a5[a5>1]<-1
      sigma[1]<-(a2+a4^2*m_nf*(swx[1]+swx[3])+a5^2*m_nf*swx[2])/(a3+a4*(n0[1]+n0[3])+a5*n0[2])
      a6<-abs(sigma[1]-a1);a1<-sigma[1]
      if(n_iter>20) break
    }
    sigma[2]<-sigma[3]<-sigma[1];sigma[4]<-sigma[1]/m_nf+sigma40
    sigma[6]<-sigma[4];sigma[5]<-sigma[4]+g_a[1]
    ##############critera for iterations to stop######################
    pi<-m[c(4,5,6)]; gh<-sigma[c(4,5,6)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*5
  ################second order genetic parameter process################
  jj<-sigmaf3 - sigma[4];
  gg<-sigma[4]-sigma[1]/m_nf;
  jj[jj<0]<-0
  gg[gg<0]<-0
  ll<-jj/sigmaf3
  rr<-gg/sigmaf3
  #############hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma[1]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < m_sam[1]){P1bmw <- P1bmw+runif(m_sam[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:m_sam[1])) - 0.5)/m_sam[1])^2)
  P1u<- as.matrix(c(12*m_sam[1]*((P1dd[1]/m_sam[1]-0.5)^2),((45*m_sam[1])/4)*((P1dd[2]/m_sam[1]-1/3)^2),180*m_sam[1]*((P1dd[3]/m_sam[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,m_sam[1]))))

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*m_sam[2])
  F1bmw <- matrix(0,m_sam[2],1)
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma[2]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < m_sam[2]){F1bmw <- F1bmw+runif(m_sam[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:m_sam[2])) - 0.5)/m_sam[2])^2)
  F1u<- as.matrix(c(12*m_sam[2]*((F1dd[1]/m_sam[2]-0.5)^2),((45*m_sam[2])/4)*((F1dd[2]/m_sam[2]-1/3)^2),180*m_sam[2]*((F1dd[3]/m_sam[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,m_sam[2]))))

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*m_sam[3])
  P2bmw <- matrix(0,m_sam[3],1)
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma[3]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < m_sam[3]){P2bmw <- P2bmw+runif(m_sam[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:m_sam[3])) - 0.5)/m_sam[3])^2)
  P2u<- as.matrix(c(12*m_sam[3]*((P2dd[1]/m_sam[3]-0.5)^2),((45*m_sam[3])/4)*((P2dd[2]/m_sam[3]-1/3)^2),180*m_sam[3]*((P2dd[3]/m_sam[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,m_sam[3]))))

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3);
  F3w1<-1/(12*m_sam[4])
  F3bmw <- matrix(0,m_sam[4],1); F3bmwsl <- matrix(0,m_sam[4],d1)
  for(i in 1:d1){
    F3gg <- (dataF3 - pi[i])/sqrt(gh[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mix_pi[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[4]){F3P2 <- F3P2+runif(m_sam[4])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[4]) + sum((F3P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F3u <- as.matrix(c(12*m_sam[4]*((F3dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F3dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F3dd[3]/m_sam[4]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[4]))))

  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX1-NCD-AD",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," "," ",round(sigma[1],4),round(sigma[4],4),round(sigma[5],4),round(sigma[6],4)," "," "," "," "," "," ",
                       round(t(mix_pi),4)," "," "," "," "," "," ",round(B16[1],4)," "," "," ",round(B16[2],4)," "," "," "," "," "," "," ",round(B16[3],4),round(B16[4],4),round(jj,4),round(ll*100,4),round(gg,4),round(rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
####################MX2-ADI-ADI(E-0)###########################################################################
G4F3ModelFun[[18]] <- function(K1,logL,df11,df21,df31,df41,G4F3text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF3 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF3)))
  sigmaf3<-var(dataF3);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a9<-sqrt(sigmaf3)
  m_nf <- as.numeric(G4F3text2)
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF3)))
  ##################procedure start#############################################
  d3<-9
  mi<-as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2))
  m[c(1:4)]<-m
  if (m[1]<m[3]) {a9<--a9}
  m[5]<-m[4]+2.4*a9;m[6]<-m[4]+1.6*a9;m[7]<-m[4]+0.8*a9
  m[8]<-m[4];m[9]<-m[4]-0.8*a9;m[10]<-m[4]-1.6*a9
  m[11]<-m[4]-2.4*a9;m[12]<-m[4]-3.2*a9;m[4]<-m[4]+3.2*a9
  pi<-m[c(4:12)];gh<-sigma[c(4:12)]
  #m[c(1:3)] <- c(0.17,0.37,0.45)
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mi)))
  ################iteration process##########################
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,9,m_sam[4]); swx <- matrix(0,9,1)
  g_a<-matrix(0,9,1);s0<-matrix(0,2,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:12)];ssigma<-sigma[c(4:12)]
    for(i in 1:d3) { WW[i,] <- mi[i]*dnorm(dataF3,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF3,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF3
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    m[c(1,2,3)]<-sumx[c(1,2,3)]/m_sam[c(1,2,3)]
    m[c(4:12)]<-sumwx[c(1:9)]/n0[1:9]
    ###########first order genetic parameter process####################
    hh17<-matrix(c(1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,
                   0,0,0,1,1,1,1,1,1,1,1,1,1,0,-1,1,1,1,0,0,0,-1,-1,-1,1,0,-1,1,0,-1,1,0,-1,1,0,-1,
                   0,1,0,0,0,0,0.5,0.5,0.5,0,0,0,0,1,0,0,0.5,0,0,0.5,0,0,0.5,0,1,0,1,1,0,-1,0,0,0,-1,0,1,
                   0,0,0,0,0.5,0,0,0,0,0,-0.5,0,0,0,0,0,0,0,0.5,0,-0.5,0,0,0,0,1,0,0,0,0,0,0.25,0,0,0,0),12,12)
    B17 <- solve(hh17,m)
    g_a[c(1,3,7,9)]<-0
    g_a[2]<-(0.5*(B17[6]+B17[9])^2+0.25*(B17[8]+B17[10])^2)/m_nf
    g_a[4]<-(0.5*(B17[5]+B17[9])^2+0.25*(B17[7]+B17[11])^2)/m_nf
    g_a[5]<-0.25*((B17[5]+B17[10])^2+(B17[6]+B17[11])^2+(B17[7]+0.5*B17[12])^2+(B17[8]+0.5*B17[12])^2+0.25*B17[12]^2)/m_nf
    g_a[6]<-(0.5*(B17[5]-B17[9])^2+0.25*(B17[7]-B17[11])^2)/m_nf
    g_a[8]<-(0.5*(B17[6]-B17[9])^2+0.25*(B17[8]-B17[10])^2)/m_nf
    #################obtain variance#######################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d3) {swx[i] <- WW[i,]%*%(dataF3-m[i+3])^2 }
    ############polygene variance process###########
    a1<-sigma[4];n_iter<-0;a2<-1000
    while(a2>0.0001){
      n_iter<-n_iter+1
      h1<-sigma[4]/(sigma[4]+g_a)
      s0[1]<-sum(swx*h1^2)
      s0[2]<-sum(n0*h1)
      sigma[4]<-s0[1]/s0[2]
      a2<-abs(sigma[4]-a1)
      a1<-sigma[4]
      if(n_iter>20) break
    }

    sigma40<-sigma[4]-sigma[1]/m_nf

    if (sigma40<0){sigma40<-0;sigma[4]<-sigma[1]/m_nf}
    sigma[c(5:12)]<-sigma[4]+g_a[c(2:9)]
    ##############error process#################
    a1<-sigma[1];n_iter<-0;a2<-sigma40;a2[a2<0]<-0;a4<-1000
    while(a4>0.0001){
      n_iter<-n_iter+1
      h1<-(sigma[1]/m_nf)/(a2+sigma[1]/m_nf+g_a)
      s0[1]<-sum(swx*h1^2*m_nf)
      s0[2]<-sum(h1*n0)
      s0[1]<-s0[1]+ss1+ss2+ss3
      s0[2]<-s0[2]+m_sam[1]+m_sam[2]+m_sam[3]
      sigma[1]<-s0[1]/s0[2]
      a4<-abs(sigma[1]-a1)
      a1<-sigma[1]
      if(n_iter>20) break
    }
    sigma[2]<-sigma[3]<-sigma[1]
    sigma[c(4:12)]<-sigma[1]/m_nf+sigma40+g_a[c(1:9)]
    ##############criteria for iterations to stop#############
    pi<-m[c(4:12)] ;gh<-sigma[c(4:12)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*14
  #################second order genetic parameter process###################
  jj<-sigmaf3 - sigma[4]
  gg<-sigma[4]-sigma[1]/m_nf
  jj[jj<0]<-0
  gg[gg<0]<-0
  ll<-jj/sigmaf3
  rr<-gg/sigmaf3
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma[1]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < m_sam[1]){P1bmw <- P1bmw+runif(m_sam[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:m_sam[1])) - 0.5)/m_sam[1])^2)
  P1u<- as.matrix(c(12*m_sam[1]*((P1dd[1]/m_sam[1]-0.5)^2),((45*m_sam[1])/4)*((P1dd[2]/m_sam[1]-1/3)^2),180*m_sam[1]*((P1dd[3]/m_sam[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,m_sam[1]))))

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*m_sam[2])
  F1bmw <- matrix(0,m_sam[2],1)
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma[2]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < m_sam[2]){F1bmw <- F1bmw+runif(m_sam[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:m_sam[2])) - 0.5)/m_sam[2])^2)
  F1u<- as.matrix(c(12*m_sam[2]*((F1dd[1]/m_sam[2]-0.5)^2),((45*m_sam[2])/4)*((F1dd[2]/m_sam[2]-1/3)^2),180*m_sam[2]*((F1dd[3]/m_sam[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,m_sam[2]))))

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*m_sam[3])
  P2bmw <- matrix(0,m_sam[3],1)
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma[3]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < m_sam[3]){P2bmw <- P2bmw+runif(m_sam[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:m_sam[3])) - 0.5)/m_sam[3])^2)
  P2u<- as.matrix(c(12*m_sam[3]*((P2dd[1]/m_sam[3]-0.5)^2),((45*m_sam[3])/4)*((P2dd[2]/m_sam[3]-1/3)^2),180*m_sam[3]*((P2dd[3]/m_sam[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,m_sam[3]))))

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3);
  F3w1<-1/(12*m_sam[4])
  F3bmw <- matrix(0,m_sam[4],1); F3bmwsl <- matrix(0,m_sam[4],d3)
  for(i in 1:d3){
    F3gg <- (dataF3 - pi[i])/sqrt(gh[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mix_pi[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[4]){F3P2 <- F3P2+runif(m_sam[4])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[4]) + sum((F3P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F3u <- as.matrix(c(12*m_sam[4]*((F3dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F3dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F3dd[3]/m_sam[4]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[4]))))

  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX2-ADI-ADI",round(abc,4),round(AIC,4),round(t(m),4),round(sigma[1],4),round(sigma[4],4),round(sigma[5],4),round(sigma[6],4),round(sigma[7],4),round(sigma[8],4),round(sigma[9],4),round(sigma[10],4),round(sigma[11],4),round(sigma[12],4),
                       round(t(mix_pi),4),round(B17[1],4),round(B17[2],4),round(B17[3],4),round(B17[4],4),round(B17[5],4),round(B17[6],4),round(B17[7],4),round(B17[8],4),round(B17[9],4),round(B17[10],4),round(B17[11],4),round(B17[12],4)," "," ",round(jj,4),round(ll*100,4),round(gg,4),round(rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
########################################MX2-ADI-AD(E-1)##############################################################
G4F3ModelFun[[19]] <- function(K1,logL,df11,df21,df31,df41,G4F3text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF3 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF3)))
  sigmaf3<-var(dataF3);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a9<-sqrt(sigmaf3)
  m_nf <- as.numeric(G4F3text2)
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF3)))
  ##################procedure start#############################################
  d3<-9
  mi<-as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2))
  m[c(1:4)]<-m
  if(m[1]<m[3]) {a9<--a9}
  m[5]<-m[4]+2.4*a9;m[6]<-m[4]+1.6*a9;m[7]<-m[4]+0.8*a9
  m[8]<-m[4];m[9]<-m[4]-0.8*a9;m[10]<-m[4]-1.6*a9
  m[11]<-m[4]-2.4*a9;m[12]<-m[4]-3.2*a9;m[4]<-m[4]+3.2*a9
  pi<-m[c(4:12)] ;gh<-sigma[c(4:12)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mi)))
  #####################iteration process##########################################
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,9,m_sam[4]); swx <- matrix(0,9,1)
  g_a<-matrix(0,9,1); s0<-matrix(0,2,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:12)];ssigma<-sigma[c(4:12)]
    for(i in 1:d3) { WW[i,] <- mi[i]*dnorm(dataF3,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF3,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF3
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    #################restrictions#################################
    aa1<-6*sumx[1]/m_sam[1]+4*sumx[2]/m_sam[2]+6*sumx[3]/m_sam[3]-7*(sumwx[1]/n0[1]+sumwx[9]/n0[9])-(sumwx[3]/n0[3]+sumwx[7]/n0[7])-16*sumwx[5]/n0[5]+4.0*(sumwx[2]/n0[2]+sumwx[4]/n0[4]+sumwx[6]/n0[6]+sumwx[8]/n0[8])
    aa2<-36*(sigma[1]/m_sam[1]+sigma[1]/m_sam[3])+16*sigma[1]/m_sam[2]+49*(sigma[4]/n0[1]+sigma[12]/n0[9])+(sigma[6]/n0[3]+sigma[10]/n0[7])+256*sigma[8]/n0[5]+16*(sigma[5]/n0[2]+sigma[7]/n0[4]+sigma[9]/n0[6]+sigma[11]/n0[8])
    aa3<-aa1/aa2
    m[1]<-(sumx[1]-sigma[1]*aa3*6)/m_sam[1];m[2]<-(sumx[2]-sigma[1]*aa3*4)/m_sam[2]
    m[3]<-(sumx[3]-sigma[1]*aa3*6)/m_sam[3];m[4]<-(sumwx[1]+sigma[4]*aa3*7)/n0[1]
    m[5]<-(sumwx[2]-sigma[5]*aa3*4)/n0[2];m[6]<-(sumwx[3]+sigma[6]*aa3)/n0[3]
    m[7]<-(sumwx[4]-sigma[7]*aa3*4)/n0[4];m[8]<-(sumwx[5]+sigma[8]*aa3*16)/n0[5]
    m[9]<-(sumwx[6]-sigma[9]*aa3*4)/n0[6];m[10]<-(sumwx[7]+sigma[10]*aa3)/n0[7]
    m[11]<-(sumwx[8]-sigma[11]*aa3*4)/n0[8];m[12]<-(sumwx[9]+sigma[12]*aa3*7)/n0[9]
    #############first order genetic parameter process#######################
    hh18<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,0,-1,1,1,1,0,0,0,-1,-1,-1,1,0,-1,1,0,-1,1,0,-1,1,0,-1,
                   0,1,0,0,0,0,0.5,0.5,0.5,0,0,0,0,1,0,0,0.5,0,0,0.5,0,0,0.5,0,1,0,1,1,0,-1,0,0,0,-1,0,1,
                   0,0,0,0,0.5,0,0,0,0,0,-0.5,0,0,0,0,0,0,0,0.5,0,-0.5,0,0,0,0,1,0,0,0,0,0,0.25,0,0,0,0,
                   1,0,-1,0,0,0,0,0,0,0,0,0,0,1,0,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25),12,11)
    B18<-solve(crossprod(hh18,hh18))%*%crossprod(hh18,m)
    g_a[c(1,3,7,9)]<-0
    g_a[2]<-(0.5*(B18[3]+B18[6])^2+0.25*(B18[5]+B18[7])^2)/m_nf;
    g_a[4]<-(0.5*(B18[2]+B18[6])^2+0.25*(B18[4]+B18[8])^2)/m_nf;
    g_a[5]<-0.25*((B18[2]+B18[7])^2+(B18[3]+B18[8])^2+(B18[4]+0.5*B18[9])^2+(B18[5]+0.5*B18[9])^2+0.25*B18[9]^2)/m_nf;
    g_a[6]<-(0.5*(B18[2]-B18[6])^2+0.25*(B18[4]-B18[8])^2)/m_nf;
    g_a[8]<-(0.5*(B18[3]-B18[6])^2+0.25*(B18[5]-B18[7])^2)/m_nf;
    ################obtain variance########################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d3) {swx[i] <- WW[i,]%*%(dataF3-m[i+3])^2 }
    ############polygene variance####################
    a1<-sigma[4];n_iter<-0;a2<-1000
    while(a2>0.0001){
      n_iter<-n_iter+1
      h1<-sigma[4]/(sigma[4]+g_a)
      s0[1]<-sum(swx*h1^2)
      s0[2]<-sum(n0*h1)
      sigma[4]<-s0[1]/s0[2]
      a2<-abs(sigma[4]-a1)
      a1<-sigma[4]
      if(n_iter>20) break
    }
    sigma[2]<-sigma[3]<-sigma[1]

    sigma40<-sigma[4]-sigma[1]/m_nf

    if (sigma40<0){sigma40<-0;sigma[4]<-sigma[1]/m_nf}
    sigma[c(5:12)]<-sigma[4]+g_a[c(2:9)]
    ################error process######################################
    a1<-sigma[1];n_iter<-0;a2<-sigma40;a2[a2<0]<-0;a4<-1000
    while(a4>0.0001){
      n_iter<-n_iter+1
      h1<-(sigma[1]/m_nf)/(a2+sigma[1]/m_nf+g_a)
      s0[1]<-sum(swx*h1^2*m_nf)
      s0[2]<-sum(h1*n0)
      s0[1]<-s0[1]+ss1+ss2+ss3
      s0[2]<-s0[2]+m_sam[1]+m_sam[2]+m_sam[3]
      sigma[1]<-s0[1]/s0[2]
      a4<-abs(sigma[1]-a1)
      a1<-sigma[1]
      if(n_iter>20) break
    }
    sigma[c(4:12)]<-sigma[1]/m_nf+sigma40+g_a[c(1:9)]
    ###############criteria for iterations to stop###################
    pi<-m[c(4:12)] ;gh<-sigma[c(4:12)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*13
  ################second order genetic parameter process################
  jj<-sigmaf3 - sigma[4];
  gg<-sigma[4]-sigma[1]/m_nf;
  jj[jj<0]<-0
  gg[gg<0]<-0
  ll<-jj/sigmaf3
  rr<-gg/sigmaf3
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma[1]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < m_sam[1]){P1bmw <- P1bmw+runif(m_sam[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:m_sam[1])) - 0.5)/m_sam[1])^2)
  P1u<- as.matrix(c(12*m_sam[1]*((P1dd[1]/m_sam[1]-0.5)^2),((45*m_sam[1])/4)*((P1dd[2]/m_sam[1]-1/3)^2),180*m_sam[1]*((P1dd[3]/m_sam[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,m_sam[1]))))

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*m_sam[2])
  F1bmw <- matrix(0,m_sam[2],1)
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma[2]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < m_sam[2]){F1bmw <- F1bmw+runif(m_sam[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:m_sam[2])) - 0.5)/m_sam[2])^2)
  F1u<- as.matrix(c(12*m_sam[2]*((F1dd[1]/m_sam[2]-0.5)^2),((45*m_sam[2])/4)*((F1dd[2]/m_sam[2]-1/3)^2),180*m_sam[2]*((F1dd[3]/m_sam[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,m_sam[2]))))

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*m_sam[3])
  P2bmw <- matrix(0,m_sam[3],1)
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma[3]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < m_sam[3]){P2bmw <- P2bmw+runif(m_sam[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:m_sam[3])) - 0.5)/m_sam[3])^2)
  P2u<- as.matrix(c(12*m_sam[3]*((P2dd[1]/m_sam[3]-0.5)^2),((45*m_sam[3])/4)*((P2dd[2]/m_sam[3]-1/3)^2),180*m_sam[3]*((P2dd[3]/m_sam[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,m_sam[3]))))

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3);
  F3w1<-1/(12*m_sam[4])
  F3bmw <- matrix(0,m_sam[4],1); F3bmwsl <- matrix(0,m_sam[4],d3)
  for(i in 1:d3){
    F3gg <- (dataF3 - pi[i])/sqrt(gh[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mix_pi[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[4]){F3P2 <- F3P2+runif(m_sam[4])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[4]) + sum((F3P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F3u <- as.matrix(c(12*m_sam[4]*((F3dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F3dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F3dd[3]/m_sam[4]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[4]))))

  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX2-ADI-AD",round(abc,4),round(AIC,4),round(t(m),4),round(sigma[1],4),round(sigma[4],4),round(sigma[5],4),round(sigma[6],4),round(sigma[7],4),round(sigma[8],4),round(sigma[9],4),round(sigma[10],4),round(sigma[11],4),round(sigma[12],4),
                       round(t(mix_pi),4),round(B18[1],4)," "," "," ",round(B18[2],4),round(B18[3],4),round(B18[4],4),round(B18[5],4),round(B18[6],4),round(B18[7],4),round(B18[8],4),round(B18[9],4),round(B18[10],4),round(B18[11],4),round(jj,4),round(ll*100,4),round(gg,4),round(rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
############################MX2-AD-AD(E-2)##########################################################
G4F3ModelFun[[20]] <- function(K1,logL,df11,df21,df31,df41,G4F3text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF3 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF3)))
  sigmaf3<-var(dataF3);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a9<-sqrt(sigmaf3/(m_sam[4]*(m_sam[4]-1)))
  m_nf <- as.numeric(G4F3text2)
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF3)))
  ##################procedure start#############################################
  d3<-9
  mi<-as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2))
  m[c(1:4)]<-m
  if (m[1]<m[3]) {a9<--a9}
  m[5]<-m[4]+2.4*a9;m[6]<-m[4]+1.6*a9;m[7]<-m[4]+0.8*a9
  m[8]<-m[4];m[9]<-m[4]-0.8*a9;m[10]<-m[4]-1.6*a9
  m[11]<-m[4]-2.4*a9;m[12]<-m[4]-3.2*a9;m[4]<-m[4]+3.2*a9
  pi<-m[c(4:12)];gh<-sigma[c(4:12)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mi)))
  ####################iteration process#############################
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,9,m_sam[4]); swx <- matrix(0,9,1)
  hh<-matrix(0,5,5); b_line<-matrix(0,5,1)
  g_a<-matrix(0,9,1);s0<-matrix(0,2,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:12)]; ssigma<-sigma[c(4:12)]
    for(i in 1:d3) { WW[i,] <- mi[i]*dnorm(dataF3,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF3,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF3
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    m[c(1,2,3)]<-sumx[c(1,2,3)]/m_sam[c(1,2,3)]
    m[c(4:12)]<-sumwx[c(1:9)]/n0[1:9]
    ##########################################################################
    hh[1,1]<-25.0*sigma[1]/m_sam[1]+4.0*sigma[1]/m_sam[2]+25.0*sigma[1]/m_sam[3]+4.0*sigma[4]/n0[1]+4.0*sigma[6]/n0[3]+16.0*sigma[8]/n0[5]+4.0*sigma[10]/n0[7]+4.0*sigma[12]/n0[9]
    hh[1,2]<--2.0*sigma[4]/n0[1]-4.0*sigma[8]/n0[5]
    hh[1,3]<--2.0*sigma[4]/n0[1]+2.0*sigma[10]/n0[7]
    hh[1,4]<--2.0*sigma[4]/n0[1]+2.0*sigma[6]/n0[3]+2.0*sigma[10]/ n0[7]-2.0*sigma[12]/n0[9]
    hh[1,5]<--4.0*sigma[8]/n0[5]-2.0*sigma[12]/n0[9]
    hh[2,2]<-sigma[4]/n0[1]+sigma[5]/n0[2]+sigma[7]/n0[4]+sigma[8]/n0[5]
    hh[2,3]<-sigma[4]/n0[1]+sigma[5]/n0[2]
    hh[2,4]<-sigma[4]/n0[1]
    hh[2,5]<-sigma[8]/n0[5]
    hh[3,3]<-sigma[4]/n0[1]+sigma[5]/n0[2]+sigma[10]/n0[7]+sigma[11]/n0[8]
    hh[3,4]<-sigma[4]/n0[1]+sigma[10]/n0[7]
    hh[3,5]<--sigma[11]/n0[8]
    hh[4,4]<-sigma[4]/n0[1]+sigma[6]/n0[3]+sigma[10]/n0[7]+sigma[12]/n0[9]
    hh[4,5]<-sigma[12]/n0[9]
    hh[5,5]<-sigma[8]/n0[5]+sigma[9]/n0[6]+sigma[11]/n0[8]+sigma[12]/n0[9]
    for(i in 2:5)
    {
      for(j in 1:(i-1))
      {
        hh[i,j]<-hh[j,i]
      }
    }
    ###################################################################################
    b_line[1]<-5*sumx[1]/m_sam[1]+2*sumx[2]/m_sam[2]+5*sumx[3]/m_sam[3]-2*(sumwx[1]/n0[1]+sumwx[3]/n0[3]+sumwx[7]/n0[7]+sumwx[9]/n0[9])-4*sumwx[5]/n0[5]
    b_line[2]<-sumwx[1]/n0[1]-sumwx[2]/n0[2]-sumwx[4]/n0[4]+sumwx[5]/n0[5]
    b_line[3]<-sumwx[1]/n0[1]-sumwx[2]/n0[2]-sumwx[7]/n0[7]+sumwx[8]/n0[8]
    b_line[4]<-sumwx[1]/n0[1]-sumwx[3]/n0[3]-sumwx[7]/n0[7]+sumwx[9]/n0[9]
    b_line[5]<-sumwx[5]/n0[5]-sumwx[6]/n0[6]-sumwx[8]/n0[8]+sumwx[9]/n0[9]
    B19<-solve(hh,b_line)
    ####################################################################
    m[1]<-(sumx[1]-sigma[1]*B19[1]*5)/m_sam[1]
    m[2]<-(sumx[2]-sigma[1]*B19[1]*2)/m_sam[2]
    m[3]<-(sumx[3]-sigma[1]*B19[1]*5)/m_sam[3]
    m[4]<-(sumwx[1]+sigma[4]*(2*B19[1]-B19[2]-B19[3]-B19[4]))/n0[1]
    m[5]<-(sumwx[2]+sigma[5]*(B19[2]+B19[3]))/n0[2]
    m[6]<-(sumwx[3]+sigma[6]*(2*B19[1]+B19[4]))/n0[3]
    m[7]<-(sumwx[4]+sigma[7]*B19[2])/n0[4]
    m[8]<-(sumwx[5]+sigma[8]*(4*B19[1]-B19[2]-B19[5]))/n0[5]
    m[9]<-(sumwx[6]+sigma[9]*B19[5])/n0[6]
    m[10]<-(sumwx[7]+sigma[10]*(2*B19[1]+B19[3]+B19[4]))/n0[7]
    m[11]<-(sumwx[8]+sigma[11]*(-B19[3]+B19[5]))/n0[8]
    m[12]<-(sumwx[9]+sigma[12]*(2*B19[1]-B19[4]-B19[5]))/n0[9]
    #########first order genetic parameter process#################
    hh191<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,0,-1,1,1,1,0,0,0,-1,-1,-1,1,0,-1,1,0,-1,1,0,-1,1,0,-1,
                    0,1,0,0,0,0,0.5,0.5,0.5,0,0,0,0,1,0,0,0.5,0,0,0.5,0,0,0.5,0,
                    1,0,-1,0,0,0,0,0,0,0,0,0,0,1,0,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25),12,7)
    B191<-solve(crossprod(hh191,hh191))%*%crossprod(hh191,m)
    g_a[c(1,3,7,9)]<-0
    g_a[2]<-(0.5*B191[3]^2+0.25*B191[5]^2)/m_nf;g_a[4]<-(0.5*B191[2]^2+0.25*B191[4]^2)/m_nf
    g_a[5]<-(0.5*B191[2]^2+0.5*B191[3]^2+0.25*B191[4]^2+0.25*B191[5]^2)/m_nf
    g_a[6]<-(0.5*B191[2]^2+0.25*B191[4]^2)/m_nf;g_a[8]<-(0.5*B191[3]^2+0.25*B191[5]^2)/m_nf
    ################obtain variance##################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d3) {swx[i] <- WW[i,]%*%(dataF3-m[i+3])^2 }
    ###########polygene variance process###############
    a1<-sigma[4];n_iter<-0;a2<-1000
    while(a2>0.0001){
      n_iter<-n_iter+1
      h1<-sigma[4]/(sigma[4]+g_a)
      s0[1]<-sum(swx*h1^2)
      s0[2]<-sum(n0*h1)
      sigma[4]<-s0[1]/s0[2]
      a2<-abs(sigma[4]-a1)
      a1<-sigma[4]
      if(n_iter>20) break
    }

    sigma40<-sigma[4]-sigma[1]/m_nf

    if (sigma40<0){sigma40<-0;sigma[4]<-sigma[1]/m_nf}
    sigma[c(5:12)]<-sigma[4]+g_a[c(2:9)]
    ################error process#####################
    a1<-sigma[1];n_iter<-0;a2<-sigma40;a2[a2<0]<-0;a4<-1000
    while(a4>0.0001){
      n_iter<-n_iter+1
      h1<-(sigma[1]/m_nf)/(a2+sigma[1]/m_nf+g_a)
      s0[1]<-sum(swx*h1^2*m_nf)
      s0[2]<-sum(h1*n0)
      s0[1]<-s0[1]+ss1+ss2+ss3
      s0[2]<-s0[2]+m_sam[1]+m_sam[2]+m_sam[3]
      sigma[1]<-s0[1]/s0[2]
      a4<-abs(sigma[1]-a1)
      a1<-sigma[1]
      if(n_iter>20) break
    }
    sigma[2]<-sigma[3]<-sigma[1]
    sigma[c(4:12)]<-sigma[1]/m_nf+sigma40+g_a[c(1:9)]
    ############criteria for iterations to stop###########
    pi<-m[c(4:12)];gh<-sigma[c(4:12)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*9
  ############second order genetic parameter process#############################################
  jj<-sigmaf3 - sigma[4];
  gg<-sigma[4]-sigma[1]/m_nf;
  jj[jj<0]<-0
  gg[gg<0]<-0
  ll<-jj/sigmaf3
  rr<-gg/sigmaf3
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma[1]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < m_sam[1]){P1bmw <- P1bmw+runif(m_sam[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:m_sam[1])) - 0.5)/m_sam[1])^2)
  P1u<- as.matrix(c(12*m_sam[1]*((P1dd[1]/m_sam[1]-0.5)^2),((45*m_sam[1])/4)*((P1dd[2]/m_sam[1]-1/3)^2),180*m_sam[1]*((P1dd[3]/m_sam[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,m_sam[1]))))

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*m_sam[2])
  F1bmw <- matrix(0,m_sam[2],1)
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma[2]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < m_sam[2]){F1bmw <- F1bmw+runif(m_sam[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:m_sam[2])) - 0.5)/m_sam[2])^2)
  F1u<- as.matrix(c(12*m_sam[2]*((F1dd[1]/m_sam[2]-0.5)^2),((45*m_sam[2])/4)*((F1dd[2]/m_sam[2]-1/3)^2),180*m_sam[2]*((F1dd[3]/m_sam[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,m_sam[2]))))

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*m_sam[3])
  P2bmw <- matrix(0,m_sam[3],1)
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma[3]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < m_sam[3]){P2bmw <- P2bmw+runif(m_sam[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:m_sam[3])) - 0.5)/m_sam[3])^2)
  P2u<- as.matrix(c(12*m_sam[3]*((P2dd[1]/m_sam[3]-0.5)^2),((45*m_sam[3])/4)*((P2dd[2]/m_sam[3]-1/3)^2),180*m_sam[3]*((P2dd[3]/m_sam[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,m_sam[3]))))

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3);
  F3w1<-1/(12*m_sam[4])
  F3bmw <- matrix(0,m_sam[4],1); F3bmwsl <- matrix(0,m_sam[4],d3)
  for(i in 1:d3){
    F3gg <- (dataF3 - pi[i])/sqrt(gh[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mix_pi[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[4]){F3P2 <- F3P2+runif(m_sam[4])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[4]) + sum((F3P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F3u <- as.matrix(c(12*m_sam[4]*((F3dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F3dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F3dd[3]/m_sam[4]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[4]))))

  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX2-AD-AD",round(abc,4),round(AIC,4),round(t(m),4),round(sigma[1],4),round(sigma[4],4),round(sigma[5],4),round(sigma[6],4),round(sigma[7],4),round(sigma[8],4),round(sigma[9],4),round(sigma[10],4),round(sigma[11],4),round(sigma[12],4),
                       round(t(mix_pi),4),round(B191[1],4)," "," "," ",round(B191[2],4),round(B191[3],4),round(B191[4],4),round(B191[5],4)," "," "," "," ",round(B191[6],4),round(B191[7],4),round(jj,4),round(ll*100,4),round(gg,4),round(rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
#####################################MX2-A-AD(E-3)####################################################################
G4F3ModelFun[[21]] <- function(K1,logL,df11,df21,df31,df41,G4F3text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF3 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF3)))
  sigmaf3<-var(dataF3)
  sigmaP1<-var(dataP1)
  sigmaP2<-var(dataP2)
  sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a9<-sqrt(sigmaf3)
  m_nf <- as.numeric(G4F3text2)
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF3)))
  ##################procedure start#############################################
  d3<-9
  mi<-as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2))
  m[c(1:4)]<-m
  if (m[1]<m[3]) {a9<--a9};
  m[5]<-m[4]+2.4*a9;m[6]<-m[4]+1.6*a9;m[7]<-m[4]+0.8*a9;
  m[8]<-m[4];m[9]<-m[4]-0.8*a9;m[10]<-m[4]-1.6*a9;
  m[11]<-m[4]-2.4*a9;m[12]<-m[4]-3.2*a9;m[4]<-m[4]+3.2*a9;
  pi<-m[c(4:12)];gh<-sigma[c(4:12)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mi)))
  #####################iteration process##############################
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,9,m_sam[4]); swx <- matrix(0,9,1)
  hh<-matrix(0,7,7);b_line<-matrix(0,7,1)
  g_a<-matrix(0,9,1);s0<-matrix(0,2,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:12)];ssigma<-sigma[c(4:12)]
    for(i in 1:d3) { WW[i,] <- mi[i]*dnorm(dataF3,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF3,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF3
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    ################################################################
    hh[1,1]<-9*sigma[1]/m_sam[1]+4*sigma[1]/m_sam[2]+9*sigma[1]/m_sam[3]+64*sigma[8]/n0[5]
    hh[1,2]<-0
    hh[1,3]<-16*sigma[8]/n0[5]
    hh[1,4]<-0
    hh[1,5]<-16*sigma[8]/n0[5]
    hh[1,6]<-16*sigma[8]/n0[5]
    hh[1,7]<-0
    hh[2,2]<-sigma[4]/n0[1]+4*sigma[5]/n0[2]+sigma[6]/n0[3]
    hh[2,3]<-sigma[4]/n0[1]
    hh[2,4]<--2*sigma[5]/n0[2]-sigma[6]/n0[3]
    hh[2,5]<-sigma[6]/n0[3]
    hh[2,6]<-0
    hh[2,7]<-0
    hh[3,3]<-sigma[4]/n0[1]+4*sigma[8]/n0[5]+sigma[12]/n0[9]
    hh[3,4]<-sigma[12]/n0[9]
    hh[3,5]<-4*sigma[8]/n0[5]
    hh[3,6]<-4*sigma[8]/n0[5]
    hh[3,7]<-2*sigma[12]/n0[9]
    hh[4,4]<-sigma[5]/n0[2]+sigma[6]/n0[3]+sigma[11]/n0[8]+sigma[12]/n0[9]
    hh[4,5]<--sigma[6]/n0[3]
    hh[4,6]<-0
    hh[4,7]<-2*(sigma[11]/n0[8]+sigma[12]/n0[9])
    hh[5,5]<-sigma[6]/n0[3]+4.0*sigma[8]/n0[5]+sigma[10]/n0[7]
    hh[5,6]<-4*sigma[8]/n0[5]
    hh[5,7]<-0
    hh[6,6]<-sigma[7]/n0[4]+4*sigma[8]/n0[5]+sigma[9]/n0[6]
    hh[6,7]<-sigma[7]/n0[4]-sigma[9]/n0[6]
    hh[7,7]<-sigma[7]/n0[4]+sigma[9]/n0[6]+4*sigma[11]/n0[8]+4*sigma[12]/n0[9]
    for(i in 2:7)
    {
      for(j in 1:(i-1))
      {
        hh[i,j]<-hh[j,i]
      }
    }
    ##############################################################################
    b_line[1]<-3*sumx[1]/m_sam[1]+2*sumx[2]/m_sam[2]+3*sumx[3]/m_sam[3]-8*sumwx[5]/n0[5]
    b_line[2]<-sumwx[1]/n0[1]-2*sumwx[2]/n0[2]+sumwx[3]/n0[3]
    b_line[3]<-sumwx[1]/n0[1]-2*sumwx[5]/n0[5]+sumwx[9]/n0[9]
    b_line[4]<-sumwx[2]/n0[2]-sumwx[3]/n0[3]-sumwx[8]/n0[8]+sumwx[9]/n0[9]
    b_line[5]<-sumwx[3]/n0[3]-2*sumwx[5]/n0[5]+sumwx[7]/n0[7]
    b_line[6]<-sumwx[4]/n0[4]-2*sumwx[5]/n0[5]+sumwx[6]/n0[6]
    b_line[7]<-sumwx[4]/n0[4]-sumwx[6]/n0[6]-2*sumwx[8]/n0[8]+2*sumwx[9]/n0[9]
    B20<-solve(hh,b_line)
    ###############################################################################
    m[1]<-(sumx[1]-sigma[1]*B20[1]*3)/m_sam[1]
    m[2]<-(sumx[2]-sigma[1]*B20[1]*2)/m_sam[2]
    m[3]<-(sumx[3]-sigma[1]*B20[1]*3)/m_sam[3]
    m[4]<-(sumwx[1]-sigma[4]*(B20[2]+B20[3]))/n0[1]
    m[5]<-(sumwx[2]+sigma[5]*(2*B20[2]-B20[4]))/n0[2]
    m[6]<-(sumwx[3]+sigma[6]*(-B20[2]+B20[4]-B20[5]))/n0[3]
    m[7]<-(sumwx[4]+sigma[7]*(-B20[6]-B20[7]))/n0[4]
    m[8]<-(sumwx[5]+sigma[8]*(8*B20[1]+2*B20[3]+2*B20[5]+2.0*B20[6]))/n0[5]
    m[9]<-(sumwx[6]+sigma[9]*(-B20[6]+B20[7]))/n0[6]
    m[10]<-(sumwx[7]-sigma[10]*B20[5])/n0[7]
    m[11]<-(sumwx[8]+sigma[11]*(B20[4]+2*B20[7]))/n0[8]
    m[12]<-(sumwx[9]-sigma[12]*(B20[3]+B20[4]+2*B20[7]))/n0[9]
    #################first order genetic parameter process###############
    hh201<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,0,-1,1,1,1,0,0,0,-1,-1,-1,1,0,-1,1,0,-1,1,0,-1,1,0,-1,1,
                    0,-1,0,0,0,0,0,0,0,0,0,0,1,0,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25),12,5)
    B201<-solve(crossprod(hh201,hh201))%*%crossprod(hh201,m)
    g_a[c(1,3,7,9)]<-0
    g_a[2]<-0.5*B201[3]^2/m_nf
    g_a[4]<-0.5*B201[2]^2/m_nf
    g_a[5]<-(0.5*B201[2]^2+0.5*B201[3]^2)/m_nf
    g_a[6]<-0.5*B201[2]^2/m_nf
    g_a[8]<-0.5*B201[3]^2/m_nf
    ###########obtain variance##############################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d3) {swx[i] <- WW[i,]%*%(dataF3-m[i+3])^2 }
    ##########polygene variance process######################
    a1<-sigma[4];n_iter<-0;a2<-1000
    while(a2>0.0001){
      n_iter<-n_iter+1
      h1<-sigma[4]/(sigma[4]+g_a)
      s0[1]<-sum(swx*h1^2)
      s0[2]<-sum(n0*h1)
      sigma[4]<-s0[1]/s0[2]
      a2<-abs(sigma[4]-a1)
      a1<-sigma[4]
      if(n_iter>20) break
    }

    sigma40<-sigma[4]-sigma[1]/m_nf

    if (sigma40<0){sigma40<-0;sigma[4]<-sigma[1]/m_nf}
    sigma[c(5:12)]<-sigma[4]+g_a[c(2:9)]
    ############error process###########################
    a1<-sigma[1];n_iter<-0;a2<-sigma40;a2[a2<0]<-0;a4<-1000
    while(a4>0.0001){
      n_iter<-n_iter+1
      h1<-(sigma[1]/m_nf)/(a2+sigma[1]/m_nf+g_a)
      s0[1]<-sum(swx*h1^2*m_nf)
      s0[2]<-sum(h1*n0)
      s0[1]<-s0[1]+ss1+ss2+ss3
      s0[2]<-s0[2]+m_sam[1]+m_sam[2]+m_sam[3]
      sigma[1]<-s0[1]/s0[2]
      a4<-abs(sigma[1]-a1)
      a1<-sigma[1]
      if(n_iter>20) break
    }
    sigma[2]<-sigma[3]<-sigma[1]
    sigma[c(4:12)]<-sigma[1]/m_nf+sigma40+g_a[c(1:9)]
    ##########criteria for iterations to stop####################
    pi<-m[c(4:12)];gh<-sigma[c(4:12)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*9
  ############second order genetic parameter process#######################
  jj<-sigmaf3 - sigma[4];
  gg<-sigma[4]-sigma[1]/m_nf;
  jj[jj<0]<-0
  gg[gg<0]<-0
  ll<-jj/sigmaf3
  rr<-gg/sigmaf3
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma[1]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < m_sam[1]){P1bmw <- P1bmw+runif(m_sam[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:m_sam[1])) - 0.5)/m_sam[1])^2)
  P1u<- as.matrix(c(12*m_sam[1]*((P1dd[1]/m_sam[1]-0.5)^2),((45*m_sam[1])/4)*((P1dd[2]/m_sam[1]-1/3)^2),180*m_sam[1]*((P1dd[3]/m_sam[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,m_sam[1]))))

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*m_sam[2])
  F1bmw <- matrix(0,m_sam[2],1)
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma[2]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < m_sam[2]){F1bmw <- F1bmw+runif(m_sam[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:m_sam[2])) - 0.5)/m_sam[2])^2)
  F1u<- as.matrix(c(12*m_sam[2]*((F1dd[1]/m_sam[2]-0.5)^2),((45*m_sam[2])/4)*((F1dd[2]/m_sam[2]-1/3)^2),180*m_sam[2]*((F1dd[3]/m_sam[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,m_sam[2]))))

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*m_sam[3])
  P2bmw <- matrix(0,m_sam[3],1)
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma[3]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < m_sam[3]){P2bmw <- P2bmw+runif(m_sam[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:m_sam[3])) - 0.5)/m_sam[3])^2)
  P2u<- as.matrix(c(12*m_sam[3]*((P2dd[1]/m_sam[3]-0.5)^2),((45*m_sam[3])/4)*((P2dd[2]/m_sam[3]-1/3)^2),180*m_sam[3]*((P2dd[3]/m_sam[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,m_sam[3]))))

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3);
  F3w1<-1/(12*m_sam[4])
  F3bmw <- matrix(0,m_sam[4],1); F3bmwsl <- matrix(0,m_sam[4],d3)
  for(i in 1:d3){
    F3gg <- (dataF3 - pi[i])/sqrt(gh[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mix_pi[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[4]){F3P2 <- F3P2+runif(m_sam[4])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[4]) + sum((F3P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F3u <- as.matrix(c(12*m_sam[4]*((F3dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F3dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F3dd[3]/m_sam[4]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[4]))))

  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX2-A-AD",round(abc,4),round(AIC,4),round(t(m),4),round(sigma[1],4),round(sigma[4],4),round(sigma[5],4),round(sigma[6],4),round(sigma[7],4),round(sigma[8],4),round(sigma[9],4),round(sigma[10],4),round(sigma[11],4),round(sigma[12],4),
                       round(t(mix_pi),4),round(B201[1],4)," "," "," ",round(B201[2],4),round(B201[3],4)," "," ", " "," "," "," ",round(B201[4],4),round(B201[5],4),round(jj,4),round(ll*100,4),round(gg,4),round(rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
#################################MX2-EA-AD(E-4)####################################################################
G4F3ModelFun[[22]] <- function(K1,logL,df11,df21,df31,df41,G4F3text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF3 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF3)))
  sigmaf3<-var(dataF3);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a9<-sqrt(sigmaf3)
  m_nf <- as.numeric(G4F3text2)
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF3)))
  ##################procedure start#############################################
  d2<-6
  mi<-as.matrix(c(0.0625,0.25,0.125,0.25,0.25,0.0625))
  sigma<-as.matrix(c(10,10,10,10,172,10,172,172,10))
  m[c(1:4)]<-m
  if (m[1]<m[3]) {a9<--a9}
  m[5]<-m[4]+1.5*a9;m[6]<-m[4];m[7]<-sumx[4]/m_sam[4]
  m[8]<-m[4]-1.5*a9;m[9]<-m[4]-3*a9;m[4]<-m[4]+3*a9
  pi<-m[c(4:9)];gh<-sigma[c(4:9)]
  m[c(4:9)] <- c(138,120,102,84,54,66.5)
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mi)))
  ####################iteration process###############################
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,6,m_sam[4]); swx <- matrix(0,6,1)
  hh<-matrix(0,4,4);b_line<-matrix(0,4,1)
  g_a<-matrix(0,6,1);s0<-matrix(0,2,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:9)];ssigma<-sigma[c(4:9)]
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataF3,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF3,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF3
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001

    ##############restrictions#################################################
    a1<-sumwx[3]*sigma[7]+sumwx[4]*sigma[6]
    a2<-n0[3]*sigma[7]+n0[4]*sigma[6]
    a3<-sigma[6]*sigma[7]
    hh[1,1]<-9*sigma[1]/m_sam[1]+4*sigma[1]/m_sam[2]+9*sigma[1]/m_sam[3]+64*a3/a2
    hh[1,2]<--8*a3/a2
    hh[1,3]<-16*a3/a2
    hh[1,4]<--8*a3/a2
    hh[2,2]<-sigma[4]/n0[1]+4*sigma[5]/n0[2]+a3/a2
    hh[2,3]<--2*sigma[5]/n0[2]-2*a3/a2
    hh[2,4]<-a3/a2
    hh[3,3]<-sigma[5]/n0[2]+4*a3/a2+sigma[8]/n0[5]
    hh[3,4]<--2*a3/a2-2*sigma[8]/n0[5]
    hh[4,4]<-a3/a2+4*sigma[8]/n0[5]+sigma[9]/n0[6]
    #######################################################################
    b_line[1]<-3*sumx[1]/m_sam[1]+2*sumx[2]/m_sam[2]+3*sumx[3]/m_sam[3]-8*a1/a2
    b_line[2]<-sumwx[1]/n0[1]-2*sumwx[2]/n0[2]+a1/a2
    b_line[3]<-sumwx[2]/n0[2]-2*a1/a2+sumwx[5]/n0[5]
    b_line[4]<-a1/a2-2*sumwx[5]/n0[5]+sumwx[6]/n0[6]
    B21<-solve(hh,b_line)

    # b_line[1]<-3*sumx[1]/m_sam[1]+2*sumx[2]/m_sam[2]+3*sumx[3]/m_sam[3]-8*sumwx[3]/n0[3]
    # b_line[2]<-sumwx[1]/n0[1]-2*sumwx[2]/n0[2]+sumwx[3]/n0[3]
    # b_line[3]<-sumwx[2]/n0[2]-2*a1/a2+sumwx[5]/n0[5]
    # b_line[4]<-a1/a2-2*sumwx[5]/n0[5]+sumwx[6]/n0[6]
    # B21<-solve(hh,b_line)
    # ##########################################################################
    m[1]<-(sumx[1]-sigma[1]*B21[1]*3)/m_sam[1]
    m[2]<-(sumx[2]-sigma[1]*B21[1]*2)/m_sam[2]
    m[3]<-(sumx[3]-sigma[1]*B21[1]*3)/m_sam[3]
    m[4]<-(sumwx[1]-sigma[4]*B21[2])/n0[1]
    m[5]<-(sumwx[2]+sigma[5]*(2*B21[2]-B21[3]))/n0[2]
    m[6]<-(sumwx[3]*sigma[7]+sumwx[4]*sigma[6]+sigma[6]*sigma[7]*(8*B21[1]-B21[2]+2.0*B21[3]-B21[4]))/(n0[3]*sigma[7]+n0[4]*sigma[6])
    m[7]<-m[3]
    m[8]<-(sumwx[5]+sigma[8]*(-B21[3]+2*B21[4]))/n0[5]
    m[9]<-(sumwx[6]-sigma[9]*B21[4])/n0[6]
    ###############first order genetic parameter process#######################
    hh211<-matrix(c(1,1,1,1,1,1,1,1,2,0,-2,2,1,0,-1,-2,1,0,-1,0,0,0,0,0,0,1,0,0.25,0.25,0.25,0.25,0.25),8,4)
    B211<-solve(crossprod(hh211,hh211))%*%crossprod(hh211,m[c(1,2,3,4,5,6,8,9)])
    g_a[c(1,3,6)]<-0
    g_a[2]<-0.5*B211[2]^2/m_nf
    g_a[4]<-B211[2]^2/m_nf
    g_a[5]<-0.5*B211[2]^2/m_nf
    ##########obtain variance#########################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d2) {swx[i] <- WW[i,]%*%(dataF3-m[i+3])^2 }
    #########polygene variance process###################
    a1<-sigma[4];n_iter<-0;a2<-1000
    while(a2>0.0001){
      n_iter<-n_iter+1
      h1<-sigma[4]/(sigma[4]+g_a)
      s0[1]<-sum(swx*h1^2)
      s0[2]<-sum(n0*h1)
      sigma[4]<-s0[1]/s0[2]
      a2<-abs(sigma[4]-a1)
      a1<-sigma[4]
      if(n_iter>20) break
    }

    sigma40<-sigma[4]-sigma[1]/m_nf

    if (sigma40<0){sigma40<-0;sigma[4]<-sigma[1]/m_nf}
    sigma[c(5:9)]<-sigma[4]+g_a[c(2:6)]
    ###########error process#######################
    a1<-sigma[1]; n_iter<-0;a2<-sigma40;a2[a2<0]<-0;a4<-1000
    while(a4>0.0001){
      n_iter<-n_iter+1
      h1<-(sigma[1]/m_nf)/(a2+sigma[1]/m_nf+g_a)
      s0[1]<-sum(swx*h1^2*m_nf)
      s0[2]<-sum(h1*n0)
      s0[1]<-s0[1]+ss1+ss2+ss3
      s0[2]<-s0[2]+m_sam[1]+m_sam[2]+m_sam[3]
      sigma[1]<-s0[1]/s0[2]
      a4<-abs(sigma[1]-a1)
      a1<-sigma[1]
      if(n_iter>20) break
    }
    sigma[2]<-sigma[3]<-sigma[1]
    sigma[c(4:9)]<-sigma[1]/m_nf+sigma40+g_a[c(1:6)]
    #############criteria for iterations to stop########
    pi<-m[c(4:9)];gh<-sigma[c(4:9)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*6
  ###########second order genetic parameter process#################
  jj<-sigmaf3 - sigma[4];
  gg<-sigma[4]-sigma[1]/m_nf;
  jj[jj<0]<-0
  gg[gg<0]<-0
  ll<-jj/sigmaf3
  rr<-gg/sigmaf3
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma[1]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < m_sam[1]){P1bmw <- P1bmw+runif(m_sam[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:m_sam[1])) - 0.5)/m_sam[1])^2)
  P1u<- as.matrix(c(12*m_sam[1]*((P1dd[1]/m_sam[1]-0.5)^2),((45*m_sam[1])/4)*((P1dd[2]/m_sam[1]-1/3)^2),180*m_sam[1]*((P1dd[3]/m_sam[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,m_sam[1]))))

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*m_sam[2])
  F1bmw <- matrix(0,m_sam[2],1)
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma[2]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < m_sam[2]){F1bmw <- F1bmw+runif(m_sam[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:m_sam[2])) - 0.5)/m_sam[2])^2)
  F1u<- as.matrix(c(12*m_sam[2]*((F1dd[1]/m_sam[2]-0.5)^2),((45*m_sam[2])/4)*((F1dd[2]/m_sam[2]-1/3)^2),180*m_sam[2]*((F1dd[3]/m_sam[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,m_sam[2]))))

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*m_sam[3])
  P2bmw <- matrix(0,m_sam[3],1)
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma[3]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < m_sam[3]){P2bmw <- P2bmw+runif(m_sam[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:m_sam[3])) - 0.5)/m_sam[3])^2)
  P2u<- as.matrix(c(12*m_sam[3]*((P2dd[1]/m_sam[3]-0.5)^2),((45*m_sam[3])/4)*((P2dd[2]/m_sam[3]-1/3)^2),180*m_sam[3]*((P2dd[3]/m_sam[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,m_sam[3]))))

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3);
  F3w1<-1/(12*m_sam[4])
  F3bmw <- matrix(0,m_sam[4],1); F3bmwsl <- matrix(0,m_sam[4],d2)
  for(i in 1:d2){
    F3gg <- (dataF3 - pi[i])/sqrt(gh[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mix_pi[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[4]){F3P2 <- F3P2+runif(m_sam[4])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[4]) + sum((F3P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F3u <- as.matrix(c(12*m_sam[4]*((F3dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F3dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F3dd[3]/m_sam[4]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[4]))))

  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX2-EA-AD",round(abc,4),round(AIC,4),round(t(m),4)," "," "," ",round(sigma[1],4),round(sigma[4],4),round(sigma[5],4),round(sigma[6],4),round(sigma[7],4),round(sigma[8],4),round(sigma[9],4)," "," "," ",
                       round(t(mix_pi),4)," "," "," ",round(B211[1],4)," "," "," ",round(B211[2],4)," "," "," "," "," "," "," ",round(B211[3],4),round(B211[4],4),round(jj,4),round(ll*100,4),round(gg,4),round(rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
###################################MX2-CD-AD(E-5)##################################################################
G4F3ModelFun[[23]] <- function(K1,logL,df11,df21,df31,df41,G4F3text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF3 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF3)))
  sigmaf3<-var(dataF3);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a9<-sqrt(sigmaf3)
  m_nf <- as.numeric(G4F3text2)
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF3)))
  ##################procedure start#############################################
  d3<-9
  mi<-as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2))
  m[c(1:4)]<-m
  if (m[1]<m[3]) {a9<--a9}
  m[5]<-m[4]+2.4*a9;m[6]<-m[4]+1.6*a9;m[7]<-m[4]+0.8*a9
  m[8]<-m[4];m[9]<-m[4]-0.8*a9;m[10]<-m[4]-1.6*a9
  m[11]<-m[4]-2.4*a9;m[12]<-m[4]-3.2*a9; m[4]<-m[4]+3.2*a9
  pi<-m[c(4:12)];gh<-sigma[c(4:12)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mi)))
  ##################iteration process#############################
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,9,m_sam[4]); swx <- matrix(0,9,1)
  hh<-matrix(0,7,7);b_line<-matrix(0,7,1)
  g_a<-matrix(0,9,1);s0<-matrix(0,2,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:12)];ssigma<-sigma[c(4:12)]
    for(i in 1:d3) { WW[i,] <- mi[i]*dnorm(dataF3,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF3,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF3
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    #######################restrictions################################################
    hh[1,1]<-9.0*sigma[1]/m_sam[1]+4.0*sigma[1]/m_sam[2]+9.0*sigma[1]/m_sam[3]+4.0*sigma[4]/n0[1]+16.0*sigma[8]/n0[5]+4.0*sigma[12]/n0[9]
    hh[1,2]<--6.0*sigma[4]/n0[1]
    hh[1,3]<-16.0*sigma[8]/n0[5]
    hh[1,4]<-16.0*sigma[8]/n0[5]
    hh[1,5]<--2.0*sigma[4]/n0[1]
    hh[1,6]<--2.0*sigma[4]/n0[1]
    hh[1,7]<--4.0*sigma[8]/n0[5]-2.0*sigma[12]/n0[9]
    hh[2,2]<-9.0*sigma[4]/n0[1]+16.0*sigma[5]/n0[2]+sigma[6]/n0[3]
    hh[2,3]<-0.0
    hh[2,4]<--12.0*sigma[5]/n0[2]
    hh[2,5]<-3.0*sigma[4]/n0[1]+4.0*sigma[5]/n0[2]
    hh[2,6]<-3.0*sigma[4]/n0[1]-sigma[6]/n0[3]
    hh[2,7]<-0.0
    hh[3,3]<-9.0*sigma[7]/n0[4]+16.0*sigma[8]/n0[5]+sigma[9]/n0[6]
    hh[3,4]<-16.0*sigma[8]/n0[5]
    hh[3,5]<-0.0
    hh[3,6]<--3.0*sigma[7]/n0[4]+sigma[9]/n0[6]
    hh[3,7]<--4.0*sigma[8]/n0[5]-sigma[9]/n0[6]
    hh[4,4]<-9.0*sigma[5]/n0[2]+16.0*sigma[8]/n0[5]+sigma[11]/n0[8]
    hh[4,5]<--3.0*sigma[5]/n0[2]+sigma[11]/n0[8]
    hh[4,6]<-0.0
    hh[4,7]<--4.0*sigma[8]/n0[5]-sigma[11]/n0[8]
    hh[5,5]<-sigma[4]/n0[1]+sigma[5]/n0[2]+sigma[10]/n0[7]+sigma[11]/n0[8]
    hh[5,6]<-sigma[4]/n0[1]
    hh[5,7]<--sigma[11]/n0[8]
    hh[6,6]<-sigma[4]/n0[1]+sigma[6]/n0[3]+sigma[7]/n0[4]+sigma[9]/n0[6]
    hh[6,7]<--sigma[9]/n0[6]
    hh[7,7]<-sigma[8]/n0[5]+sigma[9]/n0[6]+sigma[11]/n0[8]+sigma[12]/n0[9]
    for(i in 2:7)
    {
      for(j in 1:(i-1))
      {
        hh[i,j]<-hh[j,i]
      }
    }
    ##########################################################################################################
    b_line[1]<-3*sumx[1]/m_sam[1]+2*sumx[2]/m_sam[2]+3*sumx[3]/m_sam[3]-2*sumwx[1]/n0[1]-4*sumwx[5]/n0[5]-2*sumwx[9]/n0[9]
    b_line[2]<-3*sumwx[1]/n0[1]-4*sumwx[2]/n0[2]+sumwx[3]/n0[3]
    b_line[3]<-3*sumwx[4]/n0[4]-4*sumwx[5]/n0[5]+sumwx[6]/n0[6]
    b_line[4]<-3*sumwx[2]/n0[2]-4*sumwx[5]/n0[5]+sumwx[8]/n0[8]
    b_line[5]<-sumwx[1]/n0[1]-sumwx[2]/n0[2]-sumwx[7]/n0[7]+sumwx[8]/n0[8]
    b_line[6]<-sumwx[1]/n0[1]-sumwx[3]/n0[3]-sumwx[4]/n0[4]+sumwx[6]/n0[6]
    b_line[7]<-sumwx[5]/n0[5]-sumwx[6]/n0[6]-sumwx[8]/n0[8]+sumwx[9]/n0[9]
    B22<-solve(hh,b_line)
    ##################################################################################################
    m[1]<-(sumx[1]-sigma[1]*B22[1]*3)/m_sam[1]
    m[2]<-(sumx[2]-sigma[1]*B22[1]*2)/m_sam[2]
    m[3]<-(sumx[3]-sigma[1]*B22[1]*3)/m_sam[3]
    m[4]<-(sumwx[1]+sigma[4]*(2*B22[1]-3*B22[2]-B22[5]-B22[6]))/n0[1]
    m[5]<-(sumwx[2]+sigma[5]*(4*B22[2]-3*B22[4]+B22[5]))/n0[2]
    m[6]<-(sumwx[3]+sigma[6]*(-B22[2]+B22[6]))/n0[3]
    m[7]<-(sumwx[4]+sigma[7]*(-3*B22[3]+B22[6]))/n0[4]
    m[8]<-(sumwx[5]+sigma[8]*(4*B22[1]+4*B22[3]+4*B22[4]-B22[7]))/n0[5]
    m[9]<-(sumwx[6]+sigma[9]*(-B22[3]-B22[6]+B22[7]))/n0[6]
    m[10]<-(sumwx[7]+sigma[10]*B22[5])/n0[7]
    m[11]<-(sumwx[8]+sigma[11]*(-B22[4]-B22[5]+B22[7]))/n0[8]
    m[12]<-(sumwx[9]+sigma[12]*(2*B22[1]-B22[7]))/n0[9]
    #######################first order genetic parameter process#########################
    hh221<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,-1,1,1,1,0.5,0.5,0.5,-1,-1,-1,
                    1,1,-1,1,0.5,-1,1,0.5,-1,1,0.5,-1,1,0,-1,0,0,0,0,0,0,0,0,0,
                    0,1,0,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25),12,5)
    B221<-solve(crossprod(hh221,hh221))%*%crossprod(hh221,m)
    g_a[c(1,3,7,9)]<-0
    g_a[2]<-0.75*B221[3]^2/m_nf
    g_a[4]<-0.75*B221[2]^2/m_nf
    g_a[5]<-(0.75*B221[2]^2+0.75*B221[3]^2)/m_nf
    g_a[6]<-0.75*B221[2]^2/m_nf
    g_a[8]<-0.75*B221[3]^2/m_nf
    ###################obtain variance##############################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d3) {swx[i] <- WW[i,]%*%(dataF3-m[i+3])^2 }
    ###########polygene variance process#####################
    a1<-sigma[4];n_iter<-0;a2<-1000
    while(a2>0.0001){
      n_iter<-n_iter+1
      h1<-sigma[4]/(sigma[4]+g_a)
      s0[1]<-sum(swx*h1^2)
      s0[2]<-sum(n0*h1)
      sigma[4]<-s0[1]/s0[2]
      a2<-abs(sigma[4]-a1)
      a1<-sigma[4]
      if(n_iter>20) break
    }

    sigma40<-sigma[4]-sigma[1]/m_nf

    if (sigma40<0){sigma40<-0;sigma[4]<-sigma[1]/m_nf}
    sigma[c(5:12)]<-sigma[4]+g_a[c(2:9)]
    ##########error process#####################
    a1<-sigma[1];n_iter<-0;a2<-sigma40;a2[a2<0]<-0;a4<-1000
    while(a4>0.0001){
      n_iter<-n_iter+1
      h1<-(sigma[1]/m_nf)/(a2+sigma[1]/m_nf+g_a)
      s0[1]<-sum(swx*h1^2*m_nf)
      s0[2]<-sum(h1*n0)
      s0[1]<-s0[1]+ss1+ss2+ss3
      s0[2]<-s0[2]+m_sam[1]+m_sam[2]+m_sam[3]
      sigma[1]<-s0[1]/s0[2]
      a4<-abs(sigma[1]-a1)
      a1<-sigma[1]
      if(n_iter>20) break
    }
    sigma[2]<-sigma[3]<-sigma[1]
    sigma[c(4:12)]<-sigma[1]/m_nf+sigma40+g_a[c(1:9)]
    ############criteria for iterations to stop################
    pi<-m[c(4:12)];gh<-sigma[c(4:12)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*7
  ###########second order genetic parameter process############
  jj<-sigmaf3 - sigma[4]
  gg<-sigma[4]-sigma[1]/m_nf
  jj[jj<0]<-0
  gg[gg<0]<-0
  ll<-jj/sigmaf3
  rr<-gg/sigmaf3
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma[1]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < m_sam[1]){P1bmw <- P1bmw+runif(m_sam[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:m_sam[1])) - 0.5)/m_sam[1])^2)
  P1u<- as.matrix(c(12*m_sam[1]*((P1dd[1]/m_sam[1]-0.5)^2),((45*m_sam[1])/4)*((P1dd[2]/m_sam[1]-1/3)^2),180*m_sam[1]*((P1dd[3]/m_sam[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,m_sam[1]))))

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ###################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*m_sam[2])
  F1bmw <- matrix(0,m_sam[2],1)
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma[2]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < m_sam[2]){F1bmw <- F1bmw+runif(m_sam[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:m_sam[2])) - 0.5)/m_sam[2])^2)
  F1u<- as.matrix(c(12*m_sam[2]*((F1dd[1]/m_sam[2]-0.5)^2),((45*m_sam[2])/4)*((F1dd[2]/m_sam[2]-1/3)^2),180*m_sam[2]*((F1dd[3]/m_sam[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,m_sam[2]))))

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*m_sam[3])
  P2bmw <- matrix(0,m_sam[3],1)
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma[3]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < m_sam[3]){P2bmw <- P2bmw+runif(m_sam[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:m_sam[3])) - 0.5)/m_sam[3])^2)
  P2u<- as.matrix(c(12*m_sam[3]*((P2dd[1]/m_sam[3]-0.5)^2),((45*m_sam[3])/4)*((P2dd[2]/m_sam[3]-1/3)^2),180*m_sam[3]*((P2dd[3]/m_sam[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,m_sam[3]))))

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3);
  F3w1<-1/(12*m_sam[4])
  F3bmw <- matrix(0,m_sam[4],1); F3bmwsl <- matrix(0,m_sam[4],d3)
  for(i in 1:d3){
    F3gg <- (dataF3 - pi[i])/sqrt(gh[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mix_pi[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[4]){F3P2 <- F3P2+runif(m_sam[4])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[4]) + sum((F3P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F3u <- as.matrix(c(12*m_sam[4]*((F3dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F3dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F3dd[3]/m_sam[4]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[4]))))

  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX2-CD-AD",round(abc,4),round(AIC,4),round(t(m),4),round(sigma[1],4),round(sigma[4],4),round(sigma[5],4),round(sigma[6],4),round(sigma[7],4),round(sigma[8],4),round(sigma[9],4),round(sigma[10],4),round(sigma[11],4),round(sigma[12],4),
                       round(t(mix_pi),4),round(B221[1],4)," "," "," ",round(B221[2],4),round(B221[3],4)," "," ", " "," "," "," ",round(B221[4],4),round(B221[5],4),round(jj,4),round(ll*100,4),round(gg,4),round(rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}
#################################MX2-EAD-AD(E-6)###############################################################
G4F3ModelFun[[24]] <- function(K1,logL,df11,df21,df31,df41,G4F3text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF3 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF3)))
  sigmaf3<-var(dataF3);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a9<-sqrt(sigmaf3)
  m_nf <- as.numeric(G4F3text2)
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF3)))
  ##################procedure start#############################################
  d2<-6
  mi<-as.matrix(c(0.0625,0.25,0.25,0.125,0.25,0.0625))
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2,sigmaf3/2))
  m[c(1:4)]<-m
  if(m[1]<m[3]) {a9<--a9}
  m[5]<-m[4]+2.1*a9;m[6]<-m[4]+1.2*a9;m[7]<-m[4];
  m[8]<-m[4]-0.7*a9;m[9]<-m[4]-3.0*a9;m[4]<-m[4]+3.0*a9;
  pi<-m[c(4:9)];gh<-sigma[c(4:9)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mi)))
  #################iteration process#########################################
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,6,m_sam[4]); swx <- matrix(0,6,1)
  hh<-matrix(0,5,5);b_line<-matrix(0,5,1)
  g_a<-matrix(0,6,1);s0<-matrix(0,2,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:9)];ssigma<-sigma[c(4:9)]
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataF3,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF3,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF3
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    #################################restrictions#########################################
    hh[1,1]<-9.0*sigma[1]/m_sam[1]+4.0*sigma[1]/m_sam[2]+9.0*sigma[1]/m_sam[3]+4.0*sigma[4]/n0[1]+16.0*sigma[6]/n0[3]+4.0*sigma[9]/n0[6]
    hh[1,2]<--2.0*sigma[4]/n0[1]-4.0*sigma[6]/n0[3]
    hh[1,3]<--2.0*sigma[4]/n0[1]-2.0*sigma[9]/n0[6]
    hh[1,4]<--2.0*sigma[4]/n0[1]+8.0*sigma[6]/n0[3]
    hh[1,5]<--2.0*sigma[9]/n0[6]
    hh[2,2]<-sigma[4]/n0[1]+4.0*sigma[5]/n0[2]+sigma[6]/n0[3]
    hh[2,3]<-sigma[4]/n0[1]
    hh[2,4]<-sigma[4]/n0[1]-2.0*sigma[6]/n0[3]
    hh[2,5]<-0.0
    hh[3,3]<-sigma[4]/n0[1]+4.0*sigma[7]/n0[4]+sigma[9]/n0[6]
    hh[3,4]<-sigma[4]/n0[1]-2.0*sigma[7]/n0[4]
    hh[3,5]<--6.0*sigma[7]/n0[4]+sigma[9]/n0[6]
    hh[4,4]<-sigma[4]/n0[1]+4.0*sigma[6]/n0[3]+sigma[7]/n0[4]
    hh[4,5]<-3.0*sigma[7]/n0[4]
    hh[5,5]<-9.0*sigma[7]/n0[4]+16.0*sigma[8]/n0[5]+sigma[9]/n0[6]
    for(i in 2:5)
    {
      for(j in 1:(i-1))
      {
        hh[i,j]<-hh[j,i]
      }
    }
    ############################################################################
    b_line[1]<-3.0*sumx[1]/m_sam[1]+2.0*sumx[2]/m_sam[2]+3.0*sumx[3]/m_sam[3]-2.0*(sumwx[1]/n0[1]+2.0*sumwx[3]/n0[3]+sumwx[6]/n0[6])
    b_line[2]<-sumwx[1]/n0[1]-2.0*sumwx[2]/n0[2]+sumwx[3]/n0[3]
    b_line[3]<-sumwx[1]/n0[1]-2.0*sumwx[4]/n0[4]+sumwx[6]/n0[6]
    b_line[4]<-sumwx[1]/n0[1]-2.0*sumwx[3]/n0[3]+sumwx[4]/n0[4]
    b_line[5]<-3.0*sumwx[4]/n0[4]-4.0*sumwx[5]/n0[5]+sumwx[6]/n0[6]
    B231<-solve(hh,b_line)
    #########################################################################################
    m[1]<-(sumx[1]-sigma[1]*B231[1]*3.0)/m_sam[1]
    m[2]<-(sumx[2]-sigma[1]*B231[1]*2.0)/m_sam[2]
    m[3]<-(sumx[3]-sigma[1]*B231[1]*3.0)/m_sam[3]
    m[4]<-(sumwx[1]+sigma[4]*(2.0*B231[1]-B231[2]-B231[3]-B231[4]))/n0[1]
    m[5]<-(sumwx[2]+sigma[5]*2.0*B231[2])/n0[2]
    m[6]<-(sumwx[3]+sigma[6]*(4.0*B231[1]-B231[2]+2.0*B231[4]))/n0[3]
    m[7]<-(sumwx[4]+sigma[7]*(2.0*B231[3]-B231[4]-3.0*B231[5]))/n0[4]
    m[8]<-(sumwx[5]+sigma[8]*4.0*B231[5])/n0[5]
    m[9]<-(sumwx[6]+sigma[9]*(2.0*B231[1]-B231[3]-B231[5]))/n0[6]
    ######################first order genetic parameter process#######################
    hh231<-matrix(c(1,1,1,1,1,1,1,1,1,2,2,-2,2,1.5,1,0,-0.5,-2,1,0,-1,0,0,0,0,0,0,0,1,0,0.25,0.25,0.25,0.25,0.25,0.25),9,4)
    B231<-solve(crossprod(hh231,hh231))%*%crossprod(hh231,m)
    g_a[c(1,4,6)]<-0
    g_a[2]<-0.75*B231[2]^2/m_nf;
    g_a[3]<-1.50*B231[2]^2/m_nf;
    g_a[5]<-0.75*B231[2]^2/m_nf;
    ##############obtain variance###########################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d2) {swx[i] <- WW[i,]%*%(dataF3-m[i+3])^2 }
    #########polygene variance process########################
    a1<-sigma[4];n_iter<-0;a2<-1000
    while(a2>0.0001){
      n_iter<-n_iter+1
      h1<-sigma[4]/(sigma[4]+g_a)
      s0[1]<-sum(swx*h1^2)
      s0[2]<-sum(n0*h1)
      sigma[4]<-s0[1]/s0[2]
      a2<-abs(sigma[4]-a1)
      a1<-sigma[4]
      if(n_iter>20) break

    }

    sigma40<-sigma[4]-sigma[1]/m_nf

    if (sigma40<0){sigma40<-0;sigma[4]<-sigma[1]/m_nf}
    sigma[c(5:9)]<-sigma[4]+g_a[c(2:6)]
    #######################error process###############
    a1<-sigma[1];n_iter<-0.0;a2<-sigma40;a2[a2<0]<-0;a4<-1000
    while(a4>0.0001){
      n_iter<-n_iter+1
      h1<-(sigma[1]/m_nf)/(a2+sigma[1]/m_nf+g_a)
      s0[1]<-sum(swx*h1^2*m_nf)
      s0[2]<-sum(h1*n0)
      s0[1]<-s0[1]+ss1+ss2+ss3
      s0[2]<-s0[2]+m_sam[1]+m_sam[2]+m_sam[3]
      sigma[1]<-s0[1]/s0[2]
      a4<-abs(sigma[1]-a1)
      a1<-sigma[1]
      if(n_iter>20) break
    }
    sigma[2]<-sigma[3]<-sigma[1]
    sigma[c(4:9)]<-sigma[1]/m_nf+sigma40+g_a[c(1:6)]
    ###############criteria for iterations to stop#################
    pi<-m[c(4:9)];gh<-sigma[c(4:9)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF3,pi,sqrt(gh),mix_pi)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*6
  ###########second order genetic parameter process############
  jj<-sigmaf3 - sigma[4]
  gg<-sigma[4]-sigma[1]/m_nf
  jj[jj<0]<-0
  gg[gg<0]<-0
  ll<-jj/sigmaf3
  rr<-gg/sigmaf3
  #####################################hypothesis testing for P1########################################
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma[1]))
  P1bmw[which(P1gg>=0)] <- pnorm(P1gg[P1gg>=0])
  P1bmw[which(P1gg<0)] <- 1 - pnorm(abs(P1gg[P1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P1bmw)))[1]
  if(nn < m_sam[1]){P1bmw <- P1bmw+runif(m_sam[1])/1e4}
  #########################################################
  P1dd<-c((sum(P1bmw)),(sum(P1bmw^2)),sum((P1bmw-0.5)^2))
  P1w<-P1w1+sum((P1bmw - (as.matrix(c(1:m_sam[1])) - 0.5)/m_sam[1])^2)
  P1u<- as.matrix(c(12*m_sam[1]*((P1dd[1]/m_sam[1]-0.5)^2),((45*m_sam[1])/4)*((P1dd[2]/m_sam[1]-1/3)^2),180*m_sam[1]*((P1dd[3]/m_sam[1]-1/12)^2)))
  P1D<-as.numeric(ks.test(P1bmw,"punif")[[1]][1])
  P1tt <- as.matrix(c((1 - pchisq(P1u[1],1)),(1 - pchisq(P1u[2],1)),(1 - pchisq(P1u[3],1)),K1(P1w),(1-pkolm(P1D,m_sam[1]))))

  P1tt[which( P1tt>=10e-4)]<-round(P1tt[which(P1tt>=10e-4)],4);P1tt[which(P1tt<10e-4)]<-format(P1tt[which(P1tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F1#########################################
  dataF1<-sort(dataF1)
  F1w1<-1/(12*m_sam[2])
  F1bmw <- matrix(0,m_sam[2],1)
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma[2]))
  F1bmw[which(F1gg>=0)] <- pnorm(F1gg[F1gg>=0])
  F1bmw[which(F1gg<0)] <- 1 - pnorm(abs(F1gg[F1gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F1bmw)))[1]
  if(nn < m_sam[2]){F1bmw <- F1bmw+runif(m_sam[2])/1e4}
  ##########################################################
  F1dd<-c((sum(F1bmw)),(sum(F1bmw^2)),sum((F1bmw-0.5)^2))
  F1w<-F1w1+sum((F1bmw - (as.matrix(c(1:m_sam[2])) - 0.5)/m_sam[2])^2)
  F1u<- as.matrix(c(12*m_sam[2]*((F1dd[1]/m_sam[2]-0.5)^2),((45*m_sam[2])/4)*((F1dd[2]/m_sam[2]-1/3)^2),180*m_sam[2]*((F1dd[3]/m_sam[2]-1/12)^2)))
  F1D<-as.numeric(ks.test(F1bmw,"punif")[[1]][1])
  F1tt <- as.matrix(c((1 - pchisq(F1u[1],1)),(1 - pchisq(F1u[2],1)),(1 - pchisq(F1u[3],1)),K1(F1w),(1-pkolm(F1D,m_sam[2]))))

  F1tt[which(F1tt>=10e-4)]<-round(F1tt[which(F1tt>=10e-4)],4);F1tt[which(F1tt<10e-4)]<-format(F1tt[which(F1tt<10e-4)],scientific=TRUE,digit=4)

  #################################hypothesis testing for P2#############################################
  dataP2<-sort(dataP2)
  P2w1<-1/(12*m_sam[3])
  P2bmw <- matrix(0,m_sam[3],1)
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma[3]))
  P2bmw[which(P2gg>=0)] <- pnorm(P2gg[P2gg>=0])
  P2bmw[which(P2gg<0)] <- 1 - pnorm(abs(P2gg[P2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(P2bmw)))[1]
  if(nn < m_sam[3]){P2bmw <- P2bmw+runif(m_sam[3])/1e4}
  ##########################################################
  P2dd<-c((sum(P2bmw)),(sum(P2bmw^2)),sum((P2bmw-0.5)^2))
  P2w<-P2w1+sum((P2bmw - (as.matrix(c(1:m_sam[3])) - 0.5)/m_sam[3])^2)
  P2u<- as.matrix(c(12*m_sam[3]*((P2dd[1]/m_sam[3]-0.5)^2),((45*m_sam[3])/4)*((P2dd[2]/m_sam[3]-1/3)^2),180*m_sam[3]*((P2dd[3]/m_sam[3]-1/12)^2)))
  P2D<-as.numeric(ks.test(P2bmw,"punif")[[1]][1])
  P2tt <- as.matrix(c((1 - pchisq(P2u[1],1)),(1 - pchisq(P2u[2],1)),(1 - pchisq(P2u[3],1)),K1(P2w),(1-pkolm(P2D,m_sam[3]))))

  P2tt[which(P2tt>=10e-4)]<-round(P2tt[which(P2tt>=10e-4)],4);P2tt[which(P2tt<10e-4)]<-format(P2tt[which(P2tt<10e-4)],scientific=TRUE,digit=4)

  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3);
  F3w1<-1/(12*m_sam[4])
  F3bmw <- matrix(0,m_sam[4],1); F3bmwsl <- matrix(0,m_sam[4],d2)
  for(i in 1:d2){
    F3gg <- (dataF3 - pi[i])/sqrt(gh[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mix_pi[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[4]){F3P2 <- F3P2+runif(m_sam[4])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[4]) + sum((F3P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F3u <- as.matrix(c(12*m_sam[4]*((F3dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F3dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F3dd[3]/m_sam[4]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[4]))))

  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)

  output <- data.frame("MX2-EAD-AD",round(abc,4),round(AIC,4),round(t(m),4)," "," "," ",round(sigma[1],4),round(sigma[4],4),round(sigma[5],4),round(sigma[6],4),round(sigma[7],4),round(sigma[8],4),round(sigma[9],4)," "," "," ",
                       round(t(mix_pi),4)," "," "," ",round(B231[1],4)," "," "," ",round(B231[2],4)," "," "," "," "," "," "," ",round(B231[3],4),round(B231[4],4),round(jj,4),round(ll*100,4),round(gg,4),round(rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}


K1G4F3 <- function(x){
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

logLG4F3 <- function(nm,nng,mi,mn,s,d1) { sum2 <- sum(log(dmixnorm(d1,mn,sqrt(s),mi)));return (sum2) }



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
    G4F3ModelFun[[i]](K1G4F3,logLG4F3,df11,df21,df31,df41,G4F3text2)[[1]]
  }
  stopCluster(cl)
  mi<-NULL

}else{

  allresultq=switch(model,"1MG-AD" = G4F3ModelFun[[1]](K1G4F3,logLG4F3,df11,df21,df31,df41,G4F3text2),"1MG-A"=G4F3ModelFun[[2]](K1G4F3,logLG4F3,df11,df21,df31,df41,G4F3text2),"1MG-EAD"=G4F3ModelFun[[3]](K1G4F3,logLG4F3,df11,df21,df31,df41,G4F3text2),"1MG-NCD"=G4F3ModelFun[[4]](K1G4F3,logLG4F3,df11,df21,df31,df41,G4F3text2),"2MG-ADI"=G4F3ModelFun[[5]](K1G4F3,logLG4F3,df11,df21,df31,df41,G4F3text2),
                    "2MG-AD"=G4F3ModelFun[[6]](K1G4F3,logLG4F3,df11,df21,df31,df41,G4F3text2),"2MG-A"=G4F3ModelFun[[7]](K1G4F3,logLG4F3,df11,df21,df31,df41,G4F3text2),"2MG-EA"=G4F3ModelFun[[8]](K1G4F3,logLG4F3,df11,df21,df31,df41,G4F3text2),"2MG-CD"=G4F3ModelFun[[9]](K1G4F3,logLG4F3,df11,df21,df31,df41,G4F3text2),"2MG-EAD"=G4F3ModelFun[[10]](K1G4F3,logLG4F3,df11,df21,df31,df41,G4F3text2),
                    "PG-ADI"=G4F3ModelFun[[11]](K1G4F3,logLG4F3,df11,df21,df31,df41,G4F3text2),"PG-AD"=G4F3ModelFun[[12]](K1G4F3,logLG4F3,df11,df21,df31,df41,G4F3text2),"MX1-AD-ADI"=G4F3ModelFun[[13]](K1G4F3,logLG4F3,df11,df21,df31,df41,G4F3text2),"MX1-AD-AD"=G4F3ModelFun[[14]](K1G4F3,logLG4F3,df11,df21,df31,df41,G4F3text2),"MX1-A-AD"=G4F3ModelFun[[15]](K1G4F3,logLG4F3,df11,df21,df31,df41,G4F3text2),
                    "MX1-EAD-AD"=G4F3ModelFun[[16]](K1G4F3,logLG4F3,df11,df21,df31,df41,G4F3text2),"MX1-NCD-AD"=G4F3ModelFun[[17]](K1G4F3,logLG4F3,df11,df21,df31,df41,G4F3text2),"MX2-ADI-ADI"=G4F3ModelFun[[18]](K1G4F3,logLG4F3,df11,df21,df31,df41,G4F3text2),"MX2-ADI-AD"=G4F3ModelFun[[19]](K1G4F3,logLG4F3,df11,df21,df31,df41,G4F3text2),"MX2-AD-AD"=G4F3ModelFun[[20]](K1G4F3,logLG4F3,df11,df21,df31,df41,G4F3text2),
                    "MX2-A-AD"=G4F3ModelFun[[21]](K1G4F3,logLG4F3,df11,df21,df31,df41,G4F3text2),"MX2-EA-AD"=G4F3ModelFun[[22]](K1G4F3,logLG4F3,df11,df21,df31,df41,G4F3text2),"MX2-CD-AD"=G4F3ModelFun[[23]](K1G4F3,logLG4F3,df11,df21,df31,df41,G4F3text2),"MX2-EAD-AD"=G4F3ModelFun[[24]](K1G4F3,logLG4F3,df11,df21,df31,df41,G4F3text2))

  allresult<-allresultq[[1]]
  if(model=="PG-AD"||model=="PG-ADI"){
    mi<-NULL
  }else{
    mi<-allresultq[[2]]
  }
}
colnames(allresult) <- G4F3colname
out<-list(allresult,mi)
return(out)
}




