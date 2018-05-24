
G4F2Fun<-function(df,model){

data<-sapply(df,as.character)
dP1<-data[-1,which(data[1,]=="P1")];P1<-as.numeric(dP1[which(is.na(as.numeric(dP1))==FALSE)]);df11<-as.data.frame(P1)
dF1<-data[-1,which(data[1,]=="F1")];F1<-as.numeric(dF1[which(is.na(as.numeric(dF1))==FALSE)]);df21<-as.data.frame(F1)
dP2<-data[-1,which(data[1,]=="P2")];P2<-as.numeric(dP2[which(is.na(as.numeric(dP2))==FALSE)]);df31<-as.data.frame(P2)
dF2<-data[-1,which(data[1,]=="F2")];F2<-as.numeric(dF2[which(is.na(as.numeric(dF2))==FALSE)]);df41<-as.data.frame(F2)
##################################################

G4F2colname <- c("Model","Log_Max_likelihood_Value","AIC","meanP1","meanF1","meanP2","mean[1]","mean[2]","mean[3]","mean[4]","mean[5]","mean[6]","mean[7]","mean[8]","mean[9]",
                 "Var(Residual)","Var(Component)","Proportion[1]","Proportion[2]","Proportion[3]","Proportion[4]","Proportion[5]","Proportion[6]","Proportion[7]","Proportion[8]","Proportion[9]",
                 "m1(m)","m2","m3","m4","da(d)","db","ha(h)","hb","i","jab","jba","l","[d]","[h]","Major-Gene Var","Heritability(Major-Gene)(%)","Poly-Gene Var","Heritability(Poly-Gene)(%)","U1 square(P1)","P(U1 square(P1))","U2 square(P1)","P(U2 square(P1))","U3 square(P1)","P(U3 square(P1))","nW(P1)","P(nW(P1))","Dn(P1)","P(Dn(P1))",
                 "U1 square(F1)","P(U1 square(F1))","U2 square(F1)","P(U2 square(F1))","U3 square(F1)","P(U3 square(F1))","nW(F1)","P(nW(F1))","Dn(F1)","P(Dn(F1))","U1 square(P2)","P(U1 square(P2))","U2 square(P2)","P(U2 square(P2))","U3 square(P2)","P(U3 square(P2))","nW(P2)","P(nW(P2))","Dn(P2)","P(Dn(P2))",
                 "U1 square(F2)","P(U1 square(F2))","U2 square(F2)","P(U2 square(F2))","U3 square(F2)","P(U3 square(F2))","nW(F2)","P(nW(F2))","Dn(F2)","P(Dn(F2))")
G4F2ModelFun<-list(NA)
###################define each model function############################
##########################1MG_AD(A_1)#################################
G4F2ModelFun[[1]] <- function(K1,logL,df11,df21,df31,df41){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF2 <- as.matrix(as.numeric(df41[,1]))
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2)))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1])
  m_esp<-0.0001
  sigmaF2<-var(dataF2);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a1<-sqrt(sigmaF2/m_sam[4])
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2)))
  ###############procedure start##############################
  d1<-3
  mi<-as.matrix(c(0.25,0.5,0.25))
  sigma<-matrix((sigmaF2/2),6,1)
  m[c(1:4)]<-m
  if(m[1]<m[3]){a1<--a1}
  m[5]<-m[4];m[4]<-m[5]+2*a1;m[6]<-m[5]-2*a1
  pi<-m[c(4:6)];gh<-sigma[c(4:6)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mi)))
  ##############iteration process###############
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,3,m_sam[4]); swx <- matrix(0,3,1);s0<-matrix(0,1,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:6)];ssigma<-sigma[c(4:6)]
    for(i in 1:d1) { WW[i,] <- mi[i]*dnorm(dataF2,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF2,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF2
    sigma[4]<-sigma[5]<- sigma[6]<-sigma[1]
    ##############obtain mean###################
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    m[c(1:3)]<-(sumx[c(1:3)]+sumwx[c(1:3)])/(m_sam[c(1:3)]+n0[c(1:3)])
    m[c(4:6)]<-m[c(1:3)]
    #############obtain variance#################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d1) {swx[i] <- WW[i,]%*%(dataF2-m[i+3])^2 }
    s0[1]<-ss1+ss2+ss3+swx[1]+swx[2]+swx[3]
    sigma[1]<-s0[1]/(m_sam[1]+m_sam[2]+m_sam[3]+m_sam[4])
    sigma[2]<-sigma[3]<-sigma[4]<-sigma[5]<-sigma[6]<-sigma[1]
    ########criteria for iterations to stop####################
    pi<-m[c(4:6)];gh<-sigma[c(4:6)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mix_pi))) 
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0;
  AIC<--2*abc+2*6
  #########first order genetic parameter process##########
  hh1 <- matrix(c(1,1,1,1,0,-1,0,1,0),3,3)
  B1<-solve(hh1,m[c(4:6)])
  #################second oder genetic parameter process#############
  jj<-sigmaF2-sigma[1]
  if(jj<0){jj<-0}
  ll<-jj/sigmaF2
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
  
  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2); 
  F2w1<-1/(12*m_sam[4])
  F2bmw <- matrix(0,m_sam[4],1); F2bmwsl <- matrix(0,m_sam[4],d1)
  for(i in 1:d1){
    F2gg <- (dataF2 - pi[i])/sqrt(gh[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < m_sam[4]){F2P2 <- F2P2+runif(m_sam[4])/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*m_sam[4]) + sum((F2P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F2u <- as.matrix(c(12*m_sam[4]*((F2dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F2dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F2dd[3]/m_sam[4]-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),(1-pkolm(F2D,m_sam[4]))))
  
  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("1MG-AD",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," "," ",round(sigma[1],4),round(sigma[4],4),
                       round(t(mix_pi),4)," "," "," "," "," "," ",round(B1[1],4)," "," "," ",round(B1[2],4)," ",round(B1[3],4)," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}

################1MG-A(A_2)#########################     
G4F2ModelFun[[2]] <- function(K1,logL,df11,df21,df31,df41){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF2 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2)))
  m_esp<-0.0001
  sigmaF2<-var(dataF2);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a1<-sqrt(sigmaF2/m_sam[4])
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2)))
  ###############procedure start##############################    
  mi<-as.matrix(c(0.25,0.5,0.25))
  sigma<-matrix((sigmaF2/2),6,1)
  m[c(1:4)]<-m
  if(m[1]<m[3]){a1<--a1}
  m[5]<-m[4] ;m[4]<-m[5]+2*a1;m[6]<-m[5]-2*a1
  pi<-m[c(4:6)];gh<-sigma[c(4:6)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mi)))
  ########iteration process############
  d1<-3
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,3,m_sam[4]); swx <- matrix(0,3,1); s0<-matrix(0,1,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:6)];ssigma<-sigma[c(4:6)]
    for(i in 1:d1) { WW[i,] <- mi[i]*dnorm(dataF2,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF2,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF2
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    #########restriction##########
    aa1<-(sumx[1]+sumwx[1])/(m_sam[1]+n0[1])-2*(sumx[2]+sumwx[2])/(m_sam[2]+n0[2])+(sumx[3]+sumwx[3])/(m_sam[3]+n0[3])
    aa2<-sigma[1]/(m_sam[1]+n0[1])+4*sigma[1]/(m_sam[2]+n0[2])+sigma[1]/(m_sam[3]+n0[3])
    aa3<-aa1/aa2
    #######obtain mean#####################
    m[1]<-(sumx[1]+sumwx[1]-sigma[1]*aa3)/(m_sam[1]+n0[1])
    m[2]<-(sumx[2]+sumwx[2]+sigma[2]*aa3*2)/(m_sam[2]+n0[2])
    m[3]<-(sumx[3]+sumwx[3]-sigma[3]*aa3)/(m_sam[3]+n0[3])
    m[c(4:6)]<-m[c(1:3)]
    ######obtain variance##############
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2); ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d1) {swx[i] <- WW[i,]%*%(dataF2-m[i+3])^2 }
    s0[1]<-ss1+ss2+ss3+ swx[1]+ swx[2]+ swx[3]
    sigma[1]<-s0[1]/(m_sam[1]+m_sam[2]+m_sam[3]+m_sam[4])
    sigma[c(2:6)]<-sigma[1]
    ########criteria for iterations to stop####################
    pi<-m[c(4:6)];gh<-sigma[c(4:6)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mix_pi))) 
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*5  
  #########first order genetic parameter process##########
  hh2<-matrix(c(1,1,1,1,0,-1),3,2)
  B2<-solve(crossprod(hh2,hh2))%*%crossprod(hh2,m[c(4:6)])
  #################second oder genetic parameter process#############
  jj<-sigmaF2-sigma[1]
  if(jj<0){jj<-0}
  ll<-jj/sigmaF2
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
  
  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2); 
  F2w1<-1/(12*m_sam[4])
  F2bmw <- matrix(0,m_sam[4],1); F2bmwsl <- matrix(0,m_sam[4],d1)
  for(i in 1:d1){
    F2gg <- (dataF2 - pi[i])/sqrt(gh[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < m_sam[4]){F2P2 <- F2P2+runif(m_sam[4])/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*m_sam[4]) + sum((F2P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F2u <- as.matrix(c(12*m_sam[4]*((F2dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F2dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F2dd[3]/m_sam[4]-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),(1-pkolm(F2D,m_sam[4]))))
  
  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("1MG-A",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," "," ",round(sigma[1],4),round(sigma[4],4),
                       round(t(mix_pi),4)," "," "," "," "," "," ",round(B2[1],4)," "," "," ",round(B2[2],4)," "," "," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}

##################1MG-EAD(A_3)#####################################
G4F2ModelFun[[3]] <- function(K1,logL,df11,df21,df31,df41){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF2 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2)))
  m_esp<-0.0001
  sigmaF2<-var(dataF2);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a1<-sqrt(sigmaF2/m_sam[4])
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2)))
  ###############procedure start##############################   
  d2<-2
  mi<-as.matrix(c(0.75,0.25))
  m[c(1:4)]<-m
  if(m[1]<m[3]){a1<--a1}
  m[5]<-m[4]-1.5*a1 ;m[4]<-m[4]+1.5*a1
  sigma<-matrix((sigmaF2/2),5,1)
  pi<-m[c(4:5)];gh<-sigma[c(4:5)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mi)))
  ########iteration process############
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,2,m_sam[4]); swx <- matrix(0,2,1);  s0<-matrix(0,1,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4,5)];ssigma<-sigma[c(4,5)]
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataF2,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF2,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF2
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    ###########obtain mean##############################
    m[1]<-(sumx[1]+sumx[2]+sumwx[1])/(m_sam[1]+m_sam[2]+n0[1])
    m[3]<-(sumx[3]+sumwx[2])/(m_sam[3]+n0[2]);
    m[2]<-m[1];m[4]<-m[1];m[5]<-m[3]
    ######obtain variance##############
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d2) {swx[i] <- WW[i,]%*%(dataF2-m[i+3])^2 }
    ########criteria for iterations to stop####################
    s0[1]<-ss1+ss2+ss3+ swx[1]+ swx[2]
    sigma[1]<-s0[1]/(m_sam[1]+m_sam[2]+m_sam[3]+m_sam[4])
    sigma[c(2:5)]<-sigma[1]
    ########criteria for iterations to stop####################
    pi<-m[c(4:5)];gh<-sigma[c(4:5)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mix_pi))) 
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*5  
  #########first order genetic parameter process##########
  hh3<-matrix(c(1,1,1,-1),2,2)
  B3<-solve(hh3,m[c(4,5)])
  #################second oder genetic parameter process#############
  jj<-sigmaF2-sigma[1]
  if(jj<0){jj<-0}
  ll<-jj/sigmaF2
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
  
  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2); 
  F2w1<-1/(12*m_sam[4])
  F2bmw <- matrix(0,m_sam[4],1); F2bmwsl <- matrix(0,m_sam[4],d2)
  for(i in 1:d2){
    F2gg <- (dataF2 - pi[i])/sqrt(gh[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < m_sam[4]){F2P2 <- F2P2+runif(m_sam[4])/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*m_sam[4]) + sum((F2P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F2u <- as.matrix(c(12*m_sam[4]*((F2dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F2dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F2dd[3]/m_sam[4]-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),(1-pkolm(F2D,m_sam[4]))))
  
  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("1MG-EAD",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," "," "," ",round(sigma[1],4),round(sigma[4],4),
                       round(t(mix_pi),4)," "," "," "," "," "," "," ",round(B3[1],4)," "," "," ",round(B3[2],4)," "," "," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}

###################1MG-NCD(A_4)#################################
G4F2ModelFun[[4]] <- function(K1,logL,df11,df21,df31,df41){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF2 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2)))
  m_esp<-0.0001
  sigmaF2<-var(dataF2);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a1<-sqrt(sigmaF2/m_sam[4])
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2)))
  ###############procedure start##############################   
  d2<-2
  mi<-as.matrix(c(0.25,0.75))
  m[c(1:4)]<-m
  if(m[1]<m[3]){a1<--a1}
  m[5]<-m[4]-1.5*a1 ;m[4]<-m[4]+1.5*a1
  sigma<-matrix((sigmaF2/2),5,1)
  pi<-m[c(4:5)];gh<-sigma[c(4:5)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mi)))
  ########iteration process############
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,2,m_sam[4]); swx <- matrix(0,2,1); s0<-matrix(0,1,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4,5)];ssigma<-sigma[c(4,5)]
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataF2,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF2,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF2
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    ##########obtain mean####################
    m[1]<-(sumx[1]+sumwx[1])/(m_sam[1]+n0[1])
    m[2]<-(sumx[2]+sumx[3]+sumwx[2])/(m_sam[2]+m_sam[3]+n0[2])
    m[3]<-m[2];m[4]<-m[1];m[5]<-m[2]
    ######obtain variance##############
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d2) {swx[i] <- WW[i,]%*%(dataF2-m[i+3])^2 }
    s0[1]<-ss1+ss2+ss3+swx[1]+swx[2]
    sigma[1]<-s0[1]/(m_sam[1]+m_sam[2]+m_sam[3]+m_sam[4])
    sigma[c(2:5)]<-sigma[1]
    ########criteria for iterations to stop####################
    pi<-m[c(4:5)];gh<-sigma[c(4:5)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mix_pi))) 
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*5  
  #########first order genetic parameter process##########
  hh4<-matrix(c(1,1,1,-1),2,2)
  B4<-solve(hh4,m[c(4,5)])
  #################second oder genetic parameter process#############
  jj<-sigmaF2-sigma[1]
  if(jj<0){jj<-0}
  ll<-jj/sigmaF2
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
  
  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2); 
  F2w1<-1/(12*m_sam[4])
  F2bmw <- matrix(0,m_sam[4],1); F2bmwsl <- matrix(0,m_sam[4],d2)
  for(i in 1:d2){
    F2gg <- (dataF2 - pi[i])/sqrt(gh[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < m_sam[4]){F2P2 <- F2P2+runif(m_sam[4])/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*m_sam[4]) + sum((F2P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F2u <- as.matrix(c(12*m_sam[4]*((F2dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F2dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F2dd[3]/m_sam[4]-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),(1-pkolm(F2D,m_sam[4])))) 
  
  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("1MG-NCD",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," "," "," ",round(sigma[1],4),round(sigma[4],4),
                       round(t(mix_pi),4)," "," "," "," "," "," "," ",round(B4[1],4)," "," "," ",round(B4[2],4)," "," "," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}

################2MG-ADI(B-1)####################################### 
G4F2ModelFun[[5]] <- function(K1,logL,df11,df21,df31,df41){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF2 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2)))
  m_esp<-0.0001
  sigmaF2<-var(dataF2);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a2<-sqrt(sigmaF2/m_sam[4]*(m_sam[4]-1))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2)))
  ###############procedure start##############################  
  d3<-9
  mi<-as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  sigma<-matrix((sigma0),12,1)
  m[c(1:3)]<-m[c(1:3)]
  if (m[1]<m[3]){a2<--a2}
  m[5]<-m[4]+2.4*a2;m[6]<-m[4]+1.6*a2;m[7]<-m[4]+0.8*a2
  m[8]<-m[4];m[9]<-m[4]-0.8*a2;m[10]<-m[4]-1.6*a2
  m[11]<-m[4]-2.4*a2;m[12]<-m[4]-3.2*a2;m[4]<-m[4]+3.2*a2
  pi<-m[c(4:12)];gh<-sigma[c(4:12)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mi)))
  ########iteration process############
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,9,m_sam[4]);swx <- matrix(0,9,1); s0<-matrix(0,1,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:12)];ssigma<-sigma[c(4:12)]
    for(i in 1:d3) { WW[i,] <- mi[i]*dnorm(dataF2,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF2,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF2
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    #########obtain mean#############
    m[1]<-(sumx[1]+sumwx[1])/(m_sam[1]+n0[1]);m[2]<-(sumx[2]+sumwx[5])/(m_sam[2]+n0[5])
    m[3]<-(sumx[3]+sumwx[9])/(m_sam[3]+n0[9]);m[5]<-sumwx[2]/n0[2];m[6]<-sumwx[3]/n0[3]
    m[7]<-sumwx[4]/n0[4];m[9]<-sumwx[6]/n0[6];m[10]<-sumwx[7]/n0[7];m[11]<-sumwx[8]/n0[8]
    m[4]<-m[1];m[8]<-m[2];m[12]<-m[3]
    ######obtain variance##############
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d3) {swx[i] <- WW[i,]%*%(dataF2-m[i+3])^2 }
    s0[1]<-sum(swx);s0[1]<-s0[1]+ss1+ss2+ss3
    sigma[1]<-s0[1]/(m_sam[1]+m_sam[2]+m_sam[3]+m_sam[4])
    sigma[c(2:12)]<-sigma[1]
    ########criteria for iterations to stop####################
    pi<-m[c(4:12)];gh<-sigma[c(4:12)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mix_pi))) 
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*11
  #############first order genetic parameter process##############
  hh5<-matrix(c(1,1,1,1,1,1,1,1,1,1,0,-1,1,1,0,0,-1,-1,1,0,-1,0,-1,1,-1,1,0,0,1,0,0,0,1,1,0,0,
                0,1,0,1,0,0,0,0,1,1,0,1,0,-1,0,0,-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,0,1,-1,0,0,0,
                1,0,0,0,0,0,0,0),9,9)
  B5<-solve(hh5,m[c(1,2,3,5,6,7,9,10,11)])
  #################second oder genetic parameter process#############
  jj<-sigmaF2-sigma[1]
  if(jj<0){jj<-0}
  ll<-jj/sigmaF2
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
  
  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2); 
  F2w1<-1/(12*m_sam[4])
  F2bmw <- matrix(0,m_sam[4],1); F2bmwsl <- matrix(0,m_sam[4],d3)
  for(i in 1:d3){
    F2gg <- (dataF2 - pi[i])/sqrt(gh[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < m_sam[4]){F2P2 <- F2P2+runif(m_sam[4])/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*m_sam[4]) + sum((F2P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F2u <- as.matrix(c(12*m_sam[4]*((F2dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F2dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F2dd[3]/m_sam[4]-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),(1-pkolm(F2D,m_sam[4]))))
  
  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("2MG-ADI",round(abc,4),round(AIC,4),round(t(m),4),round(sigma[1],4),round(sigma[4],4),
                       round(t(mix_pi),4),round(B5[1],4)," "," "," ",round(B5[2],4),round(B5[3],4),round(B5[4],4),round(B5[5],4),round(B5[6],4),round(B5[7],4),round(B5[8],4),round(B5[9],4)," "," ",round(jj,4),round(ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}

####################2MG-AD(B-2)#################################### 
G4F2ModelFun[[6]] <- function(K1,logL,df11,df21,df31,df41){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF2 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2)))
  m_esp<-0.0001
  sigmaF2<-var(dataF2);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a2<-sqrt(sigmaF2/m_sam[4]*(m_sam[4]-1))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2)))
  ###############procedure start##############################   
  d3<-9
  mi<-as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  sigma<-matrix((sigma0),12,1)
  m[c(1:3)]<-m[c(1:3)]
  if(m[1]<m[3]){a2<--a2}
  m[5]<-m[4]+2.4*a2;m[6]<-m[4]+1.6*a2;m[7]<-m[4]+0.8*a2
  m[8]<-m[4];m[9]<-m[4]-0.8*a2;m[10]<-m[4]-1.6*a2
  m[11]<-m[4]-2.4*a2;m[12]<-m[4]-3.2*a2;m[4]<-m[4]+3.2*a2
  pi<-m[c(4:12)];gh<-sigma[c(4:12)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mi)))
  ########iteration process############
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,9,m_sam[4]);swx <- matrix(0,9,1);s0<-matrix(0,1,1)
  hh6<-matrix(0,4,4);b_line6<-matrix(0,4,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:12)];ssigma<-sigma[c(4:12)]
    for(i in 1:d3) { WW[i,] <- mi[i]*dnorm(dataF2,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF2,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF2
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    #####################################################
    hh6[1,1]<-sigma[1]/(m_sam[1]+n0[1])+sigma[6]/n0[3]+sigma[7]/n0[4]+sigma[9]/n0[6]
    hh6[1,2]<-sigma[6]/n0[3]
    hh6[1,3]<-sigma[1]/(m_sam[1]+n0[1])+sigma[6]/n0[3]
    hh6[1,4]<--sigma[1]/(m_sam[1]+n0[1])-sigma[7]/n0[4]+sigma[9]/n0[6]
    hh6[2,2]<-sigma[5]/n0[2]+sigma[6]/n0[3]+sigma[11]/n0[8]+sigma[3]/(m_sam[3]+n0[9])
    hh6[2,3]<-sigma[3]/n0[3]+sigma[3]/(m_sam[3]+n0[9])
    hh6[2,4]<-sigma[5]/n0[2]-sigma[11]/n0[8]-sigma[3]/(m_sam[3]+n0[9])
    hh6[3,3]<-sigma[1]/(m_sam[1]+n0[1])+sigma[3]/n0[3]+sigma[7]/n0[7]+sigma[9]/(m_sam[3]+n0[9])
    hh6[3,4]<--sigma[1]/(m_sam[1]+n0[1])-sigma[9]/(m_sam[3]+n0[9])
    hh6[4,4]<-sigma[1]/(m_sam[1]+n0[1])+sigma[5]/n0[2]+sigma[7]/n0[4]+4.0*sigma[2]/(m_sam[2]+n0[5])+sigma[9]/n0[6]+sigma[11]/n0[8]+sigma[3]/(m_sam[3]+n0[9])
    for(i in 2:4)
    {
      for(j in 1:(i-1))
      {
        hh6[i,j]<-hh6[j,i]
      }
    }
    ###########################################################
    b_line6[1]<-(sumx[1]+sumwx[1])/(m_sam[1]+n0[1])-sumwx[3]/n0[3]-sumwx[4]/n0[4]+sumwx[6]/n0[6]
    b_line6[2]<-sumwx[2]/n0[2]-sumwx[3]/n0[3]-sumwx[8]/n0[8]+(sumx[3]+sumwx[9])/(m_sam[3]+n0[9])
    b_line6[3]<-(sumx[1]+sumwx[1])/(m_sam[1]+n0[1])-sumwx[3]/n0[3]-sumwx[7]/n0[7]+(sumx[3]+sumwx[9])/(m_sam[3]+n0[9])
    b_line6[4]<-sumwx[2]/n0[2]+sumwx[4]/n0[4]+sumwx[6]/n0[6]+sumwx[8]/n0[8]-(sumx[1]+sumwx[1])/(m_sam[1]+n0[1])-2*(sumx[2]+sumwx[5])/(m_sam[2]+n0[5])-(sumx[3]+sumwx[9])/(m_sam[3]+n0[9]) 
    B6<-solve(hh6,b_line6)
    #################obtain mean##############################
    m[1]<-(sumx[1]+sumwx[1]-sigma[1]*(B6[1]+B6[3]-B6[4]))/(m_sam[1]+n0[1])
    m[2]<-(sumx[2]+sumwx[5]+sigma[2]*2*B6[4])/(m_sam[2]+n0[5])
    m[3]<-(sumx[3]+sumwx[9]-sigma[3]*(B6[2]+B6[3]-B6[4]))/(m_sam[3]+n0[9])
    m[5]<-(sumwx[2]-sigma[5]*(B6[2]+B6[4]))/n0[2]
    m[6]<-(sumwx[3]+sigma[6]*(B6[1]+B6[2]+B6[3]))/n0[3]
    m[7]<-(sumwx[4]+sigma[7]*(B6[1]-B6[4]))/n0[4]
    m[9]<-(sumwx[6]-sigma[9]*(B6[1]+B6[4]))/n0[6]
    m[10]<-(sumwx[7]+sigma[10]*B6[3])/n0[7]
    m[11]<-(sumwx[8]+sigma[11]*(B6[2]-B6[4]))/n0[8]
    m[4]<-m[1]; m[8]<-m[2];m[12]<-m[3]
    ################obtain variance##########################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d3) {swx[i] <- WW[i,]%*%(dataF2-m[i+3])^2 }
    s0[1]<-sum(swx);s0[1]<-s0[1]+ss1+ss2+ss3
    sigma[1]<-s0[1]/(m_sam[1]+m_sam[2]+m_sam[3]+m_sam[4])
    sigma[c(2:12)]<-sigma[1]
    ########criteria for iterations to stop####################
    pi<-m[c(4:12)];gh<-sigma[c(4:12)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mix_pi))) 
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*7
  #############first order genetic parameter process##############
  hh61<-matrix(c(1,1,1,1,1,1,1,1,1,1,0,-1,1,1,0,0,-1,-1,1,0,-1,0,-1,1,-1,1,0,0,1,0,0,0,1,1,0,0,0,1,0,1,0,0,0,0,1),9,5)
  B61<-solve(crossprod(hh61,hh61))%*%crossprod(hh61,m[c(1,2,3,5,6,7,9,10,11)])
  #################second oder genetic parameter process#############
  jj<-sigmaF2-sigma[1]
  if(jj<0){jj<-0}
  ll<-jj/sigmaF2
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
  
  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2); 
  F2w1<-1/(12*m_sam[4])
  F2bmw <- matrix(0,m_sam[4],1); F2bmwsl <- matrix(0,m_sam[4],d3)
  for(i in 1:d3){
    F2gg <- (dataF2 - pi[i])/sqrt(gh[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < m_sam[4]){F2P2 <- F2P2+runif(m_sam[4])/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*m_sam[4]) + sum((F2P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F2u <- as.matrix(c(12*m_sam[4]*((F2dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F2dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F2dd[3]/m_sam[4]-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),(1-pkolm(F2D,m_sam[4]))))
  
  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("2MG-AD",round(abc,4),round(AIC,4),round(t(m),4),round(sigma[1],4),round(sigma[4],4),
                       round(t(mix_pi),4),round(B61[1],4)," "," "," ",round(B61[2],4),round(B61[3],4),round(B61[4],4),round(B61[5],4)," "," "," "," "," "," ",round(jj,4),round(ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}

##############2MG-A(B-3)#################################### 
G4F2ModelFun[[7]] <- function(K1,logL,df11,df21,df31,df41){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF2 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2)))
  m_esp<-0.0001
  sigmaF2<-var(dataF2);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a2<-sqrt(sigmaF2/m_sam[4]*(m_sam[4]-1))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2)))
  ###############procedure start##############################   
  d3<-9
  mi<-as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  sigma<-matrix((sigma0),12,1)
  m[c(1:3)]<-m[c(1:3)]
  if(m[1]<m[3]){a2<--a2}
  m[5]<-m[4]+2.4*a2;m[6]<-m[4]+1.6*a2;m[7]<-m[4]+0.8*a2
  m[8]<-m[4];m[9]<-m[4]-0.8*a2;m[10]<-m[4]-1.6*a2
  m[11]<-m[4]-2.4*a2;m[12]<-m[4]-3.2*a2;m[4]<-m[4]+3.2*a2
  pi<-m[c(4:12)];gh<-sigma[c(4:12)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mi)))
  ########iteration process############
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,9,m_sam[4]);swx <- matrix(0,9,1)
  hh7<-matrix(0,6,6);b_line7<-matrix(0,6,1);s0<-matrix(0,1,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:12)];ssigma<-sigma[c(4:12)]
    for(i in 1:d3) { WW[i,] <- mi[i]*dnorm(dataF2,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF2,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF2
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    ########################################
    hh7[1,1]<-sigma[1]/(m_sam[1]+n0[1])+4*sigma[2]/n0[2]+sigma[3]/n0[3]
    hh7[1,2]<-sigma[1]/(m_sam[1]+n0[1])
    hh7[1,3]<--2*sigma[2]/n0[2]-sigma[3]/n0[3]
    hh7[1,4]<-sigma[3]/n0[3]
    hh7[1,5]<-0;hh7[1,6]<-0
    hh7[2,2]<-sigma[1]/(m_sam[1]+n0[1])+4*sigma[5]/(m_sam[2]+n0[5])+sigma[9]/(m_sam[3]+n0[9])
    hh7[2,3]<-sigma[9]/(m_sam[3]+n0[9])
    hh7[2,4]<-4*sigma[5]/(m_sam[2]+n0[5])
    hh7[2,5]<-4*sigma[5]/(m_sam[2]+n0[5])
    hh7[2,6]<-2*sigma[9]/(m_sam[3]+n0[9])
    hh7[3,3]<-sigma[2]/n0[2]+sigma[3]/n0[3]+sigma[8]/n0[8]+sigma[9]/(m_sam[3]+n0[9])
    hh7[3,4]<--sigma[3]/n0[3]
    hh7[3,5]<-0
    hh7[3,6]<-2*sigma[8]/n0[8]+2*sigma[9]/(m_sam[3]+n0[9])
    hh7[4,4]<-sigma[3]/n0[3]+4*sigma[5]/(m_sam[2]+n0[5])+sigma[7]/n0[7]
    hh7[4,5]<-4*sigma[5]/(m_sam[2]+n0[5])
    hh7[4,6]<-0
    hh7[5,5]<-sigma[4]/n0[4]+4*sigma[5]/(m_sam[2]+n0[5])+sigma[6]/n0[6]
    hh7[5,6]<-sigma[4]/n0[4]-sigma[6]/n0[6]
    hh7[6,6]<-sigma[4]/n0[4]+sigma[6]/n0[6]+4*sigma[8]/n0[8]+4*sigma[9]/(m_sam[3]+n0[9])
    for(i in 2:6)
    {
      for(j in 1:(i-1))
      {
        hh7[i,j]<-hh7[j,i]
      }
    }
    ############################################################## 
    b_line7[1]<-(sumx[1]+sumwx[1])/(m_sam[1]+n0[1])-2*sumwx[2]/n0[2]+sumwx[3]/n0[3]
    b_line7[2]<-(sumx[1]+sumwx[1])/(m_sam[1]+n0[1])-2*(sumx[2]+sumwx[5])/(m_sam[2]+n0[5])+(sumx[3]+sumwx[9])/(m_sam[3]+n0[9])
    b_line7[3]<-sumwx[2]/n0[2]-sumwx[3]/n0[3]-sumwx[8]/n0[8]+(sumx[3]+sumwx[9])/(m_sam[3]+n0[9])
    b_line7[4]<-sumwx[3]/n0[3]-2*(sumx[2]+sumwx[5])/(m_sam[2]+n0[5])+sumwx[7]/n0[7]
    b_line7[5]<-sumwx[4]/n0[4]-2*(sumx[2]+sumwx[5])/(m_sam[2]+n0[5])+sumwx[6]/n0[6]
    b_line7[6]<-sumwx[4]/n0[4]-sumwx[6]/n0[6]-2*sumwx[8]/n0[8]+2*(sumx[3]+sumwx[9])/(m_sam[3]+n0[9]) 
    B7<-solve(hh7,b_line7)
    ################obtain mean##############################
    m[1]<-(sumx[1]+sumwx[1]-sigma[1]*(B7[1]+B7[2]))/(m_sam[1]+n0[1])
    m[2]<-(sumx[2]+sumwx[5]+sigma[2]*(2*B7[2]+2*B7[4]+2*B7[5]))/(m_sam[2]+n0[5])
    m[3]<-(sumx[3]+sumwx[9]-sigma[3]*(B7[2]+B7[3]+2*B7[6]))/(m_sam[3]+n0[9])
    m[5]<-(sumwx[2]+sigma[2]*(2*B7[1]-B7[3]))/n0[2]
    m[6]<-(sumwx[3]+sigma[3]*(-B7[1]+B7[3]-B7[4]))/n0[3]
    m[7]<-(sumwx[4]+sigma[4]*(-B7[5]-B7[6]))/n0[4]
    m[9]<-(sumwx[6]+sigma[6]*(-B7[5]+B7[6]))/n0[6]
    m[10]<-(sumwx[7]-sigma[7]*B7[4])/n0[7]
    m[11]<-(sumwx[8]+sigma[8]*(B7[3]+2*B7[6]))/n0[8]
    m[4]<-m[1]; m[8]<-m[2];m[12]<-m[3]
    ##############obtain variance#############################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d3) {swx[i] <- WW[i,]%*%(dataF2-m[i+3])^2 }
    s0[1]<-sum(swx);s0[1]<-s0[1]+ss1+ss2+ss3
    sigma[1]<-s0[1]/(m_sam[1]+m_sam[2]+m_sam[3]+m_sam[4])
    sigma[c(2:12)]<-sigma[1]
    ########criteria for iterations to stop####################
    pi<-m[c(4:12)];gh<-sigma[c(4:12)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mix_pi))) 
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*5
  #############first order genetic parameter process##############
  hh71<-matrix(c(1,1,1,1,1,1,1,1,1,1,0,-1,1,1,0,0,-1,-1,1,0,-1,0,-1,1,-1,1,0),9,3)
  B71<-solve(crossprod(hh71,hh71))%*%crossprod(hh71,m[c(1,2,3,5,6,7,9,10,11)])
  #################second oder genetic parameter process#############
  jj<-sigmaF2-sigma[1]
  if(jj<0){jj<-0}
  ll<-jj/sigmaF2
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
  
  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2); 
  F2w1<-1/(12*m_sam[4])
  F2bmw <- matrix(0,m_sam[4],1); F2bmwsl <- matrix(0,m_sam[4],d3)
  for(i in 1:d3){
    F2gg <- (dataF2 - pi[i])/sqrt(gh[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < m_sam[4]){F2P2 <- F2P2+runif(m_sam[4])/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*m_sam[4]) + sum((F2P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F2u <- as.matrix(c(12*m_sam[4]*((F2dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F2dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F2dd[3]/m_sam[4]-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),(1-pkolm(F2D,m_sam[4]))))
  
  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("2MG-A",round(abc,4),round(AIC,4),round(t(m),4),round(sigma[1],4),round(sigma[4],4),
                       round(t(mix_pi),4),round(B71[1],4)," "," "," ",round(B71[2],4),round(B71[3],4)," "," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}

###################2MG-EA(B-4)########################### 
G4F2ModelFun[[8]] <- function(K1,logL,df11,df21,df31,df41){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF2 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2)))
  m_esp<-0.0001
  sigmaF2<-var(dataF2);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a2<-sqrt(sigmaF2/m_sam[4]*(m_sam[4]-1))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2)))
  ###############procedure start##############################   
  d4<-5
  mi<-matrix(c(0.0625,0.25,0.375,0.25,0.0625))
  sigma<-matrix((sigma0),8,1)
  m[c(1:3)]<-m[c(1:3)]
  if (m[1]<m[3]) {a2<--a2}
  m[5]<-m[4]+1.5*a2;m[6]<-m[4];m[7]<-m[4]-1.5*a2
  m[8]<-m[4]-3*a2;m[4]<-m[4]+3*a2
  pi<-m[c(4:8)];gh<-sigma[c(4:8)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mi)))
  ########iteration process############
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,5,m_sam[4]);swx <- matrix(0,5,1)
  hh8<-matrix(0,3,3);b_line8<-matrix(0,3,1);s0<-matrix(0,1,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:8)];ssigma<-sigma[c(4:8)]
    for(i in 1:d4) { WW[i,] <- mi[i]*dnorm(dataF2,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF2,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF2
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    #################################################
    hh8[1,1]<-sigma[1]/(m_sam[1]+n0[1])+4*sigma[5]/n0[2]+sigma[2]/(m_sam[2]+n0[3])
    hh8[1,2]<--2*sigma[5]/n0[2]-2*sigma[2]/(m_sam[2]+n0[3])
    hh8[1,3]<-sigma[2]/(m_sam[2]+n0[3])
    hh8[2,2]<-sigma[5]/n0[2]+4*sigma[2]/(m_sam[2]+n0[3])+sigma[7]/n0[4]
    hh8[2,3]<--2*sigma[2]/(m_sam[2]+n0[3])-2*sigma[7]/n0[4]
    hh8[3,3]<-sigma[2]/(m_sam[2]+n0[3])+4*sigma[7]/n0[4]+sigma[3]/(m_sam[3]+n0[5])
    for(i in 2:3)
    {
      for(j in 1:(i-1))
      {
        hh8[i,j]<-hh8[j,i]
      }
    }
    ####################################################
    b_line8[1]<-(sumx[1]+sumwx[1])/(m_sam[1]+n0[1])-2*sumwx[2]/n0[2]+(sumx[2]+sumwx[3])/(m_sam[2]+n0[3])
    b_line8[2]<-sumwx[2]/n0[2]-2*(sumx[2]+sumwx[3])/(m_sam[2]+n0[3])+sumwx[4]/n0[4]
    b_line8[3]<-(sumx[2]+sumwx[3])/(m_sam[2]+n0[3])-2*sumwx[4]/n0[4]+(sumx[3]+sumwx[5])/(m_sam[3]+n0[5])
    B8<-solve(hh8,b_line8)
    ###########obtain varianc######################
    m[1]<-(sumx[1]+sumwx[1]-sigma[1]*B8[1])/(m_sam[1]+n0[1])
    m[2]<-(sumx[2]+sumwx[3]+sigma[2]*(-B8[1]+2*B8[2]-B8[3]))/(m_sam[2]+n0[3])
    m[3]<-(sumx[3]+sumwx[5]-sigma[3]*B8[3])/(m_sam[3]+n0[5])
    m[5]<-(sumwx[2]+sigma[5]*(2*B8[1]-B8[2]))/n0[2]
    m[7]<-(sumwx[4]+sigma[7]*(-B8[2]+2*B8[3]))/n0[4]
    m[4]<-m[1];m[6]<-m[2];m[8]<-m[3]
    #############obtain variance################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d4) {swx[i] <- WW[i,]%*%(dataF2-m[i+3])^2 }
    s0[1]<-sum(swx);s0[1]<-s0[1]+ss1+ss2+ss3
    sigma[1]<-s0[1]/(m_sam[1]+m_sam[2]+m_sam[3]+m_sam[4])
    sigma[c(2:8)]<-sigma[1]
    ########criteria for iterations to stop####################
    pi<-m[c(4:8)];gh<-sigma[c(4:8)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mix_pi))) 
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*4
  #############first order genetic parameter process##############
  hh81<-matrix(c(1,1,1,1,1,2,0,-2,1,-1),5,2)
  B81<-solve(crossprod(hh81,hh81))%*%crossprod(hh81,m[c(1,2,3,5,7)])
  #################second oder genetic parameter process#############
  jj<-sigmaF2-sigma[1]
  if(jj<0){jj<-0}
  ll<-jj/sigmaF2   
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
  
  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2); 
  F2w1<-1/(12*m_sam[4])
  F2bmw <- matrix(0,m_sam[4],1); F2bmwsl <- matrix(0,m_sam[4],d4)
  for(i in 1:d4){
    F2gg <- (dataF2 - pi[i])/sqrt(gh[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < m_sam[4]){F2P2 <- F2P2+runif(m_sam[4])/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*m_sam[4]) + sum((F2P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F2u <- as.matrix(c(12*m_sam[4]*((F2dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F2dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F2dd[3]/m_sam[4]-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),(1-pkolm(F2D,m_sam[4]))))
  
  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("2MG-EA",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," ",round(sigma[1],4),round(sigma[4],4),
                       round(t(mix_pi),4)," "," "," "," ",round(B81[1],4)," "," "," ",round(B81[2],4)," "," "," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}

###############2MG-CD(B-5)######################
G4F2ModelFun[[9]] <- function(K1,logL,df11,df21,df31,df41){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF2 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2)))
  m_esp<-0.0001
  sigmaF2<-var(dataF2);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a2<-sqrt(sigmaF2/m_sam[4]*(m_sam[4]-1))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2)))
  ###############procedure start############################## 
  d5<-4
  mi<-matrix(c(0.5625,0.1875,0.1875,0.0625))
  sigma<-matrix((sigma0),7,1)
  m[c(1:3)]<-m[c(1:3)]
  if (m[1]<m[3]) {a2<--a2}
  m[5]<-m[4]+0.5*a2;m[6]<-m[4]-0.5*a2;m[7]<-m[4]-3*a2;m[4]<-m[4]+3*a2
  pi<-m[c(4:7)];gh<-sigma[c(4:7)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mi)))
  ########iteration process############
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,4,m_sam[4]);swx <- matrix(0,4,1);s0<-matrix(0,1,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:7)];ssigma<-sigma[c(4:7)]
    for(i in 1:d5) { WW[i,] <- mi[i]*dnorm(dataF2,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF2,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF2
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    aa1<-(sumx[1]+sumx[2]+sumwx[1])/(m_sam[1]+m_sam[2]+n0[1])-sumwx[2]/n0[2]-sumwx[3]/n0[3]+(sumx[3]+sumwx[4])/(m_sam[3]+n0[4])
    aa2<-sigma[1]/(m_sam[1]+m_sam[2]+n0[1])+sigma[5]/n0[2]+sigma[6]/n0[3]+sigma[3]/(m_sam[3]+n0[4])
    aa3<-aa1/aa2
    ###############obtain mean################
    m[1]<-(sumx[1]+sumx[2]+sumwx[1]-sigma[1]*aa3)/(m_sam[1]+m_sam[2]+n0[1])
    m[3]<-(sumx[3]+sumwx[4]-sigma[3]*aa3)/(m_sam[3]+n0[4])
    m[5]<-(sumwx[2]+sigma[5]*aa3)/n0[2]
    m[6]<-(sumwx[3]+sigma[6]*aa3)/n0[3]
    m[2]<-m[1]; m[4]<-m[1];m[7]<-m[3]
    ###############obtain variance##############################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d5) {swx[i] <- WW[i,]%*%(dataF2-m[i+3])^2 }
    s0[1]<-sum(swx);s0[1]<-s0[1]+ss1+ss2+ss3
    sigma[1]<-s0[1]/(m_sam[1]+m_sam[2]+m_sam[3]+m_sam[4])
    sigma[c(2:7)]<-sigma[1]
    ########criteria for iterations to stop####################
    pi<-m[c(4:7)];gh<-sigma[c(4:7)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mix_pi))) 
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*5
  #################first oder genetic parameter process#############
  hh9<-matrix(c(1,1,1,1,1,-1,1,-1,1,-1,-1,1),4,3)
  B9<-solve(crossprod(hh9,hh9))%*%crossprod(hh9,m[c(1,3,5,6)])
  #################second oder genetic parameter process#############
  jj<-sigmaF2-sigma[1]
  if(jj<0){jj<-0}
  ll<-jj/sigmaF2   
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
  
  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2); 
  F2w1<-1/(12*m_sam[4])
  F2bmw <- matrix(0,m_sam[4],1); F2bmwsl <- matrix(0,m_sam[4],d5)
  for(i in 1:d5){
    F2gg <- (dataF2 - pi[i])/sqrt(gh[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < m_sam[4]){F2P2 <- F2P2+runif(m_sam[4])/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*m_sam[4]) + sum((F2P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F2u <- as.matrix(c(12*m_sam[4]*((F2dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F2dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F2dd[3]/m_sam[4]-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),(1-pkolm(F2D,m_sam[4]))))
  
  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("2MG-CD",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," ",round(sigma[1],4),round(sigma[4],4),
                       round(t(mix_pi),4)," "," "," "," "," ",round(B9[1],4)," "," "," ",round(B9[2],4),round(B9[3],4)," "," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}

########################2MG-EAD(B-6)############################################
G4F2ModelFun[[10]] <- function(K1,logL,df11,df21,df31,df41){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF2 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2)))
  m_esp<-0.0001
  sigmaF2<-var(dataF2);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a2<-sqrt(sigmaF2/m_sam[4]*(m_sam[4]-1))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2)))
  ###############procedure start############################## 
  d1<-3
  mi<-as.matrix(c(0.5625,0.375,0.0625))
  sigma<-matrix((sigma0),6,1)
  m[c(1:3)]<-m[c(1:3)]
  if (m[1]<m[3]) {a2<--a2}
  m[5]<-m[4];m[6]<-m[4]-3*a2;m[4]<-m[4]+3*a2
  pi<-m[c(4:6)];gh<-sigma[c(4:6)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mi)))
  ########iteration process############
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,3,m_sam[4]);swx <- matrix(0,3,1);s0<-matrix(0,1,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:6)];ssigma<-sigma[c(4:6)]
    for(i in 1:d1) { WW[i,] <- mi[i]*dnorm(dataF2,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF2,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF2
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    aa1<-(sumx[1]+sumx[2]+sumwx[1])/(m_sam[1]+m_sam[2]+n0[1])-2*sumwx[2]/n0[2]+(sumx[3]+sumwx[3])/(m_sam[3]+n0[3])
    aa2<-sigma[1]/(m_sam[1]+m_sam[2]+n0[1])+4*sigma[5]/n0[2]+sigma[3]/(m_sam[3]+n0[3])
    aa3<-aa1/aa2
    #########obtain mean##########################
    m[1]<-(sumx[1]+sumx[2]+sumwx[1]-sigma[1]*aa3)/(m_sam[1]+m_sam[2]+n0[1])
    m[3]<-(sumx[3]+sumwx[3]-sigma[3]*aa3)/(m_sam[3]+n0[3])
    m[5]<-(sumwx[2]+sigma[2]*2*aa3)/n0[2]
    m[4]<-m[2]<-m[1]; m[6]<-m[3]
    ############obtain variance######################################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d1) {swx[i] <- WW[i,]%*%(dataF2-m[i+3])^2 }
    s0[1]<-sum(swx);s0[1]<-s0[1]+ss1+ss2+ss3
    sigma[1]<-s0[1]/(m_sam[1]+m_sam[2]+m_sam[3]+m_sam[4])
    sigma[c(2:6)]<-sigma[1]
    ########criteria for iterations to stop####################
    pi<-m[c(4:6)];gh<-sigma[c(4:6)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mix_pi))) 
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*4
  #############first order genetic parameter process##############
  hh10<-matrix(c(1,1,1,2,-2,0),3,2)
  B10<-solve(crossprod(hh10,hh10))%*%crossprod(hh10,m[c(1,3,5)])
  #################second oder genetic parameter process#############
  jj<-sigmaF2-sigma[1]
  if(jj<0){jj<-0}
  ll<-jj/sigmaF2
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
  
  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2); 
  F2w1<-1/(12*m_sam[4])
  F2bmw <- matrix(0,m_sam[4],1); F2bmwsl <- matrix(0,m_sam[4],d1)
  for(i in 1:d1){
    F2gg <- (dataF2 - pi[i])/sqrt(gh[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < m_sam[4]){F2P2 <- F2P2+runif(m_sam[4])/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*m_sam[4]) + sum((F2P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F2u <- as.matrix(c(12*m_sam[4]*((F2dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F2dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F2dd[3]/m_sam[4]-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),(1-pkolm(F2D,m_sam[4]))))
  
  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("2MG-EAD",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," "," ",round(sigma[1],4),round(sigma[4],4),
                       round(t(mix_pi),4)," "," "," "," "," "," ",round(B10[1],4)," "," "," ",round(B10[2],4)," "," "," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}

###############PG-ADI(C-0)################################
G4F2ModelFun[[11]] <- function(K1,logL,df11,df21,df31,df41){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF2 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2)))
  m_esp<-0.0001
  sigmaF2<-var(dataF2);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2)))
  ###############procedure start############################## 
  d6<-1;mi<-1
  sigma<-matrix(c(sigma0,sigma0,sigma0,2*2*sigma0))
  m[c(1:4)]<-m
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dnorm(dataF2,m[4],sqrt(sigma[4]))))
  ###########################iteration process#########################################
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,1,m_sam[4]);swx <- matrix(0,1,1);s0<-matrix(0,1,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    WW[1,] <- mi*dnorm(dataF2,m[4],sqrt(sigma[4]))/dmixnorm(dataF2,m[4],sqrt(sigma[4]),mi)
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF2
    ############obtain mean#####################
    mix_pi[1]<-1;m[c(1:4)]<-m[c(1:4)]
    ############obtain variance######################################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    swx[1] <- WW[1,]%*%(dataF2-m[4])^2 
    s0[1]<-swx[1];sigma[4]<-s0[1]/m_sam[4]
    
    sigma00<-sigma[4]-sigma[1]
    if (sigma00<0) {sigma00<-0;sigma[4]<-sigma[1]}
    aa0<-sigma[1];aa2<-ss1+ss2+ss3;aa3<-m_sam[1]+m_sam[2]+m_sam[3]
    n_iter<-0
    aa4<-1000
    while(aa4>0.0001){
      n_iter<-n_iter+1
      aa1<-sigma[1]/(sigma[1]+sigma00)
      sigma[1]<-(aa2+aa1^2*s0[1])/(aa3+aa1*m_sam[4])
      aa4<-abs(sigma[1]-aa0)
      aa0<-sigma[1]
      if(n_iter>20)break
    }
    sigma[2]<-sigma[3]<-sigma[1];sigma[4]<-sigma[1]+sigma00
    ###################criteria for iterations to stop#######################################
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+ sum(log(dnorm(dataF2,m[4],sqrt(sigma[4]))))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*6
  #######first order genetic parameter process############  
  ma1<-m[1];ma2<-m[2];ma3<-m[3];ma4<-m[4]
  #######second order genetic parameter process############  
  gg <- sigmaF2-sigma[1]
  if(gg<0 || gg>sigmaF2) {gg<-0}
  rr <- gg/sigmaF2
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
  
  ####################################hypothesis testing for F2#########################################
  dataF2<-sort(dataF2)
  F2w1<-1/(12*m_sam[4])
  F2bmw <- matrix(0,m_sam[4],1)
  F2gg <- (dataF2 - m[4])/sqrt(as.vector(sigma[4]))
  F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
  F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2bmw)))[1]
  if(nn < m_sam[4]){F2bmw <- F2bmw+runif(m_sam[4])/1e4}
  ##########################################################
  F2dd<-c((sum(F2bmw)),(sum(F2bmw^2)),sum((F2bmw-0.5)^2))
  F2w<-F2w1+sum((F2bmw - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F2u<- as.matrix(c(12*m_sam[4]*((F2dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F2dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F2dd[3]/m_sam[4]-1/12)^2)))
  F2D<-as.numeric(ks.test(F2bmw,"punif")[[1]][1])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2w),(1-pkolm(F2D,m_sam[4]))))
  
  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("PG-ADI",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," "," "," "," ",round(sigma[1],4),round(sigma[4],4),
                       round(t(mix_pi),4)," "," "," "," "," "," "," "," ",round(ma1,4),round(ma2,4),round(ma3,4),round(ma4,4)," "," "," "," "," "," "," "," "," "," "," "," ",round(gg,4),round(rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2w,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output)
  return(OUTPUT)
}

###################PG-AD(C-1)###############################  
G4F2ModelFun[[12]] <- function(K1,logL,df11,df21,df31,df41){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF2 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2)))
  m_esp<-0.0001
  sigmaF2<-var(dataF2);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2)))
  ###############procedure start##############################   
  d6<-1;mi<-1
  sigma<-matrix(c(sigma0,sigma0,sigma0,2*2*sigma0))
  m[c(1:4)]<-m
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dnorm(dataF2,m[4],sqrt(sigma[4]))))
  ###########################iteration process#########################################
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,1,m_sam[4]);swx <- matrix(0,1,1); s0<-matrix(0,1,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    WW[1,] <- mi*dnorm(dataF2,m[4],sqrt(sigma[4]))/dmixnorm(dataF2,m[4],sqrt(sigma[4]),mi)
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF2
    ############obtain mean#####################
    mix_pi[1]<-1
    aa1<-sumx[1]/m_sam[1]+2*sumx[2]/m_sam[2]+sumx[3]/m_sam[3]-4*m[4]
    aa2<-sigma[1]/m_sam[1]+4*sigma[2]/m_sam[2]+sigma[3]/m_sam[3]+16*sigma[4]/m_sam[4]
    aa3<-aa1/aa2
    m[1]<-(sumx[1]-sigma[1]*aa3)/m_sam[1];m[2]<-(sumx[2]-sigma[2]*aa3*2)/m_sam[2]
    m[3]<-(sumx[3]-sigma[3]*aa3)/m_sam[3];m[4]<-(m[4]*m_sam[4]+sigma[4]*aa3*4)/m_sam[4]
    ############obtain variance######################################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    swx[1] <- WW[1,]%*%(dataF2-m[4])^2 
    s0[1]<-swx[1];sigma[4]<-s0[1]/m_sam[4]
    
    sigma00<-sigma[4]-sigma[1]
    if (sigma00<0) {sigma00<-0;sigma[4]<-sigma[1]}
    aa0<-sigma[1]
    aa2<-ss1+ss2+ss3;aa3<-m_sam[1]+m_sam[2]+m_sam[3]
    n_iter<-0
    aa4<-1000
    while(aa4>0.0001){
      n_iter<-n_iter+1
      aa1<-sigma[1]/(sigma[1]+sigma00)
      sigma[1]<-(aa2+aa1^2*s0[1])/(aa3+aa1*m_sam[4])
      aa4<-abs(sigma[1]-aa0)
      aa0<-sigma[1]
      if(n_iter>20)break
    }
    sigma[2]<-sigma[3]<-sigma[1];sigma[4]<-sigma[1]+sigma00
    ###################criteria for iterations to stop#######################################
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+ sum(log(dnorm(dataF2,m[4],sqrt(sigma[4]))))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*5
  #############first order genetic parameter process##############
  hh11<-matrix(c(1,1,1,1,1,0,-1,0,0,1,0,0.5),4,3)
  B11<-solve(crossprod(hh11,hh11))%*%crossprod(hh11,m)
  #################second oder genetic parameter process#############
  gg<-sigmaF2-sigma[1]
  if(gg<0 || gg>sigmaF2) {gg<-0}
  rr<-gg/sigmaF2
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
  
  ####################################hypothesis testing for F2#########################################
  dataF2<-sort(dataF2)
  F2w1<-1/(12*m_sam[4])
  F2bmw <- matrix(0,m_sam[4],1)
  F2gg <- (dataF2 - m[4])/sqrt(as.vector(sigma[4]))
  F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
  F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2bmw)))[1]
  if(nn < m_sam[4]){F2bmw <- F2bmw+runif(m_sam[4])/1e4}
  ##########################################################
  F2dd<-c((sum(F2bmw)),(sum(F2bmw^2)),sum((F2bmw-0.5)^2))
  F2w<-F2w1+sum((F2bmw - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F2u<- as.matrix(c(12*m_sam[4]*((F2dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F2dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F2dd[3]/m_sam[4]-1/12)^2)))
  F2D<-as.numeric(ks.test(F2bmw,"punif")[[1]][1])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2w),(1-pkolm(F2D,m_sam[4]))))
  
  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("PG-AD",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," "," "," "," ",round(sigma[1],4),round(sigma[4],4),
                       round(t(mix_pi),4)," "," "," "," "," "," "," "," ",round(B11[1],4)," "," "," "," "," "," "," "," "," "," "," ",round(B11[2],4),round(B11[3],4)," "," ",round(gg,4),round(rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2w,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output)
  return(OUTPUT)
}

###################MX1-AD-ADI(D-0)####################################  
G4F2ModelFun[[13]] <- function(K1,logL,df11,df21,df31,df41){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF2 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2)))
  m_esp<-0.0001
  sigmaF2<-var(dataF2);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a1<-sqrt(sigmaF2/m_sam[4])
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2)))
  ###############procedure start##############################     
  mi<-as.matrix(c(0.25,0.5,0.25))
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaF2/2,sigmaF2/2,sigmaF2/2))
  m[c(1:4)]<-m
  if(m[1]<m[3]) {a1<--a1}
  m[5]<-m[4];m[4]<-m[5]+2*a1;m[6]<-m[5]-2*a1
  pi<-m[c(4:6)];gh<-sigma[c(4:6)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mi)))
  ########iteration process############
  d1<-3
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,3,m_sam[4]); swx <- matrix(0,3,1);s0<-matrix(0,1,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4,5,6)];ssigma<-sigma[c(4,5,6)]
    for(i in 1:d1) { WW[i,] <- mi[i]*dnorm(dataF2,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF2,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF2
    ###########obtain mean#########################
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    m[c(1:3)]<-sumx[c(1:3)]/m_sam[c(1:3)]
    m[c(4:6)]<-sumwx[1:3]/n0[1:3]
    ############obtain variance##############################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d1) {swx[i] <- WW[i,]%*%(dataF2-m[i+3])^2 }
    s0[1]<-swx[1]+swx[2]+swx[3];sigma[4]<-s0[1]/m_sam[4]
    
    sigma00<-sigma[4]-sigma[1]
    if (sigma00<0) {sigma00<-0;sigma[4]<-sigma[1]}
    sigma[5]<-sigma[6]<-sigma[4]
    aa0<-sigma[1]
    aa2<-ss1+ss2+ss3;aa3<-m_sam[1]+m_sam[2]+m_sam[3];
    n_iter<-0;
    aa4<-1000
    while(aa4>0.0001){
      n_iter<-n_iter+1
      aa1<-sigma[1]/(sigma[1]+sigma00)
      sigma[1]<-(aa2+aa1^2*s0[1])/(aa3+aa1*m_sam[4])
      aa4<-abs(sigma[1]-aa0)
      aa0<-sigma[1]
      if(n_iter>20)break
    }
    sigma[2]<-sigma[3]<-sigma[1]
    sigma[4]<-sigma[1]+sigma00;sigma[5]<-sigma[4];sigma[6]<-sigma[4];
    ########criteria for iterations to stop####################
    pi<-m[c(4:6)];gh<-sigma[c(4:6)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mix_pi))) 
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*8
  #########first order genetic parameter process##########
  hh12 <- matrix(c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,1,1,1,0,-1,1,0,-1,0,1,0,0,1,0),6,6)
  B12 <- solve(hh12,m)
  #################second oder genetic parameter process#############
  jj<-sigmaF2 - sigma[4]
  gg<-sigma[4]-sigma[1]
  if(jj<0){jj<-0}
  if(gg<0){gg<-0}
  ll<-jj/sigmaF2
  rr<-gg/sigmaF2
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
  
  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2); 
  F2w1<-1/(12*m_sam[4])
  F2bmw <- matrix(0,m_sam[4],1); F2bmwsl <- matrix(0,m_sam[4],d1)
  for(i in 1:d1){
    F2gg <- (dataF2 - pi[i])/sqrt(gh[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < m_sam[4]){F2P2 <- F2P2+runif(m_sam[4])/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*m_sam[4]) + sum((F2P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F2u <- as.matrix(c(12*m_sam[4]*((F2dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F2dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F2dd[3]/m_sam[4]-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),(1-pkolm(F2D,m_sam[4]))))
  
  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("MX1-AD-ADI",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," "," ",round(sigma[1],4),round(sigma[4],4),
                       round(t(mix_pi),4)," "," "," "," "," "," ",round(B12[1],4),round(B12[2],4),round(B12[3],4),round(B12[4],4),round(B12[5],4)," ",round(B12[6],4)," "," "," "," "," "," "," ",round(jj,4),round(ll*100,4),round(gg,4),round(rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}

############MX1-AD-AD(D-1)#############################################   
G4F2ModelFun[[14]] <- function(K1,logL,df11,df21,df31,df41){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF2 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2)))
  m_esp<-0.0001
  sigmaF2<-var(dataF2);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a1<-sqrt(sigmaF2/m_sam[4])
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2)))
  ###############procedure start##############################        
  mi<-as.matrix(c(0.25,0.5,0.25))
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaF2/2,sigmaF2/2,sigmaF2/2))
  m[c(1:4)]<-m
  if(m[1]<m[3]) {a1<--a1}
  m[5]<-m[4];m[4]<-m[5]+2*a1;m[6]<-m[5]-2*a1
  pi<-m[c(4:6)];gh<-sigma[c(4:6)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mi)))
  ########iteration process############
  d1<-3
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,3,m_sam[4]); swx <- matrix(0,3,1); s0<-matrix(0,1,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:6)];ssigma<-sigma[c(4:6)]
    for(i in 1:d1) { WW[i,] <- mi[i]*dnorm(dataF2,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF2,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF2
    ###########obtain mean#########################
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    #############restriction#######################
    aa1<-sumx[1]/m_sam[1]+2*sumx[2]/m_sam[2]+sumx[3]/m_sam[3]-sumwx[1]/n0[1]-2*sumwx[2]/n0[2]-sumwx[3]/n0[3]
    aa2<-sigma[1]/m_sam[1]+4*sigma[2]/m_sam[2]+sigma[3]/m_sam[3]+sigma[4]/n0[1]+4.0*sigma[5]/n0[2]+sigma[6]/n0[3]
    aa3<-aa1/aa2
    m[1]<-(sumx[1]-sigma[1]*aa3)/m_sam[1];m[2]<-(sumx[2]-sigma[2]*aa3*2)/m_sam[2]
    m[3]<-(sumx[3]-sigma[3]*aa3)/m_sam[3];m[4]<-(sumwx[1]+sigma[4]*aa3)/n0[1]
    m[5]<-(sumwx[2]+sigma[5]*aa3*2)/n0[2];m[6]<-(sumwx[3]+sigma[6]*aa3)/n0[3]
    ##########obtain variance####################################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d1) {swx[i] <- WW[i,]%*%(dataF2-m[i+3])^2 }
    s0[1]<-swx[1]+swx[2]+swx[3]
    sigma[4]<-s0[1]/m_sam[4]
    
    sigma00<-sigma[4]-sigma[1]
    if (sigma00<0){sigma00<-0;sigma[4]<-sigma[1]}
    sigma[5]<-sigma[6]<-sigma[4]
    aa0<-sigma[1]
    aa2<-ss1+ss2+ss3;aa3<-m_sam[1]+m_sam[2]+m_sam[3]
    n_iter<-0
    aa4<-1000
    while(aa4>0.0001){
      n_iter<-n_iter+1
      aa1<-sigma[1]/(sigma[1]+sigma00)
      if (aa1>1) {aa1<-1}
      sigma[1]<-(aa2+aa1^2*s0[1])/(aa3+aa1*m_sam[4])
      aa4<-abs(sigma[1]-aa0)
      aa0<-sigma[1]
      if(n_iter>20)break
    }
    sigma[2]<-sigma[3]<-sigma[1]
    sigma[4]<-sigma[1]+sigma00;sigma[5]<-sigma[4];sigma[6]<-sigma[4]
    ########criteria for iterations to stop####################
    pi<-m[c(4:6)];gh<-sigma[c(4:6)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mix_pi))) 
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*7 
  #########first order genetic parameter process##########
  hh13 <- matrix(c(1,1,1,1,1,1,1,0,-1,1,0,-1,0,1,0,0,1,0,1,0,-1,0,0,0,0,1,0,0.5,0.5,0.5),6,5)
  B13<-solve(crossprod(hh13,hh13))%*%crossprod(hh13,m)
  #################second oder genetic parameter process#############
  jj<-sigmaF2 - sigma[4]
  gg<-sigma[4]-sigma[1]
  if(jj<0){jj<-0}
  if(gg<0){gg<-0}
  ll<-jj/sigmaF2
  rr<-gg/sigmaF2
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
  
  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2); 
  F2w1<-1/(12*m_sam[4])
  F2bmw <- matrix(0,m_sam[4],1); F2bmwsl <- matrix(0,m_sam[4],d1)
  for(i in 1:d1){
    F2gg <- (dataF2 - pi[i])/sqrt(gh[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < m_sam[4]){F2P2 <- F2P2+runif(m_sam[4])/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*m_sam[4]) + sum((F2P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F2u <- as.matrix(c(12*m_sam[4]*((F2dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F2dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F2dd[3]/m_sam[4]-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),(1-pkolm(F2D,m_sam[4]))))
  
  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("MX1-AD-AD",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," "," ",round(sigma[1],4),round(sigma[4],4),
                       round(t(mix_pi),4)," "," "," "," "," "," ",round(B13[1],4)," "," "," ",round(B13[2],4)," ",round(B13[3],4)," "," "," "," "," ",round(B13[4],4),round(B13[5],4),round(jj,4),round(ll*100,4),round(gg,4),round(rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}

##############MX1-A-AD(D-2)############################ 
G4F2ModelFun[[15]] <- function(K1,logL,df11,df21,df31,df41){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF2 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2)))
  m_esp<-0.0001
  sigmaF2<-var(dataF2);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a1<-sqrt(sigmaF2/m_sam[4])
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2)))
  ###############procedure start##############################        
  mi<-as.matrix(c(0.25,0.5,0.25))
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaF2/2,sigmaF2/2,sigmaF2/2))
  m[c(1:4)]<-m
  if(m[1]<m[3]) {a1<--a1}
  m[5]<-m[4];m[4]<-m[5]+2*a1;m[6]<-m[5]-2*a1
  pi<-m[c(4:6)];gh<-sigma[c(4:6)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mi)))
  ########iteration process############
  d1<-3
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,3,m_sam[4]); swx <- matrix(0,3,1);rr<-matrix(0,1,1);s0<-matrix(0,1,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:6)];ssigma<-sigma[c(4:6)]
    for(i in 1:d1) { WW[i,] <- mi[i]*dnorm(dataF2,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF2,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF2
    ###########obtain mean#########################
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    #############restriction#######################
    aa1<-sumx[1]/m_sam[1]+2*sumx[2]/m_sam[2]+sumx[3]/m_sam[3]-4*sumwx[2]/n0[2]
    aa2<-sumwx[1]/n0[1]-2*sumwx[2]/n0[2]+sumwx[3]/n0[3]
    aa3<-sigma[1]/m_sam[1]+4*sigma[2]/m_sam[2]+sigma[3]/m_sam[3]+16*sigma[5]/n0[2]
    aa4<-8*sigma[5]/n0[2]
    aa5<-sigma[4]/n0[1]+4*sigma[5]/n0[2]+sigma[6]/n0[3]
    aa6<-aa3*aa5-aa4^2;aa7<-aa1*aa5-aa2*aa4;aa8<-aa2*aa3-aa1*aa4
    rr[1]<-aa7/aa6;rr[2]<-aa8/aa6
    m[1]<-(sumx[1]-sigma[1]*rr[1])/m_sam[1];m[2]<-(sumx[2]-sigma[2]*rr[1]*2)/m_sam[2]
    m[3]<-(sumx[3]-sigma[3]*rr[1])/m_sam[3];m[4]<-(sumwx[1]-sigma[4]*rr[2])/n0[1]
    m[5]<-(sumwx[2]+sigma[5]*(4*rr[1]+2*rr[2]))/n0[2];m[6]<-(sumwx[3]-sigma[6]*rr[2])/n0[3]
    ##########obtain variance####################################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d1) {swx[i] <- WW[i,]%*%(dataF2-m[i+3])^2 }
    s0[1]<-swx[1]+swx[2]+swx[3]
    sigma[4]<-s0[1]/m_sam[4]
    
    sigma00<-sigma[4]-sigma[1]
    if (sigma00<0) {sigma00<-0;sigma[4]<-sigma[1]}
    sigma[5]<-sigma[6]<-sigma[4]
    aa0<-sigma[1]
    aa2<-ss1+ss2+ss3;aa3<-m_sam[1]+m_sam[2]+m_sam[3]
    n_iter<-0
    aa4<-1000
    while(aa4>0.0001){
      n_iter<-n_iter+1 
      aa1<-sigma[1]/(sigma[1]+sigma00)
      if (aa1>1) aa1<-1
      sigma[1]<-(aa2+aa1^2*s0[1])/(aa3+aa1*m_sam[4])
      aa4<-abs(sigma[1]-aa0)
      aa0<-sigma[1]
      if(n_iter>20)break
    }
    sigma[2]<-sigma[3]<-sigma[1]
    sigma[4]<-sigma[1]+sigma00;sigma[5]<-sigma[4];sigma[6]<-sigma[4]
    ########criteria for iterations to stop####################
    pi<-m[c(4:6)];gh<-sigma[c(4:6)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mix_pi))) 
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*6
  #########first order genetic parameter process##########
  hh14 <- matrix(c(1,1,1,1,1,1,1,0,-1,1,0,-1,1,0,-1,0,0,0,0,1,0,0.5,0.5,0.5),6,4)
  B14<-solve(crossprod(hh14,hh14))%*%crossprod(hh14,m)
  #################second oder genetic parameter process#############
  jj<-sigmaF2 - sigma[4]
  gg<-sigma[4]-sigma[1]
  if(jj<0){jj<-0}
  if(gg<0){gg<-0}
  ll<-jj/sigmaF2
  rr<-gg/sigmaF2
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
  
  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2); 
  F2w1<-1/(12*m_sam[4])
  F2bmw <- matrix(0,m_sam[4],1); F2bmwsl <- matrix(0,m_sam[4],d1)
  for(i in 1:d1){
    F2gg <- (dataF2 - pi[i])/sqrt(gh[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < m_sam[4]){F2P2 <- F2P2+runif(m_sam[4])/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*m_sam[4]) + sum((F2P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F2u <- as.matrix(c(12*m_sam[4]*((F2dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F2dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F2dd[3]/m_sam[4]-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),(1-pkolm(F2D,m_sam[4]))))
  
  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("MX1-A-AD",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," "," ",round(sigma[1],4),round(sigma[4],4),
                       round(t(mix_pi),4)," "," "," "," "," "," ",round(B14[1],4)," "," "," ",round(B14[2],4)," "," "," "," "," "," "," ",round(B14[3],4),round(B14[4],4),round(jj,4),round(ll*100,4),round(gg,4),round(rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}


###############MX1-EAD-AD(D-3)#####################   
G4F2ModelFun[[16]] <- function(K1,logL,df11,df21,df31,df41){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF2 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2)))
  m_esp<-0.0001
  sigmaF2<-var(dataF2);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a1<-sqrt(sigmaF2/m_sam[4])
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2)))
  ###############procedure start##############################           
  d2<-2
  mi<-as.matrix(c(0.75,0.25))
  m[c(1:4)]<-m
  if(m[1]<m[3]){a1<-a1}
  m[5]<-m[4]-1.5*a1;m[4]<-m[4]+1.5*a1
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaF2/2,sigmaF2/2))
  pi<-m[c(4:5)];gh<-sigma[c(4:5)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mi)))
  ########iteration process############
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,2,m_sam[4]); swx <- matrix(0,2,1);  s0<-matrix(0,1,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:5)];ssigma<-sigma[c(4:5)]
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataF2,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF2,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF2
    ###########obtain mean#########################
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    ##############restriction#####################
    aa1<-sumx[1]/m_sam[1]+2*sumx[2]/m_sam[2]+sumx[3]/m_sam[3]-3*sumwx[1]/n0[1]-sumwx[2]/n0[2]
    aa2<-sigma[1]/m_sam[1]+4*sigma[2]/m_sam[2]+sigma[3]/m_sam[3]+9*sigma[4]/n0[1]+sigma[5]/n0[2]
    aa3<-aa1/aa2
    m[1]<-(sumx[1]-sigma[1]*aa3)/m_sam[1];m[2]<-(sumx[2]-sigma[2]*aa3*2)/m_sam[2]
    m[3]<-(sumx[3]-sigma[3]*aa3)/m_sam[3];m[4]<-(sumwx[1]+sigma[4]*aa3*3)/n0[1]
    m[5]<-(sumwx[2]+sigma[5]*aa3)/n0[2] 
    ##########obtain variance####################################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d2) {swx[i] <- WW[i,]%*%(dataF2-m[i+3])^2 }
    s0[1]<-swx[1]+swx[2]
    sigma[4]<-s0[1]/m_sam[4]
    
    sigma00<-sigma[4]-sigma[1]
    if (sigma00<0) {sigma00<-0;sigma[4]<-sigma[1]}
    sigma[5]<-sigma[4]
    aa0<-sigma[1]
    aa2<-ss1+ss2+ss3;aa3<-m_sam[1]+m_sam[2]+m_sam[3]
    n_iter<-0
    aa4<-1000
    while(aa4>0.0001){
      n_iter<-n_iter+1
      aa1<-sigma[1]/(sigma[1]+sigma00)
      if (aa1>1) {aa1<-1}
      sigma[1]<-(aa2+aa1^2*s0[1])/(aa3+aa1*m_sam[4])
      aa4<-abs(sigma[1]-aa0)
      aa0<-sigma[1]
      if(n_iter>20)break
    }
    sigma[2]<-sigma[3]<-sigma[1]
    sigma[4]<-sigma[1]+sigma00;sigma[5]<-sigma[4]
    ########criteria for iterations to stop####################
    pi<-m[c(4:5)];gh<-sigma[c(4:5)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mix_pi))) 
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*6 
  #########first order genetic parameter process##########
  hh15 <- matrix(c(1,1,1,1,1,1,1,-1,1,-1,1,0,-1,0,0,0,1,0,0.5,0.5),5,4)
  B15<-solve(crossprod(hh15,hh15))%*%crossprod(hh15,m)
  #################second oder genetic parameter process#############
  jj<-sigmaF2 - sigma[4]
  gg<-sigma[4]-sigma[1]
  if(jj<0){jj<-0}
  if(gg<0){gg<-0}
  ll<-jj/sigmaF2
  rr<-gg/sigmaF2
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
  
  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2); 
  F2w1<-1/(12*m_sam[4])
  F2bmw <- matrix(0,m_sam[4],1); F2bmwsl <- matrix(0,m_sam[4],d2)
  for(i in 1:d2){
    F2gg <- (dataF2 - pi[i])/sqrt(gh[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < m_sam[4]){F2P2 <- F2P2+runif(m_sam[4])/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*m_sam[4]) + sum((F2P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F2u <- as.matrix(c(12*m_sam[4]*((F2dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F2dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F2dd[3]/m_sam[4]-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),(1-pkolm(F2D,m_sam[4]))))
  
  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("MX1-EAD-AD",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," "," "," ",round(sigma[1],4),round(sigma[4],4),
                       round(t(mix_pi),4)," "," "," "," "," "," "," ",round(B15[1],4)," "," "," ",round(B15[2],4)," "," "," "," "," "," "," ",round(B15[3],4),round(B15[4],4),round(jj,4),round(ll*100,4),round(gg,4),round(rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}

############MX1-NCD-AD(D-4)##################################  
G4F2ModelFun[[17]] <- function(K1,logL,df11,df21,df31,df41){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF2 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2)))
  m_esp<-0.0001
  sigmaF2<-var(dataF2);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a1<-sqrt(sigmaF2/m_sam[4])
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2)))
  ###############procedure start##############################          
  d2<-2
  mi<-as.matrix(c(0.25,0.75))
  m[c(1:4)]<-m
  if(m[1]<m[3]){a1<-a1}
  m[5]<-m[4]-1.5*a1;m[4]<-m[4]+1.5*a1
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaF2/2,sigmaF2/2))
  pi<-m[c(4:5)];gh<-sigma[c(4:5)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mi)))
  ########iteration process############
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,2,m_sam[4]); swx <- matrix(0,2,1);s0<-matrix(0,1,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:5)];ssigma<-sigma[c(4:5)]
    for(i in 1:d2) { WW[i,] <- mi[i]*dnorm(dataF2,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF2,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF2
    ###########obtain mean#########################
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    ##############restriction#####################
    aa1<-sumx[1]/m_sam[1]+2*sumx[2]/m_sam[2]+sumx[3]/m_sam[3]-sumwx[1]/n0[1]-3*sumwx[2]/n0[2]
    aa2<-sigma[1]/m_sam[1]+4*sigma[2]/m_sam[2]+sigma[3]/m_sam[3]+sigma[4]/n0[1]+9*sigma[5]/n0[2]
    aa3<-aa1/aa2
    m[1]<-(sumx[1]-sigma[1]*aa3)/m_sam[1];m[2]<-(sumx[2]-sigma[2]*aa3*2)/m_sam[2]
    m[3]<-(sumx[3]-sigma[3]*aa3)/m_sam[3];m[4]<-(sumwx[1]+sigma[4]*aa3)/n0[1]
    m[5]<-(sumwx[2]+sigma[5]*aa3*3)/n0[2]
    ##########obtain variance####################################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d2) {swx[i] <- WW[i,]%*%(dataF2-m[i+3])^2 }
    s0[1]<-swx[1]+swx[2]
    sigma[4]<-s0[1]/m_sam[4]
    
    sigma00<-sigma[4]-sigma[1]
    if (sigma00<0) {sigma00<-0;sigma[4]<-sigma[1]}
    sigma[5]<-sigma[4]
    aa0<-sigma[1]
    aa2<-ss1+ss2+ss3;aa3<-m_sam[1]+m_sam[2]+m_sam[3]
    n_iter<-0
    aa4<-1000
    while(aa4>0.0001){
      n_iter<-n_iter+1
      aa1<-sigma[1]/(sigma[1]+sigma00)
      if (aa1>1) {aa1<-1}
      sigma[1]<-(aa2+aa1^2*s0[1])/(aa3+aa1*m_sam[4])
      aa4<-abs(sigma[1]-aa0)
      aa0<-sigma[1]  
      if(n_iter>20)break
    }
    sigma[2]<-sigma[3]<-sigma[1]
    sigma[4]<-sigma[1]+sigma00;sigma[5]<-sigma[4]  
    ########criteria for iterations to stop####################
    pi<-m[c(4:5)];gh<-sigma[c(4:5)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mix_pi))) 
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*6 
  #########first order genetic parameter process##########
  hh16 <- matrix(c(1,1,1,1,1,1,-1,-1,1,-1,1,0,-1,0,0,0,1,0,0.5,0.5),5,4)
  B16<-solve(crossprod(hh16,hh16))%*%crossprod(hh16,m)
  #################second oder genetic parameter process#############
  jj<-sigmaF2 - sigma[4]
  gg<-sigma[4]-sigma[1]
  if(jj<0){jj<-0}
  if(gg<0){gg<-0}
  ll<-jj/sigmaF2
  rr<-gg/sigmaF2   
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
  
  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2); 
  F2w1<-1/(12*m_sam[4])
  F2bmw <- matrix(0,m_sam[4],1); F2bmwsl <- matrix(0,m_sam[4],d2)
  for(i in 1:d2){
    F2gg <- (dataF2 - pi[i])/sqrt(gh[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < m_sam[4]){F2P2 <- F2P2+runif(m_sam[4])/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*m_sam[4]) + sum((F2P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F2u <- as.matrix(c(12*m_sam[4]*((F2dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F2dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F2dd[3]/m_sam[4]-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),(1-pkolm(F2D,m_sam[4]))))
  
  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("MX1-NCD-AD",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," "," "," ",round(sigma[1],4),round(sigma[4],4),
                       round(t(mix_pi),4)," "," "," "," "," "," "," ",round(B16[1],4)," "," "," ",round(B16[2],4)," "," "," "," "," "," "," ",round(B16[3],4),round(B16[4],4),round(jj,4),round(ll*100,4),round(gg,4),round(rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}

#############MX2-ADI-ADI(E-0)##########################
G4F2ModelFun[[18]] <- function(K1,logL,df11,df21,df31,df41){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF2 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2)))
  m_esp<-0.0001
  sigmaF2<-var(dataF2);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a2<-sqrt(sigmaF2/m_sam[4]*(m_sam[4]-1))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2)))
  ###############procedure start##############################      
  d3<-9
  mi<-as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaF2/2,sigmaF2/2,sigmaF2/2,sigmaF2/2,sigmaF2/2,sigmaF2/2,sigmaF2/2,sigmaF2/2,sigmaF2/2))
  m[c(1:3)]<-m[c(1:3)]
  if (m[1]<m[3]) {a2<--a2}
  m[5]<-m[4]+2.4*a2;m[6]<-m[4]+1.6*a2;m[7]<-m[4]+0.8*a2
  m[8]<-m[4];m[9]<-m[4]-0.8*a2;m[10]<-m[4]-1.6*a2
  m[11]<-m[4]-2.4*a2;m[12]<-m[4]-3.2*a2;m[4]<-m[4]+3.2*a2
  pi<-m[c(4:12)];gh<-sigma[c(4:12)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mi)))
  ########iteration process############
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,9,m_sam[4]);swx <- matrix(0,9,1); s0<-matrix(0,1,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:12)];ssigma<-sigma[c(4:12)]
    for(i in 1:d3) { WW[i,] <- mi[i]*dnorm(dataF2,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF2,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF2
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    m[c(1:3)]<-sumx[c(1:3)]/m_sam[c(1:3)];m[c(4:12)]<-sumwx[c(1:9)]/n0[c(1:9)]
    ######obtain variance##############
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d3) {swx[i] <- WW[i,]%*%(dataF2-m[i+3])^2 }
    s0[1]<-sum(swx)
    sigma[4]<-s0[1]/m_sam[4]
    
    sigma00<-sigma[4]-sigma[1]
    if (sigma00<0){sigma00<-0;sigma[4]<-sigma[1]}
    sigma[c(5:12)]<-sigma[4]
    aa1<-ss1+ss2+ss3;aa2<-m_sam[1]+m_sam[2]+m_sam[3]
    n_iter<-0
    aa0<-sigma[1]
    aa4<-1000
    while(aa4>0.0001){
      n_iter<-n_iter+1
      aa3<-sigma[1]/(sigma[1]+sigma00)
      sigma[1]<-(aa1+aa3^2*s0[1])/(aa2+aa3*m_sam[4])
      aa4<-abs(sigma[1]-aa0)
      aa0<-sigma[1]
      if (n_iter>20) break 
    }
    sigma[2]<-sigma[3]<-sigma[1];sigma[c(4:12)]<-sigma[1]+sigma00
    ########criteria for iterations to stop####################
    pi<-m[c(4:12)];gh<-sigma[c(4:12)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mix_pi))) 
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*12
  #########first order genetic parameter process##########
  hh17 <- matrix(c(1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,
                   0,-1,1,1,1,0,0,0,-1,-1,-1,1,0,-1,1,0,-1,1,0,-1,1,0,-1,0,1,0,0,0,0,1,1,1,0,0,0,0,1,0,0,1,0,0,1,0,0,1,0,1,0,1,1,0,-1,
                   0,0,0,-1,0,1,0,0,0,0,1,0,0,0,0,0,-1,0,0,0,0,0,0,0,1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0),12,12)
  B17<-solve(hh17,m)
  #################second oder genetic parameter process#############
  jj<-sigmaF2 - sigma[4]
  gg<-sigma[4]-sigma[1]
  if(jj<0){jj<-0}
  if(gg<0){gg<-0}
  ll<-jj/sigmaF2
  rr<-gg/sigmaF2   
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
  
  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2); 
  F2w1<-1/(12*m_sam[4])
  F2bmw <- matrix(0,m_sam[4],1); F2bmwsl <- matrix(0,m_sam[4],d3)
  for(i in 1:d3){
    F2gg <- (dataF2 - pi[i])/sqrt(gh[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < m_sam[4]){F2P2 <- F2P2+runif(m_sam[4])/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*m_sam[4]) + sum((F2P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F2u <- as.matrix(c(12*m_sam[4]*((F2dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F2dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F2dd[3]/m_sam[4]-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),(1-pkolm(F2D,m_sam[4]))))
  
  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("MX2-ADI-ADI",round(abc,4),round(AIC,4),round(t(m),4),round(sigma[1],4),round(sigma[4],4),
                       round(t(mix_pi),4),round(B17[1],4),round(B17[2],4),round(B17[3],4),round(B17[4],4),round(B17[5],4),round(B17[6],4),round(B17[7],4),round(B17[8],4),round(B17[9],4),round(B17[10],4),round(B17[11],4),round(B17[12],4)," "," ",round(jj,4),round(ll*100,4),round(gg,4),round(rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}

#################MX2-ADI-AD(E-1)###########################
G4F2ModelFun[[19]] <- function(K1,logL,df11,df21,df31,df41){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF2 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2)))
  m_esp<-0.0001
  sigmaF2<-var(dataF2);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a2<-sqrt(sigmaF2/m_sam[4]*(m_sam[4]-1))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2)))  
  ###############procedure start##############################   
  d3<-9
  mi<-as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaF2/2,sigmaF2/2,sigmaF2/2,sigmaF2/2,sigmaF2/2,sigmaF2/2,sigmaF2/2,sigmaF2/2,sigmaF2/2))
  m[c(1:3)]<-m[c(1:3)]
  if (m[1]<m[3]) {a2<--a2}
  m[5]<-m[4]+2.4*a2;m[6]<-m[4]+1.6*a2;m[7]<-m[4]+0.8*a2
  m[8]<-m[4];m[9]<-m[4]-0.8*a2;m[10]<-m[4]-1.6*a2
  m[11]<-m[4]-2.4*a2;m[12]<-m[4]-3.2*a2;m[4]<-m[4]+3.2*a2
  pi<-m[c(4:12)];gh<-sigma[c(4:12)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mi)))
  ########iteration process############
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,9,m_sam[4]); swx <- matrix(0,9,1); s0<-matrix(0,1,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:12)];ssigma<-sigma[c(4:12)]
    for(i in 1:d3) { WW[i,] <- mi[i]*dnorm(dataF2,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF2,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF2
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    aa1<-sumx[1]/m_sam[1]+2*sumx[2]/m_sam[2]+sumx[3]/m_sam[3]-sumwx[1]/n0[1]-2*sumwx[5]/n0[5]-sumwx[9]/n0[9]
    aa2<-sigma[1]/m_sam[1]+4*sigma[2]/m_sam[2]+sigma[3]/m_sam[3]+sigma[4]/n0[1]+4*sigma[8]/n0[5]+sigma[12]/n0[9]
    aa3<-aa1/aa2   
    ###########obtain mean#######################
    m[1]<-(sumx[1]-sigma[1]*aa3)/m_sam[1];m[2]<-(sumx[2]-sigma[2]*aa3*2)/m_sam[2];m[3]<-(sumx[3]-sigma[3]*aa3)/m_sam[3]       
    m[4]<-(sumwx[1]+sigma[4]*aa3)/n0[1];m[5]<-sumwx[2]/n0[2];m[6]<-sumwx[3]/n0[3]
    m[7]<-sumwx[4]/n0[4];m[8]<-(sumwx[5]+sigma[8]*aa3*2)/n0[5];m[9]<-sumwx[6]/n0[6]
    m[10]<-sumwx[7]/n0[7];m[11]<-sumwx[8]/n0[8];m[12]<-(sumwx[9]+sigma[12]*aa3)/n0[9]
    ######obtain variance##############
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d3) {swx[i] <- WW[i,]%*%(dataF2-m[i+3])^2 }
    s0[1]<-sum(swx);sigma[4]<-s0[1]/m_sam[4]
    
    sigma00<-sigma[4]-sigma[1]
    if (sigma00<0){sigma00<-0;sigma[4]<-sigma[1]}
    sigma[c(5:12)]<-sigma[4]
    aa1<-ss1+ss2+ss3;aa2<-m_sam[1]+m_sam[2]+m_sam[3]
    n_iter<-0
    aa0<-sigma[1]
    aa4<-1000
    while(aa4>0.0001){
      n_iter<-n_iter+1
      aa3<-sigma[1]/(sigma[1]+sigma00)
      sigma[1]<-(aa1+aa3^2*s0[1])/(aa2+aa3*m_sam[4])
      aa4<-abs(sigma[1]-aa0)
      aa0<-sigma[1]
      if (n_iter>20) break 
    }
    sigma[2]<-sigma[3]<-sigma[1]
    sigma[c(4:12)]<-sigma[1]+sigma00
    ########criteria for iterations to stop####################
    pi<-m[c(4:12)];gh<-sigma[c(4:12)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mix_pi))) 
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*9
  #########first order genetic parameter process##########
  hh18 <- matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,0,-1,1,1,1,0,0,0,-1,-1,-1,1,0,-1,1,0,-1,1,0,-1,1,0,-1,
                   0,1,0,0,0,0,1,1,1,0,0,0,0,1,0,0,1,0,0,1,0,0,1,0,1,0,1,1,0,-1,0,0,0,-1,0,1,0,0,0,0,
                   1,0,0,0,0,0,-1,0,0,0,0,0,0,0,1,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,-1,0,0,0,0,0,
                   0,0,0,0,0,1,0,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5),12,11)
  B18<-solve(crossprod(hh18,hh18))%*%crossprod(hh18,m)
  #################second oder genetic parameter process#############
  jj<-sigmaF2 - sigma[4]
  gg<-sigma[4]-sigma[1]
  if(jj<0){jj<-0}
  if(gg<0){gg<-0}
  ll<-jj/sigmaF2
  rr<-gg/sigmaF2   
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
  
  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2); 
  F2w1<-1/(12*m_sam[4])
  F2bmw <- matrix(0,m_sam[4],1); F2bmwsl <- matrix(0,m_sam[4],d3)
  for(i in 1:d3){
    F2gg <- (dataF2 - pi[i])/sqrt(gh[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < m_sam[4]){F2P2 <- F2P2+runif(m_sam[4])/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*m_sam[4]) + sum((F2P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F2u <- as.matrix(c(12*m_sam[4]*((F2dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F2dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F2dd[3]/m_sam[4]-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),(1-pkolm(F2D,m_sam[4]))))
  
  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("MX2-ADI-AD",round(abc,4),round(AIC,4),round(t(m),4),round(sigma[1],4),round(sigma[4],4),
                       round(t(mix_pi),4),round(B18[1],4)," "," "," ",round(B18[2],4),round(B18[3],4),round(B18[4],4),round(B18[5],4),round(B18[6],4),round(B18[7],4),round(B18[8],4),round(B18[9],4),round(B18[10],4),round(B18[11],4),round(jj,4),round(ll*100,4),round(gg,4),round(rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}

###########MX2-AD-AD(E-2)####################################
G4F2ModelFun[[20]] <- function(K1,logL,df11,df21,df31,df41){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF2 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2)))
  m_esp<-0.0001
  sigmaF2<-var(dataF2);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a2<-sqrt(sigmaF2/m_sam[4]*(m_sam[4]-1))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2)))
  ###############procedure start##############################      
  d3<-9
  mi<-as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaF2/2,sigmaF2/2,sigmaF2/2,sigmaF2/2,sigmaF2/2,sigmaF2/2,sigmaF2/2,sigmaF2/2,sigmaF2/2))
  m[c(1:3)]<-m[c(1:3)]
  if (m[1]<m[3]) {a2<--a2}
  m[5]<-m[4]+2.4*a2;m[6]<-m[4]+1.6*a2;m[7]<-m[4]+0.8*a2
  m[8]<-m[4];m[9]<-m[4]-0.8*a2;m[10]<-m[4]-1.6*a2
  m[11]<-m[4]-2.4*a2;m[12]<-m[4]-3.2*a2;m[4]<-m[4]+3.2*a2
  pi<-m[c(4:12)];gh<-sigma[c(4:12)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mi)))
  ########iteration process############
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,9,m_sam[4]); swx <- matrix(0,9,1); s0<-matrix(0,1,1)
  hh19<-matrix(0,5,5);b_line19<-matrix(0,5,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:12)];ssigma<-sigma[c(4:12)]
    for(i in 1:d3) { WW[i,] <- mi[i]*dnorm(dataF2,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF2,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF2
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    ############################################################
    hh19[1,1]<-sigma[1]/m_sam[1]+4*sigma[2]/m_sam[2]+sigma[3]/m_sam[3]+sigma[4]/n0[1]+4.0*sigma[4]/n0[5]+sigma[4]/n0[9]
    hh19[1,2]<--sigma[4]/n0[1]-2*sigma[4]/n0[5]
    hh19[1,3]<--sigma[4]/n0[1]
    hh19[1,4]<--sigma[4]/n0[1]-sigma[4]/n0[9]
    hh19[1,5]<--2*sigma[4]/n0[5]-sigma[4]/n0[9]
    hh19[2,2]<-sigma[4]/n0[1]+sigma[4]/n0[2]+sigma[4]/n0[4]+sigma[4]/n0[5]
    hh19[2,3]<-sigma[4]/n0[1]+sigma[4]/n0[2]
    hh19[2,4]<-sigma[4]/n0[1]
    hh19[2,5]<-sigma[4]/n0[5]
    hh19[3,3]<-sigma[4]/n0[1]+sigma[4]/n0[2]+sigma[4]/n0[7]+sigma[4]/n0[8]
    hh19[3,4]<-sigma[4]/n0[1]+sigma[4]/n0[7]
    hh19[3,5]<--sigma[4]/n0[8]
    hh19[4,4]<-sigma[4]/n0[1]+sigma[4]/n0[3]+sigma[4]/n0[7]+sigma[4]/n0[9]
    hh19[4,5]<-sigma[4]/n0[9]
    hh19[5,5]<-sigma[4]/n0[5]+sigma[4]/n0[6]+sigma[4]/n0[8]+sigma[4]/n0[9]
    for(i in 2:5)
    {
      for(j in 1:(i-1))
      {
        hh19[i,j]<-hh19[j,i]
      }
    }
    ##############################################################
    b_line19[1]<-sumx[1]/m_sam[1]+2*sumx[2]/m_sam[2]+sumx[3]/m_sam[3]-sumwx[1]/n0[1]-2*sumwx[5]/n0[5]-sumwx[9]/n0[9]
    b_line19[2]<-sumwx[1]/n0[1]-sumwx[2]/n0[2]-sumwx[4]/n0[4]+sumwx[5]/n0[5]
    b_line19[3]<-sumwx[1]/n0[1]-sumwx[2]/n0[2]-sumwx[7]/n0[7]+sumwx[8]/n0[8]
    b_line19[4]<-sumwx[1]/n0[1]-sumwx[3]/n0[3]-sumwx[7]/n0[7]+sumwx[9]/n0[9]
    b_line19[5]<-sumwx[5]/n0[5]-sumwx[6]/n0[6]-sumwx[8]/n0[8]+sumwx[9]/n0[9]
    B19<-solve(hh19,b_line19)
    ###############obtain mean#####################
    m[1]<-(sumx[1]-sigma[1]*B19[1])/m_sam[1]
    m[2]<-(sumx[2]-sigma[2]*B19[1]*2)/m_sam[2]
    m[3]<-(sumx[3]-sigma[3]*B19[1])/m_sam[3]
    m[4]<-(sumwx[1]+sigma[4]*(B19[1]-B19[2]-B19[3]-B19[4]))/n0[1]
    m[5]<-(sumwx[2]+sigma[5]*(B19[2]+B19[3]))/n0[2]
    m[6]<-(sumwx[3]+sigma[6]*B19[4])/n0[3]
    m[7]<-(sumwx[4]+sigma[7]*B19[2])/n0[4]
    m[8]<-(sumwx[5]+sigma[8]*(2*B19[1]-B19[2]-B19[5]))/n0[5]
    m[9]<-(sumwx[6]+sigma[9]*B19[5])/n0[6]
    m[10]<-(sumwx[7]+sigma[10]*(B19[3]+B19[4]))/n0[7]
    m[11]<-(sumwx[8]+sigma[11]*(-B19[3]+B19[5]))/n0[8]
    m[12]<-(sumwx[9]+sigma[12]*(B19[1]-B19[4]-B19[5]))/n0[9]
    ######obtain variance##############
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d3) {swx[i] <- WW[i,]%*%(dataF2-m[i+3])^2 }
    s0[1]<-sum(swx);sigma[4]<-s0[1]/m_sam[4]
    
    sigma00<-sigma[4]-sigma[1]
    if (sigma00<0){sigma00<-0;sigma[4]<-sigma[1]}
    sigma[c(5:12)]<-sigma[4]
    aa1<-ss1+ss2+ss3;aa2<-m_sam[1]+m_sam[2]+m_sam[3]
    n_iter<-0
    aa0<-sigma[1]
    aa4<-1000
    while(aa4>0.0001){
      n_iter<-n_iter+1
      aa3<-sigma[1]/(sigma[1]+sigma00)
      sigma[1]<-(aa1+aa3^2*s0[1])/(aa2+aa3*m_sam[4])
      aa4<-abs(sigma[1]-aa0)
      aa0<-sigma[1]
      if (n_iter>20) break 
    }
    sigma[2]<-sigma[3]<-sigma[1]
    sigma[c(4:12)]<-sigma[1]+sigma00
    ########criteria for iterations to stop####################
    pi<-m[c(4:12)];gh<-sigma[c(4:12)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mix_pi))) 
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*5
  #########first order genetic parameter process##########
  hh191 <- matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,0,-1,1,1,1,0,0,0,-1,-1,-1,1,0,-1,1,0,-1,1,
                    0,-1,1,0,-1,0,1,0,0,0,0,1,1,1,0,0,0,0,1,0,0,1,0,0,1,0,0,1,0,1,0,-1,0,0,
                    0,0,0,0,0,0,0,0,1,0,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5),12,7)
  B191<-solve(crossprod(hh191,hh191))%*%crossprod(hh191,m)
  #################second oder genetic parameter process#############
  jj<-sigmaF2 - sigma[4]
  gg<-sigma[4]-sigma[1]
  if(jj<0){jj<-0}
  if(gg<0){gg<-0}
  ll<-jj/sigmaF2
  rr<-gg/sigmaF2   
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
  
  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2); 
  F2w1<-1/(12*m_sam[4])
  F2bmw <- matrix(0,m_sam[4],1); F2bmwsl <- matrix(0,m_sam[4],d3)
  for(i in 1:d3){
    F2gg <- (dataF2 - pi[i])/sqrt(gh[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < m_sam[4]){F2P2 <- F2P2+runif(m_sam[4])/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*m_sam[4]) + sum((F2P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F2u <- as.matrix(c(12*m_sam[4]*((F2dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F2dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F2dd[3]/m_sam[4]-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),(1-pkolm(F2D,m_sam[4]))))
  
  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("MX2-AD-AD",round(abc,4),round(AIC,4),round(t(m),4),round(sigma[1],4),round(sigma[4],4),
                       round(t(mix_pi),4),round(B191[1],4)," "," "," ",round(B191[2],4),round(B191[3],4),round(B191[4],4),round(B191[5],4)," "," "," "," ",round(B191[6],4),round(B191[7],4),round(jj,4),round(ll*100,4),round(gg,4),round(rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}

#####################MX2-A-AD(E-3)##########################
G4F2ModelFun[[21]] <- function(K1,logL,df11,df21,df31,df41){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF2 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2)))
  m_esp<-0.0001
  sigmaF2<-var(dataF2);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a2<-sqrt(sigmaF2/m_sam[4]*(m_sam[4]-1))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2)))   
  ###############procedure start##############################    
  d3<-9
  mi<-as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaF2/2,sigmaF2/2,sigmaF2/2,sigmaF2/2,sigmaF2/2,sigmaF2/2,sigmaF2/2,sigmaF2/2,sigmaF2/2))
  m[c(1:3)]<-m[c(1:3)]
  if (m[1]<m[3]) {a2<--a2}
  m[5]<-m[4]+2.4*a2;m[6]<-m[4]+1.6*a2;m[7]<-m[4]+0.8*a2
  m[8]<-m[4];m[9]<-m[4]-0.8*a2;m[10]<-m[4]-1.6*a2
  m[11]<-m[4]-2.4*a2;m[12]<-m[4]-3.2*a2;m[4]<-m[4]+3.2*a2
  pi<-m[c(4:12)];gh<-sigma[c(4:12)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mi)))
  ########iteration process############
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,9,m_sam[4]); swx <- matrix(0,9,1);s0<-matrix(0,1,1)
  hh20<-matrix(0,7,7);b_line20<-matrix(0,7,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:12)];ssigma<-sigma[c(4:12)]
    for(i in 1:d3) { WW[i,] <- mi[i]*dnorm(dataF2,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF2,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF2
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    ########################################  
    hh20[1,1]<-sigma[1]/m_sam[1]+4*sigma[2]/m_sam[2]+sigma[3]/m_sam[3]+16*sigma[4]/n0[5]
    hh20[1,2]<-0
    hh20[1,3]<-0
    hh20[1,4]<-8*sigma[4]/n0[5]
    hh20[1,5]<-8*sigma[4]/n0[5]
    hh20[1,6]<-0;hh20[1,7]<-0;
    hh20[2,2]<-sigma[4]/n0[1]+4*sigma[5]/n0[2]+sigma[6]/n0[3]
    hh20[2,3]<-sigma[4]/n0[3]
    hh20[2,4]<-0
    hh20[2,5]<--2*sigma[4]/n0[2]
    hh20[2,6]<-0
    hh20[2,7]<-sigma[4]/n0[1]-sigma[4]/n0[3]
    hh20[3,3]<-sigma[4]/n0[3]+4*sigma[4]/n0[6]+sigma[4]/n0[9]
    hh20[3,4]<--2*sigma[4]/n0[6]
    hh20[3,5]<-0
    hh20[3,6]<-sigma[4]/n0[9]
    hh20[3,7]<--sigma[4]/n0[3]
    hh20[4,4]<-sigma[4]/n0[4]+4*sigma[4]/n0[5]+sigma[4]/n0[6]
    hh20[4,5]<-4*sigma[4]/n0[5]
    hh20[4,6]<-0
    hh20[4,7]<-0
    hh20[5,5]<-sigma[4]/n0[2]+4*sigma[4]/n0[5]+sigma[4]/n0[8]
    hh20[5,6]<--2*sigma[4]/n0[8]
    hh20[5,7]<-2*sigma[4]/n0[8]
    hh20[6,6]<-sigma[4]/n0[7]+4*sigma[4]/n0[8]+sigma[4]/n0[9]
    hh20[6,7]<--2*sigma[4]/n0[7]-4*sigma[4]/n0[8]
    hh20[7,7]<-sigma[4]/n0[1]+sigma[4]/n0[3]+4*sigma[4]/n0[7]+4*sigma[4]/n0[8]
    for(i in 2:7)
    {
      for(j in 1:(i-1))
      {
        hh20[i,j]<-hh20[j,i]
      }
    }
    ##########################################################################
    b_line20[1]<-sumx[1]/m_sam[1]+2*sumx[2]/m_sam[2]+sumx[3]/m_sam[3]-4*sumwx[5]/n0[5]
    b_line20[2]<-sumwx[1]/n0[1]-2*sumwx[2]/n0[2]+sumwx[3]/n0[3]
    b_line20[3]<-sumwx[3]/n0[3]-2*sumwx[6]/n0[6]+sumwx[9]/n0[9]
    b_line20[4]<-sumwx[4]/n0[4]-2*sumwx[5]/n0[5]+sumwx[6]/n0[6]
    b_line20[5]<-sumwx[2]/n0[2]-2*sumwx[5]/n0[5]+sumwx[8]/n0[8]
    b_line20[6]<-sumwx[7]/n0[7]-2*sumwx[8]/n0[8]+sumwx[9]/n0[9]
    b_line20[7]<-sumwx[1]/n0[1]-sumwx[3]/n0[3]-2*sumwx[7]/n0[7]+2*sumwx[8]/n0[8]
    B20<-solve(hh20,b_line20)
    ##################obtain mean######################
    m[1]<-(sumx[1]-sigma[1]*B20[1])/m_sam[1]
    m[2]<-(sumx[2]-sigma[1]*B20[1]*2)/m_sam[2]
    m[3]<-(sumx[3]-sigma[1]*B20[1])/m_sam[3]
    m[4]<-(sumwx[1]-sigma[4]*(B20[2]+B20[7]))/n0[1]
    m[5]<-(sumwx[2]+sigma[4]*(2*B20[2]-B20[5]))/n0[2]
    m[6]<-(sumwx[3]-sigma[4]*(B20[2]+B20[3]-B20[7]))/n0[3]
    m[7]<-(sumwx[4]-sigma[4]*B20[4])/n0[4]
    m[8]<-(sumwx[5]+sigma[4]*(4*B20[1]+2.0*B20[4]+2*B20[5]))/n0[5]
    m[9]<-(sumwx[6]+sigma[4]*(2*B20[3]-B20[4]))/n0[6]
    m[10]<-(sumwx[7]-sigma[4]*(B20[6]-2*B20[7]))/n0[7]
    m[11]<-(sumwx[8]-sigma[4]*(B20[5]-2*B20[6]+2*B20[7]))/n0[8]
    m[12]<-(sumwx[9]-sigma[4]*(B20[3]+B20[6]))/n0[9]
    ######obtain variance##############
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d3) {swx[i] <- WW[i,]%*%(dataF2-m[i+3])^2 }
    s0[1]<-sum(swx);sigma[4]<-s0[1]/m_sam[4]
    
    sigma00<-sigma[4]-sigma[1]
    if (sigma00<0){sigma00<-0;sigma[4]<-sigma[1]}
    sigma[c(5:12)]<-sigma[4]
    aa1<-ss1+ss2+ss3;aa2<-m_sam[1]+m_sam[2]+m_sam[3]
    n_iter<-0
    aa0<-sigma[1]
    aa4<-1000
    while(aa4>0.0001){
      n_iter<-n_iter+1
      aa3<-sigma[1]/(sigma[1]+sigma00)
      sigma[1]<-(aa1+aa3^2*s0[1])/(aa2+aa3*m_sam[4])
      aa4<-abs(sigma[1]-aa0)
      aa0<-sigma[1]
      if (n_iter>20) break 
    }
    sigma[2]<-sigma[3]<-sigma[1]
    sigma[c(4:12)]<-sigma[1]+sigma00
    #######criteria for iterations to stop####################
    pi<-m[c(4:12)];gh<-sigma[c(4:12)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mix_pi))) 
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*3
  #########first order genetic parameter process##########
  hh201 <- matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,0,-1,1,1,1,0,0,0,-1,-1,-1,1,0,-1,1,0,-1,1,0,-1,1,0,-1,1,
                    0,-1,0,0,0,0,0,0,0,0,0,0,1,0,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5),12,5)
  B201<-solve(crossprod(hh201,hh201))%*%crossprod(hh201,m)
  #################second oder genetic parameter process#############
  jj<-sigmaF2 - sigma[4]
  gg<-sigma[4]-sigma[1]
  if(jj<0){jj<-0}
  if(gg<0){gg<-0}
  ll<-jj/sigmaF2
  rr<-gg/sigmaF2   
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
  
  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2); 
  F2w1<-1/(12*m_sam[4])
  F2bmw <- matrix(0,m_sam[4],1); F2bmwsl <- matrix(0,m_sam[4],d3)
  for(i in 1:d3){
    F2gg <- (dataF2 - pi[i])/sqrt(gh[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < m_sam[4]){F2P2 <- F2P2+runif(m_sam[4])/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*m_sam[4]) + sum((F2P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F2u <- as.matrix(c(12*m_sam[4]*((F2dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F2dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F2dd[3]/m_sam[4]-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),(1-pkolm(F2D,m_sam[4]))))
  
  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("MX2-A-AD",round(abc,4),round(AIC,4),round(t(m),4),round(sigma[1],4),round(sigma[4],4),
                       round(t(mix_pi),4),round(B201[1],4)," "," "," ",round(B201[2],4),round(B201[3],4)," "," ", " "," "," "," ",round(B201[4],4),round(B201[5],4),round(jj,4),round(ll*100,4),round(gg,4),round(rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}

#####################MX2-EA-AD(E-4)#####################
G4F2ModelFun[[22]] <- function(K1,logL,df11,df21,df31,df41){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF2 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2)))
  m_esp<-0.0001
  sigmaF2<-var(dataF2);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a2<-sqrt(sigmaF2/m_sam[4]*(m_sam[4]-1))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2)))
  ###############procedure start##############################    
  d4<-5
  mi<-as.matrix(c(0.0625,0.25,0.375,0.25,0.0625))
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaF2/2,sigmaF2/2,sigmaF2/2,sigmaF2/2,sigmaF2/2))
  m[c(1:3)]<-m[c(1:3)]
  if (m[1]<m[3]) {a2<--a2}
  m[5]<-m[4]+1.5*a2;m[6]<-m[4];m[7]<-m[4]-1.5*a2
  m[8]<-m[4]-3*a2;m[4]<-m[4]+3*a2
  pi<-m[c(4:8)];gh<-sigma[c(4:8)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mi)))
  ########iteration process############
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,5,m_sam[4]); swx <- matrix(0,5,1); s0<-matrix(0,1,1)
  hh21<-matrix(0,4,4);b_line21<-matrix(0,4,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:8)]
    ssigma<-sigma[c(4:8)]
    for(i in 1:d4) { WW[i,] <- mi[i]*dnorm(dataF2,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF2,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF2
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001
    ####################################################
    hh21[1,1]<-sigma[1]/m_sam[1]+4*sigma[2]/m_sam[2]+sigma[3]/m_sam[3]+16*sigma[4]/n0[3]
    hh21[1,2]<--4*sigma[4]/n0[3]
    hh21[1,3]<-8*sigma[4]/n0[3]
    hh21[1,4]<--4*sigma[4]/n0[3]
    hh21[2,2]<-sigma[4]/n0[1]+4*sigma[4]/n0[2]+sigma[4]/n0[3]
    hh21[2,3]<--2*sigma[4]/n0[2]-2*sigma[4]/n0[3]
    hh21[2,4]<-sigma[4]/n0[3]
    hh21[3,3]<-sigma[4]/n0[2]+4*sigma[4]/n0[3]+sigma[4]/n0[4]
    hh21[3,4]<--2*sigma[4]/n0[3]-2*sigma[4]/n0[4]
    hh21[4,4]<-sigma[4]/n0[3]+4*sigma[4]/n0[4]+sigma[4]/n0[5]
    for(i in 2:4)
    {
      for(j in 1:(i-1))
      {
        hh21[i,j]<-hh21[j,i]
      }
    }
    #######################################################
    b_line21[1]<-sumx[1]/m_sam[1]+2*sumx[2]/m_sam[2]+sumx[3]/m_sam[3]-4*sumwx[3]/n0[3]
    b_line21[2]<-sumwx[1]/n0[1]-2*sumwx[2]/n0[2]+sumwx[3]/n0[3]
    b_line21[3]<-sumwx[2]/n0[2]-2*sumwx[3]/n0[3]+sumwx[4]/n0[4]
    b_line21[4]<-sumwx[3]/n0[3]-2*sumwx[4]/n0[4]+sumwx[5]/n0[5]
    B21<-solve(hh21,b_line21)
    #############obtain mean#########################
    m[1]<-(sumx[1]-sigma[1]*B21[1])/m_sam[1]
    m[2]<-(sumx[2]-sigma[2]*B21[1]*2)/m_sam[2]
    m[3]<-(sumx[3]-sigma[3]*B21[1])/m_sam[3]
    m[4]<-(sumwx[1]-sigma[4]*B21[2])/n0[1]
    m[5]<-(sumwx[2]+sigma[5]*(2*B21[2]-B21[3]))/n0[2]
    m[6]<-(sumwx[3]+sigma[6]*(4*B21[1]-B21[2]+2*B21[3]-B21[4]))/n0[3]
    m[7]<-(sumwx[4]+sigma[7]*(-B21[3]+2*B21[4]))/n0[4]
    m[8]<-(sumwx[5]-sigma[8]*B21[4])/n0[5]
    ######obtain variance##############
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d4) {swx[i] <- WW[i,]%*%(dataF2-m[i+3])^2 }
    s0[1]<-sum(swx)
    sigma[4]<-s0[1]/m_sam[4]
    
    sigma00<-sigma[4]-sigma[1]
    if (sigma00<0){sigma00<-0;sigma[4]<-sigma[1]}
    sigma[c(5:8)]<-sigma[4]
    aa1<-ss1+ss2+ss3;aa2<-m_sam[1]+m_sam[2]+m_sam[3]
    n_iter<-0
    aa0<-sigma[1]
    aa4<-1000
    while(aa4>0.0001){
      n_iter<-n_iter+1
      aa3<-sigma[1]/(sigma[1]+sigma00)
      sigma[1]<-(aa1+aa3^2*s0[1])/(aa2+aa3*m_sam[4])
      aa4<-abs(sigma[1]-aa0)
      aa0<-sigma[1]
      if (n_iter>20) break 
    }
    sigma[2]<-sigma[3]<-sigma[1]
    sigma[c(4:8)]<-sigma[1]+sigma00
    ########criteria for iterations to stop####################
    pi<-m[c(4:8)];gh<-sigma[c(4:8)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mix_pi))) 
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*2
  #########first order genetic parameter process##########  
  hh211<- matrix(c(1,1,1,1,1,1,1,1,2,0,-2,2,1,0,-1,-2,1,0,-1,0,0,0,0,0,0,1,0,0.5,0.5,0.5,0.5,0.5),8,4)
  B211<-solve(crossprod(hh211,hh211))%*%crossprod(hh211,m)
  #################second oder genetic parameter process#############
  jj<-sigmaF2 - sigma[4]
  gg<-sigma[4]-sigma[1]
  if(jj<0){jj<-0}
  if(gg<0){gg<-0}
  ll<-jj/sigmaF2
  rr<-gg/sigmaF2     
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
  
  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2); 
  F2w1<-1/(12*m_sam[4])
  F2bmw <- matrix(0,m_sam[4],1); F2bmwsl <- matrix(0,m_sam[4],d4)
  for(i in 1:d4){
    F2gg <- (dataF2 - pi[i])/sqrt(gh[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < m_sam[4]){F2P2 <- F2P2+runif(m_sam[4])/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*m_sam[4]) + sum((F2P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F2u <- as.matrix(c(12*m_sam[4]*((F2dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F2dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F2dd[3]/m_sam[4]-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),(1-pkolm(F2D,m_sam[4]))))
  
  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("MX2-EA-AD",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," ",round(sigma[1],4),round(sigma[4],4),
                       round(t(mix_pi),4)," "," "," "," ",round(B211[1],4)," "," "," ",round(B211[2],4)," "," "," "," "," "," "," ",round(B211[3],4),round(B211[4],4),round(jj,4),round(ll*100,4),round(gg,4),round(rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}

########################MX2-CD-AD(E-5)#############################
G4F2ModelFun[[23]] <- function(K1,logL,df11,df21,df31,df41){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF2 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2)))
  m_esp<-0.0001
  sigmaF2<-var(dataF2);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a2<-sqrt(sigmaF2/m_sam[4]*(m_sam[4]-1))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2)))
  ###############procedure start##############################   
  d5<-4
  mi<-as.matrix(c(0.5625,0.1875,0.1875,0.0625))
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaF2/2,sigmaF2/2,sigmaF2/2,sigmaF2/2))
  m[c(1:3)]<-m[c(1:3)]
  if (m[1]<m[3]) {a2<--a2}
  m[5]<-m[4]+0.5*a2;m[6]<-m[4]-0.5*a2 ;m[7]<-m[4]-3*a2;m[4]<-m[4]+3*a2
  pi<-m[c(4:7)];gh<-sigma[c(4:7)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mi))) 
  ########iteration process############
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,4,m_sam[4]); swx <- matrix(0,4,1);rr7<-matrix(0,2,1);s0<-matrix(0,1,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:7)];ssigma<-sigma[c(4:7)]
    for(i in 1:d5) { WW[i,] <- mi[i]*dnorm(dataF2,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF2,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF2
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001  
    aa1<-sumx[1]/m_sam[1]+2*sumx[2]/m_sam[2]+sumx[3]/m_sam[3]-3*sumwx[1]/n0[1]-sumwx[4]/n0[4]
    aa2<-sumwx[1]/n0[1]-sumwx[2]/n0[2]-sumwx[3]/n0[3]+sumwx[4]/n0[4]
    aa3<-sigma[1]/m_sam[1]+4*sigma[2]/m_sam[2]+sigma[3]/m_sam[3]+9*sigma[4]/n0[1]+sigma[7]/n0[4]
    aa4<--3*sigma[4]/n0[1]-sigma[7]/n0[4]
    aa5<-sigma[4]/n0[1]+sigma[5]/n0[2]+sigma[6]/n0[3]+sigma[7]/n0[4]
    aa6<-aa3*aa5-aa4^2;aa7<-aa1*aa5-aa2*aa4;aa8<-aa2*aa3-aa1*aa4
    rr7[1]<-aa7/aa6 ;rr7[2]<-aa8/aa6
    ############obtain mean###########################
    m[1]<-(sumx[1]-sigma[1]*rr7[1])/m_sam[1]
    m[2]<-(sumx[2]-sigma[2]*rr7[1]*2)/m_sam[2]
    m[3]<-(sumx[3]-sigma[3]*rr7[1])/m_sam[3]
    m[4]<-(sumwx[1]+sigma[4]*(3*rr7[1]+rr7[2]))/n0[1]
    m[5]<-(sumwx[2]+sigma[5]*rr7[2])/n0[2]
    m[6]<-(sumwx[3]+sigma[6]*rr7[2])/n0[3]
    m[7]<-(sumwx[4]+sigma[7]*(rr7[1]-rr7[2]))/n0[4]
    ######obtain variance##############
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d5) {swx[i] <- WW[i,]%*%(dataF2-m[i+3])^2 }
    s0[1]<-sum(swx)
    sigma[4]<-s0[1]/m_sam[4]
    
    sigma00<-sigma[4]-sigma[1]
    if (sigma00<0){sigma00<-0;sigma[4]<-sigma[1]}
    sigma[c(5:7)]<-sigma[4]
    aa1<-ss1+ss2+ss3;aa2<-m_sam[1]+m_sam[2]+m_sam[3]
    n_iter<-0
    aa0<-sigma[1]
    aa4<-1000
    while(aa4>0.0001){
      n_iter<-n_iter+1
      aa3<-sigma[1]/(sigma[1]+sigma00)
      sigma[1]<-(aa1+aa3^2*s0[1])/(aa2+aa3*m_sam[4])
      aa4<-abs(sigma[1]-aa0)
      aa0<-sigma[1]
      if (n_iter>20) break 
    }
    sigma[2]<-sigma[3]<-sigma[1]
    sigma[c(4:7)]<-sigma[1]+sigma00
    ########criteria for iterations to stop####################
    pi<-m[c(4:7)];gh<-sigma[c(4:7)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mix_pi))) 
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*3 
  #########first order genetic parameter process##########  
  hh22<- matrix(c(1,1,1,1,1,1,1,1,1,-1,1,1,-1,-1,1,1,-1,1,-1,1,-1,1,
                  0,-1,0,0,0,0,0,1,0,0.5,0.5,0.5,0.5),7,5)
  B22<-solve(crossprod(hh22,hh22))%*%crossprod(hh22,m)
  #################second oder genetic parameter process#############
  jj<-sigmaF2 - sigma[4]
  gg<-sigma[4]-sigma[1]
  if(jj<0){jj<-0}
  if(gg<0){gg<-0}
  ll<-jj/sigmaF2
  rr<-gg/sigmaF2     
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
  
  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2); 
  F2w1<-1/(12*m_sam[4])
  F2bmw <- matrix(0,m_sam[4],1); F2bmwsl <- matrix(0,m_sam[4],d5)
  for(i in 1:d5){
    F2gg <- (dataF2 - pi[i])/sqrt(gh[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < m_sam[4]){F2P2 <- F2P2+runif(m_sam[4])/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*m_sam[4]) + sum((F2P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F2u <- as.matrix(c(12*m_sam[4]*((F2dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F2dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F2dd[3]/m_sam[4]-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),(1-pkolm(F2D,m_sam[4]))))
  
  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("MX2-CD-AD",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," ",round(sigma[1],4),round(sigma[4],4),
                       round(t(mix_pi),4)," "," "," "," "," ",round(B22[1],4)," "," "," ",round(B22[2],4),round(B22[3],4)," "," ", " "," "," "," ",round(B22[4],4),round(B22[5],4),round(jj,4),round(ll*100,4),round(gg,4),round(rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}

####################MX2-EAD_AD(E-6)###########################
G4F2ModelFun[[24]] <- function(K1,logL,df11,df21,df31,df41){
  dataP1 <- as.matrix(as.numeric(df11[,1]))
  dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]))
  dataF2 <- as.matrix(as.numeric(df41[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1])
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2)))
  m_esp<-0.0001
  sigmaF2<-var(dataF2);sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1)
  sigma0<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  a2<-sqrt(sigmaF2/m_sam[4]*(m_sam[4]-1))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2)))
  ###############procedure start##############################   
  d1<-3
  mi<-as.matrix(c(0.5625,0.375,0.0625))
  sigma<-as.matrix(c(sigma0,sigma0,sigma0,sigmaF2/2,sigmaF2/2,sigmaF2/2))
  m[c(1:3)]<-m[c(1:3)]
  if(m[1]<m[3]){a2<--a2}
  m[5]<-m[4];m[6]<-m[4]-3*a2;m[4]<-m[4]+3*a2
  pi<-m[c(4:6)];gh<-sigma[c(4:6)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mi))) 
  ########iteration process############
  iteration <- 0; stopa <- 1000
  WW <- matrix(0,3,m_sam[4]); swx <- matrix(0,3,1);rr8<-matrix(0,2,1);s0<-matrix(0,1,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    ############E-step#############
    mm<-m[c(4:6)];ssigma<-sigma[c(4:6)]
    for(i in 1:d1) { WW[i,] <- mi[i]*dnorm(dataF2,m[i+3],sqrt(sigma[i+3]))/dmixnorm(dataF2,mm,sqrt(ssigma),mi) }
    mix_pi <- as.matrix(rowSums(WW)/m_sam[4])
    sumwx <- WW%*%dataF2
    n0 <- m_sam[4]*mix_pi
    n0[abs(n0)<0.000001] <- 0.000001   
    ##########obtain mean###########################
    aa1<-sumx[1]/m_sam[1]+2*sumx[2]/m_sam[2]+sumx[3]/m_sam[3]-3*sumwx[1]/n0[1]-sumwx[3]/n0[3]
    aa2<-sumwx[1]/n0[1]-2*sumwx[2]/n0[2]+sumwx[3]/n0[3]
    aa3<-sigma[1]/m_sam[1]+4*sigma[2]/m_sam[2]+sigma[3]/m_sam[3]+9*sigma[4]/n0[1]+sigma[6]/n0[3]
    aa4<--3*sigma[4]/n0[1]-sigma[6]/n0[3]
    aa5<-sigma[4]/n0[1]+4*sigma[5]/n0[2]+sigma[6]/n0[3]
    aa6<-aa3*aa5-aa4^2;aa7<-aa1*aa5-aa2*aa4;aa8<-aa2*aa3-aa1*aa4
    rr8[1]<-aa7/aa6;rr8[2]<-aa8/aa6
    m[1]<-(sumx[1]-sigma[1]*rr8[1])/m_sam[1]
    m[2]<-(sumx[2]-sigma[2]*rr8[1]*2)/m_sam[2]
    m[3]<-(sumx[3]-sigma[3]*rr8[1])/m_sam[3]
    m[4]<-(sumwx[1]+sigma[4]*(3*rr8[1]-rr8[2]))/n0[1]
    m[5]<-(sumwx[2]+sigma[5]*2*rr8[2])/n0[2]
    m[6]<-(sumwx[3]+sigma[6]*(rr8[1]-rr8[2]))/n0[3]
    #####obtain variance##############
    ss1<-sum((dataP1-m[1])^2) ;ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d1) {swx[i] <- WW[i,]%*%(dataF2-m[i+3])^2 }
    s0[1]<-sum(swx);sigma[4]<-s0[1]/m_sam[4]
    
    sigma00<-sigma[4]-sigma[1]
    if (sigma00<0){sigma00<-0;sigma[4]<-sigma[1]}
    sigma[c(5,6)]<-sigma[4]
    aa1<-ss1+ss2+ss3;aa2<-m_sam[1]+m_sam[2]+m_sam[3]
    n_iter<-0
    aa0<-sigma[1]
    aa4<-1000
    while(aa4>0.0001){
      n_iter<-n_iter+1
      aa3<-sigma[1]/(sigma[1]+sigma00)
      sigma[1]<-(aa1+aa3*aa3*s0[1])/(aa2+aa3*m_sam[4])
      aa4<-abs(sigma[1]-aa0)
      aa0<-sigma[1]
      if (n_iter>20) break 
    }
    sigma[2]<-sigma[3]<-sigma[1]
    sigma[c(4:6)]<-sigma[1]+sigma00
    ########criteria for iterations to stop####################
    pi<-m[c(4:6)];gh<-sigma[c(4:6)]
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma[1]))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma[2]))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma[3]))))+sum(log(dmixnorm(dataF2,pi,sqrt(gh),mix_pi))) 
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*2
  #########first order genetic parameter process##########  
  hh23<- matrix(c(1,1,1,1,1,1,2,2,-2,2,0,-2,1,0,-1,0,0,0,0,1,0,0.5,0.5,0.5),6,4)
  B23<-solve(crossprod(hh23,hh23))%*%crossprod(hh23,m)
  #################second oder genetic parameter process#############
  jj<-sigmaF2 - sigma[4]
  gg<-sigma[4]-sigma[1]
  if(jj<0){jj<-0}
  if(gg<0){gg<-0}
  ll<-jj/sigmaF2
  rr<-gg/sigmaF2       
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
  
  ####################################hypothesis testing for F2#########################################
  dataF2 <- sort(dataF2); 
  F2w1<-1/(12*m_sam[4])
  F2bmw <- matrix(0,m_sam[4],1); F2bmwsl <- matrix(0,m_sam[4],d1)
  for(i in 1:d1){
    F2gg <- (dataF2 - pi[i])/sqrt(gh[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mix_pi[i]
  }
  F2P2 <- rowSums(F2bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F2P2)))[1]
  if(nn < m_sam[4]){F2P2 <- F2P2+runif(m_sam[4])/1e4}
  ##########################################################
  F2dd <- as.matrix(c(sum(F2P2),sum(F2P2^2),sum((F2P2-0.5)^2)))
  F2WW2 <- 1/(12*m_sam[4]) + sum((F2P2 - (as.matrix(c(1:m_sam[4])) - 0.5)/m_sam[4])^2)
  F2u <- as.matrix(c(12*m_sam[4]*((F2dd[1]/m_sam[4]-0.5)^2),((45*m_sam[4])/4)*((F2dd[2]/m_sam[4]-1/3)^2),180*m_sam[4]*((F2dd[3]/m_sam[4]-1/12)^2)))
  F2D <- as.numeric(ks.test(F2P2,"punif")[[1]][1])
  F2tt <- as.matrix(c((1 - pchisq(F2u[1],1)),(1 - pchisq(F2u[2],1)),(1 - pchisq(F2u[3],1)),K1(F2WW2),(1-pkolm(F2D,m_sam[4]))))
  
  F2tt[which( F2tt>=10e-4)]<-round(F2tt[which(F2tt>=10e-4)],4);F2tt[which(F2tt<10e-4)]<-format(F2tt[which(F2tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("MX2-EAD-AD",round(abc,4),round(AIC,4),round(t(m),4)," "," "," "," "," "," ",round(sigma[1],4),round(sigma[4],4),
                       round(t(mix_pi),4)," "," "," "," "," "," ",round(B23[1],4)," "," "," ",round(B23[2],4)," "," "," "," "," "," "," ",round(B23[3],4),round(B23[4],4),round(jj,4),round(ll*100,4),round(gg,4),round(rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mi)
  return(OUTPUT)
}


K1G4F2 <- function(x){
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

logLG4F2 <- function(nm,nng,mi,mn,s,d1) { sum2 <- sum(log(dmixnorm(d1,mn,sqrt(s),mi)));return (sum2) } 



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
    G4F2ModelFun[[i]](K1G4F2,logLG4F2,df11,df21,df31,df41)[[1]]
  }
  stopCluster(cl)
  mi<-NULL
  
}else{
  
allresultq=switch(model,"1MG-AD" = G4F2ModelFun[[1]](K1G4F2,logLG4F2,df11,df21,df31,df41),"1MG-A"=G4F2ModelFun[[2]](K1G4F2,logLG4F2,df11,df21,df31,df41),"1MG-EAD"=G4F2ModelFun[[3]](K1G4F2,logLG4F2,df11,df21,df31,df41),"1MG-NCD"=G4F2ModelFun[[4]](K1G4F2,logLG4F2,df11,df21,df31,df41),"2MG-ADI"=G4F2ModelFun[[5]](K1G4F2,logLG4F2,df11,df21,df31,df41),
         "2MG-AD"=G4F2ModelFun[[6]](K1G4F2,logLG4F2,df11,df21,df31,df41),"2MG-A"=G4F2ModelFun[[7]](K1G4F2,logLG4F2,df11,df21,df31,df41),"2MG-EA"=G4F2ModelFun[[8]](K1G4F2,logLG4F2,df11,df21,df31,df41),"2MG-CD"=G4F2ModelFun[[9]](K1G4F2,logLG4F2,df11,df21,df31,df41),"2MG-EAD"=G4F2ModelFun[[10]](K1G4F2,logLG4F2,df11,df21,df31,df41),
         "PG-ADI"=G4F2ModelFun[[11]](K1G4F2,logLG4F2,df11,df21,df31,df41),"PG-AD"=G4F2ModelFun[[12]](K1G4F2,logLG4F2,df11,df21,df31,df41),"MX1-AD-ADI"=G4F2ModelFun[[13]](K1G4F2,logLG4F2,df11,df21,df31,df41),"MX1-AD-AD"=G4F2ModelFun[[14]](K1G4F2,logLG4F2,df11,df21,df31,df41),"MX1-A-AD"=G4F2ModelFun[[15]](K1G4F2,logLG4F2,df11,df21,df31,df41),
         "MX1-EAD-AD"=G4F2ModelFun[[16]](K1G4F2,logLG4F2,df11,df21,df31,df41),"MX1-NCD-AD"=G4F2ModelFun[[17]](K1G4F2,logLG4F2,df11,df21,df31,df41),"MX2-ADI-ADI"=G4F2ModelFun[[18]](K1G4F2,logLG4F2,df11,df21,df31,df41),"MX2-ADI-AD"=G4F2ModelFun[[19]](K1G4F2,logLG4F2,df11,df21,df31,df41),"MX2-AD-AD"=G4F2ModelFun[[20]](K1G4F2,logLG4F2,df11,df21,df31,df41),
         "MX2-A-AD"=G4F2ModelFun[[21]](K1G4F2,logLG4F2,df11,df21,df31,df41),"MX2-EA-AD"=G4F2ModelFun[[22]](K1G4F2,logLG4F2,df11,df21,df31,df41),"MX2-CD-AD"=G4F2ModelFun[[23]](K1G4F2,logLG4F2,df11,df21,df31,df41),"MX2-EAD-AD"=G4F2ModelFun[[24]](K1G4F2,logLG4F2,df11,df21,df31,df41))
  
allresult<-allresultq[[1]]
if(model=="PG-ADI"||model=="PG-AD"){
  mi<-NULL  
}else{
  mi<-allresultq[[2]]
} 
}
colnames(allresult) <- G4F2colname
out<-list(allresult,mi)
return(out)
}




