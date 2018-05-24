G5Fun<-function(df,model,G5text2){

data<-sapply(df,as.character)

dP1<-data[-1,which(data[1,]=="P1")];P1<-as.numeric(dP1[which(is.na(as.numeric(dP1))==FALSE)]);df11<-as.data.frame(P1)
dF1<-data[-1,which(data[1,]=="F1")];F1<-as.numeric(dF1[which(is.na(as.numeric(dF1))==FALSE)]);df21<-as.data.frame(F1)
dP2<-data[-1,which(data[1,]=="P2")];P2<-as.numeric(dP2[which(is.na(as.numeric(dP2))==FALSE)]);df31<-as.data.frame(P2)
dF2<-data[-1,which(data[1,]=="F2")];F2<-as.numeric(dF2[which(is.na(as.numeric(dF2))==FALSE)]);df41<-as.data.frame(F2)
dF23<-data[-1,which(data[1,]=="F23")];F23<-as.numeric(dF23[which(is.na(as.numeric(dF23))==FALSE)]);df51<-as.data.frame(F23)

G5colname<-c("Model","Log_Max_likelihood_Value","AIC","meanP1","meanF1","meanP2","meanF2[1]","meanF2[2]","meanF2[3]","meanF2[4]","meanF2[5]","meanF2[6]","meanF2[7]","meanF2[8]","meanF2[9]","VarF2(Component)",	
             "Proportion(F2)[1]","Proportion(F2)[2]","Proportion(F2)[3]","Proportion(F2)[4]","Proportion(F2)[5]","Proportion(F2)[6]","Proportion(F2)[7]","Proportion(F2)[8]","Proportion(F2)[9]",	
             "meanF2:3[1]","meanF2:3[2]","meanF2:3[3]","meanF2:3[4]","meanF2:3[5]","meanF2:3[6]","meanF2:3[7]","meanF2:3[8]","meanF2:3[9]","VarF2:3[1]","VarF2:3[2]","VarF2:3[3]","VarF2:3[4]","VarF2:3[5]","VarF2:3[6]","VarF2:3[7]","VarF2:3[8]","VarF2:3[9]",
             "Proportion(F2:3)[1]","Proportion(F2:3)[2]","Proportion(F2:3)[3]","Proportion(F2:3)[4]","Proportion(F2:3)[5]","Proportion(F2:3)[6]","Proportion(F2:3)[7]","Proportion(F2:3)[8]","Proportion(F2:3)[9]","Var(Residual)",
             "m1(m)","m2","m3","m4","m5","da(d)","db","ha(h)","hb","i","jab","jba","l","[d]","[h]","Major-Gene Var(F2)","Heritability(Major-Gene(F2))(%)","Poly-Gene Var(F2)","Heritability(Poly-Gene(F2))(%)","Major-Gene Var(F2:3)","Heritability(Major-Gene(F2:3))(%)","Poly-Gene Var(F2:3)",	"Heritability(Poly-Gene(F2:3))(%)",	
             "U1 square(P1)","P(U1 square(P1))","U2 square(P1)","P(U2 square(P1))","U3 square(P1)","P(U3 square(P1))","nW square(P1)","P(nW square(P1))","Dn(P1)","P(Dn(P1))","U1 square(F1)","P(U1 square(F1))","U2 square(F1)","P(U2 square(F1))","U3 square(F1)","P(U3 square(F1))","nW square(F1)","P(nW square(F1))","Dn(F1)","P(Dn(F1))",	
             "U1 square(P2)","P(U1 square(P2))","U2 square(P2)","P(U2 square(P2))","U3 square(P2)","P(U3 square(P2))","nW square(P2)","P(nW square(P2))","Dn(P2)","P(Dn(P2))","U1 square(F2)","P(U1 square(F2))","U2 square(F2)","P(U2 square(F2))","U3 square(F2)","P(U3 square(F2))","nW square(F2)","P(nW square(F2))","Dn(F2)","P(Dn(F2))",	
             "U1 square(F2:3)","P(U1 square(F2:3))","U2 square(F2:3)","P(U2 square(F2:3))","U3 square(F2:3)","P(U3 square(F2:3))","nW square(F2:3)","P(nW square(F2:3))","Dn(F2:3)","P(Dn(F2:3))")  


G5ModelFun<-list(NA)
###################define each model function############################
##########################1MG_AD(A_1)#################################
G5ModelFun[[1]] <- function(K1,logL,df11,df21,df31,df41,df51,G5text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataF2 <- as.matrix(as.numeric(df41[,1]))
  dataF3 <- as.matrix(as.numeric(df51[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2),sum(dataF3)))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2),mean(dataF3)))
  sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1);sigmaF2<-var(dataF2);sigmaF3<-var(dataF3)
  sumx6<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  sigma<-sumx6;sigma0<-sigma;sigma40<-sigmaF2;sigma400<-sigma40;sigma50<-sigmaF3;sigma500<-sigma50
  m_nf <- as.numeric(G5text2)
  ###############procedure start###########################
  d1<-3
  mix_pi4<-as.matrix(c(0.25,0.5,0.25));mix_pi5<-as.matrix(c(0.25,0.5,0.25))
  sigma4<-matrix(0,3,1);sigma5<-matrix(0,3,1)
  sigma4[c(1:3)]<-sigma;sigma5[c(1,3)]<-sigma/m_nf
  m4<-m[c(1:3)]
  m5<-as.matrix(c(m[1],0.25*m[1]+0.5*m[2]+0.25*m[3],m[3]))
  ########first order genetic parameters###############
  hh1<-matrix(c(1,1,1,1,1,0,-1,0,0,1,0,0.5),4,3)
  mm<-as.matrix(c(m[1:3],m5[2]))
  B1<-solve(crossprod(hh1,hh1))%*%crossprod(hh1,mm)
  a1<-B1[2];a2<-B1[3]
  sigma5[2]<-sigma5[1]+(0.5*a1^2+0.25*a2^2)/m_nf
  mi4<-mix_pi4[c(1:3)];mi5<-mix_pi5[c(1:3)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
    sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
  ##########iteration process###########
  ############E-step###################
  ############F2 process###############
  iteration <- 0; stopa <- 1000
  W4 <- matrix(0,3,m_sam[4]); swx24 <- matrix(0,3,1)
  W5 <- matrix(0,3,m_sam[5]); swx25 <- matrix(0,3,1)
  s0<-matrix(0,5,1);n0<-matrix(0,6,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:d1) { W4[i,] <- mix_pi4[i]*dnorm(dataF2,m4[i],sqrt(sigma4[i]))/dmixnorm(dataF2,m4,sqrt(sigma4),mix_pi4) }
    mi4 <- as.matrix(rowSums(W4)/m_sam[4])
    sumwx4 <- W4%*%dataF2
    #############F3 process############
    for(i in 1:d1) { W5[i,] <- mix_pi5[i]*dnorm(dataF3,m5[i],sqrt(sigma5[i]))/dmixnorm(dataF3,m5,sqrt(sigma5),mix_pi5) }
    mi5 <- as.matrix(rowSums(W5)/m_sam[5])
    sumwx5 <- W5%*%dataF3
    aaa0<-0
    s0[1]<-sumx[1]+sumwx4[1]+m_nf*sumwx5[1];s0[2]<-sumx[2]+sumwx4[2]
    s0[3]<-sumx[3]+sumwx4[3]+m_nf*sumwx5[3];s0[4]<-sumwx5[2]
    n0[1]<-m_sam[1]+mi4[1]*m_sam[4]+mi5[1]*m_sam[5]*m_nf
    n0[2]<-m_sam[2]+mi4[2]*m_sam[4]
    n0[3]<-m_sam[3]+mi4[3]*m_sam[4]+mi5[3]*m_sam[5]*m_nf
    n0[4]<-mi5[2]*m_sam[5]
    n0[c(1:4)][abs(n0[c(1:4)])<0.000001] <- 0.000001
    aa3<-s0[1]/n0[1]+s0[3]/n0[3]+2*s0[2]/n0[2]-4*s0[4]/n0[4]
    aa4<-sigma*(1/n0[1]+1/n0[3]+4/n0[2])   
    aa1<-1000
    n_iter<-0
    while(aa1>0.0001){
      n_iter<-n_iter+1 
      ########first order genetic parameters###############
      mm<-as.matrix(c(m[1:3],m5[2]))
      B11<-solve(crossprod(hh1,hh1))%*%crossprod(hh1,mm)
      aa1<-B11[2]
      aa2<-B11[3]
      sigma5[2]<-(sigma+0.5*aa1^2+0.25*aa2^2)/m_nf
      aa2<-aa4+16*sigma5[2]/n0[4]
      aaa1<-aa3/aa2
      m[1]<-(s0[1]-aaa1*sigma)/n0[1]
      m[2]<-(s0[2]-2*aaa1*sigma)/n0[2]
      m[3]<-(s0[3]-aaa1*sigma)/n0[3]
      m5[2]<-(s0[4]+4*aaa1*sigma5[2])/n0[4]
      aa1<-abs(aaa1-aaa0)
      aaa0<-aaa1
      if(n_iter>20) break 
    }
    m4[c(1:3)]<-m[c(1:3)];m5[c(1,3)]<-m[c(1,3)] 
    ########obtain variance#######################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d1) {swx24[i] <- W4[i,]%*%(dataF2-m4[i])^2 } ;for(i in 1:d1) {swx25[i] <- W5[i,]%*%(dataF3-m5[i])^2 }      
    s0[5]<-ss1+ss2+ss3+swx24[1]+swx24[2]+swx24[3]+(swx25[1]+swx25[3])*m_nf
    n0[5]<-m_sam[1]+m_sam[2]+m_sam[3]+m_sam[4]+(mi5[1]+mi5[3])*m_sam[5]
    n0[6]<-mi5[2]*m_sam[5]
    aaa0<-sigma
    ########first order genetic parameters###############
    mm<-as.matrix(c(m[1:3],m5[2]))
    B111<-solve(crossprod(hh1,hh1))%*%crossprod(hh1,mm)
    aa1<-B111[2];aa2<-B111[3]
    aa1<-0.5*aa1^2+0.25*aa2^2
    n_iter<-0;aa3<-1000
    while(aa3>0.0001){
      n_iter<-n_iter+1
      aa2<-sigma/(sigma+aa1)
      sigma<-(s0[5]+aa2^2*swx25[2])/(n0[5]+aa2*n0[6])
      aa3<-abs(sigma-aaa0)
      aaa0<-sigma
      if(n_iter>20) break 
    }
    if(sigma<sigma0){sigma<-sigma0}
    sigma4[c(1,2,3)]<-sigma
    sigma5[1]<-sigma/m_nf;sigma5[2]<-(sigma+aa1)/m_nf;sigma5[3]<-sigma/m_nf
    ##############################criteria for iterations to stop######################## 
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
      sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*4
  mm<-as.matrix(c(m[1:3],m5[2]))
  B1111<-solve(crossprod(hh1,hh1))%*%crossprod(hh1,mm)
  ########second order genetic parameters###############   
  F2jj<-sigmaF2-sigma
  if(F2jj<0){F2jj<-0}
  F2ll<-F2jj/sigmaF2
  F3jj<-sigmaF3-sigma/m_nf
  if(F3jj<0) {F3jj<-0}
  F3ll<-F3jj/sigmaF3   
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma))
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
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma))
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
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma))
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
    F2gg <- (dataF2 - m4[i])/sqrt(sigma4[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mi4[i]
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
  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3); 
  F3w1<-1/(12*m_sam[5])
  F3bmw <- matrix(0,m_sam[5],1); F3bmwsl <- matrix(0,m_sam[5],d1)
  for(i in 1:d1){
    F3gg <- (dataF3 - m5[i])/sqrt(sigma5[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mi5[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[5]){F3P2 <- F3P2+runif(m_sam[5])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[5]) + sum((F3P2 - (as.matrix(c(1:m_sam[5])) - 0.5)/m_sam[5])^2)
  F3u <- as.matrix(c(12*m_sam[5]*((F3dd[1]/m_sam[5]-0.5)^2),((45*m_sam[5])/4)*((F3dd[2]/m_sam[5]-1/3)^2),180*m_sam[5]*((F3dd[3]/m_sam[5]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[5]))))
  
  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("1MG-AD",round(abc,4),round(AIC,4),round(m[1],4),round(m[2],4),round(m[3],4), round(t(m4),4)," "," "," "," "," "," ",round(sigma4[1],4),
                       round(t(mi4),4)," "," "," "," "," "," ",round(t(m5),4)," "," "," "," "," "," ",round(t(sigma5),4)," "," "," "," "," "," ",round(t(mi5),4)," "," "," "," "," "," ",round(sigma,4),          
                       round(B1111[1],4)," "," "," "," ",round(B1111[2],4)," ",round(B1111[3],4)," "," "," "," "," "," "," ",round(F2jj,4),round(F2ll*100,4)," "," ",round(F3jj,4),round(F3ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mix_pi4,mix_pi5)
  return(OUTPUT)
}
##################1MG-A(A-2)#########################
G5ModelFun[[2]] <- function(K1,logL,df11,df21,df31,df41,df51,G5text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataF2 <- as.matrix(as.numeric(df41[,1]))
  dataF3 <- as.matrix(as.numeric(df51[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2),sum(dataF3)))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2),mean(dataF3)))
  sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1);sigmaF2<-var(dataF2);sigmaF3<-var(dataF3)
  sumx6<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  sigma<-sumx6;sigma0<-sigma;sigma40<-sigmaF2;sigma400<-sigma40;sigma50<-sigmaF3;sigma500<-sigma50
  m_nf <- as.numeric(G5text2)
  ###############procedure start###########################
  d1<-3
  mix_pi4<-as.matrix(c(0.25,0.5,0.25));mix_pi5<-as.matrix(c(0.25,0.5,0.25))
  sigma4<-matrix(0,3,1);sigma5<-matrix(0,3,1)
  sigma4[c(1:3)]<-sigma;sigma5[c(1,3)]<-sigma/m_nf
  #######first order parameters############
  hh2<-matrix(c(1,1,1,1,0,-1),3,2)
  B2<-solve(crossprod(hh2,hh2))%*%crossprod(hh2,m[c(1:3)])
  a1<-B2[2]
  m4<-m[c(1:3)];m5<-m[c(1:3)]
  sigma5[2]<-sigma5[1]+(0.5*a1^2)/m_nf
  mi4<-mix_pi4[c(1:3)];mi5<-mix_pi5[c(1:3)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
    sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
  ##########iteration process###########
  ############E-step###################
  ############F2 process###############
  iteration <- 0; stopa <- 1000
  W4 <- matrix(0,3,m_sam[4]); swx24 <- matrix(0,3,1)
  W5 <- matrix(0,3,m_sam[5]); swx25 <- matrix(0,3,1)
  s0<-matrix(0,5,1);n0<-matrix(0,5,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:d1) { W4[i,] <- mix_pi4[i]*dnorm(dataF2,m4[i],sqrt(sigma4[i]))/dmixnorm(dataF2,m4,sqrt(sigma4),mix_pi4) }
    mi4 <- as.matrix(rowSums(W4)/m_sam[4])
    sumwx4 <- W4%*%dataF2
    #############F3 process############
    for(i in 1:d1) { W5[i,] <- mix_pi5[i]*dnorm(dataF3,m5[i],sqrt(sigma5[i]))/dmixnorm(dataF3,m5,sqrt(sigma5),mix_pi5) }
    mi5 <- as.matrix(rowSums(W5)/m_sam[5])
    sumwx5 <- W5%*%dataF3
    aaa0<-0
    s0[1]<-sumx[1]+sumwx4[1]+m_nf*sumwx5[1];s0[3]<-sumx[3]+sumwx4[3]+m_nf*sumwx5[3]
    n0[1]<-m_sam[1]+mi4[1]*m_sam[4]+mi5[1]*m_sam[5]*m_nf;n0[3]<-m_sam[3]+mi4[3]*m_sam[4]+mi5[3]*m_sam[5]*m_nf
    n0[c(1:3)][abs(n0[c(1:3)])<0.000001] <- 0.000001
    aa1<-s0[1]/n0[1]+s0[3]/n0[3];aa2<-sigma*(1/n0[1]+1/n0[3])   
    aa6<-1000;n_iter<-0
    while(aa6>0.0001){
      n_iter<-n_iter+1 
      ########first order genetic parameters###############
      B22<-solve(crossprod(hh2,hh2))%*%crossprod(hh2,m[c(1:3)])
      aa3<-B22[2]
      sigma5[2]<-(sigma+0.5*aa3^2)/m_nf
      s0[2]<-sumx[2]+sumwx4[2]+sigma*sumwx5[2]/sigma5[2]
      n0[2]<-m_sam[2]+mi4[2]*m_sam[4]+sigma*mi5[2]*m_sam[5]/sigma5[2]
      aa3<-aa1-2*s0[2]/n0[2]
      aa4<-aa2+4*sigma/n0[2]
      aaa1<-aa3/aa4
      m[1]<-(s0[1]-aaa1*sigma)/n0[1]
      m[2]<-(s0[2]+2*aaa1*sigma)/n0[2]
      m[3]<-(s0[3]-aaa1*sigma)/n0[3]
      aa6<-abs(aaa1-aaa0)
      aaa0<-aaa1
      if(n_iter>20) break 
    }
    m4[c(1:3)]<-m[c(1:3)];m5[c(1:3)]<-m[c(1:3)] 
    ##########obtain variance###############
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d1) {swx24[i] <- W4[i,]%*%(dataF2-m4[i])^2 } ;for(i in 1:d1) {swx25[i] <- W5[i,]%*%(dataF3-m5[i])^2 }      
    s0[5]<-ss1+ss2+ss3+swx24[1]+swx24[2]+swx24[3]+(swx25[1]+swx25[3])*m_nf
    n0[5]<-m_sam[1]+m_sam[2]+m_sam[3]+m_sam[4]+(mi5[1]+mi5[3])*m_sam[5]
    aaa0<-sigma
    ########first order genetic parameters##############
    B222<-solve(crossprod(hh2,hh2))%*%crossprod(hh2,m[c(1:3)])
    aa1<-B222[2];aa1<-0.5*aa1^2
    n_iter<-0;aa3<-1000
    while(aa3>0.0001){
      n_iter<-n_iter+1
      aa2<-sigma/(sigma+aa1)
      sigma<-(s0[5]+aa2^2*m_nf*swx25[2])/(n0[5]+aa2*mi5[2]*m_sam[5])
      aa3<-abs(sigma-aaa0)
      aaa0<-sigma
      if(n_iter>20) break 
    }
    
    if(sigma<sigma0){sigma<-sigma0}
    sigma4[c(1,2,3)]<-sigma;sigma5[1]<-sigma/m_nf;sigma5[2]<-(sigma+aa1)/m_nf;sigma5[3]<-sigma/m_nf
    ##############################criteria for iterations to stop######################## 
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
      sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*3
  ########first order genetic parameters###############  
  B2222<-solve(crossprod(hh2,hh2))%*%crossprod(hh2,m[c(1:3)])
  ########second order genetic parameters###############   
  F2jj<-sigmaF2-sigma
  if(F2jj<0){F2jj<-0}
  F2ll<-F2jj/sigmaF2
  F3jj<-sigmaF3-sigma/m_nf
  if(F3jj<0) {F3jj<-0}
  F3ll<-F3jj/sigmaF3   
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma))
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
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma))
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
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma))
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
    F2gg <- (dataF2 - m4[i])/sqrt(sigma4[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mi4[i]
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
  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3); 
  F3w1<-1/(12*m_sam[5])
  F3bmw <- matrix(0,m_sam[5],1); F3bmwsl <- matrix(0,m_sam[5],d1)
  for(i in 1:d1){
    F3gg <- (dataF3 - m5[i])/sqrt(sigma5[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mi5[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[5]){F3P2 <- F3P2+runif(m_sam[5])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[5]) + sum((F3P2 - (as.matrix(c(1:m_sam[5])) - 0.5)/m_sam[5])^2)
  F3u <- as.matrix(c(12*m_sam[5]*((F3dd[1]/m_sam[5]-0.5)^2),((45*m_sam[5])/4)*((F3dd[2]/m_sam[5]-1/3)^2),180*m_sam[5]*((F3dd[3]/m_sam[5]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[5]))))
  
  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("1MG-A",round(abc,4),round(AIC,4),round(m[1],4),round(m[2],4),round(m[3],4), round(t(m4),4)," "," "," "," "," "," ",round(sigma4[1],4),
                       round(t(mi4),4)," "," "," "," "," "," ",round(t(m5),4)," "," "," "," "," "," ",round(t(sigma5),4)," "," "," "," "," "," ",round(t(mi5),4)," "," "," "," "," "," ",round(sigma,4),          
                       round(B2222[1],4)," "," "," "," ",round(B2222[2],4)," "," "," "," "," "," "," "," "," ",round(F2jj,4),round(F2ll*100,4)," "," ",round(F3jj,4),round(F3ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
  output<-as.matrix(output)
  OUTPUT<-list(output,mix_pi4,mix_pi5)
  return(OUTPUT)
}
#######################1MG-EAD(A-3)##########################
G5ModelFun[[3]] <- function(K1,logL,df11,df21,df31,df41,df51,G5text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataF2 <- as.matrix(as.numeric(df41[,1]))
  dataF3 <- as.matrix(as.numeric(df51[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2),sum(dataF3)))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2),mean(dataF3)))
  sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1);sigmaF2<-var(dataF2);sigmaF3<-var(dataF3)
  sumx6<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  sigma<-sumx6;sigma0<-sigma;sigma40<-sigmaF2;sigma400<-sigma40;sigma50<-sigmaF3;sigma500<-sigma50
  m_nf <- as.numeric(G5text2)
  ###############procedure start###########################
  d1<-3;d2<-2
  mix_pi4<-as.matrix(c(0.75,0.25));mix_pi5<-as.matrix(c(0.25,0.5,0.25))
  sigma4<-matrix(0,2,1);sigma5<-matrix(0,3,1)
  sigma4[c(1:2)]<-sigma;sigma5[c(1,3)]<-sigma/m_nf
  m4<-as.matrix(c(m[1],m[3]))
  m5<-as.matrix(c(m[1],0.75*m[1]+0.25*m[3],m[3]))
  #######first order parameters############
  hh3<-matrix(c(1,1,1,1,-1,0.5),3,2)
  mm<-as.matrix(c(m[1],m[3],m5[2]))
  B3<-solve(crossprod(hh3,hh3))%*%crossprod(hh3,mm)
  a1<-B3[2]
  sigma5[2]<-sigma5[1]+(0.75*a1^2)/m_nf
  mi4<-mix_pi4[c(1,2)];mi5<-mix_pi5[c(1:3)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
    sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
  ##########iteration process###########
  ############E-step###################
  ############F2 process###############
  iteration <- 0; stopa <- 1000
  W4 <- matrix(0,2,m_sam[4]); swx24 <- matrix(0,2,1)
  W5 <- matrix(0,3,m_sam[5]); swx25 <- matrix(0,3,1)
  s0<-matrix(0,4,1);n0<-matrix(0,4,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:d2) { W4[i,] <- mix_pi4[i]*dnorm(dataF2,m4[i],sqrt(sigma4[i]))/dmixnorm(dataF2,m4,sqrt(sigma4),mix_pi4) }
    mi4 <- as.matrix(rowSums(W4)/m_sam[4])
    sumwx4 <- W4%*%dataF2
    #############F3 process############
    for(i in 1:d1) { W5[i,] <- mix_pi5[i]*dnorm(dataF3,m5[i],sqrt(sigma5[i]))/dmixnorm(dataF3,m5,sqrt(sigma5),mix_pi5) }
    mi5 <- as.matrix(rowSums(W5)/m_sam[5])
    sumwx5 <- W5%*%dataF3
    aaa0<-0
    s0[1]<-sumx[1]+sumx[2]+sumwx4[1]+m_nf*sumwx5[1];s0[2]<-sumx[3]+sumwx4[2]+m_nf*sumwx5[3];s0[3]<-sumwx5[2]
    n0[1]<-m_sam[1]+m_sam[2]+mi4[1]*m_sam[4]+mi5[1]*m_sam[5]*m_nf
    n0[2]<-m_sam[3]+mi4[2]*m_sam[4]+mi5[3]*m_sam[5]*m_nf
    n0[3]<-mi5[2]*m_sam[5]
    n0[c(1:3)][abs(n0[c(1:3)])<0.000001] <- 0.000001
    aa2<-3*s0[1]/n0[1]+s0[2]/n0[2]-4*s0[3]/n0[3]
    aa3<-9*sigma/n0[1]+sigma/n0[2]   
    aa1<-1000;n_iter<-0
    while(aa1>0.0001){
      n_iter<-n_iter+1 
      ########first order genetic parameters###############
      mm<-as.matrix(c(m[1],m[3],m5[2]))
      B31<-solve(crossprod(hh3,hh3))%*%crossprod(hh3,mm) 
      aa1<-B31[2];aa1<-0.75*aa1^2
      sigma5[2]<-(sigma+aa1)/m_nf
      ############restrictions#######################
      aaa1<-aa2/(aa3+16*sigma5[2]/n0[3])
      m[1]<-(s0[1]-3*aaa1*sigma)/n0[1]
      m[3]<-(s0[2]-aaa1*sigma)/n0[2]
      m5[2]<-(s0[3]+4*aaa1*sigma5[2])/n0[3]
      aa1<-abs(aaa1-aaa0)
      aaa0<-aaa1
      if(n_iter>20) break 
    }
    m[2]<-m4[1]<-m5[1]<-m[1];m4[2]<-m5[3]<-m[3]
    ##########obtain variance###############
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d2) {swx24[i] <- W4[i,]%*%(dataF2-m4[i])^2 } ;for(i in 1:d1) {swx25[i] <- W5[i,]%*%(dataF3-m5[i])^2 }      
    s0[4]<-ss1+ss2+ss3+swx24[1]+swx24[2]+m_nf*(swx25[1]+swx25[3])
    n0[4]<-m_sam[1]+m_sam[2]+m_sam[3]+m_sam[4]+(mi5[1]+mi5[3])*m_sam[5]
    aaa0<-sigma
    ########first order genetic parameters##############
    mm<-as.matrix(c(m[1],m[3],m5[2]))
    B311<-solve(crossprod(hh3,hh3))%*%crossprod(hh3,mm) 
    aa1<-B311[2];aa1<-0.75*aa1^2
    n_iter<-0;aa3<-1000
    while(aa3>0.0001){
      n_iter<-n_iter+1
      aa2<-sigma/(sigma+aa1)
      sigma<-(s0[4]+aa2^2*m_nf*swx25[2])/(n0[4]+aa2*n0[3])
      aa3<-abs(sigma-aaa0)
      aaa0<-sigma
      if(n_iter>20) break 
    }
    
    if(sigma<sigma0){sigma<-sigma0}
    sigma4[c(1,2)]<-sigma;sigma5[1]<-sigma/m_nf;sigma5[2]<-(sigma+aa1)/m_nf;sigma5[3]<-sigma/m_nf
    ##############################criteria for iterations to stop######################## 
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
      sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*3
  ########first order genetic parameters###############  
  mm<-as.matrix(c(m[1],m[3],m5[2]))
  B3111<-solve(crossprod(hh3,hh3))%*%crossprod(hh3,mm)
  ########second order genetic parameters###############   
  F2jj<-sigmaF2-sigma
  if(F2jj<0){F2jj<-0}
  F2ll<-F2jj/sigmaF2
  F3jj<-sigmaF3-sigma/m_nf
  if(F3jj<0) {F3jj<-0}
  F3ll<-F3jj/sigmaF3   
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma))
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
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma))
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
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma))
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
    F2gg <- (dataF2 - m4[i])/sqrt(sigma4[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mi4[i]
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
  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3); 
  F3w1<-1/(12*m_sam[5])
  F3bmw <- matrix(0,m_sam[5],1); F3bmwsl <- matrix(0,m_sam[5],d1)
  for(i in 1:d1){
    F3gg <- (dataF3 - m5[i])/sqrt(sigma5[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mi5[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[5]){F3P2 <- F3P2+runif(m_sam[5])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[5]) + sum((F3P2 - (as.matrix(c(1:m_sam[5])) - 0.5)/m_sam[5])^2)
  F3u <- as.matrix(c(12*m_sam[5]*((F3dd[1]/m_sam[5]-0.5)^2),((45*m_sam[5])/4)*((F3dd[2]/m_sam[5]-1/3)^2),180*m_sam[5]*((F3dd[3]/m_sam[5]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[5]))))
  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("1MG-EAD",round(abc,4),round(AIC,4),round(m[1],4),round(m[2],4),round(m[3],4), round(t(m4),4)," "," "," "," "," "," "," ",round(sigma4[1],4),
                       round(t(mi4),4)," "," "," "," "," "," "," ",round(t(m5),4)," "," "," "," "," "," ",round(t(sigma5),4)," "," "," "," "," "," ",round(t(mi5),4)," "," "," "," "," "," ",round(sigma,4),          
                       round(B3111[1],4)," "," "," "," ",round(B3111[2],4)," "," "," "," "," "," "," "," "," ",round(F2jj,4),round(F2ll*100,4)," "," ",round(F3jj,4),round(F3ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
 output<-as.matrix(output)
  OUTPUT<-list(output,mix_pi4,mix_pi5)
  return(OUTPUT)
}
##################1MG-NCD(A-4)######################
G5ModelFun[[4]] <- function(K1,logL,df11,df21,df31,df41,df51,G5text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataF2 <- as.matrix(as.numeric(df41[,1]))
  dataF3 <- as.matrix(as.numeric(df51[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2),sum(dataF3)))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2),mean(dataF3)))
  sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1);sigmaF2<-var(dataF2);sigmaF3<-var(dataF3)
  sumx6<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  sigma<-sumx6;sigma0<-sigma;sigma40<-sigmaF2;sigma400<-sigma40;sigma50<-sigmaF3;sigma500<-sigma50
  m_nf <- as.numeric(G5text2)
  ###############procedure start###########################
  d1<-3;d2<-2
  mix_pi4<-as.matrix(c(0.25,0.75));mix_pi5<-as.matrix(c(0.25,0.5,0.25))
  sigma4<-matrix(0,2,1);sigma5<-matrix(0,3,1)
  sigma4[c(1:2)]<-sigma;sigma5[c(1,3)]<-sigma/m_nf
  m4<-as.matrix(c(m[1],m[3]))
  m5<-as.matrix(c(m[1],0.25*m[1]+0.75*m[3],m[3]))
  #######first order parameters############
  hh4<-matrix(c(1,1,1,1,-1,-0.5),3,2)
  mm<-as.matrix(c(m[1],m[3],m5[2]))
  B4<-solve(crossprod(hh4,hh4))%*%crossprod(hh4,mm)
  a1<-B4[2]
  sigma5[2]<-sigma5[1]+(0.75*a1^2)/m_nf
  mi4<-mix_pi4[c(1,2)];mi5<-mix_pi5[c(1:3)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
    sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
  ##########iteration process###########
  ############E-step###################
  ############F2 process###############
  iteration <- 0; stopa <- 1000
  W4 <- matrix(0,2,m_sam[4]); swx24 <- matrix(0,2,1)
  W5 <- matrix(0,3,m_sam[5]); swx25 <- matrix(0,3,1)
  s0<-matrix(0,4,1);n0<-matrix(0,4,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:d2) { W4[i,] <- mix_pi4[i]*dnorm(dataF2,m4[i],sqrt(sigma4[i]))/dmixnorm(dataF2,m4,sqrt(sigma4),mix_pi4) }
    mi4 <- as.matrix(rowSums(W4)/m_sam[4])
    sumwx4 <- W4%*%dataF2
    #############F3 process############
    for(i in 1:d1) { W5[i,] <- mix_pi5[i]*dnorm(dataF3,m5[i],sqrt(sigma5[i]))/dmixnorm(dataF3,m5,sqrt(sigma5),mix_pi5) }
    mi5 <- as.matrix(rowSums(W5)/m_sam[5])
    sumwx5 <- W5%*%dataF3
    aaa0<-0
    s0[1]<-sumx[1]+sumwx4[1]+m_nf*sumwx5[1];s0[2]<-sumx[2]+sumx[3]+sumwx4[2]+m_nf*sumwx5[3];s0[3]<-sumwx5[2]
    n0[1]<-m_sam[1]+mi4[1]*m_sam[4]+mi5[1]*m_sam[5]*m_nf
    n0[2]<-m_sam[2]+m_sam[3]+mi4[2]*m_sam[4]+mi5[3]*m_sam[5]*m_nf
    n0[3]<-mi5[2]*m_sam[5]
    n0[c(1:3)][abs(n0[c(1:3)])<0.000001] <- 0.000001
    aa2<-s0[1]/n0[1]+3*s0[2]/n0[2]-4*s0[3]/n0[3]
    aa3<-sigma*(1/n0[1]+9/n0[2])   
    aa1<-1000;n_iter<-0
    while(aa1>0.0001){
      n_iter<-n_iter+1 
      ########first order genetic parameters###############
      mm<-as.matrix(c(m[1],m[3],m5[2]))
      B41<-solve(crossprod(hh4,hh4))%*%crossprod(hh4,mm) 
      aa1<-B41[2];aa1<-0.75*aa1^2
      sigma5[2]<-(sigma+aa1)/m_nf
      ############restrictions#######################
      aaa1<-aa2/(aa3+16*sigma5[2]/n0[3])
      m[1]<-(s0[1]-aaa1*sigma)/n0[1]
      m[3]<-(s0[2]-3*aaa1*sigma)/n0[2]
      m5[2]<-(s0[3]+4*aaa1*sigma5[2])/n0[3]
      aa1<-abs(aaa1-aaa0)
      aaa0<-aaa1
      if(n_iter>20) break 
    }
    m[2]<-m4[2]<-m5[3]<-m[3]
    m4[1]<-m5[1]<-m[1]
    ##########obtain variance###############
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d2) {swx24[i] <- W4[i,]%*%(dataF2-m4[i])^2 } ;for(i in 1:d1) {swx25[i] <- W5[i,]%*%(dataF3-m5[i])^2 }      
    s0[4]<-ss1+ss2+ss3+swx24[1]+swx24[2]+m_nf*(swx25[1]+swx25[3])
    n0[4]<-m_sam[1]+m_sam[2]+m_sam[3]+m_sam[4]+(mi5[1]+mi5[3])*m_sam[5]
    aaa0<-sigma
    ########first order genetic parameters##############
    mm<-as.matrix(c(m[1],m[3],m5[2]))
    B411<-solve(crossprod(hh4,hh4))%*%crossprod(hh4,mm) 
    aa1<-B411[2];aa1<-0.75*aa1^2
    n_iter<-0;aa3<-1000
    while(aa3>0.0001){
      n_iter<-n_iter+1
      aa2<-sigma/(sigma+aa1)
      sigma<-(s0[4]+aa2^2*m_nf*swx25[2])/(n0[4]+aa2*mi5[2]*m_sam[5])
      aa3<-abs(sigma-aaa0)
      aaa0<-sigma
      if(n_iter>20) break 
    }
    if(sigma<sigma0){sigma<-sigma0}
    sigma4[c(1,2)]<-sigma;sigma5[1]<-sigma/m_nf;sigma5[2]<-(sigma+aa1)/m_nf;sigma5[3]<-sigma/m_nf
    ##############################criteria for iterations to stop######################## 
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
      sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*3
  ########first order genetic parameters###############    
  mm<-as.matrix(c(m[1],m[3],m5[2]))
  B4111<-solve(crossprod(hh4,hh4))%*%crossprod(hh4,mm)
  ########second order genetic parameters###############   
  F2jj<-sigmaF2-sigma
  if(F2jj<0){F2jj<-0}
  F2ll<-F2jj/sigmaF2
  F3jj<-sigmaF3-sigma/m_nf
  if(F3jj<0) {F3jj<-0}
  F3ll<-F3jj/sigmaF3   
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma))
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
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma))
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
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma))
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
    F2gg <- (dataF2 - m4[i])/sqrt(sigma4[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mi4[i]
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
  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3); 
  F3w1<-1/(12*m_sam[5])
  F3bmw <- matrix(0,m_sam[5],1); F3bmwsl <- matrix(0,m_sam[5],d1)
  for(i in 1:d1){
    F3gg <- (dataF3 - m5[i])/sqrt(sigma5[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mi5[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[5]){F3P2 <- F3P2+runif(m_sam[5])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[5]) + sum((F3P2 - (as.matrix(c(1:m_sam[5])) - 0.5)/m_sam[5])^2)
  F3u <- as.matrix(c(12*m_sam[5]*((F3dd[1]/m_sam[5]-0.5)^2),((45*m_sam[5])/4)*((F3dd[2]/m_sam[5]-1/3)^2),180*m_sam[5]*((F3dd[3]/m_sam[5]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[5]))))
  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("1MG-NCD",round(abc,4),round(AIC,4),round(m[1],4),round(m[2],4),round(m[3],4), round(t(m4),4)," "," "," "," "," "," "," ",round(sigma4[1],4),
                       round(t(mi4),4)," "," "," "," "," "," "," ",round(t(m5),4)," "," "," "," "," "," ",round(t(sigma5),4)," "," "," "," "," "," ",round(t(mi5),4)," "," "," "," "," "," ",round(sigma,4),          
                       round(B4111[1],4)," "," "," "," ",round(B4111[2],4)," "," "," "," "," "," "," "," "," ",round(F2jj,4),round(F2ll*100,4)," "," ",round(F3jj,4),round(F3ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
 output<-as.matrix(output)
  OUTPUT<-list(output,mix_pi4,mix_pi5)
  return(OUTPUT)
}
###############2MG-ADI(B-1)##########################################
G5ModelFun[[5]] <- function(K1,logL,df11,df21,df31,df41,df51,G5text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataF2 <- as.matrix(as.numeric(df41[,1]))
  dataF3 <- as.matrix(as.numeric(df51[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2),sum(dataF3)))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2),mean(dataF3)))
  sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1);sigmaF2<-var(dataF2);sigmaF3<-var(dataF3)
  sumx6<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  sigma<-sumx6;sigma0<-sigma;sigma40<-sigmaF2;sigma400<-sigma40;sigma50<-sigmaF3;sigma500<-sigma50
  m_nf <- as.numeric(G5text2)
  ###############procedure start###########################
  mix_pi4<-as.matrix(c(0.0625,0.125,0.0625,0.1250,0.250,0.125,0.0625,0.125,0.0625))
  mix_pi5<-as.matrix(c(0.0625,0.125,0.0625,0.1250,0.250,0.125,0.0625,0.125,0.0625))
  a1<-sqrt(sigma40/m_sam[4])
  if (m[1]<m[3]) a1<--a1
  m4<-as.matrix(c(m[4]+2.8*a1,m[4]+2.1*a1,m[4]+1.4*a1,m[4]+0.7*a1,m[4],m[4]-0.7*a1,m[4]-1.4*a1,m[4]-2.1*a1,m[4]-2.8*a1))
  a1<-sqrt(sigma50/m_sam[5])
  if (m[1]<m[3]) {a1<--a1}
  m5<-as.matrix(c(m[5]+2.8*a1,m[5]+2.1*a1,m[5]+1.4*a1,m[5]+0.7*a1,m[5],m[5]-0.7*a1,m[5]-1.4*a1,m[5]-2.1*a1,m[5]-2.8*a1))
  sigma<-sigma400/2
  #######first order parameters############
  hh5<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,-1,1,1,0,0,-1,-1,1,0,0,0,-1,
                1,0,-1,1,-1,1,-1,1,0,0,1,0,-1,0,0,1,0,0,0,1,1,0,0,0,0.5,0.5,0.5,0,
                0,1,0,0,0,0,0,0,1,0.5,0,0.5,0,0.5,1,0,1,0,-1,0,0,-1,0,0,0,0,0,0,
                0,0,0,1,0,0,0,0,-1,0.5,0,0,0,-0.5,0,0,0,0,0,1,-1,0,0,0,0.5,0,-0.5,0,
                0,1,0,0,0,0,0,0,0,0,0,0.25,0,0),14,9)
  mm<-as.matrix(c(m[1],m[2],m[3],m4[2],m4[3],m4[4],m4[6],m4[7],m4[8],m5[2],m5[4],m5[5],m5[6],m5[8]))
  B5<-solve(crossprod(hh5,hh5))%*%crossprod(hh5,mm)
  gs<-B5[c(2:9)]
  sigma4<-matrix(0,9,1);sigma5<-matrix(0,9,1)
  sigma4[c(1:9)]<-sigma;sigma5[c(1,3,7,9)]<-sigma/m_nf
  g_aa1<-0.5*(gs[2]+gs[5])^2+0.25*(gs[4]+gs[6])^2   #  0.5(db+i)**2+0.25(hb+jab)**2.
  g_aa2<-0.5*(gs[1]+gs[5])^2+0.25*(gs[3]+gs[7])^2   #  0.5(da+i)**2+0.25(ha+jba)**2.
  g_aa3<-0.5*(gs[1]-gs[5])^2+0.25*(gs[3]-gs[7])^2   #  0.5(da-i)**2+0.25(ha-jba)**2.
  g_aa4<-0.5*(gs[2]-gs[5])^2+0.25*(gs[4]-gs[6])^2   #  0.5(db-i)**2+0.25(hb-jab)**2.
  g_aa5<-0.25*(gs[1]^2+gs[2]^2+gs[5]^2+(gs[1]+gs[6])^2+(gs[2]+gs[7])^2+(gs[3]+gs[8]/2)^2+(gs[4]+gs[8]/2)^2+gs[8]^2/4)
  sigma5[2]<-(sigma+g_aa1)/m_nf;sigma5[4]<-(sigma+g_aa2)/m_nf
  sigma5[5]<-(sigma+g_aa5)/m_nf;sigma5[6]<-(sigma+g_aa3)/m_nf
  sigma5[8]<-(sigma+g_aa4)/m_nf
  mi4<-mix_pi4[c(1:9)];mi5<-mix_pi5[c(1:9)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
    sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
  ########################################################
  ##########iteration process###########
  ############E-step###################
  ############F2 process###############
  d3<-9
  iteration <- 0; stopa <- 1000
  W4 <- matrix(0,9,m_sam[4]); swx24 <- matrix(0,9,1)
  W5 <- matrix(0,9,m_sam[5]); swx25 <- matrix(0,9,1)
  hh<-matrix(0,5,5);b_line<-matrix(0,5,1)
  s0<-matrix(0,20,1);n0<-matrix(0,20,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:d3) { W4[i,] <- mix_pi4[i]*dnorm(dataF2,m4[i],sqrt(sigma4[i]))/dmixnorm(dataF2,m4,sqrt(sigma4),mix_pi4) }
    mi4 <- as.matrix(rowSums(W4)/m_sam[4])
    sumwx4 <- W4%*%dataF2
    #############F3 process############
    for(i in 1:d3) { W5[i,] <- mix_pi5[i]*dnorm(dataF3,m5[i],sqrt(sigma5[i]))/dmixnorm(dataF3,m5,sqrt(sigma5),mix_pi5) }
    mi5 <- as.matrix(rowSums(W5)/m_sam[5])
    sumwx5 <- W5%*%dataF3
    s0[1]<-sumx[1]+sumwx4[1]+m_nf*sumwx5[1];s0[2]<-sumx[2]+sumwx4[5]
    s0[3]<-sumx[3]+sumwx4[9]+m_nf*sumwx5[9];s0[4]<-sumwx4[2]
    s0[5]<-sumwx4[3]+m_nf*sumwx5[3];s0[6]<-sumwx4[4]
    s0[7]<-sumwx4[6];s0[8]<-sumwx4[7]+m_nf*sumwx5[7]
    s0[9]<-sumwx4[8];s0[10]<-sumwx5[2]
    s0[11]<-sumwx5[4];s0[12]<-sumwx5[5]
    s0[13]<-sumwx5[6];s0[14]<-sumwx5[8]
    nn<-m_sam
    n0[1]<-nn[1]+mi4[1]*nn[4]+m_nf*mi5[1]*nn[5]
    n0[2]<-nn[2]+mi4[5]*nn[4]
    n0[3]<-nn[3]+mi4[9]*nn[4]+m_nf*mi5[9]*nn[5]
    n0[4]<-mi4[2]*nn[4]
    n0[5]<-mi4[3]*nn[4]+m_nf*mi5[3]*nn[5]
    n0[6]<-mi4[4]*nn[4];n0[7]<-mi4[6]*nn[4]
    n0[8]<-mi4[7]*nn[4]+m_nf*mi5[7]*nn[5]
    n0[9]<-mi4[8]*nn[4];n0[10]<-mi5[2]*nn[5]
    n0[11]<-mi5[4]*nn[5];n0[12]<-mi5[5]*nn[5]
    n0[13]<-mi5[6]*nn[5];n0[14]<-mi5[8]*nn[5]
    n0[c(1:14)][abs(n0[c(1:14)])<0.000001]<-0.000001
    aaa1<-1000;n_iter<-0;AA<-matrix(0,5,1)
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      ##########first order parameters###################
      mm<-as.matrix(c(m[1],m[2],m[3],m4[2],m4[3],m4[4],m4[6],m4[7],m4[8],m5[2],m5[4],m5[5],m5[6],m5[8]))
      B51<-solve(crossprod(hh5,hh5))%*%crossprod(hh5,mm)
      gs<-B51[c(2:9)]
      sigma4[c(1:9)]<-sigma;sigma5[c(1,3,7,9)]<-sigma/m_nf
      g_aa1<-0.5*(gs[2]+gs[5])^2+0.25*(gs[4]+gs[6])^2    #   0.5(db+i)**2+0.25(hb+jab)**2.
      g_aa2<-0.5*(gs[1]+gs[5])^2+0.25*(gs[3]+gs[7])^2    #   0.5(da+i)**2+0.25(ha+jba)**2.
      g_aa3<-0.5*(gs[1]-gs[5])^2+0.25*(gs[3]-gs[7])^2    #   0.5(da-i)**2+0.25(ha-jba)**2.
      g_aa4<-0.5*(gs[2]-gs[5])^2+0.25*(gs[4]-gs[6])^2    #   0.5(db-i)**2+0.25(hb-jab)**2.
      g_aa5<-0.25*(gs[1]^2+gs[2]^2+gs[5]^2+(gs[1]+gs[6])^2+(gs[2]+gs[7])^2+(gs[3]+gs[8]/2)^2+(gs[4]+gs[8]/2)^2+gs[8]^2/4)
      sigma5[2]<-(sigma+g_aa1)/m_nf;sigma5[4]<-(sigma+g_aa2)/m_nf
      sigma5[5]<-(sigma+g_aa5)/m_nf;sigma5[6]<-(sigma+g_aa3)/m_nf
      sigma5[8]<-(sigma+g_aa4)/m_nf
      #################restrictions########################################
      hh[1,1]<-sigma*(1/n0[1]+16/n0[2]+1/n0[3]+4/n0[4]+1/n0[5]+4/n0[6]+4/n0[7]+1/n0[8]+4/n0[9])+256*sigma5[5]/n0[12]
      hh[1,2]<-sigma*(1/n0[1]+4/n0[6]+1/n0[8])
      hh[1,3]<-sigma*(1/n0[1]+4/n0[4]+1/n0[5])
      hh[1,4]<-sigma*(1/n0[3]+1/n0[5]+4/n0[7])
      hh[1,5]<-sigma*(1/n0[3]+1/n0[8]+4/n0[9])
      hh[2,2]<-sigma*(1/n0[1]+4/n0[6]+1/n0[8])+16*sigma5[4]/n0[11]
      hh[2,3]<-sigma/n0[1];hh[2,4]<-0
      hh[2,5]<-sigma/n0[8]
      hh[3,3]<-sigma*(1/n0[1]+4/n0[4]+1/n0[5])+16*sigma5[2]/n0[10]
      hh[3,4]<-sigma/n0[5];hh[3,5]<-0;
      hh[4,3]<-hh[3,4]
      hh[4,4]<-sigma*(1/n0[3]+1/n0[5]+4/n0[7])+16*sigma5[6]/n0[13]
      hh[4,5]<-sigma/n0[3]
      hh[5,5]<-sigma*(1/n0[3]+1/n0[8]+4/n0[9])+16*sigma5[8]/n0[14]
      for(i in 2:5)
      {
        for(j in 1:(i-1))
        {
          hh[i,j]<-hh[j,i]
        }
      }
      ########################################################
      b_line[1]<-s0[1]/n0[1]+4*s0[2]/n0[2]+s0[3]/n0[3]+2*s0[4]/n0[4]+s0[5]/n0[5]+2*s0[6]/n0[6]+2*s0[7]/n0[7]+s0[8]/n0[8]+2*s0[9]/n0[9]-16*s0[12]/n0[12]
      b_line[2]<-s0[1]/n0[1]+2*s0[6]/n0[6]+s0[8]/n0[8]-4*s0[11]/n0[11]
      b_line[3]<-s0[1]/n0[1]+2*s0[4]/n0[4]+s0[5]/n0[5]-4*s0[10]/n0[10]
      b_line[4]<-s0[3]/n0[3]+s0[5]/n0[5]+2*s0[7]/n0[7]-4*s0[13]/n0[13]
      b_line[5]<-s0[3]/n0[3]+s0[8]/n0[8]+2*s0[9]/n0[9]-4*s0[14]/n0[14]
      B511<-solve(hh,b_line)
      ########################################################
      m[1]<-(s0[1]-sigma*(B511[1]+B511[2]+B511[3]))/n0[1]
      m[2]<-(s0[2]-4*B511[1]*sigma)/n0[2]
      m[3]<-(s0[3]-sigma*(B511[1]+B511[4]+B511[5]))/n0[3]
      m4[2]<-(s0[4]-2*sigma*(B511[1]+B511[3]))/n0[4]
      m4[3]<-(s0[5]-sigma*(B511[1]+B511[3]+B511[4]))/n0[5]
      m4[4]<-(s0[6]-sigma*2*(B511[1]+B511[2]))/n0[6]
      m4[6]<-(s0[7]-sigma*2*(B511[1]+B511[4]))/n0[7]
      m4[7]<-(s0[8]-sigma*(B511[1]+B511[2]+B511[5]))/n0[8]
      m4[8]<-(s0[9]-sigma*2*(B511[1]+B511[5]))/n0[9]
      m5[2]<-(s0[10]+4*B511[3]*sigma5[2])/n0[10]
      m5[4]<-(s0[11]+4*B511[2]*sigma5[4])/n0[11]
      m5[5]<-(s0[12]+16*B511[1]*sigma5[5])/n0[12]
      m5[6]<-(s0[13]+4*B511[4]*sigma5[6])/n0[13]
      m5[8]<-(s0[14]+4*B511[5]*sigma5[8])/n0[14]
      m4[1]<-m[1];m4[5]<-m[2];m4[9]<-m[3]
      m5[1]<-m[1];m5[3]<-m4[3];m5[7]<-m4[7]
      m5[9]<-m[3]
      aaa1<-max(abs(B511-AA))
      AA<-B511
      if (n_iter>20) break
    } 
    ##########obtain variance###############
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d3) {swx24[i] <- W4[i,]%*%(dataF2-m4[i])^2 } ;for(i in 1:d3) {swx25[i] <- W5[i,]%*%(dataF3-m5[i])^2 }      
    s0[15]<-sum(swx24)+ss1+ss2+ss3
    n0[15]<-nn[1]+nn[2]+nn[3]+nn[4]+nn[5]*(mi5[1]+mi5[3]+mi5[7]+mi5[9])
    n0[16]<-mi5[2]*nn[5];n0[17]<-mi5[4]*nn[5];n0[18]<-mi5[5]*nn[5]
    n0[19]<-mi5[6]*nn[5];n0[20]<-mi5[8]*nn[5]
    mm<-as.matrix(c(m[1],m[2],m[3],m4[2],m4[3],m4[4],m4[6],m4[7],m4[8],m5[2],m5[4],m5[5],m5[6],m5[8]))
    B5111<-solve(crossprod(hh5,hh5))%*%crossprod(hh5,mm)
    gs<-B5111[c(2:9)]
    g_aa1<-0.5*(gs[2]+gs[5])^2+0.25*(gs[4]+gs[6])^2   # 0.5(db+i)**2+0.25(hb+jab)**2.
    g_aa2<-0.5*(gs[1]+gs[5])^2+0.25*(gs[3]+gs[7])^2   # 0.5(da+i)**2+0.25(ha+jba)**2.
    g_aa3<-0.5*(gs[1]-gs[5])^2+0.25*(gs[3]-gs[7])^2   # 0.5(da-i)**2+0.25(ha-jba)**2.
    g_aa4<-0.5*(gs[2]-gs[5])^2+0.25*(gs[4]-gs[6])^2   # 0.5(db-i)**2+0.25(hb-jab)**2.
    g_aa5<-0.25*(gs[1]^2+gs[2]^2+gs[5]^2+(gs[1]+gs[6])^2+(gs[2]+gs[7])^2+(gs[3]+gs[8]/2)^2+(gs[4]+gs[8]/2)^2+gs[8]^2/4.0)
    aaa0<-sigma;n_iter<-0;aaa1<-1000
    while (aaa1>0.0001){
      n_iter<-n_iter+1
      aa1<-sigma/(sigma+g_aa1);aa2<-sigma/(sigma+g_aa2)
      aa3<-sigma/(sigma+g_aa5);aa4<-sigma/(sigma+g_aa3)
      aa5<-sigma/(sigma+g_aa4)
      sigma<-(s0[15]+m_nf*(swx25[1]+swx25[3]+swx25[7]+swx25[9]+aa1^2*swx25[2]+aa2^2*swx25[4]+aa3^2*swx25[5]+aa4^2*swx25[6]+aa5^2*swx25[8]))/(n0[15]+aa1*n0[16]+aa2*n0[17]+aa3*n0[18]+aa4*n0[19]+aa5*n0[20])
      aaa1<-abs(sigma-aaa0)
      aaa0<-sigma
      if (n_iter>20) break 
    }
    if(sigma<sigma0){sigma<-sigma0}
    sigma4[c(1:9)]<-sigma
    sigma5[c(1,3,7,9)]<-sigma/m_nf;sigma5[2]<-(sigma+g_aa1)/m_nf;sigma5[4]<-(sigma+g_aa2)/m_nf
    sigma5[5]<-(sigma+g_aa5)/m_nf;sigma5[6]<-(sigma+g_aa3)/m_nf;sigma5[8]<-(sigma+g_aa4)/m_nf
    ##############################criteria for iterations to stop######################## 
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
      sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*10
  #######first order parameters############
  mm<-as.matrix(c(m[1],m[2],m[3],m4[2],m4[3],m4[4],m4[6],m4[7],m4[8],m5[2],m5[4],m5[5],m5[6],m5[8]))
  B51111<-solve(crossprod(hh5,hh5))%*%crossprod(hh5,mm)
  ########second order genetic parameters###############   
  F2jj<-sigmaF2-sigma
  if(F2jj<0){F2jj<-0}
  F2ll<-F2jj/sigmaF2
  F3jj<-sigmaF3-sigma/m_nf
  if(F3jj<0) {F3jj<-0}
  F3ll<-F3jj/sigmaF3   
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma))
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
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma))
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
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma))
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
    F2gg <- (dataF2 - m4[i])/sqrt(sigma4[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mi4[i]
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
  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3); 
  F3w1<-1/(12*m_sam[5])
  F3bmw <- matrix(0,m_sam[5],1); F3bmwsl <- matrix(0,m_sam[5],d3)
  for(i in 1:d3){
    F3gg <- (dataF3 - m5[i])/sqrt(sigma5[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mi5[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[5]){F3P2 <- F3P2+runif(m_sam[5])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[5]) + sum((F3P2 - (as.matrix(c(1:m_sam[5])) - 0.5)/m_sam[5])^2)
  F3u <- as.matrix(c(12*m_sam[5]*((F3dd[1]/m_sam[5]-0.5)^2),((45*m_sam[5])/4)*((F3dd[2]/m_sam[5]-1/3)^2),180*m_sam[5]*((F3dd[3]/m_sam[5]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[5]))))
  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("2MG-ADI",round(abc,4),round(AIC,4),round(m[1],4),round(m[2],4),round(m[3],4), round(t(m4),4),round(sigma4[1],4),
                       round(t(mi4),4),round(t(m5),4),round(t(sigma5),4),round(t(mi5),4),round(sigma,4),          
                       round(B51111[1],4)," "," "," "," ",round(B51111[2],4),round(B51111[3],4),round(B51111[4],4),round(B51111[5],4),round(B51111[6],4),round(B51111[7],4),round(B51111[8],4),round(B51111[9],4)," "," ",round(F2jj,4),round(F2ll*100,4)," "," ",round(F3jj,4),round(F3ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
 output<-as.matrix(output)
  OUTPUT<-list(output,mix_pi4,mix_pi5)
  return(OUTPUT)
}

#######################2MG-AD(B-2)#############################
G5ModelFun[[6]] <- function(K1,logL,df11,df21,df31,df41,df51,G5text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataF2 <- as.matrix(as.numeric(df41[,1]))
  dataF3 <- as.matrix(as.numeric(df51[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2),sum(dataF3)))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2),mean(dataF3)))
  sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1);sigmaF2<-var(dataF2);sigmaF3<-var(dataF3)
  sumx6<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  sigma<-sumx6;sigma0<-sigma;sigma40<-sigmaF2;sigma400<-sigma40;sigma50<-sigmaF3;sigma500<-sigma50
  m_nf <- as.numeric(G5text2)
  ###############procedure start###########################
  mix_pi4<-as.matrix(c(0.0625,0.125,0.0625,0.1250,0.250,0.125,0.0625,0.125,0.0625))
  mix_pi5<-as.matrix(c(0.0625,0.125,0.0625,0.1250,0.250,0.125,0.0625,0.125,0.0625))
  a1<-sqrt(sigma40/m_sam[4])
  if (m[1]<m[3]) a1<--a1
  m4<-as.matrix(c(m[4]+2.8*a1,m[4]+2.1*a1,m[4]+1.4*a1,m[4]+0.7*a1,m[4],m[4]-0.7*a1,m[4]-1.4*a1,m[4]-2.1*a1,m[4]-2.8*a1))
  a1<-sqrt(sigma50/m_sam[5])
  if (m[1]<m[3]) {a1<--a1}
  m5<-as.matrix(c(m[5]+2.8*a1,m[5]+2.1*a1,m[5]+1.4*a1,m[5]+0.7*a1,m[5],m[5]-0.7*a1,m[5]-1.4*a1,m[5]-2.1*a1,m[5]-2.8*a1))
  sigma<-sigma400/2
  #######first order parameters############
  hh6<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,-1,1,1,0,0,-1,-1,1,0,0,0,-1,
                1,0,-1,1,-1,1,-1,1,0,0,1,0,-1,0,0,1,0,0,0,1,1,0,0,0,0.5,0.5,0.5,0,
                0,1,0,0,0,0,0,0,1,0.5,0,0.5,0,0.5),14,5)
  mm<-as.matrix(c(m[1],m[2],m[3],m4[2],m4[3],m4[4],m4[6],m4[7],m4[8],m5[2],m5[4],m5[5],m5[6],m5[8]))
  B6<-solve(crossprod(hh6,hh6))%*%crossprod(hh6,mm)
  gs<-B6[c(2:5)]
  sigma4<-matrix(0,9,1);sigma5<-matrix(0,9,1)
  sigma4[c(1:9)]<-sigma;sigma5[c(1,3,7,9)]<-sigma/m_nf
  g_aa1<-0.5*gs[2]^2+0.25*gs[4]^2 #   0.5*db**2+0.25*hb**2.
  g_aa2<-0.5*gs[1]^2+0.25*gs[3]^2 #   0.5*da**2+0.25*ha**2.
  g_aa3<-g_aa1+g_aa2 #  0.5(da**2+db**2)+0.25(ha**2+hb**2).
  sigma5[2]<-(sigma+g_aa1)/m_nf;sigma5[4]<-(sigma+g_aa2)/m_nf
  sigma5[5]<-(sigma+g_aa3)/m_nf;sigma5[6]<-sigma5[4]
  sigma5[8]<-sigma5[2]
  mi4<-mix_pi4[c(1:9)];mi5<-mix_pi5[c(1:9)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
    sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
  ########################################################
  ##########iteration process###########
  ############E-step###################
  ############F2 process###############
  d3<-9
  iteration <- 0; stopa <- 1000
  W4 <- matrix(0,9,m_sam[4]); swx24 <- matrix(0,9,1)
  W5 <- matrix(0,9,m_sam[5]); swx25 <- matrix(0,9,1)
  hh<-matrix(0,9,9);b_line<-matrix(0,9,1)
  s0<-matrix(0,20,1); n0<-matrix(0,20,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:d3) { W4[i,] <- mix_pi4[i]*dnorm(dataF2,m4[i],sqrt(sigma4[i]))/dmixnorm(dataF2,m4,sqrt(sigma4),mix_pi4) }
    mi4 <- as.matrix(rowSums(W4)/m_sam[4])
    sumwx4 <- W4%*%dataF2
    #############F3 process############
    for(i in 1:d3) { W5[i,] <- mix_pi5[i]*dnorm(dataF3,m5[i],sqrt(sigma5[i]))/dmixnorm(dataF3,m5,sqrt(sigma5),mix_pi5) }
    mi5 <- as.matrix(rowSums(W5)/m_sam[5])
    sumwx5 <- W5%*%dataF3
    aaa0<-0;nn<-m_sam
    s0[1]<-sumx[1]+sumwx4[1]+m_nf*sumwx5[1];s0[2]<-sumx[2]+sumwx4[5]
    s0[3]<-sumx[3]+sumwx4[9]+m_nf*sumwx5[9];s0[4]<-sumwx4[2]
    s0[5]<-sumwx4[3]+m_nf*sumwx5[3];s0[6]<-sumwx4[4]
    s0[7]<-sumwx4[6];s0[8]<-sumwx4[7]+m_nf*sumwx5[7]
    s0[9]<-sumwx4[8];s0[10]<-sumwx5[2]
    s0[11]<-sumwx5[4];s0[12]<-sumwx5[5]
    s0[13]<-sumwx5[6];s0[14]<-sumwx5[8]
    ###################################
    n0[1]<-nn[1]+mi4[1]*nn[4]+m_nf*mi5[1]*nn[5];n0[2]<-nn[2]+mi4[5]*nn[4]
    n0[3]<-nn[3]+mi4[9]*nn[4]+m_nf*mi5[9]*nn[5];n0[4]<-mi4[2]*nn[4]
    n0[5]<-mi4[3]*nn[4]+m_nf*mi5[3]*nn[5];n0[6]<-mi4[4]*nn[4]
    n0[7]<-mi4[6]*nn[4];n0[8]<-mi4[7]*nn[4]+m_nf*mi5[7]*nn[5]
    n0[9]<-mi4[8]*nn[4];n0[10]<-mi5[2]*nn[5]
    n0[11]<-mi5[4]*nn[5];n0[12]<-mi5[5]*nn[5]
    n0[13]<-mi5[6]*nn[5];n0[14]<-mi5[8]*nn[5]
    n0[c(1:14)][abs(n0[c(1:14)])<0.000001]<-0.000001
    aaa1<-1000;n_iter<-0;AA<-matrix(0,9,1)
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      ##########first order parameters###################
      mm<-as.matrix(c(m[1],m[2],m[3],m4[2],m4[3],m4[4],m4[6],m4[7],m4[8],m5[2],m5[4],m5[5],m5[6],m5[8]))
      B61<-solve(crossprod(hh6,hh6))%*%crossprod(hh6,mm)
      gs<-B61[c(2:5)]
      sigma4[c(1:9)]<-sigma;sigma5[c(1,3,7,9)]<-sigma/m_nf
      g_aa1<-0.5*gs[2]^2+0.25*gs[4]^2 #   0.5*db**2+0.25*hb**2.
      g_aa2<-0.5*gs[1]^2+0.25*gs[3]^2 #   0.5*da**2+0.25*ha**2.
      g_aa3<-g_aa1+g_aa2 #   0.5(da**2+db**2)+0.25(ha**2+hb**2).
      sigma5[2]<-(sigma+g_aa1)/m_nf;sigma5[4]<-(sigma+g_aa2)/m_nf
      sigma5[5]<-(sigma+g_aa3)/m_nf;sigma5[6]<-sigma5[4]
      sigma5[8]<-sigma5[2]
      #################restrictions########################################
      hh[1,1]<-sigma*(1/n0[1]+4/n0[2]+1/n0[3])+16*sigma5[5]/n0[12]
      hh[1,2]<-sigma/n0[1]
      hh[1,3]<-sigma/n0[1]
      hh[1,4]<-sigma*(1/n0[1]+1/n0[3])
      hh[1,5]<-sigma*3/n0[1]
      hh[1,6]<-sigma*3/n0[1]
      hh[1,7]<-sigma/n0[3]
      hh[1,8]<-sigma/n0[3]
      hh[1,9]<-16*sigma5[5]/n0[12]
      hh[2,2]<-sigma/n0[1]+4*sigma4[2]/n0[4]+sigma4[3]/n0[5]+16*sigma5[2]/n0[10]
      hh[2,3]<-sigma/n0[1]
      hh[2,4]<-sigma/n0[1]-sigma4[3]/n0[5]
      hh[2,5]<-3*sigma/n0[1]-10*sigma4[2]/n0[4]+3*sigma4[3]/n0[5]+4*sigma5[2]/n0[10];
      hh[2,6]<-3*sigma/n0[1]
      hh[2,7]<-sigma4[3]/n0[5]
      hh[2,8]<-0
      hh[2,9]<-2*sigma4[2]/n0[4]
      hh[3,3]<-sigma/n0[1]+4*sigma4[4]/n0[6]+sigma4[7]/n0[8]+16*sigma5[4]/n0[11]
      hh[3,4]<-sigma/n0[1]-sigma4[7]/n0[8]
      hh[3,5]<-sigma*3/n0[1]
      hh[3,6]<-sigma*3/n0[1]-10*sigma4[4]/n0[6]+3*sigma4[7]/n0[8]+4*sigma5[4]/n0[11]
      hh[3,7]<-0
      hh[3,8]<-sigma4[7]/n0[8]
      hh[3,9]<-2*sigma4[4]/n0[6]
      hh[4,4]<-sigma*(1/n0[1]+1/n0[3])+sigma4[3]/n0[5]+sigma4[7]/n0[8]
      hh[4,5]<-sigma*3/n0[1]-3*sigma4[3]/n0[5]
      hh[4,6]<-sigma*3/n0[1]-3*sigma4[7]/n0[8]
      hh[4,7]<-sigma/n0[3]-sigma4[3]/n0[5]
      hh[4,8]<-sigma/n0[3]-sigma4[7]/n0[8]
      hh[4,9]<-0
      hh[5,5]<-sigma*9/n0[1]+25*sigma4[2]/n0[4]+9*sigma4[3]/n0[5]+121*sigma4[8]/n0[9]+sigma5[2]/n0[10]+121*sigma5[8]/n0[14]
      hh[5,6]<-sigma*9/n0[1]
      hh[5,7]<-sigma4[3]*3/n0[5]
      hh[5,8]<-sigma4[8]*22/n0[9]+44*sigma5[8]/n0[14]
      hh[5,9]<--sigma4[2]*5/n0[4]+11*sigma4[8]/n0[9]
      hh[6,6]<-sigma*9/n0[1]+25*sigma4[4]/n0[6]+121*sigma4[6]/n0[7]+9*sigma4[7]/n0[8]+sigma5[4]/n0[11]+121*sigma5[6]/n0[13]
      hh[6,7]<-22*sigma4[6]/n0[7]+44*sigma5[6]/n0[13]
      hh[6,8]<-3*sigma4[7]/n0[8]
      hh[6,9]<--5*sigma4[4]/n0[6]+11*sigma4[6]/n0[7]
      hh[7,7]<-sigma/n0[3]+sigma4[3]/n0[5]+4*sigma4[6]/n0[7]+16*sigma5[6]/n0[13]
      hh[7,8]<-sigma/n0[3]
      hh[7,9]<-2*sigma4[6]/n0[7]
      hh[8,8]<-sigma/n0[3]+sigma4[7]/n0[8]+4*sigma4[8]/n0[9]+16*sigma5[8]/n0[14]
      hh[8,9]<-2*sigma4[8]/n0[9]
      hh[9,9]<-sigma4[2]/n0[4]+sigma4[4]/n0[6]+sigma4[6]/n0[7]+sigma4[8]/n0[9]+16*sigma5[5]/n0[12]
      for(i in 2:9)
      {
        for(j in 1:(i-1))
        {
          hh[i,j]<-hh[j,i]
        }
      }
      #################################################################### 
      b_line[1]<-s0[1]/n0[1]+2*s0[2]/n0[2]+s0[3]/n0[3]-4*s0[12]/n0[12]
      b_line[2]<-s0[1]/n0[1]+2*s0[4]/n0[4]+s0[5]/n0[5]-4*s0[10]/n0[10]
      b_line[3]<-s0[1]/n0[1]+2*s0[6]/n0[6]+s0[8]/n0[8]-4*s0[11]/n0[11]
      b_line[4]<-s0[1]/n0[1]+s0[3]/n0[3]-s0[5]/n0[5]-s0[8]/n0[8]
      b_line[5]<-3*s0[1]/n0[1]-5*s0[4]/n0[4]+3*s0[5]/n0[5]+11*s0[9]/n0[9]-s0[10]/n0[10]-11*s0[14]/n0[14]
      b_line[6]<-3*s0[1]/n0[1]-5*s0[6]/n0[6]+11*s0[7]/n0[7]+3*s0[8]/n0[8]-s0[11]/n0[11]-11*s0[13]/n0[13]
      b_line[7]<-s0[3]/n0[3]+s0[5]/n0[5]+2*s0[7]/n0[7]-4*s0[13]/n0[13]
      b_line[8]<-s0[3]/n0[3]+s0[8]/n0[8]+2*s0[9]/n0[9]-4*s0[14]/n0[14]
      b_line[9]<-s0[4]/n0[4]+s0[6]/n0[6]+s0[7]/n0[7]+s0[9]/n0[9]-4*s0[12]/n0[12]
      B611<-solve(hh,b_line)
      ################################################################
      m[1]<-(s0[1]-sigma*(B611[1]+B611[2]+B611[3]+B611[4]+3*B611[5]+3*B611[6]))/n0[1]
      m[2]<-(s0[2]-sigma*2*B611[1])/n0[2]
      m[3]<-(s0[3]-sigma*(B611[1]+B611[4]+B611[7]+B611[8]))/n0[3]
      m4[2]<-(s0[4]-sigma4[2]*(2*B611[2]-5.0*B611[5]+B611[9]))/n0[4]
      m4[3]<-(s0[5]-sigma4[3]*(B611[2]-B611[4]+3*B611[5]+B611[7]))/n0[5]
      m4[4]<-(s0[6]-sigma4[4]*(2*B611[3]-5*B611[6]+B611[9]))/n0[6]
      m4[6]<-(s0[7]-sigma*(11*B611[6]+2*B611[7]+B611[9]))/n0[7]
      m4[7]<-(s0[8]-sigma*(B611[3]-B611[4]+3*B611[6]+B611[8]))/n0[8]
      m4[8]<-(s0[9]-sigma*(11*B611[5]+2*B611[8]+B611[9]))/n0[9]
      m5[2]<-(s0[10]+sigma5[2]*(4*B611[2]+B611[5]))/n0[10]
      m5[4]<-(s0[11]+sigma5[4]*(4*B611[3]+B611[6]))/n0[11]
      m5[5]<-(s0[12]+sigma5[5]*(4*B611[1]+4*B611[9]))/n0[12]
      m5[6]<-(s0[13]+sigma5[6]*(11*B611[6]+4*B611[7]))/n0[13]
      m5[8]<-(s0[14]+sigma5[8]*(11*B611[5]+4*B611[8]))/n0[14]
      m4[1]<-m[1];m4[5]<-m[2];m4[9]<-m[3]
      m5[1]<-m[1];m5[3]<-m4[3];m5[7]<-m4[7];m5[9]<-m[3]
      aaa1<-max(abs(B611-AA))
      AA<-B611
      if (n_iter>20) break
    }  
    ##########obtain variance###############
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d3) {swx24[i] <- W4[i,]%*%(dataF2-m4[i])^2 } ;for(i in 1:d3) {swx25[i] <- W5[i,]%*%(dataF3-m5[i])^2 }      
    s0[15]<-sum(swx24)+ss1+ss2+ss3+(swx25[1]+swx25[3]+swx25[7]+swx25[9])*m_nf
    n0[15]<-nn[1]+nn[2]+nn[3]+nn[4]+nn[5]*(mi5[1]+mi5[3]+mi5[7]+mi5[9])
    mm<-as.matrix(c(m[1],m[2],m[3],m4[2],m4[3],m4[4],m4[6],m4[7],m4[8],m5[2],m5[4],m5[5],m5[6],m5[8]))
    B6111<-solve(crossprod(hh6,hh6))%*%crossprod(hh6,mm)
    gs<-B6111[c(2:5)]
    g_aa1<-0.5*gs[2]^2+0.25*gs[4]^2    #   0.5*db**2+0.25*hb**2.
    g_aa2<-0.5*gs[1]^2+0.25*gs[3]^2    #   0.5*da**2+0.25*ha**2.
    g_aa3<-g_aa1+g_aa2     #   0.5(da**2+db**2)+0.25(ha**2+hb**2).
    aaa0<-sigma;n_iter<-0;aaa6<-1000
    while (aaa6>0.0001){
      n_iter<-n_iter+1
      aa1<-sigma/(sigma+g_aa1)
      aa2<-sigma/(sigma+g_aa2)
      aa3<-sigma/(sigma+g_aa3)
      aa4<-((swx25[2]+swx25[8])*aa1^2+(swx25[4]+swx25[6])*aa2^2+swx25[5]*aa3^2)*m_nf
      aa5<-(mi5[2]+mi5[8])*nn[5]*aa1+(mi5[4]+mi5[6])*nn[5]*aa2+mi5[5]*nn[5]*aa3
      sigma<-(s0[15]+aa4)/(n0[15]+aa5)
      aa6<-abs(sigma-aaa0)
      aaa0<-sigma
      if (n_iter>20) break
    } 
    if (sigma<sigma0) {sigma<-sigma0}
    sigma4[c(1:9)]<-sigma;sigma5[c(1,3,7,9)]<-sigma/m_nf
    sigma5[2]<-(sigma+g_aa1)/m_nf;sigma5[4]<-(sigma+g_aa2)/m_nf
    sigma5[5]<-(sigma+g_aa3)/m_nf;sigma5[6]<-sigma5[4];sigma5[8]<-sigma5[2]
    ##############################criteria for iterations to stop######################## 
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
      sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*6
  #######first order parameters############
  mm<-as.matrix(c(m[1],m[2],m[3],m4[2],m4[3],m4[4],m4[6],m4[7],m4[8],m5[2],m5[4],m5[5],m5[6],m5[8]))
  B61111<-solve(crossprod(hh6,hh6))%*%crossprod(hh6,mm)
  ########second order genetic parameters###############   
  F2jj<-sigmaF2-sigma
  if(F2jj<0){F2jj<-0}
  F2ll<-F2jj/sigmaF2
  F3jj<-sigmaF3-sigma/m_nf
  if(F3jj<0) {F3jj<-0}
  F3ll<-F3jj/sigmaF3   
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma))
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
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma))
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
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma))
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
    F2gg <- (dataF2 - m4[i])/sqrt(sigma4[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mi4[i]
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
  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3); 
  F3w1<-1/(12*m_sam[5])
  F3bmw <- matrix(0,m_sam[5],1); F3bmwsl <- matrix(0,m_sam[5],d3)
  for(i in 1:d3){
    F3gg <- (dataF3 - m5[i])/sqrt(sigma5[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mi5[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[5]){F3P2 <- F3P2+runif(m_sam[5])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[5]) + sum((F3P2 - (as.matrix(c(1:m_sam[5])) - 0.5)/m_sam[5])^2)
  F3u <- as.matrix(c(12*m_sam[5]*((F3dd[1]/m_sam[5]-0.5)^2),((45*m_sam[5])/4)*((F3dd[2]/m_sam[5]-1/3)^2),180*m_sam[5]*((F3dd[3]/m_sam[5]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[5]))))
  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("2MG-AD",round(abc,4),round(AIC,4),round(m[1],4),round(m[2],4),round(m[3],4), round(t(m4),4),round(sigma4[1],4),
                       round(t(mi4),4),round(t(m5),4),round(t(sigma5),4),round(t(mi5),4),round(sigma,4),          
                       round(B61111[1],4)," "," "," "," ",round(B61111[2],4),round(B61111[3],4),round(B61111[4],4),round(B61111[5],4)," "," "," "," "," "," ",round(F2jj,4),round(F2ll*100,4)," "," ",round(F3jj,4),round(F3ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
 output<-as.matrix(output)
  OUTPUT<-list(output,mix_pi4,mix_pi5)
  return(OUTPUT)
}

##################2MG-A(B-3)############################
G5ModelFun[[7]] <- function(K1,logL,df11,df21,df31,df41,df51,G5text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataF2 <- as.matrix(as.numeric(df41[,1]))
  dataF3 <- as.matrix(as.numeric(df51[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2),sum(dataF3)))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2),mean(dataF3)))
  sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1);sigmaF2<-var(dataF2);sigmaF3<-var(dataF3)
  sumx6<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  sigma<-sumx6;sigma0<-sigma;sigma40<-sigmaF2;sigma400<-sigma40;sigma50<-sigmaF3;sigma500<-sigma50
  m_nf <- as.numeric(G5text2)
  ###############procedure start###########################
  mix_pi4<-as.matrix(c(0.0625,0.125,0.0625,0.1250,0.250,0.125,0.0625,0.125,0.0625))
  mix_pi5<-as.matrix(c(0.0625,0.125,0.0625,0.1250,0.250,0.125,0.0625,0.125,0.0625))
  a1<-sqrt(sigma40/m_sam[4])
  if (m[1]<m[3]) a1<--a1
  m4<-as.matrix(c(m[4]+2.8*a1,m[4]+2.1*a1,m[4]+1.4*a1,m[4]+0.7*a1,m[4],m[4]-0.7*a1,m[4]-1.4*a1,m[4]-2.1*a1,m[4]-2.8*a1))
  a1<-sqrt(sigma50/m_sam[5])
  if (m[1]<m[3]) {a1<--a1}
  m5<-as.matrix(c(m[5]+2.8*a1,m[5]+2.1*a1,m[5]+1.4*a1,m[5]+0.7*a1,m[5],m[5]-0.7*a1,m[5]-1.4*a1,m[5]-2.1*a1,m[5]-2.8*a1))
  sigma<-sigma400/2
  #######first order parameters############
  hh7<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,-1,1,1,0,0,-1,-1,1,0,0,0,-1,
                1,0,-1,1,-1,1,-1,1,0,0,1,0,-1,0),14,3)
  mm<-as.matrix(c(m[1],m[2],m[3],m4[2],m4[3],m4[4],m4[6],m4[7],m4[8],m5[2],m5[4],m5[5],m5[6],m5[8]))
  B7<-solve(crossprod(hh7,hh7))%*%crossprod(hh7,mm)
  gs<-B7[c(2:3)]
  sigma4<-matrix(0,9,1);sigma5<-matrix(0,9,1)
  sigma4[c(1:9)]<-sigma;sigma5[c(1,3,7,9)]<-sigma/m_nf
  g_aa1<-0.5*gs[2]^2   # 0.5*db**2.
  g_aa2<-0.5*gs[1]^2   # 0.5*da**2.
  g_aa3<-g_aa1+g_aa2  #  0.5(da**2+db**2).
  sigma5[2]<-(sigma+g_aa1)/m_nf;sigma5[4]<-(sigma+g_aa2)/m_nf
  sigma5[5]<-(sigma+g_aa3)/m_nf;sigma5[6]<-sigma5[4]
  sigma5[8]<-sigma5[2]
  mi4<-mix_pi4[c(1:9)];mi5<-mix_pi5[c(1:9)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
    sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
  ########################################################
  ##########iteration process###########
  ############E-step###################
  ############F2 process###############
  d3<-9
  iteration <- 0; stopa <- 1000
  W4 <- matrix(0,9,m_sam[4]); swx24 <- matrix(0,9,1)
  W5 <- matrix(0,9,m_sam[5]); swx25 <- matrix(0,9,1)
  hh<-matrix(0,6,6);b_line<-matrix(0,6,1)
  s0<-matrix(0,20,1);n0<-matrix(0,20,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:d3) { W4[i,] <- mix_pi4[i]*dnorm(dataF2,m4[i],sqrt(sigma4[i]))/dmixnorm(dataF2,m4,sqrt(sigma4),mix_pi4) }
    mi4 <- as.matrix(rowSums(W4)/m_sam[4])
    sumwx4 <- W4%*%dataF2
    #############F3 process############
    for(i in 1:d3) { W5[i,] <- mix_pi5[i]*dnorm(dataF3,m5[i],sqrt(sigma5[i]))/dmixnorm(dataF3,m5,sqrt(sigma5),mix_pi5) }
    mi5 <- as.matrix(rowSums(W5)/m_sam[5])
    sumwx5 <- W5%*%dataF3
    aaa0<-0;nn<-m_sam
    s0[1]<-sumx[1]+sumwx4[1]+m_nf*sumwx5[1];s0[2]<-sumx[2]+sumwx4[5]+sigma*sumwx5[5]/sigma5[5]
    s0[3]<-sumx[3]+sumwx4[9]+m_nf*sumwx5[9];s0[4]<-sumwx4[2]+sigma*sumwx5[2]/sigma5[2]
    s0[5]<-sumwx4[3]+m_nf*sumwx5[3];s0[6]<-sumwx4[4]+sumwx5[4]*sigma/sigma5[4]
    s0[7]<-sumwx4[6]+sumwx5[6]*sigma/sigma5[6];s0[8]<-sumwx4[7]+m_nf*sumwx5[7]
    s0[9]<-sumwx4[8]+sumwx5[8]*sigma/sigma5[8]
    n0[1]<-nn[1]+mi4[1]*nn[4]+m_nf*mi5[1]*nn[5];n0[2]<-nn[2]+mi4[5]*nn[4]+mi5[5]*nn[5]*sigma/sigma5[5]
    n0[3]<-nn[3]+mi4[9]*nn[4]+m_nf*mi5[9]*nn[5];n0[4]<-mi4[2]*nn[4]+mi5[2]*nn[5]*sigma/sigma5[2]
    n0[5]<-mi4[3]*nn[4]+m_nf*mi5[3]*nn[5];n0[6]<-mi4[4]*nn[4]+mi5[4]*nn[5]*sigma/sigma5[4]
    n0[7]<-mi4[6]*nn[4]+mi5[6]*nn[5]*sigma/sigma5[6];n0[8]<-mi4[7]*nn[4]+m_nf*mi5[7]*nn[5]
    n0[9]<-mi4[8]*nn[4]+mi5[8]*nn[5]*sigma/sigma5[8]
    n0[c(1:9)][abs(n0[c(1:9)])<0.000001] <- 0.000001
    aa1<-1000; n_iter<-0;AA<-matrix(0,6,1)
    while (aa1>0.0001)
    {
      n_iter<-n_iter+1
      ##########first order parameters###################
      mm<-as.matrix(c(m[1],m[2],m[3],m4[2],m4[3],m4[4],m4[6],m4[7],m4[8],m5[2],m5[4],m5[5],m5[6],m5[8]))
      B71<-solve(crossprod(hh7,hh7))%*%crossprod(hh7,mm)
      gs<-B71[c(2:3)]
      sigma4[c(1:9)]<-sigma;sigma5[c(1,3,7,9)]<-sigma/m_nf
      g_aa1<-0.5*gs[2]^2   #   0.5*db**2.
      g_aa2<-0.5*gs[1]^2   #   0.5*da**2.
      g_aa3<-g_aa1+g_aa2       #   0.5(da**2+db**2).
      sigma5[2]<-(sigma+g_aa1)/m_nf;sigma5[4]<-(sigma+g_aa2)/m_nf
      sigma5[5]<-(sigma+g_aa3)/m_nf;sigma5[6]<-sigma5[4];sigma5[8]<-sigma5[2]
      #################restrictions########################################
      hh[1,1]<-sigma*(1/n0[1]+1/n0[5]+4/n0[4])
      hh[1,2]<-sigma/n0[1]
      hh[1,3]<-sigma*2/n0[4]
      hh[1,4]<--sigma/n0[5]
      hh[1,5]<-sigma/n0[5];hh[1,6]<-0
      hh[2,2]<-sigma*(1/n0[1]+4/n0[6]+1/n0[8]);hh[2,3]<-0
      hh[2,4]<--sigma/n0[8];hh[2,5]<-0
      hh[2,6]<-sigma/n0[8]
      hh[3,3]<-sigma*(4/n0[2]+1/n0[4]+1/n0[9])
      hh[3,4]<-4*sigma/n0[2];hh[3,5]<-0
      hh[3,6]<-sigma*2/n0[9]
      hh[4,4]<-sigma*(4/n0[2]+1/n0[5]+1/n0[8])
      hh[4,5]<--sigma/n0[5]
      hh[4,6]<--sigma/n0[8]
      hh[5,5]<-sigma*(1/n0[3]+1/n0[5]+4/n0[7])
      hh[5,6]<-sigma/n0[3]
      hh[6,6]<-sigma*(1/n0[3]+1/n0[8]+4/n0[9])
      for(i in 2:6)
      {
        for(j in 1:(i-1))
        {
          hh[i,j]<-hh[j,i]
        }
      }
      ##########################################################
      b_line[1]<-s0[1]/n0[1]+s0[5]/n0[5]-2*s0[4]/n0[4]
      b_line[2]<-s0[1]/n0[1]+s0[8]/n0[8]-2*s0[6]/n0[6]
      b_line[3]<-2*s0[2]/n0[2]-s0[4]/n0[4]-s0[9]/n0[9]
      b_line[4]<-2*s0[2]/n0[2]-s0[5]/n0[5]-s0[8]/n0[8]
      b_line[5]<-s0[3]/n0[3]+s0[5]/n0[5]-2*s0[7]/n0[7]
      b_line[6]<-s0[3]/n0[3]+s0[8]/n0[8]-2*s0[9]/n0[9]
      B711<-solve(hh,b_line)
      ############################################################
      m[1]<-(s0[1]-sigma*(B711[1]+B711[2]))/n0[1]
      m[2]<-(s0[2]-sigma*(2*B711[3]+2*B711[4]))/n0[2]
      m[3]<-(s0[3]-sigma*(B711[5]+B711[6]))/n0[3]
      m4[2]<-(s0[4]+(2*B711[1]+B711[3])*sigma)/n0[4]
      m4[3]<-(s0[5]-(B711[1]-B711[4]+B711[5])*sigma)/n0[5]
      m4[4]<-(s0[6]+sigma*2*B711[2])/n0[6]
      m4[6]<-(s0[7]+sigma*2*B711[5])/n0[7]
      m4[7]<-(s0[8]+sigma*(-B711[2]+B711[4]-B711[6]))/n0[8]
      m4[8]<-(s0[9]+sigma*(B711[3]+2*B711[6]))/n0[9]
      m4[1]<-m[1];m4[5]<-m[2];m4[9]<-m[3]
      m5[1]<-m[1];m5[2]<-m4[2];m5[3]<-m4[3]
      m5[4]<-m4[4];m5[5]<-m4[5];m5[6]<-m4[6]
      m5[7]<-m4[7];m5[8]<-m4[8];m5[9]<-m[3]
      aa1<-max(abs(B711-AA))
      AA<-B711
      if (n_iter>20) break
    }  
    ##########obtain variance###############
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d3) {swx24[i] <- W4[i,]%*%(dataF2-m4[i])^2 } ;for(i in 1:d3) {swx25[i] <- W5[i,]%*%(dataF3-m5[i])^2 }      
    s0[11]<-sum(swx24)+ss1+ss2+ss3+(swx25[1]+swx25[3]+swx25[7]+swx25[9])*m_nf
    n0[11]<-nn[1]+nn[2]+nn[3]+nn[4]+(mi5[1]+mi5[3]+mi5[7]+mi5[9])*nn[5]
    n0[12]<-mi5[2]*nn[5];n0[13]<-mi5[4]*nn[5];n0[14]<-mi5[5]*nn[5]
    n0[15]<-mi5[6]*nn[5];n0[16]<-mi5[8]*nn[5]
    mm<-as.matrix(c(m[1],m[2],m[3],m4[2],m4[3],m4[4],m4[6],m4[7],m4[8],m5[2],m5[4],m5[5],m5[6],m5[8]))
    B7111<-solve(crossprod(hh7,hh7))%*%crossprod(hh7,mm)
    gs<-B7111[c(2:3)]
    g_aa1<-0.5*gs[2]^2   #   0.5*db**2.
    g_aa2<-0.5*gs[1]^2   #  0.5*da**2.
    g_aa3<-g_aa1+g_aa2       #  0.5(da**2+db**2).
    aaa0<-sigma;n_iter<-0;aa4<-1000
    while (aa4>0.0001){
      n_iter<-n_iter+1
      aa1<-sigma/(sigma+g_aa1)
      aa2<-sigma/(sigma+g_aa2)
      aa3<-sigma/(sigma+g_aa3)
      sigma<-(s0[11]+aa1^2*(swx25[2]+swx25[8])+aa2^2*(swx25[4]+swx25[6])+aa3^2*swx25[5])/(n0[11]+aa1*(n0[12]+n0[16])+aa2*(n0[13]+n0[15])+aa3*n0[14])
      aa4<-abs(sigma-aaa0)
      aaa0<-sigma
      if (n_iter>20) break
    } 
    if (sigma<sigma0) {sigma<-sigma0}
    sigma4[c(1:9)]<-sigma;sigma5[c(1,3,7,9)]<-sigma/m_nf
    sigma5[2]<-(sigma+g_aa1)/m_nf;sigma5[4]<-(sigma+g_aa2)/m_nf
    sigma5[5]<-(sigma+g_aa3)/m_nf;sigma5[6]<-sigma5[4];sigma5[8]<-sigma5[2]
    ##############################criteria for iterations to stop######################## 
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
      sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*4
  #######first order parameters############
  mm<-as.matrix(c(m[1],m[2],m[3],m4[2],m4[3],m4[4],m4[6],m4[7],m4[8],m5[2],m5[4],m5[5],m5[6],m5[8]))
  B71111<-solve(crossprod(hh7,hh7))%*%crossprod(hh7,mm)
  ########second order genetic parameters###############   
  F2jj<-sigmaF2-sigma
  if(F2jj<0){F2jj<-0}
  F2ll<-F2jj/sigmaF2
  F3jj<-sigmaF3-sigma/m_nf
  if(F3jj<0) {F3jj<-0}
  F3ll<-F3jj/sigmaF3   
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma))
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
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma))
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
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma))
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
    F2gg <- (dataF2 - m4[i])/sqrt(sigma4[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mi4[i]
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
  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3); 
  F3w1<-1/(12*m_sam[5])
  F3bmw <- matrix(0,m_sam[5],1); F3bmwsl <- matrix(0,m_sam[5],d3)
  for(i in 1:d3){
    F3gg <- (dataF3 - m5[i])/sqrt(sigma5[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mi5[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[5]){F3P2 <- F3P2+runif(m_sam[5])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[5]) + sum((F3P2 - (as.matrix(c(1:m_sam[5])) - 0.5)/m_sam[5])^2)
  F3u <- as.matrix(c(12*m_sam[5]*((F3dd[1]/m_sam[5]-0.5)^2),((45*m_sam[5])/4)*((F3dd[2]/m_sam[5]-1/3)^2),180*m_sam[5]*((F3dd[3]/m_sam[5]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[5]))))
  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("2MG-A",round(abc,4),round(AIC,4),round(m[1],4),round(m[2],4),round(m[3],4), round(t(m4),4),round(sigma4[1],4),
                       round(t(mi4),4),round(t(m5),4),round(t(sigma5),4),round(t(mi5),4),round(sigma,4),          
                       round(B71111[1],4)," "," "," "," ",round(B71111[2],4),round(B71111[3],4)," "," "," "," "," "," "," "," ",round(F2jj,4),round(F2ll*100,4)," "," ",round(F3jj,4),round(F3ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
 output<-as.matrix(output)
  OUTPUT<-list(output,mix_pi4,mix_pi5)
  return(OUTPUT)
}
##############2MG-EA(B-4)##############################
G5ModelFun[[8]] <- function(K1,logL,df11,df21,df31,df41,df51,G5text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataF2 <- as.matrix(as.numeric(df41[,1]))
  dataF3 <- as.matrix(as.numeric(df51[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2),sum(dataF3)))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2),mean(dataF3)))
  sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1);sigmaF2<-var(dataF2);sigmaF3<-var(dataF3)
  sumx6<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  sigma<-sumx6;sigma0<-sigma;sigma40<-sigmaF2;sigma400<-sigma40;sigma50<-sigmaF3;sigma500<-sigma50
  m_nf <- as.numeric(G5text2)
  ###############procedure start###########################
  sigma4<-matrix(0,5,1);sigma5<-matrix(0,6,1)
  sigma4[c(1:5)]<-sigma;sigma5[c(1,3)]<-sigma/m_nf
  mix_pi4<-as.matrix(c(0.0625,0.25,0.375,0.25,0.0625))
  mix_pi5<-as.matrix(c(0.0625,0.25,0.125,0.25,0.25,0.0625))
  m4<-as.matrix(c(m[1],(m[1]+m[2])/2,m[2],(m[2]+m[3])/2,m[3]))
  m5<-as.matrix(c(m[1],m4[2],m[2],m[2],m4[4],m[3]))
  #######first order parameters############
  hh8<-matrix(c(1,1,1,1,1,2,0,-2,1,-1),5,2)
  mm<-as.matrix(c(m[1],m[2],m[3],m4[2],m4[4]))
  B8<-solve(crossprod(hh8,hh8))%*%crossprod(hh8,mm)
  a1<-B8[2];a1<-a1^2
  sigma5[2]<-(sigma+0.5*a1)/m_nf;sigma5[4]<-(sigma+a1)/m_nf
  sigma5[5]<-sigma5[2];sigma5[6]<-sigma/m_nf
  mi4<-mix_pi4[c(1:5)];mi5<-mix_pi5[c(1:6)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
    sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
  ########################################################
  ##########iteration process###########
  ############E-step###################
  ############F2 process###############
  d5<-5;d6<-6
  iteration <- 0; stopa <- 1000
  W4 <- matrix(0,5,m_sam[4]); swx24 <- matrix(0,5,1)
  W5 <- matrix(0,6,m_sam[5]); swx25 <- matrix(0,6,1)
  hh<-matrix(0,3,3);b_line<-matrix(0,3,1)
  s0<-matrix(0,20,1);n0<-matrix(0,20,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:d5) { W4[i,] <- mix_pi4[i]*dnorm(dataF2,m4[i],sqrt(sigma4[i]))/dmixnorm(dataF2,m4,sqrt(sigma4),mix_pi4) }
    mi4 <- as.matrix(rowSums(W4)/m_sam[4])
    sumwx4 <- W4%*%dataF2
    #############F3 process############
    for(i in 1:d6) { W5[i,] <- mix_pi5[i]*dnorm(dataF3,m5[i],sqrt(sigma5[i]))/dmixnorm(dataF3,m5,sqrt(sigma5),mix_pi5) }
    mi5 <- as.matrix(rowSums(W5)/m_sam[5])
    sumwx5 <- W5%*%dataF3
    aaa0<-0;nn<-m_sam
    s0[1]<-sumx[1]+sumwx4[1]+m_nf*sumwx5[1];s0[2]<-sumx[2]+sumwx4[3]+m_nf*sumwx5[3]+sumwx5[4]*sigma/sigma5[4]
    s0[3]<-sumx[3]+sumwx4[5]+m_nf*sumwx5[6];s0[4]<-sumwx4[2]+sumwx5[2]*sigma/sigma5[2]
    s0[5]<-sumwx4[4]+sumwx5[5]*sigma/sigma5[5]
    n0[1]<-nn[1]+mi4[1]*nn[4]+mi5[1]*nn[5]*m_nf;n0[2]<-nn[2]+mi4[3]*nn[4]+mi5[3]*nn[5]*m_nf+mi5[4]*nn[5]*sigma/sigma5[4]
    n0[3]<-nn[3]+mi4[5]*nn[4]+mi5[6]*nn[5]*m_nf;n0[4]<-mi4[2]*nn[4]+mi5[2]*nn[5]*sigma/sigma5[2]
    n0[5]<-mi4[4]*nn[4]+mi5[5]*nn[5]*sigma/sigma5[5]
    n0[c(1:5)][abs(n0[c(1:5)])<0.000001]<-0.000001
    aaa1<-1000;n_iter<-0;AA<-matrix(0,3,1)
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      ##########first order parameters###################
      mm<-as.matrix(c(m[1],m[2],m[3],m4[2],m4[4]))
      B81<-solve(crossprod(hh8,hh8))%*%crossprod(hh8,mm)
      aa6<-B81[2];aa6<-aa6^2
      sigma5[2]<-(sigma+0.5*aa6)/m_nf;sigma5[4]<-(sigma+aa6)/m_nf;sigma5[5]<-sigma5[2]
      #################restrictions########################################
      hh[1,1]<-sigma/n0[1]+sigma/n0[2]+4*sigma/n0[4]
      hh[1,2]<-sigma*(2/n0[2]+2/n0[4])
      hh[1,3]<-sigma/n0[2]
      hh[2,2]<-sigma*(4/n0[2]+1/n0[4]+1/n0[5])
      hh[2,3]<-sigma*(2/n0[2]+2/n0[5])
      hh[3,3]<-sigma*(1/n0[3]+1/n0[2]+4/n0[5])
      for(i in 2:3)
      {
        for(j in 1:(i-1))
        {
          hh[i,j]<-hh[j,i]
        }
      }
      #################################################################### 
      b_line[1]<-s0[1]/n0[1]+s0[2]/n0[2]-2*s0[4]/n0[4]
      b_line[2]<-2*s0[2]/n0[2]-s0[4]/n0[4]-s0[5]/n0[5]
      b_line[3]<-s0[2]/n0[2]+s0[3]/n0[3]-2*s0[5]/n0[5]
      B811<-solve(hh,b_line)
      m[1]<-(s0[1]-sigma*B811[1])/n0[1]
      m[2]<-(s0[2]-sigma*(B811[1]+2*B811[2]+B811[3]))/n0[2]
      m[3]<-(s0[3]-sigma*B811[3])/n0[3]
      m4[2]<-(s0[4]+sigma*(2*B811[1]+B811[2]))/n0[4]
      m4[4]<-(s0[5]+sigma*(B811[2]+2*B811[3]))/n0[5]
      aaa1<-max(abs(B811-AA))
      AA<-B811
      if (n_iter>20) break
    }  
    m4[1]<-m[1];m4[3]<-m[2];m4[5]<-m[3]
    m5[1]<-m[1];m5[2]<-m4[2];m5[3]<-m[2]
    m5[4]<-m[2];m5[5]<-m4[4];m5[6]<-m[3]
    #######first order parameters############
    mm<-as.matrix(c(m[1],m[2],m[3],m4[2],m4[4]))
    B8111<-solve(crossprod(hh8,hh8))%*%crossprod(hh8,mm)
    aa6<-B8111[2]
    aa6<-aa6^2
    ##########obtain variance###############
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d5) {swx24[i] <- W4[i,]%*%(dataF2-m4[i])^2 } ;for(i in 1:d6) {swx25[i] <- W5[i,]%*%(dataF3-m5[i])^2 }      
    s0[6]<-sum(swx24)+ss1+ss2+ss3+(swx25[1]+swx25[3]+swx25[6])*m_nf
    n0[6]<-nn[1]+nn[2]+nn[3]+nn[4]+(mi5[1]+mi5[3]+mi5[6])*nn[5]
    n0[7]<-(mi5[2]+mi5[5])*nn[5]
    n0[8]<-mi5[4]*nn[5]
    aaa0<-sigma;n_iter<-0;aaa1<-1000
    while (aaa1>0.0001){
      n_iter<-n_iter+1
      aa1<-sigma/(sigma+0.5*aa6)
      aa2<-sigma/(sigma+aa6)
      sigma<-(s0[6]+m_nf*(aa1^2*(swx25[2]+swx25[5])+aa2^2*swx25[4]))/(n0[6]+aa1*n0[7]+aa2*n0[8])
      aaa1<-abs(sigma-aaa0)
      aaa0<-sigma
      if (n_iter>20) break
    } 
    if (sigma<sigma0) {sigma<-sigma0}
    sigma4[c(1:5)]<-sigma
    #######first order parameters############
    sigma5[1]<-sigma/m_nf;sigma5[2]<-(sigma+0.5*aa6)/m_nf
    sigma5[3]<-sigma/m_nf;sigma5[4]<-(sigma+aa6)/m_nf
    sigma5[5]<-sigma5[2];sigma5[6]<-sigma/m_nf
    ##############################criteria for iterations to stop######################## 
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
      sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*3
  #######first order parameters############
  mm<-as.matrix(c(m[1],m[2],m[3],m4[2],m4[4]))
  B81111<-solve(crossprod(hh8,hh8))%*%crossprod(hh8,mm)
  ########second order genetic parameters###############   
  F2jj<-sigmaF2-sigma
  if(F2jj<0){F2jj<-0}
  F2ll<-F2jj/sigmaF2
  F3jj<-sigmaF3-sigma/m_nf
  if(F3jj<0) {F3jj<-0}
  F3ll<-F3jj/sigmaF3   
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma))
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
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma))
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
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma))
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
    F2gg <- (dataF2 - m4[i])/sqrt(sigma4[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mi4[i]
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
  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3); 
  F3w1<-1/(12*m_sam[5])
  F3bmw <- matrix(0,m_sam[5],1); F3bmwsl <- matrix(0,m_sam[5],d6)
  for(i in 1:d6){
    F3gg <- (dataF3 - m5[i])/sqrt(sigma5[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mi5[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[5]){F3P2 <- F3P2+runif(m_sam[5])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[5]) + sum((F3P2 - (as.matrix(c(1:m_sam[5])) - 0.5)/m_sam[5])^2)
  F3u <- as.matrix(c(12*m_sam[5]*((F3dd[1]/m_sam[5]-0.5)^2),((45*m_sam[5])/4)*((F3dd[2]/m_sam[5]-1/3)^2),180*m_sam[5]*((F3dd[3]/m_sam[5]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[5]))))
  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("2MG-EA",round(abc,4),round(AIC,4),round(m[1],4),round(m[2],4),round(m[3],4), round(t(m4),4)," "," "," "," ",round(sigma4[1],4),
                       round(t(mi4),4)," "," "," "," ",round(t(m5),4)," "," "," ",round(t(sigma5),4)," "," "," ",round(t(mi5),4)," "," "," ",round(sigma,4),          
                       round(B81111[1],4)," "," "," "," ",round(B81111[2],4)," "," "," "," "," "," "," "," "," ",round(F2jj,4),round(F2ll*100,4)," "," ",round(F3jj,4),round(F3ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
 output<-as.matrix(output)
  OUTPUT<-list(output,mix_pi4,mix_pi5)
  return(OUTPUT)
}
###############2MG-CD(B-5)##########################
G5ModelFun[[9]] <- function(K1,logL,df11,df21,df31,df41,df51,G5text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataF2 <- as.matrix(as.numeric(df41[,1]))
  dataF3 <- as.matrix(as.numeric(df51[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2),sum(dataF3)))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2),mean(dataF3)))
  sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1);sigmaF2<-var(dataF2);sigmaF3<-var(dataF3)
  sumx6<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  sigma<-sumx6;sigma0<-sigma;sigma40<-sigmaF2;sigma400<-sigma40;sigma50<-sigmaF3;sigma500<-sigma50
  m_nf <- as.numeric(G5text2)
  ###############procedure start##########################
  sigma4<-matrix(0,4,1);sigma5<-matrix(0,9,1)
  sigma4[c(1:4)]<-sigma;sigma5[c(1,3,6,7,9)]<-sigma/m_nf
  mix_pi4<-as.matrix(c(0.5625,0.1875,0.1875,0.0625))
  mix_pi5<-as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  nn<-m_sam
  a1<-sqrt(sigma40/(nn[4]-1))
  if (m[1]<m[3]) a1<--a1
  m4<-as.matrix(c(m[1],m[4]+a1,m[4]-a1,m[3]))
  a1<-sqrt(sigma50/(nn[5]-1))
  if (m[1]<m[3]) {a1<--a1}
  m5<-as.matrix(c(m[1],m[5]+2.4*a1,m4[2],m[5]+0.8*a1,0.75*m[1]+0.25*m[3],m[5]-0.8*a1,m4[3],m[5]-2.4*a1,m[3]))
  #######first order parameters############
  hh9<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,
                1,-1,1,-1,1,1,0.5,0.5,0.5,-1,-1,
                1,-1,-1,1,0.5,-1,1,0.5,-1,1,0.5),11,3)
  mm<-as.matrix(c(m[1],m[3],m4[2],m4[3],m5[2],m5[3],m5[4],m5[5],m5[6],m5[7],m5[8]))
  B9<-solve(crossprod(hh9,hh9))%*%crossprod(hh9,mm)
  a1<-B9[2]	   #da.
  a2<-B9[3]   #db.
  a1<-a1^2;a2<-a2^2
  sigma5[2]<-(sigma+0.75*a2)/m_nf;sigma5[4]<-(sigma+0.75*a1)/m_nf
  sigma5[5]<-(sigma+0.75*(a1+a2))/m_nf;sigma5[6]<-sigma5[4]
  sigma5[7]<-sigma/m_nf;sigma5[8]<-sigma5[2];sigma5[9]<-sigma/m_nf	 
  mi4<-mix_pi4[c(1:4)];mi5<-mix_pi5[c(1:9)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
    sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
  ########################################################
  ##########iteration process###########
  ############E-step###################
  ############F2 process###############
  d4<-4;d3<-9
  iteration <- 0; stopa <- 1000
  W4 <- matrix(0,4,m_sam[4]); swx24 <- matrix(0,4,1)
  W5 <- matrix(0,9,m_sam[5]); swx25 <- matrix(0,9,1)
  hh<-matrix(0,6,6);b_line<-matrix(0,6,1)
  s0<-matrix(0,20,1);n0<-matrix(0,20,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:d4) { W4[i,] <- mix_pi4[i]*dnorm(dataF2,m4[i],sqrt(sigma4[i]))/dmixnorm(dataF2,m4,sqrt(sigma4),mix_pi4) }
    mi4 <- as.matrix(rowSums(W4)/m_sam[4])
    sumwx4 <- W4%*%dataF2
    #############F3 process############
    for(i in 1:d3) { W5[i,] <- mix_pi5[i]*dnorm(dataF3,m5[i],sqrt(sigma5[i]))/dmixnorm(dataF3,m5,sqrt(sigma5),mix_pi5) }
    mi5 <- as.matrix(rowSums(W5)/m_sam[5])
    sumwx5 <- W5%*%dataF3
    aaa0<-0
    s0[1]<-sumx[1]+sumx[2]+sumwx4[1]+m_nf*sumwx5[1];s0[2]<-sumx[3]+sumwx4[4]+m_nf*sumwx5[9]
    s0[3]<-sumwx4[2]+m_nf*sumwx5[3];s0[4]<-sumwx4[3]+m_nf*sumwx5[7]
    s0[5]<-sumwx5[2];s0[6]<-sumwx5[4];s0[7]<-sumwx5[5]
    s0[8]<-sumwx5[6];s0[9]<-sumwx5[8]
    n0[1]<-nn[1]+nn[2]+mi4[1]*nn[4]+m_nf*mi5[1]*nn[5]
    n0[2]<-nn[3]+mi4[4]*nn[4]+m_nf*mi5[9]*nn[5]
    n0[3]<-mi4[2]*nn[4]+m_nf*mi5[3]*nn[5]
    n0[4]<-mi4[3]*nn[4]+m_nf*mi5[7]*nn[5]
    n0[5]<-mi5[2]*nn[5];n0[6]<-mi5[4]*nn[5]
    n0[7]<-mi5[5]*nn[5];n0[8]<-mi5[6]*nn[5]
    n0[9]<-mi5[8]*nn[5]
    n0[c(1:9)][abs(n0[c(1:9)])<0.000001]<-0.000001
    aaa1<-1000;n_iter<-0;AA<-matrix(0,6,1)
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      ##########first order parameters###################
      mm<-as.matrix(c(m[1],m[3],m4[2],m4[3],m5[2],m5[3],m5[4],m5[5],m5[6],m5[7],m5[8]))
      B91<-solve(crossprod(hh9,hh9))%*%crossprod(hh9,mm)
      g_aa1<-0.75*B91[3]^2               #  gs[0]: da.
      #   0.75*db**2.                  #  gs[1]: db.
      g_aa2<-0.75*B91[2]^2
      #   0.75*da**2.
      g_aa3<-g_aa1+g_aa2
      #   0.75*(da**2+db**2).
      sigma5[c(1,3,7,9)]<-sigma/m_nf;sigma5[2]<-(sigma+g_aa1)/m_nf
      sigma5[4]<-(sigma+g_aa2)/m_nf;sigma5[5]<-(sigma+g_aa3)/m_nf
      sigma5[6]<-sigma5[4];sigma5[8]<-sigma5[2]
      #################restrictions########################################
      hh[1,1]<-sigma*(1/n0[1]+1/n0[2]+1/n0[3]+1/n0[4])
      hh[1,2]<-sigma*(3/n0[1]+1/n0[2])
      hh[1,3]<-sigma*(3/n0[1]-1/n0[3])
      hh[1,4]<-sigma*(3/n0[1]-1/n0[4])
      hh[1,5]<-sigma*(1/n0[2]-3/n0[3])
      hh[1,6]<-sigma*(1/n0[2]-3/n0[4])
      hh[2,2]<-sigma*(9/n0[1]+1/n0[2])+16*sigma5[5]/n0[7]
      hh[2,3]<-sigma*9/n0[1]
      hh[2,4]<-sigma*9/n0[1]
      hh[2,5]<-sigma/n0[2]
      hh[2,6]<-sigma/n0[2]
      hh[3,3]<-sigma*(9/n0[1]+1/n0[3])+16*sigma5[2]/n0[5]
      hh[3,4]<-sigma*9/n0[1]
      hh[3,5]<-sigma*3/n0[3]
      hh[3,6]<-0
      hh[4,4]<-sigma*(9/n0[1]+1/n0[4])+16*sigma5[4]/n0[6]
      hh[4,5]<-0
      hh[4,6]<-3*sigma/n0[4]
      hh[5,5]<-sigma*(1/n0[2]+9/n0[3])+16*sigma5[6]/n0[8]
      hh[5,6]<-sigma/n0[2]
      hh[6,6]<-sigma*(1/n0[2]+9/n0[4])+16*sigma5[8]/n0[9]
      for(i in 2:6)
      {
        for(j in 1:(i-1))
        {
          hh[i,j]<-hh[j,i]
        }
      }
      ###########################################################
      b_line[1]<-s0[1]/n0[1]+s0[2]/n0[2]-s0[3]/n0[3]-s0[4]/n0[4]
      b_line[2]<-3*s0[1]/n0[1]+s0[2]/n0[2]-4*s0[7]/n0[7]
      b_line[3]<-3*s0[1]/n0[1]+s0[3]/n0[3]-4*s0[5]/n0[5]
      b_line[4]<-3*s0[1]/n0[1]+s0[4]/n0[4]-4*s0[6]/n0[6]
      b_line[5]<-s0[2]/n0[2]+3*s0[3]/n0[3]-4*s0[8]/n0[8]
      b_line[6]<-s0[2]/n0[2]+3*s0[4]/n0[4]-4*s0[9]/n0[9]
      B911<-solve(hh,b_line)
      ###########################################################
      m[1]<-(s0[1]-sigma*(B911[1]+3*B911[2]+3*B911[3]+3*B911[4]))/n0[1]
      m[3]<-(s0[2]-sigma*(B911[1]+B911[2]+B911[5]+B911[6]))/n0[2]
      m4[2]<-(s0[3]+sigma*(B911[1]-B911[3]-3*B911[5]))/n0[3]
      m4[3]<-(s0[4]+sigma*(B911[1]-B911[4]-3*B911[6]))/n0[4]
      m5[2]<-(s0[5]+sigma5[2]*4*B911[3])/n0[5]
      m5[4]<-(s0[6]+sigma5[4]*4*B911[4])/n0[6]
      m5[5]<-(s0[7]+sigma5[5]*4*B911[2])/n0[7]
      m5[6]<-(s0[8]+sigma5[6]*4*B911[5])/n0[8]
      m5[8]<-(s0[9]+sigma5[8]*4*B911[6])/n0[9]
      m[2]<-m[1];m4[1]<-m[1];m5[1]<-m[1]
      m4[4]<-m[3];m5[9]<-m[3];m5[3]<-m4[2];m5[7]<-m4[3]
      aaa1<-max(abs(B911-AA))
      AA<-B911
      if (n_iter>20) break
    }  
    ##########obtain variance###############
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d4) {swx24[i] <- W4[i,]%*%(dataF2-m4[i])^2 } ;for(i in 1:d3) {swx25[i] <- W5[i,]%*%(dataF3-m5[i])^2 }      
    s0[10]<-sum(swx24)+ss1+ss2+ss3+(swx25[1]+swx25[3]+swx25[7]+swx25[9])*m_nf
    n0[10]<-nn[1]+nn[2]+nn[3]+nn[4]+(mi5[1]+mi5[3]+mi5[7]+mi5[9])*nn[5]
    n0[11]<-(mi5[2]+mi5[8])*nn[5];n0[12]<-(mi5[4]+mi5[6])*nn[5]
    n0[13]<-mi5[5]*nn[5]
    mm<-as.matrix(c(m[1],m[3],m4[2],m4[3],m5[2],m5[3],m5[4],m5[5],m5[6],m5[7],m5[8]))
    B9111<-solve(crossprod(hh9,hh9))%*%crossprod(hh9,mm)
    g_aa1<-0.75*B9111[3]^2               #  gs[0]: da.
    #   0.75*db**2.                  #  gs[1]: db.
    g_aa2<-0.75*B9111[2]^2
    #   0.75*da**2.
    g_aa3<-g_aa1+g_aa2
    #   0.75*(da**2+db**2).
    aaa0<-sigma;n_iter<-0;ab4<-1000
    while (ab4>0.001){
      n_iter<-n_iter+1
      ab1<-sigma/(sigma+g_aa1)  # g_aa1:0.75db*db.
      ab2<-sigma/(sigma+g_aa2)  # g_aa2:0.75da*da.
      ab3<-sigma/(sigma+g_aa3)
      sigma<-(s0[10]+m_nf*(ab1^2*(swx25[2]+swx25[8])+ab2^2*(swx25[4]+swx25[6])+ab3^2*swx25[5]))/(n0[10]+ab1*n0[11]+ab2*n0[12]+ab3*n0[13])
      ab4<-abs(sigma-aaa0)
      aaa0<-sigma
      if (n_iter>20) break
    } 
    if (sigma<sigma0) {sigma<-sigma0}
    sigma5[c(1,3,7,9)]<-sigma/m_nf;sigma5[2]<-(sigma+g_aa1)/m_nf
    sigma5[4]<-(sigma+g_aa2)/m_nf;sigma5[5]<-(sigma+g_aa3)/m_nf
    sigma5[6]<-sigma5[4];sigma5[8]<-sigma5[2]
    ##############################criteria for iterations to stop######################## 
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
      sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*4
  #######first order parameters############
  mm<-as.matrix(c(m[1],m[3],m4[2],m4[3],m5[2],m5[3],m5[4],m5[5],m5[6],m5[7],m5[8]))
  B91111<-solve(crossprod(hh9,hh9))%*%crossprod(hh9,mm)
  ########second order genetic parameters###############   
  F2jj<-sigmaF2-sigma
  if(F2jj<0){F2jj<-0}
  F2ll<-F2jj/sigmaF2
  F3jj<-sigmaF3-sigma/m_nf
  if(F3jj<0) {F3jj<-0}
  F3ll<-F3jj/sigmaF3   
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma))
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
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma))
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
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma))
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
    F2gg <- (dataF2 - m4[i])/sqrt(sigma4[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mi4[i]
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
  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3); 
  F3w1<-1/(12*m_sam[5])
  F3bmw <- matrix(0,m_sam[5],1); F3bmwsl <- matrix(0,m_sam[5],d3)
  for(i in 1:d3){
    F3gg <- (dataF3 - m5[i])/sqrt(sigma5[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mi5[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[5]){F3P2 <- F3P2+runif(m_sam[5])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[5]) + sum((F3P2 - (as.matrix(c(1:m_sam[5])) - 0.5)/m_sam[5])^2)
  F3u <- as.matrix(c(12*m_sam[5]*((F3dd[1]/m_sam[5]-0.5)^2),((45*m_sam[5])/4)*((F3dd[2]/m_sam[5]-1/3)^2),180*m_sam[5]*((F3dd[3]/m_sam[5]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[5]))))
  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("2MG-CD",round(abc,4),round(AIC,4),round(m[1],4),round(m[2],4),round(m[3],4), round(t(m4),4)," "," "," "," "," ",round(sigma4[1],4),
                       round(t(mi4),4)," "," "," "," "," ",round(t(m5),4),round(t(sigma5),4),round(t(mi5),4),round(sigma,4),          
                       round(B91111[1],4)," "," "," "," ",round(B91111[2],4),round(B91111[3],4)," "," "," "," "," "," "," "," ",round(F2jj,4),round(F2ll*100,4)," "," ",round(F3jj,4),round(F3ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
 output<-as.matrix(output)
  OUTPUT<-list(output,mix_pi4,mix_pi5)
  return(OUTPUT)
}
##############2MG-EAD(B-6)###################################
G5ModelFun[[10]] <- function(K1,logL,df11,df21,df31,df41,df51,G5text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataF2 <- as.matrix(as.numeric(df41[,1]))
  dataF3 <- as.matrix(as.numeric(df51[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2),sum(dataF3)))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2),mean(dataF3)))
  sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1);sigmaF2<-var(dataF2);sigmaF3<-var(dataF3)
  sumx6<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  sigma<-sumx6;sigma0<-sigma;sigma40<-sigmaF2;sigma400<-sigma40;sigma50<-sigmaF3;sigma500<-sigma50
  m_nf <- as.numeric(G5text2)
  ###############procedure start###########################
  sigma4<-matrix(0,3,1);sigma5<-matrix(0,6,1)
  sigma4[c(1:3)]<-sigma;sigma5[1]<-sigma/m_nf
  mix_pi4<-as.matrix(c(0.5625,0.375,0.0625))
  mix_pi5<-as.matrix(c(0.0625,0.25,0.125,0.25,0.25,0.0625))
  nn<-m_sam
  m4<-as.matrix(c(m[1],(m[1]+m[3])/2,m[3]))
  m5<-matrix(0,6,1)
  m5[1]<-m[1];m5[6]<-m[3];m5[4]<-m4[2]
  m5[3]<-(m5[1]+m5[4])/2;m5[2]<-(m5[1]+m5[3])/2
  m5[5]<-0.625*m[1]+0.375*m[3]
  #######first order parameters############
  hh10<-matrix(c(1,1,1,1,1,1,1,2,-2,0,1.5,1,0,-0.5),7,2)
  mm<-as.matrix(c(m[1],m[3],m4[2],m5[2],m5[3],m5[4],m5[5]))
  B10<-solve(crossprod(hh10,hh10))%*%crossprod(hh10,mm)
  a1<-B10[2]	   
  a1<-a1^2
  sigma5[2]<-(sigma+0.75*a1)/m_nf;sigma5[3]<-(sigma+1.5*a1)/m_nf
  sigma5[4]<-sigma/m_nf;sigma5[5]<-(sigma+0.75*a1)/m_nf;sigma5[6]<-sigma/m_nf
  mi4<-mix_pi4[c(1:3)];mi5<-mix_pi5[c(1:6)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
    sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
  ########################################################
  ##########iteration process###########
  ############E-step###################
  ############F2 process###############
  d1<-3
  d6<-6
  iteration <- 0; stopa <- 1000
  W4 <- matrix(0,3,m_sam[4]); swx24 <- matrix(0,3,1)
  W5 <- matrix(0,6,m_sam[5]); swx25 <- matrix(0,6,1)
  hh<-matrix(0,4,4);b_line<-matrix(0,4,1)
  s0<-matrix(0,20,1);n0<-matrix(0,20,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:d1) { W4[i,] <- mix_pi4[i]*dnorm(dataF2,m4[i],sqrt(sigma4[i]))/dmixnorm(dataF2,m4,sqrt(sigma4),mix_pi4) }
    mi4 <- as.matrix(rowSums(W4)/m_sam[4])
    sumwx4 <- W4%*%dataF2
    #############F3 process############
    for(i in 1:d6) { W5[i,] <- mix_pi5[i]*dnorm(dataF3,m5[i],sqrt(sigma5[i]))/dmixnorm(dataF3,m5,sqrt(sigma5),mix_pi5) }
    mi5 <- as.matrix(rowSums(W5)/m_sam[5])
    sumwx5 <- W5%*%dataF3
    aaa0<-0
    s0[1]<-sumx[1]+sumx[2]+sumwx4[1]+m_nf*sumwx5[1]
    s0[2]<-sumx[3]+sumwx4[3]+m_nf*sumwx5[6]
    s0[3]<-sumwx4[2]+m_nf*sumwx5[4]
    s0[4]<-sumwx5[2];s0[5]<-sumwx5[3];s0[6]<-sumwx5[5]
    n0[1]<-nn[1]+nn[2]+mi4[1]*nn[4]+m_nf*mi5[1]*nn[5]
    n0[2]<-nn[3]+mi4[3]*nn[4]+m_nf*mi5[6]*nn[5]
    n0[3]<-mi4[2]*nn[4]+m_nf*mi5[4]*nn[5]
    n0[4]<-mi5[2]*nn[5];n0[5]<-mi5[3]*nn[5];n0[6]<-mi5[5]*nn[5]
    n0[c(1:6)][abs(n0[c(1:6)])<0.000001] <- 0.000001
    aa1<-1000;n_iter<-0;AA<-matrix(0,4,1)
    while (aa1>0.0001)
    {
      n_iter<-n_iter+1
      ##########first order parameters###################
      mm<-as.matrix(c(m[1],m[3],m4[2],m5[2],m5[3],m5[4],m5[5]))
      B101<-solve(crossprod(hh10,hh10))%*%crossprod(hh10,mm)
      g_aa1<-0.75*B101[2]^2
      #   0.75*d*d.
      g_aa2<-1.5*B101[2]^2
      #   1.5*d*d.
      sigma5[1]<-sigma/m_nf;sigma5[2]<-(sigma+g_aa1)/m_nf;sigma5[3]<-(sigma+g_aa2)/m_nf
      sigma5[4]<-sigma/m_nf;sigma5[5]<-sigma5[2];sigma5[6]<-sigma/m_nf
      #################restrictions########################################
      hh[1,1]<-sigma*(1/n0[1]+1/n0[2]+4/n0[3])
      hh[1,2]<-sigma*(3/n0[1]+1/n0[2])
      hh[1,3]<-sigma*(3/n0[1]-2/n0[3])
      hh[1,4]<-sigma*(1/n0[2]-6/n0[3])
      hh[2,2]<-sigma*(9/n0[1]+1/n0[2])+16*sigma5[3]/n0[5]
      hh[2,3]<-9*sigma/n0[1]
      hh[2,4]<-sigma/n0[2]
      hh[3,3]<-sigma*(9/n0[1]+1/n0[3])+16*sigma5[2]/n0[4]
      hh[3,4]<-sigma*3/n0[3]
      hh[4,4]<-sigma*(9/n0[3]+1/n0[2])+16*sigma5[5]/n0[6]
      for(i in 2:4)
      {
        for(j in 1:(i-1))
        {
          hh[i,j]<-hh[j,i]
        }
      }
      #################################################################### 
      b_line[1]<-s0[1]/n0[1]+s0[2]/n0[2]-2*s0[3]/n0[3]
      b_line[2]<-3*s0[1]/n0[1]+s0[2]/n0[2]-4*s0[5]/n0[5]
      b_line[3]<-3*s0[1]/n0[1]+s0[3]/n0[3]-4*s0[4]/n0[4]
      b_line[4]<-s0[2]/n0[2]+3*s0[3]/n0[3]-4*s0[6]/n0[6]
      B1011<-solve(hh,b_line)
      ###################################################
      m[1]<-(s0[1]-sigma*(B1011[1]+3*B1011[2]+3*B1011[3]))/n0[1]
      m[3]<-(s0[2]-sigma*(B1011[1]+B1011[2]+B1011[4]))/n0[2]
      m4[2]<-(s0[3]+sigma*(2*B1011[1]-B1011[3]-3*B1011[4]))/n0[3]
      m5[2]<-(s0[4]+sigma5[2]*4*B1011[3])/n0[4]
      m5[3]<-(s0[5]+sigma5[3]*4*B1011[2])/n0[5]
      m5[5]<-(s0[6]+sigma5[5]*4*B1011[4])/n0[6]
      m[2]<-m[1];m4[1]<-m[1];m5[1]<-m[1]
      m4[3]<-m[3];m5[6]<-m[3];m5[4]<-m4[2]
      aa1<-max(abs(B1011-AA))
      AA<-B1011
      if (n_iter>20) break
    }  
    ##########obtain variance###############
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d1) {swx24[i] <- W4[i,]%*%(dataF2-m4[i])^2 } ; for(i in 1:d6) {swx25[i] <- W5[i,]%*%(dataF3-m5[i])^2 }      
    s0[7]<-ss1+ss2+ss3+swx24[1]+swx24[2]+swx24[3]+m_nf*(swx25[1]+swx25[4]+swx25[6])
    n0[7]<-nn[1]+nn[2]+nn[3]+nn[4]+(mi5[1]+mi5[4]+mi5[6])*nn[5]
    n0[8]<-(mi5[2]+mi5[5])*nn[5];n0[9]<-mi5[3]*nn[5]
    mm<-as.matrix(c(m[1],m[3],m4[2],m5[2],m5[3],m5[4],m5[5]))
    B10111<-solve(crossprod(hh10,hh10))%*%crossprod(hh10,mm)
    g_aa1<-0.75*B10111[2]^2
    #   0.75*d*d.
    g_aa2<-1.5*B10111[2]^2
    #   1.5*d*d.
    aaa0<-sigma;ab5<-0;n_iter<-0;aa4<-1000
    while (aa4>0.0001){
      ab5<-ab5+1
      aa1<-sigma/(sigma+g_aa1)
      aa2<-sigma/(sigma+g_aa2)
      sigma<-(s0[7]+m_nf*(aa1^2*(swx25[2]+swx25[5])+aa2^2*swx25[3]))/(n0[7]+aa1*n0[8]+aa2*n0[9])
      aa4<-abs(sigma-aaa0)
      aaa0<-sigma
      if (ab5>20) break
    } 
    if (sigma<sigma0) {sigma<-sigma0}
    sigma5[1]<-sigma/m_nf;sigma5[2]<-(sigma+g_aa1)/m_nf
    sigma5[3]<-(sigma+g_aa2)/m_nf;sigma5[4]<-sigma/m_nf
    sigma5[5]<-sigma5[2];sigma5[6]<-sigma/m_nf
    ##############################criteria for iterations to stop######################## 
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
      sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*3
  #######first order parameters############
  mm<-as.matrix(c(m[1],m[3],m4[2],m5[2],m5[3],m5[4],m5[5]))
  B101111<-solve(crossprod(hh10,hh10))%*%crossprod(hh10,mm)
  ########second order genetic parameters###############   
  F2jj<-sigmaF2-sigma
  if(F2jj<0){F2jj<-0}
  F2ll<-F2jj/sigmaF2
  F3jj<-sigmaF3-sigma/m_nf
  if(F3jj<0) {F3jj<-0}
  F3ll<-F3jj/sigmaF3   
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma))
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
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma))
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
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma))
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
    F2gg <- (dataF2 - m4[i])/sqrt(sigma4[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mi4[i]
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
  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3); 
  F3w1<-1/(12*m_sam[5])
  F3bmw <- matrix(0,m_sam[5],1); F3bmwsl <- matrix(0,m_sam[5],d6)
  for(i in 1:d6){
    F3gg <- (dataF3 - m5[i])/sqrt(sigma5[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mi5[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[5]){F3P2 <- F3P2+runif(m_sam[5])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[5]) + sum((F3P2 - (as.matrix(c(1:m_sam[5])) - 0.5)/m_sam[5])^2)
  F3u <- as.matrix(c(12*m_sam[5]*((F3dd[1]/m_sam[5]-0.5)^2),((45*m_sam[5])/4)*((F3dd[2]/m_sam[5]-1/3)^2),180*m_sam[5]*((F3dd[3]/m_sam[5]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[5]))))
  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("2MG-EAD",round(abc,4),round(AIC,4),round(m[1],4),round(m[2],4),round(m[3],4), round(t(m4),4)," "," "," "," "," "," ",round(sigma4[1],4),
                       round(t(mi4),4)," "," "," "," "," "," ",round(t(m5),4)," "," "," ",round(t(sigma5),4)," "," "," ",round(t(mi5),4)," "," "," ",round(sigma,4),          
                       round(B101111[1],4)," "," "," "," ",round(B101111[2],4)," "," "," "," "," "," "," "," "," ",round(F2jj,4),round(F2ll*100,4)," "," ",round(F3jj,4),round(F3ll*100,4)," "," ",
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
 output<-as.matrix(output)
  OUTPUT<-list(output,mix_pi4,mix_pi5)
  return(OUTPUT)
}
###############PG-ADI(C-0)########################
G5ModelFun[[11]] <- function(K1,logL,df11,df21,df31,df41,df51,G5text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataF2 <- as.matrix(as.numeric(df41[,1]))
  dataF3 <- as.matrix(as.numeric(df51[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2),sum(dataF3)))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2),mean(dataF3)))
  sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1);sigmaF2<-var(dataF2);sigmaF3<-var(dataF3)
  sumx6<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  sigma<-sumx6;sigma0<-sigma;sigma40<-sigmaF2;sigma400<-sigma40;sigma50<-sigmaF3;sigma500<-sigma50
  m_nf <- as.numeric(G5text2)
  ###############procedure start###########################
  nn<-m_sam
  mi4<-1;mi5<-1
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dnorm(dataF2,m[4],sqrt(sigma40))))+
    sum(log(dnorm(dataF3,m[5],sqrt(sigma50)))) 
  ########################################################
  ##########iteration process###########
  ############E-step###################
  ############F2 process###############
  iteration <- 0; stopa <- 1000
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    m[c(1:5)]<-sumx/m_sam
    sigma4<-matrix(0,1,1);sigma5<-matrix(0,1,1)
    D2<-sum((dataF2-m[4])^2);sigma4[1]<-D2/nn[4]
    D3<-sum((dataF3-m[5])^2);sigma5[1]<-D3/nn[5]
    ##########obtain variance###############
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2)
    ss2<-sum((dataF1-m[2])^2);ss4<-sum((dataF2-m[4])^2)
    ss5<-sum((dataF3-m[5])^2)
    abc1<-ss1+ss2+ss3;abc2<-nn[1]+nn[2]+nn[3]
    aaa0<-sigma;aaa1<-1000;n_iter<-0
    while (aaa1>0.0001){
      n_iter<-n_iter+1
      aa1<-sigma/sigma4[1]
      aa2<-sigma/sigma5[1]
      sigma<-(abc1+aa1*aa1*ss4+aa2*aa2*ss5)/(abc2+aa1*nn[4]+aa2*nn[5])
      aaa1<-abs(sigma-aaa0)
      aaa0<-sigma
      if(n_iter>20)break
    }
    if (sigma<sigma0) {sigma<-sigma0}
    sigma40<-sigma4[1]-sigma;sigma50<-sigma5[1]-sigma/m_nf
    ##############################criteria for iterations to stop######################## 
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dnorm(dataF2,m[4],sqrt(sigma4))))+
      sum(log(dnorm(dataF3,m[5],sqrt(sigma5))))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*10
  #######first order parameters############
  ma1<-m[1];ma2<-m[2];ma3<-m[3];ma4<-m[4];ma5<-m[5]
  ########second order genetic parameters###############   
  F2_gg <- sigmaF2-sigma
  if(F2_gg<0) {F2_gg<-0}
  F2_rr <- F2_gg/sigmaF2
  F3_gg <- sigmaF3-sigma/m_nf
  if(F3_gg<0) {F3_gg<-0}
  F3_rr <- F3_gg/sigmaF3
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma))
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
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma))
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
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma))
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
  F2gg <- (dataF2 - m[4])/sqrt(as.vector(sigma))
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
  ####################################hypothesis testing for F3#########################################
  dataF3<-sort(dataF3)
  F3w1<-1/(12*m_sam[5])
  F3bmw <- matrix(0,m_sam[5],1)
  F3gg <- (dataF3 - m[5])/sqrt(as.vector(sigma))
  F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
  F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3bmw)))[1]
  if(nn < m_sam[5]){F3bmw <- F3bmw+runif(m_sam[5])/1e4}
  ##########################################################
  F3dd<-c((sum(F3bmw)),(sum(F3bmw^2)),sum((F3bmw-0.5)^2))
  F3w<-F3w1+sum((F3bmw - (as.matrix(c(1:m_sam[5])) - 0.5)/m_sam[5])^2)
  F3u<- as.matrix(c(12*m_sam[5]*((F3dd[1]/m_sam[5]-0.5)^2),((45*m_sam[5])/4)*((F3dd[2]/m_sam[5]-1/3)^2),180*m_sam[5]*((F3dd[3]/m_sam[5]-1/12)^2)))
  F3D<-as.numeric(ks.test(F3bmw,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3w),(1-pkolm(F3D,m_sam[5]))))
  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("PG-ADI",round(abc,4),round(AIC,4),round(m[1],4),round(m[2],4),round(m[3],4), round(m[4],4)," "," "," "," "," "," "," "," ",round(sigma4[1],4),
                       round(mi4,4)," "," "," "," "," "," "," "," ",round(m[5],4)," "," "," "," "," "," "," "," ",round(t(sigma5),4)," "," "," "," "," "," "," "," ",round(mi5,4)," "," "," "," "," "," "," "," ",round(sigma,4),          
                       round(ma1,4),round(ma2,4),round(ma3,4),round(ma4,4),round(ma5,4)," "," "," "," "," "," "," "," "," "," "," "," ",round(F2_gg,4),round(F2_rr*100,4)," "," ",round(F3_gg,4),round(F3_rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2w,4),F2tt[4],round(F2D,4),F2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3w,4),F3tt[4],round(F3D,4),F3tt[5])
 output<-as.matrix(output)
  OUTPUT<-list(output)
  return(OUTPUT)
}
############PG-AD(C-1)###########################
G5ModelFun[[12]] <- function(K1,logL,df11,df21,df31,df41,df51,G5text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataF2 <- as.matrix(as.numeric(df41[,1]))
  dataF3 <- as.matrix(as.numeric(df51[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2),sum(dataF3)))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2),mean(dataF3)))
  sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1);sigmaF2<-var(dataF2);sigmaF3<-var(dataF3)
  sumx6<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  sigma<-sumx6;sigma0<-sigma;sigma40<-sigmaF2;sigma400<-sigma40;sigma50<-sigmaF3;sigma500<-sigma50
  m_nf <- as.numeric(G5text2)
  ###############procedure start###########################
  nn<-m_sam
  mi4<-1;mi5<-1
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dnorm(dataF2,m[4],sqrt(sigma40))))+
    sum(log(dnorm(dataF3,m[5],sqrt(sigma50)))) 
  ########################################################
  ##########iteration process###########
  ############E-step###################
  ############F2 process###############
  iteration <- 0; stopa <- 1000; rr<-matrix(0,2,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    n_iter<-0;aaa1<-1000;AA<-matrix(0,2,1)
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      sigma40<-sum((dataF2-m[4])^2);sigma40<-sigma40/nn[4]
      sigma50<-sum((dataF3-m[5])^2);sigma50<-sigma50/nn[5]
      aa1<-sigma*(1/nn[1]+4/nn[2]+1/nn[3])+16*sigma40/nn[4]
      aa2<-2*sigma/nn[2]+12*sigma40/nn[4]
      aa3<-sigma/nn[2]+9*sigma40/nn[4]+4*sigma50/nn[5]
      aa4<-sumx[1]/nn[1]+2*sumx[2]/nn[2]+sumx[3]/nn[3]-4*sumx[4]/nn[4]
      aa5<-sumx[2]/nn[2]-3*sumx[4]/nn[4]+2*sumx[5]/nn[5]
      aa6<-aa1*aa3-aa2^2
      rr[1]<-(aa3*aa4-aa2*aa5)/aa6;rr[2]<-(aa1*aa5-aa2*aa4)/aa6
      m[1]<-(sumx[1]-sigma*rr[1])/nn[1]
      m[2]<-(sumx[2]-sigma*(2*rr[1]+rr[2]))/nn[2]
      m[3]<-(sumx[3]-sigma*rr[1])/nn[3]
      m[4]<-(sumx[4]+sigma40*(4*rr[1]+3*rr[2]))/nn[4]
      m[5]<-(sumx[5]-2*rr[2]*sigma50)/nn[5]
      aaa1<-max(abs(rr-AA))
      AA<-rr
      if (n_iter>20) break
    }
    aaa0<-sigma
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    abc1<-ss1+ss2+ss3;abc2<-nn[1]+nn[2]+nn[3]
    ss4<-sum((dataF2-m[4])^2);ss5<-sum((dataF3-m[5])^2)
    sigma40<-ss4/nn[4];sigma_4<-sigma40-sigma
    sigma50<-ss5/nn[5];sigma_5<-sigma50-sigma/m_nf
    aa3<-1000;n_iter<-0
    while (aa3>0.0001){
      n_iter<-n_iter+1
      aa1<-sigma/sigma40
      if (aa1>=1) {aa1<-1}
      aa2<-sigma/sigma50
      if (aa2>=1) {aa2<-1}
      aa4<-abc1+aa1^2*ss4+aa2^2*ss5
      aa5<-abc2+aa1*nn[4]+aa2*nn[5]
      sigma<-aa4/aa5
      aa3<-abs(sigma-aaa0)
      aaa0<-sigma
      if(n_iter>20)break
    }
    if (sigma<sigma0) {sigma<-sigma0}
    sigma40<-sigma_4+sigma;sigma50<-sigma_5+sigma/m_nf
    ##############################criteria for iterations to stop######################## 
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dnorm(dataF2,m[4],sqrt(sigma40))))+
      sum(log(dnorm(dataF3,m[5],sqrt(sigma50))))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*7
  #######first order parameters############
  hh11<-matrix(c(1,1,1,1,1,1,0,-1,0,0,0,1,0,0.5,0.25),5,3)
  B11<-solve(crossprod(hh11,hh11))%*%crossprod(hh11,m)
  ########second order genetic parameters###############   
  F2_gg <- sigmaF2-sigma
  if(F2_gg<0) {F2_gg<-0}
  F2_rr <- F2_gg/sigmaF2
  F3_gg <- sigmaF3-sigma/m_nf
  if(F3_gg<0) {F3_gg<-0}
  F3_rr <- F3_gg/sigmaF3
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma))
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
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma))
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
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma))
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
  F2gg <- (dataF2 - m[4])/sqrt(as.vector(sigma))
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
  ####################################hypothesis testing for F3#########################################
  dataF3<-sort(dataF3)
  F3w1<-1/(12*m_sam[5])
  F3bmw <- matrix(0,m_sam[5],1)
  F3gg <- (dataF3 - m[5])/sqrt(as.vector(sigma))
  F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
  F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3bmw)))[1]
  if(nn < m_sam[5]){F3bmw <- F3bmw+runif(m_sam[5])/1e4}
  ##########################################################
  F3dd<-c((sum(F3bmw)),(sum(F3bmw^2)),sum((F3bmw-0.5)^2))
  F3w<-F3w1+sum((F3bmw - (as.matrix(c(1:m_sam[5])) - 0.5)/m_sam[5])^2)
  F3u<- as.matrix(c(12*m_sam[5]*((F3dd[1]/m_sam[5]-0.5)^2),((45*m_sam[5])/4)*((F3dd[2]/m_sam[5]-1/3)^2),180*m_sam[5]*((F3dd[3]/m_sam[5]-1/12)^2)))
  F3D<-as.numeric(ks.test(F3bmw,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3w),(1-pkolm(F3D,m_sam[5]))))
  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("PG-AD",round(abc,4),round(AIC,4),round(m[1],4),round(m[2],4),round(m[3],4), round(m[4],4)," "," "," "," "," "," "," "," ",round(sigma40,4),
                       round(mi4,4)," "," "," "," "," "," "," "," ",round(m[5],4)," "," "," "," "," "," "," "," ",round(sigma50,4)," "," "," "," "," "," "," "," ",round(mi5,4)," "," "," "," "," "," "," "," ",round(sigma,4),          
                       round(B11[1],4)," "," "," "," "," "," "," "," "," "," "," "," ",round(B11[2],4),round(B11[3],4)," "," ",round(F2_gg,4),round(F2_rr*100,4)," "," ",round(F3_gg,4),round(F3_rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2w,4),F2tt[4],round(F2D,4),F2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3w,4),F3tt[4],round(F3D,4),F3tt[5])
 output<-as.matrix(output)
  OUTPUT<-list(output)
  return(OUTPUT)
}
##################MX1-AD-ADI(D-0)############################
G5ModelFun[[13]] <- function(K1,logL,df11,df21,df31,df41,df51,G5text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataF2 <- as.matrix(as.numeric(df41[,1]))
  dataF3 <- as.matrix(as.numeric(df51[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2),sum(dataF3)))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2),mean(dataF3)))
  sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1);sigmaF2<-var(dataF2);sigmaF3<-var(dataF3)
  sumx6<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  sigma<-sumx6;sigma0<-sigma;sigma40<-sigmaF2;sigma400<-sigma40;sigma50<-sigmaF3;sigma500<-sigma50
  m_nf <- as.numeric(G5text2)
  ###############procedure start###########################
  d1<-3
  mix_pi4<-as.matrix(c(0.25,0.5,0.25));mix_pi5<-as.matrix(c(0.25,0.5,0.25))
  nn<-m_sam
  a1<-sqrt(sigma400/nn[4])
  if (m[1]<m[3]) {a1<--a1}
  m4<-as.matrix(c(m[4]+2*a1,m[4],m[4]-2*a1))
  a1<-sqrt(sigma500/nn[5])
  if (m[1]<m[3]) {a1<--a1}
  m5<-as.matrix(c(m[5]+2*a1,m[5],m[5]-2*a1))
  ########first order genetic parameters###############
  hh12<-matrix(c(1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,
                 0,0,0,0,0,0,1,1,1,1,0,-1,1,0,-1,1,0,-1,0,1,0,0,1,0,0,0.5,0),9,7)
  mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m5[1],m5[2],m5[3]))
  B12<-solve(crossprod(hh12,hh12))%*%crossprod(hh12,mm)
  a1<-B12[6];a2<-B12[7]
  a1<-(0.5*a1^2+0.25*a2^2)/m_nf
  sigma4<-matrix(0,3,1);sigma5<-matrix(0,3,1)
  sigma4[1]<-sigma400/2;sigma4[2]<-sigma4[1];sigma4[3]<-sigma4[1]
  sigma5[1]<-sigma500/2;sigma5[3]<-sigma5[1];sigma5[2]<-sigma5[1]+a1
  mi4<-mix_pi4[c(1:3)];mi5<-mix_pi5[c(1:3)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
    sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
  ##########iteration process###########
  ############E-step###################
  ############F2 process###############
  iteration <- 0; stopa <- 1000
  W4 <- matrix(0,3,m_sam[4]); swx24 <- matrix(0,3,1)
  W5 <- matrix(0,3,m_sam[5]); swx25 <- matrix(0,3,1)
  n0<-matrix(0,6,1);s0<-matrix(0,6,1)
  rr<-matrix(0,2,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:d1) { W4[i,] <- mix_pi4[i]*dnorm(dataF2,m4[i],sqrt(sigma4[i]))/dmixnorm(dataF2,m4,sqrt(sigma4),mix_pi4) }
    mi4 <- as.matrix(rowSums(W4)/m_sam[4])
    n0[c(1:3)]<-as.matrix(rowSums(W4))
    sumwx4 <- W4%*%dataF2
    #############F3 process############
    for(i in 1:d1) { W5[i,] <- mix_pi5[i]*dnorm(dataF3,m5[i],sqrt(sigma5[i]))/dmixnorm(dataF3,m5,sqrt(sigma5),mix_pi5) }
    mi5 <- as.matrix(rowSums(W5)/m_sam[5])
    n0[c(4:6)]<-as.matrix(rowSums(W5));n0[c(1:6)][abs(n0[c(1:6)])<0.000001]<-0.000001
    sumwx5 <- W5%*%dataF3
    aaa0<-0
    m[1]<-sumx[1]/nn[1];m[2]<-sumx[2]/nn[2];m[3]<-sumx[3]/nn[3]
    s0[1]<-sumwx4[1]/n0[1]-sumwx4[3]/n0[3]-sumwx5[1]/n0[4]+sumwx5[3]/n0[6]
    s0[2]<-2*sumwx4[2]/n0[2]-2*sumwx4[3]/n0[3]+sumwx5[1]/n0[4]-4*sumwx5[2]/n0[5]+3*sumwx5[3]/n0[6]
    n_iter<-0;aaa1<-1000;AA<-matrix(0,2,1)
    while(aaa1>0.0001){
      n_iter<-n_iter+1 
      ########first order genetic parameters###############
      mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m5[1],m5[2],m5[3]))
      B121<-solve(crossprod(hh12,hh12))%*%crossprod(hh12,mm) 
      aa1<-B121[6];aa2<-B121[7];aa1<-(0.5*aa1^2+0.25*aa2^2)/m_nf
      sigma5[2]<-sigma5[1]+aa1
      abc1<-sigma4[1]/n0[1]+sigma4[3]/n0[3]+sigma5[1]/n0[4]+sigma5[3]/n0[6]
      abc2<-2*sigma4[3]/n0[3]-sigma5[1]/n0[4]+3*sigma5[3]/n0[6]
      abc3<-4*sigma4[2]/n0[2]+4*sigma4[3]/n0[3]+sigma5[1]/n0[4]+16*sigma5[2]/n0[5]+9*sigma5[3]/n0[6]
      aa2<-abc1*abc3-abc2*abc2
      aa3<-s0[1]*abc3-s0[2]*abc2
      aa4<-s0[2]*abc1-s0[1]*abc2
      rr[1]<-aa3/aa2;rr[2]<-aa4/aa2
      m4[1]<-(sumwx4[1]-rr[1]*sigma4[1])/n0[1]
      m4[2]<-(sumwx4[2]-2*rr[2]*sigma4[2])/n0[2]
      m4[3]<-(sumwx4[3]+sigma4[3]*(rr[1]+2*rr[2]))/n0[3]
      m5[1]<-(sumwx5[1]+(rr[1]-rr[2])*sigma5[1])/n0[4]
      m5[2]<-(sumwx5[2]+4*rr[2]*sigma5[2])/n0[5]
      m5[3]<-(sumwx5[3]-(rr[1]+3*rr[2])*sigma5[3])/n0[6]
      aaa1<-max(abs(rr-AA))
      AA<-rr
      if (n_iter>20) break
    }
    ########obtain variance#######################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d1) {swx24[i] <- W4[i,]%*%(dataF2-m4[i])^2 } ;for(i in 1:d1) {swx25[i] <- W5[i,]%*%(dataF3-m5[i])^2 }      
    sigma4[1]<-(swx24[1]+swx24[2]+swx24[3])/nn[4]
    sigma40<-sigma4[1]-sigma;
    if (sigma40<0) {sigma40<-0;sigma4[1]<-sigma}
    sigma4[1]<-sigma40+sigma;sigma4[2]<-sigma4[1];sigma4[3]<-sigma4[1]
    aaa0<-sigma5[1]
    mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m5[1],m5[2],m5[3]))
    B1211<-solve(crossprod(hh12,hh12))%*%crossprod(hh12,mm) 
    aa1<-B1211[6];aa2<-B1211[7];aa1<-(0.5*aa1^2+0.25*aa2^2)/m_nf
    n_iter<-0;aa3<-1000
    while (aa3>0.0001){
      n_iter<-n_iter+1
      ab3<-sigma5[1]/(sigma5[1]+aa1)
      sigma5[1]<-(swx25[1]+swx25[3]+ab3^2*swx25[2])/(n0[4]+n0[6]+ab3*n0[5])
      aa3<-abs(sigma5[1]-aaa0)
      aaa0<-sigma5[1]
      if (n_iter>20) break
    }
    sigma50<-sigma5[1]-sigma/m_nf
    if (sigma50<0) {sigma50<-0;sigma5[1]<-sigma/m_nf}
    sigma5[1]<-sigma50+sigma;sigma5[2]<-sigma5[1]+aa1;sigma5[3]<-sigma5[1]
    s0[1]<-ss1+ss2+ss3;s0[2]<-nn[1]+nn[2]+nn[3]
    aaa0<-0;n_iter<-0;aa3<-1000
    while (aa3>0.0001){
      n_iter<-n_iter+1
      abc1<-sigma/(sigma+sigma40)
      abc2<-(sigma/m_nf)/(sigma/m_nf+sigma50)
      abc3<-(sigma/m_nf)/(sigma/m_nf+sigma50+aa1)
      aa4<-s0[1]+abc1^2*(swx24[1]+swx24[2]+swx24[3])+m_nf*abc2^2*(swx25[1]+swx25[3])+m_nf*abc3^2*swx25[2]
      aa5<-s0[2]+abc1*nn[4]+abc2*(n0[4]+n0[6])+abc3*n0[5]
      sigma<-aa4/aa5
      aa3<-abs(sigma-aaa0)
      aaa0<-sigma
      if (n_iter>20) break 
    }
    aa4<-0.5*sigma0
    if (sigma<aa4) {sigma<-aa4}
    sigma4[1]<-sigma+sigma40;sigma4[2]<-sigma4[1]
    sigma4[3]<-sigma4[1];sigma5[1]<-sigma/m_nf+sigma50
    sigma5[3]<-sigma5[1];sigma5[2]<-sigma/m_nf+sigma50+aa1
    ##############################criteria for iterations to stop######################## 
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
      sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*10
  ########first order genetic parameters###############  
  mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m5[1],m5[2],m5[3]))
  B12111<-solve(crossprod(hh12,hh12))%*%crossprod(hh12,mm)     
  ########second order genetic parameters###############   
  F2_jj <- sigmaF2 - sigma4[1]
  F2_gg <- sigma4[1]-sigma
  if(F2_jj<0) {F2_jj<-0}
  if(F2_gg<0 || F2_gg>sigmaF2) {F2_gg<-0}
  F2_ll <- F2_jj/sigmaF2
  F2_rr <- F2_gg/sigmaF2
  F3_jj <- sigmaF3 - sigma5[1]
  F3_gg <- sigma5[1]-sigma/m_nf
  if(F3_jj<0) {F3_jj<-0}
  if(F3_gg<0 || F3_gg>sigmaF3) {F3_gg<-0}
  F3_ll <- F3_jj/sigmaF3
  F3_rr <- F3_gg/sigmaF3
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma))
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
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma))
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
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma))
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
    F2gg <- (dataF2 - m4[i])/sqrt(sigma4[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mi4[i]
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
  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3); 
  F3w1<-1/(12*m_sam[5])
  F3bmw <- matrix(0,m_sam[5],1); F3bmwsl <- matrix(0,m_sam[5],d1)
  for(i in 1:d1){
    F3gg <- (dataF3 - m5[i])/sqrt(sigma5[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mi5[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[5]){F3P2 <- F3P2+runif(m_sam[5])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[5]) + sum((F3P2 - (as.matrix(c(1:m_sam[5])) - 0.5)/m_sam[5])^2)
  F3u <- as.matrix(c(12*m_sam[5]*((F3dd[1]/m_sam[5]-0.5)^2),((45*m_sam[5])/4)*((F3dd[2]/m_sam[5]-1/3)^2),180*m_sam[5]*((F3dd[3]/m_sam[5]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[5]))))
  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("MX1-AD-ADI",round(abc,4),round(AIC,4),round(m[1],4),round(m[2],4),round(m[3],4), round(t(m4),4)," "," "," "," "," "," ",round(sigma4[1],4),
                       round(t(mi4),4)," "," "," "," "," "," ",round(t(m5),4)," "," "," "," "," "," ",round(t(sigma5),4)," "," "," "," "," "," ",round(t(mi5),4)," "," "," "," "," "," ",round(sigma,4),          
                       round(B12111[1],4),round(B12111[2],4),round(B12111[3],4),round(B12111[4],4),round(B12111[5],4),round(B12111[6],4)," ",round(B12111[7],4)," "," "," "," "," "," "," ",round(F2_jj,4),round(F2_ll*100,4),round(F2_gg,4),round(F2_rr*100,4),round(F3_jj,4),round(F3_ll*100,4),round(F3_gg,4),round(F3_rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
 output<-as.matrix(output)
  OUTPUT<-list(output,mix_pi4,mix_pi5)
  return(OUTPUT)
}
###############MX1-AD-AD(D-1)###############################################
G5ModelFun[[14]] <- function(K1,logL,df11,df21,df31,df41,df51,G5text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataF2 <- as.matrix(as.numeric(df41[,1]))
  dataF3 <- as.matrix(as.numeric(df51[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2),sum(dataF3)))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2),mean(dataF3)))
  sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1);sigmaF2<-var(dataF2);sigmaF3<-var(dataF3)
  sumx6<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  sigma<-sumx6;sigma0<-sigma;sigma40<-sigmaF2;sigma400<-sigma40;sigma50<-sigmaF3;sigma500<-sigma50
  m_nf <- as.numeric(G5text2)
  ###############procedure start###########################
  d1<-3
  mix_pi4<-as.matrix(c(0.25,0.5,0.25));mix_pi5<-as.matrix(c(0.25,0.5,0.25))
  nn<-m_sam
  a1<-sqrt(sigma400/nn[4])
  if (m[1]<m[3]) {a1<--a1}
  m4<-as.matrix(c(m[4]+2*a1,m[4],m[4]-2*a1))
  a1<-sqrt(sigma500/nn[5])
  if (m[1]<m[3]) {a1<--a1}
  m5<-as.matrix(c(m[5]+2*a1,m[5],m[5]-2*a1))
  ########first order genetic parameters###############
  hh13<-matrix(c(1,1,1,1,1,1,1,1,1,1,0,-1,1,0,-1,1,0,-1,
                 0,1,0,0,1,0,0,0.5,0,1,0,-1,0,0,0,0,0,0,
                 0,1,0,0.5,0.5,0.5,0.25,0.25,0.25),9,5)
  mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m5[1],m5[2],m5[3]))
  B13<-solve(crossprod(hh13,hh13))%*%crossprod(hh13,mm)
  a1<-B13[2];a2<-B13[3]
  a1<-(0.5*a1^2+0.25*a2^2)/m_nf
  sigma4<-matrix(0,3,1);sigma5<-matrix(0,3,1)
  sigma4[1]<-sigma400/2;sigma4[2]<-sigma4[1]
  sigma4[3]<-sigma4[1];sigma5[1]<-sigma500/2
  sigma5[3]<-sigma5[1];sigma5[2]<-sigma5[1]+a1
  mi4<-mix_pi4[c(1:3)];mi5<-mix_pi5[c(1:3)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
    sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
  ##########iteration process###########
  ############E-step###################
  ############F2 process###############
  iteration <- 0; stopa <- 1000
  W4 <- matrix(0,3,m_sam[4]); swx24 <- matrix(0,3,1)
  W5 <- matrix(0,3,m_sam[5]); swx25 <- matrix(0,3,1)
  hh<-matrix(0,4,4);b_line<-matrix(0,4,1)
  n0<-matrix(0,6,1);s0<-matrix(0,2,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:d1) { W4[i,] <- mix_pi4[i]*dnorm(dataF2,m4[i],sqrt(sigma4[i]))/dmixnorm(dataF2,m4,sqrt(sigma4),mix_pi4) }
    mi4 <- as.matrix(rowSums(W4)/m_sam[4])
    n0[c(1:3)]<-as.matrix(rowSums(W4))
    sumwx4 <- W4%*%dataF2
    #############F3 process############
    for(i in 1:d1) { W5[i,] <- mix_pi5[i]*dnorm(dataF3,m5[i],sqrt(sigma5[i]))/dmixnorm(dataF3,m5,sqrt(sigma5),mix_pi5) }
    mi5 <- as.matrix(rowSums(W5)/m_sam[5])
    n0[c(4:6)]<-as.matrix(rowSums(W5));n0[c(1:6)][abs(n0[c(1:6)])<0.000001]<-0.000001
    sumwx5 <- W5%*%dataF3
    aaa0<-0;aaa1<-1000;n_iter<-0;AA<-matrix(0,4,1)
    while(aaa1>0.0001){
      n_iter<-n_iter+1 
      ########first order genetic parameters###############
      mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m5[1],m5[2],m5[3]))
      B131<-solve(crossprod(hh13,hh13))%*%crossprod(hh13,mm) 
      aa1<-B131[2];aa2<-B131[3];aa1<-(0.5*aa1^2+0.25*aa2^2)/m_nf
      sigma5[2]<-sigma5[1]+aa1
      ############restrictions######################
      hh[1,1]<-sigma*(1/nn[1]+4/nn[2]+1/nn[3])+sigma4[1]/n0[1]+4*sigma4[2]/n0[2]+sigma4[3]/n0[3]
      hh[1,2]<-sigma*2/nn[2]+2*sigma4[1]/n0[1]+2*sigma4[2]/n0[2]
      hh[1,3]<-sigma*(1/nn[1]+1/nn[3])-4*sigma4[2]/n0[2]
      hh[1,4]<--sigma4[1]/n0[1]+sigma4[3]/n0[3]
      hh[2,2]<-sigma/nn[2]+4*sigma4[1]/n0[1]+sigma4[2]/n0[2]+4*sigma5[1]/n0[4]
      hh[2,3]<--2*sigma4[2]/n0[2]
      hh[2,4]<--2*sigma4[1]/n0[1]-2*sigma5[1]/n0[4]
      hh[3,3]<-sigma*(1/nn[1]+1/nn[3])+4*sigma4[2]/n0[2]+16*sigma5[2]/n0[5]
      hh[3,4]<-0
      hh[4,4]<-sigma4[1]/n0[1]+sigma4[3]/n0[3]+sigma5[1]/n0[4]+sigma5[3]/n0[6]
      for(i in 2:4)
      {
        for(j in 1:(i-1))
        {
          hh[i,j]<-hh[j,i]
        }
      }
      ##########################################################
      b_line[1]<-sumx[1]/nn[1]+2*sumx[2]/nn[2]+sumx[3]/nn[3]-sumwx4[1]/n0[1]-2*sumwx4[2]/n0[2]-sumwx4[3]/n0[3]
      b_line[2]<-sumx[2]/nn[2]-2*sumwx4[1]/n0[1]-sumwx4[2]/n0[2]+2*sumwx5[1]/n0[4]
      b_line[3]<-sumx[1]/nn[1]+sumx[3]/nn[3]+2*sumwx4[2]/n0[2]-4*sumwx5[2]/n0[5]
      b_line[4]<-sumwx4[1]/n0[1]-sumwx4[3]/n0[3]-sumwx5[1]/n0[4]+sumwx5[3]/n0[6]
      B1311<-solve(hh,b_line)
      ###################################################
      m[1]<-(sumx[1]-sigma*(B1311[1]+B1311[3]))/nn[1]
      m[2]<-(sumx[2]-sigma*(2*B1311[1]+B1311[2]))/nn[2]
      m[3]<-(sumx[3]-sigma*(B1311[1]+B1311[3]))/nn[3]
      m4[1]<-(sumwx4[1]+sigma4[1]*(B1311[1]+2*B1311[2]-B1311[4]))/n0[1]
      m4[2]<-(sumwx4[2]+sigma4[2]*(2*B1311[1]+B1311[2]-2*B1311[3]))/n0[2]
      m4[3]<-(sumwx4[3]+sigma4[3]*(B1311[1]+B1311[4]))/n0[3]
      m5[1]<-(sumwx5[1]+sigma5[1]*(-2*B1311[2]+B1311[4]))/n0[4]
      m5[2]<-(sumwx5[2]+sigma5[2]*4*B1311[3])/n0[5]
      m5[3]<-(sumwx5[3]-sigma5[3]*B1311[4])/n0[6]
      aaa1<-max(abs(B1311-AA))
      AA<-B1311
      if (n_iter>20) break
    }
    ########obtain variance#######################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d1) {swx24[i] <- W4[i,]%*%(dataF2-m4[i])^2 } ;for(i in 1:d1) {swx25[i] <- W5[i,]%*%(dataF3-m5[i])^2 }      
    sigma4[1]<-(swx24[1]+swx24[2]+swx24[3])/nn[4];sigma40<-sigma4[1]-sigma
    if (sigma40<0) {sigma40<-0;sigma4[1]<-sigma}
    sigma4[1]<-sigma40+sigma;sigma4[2]<-sigma4[1];sigma4[3]<-sigma4[1]
    aaa0<-sigma5[1]
    mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m5[1],m5[2],m5[3]))
    B13111<-solve(crossprod(hh13,hh13))%*%crossprod(hh13,mm) 
    aa1<-B13111[2];aa2<-B13111[3];aa1<-(0.5*aa1^2+0.25*aa2^2)/m_nf
    n_iter<-0;aa3<-1000
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      ab3<-sigma5[1]/(sigma5[1]+aa1)
      sigma5[1]<-(swx25[1]+swx25[3]+ab3^2*swx25[2])/(n0[4]+n0[6]+ab3*n0[5])
      aa3<-abs(sigma5[1]-aaa0)
      aaa0<-sigma5[1]
      if (n_iter>20) break
    } 
    sigma50<-sigma5[1]-sigma/m_nf
    if (sigma50<0) {sigma50<-0;sigma5[1]-sigma/m_nf}
    sigma5[1]<-sigma50+sigma/m_nf;sigma5[2]<-sigma5[1]+aa1
    sigma5[3]<-sigma5[1]
    s0[1]<-ss1+ss2+ss3;s0[2]<-nn[1]+nn[2]+nn[3]
    aaa0<-0;n_iter<-0
    abb1<-sigma40;abb2<-sigma50;aa3<-1000
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      if (abb1<0) {abb1<-0}
      if (abb2<0) {abb2<-0}
      abc1<-sigma/(sigma+abb1)
      abc2<-(sigma/m_nf)/(sigma/m_nf+abb2)
      abc3<-(sigma/m_nf)/(sigma/m_nf+abb2+aa1)
      aa4<-s0[1]+abc1^2*(swx24[1]+swx24[2]+swx24[3])+m_nf*(abc2^2*(swx25[1]+swx25[3])+abc3^2*swx25[2])
      aa5<-s0[2]+abc1*nn[4]+abc2*(n0[4]+n0[6])+abc3*n0[5]
      sigma<-aa4/aa5
      aa3<-abs(sigma-aaa0)
      aaa0<-sigma
      if (n_iter>20) break
    } 
    if (sigma<sigma0) {sigma<-sigma0}
    sigma4[1]<-sigma+sigma40;sigma4[2]<-sigma4[1];sigma4[3]<-sigma4[1]
    sigma5[1]<-sigma/m_nf+sigma50;sigma5[3]<-sigma5[1];sigma5[2]<-sigma/m_nf+sigma50+aa1
    ##############################criteria for iterations to stop######################## 
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
      sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*8
  ########first order genetic parameters###############  
  mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m5[1],m5[2],m5[3]))
  B131111<-solve(crossprod(hh13,hh13))%*%crossprod(hh13,mm)     
  ########second order genetic parameters###############   
  F2_jj <- sigmaF2 - sigma4[1]
  F2_gg <- sigma4[1]-sigma
  if(F2_jj<0) {F2_jj<-0}
  if(F2_gg<0 || F2_gg>sigmaF2) {F2_gg<-0}
  F2_ll <- F2_jj/sigmaF2
  F2_rr <- F2_gg/sigmaF2
  F3_jj <- sigmaF3 - sigma5[1]
  F3_gg <- sigma5[1]-sigma/m_nf
  if(F3_jj<0) {F3_jj<-0}
  if(F3_gg<0 || F3_gg>sigmaF3) {F3_gg<-0}
  F3_ll <- F3_jj/sigmaF3
  F3_rr <- F3_gg/sigmaF3
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma))
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
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma))
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
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma))
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
    F2gg <- (dataF2 - m4[i])/sqrt(sigma4[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mi4[i]
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
  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3); 
  F3w1<-1/(12*m_sam[5])
  F3bmw <- matrix(0,m_sam[5],1); F3bmwsl <- matrix(0,m_sam[5],d1)
  for(i in 1:d1){
    F3gg <- (dataF3 - m5[i])/sqrt(sigma5[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mi5[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[5]){F3P2 <- F3P2+runif(m_sam[5])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[5]) + sum((F3P2 - (as.matrix(c(1:m_sam[5])) - 0.5)/m_sam[5])^2)
  F3u <- as.matrix(c(12*m_sam[5]*((F3dd[1]/m_sam[5]-0.5)^2),((45*m_sam[5])/4)*((F3dd[2]/m_sam[5]-1/3)^2),180*m_sam[5]*((F3dd[3]/m_sam[5]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[5]))))
  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("MX1-AD-AD",round(abc,4),round(AIC,4),round(m[1],4),round(m[2],4),round(m[3],4), round(t(m4),4)," "," "," "," "," "," ",round(sigma4[1],4),
                       round(t(mi4),4)," "," "," "," "," "," ",round(t(m5),4)," "," "," "," "," "," ",round(t(sigma5),4)," "," "," "," "," "," ",round(t(mi5),4)," "," "," "," "," "," ",round(sigma,4),          
                       round(B131111[1],4)," "," "," "," ",round(B131111[2],4)," ",round(B131111[3],4)," "," "," "," "," ",round(B131111[4],4),round(B131111[5],4),round(F2_jj,4),round(F2_ll*100,4),round(F2_gg,4),round(F2_rr*100,4),round(F3_jj,4),round(F3_ll*100,4),round(F3_gg,4),round(F3_rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
 output<-as.matrix(output)
  OUTPUT<-list(output,mix_pi4,mix_pi5)
  return(OUTPUT)
}
####################MX1-A-AD(D-2)###############################################
G5ModelFun[[15]] <- function(K1,logL,df11,df21,df31,df41,df51,G5text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataF2 <- as.matrix(as.numeric(df41[,1]))
  dataF3 <- as.matrix(as.numeric(df51[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2),sum(dataF3)))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2),mean(dataF3)))
  sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1);sigmaF2<-var(dataF2);sigmaF3<-var(dataF3)
  sumx6<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  sigma<-sumx6;sigma0<-sigma;sigma40<-sigmaF2;sigma400<-sigma40;sigma50<-sigmaF3;sigma500<-sigma50
  m_nf <- as.numeric(G5text2)
  ###############procedure start###########################
  d1<-3
  mix_pi4<-as.matrix(c(0.25,0.5,0.25));mix_pi5<-as.matrix(c(0.25,0.5,0.25))
  nn<-m_sam
  a1<-sqrt(sigma400/nn[4])
  if (m[1]<m[3]) {a1<--a1}
  m4<-as.matrix(c(m[4]+2*a1,m[4],m[4]-2*a1))
  a1<-sqrt(sigma500/nn[5])
  if (m[1]<m[3]) {a1<--a1}
  m5<-as.matrix(c(m[5]+2*a1,m[5],m[5]-2*a1))
  ########first order genetic parameters###############
  hh14<-matrix(c(1,1,1,1,1,1,1,1,1,1,0,-1,1,0,-1,1,0,-1,
                 1,0,-1,0,0,0,0,0,0,0,1,0,0.5,0.5,0.5,0.25,0.25,0.25),9,4)
  mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m5[1],m5[2],m5[3]))
  B14<-solve(crossprod(hh14,hh14))%*%crossprod(hh14,mm)
  a1<-B14[2];a1<-(0.5*a1^2)/m_nf
  sigma4<-matrix(0,3,1);sigma5<-matrix(0,3,1)
  sigma4[1]<-sigma400/2;sigma4[2]<-sigma4[1];sigma4[3]<-sigma4[1]
  sigma5[1]<-sigma500/2;sigma5[3]<-sigma5[1];sigma5[2]<-sigma5[1]+a1
  mi4<-mix_pi4[c(1:3)];mi5<-mix_pi5[c(1:3)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
    sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
  ##########iteration process###########
  ############E-step###################
  ############F2 process###############
  iteration <- 0; stopa <- 1000
  W4 <- matrix(0,3,m_sam[4]); swx24 <- matrix(0,3,1)
  W5 <- matrix(0,3,m_sam[5]); swx25 <- matrix(0,3,1)
  hh<-matrix(0,5,5);b_line<-matrix(0,5,1)
  n0<-matrix(0,6,1); s0<-matrix(0,2,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:d1) { W4[i,] <- mix_pi4[i]*dnorm(dataF2,m4[i],sqrt(sigma4[i]))/dmixnorm(dataF2,m4,sqrt(sigma4),mix_pi4) }
    mi4 <- as.matrix(rowSums(W4)/m_sam[4])
    n0[c(1:3)]<-as.matrix(rowSums(W4))
    sumwx4 <- W4%*%dataF2
    #############F3 process############
    for(i in 1:d1) { W5[i,] <- mix_pi5[i]*dnorm(dataF3,m5[i],sqrt(sigma5[i]))/dmixnorm(dataF3,m5,sqrt(sigma5),mix_pi5) }
    mi5 <- as.matrix(rowSums(W5)/m_sam[5])
    n0[c(4:6)]<-as.matrix(rowSums(W5)); n0[c(1:6)][abs(n0[c(1:6)])<0.000001]<-0.000001
    sumwx5 <- W5%*%dataF3
    aaa0<-0;aaa1<-1000; n_iter<-0;AA<-matrix(0,5,1)
    while(aaa1>0.0001){
      n_iter<-n_iter+1 
      ########first order genetic parameters###############
      mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m5[1],m5[2],m5[3]))
      B141<-solve(crossprod(hh14,hh14))%*%crossprod(hh14,mm) 
      aa1<-B141[2];aa1<-(0.5*aa1^2)/m_nf
      sigma5[2]<-sigma5[1]+aa1
      ############restrictions######################
      hh[1,1]<-sigma*(1/nn[1]+4/nn[2]+1/nn[3])+16*sigma4[2]/n0[2]
      hh[1,2]<-sigma*2/nn[2]+4*sigma4[2]/n0[2]
      hh[1,3]<-sigma*2/nn[2]+12*sigma4[2]/n0[2]
      hh[1,4]<-sigma*2/nn[2]+20*sigma4[2]/n0[2]
      hh[1,5]<-8*sigma4[2]/n0[2]
      hh[2,2]<-sigma*1/nn[2]+4*sigma4[1]/n0[1]+sigma4[2]/n0[2]+4*sigma5[1]/n0[4]
      hh[2,3]<-sigma/nn[2]+3*sigma4[2]/n0[2]
      hh[2,4]<-sigma/nn[2]-4*sigma4[1]/n0[1]+5*sigma4[2]/n0[2]
      hh[2,5]<-sigma4[1]*(-2/n0[1]+2/n0[2])
      hh[3,3]<-sigma/nn[2]+9*sigma4[2]/n0[2]+4*sigma5[2]/n0[5]
      hh[3,4]<-sigma/nn[2]+15*sigma4[2]/n0[2]
      hh[3,5]<-6*sigma4[2]/n0[2]
      hh[4,4]<-sigma/nn[2]+4*sigma4[1]/n0[1]+25*sigma4[2]/n0[2]+4*sigma5[3]/n0[6]
      hh[4,5]<-sigma4[1]*2/n0[1]+10*sigma4[2]/n0[2]
      hh[5,5]<-sigma4[1]*(1/n0[1]+4/n0[2]+1/n0[3])
      for(i in 2:5)
      {
        for(j in 1:(i-1))
        {
          hh[i,j]<-hh[j,i]
        }
      }
      ##########################################################
      b_line[1]<-sumx[1]/nn[1]+2*sumx[2]/nn[2]+sumx[3]/nn[3]-4*sumwx4[2]/n0[2]
      b_line[2]<-sumx[2]/nn[2]-2*sumwx4[1]/n0[1]-sumwx4[2]/n0[2]+2*sumwx5[1]/n0[4]
      b_line[3]<-sumx[2]/nn[2]-3*sumwx4[2]/n0[2]+2*sumwx5[2]/n0[5]
      b_line[4]<-sumx[2]/nn[2]+2*sumwx4[1]/n0[1]-5*sumwx4[2]/n0[2]+2*sumwx5[3]/n0[6]
      b_line[5]<-sumwx4[1]/n0[1]-2*sumwx4[2]/n0[2]+sumwx4[3]/n0[3]
      B1411<-solve(hh,b_line)
      ##########################################################
      m[1]<-(sumx[1]-sigma*B1411[1])/nn[1]
      m[2]<-(sumx[2]-sigma*(2*B1411[1]+B1411[2]+B1411[3]+B1411[4]))/nn[2]
      m[3]<-(sumx[3]-sigma*B1411[1])/nn[3]
      m4[1]<-(sumwx4[1]+sigma4[1]*(2*B1411[2]-2*B1411[4]-B1411[5]))/n0[1]
      m4[2]<-(sumwx4[2]+sigma4[2]*(4*B1411[1]+B1411[2]+3*B1411[3]+5*B1411[4]+2*B1411[5]))/n0[2]
      m4[3]<-(sumwx4[3]-sigma4[3]*B1411[5])/n0[3]
      m5[1]<-(sumwx5[1]-sigma5[1]*2*B1411[2])/n0[4]
      m5[2]<-(sumwx5[2]-sigma5[2]*2*B1411[3])/n0[5]
      m5[3]<-(sumwx5[3]-sigma5[3]*2*B1411[4])/n0[6]
      aaa1<-max(abs(B1411-AA))
      AA<-B1411
      if (n_iter>20) break
    }
    ########obtain variance#######################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d1) {swx24[i] <- W4[i,]%*%(dataF2-m4[i])^2 } ;for(i in 1:d1) {swx25[i] <- W5[i,]%*%(dataF3-m5[i])^2 }      
    sigma4[1]<-(swx24[1]+swx24[2]+swx24[3])/nn[4]
    sigma40<-sigma4[1]-sigma
    if (sigma40<0) {sigma40<-0;sigma4[1]<-sigma}
    sigma4[1]<-sigma40+sigma;sigma4[2]<-sigma4[1];sigma4[3]<-sigma4[1]
    aaa0<-sigma5[1]
    mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m5[1],m5[2],m5[3]))
    B14111<-solve(crossprod(hh14,hh14))%*%crossprod(hh14,mm) 
    aa1<-B14111[2];aa1<-(0.5*aa1^2)/m_nf
    n_iter<-0;aa3<-1000
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      ab3<-sigma5[1]/(sigma5[1]+aa1)
      sigma5[1]<-(swx25[1]+swx25[3]+ab3^2*swx25[2])/(n0[4]+n0[6]+ab3*n0[5])
      aa3<-abs(sigma5[1]-aaa0)
      aaa0<-sigma5[1]
      if (n_iter>20) break
    } 
    sigma50<-sigma5[1]-sigma/m_nf
    if (sigma50<0) {sigma50<-0;sigma5[1]-sigma/m_nf}
    sigma5[1]<-sigma50+sigma/m_nf;sigma5[2]<-sigma5[1]+aa1;sigma5[3]<-sigma5[1]
    s0[1]<-ss1+ss2+ss3;s0[2]<-nn[1]+nn[2]+nn[3]
    aaa0<-0;n_iter<-0
    abb1<-sigma40;abb2<-sigma50
    aa3<-1000
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      abc1<-sigma/(sigma+sigma40)
      abc2<-(sigma/m_nf)/(sigma/m_nf+sigma50)
      abc3<-(sigma/m_nf)/(sigma/m_nf+sigma50+aa1)
      aa4<-s0[1]+abc1^2*(swx24[1]+swx24[2]+swx24[3])+m_nf*(abc2^2*(swx25[1]+swx25[3])+abc3^2*swx25[2])
      aa5<-s0[2]+abc1*nn[4]+abc2*(n0[4]+n0[6])+abc3*n0[5]
      sigma<-aa4/aa5
      aa3<-abs(sigma-aaa0)
      aaa0<-sigma
      if (n_iter>20) break
    } 
    if (sigma<sigma0) {sigma<-sigma0}
    sigma4[1]<-sigma+sigma40;sigma4[2]<-sigma4[1];sigma4[3]<-sigma4[1]
    sigma5[1]<-sigma/m_nf+sigma50;sigma5[3]<-sigma5[1];sigma5[2]<-sigma/m_nf+sigma50+aa1
    ##############################criteria for iterations to stop######################## 
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
      sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*7
  ########first order genetic parameters###############  
  mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m5[1],m5[2],m5[3]))
  B141111<-solve(crossprod(hh14,hh14))%*%crossprod(hh14,mm)     
  ########second order genetic parameters###############   
  F2_jj <- sigmaF2 - sigma4[1]
  F2_gg <- sigma4[1]-sigma
  if(F2_jj<0) {F2_jj<-0}
  if(F2_gg<0 || F2_gg>sigmaF2) {F2_gg<-0}
  F2_ll <- F2_jj/sigmaF2
  F2_rr <- F2_gg/sigmaF2
  F3_jj <- sigmaF3 - sigma5[1]
  F3_gg <- sigma5[1]-sigma/m_nf
  if(F3_jj<0) {F3_jj<-0}
  if(F3_gg<0 || F3_gg>sigmaF3) {F3_gg<-0}
  F3_ll <- F3_jj/sigmaF3
  F3_rr <- F3_gg/sigmaF3
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma))
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
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma))
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
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma))
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
    F2gg <- (dataF2 - m4[i])/sqrt(sigma4[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mi4[i]
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
  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3); 
  F3w1<-1/(12*m_sam[5])
  F3bmw <- matrix(0,m_sam[5],1); F3bmwsl <- matrix(0,m_sam[5],d1)
  for(i in 1:d1){
    F3gg <- (dataF3 - m5[i])/sqrt(sigma5[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mi5[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[5]){F3P2 <- F3P2+runif(m_sam[5])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[5]) + sum((F3P2 - (as.matrix(c(1:m_sam[5])) - 0.5)/m_sam[5])^2)
  F3u <- as.matrix(c(12*m_sam[5]*((F3dd[1]/m_sam[5]-0.5)^2),((45*m_sam[5])/4)*((F3dd[2]/m_sam[5]-1/3)^2),180*m_sam[5]*((F3dd[3]/m_sam[5]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[5]))))
  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("MX1-A-AD",round(abc,4),round(AIC,4),round(m[1],4),round(m[2],4),round(m[3],4), round(t(m4),4)," "," "," "," "," "," ",round(sigma4[1],4),
                       round(t(mi4),4)," "," "," "," "," "," ",round(t(m5),4)," "," "," "," "," "," ",round(t(sigma5),4)," "," "," "," "," "," ",round(t(mi5),4)," "," "," "," "," "," ",round(sigma,4),          
                       round(B141111[1],4)," "," "," "," ",round(B141111[2],4)," "," "," "," "," "," "," ",round(B141111[3],4),round(B141111[4],4),round(F2_jj,4),round(F2_ll*100,4),round(F2_gg,4),round(F2_rr*100,4),round(F3_jj,4),round(F3_ll*100,4),round(F3_gg,4),round(F3_rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
 output<-as.matrix(output)
  OUTPUT<-list(output,mix_pi4,mix_pi5)
  return(OUTPUT)
}
####################MX1-EAD-AD(D-3)###############################################
G5ModelFun[[16]] <- function(K1,logL,df11,df21,df31,df41,df51,G5text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataF2 <- as.matrix(as.numeric(df41[,1]))
  dataF3 <- as.matrix(as.numeric(df51[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2),sum(dataF3)))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2),mean(dataF3)))
  sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1);sigmaF2<-var(dataF2);sigmaF3<-var(dataF3)
  sumx6<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  sigma<-sumx6;sigma0<-sigma;sigma40<-sigmaF2;sigma400<-sigma40;sigma50<-sigmaF3;sigma500<-sigma50
  m_nf <- as.numeric(G5text2)
  ###############procedure start###########################
  d1<-3;d2<-2
  mix_pi4<-as.matrix(c(0.75,0.25));mix_pi5<-as.matrix(c(0.25,0.5,0.25))
  sigma4<-matrix(0,2,1);sigma5<-matrix(0,3,1)
  abb<-sigma400/(6*sigma0)
  sigma4[1]<-sigma400/abb;sigma4[2]<-sigma4[1]
  sigma5[1]<-sigma500/abb;sigma5[3]<-sigma5[1]
  nn<-m_sam
  a1<-sqrt(sigma400/nn[4])
  if (m[1]<m[3]) {a1<--a1}
  m4<-as.matrix(c(m[4]+a1,m[4]-a1))
  a1<-sqrt(sigma500/nn[5])
  if (m[1]<m[3]) {a1<--a1}
  m5<-as.matrix(c(m[5]+2.5*a1,m[5]+0.5*a1,m[5]-2.5*a1))
  ########first order genetic parameters###############
  hh15<-matrix(c(1,1,1,1,1,1,1,1,1,1,-1,1,-1,1,0.5,-1,
                 1,0,-1,0,0,0,0,0,0,1,0,0.5,0.5,0.25,0.25,0.25),8,4)
  mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m5[1],m5[2],m5[3]))
  B15<-solve(crossprod(hh15,hh15))%*%crossprod(hh15,mm)
  a1<-B15[2];a1<-0.75*a1^2
  sigma5[2]<-sigma5[1]+a1/m_nf
  mi4<-mix_pi4[c(1:2)];mi5<-mix_pi5[c(1:3)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
    sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
  ##########iteration process###########
  ############E-step###################
  ############F2 process###############
  iteration <- 0; stopa <- 1000
  W4 <- matrix(0,2,m_sam[4]); swx24 <- matrix(0,2,1)
  W5 <- matrix(0,3,m_sam[5]); swx25 <- matrix(0,3,1)
  hh<-matrix(0,4,4);b_line<-matrix(0,4,1)
  n0<-matrix(0,6,1);s0<-matrix(0,2,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:d2) { W4[i,] <- mix_pi4[i]*dnorm(dataF2,m4[i],sqrt(sigma4[i]))/dmixnorm(dataF2,m4,sqrt(sigma4),mix_pi4) }
    mi4 <- as.matrix(rowSums(W4)/m_sam[4]);n0[c(1:2)]<-as.matrix(rowSums(W4));sumwx4 <- W4%*%dataF2
    #############F3 process############
    for(i in 1:d1) { W5[i,] <- mix_pi5[i]*dnorm(dataF3,m5[i],sqrt(sigma5[i]))/dmixnorm(dataF3,m5,sqrt(sigma5),mix_pi5) }
    mi5 <- as.matrix(rowSums(W5)/m_sam[5]);n0[c(3:5)]<-as.matrix(rowSums(W5))
    n0[c(1:5)][abs(n0[c(1:5)])<0.000001]<-0.000001;sumwx5 <- W5%*%dataF3
    aaa0<-0
    aaa1<-1000;n_iter<-0;AA<-matrix(0,4,1)
    while(aaa1>0.0001){
      n_iter<-n_iter+1 
      ########first order genetic parameters###############
      mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m5[1],m5[2],m5[3]))
      B151<-solve(crossprod(hh15,hh15))%*%crossprod(hh15,mm) 
      aa1<-B151[2]
      aa1<-(0.75*aa1^2)/m_nf
      sigma5[2]<-sigma5[1]+aa1
      ############restrictions######################
      hh[1,1]<-sigma*(1/nn[1]+4/nn[2]+1/nn[3])+9*sigma4[1]/n0[1]+sigma4[2]/n0[2]
      hh[1,2]<-sigma*2/nn[2]+9*sigma4[1]/n0[1]
      hh[1,3]<-sigma*(1/nn[1]+1/nn[3])-6*sigma4[1]/n0[1]
      hh[1,4]<-sigma*(2/nn[1]+6/nn[2]+2/nn[3])+15*sigma4[1]/n0[1]
      hh[2,2]<-sigma/nn[2]+9*sigma4[1]/n0[1]+4*sigma5[1]/n0[3]
      hh[2,3]<--6*sigma4[1]/n0[1]
      hh[2,4]<-sigma*3/nn[2]+15*sigma4[1]/n0[1]
      hh[3,3]<-sigma*(1/nn[1]+1/nn[3])+4*sigma4[1]/n0[1]+16*sigma5[2]/n0[4]
      hh[3,4]<-sigma*(2/nn[1]+2/nn[3])-10*sigma4[1]/n0[1]
      hh[4,4]<-sigma*(4/nn[1]+9/nn[2]+4/nn[3])+25*sigma4[1]/n0[1]+4*sigma5[3]/n0[5]
      for(i in 2:4)
      {
        for(j in 1:(i-1))
        {
          hh[i,j]<-hh[j,i]
        }
      }
      ##########################################################
      b_line[1]<-sumx[1]/nn[1]+2*sumx[2]/nn[2]+sumx[3]/nn[3]-3*sumwx4[1]/n0[1]-sumwx4[2]/n0[2]
      b_line[2]<-sumx[2]/nn[2]-3*sumwx4[1]/n0[1]+2*sumwx5[1]/n0[3]
      b_line[3]<-sumx[1]/nn[1]+sumx[3]/nn[3]+2*sumwx4[1]/n0[1]-4*sumwx5[2]/n0[4]
      b_line[4]<-2*sumx[1]/nn[1]+3*sumx[2]/nn[2]+2*sumx[3]/nn[3]-5*sumwx4[1]/n0[1]-2*sumwx5[3]/n0[5]
      B1511<-solve(hh,b_line)
      ###################################################
      m[1]<-(sumx[1]-sigma*(B1511[1]+B1511[3]+2*B1511[4]))/nn[1]
      m[2]<-(sumx[2]-sigma*(2*B1511[1]+B1511[2]+3*B1511[4]))/nn[2]
      m[3]<-(sumx[3]-sigma*(B1511[1]+B1511[3]+2*B1511[4]))/nn[3]
      m4[1]<-(sumwx4[1]+sigma4[1]*(3*B1511[1]+3*B1511[2]-2*B1511[3]+5*B1511[4]))/n0[1]
      m4[2]<-(sumwx4[2]+sigma4[2]*B1511[1])/n0[2]
      m5[1]<-(sumwx5[1]-sigma5[1]*2*B1511[2])/n0[3]
      m5[2]<-(sumwx5[2]+sigma5[2]*4*B1511[3])/n0[4]
      m5[3]<-(sumwx5[3]+sigma5[3]*2*B1511[4])/n0[5]
      aaa1<-max(abs(B1511-AA))
      AA<-B1511
      if (n_iter>20) break
    }
    ########obtain variance#######################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d2) {swx24[i] <- W4[i,]%*%(dataF2-m4[i])^2 } ;for(i in 1:d1) {swx25[i] <- W5[i,]%*%(dataF3-m5[i])^2 }      
    sigma4[1]<-(swx24[1]+swx24[2])/nn[4];sigma40<-sigma4[1]-sigma
    if (sigma40<0) {sigma40<-0}
    sigma4[1]<-sigma+sigma40;sigma4[2]<-sigma4[1]
    mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m5[1],m5[2],m5[3]))
    B15111<-solve(crossprod(hh15,hh15))%*%crossprod(hh15,mm) 
    aa1<-B15111[2];aa1<-(0.75*aa1^2)/m_nf
    aaa0<-sigma5[1];n_iter<-0
    aa3<-1000
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      ab3<-sigma5[1]/(sigma5[1]+aa1)
      sigma5[1]<-(swx25[1]+swx25[3]+ab3^2*swx25[2])/(n0[3]+n0[5]+ab3*n0[4])
      aa3<-abs(sigma5[1]-aaa0)
      aaa0<-sigma5[1]
      if (n_iter>20) break
    } 
    sigma50<-sigma5[1]-sigma/m_nf
    if (sigma50<0) {sigma50<-0;sigma5[1]<-sigma/m_nf}
    sigma5[1]<-sigma/m_nf+sigma50;sigma5[2]<-sigma5[1]+aa1;sigma5[3]<-sigma5[1]
    s0[1]<-ss1+ss2+ss3;s0[2]<-nn[1]+nn[2]+nn[3]
    aaa0<-0;n_iter<-0;aa3<-1000
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      abc1<-sigma/(sigma+sigma40)
      abc2<-(sigma/m_nf)/(sigma/m_nf+sigma50)
      abc3<-(sigma/m_nf)/(sigma/m_nf+sigma50+aa1)
      aa4<-s0[1]+abc1^2*(swx24[1]+swx24[2])+m_nf*(abc2^2*(swx25[1]+swx25[3])+abc3^2*swx25[2])
      aa5<-s0[2]+abc1*nn[4]+abc2*(n0[3]+n0[5])+abc3*n0[4]
      sigma<-aa4/aa5
      aa3<-abs(sigma-aaa0)
      aaa0<-sigma
      if (n_iter>20) break
    } 
    if (sigma<sigma0) {sigma<-sigma0}
    sigma4[1]<-sigma+sigma40;sigma4[2]<-sigma4[1];sigma5[1]<-sigma/m_nf+sigma50
    sigma5[3]<-sigma5[1];sigma5[2]<-sigma5[1]+aa1
    ##############################criteria for iterations to stop######################## 
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
      sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*7
  ########first order genetic parameters###############  
  mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m5[1],m5[2],m5[3]))
  B151111<-solve(crossprod(hh15,hh15))%*%crossprod(hh15,mm)     
  ########second order genetic parameters###############   
  F2_jj <- sigmaF2 - sigma4[1]
  F2_gg <- sigma4[1]-sigma
  if(F2_jj<0) {F2_jj<-0}
  if(F2_gg<0 || F2_gg>sigmaF2) {F2_gg<-0}
  F2_ll <- F2_jj/sigmaF2
  F2_rr <- F2_gg/sigmaF2
  F3_jj <- sigmaF3 - sigma5[1]
  F3_gg <- sigma5[1]-sigma/m_nf
  if(F3_jj<0) {F3_jj<-0}
  if(F3_gg<0 || F3_gg>sigmaF3) {F3_gg<-0}
  F3_ll <- F3_jj/sigmaF3
  F3_rr <- F3_gg/sigmaF3
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma))
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
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma))
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
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma))
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
    F2gg <- (dataF2 - m4[i])/sqrt(sigma4[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mi4[i]
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
  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3); 
  F3w1<-1/(12*m_sam[5])
  F3bmw <- matrix(0,m_sam[5],1); F3bmwsl <- matrix(0,m_sam[5],d1)
  for(i in 1:d1){
    F3gg <- (dataF3 - m5[i])/sqrt(sigma5[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mi5[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[5]){F3P2 <- F3P2+runif(m_sam[5])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[5]) + sum((F3P2 - (as.matrix(c(1:m_sam[5])) - 0.5)/m_sam[5])^2)
  F3u <- as.matrix(c(12*m_sam[5]*((F3dd[1]/m_sam[5]-0.5)^2),((45*m_sam[5])/4)*((F3dd[2]/m_sam[5]-1/3)^2),180*m_sam[5]*((F3dd[3]/m_sam[5]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[5]))))
  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("MX1-EAD-AD",round(abc,4),round(AIC,4),round(m[1],4),round(m[2],4),round(m[3],4), round(t(m4),4)," "," "," "," "," "," "," ",round(sigma4[1],4),
                       round(t(mi4),4)," "," "," "," "," "," "," ",round(t(m5),4)," "," "," "," "," "," ",round(t(sigma5),4)," "," "," "," "," "," ",round(t(mi5),4)," "," "," "," "," "," ",round(sigma,4),          
                       round(B151111[1],4)," "," "," "," ",round(B151111[2],4)," "," "," "," "," "," "," ",round(B151111[3],4),round(B151111[4],4),round(F2_jj,4),round(F2_ll*100,4),round(F2_gg,4),round(F2_rr*100,4),round(F3_jj,4),round(F3_ll*100,4),round(F3_gg,4),round(F3_rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
 output<-as.matrix(output)
  OUTPUT<-list(output,mix_pi4,mix_pi5)
  return(OUTPUT)
}
########################MX1-NCD-AD(D-4)#############################
G5ModelFun[[17]] <- function(K1,logL,df11,df21,df31,df41,df51,G5text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataF2 <- as.matrix(as.numeric(df41[,1]))
  dataF3 <- as.matrix(as.numeric(df51[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2),sum(dataF3)))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2),mean(dataF3)))
  sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1);sigmaF2<-var(dataF2);sigmaF3<-var(dataF3)
  sumx6<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  sigma<-sumx6;sigma0<-sigma;sigma40<-sigmaF2;sigma400<-sigma40;sigma50<-sigmaF3;sigma500<-sigma50
  m_nf <- as.numeric(G5text2)
  ###############procedure start###########################
  d1<-3;d2<-2
  mix_pi4<-as.matrix(c(0.25,0.75));mix_pi5<-as.matrix(c(0.25,0.5,0.25))
  sigma4<-matrix(0,2,1);sigma5<-matrix(0,3,1)
  abb<-sigma400/(6*sigma0)
  sigma4[1]<-sigma400/abb;sigma4[2]<-sigma4[1];sigma5[1]<-sigma500/abb;sigma5[3]<-sigma5[1]
  nn<-m_sam
  a1<-sqrt(sigma400/nn[4])
  if (m[1]<m[3]) {a1<--a1}
  m4<-as.matrix(c(m[4]+a1,m[4]-a1))
  a1<-sqrt(sigma500/nn[5])
  if (m[1]<m[3]) {a1<--a1}
  m5<-as.matrix(c(m[5]+2.5*a1,m[5]+0.5*a1,m[5]-2.5*a1))
  ########first order genetic parameters###############
  hh16<-matrix(c(1,1,1,1,1,1,1,1,1,-1,-1,1,-1,1,-0.5,-1,
                 1,0,-1,0,0,0,0,0,0,1,0,0.5,0.5,0.25,0.25,0.25),8,4)
  mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m5[1],m5[2],m5[3]))
  B16<-solve(crossprod(hh16,hh16))%*%crossprod(hh16,mm)
  a1<-B16[2];a1<-0.75*a1^2
  sigma5[2]<-sigma5[1]+a1/m_nf
  mi4<-mix_pi4[c(1:2)];mi5<-mix_pi5[c(1:3)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
    sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
  ##########iteration process###########
  ############E-step###################
  ############F2 process###############
  iteration <- 0; stopa <- 1000
  W4 <- matrix(0,2,m_sam[4]); swx24 <- matrix(0,2,1)
  W5 <- matrix(0,3,m_sam[5]); swx25 <- matrix(0,3,1)
  hh<-matrix(0,4,4);b_line<-matrix(0,4,1)
  n0<-matrix(0,6,1);s0<-matrix(0,2,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:d2) { W4[i,] <- mix_pi4[i]*dnorm(dataF2,m4[i],sqrt(sigma4[i]))/dmixnorm(dataF2,m4,sqrt(sigma4),mix_pi4) }
    mi4 <- as.matrix(rowSums(W4)/m_sam[4])
    n0[c(1:2)]<-as.matrix(rowSums(W4));sumwx4 <- W4%*%dataF2
    #############F3 process############
    for(i in 1:d1) { W5[i,] <- mix_pi5[i]*dnorm(dataF3,m5[i],sqrt(sigma5[i]))/dmixnorm(dataF3,m5,sqrt(sigma5),mix_pi5) }
    mi5 <- as.matrix(rowSums(W5)/m_sam[5])
    n0[c(3:5)]<-as.matrix(rowSums(W5))
    n0[c(1:5)][abs(n0[c(1:5)])<0.000001]<-0.000001
    sumwx5 <- W5%*%dataF3
    aaa0<-0
    aaa1<-1000;n_iter<-0;AA<-matrix(0,4,1)
    while(aaa1>0.0001){
      n_iter<-n_iter+1 
      ########first order genetic parameters###############
      mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m5[1],m5[2],m5[3]))
      B161<-solve(crossprod(hh16,hh16))%*%crossprod(hh16,mm) 
      aa1<-B161[2];aa1<-(0.75*aa1^2)/m_nf
      sigma5[2]<-sigma5[1]+aa1
      ############restrictions######################
      hh[1,1]<-sigma*(1/nn[1]+4/nn[2]+1/nn[3])+sigma4[1]/n0[1]+9*sigma4[2]/n0[2]
      hh[1,2]<-sigma*2/nn[2]+2*sigma4[1]/n0[1]+3*sigma4[2]/n0[2]
      hh[1,3]<-sigma*(1/nn[1]+1/nn[3])-6*sigma4[2]/n0[2]
      hh[1,4]<-sigma*(2/nn[1]+6/nn[2]+2/nn[3])+2*sigma4[1]/n0[1]+9*sigma4[2]/n0[2]
      hh[2,2]<-sigma/nn[2]+4*sigma4[1]/n0[1]+sigma4[2]/n0[2]+4*sigma5[1]/n0[3]
      hh[2,3]<--2*sigma4[2]/n0[2]
      hh[2,4]<-sigma*3/nn[2]+4*sigma4[1]/n0[1]+3*sigma4[2]/n0[2]
      hh[3,3]<-sigma*(1/nn[1]+1/nn[3])+4*sigma4[2]/n0[2]+16*sigma5[2]/n0[4]
      hh[3,4]<-sigma*(2/nn[1]+2/nn[3])-6*sigma4[2]/n0[2]
      hh[4,4]<-sigma*(4/nn[1]+9/nn[2]+4/nn[3])+4*sigma4[1]/n0[1]+9*sigma4[2]/n0[2]+4*sigma5[3]/n0[5]
      for(i in 2:4)
      {
        for(j in 1:(i-1))
        {
          hh[i,j]<-hh[j,i]
        }
      }
      ##########################################################
      b_line[1]<-sumx[1]/nn[1]+2*sumx[2]/nn[2]+sumx[3]/nn[3]-sumwx4[1]/n0[1]-3*sumwx4[2]/n0[2]
      b_line[2]<-sumx[2]/nn[2]-2*sumwx4[1]/n0[1]-sumwx4[2]/n0[2]+2*sumwx5[1]/n0[3]
      b_line[3]<-sumx[1]/nn[1]+sumx[3]/nn[3]+2*sumwx4[2]/n0[2]-4*sumwx5[2]/n0[4]
      b_line[4]<-2*sumx[1]/nn[1]+3*sumx[2]/nn[2]+2*sumx[3]/nn[3]-2*sumwx4[1]/n0[1]-3*sumwx4[2]/n0[2]-2*sumwx5[3]/n0[5]
      B1611<-solve(hh,b_line)
      ##########################################################
      m[1]<-(sumx[1]-sigma*(B1611[1]+B1611[3]+2*B1611[4]))/nn[1]
      m[2]<-(sumx[2]-sigma*(2*B1611[1]+B1611[2]+3*B1611[4]))/nn[2]
      m[3]<-(sumx[3]-sigma*(B1611[1]+B1611[3]+2*B1611[4]))/nn[3]
      m4[1]<-(sumwx4[1]+sigma4[1]*(B1611[1]+2*B1611[2]+2*B1611[4]))/n0[1]
      m4[2]<-(sumwx4[2]+sigma4[2]*(3*B1611[1]+B1611[2]-2*B1611[3]+3*B1611[4]))/n0[2]
      m5[1]<-(sumwx5[1]-sigma5[1]*2*B1611[2])/n0[3]
      m5[2]<-(sumwx5[2]+sigma5[2]*4*B1611[3])/n0[4]
      m5[3]<-(sumwx5[3]+sigma5[3]*2*B1611[4])/n0[5]
      aaa1<-max(abs(B1611-AA))
      AA<-B1611
      if (n_iter>20) break
    }
    ########obtain variance#######################
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d2) {swx24[i] <- W4[i,]%*%(dataF2-m4[i])^2 } ;for(i in 1:d1) {swx25[i] <- W5[i,]%*%(dataF3-m5[i])^2 }      
    sigma4[1]<-(swx24[1]+swx24[2])/nn[4];sigma40<-sigma4[1]-sigma
    if (sigma40<0) {sigma40<-0;sigma4[1]<-sigma}
    sigma4[1]<-sigma+sigma40;sigma4[2]<-sigma4[1]
    aaa0<-sigma5[1]
    mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m5[1],m5[2],m5[3]))
    B16111<-solve(crossprod(hh16,hh16))%*%crossprod(hh16,mm) 
    aa1<-B16111[2];aa1<-(0.75*aa1^2)/m_nf
    n_iter<-0;aa3<-1000
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      ab3<-sigma5[1]/(sigma5[1]+aa1)
      sigma5[1]<-(swx25[1]+swx25[3]+ab3^2*swx25[2])/(n0[3]+n0[5]+ab3*n0[4])
      aa3<-abs(sigma5[1]-aaa0)
      aaa0<-sigma5[1]
      if (n_iter>20) break
    } 
    sigma50<-sigma5[1]-sigma/m_nf
    if (sigma50<0) {sigma50<-0;sigma5[1]-sigma/m_nf}
    sigma5[1]<-sigma/m_nf+sigma50;sigma5[2]<-sigma5[1]+aa1;sigma5[3]<-sigma5[1];
    s0[1]<-ss1+ss2+ss3;s0[2]<-nn[1]+nn[2]+nn[3]
    aaa0<-0
    n_iter<-0;aa3<-1000
    while (aa3>0.0001)
    {
      n_iter<-n_iter+1
      abc1<-sigma/(sigma+sigma40)
      abc2<-(sigma/m_nf)/(sigma/m_nf+sigma50)
      abc3<-(sigma/m_nf)/(sigma/m_nf+sigma50+aa1)
      aa4<-s0[1]+abc1^2*(swx24[1]+swx24[2])+m_nf*(abc2^2*(swx25[1]+swx25[3])+abc3^2*swx25[2])
      aa5<-s0[2]+abc1*nn[4]+abc2*(n0[3]+n0[5])+abc3*n0[4]
      sigma<-aa4/aa5
      aa3<-abs(sigma-aaa0)
      aaa0<-sigma
      if (n_iter>20) break
    } 
    if (sigma<sigma0) {sigma<-sigma0}
    sigma4[1]<-sigma+sigma40;sigma4[2]<-sigma4[1]
    sigma5[1]<-sigma/m_nf+sigma50;sigma5[3]<-sigma5[1];sigma5[2]<-sigma5[1]+aa1
    ##############################criteria for iterations to stop######################## 
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
      sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*7
  ########first order genetic parameters###############  
  mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m5[1],m5[2],m5[3]))
  B161111<-solve(crossprod(hh16,hh16))%*%crossprod(hh16,mm)     
  ########second order genetic parameters###############   
  F2_jj <- sigmaF2 - sigma4[1]
  F2_gg <- sigma4[1]-sigma
  if(F2_jj<0) {F2_jj<-0}
  if(F2_gg<0 || F2_gg>sigmaF2) {F2_gg<-0}
  F2_ll <- F2_jj/sigmaF2
  F2_rr <- F2_gg/sigmaF2
  F3_jj <- sigmaF3 - sigma5[1]
  F3_gg <- sigma5[1]-sigma/m_nf
  if(F3_jj<0) {F3_jj<-0}
  if(F3_gg<0 || F3_gg>sigmaF3) {F3_gg<-0}
  F3_ll <- F3_jj/sigmaF3
  F3_rr <- F3_gg/sigmaF3
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma))
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
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma))
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
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma))
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
    F2gg <- (dataF2 - m4[i])/sqrt(sigma4[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mi4[i]
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
  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3); 
  F3w1<-1/(12*m_sam[5])
  F3bmw <- matrix(0,m_sam[5],1); F3bmwsl <- matrix(0,m_sam[5],d1)
  for(i in 1:d1){
    F3gg <- (dataF3 - m5[i])/sqrt(sigma5[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mi5[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[5]){F3P2 <- F3P2+runif(m_sam[5])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[5]) + sum((F3P2 - (as.matrix(c(1:m_sam[5])) - 0.5)/m_sam[5])^2)
  F3u <- as.matrix(c(12*m_sam[5]*((F3dd[1]/m_sam[5]-0.5)^2),((45*m_sam[5])/4)*((F3dd[2]/m_sam[5]-1/3)^2),180*m_sam[5]*((F3dd[3]/m_sam[5]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[5]))))
  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("MX1-NCD-AD",round(abc,4),round(AIC,4),round(m[1],4),round(m[2],4),round(m[3],4), round(t(m4),4)," "," "," "," "," "," "," ",round(sigma4[1],4),
                       round(t(mi4),4)," "," "," "," "," "," "," ",round(t(m5),4)," "," "," "," "," "," ",round(t(sigma5),4)," "," "," "," "," "," ",round(t(mi5),4)," "," "," "," "," "," ",round(sigma,4),          
                       round(B161111[1],4)," "," "," "," ",round(B161111[2],4)," "," "," "," "," "," "," ",round(B161111[3],4),round(B161111[4],4),round(F2_jj,4),round(F2_ll*100,4),round(F2_gg,4),round(F2_rr*100,4),round(F3_jj,4),round(F3_ll*100,4),round(F3_gg,4),round(F3_rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
 output<-as.matrix(output)
  OUTPUT<-list(output,mix_pi4,mix_pi5)
  return(OUTPUT)
}

######################MX2-ADI-ADI(E-0)##############################
G5ModelFun[[18]] <- function(K1,logL,df11,df21,df31,df41,df51,G5text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataF2 <- as.matrix(as.numeric(df41[,1]))
  dataF3 <- as.matrix(as.numeric(df51[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2),sum(dataF3)))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2),mean(dataF3)))
  sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1);sigmaF2<-var(dataF2);sigmaF3<-var(dataF3)
  sumx6<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  sigma<-sumx6;sigma0<-sigma;sigma40<-sigmaF2;sigma400<-sigma40;sigma50<-sigmaF3;sigma500<-sigma50
  m_nf <- as.numeric(G5text2)
  ###############procedure start###########################
  mix_pi4<-as.matrix(c(0.0625,0.125,0.0625,0.1250,0.250,0.125,0.0625,0.125,0.0625))
  mix_pi5<-as.matrix(c(0.0625,0.125,0.0625,0.1250,0.250,0.125,0.0625,0.125,0.0625))
  a1<-sqrt(sigma400/m_sam[4])
  if (m[1]<m[3]) a1<--a1
  m4<-as.matrix(c(m[4]+2.8*a1,m[4]+2.1*a1,m[4]+1.4*a1,m[4]+0.7*a1,m[4],m[4]-0.7*a1,m[4]-1.4*a1,m[4]-2.1*a1,m[4]-2.8*a1))
  a1<-sqrt(sigma500/m_sam[5])
  if (m[1]<m[3]) {a1<--a1}
  m5<-as.matrix(c(m[5]+2.8*a1,m[5]+2.1*a1,m[5]+1.4*a1,m[5]+0.7*a1,m[5],m[5]-0.7*a1,m[5]-1.4*a1,m[5]-2.1*a1,m[5]-2.8*a1))
  nn<-m_sam
  abb<-sigma400/(sigma*5)
  sigma4<-matrix(0,9,1);sigma5<-matrix(0,9,1)
  sigma4[1]<-sigma400/abb;sigma5[1]<-sigma500/abb
  #######first order parameters############
  hh17<-matrix(c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                 0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,-1,1,1,1,0,0,0,-1,-1,-1,1,1,1,0,0,0,-1,-1,-1,
                 1,0,-1,1,0,-1,1,0,-1,1,0,-1,1,0,-1,1,0,-1,1,0,-1,0,1,0,0,0,0,1,1,1,0,0,0,0,0,0,0.5,0.5,0.5,0,0,0,
                 0,1,0,0,1,0,0,1,0,0,1,0,0,0.5,0,0,0.5,0,0,0.5,0,1,0,1,1,0,-1,0,0,0,-1,0,1,1,0,-1,0,0,0,-1,0,1,
                 0,0,0,0,1,0,0,0,0,0,-1,0,0,0.5,0,0,0,0,0,-0.5,0,0,0,0,0,0,0,1,0,-1,0,0,0,0,0,0,0.5,0,-0.5,0,0,0,
                 0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0.25,0,0,0,0),21,13)
  mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m4[4],m4[5],m4[6],m4[7],m4[8],m4[9],m5[1],m5[2],m5[3],m5[4],m5[5],m5[6],m5[7],m5[8],m5[9]))
  B17<-solve(crossprod(hh17,hh17))%*%crossprod(hh17,mm)
  g_aa1<-(0.5*(B17[7]+B17[10])^2+0.25*(B17[9]+B17[11])^2)/m_nf            
  #   0.5(db+i)**2+0.25(hb+jab)**2.
  g_aa2<-(0.5*(B17[6]+B17[10])^2+0.25*(B17[8]+B17[12])^2)/m_nf        
  #   0.5(da+i)**2+0.25(ha+jba)**2.
  g_aa3<-(0.5*(B17[6]-B17[10])^2+0.25*(B17[8]-B17[12])^2)/m_nf          
  #   0.5(da-i)**2+0.25(ha-jba)**2.
  g_aa4<-(0.5*(B17[7]-B17[10])^2+0.25*(B17[9]-B17[11])^2)/m_nf           
  #   0.5(db-i)**2+0.25(hb-jab)**2.
  g_aa5<-0.25*(B17[6]^2+B17[7]^2+B17[10]^2+(B17[6]+B17[11])^2+(B17[7]+B17[12])^2+(B17[8]+B17[13]/2)^2+(B17[9]+B17[13]/2)^2+B17[13]^2/4)/m_nf
  sigma4[c(2:9)]<-sigma4[1];sigma5[c(3,7,9)]<-sigma5[1]
  sigma5[2]<-sigma5[1]+g_aa1;sigma5[4]<-sigma5[1]+g_aa2
  sigma5[5]<-sigma5[1]+g_aa5;sigma5[6]<-sigma5[1]+g_aa3;sigma5[8]<-sigma5[1]+g_aa4
  mi4<-mix_pi4[c(1:9)];mi5<-mix_pi5[c(1:9)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
    sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
  ########################################################
  ##########iteration process###########
  ############E-step###################
  ############F2 process###############
  d3<-9
  iteration <- 0; stopa <- 1000
  W4 <- matrix(0,9,m_sam[4]); swx24 <- matrix(0,9,1)
  W5 <- matrix(0,9,m_sam[5]); swx25 <- matrix(0,9,1)
  hh<-matrix(0,8,8);b_line<-matrix(0,8,1)
  n0<-matrix(0,18,1);s0<-matrix(0,9,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:d3) { W4[i,] <- mix_pi4[i]*dnorm(dataF2,m4[i],sqrt(sigma4[i]))/dmixnorm(dataF2,m4,sqrt(sigma4),mix_pi4) }
    mi4 <- as.matrix(rowSums(W4)/m_sam[4]);n0[c(1:9)]<-as.matrix(rowSums(W4));sumwx4 <- W4%*%dataF2
    #############F3 process############
    for(i in 1:d3) { W5[i,] <- mix_pi5[i]*dnorm(dataF3,m5[i],sqrt(sigma5[i]))/dmixnorm(dataF3,m5,sqrt(sigma5),mix_pi5) }
    mi5 <- as.matrix(rowSums(W5)/m_sam[5])
    n0[c(10:18)]<-as.matrix(rowSums(W5));n0[c(1:18)][abs(n0[c(1:18)])<0.000001]<-0.000001
    sumwx5 <- W5%*%dataF3
    aaa0<-0
    aaa1<-1000;n_iter<-0;AA<-matrix(0,8,1)
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      ##########first order parameters###################
      mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m4[4],m4[5],m4[6],m4[7],m4[8],m4[9],m5[1],m5[2],m5[3],m5[4],m5[5],m5[6],m5[7],m5[8],m5[9]))
      B171<-solve(crossprod(hh17,hh17))%*%crossprod(hh17,mm)
      g_aa1<-(0.5*(B171[7]+B171[10])^2+0.25*(B171[9]+B171[11])^2)/m_nf            
      #   0.5(db+i)**2+0.25(hb+jab)**2.
      g_aa2<-(0.5*(B171[6]+B171[10])^2+0.25*(B171[8]+B171[12])^2)/m_nf        
      #   0.5(da+i)**2+0.25(ha+jba)**2.
      g_aa3<-(0.5*(B171[6]-B171[10])^2+0.25*(B171[8]-B171[12])^2)/m_nf          
      #   0.5(da-i)**2+0.25(ha-jba)**2.
      g_aa4<-(0.5*(B171[7]-B171[10])^2+0.25*(B171[9]-B171[11])^2)/m_nf           
      #   0.5(db-i)**2+0.25(hb-jab)**2.
      g_aa5<-0.25*(B171[6]^2+B171[7]^2+B171[10]^2+(B171[6]+B171[11])^2+(B171[7]+B171[12])^2+(B171[8]+B171[13]/2)^2+(B171[9]+B171[13]/2)^2+B171[13]^2/4)/m_nf
      sigma4[c(2:9)]<-sigma4[1];sigma5[c(3,7,9)]<-sigma5[1]
      sigma5[2]<-sigma5[1]+g_aa1;sigma5[4]<-sigma5[1]+g_aa2
      sigma5[5]<-sigma5[1]+g_aa5;sigma5[6]<-sigma5[1]+g_aa3
      sigma5[8]<-sigma5[1]+g_aa4
      #################restrictions########################################
      hh[1,1]<-sigma4[1]/n0[1]+sigma4[3]/n0[3]+sigma5[1]/n0[10]+sigma5[3]/n0[12]
      hh[1,2]<-sigma4[1]/n0[1]+sigma5[1]/n0[10]
      hh[1,3]<-sigma4[1]/n0[1]+sigma5[1]/n0[10]
      hh[1,4]<-sigma4[1]/n0[1]-sigma4[3]/n0[3]+2*sigma5[1]/n0[10]-2*sigma5[3]/n0[12]
      hh[1,5]<-0;
      hh[1,6]<-sigma4[1]/n0[1]-sigma4[3]/n0[3]
      hh[1,7]<-sigma4[1]/n0[1]+sigma4[3]/n0[3]
      hh[1,8]<--sigma4[1]/n0[1]+sigma4[3]/n0[3]
      hh[2,2]<-sigma4[1]/n0[1]+sigma4[7]/n0[7]+sigma5[1]/n0[10]+sigma5[7]/n0[16]
      hh[2,3]<-sigma4[1]/n0[1]+sigma5[1]/n0[10]
      hh[2,4]<-sigma4[1]/n0[1]-sigma4[7]/n0[7]+2*sigma5[1]/n0[10]-2*sigma5[7]/n0[16]
      hh[2,5]<-0
      hh[2,6]<-sigma4[1]/n0[1]+sigma4[7]/n0[7]
      hh[2,7]<-sigma4[1]/n0[1]-sigma4[7]/n0[7]
      hh[2,8]<--sigma4[1]/n0[1]+sigma4[7]/n0[7]
      hh[3,3]<-sigma4[1]/n0[1]+sigma4[9]/n0[9]+sigma5[1]/n0[10]+sigma5[9]/n0[18]
      hh[3,4]<-sigma4[1]/n0[1]-sigma4[9]/n0[9]+2*sigma5[1]/n0[10]-2*sigma5[9]/n0[18]
      hh[3,5]<-0
      hh[3,6]<-sigma4[1]/n0[1]+sigma4[9]/n0[9]
      hh[3,7]<-sigma4[1]/n0[1]+sigma4[9]/n0[9]
      hh[3,8]<--sigma4[1]/n0[1]+sigma4[9]/n0[9]
      hh[4,4]<-sigma4[1]/n0[1]+4*sigma4[2]/n0[2]+sigma4[3]/n0[3]+sigma4[7]/n0[7]+4*sigma4[8]/n0[8]+sigma4[9]/n0[9]+4*sigma5[1]/n0[10]+16*sigma5[2]/n0[11]+4*sigma5[3]/n0[12]+4*sigma5[7]/n0[16]+16*sigma5[8]/n0[17]+4*sigma5[9]/n0[18]
      hh[4,5]<--2*sigma4[2]/n0[2]-2*sigma4[8]/n0[8]-8*sigma5[2]/n0[11]-8*sigma5[8]/n0[17]
      hh[4,6]<-sigma4[1]/n0[1]-4*sigma4[2]/n0[2]+sigma4[3]/n0[3]-sigma4[7]/n0[7]+4*sigma4[8]/n0[8]-sigma4[9]/n0[9]-16*sigma5[2]/n0[11]+16*sigma5[8]/n0[17]
      hh[4,7]<-sigma4[1]/n0[1]-sigma4[3]/n0[3]+sigma4[7]/n0[7]-sigma4[9]/n0[9]
      hh[4,8]<--sigma4[1]/n0[1]-sigma4[3]/n0[3]-sigma4[7]/n0[7]-sigma4[9]/n0[9]+16*sigma5[2]/n0[11]+16*sigma5[8]/n0[17]
      hh[5,5]<-sigma4[2]/n0[2]+sigma4[4]/n0[4]+sigma4[6]/n0[6]+sigma4[8]/n0[8]+4*sigma5[2]/n0[11]+4*sigma5[4]/n0[13]+4*sigma5[6]/n0[15]+4*sigma5[8]/n0[17]
      hh[5,6]<-2*sigma4[2]/n0[2]-2*sigma4[8]/n0[8]+8*sigma5[2]/n0[11]-8*sigma5[8]/n0[17]
      hh[5,7]<--2*sigma4[4]/n0[4]+2*sigma4[6]/n0[6]-8*sigma5[4]/n0[13]+8*sigma5[6]/n0[15]
      hh[5,8]<--8*sigma5[2]/n0[11]+8*sigma5[4]/n0[13]+8*sigma5[6]/n0[15]-8*sigma5[8]/n0[17]
      hh[6,6]<-sigma4[1]/n0[1]+4*sigma4[2]/n0[2]+sigma4[3]/n0[3]+sigma4[7]/n0[7]+4*sigma4[8]/n0[8]+sigma4[9]/n0[9]+16*sigma5[2]/n0[11]+16*sigma5[8]/n0[17]
      hh[6,7]<-sigma4[1]/n0[1]-sigma4[3]/n0[3]-sigma4[7]/n0[7]+sigma4[9]/n0[9]
      hh[6,8]<--sigma4[1]/n0[1]-sigma4[3]/n0[3]+sigma4[7]/n0[7]+sigma4[9]/n0[9]-16*sigma5[2]/n0[11]+16*sigma5[8]/n0[17]
      hh[7,7]<-sigma4[1]/n0[1]+sigma4[3]/n0[3]+4*sigma4[4]/n0[4]+4*sigma4[6]/n0[6]+sigma4[7]/n0[7]+sigma4[9]/n0[9]+16*sigma5[4]/n0[13]+16*sigma5[6]/n0[15]
      hh[7,8]<--sigma4[1]/n0[1]+sigma4[3]/n0[3]-sigma4[7]/n0[7]+sigma4[9]/n0[9]-16*sigma5[4]/n0[13]+16*sigma5[6]/n0[15]
      hh[8,8]<-sigma4[1]/n0[1]+sigma4[3]/n0[3]+16*sigma4[5]/n0[5]+sigma4[7]/n0[7]+sigma4[9]/n0[9]+16*sigma5[2]/n0[11]+16*sigma5[4]/n0[13]+16*sigma5[6]/n0[15]+16*sigma5[8]/n0[17]+256*sigma5[5]/n0[14]
      for(i in 2:8)
      {
        for(j in 1:(i-1))
        {
          hh[i,j]<-hh[j,i]
        }
      }
      ########################################################
      b_line[1]<-sumwx4[1]/n0[1]-sumwx4[3]/n0[3]-sumwx5[1]/n0[10]+sumwx5[3]/n0[12]
      b_line[2]<-sumwx4[1]/n0[1]-sumwx4[7]/n0[7]-sumwx5[1]/n0[10]+sumwx5[7]/n0[16]
      b_line[3]<-sumwx4[1]/n0[1]-sumwx4[9]/n0[9]-sumwx5[1]/n0[10]+sumwx5[9]/n0[18]
      b_line[4]<-sumwx4[1]/n0[1]-2*sumwx4[2]/n0[2]+sumwx4[3]/n0[3]+sumwx4[7]/n0[7]-2*sumwx4[8]/n0[8]+sumwx4[9]/n0[9]-2*sumwx5[2]/n0[11]+4*sumwx5[2]/n0[11]-2*sumwx5[3]/n0[12]-2*sumwx5[7]/n0[16]+4*sumwx5[8]/n0[17]-2*sumwx5[9]/n0[18]
      b_line[5]<-sumwx4[2]/n0[2]-sumwx4[4]/n0[4]-sumwx4[6]/n0[6]+sumwx4[8]/n0[8]-2*sumwx5[2]/n0[11]+2*sumwx5[4]/n0[13]+2*sumwx5[6]/n0[15]-2*sumwx5[8]/n0[17]
      b_line[6]<-sumwx4[1]/n0[1]+2*sumwx4[2]/n0[2]+sumwx4[3]/n0[3]-sumwx4[7]/n0[7]-2*sumwx4[8]/n0[8]-sumwx4[9]/n0[9]-4*sumwx5[2]/n0[11]+4*sumwx5[8]/n0[17]
      b_line[7]<-sumwx4[1]/n0[1]-sumwx4[3]/n0[3]+2*sumwx4[4]/n0[4]-2*sumwx4[6]/n0[6]+sumwx4[7]/n0[7]-sumwx4[9]/n0[9]-4*sumwx5[4]/n0[13]+4*sumwx5[6]/n0[15]
      b_line[8]<-4*sumwx4[5]/n0[5]-sumwx4[1]/n0[1]-sumwx4[3]/n0[3]-sumwx4[7]/n0[7]-sumwx4[9]/n0[9]+4*sumwx5[2]/n0[11]+4*sumwx5[4]/n0[13]+4*sumwx5[6]/n0[15]+4*sumwx5[8]/n0[17]-16*sumwx5[5]/n0[14]
      B1711<-solve(hh,b_line)
      ##########################################################
      m4[1]<-(sumwx4[1]-sigma4[1]*(B1711[1]+B1711[2]+B1711[3]+B1711[4]+B1711[6]+B1711[7]-B1711[8]))/n0[1]
      m4[2]<-(sumwx4[2]+sigma4[1]*(2*B1711[4]-B1711[5]-2*B1711[6]))/n0[2]
      m4[3]<-(sumwx4[3]+sigma4[1]*(B1711[1]-B1711[4]-B1711[6]+B1711[7]+B1711[8]))/n0[3]
      m4[4]<-(sumwx4[4]+sigma4[1]*(B1711[5]-2*B1711[7]))/n0[4]
      m4[5]<-(sumwx4[5]-sigma4[1]*4*B1711[8])/n0[5]
      m4[6]<-(sumwx4[6]+sigma4[1]*(B1711[5]+2*B1711[7]))/n0[6]
      m4[7]<-(sumwx4[7]+sigma4[1]*(B1711[2]-B1711[4]+B1711[6]-B1711[7]+B1711[8]))/n0[7]
      m4[8]<-(sumwx4[8]+sigma4[1]*(2*B1711[4]-B1711[5]+2*B1711[6]))/n0[8]
      m4[9]<-(sumwx4[9]+sigma4[1]*(B1711[3]-B1711[4]+B1711[6]+B1711[7]+B1711[8]))/n0[9]
      m5[1]<-(sumwx5[1]+sigma5[1]*(B1711[1]+B1711[2]+B1711[3]+2*B1711[4]))/n0[10]
      m5[2]<-(sumwx5[2]+sigma5[2]*(-4*B1711[4]+2*B1711[5]+4*B1711[6]-4*B1711[8]))/n0[11]
      m5[3]<-(sumwx5[3]+sigma5[3]*(-B1711[1]+2*B1711[4]))/n0[12]
      m5[4]<-(sumwx5[4]+sigma5[4]*(-2*B1711[5]+4*B1711[7]-4*B1711[8]))/n0[13]
      m5[5]<-(sumwx5[5]+sigma5[5]*16*B1711[8])/n0[14]
      m5[6]<-(sumwx5[6]+sigma5[6]*(-2*B1711[5]-4*B1711[7]-4*B1711[8]))/n0[15]
      m5[7]<-(sumwx5[7]+sigma5[7]*(-B1711[2]+2*B1711[4]))/n0[16]
      m5[8]<-(sumwx5[8]+sigma5[8]*(-4*B1711[4]+2*B1711[5]-4*B1711[6]-4*B1711[8]))/n0[17]
      m5[9]<-(sumwx5[9]+sigma5[9]*(-B1711[3]+2*B1711[4]))/n0[18]
      aaa1<-max(abs(B1711-AA))
      AA<-B1711
      if (n_iter>20) break
    } 
    ##########obtain variance###############
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d3) {swx24[i] <- W4[i,]%*%(dataF2-m4[i])^2 } ;for(i in 1:d3) {swx25[i] <- W5[i,]%*%(dataF3-m5[i])^2 } 
    ab1<-sum(swx24);sigma4[1]<-ab1/nn[4];sigma40<-sigma4[1]-sigma
    if (sigma40<0) {sigma40<-0;sigma4[1]<-sigma}
    sigma4[c(2:9)]<-sigma4[1];aaa0<-sigma5[1]
    ab2<-swx25[1]+swx25[3]+swx25[7]+swx25[9];ab3<-n0[10]+n0[12]+n0[16]+n0[18]
    mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m4[4],m4[5],m4[6],m4[7],m4[8],m4[9],m5[1],m5[2],m5[3],m5[4],m5[5],m5[6],m5[7],m5[8],m5[9]))
    B17111<-solve(crossprod(hh17,hh17))%*%crossprod(hh17,mm)
    g_aa1<-(0.5*(B17111[7]+B17111[10])^2+0.25*(B17111[9]+B17111[11])^2)/m_nf            
    #   0.5(db+i)**2+0.25(hb+jab)**2.
    g_aa2<-(0.5*(B17111[6]+B17111[10])^2+0.25*(B17111[8]+B17111[12])^2)/m_nf        
    #   0.5(da+i)**2+0.25(ha+jba)**2.
    g_aa3<-(0.5*(B17111[6]-B17111[10])^2+0.25*(B17111[8]-B17111[12])^2)/m_nf          
    #   0.5(da-i)**2+0.25(ha-jba)**2.
    g_aa4<-(0.5*(B17111[7]-B17111[10])^2+0.25*(B17111[9]-B17111[11])^2)/m_nf           
    #   0.5(db-i)**2+0.25(hb-jab)**2.
    g_aa5<-0.25*(B17111[6]^2+B17111[7]^2+B17111[10]^2+(B17111[6]+B17111[11])^2+(B17111[7]+B17111[12])^2+(B17111[8]+B17111[13]/2)^2+(B17111[9]+B17111[13]/2)^2+B17111[13]^2/4)/m_nf
    n_iter<-0;aaa1<-1000
    while (aaa1>0.0001){
      n_iter<-n_iter+1
      aa1<-sigma5[1]/(sigma5[1]+g_aa1);aa2<-sigma5[1]/(sigma5[1]+g_aa2)
      aa3<-sigma5[1]/(sigma5[1]+g_aa3);aa4<-sigma5[1]/(sigma5[1]+g_aa4)
      aa5<-sigma5[1]/(sigma5[1]+g_aa5)
      as3<-ab2+aa1^2*swx25[2]+aa2^2*swx25[4]+aa3^2*swx25[6]+aa4^2*swx25[8]+aa5^2*swx25[5]
      as4<-ab3+aa1*n0[11]+aa2*n0[13]+aa3*n0[15]+aa4*n0[17]+aa5*n0[14]
      sigma5[1]<-as3/as4
      aaa1<-abs(sigma5[1]-aaa0)
      aaa0<-sigma5[1]
      if (n_iter>20) break 
    }
    sigma50<-sigma5[1]-sigma/m_nf
    if (sigma50<0) {sigma50<-0;sigma5[1]<-sigma/m_nf}
    sigma5[c(3,7,9)]<-sigma5[1];sigma5[2]<-sigma5[1]+g_aa1
    sigma5[4]<-sigma5[1]+g_aa2;sigma5[5]<-sigma5[1]+g_aa5
    sigma5[6]<-sigma5[1]+g_aa3;sigma5[8]<-sigma5[1]+g_aa4
    ab2<-ss1+ss2+ss3;ab3<-nn[1]+nn[2]+nn[3]
    n_iter<-0;aaa0<-sigma;aaa1<-1000
    while (aaa1>0.0001){
      n_iter<-n_iter+1  
      ab4<-sigma/(sigma40+sigma)
      s0[1]<-(sigma/m_nf)/(sigma/m_nf+sigma50);s0[2]<-(sigma/m_nf)/(sigma/m_nf+sigma50+g_aa1)
      s0[3]<-(sigma/m_nf)/(sigma/m_nf+sigma50);s0[4]<-(sigma/m_nf)/(sigma/m_nf+sigma50+g_aa2)
      s0[5]<-(sigma/m_nf)/(sigma/m_nf+sigma50+g_aa5);s0[6]<-(sigma/m_nf)/(sigma/m_nf+sigma50+g_aa3)
      s0[7]<-(sigma/m_nf)/(sigma/m_nf+sigma50);s0[8]<-(sigma/m_nf)/(sigma/m_nf+sigma50+g_aa4)
      s0[9]<-(sigma/m_nf)/(sigma/m_nf+sigma50)
      as1<-sum(s0[c(1:9)]^2*swx25[c(1:9)])
      as2<-sum(s0[c(1:9)]*n0[c(10:18)])
      as3<-sum(swx24[c(1:9)])
      sigma<-(ab2+ab4^2*as3+as1*m_nf)/(ab3+ab4*nn[4]+as2)
      aaa1<-abs(sigma-aaa0)
      aaa0<-sigma
      if (n_iter>20) break
    }
    if (sigma<sigma0) {sigma<-sigma0}
    sigma4[1]<-sigma+sigma40;sigma4[c(2:9)]<-sigma4[1]
    sigma5[1]<-sigma/m_nf+sigma50;sigma5[c(3,7,9)]<-sigma5[1]
    sigma5[2]<-sigma5[1]+g_aa1;sigma5[4]<-sigma5[1]+g_aa2
    sigma5[5]<-sigma5[1]+g_aa5;sigma5[6]<-sigma5[1]+g_aa3
    sigma5[8]<-sigma5[1]+g_aa4
    ##############################criteria for iterations to stop######################## 
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
      sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*16
  #######first order parameters############
  mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m4[4],m4[5],m4[6],m4[7],m4[8],m4[9],m5[1],m5[2],m5[3],m5[4],m5[5],m5[6],m5[7],m5[8],m5[9]))
  B171111<-solve(crossprod(hh17,hh17))%*%crossprod(hh17,mm)
  ########second order genetic parameters###############   
  F2_jj <- sigmaF2 - sigma4[1]
  F2_gg <- sigma4[1]-sigma
  if(F2_jj<0) {F2_jj<-0}
  if(F2_gg<0 || F2_gg>sigmaF2) {F2_gg<-0}
  F2_ll <- F2_jj/sigmaF2
  F2_rr <- F2_gg/sigmaF2
  F3_jj <- sigmaF3 - sigma5[1]
  F3_gg <- sigma5[1]-sigma/m_nf
  if(F3_jj<0) {F3_jj<-0}
  if(F3_gg<0 || F3_gg>sigmaF3) {F3_gg<-0}
  F3_ll <- F3_jj/sigmaF3
  F3_rr <- F3_gg/sigmaF3
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma))
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
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma))
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
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma))
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
    F2gg <- (dataF2 - m4[i])/sqrt(sigma4[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mi4[i]
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
  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3); 
  F3w1<-1/(12*m_sam[5])
  F3bmw <- matrix(0,m_sam[5],1); F3bmwsl <- matrix(0,m_sam[5],d3)
  for(i in 1:d3){
    F3gg <- (dataF3 - m5[i])/sqrt(sigma5[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mi5[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[5]){F3P2 <- F3P2+runif(m_sam[5])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[5]) + sum((F3P2 - (as.matrix(c(1:m_sam[5])) - 0.5)/m_sam[5])^2)
  F3u <- as.matrix(c(12*m_sam[5]*((F3dd[1]/m_sam[5]-0.5)^2),((45*m_sam[5])/4)*((F3dd[2]/m_sam[5]-1/3)^2),180*m_sam[5]*((F3dd[3]/m_sam[5]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[5]))))
  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("MX2-ADI-ADI",round(abc,4),round(AIC,4),round(m[1],4),round(m[2],4),round(m[3],4), round(t(m4),4),round(sigma4[1],4),
                       round(t(mi4),4),round(t(m5),4),round(t(sigma5),4),round(t(mi5),4),round(sigma,4),          
                       round(B171111[1],4),round(B171111[2],4),round(B171111[3],4),round(B171111[4],4),round(B171111[5],4),round(B171111[6],4),round(B171111[7],4),round(B171111[8],4),round(B171111[9],4),round(B171111[10],4),round(B171111[11],4),round(B171111[12],4),round(B171111[13],4)," "," ",round(F2_jj,4),round(F2_ll*100,4),round(F2_gg,4),round(F2_rr*100,4),round(F3_jj,4),round(F3_ll*100,4),round(F3_gg,4),round(F3_rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
 output<-as.matrix(output)
  OUTPUT<-list(output,mix_pi4,mix_pi5)
  return(OUTPUT)
}
##################MX2-ADI-AD(E-1)#############################
G5ModelFun[[19]] <- function(K1,logL,df11,df21,df31,df41,df51,G5text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataF2 <- as.matrix(as.numeric(df41[,1]))
  dataF3 <- as.matrix(as.numeric(df51[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2),sum(dataF3)))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2),mean(dataF3)))
  sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1);sigmaF2<-var(dataF2);sigmaF3<-var(dataF3)
  sumx6<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  sigma<-sumx6;sigma0<-sigma;sigma40<-sigmaF2;sigma400<-sigma40;sigma50<-sigmaF3;sigma500<-sigma50
  m_nf <- as.numeric(G5text2)
  ###############procedure start###########################
  mix_pi4<-as.matrix(c(0.0625,0.125,0.0625,0.1250,0.250,0.125,0.0625,0.125,0.0625))
  mix_pi5<-as.matrix(c(0.0625,0.125,0.0625,0.1250,0.250,0.125,0.0625,0.125,0.0625))
  a1<-sqrt(sigma400/m_sam[4])
  if (m[1]<m[3]) a1<--a1
  m4<-as.matrix(c(m[4]+2.8*a1,m[4]+2.1*a1,m[4]+1.4*a1,m[4]+0.7*a1,m[4],m[4]-0.7*a1,m[4]-1.4*a1,m[4]-2.1*a1,m[4]-2.8*a1))
  a1<-sqrt(sigma500/m_sam[5])
  if (m[1]<m[3]) {a1<--a1}
  m5<-as.matrix(c(m[5]+2.8*a1,m[5]+2.1*a1,m[5]+1.4*a1,m[5]+0.7*a1,m[5],m[5]-0.7*a1,m[5]-1.4*a1,m[5]-2.1*a1,m[5]-2.8*a1))
  nn<-m_sam
  abb<-sigma400/(sigma*5)
  sigma4<-matrix(0,9,1);sigma5<-matrix(0,9,1)
  sigma4[1]<-sigma400/abb;sigma5[1]<-sigma500/abb
  #######first order parameters############
  hh18<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,-1,1,1,1,0,0,0,-1,-1,-1,1,1,1,0,0,0,-1,-1,-1,
                 1,0,-1,1,0,-1,1,0,-1,1,0,-1,1,0,-1,1,0,-1,1,0,-1,0,1,0,0,0,0,1,1,1,0,0,0,0,0,0,0.5,0.5,0.5,0,0,0,
                 0,1,0,0,1,0,0,1,0,0,1,0,0,0.5,0,0,0.5,0,0,0.5,0,1,0,1,1,0,-1,0,0,0,-1,0,1,1,0,-1,0,0,0,-1,0,1,
                 0,0,0,0,1,0,0,0,0,0,-1,0,0,0.5,0,0,0,0,0,-0.5,0,0,0,0,0,0,0,1,0,-1,0,0,0,0,0,0,0.5,0,-0.5,0,0,0,
                 0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0.25,0,0,0,0,1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                 0,1,0,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25),21,11)
  mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m4[4],m4[5],m4[6],m4[7],m4[8],m4[9],m5[1],m5[2],m5[3],m5[4],m5[5],m5[6],m5[7],m5[8],m5[9]))
  B18<-solve(crossprod(hh18,hh18))%*%crossprod(hh18,mm)
  g_aa1<-(0.5*(B18[3]+B18[6])^2+0.25*(B18[5]+B18[7])^2)/m_nf         
  #   0.5(db+i)**2+0.25(hb+jab)**2.
  g_aa2<-(0.5*(B18[2]+B18[6])^2+0.25*(B18[4]+B18[8])^2)/m_nf           
  #   0.5(da+i)**2+0.25(ha+jba)**2.
  g_aa3<-(0.5*(B18[2]-B18[6])^2+0.25*(B18[4]-B18[8])^2)/m_nf     
  #   0.5(da-i)**2+0.25(ha-jba)**2.
  g_aa4<-(0.5*(B18[3]-B18[6])^2+0.25*(B18[5]-B18[7])^2)/m_nf        
  #   0.5(db-i)**2+0.25(hb-jab)**2.
  g_aa5<-0.25*(B18[2]^2+B18[3]^2+B18[6]^2+(B18[2]+B18[7])^2+(B18[3]+B18[8])^2+(B18[4]+B18[9]/2)^2+(B18[5]+B18[9]/2)^2+B18[9]^2/4)/m_nf
  sigma4[c(2:9)]<-sigma4[1];sigma5[c(3,7,9)]<-sigma5[1]
  sigma5[2]<-sigma5[1]+g_aa1;sigma5[4]<-sigma5[1]+g_aa2
  sigma5[5]<-sigma5[1]+g_aa5;sigma5[6]<-sigma5[1]+g_aa3
  sigma5[8]<-sigma5[1]+g_aa4
  mi4<-mix_pi4[c(1:9)];mi5<-mix_pi5[c(1:9)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
    sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
  ########################################################
  ##########iteration process###########
  ############E-step###################
  ############F2 process###############
  d3<-9
  iteration <- 0; stopa <- 1000
  W4 <- matrix(0,9,m_sam[4]); swx24 <- matrix(0,9,1)
  W5 <- matrix(0,9,m_sam[5]); swx25 <- matrix(0,9,1)
  hh<-matrix(0,10,10);b_line<-matrix(0,10,1)
  n0<-matrix(0,18,1);s0<-matrix(0,9,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:d3) { W4[i,] <- mix_pi4[i]*dnorm(dataF2,m4[i],sqrt(sigma4[i]))/dmixnorm(dataF2,m4,sqrt(sigma4),mix_pi4) }
    mi4 <- as.matrix(rowSums(W4)/m_sam[4]);n0[c(1:9)]<-as.matrix(rowSums(W4));sumwx4 <- W4%*%dataF2
    #############F3 process############
    for(i in 1:d3) { W5[i,] <- mix_pi5[i]*dnorm(dataF3,m5[i],sqrt(sigma5[i]))/dmixnorm(dataF3,m5,sqrt(sigma5),mix_pi5) }
    mi5 <- as.matrix(rowSums(W5)/m_sam[5]);n0[c(10:18)]<-as.matrix(rowSums(W5))
    n0[c(1:18)][abs(n0[c(1:18)])<0.000001]<-0.000001;sumwx5 <- W5%*%dataF3
    aaa0<-0
    aaa1<-1000; n_iter<-0;AA<-matrix(0,10,1)
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      ##########first order parameters###################
      mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m4[4],m4[5],m4[6],m4[7],m4[8],m4[9],
                      m5[1],m5[2],m5[3],m5[4],m5[5],m5[6],m5[7],m5[8],m5[9]))
      B181<-solve(crossprod(hh18,hh18))%*%crossprod(hh18,mm)
      g_aa1<-(0.5*(B181[3]+B181[6])^2+0.25*(B181[5]+B181[7])^2)/m_nf         
      #   0.5(db+i)**2+0.25(hb+jab)**2.
      g_aa2<-(0.5*(B181[2]+B181[6])^2+0.25*(B181[4]+B181[8])^2)/m_nf           
      #   0.5(da+i)**2+0.25(ha+jba)**2.
      g_aa3<-(0.5*(B181[2]-B181[6])^2+0.25*(B181[4]-B181[8])^2)/m_nf     
      #   0.5(da-i)**2+0.25(ha-jba)**2.
      g_aa4<-(0.5*(B181[3]-B181[6])^2+0.25*(B181[5]-B181[7])^2)/m_nf        
      #   0.5(db-i)**2+0.25(hb-jab)**2.
      g_aa5<-0.25*(B181[2]^2+B181[3]^2+B181[6]^2+(B181[2]+B181[7])^2+(B181[3]+B181[8])^2+(B181[4]+B181[9]/2)^2+(B181[5]+B181[9]/2)^2+B181[9]^2/4)/m_nf
      sigma4[c(2:9)]<-sigma4[1];sigma5[c(3,7,9)]<-sigma5[1]
      sigma5[2]<-sigma5[1]+g_aa1;sigma5[4]<-sigma5[1]+g_aa2
      sigma5[5]<-sigma5[1]+g_aa5;sigma5[6]<-sigma5[1]+g_aa3
      sigma5[8]<-sigma5[1]+g_aa4
      #################restrictions########################################
      hh[1,1]<-sigma*(1/nn[1]+4/nn[2]+1/nn[3])+sigma4[1]/n0[1]+4*sigma4[5]/n0[5]+sigma4[9]/n0[9]
      hh[1,2]<-sigma*2/nn[2]+2*sigma4[5]/n0[5]+2*sigma4[9]/n0[9]
      hh[1,3]<--sigma4[1]/n0[1]
      hh[1,4]<--sigma4[1]/n0[1]
      hh[1,5]<--sigma4[1]/n0[1]+sigma4[9]/n0[9]
      hh[1,6]<--sigma4[1]/n0[1]-sigma4[9]/n0[9]
      hh[1,7]<-0
      hh[1,8]<--sigma4[1]/n0[1]+sigma4[9]/n0[9]
      hh[1,9]<--sigma4[1]/n0[1]+sigma4[9]/n0[9]
      hh[1,10]<-sigma4[1]/n0[1]-8*sigma4[5]/n0[5]+sigma4[9]/n0[9]
      hh[2,2]<-sigma/nn[2]+sigma4[5]/n0[5]+4*sigma4[9]/n0[9]+4*sigma5[9]/n0[18]
      hh[2,3]<-0
      hh[2,4]<-0
      hh[2,5]<-2*sigma4[9]/n0[9]+2*sigma5[9]/n0[18]
      hh[2,6]<--2*sigma4[9]/n0[9]-4*sigma5[9]/n0[18]
      hh[2,7]<-0
      hh[2,8]<-2*sigma4[9]/n0[9]
      hh[2,9]<-2*sigma4[9]/n0[9]
      hh[2,10]<--4*sigma4[5]/n0[5]+2*sigma4[9]/n0[9]
      hh[3,3]<-sigma4[1]/n0[1]+sigma4[3]/n0[3]+sigma5[1]/n0[10]+sigma5[3]/n0[12]
      hh[3,4]<-sigma4[1]/n0[1]+sigma5[1]/n0[10]
      hh[3,5]<-sigma4[1]/n0[1]+sigma5[1]/n0[10]
      hh[3,6]<-sigma4[1]/n0[1]-sigma4[3]/n0[3]+2*sigma5[1]/n0[10]-2*sigma5[3]/n0[12]
      hh[3,7]<-0
      hh[3,8]<-sigma4[1]/n0[1]-sigma4[3]/n0[3]
      hh[3,9]<-sigma4[1]/n0[1]+sigma4[3]/n0[3]
      hh[3,10]<--sigma4[1]/n0[1]+sigma4[3]/n0[3]
      hh[4,4]<-sigma4[1]/n0[1]+sigma4[7]/n0[7]+sigma5[1]/n0[10]+sigma5[7]/n0[16]
      hh[4,5]<-sigma4[1]/n0[1]+sigma5[1]/n0[10]
      hh[4,6]<-sigma4[1]/n0[1]-sigma4[7]/n0[7]+2*sigma5[1]/n0[10]-2*sigma5[7]/n0[16]
      hh[4,7]<-0
      hh[4,8]<-sigma4[1]/n0[1]+sigma4[7]/n0[7]
      hh[4,9]<-sigma4[1]/n0[1]-sigma4[7]/n0[7]
      hh[4,10]<--sigma4[1]/n0[1]+sigma4[7]/n0[7]
      hh[5,5]<-sigma4[1]/n0[1]+sigma4[9]/n0[9]+sigma5[1]/n0[10]+sigma5[9]/n0[18]
      hh[5,6]<-sigma4[1]/n0[1]-sigma4[9]/n0[9]+2*sigma5[1]/n0[10]-2*sigma5[9]/n0[18]
      hh[5,7]<-0
      hh[5,8]<-sigma4[1]/n0[1]+sigma4[9]/n0[9]
      hh[5,9]<-sigma4[1]/n0[1]+sigma4[9]/n0[9]
      hh[5,10]<--sigma4[1]/n0[1]+sigma4[9]/n0[9]
      hh[6,6]<-sigma4[1]/n0[1]+4*sigma4[2]/n0[2]+sigma4[3]/n0[3]+sigma4[7]/n0[7]+4*sigma4[8]/n0[8]+sigma4[9]/n0[9]+4*sigma5[1]/n0[10]+16*sigma5[2]/n0[11]+4*sigma5[3]/n0[12]+4*sigma5[7]/n0[16]+16*sigma5[8]/n0[17]+4*sigma5[9]/n0[18]
      hh[6,7]<--2*sigma4[2]/n0[2]-2*sigma4[8]/n0[8]-8*sigma5[2]/n0[11]-8*sigma5[8]/n0[17]
      hh[6,8]<-sigma4[1]/n0[1]-4*sigma4[2]/n0[2]+sigma4[3]/n0[3]-sigma4[7]/n0[7]+4*sigma4[8]/n0[8]-sigma4[9]/n0[9]-16*sigma5[2]/n0[11]+16*sigma5[8]/n0[17]
      hh[6,9]<-sigma4[1]/n0[1]-sigma4[3]/n0[3]+sigma4[7]/n0[7]-sigma4[9]/n0[9]
      hh[6,10]<--sigma4[1]/n0[1]-sigma4[3]/n0[3]-sigma4[7]/n0[7]-sigma4[9]/n0[9]+16*sigma5[2]/n0[11]+16*sigma5[8]/n0[17]
      hh[7,7]<-sigma4[2]/n0[2]+sigma4[4]/n0[4]+sigma4[6]/n0[6]+sigma4[8]/n0[8]+4*sigma5[2]/n0[11]+4*sigma5[4]/n0[13]+4*sigma5[6]/n0[15]+4*sigma5[8]/n0[17]
      hh[7,8]<-2*sigma4[2]/n0[2]-2*sigma4[8]/n0[8]+8*sigma5[2]/n0[11]-8*sigma5[8]/n0[17] 
      hh[7,9]<--2*sigma4[4]/n0[4]+2*sigma4[6]/n0[6]-8*sigma5[4]/n0[13]+8*sigma5[6]/n0[15]
      hh[7,10]<--8*sigma5[2]/n0[11]+8*sigma5[4]/n0[13]+8*sigma5[6]/n0[15]-8*sigma5[8]/n0[17]
      hh[8,8]<-sigma4[1]/n0[1]+4*sigma4[2]/n0[2]+sigma4[3]/n0[3]+sigma4[7]/n0[7]+4*sigma4[8]/n0[8]+sigma4[9]/n0[9]+16*sigma5[2]/n0[11]+16*sigma5[8]/n0[17]
      hh[8,9]<-sigma4[1]/n0[1]-sigma4[3]/n0[3]-sigma4[7]/n0[7]+sigma4[9]/n0[9]
      hh[8,10]<--sigma4[1]/n0[1]-sigma4[3]/n0[3]+sigma4[7]/n0[7]+sigma4[9]/n0[9]-16*sigma5[2]/n0[11]+16*sigma5[8]/n0[17]
      hh[9,9]<-sigma4[1]/n0[1]+sigma4[3]/n0[3]+4*sigma4[4]/n0[4]+4*sigma4[6]/n0[6]+sigma4[7]/n0[7]+sigma4[9]/n0[9]+16*sigma5[4]/n0[13]+16*sigma5[6]/n0[15]
      hh[9,10]<--sigma4[1]/n0[1]+sigma4[3]/n0[3]-sigma4[7]/n0[7]+sigma4[9]/n0[9]-16*sigma5[4]/n0[13]+16*sigma5[6]/n0[15]
      hh[10,10]<-sigma4[1]/n0[1]+sigma4[3]/n0[3]+16*sigma4[5]/n0[5]+sigma4[7]/n0[7]+sigma4[9]/n0[9]+16*sigma5[2]/n0[11]+16*sigma5[4]/n0[13]+16*sigma5[6]/n0[15]+16*sigma5[8]/n0[17]+256*sigma5[5]/n0[14]
      for(i in 2:10)
      {
        for(j in 1:(i-1))
        {
          hh[i,j]<-hh[j,i]
        }
      }
      ########################################################
      b_line[1]<-sumx[1]/nn[1]+2*sumx[2]/nn[2]+sumx[3]/nn[3]-sumwx4[1]/n0[1]-2*sumwx4[5]/n0[5]-sumwx4[9]/n0[9]
      b_line[2]<-sumx[2]/nn[2]-sumwx4[5]/n0[5]-2*sumwx4[9]/n0[9]+2*sumwx5[9]/n0[18]
      b_line[3]<-sumwx4[1]/n0[1]-sumwx4[3]/n0[3]-sumwx5[1]/n0[10]+sumwx5[3]/n0[12]
      b_line[4]<-sumwx4[1]/n0[1]-sumwx4[7]/n0[7]-sumwx5[1]/n0[10]+sumwx5[7]/n0[16]
      b_line[5]<-sumwx4[1]/n0[1]-sumwx4[9]/n0[9]-sumwx5[1]/n0[10]+sumwx5[9]/n0[18]
      b_line[6]<-sumwx4[1]/n0[1]-2*sumwx4[2]/n0[2]+sumwx4[3]/n0[3]+sumwx4[7]/n0[7]-2*sumwx4[8]/n0[8]+sumwx4[9]/n0[9]-2*sumwx5[2]/n0[11]+4*sumwx5[2]/n0[11]-2*sumwx5[3]/n0[12]-2*sumwx5[7]/n0[16]+4*sumwx5[8]/n0[17]-2*sumwx5[9]/n0[18]
      b_line[7]<-sumwx4[2]/n0[2]-sumwx4[4]/n0[4]-sumwx4[6]/n0[6]+sumwx4[8]/n0[8]-2*sumwx5[2]/n0[11]+2*sumwx5[4]/n0[13]+2*sumwx5[6]/n0[15]-2*sumwx5[8]/n0[17]
      b_line[8]<-sumwx4[1]/n0[1]+2*sumwx4[2]/n0[2]+sumwx4[3]/n0[3]-sumwx4[7]/n0[7]-2*sumwx4[8]/n0[8]-sumwx4[9]/n0[9]-4*sumwx5[2]/n0[11]+4*sumwx5[8]/n0[17]
      b_line[9]<-sumwx4[1]/n0[1]-sumwx4[3]/n0[3]+2*sumwx4[4]/n0[4]-2*sumwx4[6]/n0[6]+sumwx4[7]/n0[7]-sumwx4[9]/n0[9]-4*sumwx5[4]/n0[13]+4*sumwx5[6]/n0[15]
      b_line[10]<-4*sumwx4[5]/n0[5]-sumwx4[1]/n0[1]-sumwx4[3]/n0[3]-sumwx4[7]/n0[7]-sumwx4[9]/n0[9]+4*sumwx5[2]/n0[11]+4*sumwx5[4]/n0[13]+4*sumwx5[6]/n0[15]+4*sumwx5[8]/n0[17]-16*sumwx5[5]/n0[14]
      B1811<-solve(hh,b_line)
      ##########################################################
      m[1]<-(sumx[1]-sigma*B1811[1])/nn[1]
      m[2]<-(sumx[2]-sigma*(2*B1811[1]+B1811[2]))/nn[2]
      m[3]<-(sumx[3]-sigma*B1811[1])/nn[3]
      m4[1]<-(sumwx4[1]+sigma4[1]*(B1811[1]-B1811[3]-B1811[4]-B1811[5]-B1811[6]-B1811[8]-B1811[9]+B1811[10]))/n0[1]
      m4[2]<-(sumwx4[2]+sigma4[1]*(2*B1811[6]-B1811[7]-2*B1811[8]))/n0[2]
      m4[3]<-(sumwx4[3]+sigma4[1]*(B1811[3]-B1811[6]-B1811[8]+B1811[9]+B1811[10]))/n0[3]
      m4[4]<-(sumwx4[4]+sigma4[1]*(B1811[7]-2*B1811[9]))/n0[4]
      m4[5]<-(sumwx4[5]+sigma4[1]*(2*B1811[1]+B1811[2]-4*B1811[10]))/n0[5]
      m4[6]<-(sumwx4[6]+sigma4[1]*(B1811[7]+2*B1811[9]))/n0[6]
      m4[7]<-(sumwx4[7]+sigma4[1]*(B1811[4]-B1811[6]+B1811[8]-B1811[9]+B1811[10]))/n0[7]
      m4[8]<-(sumwx4[8]+sigma4[1]*(2*B1811[6]-B1811[7]+2*B1811[8]))/n0[8]
      m4[9]<-(sumwx4[9]+sigma4[1]*(B1811[1]+2*B1811[2]+B1811[5]-B1811[6]+B1811[8]+B1811[9]+B1811[10]))/n0[9]
      m5[1]<-(sumwx5[1]+sigma5[1]*(B1811[3]+B1811[4]+B1811[5]+2*B1811[6]))/n0[10]
      m5[2]<-(sumwx5[2]+sigma5[2]*(-4*B1811[6]+2*B1811[7]+4*B1811[8]-4*B1811[10]))/n0[11]
      m5[3]<-(sumwx5[3]+sigma5[3]*(-B1811[3]+2*B1811[6]))/n0[12]
      m5[4]<-(sumwx5[4]+sigma5[4]*(-2*B1811[7]+4*B1811[9]-4*B1811[10]))/n0[13]
      m5[5]<-(sumwx5[5]+sigma5[5]*16*B1811[10])/n0[14]
      m5[6]<-(sumwx5[6]-sigma5[6]*(2*B1811[7]+4*B1811[9]+4*B1811[10]))/n0[15]
      m5[7]<-(sumwx5[7]+sigma5[7]*(-B1811[4]+2*B1811[6]))/n0[16]
      m5[8]<-(sumwx5[8]+sigma5[8]*(-4*B1811[6]+2*B1811[7]-4*B1811[8]-4*B1811[10]))/n0[17]
      m5[9]<-(sumwx5[9]+sigma5[9]*(-2*B1811[2]-B1811[5]+2*B1811[6]))/n0[18]
      aaa1<-max(abs(B1811-AA))
      AA<-B1811
      if (n_iter>20) break
    } 
    ##########obtain variance###############
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d3) {swx24[i] <- W4[i,]%*%(dataF2-m4[i])^2 } ;for(i in 1:d3) {swx25[i] <- W5[i,]%*%(dataF3-m5[i])^2 } 
    ab1<-sum(swx24);sigma4[1]<-ab1/nn[4];sigma40<-sigma4[1]-sigma
    if (sigma40<0) {sigma40<-0;sigma4[1]<-sigma}
    sigma4[c(2:9)]<-sigma4[1];aaa0<-sigma5[1]
    ab2<-swx25[1]+swx25[3]+swx25[7]+swx25[9];ab3<-n0[10]+n0[12]+n0[16]+n0[18]
    n_iter<-0;aaa1<-1000
    mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m4[4],m4[5],m4[6],m4[7],m4[8],m4[9],
                    m5[1],m5[2],m5[3],m5[4],m5[5],m5[6],m5[7],m5[8],m5[9]))
    B18111<-solve(crossprod(hh18,hh18))%*%crossprod(hh18,mm)
    g_aa1<-(0.5*(B18111[3]+B18111[6])^2+0.25*(B18111[5]+B18111[7])^2)/m_nf         
    #   0.5(db+i)**2+0.25(hb+jab)**2.
    g_aa2<-(0.5*(B18111[2]+B18111[6])^2+0.25*(B18111[4]+B18111[8])^2)/m_nf           
    #   0.5(da+i)**2+0.25(ha+jba)**2.
    g_aa3<-(0.5*(B18111[2]-B18111[6])^2+0.25*(B18111[4]-B18111[8])^2)/m_nf     
    #   0.5(da-i)**2+0.25(ha-jba)**2.
    g_aa4<-(0.5*(B18111[3]-B18111[6])^2+0.25*(B18111[5]-B18111[7])^2)/m_nf        
    #   0.5(db-i)**2+0.25(hb-jab)**2.
    g_aa5<-0.25*(B18111[2]^2+B18111[3]^2+B18111[6]^2+(B18111[2]+B18111[7])^2+(B18111[3]+B18111[8])^2+(B18111[4]+B18111[9]/2)^2+(B18111[5]+B18111[9]/2)^2+B18111[9]^2/4)/m_nf
    while (aaa1>0.0001){
      n_iter<-n_iter+1
      aa1<-sigma5[1]/(sigma5[1]+g_aa1);aa2<-sigma5[1]/(sigma5[1]+g_aa2)
      aa3<-sigma5[1]/(sigma5[1]+g_aa3);aa4<-sigma5[1]/(sigma5[1]+g_aa4)
      aa5<-sigma5[1]/(sigma5[1]+g_aa5)
      as3<-ab2+aa1^2*swx25[2]+aa2^2*swx25[4]+aa3^2*swx25[6]+aa4^2*swx25[8]+aa5^2*swx25[5]
      as4<-ab3+aa1*n0[11]+aa2*n0[13]+aa3*n0[15]+aa4*n0[17]+aa5*n0[14]
      sigma5[1]<-as3/as4
      aaa1<-abs(sigma5[1]-aaa0)
      aaa0<-sigma5[1]
      if (n_iter>20) break 
    }
    sigma50<-sigma5[1]-sigma/m_nf
    if (sigma50<0) {sigma50<-0;sigma5[1]<-sigma/m_nf}
    sigma5[c(3,7,9)]<-sigma5[1]
    sigma5[2]<-sigma5[1]+g_aa1;sigma5[4]<-sigma5[1]+g_aa2
    sigma5[5]<-sigma5[1]+g_aa5;sigma5[6]<-sigma5[1]+g_aa3
    sigma5[8]<-sigma5[1]+g_aa4
    ab2<-ss1+ss2+ss3;ab3<-nn[1]+nn[2]+nn[3]
    n_iter<-0;aaa0<-sigma
    aaa1<-1000
    while (aaa1>0.0001){
      n_iter<-n_iter+1  
      ab4<-sigma/(sigma40+sigma)
      s0[1]<-(sigma/m_nf)/(sigma/m_nf+sigma50)
      s0[2]<-(sigma/m_nf)/(sigma/m_nf+sigma50+g_aa1)
      s0[3]<-(sigma/m_nf)/(sigma/m_nf+sigma50)
      s0[4]<-(sigma/m_nf)/(sigma/m_nf+sigma50+g_aa2)
      s0[5]<-(sigma/m_nf)/(sigma/m_nf+sigma50+g_aa5)
      s0[6]<-(sigma/m_nf)/(sigma/m_nf+sigma50+g_aa3)
      s0[7]<-(sigma/m_nf)/(sigma/m_nf+sigma50)
      s0[8]<-(sigma/m_nf)/(sigma/m_nf+sigma50+g_aa4)
      s0[9]<-(sigma/m_nf)/(sigma/m_nf+sigma50)
      as1<-sum(s0[c(1:9)]^2*swx25[c(1:9)])
      as2<-sum(s0[c(1:9)]*n0[c(10:18)])
      as3<-sum(swx24[c(1:9)])
      sigma<-(ab2+ab4^2*as3+as1*m_nf)/(ab3+ab4*nn[4]+as2)
      aaa1<-abs(sigma-aaa0)
      aaa0<-sigma
      if (n_iter>20) break
    }
    
    if (sigma<sigma0) {sigma<-sigma0}
    sigma4[1]<-sigma+sigma40;sigma4[c(2:9)]<-sigma4[1]
    sigma5[1]<-sigma/m_nf+sigma50;sigma5[c(3,7,9)]<-sigma5[1]
    sigma5[2]<-sigma5[1]+g_aa1;sigma5[4]<-sigma5[1]+g_aa2
    sigma5[5]<-sigma5[1]+g_aa5;sigma5[6]<-sigma5[1]+g_aa3
    sigma5[8]<-sigma5[1]+g_aa4
    ##############################criteria for iterations to stop######################## 
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
      sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*14
  #######first order parameters############
  mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m4[4],m4[5],m4[6],m4[7],m4[8],m4[9],m5[1],m5[2],m5[3],m5[4],m5[5],m5[6],m5[7],m5[8],m5[9]))
  B181111<-solve(crossprod(hh18,hh18))%*%crossprod(hh18,mm)
  ########second order genetic parameters###############   
  F2_jj <- sigmaF2 - sigma4[1]
  F2_gg <- sigma4[1]-sigma
  if(F2_jj<0) {F2_jj<-0}
  if(F2_gg<0 || F2_gg>sigmaF2) {F2_gg<-0}
  F2_ll <- F2_jj/sigmaF2
  F2_rr <- F2_gg/sigmaF2
  F3_jj <- sigmaF3 - sigma5[1]
  F3_gg <- sigma5[1]-sigma/m_nf
  if(F3_jj<0) {F3_jj<-0}
  if(F3_gg<0 || F3_gg>sigmaF3) {F3_gg<-0}
  F3_ll <- F3_jj/sigmaF3
  F3_rr <- F3_gg/sigmaF3
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma))
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
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma))
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
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma))
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
    F2gg <- (dataF2 - m4[i])/sqrt(sigma4[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mi4[i]
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
  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3); 
  F3w1<-1/(12*m_sam[5])
  F3bmw <- matrix(0,m_sam[5],1); F3bmwsl <- matrix(0,m_sam[5],d3)
  for(i in 1:d3){
    F3gg <- (dataF3 - m5[i])/sqrt(sigma5[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mi5[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[5]){F3P2 <- F3P2+runif(m_sam[5])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[5]) + sum((F3P2 - (as.matrix(c(1:m_sam[5])) - 0.5)/m_sam[5])^2)
  F3u <- as.matrix(c(12*m_sam[5]*((F3dd[1]/m_sam[5]-0.5)^2),((45*m_sam[5])/4)*((F3dd[2]/m_sam[5]-1/3)^2),180*m_sam[5]*((F3dd[3]/m_sam[5]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[5]))))
  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("MX2-ADI-AD",round(abc,4),round(AIC,4),round(m[1],4),round(m[2],4),round(m[3],4), round(t(m4),4),round(sigma4[1],4),
                       round(t(mi4),4),round(t(m5),4),round(t(sigma5),4),round(t(mi5),4),round(sigma,4),          
                       round(B181111[1],4)," "," "," "," ",round(B181111[2],4),round(B181111[3],4),round(B181111[4],4),round(B181111[5],4),round(B181111[6],4),round(B181111[7],4),round(B181111[8],4),round(B181111[9],4),round(B181111[10],4),round(B181111[11],4),round(F2_jj,4),round(F2_ll*100,4),round(F2_gg,4),round(F2_rr*100,4),round(F3_jj,4),round(F3_ll*100,4),round(F3_gg,4),round(F3_rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
 output<-as.matrix(output)
  OUTPUT<-list(output,mix_pi4,mix_pi5)
  return(OUTPUT)
}
##################MX2-AD-AD(E-2)#############################
G5ModelFun[[20]] <- function(K1,logL,df11,df21,df31,df41,df51,G5text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataF2 <- as.matrix(as.numeric(df41[,1]))
  dataF3 <- as.matrix(as.numeric(df51[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2),sum(dataF3)))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2),mean(dataF3)))
  sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1);sigmaF2<-var(dataF2);sigmaF3<-var(dataF3)
  sumx6<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  sigma<-sumx6;sigma0<-sigma;sigma40<-sigmaF2;sigma400<-sigma40;sigma50<-sigmaF3;sigma500<-sigma50
  m_nf <- as.numeric(G5text2)
  ###############procedure start###########################
  mix_pi4<-as.matrix(c(0.0625,0.125,0.0625,0.1250,0.250,0.125,0.0625,0.125,0.0625))
  mix_pi5<-as.matrix(c(0.0625,0.125,0.0625,0.1250,0.250,0.125,0.0625,0.125,0.0625))
  a1<-sqrt(sigma400/m_sam[4])
  if (m[1]<m[3]) a1<--a1
  m4<-as.matrix(c(m[1],m[4]+2*a1,m[4]+1.2*a1,m[4]+0.6*a1,m[2],m[4]-0.6*a1,m[4]-1.2*a1,m[4]-2*a1,m[3]))
  a1<-sqrt(sigma500/m_sam[5])
  if (m[1]<m[3]) {a1<--a1}
  m5<-as.matrix(c(m[1],m[5]+2*a1,m[5]+1.2*a1,m[5]+0.6*a1,m[5],m[5]-0.6*a1,m[5]-1.2*a1,m[5]-2*a1,m[3]))
  nn<-m_sam
  abb<-sigma400/(sigma*5)
  sigma4<-matrix(0,9,1);sigma5<-matrix(0,9,1)
  sigma4[1]<-sigma400/abb
  if (sigma4[1]>sigma400) {sigma4[1]<-sigma400/1.25}
  sigma5[1]<-sigma500/abb
  if (sigma5[1]>sigma500) {sigma5[1]<-sigma500/1.25}
  #######first order parameters############
  hh19<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,-1,1,1,1,0,0,0,-1,-1,-1,1,1,1,0,0,0,-1,-1,-1,
                 1,0,-1,1,0,-1,1,0,-1,1,0,-1,1,0,-1,1,0,-1,1,0,-1,0,1,0,0,0,0,1,1,1,0,0,0,0,0,0,0.5,0.5,0.5,0,0,0,
                 0,1,0,0,1,0,0,1,0,0,1,0,0,0.5,0,0,0.5,0,0,0.5,0,1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                 0,1,0,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25),21,7)
  mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m4[4],m4[5],m4[6],m4[7],m4[8],m4[9],m5[1],m5[2],m5[3],m5[4],m5[5],m5[6],m5[7],m5[8],m5[9]))
  B19<-solve(crossprod(hh19,hh19))%*%crossprod(hh19,mm)
  g_aa1<-(0.5*B19[3]^2+0.25*B19[5]^2)/m_nf
  #   0.5*db**2+0.25*hb**2.
  g_aa2<-(0.5*B19[2]^2+0.25*B19[4]^2)/m_nf
  #   0.5*da**2+0.25*ha**2.
  g_aa3<-g_aa1+g_aa2
  sigma4[c(2:9)]<-sigma4[1];sigma5[c(3,7,9)]<-sigma5[1]
  sigma5[2]<-sigma5[1]+g_aa1;sigma5[4]<-sigma5[1]+g_aa2
  sigma5[5]<-sigma5[1]+g_aa3;sigma5[6]<-sigma5[1]+g_aa2
  sigma5[8]<-sigma5[1]+g_aa1
  mi4<-mix_pi4[c(1:9)];mi5<-mix_pi5[c(1:9)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
    sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
  ########################################################
  ##########iteration process###########
  ############E-step###################
  ############F2 process###############
  d3<-9
  iteration <- 0; stopa <- 1000
  W4 <- matrix(0,9,m_sam[4]); swx24 <- matrix(0,9,1)
  W5 <- matrix(0,9,m_sam[5]); swx25 <- matrix(0,9,1)
  hh<-matrix(0,14,14);b_line<-matrix(0,14,1)
  n0<-matrix(0,18,1); s0<-matrix(0,9,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:d3) { W4[i,] <- mix_pi4[i]*dnorm(dataF2,m4[i],sqrt(sigma4[i]))/dmixnorm(dataF2,m4,sqrt(sigma4),mix_pi4) }
    mi4 <- as.matrix(rowSums(W4)/m_sam[4]);n0[c(1:9)]<-as.matrix(rowSums(W4));sumwx4 <- W4%*%dataF2
    #############F3 process############
    for(i in 1:d3) { W5[i,] <- mix_pi5[i]*dnorm(dataF3,m5[i],sqrt(sigma5[i]))/dmixnorm(dataF3,m5,sqrt(sigma5),mix_pi5) }
    mi5 <- as.matrix(rowSums(W5)/m_sam[5]);n0[c(10:18)]<-as.matrix(rowSums(W5))
    n0[c(1:18)][abs(n0[c(1:18)])<0.000001]<-0.000001;sumwx5 <- W5%*%dataF3
    aaa0<-0
    aaa1<-1000;n_iter<-0;AA<-matrix(0,14,1)
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      ##########first order parameters###################
      mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m4[4],m4[5],m4[6],m4[7],m4[8],m4[9],m5[1],m5[2],m5[3],m5[4],m5[5],m5[6],m5[7],m5[8],m5[9]))
      B191<-solve(crossprod(hh19,hh19))%*%crossprod(hh19,mm)
      g_aa1<-(0.5*B191[3]^2+0.25*B191[5]^2)/m_nf
      #   0.5*db**2+0.25*hb**2.
      g_aa2<-(0.5*B191[2]^2+0.25*B191[4]^2)/m_nf
      #   0.5*da**2+0.25*ha**2.
      g_aa3<-g_aa1+g_aa2
      sigma4[c(2:9)]<-sigma4[1];sigma5[c(3,7,9)]<-sigma5[1]
      sigma5[2]<-sigma5[1]+g_aa1;sigma5[4]<-sigma5[1]+g_aa2
      sigma5[5]<-sigma5[1]+g_aa3;sigma5[6]<-sigma5[1]+g_aa2
      sigma5[8]<-sigma5[1]+g_aa1
      #################restrictions########################################
      hh[1,1]<-sigma*(1/nn[1]+4/nn[2]+1/nn[3])+sigma4[1]/n0[1]+4*sigma4[5]/n0[5]+sigma4[9]/n0[9]
      hh[1,2]<-sigma*2/nn[2]+2*sigma4[5]/n0[5]+2*sigma4[9]/n0[9]
      hh[1,3]<--sigma4[1]/n0[1]
      hh[1,4]<-0
      hh[1,5]<-sigma4[9]/n0[9]
      hh[1,6]<-0;hh[1,7]<-0
      hh[1,8]<-2*sigma4[5]/n0[5]
      hh[1,9]<--sigma4[1]/n0[1]
      hh[1,10]<-0
      hh[1,11]<-sigma*(2/nn[1]+4/nn[2]+2/nn[3])+4*sigma4[5]/n0[5]
      hh[1,12]<-sigma*(1/nn[1]+4/nn[2]+1/nn[3])
      hh[1,13]<-sigma/nn[1]+sigma/nn[3]-4*sigma4[5]/n0[5]
      hh[1,14]<-2*sigma/nn[2]+6*sigma4[5]/n0[5]
      hh[2,2]<-sigma/nn[2]+sigma4[5]/n0[5]+4*sigma4[9]/n0[9]+4*sigma5[9]/n0[18]
      hh[2,3]<-0;hh[2,4]<-0
      hh[2,5]<-2*sigma4[9]/n0[9]+2*sigma5[9]/n0[18]
      hh[2,6]<-0;hh[2,7]<-0
      hh[2,8]<-sigma4[5]/n0[5]
      hh[2,9]<--4*sigma5[9]/n0[18]
      hh[2,10]<-2*sigma5[9]/n0[18]
      hh[2,11]<-2*sigma/nn[2]+2*sigma4[5]/n0[5]-2*sigma5[9]/n0[18]
      hh[2,12]<-2*sigma/nn[2]
      hh[2,13]<--2*sigma4[5]/n0[5]
      hh[2,14]<-sigma/nn[2]+3*sigma4[5]/n0[5]
      hh[3,3]<-sigma4[1]/n0[1]+sigma4[3]/n0[3]+sigma5[1]/n0[10]+sigma5[3]/n0[12]
      hh[3,4]<--sigma4[3]/n0[3]-sigma5[3]/n0[12]
      hh[3,5]<-0;hh[3,6]<-0;hh[3,7]<-0;hh[3,8]<-0
      hh[3,9]<-sigma4[1]/n0[1]-2*sigma5[3]/n0[12]
      hh[3,10]<--sigma5[3]/n0[12]
      hh[3,11]<-sigma5[1]/n0[10]-sigma5[3]/n0[12]
      hh[3,12]<-0;hh[3,13]<-0;hh[3,14]<-0
      hh[4,4]<-sigma4[3]/n0[3]+sigma4[7]/n0[7]+sigma5[3]/n0[12]+sigma5[7]/n0[16]
      hh[4,5]<--sigma4[7]/n0[7]-sigma5[7]/n0[16]
      hh[4,6]<-0;hh[4,7]<-0
      hh[4,8]<-sigma4[7]/n0[7]
      hh[4,9]<--sigma4[7]/n0[7]+2*sigma5[3]/n0[12]
      hh[4,10]<-sigma5[3]/n0[12]
      hh[4,11]<-sigma5[3]/n0[12]-sigma5[7]/n0[16]
      hh[4,12]<-0;hh[4,13]<-0;hh[4,14]<-0
      hh[5,5]<-sigma4[7]/n0[7]+sigma4[9]/n0[9]+sigma5[7]/n0[16]+sigma5[9]/n0[18]
      hh[5,6]<-0;hh[5,7]<-0
      hh[5,8]<--sigma4[7]/n0[7]
      hh[5,9]<-sigma4[7]/n0[7]-2*sigma5[9]/n0[18]
      hh[5,10]<-sigma5[9]/n0[18]
      hh[5,11]<-sigma5[7]/n0[16]-sigma5[9]/n0[18]
      hh[5,12]<-0;hh[5,13]<-0;hh[5,14]<-0
      hh[6,6]<-sigma4[2]/n0[2]+sigma4[8]/n0[8]+sigma5[2]/n0[11]+sigma5[8]/n0[17]
      hh[6,7]<-0
      hh[6,8]<--sigma4[8]/n0[8]
      hh[6,9]<-0
      hh[6,10]<--sigma5[2]/n0[11]-sigma5[8]/n0[17]
      hh[6,11]<-0
      hh[6,12]<--sigma4[2]/n0[2]+sigma4[8]/n0[8]
      hh[6,13]<-0
      hh[6,14]<-sigma5[2]/n0[11]-sigma5[8]/n0[17]
      hh[7,7]<-sigma4[4]/n0[4]+sigma4[6]/n0[6]+sigma5[4]/n0[13]+sigma5[6]/n0[15]
      hh[7,8]<-sigma4[4]/n0[4]
      hh[7,9]<--2*sigma4[4]/n0[4]+4*sigma5[6]/n0[15]
      hh[7,10]<-0;hh[7,11]<-0
      hh[7,12]<--sigma4[4]/n0[4]+sigma4[6]/n0[6]
      hh[7,13]<-0
      hh[7,14]<-sigma5[4]/n0[13]-sigma5[6]/n0[15]
      hh[8,8]<-sigma4[4]/n0[4]+sigma4[5]/n0[5]+sigma4[7]/n0[7]+sigma4[8]/n0[8]
      hh[8,9]<--2*sigma4[4]/n0[4]-sigma4[7]/n0[7]
      hh[8,10]<-0
      hh[8,11]<-2*sigma4[5]/n0[5]
      hh[8,12]<--sigma4[4]/n0[4]-sigma4[8]/n0[8]
      hh[8,13]<--2*sigma4[5]/n0[5]
      hh[8,14]<-3*sigma4[5]/n0[5]
      hh[9,9]<-sigma4[1]/n0[1]+4*sigma4[4]/n0[4]+sigma4[7]/n0[7]+4*sigma5[3]/n0[12]+16*sigma5[6]/n0[15]+4*sigma5[9]/n0[18]
      hh[9,10]<-2*sigma5[3]/n0[12]-2*sigma5[9]/n0[18]
      hh[9,11]<-2*sigma5[3]/n0[12]+2*sigma5[9]/n0[18]
      hh[9,12]<-2*sigma4[4]/n0[4]
      hh[9,13]<-0
      hh[9,14]<--4*sigma5[6]/n0[15]
      hh[10,10]<-sigma5[2]/n0[11]+sigma5[3]/n0[12]+sigma5[8]/n0[17]+sigma5[9]/n0[18]
      hh[10,11]<-sigma5[3]/n0[12]-sigma5[9]/n0[18]
      hh[10,12]<-0;hh[10,13]<-0
      hh[10,14]<--sigma5[2]/n0[11]+sigma5[8]/n0[17]
      hh[11,11]<-sigma*(4/nn[1]+4/nn[2]+4/nn[3])+4*sigma4[5]/n0[5]+sigma5[1]/n0[10]+sigma5[3]/n0[12]+sigma5[7]/n0[16]+sigma5[9]/n0[18]
      hh[11,12]<-sigma*(2/nn[1]+4/nn[2]+2/nn[3])
      hh[11,13]<-sigma*(2/nn[1]+2/nn[3])-4*sigma4[5]/n0[5]
      hh[11,14]<-2*sigma/nn[2]+6*sigma4[5]/n0[5]
      hh[12,12]<-sigma*(1/nn[1]+4/nn[2]+1/nn[3])+sigma4[2]/n0[2]+sigma4[4]/n0[4]+sigma4[6]/n0[6]+sigma4[8]/n0[8]
      hh[12,13]<-sigma*(1/nn[1]+1/nn[3])
      hh[12,14]<-2*sigma/nn[2]
      hh[13,13]<-sigma/nn[1]+sigma/nn[3]+4*sigma4[5]/n0[5]+16*sigma5[5]/n0[14]
      hh[13,14]<--6*sigma4[5]/n0[5]-24*sigma5[5]/n0[14]
      hh[14,14]<-sigma/nn[2]+9*sigma4[5]/n0[5]+sigma5[2]/n0[11]+sigma5[4]/n0[13]+36*sigma5[5]/n0[14]+sigma5[6]/n0[15]+sigma5[8]/n0[17]
      for(i in 2:14)
      {
        for(j in 1:(i-1))
        {
          hh[i,j]<-hh[j,i]
        }
      }
      ########################################################
      b_line[1]<-sumx[1]/nn[1]+2*sumx[2]/nn[2]+sumx[3]/nn[3]-sumwx4[1]/n0[1]-2*sumwx4[5]/n0[5]-sumwx4[9]/n0[9]
      b_line[2]<-sumx[2]/nn[2]-sumwx4[5]/n0[5]-2*sumwx4[9]/n0[9]+2*sumwx5[9]/n0[18]
      b_line[3]<-sumwx4[1]/n0[1]-sumwx4[3]/n0[3]-sumwx5[1]/n0[10]+sumwx5[3]/n0[12]
      b_line[4]<-sumwx4[3]/n0[3]-sumwx4[7]/n0[7]-sumwx5[3]/n0[12]+sumwx5[7]/n0[16]
      b_line[5]<-sumwx4[7]/n0[7]-sumwx4[9]/n0[9]-sumwx5[7]/n0[16]+sumwx5[9]/n0[18]
      b_line[6]<-sumwx4[2]/n0[2]-sumwx4[8]/n0[8]-sumwx5[2]/n0[11]+sumwx5[8]/n0[17]
      b_line[7]<-sumwx4[4]/n0[4]-sumwx4[6]/n0[6]-sumwx5[4]/n0[13]+sumwx5[6]/n0[15]
      b_line[8]<-sumwx4[4]/n0[4]-sumwx4[5]/n0[5]-sumwx4[7]/n0[7]+sumwx4[8]/n0[8]
      b_line[9]<-sumwx4[1]/n0[1]-2*sumwx4[4]/n0[4]+sumwx4[7]/n0[7]-2*sumwx5[3]/n0[12]+4*sumwx5[6]/n0[15]-2*sumwx5[9]/n0[18]
      b_line[10]<-sumwx5[2]/n0[11]-sumwx5[3]/n0[12]-sumwx5[8]/n0[17]+sumwx5[9]/n0[18]
      b_line[11]<-2*sumx[1]/nn[1]+2*sumx[2]/nn[2]+2*sumx[3]/nn[3]-2*sumwx4[5]/n0[5]-sumwx5[1]/n0[10]-sumwx5[3]/n0[12]-sumwx5[7]/n0[16]-sumwx5[9]/n0[18]
      b_line[12]<-sumx[1]/nn[1]+2*sumx[2]/nn[2]+sumx[3]/nn[3]-sumwx4[2]/n0[2]-sumwx4[4]/n0[4]-sumwx4[6]/n0[6]-sumwx4[8]/n0[8]
      b_line[13]<-sumx[1]/nn[1]+sumx[3]/nn[3]+2*sumwx4[5]/n0[5]-4*sumwx5[5]/n0[14]
      b_line[14]<-sumx[2]/nn[2]-3*sumwx4[5]/n0[5]-sumwx5[2]/n0[11]-sumwx5[4]/n0[13]-sumwx5[6]/n0[15]-sumwx5[8]/n0[17]+6*sumwx5[5]/n0[14]
      B1911<-solve(hh,b_line)
      ##########################################################
      m[1]<-(sumx[1]-sigma*(B1911[1]+2*B1911[11]+B1911[12]+B1911[13]))/nn[1]
      m[2]<-(sumx[2]-sigma*(2*B1911[1]+B1911[2]+2*B1911[11]+2*B1911[12]+B1911[14]))/nn[2]
      m[3]<-(sumx[3]-sigma*(B1911[1]+2*B1911[11]+B1911[12]+B1911[13]))/nn[3]
      m4[1]<-(sumwx4[1]+sigma4[1]*(B1911[1]-B1911[3]-B1911[9]))/n0[1]
      m4[2]<-(sumwx4[2]+sigma4[2]*(-B1911[6]+B1911[12]))/n0[2]
      m4[3]<-(sumwx4[3]+sigma4[3]*(B1911[3]-B1911[4]))/n0[3]
      m4[4]<-(sumwx4[4]+sigma4[4]*(-B1911[7]-B1911[8]+2*B1911[9]+B1911[12]))/n0[4]
      m4[5]<-(sumwx4[5]+sigma4[5]*(2*B1911[1]+B1911[2]+B1911[8]+2*B1911[11]-2*B1911[13]+3*B1911[14]))/n0[5]
      m4[6]<-(sumwx4[6]+sigma4[6]*(B1911[7]+B1911[12]))/n0[6]
      m4[7]<-(sumwx4[7]+sigma4[7]*(B1911[4]-B1911[5]+B1911[8]-B1911[9]))/n0[7]
      m4[8]<-(sumwx4[8]+sigma4[8]*(B1911[6]-B1911[8]+B1911[12]))/n0[8]
      m4[9]<-(sumwx4[9]+sigma4[9]*(B1911[1]+2*B1911[2]+B1911[5]))/n0[9]
      m5[1]<-(sumwx5[1]+sigma5[1]*(B1911[3]+B1911[11]))/n0[10]
      m5[2]<-(sumwx5[2]+sigma5[2]*(B1911[6]-B1911[10]+B1911[14]))/n0[11]
      m5[3]<-(sumwx5[3]+sigma5[3]*(-B1911[3]+B1911[4]+2*B1911[9]+B1911[10]+B1911[11]))/n0[12]
      m5[4]<-(sumwx5[4]+sigma5[4]*(B1911[7]+B1911[14]))/n0[13]
      m5[5]<-(sumwx5[5]+sigma5[5]*(4*B1911[13]-6*B1911[14]))/n0[14]
      m5[6]<-(sumwx5[6]+sigma5[6]*(-B1911[7]-4*B1911[9]+B1911[14]))/n0[15]
      m5[7]<-(sumwx5[7]+sigma5[7]*(-B1911[4]+B1911[5]+B1911[11]))/n0[16]
      m5[8]<-(sumwx5[8]+sigma5[8]*(-B1911[6]+B1911[10]+B1911[14]))/n0[17]
      m5[9]<-(sumwx5[9]-sigma5[9]*(2*B1911[2]+B1911[5]-2*B1911[9]+B1911[10]-B1911[11]))/n0[18]
      aaa1<-max(abs(B1911-AA))
      AA<-B1911
      if (n_iter>20) break
    } 
    ##########obtain variance###############
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d3) {swx24[i] <- W4[i,]%*%(dataF2-m4[i])^2 } ;for(i in 1:d3) {swx25[i] <- W5[i,]%*%(dataF3-m5[i])^2 } 
    ab1<-sum(swx24)
    sigma4[1]<-ab1/nn[4];sigma40<-sigma4[1]-sigma
    if (sigma40<0) {sigma40<-0;sigma4[1]<-sigma}
    sigma4[c(2:9)]<-sigma4[1]
    aaa0<-sigma5[1]
    ab2<-swx25[1]+swx25[3]+swx25[7]+swx25[9];ab3<-n0[10]+n0[12]+n0[16]+n0[18]
    mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m4[4],m4[5],m4[6],m4[7],m4[8],m4[9],m5[1],m5[2],m5[3],m5[4],m5[5],m5[6],m5[7],m5[8],m5[9]))
    B19111<-solve(crossprod(hh19,hh19))%*%crossprod(hh19,mm)
    g_aa1<-(0.5*B19111[3]^2+0.25*B19111[5]^2)/m_nf
    #   0.5*db**2+0.25*hb**2.
    g_aa2<-(0.5*B19111[2]^2+0.25*B19111[4]^2)/m_nf
    #   0.5*da**2+0.25*ha**2.
    g_aa3<-g_aa1+g_aa2
    n_iter<-0;aaa1<-1000
    while (aaa1>0.0001){
      n_iter<-n_iter+1
      aa1<-sigma5[1]/(sigma5[1]+g_aa1);aa2<-sigma5[1]/(sigma5[1]+g_aa2);aa3<-sigma5[1]/(sigma5[1]+g_aa3)
      as3<-ab2+aa1^2*(swx25[2]+swx25[8])+aa2^2*(swx25[4]+swx25[6])+aa3^2*swx25[5]
      as4<-ab3+aa1*(n0[11]+n0[17])+aa2*(n0[13]+n0[15])+aa3*n0[14]
      sigma5[1]<-as3/as4
      aaa1<-abs(sigma5[1]-aaa0)
      aaa0<-sigma5[1]
      if (n_iter>20) break 
    }
    sigma50<-sigma5[1]-sigma/m_nf
    if (sigma50<0) {sigma50<-0;sigma5[1]<-sigma/m_nf}
    sigma5[c(3,7,9)]<-sigma5[1]
    sigma5[2]<-sigma5[1]+g_aa1;sigma5[4]<-sigma5[1]+g_aa2
    sigma5[5]<-sigma5[1]+g_aa3;sigma5[6]<-sigma5[1]+g_aa2
    sigma5[8]<-sigma5[1]+g_aa1
    ab2<-ss1+ss2+ss3;ab3<-nn[1]+nn[2]+nn[3]
    n_iter<-0;aaa0<-sigma
    an1<-sigma40;an2<-sigma50
    aaa1<-1000
    while (aaa1>0.01){
      n_iter<-n_iter+1  
      if (an1<0) {an1<-0};if (an2<0) {an2<-0}
      ab4<-sigma/(sigma40+sigma)
      s0[1]<-(sigma/m_nf)/(sigma/m_nf+an2)
      s0[2]<-(sigma/m_nf)/(sigma/m_nf+an2+g_aa1)
      s0[3]<-(sigma/m_nf)/(sigma/m_nf+an2)
      s0[4]<-(sigma/m_nf)/(sigma/m_nf+an2+g_aa2)
      s0[5]<-(sigma/m_nf)/(sigma/m_nf+an2+g_aa3)
      s0[6]<-(sigma/m_nf)/(sigma/m_nf+an2+g_aa2)
      s0[7]<-(sigma/m_nf)/(sigma/m_nf+an2)
      s0[8]<-(sigma/m_nf)/(sigma/m_nf+an2+g_aa1)
      s0[9]<-(sigma/m_nf)/(sigma/m_nf+an2)
      as1<-sum(s0[c(1:9)]^2*swx25[c(1:9)])
      as2<-sum(s0[c(1:9)]*n0[c(10:18)])
      as3<-sum(swx24[c(1:9)])
      sigma<-(ab2+ab4*ab4*as3+as1*m_nf)/(ab3+ab4*nn[4]+as2)
      aaa1<-abs(sigma-aaa0)
      aaa0<-sigma
      if (n_iter>20) break
    }
    if (sigma<sigma0) {sigma<-sigma0}
    sigma4[1]<-sigma+sigma40;sigma4[c(2:9)]<-sigma4[1]
    sigma5[1]<-sigma/m_nf+sigma50;sigma5[c(3,7,9)]<-sigma5[1]
    sigma5[2]<-sigma5[1]+g_aa1;sigma5[4]<-sigma5[1]+g_aa2
    sigma5[5]<-sigma5[1]+g_aa3;sigma5[6]<-sigma5[1]+g_aa2
    sigma5[8]<-sigma5[1]+g_aa1
    ##############################criteria for iterations to stop######################## 
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
      sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*10
  #######first order parameters############
  mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m4[4],m4[5],m4[6],m4[7],m4[8],m4[9],m5[1],m5[2],m5[3],m5[4],m5[5],m5[6],m5[7],m5[8],m5[9]))
  B191111<-solve(crossprod(hh19,hh19))%*%crossprod(hh19,mm)
  ########second order genetic parameters###############   
  F2_jj <- sigmaF2 - sigma4[1]
  F2_gg <- sigma4[1]-sigma
  if(F2_jj<0) {F2_jj<-0}
  if(F2_gg<0 || F2_gg>sigmaF2) {F2_gg<-0}
  F2_ll <- F2_jj/sigmaF2
  F2_rr <- F2_gg/sigmaF2
  F3_jj <- sigmaF3 - sigma5[1]
  F3_gg <- sigma5[1]-sigma/m_nf
  if(F3_jj<0) {F3_jj<-0}
  if(F3_gg<0 || F3_gg>sigmaF3) {F3_gg<-0}
  F3_ll <- F3_jj/sigmaF3
  F3_rr <- F3_gg/sigmaF3
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma))
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
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma))
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
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma))
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
    F2gg <- (dataF2 - m4[i])/sqrt(sigma4[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mi4[i]
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
  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3); 
  F3w1<-1/(12*m_sam[5])
  F3bmw <- matrix(0,m_sam[5],1); F3bmwsl <- matrix(0,m_sam[5],d3)
  for(i in 1:d3){
    F3gg <- (dataF3 - m5[i])/sqrt(sigma5[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mi5[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[5]){F3P2 <- F3P2+runif(m_sam[5])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[5]) + sum((F3P2 - (as.matrix(c(1:m_sam[5])) - 0.5)/m_sam[5])^2)
  F3u <- as.matrix(c(12*m_sam[5]*((F3dd[1]/m_sam[5]-0.5)^2),((45*m_sam[5])/4)*((F3dd[2]/m_sam[5]-1/3)^2),180*m_sam[5]*((F3dd[3]/m_sam[5]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[5]))))
  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("MX2-AD-AD",round(abc,4),round(AIC,4),round(m[1],4),round(m[2],4),round(m[3],4), round(t(m4),4),round(sigma4[1],4),
                       round(t(mi4),4),round(t(m5),4),round(t(sigma5),4),round(t(mi5),4),round(sigma,4),          
                       round(B191111[1],4)," "," "," "," ",round(B191111[2],4),round(B191111[3],4),round(B191111[4],4),round(B191111[5],4)," "," "," "," ",round(B191111[6],4),round(B191111[7],4),round(F2_jj,4),round(F2_ll*100,4),round(F2_gg,4),round(F2_rr*100,4),round(F3_jj,4),round(F3_ll*100,4),round(F3_gg,4),round(F3_rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
 output<-as.matrix(output)
  OUTPUT<-list(output,mix_pi4,mix_pi5)
  return(OUTPUT)
}
########################MX2-A-AD(E-3)##################################
G5ModelFun[[21]] <- function(K1,logL,df11,df21,df31,df41,df51,G5text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataF2 <- as.matrix(as.numeric(df41[,1]))
  dataF3 <- as.matrix(as.numeric(df51[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2),sum(dataF3)))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2),mean(dataF3)))
  sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1);sigmaF2<-var(dataF2);sigmaF3<-var(dataF3)
  sumx6<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  sigma<-sumx6;sigma0<-sigma;sigma40<-sigmaF2;sigma400<-sigma40;sigma50<-sigmaF3;sigma500<-sigma50
  m_nf <- as.numeric(G5text2)
  ###############procedure start###########################
  mix_pi4<-as.matrix(c(0.0625,0.125,0.0625,0.1250,0.250,0.125,0.0625,0.125,0.0625))
  mix_pi5<-as.matrix(c(0.0625,0.125,0.0625,0.1250,0.250,0.125,0.0625,0.125,0.0625))
  a1<-sqrt(sigma400/m_sam[4])
  if (m[1]<m[3]) a1<--a1
  m4<-as.matrix(c(m[1],m[4]+2*a1,m[4]+1.2*a1,m[4]+0.6*a1,m[2],m[4]-0.6*a1,m[4]-1.2*a1,m[4]-2*a1,m[3]))
  a1<-sqrt(sigma500/m_sam[5])
  if (m[1]<m[3]) {a1<--a1}
  m5<-as.matrix(c(m[1],m[5]+2*a1,m[5]+1.2*a1,m[5]+0.6*a1,m[2],m[5]-0.6*a1,m[5]-1.2*a1,m[5]-2*a1,m[3]))
  nn<-m_sam
  abb<-sigma400/(sigma*5)
  sigma4<-matrix(0,9,1);sigma5<-matrix(0,9,1)
  sigma4[1]<-sigma400/abb
  if (sigma4[1]>sigma400) {sigma4[1]<-sigma400/2}
  sigma5[1]<-sigma500/abb
  if (sigma5[1]>sigma500) {sigma5[1]<-sigma500/2}
  #######first order parameters############
  hh20<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,-1,1,1,1,0,0,0,-1,-1,-1,1,1,1,0,0,0,-1,-1,-1,
                 1,0,-1,1,0,-1,1,0,-1,1,0,-1,1,0,-1,1,0,-1,1,0,-1,1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                 0,1,0,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25),21,5)
  mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m4[4],m4[5],m4[6],m4[7],m4[8],m4[9],m5[1],m5[2],m5[3],m5[4],m5[5],m5[6],m5[7],m5[8],m5[9]))
  B20<-solve(crossprod(hh20,hh20))%*%crossprod(hh20,mm)
  g_aa1<-0.5*B20[3]^2/m_nf        #  0.5*db**2.
  g_aa2<-0.5*B20[2]^2/m_nf        #   0.5*da**2.
  g_aa3<-g_aa1+g_aa2
  sigma4[c(2:9)]<-sigma4[1];sigma5[c(3,7,9)]<-sigma5[1]
  sigma5[2]<-sigma5[1]+g_aa1;sigma5[4]<-sigma5[1]+g_aa2
  sigma5[5]<-sigma5[1]+g_aa3;sigma5[6]<-sigma5[1]+g_aa2
  sigma5[8]<-sigma5[1]+g_aa1
  mi4<-mix_pi4[c(1:9)];mi5<-mix_pi5[c(1:9)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
    sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
  ########################################################
  ##########iteration process###########
  ############E-step###################
  ############F2 process###############
  d3<-9
  iteration <- 0; stopa <- 1000
  W4 <- matrix(0,9,m_sam[4]); swx24 <- matrix(0,9,1)
  W5 <- matrix(0,9,m_sam[5]); swx25 <- matrix(0,9,1)
  hh<-matrix(0,16,16);b_line<-matrix(0,16,1)
  n0<-matrix(0,18,1);s0<-matrix(0,9,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:d3) { W4[i,] <- mix_pi4[i]*dnorm(dataF2,m4[i],sqrt(sigma4[i]))/dmixnorm(dataF2,m4,sqrt(sigma4),mix_pi4) }
    mi4 <- as.matrix(rowSums(W4)/m_sam[4]);n0[c(1:9)]<-as.matrix(rowSums(W4));sumwx4 <- W4%*%dataF2
    #############F3 process############
    for(i in 1:d3) { W5[i,] <- mix_pi5[i]*dnorm(dataF3,m5[i],sqrt(sigma5[i]))/dmixnorm(dataF3,m5,sqrt(sigma5),mix_pi5) }
    mi5 <- as.matrix(rowSums(W5)/m_sam[5]);n0[c(10:18)]<-as.matrix(rowSums(W5))
    n0[c(1:18)][abs(n0[c(1:18)])<0.000001]<-0.000001;sumwx5 <- W5%*%dataF3
    aaa0<-0
    aaa1<-1000;n_iter<-0;AA<-matrix(0,16,1)
    while (aaa1>0.0001)
    {
      n_iter<-n_iter+1
      ##########first order parameters###################
      mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m4[4],m4[5],m4[6],m4[7],m4[8],m4[9],m5[1],m5[2],m5[3],m5[4],m5[5],m5[6],m5[7],m5[8],m5[9]))
      B201<-solve(crossprod(hh20,hh20))%*%crossprod(hh20,mm)
      g_aa1<-0.5*B201[3]^2/m_nf        #   0.5*db**2.
      g_aa2<-0.5*B201[2]^2/m_nf        #   0.5*da**2.
      g_aa3<-g_aa1+g_aa2
      sigma4[c(2:9)]<-sigma4[1];sigma5[c(3,7,9)]<-sigma5[1]
      sigma5[2]<-sigma5[1]+g_aa1;sigma5[4]<-sigma5[1]+g_aa2
      sigma5[5]<-sigma5[1]+g_aa3;sigma5[6]<-sigma5[1]+g_aa2
      sigma5[8]<-sigma5[1]+g_aa1
      ################restrictions########################################
      hh[1,1]<-sigma4[1]/n0[1]+4*sigma4[2]/n0[2]+sigma4[3]/n0[3]
      hh[1,2]<-0;hh[1,3]<-0;hh[1,4]<-0;hh[1,5]<-0;hh[1,6]<-0
      hh[1,7]<-sigma4[1]/n0[1]
      hh[1,8]<--2*sigma4[2]/n0[2]
      hh[1,9]<-0;hh[1,10]<-0;hh[1,11]<-0;hh[1,12]<-0
      hh[1,13]<-sigma4[1]/n0[1]-sigma4[3]/n0[3]
      hh[1,14]<--2*sigma4[2]/n0[2]
      hh[1,15]<-0
      hh[1,16]<--8*sigma4[2]/n0[2]
      hh[2,2]<-sigma4[4]/n0[4]+4*sigma4[5]/n0[5]+sigma4[6]/n0[6]
      hh[2,3]<-0;hh[2,4]<-0;hh[2,5]<-0;hh[2,6]<-0
      hh[2,7]<--2*sigma4[4]/n0[4]
      hh[2,8]<-4*sigma4[5]/n0[5]
      hh[2,9]<-0;hh[2,10]<-0
      hh[2,11]<-8*sigma4[5]/n0[5]
      hh[2,12]<-6*sigma4[5]/n0[5]
      hh[2,13]<-0;hh[2,14]<-0
      hh[2,15]<-sigma4[4]/n0[4]-sigma4[6]/n0[6]
      hh[2,16]<-0
      hh[3,3]<-sigma4[7]/n0[7]+4*sigma4[8]/n0[8]+sigma4[9]/n0[9]
      hh[3,4]<-0;hh[3,5]<-0;hh[3,6]<-0
      hh[3,7]<-sigma4[7]/n0[7]
      hh[3,8]<--2*sigma4[8]/n0[8]
      hh[3,9]<-0;hh[3,10]<-0;hh[3,11]<-0;hh[3,12]<-0;hh[3,13]<-0
      hh[3,14]<-2*sigma4[8]/n0[8]
      hh[3,15]<-0;hh[3,16]<-0
      hh[4,4]<-sigma5[1]/n0[10]+4*sigma5[2]/n0[11]+sigma5[3]/n0[12]
      hh[4,5]<-0;hh[4,6]<-0;hh[4,7]<-0;hh[4,8]<-0
      hh[4,9]<-sigma5[1]/n0[10]
      hh[4,10]<--2*sigma5[2]/n0[11]
      hh[4,11]<-0;hh[4,12]<-0;hh[4,13]<-0
      hh[4,14]<-2*sigma5[2]/n0[11]
      hh[4,15]<-0
      hh[4,16]<-10*sigma5[2]/n0[11]
      hh[5,5]<-sigma5[4]/n0[13]+4*sigma5[5]/n0[14]+sigma5[6]/n0[15]
      hh[5,6]<-0;hh[5,7]<-0;hh[5,8]<-0
      hh[5,9]<--2*sigma5[4]/n0[13]
      hh[5,10]<-4*sigma5[5]/n0[14]
      hh[5,11]<-0
      hh[5,12]<--4*sigma5[5]/n0[14]
      hh[5,13]<-0;hh[5,14]<-0
      hh[5,15]<--sigma5[4]/n0[13]+sigma5[6]/n0[15]
      hh[5,16]<--sigma5[4]/n0[13]-sigma5[6]/n0[15]
      hh[6,6]<-sigma5[7]/n0[16]+4*sigma5[8]/n0[17]+sigma5[9]/n0[18]
      hh[6,7]<-0;hh[6,8]<-0
      hh[6,9]<-sigma5[7]/n0[16]
      hh[6,10]<--2*sigma5[8]/n0[17]
      hh[6,11]<-0;hh[6,12]<-0
      hh[6,13]<--sigma5[7]/n0[16]+sigma5[9]/n0[18]
      hh[6,14]<--2*sigma5[8]/n0[17]
      hh[6,15]<-0
      hh[6,16]<-2*sigma5[8]/n0[17]
      hh[7,7]<-sigma4[1]/n0[1]+4*sigma4[4]/n0[4]+sigma4[7]/n0[7]
      hh[7,8]<-0;hh[7,9]<-0;hh[7,10]<-0;hh[7,11]<-0;hh[7,12]<-0
      hh[7,13]<-sigma4[1]/n0[1]
      hh[7,14]<-0
      hh[7,15]<--2*sigma4[4]/n0[4]
      hh[7,16]<-0
      hh[8,8]<-sigma4[2]/n0[2]+4*sigma4[5]/n0[5]+sigma4[8]/n0[8]
      hh[8,9]<-0;hh[8,10]<-0
      hh[8,11]<-8*sigma4[5]/n0[5]
      hh[8,12]<-6*sigma4[5]/n0[5]
      hh[8,13]<-0
      hh[8,14]<-sigma4[2]/n0[2]-sigma4[8]/n0[8]
      hh[8,15]<-0
      hh[8,16]<-4*sigma4[2]/n0[2]
      hh[9,9]<-sigma5[1]/n0[10]+4*sigma5[4]/n0[13]+sigma5[7]/n0[16]
      hh[9,10]<-0;hh[9,11]<-0;hh[9,12]<-0
      hh[9,13]<--sigma5[7]/n0[16]
      hh[9,14]<-0
      hh[9,15]<-2*sigma5[4]/n0[13]
      hh[9,16]<-2*sigma5[4]/n0[13]
      hh[10,10]<-sigma5[2]/n0[11]+4*sigma5[5]/n0[14]+sigma5[8]/n0[17]
      hh[10,11]<-0
      hh[10,12]<--4*sigma5[5]/n0[14]
      hh[10,13]<-0
      hh[10,14]<--sigma5[2]/n0[11]+sigma5[8]/n0[17]
      hh[10,15]<-0
      hh[10,16]<--5*sigma5[2]/n0[11]-sigma5[8]/n0[17]
      hh[11,11]<-sigma/nn[1]+4*sigma/nn[2]+sigma/nn[3]+16*sigma4[5]/n0[5]
      hh[11,12]<-sigma/nn[2]+12*sigma4[5]/n0[5]
      hh[11,13]<-0;hh[11,14]<-0;hh[11,15]<-0
      hh[11,16]<-2*sigma/nn[1]+2*sigma/nn[3]
      hh[12,12]<-sigma/nn[2]+9*sigma4[5]/n0[5]+4*sigma5[5]/n0[14]
      hh[12,13]<-0;hh[12,14]<-0;hh[12,15]<-0
      hh[12,16]<-0
      hh[13,13]<-sigma4[1]/n0[1]+sigma4[3]/n0[3]+sigma5[7]/n0[16]+sigma5[9]/n0[18]
      hh[13,14]<-0;hh[13,15]<-0
      hh[13,16]<-sigma5[7]/n0[16]
      hh[14,14]<-sigma4[2]/n0[2]+sigma4[8]/n0[8]+sigma5[2]/n0[11]+sigma5[8]/n0[17]
      hh[14,15]<-0
      hh[14,16]<-4*sigma4[2]/n0[2]+5*sigma5[2]/n0[11]-sigma5[8]/n0[17]
      hh[15,15]<-sigma4[4]/n0[4]+sigma4[6]/n0[6]+sigma5[4]/n0[13]+sigma5[6]/n0[15]
      hh[15,16]<-sigma5[4]/n0[13]-sigma5[6]/n0[15]
      hh[16,16]<-4*sigma/nn[1]+4*sigma/nn[3]+16*sigma4[2]/n0[2]+25*sigma5[2]/n0[11]+sigma5[4]/n0[13]+sigma5[6]/n0[15]+sigma5[8]/n0[17]
      for(i in 2:16)
      {
        for(j in 1:(i-1))
        {
          hh[i,j]<-hh[j,i]
        }
      }
      ########################################################
      b_line[1]<-sumwx4[1]/n0[1]-2*sumwx4[2]/n0[2]+sumwx4[3]/n0[3]
      b_line[2]<-sumwx4[4]/n0[4]-2*sumwx4[5]/n0[5]+sumwx4[6]/n0[6]
      b_line[3]<-sumwx4[7]/n0[7]-2*sumwx4[8]/n0[8]+sumwx4[9]/n0[9]
      b_line[4]<-sumwx5[1]/n0[10]-2*sumwx5[2]/n0[11]+sumwx5[3]/n0[12]
      b_line[5]<-sumwx5[4]/n0[13]-2*sumwx5[5]/n0[14]+sumwx5[6]/n0[15]
      b_line[6]<-sumwx5[7]/n0[16]-2*sumwx5[8]/n0[17]+sumwx5[9]/n0[18]
      b_line[7]<-sumwx4[1]/n0[1]-2*sumwx4[4]/n0[4]+sumwx4[7]/n0[7]
      b_line[8]<-sumwx4[2]/n0[2]-2*sumwx4[5]/n0[5]+sumwx4[8]/n0[8]
      b_line[9]<-sumwx5[1]/n0[10]-2*sumwx5[4]/n0[13]+sumwx5[7]/n0[16]
      b_line[10]<-sumwx5[2]/n0[11]-2*sumwx5[5]/n0[14]+sumwx5[8]/n0[17]
      b_line[11]<-sumx[1]/nn[1]+2*sumx[2]/nn[2]+sumx[3]/nn[3]-4*sumwx4[5]/n0[5]
      b_line[12]<-sumx[2]/nn[2]-3*sumwx4[5]/n0[5]+2*sumwx5[5]/n0[14]
      b_line[13]<-sumwx4[1]/n0[1]-sumwx4[3]/n0[3]-sumwx5[7]/n0[16]+sumwx5[9]/n0[18]
      b_line[14]<-sumwx4[2]/n0[2]-sumwx4[8]/n0[8]-sumwx5[2]/n0[11]+sumwx5[8]/n0[17]
      b_line[15]<-sumwx4[4]/n0[4]-sumwx4[6]/n0[6]-sumwx5[4]/n0[13]+sumwx5[6]/n0[15]
      b_line[16]<-2*sumx[1]/nn[1]+2*sumx[3]/nn[3]+4*sumwx4[2]/n0[2]-5*sumwx5[2]/n0[11]-sumwx5[4]/n0[13]-sumwx5[6]/n0[15]-sumwx5[8]/n0[17]
      B2011<-solve(hh,b_line)
      ##########################################################
      m[1]<-(sumx[1]+sigma*(-B2011[11]-2*B2011[16]))/nn[1]
      m[2]<-(sumx[2]-sigma*(2*B2011[11]+B2011[12]))/nn[2]
      m[3]<-(sumx[3]+sigma*(-B2011[11]-2*B2011[16]))/nn[3]
      m4[1]<-(sumwx4[1]-sigma4[1]*(B2011[1]+B2011[7]+B2011[13]))/n0[1]
      m4[2]<-(sumwx4[2]+sigma4[1]*(2*B2011[1]-B2011[8]-B2011[14]-4*B2011[16]))/n0[2]
      m4[3]<-(sumwx4[3]+sigma4[1]*(-B2011[1]+B2011[13]))/n0[3]
      m4[4]<-(sumwx4[4]+sigma4[1]*(-B2011[2]+2*B2011[7]-B2011[15]))/n0[4]
      m4[5]<-(sumwx4[5]+sigma4[1]*(2*B2011[2]+2*B2011[8]+4*B2011[11]+3*B2011[12]))/n0[5]
      m4[6]<-(sumwx4[6]+sigma4[1]*(-B2011[2]+B2011[15]))/n0[6]
      m4[7]<-(sumwx4[7]-sigma4[1]*(B2011[3]+B2011[7]))/n0[7]
      m4[8]<-(sumwx4[8]+sigma4[1]*(2*B2011[3]-B2011[8]+B2011[14]))/n0[8]
      m4[9]<-(sumwx4[9]-sigma4[1]*B2011[3])/n0[9]
      m5[1]<-(sumwx5[1]-sigma5[1]*(B2011[4]+B2011[9]))/n0[10]
      m5[2]<-(sumwx5[2]+sigma5[2]*(2*B2011[4]-B2011[10]+B2011[14]+5*B2011[16]))/n0[11]
      m5[3]<-(sumwx5[3]-sigma5[3]*B2011[4])/n0[12]
      m5[4]<-(sumwx5[4]+sigma5[4]*(-B2011[5]+2*B2011[9]+B2011[15]+B2011[16]))/n0[13]
      m5[5]<-(sumwx5[5]+sigma5[5]*(2*B2011[5]+2*B2011[10]-2*B2011[12]))/n0[14]
      m5[6]<-(sumwx5[6]-sigma5[6]*(B2011[5]+B2011[15]-B2011[16]))/n0[15]
      m5[7]<-(sumwx5[7]+sigma5[7]*(-B2011[6]-B2011[9]+B2011[13]))/n0[16]
      m5[8]<-(sumwx5[8]+sigma5[8]*(2*B2011[6]-B2011[10]-B2011[14]+B2011[16]))/n0[17]
      m5[9]<-(sumwx5[9]-sigma5[9]*(B2011[6]+B2011[13]))/n0[18]
      aaa1<-max(abs(B2011-AA))
      AA<-B2011
      if (n_iter>20) break
    } 
    ##########obtain variance###############
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d3) {swx24[i] <- W4[i,]%*%(dataF2-m4[i])^2 } ;for(i in 1:d3) {swx25[i] <- W5[i,]%*%(dataF3-m5[i])^2 } 
    ab1<-sum(swx24);sigma4[1]<-ab1/nn[4];sigma40<-sigma4[1]-sigma
    if (sigma40<0) {sigma40<-0;sigma4[1]<-sigma};sigma4[c(2:9)]<-sigma4[1]
    aaa0<-sigma5[1]
    ab2<-swx25[1]+swx25[3]+swx25[7]+swx25[9];ab3<-n0[10]+n0[12]+n0[16]+n0[18]
    mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m4[4],m4[5],m4[6],m4[7],m4[8],m4[9],m5[1],m5[2],m5[3],m5[4],m5[5],m5[6],m5[7],m5[8],m5[9]))
    B20111<-solve(crossprod(hh20,hh20))%*%crossprod(hh20,mm)
    g_aa1<-0.5*B20111[3]^2/m_nf        #   0.5*db**2.
    g_aa2<-0.5*B20111[2]^2/m_nf        #   0.5*da**2.
    g_aa3<-g_aa1+g_aa2
    n_iter<-0;aaa1<-1000
    while (aaa1>0.0001){
      n_iter<-n_iter+1
      aa1<-sigma5[1]/(sigma5[1]+g_aa1);aa2<-sigma5[1]/(sigma5[1]+g_aa2);aa3<-sigma5[1]/(sigma5[1]+g_aa3)
      as3<-ab2+aa1^2*(swx25[2]+swx25[8])+aa2^2*(swx25[4]+swx25[6])+aa3^2*swx25[5]
      as4<-ab3+aa1*(n0[11]+n0[17])+aa2*(n0[13]+n0[15])+aa3*n0[14]
      sigma5[1]<-as3/as4
      aaa1<-abs(sigma5[1]-aaa0)
      aaa0<-sigma5[1]
      if (n_iter>20) break 
    }
    sigma50<-sigma5[1]-sigma/m_nf
    if (sigma50<0) {sigma50<-0;sigma5[1]<-sigma/m_nf}
    sigma5[c(3,7,9)]<-sigma5[1];sigma5[2]<-sigma5[1]+g_aa1
    sigma5[4]<-sigma5[1]+g_aa2;sigma5[5]<-sigma5[1]+g_aa3
    sigma5[6]<-sigma5[1]+g_aa2;sigma5[8]<-sigma5[1]+g_aa1
    ab2<-ss1+ss2+ss3;ab3<-nn[1]+nn[2]+nn[3]
    n_iter<-0;aaa0<-sigma;an1<-sigma40;an2<-sigma50
    aaa1<-1000
    while (aaa1>0.01){
      n_iter<-n_iter+1  
      if (an1<0) {an1<-0};if (an2<0) {an2<-0}
      ab4<-sigma/(an1+sigma)
      s0[1]<-(sigma/m_nf)/(sigma/m_nf+an2)
      s0[2]<-(sigma/m_nf)/(sigma/m_nf+an2+g_aa1)
      s0[3]<-(sigma/m_nf)/(sigma/m_nf+an2)
      s0[4]<-(sigma/m_nf)/(sigma/m_nf+an2+g_aa2)
      s0[5]<-(sigma/m_nf)/(sigma/m_nf+an2+g_aa3)
      s0[6]<-(sigma/m_nf)/(sigma/m_nf+an2+g_aa2)
      s0[7]<-(sigma/m_nf)/(sigma/m_nf+an2)
      s0[8]<-(sigma/m_nf)/(sigma/m_nf+an2+g_aa1)
      s0[9]<-(sigma/m_nf)/(sigma/m_nf+an2)
      as1<-sum(s0[c(1:9)]^2*swx25[c(1:9)])
      as2<-sum(s0[c(1:9)]*n0[c(10:18)])
      as3<-sum(swx24[c(1:9)])
      sigma<-(ab2+ab4*ab4*as3+as1*m_nf)/(ab3+ab4*nn[4]+as2)
      aaa1<-abs(sigma-aaa0)
      aaa0<-sigma
      if (n_iter>20) break
    }
    
    if (sigma<sigma0) {sigma<-sigma0}
    sigma4[1]<-sigma+sigma40;sigma4[c(2:9)]<-sigma4[1]
    sigma5[1]<-sigma/m_nf+sigma50;sigma5[c(3,7,9)]<-sigma5[1]
    sigma5[2]<-sigma5[1]+g_aa1;sigma5[4]<-sigma5[1]+g_aa2
    sigma5[5]<-sigma5[1]+g_aa3;sigma5[6]<-sigma5[1]+g_aa2
    sigma5[8]<-sigma5[1]+g_aa1
    ##############################criteria for iterations to stop######################## 
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
      sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*8
  #######first order parameters############
  mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m4[4],m4[5],m4[6],m4[7],m4[8],m4[9],m5[1],m5[2],m5[3],m5[4],m5[5],m5[6],m5[7],m5[8],m5[9]))
  B201111<-solve(crossprod(hh20,hh20))%*%crossprod(hh20,mm)
  ########second order genetic parameters###############   
  F2_jj <- sigmaF2 - sigma4[1]
  F2_gg <- sigma4[1]-sigma
  if(F2_jj<0) {F2_jj<-0}
  if(F2_gg<0 || F2_gg>sigmaF2) {F2_gg<-0}
  F2_ll <- F2_jj/sigmaF2
  F2_rr <- F2_gg/sigmaF2
  F3_jj <- sigmaF3 - sigma5[1]
  F3_gg <- sigma5[1]-sigma/m_nf
  if(F3_jj<0) {F3_jj<-0}
  if(F3_gg<0 || F3_gg>sigmaF3) {F3_gg<-0}
  F3_ll <- F3_jj/sigmaF3
  F3_rr <- F3_gg/sigmaF3
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma))
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
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma))
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
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma))
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
    F2gg <- (dataF2 - m4[i])/sqrt(sigma4[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mi4[i]
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
  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3); 
  F3w1<-1/(12*m_sam[5])
  F3bmw <- matrix(0,m_sam[5],1); F3bmwsl <- matrix(0,m_sam[5],d3)
  for(i in 1:d3){
    F3gg <- (dataF3 - m5[i])/sqrt(sigma5[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mi5[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[5]){F3P2 <- F3P2+runif(m_sam[5])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[5]) + sum((F3P2 - (as.matrix(c(1:m_sam[5])) - 0.5)/m_sam[5])^2)
  F3u <- as.matrix(c(12*m_sam[5]*((F3dd[1]/m_sam[5]-0.5)^2),((45*m_sam[5])/4)*((F3dd[2]/m_sam[5]-1/3)^2),180*m_sam[5]*((F3dd[3]/m_sam[5]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[5]))))
  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("MX2-A-AD",round(abc,4),round(AIC,4),round(m[1],4),round(m[2],4),round(m[3],4), round(t(m4),4),round(sigma4[1],4),
                       round(t(mi4),4),round(t(m5),4),round(t(sigma5),4),round(t(mi5),4),round(sigma,4),          
                       round(B201111[1],4)," "," "," "," ",round(B201111[2],4),round(B201111[3],4)," "," "," "," "," "," ",round(B201111[4],4),round(B201111[5],4),round(F2_jj,4),round(F2_ll*100,4),round(F2_gg,4),round(F2_rr*100,4),round(F3_jj,4),round(F3_ll*100,4),round(F3_gg,4),round(F3_rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
 output<-as.matrix(output)
  OUTPUT<-list(output,mix_pi4,mix_pi5)
  return(OUTPUT)
}
########################MX2-EA-AD(E-4)##################################
G5ModelFun[[22]] <- function(K1,logL,df11,df21,df31,df41,df51,G5text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataF2 <- as.matrix(as.numeric(df41[,1]))
  dataF3 <- as.matrix(as.numeric(df51[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2),sum(dataF3)))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2),mean(dataF3)))
  sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1);sigmaF2<-var(dataF2);sigmaF3<-var(dataF3)
  sumx6<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  sigma<-sumx6;sigma0<-sigma;sigma40<-sigmaF2;sigma400<-sigma40;sigma50<-sigmaF3;sigma500<-sigma50
  m_nf <- as.numeric(G5text2)
  ###############procedure start###########################
  mix_pi4<-as.matrix(c(0.0625,0.250,0.375,0.25,0.0625))
  mix_pi5<-as.matrix(c(0.0625,0.250,0.125,0.25,0.25,0.0625))
  a1<-sqrt(sigma400/m_sam[4])
  if (m[1]<m[3]) a1<--a1
  m4<-as.matrix(c(m[4]+2*a1,m[4]+a1,m[2],m[4]-a1,m[4]-2*a1))
  a1<-sqrt(sigma500/m_sam[5])
  if (m[1]<m[3]) {a1<--a1}
  m5<-as.matrix(c(m[5]+2*a1,m[5]+a1,m[2],m[2],m[5]-a1,m[3]-2*a1))
  nn<-m_sam
  abb<-sigma400/(sigma*5)
  sigma4<-matrix(0,5,1);sigma5<-matrix(0,6,1)
  sigma4[1]<-sigma400/abb;if (sigma4[1]>sigma400) {sigma4[1]<-sigma400/2}
  sigma5[1]<-sigma500/abb;if (sigma5[1]>sigma500) {sigma5[1]<-sigma500/2}
  #######first order parameters############
  hh21<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,-2,2,1,0,-1,-2,2,1,0,-1,-2,
                 1,0,-1,0,0,0,0,0,0,0,0,0,0,0,1,0,0.5,0.5,0.5,0.5,0.5,0.25,0.25,0.25,0.25,0.25),13,4)
  mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m4[4],m4[5],m5[1],m5[2],m5[3],m5[5],m5[6]))
  B21<-solve(crossprod(hh21,hh21))%*%crossprod(hh21,mm)
  g_aa1<-0.5*B21[2]^2/m_nf        #   0.5*d**2.
  g_aa2<-2*g_aa1
  sigma4[c(2:5)]<-sigma4[1];sigma5[2]<-sigma5[1]+g_aa1
  sigma5[3]<-sigma5[1];sigma5[4]<-sigma5[1]+g_aa2
  sigma5[5]<-sigma5[2];sigma5[6]<-sigma5[1]
  mi4<-mix_pi4[c(1:5)];mi5<-mix_pi5[c(1:6)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
    sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
  ########################################################
  ##########iteration process###########
  ############E-step###################
  ############F2 process###############
  d5<-5;d6<-6
  iteration <- 0; stopa <- 1000
  W4 <- matrix(0,5,m_sam[4]); swx24 <- matrix(0,5,1)
  W5 <- matrix(0,6,m_sam[5]); swx25 <- matrix(0,6,1)
  hh<-matrix(0,9,9);b_line<-matrix(0,9,1)
  n0<-matrix(0,11,1);s0<-matrix(0,6,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:d5) { W4[i,] <- mix_pi4[i]*dnorm(dataF2,m4[i],sqrt(sigma4[i]))/dmixnorm(dataF2,m4,sqrt(sigma4),mix_pi4) }
    mi4 <- as.matrix(rowSums(W4)/m_sam[4]);n0[c(1:5)]<-as.matrix(rowSums(W4));sumwx4 <- W4%*%dataF2
    #############F3 process############
    for(i in 1:d6) { W5[i,] <- mix_pi5[i]*dnorm(dataF3,m5[i],sqrt(sigma5[i]))/dmixnorm(dataF3,m5,sqrt(sigma5),mix_pi5) }
    mi5 <- as.matrix(rowSums(W5)/m_sam[5]);n0[c(6:11)]<-as.matrix(rowSums(W5))
    n0[c(1:11)][abs(n0[c(1:11)])<0.000001]<-0.000001;sumwx5 <- W5%*%dataF3
    aaa0<-0
    aaa1<-1000;n_iter<-0;AA<-matrix(0,9,1)
    while (aaa1>0.001)
    {
      n_iter<-n_iter+1
      ##########first order parameters###################
      mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m4[4],m4[5],m5[1],m5[2],m5[3],m5[5],m5[6]))
      B211<-solve(crossprod(hh21,hh21))%*%crossprod(hh21,mm)
      g_aa1<-0.5*B211[2]^2/m_nf        #   0.5*d**2.
      g_aa2<-2*g_aa1
      sigma4[c(2:5)]<-sigma4[1];sigma5[2]<-sigma5[1]+g_aa1
      sigma5[3]<-sigma5[1];sigma5[4]<-sigma5[1]+g_aa2
      sigma5[5]<-sigma5[2];sigma5[6]<-sigma5[1]
      #################restrictions########################################
      aa1<-sumwx5[3]*sigma5[4]+sigma5[3]*sumwx5[4];aa2<-n0[8]*sigma5[4]+sigma5[3]*n0[9]
      aa3<-aa1/aa2;aa4<-sigma5[3]*sigma5[4]/aa2
      hh[1,1]<-sigma*(1/nn[1]+4/nn[2]+1/nn[3])+16*sigma4[3]/n0[3]
      hh[1,2]<-sigma*2/nn[2]+12*sigma4[3]/n0[3]
      hh[1,3]<--4*sigma4[3]/n0[3]
      hh[1,4]<--4*sigma4[3]/n0[3]
      hh[1,5]<-0;hh[1,6]<-0;hh[1,7]<-0;hh[1,8]<-0
      hh[1,9]<-8*sigma4[3]/n0[3]
      hh[2,2]<-sigma/nn[2]+9*sigma4[3]/n0[3]+4*aa4
      hh[2,3]<--3*sigma4[3]/n0[3]
      hh[2,4]<--3*sigma4[3]/n0[3]
      hh[2,5]<-2*aa4
      hh[2,6]<-2*aa4
      hh[2,7]<-0;hh[2,8]<-0
      hh[2,9]<-6*sigma4[3]/n0[3]
      hh[3,3]<-sigma4[1]/n0[1]+4*sigma4[2]/n0[2]+sigma4[3]/n0[3]
      hh[3,4]<-sigma4[3]/n0[3]
      hh[3,5]<-0;hh[3,6]<-0
      hh[3,7]<-sigma4[1]/n0[1]+2*sigma4[2]/n0[2]
      hh[3,8]<-0
      hh[3,9]<-sigma4[1]/n0[1]-2*sigma4[3]/n0[3]
      hh[4,4]<-sigma4[3]/n0[3]+4*sigma4[4]/n0[4]+sigma4[5]/n0[5]
      hh[4,5]<-0;hh[4,6]<-0;hh[4,7]<-0
      hh[4,8]<--2*sigma4[4]/n0[4]-sigma4[5]/n0[5]
      hh[4,9]<--2*sigma4[3]/n0[3]+sigma4[5]/n0[5]
      hh[5,5]<-sigma5[1]/n0[6]+4*sigma5[2]/n0[7]+aa4
      hh[5,6]<-aa4
      hh[5,7]<--sigma5[1]/n0[6]-2*sigma5[2]/n0[7]
      hh[5,8]<-0;hh[5,9]<-0
      hh[6,6]<-aa4+4*sigma5[5]/n0[10]+sigma5[6]/n0[11];
      hh[6,7]<-0
      hh[6,8]<-2*sigma5[5]/n0[10]+sigma5[6]/n0[11]
      hh[6,9]<-0
      hh[7,7]<-sigma4[1]/n0[1]+sigma4[2]/n0[2]+sigma5[1]/n0[6]+sigma5[2]/n0[7]
      hh[7,8]<-0
      hh[7,9]<-sigma4[1]/n0[1]
      hh[8,8]<-sigma4[4]/n0[4]+sigma4[5]/n0[5]+sigma5[5]/n0[10]+sigma5[6]/n0[11]
      hh[8,9]<--sigma4[5]/n0[5]
      hh[9,9]<-sigma4[1]/n0[1]+4*sigma4[3]/n0[3]+sigma4[5]/n0[5]
      for(i in 2:9)
      {
        for(j in 1:(i-1))
        {
          hh[i,j]<-hh[j,i]
        }
      }
      ########################################################
      b_line[1]<-sumx[1]/nn[1]+2*sumx[2]/nn[2]+sumx[3]/nn[3]-4*sumwx4[3]/n0[3]
      b_line[2]<-sumx[2]/nn[2]-3*sumwx4[3]/n0[3]+2*aa3
      b_line[3]<-sumwx4[1]/n0[1]-2*sumwx4[2]/n0[2]+sumwx4[3]/n0[3]
      b_line[4]<-sumwx4[3]/n0[3]-2*sumwx4[4]/n0[4]+sumwx4[5]/n0[5]
      b_line[5]<-sumwx5[1]/n0[6]-2*sumwx5[2]/n0[7]+aa3
      b_line[6]<-aa3-2*sumwx5[5]/n0[10]+sumwx5[6]/n0[11]
      b_line[7]<-sumwx4[1]/n0[1]-sumwx4[2]/n0[2]-sumwx5[1]/n0[6]+sumwx5[2]/n0[7]
      b_line[8]<-sumwx4[4]/n0[4]-sumwx4[5]/n0[5]-sumwx5[5]/n0[10]+sumwx5[6]/n0[11]
      b_line[9]<-sumwx4[1]/n0[1]-2*sumwx4[3]/n0[3]+sumwx4[5]/n0[5]
      B2111<-solve(hh,b_line)
      ##########################################################
      m[1]<-(sumx[1]-sigma*B2111[1])/nn[1]
      m[2]<-(sumx[2]-sigma*(2*B2111[1]+B2111[2]))/nn[2]
      m[3]<-(sumx[3]-sigma*B2111[1])/nn[3]
      m4[1]<-(sumwx4[1]-sigma4[1]*(B2111[3]+B2111[7]+B2111[9]))/n0[1]
      m4[2]<-(sumwx4[2]+sigma4[1]*(2*B2111[3]+B2111[7]))/n0[2]
      m4[3]<-(sumwx4[3]+sigma4[1]*(4*B2111[1]+3*B2111[2]-B2111[3]-B2111[4]+2*B2111[9]))/n0[3]
      m4[4]<-(sumwx4[4]+sigma4[1]*(2*B2111[4]-B2111[8]))/n0[4]
      m4[5]<-(sumwx4[5]+sigma4[1]*(-B2111[4]+B2111[8]-B2111[9]))/n0[5]
      m5[1]<-(sumwx5[1]+sigma5[1]*(-B2111[5]+B2111[7]))/n0[6]
      m5[2]<-(sumwx5[2]+sigma5[2]*(2*B2111[5]-B2111[7]))/n0[7]
      ak10<-sumwx5[3]*sigma5[4]+sigma5[3]*sumwx5[4]
      ak11<-n0[8]*sigma5[4]+sigma5[3]*n0[9]
      m5[3]<-(ak10-sigma5[3]*sigma5[4]*(2*B2111[2]+B2111[5]+B2111[6]))/ak11
      m5[4]<-m5[3]
      m5[5]<-(sumwx5[5]+sigma5[5]*(2*B2111[6]+B2111[8]))/n0[10]
      m5[6]<-(sumwx5[6]-sigma5[6]*(B2111[6]+B2111[8]))/n0[11]
      aaa1<-max(abs(B2111-AA))
      AA<-B2111
      if (n_iter>20) break
    } 
    ##########obtain variance###############
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d5) {swx24[i] <- W4[i,]%*%(dataF2-m4[i])^2 } ;for(i in 1:d6) {swx25[i] <- W5[i,]%*%(dataF3-m5[i])^2 } 
    ab1<-sum(swx24);sigma4[1]<-ab1/nn[4];sigma40<-sigma4[1]-sigma
    if (sigma40<0) {sigma40<-0;sigma4[1]<-sigma}
    sigma4[c(2:5)]<-sigma4[1];aaa0<-sigma5[1]
    ab2<-swx25[1]+swx25[3]+swx25[6];ab3<-n0[6]+n0[8]+n0[11]
    #######first order parameters############
    mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m4[4],m4[5],m5[1],m5[2],m5[3],m5[5],m5[6]))
    B21111<-solve(crossprod(hh21,hh21))%*%crossprod(hh21,mm)
    n_iter<-0;aaa1<-1000
    aa12<-B21111[2];aa12<-0.5*aa12^2/m_nf
    while (aaa1>0.001)
    {
      n_iter<-n_iter+1
      aa1<-sigma5[1]/(sigma5[1]+aa12)
      aa2<-sigma5[1]/(sigma5[1]+2*aa12)
      as3<-ab2+aa1^2*(swx25[2]+swx25[5])+aa2^2*swx25[4]
      as4<-ab3+aa1*(n0[7]+n0[10])+aa2*n0[9]
      sigma5[1]<-as3/as4
      aaa1<-abs(sigma5[1]-aaa0)
      aaa0<-sigma5[1]
      if (n_iter>20) break
    } 
    sigma50<-sigma5[1]-sigma/m_nf;
    if (sigma50<0) {sigma50<-0;sigma5[1]<-sigma/m_nf}
    sigma5[2]<-sigma5[1]+aa12;sigma5[3]<-sigma5[1]
    sigma5[4]<-sigma5[1]+2*aa12;sigma5[5]<-sigma5[1]+aa12
    sigma5[6]<-sigma5[1]
    ab2<-ss1+ss2+ss3;ab3<-nn[1]+nn[2]+nn[3]
    n_iter<-0;aaa0<-sigma;an1<-sigma40;an2<-sigma50
    aaa1<-1000
    while (aaa1>0.001)
    {
      n_iter<-n_iter+1
      if (an1<0) {an1<-0};if (an2<0) {an2<-0}
      ab4<-sigma/(an1+sigma)
      s0[1]<-(sigma/m_nf)/(sigma/m_nf+an2);s0[2]<-(sigma/m_nf)/(sigma/m_nf+an2+aa12)
      s0[3]<-(sigma/m_nf)/(sigma/m_nf+an2);s0[4]<-(sigma/m_nf)/(sigma/m_nf+an2+2*aa12)
      s0[5]<-(sigma/m_nf)/(sigma/m_nf+an2+aa12);s0[6]<-(sigma/m_nf)/(sigma/m_nf+an2)
      as1<-sum(s0[c(1:6)]^2*swx25[c(1:6)])
      as2<-sum(s0[c(1:6)]*n0[c(6:11)])
      as3<-sum(swx24)
      sigma<-(ab2+ab4*ab4*as3+as1*m_nf)/(ab3+ab4*nn[4]+as2)
      aaa1<-abs(sigma-aaa0)
      aaa0<-sigma
      if (n_iter>20) break
    } 
    if (sigma<sigma0) {sigma<-sigma0}
    sigma4[1]<-sigma+sigma40;sigma4[c(2:5)]<-sigma4[1]
    sigma5[1]<-sigma/m_nf+sigma50;sigma5[2]<-sigma5[1]+aa12
    sigma5[3]<-sigma5[1];sigma5[4]<-sigma5[1]+2*aa12
    sigma5[5]<-sigma5[1]+aa12;sigma5[6]<-sigma5[1]
    ##############################criteria for iterations to stop######################## 
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
      sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*7
  #######first order parameters############
  mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m4[4],m4[5], m5[1],m5[2],m5[3],m5[5],m5[6]))
  B211111<-solve(crossprod(hh21,hh21))%*%crossprod(hh21,mm)
  ########second order genetic parameters###############   
  F2_jj <- sigmaF2 - sigma4[1]
  F2_gg <- sigma4[1]-sigma
  if(F2_jj<0) {F2_jj<-0}
  if(F2_gg<0 || F2_gg>sigmaF2) {F2_gg<-0}
  F2_ll <- F2_jj/sigmaF2
  F2_rr <- F2_gg/sigmaF2
  F3_jj <- sigmaF3 - sigma5[1]
  F3_gg <- sigma5[1]-sigma/m_nf
  if(F3_jj<0) {F3_jj<-0}
  if(F3_gg<0 || F3_gg>sigmaF3) {F3_gg<-0}
  F3_ll <- F3_jj/sigmaF3
  F3_rr <- F3_gg/sigmaF3
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma))
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
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma))
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
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma))
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
    F2gg <- (dataF2 - m4[i])/sqrt(sigma4[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mi4[i]
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
  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3); 
  F3w1<-1/(12*m_sam[5])
  F3bmw <- matrix(0,m_sam[5],1); F3bmwsl <- matrix(0,m_sam[5],d6)
  for(i in 1:d6){
    F3gg <- (dataF3 - m5[i])/sqrt(sigma5[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mi5[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[5]){F3P2 <- F3P2+runif(m_sam[5])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[5]) + sum((F3P2 - (as.matrix(c(1:m_sam[5])) - 0.5)/m_sam[5])^2)
  F3u <- as.matrix(c(12*m_sam[5]*((F3dd[1]/m_sam[5]-0.5)^2),((45*m_sam[5])/4)*((F3dd[2]/m_sam[5]-1/3)^2),180*m_sam[5]*((F3dd[3]/m_sam[5]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[5]))))
  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("MX2-EA-AD",round(abc,4),round(AIC,4),round(m[1],4),round(m[2],4),round(m[3],4), round(t(m4),4)," "," "," "," ",round(sigma4[1],4),
                       round(t(mi4),4)," "," "," "," ",round(t(m5),4)," "," "," ",round(t(sigma5),4)," "," "," ",round(t(mi5),4)," "," "," ",round(sigma,4),          
                       round(B211111[1],4)," "," "," "," ",round(B211111[2],4)," "," "," "," "," "," "," ",round(B211111[3],4),round(B211111[4],4),round(F2_jj,4),round(F2_ll*100,4),round(F2_gg,4),round(F2_rr*100,4),round(F3_jj,4),round(F3_ll*100,4),round(F3_gg,4),round(F3_rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
 output<-as.matrix(output)
  OUTPUT<-list(output,mix_pi4,mix_pi5)
  return(OUTPUT)
}
#########################MX2-CD-AD(E-5)###################################
G5ModelFun[[23]] <- function(K1,logL,df11,df21,df31,df41,df51,G5text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataF2 <- as.matrix(as.numeric(df41[,1]))
  dataF3 <- as.matrix(as.numeric(df51[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2),sum(dataF3)))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2),mean(dataF3)))
  sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1);sigmaF2<-var(dataF2);sigmaF3<-var(dataF3)
  sumx6<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  sigma<-sumx6;sigma0<-sigma;sigma40<-sigmaF2;sigma400<-sigma40;sigma50<-sigmaF3;sigma500<-sigma50
  m_nf <- as.numeric(G5text2)
  ###############procedure start###########################
  mix_pi4<-as.matrix(c(0.5625,0.1875,0.1875,0.0625))
  mix_pi5<-as.matrix(c(0.0625,0.125,0.0625,0.125,0.25,0.125,0.0625,0.125,0.0625))
  a1<-sqrt(sigma400/m_sam[4])
  if (m[1]<m[3]) a1<--a1
  m4<-as.matrix(c(m[1],m[4]+1.5*a1,m[4]-1.5*a1,m[3]))
  a1<-sqrt(sigma500/m_sam[5])
  if (m[1]<m[3]) {a1<--a1}
  m5<-as.matrix(c(m[5]+2.4*a1,m[5]+1.8*a1,m[5]+1.2*a1,m[5]+0.6*a1,m[5],m[5]-0.6*a1,m[5]-1.2*a1,m[5]-1.8*a1,m[5]-2.4*a1))
  nn<-m_sam
  abb<-sigma400/(sigma*5)
  sigma4<-matrix(0,4,1);sigma5<-matrix(0,9,1)
  sigma4[1]<-sigma400/abb;if (sigma4[1]>sigma400) {sigma4[1]<-sigma400/2}
  sigma5[1]<-sigma500/abb;if (sigma5[1]>sigma500) {sigma5[1]<-sigma500/2}
  #######first order parameters############
  hh22<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,-1,1,1,-1,-1,1,1,1,0.5,0.5,0.5,-1,-1,-1,
                 1,1,-1,1,-1,1,-1,1,0.5,-1,1,0.5,-1,1,0.5,-1,1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,
                 0,1,0,0.5,0.5,0.5,0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25),16,5)
  mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m4[4],m5[1],m5[2],m5[3],m5[4],m5[5],m5[6],m5[7],m5[8],m5[9]))
  B22<-solve(crossprod(hh22,hh22))%*%crossprod(hh22,mm)
  g_aa1<-0.5*B22[3]^2/m_nf        #   0.5*db**2.
  g_aa2<-0.5*B22[2]^2/m_nf        #   0.5*da**2.
  g_aa3<-g_aa1+g_aa2
  sigma4[c(2:4)]<-sigma4[1];sigma5[c(3,7,9)]<-sigma5[1]
  sigma5[2]<-sigma5[1]+g_aa1;sigma5[4]<-sigma5[1]+g_aa2
  sigma5[5]<-sigma5[1]+g_aa3;sigma5[6]<-sigma5[4];sigma5[8]<-sigma5[2]
  mi4<-mix_pi4[c(1:4)];mi5<-mix_pi5[c(1:9)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
    sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
  ########################################################
  ##########iteration process###########
  ############E-step###################
  ############F2 process###############
  d4<-4;d3<-9
  iteration <- 0; stopa <- 1000
  W4 <- matrix(0,4,m_sam[4]); swx24 <- matrix(0,4,1)
  W5 <- matrix(0,9,m_sam[5]); swx25 <- matrix(0,9,1)
  hh<-matrix(0,11,11);b_line<-matrix(0,11,1)
  n0<-matrix(0,13,1);s0<-matrix(0,9,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:d4) { W4[i,] <- mix_pi4[i]*dnorm(dataF2,m4[i],sqrt(sigma4[i]))/dmixnorm(dataF2,m4,sqrt(sigma4),mix_pi4) }
    mi4 <- as.matrix(rowSums(W4)/m_sam[4]);n0[c(1:4)]<-as.matrix(rowSums(W4));sumwx4 <- W4%*%dataF2
    #############F3 process############
    for(i in 1:d3) { W5[i,] <- mix_pi5[i]*dnorm(dataF3,m5[i],sqrt(sigma5[i]))/dmixnorm(dataF3,m5,sqrt(sigma5),mix_pi5) }
    mi5 <- as.matrix(rowSums(W5)/m_sam[5]);n0[c(5:13)]<-as.matrix(rowSums(W5))
    n0[c(1:13)][abs(n0[c(1:13)])<0.000001]<-0.000001;sumwx5 <- W5%*%dataF3
    aaa0<-0
    aaa1<-1000;n_iter<-0;AA<-matrix(0,11,1)
    while (aaa1>0.001)
    {
      n_iter<-n_iter+1
      ##########first order parameters###################
      mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m4[4],m5[1],m5[2],m5[3],m5[4],m5[5],m5[6],m5[7],m5[8],m5[9]))
      B221<-solve(crossprod(hh22,hh22))%*%crossprod(hh22,mm)
      g_aa1<-0.5*B221[3]^2/m_nf        #   0.5*d**2.
      g_aa2<-0.5*B221[2]^2/m_nf
      g_aa3<-g_aa1+g_aa2
      sigma4[c(2:4)]<-sigma4[1];sigma5[c(3,7,9)]<-sigma5[1];sigma5[2]<-sigma5[1]+g_aa1
      sigma5[4]<-sigma5[1]+g_aa2;sigma5[5]<-sigma5[1]+g_aa3
      sigma5[6]<-sigma5[4];sigma5[8]<-sigma5[2]
      #################restrictions########################################
      hh[1,1]<-sigma*(1/nn[1]+4/nn[2]+1/nn[3])+9*sigma4[1]/n0[1]+sigma4[4]/n0[4]
      hh[1,2]<-sigma*2/nn[2]+6*sigma4[1]/n0[1]+sigma4[4]/n0[4]
      hh[1,3]<--3*sigma4[1]/n0[1]-sigma4[4]/n0[4]
      hh[1,4]<-sigma/nn[1]+2*sigma/nn[2]+sigma/nn[3]+3*sigma4[1]/n0[1]
      hh[1,5]<-0;hh[1,6]<-0;hh[1,7]<-0
      hh[1,8]<-3*sigma/nn[1]+4*sigma/nn[2]+3*sigma/nn[3]
      hh[1,9]<--3*sigma4[1]/n0[1]+sigma4[4]/n0[4]
      hh[1,10]<-0
      hh[1,11]<-0
      hh[2,2]<-sigma/nn[2]+4*sigma4[1]/n0[1]+sigma4[4]/n0[4]+sigma5[1]/n0[5]+sigma5[9]/n0[13]
      hh[2,3]<--2*sigma4[1]/n0[1]-sigma4[4]/n0[4]
      hh[2,4]<-sigma/nn[2]+2*sigma4[1]/n0[1]
      hh[2,5]<-0
      hh[2,6]<--sigma5[9]/n0[13]
      hh[2,7]<-sigma5[1]/n0[5]
      hh[2,8]<-2*sigma/nn[2]
      hh[2,9]<--2*sigma4[1]/n0[1]+sigma4[4]/n0[4]-sigma5[1]/n0[5]+sigma5[9]/n0[13]
      hh[2,10]<-0
      hh[2,11]<-0
      hh[3,3]<-sigma4[1]/n0[1]+sigma4[2]/n0[2]+sigma4[3]/n0[3]+sigma4[4]/n0[4]
      hh[3,4]<--sigma4[1]/n0[1]-sigma4[2]/n0[2]+sigma4[3]/n0[3]
      hh[3,5]<--sigma4[2]/n0[2]+sigma4[3]/n0[3]
      hh[3,6]<-0;hh[3,7]<-0;hh[3,8]<-0
      hh[3,9]<-sigma4[1]/n0[1]-sigma4[4]/n0[4]
      hh[3,10]<-0
      hh[3,11]<-0
      hh[4,4]<-sigma/nn[1]+sigma/nn[2]+sigma/nn[3]+sigma4[1]/n0[1]+sigma4[2]/n0[2]+sigma4[3]/n0[3]+4*sigma5[3]/n0[7]
      hh[4,5]<-sigma4[2]/n0[2]+sigma4[3]/n0[3]+2*sigma5[3]/n0[7]
      hh[4,6]<-0
      hh[4,7]<-0
      hh[4,8]<-3*sigma/nn[1]+2*sigma/nn[2]+3*sigma/nn[3]
      hh[4,9]<--sigma4[1]/n0[1]
      hh[4,10]<-0
      hh[4,11]<-0
      hh[5,5]<-sigma4[2]/n0[2]+sigma4[3]/n0[3]+sigma5[3]/n0[7]+sigma5[7]/n0[11]
      hh[5,6]<--sigma5[7]/n0[11]
      hh[5,7]<-sigma5[7]/n0[11]
      hh[5,8]<-0
      hh[5,9]<-0
      hh[5,10]<-0
      hh[5,11]<-0
      hh[6,6]<-sigma5[4]/n0[8]+4*sigma5[5]/n0[9]+sigma5[6]/n0[10]+sigma5[7]/n0[11]+4*sigma5[8]/n0[12]+sigma5[9]/n0[13]
      hh[6,7]<--2*sigma5[4]/n0[8]-4*sigma5[5]/n0[9]-sigma5[7]/n0[11]-2*sigma5[8]/n0[12]
      hh[6,8]<--2*sigma5[4]/n0[8]-2*sigma5[6]/n0[10]-4*sigma5[8]/n0[12]
      hh[6,9]<--sigma5[9]/n0[13]
      hh[6,10]<-8*sigma5[5]/n0[9]+2*sigma5[8]/n0[12]
      hh[6,11]<-3*sigma5[4]/n0[8]+8*sigma5[5]/n0[9]+sigma5[6]/n0[10]
      hh[7,7]<-sigma5[1]/n0[5]+sigma5[2]/n0[6]+4*sigma5[4]/n0[8]+4*sigma5[5]/n0[9]+sigma5[7]/n0[11]+sigma5[8]/n0[12]
      hh[7,8]<-2*sigma5[2]/n0[6]+4*sigma5[4]/n0[8]+2*sigma5[8]/n0[12]
      hh[7,9]<--sigma5[1]/n0[5]
      hh[7,10]<--3*sigma5[2]/n0[6]-8*sigma5[5]/n0[9]-sigma5[8]/n0[12]
      hh[7,11]<--6*sigma5[4]/n0[8]-8*sigma5[5]/n0[9]
      hh[8,8]<-9*sigma/nn[1]+4*sigma/nn[2]+9*sigma/nn[3]+4*sigma5[2]/n0[6]+4*sigma5[4]/n0[8]+4*sigma5[6]/n0[10]+4*sigma5[8]/n0[12]
      hh[8,9]<-0
      hh[8,10]<--6*sigma5[2]/n0[6]-2*sigma5[8]/n0[12]
      hh[8,11]<--6*sigma5[4]/n0[8]-2*sigma5[6]/n0[10]
      hh[9,9]<-sigma4[1]/n0[1]+sigma4[4]/n0[4]+sigma5[1]/n0[5]+sigma5[9]/n0[13]
      hh[9,10]<-0
      hh[9,11]<-0
      hh[10,10]<-9*sigma5[2]/n0[6]+16*sigma5[5]/n0[9]+sigma5[8]/n0[12]
      hh[10,11]<-16*sigma5[5]/n0[9]
      hh[11,11]<-9*sigma5[4]/n0[8]+16*sigma5[5]/n0[9]+sigma5[6]/n0[10]
      for(i in 2:11)
      {
        for(j in 1:(i-1))
        {
          hh[i,j]<-hh[j,i]
        }
      }
      ##################################################################
      b_line[1]<-sumx[1]/nn[1]+2*sumx[2]/nn[2]+sumx[3]/nn[3]-3*sumwx4[1]/n0[1]-sumwx4[4]/n0[4]
      b_line[2]<-sumx[2]/nn[2]-2*sumwx4[1]/n0[1]-sumwx4[4]/n0[4]+sumwx5[1]/n0[5]+sumwx5[9]/n0[13]
      b_line[3]<-sumwx4[1]/n0[1]-sumwx4[2]/n0[2]-sumwx4[3]/n0[3]+sumwx4[4]/n0[4]
      b_line[4]<-sumx[1]/nn[1]+sumx[2]/nn[2]+sumx[3]/nn[3]-sumwx4[1]/n0[1]+sumwx4[2]/n0[2]-sumwx4[3]/n0[3]-2*sumwx5[3]/n0[7]
      b_line[5]<-sumwx4[2]/n0[2]-sumwx4[3]/n0[3]-sumwx5[3]/n0[7]+sumwx5[7]/n0[11]
      b_line[6]<-sumwx5[4]/n0[8]-2*sumwx5[5]/n0[9]+sumwx5[6]/n0[10]-sumwx5[7]/n0[11]+2*sumwx5[8]/n0[12]-sumwx5[9]/n0[13]
      b_line[7]<-sumwx5[1]/n0[5]-2*sumwx5[4]/n0[8]+sumwx5[7]/n0[11]-sumwx5[2]/n0[6]+2*sumwx5[5]/n0[9]-sumwx5[8]/n0[12]
      b_line[8]<-3*sumx[1]/nn[1]+2*sumx[2]/nn[2]+3*sumx[3]/nn[3]-2*sumwx5[2]/n0[6]-2*sumwx5[4]/n0[8]-2*sumwx5[6]/n0[10]-2*sumwx5[8]/n0[12]
      b_line[9]<-sumwx4[1]/n0[1]-sumwx4[4]/n0[4]-sumwx5[1]/n0[5]+sumwx5[9]/n0[13]
      b_line[10]<-3*sumwx5[2]/n0[6]-4*sumwx5[5]/n0[9]+sumwx5[8]/n0[12]
      b_line[11]<-3*sumwx5[4]/n0[8]-4*sumwx5[5]/n0[9]+sumwx5[6]/n0[10]
      B2211<-solve(hh,b_line)
      ##########################################################
      m[1]<-(sumx[1]-sigma*(B2211[1]+3*B2211[8]+B2211[4]))/nn[1]
      m[2]<-(sumx[2]-sigma*(2*B2211[1]+B2211[2]+2*B2211[8]+B2211[4]))/nn[2]
      m[3]<-(sumx[3]-sigma*(B2211[1]+3*B2211[8]+B2211[4]))/nn[3]
      m4[1]<-(sumwx4[1]+sigma4[1]*(3*B2211[1]+2*B2211[2]-B2211[3]+B2211[4]-B2211[9]))/n0[1]
      m4[2]<-(sumwx4[2]+sigma4[2]*(B2211[3]-B2211[4]-B2211[5]))/n0[2]
      m4[3]<-(sumwx4[3]+sigma4[3]*(B2211[3]+B2211[4]+B2211[5]))/n0[3]
      m4[4]<-(sumwx4[4]+sigma4[4]*(B2211[1]+B2211[2]-B2211[3]+B2211[9]))/n0[4]
      m5[1]<-(sumwx5[1]+sigma5[1]*(-B2211[2]-B2211[7]+B2211[9]))/n0[5]
      m5[2]<-(sumwx5[2]+sigma5[2]*(B2211[7]+2*B2211[8]-3*B2211[10]))/n0[6]
      m5[3]<-(sumwx5[3]+sigma5[3]*(2*B2211[4]+B2211[5]))/n0[7]
      m5[4]<-(sumwx5[4]+sigma5[4]*(-B2211[6]+2*B2211[7]+2*B2211[8]-3*B2211[11]))/n0[8]
      m5[5]<-(sumwx5[5]+sigma5[5]*(2*B2211[6]-2*B2211[7]+4*B2211[10]+4*B2211[11]))/n0[9]
      m5[6]<-(sumwx5[6]+sigma5[6]*(-B2211[6]+2*B2211[8]-B2211[11]))/n0[10]
      m5[7]<-(sumwx5[7]+sigma5[7]*(-B2211[5]+B2211[6]-B2211[7]))/n0[11]
      m5[8]<-(sumwx5[8]+sigma5[8]*(-2*B2211[6]+B2211[7]+2*B2211[8]-B2211[10]))/n0[12]
      m5[9]<-(sumwx5[9]+sigma5[9]*(-B2211[2]+B2211[6]-B2211[9]))/n0[13]
      aaa1<-max(abs(B2211-AA))
      AA<-B2211
      if (n_iter>20) break
    } 
    ##########obtain variance###############
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d4) {swx24[i] <- W4[i,]%*%(dataF2-m4[i])^2 } ;for(i in 1:d3) {swx25[i] <- W5[i,]%*%(dataF3-m5[i])^2 } 
    ab1<-sum(swx24);sigma4[1]<-ab1/nn[4];sigma40<-sigma4[1]-sigma
    if (sigma40<0) {sigma40<-0;sigma4[1]<-sigma};sigma4[c(2:4)]<-sigma4[1]
    aaa0<-sigma5[1]
    n_iter<-0
    ab2<-swx25[1]+swx25[3]+swx25[7]+swx25[9];ab3<-n0[5]+n0[7]+n0[11]+n0[13]
    #######first order parameters############
    mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m4[4],m5[1],m5[2],m5[3],m5[4],m5[5],m5[6],m5[7],m5[8],m5[9]))
    B22111<-solve(crossprod(hh22,hh22))%*%crossprod(hh22,mm)
    g_aa1<-0.5*B22111[3]^2/m_nf        #   0.5*db**2.
    g_aa2<-0.5*B22111[2]^2/m_nf        #   0.5*da**2.
    g_aa3<-g_aa1+g_aa2
    n_iter<-0;aaa1<-1000
    while (aaa1>0.001)
    {
      n_iter<-n_iter+1
      aa1<-sigma5[1]/(sigma5[1]+g_aa1)
      aa2<-sigma5[1]/(sigma5[1]+g_aa2)
      aa3<-sigma5[1]/(sigma5[1]+g_aa3)
      as3<-ab2+aa1^2*(swx25[2]+swx25[8])+aa2^2*(swx25[4]+swx25[6])+aa3^2*swx25[5]
      as4<-ab3+aa1*(n0[6]+n0[12])+aa2*(n0[8]+n0[10])+aa3*n0[9]
      sigma5[1]<-as3/as4
      aaa1<-abs(sigma5[1]-aaa0)
      aaa0<-sigma5[1]
      if (n_iter>20) break
    } 
    sigma50<-sigma5[1]-sigma/m_nf;
    if (sigma50<0) {sigma50<-0;sigma5[1]<-sigma/m_nf}
    sigma5[c(3,7,9)]<-sigma5[1]
    sigma5[2]<-sigma5[1]+g_aa1;sigma5[4]<-sigma5[1]+g_aa2
    sigma5[5]<-sigma5[1]+g_aa3;sigma5[6]<-sigma5[4];sigma5[8]<-sigma5[2]
    ab2<-ss1+ss2+ss3;ab3<-nn[1]+nn[2]+nn[3]
    n_iter<-0;aaa0<-sigma
    an1<-sigma40;an2<-sigma50
    while (aaa1>0.001)
    {
      n_iter<-n_iter+1
      if (an1<0) {an1<-0};if (an2<0) {an2<-0}
      ab4<-sigma/(an1+sigma)
      s0[1]<-(sigma/m_nf)/(sigma/m_nf+an2)
      s0[2]<-(sigma/m_nf)/(sigma/m_nf+an2+g_aa1)
      s0[3]<-(sigma/m_nf)/(sigma/m_nf+an2)
      s0[4]<-(sigma/m_nf)/(sigma/m_nf+an2+g_aa2)
      s0[5]<-(sigma/m_nf)/(sigma/m_nf+an2+g_aa3)
      s0[6]<-(sigma/m_nf)/(sigma/m_nf+an2+g_aa2)
      s0[7]<-(sigma/m_nf)/(sigma/m_nf+an2)
      s0[8]<-(sigma/m_nf)/(sigma/m_nf+an2+g_aa1)
      s0[9]<-(sigma/m_nf)/(sigma/m_nf+an2)
      as1<-sum(s0[c(1:9)]^2*swx25[c(1:9)])
      as2<-sum(s0[c(1:9)]*n0[c(5:13)])
      as3<-sum(swx24)
      sigma<-(ab2+ab4^2*as3+as1*m_nf)/(ab3+ab4*nn[4]+as2)
      aaa1<-abs(sigma-aaa0)
      aaa0<-sigma
      if (n_iter>20) break
    } 
    if (sigma<sigma0) {sigma<-sigma0}
    sigma4[1]<-sigma+sigma40;sigma4[c(2:4)]<-sigma4[1]
    sigma5[1]<-sigma/m_nf+sigma50;sigma5[c(3,7,9)]<-sigma5[1]
    sigma5[2]<-sigma5[1]+g_aa1;sigma5[4]<-sigma5[1]+g_aa2;sigma5[5]<-sigma5[1]+g_aa3
    sigma5[6]<-sigma5[4];sigma5[8]<-sigma5[2]
    ##############################criteria for iterations to stop######################## 
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
      sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*8
  #######first order parameters############
  mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m4[4],m5[1],m5[2],m5[3],m5[4],m5[5],m5[6],m5[7],m5[8],m5[9]))
  B221111<-solve(crossprod(hh22,hh22))%*%crossprod(hh22,mm)
  ########second order genetic parameters###############   
  F2_jj <- sigmaF2 - sigma4[1]
  F2_gg <- sigma4[1]-sigma
  if(F2_jj<0) {F2_jj<-0}
  if(F2_gg<0 || F2_gg>sigmaF2) {F2_gg<-0}
  F2_ll <- F2_jj/sigmaF2
  F2_rr <- F2_gg/sigmaF2
  F3_jj <- sigmaF3 - sigma5[1]
  F3_gg <- sigma5[1]-sigma/m_nf
  if(F3_jj<0) {F3_jj<-0}
  if(F3_gg<0 || F3_gg>sigmaF3) {F3_gg<-0}
  F3_ll <- F3_jj/sigmaF3
  F3_rr <- F3_gg/sigmaF3
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma))
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
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma))
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
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma))
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
    F2gg <- (dataF2 - m4[i])/sqrt(sigma4[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mi4[i]
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
  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3); 
  F3w1<-1/(12*m_sam[5])
  F3bmw <- matrix(0,m_sam[5],1); F3bmwsl <- matrix(0,m_sam[5],d3)
  for(i in 1:d3){
    F3gg <- (dataF3 - m5[i])/sqrt(sigma5[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mi5[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[5]){F3P2 <- F3P2+runif(m_sam[5])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[5]) + sum((F3P2 - (as.matrix(c(1:m_sam[5])) - 0.5)/m_sam[5])^2)
  F3u <- as.matrix(c(12*m_sam[5]*((F3dd[1]/m_sam[5]-0.5)^2),((45*m_sam[5])/4)*((F3dd[2]/m_sam[5]-1/3)^2),180*m_sam[5]*((F3dd[3]/m_sam[5]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[5]))))
  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("MX2-CD-AD",round(abc,4),round(AIC,4),round(m[1],4),round(m[2],4),round(m[3],4), round(t(m4),4)," "," "," "," "," ",round(sigma4[1],4),
                       round(t(mi4),4)," "," "," "," "," ",round(t(m5),4),round(t(sigma5),4),round(t(mi5),4),round(sigma,4),          
                       round(B221111[1],4)," "," "," "," ",round(B221111[2],4),round(B221111[3],4)," "," "," "," "," "," ",round(B221111[4],4),round(B221111[5],4),round(F2_jj,4),round(F2_ll*100,4),round(F2_gg,4),round(F2_rr*100,4),round(F3_jj,4),round(F3_ll*100,4),round(F3_gg,4),round(F3_rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
 output<-as.matrix(output)
  OUTPUT<-list(output,mix_pi4,mix_pi5)
  return(OUTPUT)
}

#####################MX2-EAD-AD(E-6)####################################
G5ModelFun[[24]] <- function(K1,logL,df11,df21,df31,df41,df51,G5text2){
  dataP1 <- as.matrix(as.numeric(df11[,1]));dataF1 <- as.matrix(as.numeric(df21[,1]))
  dataP2 <- as.matrix(as.numeric(df31[,1]));dataF2 <- as.matrix(as.numeric(df41[,1]))
  dataF3 <- as.matrix(as.numeric(df51[,1]))
  m_sam<-c(dim(dataP1)[1],dim(dataF1)[1],dim(dataP2)[1],dim(dataF2)[1],dim(dataF3)[1])
  m_esp <- 0.0001
  sumx<-as.matrix(c(sum(dataP1),sum(dataF1),sum(dataP2),sum(dataF2),sum(dataF3)))
  m<-as.matrix(c(mean(dataP1),mean(dataF1),mean(dataP2),mean(dataF2),mean(dataF3)))
  sigmaP1<-var(dataP1);sigmaP2<-var(dataP2);sigmaF1<-var(dataF1);sigmaF2<-var(dataF2);sigmaF3<-var(dataF3)
  sumx6<-((m_sam[1]-1)*sigmaP1+(m_sam[3]-1)*sigmaP2+(m_sam[2]-1)*sigmaF1)/((m_sam[1]+m_sam[3]+m_sam[2])-3)
  sigma<-sumx6;sigma0<-sigma;sigma40<-sigmaF2;sigma400<-sigma40;sigma50<-sigmaF3;sigma500<-sigma50
  m_nf <- as.numeric(G5text2)
  ###############procedure start###########################
  mix_pi4<-as.matrix(c(0.5625,0.375,0.0625))
  mix_pi5<-as.matrix(c(0.0625,0.25,0.25,0.125,0.25,0.0625))
  a1<-sqrt(sigma400/m_sam[4])
  if (m[1]<m[3]) a1<--a1
  m4<-as.matrix(c(m[1],m[4],m[3]))
  a1<-sqrt(sigma500/m_sam[5])
  if (m[1]<m[3]) {a1<--a1}
  m5<-as.matrix(c(m[5]+2.4*a1,m[5]+1.6*a1,m[5]+1.2*a1,m[5],m[5]-0.5*a1,m[5]-2*a1))
  nn<-m_sam
  abb<-sigma400/(sigma*5)
  sigma4<-matrix(0,3,1);sigma5<-matrix(0,6,1)
  sigma4[1]<-sigma400/abb;if (sigma4[1]>sigma400) {sigma4[1]<-sigma400/2}
  sigma5[1]<-sigma500/abb;if (sigma5[1]>sigma500) {sigma5[1]<-sigma500/2}
  #######first order parameters############
  hh23<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,2,2,-2,2,0,-2,2,1.5,1,0,-0.5,-0.5,
                 1,0,-1,0,0,0,0,0,0,0,0,0,0,1,0,0.5,0.5,0.5,0.25,0.25,0.25,0.25,0.25,0.25),12,4)
  mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m5[1],m5[2],m5[3],m5[4],m5[5],m5[6]))
  B23<-solve(crossprod(hh23,hh23))%*%crossprod(hh23,mm)
  g_aa1<-0.75*B23[2]^2/m_nf      #   0.75*d**2.
  g_aa2<-1.5*B23[2]^2/m_nf        #   1.5*d**2.
  sigma4[c(2,3)]<-sigma4[1]
  sigma5[2]<-sigma5[1]+g_aa1;sigma5[3]<-sigma5[1]+g_aa2
  sigma5[4]<-sigma5[1];sigma5[5]<-sigma5[2];sigma5[6]<-sigma5[1]
  mi4<-mix_pi4[c(1:3)];mi5<-mix_pi5[c(1:6)]
  L0<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
    sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
    sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
  ########################################################
  ##########iteration process###########
  ############E-step###################
  ############F2 process###############
  d1<-3;d6<-6
  iteration <- 0; stopa <- 1000
  W4 <- matrix(0,3,m_sam[4]); swx24 <- matrix(0,3,1)
  W5 <- matrix(0,6,m_sam[5]); swx25 <- matrix(0,6,1)
  hh<-matrix(0,8,8);b_line<-matrix(0,8,1)
  n0<-matrix(0,9,1);s0<-matrix(0,6,1)
  while(stopa > m_esp&&iteration<=1000){
    iteration <- iteration + 1
    for(i in 1:d1) { W4[i,] <- mix_pi4[i]*dnorm(dataF2,m4[i],sqrt(sigma4[i]))/dmixnorm(dataF2,m4,sqrt(sigma4),mix_pi4) }
    mi4 <- as.matrix(rowSums(W4)/m_sam[4]);n0[c(1:3)]<-as.matrix(rowSums(W4));sumwx4 <- W4%*%dataF2
    #############F3 process############
    for(i in 1:d6) { W5[i,] <- mix_pi5[i]*dnorm(dataF3,m5[i],sqrt(sigma5[i]))/dmixnorm(dataF3,m5,sqrt(sigma5),mix_pi5) }
    mi5 <- as.matrix(rowSums(W5)/m_sam[5]);n0[c(4:9)]<-as.matrix(rowSums(W5))
    n0[c(1:9)][abs(n0[c(1:9)])<0.000001]<-0.000001;sumwx5 <- W5%*%dataF3
    aaa0<-0
    aaa1<-1000;n_iter<-0;AA<-matrix(0,8,1)
    while (aaa1>0.01)
    {
      n_iter<-n_iter+1
      ##########first order parameters###################
      mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m5[1],m5[2],m5[3],m5[4],m5[5],m5[6]))
      B231<-solve(crossprod(hh23,hh23))%*%crossprod(hh23,mm)
      g_aa1<-0.75*B231[2]^2/m_nf        #   0.5*db**2.
      g_aa2<-1.5*B231[2]^2/m_nf        #   0.5*da**2.
      sigma4[c(2:3)]<-sigma4[1];sigma5[2]<-sigma5[1]+g_aa1
      sigma5[3]<-sigma5[1]+g_aa2;sigma5[4]<-sigma5[1]
      sigma5[5]<-sigma5[2];sigma5[6]<-sigma5[1]
      #################restrictions########################################
      hh[1,1]<-sigma*(1/nn[1]+4/nn[2]+1/nn[3])+9*sigma4[1]/n0[1]+sigma4[3]/n0[3]
      hh[1,2]<-sigma*2/nn[2]+6*sigma4[1]/n0[1]+sigma4[3]/n0[3]
      hh[1,3]<--3*sigma4[1]/n0[1]
      hh[1,4]<--3*sigma4[1]/n0[1]+sigma4[3]/n0[3]
      hh[1,5]<--9*sigma4[1]/n0[1]+sigma4[3]/n0[3]
      hh[1,6]<-0;hh[1,7]<-0;hh[1,8]<-0
      hh[2,2]<-sigma/nn[2]+4*sigma4[1]/n0[1]+sigma4[3]/n0[3]+sigma5[1]/n0[4]+ sigma5[6]/n0[9]
      hh[2,3]<--2*sigma4[1]/n0[1]-sigma5[1]/n0[4]
      hh[2,4]<--2*sigma4[1]/n0[1]+sigma4[3]/n0[3]-sigma5[1]/n0[4]+sigma5[6]/n0[9]
      hh[2,5]<--6*sigma4[1]/n0[1]+sigma4[3]/n0[3]
      hh[2,6]<-3*sigma5[1]/n0[4]+sigma5[6]/n0[9]
      hh[2,7]<-sigma5[1]/n0[4]+3*sigma5[6]/n0[9]
      hh[2,8]<-0
      hh[3,3]<-sigma4[1]/n0[1]+sigma4[2]/n0[2]+sigma5[1]/n0[4]+sigma5[4]/n0[7]
      hh[3,4]<-sigma4[1]/n0[1]+sigma5[1]/n0[4]
      hh[3,5]<-3*sigma4[1]/n0[1]+2*sigma4[2]/n0[2]
      hh[3,6]<--3*sigma5[1]/n0[4]
      hh[3,7]<--sigma5[1]/n0[4]-8*sigma5[4]/n0[7]
      hh[3,8]<--sigma5[4]/n0[7]
      hh[4,4]<-sigma4[1]/n0[1]+sigma4[3]/n0[3]+sigma5[1]/n0[4]+sigma5[6]/n0[9]
      hh[4,5]<-3*sigma4[1]/n0[1]+sigma4[3]/n0[3]
      hh[4,6]<--3*sigma5[1]/n0[4]+sigma5[6]/n0[9]
      hh[4,7]<--sigma5[1]/n0[4]+3*sigma5[6]/n0[9]
      hh[4,8]<-0
      hh[5,5]<-9*sigma4[1]/n0[1]+4*sigma4[2]/n0[2]+sigma4[3]/n0[3]+16*sigma5[2]/n0[5]+16*sigma5[5]/n0[8]
      hh[5,6]<-0;hh[5,7]<-0
      hh[5,8]<--4*sigma5[2]/n0[5]+4*sigma5[5]/n0[8]
      hh[6,6]<-9*sigma5[1]/n0[4]+16*sigma5[3]/n0[6]+sigma5[6]/n0[9]
      hh[6,7]<-3*sigma5[1]/n0[4]-16*sigma5[3]/n0[6]+3*sigma5[6]/n0[9]
      hh[6,8]<-4*sigma5[3]/n0[6]
      hh[7,7]<-sigma5[1]/n0[4]+16*sigma5[3]/n0[6]+64*sigma5[4]/n0[7]+9*sigma5[6]/n0[9]
      hh[7,8]<--4*sigma5[3]/n0[6]+8*sigma5[4]/n0[7]
      hh[8,8]<-sigma5[2]/n0[5]+sigma5[3]/n0[6]+sigma5[4]/n0[7]+sigma5[5]/n0[8]
      for(i in 2:8)
      {
        for(j in 1:(i-1))
        {
          hh[i,j]<-hh[j,i]
        }
      }
      ##################################################################### 
      b_line[1]<-sumx[1]/nn[1]+2*sumx[2]/nn[2]+sumx[3]/nn[3]-3*sumwx4[1]/n0[1]-sumwx4[3]/n0[3]
      b_line[2]<-sumx[2]/nn[2]-2*sumwx4[1]/n0[1]-sumwx4[3]/n0[3]+sumwx5[1]/n0[4]+sumwx5[6]/n0[9]
      b_line[3]<-sumwx4[1]/n0[1]-sumwx4[2]/n0[2]-sumwx5[1]/n0[4]+sumwx5[4]/n0[7]
      b_line[4]<-sumwx4[1]/n0[1]-sumwx4[3]/n0[3]-sumwx5[1]/n0[4]+sumwx5[6]/n0[9]
      b_line[5]<-3*sumwx4[1]/n0[1]-2*sumwx4[2]/n0[2]-sumwx4[3]/n0[3]-4*sumwx5[2]/n0[5]+4*sumwx5[5]/n0[8]
      b_line[6]<-3*sumwx5[1]/n0[4]-4*sumwx5[3]/n0[6]+sumwx5[6]/n0[9]
      b_line[7]<-sumwx5[1]/n0[4]+4*sumwx5[3]/n0[6]-8*sumwx5[4]/n0[7]+3*sumwx5[6]/n0[9]
      b_line[8]<-sumwx5[2]/n0[5]-sumwx5[3]/n0[6]-sumwx5[4]/n0[7]+sumwx5[5]/n0[8]
      B2311<-solve(hh,b_line)
      ##########################################################
      m[1]<-(sumx[1]-sigma*B2311[1])/nn[1]
      m[2]<-(sumx[2]-sigma*(2*B2311[1]+B2311[2]))/nn[2]
      m[3]<-(sumx[3]-sigma*B2311[1])/nn[3]
      m4[1]<-(sumwx4[1]+sigma4[1]*(3*B2311[1]+2*B2311[2]-B2311[3]-B2311[4]-3*B2311[5]))/n0[1]
      m4[2]<-(sumwx4[2]+sigma4[1]*(B2311[3]+2*B2311[5]))/n0[2]
      m4[3]<-(sumwx4[3]+sigma4[1]*(B2311[1]+B2311[2]+B2311[4]+B2311[5]))/n0[3]
      m5[1]<-(sumwx5[1]+sigma5[1]*(-B2311[2]+B2311[3]+B2311[4]-3*B2311[6]-B2311[7]))/n0[4]
      m5[2]<-(sumwx5[2]+sigma5[2]*(4*B2311[5]-B2311[8]))/n0[5]
      m5[3]<-(sumwx5[3]+sigma5[3]*(4*B2311[6]-4*B2311[7]+B2311[8]))/n0[6]
      m5[4]<-(sumwx5[4]+sigma5[4]*(-B2311[3]+8*B2311[7]+B2311[8]))/n0[7]
      m5[5]<-(sumwx5[5]-sigma5[5]*(4*B2311[5]+B2311[8]))/n0[8]
      m5[6]<-(sumwx5[6]-sigma5[6]*(B2311[2]+B2311[4]+B2311[6]+3*B2311[7]))/n0[9]
      aaa1<-max(abs(B2311-AA))
      AA<-B2311
      if (n_iter>20) break
    } 
    ##########obtain variance###############
    ss1<-sum((dataP1-m[1])^2);ss3<-sum((dataP2-m[3])^2);ss2<-sum((dataF1-m[2])^2)
    for(i in 1:d1) {swx24[i] <- W4[i,]%*%(dataF2-m4[i])^2 } ;for(i in 1:d6) {swx25[i] <- W5[i,]%*%(dataF3-m5[i])^2 } 
    ab1<-sum(swx24);sigma4[1]<-ab1/nn[4];sigma40<-sigma4[1]-sigma
    if (sigma40<0) {sigma40<-0;sigma4[1]<-sigma};sigma4[c(2:3)]<-sigma4[1];aaa0<-sigma5[1]
    n_iter<-0
    ab2<-swx25[1]+swx25[4]+swx25[6];ab3<-n0[4]+n0[7]+n0[9]
    mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m5[1],m5[2],m5[3],m5[4],m5[5],m5[6]))
    B23111<-solve(crossprod(hh23,hh23))%*%crossprod(hh23,mm)
    g_aa1<-0.75*B23111[2]^2/m_nf        #   0.75*d**2.
    g_aa2<-1.5*B23111[2]^2/m_nf        #   1.5*d**2.
    aaa1<-1000;n_iter<-0
    while (aaa1>0.001){
      n_iter<-n_iter+1
      aa1<-sigma5[1]/(sigma5[1]+g_aa1);
      aa2<-sigma5[1]/(sigma5[1]+g_aa2);
      as3<-ab2+aa1^2*(swx25[2]+swx25[5])+aa2^2*swx25[3];
      as4<-ab3+aa1*(n0[5]+n0[8])+aa2*n0[6];
      sigma5[1]<-as3/as4;
      aaa1<-abs(sigma5[1]-aaa0);
      aaa0<-sigma5[1];
      if (n_iter>20) break 
    }
    sigma50<-sigma5[1]-sigma/m_nf
    if (sigma50<0) {sigma50<-0;sigma5[1]<-sigma/m_nf}
    sigma5[2]<-sigma5[1]+g_aa1;sigma5[3]<-sigma5[1]+g_aa2
    sigma5[4]<-sigma5[1];sigma5[5]<-sigma5[1]+g_aa1
    sigma5[6]<-sigma5[1]
    ab2<-ss1+ss2+ss3;ab3<-nn[1]+nn[2]+nn[3]
    n_iter<-0;aaa0<-sigma;aaa1<-1000
    while (aaa1>0.0001){
      n_iter<-n_iter+1
      ab4<-sigma/(sigma40+sigma)
      s0[1]<-(sigma/m_nf)/(sigma/m_nf+sigma50)
      s0[2]<-(sigma/m_nf)/(sigma/m_nf+sigma50+g_aa1)
      s0[3]<-(sigma/m_nf)/(sigma/m_nf+sigma50+g_aa2)
      s0[4]<-(sigma/m_nf)/(sigma/m_nf+sigma50)
      s0[5]<-(sigma/m_nf)/(sigma/m_nf+sigma50+g_aa1)
      s0[6]<-(sigma/m_nf)/(sigma/m_nf+sigma50)
      as1<-sum(s0[c(1:6)]^2*swx25[c(1:6)])
      as2<-sum(s0[c(1:6)]*n0[c(4:9)])
      as3<-sum(swx24)
      sigma<-(ab2+ab4*ab4*as3+as1*m_nf)/(ab3+ab4*nn[4]+as2)
      aaa1<-abs(sigma-aaa0)
      aaa0<-sigma
      if (n_iter>20) break
    } 
    if (sigma<sigma0) {sigma<-sigma0}
    sigma4[1]<-sigma+sigma40;sigma4[c(2,3)]<-sigma4[1]
    sigma5[1]<-sigma/m_nf+sigma50;sigma5[2]<-sigma5[1]+g_aa1
    sigma5[3]<-sigma5[1]+g_aa2;sigma5[4]<-sigma5[1]
    sigma5[5]<-sigma5[1]+g_aa1;sigma5[6]<-sigma5[1]
    ##############################criteria for iterations to stop######################## 
    L1<-sum(log(dnorm(dataP1,m[1],sqrt(sigma))))+sum(log(dnorm(dataF1,m[2],sqrt(sigma))))+
      sum(log(dnorm(dataP2,m[3],sqrt(sigma))))+sum(log(dmixnorm(dataF2,m4,sqrt(sigma4),mi4)))+
      sum(log(dmixnorm(dataF3,m5,sqrt(sigma5),mi5)))
    stopa <- L1 - L0
    L0 <- L1
    if(stopa < 0) {stopa <- -stopa}
  }
  abc<-L0
  AIC<--2*abc+2*7
  #######first order parameters############
  mm<-as.matrix(c(m[1],m[2],m[3],m4[1],m4[2],m4[3],m5[1],m5[2],m5[3],m5[4],m5[5],m5[6]))
  B231111<-solve(crossprod(hh23,hh23))%*%crossprod(hh23,mm)
  ########second order genetic parameters###############   
  F2_jj <- sigmaF2 - sigma4[1]
  F2_gg <- sigma4[1]-sigma
  if(F2_jj<0) {F2_jj<-0}
  if(F2_gg<0 || F2_gg>sigmaF2) {F2_gg<-0}
  F2_ll <- F2_jj/sigmaF2
  F2_rr <- F2_gg/sigmaF2
  F3_jj <- sigmaF3 - sigma5[1]
  F3_gg <- sigma5[1]-sigma/m_nf
  if(F3_jj<0) {F3_jj<-0}
  if(F3_gg<0 || F3_gg>sigmaF3) {F3_gg<-0}
  F3_ll <- F3_jj/sigmaF3
  F3_rr <- F3_gg/sigmaF3
  #####################################hypothesis testing for P1########################################  
  dataP1<-sort(dataP1)
  P1w1<-1/(12*m_sam[1])
  P1bmw <- matrix(0,m_sam[1],1)
  P1gg <- (dataP1 - m[1])/sqrt(as.vector(sigma))
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
  F1gg <- (dataF1 - m[2])/sqrt(as.vector(sigma))
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
  P2gg <- (dataP2 - m[3])/sqrt(as.vector(sigma))
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
    F2gg <- (dataF2 - m4[i])/sqrt(sigma4[i])
    F2bmw[which(F2gg>=0)] <- pnorm(F2gg[F2gg>=0])
    F2bmw[which(F2gg<0)] <- 1 - pnorm(abs(F2gg[F2gg<0]))
    F2bmwsl[,i] <- F2bmw*mi4[i]
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
  ####################################hypothesis testing for F3#########################################
  dataF3 <- sort(dataF3); 
  F3w1<-1/(12*m_sam[5])
  F3bmw <- matrix(0,m_sam[5],1); F3bmwsl <- matrix(0,m_sam[5],d6)
  for(i in 1:d6){
    F3gg <- (dataF3 - m5[i])/sqrt(sigma5[i])
    F3bmw[which(F3gg>=0)] <- pnorm(F3gg[F3gg>=0])
    F3bmw[which(F3gg<0)] <- 1 - pnorm(abs(F3gg[F3gg<0]))
    F3bmwsl[,i] <- F3bmw*mi5[i]
  }
  F3P2 <- rowSums(F3bmwsl)
  #############deal with ties ####################
  nn <- dim(as.matrix(unique(F3P2)))[1]
  if(nn < m_sam[5]){F3P2 <- F3P2+runif(m_sam[5])/1e4}
  ##########################################################
  F3dd <- as.matrix(c(sum(F3P2),sum(F3P2^2),sum((F3P2-0.5)^2)))
  F3WW2 <- 1/(12*m_sam[5]) + sum((F3P2 - (as.matrix(c(1:m_sam[5])) - 0.5)/m_sam[5])^2)
  F3u <- as.matrix(c(12*m_sam[5]*((F3dd[1]/m_sam[5]-0.5)^2),((45*m_sam[5])/4)*((F3dd[2]/m_sam[5]-1/3)^2),180*m_sam[5]*((F3dd[3]/m_sam[5]-1/12)^2)))
  F3D <- as.numeric(ks.test(F3P2,"punif")[[1]][1])
  F3tt <- as.matrix(c((1 - pchisq(F3u[1],1)),(1 - pchisq(F3u[2],1)),(1 - pchisq(F3u[3],1)),K1(F3WW2),(1-pkolm(F3D,m_sam[5]))))
  F3tt[which( F3tt>=10e-4)]<-round(F3tt[which(F3tt>=10e-4)],4);F3tt[which(F3tt<10e-4)]<-format(F3tt[which(F3tt<10e-4)],scientific=TRUE,digit=4)
  
  output <- data.frame("MX2-EAD-AD",round(abc,4),round(AIC,4),round(m[1],4),round(m[2],4),round(m[3],4), round(t(m4),4)," "," "," "," "," "," ",round(sigma4[1],4),
                       round(t(mi4),4)," "," "," "," "," "," ",round(t(m5),4)," "," "," ",round(t(sigma5),4)," "," "," ",round(t(mi5),4)," "," "," ",round(sigma,4),          
                       round(B231111[1],4)," "," "," "," ",round(B231111[2],4)," "," "," "," "," "," "," ",round(B231111[3],4),round(B231111[4],4),round(F2_jj,4),round(F2_ll*100,4),round(F2_gg,4),round(F2_rr*100,4),round(F3_jj,4),round(F3_ll*100,4),round(F3_gg,4),round(F3_rr*100,4),
                       round(P1u[1],4),P1tt[1],round(P1u[2],4),P1tt[2],round(P1u[3],4),P1tt[3],round(P1w,4),P1tt[4],round(P1D,4),P1tt[5],
                       round(F1u[1],4),F1tt[1],round(F1u[2],4),F1tt[2],round(F1u[3],4),F1tt[3],round(F1w,4),F1tt[4],round(F1D,4),F1tt[5],
                       round(P2u[1],4),P2tt[1],round(P2u[2],4),P2tt[2],round(P2u[3],4),P2tt[3],round(P2w,4),P2tt[4],round(P2D,4),P2tt[5],
                       round(F2u[1],4),F2tt[1],round(F2u[2],4),F2tt[2],round(F2u[3],4),F2tt[3],round(F2WW2,4),F2tt[4],round(F2D,4),F2tt[5],
                       round(F3u[1],4),F3tt[1],round(F3u[2],4),F3tt[2],round(F3u[3],4),F3tt[3],round(F3WW2,4),F3tt[4],round(F3D,4),F3tt[5])
 output<-as.matrix(output)
  OUTPUT<-list(output,mix_pi4,mix_pi5)
  return(OUTPUT)
  
}


K1G5 <- function(x){
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

logLG5 <- function(nm,nng,mi,mn,s,d1) { sum2 <- sum(log(dmixnorm(d1,mn,sqrt(s),mi)));return (sum2) } 


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
  G5ModelFun[[i]](K1G5,logLG5,df11,df21,df31,df41,df51,G5text2)[[1]]
}
stopCluster(cl)
mix_pi4<-NULL;mix_pi5<-NULL

}else{
  
  allresultq=switch(model,"1MG-AD" = G5ModelFun[[1]](K1G5,logLG5,df11,df21,df31,df41,df51,G5text2),"1MG-A"=G5ModelFun[[2]](K1G5,logLG5,df11,df21,df31,df41,df51,G5text2),"1MG-EAD"=G5ModelFun[[3]](K1G5,logLG5,df11,df21,df31,df41,df51,G5text2),"1MG-NCD"=G5ModelFun[[4]](K1G5,logLG5,df11,df21,df31,df41,df51,G5text2),"2MG-ADI"=G5ModelFun[[5]](K1G5,logLG5,df11,df21,df31,df41,df51,G5text2),
                    "2MG-AD"=G5ModelFun[[6]](K1G5,logLG5,df11,df21,df31,df41,df51,G5text2),"2MG-A"=G5ModelFun[[7]](K1G5,logLG5,df11,df21,df31,df41,df51,G5text2),"2MG-EA"=G5ModelFun[[8]](K1G5,logLG5,df11,df21,df31,df41,df51,G5text2),"2MG-CD"=G5ModelFun[[9]](K1G5,logLG5,df11,df21,df31,df41,df51,G5text2),"2MG-EAD"=G5ModelFun[[10]](K1G5,logLG5,df11,df21,df31,df41,df51,G5text2),
                    "PG-ADI"=G5ModelFun[[11]](K1G5,logLG5,df11,df21,df31,df41,df51,G5text2),"PG-AD"=G5ModelFun[[12]](K1G5,logLG5,df11,df21,df31,df41,df51,G5text2),"MX1-AD-ADI"=G5ModelFun[[13]](K1G5,logLG5,df11,df21,df31,df41,df51,G5text2),"MX1-AD-AD"=G5ModelFun[[14]](K1G5,logLG5,df11,df21,df31,df41,df51,G5text2),"MX1-A-AD"=G5ModelFun[[15]](K1G5,logLG5,df11,df21,df31,df41,df51,G5text2),
                    "MX1-EAD-AD"=G5ModelFun[[16]](K1G5,logLG5,df11,df21,df31,df41,df51,G5text2),"MX1-NCD-AD"=G5ModelFun[[17]](K1G5,logLG5,df11,df21,df31,df41,df51,G5text2),"MX2-ADI-ADI"=G5ModelFun[[18]](K1G5,logLG5,df11,df21,df31,df41,df51,G5text2),"MX2-ADI-AD"=G5ModelFun[[19]](K1G5,logLG5,df11,df21,df31,df41,df51,G5text2),"MX2-AD-AD"=G5ModelFun[[20]](K1G5,logLG5,df11,df21,df31,df41,df51,G5text2),
                    "MX2-A-AD"=G5ModelFun[[21]](K1G5,logLG5,df11,df21,df31,df41,df51,G5text2),"MX2-EA-AD"=G5ModelFun[[22]](K1G5,logLG5,df11,df21,df31,df41,df51,G5text2),"MX2-CD-AD"=G5ModelFun[[23]](K1G5,logLG5,df11,df21,df31,df41,df51,G5text2),"MX2-EAD-AD"=G5ModelFun[[24]](K1G5,logLG5,df11,df21,df31,df41,df51,G5text2))
  
  
  allresult<-allresultq[[1]]
  if(model=="PG-AD"||model=="PG-ADI"){
    mix_pi4<-NULL;mix_pi5<-NULL
  }else{
    mix_pi4<-allresultq[[2]];mix_pi5<-allresultq[[3]] 
  }
}
colnames(allresult) <- G5colname
out<-list(allresult,mix_pi4,mix_pi5)
return(out)
} 

  
  
  
  
  
