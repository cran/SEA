PosPro<-function(Population,result,data){
  
  if(Population=="F2"){
    data<-sapply(data,as.character)
    dF2<-data[-1,which(data[1,]=="F2")];F2<-as.numeric(dF2[which(is.na(as.numeric(dF2))==FALSE)]);data<-as.matrix(F2) 
    
    datac<-result[[1]];mi<- result[[2]]
    mm<-datac[4:12]
    if(length(mm[-which(mm[]==" ")])==0)
    {m<-as.matrix(as.numeric(mm))}else{m<-as.matrix(as.numeric(mm[-which(mm[]==" ")]))}
    sigma<-matrix(as.numeric(datac[13]),dim(m)[1],1)
    
    WW <- matrix(0,dim(mi)[1],dim(data)[1])
    for(i in 1:dim(mi)[1]) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    WW<-round(WW,4)
    RowBl<-matrix(" ",(9-dim(mi)[1]),dim(data)[1])
    PPW<-t(rbind(WW,RowBl))
    dimnames(PPW)<-list(c(1:dim(data)[1]),c("F2(1)","F2(2)","F2(3)","F2(4)","F2(5)","F2(6)","F2(7)","F2(8)","F2(9)"))
    
  }else if(Population=="F2:3"){
    data<-sapply(data,as.character)
    dF23<-data[-1,which(data[1,]=="F23")];F23<-as.numeric(dF23[which(is.na(as.numeric(dF23))==FALSE)]);data<-as.matrix(F23) 
   
    datac<-result[[1]];mi<- result[[2]]
    mm<-datac[4:12]
    if(length(mm[-which(mm[]==" ")])==0)
    {m<-as.matrix(as.numeric(mm))}else{m<-as.matrix(as.numeric(mm[-which(mm[]==" ")]))}
    ssigma<-datac[13:21]
    if(length(ssigma[-which(ssigma[]==" ")])==0)
    {sigma<-as.matrix(as.numeric(ssigma))}else{sigma<-as.matrix(as.numeric(ssigma[-which(ssigma[]==" ")]))}
    
    WW <- matrix(0,dim(mi)[1],dim(data)[1])
    for(i in 1:dim(mi)[1]) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    WW<-round(WW,4)
    RowBl<-matrix(" ",(9-dim(mi)[1]),dim(data)[1])
    PPW<-t(rbind(WW,RowBl))
    dimnames(PPW)<-list(c(1:dim(data)[1]),c("F2:3(1)","F2:3(2)","F2:3(3)","F2:3(4)","F2:3(5)","F2:3(6)","F2:3(7)","F2:3(8)","F2:3(9)"))
  }else if(Population=="DH"){
    data<-sapply(data,as.character)
    dDH<-data[-1,which(data[1,]=="DH")];DH<-as.numeric(dDH[which(is.na(as.numeric(dDH))==FALSE)]);data<-as.matrix(DH) 
    
    datac<-result[[1]];mi<- result[[2]]
    mm4<-datac[4:19]
    if(length(mm4[-which(mm4[]==" ")])==0)
    {m4<-as.matrix(as.numeric(mm4))}else{m4<-as.matrix(as.numeric(mm4[-which(mm4[]==" ")]))}  
    sigma4<-matrix(as.numeric(datac[20]),dim(m4)[1],1) 
    
    W4 <- matrix(0,dim(mi)[1],dim(data)[1])
    for(i in 1:dim(mi)[1]) { W4[i,] <- mi[i]*dnorm(data,m4[i],sqrt(sigma4[i]))/dmixnorm(data,m4,sqrt(sigma4),mi) } 
    PW4<-round(W4,4)
    RowBl4<-matrix(" ",(16-dim(mi)[1]),dim(data)[1])
    PPW4<-rbind(PW4,RowBl4)
    PPW<-t(PPW4)
    dimnames(PPW)<-list(c(1:dim(data)[1]),c("DH(1)","DH(2)","DH(3)","DH(4)","DH(5)","DH(6)","DH(7)","DH(8)",
                                            "DH(9)","DH(10)","DH(11)","DH(12)","DH(13)","DH(14)","DH(15)","DH(16)"))
    
  }else if(Population=="BIL"){
    data<-sapply(data,as.character)
    dBIL<-data[-1,which(data[1,]=="BIL")];BIL<-as.numeric(dBIL[which(is.na(as.numeric(dBIL))==FALSE)]);data<-as.matrix(BIL) 
    
    datac<-result[[1]];mi<- result[[2]]
    mm<-datac[4:11]
    if(length(mm[-which(mm[]==" ")])==0)
    {m<-as.matrix(as.numeric(mm))}else{m<-as.matrix(as.numeric(mm[-which(mm[]==" ")]))}  
    sigma<-matrix(as.numeric(datac[12]),dim(m)[1],1) 
    
    WW <- matrix(0,dim(mi)[1],dim(data)[1])
    for(i in 1:dim(mi)[1]) { WW[i,] <- mi[i]*dnorm(data,m[i],sqrt(sigma[i]))/dmixnorm(data,m,sqrt(sigma),mi) }
    WW<-round(WW,4)
    RowBl<-matrix(" ",(8-dim(mi)[1]),dim(data)[1])
    PPW<-t(rbind(WW,RowBl))
    dimnames(PPW)<-list(c(1:dim(data)[1]),c("BIL(1)","BIL(2)","BIL(3)","BIL(4)","BIL(5)","BIL(6)","BIL(7)","BIL(8)"))
    
  }else if(Population=="BC (B1 B2)"){
    
    data<-sapply(data,as.character)
    dB1<-data[-1,which(data[1,]=="B1")];B1<-as.numeric(dB1[which(is.na(as.numeric(dB1))==FALSE)]);dataB1<-as.matrix(B1)
    dB2<-data[-1,which(data[1,]=="B2")];B2<-as.numeric(dB2[which(is.na(as.numeric(dB2))==FALSE)]);dataB2<-as.matrix(B2)
    
    datac<-result[[1]];mi_1<- result[[2]];mi_2<- result[[3]]
    mm4<-datac[4:7]
    if(length(mm4[-which(mm4[]==" ")])==0)
    {m4<-as.matrix(as.numeric(mm4))}else{m4<-as.matrix(as.numeric(mm4[-which(mm4[]==" ")]))}  
    sigma4<-matrix(as.numeric(datac[8]),dim(m4)[1],1) 
    
    mm5<-datac[13:16]
    if(length(mm5[-which(mm5[]==" ")])==0)
    {m5<-as.matrix(as.numeric(mm5))}else{m5<-as.matrix(as.numeric(mm5[-which(mm5[]==" ")]))}  
    sigma5<-matrix(as.numeric(datac[17]),dim(m5)[1],1) 
    
    W4 <- matrix(0,dim(mi_1)[1],dim(dataB1)[1]);W5 <- matrix(0,dim(mi_2)[1],dim(dataB2)[1])
    for(i in 1:dim(mi_1)[1]) { W4[i,] <- mi_1[i]*dnorm(dataB1,m4[i],sqrt(sigma4[i]))/dmixnorm(dataB1,m4,sqrt(sigma4),mi_1) }
    for(i in 1:dim(mi_2)[1]) { W5[i,] <- mi_2[i]*dnorm(dataB2,m5[i],sqrt(sigma5[i]))/dmixnorm(dataB2,m5,sqrt(sigma5),mi_2) }
    ColBlNum<-abs(dim(dataB1)[1]-dim(dataB2)[1])
    ColBl4<-matrix(" ",dim(mi_1)[1],ColBlNum)
    ColBl5<-matrix(" ",dim(mi_2)[1],ColBlNum)
    W4<-round(W4,4);W5<-round(W5,4)
    if(dim(dataB1)[1]<dim(dataB2)[1])
    {
      PW4<-cbind(W4,ColBl4)
      PW5<-W5
    }else{
      PW4<-W4 
      PW5<-cbind(W5,ColBl5) 
    }  
    RowBl4<-matrix(" ",(4-dim(mi_1)[1]),max(dim(dataB1)[1],dim(dataB2)[1]))
    PPW4<-rbind(PW4,RowBl4)
    RowBl5<-matrix(" ",(4-dim(mi_2)[1]),max(dim(dataB1)[1],dim(dataB2)[1]))
    PPW5<-rbind(PW5,RowBl5)
    PPW<-t(rbind(PPW4,PPW5))
    dimnames(PPW)<-list(c(1:max(dim(dataB1)[1],dim(dataB2)[1])),c("B1(1)","B1(2)","B1(3)","B1(4)","B2(1)","B2(2)","B2(3)","B2(4)"))
  }else if(Population=="BCF (B1:2 B2:2)"){
    data<-sapply(data,as.character)
    dB12<-data[-1,which(data[1,]=="B12")];B12<-as.numeric(dB12[which(is.na(as.numeric(dB12))==FALSE)]);dataB1<-as.matrix(B12)
    dB22<-data[-1,which(data[1,]=="B22")];B22<-as.numeric(dB22[which(is.na(as.numeric(dB22))==FALSE)]);dataB2<-as.matrix(B22)
    
    datac<-result[[1]];mi_1<- result[[2]];mi_2<- result[[3]]
    
    mm4<-datac[4:7]
    if(length(mm4[-which(mm4[]==" ")])==0)
    {m4<-as.matrix(as.numeric(mm4))}else{m4<-as.matrix(as.numeric(mm4[-which(mm4[]==" ")]))}  
    ssigma4<-datac[8:11]
    if(length(ssigma4[-which(ssigma4[]==" ")])==0)
    {sigma4<-as.matrix(as.numeric(ssigma4))}else{sigma4<-as.matrix(as.numeric(ssigma4[-which(ssigma4[]==" ")]))}  
    
    mm5<-datac[16:19]
    if(length(mm5[-which(mm5[]==" ")])==0)
    {m5<-as.matrix(as.numeric(mm5))}else{m5<-as.matrix(as.numeric(mm5[-which(mm5[]==" ")]))}  
    ssigma5<-datac[20:23]
    if(length(ssigma5[-which(ssigma5[]==" ")])==0)
    {sigma5<-as.matrix(as.numeric(ssigma5))}else{sigma5<-as.matrix(as.numeric(ssigma5[-which(ssigma5[]==" ")]))}  
    
    W4 <- matrix(0,dim(mi_1)[1],dim(dataB1)[1]);W5 <- matrix(0,dim(mi_2)[1],dim(dataB2)[1])
    for(i in 1:dim(mi_1)[1]) { W4[i,] <- mi_1[i]*dnorm(dataB1,m4[i],sqrt(sigma4[i]))/dmixnorm(dataB1,m4,sqrt(sigma4),mi_1) }
    for(i in 1:dim(mi_2)[1]) { W5[i,] <- mi_2[i]*dnorm(dataB2,m5[i],sqrt(sigma5[i]))/dmixnorm(dataB2,m5,sqrt(sigma5),mi_2) }
    ColBlNum<-abs(dim(dataB1)[1]-dim(dataB2)[1])
    ColBl4<-matrix(" ",dim(mi_1)[1],ColBlNum)
    ColBl5<-matrix(" ",dim(mi_2)[1],ColBlNum)
    W4<-round(W4,4);W5<-round(W5,4)
    if(dim(dataB1)[1]<dim(dataB2)[1])
    {
      PW4<-cbind(W4,ColBl4)
      PW5<-W5
    }else{
      PW4<-W4 
      PW5<-cbind(W5,ColBl5) 
    }  
    RowBl4<-matrix(" ",(4-dim(mi_1)[1]),max(dim(dataB1)[1],dim(dataB2)[1]))
    PPW4<-rbind(PW4,RowBl4)
    RowBl5<-matrix(" ",(4-dim(mi_2)[1]),max(dim(dataB1)[1],dim(dataB2)[1]))
    PPW5<-rbind(PW5,RowBl5)
    PPW<-t(rbind(PPW4,PPW5))
    dimnames(PPW)<-list(c(1:max(dim(dataB1)[1],dim(dataB2)[1])),c("B1:2(1)","B1:2(2)","B1:2(3)","B1:2(4)","B2:2(1)","B2:2(2)","B2:2(3)","B2:2(4)"))
    
  }else if(Population=="G4F2 (P1 P2 F1 F2)"){
    data<-sapply(data,as.character)
    dF2<-data[-1,which(data[1,]=="F2")];F2<-as.numeric(dF2[which(is.na(as.numeric(dF2))==FALSE)]);dataF2<-as.matrix(F2)
    
    datac<-result[[1]];mi<- result[[2]]
    
    mm<-datac[7:15]
    if(length(mm[-which(mm[]==" ")])==0)
    {m<-as.matrix(as.numeric(mm))}else{m<-as.matrix(as.numeric(mm[-which(mm[]==" ")]))}  
    sigma<-matrix(as.numeric(datac[16]),dim(m)[1],1) 
    
    WW <- matrix(0,dim(mi)[1],dim(dataF2)[1]);
    for(i in 1:dim(mi)[1]) { WW[i,] <- mi[i]*dnorm(dataF2,m[i],sqrt(sigma[i]))/dmixnorm(dataF2,m,sqrt(sigma),mi) }
    WW<-round(WW,4)
    RowBl<-matrix(" ",(9-dim(mi)[1]),dim(dataF2)[1])
    PPW<-t(rbind(WW,RowBl))
    dimnames(PPW)<-list(c(1:dim(dataF2)[1]),c("F2(1)","F2(2)","F2(3)","F2(4)","F2(5)","F2(6)","F2(7)","F2(8)","F2(9)"))
    
  }else if(Population=="G4F3 (P1 P2 F1 F2:3)"){
    data<-sapply(data,as.character)
    dF23<-data[-1,which(data[1,]=="F23")];F23<-as.numeric(dF23[which(is.na(as.numeric(dF23))==FALSE)]);dataF3<-as.matrix(F23)
    
    datac<-result[[1]];mi<- result[[2]]
    mm<-datac[7:15]
    if(length(mm[-which(mm[]==" ")])==0)
    {m<-as.matrix(as.numeric(mm))}else{m<-as.matrix(as.numeric(mm[-which(mm[]==" ")]))}  
    
    ssigma<-datac[17:25]
    if(length(ssigma[-which(ssigma[]==" ")])==0)
    {sigma<-as.matrix(as.numeric(ssigma))}else{sigma<-as.matrix(as.numeric(ssigma[-which(ssigma[]==" ")]))}  
    
    WW <- matrix(0,dim(mi)[1],dim(dataF3)[1]);
    for(i in 1:dim(mi)[1]) { WW[i,] <- mi[i]*dnorm(dataF3,m[i],sqrt(sigma[i]))/dmixnorm(dataF3,m,sqrt(sigma),mi) }
    WW<-round(WW,4)
    RowBl<-matrix(" ",(9-dim(mi)[1]),dim(dataF3)[1])
    PPW<-t(rbind(WW,RowBl))
    dimnames(PPW)<-list(c(1:dim(dataF3)[1]),c("F3(1)","F3(2)","F3(3)","F3(4)","F3(5)","F3(6)","F3(7)","F3(8)","F3(9)"))
    
  }else if(Population=="G3DH (P1 P2 DH)"){
    data<-sapply(data,as.character)
    dDH<-data[-1,which(data[1,]=="DH")];DH<-as.numeric(dDH[which(is.na(as.numeric(dDH))==FALSE)]);dataDH<-as.matrix(DH)
    
    datac<-result[[1]];mi<- result[[2]]
    
    mm4<-datac[7:22]
    if(length(mm4[-which(mm4[]==" ")])==0)
    {m4<-as.matrix(as.numeric(mm4))}else{m4<-as.matrix(as.numeric(mm4[-which(mm4[]==" ")]))}  
    sigma4<-matrix(as.numeric(datac[23]),dim(m4)[1],1) 
    
    W4 <- matrix(0,dim(mi)[1],dim(dataDH)[1])
    for(i in 1:dim(mi)[1]) { W4[i,] <- mi[i]*dnorm(dataDH,m4[i],sqrt(sigma4[i]))/dmixnorm(dataDH,m4,sqrt(sigma4),mi) }
    PW4<-round(W4,4)
    RowBl4<-matrix(" ",(16-dim(mi)[1]),dim(dataDH)[1])
    PPW<-t(rbind(PW4,RowBl4))
    dimnames(PPW)<-list(c(1:dim(dataDH)[1]),c("DH(1)","DH(2)","DH(3)","DH(4)","DH(5)","DH(6)","DH(7)","DH(8)","DH(9)","DH(10)","DH(11)","DH(12)","DH(13)","DH(14)","DH(15)","DH(16)"))
    
  }else if(Population=="G5BC (P1 P2 F1 B1 B2)"){
    data<-sapply(data,as.character)
    dB1<-data[-1,which(data[1,]=="B1")];B1<-as.numeric(dB1[which(is.na(as.numeric(dB1))==FALSE)]);dataB1<-as.matrix(B1)
    dB2<-data[-1,which(data[1,]=="B2")];B2<-as.numeric(dB2[which(is.na(as.numeric(dB2))==FALSE)]);dataB2<-as.matrix(B2)
    
    datac<-result[[1]];mi_1<- result[[2]];mi_2<- result[[3]]
    
    mm4<-datac[8:11]
    if(length(mm4[-which(mm4[]==" ")])==0)
    {m4<-as.matrix(as.numeric(mm4))}else{m4<-as.matrix(as.numeric(mm4[-which(mm4[]==" ")]))}  
    sigma4<-matrix(as.numeric(datac[12]),dim(m4)[1],1) 
    mm5<-datac[17:20]
    if(length(mm5[-which(mm5[]==" ")])==0)
    {m5<-as.matrix(as.numeric(mm5))}else{m5<-as.matrix(as.numeric(mm5[-which(mm5[]==" ")]))}  
    sigma5<-matrix(as.numeric(datac[21]),dim(m5)[1],1) 
    
    W4 <- matrix(0,dim(mi_1)[1],dim(dataB1)[1]);W5 <- matrix(0,dim(mi_2)[1],dim(dataB2)[1])
    for(i in 1:dim(mi_1)[1]) { W4[i,] <- mi_1[i]*dnorm(dataB1,m4[i],sqrt(sigma4[i]))/dmixnorm(dataB1,m4,sqrt(sigma4),mi_1) }
    for(i in 1:dim(mi_2)[1]) { W5[i,] <- mi_2[i]*dnorm(dataB2,m5[i],sqrt(sigma5[i]))/dmixnorm(dataB2,m5,sqrt(sigma5),mi_2) }
    ColBlNum<-abs(dim(dataB1)[1]-dim(dataB2)[1])
    ColBl4<-matrix(" ",dim(mi_1)[1],ColBlNum)
    ColBl5<-matrix(" ",dim(mi_2)[1],ColBlNum)
    W4<-round(W4,4);W5<-round(W5,4)
    if(dim(dataB1)[1]<dim(dataB2)[1])
    {
      PW4<-cbind(W4,ColBl4)
      PW5<-W5
    }else{
      PW4<-W4 
      PW5<-cbind(W5,ColBl5) 
    }  
    RowBl4<-matrix(" ",(4-dim(mi_1)[1]),max(dim(dataB1)[1],dim(dataB2)[1]))
    PPW4<-rbind(PW4,RowBl4)
    RowBl5<-matrix(" ",(4-dim(mi_2)[1]),max(dim(dataB1)[1],dim(dataB2)[1]))
    PPW5<-rbind(PW5,RowBl5)
    PPW<-t(rbind(PPW4,PPW5))
    dimnames(PPW)<-list(c(1:max(dim(dataB1)[1],dim(dataB2)[1])),c("B1(1)","B1(2)","B1(3)","B1(4)","B2(1)","B2(2)","B2(3)","B2(4)"))
    
  }else if(Population=="G5BCF (P1 P2 F1 B1:2 B2:2)"){
    data<-sapply(data,as.character)
    dB12<-data[-1,which(data[1,]=="B12")];B12<-as.numeric(dB12[which(is.na(as.numeric(dB12))==FALSE)]);dataB1<-as.matrix(B12)
    dB22<-data[-1,which(data[1,]=="B22")];B22<-as.numeric(dB22[which(is.na(as.numeric(dB22))==FALSE)]);dataB2<-as.matrix(B22)
    
    datac<-result[[1]];mmi1<- result[[2]];mmi2<- result[[3]]
    
    mm1<-datac[7:10]
    if(length(mm1[-which(mm1[]==" ")])==0)
    {m1<-as.matrix(as.numeric(mm1))}else{m1<-as.matrix(as.numeric(mm1[-which(mm1[]==" ")]))}  
    ssigma1<-datac[11:14]
    if(length(ssigma1[-which(ssigma1[]==" ")])==0)
    {sigma1<-as.matrix(as.numeric(ssigma1))}else{sigma1<-as.matrix(as.numeric(ssigma1[-which(ssigma1[]==" ")]))} 
    mm2<-datac[19:22]
    if(length(mm2[-which(mm2[]==" ")])==0)
    {m2<-as.matrix(as.numeric(mm2))}else{m2<-as.matrix(as.numeric(mm2[-which(mm2[]==" ")]))}  
    ssigma2<-datac[23:26]
    if(length(ssigma2[-which(ssigma2[]==" ")])==0)
    {sigma2<-as.matrix(as.numeric(ssigma2))}else{sigma2<-as.matrix(as.numeric(ssigma2[-which(ssigma2[]==" ")]))} 
    
    
    W1 <- matrix(0,dim(mmi1)[1],dim(dataB1)[1]);W2 <- matrix(0,dim(mmi2)[1],dim(dataB2)[1])
    for(i in 1:dim(mmi1)[1]) { W1[i,] <- mmi1[i]*dnorm(dataB1,m1[i],sqrt(sigma1[i]))/dmixnorm(dataB1,m1,sqrt(sigma1),mmi1) }
    for(i in 1:dim(mmi2)[1]) { W2[i,] <- mmi2[i]*dnorm(dataB2,m2[i],sqrt(sigma2[i]))/dmixnorm(dataB2,m2,sqrt(sigma2),mmi2) }
    ColBlNum<-abs(dim(dataB1)[1]-dim(dataB2)[1])
    ColBl1<-matrix(" ",dim(mmi1)[1],ColBlNum)
    ColBl2<-matrix(" ",dim(mmi2)[1],ColBlNum)
    W1<-round(W1,4);W2<-round(W2,4)
    if(dim(dataB1)[1]<dim(dataB2)[1])
    {
      PW1<-cbind(W1,ColBl1)
      PW2<-W2
    }else{
      PW1<-W1 
      PW2<-cbind(W2,ColBl2) 
    }  
    RowBl1<-matrix(" ",(4-dim(mmi1)[1]),max(dim(dataB1)[1],dim(dataB2)[1]))
    PPW1<-rbind(PW1,RowBl1)
    RowBl2<-matrix(" ",(4-dim(mmi2)[1]),max(dim(dataB1)[1],dim(dataB2)[1]))
    PPW2<-rbind(PW2,RowBl2)
    PPW<-t(rbind(PPW1,PPW2))
    dimnames(PPW)<-list(c(1:max(dim(dataB1)[1],dim(dataB2)[1])),c("B1:2(1)","B1:2(2)","B1:2(3)","B1:2(4)","B2:2(1)","B2:2(2)","B2:2(3)","B2:2(4)"))
    
  }else if(Population=="G5 (P1 P2 F1 F2 F2:3)"){
    data<-sapply(data,as.character)
    dF2<-data[-1,which(data[1,]=="F2")];F2<-as.numeric(dF2[which(is.na(as.numeric(dF2))==FALSE)]);dataF2<-as.matrix(F2)
    dF23<-data[-1,which(data[1,]=="F23")];F23<-as.numeric(dF23[which(is.na(as.numeric(dF23))==FALSE)]);dataF3<-as.matrix(F23)
    
    datac<-result[[1]];mix_pi4<- result[[2]];mix_pi5<- result[[3]]
    
    mm4<-datac[7:15]
    if(length(mm4[-which(mm4[]==" ")])==0)
    {m4<-as.matrix(as.numeric(mm4))}else{m4<-as.matrix(as.numeric(mm4[-which(mm4[]==" ")]))}  
    sigma4<-matrix(as.numeric(datac[16]),dim(m4)[1],1) 
    mm5<-datac[26:34]
    if(length(mm5[-which(mm5[]==" ")])==0)
    {m5<-as.matrix(as.numeric(mm5))}else{m5<-as.matrix(as.numeric(mm5[-which(mm5[]==" ")]))}  
    ssigma5<-datac[35:43]
    if(length(ssigma5[-which(ssigma5[]==" ")])==0)
    {sigma5<-as.matrix(as.numeric(ssigma5))}else{sigma5<-as.matrix(as.numeric(ssigma5[-which(ssigma5[]==" ")]))}  
    
    
    W4 <- matrix(0,dim(mix_pi4)[1],dim(dataF2)[1]);W5 <- matrix(0,dim(mix_pi5)[1],dim(dataF3)[1])
    for(i in 1:dim(mix_pi4)[1]) { W4[i,] <- mix_pi4[i]*dnorm(dataF2,m4[i],sqrt(sigma4[i]))/dmixnorm(dataF2,m4,sqrt(sigma4),mix_pi4) }
    for(i in 1:dim(mix_pi5)[1]) { W5[i,] <- mix_pi5[i]*dnorm(dataF3,m5[i],sqrt(sigma5[i]))/dmixnorm(dataF3,m5,sqrt(sigma5),mix_pi5) }
    ColBlNum<-abs(dim(dataF2)[1]-dim(dataF3)[1])
    ColBl4<-matrix(" ",dim(mix_pi4)[1],ColBlNum)
    ColBl5<-matrix(" ",dim(mix_pi5)[1],ColBlNum)
    W4<-round(W4,4);W5<-round(W5,4)
    if(dim(dataF2)[1]<dim(dataF3)[1])
    {
      PW4<-cbind(W4,ColBl4)
      PW5<-W5
    }else{
      PW4<-W4 
      PW5<-cbind(W5,ColBl5) 
    }  
    RowBl4<-matrix(" ",(9-dim(mix_pi4)[1]),max(dim(dataF2)[1],dim(dataF3)[1]))
    PPW4<-rbind(PW4,RowBl4)
    RowBl5<-matrix(" ",(9-dim(mix_pi5)[1]),max(dim(dataF2)[1],dim(dataF3)[1]))
    PPW5<-rbind(PW5,RowBl5)
    PPW<-t(rbind(PPW4,PPW5))
    dimnames(PPW)<-list(c(1:max(dim(dataF2)[1],dim(dataF3)[1])),c("F2(1)","F2(2)","F2(3)","F2(4)","F2(5)","F2(6)","F2(7)","F2(8)","F2(9)",
                                                                  "F2:3(1)","F2:3(2)","F2:3(3)","F2:3(4)","F2:3(5)","F2:3(6)","F2:3(7)","F2:3(8)","F2:3(9)"))
  }else if(Population=="G6 (P1 P2 F1 F2 B1 B2)"){
    data<-sapply(data,as.character)
    dB1<-data[-1,which(data[1,]=="B1")];B1<-as.numeric(dB1[which(is.na(as.numeric(dB1))==FALSE)]);dataB1<-as.matrix(B1)
    dB2<-data[-1,which(data[1,]=="B2")];B2<-as.numeric(dB2[which(is.na(as.numeric(dB2))==FALSE)]);dataB2<-as.matrix(B2)
    dF2<-data[-1,which(data[1,]=="F2")];F2<-as.numeric(dF2[which(is.na(as.numeric(dF2))==FALSE)]);dataF2<-as.matrix(F2)
    
    datac<-result[[1]];mi_1<- result[[2]];mi_2<- result[[3]];mi_3<- result[[4]]
    
    mm4<-datac[8:11]
    if(length(mm4[-which(mm4[]==" ")])==0)
    {m4<-as.matrix(as.numeric(mm4))}else{m4<-as.matrix(as.numeric(mm4[-which(mm4[]==" ")]))}  
    sigma4<-matrix(as.numeric(datac[12]),dim(m4)[1],1) 
    
    mm5<-datac[17:20]
    if(length(mm5[-which(mm5[]==" ")])==0)
    {m5<-as.matrix(as.numeric(mm5))}else{m5<-as.matrix(as.numeric(mm5[-which(mm5[]==" ")]))}  
    sigma5<-matrix(as.numeric(datac[21]),dim(m5)[1],1) 
    
    mm6<-datac[26:34]
    if(length(mm6[-which(mm6[]==" ")])==0)
    {m6<-as.matrix(as.numeric(mm6))}else{m6<-as.matrix(as.numeric(mm6[-which(mm6[]==" ")]))}  
    sigma6<-matrix(as.numeric(datac[35]),dim(m6)[1],1) 
    
    W4 <- matrix(0,dim(mi_1)[1],dim(dataB1)[1]); W5 <- matrix(0,dim(mi_2)[1],dim(dataB2)[1]); W6 <- matrix(0,dim(mi_3)[1],dim(dataF2)[1])
    for(i in 1:dim(mi_1)[1]) { W4[i,] <- mi_1[i]*dnorm(dataB1,m4[i],sqrt(sigma4[i]))/dmixnorm(dataB1,m4,sqrt(sigma4),mi_1) }
    for(i in 1:dim(mi_2)[1]) { W5[i,] <- mi_2[i]*dnorm(dataB2,m5[i],sqrt(sigma5[i]))/dmixnorm(dataB2,m5,sqrt(sigma5),mi_2) }
    for(i in 1:dim(mi_3)[1]) { W6[i,] <- mi_3[i]*dnorm(dataF2,m6[i],sqrt(sigma6[i]))/dmixnorm(dataF2,m6,sqrt(sigma6),mi_3) }
    Colmax<-max(dim(dataB1)[1],dim(dataB2)[1],dim(dataF2)[1])
    ColBlNum<-c(Colmax-dim(dataB1)[1],Colmax-dim(dataB2)[1],Colmax-dim(dataF2)[1])
    ColBl4<-matrix(" ",dim(mi_1)[1],ColBlNum[1])
    ColBl5<-matrix(" ",dim(mi_2)[1],ColBlNum[2])
    ColBl6<-matrix(" ",dim(mi_3)[1],ColBlNum[3])
    W4<-round(W4,4);W5<-round(W5,4);W6<-round(W6,4)
    PW4<-cbind(W4,ColBl4); PW5<-cbind(W5,ColBl5); PW6<-cbind(W6,ColBl6) 
    RowBl4<-matrix(" ",(4-dim(mi_1)[1]),Colmax)
    PPW4<-rbind(PW4,RowBl4)
    RowBl5<-matrix(" ",(4-dim(mi_2)[1]),Colmax)
    PPW5<-rbind(PW5,RowBl5)
    RowBl6<-matrix(" ",(9-dim(mi_3)[1]),Colmax)
    PPW6<-rbind(PW6,RowBl6)
    
    PPW<-t(rbind(PPW4,PPW5,PPW6))
    dimnames(PPW)<-list(c(1:Colmax),c("B1(1)","B1(2)","B1(3)","B1(4)","B2(1)","B2(2)","B2(3)","B2(4)",
                                      "F2(1)","F2(2)","F2(3)","F2(4)","F2(5)","F2(6)","F2(7)","F2(8)","F2(9)"))
    
  }else if(Population=="G6F (P1 F1 P2 B1:2 B2:2 F2:3)"){
    data<-sapply(data,as.character)
    dB12<-data[-1,which(data[1,]=="B12")];B12<-as.numeric(dB12[which(is.na(as.numeric(dB12))==FALSE)]);dataB1<-as.matrix(B12)
    dB22<-data[-1,which(data[1,]=="B22")];B22<-as.numeric(dB22[which(is.na(as.numeric(dB22))==FALSE)]);dataB2<-as.matrix(B22)
    dF23<-data[-1,which(data[1,]=="F23")];F23<-as.numeric(dF23[which(is.na(as.numeric(dF23))==FALSE)]);dataF2<-as.matrix(F23)
    
    datac<-result[[1]];mi1<- result[[2]];mi2<- result[[3]];mi3<- result[[4]]
    
    mm1<-datac[8:11]
    if(length(mm1[-which(mm1[]==" ")])==0)
    {m1<-as.matrix(as.numeric(mm1))}else{m1<-as.matrix(as.numeric(mm1[-which(mm1[]==" ")]))}  
    ssigma1<-datac[12:15]
    if(length(ssigma1[-which(ssigma1[]==" ")])==0)
    {sigma1<-as.matrix(as.numeric(ssigma1))}else{sigma1<-as.matrix(as.numeric(ssigma1[-which(ssigma1[]==" ")]))}  
    mm2<-datac[20:23]
    if(length(mm2[-which(mm2[]==" ")])==0)
    {m2<-as.matrix(as.numeric(mm2))}else{m2<-as.matrix(as.numeric(mm2[-which(mm2[]==" ")]))}  
    ssigma2<-datac[24:27]
    if(length(ssigma2[-which(ssigma2[]==" ")])==0)
    {sigma2<-as.matrix(as.numeric(ssigma2))}else{sigma2<-as.matrix(as.numeric(ssigma2[-which(ssigma2[]==" ")]))}  
    mm3<-datac[32:40]
    if(length(mm3[-which(mm3[]==" ")])==0)
    {m3<-as.matrix(as.numeric(mm3))}else{m3<-as.matrix(as.numeric(mm3[-which(mm3[]==" ")]))}  
    ssigma3<-datac[41:49]
    if(length(ssigma3[-which(ssigma3[]==" ")])==0)
    {sigma3<-as.matrix(as.numeric(ssigma3))}else{sigma3<-as.matrix(as.numeric(ssigma3[-which(ssigma3[]==" ")]))}  
    
    W1 <- matrix(0,dim(mi1)[1],dim(dataB1)[1]);W2 <- matrix(0,dim(mi2)[1],dim(dataB2)[1]);W3 <- matrix(0,dim(mi3)[1],dim(dataF2)[1])
    for(i in 1:dim(mi1)[1]) { W1[i,] <- mi1[i]*dnorm(dataB1,m1[i],sqrt(sigma1[i]))/dmixnorm(dataB1,m1,sqrt(sigma1),mi1) }
    for(i in 1:dim(mi2)[1]) { W2[i,] <- mi2[i]*dnorm(dataB2,m2[i],sqrt(sigma2[i]))/dmixnorm(dataB2,m2,sqrt(sigma2),mi2) }
    for(i in 1:dim(mi3)[1]) { W3[i,] <- mi3[i]*dnorm(dataF2,m3[i],sqrt(sigma3[i]))/dmixnorm(dataF2,m3,sqrt(sigma3),mi3) }
    ColBl1<-matrix(" ",dim(mi1)[1],abs(dim(dataF2)[1]-dim(dataB1)[1]))
    ColBl2<-matrix(" ",dim(mi2)[1],abs(dim(dataF2)[1]-dim(dataB2)[1]))
    ColBl3<-matrix(" ",dim(mi3)[1],abs(dim(dataB2)[1]-dim(dataB1)[1]))
    W1<-round(W1,4);W2<-round(W2,4);W3<-round(W3,4)
    if(max(dim(dataB1)[1],dim(dataB2)[1],dim(dataF2)[1])==dim(dataF2)[1])
    {
      ColBl1<-matrix(" ",dim(mi1)[1],abs(dim(dataF2)[1]-dim(dataB1)[1]))
      ColBl2<-matrix(" ",dim(mi2)[1],abs(dim(dataF2)[1]-dim(dataB2)[1]))
      PW1<-cbind(W1,ColBl1)
      PW2<-cbind(W2,ColBl2)
      PW3<-W3
    }else if(max(dim(dataB1)[1],dim(dataB2)[1],dim(dataF2)[1])==dim(dataB2)[1])
      
    {
      ColBl3<-matrix(" ",dim(mi1)[1],abs(dim(dataB2)[1]-dim(dataB1)[1]))
      ColBl2<-matrix(" ",dim(mi3)[1],abs(dim(dataF2)[1]-dim(dataB2)[1]))
      PW1<-cbind(W1,ColBl3)
      PW2<-W2
      PW3<-cbind(W3,ColBl2) 
    }else if(max(dim(dataB1)[1],dim(dataB2)[1],dim(dataF2)[1])==dim(dataB1)[1])  
      
    { 
      ColBl3<-matrix(" ",dim(mi2)[1],abs(dim(dataB2)[1]-dim(dataB1)[1]))
      ColBl1<-matrix(" ",dim(mi3)[1],abs(dim(dataF2)[1]-dim(dataB1)[1]))
      PW1<-W1
      PW2<-cbind(W2,ColBl3)
      PW3<-cbind(W3,ColBl1) 
    }
    
    RowBl1<-matrix(" ",(4-dim(mi1)[1]),max(dim(dataB1)[1],dim(dataB2)[1],dim(dataF2)[1]))
    PPW1<-rbind(PW1,RowBl1)
    RowBl2<-matrix(" ",(4-dim(mi2)[1]),max(dim(dataB1)[1],dim(dataB2)[1],dim(dataF2)[1]))
    PPW2<-rbind(PW2,RowBl2)
    RowBl3<-matrix(" ",(9-dim(mi3)[1]),max(dim(dataB1)[1],dim(dataB2)[1],dim(dataF2)[1]))
    PPW3<-rbind(PW3,RowBl3)
    PPW<-t(rbind(PPW1,PPW2,PPW3))
    dimnames(PPW)<-list(c(1:max(dim(dataB1)[1],dim(dataB2)[1],dim(dataF2)[1])),c("B1:2(1)","B1:2(2)","B1:2(3)","B1:2(4)","B2:2(1)","B2:2(2)","B2:2(3)","B2:2(4)",
                                                                                 "F2:3(1)","F2:3(2)","F2:3(3)","F2:3(4)","F2:3(5)","F2:3(6)","F2:3(7)","F2:3(8)","F2:3(9)"))
    
  }
  
  return(PPW)
  
}

plotfun<-function(Population,Population2,data,result,bins,pr,colour){
  
  dplot<-function(x,m,sigma,mix_pi,label,ib,precision2,color){
    x2 <- seq(min(x), max(x), by = (max(x)-min(x))/10000)
    y2 <- dmixnorm(x2,m,sqrt(sigma),mix_pi)
    y<-matrix(0,dim(m)[1],length(x2))
    for(i in 1:dim(m)[1]){y[i,]<-mix_pi[i]*dnorm(x2,m[i],sqrt(sigma[i]))}
    par(mar=c(5,5,4,5)+0.1)
    bins <- seq(min(x), max(x), length.out = ib + 1)
    h<-hist(x,breaks = bins,xlab =label,main=" ")
    par(new=T)
    maxvalue<-round(max(y2,y),2)
    plot(x2,y2,axes=F,ylim=c(0,max(y2,y)),xlab="",ylab="",col=color,type='l',yaxt="n")
    axis(4,col=color,col.ticks=color,col.axis=color,yaxp=c(0,maxvalue,precision2))
    mtext("Density",side=4,line=3,col=color)
    a<-0
    while(a<dim(m)[1]){
      a<-a+1
      par(new=T)
      plot(x2,y[a,],axes=F,ylim=c(0,max(y2,y)),xlab="",ylab="",type='l',lty=2)
    }
  }  
  
  if(Population=="F2"){
    
    data<-sapply(data,as.character)
    dF2<-data[-1,which(data[1,]=="F2")];F2<-as.numeric(dF2[which(is.na(as.numeric(dF2))==FALSE)]);data<-as.matrix(F2) 
    
    datac<-result[[1]]  
    mmi<-datac[14:22]
    if(length(mmi[-which(mmi[]==" ")])==0)
    {mi<-as.matrix(as.numeric(mmi))}else{mi<-as.matrix(as.numeric(mmi[-which(mmi[]==" ")]))}
    
    mm<-datac[4:12]
    if(length(mm[-which(mm[]==" ")])==0)
    {m<-as.matrix(as.numeric(mm))}else{m<-as.matrix(as.numeric(mm[-which(mm[]==" ")]))}
    sigma<-matrix(as.numeric(datac[13]),dim(m)[1],1)
    
    
    dplot(data,m,sigma,mi,Population,bins,as.numeric(pr),colour)
    
  }else if(Population=="F2:3"){
    data<-sapply(data,as.character)
    dF23<-data[-1,which(data[1,]=="F23")];F23<-as.numeric(dF23[which(is.na(as.numeric(dF23))==FALSE)]);data<-as.matrix(F23) 
    
    datac<-result[[1]]
    mmi<-datac[22:30]
    if(length(mmi[-which(mmi[]==" ")])==0)
    {mi<-as.matrix(as.numeric(mmi))}else{mi<-as.matrix(as.numeric(mmi[-which(mmi[]==" ")]))} 
    mm<-datac[4:12]
    if(length(mm[-which(mm[]==" ")])==0)
    {m<-as.matrix(as.numeric(mm))}else{m<-as.matrix(as.numeric(mm[-which(mm[]==" ")]))}
    ssigma<-datac[13:21]
    if(length(ssigma[-which(ssigma[]==" ")])==0)
    {sigma<-as.matrix(as.numeric(ssigma))}else{sigma<-as.matrix(as.numeric(ssigma[-which(ssigma[]==" ")]))}
    
    
    dplot(data,m,sigma,mi,Population,bins,as.numeric(pr),colour)
    
    
  }else if(Population=="DH"){
    data<-sapply(data,as.character)
    dDH<-data[-1,which(data[1,]=="DH")];DH<-as.numeric(dDH[which(is.na(as.numeric(dDH))==FALSE)]);data<-as.matrix(DH) 
    
    datac<-result[[1]]
    mmi4<-datac[21:36]
    if(length(mmi4[-which(mmi4[]==" ")])==0)
    {mi4<-as.matrix(as.numeric(mmi4))}else{mi4<-as.matrix(as.numeric(mmi4[-which(mmi4[]==" ")]))}
    
    mm4<-datac[4:19]
    if(length(mm4[-which(mm4[]==" ")])==0)
    {m4<-as.matrix(as.numeric(mm4))}else{m4<-as.matrix(as.numeric(mm4[-which(mm4[]==" ")]))}  
    sigma4<-matrix(as.numeric(datac[20]),dim(m4)[1],1) 
    
    dplot(data,m4,sigma4,mi4,Population,bins,as.numeric(pr),colour)
    
  }else if(Population=="BIL"){
    data<-sapply(data,as.character)
    dBIL<-data[-1,which(data[1,]=="BIL")];BIL<-as.numeric(dBIL[which(is.na(as.numeric(dBIL))==FALSE)]);data<-as.matrix(BIL) 
    
    datac<-result[[1]]
    mmi<-datac[13:20]
    if(length(mmi[-which(mmi[]==" ")])==0)
    {mi<-as.matrix(as.numeric(mmi))}else{mi<-as.matrix(as.numeric(mmi[-which(mmi[]==" ")]))}
    
    mm<-datac[4:11]
    if(length(mm[-which(mm[]==" ")])==0)
    {m<-as.matrix(as.numeric(mm))}else{m<-as.matrix(as.numeric(mm[-which(mm[]==" ")]))}  
    sigma<-matrix(as.numeric(datac[12]),dim(m)[1],1) 
    
    dplot(data,m,sigma,mi,Population,bins,as.numeric(pr),colour)
    
  }else if(Population=="BC (B1 B2)"){
    
    data<-sapply(data,as.character)
    dB1<-data[-1,which(data[1,]=="B1")];B1<-as.numeric(dB1[which(is.na(as.numeric(dB1))==FALSE)]);dataB1<-as.matrix(B1)
    dB2<-data[-1,which(data[1,]=="B2")];B2<-as.numeric(dB2[which(is.na(as.numeric(dB2))==FALSE)]);dataB2<-as.matrix(B2)
    
    datac<-result[[1]]
    
    mmi4<-datac[9:12]
    if(length(mmi4[-which(mmi4[]==" ")])==0)
    {mi4<-as.matrix(as.numeric(mmi4))}else{mi4<-as.matrix(as.numeric(mmi4[-which(mmi4[]==" ")]))}
    mmi5<-datac[18:21]
    if(length(mmi5[-which(mmi5[]==" ")])==0)
    {mi5<-as.matrix(as.numeric(mmi5))}else{mi5<-as.matrix(as.numeric(mmi5[-which(mmi5[]==" ")]))}
    
    mm4<-datac[4:7]
    if(length(mm4[-which(mm4[]==" ")])==0)
    {m4<-as.matrix(as.numeric(mm4))}else{m4<-as.matrix(as.numeric(mm4[-which(mm4[]==" ")]))}  
    sigma4<-matrix(as.numeric(datac[8]),dim(m4)[1],1) 
    
    mm5<-datac[13:16]
    if(length(mm5[-which(mm5[]==" ")])==0)
    {m5<-as.matrix(as.numeric(mm5))}else{m5<-as.matrix(as.numeric(mm5[-which(mm5[]==" ")]))}  
    sigma5<-matrix(as.numeric(datac[17]),dim(m5)[1],1) 
    
    if(Population2=="B1"){
      dplot(dataB1,m4,sigma4,mi4,Population2,bins,as.numeric(pr),colour)  
    }else if(Population2=="B2"){
      dplot(dataB2,m5,sigma5,mi5,Population2,bins,as.numeric(pr),colour)
    }
    
  }else if(Population=="BCF (B1:2 B2:2)"){
    
    data<-sapply(data,as.character)
    dB12<-data[-1,which(data[1,]=="B12")];B12<-as.numeric(dB12[which(is.na(as.numeric(dB12))==FALSE)]);dataB1<-as.matrix(B12)
    dB22<-data[-1,which(data[1,]=="B22")];B22<-as.numeric(dB22[which(is.na(as.numeric(dB22))==FALSE)]);dataB2<-as.matrix(B22)
    
    datac<-result[[1]]
    
    mmi4<-datac[12:15]
    if(length(mmi4[-which(mmi4[]==" ")])==0)
    {mi4<-as.matrix(as.numeric(mmi4))}else{mi4<-as.matrix(as.numeric(mmi4[-which(mmi4[]==" ")]))}
    
    mmi5<-datac[24:27]
    if(length(mmi5[-which(mmi5[]==" ")])==0)
    {mi5<-as.matrix(as.numeric(mmi5))}else{mi5<-as.matrix(as.numeric(mmi5[-which(mmi5[]==" ")]))}
    
    mm4<-datac[4:7]
    if(length(mm4[-which(mm4[]==" ")])==0)
    {m4<-as.matrix(as.numeric(mm4))}else{m4<-as.matrix(as.numeric(mm4[-which(mm4[]==" ")]))}  
    ssigma4<-datac[8:11]
    if(length(ssigma4[-which(ssigma4[]==" ")])==0)
    {sigma4<-as.matrix(as.numeric(ssigma4))}else{sigma4<-as.matrix(as.numeric(ssigma4[-which(ssigma4[]==" ")]))}  
    
    mm5<-datac[16:19]
    if(length(mm5[-which(mm5[]==" ")])==0)
    {m5<-as.matrix(as.numeric(mm5))}else{m5<-as.matrix(as.numeric(mm5[-which(mm5[]==" ")]))}  
    ssigma5<-datac[20:23]
    if(length(ssigma5[-which(ssigma5[]==" ")])==0)
    {sigma5<-as.matrix(as.numeric(ssigma5))}else{sigma5<-as.matrix(as.numeric(ssigma5[-which(ssigma5[]==" ")]))}
    
    if(Population2=="B1:2"){
      dplot(dataB1,m4,sigma4,mi4,Population2,bins,as.numeric(pr),colour)   
    }else if(Population2=="B2:2"){
      dplot(dataB2,m5,sigma5,mi5,Population2,bins,as.numeric(pr),colour)  
    }
    
  }else if(Population=="G4F2 (P1 P2 F1 F2)"){
    
    data<-sapply(data,as.character)
    dF2<-data[-1,which(data[1,]=="F2")];F2<-as.numeric(dF2[which(is.na(as.numeric(dF2))==FALSE)]);dataF2<-as.matrix(F2)
    
    datac<-result[[1]]
    
    mmi<-datac[18:26]
    if(length(mmi[-which(mmi[]==" ")])==0)
    {mi<-as.matrix(as.numeric(mmi))}else{mi<-as.matrix(as.numeric(mmi[-which(mmi[]==" ")]))}
    
    mm<-datac[7:15]
    if(length(mm[-which(mm[]==" ")])==0)
    {m<-as.matrix(as.numeric(mm))}else{m<-as.matrix(as.numeric(mm[-which(mm[]==" ")]))}  
    sigma<-matrix(as.numeric(datac[16]),dim(m)[1],1) 
    
    dplot(dataF2,m,sigma,mi,Population,bins,as.numeric(pr),colour)  
    
  }else if(Population=="G4F3 (P1 P2 F1 F2:3)"){
    data<-sapply(data,as.character)
    dF23<-data[-1,which(data[1,]=="F23")];F23<-as.numeric(dF23[which(is.na(as.numeric(dF23))==FALSE)]);dataF3<-as.matrix(F23)
    
    
    datac<-result[[1]] 
    
    mmi<-datac[26:34]
    if(length(mmi[-which(mmi[]==" ")])==0)
    {mi<-as.matrix(as.numeric(mmi))}else{mi<-as.matrix(as.numeric(mmi[-which(mmi[]==" ")]))}
    
    mm<-datac[7:15]
    if(length(mm[-which(mm[]==" ")])==0)
    {m<-as.matrix(as.numeric(mm))}else{m<-as.matrix(as.numeric(mm[-which(mm[]==" ")]))}  
    
    ssigma<-datac[17:25]
    if(length(ssigma[-which(ssigma[]==" ")])==0)
    {sigma<-as.matrix(as.numeric(ssigma))}else{sigma<-as.matrix(as.numeric(ssigma[-which(ssigma[]==" ")]))}  
    
    dplot(dataF3,m,sigma,mi,Population,bins,as.numeric(pr),colour)  
    
    
  }else if(Population=="G3DH (P1 P2 DH)"){
    
    data<-sapply(data,as.character)
    dDH<-data[-1,which(data[1,]=="DH")];DH<-as.numeric(dDH[which(is.na(as.numeric(dDH))==FALSE)]);dataDH<-as.matrix(DH)
    
    datac<-result[[1]] 
    
    mmi4<-datac[24:39]
    if(length(mmi4[-which(mmi4[]==" ")])==0)
    {mi4<-as.matrix(as.numeric(mmi4))}else{mi4<-as.matrix(as.numeric(mmi4[-which(mmi4[]==" ")]))}
    
    mm4<-datac[7:22]
    if(length(mm4[-which(mm4[]==" ")])==0)
    {m4<-as.matrix(as.numeric(mm4))}else{m4<-as.matrix(as.numeric(mm4[-which(mm4[]==" ")]))}  
    sigma4<-matrix(as.numeric(datac[23]),dim(m4)[1],1) 
    
    dplot(dataDH,m4,sigma4,mi4,Population,bins,as.numeric(pr),colour)  
    
  }else if(Population=="G5BC (P1 P2 F1 B1 B2)"){
    
    data<-sapply(data,as.character)
    dB1<-data[-1,which(data[1,]=="B1")];B1<-as.numeric(dB1[which(is.na(as.numeric(dB1))==FALSE)]);dataB1<-as.matrix(B1)
    dB2<-data[-1,which(data[1,]=="B2")];B2<-as.numeric(dB2[which(is.na(as.numeric(dB2))==FALSE)]);dataB2<-as.matrix(B2)
    
    datac<-result[[1]] 
    
    mmi4<-datac[13:16]
    if(length(mmi4[-which(mmi4[]==" ")])==0)
    {mi4<-as.matrix(as.numeric(mmi4))}else{mi4<-as.matrix(as.numeric(mmi4[-which(mmi4[]==" ")]))}
    
    mmi5<-datac[22:25]
    if(length(mmi5[-which(mmi5[]==" ")])==0)
    {mi5<-as.matrix(as.numeric(mmi5))}else{mi5<-as.matrix(as.numeric(mmi5[-which(mmi5[]==" ")]))}
    
    
    mm4<-datac[8:11]
    if(length(mm4[-which(mm4[]==" ")])==0)
    {m4<-as.matrix(as.numeric(mm4))}else{m4<-as.matrix(as.numeric(mm4[-which(mm4[]==" ")]))}  
    sigma4<-matrix(as.numeric(datac[12]),dim(m4)[1],1) 
    mm5<-datac[17:20]
    if(length(mm5[-which(mm5[]==" ")])==0)
    {m5<-as.matrix(as.numeric(mm5))}else{m5<-as.matrix(as.numeric(mm5[-which(mm5[]==" ")]))}  
    sigma5<-matrix(as.numeric(datac[21]),dim(m5)[1],1) 
    
    if(Population2=="B1"){
      dplot(dataB1,m4,sigma4,mi4,Population2,bins,as.numeric(pr),colour)   
      
    }else if(Population2=="B2"){
      dplot(dataB2,m5,sigma5,mi5,Population2,bins,as.numeric(pr),colour)  
      
    }
    
  }else if(Population=="G5BCF (P1 P2 F1 B1:2 B2:2)"){
    data<-sapply(data,as.character)
    dB12<-data[-1,which(data[1,]=="B12")];B12<-as.numeric(dB12[which(is.na(as.numeric(dB12))==FALSE)]);dataB1<-as.matrix(B12)
    dB22<-data[-1,which(data[1,]=="B22")];B22<-as.numeric(dB22[which(is.na(as.numeric(dB22))==FALSE)]);dataB2<-as.matrix(B22)
    
    datac<-result[[1]] 
    
    mmmi1<-datac[15:18]
    if(length(mmmi1[-which(mmmi1[]==" ")])==0)
    {mmi1<-as.matrix(as.numeric(mmmi1))}else{mmi1<-as.matrix(as.numeric(mmmi1[-which(mmmi1[]==" ")]))}
    
    mmmi2<-datac[27:30]
    if(length(mmmi2[-which(mmmi2[]==" ")])==0)
    {mmi2<-as.matrix(as.numeric(mmmi2))}else{mmi2<-as.matrix(as.numeric(mmmi2[-which(mmmi2[]==" ")]))}
    
    mm1<-datac[7:10]
    if(length(mm1[-which(mm1[]==" ")])==0)
    {m1<-as.matrix(as.numeric(mm1))}else{m1<-as.matrix(as.numeric(mm1[-which(mm1[]==" ")]))}  
    ssigma1<-datac[11:14]
    if(length(ssigma1[-which(ssigma1[]==" ")])==0)
    {sigma1<-as.matrix(as.numeric(ssigma1))}else{sigma1<-as.matrix(as.numeric(ssigma1[-which(ssigma1[]==" ")]))} 
    mm2<-datac[19:22]
    if(length(mm2[-which(mm2[]==" ")])==0)
    {m2<-as.matrix(as.numeric(mm2))}else{m2<-as.matrix(as.numeric(mm2[-which(mm2[]==" ")]))}  
    ssigma2<-datac[23:26]
    if(length(ssigma2[-which(ssigma2[]==" ")])==0)
    {sigma2<-as.matrix(as.numeric(ssigma2))}else{sigma2<-as.matrix(as.numeric(ssigma2[-which(ssigma2[]==" ")]))} 
    
    if(Population2=="B1:2"){
      dplot(dataB1,m1,sigma1,mmi1,Population2,bins,as.numeric(pr),colour)   
      
    }else if(Population2=="B2:2"){
      dplot(dataB2,m2,sigma2,mmi2,Population2,bins,as.numeric(pr),colour)  
    }
  }else if(Population=="G5 (P1 P2 F1 F2 F2:3)"){
    data<-sapply(data,as.character)
    dF2<-data[-1,which(data[1,]=="F2")];F2<-as.numeric(dF2[which(is.na(as.numeric(dF2))==FALSE)]);dataF2<-as.matrix(F2)
    dF23<-data[-1,which(data[1,]=="F23")];F23<-as.numeric(dF23[which(is.na(as.numeric(dF23))==FALSE)]);dataF3<-as.matrix(F23)
    
    datac<-result[[1]]
    
    mmi4<-datac[17:25]
    if(length(mmi4[-which(mmi4[]==" ")])==0)
    {mi4<-as.matrix(as.numeric(mmi4))}else{mi4<-as.matrix(as.numeric(mmi4[-which(mmi4[]==" ")]))}
    
    mmi5<-datac[44:52]
    if(length(mmi5[-which(mmi5[]==" ")])==0)
    {mi5<-as.matrix(as.numeric(mmi5))}else{mi5<-as.matrix(as.numeric(mmi5[-which(mmi5[]==" ")]))}
    
    mm4<-datac[7:15]
    if(length(mm4[-which(mm4[]==" ")])==0)
    {m4<-as.matrix(as.numeric(mm4))}else{m4<-as.matrix(as.numeric(mm4[-which(mm4[]==" ")]))}  
    sigma4<-matrix(as.numeric(datac[16]),dim(m4)[1],1) 
    mm5<-datac[26:34]
    if(length(mm5[-which(mm5[]==" ")])==0)
    {m5<-as.matrix(as.numeric(mm5))}else{m5<-as.matrix(as.numeric(mm5[-which(mm5[]==" ")]))}  
    ssigma5<-datac[35:43]
    if(length(ssigma5[-which(ssigma5[]==" ")])==0)
    {sigma5<-as.matrix(as.numeric(ssigma5))}else{sigma5<-as.matrix(as.numeric(ssigma5[-which(ssigma5[]==" ")]))}  
    
    if(Population2=="F2"){
      dplot(dataF2,m4,sigma4,mi4,Population2,bins,as.numeric(pr),colour)   
    }else if(Population2=="F2:3"){
      dplot(dataF3,m5,sigma5,mi5,Population2,bins,as.numeric(pr),colour)  
    }
    
  }else if(Population=="G6 (P1 P2 F1 F2 B1 B2)"){
    
    data<-sapply(data,as.character)
    dB1<-data[-1,which(data[1,]=="B1")];B1<-as.numeric(dB1[which(is.na(as.numeric(dB1))==FALSE)]);dataB1<-as.matrix(B1)
    dB2<-data[-1,which(data[1,]=="B2")];B2<-as.numeric(dB2[which(is.na(as.numeric(dB2))==FALSE)]);dataB2<-as.matrix(B2)
    dF2<-data[-1,which(data[1,]=="F2")];F2<-as.numeric(dF2[which(is.na(as.numeric(dF2))==FALSE)]);dataF2<-as.matrix(F2)
    
    datac<-result[[1]]
    
    mmi4<-datac[13:16]
    if(length(mmi4[-which(mmi4[]==" ")])==0)
    {mi4<-as.matrix(as.numeric(mmi4))}else{mi4<-as.matrix(as.numeric(mmi4[-which(mmi4[]==" ")]))}
    
    mmi5<-datac[22:25]
    if(length(mmi5[-which(mmi5[]==" ")])==0)
    {mi5<-as.matrix(as.numeric(mmi5))}else{mi5<-as.matrix(as.numeric(mmi5[-which(mmi5[]==" ")]))}
    
    mmi6<-datac[36:44]
    if(length(mmi6[-which(mmi6[]==" ")])==0)
    {mi6<-as.matrix(as.numeric(mmi6))}else{mi6<-as.matrix(as.numeric(mmi6[-which(mmi6[]==" ")]))}
    
    mm4<-datac[8:11]
    if(length(mm4[-which(mm4[]==" ")])==0)
    {m4<-as.matrix(as.numeric(mm4))}else{m4<-as.matrix(as.numeric(mm4[-which(mm4[]==" ")]))}  
    sigma4<-matrix(as.numeric(datac[12]),dim(m4)[1],1) 
    
    mm5<-datac[17:20]
    if(length(mm5[-which(mm5[]==" ")])==0)
    {m5<-as.matrix(as.numeric(mm5))}else{m5<-as.matrix(as.numeric(mm5[-which(mm5[]==" ")]))}  
    sigma5<-matrix(as.numeric(datac[21]),dim(m5)[1],1) 
    
    mm6<-datac[26:34]
    if(length(mm6[-which(mm6[]==" ")])==0)
    {m6<-as.matrix(as.numeric(mm6))}else{m6<-as.matrix(as.numeric(mm6[-which(mm6[]==" ")]))}  
    sigma6<-matrix(as.numeric(datac[35]),dim(m6)[1],1) 
    
    if(Population2=="B1"){
      dplot(dataB1,m4,sigma4,mi4,Population2,bins,as.numeric(pr),colour)   
    }else if(Population2=="B2"){
      dplot(dataB2,m5,sigma5,mi5,Population2,bins,as.numeric(pr),colour)  
    }else if(Population2=="F2"){
      dplot(dataF2,m6,sigma6,mi6,Population2,bins,as.numeric(pr),colour)  
    }
    
    
  }else if(Population=="G6F (P1 F1 P2 B1:2 B2:2 F2:3)"){
    data<-sapply(data,as.character)
    dB12<-data[-1,which(data[1,]=="B12")];B12<-as.numeric(dB12[which(is.na(as.numeric(dB12))==FALSE)]);dataB1<-as.matrix(B12)
    dB22<-data[-1,which(data[1,]=="B22")];B22<-as.numeric(dB22[which(is.na(as.numeric(dB22))==FALSE)]);dataB2<-as.matrix(B22)
    dF23<-data[-1,which(data[1,]=="F23")];F23<-as.numeric(dF23[which(is.na(as.numeric(dF23))==FALSE)]);dataF2<-as.matrix(F23)
    
    datac<-result[[1]]
    
    mmi1<-datac[16:19]
    if(length(mmi1[-which(mmi1[]==" ")])==0)
    {mi1<-as.matrix(as.numeric(mmi1))}else{mi1<-as.matrix(as.numeric(mmi1[-which(mmi1[]==" ")]))}
    mmi2<-datac[28:31]
    if(length(mmi2[-which(mmi2[]==" ")])==0)
    {mi2<-as.matrix(as.numeric(mmi2))}else{mi2<-as.matrix(as.numeric(mmi2[-which(mmi2[]==" ")]))}
    mmi3<-datac[50:58]
    if(length(mmi3[-which(mmi3[]==" ")])==0)
    {mi3<-as.matrix(as.numeric(mmi3))}else{mi3<-as.matrix(as.numeric(mmi3[-which(mmi3[]==" ")]))}
    
    mm1<-datac[8:11]
    if(length(mm1[-which(mm1[]==" ")])==0)
    {m1<-as.matrix(as.numeric(mm1))}else{m1<-as.matrix(as.numeric(mm1[-which(mm1[]==" ")]))}  
    ssigma1<-datac[12:15]
    if(length(ssigma1[-which(ssigma1[]==" ")])==0)
    {sigma1<-as.matrix(as.numeric(ssigma1))}else{sigma1<-as.matrix(as.numeric(ssigma1[-which(ssigma1[]==" ")]))}  
    mm2<-datac[20:23]
    if(length(mm2[-which(mm2[]==" ")])==0)
    {m2<-as.matrix(as.numeric(mm2))}else{m2<-as.matrix(as.numeric(mm2[-which(mm2[]==" ")]))}  
    ssigma2<-datac[24:27]
    if(length(ssigma2[-which(ssigma2[]==" ")])==0)
    {sigma2<-as.matrix(as.numeric(ssigma2))}else{sigma2<-as.matrix(as.numeric(ssigma2[-which(ssigma2[]==" ")]))}  
    mm3<-datac[32:40]
    if(length(mm3[-which(mm3[]==" ")])==0)
    {m3<-as.matrix(as.numeric(mm3))}else{m3<-as.matrix(as.numeric(mm3[-which(mm3[]==" ")]))}  
    ssigma3<-datac[41:49]
    if(length(ssigma3[-which(ssigma3[]==" ")])==0)
    {sigma3<-as.matrix(as.numeric(ssigma3))}else{sigma3<-as.matrix(as.numeric(ssigma3[-which(ssigma3[]==" ")]))} 
    
    if(Population2=="B1:2"){
      dplot(dataB1,m1,sigma1,mi1,Population2,bins,as.numeric(pr),colour)   
    }else if(Population2=="B2:2"){
      dplot(dataB2,m2,sigma2,mi2,Population2,bins,as.numeric(pr),colour)  
    }else if(Population2=="F2:3"){
      dplot(dataF2,m3,sigma3,mi3,Population2,bins,as.numeric(pr),colour)  
    }
  }  
  
}

