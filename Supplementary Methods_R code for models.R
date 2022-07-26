# R code for models in paper by Irith Aloni & Amiyaal Ilany

# Content:
# Female only model - Line 15
# One mutation - Scenario I - Line 153
# One mutation - Scenario II - Line 463
# One mutation - Scenario III - Line 778
# One mutation - Scenario IV - Line 1086
# Two mutations - Scenario I - Line 1403
# Two mutations - Scenario II - Line 1740
# Two mutations - Scenario III - Line 2082
# Two mutations - Scenario IV - Line 2426
#**********************************************************************************

# Female only Model

library(Matrix)

update.network <- function(Kids,Rank,Fouder,Name,NameIn,Names) { 
  
  # STEP 1: one individual dies
  
  Die=N+Rank  #Death probability weight
  del     <- sample(1:N, 1, prob=Die)
  ndel=Name[del]
  Name=Name[-del]
  Kids=Kids[-del]
  Founder=Founder[-del]
  
  up=which(Rank>Rank[del])
  Rank[up]=Rank[up]-1
  Rank=Rank[-del]
  
  Chance=1.5*N-Rank #Reproduction probability weight
  
  
  # STEP 2: Select a female that reproduces 
  
  repro=sample(c(1:(N-1)),1, prob=Chance)
  Kids[repro]=Kids[repro]+1
  
  # Add newborn to variables
  
  down=which(Rank>Rank[repro])
  Rank[down]=Rank[down]+1
  Rank[N]=Rank[repro]+1
  
  Founder[N]=Founder[repro]
  Kids[N]=0
  NameIn=NameIn+1
  Name[N]=NameIn
  Names[zz+N,1:3]=c(NameIn,Name[repro],Founder[N])
  Names[Name[repro],"Kids"]=Kids[repro]
  Names[Name[repro],6+Kids[repro]]=Name[N]
  Names[zz+N,5]=Rank[N]
  Names[zz+N,6]=ndel
  
  return(list(Kids=Kids, Founder=Founder, Rank=Rank, Name=Name,  NameIn=NameIn, Names=Names))
}


# SIMULATION

# STEP 1: set initial parameters
N             <- 50
n.rep         <- 200 # iterations
burn.in       <- 500 # simulations

Foundout=matrix(0,n.rep,8)
colnames(Foundout)=c("Simulation","No.Founders","MainF","MainFreq","2ndF","2ndFrerq","3rdF","3rdFreq")

for (k in 1:n.rep) {                     
  
  Rank=c(1:N)
  Kids=rep(0,N)
  Founder=c(1:N)
  
  Name=c(1:N)
  NameIn=N
  
  A=matrix(0,burn.in+N,6)
  colnames(A)=c("Name","Mother","Founder","Kids","Initial rank","Died")
  B=matrix(0,burn.in+N,30)
  colnames(B)=paste0("Kid ",c(1:30))
  Names=cbind(A,B)
  
  Names[1:N,1]=Name
  Names[1:N,"Founder"]=Founder
  Names[1:N,"Initial rank"]=Founder
  
  # Start iterating
  for (zz in 1:burn.in) {
    output       <- update.network(Kids, Rank, Founder, Name,NameIn,Names)
    
    Kids=output$Kids
    Founder=output$Founder
    Rank=output$Rank
    Name=output$Name
    NameIn=output$NameIn
    Names=output$Names
  }
  
  # checking kid reproduction
  
  kidsRep=rep(0,N)
  
  for (i in 1:N) {
    KR=0
    if (Kids[i]>0)  {
      
      for (j in 1:(Kids[i]))  {
        if (Names[Names[Name[i],j+6],"Kids"]>0) {KR=KR+1}
      }
    }
    kidsRep[i]=KR
  }
  
  
  netout=matrix(0,N,7)
  colnames(netout)=c("Ind","Founder","Rank","N.kids","N.kids.repro","Initial.Rank","simulation")
  
  netout[,1]=c(1:N)
  netout[,2]=Founder
  netout[,3]=Rank
  netout[,4]=Kids
  netout[,5]=kidsRep
  netout[,6]=Names[Name,"Initial rank"]
  netout[,7]=k
  
  write.table(netout,"netOutput.csv",col.names=TRUE,row.names=F,sep=",",append=T)
  write.table(netout[1:40,],"netOutput40.csv",col.names=TRUE,row.names=F,sep=",",append=T)
  
  Foundout[k,1]=k
  Foundout[k,2]=length(unique(Founder))
  
  val=rle(sort(Founder))$values
  Freq=rle(sort(Founder))$lengths
  
  Foundout[k,3]=val[1]
  Foundout[k,4]=Freq[1]/N
  Foundout[k,5]=val[2]
  Foundout[k,6]=Freq[2]/N
  Foundout[k,7]=val[3]
  Foundout[k,8]=Freq[3]/N
  
}#end of k

write.table(Foundout,"Founder.csv",col.names=TRUE,row.names=F,sep=",", append=T)

#****************************************************************************************


# One mutation, Scenario I


library(Matrix)

update.network <- function(datin) {
  
  for (j in 1:length(datin)){
    assign(datanames[j],datin[[j]])
  }
  
  SexNew=1
  nameNewM=0
  deadClan=FALSE
  
  if(N>Male & Male>0) {
    
    # One individual dies
    
    del=ifelse(zz==1, sample(Name[-mutant],1),sample(1:N, 1, prob=Die))
    ndel=Name[del]
    Name=Name[-del]
    Die=Die[-del]
    
    if (Sex[del]==1) {
      if (N-Male==1) {
        N=Male
        Founder=NULL
        Kids=NULL
        Rank=NULL
        GenF=NULL
        deadClan=TRUE  
      }
      
      if (N-Male>1) {
        Founder=Founder[-del]
        Kids=Kids[-del]
        up=which(Rank>Rank[del])
        Rank[up]=Rank[up]-1
        Rank=Rank[-del]
        GenF=GenF[-del,]
      }
    }
    
    
    if(Sex[del]==0){
      if(Male>1) {GenM=GenM[-(del-N+Male),]}
      else {GenM=NULL}
      Male=Male-1
      if (Male==0) {N=N-1}
    }
    Sex=Sex[-del]
    
    if (deadClan==FALSE & Male>0) {
      
      Chance=N-Rank
      
      # Select male & female for reproduction
      
      repro=ifelse(zz==1,which(Name==mutant),sample(c(1:(N-1-Male)),1, prob=Chance))
      mate=sample( c(1:Male),1)
      Pa=Name[N-Male-1+mate]
      Kids[repro]=Kids[repro]+1
      
      NameIn=NameIn+1
      
      if (zz==1) {
        GenF[repro,]=c(1,0)
        Die[repro]=1
        Names[Name[repro],"Phen"]=1 #phenotype
      }
      
      # set variables for newborn
      
      GeN=c(0,0)
      GeN[1]=ifelse(length(GenF)>2,sample(GenF[repro,],1),sample(GenF,1))
      GeN[2]=ifelse(length(GenM)>2,sample(GenM[mate,],1),sample(GenM,1))
      SexNew=sample(c(0,1),1)
      
      if (SexNew==0) {
        N=N-1
        nameNewM=paste0(NameIn,clan)
        Names[NameIn,"Name"]=nameNewM
      }
      
      if(SexNew==1) {
        
        Name[(N-Male+1):N]=Name[(N-Male):(N-1)]
        Name[N-Male]=NameIn
        Names[NameIn,"Name"]=NameIn
        Sex[N-Male]=1
        if (Male>0) {Sex[N]=0}
        down=which(Rank>Rank[repro])
        Rank[down]=Rank[down]+1
        Rank[N-Male]=Rank[repro]+1
        
        Founder[N-Male]=Founder[repro]
        Kids[N-Male]=0
        GenF=rbind(GenF,GeN)
        rownames(GenF) <- NULL
        
        Die[(N-Male+1):N]=Die[(N-Male):(N-1)]
        Die[N-Male]=ifelse(sum(GeN)>0,1,2)
        
      }
      
      Names[NameIn,2:4]=c(SexNew,Name[repro],Founder[repro])
      Names[as.numeric(Name[repro]),"Kids"]=Kids[repro]
      if(Kids[repro]<31) {Names[as.numeric(Name[repro]),8+Kids[repro]]=NameIn}
      if(SexNew==1) {Names[NameIn,"Initial rank"]=Rank[N-Male]}
      Names[NameIn,"Died"]=ndel
      Names[NameIn, "Father"]=Pa
      Names[NameIn,"Phen"]=sum(GeN) 
      
    }}
  
  
  datout=list(N,Sex, Male, Die, Kids, Rank, Founder, Name,NameIn, Names ,GenM,GenF)
  
  return(list(datout=datout, SexNew=SexNew,nameNewM=nameNewM,GeN=GeN, deadClan=deadClan))
  
}

# In case the newborn is a male - add to new clan with the following function

update.male=function(datin, nameNewM, GeN) {
  
  for (j in 1:length(datin)){
    assign(datanames[j],datin[[j]])
  }
  
  N=N+1
  Male=Male+1
  Name[N]=nameNewM
  Sex[N]=0
  GenM=rbind(GenM,GeN)
  rownames(GenM)=NULL
  Die[N]=ifelse(sum(GeN)==0,2,1)
  
  datout=list(N,Sex, Male, Die, Kids, Rank, Founder, Name,NameIn, Names ,GenM,GenF)
  
  return(list(datMout=datout))
}


# SIMULATION

# Set initial parameters
N      <- 60
n.rep         <-100   # of simulations (for each mutant's rank)
burn.in       <-3000  # of iterations


out=matrix(0,n.rep*30,26)
colnames(out)=c("Mutant.Rank", "Mutant clan","NA","MaleA", "No.Mmutants A", "N.Fmutants A","NB",
                "MaleB", "No.Mmutants B","N.Fmutants B", "NC","MaleC", "No.Mmutants C","N.Fmutants C",
                "ND","MaleD", "No.Mmutants D", "N.Fmutants D", "NE", "MaleE","No.Mmutants E", 
                "N.Fmutants E", "Total N.iterattion", "simulation","No.Mutants","Prop.Mutants" )

# Start simulations

for(mutant in 1:30) {
  for (k in 1:n.rep) {
    
    N=60
    Male=30
    Sex=c(rep(1,N-Male),rep(0,Male)) #1=female 0=male
    Rank=c(1:(N-Male))
    
    # set Kids, Die, founder, Name, NameIndex, Names, Gen
    
    Die=rep(2,N)
    Kids=rep(0,(N-Male))
    Founder=c(1:(N-Male))
    
    GenF=matrix(0,N-Male,2)
    GenM=matrix(0,Male,2)
    
    Name=c(1:N)
    NameIn=N
    
    A=matrix(0,burn.in+N,8)
    colnames(A)=c("Name","Sex","Mother","Founder","Father","Kids","Initial rank","Died")
    B=matrix(0,burn.in+N,30)
    colnames(B)=paste0("Kid ",c(1:30))
    Phen=rep(0,burn.in+N)
    Names=cbind(A,B,Phen)
    
    Names[1:N,1]=Name
    Names[1:N,"Sex"]=Sex
    Names[1:(N-Male),"Founder"]=Founder
    Names[1:(N-Male),"Initial rank"]=Rank
    
    # Set clans 
    
    Clans= c("A","B","C","D","E")
    
    
    datlist=list(N,Sex, Male, Die, Kids, Rank, Founder, Name,NameIn,Names ,GenM,GenF)
    
    names(datlist)=c("Nt","Sex", "Male", "Die", "Kids", "Rank", "Founder", "Name","NameIn","Names","GenM","GenF")
    datnames=names(datlist)
    
    datanames=datnames  #for passing to function
    datanames[1]="N"
    
    
    alldat=list()
    for (i in Clans) {
      ldat=list()
      n=1
      for (j in datnames) {
        ldat[[n]]=assign(paste0(j,i),datlist[[j]])
        n=n+1
      }
      names(ldat)=paste0(datnames,i)
      assign(paste0("dat",i),ldat)
    }
    alldat=list(datA=datA,datB=datB,datC=datC,datD=datD,datE=datE)  
    
    
    # Start Iterations
    
    for (zz in 1:burn.in) {
      
      clan=sample(Clans,1)
      dat=alldat[[which(Clans==clan)]] 
      datSnames=names(dat)
      
      if(zz==1) {clanMut=clan}    
      
      output<-update.network(dat)
      
      dat=output$datout
      names(dat)=datSnames
      assign(paste0("dat",clan),dat)
      alldat[[which(Clans==clan)]]=dat
      names(alldat[[which(Clans==clan)]])=datSnames
      
      SexNew=output$SexNew
      nameNewM=output$nameNewM
      GeN=output$GeN
      deadClan=output$deadClan
      
      
      if (SexNew==0) {
        
        clanM=sample(Clans[-which(Clans==clan)],1)
        dat=alldat[[which(Clans==clanM)]]
        datSnames=names(dat)
        
        outMale=update.male(dat,nameNewM, GeN)
        
        dat=outMale$datMout
        names(dat)=datSnames
        assign(paste0("dat",clanM),dat)
        
        alldat[[which(Clans==clanM)]]=dat
        names(alldat[[which(Clans==clanM)]])=datSnames
      }
      
      
      if (deadClan==TRUE) {
        clandel=which(Clans==clan)
        alldat[[clandel]]=NULL
        Clans=Clans[-which(Clans==clan)]
      }
    } #end of a simulation
    
    # output from simulation
    
    alldat=list(datA=datA,datB=datB,datC=datC,datD=datD,datE=datE)
    
    out[n.rep*(mutant-1)+k,1]=mutant
    out[n.rep*(mutant-1)+k,2]=clanMut
    totMut=0
    
    for (i in 1:5) {
      out[n.rep*(mutant-1)+k,3+4*(i-1)]=alldat[[i]]$Nt
      out[n.rep*(mutant-1)+k,4*i]=alldat[[i]]$Male
      if(alldat[[i]]$Male>1) {
        out[n.rep*(mutant-1)+k,4*i+1]=length(which(rowSums(alldat[[i]][[11]])>0))
        totMut=totMut+length(which(rowSums(alldat[[i]][[11]])>0))
      }
      if(alldat[[i]]$Male==1){
        out[k,4*i+1]=ifelse(sum(alldat[[i]][[11]])>0,1,0)
        totMut=totMut+ceiling(sum(alldat[[i]][[11]])/2)
      }
      if(alldat[[i]]$Nt>alldat[[i]]$Male+1) { 
        out[n.rep*(mutant-1)+k,4*i+2]=length(which(rowSums(alldat[[i]][[12]])>0))
        totMut=totMut+length(which(rowSums(alldat[[i]][[12]])>0))
      }
      if(alldat[[i]]$Nt-alldat[[i]]$Male+1==1){
        out[n.rep*(mutant-1)+k,4*i+2]=ifelse(sum(alldat[[i]][[12]])>0,1,0)
        totMut=totMut+ceiling(sum(alldat[[i]][[12]])/2)      }
    }
    
    out[n.rep*(mutant-1)+k,23]=burn.in
    out[n.rep*(mutant-1)+k,24]=k 
    out[n.rep*(mutant-1)+k,25]=totMut
    out[n.rep*(mutant-1)+k,26]=totMut/(5*N)
    
  } #end of k
} #end of mutant


write.table(out,"mutation Sum.csv",col.names=TRUE, row.names=F,sep=",",append=T)

#***************************************************************************************

# One mutation - Scenario II


library(Matrix)

update.network <- function(datin) {
  
  for (j in 1:length(datin)){
    assign(datanames[j],datin[[j]])
  }
  
  SexNew=1
  nameNewM=0
  deadClan=FALSE
  
  if(N>Male & Male>0) {
    
    # One individual dies
    
    del=ifelse(zz==1,sample(c((mutant+1):N),1), sample(1:N, 1, prob=Die))
    ndel=Name[del]
    Name=Name[-del]
    Die=Die[-del]
    up=which(Rank>Rank[del])
    Rank[up]=Rank[up]-1
    Rank=Rank[-del]
    
    
    if (Sex[del]==1) {
      if (N-Male==1) {
        N=Male
        Founder=NULL
        Kids=NULL
        GenF=NULL
        deadClan=TRUE  
      }
      
      if (N-Male>1) {
        Founder=Founder[-del]
        Kids=Kids[-del]
        GenF=GenF[-del,]
      }
    }
    
    if(Sex[del]==0){
      if(Male>1) {GenM=GenM[-(del-N+Male),]}
      else {GenM=NULL}
      Male=Male-1
      if (Male==0) {N=N-1}
    }
    Sex=Sex[-del]
    
    if (deadClan==FALSE & Male>0) {
      
      Chance=N-Rank[1:(N-1-Male)]
      
      # Select female & male for reproduction
      
      repro=ifelse(zz==1,mutant,sample(c(1:(N-1-Male)),1, prob=Chance)) 
      mate=sample( c(1:Male),1)
      Pa=Name[N-Male-1+mate]
      Kids[repro]=Kids[repro]+1
      
      NameIn=NameIn+1
      
      if (zz==1) {
        GenF[repro,]=c(1,0)
        Names[repro,"Phen"]=1 #phenotype
      }
      
      # Set variables for newborn
      
      GeN=c(0,0)
      GeN[1]=ifelse(length(GenF)>2,sample(GenF[repro,],1),sample(GenF,1))
      GeN[2]=ifelse(length(GenM)>2,sample(GenM[mate,],1),sample(GenM,1))
      SexNew=sample(c(0,1),1)
      
      if (SexNew==0) {
        N=N-1
        nameNewM=paste0(NameIn,clan)
        Names[NameIn,"Name"]=nameNewM
      }
      
      if(SexNew==1) {
        
        Name[(N-Male+1):N]=Name[(N-Male):(N-1)]
        Name[N-Male]=NameIn
        Names[NameIn,"Name"]=NameIn
        Sex[N-Male]=1
        if (Male>0) {Sex[N]=0}
        down=which(Rank>Rank[repro])
        Rank[down]=Rank[down]+1
        Rank[(N-Male+1):(N)]=Rank[(N-Male):(N-1)]
        Rank[N-Male]=Rank[repro]+1
        
        Founder[N-Male]=Founder[repro]
        Kids[N-Male]=0
        
        GenF=rbind(GenF,GeN)
        rownames(GenF) <- NULL
        
      }
      
      if(N-Male>1) {
        for (j in 1:(N-Male)) {
          Die[j]=ifelse(sum(GenF[j,])>0,0.5*(N+Rank[j]),N+Rank[j])
        } }
      if(N-Male==1) {Die[1]=ifelse(sum(GenF)>0,0.5*(N+1),N+1)}
      
      
      if(Male>1) {
        for (j in (N-Male+1):N) {
          Die[j]=ifelse(sum(GenM[j-N+Male,])>0,0.5*(N+Rank[j]),N+Rank[j])
        }}
      if(Male==1) {Die[N]=ifelse(sum(GenM)>0,N,2*N)}
      
      Names[NameIn,2:4]=c(SexNew,Name[repro],Founder[repro])
      Names[as.numeric(Name[repro]),"Kids"]=Kids[repro]
      if (Kids[repro]<31) {Names[as.numeric(Name[repro]),8+Kids[repro]]=NameIn}
      if(SexNew==1) {Names[NameIn,"Initial rank"]=Rank[N-Male]}
      Names[NameIn,"Died"]=ndel
      Names[NameIn, "Father"]=Pa
      Names[NameIn,"Phen"]=sum(GeN) 
      
    }}
  
  datout=list(N,Sex, Male, Die, Kids, Rank, Founder, Name,NameIn, Names ,GenM,GenF)
  
  return(list(datout=datout, SexNew=SexNew,nameNewM=nameNewM,GeN=GeN, deadClan=deadClan))
}

# if a male newborn - add to another clan using the following function

update.male=function(datin, nameNewM, GeN) {
  
  for (j in 1:length(datin)){
    assign(datanames[j],datin[[j]])
  }
  
  N=N+1
  Male=Male+1
  Name[N]=nameNewM
  Sex[N]=0
  Rank[N]=N
  GenM=rbind(GenM,GeN)
  rownames(GenM)=NULL
  Die[N]=ifelse(sum(GeN)==0,2*N,N)
  
  datout=list(N,Sex, Male, Die, Kids, Rank, Founder, Name,NameIn, Names ,GenM,GenF)
  
  return(list(datMout=datout))
}


# SIMULATION

#  set initial parameters
N      <- 60
n.rep         <-100
burn.in       <-3000

out=matrix(0,30*n.rep,26)
colnames(out)=c("Mutant.Rank", "Mutant clan","NA","MaleA", "No.Mmutants A", "N.Fmutants A", "NB",
                "MaleB", "No.Mmutants B","N.Fmutants B", "NC","MaleC", "No.Mmutants C", "N.Fmutants C",
                "ND","MaleD", "No.Mmutants D", "N.Fmutants D", "NE", "MaleE","No.Mmutants E",
                "N.Fmutants E", "Total N.iterattion", "simulation","No.Mutants","Prop.Mutants")

for (mutant in 1:30) {
  
  for (k in 1:n.rep) {
    
    N=60
    Male=30
    Sex=c(rep(1,N-Male),rep(0,Male)) #1-female 0-male
    Rank=c(1:N)
    
    
    # set Kids, Die, founder, Name, NameIndex, Names, Gen
    
    Die=N+Rank
    Kids=rep(0,(N-Male))
    Founder=c(1:(N-Male))
    
    GenF=matrix(0,N-Male,2)
    GenM=matrix(0,Male,2)
    
    Name=c(1:N)
    NameIn=N
    
    A=matrix(0,burn.in+N,8)
    colnames(A)=c("Name","Sex","Mother","Founder","Father","Kids","Initial rank","Died")
    B=matrix(0,burn.in+N,30)
    colnames(B)=paste0("Kid ",c(1:30))
    Phen=rep(0,burn.in+N)
    Names=cbind(A,B,Phen)
    
    Names[1:N,1]=Name
    Names[1:N,"Sex"]=Sex
    Names[1:(N-Male),"Founder"]=Founder
    Names[1:N,"Initial rank"]=Rank
    
    # Set clans
    
    Clans= c("A","B","C","D","E")
    
    datlist=list(N,Sex, Male, Die, Kids, Rank, Founder, Name,NameIn,Names ,GenM,GenF)
    
    names(datlist)=c("Nt","Sex", "Male", "Die", "Kids", "Rank", "Founder", "Name","NameIn","Names","GenM","GenF")
    datnames=names(datlist)
    
    datanames=datnames  #for passing to function
    datanames[1]="N"
    
    
    alldat=list()
    for (i in Clans) {
      ldat=list()
      n=1
      for (j in datnames) {
        ldat[[n]]=assign(paste0(j,i),datlist[[j]])
        n=n+1
      }
      names(ldat)=paste0(datnames,i)
      assign(paste0("dat",i),ldat)
    }
    alldat=list(datA=datA,datB=datB,datC=datC,datD=datD,datE=datE)  
    
    
    # Iterations
    
    for (zz in 1:burn.in) {
      
      clan=sample(Clans,1)
      if (zz==1) {clanMut=clan}
      dat=alldat[[which(Clans==clan)]] 
      datSnames=names(dat)
      
      output<-update.network(dat)
      
      dat=output$datout
      names(dat)=datSnames
      assign(paste0("dat",clan),dat)
      alldat[[which(Clans==clan)]]=dat
      names(alldat[[which(Clans==clan)]])=datSnames
      
      
      SexNew=output$SexNew
      nameNewM=output$nameNewM
      GeN=output$GeN
      deadClan=output$deadClan
      
      if (SexNew==0) {
        
        clanM=sample(Clans[-which(Clans==clan)],1)
        dat=alldat[[which(Clans==clanM)]]
        datSnames=names(dat)
        
        outMale=update.male(dat,nameNewM, GeN)
        
        dat=outMale$datMout
        names(dat)=datSnames
        assign(paste0("dat",clanM),dat)
        
        alldat[[which(Clans==clanM)]]=dat
        names(alldat[[which(Clans==clanM)]])=datSnames
      }
      
      
      if (deadClan==TRUE) {
        clandel=which(Clans==clan)
        alldat[[clandel]]=NULL
        Clans=Clans[-which(Clans==clan)]
      }
    }
    
    # output from simulation
    
    alldat=list(datA=datA,datB=datB,datC=datC,datD=datD,datE=datE)
    
    out[n.rep*(mutant-1)+k,1]=mutant
    out[n.rep*(mutant-1)+k,2]=clanMut
    totMut=0
    for (i in 1:5) {
      out[n.rep*(mutant-1)+k,3+4*(i-1)]=alldat[[i]]$Nt
      out[n.rep*(mutant-1)+k,4*i]=alldat[[i]]$Male
      if(alldat[[i]]$Male>1) {
        out[n.rep*(mutant-1)+k,4*i+1]=length(which(rowSums(alldat[[i]][[11]])>0))
        totMut=totMut+length(which(rowSums(alldat[[i]][[11]])>0))
      }
      if(alldat[[i]]$Male==1){
        out[k,4*i+1]=ifelse(sum(alldat[[i]][[11]])>0,1,0)
        totMut=totMut+ceiling(sum(alldat[[i]][[11]])/2)
      }
      if(alldat[[i]]$Nt>alldat[[i]]$Male+1) { 
        out[n.rep*(mutant-1)+k,4*i+2]=length(which(rowSums(alldat[[i]][[12]])>0))
        totMut=totMut+length(which(rowSums(alldat[[i]][[12]])>0))
      }
      if(alldat[[i]]$Nt-alldat[[i]]$Male+1==1){
        out[n.rep*(mutant-1)+k,4*i+2]=ifelse(sum(alldat[[i]][[12]])>0,1,0)
        totMut=totMut+ceiling(sum(alldat[[i]][[12]])/2)      }
      
    }
    out[n.rep*(mutant-1)+k,23]=burn.in
    out[n.rep*(mutant-1)+k,24]=k 
    out[n.rep*(mutant-1)+k,25]=totMut
    out[n.rep*(mutant-1)+k,26]=totMut/(5*N)
    
    
  }#end of k
} #end of mutant

write.table(out,"mutation Sum.csv",col.names=TRUE, row.names=F,sep=",",append=T)

#************************************************************************************

# One mutation - Scenario III


library(Matrix)

update.network <- function(datin) {
  
  for (j in 1:length(datin)){
    assign(datanames[j],datin[[j]])
  }
  
  SexNew=1
  nameNewM=0
  deadClan=FALSE
  
  if(N>Male & Male>0) {
    
    # One individual dies
    
    del=ifelse(zz==1,sample(c((mutant+1):N),1),sample(1:N, 1, prob=Die))
    ndel=Name[del]
    Name=Name[-del]
    Die=Die[-del]
    up=which(Rank>Rank[del])
    Rank[up]=Rank[up]-1
    Rank=Rank[-del]
    
    if (Sex[del]==1) {
      if (N-Male==1) {
        N=Male
        Founder=NULL
        Kids=NULL
        GenF=NULL
        deadClan=TRUE  
      }
      
      if (N-Male>1) {
        Founder=Founder[-del]
        Kids=Kids[-del]
        GenF=GenF[-del,]
      }
    }
    
    if(Sex[del]==0){
      if(Male>1) {GenM=GenM[-(del-N+Male),]}
      else {GenM=NULL}
      Male=Male-1
      if (Male==0) {N=N-1}
    }
    Sex=Sex[-del]
    
    if (deadClan==FALSE & Male>0) {
      
      ChanceF=N-Rank[1:(N-1-Male)]
      ChanceM=1.5*N-Rank[(N-Male):(N-1)]
      
      #  Select female & male for reproduction
      
      repro=ifelse(zz==1,mutant,sample(c(1:(N-1-Male)),1, prob=ChanceF))
      mate=sample(c(1:Male),1,prob=ChanceM)
      Pa=Name[N-Male-1+mate]
      Kids[repro]=Kids[repro]+1
      
      NameIn=NameIn+1
      
      if (zz==1) {
        GenF[repro,]=c(1,0)
        Names[repro,"Phen"]=1 #phenotype
        
      }
      
      # Create variables for newborn
      
      GeN=c(0,0)
      GeN[1]=ifelse(length(GenF)>2,sample(GenF[repro,],1),sample(GenF,1))
      GeN[2]=ifelse(length(GenM)>2,sample(GenM[mate,],1),sample(GenM,1))
      SexNew=sample(c(0,1),1)
      
      if (SexNew==0) {
        N=N-1
        nameNewM=paste0(NameIn,clan)
        Names[NameIn,"Name"]=nameNewM
      }
      
      if(SexNew==1) {
        Name[(N-Male+1):N]=Name[(N-Male):(N-1)]
        Name[N-Male]=NameIn
        Names[NameIn,"Name"]=NameIn
        Sex[N-Male]=1
        if (Male>0) {Sex[N]=0}
        down=which(Rank>Rank[repro])
        Rank[down]=Rank[down]+1
        Rank[(N-Male+1):(N)]=Rank[(N-Male):(N-1)]
        Rank[N-Male]=Rank[repro]+1
        
        Founder[N-Male]=Founder[repro]
        Kids[N-Male]=0
        
        GenF=rbind(GenF,GeN)
        rownames(GenF) <- NULL
        
      }
      
      if(N-Male>1) {
        for (j in 1:(N-Male)) {
          Die[j]=ifelse(sum(GenF[j,])>0,0.5(N+Rank[j]),N+Rank[j])
        }}
      if(N-Male==1) {Die[1]=ifelse(sum(GenF)>0,0.5*(N+1),N+1)}
      
      
      if(Male>1) {
        for (j in (N-Male+1):N) {
          Die[j]=ifelse(sum(GenM[j-N+Male,])>0,0.5*(N+Rank[j]),N+Rank[j])
        }}
      if(Male==1) {Die[N]=ifelse(sum(GenM)>0,N,2*N)}
      
      Names[NameIn,2:4]=c(SexNew,Name[repro],Founder[repro])
      Names[as.numeric(Name[repro]),"Kids"]=Kids[repro]
      if (Kids[repro]<31) {Names[as.numeric(Name[repro]),8+Kids[repro]]=NameIn}
      if(SexNew==1) {Names[NameIn,"Initial rank"]=Rank[N-Male]}
      Names[NameIn,"Died"]=ndel
      Names[NameIn, "Father"]=Pa
      Names[NameIn,"Phen"]=sum(GeN) 
      
    }}
  
  datout=list(N,Sex, Male, Die, Kids, Rank, Founder, Name,NameIn, Names ,GenM,GenF)
  
  return(list(datout=datout, SexNew=SexNew,nameNewM=nameNewM,GeN=GeN, deadClan=deadClan))
  
}


update.male=function(datin, nameNewM, GeN) {
  
  for (j in 1:length(datin)){
    assign(datanames[j],datin[[j]])
  }
  
  N=N+1
  Male=Male+1
  Name[N]=nameNewM
  Sex[N]=0
  Rank[N]=N
  GenM=rbind(GenM,GeN)
  rownames(GenM)=NULL
  Die[N]=ifelse(sum(GeN)==0,2*N,N)
  
  datout=list(N,Sex, Male, Die, Kids, Rank, Founder, Name,NameIn, Names ,GenM,GenF)
  
  return(list(datMout=datout))
}


# SIMULATION

N      <- 60
n.rep         <-100
burn.in       <-3000

out=matrix(0,n.rep*30,26)
colnames(out)=c("Mutant.Rank", "Mutant clan","NA","MaleA", "No.Mmutants A", "N.Fmutants A", "NB",
                "MaleB","No.Mmutants B","N.Fmutants B", "NC","MaleC", "No.Mmutants C", "N.Fmutants C",
                "ND","MaleD","No.Mmutants D", "N.Fmutants D", "NE", "MaleE","No.Mmutants E","N.Fmutants E",
                "Total N.iterattion", "simulation","No.Mutants","Prop.Mutants")

for (mutant in 1:30) {
  
  for (k in 1:n.rep) {
    
    N=60
    Male=30
    Sex=c(rep(1,N-Male),rep(0,Male)) #1-female 0-male
    Rank=c(1:N)
    Die=N+Rank
    Kids=rep(0,(N-Male))
    Founder=c(1:(N-Male))
    
    GenF=matrix(0,N-Male,2)
    GenM=matrix(0,Male,2)
    
    Name=c(1:N)
    NameIn=N
    
    A=matrix(0,burn.in+N,8)
    colnames(A)=c("Name","Sex","Mother","Founder","Father","Kids","Initial rank","Died")
    B=matrix(0,burn.in+N,30)
    colnames(B)=paste0("Kid ",c(1:30))
    Phen=rep(0,burn.in+N)
    Names=cbind(A,B,Phen)
    
    Names[1:N,1]=Name
    Names[1:N,"Sex"]=Sex
    Names[1:(N-Male),"Founder"]=Founder
    Names[1:N,"Initial rank"]=Rank
    
    Clans= c("A","B","C","D","E")
    
    
    datlist=list(N,Sex, Male, Die, Kids, Rank, Founder, Name,NameIn,Names ,GenM,GenF)
    
    names(datlist)=c("Nt","Sex", "Male", "Die", "Kids", "Rank", "Founder", "Name","NameIn","Names","GenM","GenF")
    datnames=names(datlist)
    
    datanames=datnames  #for passing to function
    datanames[1]="N"
    
    
    alldat=list()
    for (i in Clans) {
      ldat=list()
      n=1
      for (j in datnames) {
        ldat[[n]]=assign(paste0(j,i),datlist[[j]])
        n=n+1
      }
      names(ldat)=paste0(datnames,i)
      assign(paste0("dat",i),ldat)
    }
    alldat=list(datA=datA,datB=datB,datC=datC,datD=datD,datE=datE)  
    
    
    # Iterations
    
    for (zz in 1:burn.in) {
      
      clan=sample(Clans,1)
      dat=alldat[[which(Clans==clan)]] 
      datSnames=names(dat)
      
      if(zz==1) {clanMut=clan}    
      
      output<-update.network(dat)
      
      dat=output$datout
      names(dat)=datSnames
      assign(paste0("dat",clan),dat)
      alldat[[which(Clans==clan)]]=dat
      names(alldat[[which(Clans==clan)]])=datSnames
      
      
      SexNew=output$SexNew
      nameNewM=output$nameNewM
      GeN=output$GeN
      deadClan=output$deadClan
      
      if (SexNew==0) {
        
        clanM=sample(Clans[-which(Clans==clan)],1)
        dat=alldat[[which(Clans==clanM)]]
        datSnames=names(dat)
        
        outMale=update.male(dat,nameNewM, GeN)
        
        dat=outMale$datMout
        names(dat)=datSnames
        assign(paste0("dat",clanM),dat)
        
        alldat[[which(Clans==clanM)]]=dat
        names(alldat[[which(Clans==clanM)]])=datSnames
      }
      
      
      if (deadClan==TRUE) {
        clandel=which(Clans==clan)
        alldat[[clandel]]=NULL
        Clans=Clans[-which(Clans==clan)]
      }
    }
    
    alldat=list(datA=datA,datB=datB,datC=datC,datD=datD,datE=datE)
    
    out[n.rep*(mutant-1)+k,1]=mutant
    out[n.rep*(mutant-1)+k,2]=clanMut
    totMut=0
    for (i in 1:5) {
      out[n.rep*(mutant-1)+k,3+4*(i-1)]=alldat[[i]]$Nt
      out[n.rep*(mutant-1)+k,4*i]=alldat[[i]]$Male
      if(alldat[[i]]$Male>1) {
        out[n.rep*(mutant-1)+k,4*i+1]=length(which(rowSums(alldat[[i]][[11]])>0))
        totMut=totMut+length(which(rowSums(alldat[[i]][[11]])>0))
      }
      if(alldat[[i]]$Male==1){
        out[k,4*i+1]=ifelse(sum(alldat[[i]][[11]])>0,1,0)
        totMut=totMut+ceiling(sum(alldat[[i]][[11]])/2)
      }
      if(alldat[[i]]$Nt>alldat[[i]]$Male+1) {
        out[n.rep*(mutant-1)+k,4*i+2]=length(which(rowSums(alldat[[i]][[12]])>0))
        totMut=totMut+length(which(rowSums(alldat[[i]][[12]])>0))
      }
      if(alldat[[i]]$Nt-alldat[[i]]$Male+1==1){
        out[n.rep*(mutant-1)+k,4*i+2]=ifelse(sum(alldat[[i]][[12]])>0,1,0)
        totMut=totMut+ceiling(sum(alldat[[i]][[12]])/2)      }
    }
    out[n.rep*(mutant-1)+k,23]=burn.in
    out[n.rep*(mutant-1)+k,24]=k 
    out[n.rep*(mutant-1)+k,25]=totMut
    out[n.rep*(mutant-1)+k,26]=totMut/(5*N)
    #write.table(netout,"netOutput.csv",col.names=TRUE, row.names=F,sep=",",append=T)
    
    
  }#end of k
}#end of mutant

write.table(out,"mutation MaleR Sum.csv",col.names=TRUE, row.names=F,sep=",",append=T)

#***************************************************************************************

# One mutation - Scenario IV


library(Matrix)

update.network <- function(datin) {
  
  for (j in 1:length(datin)){
    assign(datanames[j],datin[[j]])
  }
  
  SexNew=1
  nameNewM=0
  deadClan=FALSE
  RankM=0
  
  if(N>Male & Male>0) {
    
    # One individual dies
    
    del=ifelse(zz==1,sample(c((mutant+1):N), 1),sample(1:N,1, prob=Die))
    ndel=Name[del]
    Name=Name[-del]
    Die=Die[-del]
    up=which(Rank>Rank[del])
    Rank[up]=Rank[up]-1
    Rank=Rank[-del]
    
    if (Sex[del]==1) {
      if (N-Male==1) {
        N=Male
        Founder=NULL
        Kids=NULL
        GenF=NULL
        deadClan=TRUE  
      }
      
      if (N-Male>1) {
        Founder=Founder[-del]
        Kids=Kids[-del]
        GenF=GenF[-del,]
      }
    }
    
    
    if(Sex[del]==0){
      if(Male>1) {GenM=GenM[-(del-N+Male),]}
      else {GenM=NULL}
      Male=Male-1
      if (Male==0) {N=N-1}
    }
    Sex=Sex[-del]
    
    if (deadClan==FALSE & Male>0) {
      
      ChanceF=N-Rank[1:(N-1-Male)]
      ChanceM=1.5*N-Rank[(N-Male):(N-1)]
      
      #  Select female & male for reproduction
      
      repro=ifelse(zz==1,mutant,sample(c(1:(N-1-Male)),1, prob=ChanceF))
      mate=sample(c(1:Male),1,prob=ChanceM)
      Pa=Name[N-Male-1+mate]
      Kids[repro]=Kids[repro]+1
      
      NameIn=NameIn+1
      
      if (zz==1) {
        GenF[repro,]=c(1,0)
        Names[Name[repro],"Phen"]=1 #phenotype
      }
      
      GeN=c(0,0)
      GeN[1]=ifelse(length(GenF)>2,sample(GenF[repro,],1),sample(GenF,1))
      GeN[2]=ifelse(length(GenM)>2,sample(GenM[mate,],1),sample(GenM,1))
      SexNew=sample(c(0,1),1)
      
      if (SexNew==0) {
        N=N-1
        nameNewM=paste0(NameIn,clan)
        RankM=Rank[repro]
        Names[NameIn,"Name"]=nameNewM
        
      }
      
      if(SexNew==1) {
        Name[(N-Male+1):N]=Name[(N-Male):(N-1)]
        Name[N-Male]=NameIn
        Names[NameIn,"Name"]=NameIn
        Sex[N-Male]=1
        if (Male>0) {Sex[N]=0}
        down=which(Rank>Rank[repro])
        Rank[down]=Rank[down]+1
        Rank[(N-Male+1):(N)]=Rank[(N-Male):(N-1)]
        Rank[N-Male]=Rank[repro]+1
        
        Founder[N-Male]=Founder[repro]
        Kids[N-Male]=0
        
        GenF=rbind(GenF,GeN)
        rownames(GenF) <- NULL
        
      }
      
      if(N-Male>1) {
        for (j in 1:(N-Male)) {
          Die[j]=ifelse(sum(GenF[j,])>0,0.5*(N+Rank[j]),N+Rank[j])
        } }
      if(N-Male==1) {Die[1]=ifelse(sum(GenF)>0,0.5*(N+1),N+1)}
      
      if(Male>1) {
        for (j in (N-Male+1):N) {
          Die[j]=ifelse(sum(GenM[j-N+Male,])>0,0.5*(N+Rank[j]),N+Rank[j])
        }}
      if(Male==1) {Die[N]=ifelse(sum(GenM)>0,N,2*N)}
      
      
      Names[NameIn,2:4]=c(SexNew,Name[repro],Founder[repro])
      Names[as.numeric(Name[repro]),"Kids"]=Kids[repro]
      if (Kids[repro]<31) {Names[as.numeric(Name[repro]),8+Kids[repro]]=NameIn}
      if(SexNew==1) {Names[NameIn,"Initial rank"]=Rank[N-Male]}
      Names[NameIn,"Died"]=ndel
      Names[NameIn, "Father"]=Pa
      Names[NameIn,"Phen"]=sum(GeN) 
      
    }}
  
  
  datout=list(N,Sex, Male, Die, Kids, Rank, Founder, Name,NameIn, Names ,GenM,GenF)
  
  return(list(datout=datout, SexNew=SexNew,nameNewM=nameNewM,GeN=GeN, RankM=RankM, deadClan=deadClan))
  
}


update.male=function(datin, nameNewM, GeN, RankM) {
  
  for (j in 1:length(datin)){
    assign(datanames[j],datin[[j]])
  }
  
  N=N+1
  Male=Male+1
  Name[N]=nameNewM
  Sex[N]=0
  down=which(Rank>RankM+N-Male)
  Rank[down]=Rank[down]+1
  Rank[N]=ifelse(length(down)>0,RankM+N-Male+1,N)
  GenM=rbind(GenM,GeN)
  rownames(GenM)=NULL
  
  if (length(down)>0){
    for (j in ((RankM+1+N-Male):N)) {
      l=which(Rank==j)
      Die[l]=ifelse(sum(GenM[l-N+Male,])>0,0.5*(N+j),N+j)
    }}
  if (length(down)==0) {Die[N]=ifelse(sum(GeN)>0,N,2*N)}
  
  
  datout=list(N,Sex, Male, Die, Kids, Rank, Founder, Name,NameIn, Names ,GenM,GenF)
  
  return(list(datMout=datout))
}


# SIMULATION

N      <- 60
n.rep         <-100
burn.in       <-3000

out=matrix(0,n.rep*30,26)
colnames(out)=c("Mutant.Rank", "Mutant clan","NA","MaleA", "No.Mmutants A", "N.Fmutants A", "NB",
                "MaleB", "No.Mmutants B","N.Fmutants B", "NC","MaleC", "No.Mmutants C", "N.Fmutants C",
                "ND","MaleD", "No.Mmutants D", "N.Fmutants D", "NE", "MaleE","No.Mmutants E", 
                "N.Fmutants E", "Total N.iterattion", "simulation","No.Mutants","Prop.Mutants")

for(mutant in 1:30) {
  
  for (k in 1:n.rep) {
    
    N=60
    Male=30
    Sex=c(rep(1,N-Male),rep(0,Male)) #1=female 0=male
    Rank=c(1:N)
    Die=N+Rank
    Kids=rep(0,(N-Male))
    Founder=c(1:(N-Male))
    
    GenF=matrix(0,N-Male,2)
    GenM=matrix(0,Male,2)
    
    Name=c(1:N)
    NameIn=N
    
    A=matrix(0,burn.in+N,8)
    colnames(A)=c("Name","Sex","Mother","Founder","Father","Kids","Initial rank","Died")
    B=matrix(0,burn.in+N,30)
    colnames(B)=paste0("Kid ",c(1:30))
    Phen=rep(0,burn.in+N)
    Names=cbind(A,B,Phen)
    
    Names[1:N,1]=Name
    Names[1:N,"Sex"]=Sex
    Names[1:(N-Male),"Founder"]=Founder
    Names[1:N,"Initial rank"]=Rank
    
    Clans= c("A","B","C","D","E")
    
    datlist=list(N,Sex, Male, Die, Kids, Rank, Founder, Name,NameIn,Names ,GenM,GenF)
    
    names(datlist)=c("Nt","Sex", "Male", "Die", "Kids", "Rank", "Founder", "Name","NameIn","Names","GenM","GenF")
    datnames=names(datlist)
    
    datanames=datnames  #for passing to function
    datanames[1]="N"
    
    
    alldat=list()
    for (i in Clans) {
      ldat=list()
      n=1
      for (j in datnames) {
        ldat[[n]]=assign(paste0(j,i),datlist[[j]])
        n=n+1
      }
      names(ldat)=paste0(datnames,i)
      assign(paste0("dat",i),ldat)
    }
    alldat=list(datA=datA,datB=datB,datC=datC,datD=datD,datE=datE)  
    
    
    # Iterations
    
    for (zz in 1:burn.in) {
      
      clan=sample(Clans,1)
      dat=alldat[[which(Clans==clan)]] 
      datSnames=names(dat)
      
      if(zz==1) {clanMut=clan} 
      
      output<-update.network(dat)
      
      dat=output$datout
      names(dat)=datSnames
      assign(paste0("dat",clan),dat)
      alldat[[which(Clans==clan)]]=dat
      names(alldat[[which(Clans==clan)]])=datSnames
      
      
      SexNew=output$SexNew
      nameNewM=output$nameNewM
      GeN=output$GeN
      RankM=output$RankM
      deadClan=output$deadClan
      
      if (SexNew==0) {
        clanM=sample(Clans[-which(Clans==clan)],1)
        dat=alldat[[which(Clans==clanM)]]
        datSnames=names(dat)
        
        outMale=update.male(dat,nameNewM, GeN,RankM)
        
        dat=outMale$datMout
        names(dat)=datSnames
        assign(paste0("dat",clanM),dat)
        
        alldat[[which(Clans==clanM)]]=dat
        names(alldat[[which(Clans==clanM)]])=datSnames
      }
      
      
      if (deadClan==TRUE) {
        clandel=which(Clans==clan)
        alldat[[clandel]]=NULL
        Clans=Clans[-which(Clans==clan)]
      }
    }
    
    alldat=list(datA=datA,datB=datB,datC=datC,datD=datD,datE=datE)
    
    out[n.rep*(mutant-1)+k,1]=mutant
    out[n.rep*(mutant-1)+k,2]=clanMut
    totMut=0
    for (i in 1:5) {
      out[n.rep*(mutant-1)+k,3+4*(i-1)]=alldat[[i]]$Nt
      out[n.rep*(mutant-1)+k,4*i]=alldat[[i]]$Male
      if(alldat[[i]]$Male>1) {
        out[n.rep*(mutant-1)+k,4*i+1]=length(which(rowSums(alldat[[i]][[11]])>0))
        totMut=totMut+length(which(rowSums(alldat[[i]][[11]])>0))
      }
      if(alldat[[i]]$Male==1){
        out[n.rep*(mutant-1)+k,4*i+1]=ifelse(sum(alldat[[i]][[11]])>0,1,0)
        totMut=totMut+ceiling(sum(alldat[[i]][[11]])/2)
      }
      if(alldat[[i]]$Nt>alldat[[i]]$Male+1) { 
        out[n.rep*(mutant-1)+k,4*i+2]=length(which(rowSums(alldat[[i]][[12]])>0))
        totMut=totMut+length(which(rowSums(alldat[[i]][[12]])>0))
      }
      if(alldat[[i]]$Nt-alldat[[i]]$Male+1==1){
        out[n.rep*(mutant-1)+k,4*i+2]=ifelse(sum(alldat[[i]][[12]])>0,1,0)
        totMut=totMut+ceiling(sum(alldat[[i]][[12]])/2)      }
    }
    out[n.rep*(mutant-1)+k,23]=burn.in
    out[n.rep*(mutant-1)+k,24]=k 
    out[n.rep*(mutant-1)+k,25]=totMut
    out[n.rep*(mutant-1)+k,26]=totMut/(5*N)
    
    
  }#end of k
}#end of mutant

write.table(out,"mutation RankMa Sum.csv",col.names=TRUE, row.names=F,sep=",",append=T)

#*********************************************************************************************

# Two mutations - Scenario I


library(Matrix)

update.network <- function(datin) {
  
  for (j in 1:length(datin)){
    assign(datanames[j],datin[[j]])
  }
  
  SexNew=1
  nameNewM=0
  deadClan=FALSE
  GeN=rep(0,4)
  
  if(N>Male & Male>0) {
    
    # One individual dies
    
    del=sample(1:N, 1, prob=Die)
    ndel=Name[del]
    Name=Name[-del]
    Die=Die[-del]
    
    if (Sex[del]==1) {
      if (N-Male==1) {
        N=Male
        Founder=NULL
        Kids=NULL
        Rank=NULL
        GenF=NULL
        deadClan=TRUE  
      }
      
      if (N-Male>1) {
        Founder=Founder[-del]
        Kids=Kids[-del]
        up=which(Rank>Rank[del])
        Rank[up]=Rank[up]-1
        Rank=Rank[-del]
        GenF=GenF[-del,]
      }
    }
    
    
    if(Sex[del]==0){
      if(Male>1) {GenM=GenM[-(del-N+Male),]}
      else {GenM=NULL}
      Male=Male-1
      if (Male==0) {N=N-1}
    }
    Sex=Sex[-del]
    
    if (deadClan==FALSE & Male>0) {
      
      Chance=N-Rank
      
      # Select female & male for reproduction
      
      repro=ifelse(zz==1,sample(which(Rank %in% c(25:30)),1),
                   ifelse(zz==2,sample(which(Rank %in% c(1:6)),1),
                          sample(c(1:(N-1-Male)),1, prob=Chance)))
      mate=sample( c(1:Male),1)
      Pa=Name[N-Male-1+mate]
      Kids[repro]=Kids[repro]+1
      
      NameIn=NameIn+1
      
      if (zz==1) {    
        GenF[repro,1:2]=c(1,0) #Mutantion A; originating from low rank 
        Die[repro]=1
        Names[Name[repro],"Phen"]=1 #phenotype
        mutant=Name[repro]
      } # 
      
      if (zz==2) {
        GenF[repro,3]=1 # Mutation B; originating from high rank
        Die[repro]=1
        Names[as.numeric(Name[repro]),"Phen"]=2
        mutant=Name[repro]
      }
      
      GeN[1]=ifelse(length(GenF)>4,sample(GenF[repro,1:2],1),sample(GenF[1:2],1))
      GeN[2]=ifelse(length(GenM)>4,sample(GenM[mate,1:2],1),sample(GenM[1:2],1))
      GeN[3]=ifelse(length(GenF)>4,sample(GenF[repro,3:4],1),sample(GenF[3:4],1))
      GeN[4]=ifelse(length(GenM)>4,sample(GenM[mate,3:4],1),sample(GenM[3:4],1))
      
      SexNew=sample(c(0,1),1)
      
      if (SexNew==0) {
        N=N-1
        nameNewM=paste0(NameIn,clan)
        Names[NameIn,"Name"]=nameNewM
        
      }
      
      if(SexNew==1) {
        
        Name[(N-Male+1):N]=Name[(N-Male):(N-1)]
        Name[N-Male]=NameIn
        Names[NameIn,"Name"]=NameIn
        Sex[N-Male]=1
        if (Male>0) {Sex[N]=0}
        down=which(Rank>Rank[repro])
        Rank[down]=Rank[down]+1
        Rank[N-Male]=Rank[repro]+1
        
        Founder[N-Male]=Founder[repro]
        Kids[N-Male]=0
        
        GenF=rbind(GenF,GeN)
        rownames(GenF) <- NULL
        
        Die[(N-Male+1):N]=Die[(N-Male):(N-1)]
        Die[N-Male]=ifelse(sum(GeN)>0,1,2)
        
      }
      
      Names[NameIn,2:4]=c(SexNew,Name[repro],Founder[repro])
      Names[as.numeric(Name[repro]),"Kids"]=Kids[repro]
      if(Kids[repro]<31) {Names[as.numeric(Name[repro]),8+Kids[repro]]=NameIn}
      if(SexNew==1) {Names[NameIn,"Initial rank"]=Rank[N-Male]}
      Names[NameIn,"Died"]=ndel
      Names[NameIn, "Father"]=Pa
      
     Names[NameIn,"Phen"]=ifelse((sum(GeN[1:2]>0 & sum(GeN[3:4])==0))>1,1,  
     ifelse((sum(GeN[1:2])>0 & sum(GeN[3:4]))>0,3, ifelse((sum(GeN[1:2])==0 & 
     sum(GeN[3:4])>0),2,0)))   # phenotype: 1-mutA(lowR);2-mutB(highR);3-both;0-none
      
      
    }}
  
  
  datout=list(N,Sex, Male, Die, Kids, Rank, Founder, Name,NameIn, Names ,GenM,GenF)
  
  if (zz<3) {return(list(datout=datout,SexNew=SexNew,nameNewM=nameNewM,GeN=GeN, deadClan=FALSE,mutant=mutant))} 
  else {
    return(list(datout=datout, SexNew=SexNew,nameNewM=nameNewM,GeN=GeN, deadClan=deadClan))
  }
}


update.male=function(datin, nameNewM, GeN) {
  
  for (j in 1:length(datin)){
    assign(datanames[j],datin[[j]])
  }
  
  N=N+1
  Male=Male+1
  Name[N]=nameNewM
  Sex[N]=0
  GenM=rbind(GenM,GeN)
  rownames(GenM)=NULL
  Die[N]=ifelse(sum(GeN)>0,1,2)
  
  datout=list(N,Sex, Male, Die, Kids, Rank, Founder, Name,NameIn, Names ,GenM,GenF)
  
  return(list(datMout=datout))
}


# Simulation

N      <- 60
n.rep         <-200
burn.in       <-3000

out=matrix(0,n.rep,31)
colnames(out)=c("MutantA.Rank", "MutantA clan","MutantB.Rank", "MutantB clan","NA","MaleA",
                "MutantsA_A", "MutantsB_A", "NB","MaleB", "MutantsA_B","MutantsB_B", "NC",
                "MaleC", "MutantsA_C", "MutantsB_C","ND","MaleD", "MutantsA_D", "MutantBs_D",
                "NE", "MaleE","MutantsA_E", "MutantsB_E","Total N","N.iterattions","simulation",
                "No.MutantsA","Prop.MutantsA","No.MutantsB","Prop.MutantsB")


for (k in 1:n.rep) {
  
  N=60
  Male=30
  Sex=c(rep(1,N-Male),rep(0,Male)) #1-female 0-male
  Rank=c(1:(N-Male))
  Die=rep(2,N)
  Kids=rep(0,(N-Male))
  Founder=c(1:(N-Male))
  
  GenF=matrix(0,N-Male,4)
  GenM=matrix(0,Male,4)
  
  Name=c(1:N)
  NameIn=N
  
  A=matrix(0,burn.in+N,8)
  colnames(A)=c("Name","Sex","Mother","Founder","Father","Kids","Initial rank","Died")
  B=matrix(0,burn.in+N,30)
  colnames(B)=paste0("Kid ",c(1:30))
  Phen=rep(0,burn.in+N)
  Names=cbind(A,B,Phen)
  
  Names[1:N,1]=Name
  Names[1:N,"Sex"]=Sex
  Names[1:(N-Male),"Founder"]=Founder
  Names[1:(N-Male),"Initial rank"]=Rank
  
  Clans= c("A","B","C","D","E")
  
  datlist=list(N,Sex, Male, Die, Kids, Rank, Founder, Name,NameIn,Names ,GenM,GenF)
  
  names(datlist)=c("Nt","Sex","Male","Die","Kids","Rank","Founder","Name","NameIn","Names",
                   "GenM","GenF")
  datnames=names(datlist)
  
  datanames=datnames  #for passing to function
  datanames[1]="N"
  
  
  alldat=list()
  for (i in Clans) {
    ldat=list()
    n=1
    for (j in datnames) {
      ldat[[n]]=assign(paste0(j,i),datlist[[j]])
      n=n+1
    }
    names(ldat)=paste0(datnames,i)
    assign(paste0("dat",i),ldat)
  }
  alldat=list(datA=datA,datB=datB,datC=datC,datD=datD,datE=datE)  
  
  
  # Iterations
  
  for (zz in 1:burn.in) {
    
    clan=sample(Clans,1)
    dat=alldat[[which(Clans==clan)]] 
    datSnames=names(dat)
    
    output<-update.network(dat)
    if(zz==1) {mutantA=output$mutant
    clanMutA=clan}    
    if(zz==2) {mutantB=output$mutant
    clanMutB=clan} 
    
    dat=output$datout
    names(dat)=datSnames
    assign(paste0("dat",clan),dat)
    alldat[[which(Clans==clan)]]=dat
    names(alldat[[which(Clans==clan)]])=datSnames
    
    
    SexNew=output$SexNew
    nameNewM=output$nameNewM
    GeN=output$GeN
    deadClan=output$deadClan
    
    if (SexNew==0) {
      clanM=sample(Clans[-which(Clans==clan)],1)
      dat=alldat[[which(Clans==clanM)]]
      datSnames=names(dat)
      
      outMale=update.male(dat,nameNewM, GeN)
      
      dat=outMale$datMout
      names(dat)=datSnames
      assign(paste0("dat",clanM),dat)
      
      alldat[[which(Clans==clanM)]]=dat
      names(alldat[[which(Clans==clanM)]])=datSnames
    }
    
    
    if (deadClan==TRUE) {
      clandel=which(Clans==clan)
      alldat[[clandel]]=NULL
      Clans=Clans[-which(Clans==clan)]
    }
  }
  
  alldat=list(datA=datA,datB=datB,datC=datC,datD=datD,datE=datE)
  
  out[k,1]=mutantA
  out[k,2]=clanMutA
  out[k,3]=mutantB
  out[k,4]=clanMutB
  totMutA=0
  totMutB=0
  TotalN=0
  for (i in 1:5) {
    out[k,4*i+1]=alldat[[i]]$Nt
    TotalN=TotalN+alldat[[i]]$Nt
    out[k,4*i+2]=alldat[[i]]$Male
    mutA=0
    mutB=0
    if(alldat[[i]]$Male>1) {    
      mutA=length((which((alldat[[i]]$GenM[,1]+alldat[[i]]$GenM[,2])>0)))
      mutB=length((which((alldat[[i]]$GenM[,3]+alldat[[i]]$GenM[,4])>0)))
    }
    
    if(alldat[[i]]$Male==1){
      if (sum(alldat[[i]]$GenM[1:2])>0) {mutA=1}
      if (sum(alldat[[i]]$GenM[3:4])>0) {mutB=1} 
    }
    
    
    if(alldat[[i]]$Nt>alldat[[i]]$Male+1) {
      mutA=mutA+length((which((alldat[[i]]$GenF[,1]+alldat[[i]]$GenF[,2])>0)))
      mutB=mutB+length((which((alldat[[i]]$GenF[,3]+alldat[[i]]$GenF[,4])>0)))
    }
    
    if(alldat[[i]]$Nt-alldat[[i]]$Male+1==1){
      mutA=mutA+ifelse(sum(alldat[[i]]$GenF[1:2])>0,1,0)
      mutB=mutB+ifelse(sum(alldat[[i]]$GenF[3:4])>0,1,0)
    }
    out[k,4*i+3]=mutA
    out[k,4*i+4]=mutB
    
    totMutA=totMutA+mutA  
    totMutB=totMutB+mutB
  }
  
  out[k,25]=TotalN 
  out[k,26]=burn.in
  out[k,27]=k 
  out[k,28]=totMutA
  out[k,29]=totMutA/TotalN
  out[k,30]=totMutB
  out[k,31]=totMutB/TotalN
  
}#end of k


write.table(out,"Two mutations Sum.csv",col.names=TRUE, row.names=F,sep=",",append=T)

#*******************************************************************************************

#  Two mutations - Scenario II


library(Matrix)

update.network <- function(datin) {
  
  for (j in 1:length(datin)){
    assign(datanames[j],datin[[j]])
  }
  
  SexNew=1
  nameNewM=0
  deadClan=FALSE
  GeN=rep(0,4)
  
  if(N>Male & Male>0) {
    
    # One individual dies
    
    del=sample(1:N, 1, prob=Die)
    ndel=Name[del]
    Name=Name[-del]
    Die=Die[-del]
    up=which(Rank>Rank[del])
    Rank[up]=Rank[up]-1
    Rank=Rank[-del]
    
    if (Sex[del]==1) {
      if (N-Male==1) {
        N=Male
        Founder=NULL
        Kids=NULL
        GenF=NULL
        deadClan=TRUE  
      }
      
      if (N-Male>1) {
        Founder=Founder[-del]
        Kids=Kids[-del]
        GenF=GenF[-del,]
      }
    }
    
    if(Sex[del]==0){
      if(Male>1) {GenM=GenM[-(del-N+Male),]}
      else {GenM=NULL}
      Male=Male-1
      if (Male==0) {N=N-1}
    }
    Sex=Sex[-del]
    
    if (deadClan==FALSE & Male>0) {
      
      Chance=N-Rank[1:(N-1-Male)]
      
      # Select female & male for reproduction
      
      repro=ifelse(zz==1,sample(which(Rank %in% 25:30 & Sex==1),1),
                   ifelse(zz==2,sample(which(Rank %in% 1:6),1),sample(c(1:(N-1-Male)),1,prob=Chance)))
      mate=sample(c(1:Male),1)
      Pa=Name[N-Male-1+mate]
      Kids[repro]=Kids[repro]+1
      
      NameIn=NameIn+1
      
      if (zz==1) {
        GenF[repro,1:2]=c(1,0)  #Mutation A originates in low rank
        Names[Name[repro],"Phen"]=1 #phenotype 
        mutant=Name[repro]
      }
      
      if (zz==2) {
        GenF[repro,3]=1  # Mutation B originates in high rank
        Names[as.numeric(Name[repro]),"Phen"]=2 
        mutant=Name[repro]
      }
      
      GeN[1]=ifelse(length(GenF)>4,sample(GenF[repro,1:2],1),sample(GenF[1:2],1))
      GeN[2]=ifelse(length(GenM)>4,sample(GenM[mate,1:2],1),sample(GenM[1:2],1))
      GeN[3]=ifelse(length(GenF)>4,sample(GenF[repro,3:4],1),sample(GenF[3:4],1))
      GeN[4]=ifelse(length(GenM)>4,sample(GenM[mate,3:4],1),sample(GenM[3:4],1))
      
      SexNew=sample(c(0,1),1)
      
      if (SexNew==0) {
        N=N-1
        nameNewM=paste0(NameIn,clan)
        Names[NameIn,"Name"]=nameNewM
      }
      
      if(SexNew==1) {
        
        Name[(N-Male+1):N]=Name[(N-Male):(N-1)]
        Name[N-Male]=NameIn
        Names[NameIn,"Name"]=NameIn
        Sex[N-Male]=1
        if (Male>0) {Sex[N]=0}
        down=which(Rank>Rank[repro])
        Rank[down]=Rank[down]+1
        Rank[(N-Male+1):(N)]=Rank[(N-Male):(N-1)]
        Rank[N-Male]=Rank[repro]+1
        
        Founder[N-Male]=Founder[repro]
        Kids[N-Male]=0
        
        GenF=rbind(GenF,GeN)
        rownames(GenF) <- NULL
      }
      
      
      if(N-Male>1) {
        for (j in 1:(N-Male)) {
          Die[j]=ifelse(sum(GenF[j,])>0,0.5*(N+Rank[j]),N+Rank[j])
        }}
      if(N-Male==1) {Die[1]=ifelse(sum(GenF)>0,0.5*(N+1),N+1)}
      
      if(Male>1) {
        for (j in (N-Male+1):N) {
          Die[j]=ifelse(sum(GenM[j-N+Male,])>0,0.5(N+Rank[j]),N+Rank[j])
        }}
      if(Male==1) {Die[N]=ifelse(sum(GenM)>0,N,2*N)}
      
      
      
      Names[NameIn,2:4]=c(SexNew,Name[repro],Founder[repro])
      Names[as.numeric(Name[repro]),"Kids"]=Kids[repro]
      if(Kids[repro]<31) {Names[as.numeric(Name[repro]),8+Kids[repro]]=NameIn}
      if(SexNew==1) {Names[NameIn,"Initial rank"]=Rank[N-Male]}
      Names[NameIn,"Died"]=ndel
      Names[NameIn, "Father"]=Pa
      
      Names[NameIn,"Phen"]=ifelse((sum(GeN[1:2]>0 & sum(GeN[3:4])==0))>1,1,  
      ifelse((sum(GeN[1:2])>0 & sum(GeN[3:4]))>0,3, ifelse((sum(GeN[1:2])==0 & 
      sum(GeN[3:4])>0),2,0)))
      
    }}
  
  datout=list(N,Sex, Male, Die, Kids, Rank, Founder, Name,NameIn, Names ,GenM,GenF)
  
  if (zz<3) {return(list(datout=datout,SexNew=SexNew,nameNewM=nameNewM,GeN=GeN, mutant=mutant,deadClan=FALSE))} 
  else {
    return(list(datout=datout, SexNew=SexNew,nameNewM=nameNewM,GeN=GeN, deadClan=deadClan))
  }
}


update.male=function(datin, nameNewM, GeN) {
  
  for (j in 1:length(datin)){
    assign(datanames[j],datin[[j]])
  }
  
  N=N+1
  Male=Male+1
  Name[N]=nameNewM
  Sex[N]=0
  Rank[N]=N
  GenM=rbind(GenM,GeN)
  rownames(GenM)=NULL
  Die[N]=ifelse(sum(GeN)>0,N,2*N)
  
  datout=list(N,Sex, Male, Die, Kids, Rank, Founder, Name,NameIn, Names ,GenM,GenF)
  
  return(list(datMout=datout))
}


# Simulation

N      <- 60
n.rep         <-200
burn.in       <-3000

out=matrix(0,n.rep,31)
colnames(out)=c("MutantA.Rank","MutantA clan","MutantB.Rank","MutantB clan","NA","MaleA","MutantsA_A",
                "MutantsB_A","NB","MaleB","MutantsA_B","MutantsB_B","NC","MaleC","MutantsA_C",
                "MutantsB_C","ND","MaleD", "MutantsA_D", "MutantBs_D", "NE", "MaleE","MutantsA_E",
                "MutantsB_E","Total N","N.iterattions","simulation","No.MutantsA","Prop.MutantsA",
                "No.MutantsB","Prop.MutantsB")


for (k in 1:n.rep) {
  
  N=60
  Male=30
  Sex=c(rep(1,N-Male),rep(0,Male)) #1=female 0=male
  Rank=c(1:N)
  Die=Rank+N
  Kids=rep(0,(N-Male))
  Founder=c(1:(N-Male))
  
  GenF=matrix(0,N-Male,4)
  GenM=matrix(0,Male,4)
  
  Name=c(1:N)
  NameIn=N
  
  A=matrix(0,burn.in+N,8)
  colnames(A)=c("Name","Sex","Mother","Founder","Father","Kids","Initial rank","Died")
  B=matrix(0,burn.in+N,30)
  colnames(B)=paste0("Kid ",c(1:30))
  Phen=rep(0,burn.in+N)
  Names=cbind(A,B,Phen)
  
  Names[1:N,1]=Name
  Names[1:N,"Sex"]=Sex
  Names[1:(N-Male),"Founder"]=Founder
  Names[1:N,"Initial rank"]=Rank
  
  Clans= c("A","B","C","D","E")
  
  
  datlist=list(N,Sex, Male, Die, Kids, Rank, Founder, Name,NameIn,Names ,GenM,GenF)
  
  names(datlist)=c("Nt","Sex","Male", "Die", "Kids", "Rank", "Founder", "Name","NameIn","Names","GenM","GenF")
  datnames=names(datlist)
  
  datanames=datnames  #for passing to function
  datanames[1]="N"
  
  
  alldat=list()
  for (i in Clans) {
    ldat=list()
    n=1
    for (j in datnames) {
      ldat[[n]]=assign(paste0(j,i),datlist[[j]])
      n=n+1
    }
    names(ldat)=paste0(datnames,i)
    assign(paste0("dat",i),ldat)
  }
  alldat=list(datA=datA,datB=datB,datC=datC,datD=datD,datE=datE)  
  
  
  # Iterations
  
  for (zz in 1:burn.in) {
    
    clan=sample(Clans,1)
    dat=alldat[[which(Clans==clan)]] 
    datSnames=names(dat)
    
    output<-update.network(dat)
    if(zz==1) {mutantA=output$mutant
    clanMutA=clan}    
    if(zz==2) {mutantB=output$mutant
    clanMutB=clan} 
    
    dat=output$datout
    names(dat)=datSnames
    assign(paste0("dat",clan),dat)
    alldat[[which(Clans==clan)]]=dat
    names(alldat[[which(Clans==clan)]])=datSnames
    
    SexNew=output$SexNew
    nameNewM=output$nameNewM
    GeN=output$GeN
    deadClan=output$deadClan
    
    if (SexNew==0) {
      clanM=sample(Clans[-which(Clans==clan)],1)
      dat=alldat[[which(Clans==clanM)]]
      datSnames=names(dat)
      
      outMale=update.male(dat,nameNewM, GeN)
      
      dat=outMale$datMout
      names(dat)=datSnames
      assign(paste0("dat",clanM),dat)
      
      alldat[[which(Clans==clanM)]]=dat
      names(alldat[[which(Clans==clanM)]])=datSnames
    }
    
    
    if (deadClan==TRUE) {
      clandel=which(Clans==clan)
      alldat[[clandel]]=NULL
      Clans=Clans[-which(Clans==clan)]
    }
  }
  
  alldat=list(datA=datA,datB=datB,datC=datC,datD=datD,datE=datE)
  
  out[k,1]=mutantA
  out[k,2]=clanMutA
  out[k,3]=mutantB
  out[k,4]=clanMutB
  totMutA=0
  totMutB=0
  TotalN=0
  for (i in 1:5) {
    out[k,4*i+1]=alldat[[i]]$Nt
    TotalN=TotalN+alldat[[i]]$Nt
    out[k,4*i+2]=alldat[[i]]$Male
    mutA=0
    mutB=0
    if(alldat[[i]]$Male>1) {    
      mutA=length((which((alldat[[i]]$GenM[,1]+alldat[[i]]$GenM[,2])>0)))
      mutB=length((which((alldat[[i]]$GenM[,3]+alldat[[i]]$GenM[,4])>0)))
    }
    
    if(alldat[[i]]$Male==1){
      if (sum(alldat[[i]]$GenM[1:2])>0) {mutA=1}
      if (sum(alldat[[i]]$GenM[3:4])>0) {mutB=1} 
    }
    
    
    if(alldat[[i]]$Nt>alldat[[i]]$Male+1) {
      mutA=mutA+length((which((alldat[[i]]$GenF[,1]+alldat[[i]]$GenF[,2])>0)))
      mutB=mutB+length((which((alldat[[i]]$GenF[,3]+alldat[[i]]$GenF[,4])>0)))
    }
    
    if(alldat[[i]]$Nt-alldat[[i]]$Male+1==1){
      mutA=mutA+ifelse(sum(alldat[[i]]$GenF[1:2])>0,1,0)
      mutB=mutB+ifelse(sum(alldat[[i]]$GenF[3:4])>0,1,0)
    }
    out[k,4*i+3]=mutA
    out[k,4*i+4]=mutB
    
    totMutA=totMutA+mutA  
    totMutB=totMutB+mutB
  }
  
  out[k,25]=TotalN 
  out[k,26]=burn.in
  out[k,27]=k 
  out[k,28]=totMutA
  out[k,29]=totMutA/TotalN
  out[k,30]=totMutB
  out[k,31]=totMutB/TotalN
  
}#end of k


write.table(out,"Two mutations MaleR Sum.csv",col.names=TRUE,row.names=F,sep=",",append=T)

#*******************************************************************************************


# Two mutations = Scenario III


library(Matrix)

update.network <- function(datin) {
  
  for (j in 1:length(datin)){
    assign(datanames[j],datin[[j]])
  }
  
  SexNew=1
  nameNewM=0
  deadClan=FALSE
  GeN=rep(0,4)
  
  if(N>Male & Male>0) {
    
    # One individual dies
    
    del=sample(1:N, 1, prob=Die)
    ndel=Name[del]
    Name=Name[-del]
    Die=Die[-del]
    up=which(Rank>Rank[del])
    Rank[up]=Rank[up]-1
    Rank=Rank[-del]
    
    if (Sex[del]==1) {
      if (N-Male==1) {
        N=Male
        Founder=NULL
        Kids=NULL
        GenF=NULL
        deadClan=TRUE  
      }
      
      if (N-Male>1) {
        Founder=Founder[-del]
        Kids=Kids[-del]
        GenF=GenF[-del,]
      }
    }
    
    
    if(Sex[del]==0){
      if(Male>1) {GenM=GenM[-(del-N+Male),]}
      else {GenM=NULL}
      Male=Male-1
      if (Male==0) {N=N-1}
    }
    Sex=Sex[-del]
    
    if (deadClan==FALSE & Male>0) {
      
      ChanceF=N-Rank[1:(N-1-Male)]
      ChanceM=1.5*N-Rank[(N-Male):(N-1)]
      
      # Select female & male for reproduction
      
      repro=ifelse(zz==1,sample(which(Rank %in% 25:30 & Sex==1),1),
                   ifelse(zz==2,sample(which(Rank %in% 1:6),1),sample(c(1:(N-1-Male)),1,prob=ChanceF)))
      mate=sample(c(1:Male),1,prob=ChanceM)
      Pa=Name[N-Male-1+mate]
      Kids[repro]=Kids[repro]+1
      
      NameIn=NameIn+1
      
      if (zz==1) {
        GenF[repro,1:2]=c(1,0)
        Names[Name[repro],"Phen"]=1 #phenotype
        mutant=Name[repro]
      }
      
      if (zz==2) {
        GenF[repro,3]=1
        Names[as.numeric(Name[repro]),"Phen"]=2
        mutant=Name[repro]
      }
      
      GeN[1]=ifelse(length(GenF)>4,sample(GenF[repro,1:2],1),sample(GenF[1:2],1))
      GeN[2]=ifelse(length(GenM)>4,sample(GenM[mate,1:2],1),sample(GenM[1:2],1))
      GeN[3]=ifelse(length(GenF)>4,sample(GenF[repro,3:4],1),sample(GenF[3:4],1))
      GeN[4]=ifelse(length(GenM)>4,sample(GenM[mate,3:4],1),sample(GenM[3:4],1))
      SexNew=sample(c(0,1),1)
      
      if (SexNew==0) {
        N=N-1
        nameNewM=paste0(NameIn,clan)
        Names[NameIn,"Name"]=nameNewM
      }
      
      if(SexNew==1) {
        
        Name[(N-Male+1):N]=Name[(N-Male):(N-1)]
        Name[N-Male]=NameIn
        Names[NameIn,"Name"]=NameIn
        Sex[N-Male]=1
        if (Male>0) {Sex[N]=0}
        down=which(Rank>Rank[repro])
        Rank[down]=Rank[down]+1
        Rank[(N-Male+1):(N)]=Rank[(N-Male):(N-1)]
        Rank[N-Male]=Rank[repro]+1
        
        Founder[N-Male]=Founder[repro]
        Kids[N-Male]=0
        
        GenF=rbind(GenF,GeN)
        rownames(GenF) <- NULL
      }
      
      if(N-Male>1) {
        for (j in 1:(N-Male)) {
          Die[j]=ifelse(sum(GenF[j,])>0,0.5*(N+Rank[j]),N+Rank[j])}
      }
      if(N-Male==1) {Die[1]=ifelse(sum(GenF)>0,0.5*(N+1),N+1)}
      
      if(Male>1) {
        for (j in (N-Male+1):N) {
          Die[j]=ifelse(sum(GenM[j-N+Male,])>0,0.5*(N+Rank[j]),N+Rank[j])}
      }
      if(Male==1) {Die[N]=ifelse(sum(GenM)>0,N,2*N)}
      
      Names[NameIn,2:4]=c(SexNew,Name[repro],Founder[repro])
      Names[as.numeric(Name[repro]),"Kids"]=Kids[repro]
      if(Kids[repro]<31) {Names[as.numeric(Name[repro]),8+Kids[repro]]=NameIn}
      if(SexNew==1) {Names[NameIn,"Initial rank"]=Rank[N-Male]}
      Names[NameIn,"Died"]=ndel
      Names[NameIn, "Father"]=Pa
      
      Names[NameIn,"Phen"]=ifelse((sum(GeN[1:2]>0 & sum(GeN[3:4])==0))>1,1,  
      ifelse((sum(GeN[1:2])>0 & sum(GeN[3:4]))>0,3, ifelse((sum(GeN[1:2])==0 & 
      sum(GeN[3:4])>0),2,0)))
                  
      
    }}
  
  
  datout=list(N,Sex, Male, Die, Kids, Rank, Founder, Name,NameIn, Names ,GenM,GenF)
  
  if (zz<3) {return(list(datout=datout,SexNew=SexNew,nameNewM=nameNewM,GeN=GeN,mutant=mutant,deadClan=FALSE))} 
  else {
    return(list(datout=datout, SexNew=SexNew,nameNewM=nameNewM,GeN=GeN,deadClan=deadClan))
  }
 }


update.male=function(datin, nameNewM, GeN) {
  
  for (j in 1:length(datin)){
    assign(datanames[j],datin[[j]])
  }
  
  N=N+1
  Male=Male+1
  Name[N]=nameNewM
  Sex[N]=0
  Rank[N]=N
  GenM=rbind(GenM,GeN)
  rownames(GenM)=NULL
  Die[N]=ifelse(sum(GeN)>0,N,2*N)
  
  datout=list(N,Sex, Male, Die, Kids, Rank, Founder, Name,NameIn, Names ,GenM,GenF)
  
  return(list(datMout=datout))
}


# Simulation

N      <- 60
n.rep         <-200
burn.in       <-3000

out=matrix(0,n.rep,31)
colnames(out)=c("MutantA.Rank","MutantA clan","MutantB.Rank","MutantB clan","NA","MaleA",
                "MutantsA_A","MutantsB_A","NB","MaleB","MutantsA_B","MutantsB_B","NC",
                "MaleC","MutantsA_C","MutantsB_C","ND","MaleD","MutantsA_D","MutantBs_D",
                "NE","MaleE","MutantsA_E","MutantsB_E","Total N","N.iterattions",
                "simulation","No.MutantsA","Prop.MutantsA","No.MutantsB","Prop.MutantsB")


for (k in 1:n.rep) {
  
  N=60
  Male=30
  Sex=c(rep(1,N-Male),rep(0,Male)) #1-female 0-male
  Rank=c(1:N)
  
  Die=Rank+N
  Kids=rep(0,(N-Male))
  Founder=c(1:(N-Male))
  
  GenF=matrix(0,N-Male,4)
  GenM=matrix(0,Male,4)
  
  Name=c(1:N)
  NameIn=N
  
  A=matrix(0,burn.in+N,8)
  colnames(A)=c("Name","Sex","Mother","Founder","Father","Kids","Initial rank","Died")
  B=matrix(0,burn.in+N,30)
  colnames(B)=paste0("Kid ",c(1:30))
  Phen=rep(0,burn.in+N)
  Names=cbind(A,B,Phen)
  
  Names[1:N,1]=Name
  Names[1:N,"Sex"]=Sex
  Names[1:(N-Male),"Founder"]=Founder
  Names[1:N,"Initial rank"]=Rank
  
  Clans= c("A","B","C","D","E")
  
  
  datlist=list(N,Sex, Male, Die, Kids, Rank, Founder, Name,NameIn,Names ,GenM,GenF)
  
  names(datlist)=c("Nt","Sex","Male","Die","Kids","Rank","Founder","Name","NameIn","Names","GenM","GenF")
  datnames=names(datlist)
  
  datanames=datnames  #for passing to function
  datanames[1]="N"
  
  
  alldat=list()
  for (i in Clans) {
    ldat=list()
    n=1
    for (j in datnames) {
      ldat[[n]]=assign(paste0(j,i),datlist[[j]])
      n=n+1
    }
    names(ldat)=paste0(datnames,i)
    assign(paste0("dat",i),ldat)
  }
  alldat=list(datA=datA,datB=datB,datC=datC,datD=datD,datE=datE)  
  
  
  # Iterations
  
  for (zz in 1:burn.in) {
    
    clan=sample(Clans,1)
    dat=alldat[[which(Clans==clan)]] 
    datSnames=names(dat)
    
    output<-update.network(dat)
    if(zz==1) {mutantA=output$mutant
    clanMutA=clan}    
    if(zz==2) {mutantB=output$mutant
    clanMutB=clan} 
    
    dat=output$datout
    names(dat)=datSnames
    assign(paste0("dat",clan),dat)
    alldat[[which(Clans==clan)]]=dat
    names(alldat[[which(Clans==clan)]])=datSnames
    
    
    SexNew=output$SexNew
    nameNewM=output$nameNewM
    GeN=output$GeN
    deadClan=output$deadClan
    
    if (SexNew==0) {
      
      clanM=sample(Clans[-which(Clans==clan)],1)
      dat=alldat[[which(Clans==clanM)]]
      datSnames=names(dat)
      
      outMale=update.male(dat,nameNewM, GeN)
      
      dat=outMale$datMout
      names(dat)=datSnames
      assign(paste0("dat",clanM),dat)
      
      alldat[[which(Clans==clanM)]]=dat
      names(alldat[[which(Clans==clanM)]])=datSnames
    }
    
    
    if (deadClan==TRUE) {
      clandel=which(Clans==clan)
      alldat[[clandel]]=NULL
      Clans=Clans[-which(Clans==clan)]
    }
  }
  
  alldat=list(datA=datA,datB=datB,datC=datC,datD=datD,datE=datE)
  
  out[k,1]=mutantA
  out[k,2]=clanMutA
  out[k,3]=mutantB
  out[k,4]=clanMutB
  totMutA=0
  totMutB=0
  TotalN=0
  for (i in 1:5) {
    out[k,4*i+1]=alldat[[i]]$Nt
    TotalN=TotalN+alldat[[i]]$Nt
    out[k,4*i+2]=alldat[[i]]$Male
    mutA=0
    mutB=0
    if(alldat[[i]]$Male>1) {    
      mutA=length((which((alldat[[i]]$GenM[,1]+alldat[[i]]$GenM[,2])>0)))
      mutB=length((which((alldat[[i]]$GenM[,3]+alldat[[i]]$GenM[,4])>0)))
    }
    
    if(alldat[[i]]$Male==1){
      if (sum(alldat[[i]]$GenM[1:2])>0) {mutA=1}
      if (sum(alldat[[i]]$GenM[3:4])>0) {mutB=1} 
    }
    
    
    if(alldat[[i]]$Nt>alldat[[i]]$Male+1) {
      mutA=mutA+length((which((alldat[[i]]$GenF[,1]+alldat[[i]]$GenF[,2])>0)))
      mutB=mutB+length((which((alldat[[i]]$GenF[,3]+alldat[[i]]$GenF[,4])>0)))
    }
    
    if(alldat[[i]]$Nt-alldat[[i]]$Male+1==1){
      mutA=mutA+ifelse(sum(alldat[[i]]$GenF[1:2])>0,1,0)
      mutB=mutB+ifelse(sum(alldat[[i]]$GenF[3:4])>0,1,0)
    }
    out[k,4*i+3]=mutA
    out[k,4*i+4]=mutB
    
    totMutA=totMutA+mutA  
    totMutB=totMutB+mutB
  }
  
  out[k,25]=TotalN 
  out[k,26]=burn.in
  out[k,27]=k 
  out[k,28]=totMutA
  out[k,29]=totMutA/TotalN
  out[k,30]=totMutB
  out[k,31]=totMutB/TotalN
 
}#end of k


write.table(out,"Two mutations MaleR_Birth Sum.csv",col.names=TRUE,row.names=F,sep=",",append=T)

#************************************************************************************************

# Two mutations - Scenario IV


library(Matrix)

update.network <- function(datin) {
  
  for (j in 1:length(datin)){
    assign(datanames[j],datin[[j]])
  }
  
  SexNew=1
  nameNewM=0
  deadClan=FALSE
  GeN=rep(0,4)
  repro=0
  
  if(N>Male & Male>0) {
    
    # One individual dies
    
    del=sample(1:N, 1, prob=Die)
    ndel=Name[del]
    Name=Name[-del]
    Die=Die[-del]
    up=which(Rank>Rank[del])
    Rank[up]=Rank[up]-1
    Rank=Rank[-del]
    
    
    if (Sex[del]==1) {
      if (N-Male==1) {
        N=Male
        Founder=NULL
        Kids=NULL
        GenF=NULL
        deadClan=TRUE  
      }
      
      if (N-Male>1) {
        Founder=Founder[-del]
        Kids=Kids[-del]
        GenF=GenF[-del,]
      }
    }
    
    
    if(Sex[del]==0){
      if(Male>1) {GenM=GenM[-(del-N+Male),]}
      else {GenM=NULL}
      Male=Male-1
      if (Male==0) {N=N-1}
    }
    Sex=Sex[-del]
    
    if (deadClan==FALSE & Male>0) {
      
      ChanceF=N-Rank[1:(N-1-Male)]
      ChanceM=1.5*N-Rank[(N-Male):(N-1)]
      
      # Select female & male for reproductio
      
      repro=ifelse(zz==1,sample(which(Rank %in% 25:30 & Sex==1),1), ifelse(zz==2,sample(which(Rank %in% 1:6),1),sample(c(1:(N-1-Male)),1, prob=ChanceF)))
      mate=sample( c(1:Male),1,prob=ChanceM)
      Pa=Name[N-Male-1+mate]
      Kids[repro]=Kids[repro]+1
      
      NameIn=NameIn+1
      
      if (zz==1) {
        GenF[repro,1:2]=c(1,0)
        Names[Name[repro],"Phen"]=1 #phenotype
        mutant=Name[repro]
      }
      
      if (zz==2) {
        GenF[repro,3]=1
        Names[as.numeric(Name[repro]),"Phen"]=2
        mutant=Name[repro]
      }
      
      GeN[1]=ifelse(length(GenF)>4,sample(GenF[repro,1:2],1),sample(GenF[1:2],1))
      GeN[2]=ifelse(length(GenM)>4,sample(GenM[mate,1:2],1),sample(GenM[1:2],1))
      GeN[3]=ifelse(length(GenF)>4,sample(GenF[repro,3:4],1),sample(GenF[3:4],1))
      GeN[4]=ifelse(length(GenM)>4,sample(GenM[mate,3:4],1),sample(GenM[3:4],1))
      
      SexNew=sample(c(0,1),1)
      
      
      if (SexNew==0) {
        N=N-1
        nameNewM=paste0(NameIn,clan)
        Names[NameIn,"Name"]=nameNewM
      }
      
      if(SexNew==1) {
        Name[(N-Male+1):N]=Name[(N-Male):(N-1)]
        Name[N-Male]=NameIn
        Names[NameIn,"Name"]=NameIn
        Sex[N-Male]=1
        if (Male>0) {Sex[N]=0}
        down=which(Rank>Rank[repro])
        Rank[down]=Rank[down]+1
        Rank[(N-Male+1):(N)]=Rank[(N-Male):(N-1)]
        Rank[N-Male]=Rank[repro]+1
        
        Founder[N-Male]=Founder[repro]
        Kids[N-Male]=0
        
        GenF=rbind(GenF,GeN)
        rownames(GenF) <- NULL
      }
      
      
      if(N-Male>1) {
        for (j in 1:(N-Male)) {
          Die[j]=ifelse(sum(GenF[j,])>0,0.5*(N+Rank[j]),N+Rank[j])
        }}
      if(N-Male==1) {Die[1]=ifelse(sum(GenF)>0,0.5*(N+1),N+1)}
      
      if(Male>1) {
        for (j in (N-Male+1):N) {
          Die[j]=ifelse(sum(GenM[j-N+Male,])>0,0.5*(N+Rank[j]),N+Rank[j])
        }}
      if(Male==1) {Die[N]=ifelse(sum(GenM)>0,N,2*N)
      }
      
      Names[NameIn,2:4]=c(SexNew,Name[repro],Founder[repro])
      Names[as.numeric(Name[repro]),"Kids"]=Kids[repro]
      if(Kids[repro]<31) {Names[as.numeric(Name[repro]),8+Kids[repro]]=NameIn}
      if(SexNew==1) {Names[NameIn,"Initial rank"]=Rank[N-Male]}
      Names[NameIn,"Died"]=ndel
      Names[NameIn, "Father"]=Pa
      
      Names[NameIn,"Phen"]=ifelse((sum(GeN[1:2]>0 & sum(GeN[3:4])==0))>1,1, 
      ifelse((sum(GeN[1:2])>0 & sum(GeN[3:4]))>0,3, ifelse((sum(GeN[1:2])==0 & 
      sum(GeN[3:4])>0),2,0)))
      
      
    }}
  
  datout=list(N,Sex, Male, Die, Kids, Rank, Founder, Name,NameIn, Names ,GenM,GenF)
  
  if (zz<3) {return(list(datout=datout,SexNew=SexNew,nameNewM=nameNewM,GeN=GeN, mutant=mutant,deadClan=FALSE, RankM=Rank[repro]))} 
  else {
    return(list(datout=datout, SexNew=SexNew,nameNewM=nameNewM,GeN=GeN, deadClan=deadClan,RankM=Rank[repro]))
  }
}


update.male=function(datin, nameNewM, GeN) {
  
  for (j in 1:length(datin)){
    assign(datanames[j],datin[[j]])
  }
  
  N=N+1
  Male=Male+1
  Name[N]=nameNewM
  Sex[N]=0
  down=which(Rank>RankM+N-Male)
  Rank[down]=Rank[down]+1
  Rank[N]=ifelse(length(down)>0,RankM+N-Male+1,N)
  GenM=rbind(GenM,GeN)
  rownames(GenM)=NULL
  
  if (length(down)>0){
    for (j in ((RankM+1+N-Male):N)) {
      l=which(Rank==j)
      Die[l]=ifelse(sum(GenM[l-N+Male,])>0,0.5*(N+j),N+j)
    }}
  if (length(down)==0) {Die[N]=ifelse(sum(GeN)>0,N,2*N)}
  
  
  datout=list(N,Sex, Male, Die, Kids, Rank, Founder, Name,NameIn, Names ,GenM,GenF)
  
  return(list(datMout=datout))
}


# Simulation

N      <- 60
n.rep         <-200
burn.in       <-3000

out=matrix(0,n.rep,31)
colnames(out)=c("MutantA.Rank","MutantA clan","MutantB.Rank","MutantB clan","NA","MaleA","MutantsA_A",
                "MutantsB_A","NB","MaleB","MutantsA_B","MutantsB_B","NC","MaleC","MutantsA_C","MutantsB_C",
                "ND","MaleD","MutantsA_D","MutantBs_D","NE","MaleE","MutantsA_E","MutantsB_E","Total N",
                "N.iterattions","simulation","No.MutantsA","Prop.MutantsA","No.MutantsB","Prop.MutantsB")


for (k in 1:n.rep) {
  
  N=60
  Male=30
  Sex=c(rep(1,N-Male),rep(0,Male)) #1-female 0-male
  Rank=c(1:N)
  Die=Rank+N
  Kids=rep(0,(N-Male))
  Founder=c(1:(N-Male))
  
  GenF=matrix(0,N-Male,4)
  GenM=matrix(0,Male,4)
  
  Name=c(1:N)
  NameIn=N
  
  A=matrix(0,burn.in+N,8)
  colnames(A)=c("Name","Sex","Mother","Founder","Father","Kids","Initial rank","Died")
  B=matrix(0,burn.in+N,30)
  colnames(B)=paste0("Kid ",c(1:30))
  Phen=rep(0,burn.in+N)
  Names=cbind(A,B,Phen)
  
  Names[1:N,1]=Name
  Names[1:N,"Sex"]=Sex
  Names[1:(N-Male),"Founder"]=Founder
  Names[1:N,"Initial rank"]=Rank
  
  Clans= c("A","B","C","D","E")
  
  
  datlist=list(N,Sex, Male, Die, Kids, Rank, Founder, Name,NameIn,Names ,GenM,GenF)
  
  names(datlist)=c("Nt","Sex", "Male", "Die", "Kids", "Rank", "Founder", "Name","NameIn","Names","GenM","GenF")
  datnames=names(datlist)
  
  datanames=datnames  #for passing to function
  datanames[1]="N"
  
  
  alldat=list()
  for (i in Clans) {
    ldat=list()
    n=1
    for (j in datnames) {
      ldat[[n]]=assign(paste0(j,i),datlist[[j]])
      n=n+1
    }
    names(ldat)=paste0(datnames,i)
    assign(paste0("dat",i),ldat)
  }
  alldat=list(datA=datA,datB=datB,datC=datC,datD=datD,datE=datE)  
  
  
  # Iterations
  
  for (zz in 1:burn.in) {
    
    clan=sample(Clans,1)
    dat=alldat[[which(Clans==clan)]] 
    datSnames=names(dat)
    
    output<-update.network(dat)
    if(zz==1) {mutantA=output$mutant
    clanMutA=clan}    
    if(zz==2) {mutantB=output$mutant
    clanMutB=clan} 
    
    dat=output$datout
    names(dat)=datSnames
    assign(paste0("dat",clan),dat)
    alldat[[which(Clans==clan)]]=dat
    names(alldat[[which(Clans==clan)]])=datSnames
    
    
    SexNew=output$SexNew
    nameNewM=output$nameNewM
    GeN=output$GeN
    deadClan=output$deadClan
    RankM=output$RankM
    
    if (SexNew==0) {
      
      clanM=sample(Clans[-which(Clans==clan)],1)
      dat=alldat[[which(Clans==clanM)]]
      datSnames=names(dat)
      
      outMale=update.male(dat,nameNewM, GeN)
      
      dat=outMale$datMout
      names(dat)=datSnames
      assign(paste0("dat",clanM),dat)
      
      alldat[[which(Clans==clanM)]]=dat
      names(alldat[[which(Clans==clanM)]])=datSnames
    }
    
    
    if (deadClan==TRUE) {
      clandel=which(Clans==clan)
      alldat[[clandel]]=NULL
      Clans=Clans[-which(Clans==clan)]
    }
  }
  
  alldat=list(datA=datA,datB=datB,datC=datC,datD=datD,datE=datE)
  
  out[k,1]=mutantA
  out[k,2]=clanMutA
  out[k,3]=mutantB
  out[k,4]=clanMutB
  totMutA=0
  totMutB=0
  TotalN=0
  for (i in 1:5) {
    out[k,4*i+1]=alldat[[i]]$Nt
    TotalN=TotalN+alldat[[i]]$Nt
    out[k,4*i+2]=alldat[[i]]$Male
    mutA=0
    mutB=0
    if(alldat[[i]]$Male>1) {    
      mutA=length((which((alldat[[i]]$GenM[,1]+alldat[[i]]$GenM[,2])>0)))
      mutB=length((which((alldat[[i]]$GenM[,3]+alldat[[i]]$GenM[,4])>0)))
    }
    
    if(alldat[[i]]$Male==1){
      if (sum(alldat[[i]]$GenM[1:2])>0) {mutA=1}
      if (sum(alldat[[i]]$GenM[3:4])>0) {mutB=1} 
    }
    
    if(alldat[[i]]$Nt>alldat[[i]]$Male+1) {
      mutA=mutA+length((which((alldat[[i]]$GenF[,1]+alldat[[i]]$GenF[,2])>0)))
      mutB=mutB+length((which((alldat[[i]]$GenF[,3]+alldat[[i]]$GenF[,4])>0)))
    }
    
    if(alldat[[i]]$Nt-alldat[[i]]$Male+1==1){
      mutA=mutA+ifelse(sum(alldat[[i]]$GenF[1:2])>0,1,0)
      mutB=mutB+ifelse(sum(alldat[[i]]$GenF[3:4])>0,1,0)
    }
    out[k,4*i+3]=mutA
    out[k,4*i+4]=mutB
    
    totMutA=totMutA+mutA  
    totMutB=totMutB+mutB
  }
  
  out[k,25]=TotalN 
  out[k,26]=burn.in
  out[k,27]=k 
  out[k,28]=totMutA
  out[k,29]=totMutA/TotalN
  out[k,30]=totMutB
  out[k,31]=totMutB/TotalN
  
}#end of k


write.table(out,"Two mutations MaleR_Ma Sum.csv",col.names=TRUE, row.names=F,sep=",",append=T)


