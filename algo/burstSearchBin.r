#!/bin/env Rscript
setwd('/share2/data/liuz/smFRET/48diUb_20170228')
padz<- 1 #连续出现次数
phT<-  25
lines2Proc<- 0 #最大处理行数,小于maxWindowTime时处理所有行
threadL<-c(5.1,4.1,3.1,9.1) #delete v< threadL 单个bin中的最小光子数
rawDatafileName<-"LS150_FRET.dat" #psFRET.dat是stSum
addPaste <- function(x) Reduce("paste", x)
strThr<-addPaste(threadL)
fretfig<-paste("fretFig",strThr,padz,".pdf",sep="")

burstsFileName<-paste("LS150_65410_25_",padz,".dat",sep="")
burstsSumFileName<-paste("LS150_Sum_65410_25_",padz,".dat",sep="")
maxWindowTime<-20 #ms



blockSiz<-20000
library(compiler)
unlink(burstsFileName)
unlink(burstsSumFileName)
findruns<-function(x,k,tt){
  #print(x)
  #print(tt)
  n<-nrow(x)
  r<-ncol(x)
  runs<-NULL
  i=1
  pp<-(t(x)<tt)
  pp<-apply(pp, 2, any)
  while(i <= (n-k+1)){    
    #print(i)
    if(all( pp[i:(i+k-1)])){       
      ii<-0
      while(pp[i+k+ii]&&i+k+ii<=n){
        ii<-ii+1
      }
      runs<-rbind(runs,c(i,ii+k))
      i<-i+ii+k
    }
    else
      i<-i+1    
  }
  return(runs)
}

c_findruns<-cmpfun(findruns)

judgePhotonSum<-function(ichuck,rdata,pos0,photons,deCol){
  n<-nrow(pos0)
  burst<-NULL
  burstSum<-NULL
  for(i in ichuck){
    if(i+1<=n){
      photonSum<-sum(rdata[(pos0[i,1]+pos0[i,2]):(pos0[i+1,1]-1),deCol])
      if(photonSum>=photons){
	phoColSum<-NULL
	if(deCol-1>=2){
            cols<-2:(deCol-1)	    
            for(ic in cols){
		tphotonSum<-sum(rdata[(pos0[i,1]+pos0[i,2]):(pos0[i+1,1]-1),ic])
		phoColSum<-cbind(phoColSum,tphotonSum)
	    }
	}
        #DphotonSum<-sum(rdata[(pos0[i,1]+pos0[i,2]):(pos0[i+1,1]-1),2])
        #AphotonSum<-sum(rdata[(pos0[i,1]+pos0[i,2]):(pos0[i+1,1]-1),3])
        burst<-rbind(burst,rdata[(pos0[i,1]+pos0[i,2]):(pos0[i+1,1]-1),])
        burstSum<-rbind(burstSum,c(rdata[(pos0[i,1]+pos0[i,2]),1],rdata[pos0[i+1,1]-1,1],phoColSum,photonSum))
      }
    }
  }  
  return(list(burst,burstSum))
}
findBursts<-function(cl,bdata,br,ncols=5,deCol=5,pad0=2,photons=25,thread=c(1,0,2)){    
   lastL<-1
   photonSum<-0
   br$found<-F
   #print(bdata)
   #bdata<-matrix(bdata,ncol=ncols,byrow=T)
   
   pos0<- c_findruns (bdata[,2:ncols],pad0,thread)
   #print(pos0)
   lenprow<-nrow(pos0)
   if(lenprow>0){
     br$nextStartL<-pos0[,1][lenprow]+pos0[,2][lenprow]+br$nextStartL-1       
     if(pos0[1,1]!=1){       
       photonSum<-sum(bdata[1:pos0[1,1]-1,deCol])       
       if(photonSum>=photons){  
	 phoColSum<-NULL       
         #DphotonSum<-sum(bdata[1:pos0[1,1]-1,2])       
         #AphotonSum<-sum(bdata[1:pos0[1,1]-1,3])       
	if(deCol-1>=2){
            cols<-2:(deCol-1)	    
            for(ic in cols){
		tphotonSum<-sum(bdata[1:pos0[1,1]-1,ic])
		phoColSum<-cbind(phoColSum,tphotonSum)
	    }
	}
         br$burst<-rbind(br$burst,bdata[1:pos0[1,1]-1,]) 
         br$burstSum<-rbind(br$burstSum,c(bdata[1,1],bdata[pos0[1,1]-1,1],phoColSum,photonSum))
         br$found<-T
       }
     }
     n<-nrow(pos0)-1
     nc<-length(cl)
     options(warn=-1)
     ichunks<-split(1:n,1:nc)
     options(warn=0)
     burstA<-clusterApply(cl,ichunks,judgePhotonSum,bdata,pos0,photons,deCol)
     
     for (burst in burstA){
       if(!is.null(burst[[1]])){
         #print(burst[[1]])
         br$burst<-rbind(br$burst,burst[[1]])  
         #print("=======")
         br$burstSum<-rbind(br$burstSum,burst[[2]])
         br$found<-T
       }
     }
     
   }
   else{
     photonSum<-sum(bdata[,deCol])
     if(photonSum>=photons){
       #DphotonSum<-sum(bdata[,2])
       #AphotonSum<-(bdata[,3])
	phoColSum<-NULL
	if(deCol-1>=2){
            cols<-2:(deCol-1)	    
            for(ic in cols){
		tphotonSum<-sum(bdata[,cols])
		phoColSum<-cbind(phoColSum,tphotonSum)
	    }
	}
       br$burst<-rbind(br$burst,bdata)
       br$burstSum<-rbind(br$burstSum,c(bdata[1,1],bdata[nrow(bdata),1],phoColSum,photonSum))
       br$found<-T
     }       
     br$nextStartL<-nrow(bdata)+br$nextStartL
   }
  return(br)
  
}
c_findBursts<-cmpfun(findBursts)

startTime <- Sys.time()
br<-list(found=F,lines=0,startLine=1,burst=NULL,burstSum=NULL,nextStartL=1)  
rawdata<-NULL
nlinesR <- 0
if(lines2Proc>maxWindowTime)
  nlinesR<-lines2Proc
rawdata<-scan(rawDatafileName,nlines=nlinesR,skip=1,quiet=T,what=numeric())
rdata<-matrix(rawdata,ncol=5,byrow=T)
rdataRow<-nrow(rdata)

library(snow)
cl<-makeCluster(type="SOCK",c("localhost","localhost","localhost","localhost","localhost","localhost","localhost","localhost"))
while(br$nextStartL<rdataRow){
  trdata<-rdata[br$nextStartL:min(rdataRow,br$nextStartL+blockSiz),]
  #print(trdata)
  #print(br$nextStartL+blockSiz)
  br<-c_findBursts(cl,trdata,br,pad0=padz,thread=threadL,photons=phT)
  if(br$found){       
    write(t(br$burst),file=burstsFileName,sep="\t",ncolumns=5,append=TRUE)   
    write(t(br$burstSum),file=burstsSumFileName,sep="\t",ncolumns=6,append=TRUE)   
    br$burst<-NULL
    br$burstSum<-NULL
  }
  br$found<-F
  #br$nextStartL<-br$nextStartL+maxWindowTime*2  
}
Sys.time() - startTime
stopCluster(cl)
#E_fret<-NULL
#fE_fret<-scan(burstsSumFileName,quiet=T,what=numeric())
#fE_fret<-matrix(fE_fret,ncol=6,byrow=T)
#nE_fret<-nrow(fE_fret)
#E_fret<-seq(nE_fret)
#for(lE in seq(nE_fret)){
#  E_fret[lE]<-fE_fret[lE,4]/fE_fret[lE,6]
#}
#pdf(fretfig)
#devl<-dev.list()
#dev.set(devl["pdf"])
#hist(E_fret,500)
#dev.off()
