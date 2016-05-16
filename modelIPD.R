if(!exists('baseSplit')){
  source('readContigs.R')
  #make sure have 5 positions around
  contigs$plusMinus5<-ave(contigs$Pos,contigs$ID,contigs$Strand,FUN=function(x)lagNA(x,5,-9999)==x+5&lagNA(x,-5,-9999)==x-5)==1
  kmers<-do.call(cbind,cacheOperation('work/kmers.Rdat',lapply,1:4,function(kk,bases){
    message(kk)
    offsets<-(-kk+1):0
    out<-do.call(cbind,lapply(offsets,function(offset){
      message(offset)
      do.call(paste,c(lapply(offset:(offset+kk-1),function(lag)lagNA(bases,lag,'X')),list(sep='')))
    }))
    colnames(out)<-sprintf('%d_%d',kk,offsets)
    return(out)
  },contigs$Base))

  bases<-c('A'='A','C'='C','G'='G','T'='T')
  splitIds<-split(1:nrow(contigs),contigs$ID)
  baseSplit<-lapply(bases,function(base){
    selector<-contigs$Base==base&contigs$plusMinus5 
    info<-contigs[selector,c('ID','Pos','Strand','Base','IPDRatio')]
    thisKmers<-kmers[selector,]
    allMotifs<-unlist(lapply(colnames(thisKmers),function(x){sprintf('%s_%s',x,unique(thisKmers[,x]))}))
    motifIds<-structure(1:length(allMotifs),.Names=allMotifs)
    thisKmers[,]<-sprintf('%s_%s',rep(colnames(thisKmers),each=nrow(thisKmers)),thisKmers)
    kmerN<-apply(thisKmers,2,function(x)motifIds[x])
    return(list('info'=info,'kmer'=kmerN,'motifs'=allMotifs))
  })

  rm(contigs)
  rm(kmers)
  rm(seqs)
  gc()
}



#kmerMat nMotif x nPos
#ipdLikeByContig<-function(ipds,kmerMat,motifPs,mus,sigmas,baseMu,baseSigma){
  #pMod<-kmerMat*motifPs
  #pNotMod<-1-apply(pMod,2,sum)
  #xy<-which(pMod>0,arr.ind=TRUE)
  #indexMatrix(xy$row,xy$col,pMod)*dnorm(ipds[xy$col
  ##NEED TO DEAL WITH p(motif1)=.9 p(motif2)=.9 ...  interactions are a pain reformatting
#}

#kmerMat nMotif + 1 x nPos with first as intercept
#,sigma2s
ipdLikeByContig<-function(ipds,kmers,motifActives,mus,baseMu,baseSigma){
  activeMu<-motifActives*mus
  #activeSigmas<-motifActives*sigma2s
  #muAdds<-motifActives %*% mus
  #system.time(muAdds<-apply(kmers,1,function(x)sum(activeMu[x])))
  #system.time(muAdds<-apply(matrix(activeMu[kmers],nrow=nrow(kmers)),1,sum))
  muAdds<-apply(matrix(activeMu[kmers],nrow=nrow(kmers)),1,sum)
  #sigmaAdds<-apply(kmers,1,function(x)sum(activeSigmas[x]))
  #sigmaAdds <-motifActives %*% sigma2s
  ps<-dnorm(ipds,baseMu+muAdds,baseSigma,log=TRUE)
  return(sum(ps))
}
ipdLike<-function(ipds,kmers,splitIds,activeMotifs,...){
  lapply(1:length(splitIds),function(ii,ipds,kmers,...){
    cat('.')
    ipdLikeByContig(ipds[splitIds[[ii]]],kmers[splitIds[[ii]],],activeMotifs[ii,],...)
  },ipds,kmers,...)
}

analyzeBase<-function(info,kmer,allMotifs){
  mus<-structure(rep(0,length(allMotifs)),.Names=allMotifs)
  #kmer[,]<-as.numeric(motifIds[sprintf('%s_%s',rep(colnames(kmer),each=nrow(kmer)),kmer)])
  baseMu<-0
  baseSigma<-2
  splitIds<-split(1:nrow(info),info$ID)
  activeMotifs<-matrix(FALSE,nrow=length(splitIds),ncol=length(allMotifs),dimnames=list(names(splitIds),allMotifs))
  browser()
  system.time(ipdLike(info$IPDRatio,kmer,splitIds,activeMotifs,mus,baseMu,baseSigma))
}
analyzeBase(baseSplit[[1]]$info,baseSplit[[1]]$kmer,baseSplit[[1]]$motifs)

#need dummy value to avoid model.matrix complain about singular columns
#kmerMat<-model.matrix(~.,as.data.frame(rbind(rep('ZZZZ',ncol(baseSplit[[1]]$kmer)),baseSplit[[1]]$kmer[,])))
#kmerMat<-kmerMat[-1,!grepl('ZZZZ',colnames(kmerMat))]
#maybe dont need kmer matrix just use indexing

#emAlgo<-function(info,kmer,allMotifs){
  #motifLookup<-lapply(1:length(allMotifs),function(xx)which(kmer==xx,arr.ind=TRUE)[,'row'])
  #newMu<-jj
#}
#newMu<-function(ipd,selector,baseMu){
 #
#}

#ipdLikeByContig(baseSplit[[1]]$info$IPD,z,

#library('Matrix')

#zz<-Matrix(z,sparse=TRUE)
#object.size(z)
#object.size(zz)




