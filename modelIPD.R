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
    kmers<-kmers[selector,]
    return(list('info'=info,'kmer'=kmers))
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
ipdLikeByContig<-function(ipds,kmers,motifActives,mus,sigma2s){
  pMod<-kmerMat*motifActives
  muSum<-apply(pMod,1,function(x)sum(mus[x]))
  sigma2Sum<-apply(pMod,1,function(x)sum(sigma2s[x]))
  return(sum(dnorm(ipds,muSum,sigma2Sum)))
}


#need dummy value to avoid model.matrix complain about singular columns
#kmerMat<-model.matrix(~.,as.data.frame(rbind(rep('ZZZZ',ncol(baseSplit[[1]]$kmer)),baseSplit[[1]]$kmer[,])))
#kmerMat<-kmerMat[-1,!grepl('ZZZZ',colnames(kmerMat))]
#maybe dont need kmer matrix just use indexing


#ipdLikeByContig(baseSplit[[1]]$info$IPD,z,

#library('Matrix')

#zz<-Matrix(z,sparse=TRUE)
#object.size(z)
#object.size(zz)




