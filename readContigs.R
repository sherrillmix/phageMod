library(dnar)
library(parallel)
contigFiles<-list.files('data/','^Normal.*txt.gz$',full.names=TRUE)
contigs<-do.call(rbind,cacheOperation('work/allData.Rdat',lapply,contigFiles,read.table,header=TRUE,stringsAsFactors=FALSE))

posBaseToSeq<-function(position,base,strand){
  positionOrder<-order(position)
  position<-position[positionOrder]
  base<-base[positionOrder]
  strand<-strand[positionOrder]
  stranded<-rep(paste(rep('-',max(position)),collapse=''),2)
  names(stranded)<-c('+','-')
  isPos<-strand=='+'
  posRanges<-index2range(position[isPos])
  negRanges<-index2range(position[!isPos])
  #ordered correctly above
  concatPos<-paste(base[isPos],collapse='')
  concatNeg<-paste(base[!isPos],collapse='')
  if(nrow(posRanges)>0){
    posRanges$width<-posRanges$end-posRanges$start+1
    posRanges$seqStart<-cumsum(c(0,posRanges$width[-nrow(posRanges)]))+1
    for(ii in 1:nrow(posRanges))substring(stranded['+'],posRanges[ii,'start'],posRanges[ii,'end'])<-substring(concatPos,posRanges[ii,'seqStart'],posRanges[ii,'seqStart']+posRanges[ii,'width']-1)
  }
  if(nrow(negRanges)>0){
    negRanges$width<-negRanges$end-negRanges$start+1
    negRanges$seqStart<-cumsum(c(0,negRanges$width[-nrow(negRanges)]))+1
    for(ii in 1:nrow(negRanges))substring(stranded['-'],negRanges[ii,'start'],negRanges[ii,'end'])<-substring(concatNeg,negRanges[ii,'seqStart'],negRanges[ii,'seqStart']+negRanges[ii,'width']-1)
  }
  #fix info missing in one strand
  stranded['-']<-complementDna(stranded['-'])
  out<-paste(apply(do.call(rbind,strsplit(stranded,'')),2,function(x){
      out<-unique(x[x!='-'])
      if(length(out)==0)return('-')
      if(length(out)>1)stop(simpleError('Problem rectifying strands'))
      return(out)}
  ),collapse='')
  #too slow
  #mapply(function(x,y){if(runif(1)<.001)cat('.');substring(stranded['+'],y,y)<<-x},base[isPos],position[isPos])
  #mapply(function(x,y)substring(stranded['-'],y,y)<<-x,base[!isPos],position[!isPos])
  #if(stranded['+']!=complementDna(stranded['-']))stop(simpleError('Stranded contigs do not match'))
  return(out)
}

splitIds<-split(1:nrow(contigs),contigs$ID)
seqs<-cacheOperation('work/seqs.Rdat',lapply,splitIds,function(x){cat('.');posBaseToSeq(contigs[x,'Pos'],contigs[x,'Base'],ifelse(contigs[x,'Strand']=='Fwd','+','-'))})
write.fa(names(seqs),seqs,'work/seqs.fa')
