library(dnar)
isLeft<-function(x1,x2,x3,y1,y2,y3){
  cross<-(x2 - x1) * (y3- y1) - (y2 - y1) * (x3 - x1)
  if(is.na(cross)||abs(cross)<1e-10)return(NA)
  return(cross>0)
}
checkLeft<-function(xs,ys){
  out<-rep(NA,length(xs))
  if(length(xs)!=length(ys))stop(simpleError('xs and ys not same length'))
  if(length(xs)<3)return(out)
  for(ii in 3:length(xs)){
    out[ii]<-isLeft(xs[ii-2],xs[ii-1],xs[ii],ys[ii-2],ys[ii-1],ys[ii]) 
  }
  return(out)
}


phage<-read.csv('data/phage_reads_for_abundance_analysis.csv',stringsAsFactors=FALSE)
bac<-read.table('data/bacteria_16s_reads_for_abundance_analysis.tsv',stringsAsFactors=FALSE)


#linkFiles<-c(
#  'data/Scott_table_p.adjust_phage_bac_pairs_Family_Leveluse_annotated_phage_TRUE_Motif_occur_3_use_unclassified_bac_TRUE_use_GATCs_FALSE.csv',
#  'data/Scott_table_p.adjust_phage_bac_pairs_Genus_Leveluse_annotated_phage_TRUE_Motif_occur_3_use_unclassified_bac_TRUE_use_GATCs_FALSE.csv'
#)
#linkFiles<-list.files('data/link/','*.csv',full.names=TRUE)
linkFiles<-c(
  'select'='data/link/3known_5kmer_occur_table_p.adjust_phage_bac_pairs_Genus_Leveluse_annotated_phage_TRUE_Motif_occur_3_use_unclassified_bac_TRUE_use_GATCs_FALSE.csv',
  'family'='data/link/table_p.adjust_phage_bac_pairs_Family_Leveluse_annotated_phage_FALSE_Motif_occur_3_use_unclassified_bac_TRUE_use_GATCs_TRUE.csv',
  'all'='data/link/all_3mers_occur_table_p.adjust_phage_bac_pairs_Genus_Leveluse_annotated_phage_FALSE_Motif_occur_3_use_unclassified_bac_TRUE_use_GATCs_TRUE.csv'
)

links<-do.call(rbind,lapply(names(linkFiles),function(x){out<-read.csv(linkFiles[x],stringsAsFactors=FALSE);out$set<-x;return(out)}))
bak<-links
links<-links[!grepl('Unclassified|Uncultured|uncultured|unidentified',links$Bacteria),]
links<-links[!grepl('GATC.*\\(A',links$motif_mod_base),]
links$allMotif<-ave(links$motif_mod_base,paste(links$contig,links$Bacteria,links$set,sep='_!#'),FUN=function(x)paste(x,collapse=','))
links<-links[!duplicated(links[,c('contig','Bacteria','set')]),]
#links<-links[,colnames(links)!='X']
#links<-unique(links[!grepl('Unclassified|Uncultured',links$Bacteria),])

pCountCols<-colnames(phage)[grep('^X[0-9.]+$',colnames(phage))]
bCountCols<-colnames(bac)[grep('^X[0-9.]+$',colnames(bac))]
overlapDays<-sort(as.numeric(unique(sub('\\.[0-9]','',sub('^X','',pCountCols)))))
overlapDays<-overlapDays[overlapDays %in% as.numeric(sub('\\.[0-9]','',sub('^X','',bCountCols)))]
pDays<-lapply(overlapDays,function(x)pCountCols[grep(sprintf('^X%d(\\.[0-9]+)?$',x),pCountCols)])
bDays<-lapply(overlapDays,function(x)bCountCols[grep(sprintf('^X%d(\\.[0-9]+)?$',x),bCountCols)])
names(bDays)<-names(pDays)<-overlapDays
oneDiffs<-binary2range(diff(overlapDays)==1)
oneDiffs$end<-oneDiffs$end+1
oneDiffDays<-lapply(split(oneDiffs,1:nrow(oneDiffs)),function(x)overlapDays[x$start:x$end])


phageProp<-apply(phage[,pCountCols],2,function(x)x/sum(x))
bacProp<-apply(bac[,bCountCols],2,function(x)x/sum(x))

bacPhage<-mcmapply(function(bacId,phageId,taxa){
  cat('.')
  thisPhage<-phageProp[phage$Phage_contig==phageId,,drop=FALSE]
  thisBac<-bacProp[bac[,taxa]==bacId&!is.na(bac[,taxa]),,drop=FALSE]
  pOut<-do.call(cbind,lapply(pDays,function(x)apply(thisPhage[,x,drop=FALSE],1,mean)))
  bOut<-do.call(cbind,lapply(bDays,function(x)apply(thisBac[,x,drop=FALSE],1,mean)))
  colnames(pOut)<-colnames(bOut)<-names(bDays)
  return(list('bac'=bOut,'phage'=pOut))
},links$Bacteria,links$contig,links$Bact_Phylo_Level,SIMPLIFY=FALSE,mc.cores=10)

otuCounts<-do.call(rbind,lapply(bacPhage,sapply,nrow))

pVsB<-mclapply(bacPhage,function(x){
  lapply(oneDiffDays,function(days){
    out<-data.frame(
      'phage'=apply(x[['phage']][,as.character(days),drop=FALSE],2,sum),
      'bac'=apply(x[['bac']][,as.character(days),drop=FALSE],2,sum)
    )
    out$isLeft=checkLeft(out$bac,out$phage)
    return(out)
  })
},mc.cores=10)
nLeft<-sapply(pVsB,function(x)sum(unlist(lapply(x,function(y)y[,'isLeft'])),na.rm=TRUE))
nRight<-sapply(pVsB,function(x)sum(unlist(lapply(x,function(y)!y[,'isLeft'])),na.rm=TRUE))
table(nLeft,nRight)
plotOrder<-which(nLeft+nRight==max(nLeft+nRight))
plotOrder<-plotOrder[order(nLeft[plotOrder]-nRight[plotOrder])]
table('left'=nLeft[plotOrder],'right'=nRight[plotOrder])

plotLink<-function(ii){
  cols<-rainbow.lab(3,lightScale=0,lightMultiple=.7,alpha=.6)
    thisInfo<-links[ii,]
    thisPB<-pVsB[[ii]]
    thisCat<-do.call(rbind,thisPB)
    plot(thisCat$bac,thisCat$phage,main=sprintf('%s %s Lefts: %d Rights: %d\nMotif: %s',thisInfo$Bacteria,thisInfo$contig,nLeft[ii],nRight[ii],thisInfo$allMotif),xlab='Bacteria proportion',ylab='Phage proportion',lwd=2,pch=21,bg=cols[rep(1:length(thisPB),sapply(thisPB,nrow))],col=NA,cex=2)
    lapply(thisPB,function(x)arrows(x$bac[-nrow(x)],x$phage[-nrow(x)],x$bac[-1],x$phage[-1],length=.1))
}
indices<-plotOrder[plotOrder %in% which(links$set=='all')]
message('genus (n=',length(indices),')')
print(table('left'=nLeft[indices],'right'=nRight[indices]))
pdf('out/circle.pdf')
  lapply(indices,plotLink)
dev.off()

indices<-plotOrder[plotOrder %in% which(links$set=='select')]
message('select (n=',length(indices),')')
print(table('left'=nLeft[indices],'right'=nRight[indices]))
pdf('out/circle_select.pdf')
  lapply(indices,plotLink)
dev.off()

indices<-plotOrder[plotOrder %in% which(links$set=='family')]
message('family (n=',length(indices),')')
print(table('left'=nLeft[indices],'right'=nRight[indices]))
pdf('out/circle_family.pdf')
  lapply(indices,plotLink)
dev.off()



table(do.call(rbind,lapply(pVsB,function(x)do.call(rbind,x)))$isLeft)

