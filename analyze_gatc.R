library(vipor)
library(dnar)
library(parallel)
gatc<-read.csv('data/GATCs_per_contig_with_expected_and_mod_status.csv',stringsAsFactors=FALSE)
#pdf('test.pdf')
#hist(gatc$total_motif_occurrences-gatc$expected_GATCs,breaks=200)
#dev.off()

#mod<-glm(I(total_motif_occurrences==0)~I(data_set=='Phage_25'),offset=log(expected_GATCs),data=gatc[!grepl('T4',gatc$data_set),],family='binomial')
#print(summary(mod))
#mod<-glm(I(total_motif_occurrences==0)~I(data_set=='Phage_25')+log(expected_GATCs),data=gatc[!grepl('T4',gatc$data_set),],family='binomial')
#print(summary(mod))

#mod<-glm(total_motif_occurrences~data_set,offset=log(expected_GATCs),data=gatc[!grepl('T4',gatc$data_set),],family='poisson')
mod<-glm(total_motif_occurrences~I(data_set=='Phage_25'),offset=log(expected_GATCs),data=gatc[!grepl('T4',gatc$data_set),],family='poisson')
#print(summary(mod))

#mod<-glm(I(total_motif_occurrences==0)~data_set,offset=log(expected_GATCs),data=gatc[!grepl('T4',gatc$data_set),],family='binomial')
#print(summary(mod))

#library(pscl)
inv.logit<-plogis
logit<-qlogis

seqs<-read.fa('work/seqs.fa')
rownames(seqs)<-seqs$name
with(seqs[gatc[logit(gatc$binomProb)< -50,'contig'],],write.fa(name,seq,'work/gatcDepleted.fa'))

selector<-gatc$contig %in% rownames(seqs)
warning('Throwing out ',sum(!selector),' contigs for not being in sequences')
gatc<-gatc[selector,]
rownames(seqs)[rownames(seqs) %in% gatc$contig][1:10]
warning('Throwing out ',sum(!rownames(seqs) %in% gatc$contig),' contigs for not being in sequences')
seqs<-seqs[gatc$contig,]
gatc$nGATC<-sapply(gregexpr('(?=GATC)',seqs$seq,perl=TRUE),function(x)sum(x!=-1))
gatc$length<-nchar(degap(seqs$seq))
gatc$gappy<-grepl('-',sub('-+$','',sub('^-+','',seqs$seq)))
gatc$gc<-gcPercent(degap(seqs$seq))
# don't forget 2* for forward and reverse complement if not palindrome
gatc$prob<-(gatc$gc/2)^2*((1-gatc$gc)/2)^2 
gatc$totalProb<-gatc$contig_len*log(1-gatc$prob)
gatc$logitProb<-gatc$totalProb-log(1-exp(gatc$totalProb))
gatc$binomProb<-pbinom(gatc$total_motif_occurrences,gatc$contig_len,gatc$prob)
#mod<-glm(I(total_motif_occurrences==0)~data_set,offset=logitProb,data=gatc[!grepl('T4',gatc$data_set),],family='binomial')

pdf('out/gatc.pdf')
  par(las=2,mar=c(7,4.5,.5,.5))
  with(gatc[!grepl('T4',gatc$data_set)&gatc$total_motif_occurrences==0,],vpPlot(data_set,totalProb*log10(exp(1)),ylab='Probability of no GATC log10(p)',offsetXArgs=list(varwidth=TRUE)))
  with(gatc[!grepl('T4',gatc$data_set),],vpPlot(data_set,logit(binomProb),ylab='Probability <= GATC count logit(p)',offsetXArgs=list(varwidth=TRUE)))
  abline(h=0,lty=2,col='#FF000099')
dev.off()



dna<-c('A','C','T','G')
twoMers<-as.vector(outer(dna,dna,paste,sep=''))
fourMers<-as.vector(outer(twoMers,twoMers,paste,sep=''))
names(fourMers)<-fourMers
palin<-fourMers[fourMers==revComp(fourMers)]

palinCounts<-do.call(cbind,cacheOperation('work/palin.Rdat',mclapply,palin,function(mer){
  message(mer)
  sapply(gregexpr(sprintf('(?=%s)',mer),seqs$seq,perl=TRUE),length)
},mc.cores=8))
rownames(palinCounts)<-seqs$name
colnames(palinCounts)<-palin

fourCounts<-do.call(cbind,cacheOperation('work/4mer.Rdat',mclapply,fourMers,function(mer){
  message(mer)
  sapply(gregexpr(sprintf('(?=%s)',mer),seqs$seq,perl=TRUE),length)
},mc.cores=8))
rownames(fourCounts)<-seqs$name
colnames(fourCounts)<-fourMers

gcProb<-function(seq,gc=.5){
  nAT<-nchar(gsub('[^AT]+','',seq))
  nGC<-nchar(gsub('[^GC]+','',seq))
  p<-(gc/2)^nGC*((1-gc)/2)^nAT
  if(revComp(seq)!=seq)p<-p*2
  return(p)
}

palinP<-do.call(cbind,lapply(1:ncol(palinCounts),function(ii)pbinom(palinCounts[,ii],gatc$contig_len,gcProb(colnames(palinCounts)[ii],gatc$G.C))))
colnames(palinP)<-colnames(palinCounts)
palinP[palinP==0]<-min(palinP[palinP>0])

fourP<-do.call(cbind,lapply(1:ncol(fourCounts),function(ii)pbinom(fourCounts[,ii],gatc$contig_len,gcProb(colnames(fourCounts)[ii],gatc$G.C))))
colnames(fourP)<-colnames(fourCounts)
fourP[fourP==0]<-min(fourP[fourP>0])

fourCI<-do.call(cbind,mclapply(1:ncol(fourCounts),function(ii){
  p<-gcProb(colnames(fourCounts)[ii],gatc$G.C)
  out<-mapply(function(xx,yy,zz)logit(conservativeBoundary(binom.test(xx,yy)$conf.int,zz))-logit(zz),fourCounts[,ii],gatc$contig_len,p)
  return(out)
},mc.cores=12))
colnames(fourCI)<-colnames(fourCounts)

png('out/palin4mers.png',height=1500,width=1500,res=200)
  vpPlot(rep(colnames(palinP),each=nrow(palinP)),logit(as.vector(palinP)),las=2,col=NA,bg='#00000055',pch=21)
dev.off()

png('out/4mers.png',height=1500,width=6000,res=200)
  vpPlot(rep(colnames(fourP),each=nrow(fourP)),logit(as.vector(fourP)),las=2,col=NA,bg='#00000055',pch=21)
dev.off()

png('out/4CI.png',height=1500,width=6000,res=200)
  vpPlot(rep(colnames(fourP),each=nrow(fourP)),as.vector(fourCI),las=2,col=NA,bg='#00000055',pch=21,xaxs='i',ylab='Conservative boundary of predicted - expected (logit)',xlim=c(0,ncol(fourP)+1))
  abline(h=0,lty=2,col='#FF0000')
dev.off()

png('out/4CI_separate.png',height=8000,width=8000,res=200)
  par(mfrow=c(16,16))
  for(ii in sort(colnames(fourCI))){
    message(ii)
    vpPlot(gatc$data_set,fourCI[,ii],las=2,col=NA,bg='#00000055',pch=21,xaxs='i',ylab='Conservative boundary of predicted - expected (logit)',main=ii,ylim=range(fourCI))
    abline(h=0,lty=2,col='#FF0000')
  }
dev.off()


library(pscl)
summary(zeroinfl(total_motif_occurrences~I(data_set=='Phage_25')|I(data_set=='Phage_25'),offset=log(expected_GATCs),data=gatc[!grepl('T4',gatc$data_set),]))


