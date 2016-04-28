gatc<-read.csv('data/GATCs_per_contig_with_expected_and_mod_status.csv',stringsAsFactors=FALSE)
pdf('test.pdf')
hist(gatc$total_motif_occurrences-gatc$expected_GATCs,breaks=200)
dev.off()

mod<-glm(I(total_motif_occurrences==0)~I(data_set=='Phage_25'),offset=log(expected_GATCs),data=gatc[!grepl('T4',gatc$data_set),],family='binomial')
print(summary(mod)) mod<-glm(I(total_motif_occurrences==0)~I(data_set=='Phage_25')+log(expected_GATCs),data=gatc[!grepl('T4',gatc$data_set),],family='binomial')
print(summary(mod))

mod<-glm(total_motif_occurrences~data_set,offset=log(expected_GATCs),data=gatc[!grepl('T4',gatc$data_set),],family='poisson')
mod<-glm(total_motif_occurrences~I(data_set=='Phage_25'),offset=log(expected_GATCs),data=gatc[!grepl('T4',gatc$data_set),],family='poisson')
print(summary(mod))

mod<-glm(I(total_motif_occurrences==0)~data_set,offset=log(expected_GATCs),data=gatc[!grepl('T4',gatc$data_set),],family='binomial')
print(summary(mod))

#library(pscl)
inv.logit<-plogis
logit<-qlogis

# don't forget 2* for forward and reverse complement if not palindrome
gatc$prob<-(gatc$G.C/2)^2*((1-gatc$G.C)/2)^2 
gatc$totalProb<-gatc$contig_len*log(1-gatc$prob)
gatc$logitProb<-gatc$totalProb-log(1-exp(gatc$totalProb))
mod<-glm(I(total_motif_occurrences==0)~data_set,offset=logitProb,data=gatc[!grepl('T4',gatc$data_set),],family='binomial')

pdf('test.pdf')
  par(las=2,mar=c(10,5,.5,.5))
  with(gatc[!grepl('T4',gatc$data_set),],vpPlot(paste(ifelse(total_motif_occurrences==0,'No GATC','GATC'),data_set),log10(-totalProb*log10(exp(1))),ylab='Probability of no GATC log10(-log10(p))',offsetXArgs=list(varwidth=TRUE)))
dev.off()

pdf('test.pdf')
  par(las=2,mar=c(10,5,.5,.5))
  with(gatc[!grepl('T4',gatc$data_set)&gatc$total_motif_occurrences==0,],vpPlot(data_set,totalProb*log10(exp(1)),ylab='Probability of no GATC log10(p)',offsetXArgs=list(varwidth=TRUE)))
dev.off()
