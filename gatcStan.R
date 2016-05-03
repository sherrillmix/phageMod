library(rstan)

if(!exists('gatc'))source('analyze_gatc.R')

gatc_dat<-list(
  'nObs'=sum(!grepl('T4',gatc$data_set)),
  'counts'=gatc[!grepl('T4',gatc$data_set),'total_motif_occurrences'],
  'isPhage'=gatc[!grepl('T4',gatc$data_set),'data_set']=='Phage_25'
  'expected'=gatc[!grepl('T4',gatc$data_set),'prob']*gatc[!grepl('T4',gatc$data_set),'contig_len']
)
fit<-stan(file='poisson.stan',data=gatc_dat,iter=1000,chains=4,cores=4)
pdf('test.pdf');traceplot(fit);dev.off()

