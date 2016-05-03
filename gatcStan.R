library(rstan)

library(ggmcmc)

if(!exists('gatc'))source('analyze_gatc.R')

gatc_dat<-list(
  'nObs'=sum(!grepl('T4',gatc$data_set)),
  'counts'=gatc[!grepl('T4',gatc$data_set),'total_motif_occurrences'],
  'isPhage'=gatc[!grepl('T4',gatc$data_set),'data_set']=='Phage_25',
  'expected'=gatc[!grepl('T4',gatc$data_set),'prob']*gatc[!grepl('T4',gatc$data_set),'contig_len']
)
fit<-stan(file='poisson.stan',data=gatc_dat,iter=10000,chains=8,cores=8)
pdf('test.pdf');traceplot(fit,c('beta','theta'));dev.off()

