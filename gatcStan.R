library(rstan)

if(!exists('gatc'))source('analyze_gatc.R')

selector<-!grepl('T4',gatc$data_set)&!gatc$gappy

gatc_dat<-list(
  'nObs'=sum(selector),
  'counts'=gatc[selector,'nGATC'],
  'isPhage'=gatc[selector,'data_set']=='Phage_25',
  'expected'=gatc[selector,'prob']*gatc[selector,'length']
  #'isModified'=ifelse(is.na(gatc[selector,'unmodified_norm_mean']),FALSE,p.adjust(gatc[selector,'unmodified_norm_mean'],'fdr')<.05)
)

fit<-stan(file='poisson.stan',data=gatc_dat,iter=2000,chains=8,cores=8)
#pairs(fit,pars=c('betaIntercept','betaPhage','betaMod','betaPhageMod','theta'));
pdf('out/bayes.pdf',width=10);traceplot(fit,c('betaIntercept','betaPhage','theta'));dev.off()

