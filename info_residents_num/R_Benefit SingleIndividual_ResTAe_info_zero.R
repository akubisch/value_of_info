# Version for new Alex Data with modiefied mutation (mutating logit transform) and with explicit data on dispersal trait. 
# Estimating of emigration based on distribution statistics (mean, s.d. thus not necessary)

# extractin mean traits of resident TAe population at end of simulation (averaged from t=9100 to t=10000)
#setwd("/media/thomas/USBArbeit/Documents/Publications/09_TimingEmigration_Alex01/09_03ValueOfInformation/popsizes_I")
TAeInfo=1; TAeInfo # Infomrationslevels der residenten Population (Szenariem)
ndata=length(TAeInfo)

# we first need to extract the mean trait data

repli = 20
m_emi_res = numeric()

for (i in 1:ndata) {
  dummy = numeric()
  for (r in 1:repli) {
    filename = paste0("../info_invader_ibm/results/output_I1.0_rep",r,".txt")
    file = read.table(filename,header=T)
    dummy[r] = mean(file$emi_res)
  }
  m_emi_res[i] = mean(dummy)
}

####################################################
data_pre=matrix(nrow=1000000, ncol=ndata)
data_post=matrix(nrow=1000000, ncol=ndata)
lengthfile=numeric(ndata)

for (i in 1:ndata) # einlesen der unterschiedlichen (Szenarien) Populationsdaten in Matrix
{
  filename = paste0("../info_invader_ibm/results/popsizes_I",format(TAeInfo[i],nsmall=1),".txt")
  pops=read.table(filename, header=T) # file "popsize_distributions..." pre- and postdispersal pop. sizes
  head(pops)
  lengthfile[i]=length(pops$t)
  data_pre[,i]=pops$pre_disp[1:1000000]
  data_post[,i]=pops$post_disp[1:1000000]
  rm(pops)
}



############################################################################################################
## simulation parameters - check!!!
Cap=100
mu=0.1
Lambda=2
sigma=1


CTs=round(0.2*Cap); CTe=round(5*Cap) # range of possible threhold values
vCT=CTs:CTe # hier CT in 1-er Schritten 

# set up vectors and matrices to store results
M_resDIE_I0=numeric(ndata)   #mean emigration for resident DIE population (last generation) -> from data
M_bestVT_I0=numeric(ndata)   #find INDEX of threhold value that generates optimal fitness
M_bestgVT_I0=numeric(ndata)  #dito for geometric mean
M_bestemi_DD_I0=numeric(ndata) # ... and corresponding realized emigration prob.
M_fitMax_DD_I0=numeric(ndata) # max. arith fitness
M_gfitMax_DD_I0=numeric(ndata) # max. geom fitness
M_mfit_DI_I0=numeric(ndata) # fitness of non-informed invader
M_gfit_DI_I0=numeric(ndata) # dito geom. mean
M_geo0_DI_I0=numeric(ndata) # geom. fitness never emigrating invader
M_geo1_DI_I0=numeric(ndata) # geom. fitness always emigrating invader
M_vemip_DD_I0=matrix(nrow=length(vCT), ncol=ndata) # emigration prob for full-informed for given CTs
M_vmfit_DD_I0=matrix(nrow=length(vCT), ncol=ndata) # dito fitness
M_vgmfit_DD_I0=matrix(nrow=length(vCT), ncol=ndata)# dito geom. fitness

#select data (scenario) to be analyzed
for(i in 1:ndata)
{
  
  pre=data_pre[,i]
  post=data_post[,i]
  resDIE=m_emi_res[i] #mean emigration for resident DTe population (last generation)
  aa=(Lambda-1)/Cap # term a in Beverton-Holt
  exp_fertility=Lambda*1/(1+post*aa) # lambda*survival; expected payoff in patch j AFTER dispersal
  mfit_emi=mean(exp_fertility) # expected mean pay-off for SUCCESSFUL dispersers
  mfit_emimu=mfit_emi*(1-mu) # dito but and accounting for dispersal mortality
  # note that this is the expeted fitness for any dispersing individual, regardless of DD or DI strategy

  mfit_philo_DI=sum(exp_fertility*pre)/sum(pre)
  mfit_DI=resDIE*mfit_emimu+(1-resDIE)*mfit_philo_DI

  geof_philo_DI=(1-resDIE)*sum(log(exp_fertility)*pre) # sum logarithms of expected fertilities for philos
  geof_disp_DI=resDIE*sum(pre)*log(mfit_emimu) # geopm fitness of DI dispersers 
  gfit_DI=exp((geof_philo_DI+geof_disp_DI)/sum(pre)) # 
  vemip=c(0,1)
  vmfit_DI=vemip*mfit_emimu+(1-vemip)*mfit_philo_DI
  geo0_DI=exp(sum(log(exp_fertility)*pre)/sum(pre))
  geo1_DI=exp(log(mfit_emimu))
  
  vemip_DD=numeric(length(vCT))
  vmfit_philoDD=numeric(length(vCT))
  vmfit_DD=numeric(length(vCT))
  vgmfit_DD=numeric(length(vCT))

  for (k in 1:length(vCT))
  {
    Emigrants=sum(pre[pre>vCT[k]]) #absolut number of emigrantion events under Bang-Bang with threshold CT
    vemip_DD[k]=Emigrants/sum(pre) # virtual emigration probability of BanG-Bang if applied by all
    # arithmetic mean of fitness
    vmfit_philoDD[k]=sum(pre*(pre<=vCT[k])*exp_fertility)/sum(pre)
    vmfit_DD[k]=vemip_DD[k]*mfit_emimu+vmfit_philoDD[k]
    # geometric mean
    philos=sum(pre*(pre<=vCT[k])*log(exp_fertility))
    disps=sum(pre*(pre>vCT[k])*log(mfit_emimu))
    combined=(philos+disps)/sum(pre)
    vgmfit_DD[k]=exp(combined)

  }

  M_resDIE_I0[i]=resDIE  #mean emigration for resident DIE population (last generation)
  M_bestVT_I0[i]=which(vmfit_DD==max(vmfit_DD))[1] # find index for threhold value CT that generates optimal fitness
  M_bestgVT_I0[i]=which(vgmfit_DD==max(vgmfit_DD))[1]; # dito for geometric mean
  M_fitMax_DD_I0[i]=vmfit_DD[M_bestVT_I0[i]] # max fitness for arithmetic mean
  M_gfitMax_DD_I0[i]=vgmfit_DD[M_bestgVT_I0[i]] #max fitness  for geometric mean
  M_mfit_DI_I0[i]=mfit_DI
  M_gfit_DI_I0[i]=gfit_DI
  M_geo0_DI_I0[i]=geo0_DI
  M_geo1_DI_I0[i]=geo1_DI
  M_vemip_DD_I0[,i]=vemip_DD # threshold dependen emigration prob of CT
  M_bestemi_DD_I0[i]=vemip_DD[M_bestVT_I0[i]]
  M_vmfit_DD_I0[,i]=vmfit_DD
  M_vgmfit_DD_I0[,i]=vgmfit_DD
  
  rm(resDIE, vmfit_DD, vmfit_DI, vmfit_philoDD, vgmfit_DD, gfit_DI, geo1_DI, geo0_DI, geof_disp_DI, geof_philo_DI,
     vemip, vemip_DD, philos, disps, combined, Emigrants, exp_fertility, 
     mfit_emi, mfit_emimu, mfit_DI, mfit_philo_DI)
}

rm(i,ii, post, pre,  aa)

M_bestCT_I0=vCT[M_bestVT_I0]
M_bestgCT_I0=vCT[M_bestgVT_I0]

##################################### end data processing ##################################################

save.image(file="processed_data_I0.RDat")
