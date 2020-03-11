library(reshape)
library(plyr)
library(dplyr)
library(ggplot2)

dirMy="/home/laurence-kell/Desktop/papers/data-poor"



###########################################################################################################
## Datasets ###############################################################################################
###########################################################################################################

## FAO Taxa
load(file.path(dirMy,"data/fao_taxa.rda"))

## ICES
load(file.path(dirMy,"data/sag_2019.RData"))

## Myers database
load(file.path(dirMy,"data/myersDB.RData"))
myers=subset(timeseries,tsid%in%c("BdivBmsytouse-dimensionless","SSB-MT","TCbest-MT"))
myers=cast(myers,assessid+stockid+tsyear~tsid,value="tsvalue")
names(myers)[3:6]=c("year","bbmsy","biomass","catch")
myers=subset(myers,!(is.na(biomass)|is.na(catch)|is.na(bbmsy)))
myers=merge(myers,stock[,c("stockid","scientificname","region")],by="stockid")
names(myers)[7]="species"
length(unique(myers$assessid))

myers=ddply(myers,.(assessid), with, data.frame(
                  stockid=as.character(stockid),
                  sp     =(c(biomass[-1],NA)-biomass+catch)/mean(biomass/bbmsy,na.rm=T),
                  bbmsy  =bbmsy,
                  catch  =catch,
                  biomass=biomass,
                  year   =year,
                  stringsAsFactors=FALSE))

myers=merge(myers,stock,by="stockid")
names(myers)[9]="species"

## Fishbase
load(file.path(dirMy,"data/fb.RData"))
fb=ddply(fb, .(species),  function(x)  data.frame(
  k   =mean(x$k,   na.rm=T),
  linf=mean(x$linf,na.rm=T),
  t0  =mean(x$t0,  na.rm=T),
  l50 =mean(x$lmat,   na.rm=T),
  a   =mean(x$a,      na.rm=T),
  b   =mean(x$b,      na.rm=T),
  m   =mean(x$m,      na.rm=T)))

###############################################################################################################
## Functions ##################################################################################################  
###############################################################################################################

fnBD<-function(assessid,stockid,biomass,catch,year,species,dir,cv=NULL,p=1.01){
  
  stockid =stockid[1]
  assessid=assessid[1]
  species =species[1]
  
  print(paste(species,assessid))
  
  cat(paste(assessid,"\n"),file=file.path(dir,stockid),append=TRUE)
  
  index=biomass
  if (is.numeric(cv)){
    set.seed(1233)
    index=index*exp(rlnorm(length(biomass),0,cv))}
  
  cdp<-try(sraplus::format_driors(
    index      = index,
    index_year = year,
    catch      = catch,
    years      = year,
    
    shape_prior   =p,
    
    taxa       = species,
    use_heuristics=!TRUE))
  
  if ("try-error"%in%is(cdp)) return(NULL)
  
  fit<-try(sraplus::fit_sraplus(
    driors = cdp,      
    engine = "tmb",
    model  = "sraplus_tmb",
    estimate_proc_error = TRUE,
    estimate_qslope = FALSE))
  
  if (!("try-error"%in%is(fit)))
    save(fit,file=file.path(dir,paste(stockid,assessid,"RData",sep=".")))
  
  save(index,file=file.path(dir,paste(stockid,assessid,"u",sep=".")))
  
  "try-error"%in%is(fit)}

##SRA
fnSRA<-function(assessid,stockid,biomass,catch,year,spp,dir,cv=NULL,p=1.01){
  
  stockid =stockid[1]
  assessid=assessid[1]
  
  print(paste(spp[1],assessid))
  
  cat(paste(assessid,"\n"),file=file.path(dir,stockid),append=TRUE)
  
  cdp<-try(sraplus::format_driors(
    taxa           =spp[1],
    use_heuristics =TRUE,
    
    shape_prior    =p,
    #shape_prior_cv =0.1,
    
    catch          =catch,
    years          =year))
  
  if ("try-error"%in%is(cdp)) return(NULL)
  
  fit<-try(sraplus::fit_sraplus(
    driors = cdp,
    engine = "sir",
    adapt_delta = 0.9,
    max_treedepth = 10,
    n_keep = 4000,
    chains = 1,
    cores = 1,
    estimate_qslope = FALSE,
    estimate_shape  = FALSE,
    estimate_proc_error = TRUE))
  
  if (!("try-error"%in%is(fit)))
    save(fit,file=file.path(dir,paste(stockid,assessid,"RData",sep=".")))
  
  "try-error"%in%is(fit)}

fnDan<-function(assessid,stockid,biomass,catch,year,spp,dir,cv=NULL){
  
  stockid =stockid[1]
  assessid=assessid[1]
  
  print(paste(spp[1],assessid))
  
  cat(paste(assessid,"\n"),file=file.path(dir,stockid),append=TRUE)
  
  if (is.numeric(cv)){
    set.seed(1233)
    biomass=biomass*exp(rlnorm(length(biomass),0,cv))}
  biomass=biomass*mean(catch)*4
  
  cdp<-format_driors(taxa = spp[1],
                     b_ref_type ="b",
                     catch = catch,
                     years = year,
                     initial_state = 1,
                     initial_state_cv = .1,
                     terminal_state = NA,
                     shape_prior = 1.01,
                     shape_prior_cv = 1,
                     growth_rate_prior = 0.4,
                     growth_rate_prior_cv = 0.5,
                     sigma_ratio_prior = 1,
                     sigma_ratio_prior_cv = .1,
                     sigma_obs_prior = 0.1,
                     sigma_obs_prior_cv = 1)
  
  if ("try-error"%in%is(cdp)) return(NULL)
  
  fit <- fit_sraplus(driors = cdp,
                     engine = "tmb",
                     model = "sraplus_tmb",
                     estimate_shape = FALSE,
                     estimate_proc_error = FALSE,
                     estimate_k = TRUE,
                     learn_rate = 2e-1,
                     n_keep = 2000,
                     eps = 1e-12,
                     adapt_delta = 0.95,
                     marginalize_q = FALSE,
                     max_treedepth = 12)
  
  if (!("try-error"%in%is(fit)))
    save(fit,file=file.path(dir,paste(stockid,assessid,"RData",sep=".")))
  
  save(biomass,file=file.path(dir,paste(stockid,assessid,"u",sep=".")))
  
  "try-error"%in%is(fit)}

###########################################################################################################
## Control  ###############################################################################################
###########################################################################################################

control=myers[!duplicated(myers[,c("assessid","stockid")]),c("assessid","stockid")]
control=subset(control,!(assessid%in%c("CSIRO-TIGERFLATSE-1915-2015-MOESENEDER",
                                       "WGCSE-MEGSPPIVa-VIa-1985-2016-ICESIMP2018",
                                       "NEFSC-ATHAL5YZ-1800-2007-COL",
                                       "SEFSC-GTRIGGM-1945-2013-SISIMP2016",
                                       "NIWA-OROUGHYNZMEC-1882-2014-FU",
                                       "NIWA-OROUGHYNZMEC-1909-2011-CORDUE",
                                       "FAFRFJ-PILCHTSST-1960-2013-JPNIMP2016",
                                       "SEFSC-PINKSHRIMPGM-1984-2013-SISIMP2016",
                                       "NAFO-SC-REDFISHSPP3LN-1959-2013-WATSON",
                                       "IFOPCH-RSQLOBSTERNCH-1969-2015-PARMA",
                                       "NIWA-SCMPMB-1989-2014-FU",
                                       "CSIRO-TIGERFLATSE-1913-2006-FULTON",
                                       "SEFSC-WSHRIMPGM-1984-2014-SISIMP2016",
                                       "SEFSC-YEGROUPGM-1975-2009-HIVELY",
                                       "NEFSC-SFLOUNMATLC-1982-2014-SISIMP2016",
                                       "NIWA-NZSNAPNZ1BOP-HAGU-1899-2013-FU",
                                       "NIWA-NZSNAPNZ1ENLD-1899-2013-FU",
                                       "AFSC-TANNERCRABBSAI-1965-2015-SISIMP2016",
                                       "NEFSC-SFLOUNMATLC-1940-2012-HIVELY",
                                       "SPC-SKJCWPAC-1950-2015-PONS",
                                       "SPC-SKJCWPAC-1950-2012-PONS",
                                       "NWFSC-BLACKROCKSPCOAST-1915-2007-BRANCH",
                                       "NWFSC-CHROCKCPCOAST-1900-2015-SISIMP2016",
                                       "SWFSC-COWCODSCAL-1900-2012-SISIMP2016 COWCODSCAL",
                                       "SWFSC-CPRROCKPCOAST-1916-2012-SISIMP2016",
                                       "NWFSC-CROCKPCOAST-1916-2011-STACHURA",
                                       "NWFSC-ESOLEPCOAST-1876-2007-BRANCH",
                                       "SWFSC-GRSPROCKSCAL-1916-2010-STACHURA",
                                       "NWFSC-LSTHORNHPCOAST-1964-2012-HIVELY",
                                       "CSIRO-MORWONGSE-1913-2007-FULTON",
                                       "RSNAPSATLC-1954-2010-HIVELY")))

###########################################################################################################
## p=1.01   ###############################################################################################
###########################################################################################################

##Benchmark
for (i in rev(seq(dim(control)[1])))
     with(subset(myers,assessid==control[i,"assessid"]&stockid==control[i,"stockid"]),  
          fnBD(assessid,stockid,biomass,catch,year,species,dir=file.path(dirMy,"results/base/bd")))

##CV of 40%
for (i in seq(dim(control)[1]))
  with(subset(myers,assessid==control[i,"assessid"]&stockid==control[i,"stockid"]),  
       fnBD(assessid,stockid,biomass,catch,year,species,dir=file.path(dirMy,"results/base/bd4"),cv=0.4))

seq(dim(control)[1])[control[,1]=="SEFSC-RSNAPGM-1872-2013-SISIMP2016"]+1

##SRA
for (i in 303:dim(control)[1])
  with(subset(myers,assessid==control[i,"assessid"]&stockid==control[i,"stockid"]),  
       fnSRA(assessid,stockid,biomass,catch,year,species,dir=file.path(dirMy,"results/base/sra")))

runs=as.data.frame(rbind(
           cbind(run="bd",   assessid=list.files(file.path(dirMy,"results/base/bd" ))),
           cbind(run="bd4",  assessid=list.files(file.path(dirMy,"results/base/bd4"))),
           cbind(run="sra",  assessid=list.files(file.path(dirMy,"results/base/sra")))))
runs=runs[grep(".RData",runs$assessid),]

runs=mdply(runs, function(run,assessid) {
       print(paste(run,assessid))  
       load(file.path(file.path(dirMy,"results/base"),run,assessid))
       
       if ("try-error"%in%is(fit)) return (NULL)
       
       fit$results})
runs$assessid=gsub(".RData","",runs$assessid)
runs$assessid=laply(strsplit(runs[,2],"\\."),function(x) x[[2]])

save(runs,file=file.path(dirMy,"results/myersp2Runs.RData"))

##Dan
for (i in seq(dim(control)[1]))
  with(subset(myers,assessid==control[i,"assessid"]&stockid==control[i,"stockid"]),  
       fnDan(assessid,stockid,biomass,catch,year,species,dir=file.path(dirMy,"results/dan"),0.4))

###########################################################################################################
## p=2 ####################################################################################################
###########################################################################################################
##Benchmark
for (i in 1:dim(control)[1])
  with(subset(myers,assessid==control[i,"assessid"]&stockid==control[i,"stockid"]),  
       fnBD(assessid,stockid,biomass,catch,year,species,dir=file.path(dirMy,"results/p2/bd"),p=2))

##CV of 40%
c("PAUAPAU7")
for (i in 185:dim(control)[1])
  with(subset(myers,assessid==control[i,"assessid"]&stockid==control[i,"stockid"]),  
       fnBD(assessid,stockid,biomass,catch,year,species,dir=file.path(dirMy,"results/p2/bd4"),cv=0.4,p=2))

##SRA
c("NEFSC-ACADREDGOMGB-1913-2007-MILLER",
  "NEFSC-ACADREDGOMGB-1913-2014-SISIMP2016",
  "NWFSC-ARFLOUNDPCOAST-1916-2007-BRANCH",
  "NIWA-AUSSALMONNZ-1975-2006-JENSEN",
  "NWFSC-BLACKROCKORECOAST-1892-2015-SISIMP2016",
  "FBDAFF-CHSQUIDSA-1971-2013-GLAZER",
  "DFO-CROCKWCVANISOGQCI-1945-2009-STANTON",
  "NWFSC-PSOLEPCOAST-1876-2015-SISIMP2016",
  "NWFSC-SABLEFPCOAST-1900-2015-SISIMP2016",
  "NWFSC-SSTHORNHPCOAST-1901-2005-STANTON")
for (i in (seq(dim(control)[1])[control[,1]=="NWFSC-SSTHORNHPCOAST-1901-2005-STANTON"]+1):dim(control)[1])
  with(subset(myers,assessid==control[i,"assessid"]&stockid==control[i,"stockid"]),  
       fnSRA(assessid,stockid,biomass,catch,year,species,dir=file.path(dirMy,"results/p2/sra"),p=2))

runs=as.data.frame(rbind(
  cbind(run="bd",   assessid=list.files(file.path(dirMy,"results/p2/bd" ))),
  cbind(run="bd4",  assessid=list.files(file.path(dirMy,"results/p2/bd4"))),
  cbind(run="sra",  assessid=list.files(file.path(dirMy,"results/p2/sra")))))
runs=runs[grep(".RData",runs$assessid),]

runs=mdply(runs, function(run,assessid) {
  print(paste(run,assessid))  
  load(file.path(file.path(dirMy,"results/p2"),run,assessid))
  
  if ("try-error"%in%is(fit)) return (NULL)
  
  fit$results})
runs$assessid=gsub(".RData","",runs$assessid)
runs$assessid=laply(strsplit(runs[,2],"\\."),function(x) x[[2]])

save(runs,file=file.path(dirMy,"results/myersp2Runs.RData"))
  
###########################################################################################################
## Failures ###############################################################################################
###########################################################################################################

failures=subset(myers,(assessid%in%c("CSIRO-TIGERFLATSE-1915-2015-MOESENEDER",
                                     "WGCSE-MEGSPPIVa-VIa-1985-2016-ICESIMP2018",
                                     "NEFSC-ATHAL5YZ-1800-2007-COL",
                                     "SEFSC-GTRIGGM-1945-2013-SISIMP2016",
                                     "NIWA-OROUGHYNZMEC-1882-2014-FU",
                                     "NIWA-OROUGHYNZMEC-1909-2011-CORDUE",
                                     "FAFRFJ-PILCHTSST-1960-2013-JPNIMP2016",
                                     "SEFSC-PINKSHRIMPGM-1984-2013-SISIMP2016",
                                     "NAFO-SC-REDFISHSPP3LN-1959-2013-WATSON",
                                     "IFOPCH-RSQLOBSTERNCH-1969-2015-PARMA",
                                     "NIWA-SCMPMB-1989-2014-FU",
                                     "CSIRO-TIGERFLATSE-1913-2006-FULTON",
                                     "SEFSC-WSHRIMPGM-1984-2014-SISIMP2016",
                                     "SEFSC-YEGROUPGM-1975-2009-HIVELY",
                                     "NEFSC-SFLOUNMATLC-1982-2014-SISIMP2016",
                                     "NIWA-NZSNAPNZ1BOP-HAGU-1899-2013-FU",
                                     "NIWA-NZSNAPNZ1ENLD-1899-2013-FU",
                                     "AFSC-TANNERCRABBSAI-1965-2015-SISIMP2016",
                                     "NEFSC-SFLOUNMATLC-1940-2012-HIVELY",
                                     "SPC-SKJCWPAC-1950-2015-PONS",
                                     "SPC-SKJCWPAC-1950-2012-PONS",
                                     "NWFSC-BLACKROCKSPCOAST-1915-2007-BRANCH",
                                     "NWFSC-CHROCKCPCOAST-1900-2015-SISIMP2016",
                                     "SWFSC-COWCODSCAL-1900-2012-SISIMP2016 COWCODSCAL",
                                     "SWFSC-CPRROCKPCOAST-1916-2012-SISIMP2016",
                                     "NWFSC-CROCKPCOAST-1916-2011-STACHURA",
                                     "NWFSC-ESOLEPCOAST-1876-2007-BRANCH",
                                     "SWFSC-GRSPROCKSCAL-1916-2010-STACHURA",
                                     "NWFSC-LSTHORNHPCOAST-1964-2012-HIVELY",
                                     "CSIRO-MORWONGSE-1913-2007-FULTON",
                                     "RSNAPSATLC-1954-2010-HIVELY")))

save(failures,file="/home/laurence-kell/Desktop/papers/data-poor/results/failures.RData")
