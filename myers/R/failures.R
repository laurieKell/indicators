library(sraplus)
library(plyr)

## Set of stocks from RAM legacy database that failed
load("results/failures.RData")

## Where to save results
## Need subdir called db, db4 and SRA for the different procedures
dirMy=getwd()

## Wrapper functions
### Biomass dynamic model
fnBD<-function(assessid,stockid,biomass,catch,year,spp,dir,cv=NULL){
  
  stockid =stockid[1]
  assessid=assessid[1]
  
  print(paste(spp[1],assessid))
  
  cat(paste(assessid,"\n"),file=file.path(dir,stockid),append=TRUE)
  
  if (is.numeric(cv)){
    set.seed(1233)
    biomass=biomass*exp(rlnorm(length(biomass),0,cv))}
  biomass=biomass*mean(catch)*4
  
  cdp<-try(sraplus::format_driors(
    b_ref_type ="b",
    index      = biomass,
    index_year = year,
    catch      = catch,
    years      = year,
    taxa       = spp[1],
    use_heuristics=!TRUE))
  
  if ("try-error"%in%is(cdp)) return(NULL)
  
  fit<-try(sraplus::fit_sraplus(
    driors = cdp,      
    engine = "tmb",
    model  = "sraplus_tmb",
    adapt_delta = 0.9,
    max_treedepth = 10,
    n_keep = 4000,
    chains = 1,
    cores = 1,
    estimate_qslope = FALSE,
    estimate_proc_error = TRUE))
  
  if (!("try-error"%in%is(fit)))
    save(fit,file=file.path(dir,paste(stockid,assessid,"RData",sep=".")))
  
  save(biomass,file=file.path(dir,paste(stockid,assessid,"u",sep=".")))
  
  "try-error"%in%is(fit)}

##SRA
fnSRA<-function(assessid,stockid,biomass,catch,year,spp,dir,cv=NULL){
  
  stockid =stockid[1]
  assessid=assessid[1]
  
  print(paste(spp[1],assessid))
  
  cat(paste(assessid,"\n"),file=file.path(dir,stockid),append=TRUE)
  
  if (is.numeric(cv)){
    set.seed(1233)
    biomass=biomass*exp(rlnorm(length(biomass),0,cv))}
  
  cdp<-try(sraplus::format_driors(
    taxa           =spp[1],
    use_heuristics =TRUE,
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
    estimate_proc_error = TRUE))
  
  if (!("try-error"%in%is(fit)))
    save(fit,file=file.path(dir,paste(stockid,assessid,"RData",sep=".")))
  
  "try-error"%in%is(fit)}

control=failures[!duplicated(failures[,c("assessid","stockid")]),c("assessid","stockid")]

##Benchmark with perfect biomass index
for (i in seq(dim(control)[1]))
  with(subset(failures,assessid==control[i,"assessid"]&stockid==control[i,"stockid"]),  
       fnBD(assessid,stockid,biomass,catch,year,species,dir=file.path(dirMy,"results/bd")))

##CV of 40%
for (i in seq(dim(control)[1]))
  with(subset(failures,assessid==control[i,"assessid"]&stockid==control[i,"stockid"]),  
       fnBD(assessid,stockid,biomass,catch,year,species,dir=file.path(dirMy,"results/bd4"),cv=0.4))

##SRA
for (i in seq(control)[1])
  with(subset(failures,assessid==control[i,"assessid"]&stockid==control[i,"stockid"]),  
       fnSRA(assessid,stockid,biomass,catch,year,species,dir=file.path(dirMy,"results/sra")))
