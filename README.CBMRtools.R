#CBMGIT <- Sys.getenv('CBMGIT')
CBMGIT<-getwd()
if (CBMGIT=="") stop( "Use 'setenv CBMGIT ..' to set CBMgithub's base directory" )

setwd(paste(CBMGIT,'/CBMRtools/inst/package_make',sep=''))
source("CBMRtools.build.R") 

