#CBMGIT <- Sys.getenv('CBMGIT')
CBMGIT<-getwd()
if (CBMGIT=="") stop( "Use 'setenv CBMGIT ..' to set CBMgithub's base directory" )

<<<<<<< HEAD
setwd(paste(CBMGIT,'/CBMRtools/CBMRtools/inst/package_make',sep=''))
=======
setwd(paste(CBMGIT,'/CBMRtools/inst/package_make',sep=''))
>>>>>>> 1609597a9cf9c5dd080931f083ffb9b3a6370a98
source("CBMRtools.build.R") 

