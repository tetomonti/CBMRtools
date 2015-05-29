CBMGIT <- Sys.getenv('CBMGIT')
if (CBMGIT=="") stop( "Use 'setenv CBMGIT ..' to set CBMgithub's base directory" )

setwd(paste(CBMGIT,'/CBMRtools/CBMRtools/inst/package_make',sep=''))
source("CBMRtools.build.R") 

