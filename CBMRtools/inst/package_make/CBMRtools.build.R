## Sequence of commands to check and build the package 
## and generate html page documents

#install.packages(devtools)
#require(devtools)
#install_github('hadley/pkgdown') 

require(devtools)
require(pkgdown)

#set path to package home directory
package.dir <- normalizePath("../../../CBMRtools")

cat("Documenting...\n")
document(package.dir) # creates help pages

#cat("Checking...\n")
#check(package.dir) # checking

#cat("Loading...\n")
#load_all(package.dir) # loading

cat("Installing..\n")
install(package.dir, dependencies = TRUE) #installing

library(CBMRtools)

## the directory 'staticdocs' must exist under CBMRtools/inst/ for the
## ..help pages to be installed at CBMRtools/inst/web

STATICDIR <- '../staticdocs'
if ( is.na(file.info(STATICDIR)$isdir) ) system(paste('mkdir',STATICDIR))

cat("Making html pages...\n")
## generate html pages
## update ./inst/staticdocs/index.r for formatting of help page sections
setwd(package.dir)
build_site(pkg = package.dir, examples = TRUE)

#cat("Installing locally...\n")l
## install CBMRtools locally
#install.packages(package.dir, dependencies = TRUE, repos = NULL, type = "source")


##build tar.gz 
system("R CMD build ../CBMRtools")
