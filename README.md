
.. use and add to CBMRtools
**Install CBMRtools with install_github**
========

You need R 3.0.0 or higher. If on scc4, before starting R, type
module load R/R-3.1.1

Some dependencies you may need to install first:
install.packages("devtools")
require("devtools")
install_github("hadley/staticdocs")
install_github("andrie/ggdendro")


From within R, run the following commands:
library(devtools)
PAT <- "04fe676593e46b6bda5a5d09431156e8a500349a"

#old tag
install_github("montilab/CBMgithub/scripts/R/CBMRtools",
ref="CBMRtools_v1.2.0", # desired tag
auth_token = PAT)       # user token
	#newest version
	#to be added

require(CBMRtools)



You need to do this only once (or after the package is modified by yourself or others). Thereafter, you will only need to include the ‘require(CBMRtools)’ command in your R code. 

See all available functions and datasets:
require(“CBMRtools”)
ls("package:CBMRtools")

See usage for particular functions listed above, e.g.:
?nnAnalysis

for detailed help pages and outputs of example runs: 
open <path to CBMRtools>/inst/web/index.html

Note on Installation on PC: on installation on a PC. You may need to install the devtools and Rtools, and you may encounter some problem installing the required packages. You may need to install these packages manually and then reinstall the CBMRtools, then it will work for your PC. 

**Contributing to CBMRtools**
========

Set up git repo for CBMRtools:

git clone https://lia978@github.com/montilab/CBMRtools.git
replace lia978 with your username, will prompt for your git password

or, if you already have it, pull the latest version.

CBMRtools is located at /CBMRtools/R

Add your R function to <path to CBMRtools>/R/ **do not add standalone scripts, wrap all R code into functions**

Specify dependencies: within the file CBMRtools/DESCRIPTION, add to the section Depends the list of packages your function depends on (if any)

Add license info to the file and roxygen comments at the top of each function you want to document. You don’t need to write this from scratch, see template at
	CBMgithub/scripts/R/template.R 
and an example of how to use this template at 
	<path to CBMRtools>/R/variationFilter.R. 

You can then build the package as follows:
	setwd(“<path to CBMRtools>/inst/package_make/”)
	source(“CBMRtools.build.R”) 

this script does the package checking, compiling, and installation, also generates html page pages. The help pages are be found in CBMRtools/inst/web, see index.html and make sure your page is properly added. 

If the above steps work, sync your changes with the git repo: instructions.
Notice that these steps need to be performed only once (or after you modify functions within the package). Thereafter, you only need to include ‘require(CBMRtools)’ in your R code.

Notes: 
This will only work if you run it within the package_make folder as the compilation takes place in the CBMRtools subfolder relative to this folder.  
Make sure there is a folder under package_make called staticdocs within the inst folder or the build will return an error while trying to write the html documentation pages.
All R objects that are used to illustrate example usage must be added to the CBMRtools/R/data.R script. Avoid using rds objects that need pointing to absolute paths to load/read specific files. This way, your data object can be directly recognized and loaded using the data() function.
Any commands within @examples that calls on helper functions need to have the helper functions be exported with @export, e.g. see heatmap.ggplot.R‘s @examples, which uses hcopt.R and to.hex.R

