**Installing CBMRtools with Github repo**
========

You need R 3.0.0 or higher. If on scc4, before starting R, type:

    module load R/R-3.1.1

it is also recommended to ssh with X11 forward (ssh -X)

Some dependencies you may need to install first are specified in the [DESCRIPTION](https://github.com/montilab/CBMRtools/blob/master/CBMRtools/DESCRIPTION) file. While  install_github() tries to install package dependencies from CRAN, if you non-CRAN packages such as Bioconductor packages, you may need to install them independently. 
See this [package](https://github.com/lia978/RPackageDependenciesInstall) for a fast and easy way to install both Cran and bioconductor packages listed from the DESCRIPTION file.


Next, to install CBMRtools, run the following commands within R:

```R
    library(devtools)
	
    #install from master
    install_github("montilab/CBMRtools/CBMRtools")
    #install from specific branch
    install_github("montilab/CBMRtools/CBMRtools",ref="v1.1.3")
    require(CBMRtools)
```

# Documentations
documentation for this package is currently found at:
http://montilab.bumc.bu.edu/~montilab/CBMRtoolsHtml/web/index.html

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

Within this github, CBMRtools' code is located at /CBMRtools/R

Add your R function to <code>/CBMRtools/R/</code> **do not add standalone scripts, wrap all R code into functions**

Specify dependencies: within the file <code>/CBMRtools/DESCRIPTION</code>, add to the section Depends the list of packages your function depends on (if any)

Add license info to the file and roxygen comments at the top of each function you want to document. You don’t need to write this from scratch, see template at	<code>template.R</code> and an example of how to use this template at 
<code>/CBMRtools/R/variationFilter.R</code>.

You can then build the package as follows:

    setwd(“<path to CBMRtools>/inst/package_make/”)
    source(“CBMRtools.build.R”) 

this script does the package checking, compiling, and installation, also generates html page pages. The help pages are be found in <code>CBMRtools/inst/web</code>, see index.html and make sure your page is properly added. 

If the above steps work, sync your changes with the git repo: instructions.
Notice that these steps need to be performed only once (or after you modify functions within the package). Thereafter, you only need to include <code>require(CBMRtools)</code> in your R code.

Notes: 
- This will only work if you run it within the package_make folder as the compilation takes place in the CBMRtools subfolder relative to this folder.  
- Make sure there is a folder under <code>package_make</code> called <code>staticdocs</code> within the inst folder or the build will return an error while trying to write the html documentation pages.
- All R objects that are used to illustrate example usage must be added to the <code>CBMRtools/R/data.R</code> script. Avoid using rds objects that need pointing to absolute paths to load/read specific files. This way, your data object can be directly recognized and loaded using the data() function.
- Any commands within @examples that calls on helper functions need to have the helper functions be exported with @export, e.g. see <code>heatmap.ggplot.R</code>‘s @examples, which uses <code>hcopt.R</code> and <code>to.hex.R</code>.

