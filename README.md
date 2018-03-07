# mmgenome2

## Installation
As ADMINISTRATOR!
```r
#check for bioconductor installer
if(!require(BiocInstaller)) 
  source("https://bioconductor.org/biocLite.R")
  
#check for devtools
if(!require(devtools))
  install.packages("devtools")
  
#install mmgenome2
BiocInstaller::biocLite("kasperskytte/mmgenome2")
```