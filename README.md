[![Travis Build Status](https://travis-ci.org/KasperSkytte/mmgenome2.svg?branch=master)](https://travis-ci.org/KasperSkytte/mmgenome2) 

# Tools for extracting individual genomes from metagenomes
`mmgenome2` is an R-package designed to facilitate reproducible extraction of individual genomes from metagenomes. It is an implementation of the different binning strategies described in the [multi-metagenome](http://madsalbertsen.github.io/multi-metagenome/) project, and makes it possible to apply these to any metagenome data. In combination with the [RMarkdown](https://rmarkdown.rstudio.com/) format, the mmgenome2 package allows for reproducible step-by-step extraction of high-quality genomes from metagenomes as well as generating publication-ready figures with minimal effort. 

## Installation
First, install [R (3.4.3 or later)](https://mirrors.dotsrc.org/cran/) and [RStudio](https://www.rstudio.com/products/rstudio/download/#download). Windows users should also install [RTools](https://mirrors.dotsrc.org/cran/bin/windows/Rtools/). Then open RStudio as administrator (!) and run the following commands (just copy/paste) to install `mmgenome2` from the console:

```r
#check for remotes
if(!require(remotes))
  install.packages("remotes")
  
#install mmgenome2 using remotes
remotes::install_github("kasperskytte/mmgenome2")
```

### Installation on MAC
To install `mmgenome2` on MAC please see [this](articles\MACinstall.html) before running the above commands.

## Get started
For a brief guide about the basics of mmgenome2 go to the [Get Started](https://kasperskytte.github.io/mmgenome2/articles/mmgenome2.html) page. 