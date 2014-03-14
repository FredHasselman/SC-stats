sciCure
========
### A Toolbox for Curating Scientific Knowledge Extracted from Scientific Publications
#### Tool development for the Curate Science website

*Very BETA version - [Fred Hasselman](http://fredhasselman.com)*

All the functions that are required are available in `scicuRe_source.R`, a package will follow soon. 
Use the code below to source it directly from GitHub (the source function was found [here](http://tonybreyal.wordpress.com/2011/11/24/source_https-sourcing-an-r-script-from-github/)). It requires you to install and load the `RCurl` package

```
source_https <- function(url, ...) {
  # load the package 
  require(RCurl)

  # parse and evaluate each .R script
  sapply(c(url, ...), function(u) {
    eval(parse(text = getURL(u, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), envir = .GlobalEnv)
  })
}

# Source the scicuRe_source.R toolbox!
source_https("https://raw.github.com/FredHasselman/SC-stats/master/scicuRe_source.R")
```

**SEARCH STRATEGY**

 - Need stats labelled for each experiment in a multi-experiment study
 - Need measures of effect size
 - Need aggregate stats, i.e. based on post-hoc power and more
 - Need information to decide if proper analysis was used on stats ([`statcheck`](https://github.com/MicheleNuijten/statcheck) provides some interesting tests!)
    
*Therefore:*    
 - Search is implemented to be driven by the hierarchical structure of an empirical report
 - Structure is: Header (e.g., Experiment 1) -> Section (e.g., Results) -> [Table 1. ->] *STAT OF INTEREST* <- Page number [<- Appendix]
 - Extract as much information as possible!
    
The `statcheck` package bij M. Nuijten is available [here](https://github.com/MicheleNuijten/statcheck)   
Download, unzip and install, e.g. by running:    
```
install.packages("~/Downloads/statcheck-master/", repos = NULL, type="source")
```   

The authors of `statcheck` note:   
> The `pdftotext` program (http://www.foolabs.com/xpdf/download.html) is used to convert PDF files to plain text files. This must be installed and PATH variables must be properly set so that this program can be used from command line.
   
For the code in this package to execute correctly `pdftotxt` and `pdfinfo` need to be available on your system as well.   
Check it now:
```
if(all(file.exists(Sys.which("pdftotext")))) print("YES!") else print("NO!")
```   

