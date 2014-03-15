# [sciCure](http://fredhasselman.github.io/scicuRe/)
# ========
# ### A Toolbox for Curating Scientific Knowledge by Extracting it from Scientific Publications
# #### Tool development for the Curate Science website
# 
# *Very BETA version - [Fred Hasselman](http://fredhasselman.com)*
# 
# All the functions that are required are available in `scicuRe_source.R`, a package will follow soon. 
# Use this code to source it directly from GitHub:

source_https <- function(url, ...) {
  require(RCurl)
  # parse and evaluate each .R script
  sapply(c(url, ...), function(u) {
    eval(parse(text = getURL(u, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), envir = .GlobalEnv)
  })
}
# Source the scicuRe_source.R toolbox!
source_https("https://raw.github.com/FredHasselman/scicuRe/master/scicuRe_source.R")

# The `source_https()` function was found [here](http://tonybreyal.wordpress.com/2011/11/24/source_https-sourcing-an-r-script-from-github/)


# DEMONSTRATION -----------------------------------------------------------
#
# The code below will:
# 1. Convert all PDFfiles in a folder to text files.
# 2. Look for common header and section indicators in those files
# 3. Search for specific references to statistical tests and sample info by combining pre-defined search patterns
# 4. Locate the statistics relative to the information found in 2
# 5. Run tests, calculate ESCI and more


# PREPARE -----------------------------------------------------------------
require(plyr)
require(reshape)
setwd("/Volumes/Fred HD/Rplus/gitR/scicuRe")
PDFpath <- paste(getwd(),"PDFfolder/",sep="/")

# 1. CONVERT PDF ----------------------------------------------------------
fnames          <- list.files(PDFpath,pattern=".pdf",full.names=TRUE)
txtfiles        <- sapply(fnames,getPDFs)
names(txtfiles) <- gsub(".pdf","",basename(fnames))


# 2. BUILD SKELETON -------------------------------------------------------
skeleton        <- structureIT(fnames,txtfiles)


# 3. META PATTERNS --------------------------------------------------------
#
# Most results have a structure:
# [name of test statistic] [df1..n] [value of statistic] [p-value] [effect size] 
#
# In the neigbourhood of the reported test we may encounter:
# [Mean/Median] [SD/SE] [N] [CI] 
#
# Also:
# [one/two tailed/...] [bonferroni/FWE/...] [levene/equal variances assumed/...]
#
# These are all combinations of patterns defined in loadPatterns()
# Let's call them meta-patterns:
meta      <- vector("list",length=6)
metafound <- vector("list",length=6)

# Strategy:
# 3a. Create meta patterns
# 3b. Collect the metapatterns for all files
# 4a. Place meta txt where they "belong" by preferring stats with expected locations in the article:
#    Major Header -> Section header -> [Table 1. ->] *STAT OF INTEREST* <- Page number [<- Appendix]
# 4b. Connect sample info to stats Major Header -> Section header = *PARTICIPANTS* -> [Table 1. ->] *MEAN, SD* -> *STAT OF INTEREST* <- Page number [<- Appendix]


# 3a. CREATE META ---------------------------------------------------------
# Load predefined regex search patterns
structfind <- loadPatterns()

# Here are 6 patterns commonly encountered in the behavioural sciences:

meta[[1]] <- paste("(",structfind[["statsPatterns"]][["F"]],"[[:blank:]\\,\\;]*)(",structfind[["evidencePatterns"]][["p"]],"[[:blank:]\\,\\;]*)*([[:blank:]\\,\\;]*",structfind[["effectPatterns"]][["eta2"]],")*",sep="")
# Meta 1 will return:
# F(int,int) = real
# F(int,int) = real, p = real
# F(int,int) = real, p = real, (partial)eta(2) = real
# Fs = real
#tidyIT(matchIT(meta[[1]],txtfiles))

meta[[2]] <- paste("(",structfind[["statsPatterns"]][["t"]],"[[:blank:]\\,\\;]*)(",structfind[["evidencePatterns"]][["p"]],")*([[:blank:]\\,\\;])*(",structfind[["effectPatterns"]][["cd"]],")*",sep="")
# Meta 2 will return:
# t(int) = ±real
# t(int) = ±real, p = real
# t(int) = ±real, p = real, d = real
# ts = ±real
#tidyIT(matchIT(meta[[2]],txtfiles))

meta[[3]] <- paste("(",structfind[["statsPatterns"]][["r"]],"[[:blank:]\\,\\;]*)(",structfind[["statsPatterns"]][["t"]],")*([[:blank:]\\,\\;])*(",structfind[["evidencePatterns"]][["p"]],")*", sep="")
# Meta 3 will return:
# r = ±real
# r(int) = ±real
# r = ±real, t(int) = ±real, p = real
# rs = ±real
#tidyIT(matchIT(meta[[3]],txtfiles))

meta[[4]] <- paste("(",structfind[["statsPatterns"]][["chi2"]],"[[:blank:]\\,\\;]*)(",structfind[["evidencePatterns"]][["p"]],")*([[:blank:]\\,\\;])*(",structfind[["effectPatterns"]][["omega2"]],")*","([[:blank:]\\,\\;])*(",structfind[["effectsPatterns"]][["phi"]],")*",sep="")
# Meta 4 will return:
# χ² = real
# χ²(int) = real
# χ² = real, p = real
# χ² = real,
#tidyIT(matchIT(meta[[4]],txtfiles))

meta[[5]] <- wordXsYo(structfind[["descriptivePatterns"]][["M"]],structfind[["descriptivePatterns"]][["sde"]])
# Meta 5 will return:
# means, medians and SD or SE
#tidyIT(matchIT(meta[[5]],txtfiles))

# Meta 6 will return:
# References to sample size and subjects such as "33 female participants"
tags1 <- tagIT(gregexpr(wordXsYro(structfind[["descriptivePatterns"]][["Human1"]],structfind[["descriptivePatterns"]][["Human2"]]),txtfiles,perl=TRUE),txtfiles)
tags2 <- tagIT(gregexpr(wordXsYr(tags1,structfind[[7]][["int_pat"]]),txtfiles,perl=TRUE),txtfiles)
tags3 <- tagIT(gregexpr(structfind[["descriptivePatterns"]][["N"]],txtfiles,perl=TRUE),txtfiles)
tags  <- unlist(c(tags2,tags3))
meta[[6]] <- wordXX(tags)
#tidyIT(matchIT(meta[[6]],txtfiles))


# 3b. COLLECT --------------------------------------------------------------
# Now get the meta text and locations (will be stored in attribute "location")
metaFound <- llply(meta,matchIT,IT=txtfiles)


# 4. PLACE ------------------------------------------------------------------
# Get a datafile with stats and their location 
extracted          <- ldply(metaFound,placeIT,skeleton=skeleton)
extracted$pdf.id   <-as.numeric(as.vector(extracted$pdf.id))
extracted$location <-as.numeric(as.vector(extracted$location))

data        <- extracted[ do.call(order,extracted), ] 
data$name   <- factor(data$pdf.id,levels=seq(along=fnames),labels=names(txtfiles))
data        <- subset(data, location!=-1)
  
rm(extracted)
write.csv(data,file="autostats.csv")
  

  # # Article level stats -----------------------------------------------------
  # alpha     <- .05
  # Stats$sig <- 0
  # Stats$sig[Stats$Computed <= alpha] <- 1 
  # 
  # # Summarise article level stats
  # studysum <- summary(Stats)
  # 
  # # Number of Sig. p values extracted per article and total
  # pValues.Sig <- c(ddply(Stats,"Source",function(x) sum(x$sig, na.rm=TRUE))[,2],sum(Stats$sig, na.rm=TRUE))
  # # Average post-hoc power (1 - Computed p value) and total
  # posthoc.Pwr <- c(ddply(Stats,"Source",function(x) mean((1-x$Computed), na.rm=TRUE))[,2],mean((1-Stats$Computed), na.rm=TRUE))
  # # Incredibility Index
  # BonesPsi <- pValues.Sig / studysum$pValues
  # IC       <- BonesPsi
  # 
  # power.t.test(n=100,delta=.23)
  # 
  # write.csv(stats,file="autostats.csv")
  