# [sciCure](http://fredhasselman.github.io/scicuRe/)
# ========
# ### A Toolbox for Curating Scientific Knowledge Extracted from Scientific Publications
# #### Tool development for the Curate Science website
# 
# *Very BETA version - [Fred Hasselman](http://fredhasselman.com)*
# 
# All the functions that are required are available in `scicuRe_source.R`, a package will follow soon. 
# Use the code below to source it directly from GitHub (the source function was found [here](http://tonybreyal.wordpress.com/2011/11/24/source_https-sourcing-an-r-script-from-github/)). It requires you to install and load the `RCurl` package
# 
# ```
# source_https <- function(url, ...) {
#   # load the package 
#   require(RCurl)
# 
#   # parse and evaluate each .R script
#   sapply(c(url, ...), function(u) {
#     eval(parse(text = getURL(u, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), envir = .GlobalEnv)
#   })
# }
# 
# # Source the scicuRe_source.R toolbox!
# source_https("https://raw.github.com/FredHasselman/SC-stats/master/scicuRe_source.R")
# ```
# 
# **SEARCH STRATEGY**
# 
#  - Need stats labelled for each experiment in a multi-experiment study
#  - Need measures of effect size
#  - Need aggregate stats, i.e. based on post-hoc power and more
#  - Need information to decide if proper analysis was used on stats ([`statcheck`](https://github.com/MicheleNuijten/statcheck) provides some interesting tests!)
#     
# *Therefore:*    
#  - Search is implemented to be driven by the hierarchical structure of an empirical report
#  - Structure is: Header (e.g., Experiment 1) -> Section (e.g., Results) -> [Table 1. ->] *STAT OF INTEREST* <- Page number [<- Appendix]
#  - Extract as much information as possible!
#     
# The `statcheck` package bij M. Nuijten is available [here](https://github.com/MicheleNuijten/statcheck)   
# Download, unzip and install, e.g. by running:    
# ```
# install.packages("~/Downloads/statcheck-master/", repos = NULL, type="source")
# ```   
# 
# The authors of `statcheck` note:   
# > The `pdftotext` program (http://www.foolabs.com/xpdf/download.html) is used to convert PDF files to plain text files. This must be installed and PATH variables must be properly set so that this program can be used from command line.
#    
# For the code in this package to execute correctly `pdftotxt` and `pdfinfo` need to be available on your system as well.   
# Check it now:
# ```
# if(all(file.exists(Sys.which("pdftotext")))) print("YES!") else print("NO!")
# ```   

# PREPARE -----------------------------------------------------------------
require(plyr)
require(reshape)
setwd("/Volumes/Fred HD/Rplus/gitR/SciComStats")
PDFpath <- paste(getwd(),"PDFfolder/",sep="/")


# Get Text and PDFinfo of PDFs in PDFpath ---------------------------------
fnames          <- list.files(PDFpath,pattern=".pdf",full.names=TRUE)
txtfiles        <- sapply(fnames,getPDFs)
names(txtfiles) <- gsub(".pdf","",basename(fnames))

# Load predefined search patterns
structfind <- loadPatterns()


# Build a skeleton for each article ---------------------------------------

# 1. Get location of article structure indices
element         <- c("header","section","table","figure")
skeleton        <- ldply(seq(along=element), function(e) cbind(melt(gregexpr(structfind[["skeletonPatterns"]][[element[e]]],txtfiles,perl=T,ignore.case=T)),name="",element=element[e]))
names(skeleton) <- c("location","pdf.id","name","element")
# 2. Find header info
skeleton$header <- NA
skeleton$header[skeleton$location!=-1] <- unlist(llply(seq(along=element), function(e) matchIT(structfind[["skeletonPatterns"]][[element[e]]],txtfiles,ic=T,cln=F)))
# 3. Find/infer and add pagenumbers to the skeleton
pagenumbers     <- ldply(seq(along=txtfiles), function(f){getPageSeq(structfind[["skeletonPatterns"]][["page"]],fnames[[f]],txtfiles[[f]],f)})
skeleton        <- rbind(skeleton, pagenumbers)
# 4. Add names, tidy up, sort results within each article
skeleton$name   <- factor(skeleton$pdf.id,levels=seq(along=fnames),labels=names(txtfiles))
skeleton        <- subset(skeleton, location!=-1)
skeleton        <- skeleton[ order(skeleton$pdf.id,skeleton$location), ]


# META pattern search -----------------------------------------------------

# 1. Collect the metapatterns for all files
# 2. Find where they "belong" by preferring stats with expected locations in the article:
#    Major Header -> Section header -> [Table 1. ->] *STAT OF INTEREST* <- Page number [<- Appendix]
# 3. Connect sample info to stats Major Header -> Section header = *PARTICIPANTS* -> [Table 1. ->] *MEAN, SD* -> *STAT OF INTEREST* <- Page number [<- Appendix]
# 4. Calculate!


# Most results have a structure:
# [name of test statistic] [df1..n] [value of statistic] [p-value] [effect size] 
#
# In the neigbourhood of the reported test we may encounter:
# [Mean/Median] [SD/SE] [N] [CI] 
#
# Also:
# [one/two tailed/...] [bonferroni/FWE/...] [levene/equal variances assumed/...]

# Let's call them meta-patterns
meta      <- vector("list",length=6)
metafound <- vector("list",length=6)
# Here are 6 meta patterns commonly encountered in the behavioural sciences:

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

metafound <- llply(meta, function(m) tidyIT(matchIT(m,IT=txtfiles)))
metaloc   <- llply(meta, function(m) locIT(m,IT=txtfiles))

# Assign skeleton to meta -------------------------------------------------

for(m in seq(along=txtfiles){
  skeletor <- subset(skeleton, pdf.id==m)
           #subset(pagenumbers, pdf.id==1)$location
  
  id<-ldply(as.vector(locIT(meta[[1]],IT=txtfiles[[1]])[[1]]),function(s) max(which(subset(pagenumbers, pdf.id==1)$location<s)))
  cbind(tidyIT(matchIT(meta[[1]],IT=txtfiles[[1]]))[[1]],subset(pagenumbers, pdf.id==1)$location[id$V1[1]])
  
 >skeletor$location
  llply(metaloc,)
}


# Calculate ES based on test-statistics ------------------------------------
out <- ldply(1:nrow(Stats), function(r) getESCI(Stats[r,]))




# Find (1st) location of the raw statistics extracted by "statcheck"
for(i in seq(along=files)){
 rawstats[[i]]  <- ldply(levels(Stats$Raw)[Stats$Raw[Stats$Source==titles[i]]], function(rs) regexpr(rs,files[i],fixed=TRUE)[[1]])  
 #rawstats[[i]] <- cbind(stat=Stats$Raw[Stats$Source==titles[i]],pos=raw_pos)
}
x <- matrix(rnorm(100), nrow = 5)
dist(x,method="maximum")
c(hdr_pos[[1]], sec_pos[[1]], page_pos[[1]], tab_pos[[1]])

d-rawstats

# Suggest a location for the stats in hierarchical order:
# Header > Section > Table > Page based on minimal distance to stats label
labels <- getLabels(raw_pos,hdr_pos,page_pos,tab_pos,txt=files[i],levels(titles)[i])
  
rm(hdr_pos,page_pos,tab_pos,fig_pos,raw_pos)

labs    <- ldply(labels)
outfile <- as.data.frame(cbind(labs,out,ErrorDiag[,2:10]))

# Article level stats -----------------------------------------------------
alpha     <- .05
Stats$sig <- 0
Stats$sig[Stats$Computed <= alpha] <- 1 

# Summarise article level stats
studysum <- summary(Stats)

# Number of Sig. p values extracted per article and total
pValues.Sig <- c(ddply(Stats,"Source",function(x) sum(x$sig, na.rm=TRUE))[,2],sum(Stats$sig, na.rm=TRUE))
# Average post-hoc power (1 - Computed p value) and total
posthoc.Pwr <- c(ddply(Stats,"Source",function(x) mean((1-x$Computed), na.rm=TRUE))[,2],mean((1-Stats$Computed), na.rm=TRUE))
# Incredibility Index
BonesPsi <- pValues.Sig / studysum$pValues
IC       <- BonesPsi

power.t.test(n=100,delta=.23)

write.csv(stats,file="autostats.csv")





# 
# # Find labels such as study header and page number for each statistic 
# getLabels <- function(raw_pos,hdr_pos,page_pos,tab_pos,txt,title){
#   
#   Ori <- names(raw_pos)
#   page_label <- rep(NA,times=length(Ori))
#   hdr_label  <- rep(NA,times=length(Ori))
#   tab_label  <- rep(NA,times=length(Ori))
#   
#   if(attr(page_pos,"match.length")!=-1){
#   # Get page numbers, discard any numbers that destroy sequential order
#   Numbers  <- gsub("(^\\D*)|(\\D*$)","",substring(txt,page_pos,page_pos+attr(page_pos,"match.length")-1))
#   Pages    <- Numbers[which(abs(diff(as.numeric(Numbers)))<=2)]
#   while(any(as.numeric(Numbers)>max(as.numeric(Pages)))){
#     Numbers  <- c(Pages,Numbers[as.numeric(Numbers)>max(as.numeric(Pages))])
#     id       <- which(abs(diff(as.numeric(Numbers)))<=2)
#     Pages    <- Numbers[c(id,max(id)+1)]
#   }
#   Numbers  <- gsub("(^\\D*)|(\\D*$)","",substring(txt,page_pos,page_pos+attr(page_pos,"match.length")-1))
#   page_pos <- page_pos[which(Numbers%in%Pages)]
#   
#   # Decide on which page a statistic is printed
#   page_diff  <- as.matrix(sapply(raw_pos,function(d) d-page_pos))
#   page_label <- sapply(1:ncol(page_diff), function(d) ifelse(any(page_diff[,d]<0),{
#     Pages[which(page_diff[,d]==(max(page_diff[page_diff[,d]<0,d])))]},{
#       Pages[which(page_diff[,d]==(min(page_diff[page_diff[,d]>=0,d])))]
#     }))
#   }
#   
#   if(attr(hdr_pos,"match.length")!=-1){
#   # Decide which Experiment or Study header belongs to the statistics
#   Headers <- gsub("(^[^[:print:]]*)|([^[:print:]]*$)","",substring(txt,hdr_pos,hdr_pos+attr(hdr_pos,"match.length")-1))  
#   # Decide which section a statistic belongs to
#   hdr_diff  <- as.matrix(sapply(raw_pos,function(d) d-hdr_pos))
#   hdr_label <- sapply(1:ncol(hdr_diff), function(d)  ifelse(any(hdr_diff[,d]>=0),{
#     Headers[which(hdr_diff[,d]==min(hdr_diff[hdr_diff[,d]>=0,d]))]},{
#       Headers[which(hdr_diff[,d]==min(hdr_diff[hdr_diff[,d]<0,d]))]
#     }))
#   }
#   
#   if(attr(tab_pos,"match.length")!=-1){
#   # Decide if it is likely stats were extracted from a table
#   Tables <- gsub("(^[^[:print:]]*)|([^[:print:]]*$)","",substring(txt,tab_pos,tab_pos+attr(tab_pos,"match.length")-1))  
#   # Decide which table a statistic may have come from
#   tab_diff  <- as.matrix(sapply(raw_pos,function(d) d-tab_pos))
#   tab_label <- sapply(1:ncol(tab_diff), function(d)  ifelse(any(tab_diff[,d]>=0),{
#       Tables[which(tab_diff[,d]==min(tab_diff[tab_diff[,d]>=0,d]))]},{
#         rep(NA,times=length(Tables))
#         }))
#   }
#  
#   return(as.data.frame(cbind(title,Ori,page_label,hdr_label,tab_label)))
# }
