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

# META functions ----------------------------------------------------------

structureIT <- function(fnames,IT){
  # Load predefined search patterns
  structfind <- loadPatterns()
  
  # Build a skeleton for each article ---------------------------------------
  
  # 1. Get location of article structure indices
  element         <- c("header","section","table","figure")
  skeleton        <- ldply(seq(along=element), function(e) cbind(melt(gregexpr(structfind[["skeletonPatterns"]][[element[e]]],IT,perl=T,ignore.case=T)),name="",element=element[e]))
  names(skeleton) <- c("location","pdf.id","name","element")
  # 2. Find header info
  skeleton$header <- NA
  skeleton$header[skeleton$location!=-1] <- unlist(llply(seq(along=element), function(e) matchIT(structfind[["skeletonPatterns"]][[element[e]]],IT,ic=T,cln=F)))
  # 3. Find/infer and add pagenumbers to the skeleton
  pagenumbers     <- ldply(seq(along=IT), function(f){getPageSeq(structfind[["skeletonPatterns"]][["page"]],fnames[[f]],IT[[f]],f)})
  skeleton        <- rbind(skeleton, pagenumbers)
  # 4. Add names, tidy up, sort results within each article
  skeleton$name   <- factor(skeleton$pdf.id,levels=seq(along=fnames),labels=names(IT))
  skeleton        <- subset(skeleton, location!=-1)
  skeleton        <- skeleton[ order(skeleton$pdf.id,skeleton$location), ]
  return(skeleton)
}

IT<-metaFound[[6]]
p=1
placeIT <- function(IT,skeleton){
  out <- vector("list")
  
  cnt=0
  for(p in seq(along=IT)){
    if(length(IT[[p]])>0){
    tmp1 <- melt(IT[[p]])
    names(tmp1)[1]<- "header"
    tmp1$location <- attr(IT[[p]],"location")
    tmp1$element  <- "metastat"
    skel.tmp      <- rbind(subset(skeleton, pdf.id==p, select=c(location,element,header)))
    tmp2          <- rbind(skel.tmp,tmp1)
    tmp2$location <- as.numeric(as.vector(tmp2$location))
    tmp2          <- tmp2[ do.call(order,tmp2), ]
    metaID        <- which(tmp2$element=="metastat")
    
    for(i in metaID){
      idAbove   <- which(tmp2$location<tmp2$location[i])
      idHeader  <- idAbove[which(tmp2$element[idAbove]=="header")]
      if(length(idHeader) > 0){idHeader<-max(idHeader)} else {idHeader<-nrow(tmp2)+1}
      idSection <- idAbove[which(tmp2$element[idAbove]=="section")]
      if(length(idSection) > 0){idSection<-max(idSection)} else {idSection<-nrow(tmp2)+1}
      idTable   <- idAbove[which(tmp2$element[idAbove]=="table")]
      if(length(idTable) > 0){idTable<-max(idTable)} else {idTable<-nrow(tmp2)+1}
      idBelow  <- which(tmp2$location>tmp2$location[i])
      idPage   <- idBelow[which(tmp2$element[idBelow]=="page")]
      if(length(idPage) > 0){idPage<-min(idPage)} else {idPage<-nrow(tmp2)+1}

      tmp <- rbind(tmp2,cbind(location=-1,element="none",header="NOT FOUND"))
      cnt=cnt+1
      out[[cnt]] <- as.data.frame(cbind(pdf.id=p,location=tmp$location[i],stat=tmp$header[i],page=tmp$header[idPage],section=tmp$header[idSection],table=tmp$header[idTable],header=tmp$header[idHeader],name=names(IT)[p]),row.names=paste(p,cnt,sep="."),stringsAsFactors = FALSE)
      rm(idAbove,idHeader,idSection,idTable,idBelow,idPage,tmp)
    }
    } else {
      cnt=cnt+1
      out[[cnt]] <- as.data.frame(cbind(pdf.id=p,location=-1,stat=paste(i,sep=""),page="",section="",table="",header="",name=names(IT)[p]),row.names=paste(p,cnt,sep="."),stringsAsFactors = FALSE)
    }
  }
  return(ldply(out))
}

# Extraction functions ----------------------------------------------------

getPDFs <- function(x)
{# Based on the function "getPDF()" in the package "statscheck" (https://github.com/MicheleNuijten/statcheck)
  if(all(file.exists(Sys.which(c("pdftotext"))))){
    
    txtfiles <- character(length(x))
    
    for (i in 1:length(x))
    {
      system(paste('pdftotext -q -enc "UTF-8" "',x[i],'"',sep=""))
      if (file.exists(gsub("\\.pdf","\\.txt",x[i])))
      {
        fileName <- gsub("\\.pdf","\\.txt",x[i])
        txtfiles[i] <- readChar(fileName, file.info(fileName)$size)            
        #         info <- system(paste('pdfinfo "',x[i],'"',sep=""), intern=TRUE)
        #         txtfiles[[i]] <- list(text=text,info=info)
      } else
      {
        warning(paste("Failure in file",x[i]))
        txtfiles[i] <- ""
      }
    }}
  else {
    warning("Could not find 'pdftotext'...")
  }
  return(enc2utf8(txtfiles))
}

returnPDFdir <- function(dir, ...){
  # Function to check directory of PDFs
  # Based on the function "checkPDFdir()" from package "statchek" (https://github.com/MicheleNuijten/statcheck)
  if (missing(dir)) dir <- tk_choose.dir()
  files <- list.files(dir,pattern=".pdf",full.names=TRUE)
  txts <- character(length(files))
  message("Importing PDF files...")
  pb <- txtProgressBar(max=length(files),style=3)
  for (i in 1:length(files))
  {
    txts[i] <-  getPDFs(files[i])    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  names(txts) <- gsub(".pdf","",basename(files))
  return(txts)
}

getPDFpageinfo <- function(pdfname) {
  if(file.exists(Sys.which(c("pdfinfo")))){
    if(file.exists(pdfname)){
      x <- system(paste('pdfinfo "',pdfname,'"',sep=""), intern=TRUE)
    } else {
      warning("Could not find PDF file...")
    }
  } else {
    warning("Could not find 'pdfinfo'...")
  }
  return(as.numeric(sub("\\D+","",x[grep("(?<=Pages\\:)(\\s*)(\\d+)",x,perl=T)],perl=T)))
}

loadPatterns <- function(skelet=T,stats=T,evid=T,eff=T,desc=T,bib=T){
  # Ready to use regexpr patterns to find common sections of articles, stats, etc.
  
  # Greek   small   capital     alt.
  #
  # alpha    U+03B1
  # beta     U+03B2
  # gamma    U+03B3
  # delta    U+03B4    U+0394   
  # epsilon  U+03B5
  # zeta     U+03B6
  # eta      U+03B7
  # phi      U+03C6    U+03A6     U+03D5
  # chi      U+03C7    U+03A7     U+03BE
  # psy      U+03C8    U+03A8
  # omega    U+03C9    U+03A9
  
  # Characters encountered after 'pdftotxt' encoding for minus sign
  sign_pat<- "([-–ϪÀ±~])?"
  # Characters encountered after 'pdftotxt' encoding for (in)equality
  eq_pat  <- "([=<>≤≥≦≧≨≩≪≫≭≮≯≰≱≲≳≴~≈≉≠≡≢≣ϭ5 ])"
  # Characters encountered after 'pdftotxt' encoding for non significant
  ns_pat  <- "((n(\\.)?s(\\.)?)|(no(n|\\s+)*sig(\\.)?))"
  # Numbers
  int_pat <- paste("((\\d+","(\\,)",")*\\d+)",sep="[[:blank:]]*")
  real_pat<- paste("((\\d+(\\,))*(\\d*(\\.))*\\d+)",sep="[[:blank:]]*")
  
  skeletonPatterns <- list(
    # Ordered header indicators..."Experiment 1", "Study 4b", "Appendix B"
    header="(?<=\\n(\\f|\\n))(?<header>(Replication(s)*|Meta[-]?|Data|Exploratory|Confirmatory|Computer|Model|Network)*[[:blank:]]*(Experiment(s*)|Stud(y|ies)|Analys[ie]s|Appendi(x|ces))[[:blank:]]*(\\d|\\w)+[[:blank:]]*[[:punct:]]*[[:blank:]]*)(?<title>[[:print:]]{0,150})(?=\\n(\\f|\\n)*)",
    # Common section / paragraph indicators..."Method", "Procedure", "Results"
    section="(?<=\\n(\\f|\\n))(?<section>(Abstract|Summary|(General[[:blank:]]*)?Introduction)|Method(s*)|Procedure(s)*|(Participant|Subject|Sample|Patient|Client)(s)*|Result(s)*|(General[[:blank:]]*)?((Conclusion|Result)(s)*)?[[:blank:]]*(and|&)*[[:blank:]]*((Discussion|Result)(s)*)?|Material(s)*|Instrumen(s)*[[:print:]]{0,3})(?=\\n(\\f|\\n)*)",
    # Tables + titles
    table="(?<=\\n(\\n|\\f))(?<table>Table(s)*[[:blank:]]*(\\d|\\w)+[[:blank:]]*[[:punct:]]*[[:blank:]]*)(?<title>[[:print:]]*)(?=\\n(\\n|\\f)*)",
    # Figures + titles
    figure="(?<=\\n(\\n|\\f))(?<figure>(Figure|Box|Graph|Picture|Display|Diagram|Photo|Frame)(s)*[[:blank:]]*(\\d|\\w)+[[:blank:]]*[[:punct:]]*[[:blank:]]*)(?<title>[[:print:]]*)(?=\\n(\\n|\\f))",
    # Page numbers
    page="(?<=\\n(\\n|\\f))(?<page>\\w{0,2}\\d+\\w{0,2})(?=\\n(\\n|\\f)+)"
  )
  
  statsPatterns <- list(
    # Regexpr patterns to find statistics
    t = paste("\\b[tT]","[s]?","([(]",int_pat,"[)])?",eq_pat,sign_pat,"(",real_pat,"|",ns_pat,")",sep="[[:blank:]]*"),
    F = paste("\\b[fF]","[s]?","([(,]",int_pat,"[)])+","",eq_pat,"(",real_pat,"|",ns_pat,")",sep="[[:blank:]]*"),
    r = paste("\\b[rρΡ]","[s]?","([(]",int_pat,"[)])?",eq_pat,sign_pat,"(",real_pat,"|",ns_pat,")",sep="[[:blank:]]*"),
    Z = paste("\\b[zZ]","[s]?",eq_pat,sign_pat,real_pat,"(",real_pat,"|",ns_pat,")",sep="[[:blank:]]*"),
    W = paste("\\b[wW](ald)*",eq_pat,"(",real_pat,"|",ns_pat,")",sep="[[:blank:]]*"),
    U = paste("\\b[uU]",eq_pat,"(",real_pat,"|",ns_pat,")",sep="[[:blank:]]*"),
    chi2 = paste("\\b[χΧx]","[2²]?","[s]?","[(]",int_pat,"[)]",eq_pat,"(",real_pat,"|",ns_pat,")",sep="[[:blank:]]*")
  )
  
  evidencePatterns <- list(
    p   = paste("\\b(",ns_pat,")|([pP]","(rep)?","[s]?)",eq_pat,sign_pat,"(",real_pat,"|",ns_pat,")",sep="[[:blank:]]*"),
    ll  = paste("\\b(",sign_pat,"[2*]*(ll)*((log(\\s)*)*lik(elihood))*(d)*)",eq_pat,sign_pat,"(",real_pat,"|",ns_pat,")",sep="[[:blank:]]*"),
    inf = paste("\\b((A|B|D)IC|\\s*information\\s*criterion)",eq_pat,sign_pat,"(",real_pat,"|",ns_pat,")",sep="[[:blank:]]*")
  )
  
  effectPatterns <- list( 
    # Effect sizes
    omega2 = paste("\\b[ωΩ][2²]?",eq_pat,real_pat, sep="[[:blank:]]*"),
    eta2 = paste("\\b(partial)?","[ηΗgZ]","[pₚ]?","[2²]?","[pₚ]?",eq_pat,"(",real_pat,"|",ns_pat,")",sep="[[:blank:]]*"),
    phi = paste("\\b([cC]ramer","(['`]","[s])?)?","([φΦ]|[pP]hi)","[cC]*",eq_pat,sign_pat,real_pat, sep="[[:blank:]]*"),
    ICC = paste("\\b[iIcC]+",eq_pat,sign_pat,real_pat, sep="[[:blank:]]*"),
    psi = paste("\\b([Ψψ]|[pP]s(i|y))",eq_pat,sign_pat,real_pat, sep="[[:blank:]]*"),
    r2 = paste("\\b[rρΡ][2²]?",eq_pat,real_pat, sep="[[:blank:]]*"),
    OR = paste("\\b(l(n|og)[(]*",")*[oO]([rR]|dds)","[)]*",eq_pat,sign_pat,real_pat, sep="[[:blank:]]*"),
    RR = paste("\\b(l(n|og)[(]*",")*[rR]","[)]*",eq_pat,sign_pat,real_pat, sep="[[:blank:]]*"),
    gd = paste("\\b([gG]lass","(['`]","[s])?)?","[δΔd]",eq_pat,real_pat, sep="[[:blank:]]*"),
    hv = paste("\\b([hH]edges","(['`]","[s])?)?","[gG]",eq_pat,real_pat, sep="[[:blank:]]*"),
    cd = paste("\\b([cC]ohen","(['`]","[s])?)?","[δΔd]",eq_pat,real_pat, sep="[[:blank:]]*"),
    cf = paste("\\b([cC]ohen","(['`]","[s])?)?","[fF]","[2²]?",eq_pat,real_pat, sep="[[:blank:]]*")
    )
  
  # Regexpr patterns to find information about analyses
  descriptivePatterns <- list(
    # Sample descriptives
    sde = paste("\\b(([σΣ])|([sS](\\.)?([dD]|[eE])(\\.)?))+",eq_pat,real_pat, sep="[[:blank:]]*"),
    var = paste("\\b(([σΣ])|([sS](\\.)?[dD](\\.)?)[2²]?)+",eq_pat,real_pat, sep="[[:blank:]]*"),
    M   = paste("\\b(([μΜ])|([mM](\\.)?[nN]?(\\.)?)|([mM]ean|[aA]verage|[mM]edian))+[s]?",eq_pat,sign_pat,real_pat, sep="[[:blank:]]*"),
    N = paste("(\\b[nN]",eq_pat,")",int_pat,sep="[[:blank:]]*"),
    # Human sample info
    Human1 = "(m(e|a)n|wom(e|a)n|girl(s)*|boy(s)*|male(s)*|female(s)*|college|university|(high|pre|elementary)(-|\\s)school)*[[:blank:]]*(participant(s)*|undergraduate(s)*|graduate(s)*|subject(s)*|client(s)*|patient(s)*|student(s)*|adult(s)*|child(ren)*|kindergartner(s)*|preschooler(s)*|teen(ager)?(s)*|adolescent(s)*|toddler(s)*|infant(s)*|neonate(s)*|m(e|a)n|wom(e|a)n|girl(s)*|boy(s)*|male(s)*|female(s)*|college|university|(high|pre|elementary)(-|\\s)school)[[:blank:]]*(m(e|a)n|wom(e|a)n|girl(s)*|boy(s)*|male(s)*|female(s)*|college|university|(high|pre|elementary)(-|\\s)school)*",
    Human2 = "((American\\sIndian|Alaska\\sNative|East\\sAsian|South\\sAsian|Native\\sHawaiian|Pacific\\sIslander|Hispanic|Latino|Puerto\\sRican|Black|African\\sAmerican|White|Caucasian|European|Russian|Eastern\\sEuropean|Indian|Pakistani|Dutch|German|Swedish|British|Norwegian|Danish|Polish|French|Spanish|Italian|Greek|Portugese|South\\sAmerican|Israeli|African|Japanese|Chinese|Australian|Mixed(\\sRace)*|Arabic|Muslim|Christian|Jewish|Hindu|Buddhist|Orthodox|Catholic|Reformed|Protestant|Humanist)(s)*)",
    # Models
    models = "(hierachical|(log|non|general(ised)?)([-]|\\s*)?|linear|multi(ple|level)|stepwise|logistic|probit|poisson|categorical|polynomial|ridge|robust|Cox|survival|cluster|factor|(principal\\s)*component|path|dynamic(al)*|(multi|mono)[-]?fractal|correlation(al)*|spectral|fourier|wavelet|conjoint|(lag)?[-]sequential|)[[:blank:]]*(regression|ANOVA|ANCOVA|MANOVA|GLM|GLMM|ARMA|ARiMA|ARfiMA|GARCH|SEM|HLM|([(]?C[)]?)?RQA|(MF)?[-]?DFA|SDA|WMMT|ROC|MDS)*[[:blank:]]*(model([l]ing)|analys[ie]s|equation(s)*|simulation(s)*)",
    # Adjustments for multiple comparisons
    adjusted="\\b(Bonferroni|Holm|Hochberg|Hommel|Benjamini|Yekutieli|Sidak|Shaffer|Simes|Familywise|FWE|LSD|(False\\s*Discovery)|FDR|(Multiple\\s*Comparisons)|((Adjusted|Corrected)\\s*P[-]?\\s*value(s)*))",
    # Direction / tailed / sided -ness of tests
    directed="[[:blank:]]*[(]?[[:blank:]]*((Un)*Directed|((\\d+|one|two)?\\s*[-]?\\s*sided)|((\\d+|one|two)?\\s*[-]?\\s*tailed))[[:blank:]]*[)]?[[:blank:]]*",
    dependend="[[:blank:]]*[(]?[[:blank:]]*((((In)*Dependent|Paired|Matched)\\s*(Sample(s)*|Group(s)*|Pairs|Observation(s)*|Measurement(s)*|(\\w*[-]?test(s)*))+)|(Repeated\\s*Measure(s|ment(s)*)*)|((within|between)\\s*(subject(s)*|participant(s)*|factor(s)*|design|predictor(s)*|covariate(s)*|level(s)*|effect(s)*)))[[:blank:]]*[)]?[[:blank:]]*",
    # Randomisation
    random="[[:blank:]]*[(]?[[:blank:]]*((Assign(ed|ment)*|Select(ed|ion)*)?(At)*Random(ly)*\\s*(Assign(ed|ment)*|Select(ed|ion)*)?)[[:blank:]]*[)]?[[:blank:]]*"
    )
  
  # Regexpr patterns to find reference list and citation info
  bibPatterns <- list(
    abstract = "",
    authors  = "",
    citation = "",
    reflist  = "",
    pubdates = "",
    keywords = "",
    authorloc  = "",
    authornote = "",
    funding  = ""
  )
  
  return(list(skeletonPatterns=skeletonPatterns,
              statsPatterns=statsPatterns,
              evidencePatterns=evidencePatterns,
              effectPatterns=effectPatterns,
              descriptivePatterns=descriptivePatterns,
              bibPatterns=bibPatterns,
              c(sign_pat=sign_pat,eq_pat=eq_pat,ns_pat=ns_pat,int_pat=int_pat,real_pat=real_pat)
              )[c(skelet,stats,evid,eff,desc,bib,TRUE)]
         )
}

splitseq <- function(seq){# Split a vector of numbers into sequences of lag = 1 (if any)
  seq <- seq[seq>0]
  cutat  <- c(1,(which(diff(sort(seq),lag=1)!=1)+1))
  if(length(cutat)>0){
    cutsz   <- diff(c(cutat,(length(seq)+1)))
    splitz  <- split(sort(seq),unlist(lapply(seq(along=cutat),function(s) rep(cutat[s],times=cutsz[s]))))
    seqname <- sapply(seq(along=splitz),function(s) paste("Split sequence ",s,sep=""))
    names(splitz) <- seqname
  } else {
    splitz <- seq
    if(length(cutat)==1){
      names(splitz) <- "Ordered vector is 1 sequence!"
    }
    names(splitz) <- "No sequences in this vector"
  }
  return(splitz)
}

# SEARCHERS ---------------------------------------------------------------

cleanIT <- function(IT){# Clean up the output from matchIT()
  regpat <- loadPatterns(rep(FALSE,times=6))
  for(i in seq(along=IT)){
    # Clean the minus sign
    IT[i] <- gsub(paste("(\\s*([-–ϪÀ±~])\\s*)(?=",regpat[[1]][5],")",sep="")," -", IT[i], perl=TRUE)
    
    # Clean the stat symbols etc.
    IT[i] <- gsub("\\b[χΧx][[:blank:]]*[2²]?[[:blank:]]*[s]?[[:blank:]]*","χ²", IT[i], perl = TRUE)
    IT[i] <- gsub("\\b[ηΗgZ][[:blank:]]*[2²][[:blank:]]*","η²", IT[i], perl = TRUE)
    IT[i] <- gsub("\\b[ωΩ][[:blank:]]*[2²]?[[:blank:]]*","ω²", IT[i], perl = TRUE)
    IT[i] <- gsub(paste("\\b([cC]ramer","(['`]","[s]\\s*)?)?","([φΦ]|[pP]hi)","[cC]*", sep="[[:blank:]]*"),"\\1φ",IT[i],perl = TRUE)
    IT[i] <- gsub("\\b[fF][[:blank:]]*(?=[(])","F", IT[i], perl = TRUE)
    IT[i] <- gsub("\\b[tT][[:blank:]]*(?=[(])","t", IT[i], perl = TRUE)
    IT[i] <- gsub(paste("\\b[zZ][[:blank:]]*(?=",regpat[[1]][2],")",sep=""),"Z", IT[i], perl = TRUE)
    IT[i] <- gsub(paste("(\\b[pP][[:blank:]]*)([[:blank:]]*rep[s]*)*(?=",regpat[[1]][2],")",sep=""),"p\\2 ",IT[i], perl = TRUE)
    
    # Clean up spaces
    IT[i] <- gsub("(^[[:blank:]]+)|([[:blank:]]+$)|((?<=[(])[[:blank:]]+)","", IT[i], perl = TRUE)
    
    # Clean up equal sign
    IT[i] <- gsub("(?<=[[:graph:]])([[:blank:]]*(=|ϭ|(\\s*5\\s))[[:blank:]]*)(?=[[:digit:]\\.\\,nNsS-])"," = ",IT[i], perl=TRUE)
  }
  return(list(IT))
}

tidyIT <- function(IT,cln=FALSE){# Tidy up an IT)
  if(cln){cleantIT(IT)}
  for(i in seq(along=IT)){
    IT[[i]] <- gsub("^[[:punct:]]+|[[:punct:]]+$","",IT[[i]], perl = TRUE)
    IT[[i]] <- gsub("(^[[:blank:]]+)|([[:blank:]]+$)|((?<=[(])[[:blank:]]+)","", IT[[i]], perl = TRUE)
  }
  return(IT)
}

regIT   <- function(str,IT) {# A wrapper for grep, searches for string in source
  grep(str,IT,ignore.case=TRUE,value=TRUE)
  return(IT)
}

locIT <- function(str,IT) {# Return matching positions in source using regmatches
  n <- names(IT)
  m <- gregexpr(str,IT,perl=TRUE,ignore.case=TRUE) 
  names(m) <- n
  #print(regmatches(IT,m))
  return(m)
}

matchIT <- function(str,IT, ic=TRUE,cln=TRUE,loc=TRUE) {# Return matching strings in source using regmatches
  require(reshape)
  n <- names(IT)
  m <- gregexpr(str,IT,perl=TRUE,ignore.case=ic) 
  t <- character(length(IT))
  for(i in seq(along=IT)){
    tmp <- regmatches(IT[[i]],m[i])
    ifelse(cln,{t[i] <- cleanIT(tmp[[1]])},{t[i] <- list(tmp[[1]])})
    if(loc==TRUE){attr(t[[i]],"location") <- as.vector(melt(m[[i]])$value)}
  }
  names(t) <- n
  return(t)
}

tagIT<-function(tags,IT){# Find tags returned from grepexpr in the source and return unique matches as a regex pattern to use in subsequent searches
  tmp<-unique(sub("[[:blank:]]*$","",unlist(regmatches(IT,tags)), perl=T))
  tmp<-sub("^[[:blank:]]*","",tmp, perl=T)
  regmatches(tmp,gregexpr("(\\n|\\f|\\r)",tmp,ignore.case=T))<-rep("\\n",length(tags))
  return(tmp)
}

# subIT  <- function(IT) {# Substitute tags in long character vector / corpus
#   if(any(which(IT$tags==""))){IT$tags<-IT$tags[which(IT$tags!="",arr.ind=T)]}
#   ifelse((length(IT$tags)==0),{
#     print("No tags in IT$tags !!!")
#     return(IT$txt)
#   },{
#     m <-gregexpr(IT$tags,IT$txt,perl=TRUE)
#     subs<-regmatches(IT$txt,m)
#     regmatches(IT$txt,m)<-rep(IT$str,length(m))
#     print(paste("Changed ",length(which(unlist(m)>0))," strings into: ",IT$str))
#     return(IT$txt)
#   })
# }


# Heuristics  -------------------------------------------------------------

inferIT <- function(signPat,sTemp,nExpected,nTol,IT){# Try to infer page range from entries such as 'XXX-YYY = #Pages' in the PDF text
  require(plyr)
  tmpSeq <- NULL 
  rng <- as.numeric((nExpected-nTol):(nExpected+nTol))
  nrs <- sub(paste(signPat,sep=")"),"-",unlist(matchIT(wordXnY(mid=paste(signPat,sep=""),x=sTemp,y=sTemp),IT,cln=F)))
  if(length(nrs)!=0){
    if((abs(ldply(as.quoted(nrs),eval))+1)%in%rng){
      tmpSeq1  <- as.numeric(strsplit(nrs[which((abs(ldply(as.quoted(nrs),eval))+1)%in%rng)],signPat)[[1]])
      tmpSeq   <- seq(tmpSeq1[1],tmpSeq1[2],by=1)
    }
  }
  return(tmpSeq)
}

getPageSeq <- function(str,PDFname,IT,PDFid,closeFITdiff=3,closeFITmiss=3){
  # Lot of seting up to do...
  require(plyr)
  # Get sign regex pattern, without the ?
  signPat <- sub("\\?","",loadPatterns(rep(FALSE,times=6))[[1]]["sign_pat"])
  # Get the look ahead/behind conditions
  prepost <- strsplit(str,"?<page>\\w{0,2}\\d+\\w{0,2}",fixed=T)[[1]]
  # Get info on number of pages in the pdf from 'pdfinfo'
  nPages <- getPDFpageinfo(PDFname)
  # Define some tools (requires plyr)
  f_l <- each(length)
  f_r <- each(min,max)
  f_d <- each(diff)
  # Sequence (s) found in IT using str
  sFound  <- sort(as.numeric(matchIT(str,IT,cln=F)[[1]]))
  # Length (n) of found sequence
  nFound  <- f_l(sFound)
  # Range (r) of found sequence
  rFound  <- f_r(sFound)
  # Difference (d) between pdfinfo and found sequence length
  dFound  <- f_d(f_r(c(nFound,nPages)))
  NOcigar <- ((nPages-closeFITdiff):(nPages+closeFITdiff))
  # Location of found page nr. in units of IT size
  locatePages <- function(seqF,norm=nchar(IT)){ldply(seqF,function(n) c(page=n, loc=as.vector(locIT(wordXX(pre=prepost[1],post=prepost[2],x=n),IT)[[1]]/norm)))}
  # Min and Max values for first and last pagenumbers
  minP     <- .11
  maxP     <- .94
  avg_page <- round(nchar(IT)/nPages)
  fst_page <- round(minP*nchar(IT))
  lst_page <- round(maxP*nchar(IT))
  
  # Were multiple sequences returned in the sFound?
  splitsFound <- splitseq(sFound)
  seqL        <- f_l(splitsFound)
  locSplits   <- llply(splitsFound,locatePages)
  
  # Check true if any missing info is found
  foundD  <- FALSE
  foundR  <- FALSE
  foundS  <- FALSE
  
  while(!any(c(foundD,foundR,foundS))){
    # First try to find a range of numbers based on page locations in the PDF file
    # Is a likely first or last pagenr. in the sequence?
    if((length(firstP <- ldply(locSplits,function(s) s[which(s$loc<=minP),]))==0)){
      firstP <- integer(0)} else {firstP <- firstP[which.min(firstP$page),]}
    if((length(lastP  <- ldply(locSplits,function(s) s[which(s$loc>=maxP),]))==0)){
      lastP <- integer(0)} else {lastP <- lastP[which.max(lastP$page),]}
    
    # Evaluate results
    if((length(firstP$page)!=0)&(length(lastP$page)!=0)){
      # Range was likely found
      foundR <- TRUE
      newSeq   <- (firstP$page:lastP$page)
      tempSeqs <- c(splitsFound[firstP[[".id"]]],splitsFound[lastP[[".id"]]])
    } else {
      if(length(lastP$page)>0){
        # Guess firstP based on average expected pagesize
        firstP  <- locatePages(sFound,norm=1)[1,]
        extra   <- round((firstP$loc-round((minP/2)*nchar(IT)))/avg_page)
        newSeq  <- ((firstP$page-extra):max(sFound))
        if(extra>0){tempSeqs<- list(((firstP$page-extra):(firstP$page-1)),sFound)
        } else {
          tempSeqs<- list(sFound)
        }
      } else {
        # Guess lastP
        lastP   <- locatePages(sFound,norm=1)[nFound,]
        extra   <- round((nchar(IT)-lastP$loc)/avg_page)
        newSeq  <- (min(sFound):(lastP$page+extra))
        if(extra>0){tempSeqs<- list(sFound,(min(sFound):(lastP$page+extra)))
        } else {
          tempSeqs<- list(sFound)
        }
      }
    }
    attr(tempSeqs,"Close")   <- (f_l(unique(unlist(tempSeqs)))%in%NOcigar)
    attr(tempSeqs,"Missing") <- newSeq[!(newSeq%in%unique(unlist(sFound)))]
    
    if(attr(tempSeqs,"Close")){foundD=TRUE} 
    if((length(attr(tempSeqs,"Missing"))>0)&(length(attr(tempSeqs,"Missing"))<=closeFITmiss)){foundR=TRUE}
    
    # Check how we are doing, try some freaky solutions
    # Check 1
    if(!(foundD|foundR)){
      vault <- list(newSeq,tempSeqs)
      # It is possible nPages of the PDF file <> actual publication pages
      # See if there is a page indication in the text based on tempSeqs (e.g. 1-10) that yields [nPages-closeFITdiff,nPages+closeFITdiff]
      newSeq  <- llply(tempSeqs,function(s) inferIT(signPat=signPat,sTemp=s,nExpected=nPages,nTol=closeFITdiff,PDFIT=IT))
      if(!all(is.null(newSeq))){
        tempSeqs <- newSeq[!sapply(newSeq,is.null)]
        attr(tempSeqs,"Close") <- (f_l(unique(unlist(tempSeqs)))%in%NOcigar)
        attr(tempSeqs,"Missing") <- newSeq[!(newSeq%in%unique(unlist(sFound)))]
      }
      if(attr(tempSeqs,"Close")){foundD=TRUE} 
      if((length(attr(tempSeqs,"Missing"))>0)&(length(attr(tempSeqs,"Missing"))<=closeFITmiss)){foundR=TRUE}
      
      if(all(sapply(newSeq,is.null))){
        newSeq <- vault[[1]]
        tempSeqs <- vault[[2]]
      }
    }
    
    #     # Check 2
    #     if(!(foundD|foundR)){
    #       if(seqL > 1){
    #         vault <- list(newSeq,tempSeqs)
    #         # 1. Are certain seq combinations a close fit with expected page range: ((sum(sequence lengths) - nPages) <= closeFITdiff)?
    #         # 2. Are those combinations sensible sequences (missing number in sequence <= closeFITmiss)?
    #         combi    <- unlist(llply(2:(seqL-1), function(c) combn(seq(along=splitsFound),c,simplify=F)),recursive=F)
    #         dCombi   <- ldply(combi, function(c) nPages-f_l(unlist(splitsFound[c])))
    #         tempList <- sapply(combi[abs(dCombi)<=closeFITdiff], function(c) if(all(diff(sort(unlist(splitsFound[c])))<=closeFITmiss)){sort(unlist(splitsFound[c], use.names=F))})
    #         
    #         tempSeqs <- tempList[!sapply(tempList,is.null)]
    #         newSeq   <- min(unique(unlist(tempSeqs))):max(unique(unlist(tempSeqs)))
    #         attr(tempSeqs,"Close") <- (f_l(unique(unlist(tempSeqs)))%in%NOcigar)
    #         attr(tempSeqs,"Missing") <- newSeq[!(newSeq%in%unique(unlist(tempSeqs)))]
    #       }
    #       if(attr(tempSeqs,"Close")){foundD=TRUE} 
    #       if((length(attr(tempSeqs,"Missing"))>0)&(length(attr(tempSeqs,"Missing"))<=closeFITmiss)){foundR=TRUE}
    #      
    #      if(all(sapply(newSeq,is.null))){
    #         newSeq <- vault[[1]]
    #         tempSeqs <- vault[[2]]
    #       }
    #     }
    #     
    if((foundD&foundR)){
      # At this stage... yes
      foundS <- TRUE
    } else {
      tempSeqs<-list(sFound)
      attr(tempSeqs,"Close") <- (f_l(unique(unlist(tempSeqs)))%in%NOcigar)
      attr(tempSeqs,"Missing") <- ""
    }
  } # While any
  
  # Fill in the blanks
  pagelocs <- locatePages(unique(unlist(c(tempSeqs,attr(tempSeqs,"Missing")))),norm=1)
  
  if(pagelocs$loc[length(pagelocs$loc)]==-1){
    pagelocs$loc[length(pagelocs$loc)]<-lst_page
  }
  if(pagelocs$loc[1]==-1){
    pagelocs$loc[1]<-fst_page
  }
  pagelocs$loc[pagelocs$loc==-1]<-NA
  pagelocs$loc[which(pagelocs$page%in%attr(tempSeqs,"Missing"))] <- round(approx(x=pagelocs$page,y=pagelocs$loc,xout=attr(tempSeqs,"Missing"))$y)
  
  # Create a join-ready structure
  pageSkeleton <- data.frame(location=pagelocs$loc,pdf.id=PDFid,name="",element="page",header=pagelocs$page)
  
  return(pageSkeleton)
} # getPage


# ESCI functions ----------------------------------------------------------

# Conversion formula's from Friedman(1982) and Wolf(1986)
# Also see http://www.soph.uab.edu/Statgenetics/People/MBeasley/Courses/EffectSizeConversion.pdf
f_d <- function(f,df1,df2){2*sqrt(df1*f/df2)}
f_r <- function(f,df1,df2){sqrt((df1*f)/((df1*f) + df2))}
t_d <- function(t,df){(2*t)/sqrt(df)}
t_r <- function(t,df){sqrt(t^2/(t^2 + df))}
X_d <- function(X,N,df=1){ifelse(df>1,{sqrt(2*X/N)},{2*X/sqrt(N-X)})}
X_r <- function(X,N,df=1){ifelse(df>1,{sqrt(X/(X+N))},{sqrt(X/N)})}
r_d <- function(r){sqrt((4*(r^2))/(1-r^2))}  
d_r <- function(d){sqrt(d^2/(4+d^2))}  

getESCI <- function(data, CL=.95){
  # Function to get ESCIs based on the test statistic
  require(MBESS)
  if(is.na(data$Value)){
    res <- rep(NA,times=4)
    dr  <- as.data.frame(rbind(rep(NA, times=6))) # Need N!
    names(dr) <- c("d","d.lo","d.hi","r","r.lo","r.hi")}
  else {
    if(data$Statistic == "F"){
      res <- conf.limits.ncf(F.value=data$Value, conf.level=CL, df.1=data$df1, df.2=data$df2)
      dr  <- cbind(d=f_d(data$Value, df1=data$df1, df2=data$df2), d.lo=f_d(res[[1]], df1=data$df1, df2=data$df2), d.hi=f_d(res[[3]], df1=data$df1, df2=data$df2), r=f_r(data$Value, df1=data$df1, df2=data$df2), r.lo=f_r(res[[1]], df1=data$df1, df2=data$df2), r.hi=f_r(res[[3]], df1=data$df1, df2=data$df2) )
    }
    if(data$Statistic == "t"){
      res <- conf.limits.nct(t.value=data$Value, conf.level=CL, df=data$df1)
      dr  <- cbind(d=t_d(data$Value, df=data$df1), d.lo=t_d(res[[1]], df=data$df1), d.hi=t_d(res[[3]], df=data$df1), r=t_r(data$Value, df=data$df1), r.lo=t_r(res[[1]], df=data$df1), r.hi=t_r(res[[3]], df=data$df1) )}
    if(data$Statistic == "Chi2"){
      res <- conf.limits.nc.chisq(Chi.Square=data$Value, conf.level=CL, df=data$df1)
      dr  <- as.data.frame(rbind(rep(NA, times=6))) # Need N!
      names(dr) <- c("d","d.lo","d.hi","r","r.lo","r.hi")
    }
    if(data$Statistic == "Z"){
      res <- conf.limits.nct(t.value=data$Value, conf.level=CL, df=37.62) #Using Max. ncp for t-distribution
      dr  <- as.data.frame(rbind(rep(NA, times=6))) # Need N!
      names(dr) <- c("d","d.lo","d.hi","r","r.lo","r.hi")
    }
    if(data$Statistic == "r"){
      res <- ci.R(R=data$Value, conf.level=CL, N=(data$df1+2), K=1)
      dr  <- cbind(d=r_d(data$Value), d.lo=r_d(res[[1]]), d.hi=r_d(res[[3]]), r=data$Value, r.lo=res[[1]], r.hi=res[[3]] )}
  }
  statci <- as.data.frame(cbind(levels(data$Statistic)[data$Statistic],CL,data$Value,rbind(res[c(1,3)])))
  names(statci) <- c("stat","CL","stat.v","stat.lo","stat.hi")
  return(cbind(statci,as.data.frame(dr)))
}


# PATTERN GENERATORS  -----------------------------------------------------

wordXX   <- function(x,clps="|",pre="\\b(",post=")\\b") {# ( X1|X2 ) Example: wordXX(c("yes","no","maybe"))
  paste(pre, paste(x,collapse=clps), post, sep="")
}

# wordXXp  <- function(x,clps=")|(",pre="\\b((",post="))\\b") {# Same as wordXX + escape puntuation
#   x <- gsub("(?=[[:punct:]])","\\",x,perl=T)
#   paste(pre, paste(x,collapse=clps), post, sep="")
# }

wordXnY   <- function(x,y,clps="|",pre="\\b(",post=")\\b",mid="[[:print:]]",n=1,m=1){# ( XanyY ) Example: wordXnY("yes","no")
  paste(pre, paste(x,collapse=clps),")",mid,"{",n,",",m,"}(", paste(y,collapse=clps), post, sep="")
}

wordXnYr  <- function(x,y,clps="|",pre="\\b(",post=")\\b",mid="[[:print:]]",n=1,m=1){# ( XanyY | YanyX ) Example: wordXnYr("yes","no")
  paste("(",paste(pre,paste(x,collapse=clps),")",mid,"{",n,",",m,"}(",paste(y,collapse=clps),post,sep=""),
        ")|(",paste(pre,paste(y,collapse=clps),")",mid,"{",n,",",m,"}(",paste(x,collapse=clps),post,sep=""),")", sep="")
}

wordXnYro  <- function(x,y,clps="|",pre="\\b(",post=")?\\b",mid="[[:print:]]",n=1,m=1){# ( XanyY? | Y?anyX ) Example: wordXnYr("yes","no")
  paste("(",paste(pre,paste(x,collapse=clps),")",mid,"{",n,",",m,"}(",paste(y,collapse=clps),post,sep=""),
        ")|(",paste(pre,paste(y,collapse=clps),")",mid,"{",n,",",m,"}(",paste(x,collapse=clps),post,sep=""),")", sep="")
}

wordXsYo   <- function(x,y,clps="|",pre="\\b(",post=")?\\b") {# ( XspaceY? ) Example: wordXsY("yes","no")
  paste(pre, paste(x,collapse=clps), ")[[:blank:]]*(", paste(y,collapse=clps), post, sep="")
}

wordXsYr  <- function(x,y,clps="|",pre="\\b(",post=")\\b") {# ( XspaceY | YspaceX ) Example: wordXsYr("yes","no")
  paste("(",paste(pre,paste(x,collapse=clps),")[[:blank:]]*(",paste(y,collapse=clps),")\\b",sep=""),
        ")|(",paste("\\b(",paste(y,collapse=clps),")[[:blank:]]*(",paste(x,collapse=clps),post,sep=""),")", sep="")
}

# # ( XspaceYspaceZ | YspaceZspaceX |  ZspaceXspaceY ) Example: wordXsYsZr("yes","no","maybe")
# wordXsYsZr <- function(x,y,z) {
#   paste("(",paste("\\b(",paste(x,collapse="|"),")\\s(",paste(y,collapse="|"),")\\s(",paste(z,collapse="|"), ")\\b", sep=""),
#     ")|(",paste("\\b(",paste(y,collapse="|"),")\\s(",paste(z,collapse="|"),")\\s(",paste(x,collapse="|"), ")\\b", sep=""),
#     ")|(",paste("\\b(",paste(z,collapse="|"),")\\s(",paste(x,collapse="|"),")\\s(",paste(y,collapse="|"), ")\\b", sep=""),")")
# }

wordXsYro <- function(x,y,clps="|",pre="\\b(",post=")?\\b") {# (XspaceY)? | (YspaceY)?
  paste("(",paste(pre, paste(x,collapse=clps),")[[:blank:]]*(",paste(y,collapse=clps),post,sep=""),
        ")|(",paste(pre, paste(y,collapse=clps),")[[:blank:]]*(",paste(x,collapse=clps),post,sep=""),")",sep="")
}


# LOOK FOR NEIGHBOURS (WORDS) ---------------------------------------------
# Use these tools to get an idea about which regexpatterns will give max results 

NwN   <- function(Npre=NULL, word=NULL, Npos=NULL){
  # ( Npre X Npos ) Example: txt="Love thy neigbours"  
  # matchIT(NwN(1,"thy",0),txt)
  ifelse(all(is.null(Npre),is.null(word),is.null(Npos)),
         print("NwordN: One or more arguments missing!"),
{
  #word<-gsub("[[:punct:]]+","",word)
  regstr <-paste(c("\\b", paste(rep("(\\w+[[:graph:]]*)?",Npre),collapse="\\s?"),
                   paste(c("",paste(c("(",paste(word,collapse="|"),")"),collapse=""),""),collapse="\\s"),
                   paste(rep("(\\w+[[:graph:]]*)?",Npos),collapse="\\s?"),"\\b"),collapse="")
  print(regstr) 
  return(regstr)
})
}

NwNwN   <- function(Npre=NULL, word1=NULL, Nmid=NULL, word2=NULL, Npos=NULL){
  # ( Npre X Nmid Y Npos ) Example:  txt="Didn't you know? Love thy neigbours' neigbours, too!"  
  # matchIT(NwNwN(0,c("Love","hate","you"),2,c("neigbours","friends","family"),0),gsub("[[:punct:]]+","",txt))  
  # matchIT(NwNwN(1,"Love",3,"too",0),gsub("[[:punct:]]+","",txt))                                   
  
  ifelse(all(is.null(Npre),is.null(word1),is.null(Nmid),is.null(word2),is.null(Npos)),
         print("NwNwN: One or more arguments missing!"),
{
  #   word1<-gsub("[[:punct:]]+","",word1)
  #   word2<-gsub("[[:punct:]]+","",word2)
  regstr <-paste(c("\\b", paste(rep("(\\w+[[:graph:]]*)?",Npre),collapse="\\s?"),
                   paste(c("",paste(c("(",paste(word1,collapse="|"),")"),collapse=""),""),collapse="\\s*"),
                   paste(rep("(\\w+[[:graph:]]*)?",Nmid),collapse="\\s?"),
                   paste(c("",paste( c("(",paste(word2,collapse="|"),")"),collapse=""),""),collapse="\\s*"),
                   paste(rep("(\\w+[[:graph:]]*)?",Npos),collapse="\\s?"),"\\b"),collapse="")
  print(regstr) 
  return(regstr)
})
}


NwNwNr  <- function(Npre=NULL, word1=NULL, Nmid=NULL, word2=NULL, Npos=NULL){
  # (Npre X Nmid Y Npos) | (Npre Y Nmid X Npos) Example:    txt="Didn't you know? Your neigbours love you!"  
  # matchIT(NwNwNr(1,"know",2,"neigbours",1),gsub("[[:punct:]]+","",txt))  
  # matchIT(NwNwNr(0,"Your",2,"love",0),gsub("[[:punct:]]+","",txt))
  
  ifelse(all(is.null(Npre),is.null(word1),is.null(Nmid),is.null(word2),is.null(Npos)),
         print("NwNwN: One or more arguments missing!"),
{ #word1<-gsub("[[:punct:]]+","",word1)
  #   word2<-gsub("[[:punct:]]+","",word2)
  regstr1 <-paste(c("\\b", paste(rep("(\\w+[[:graph:]]*)?",Npre),collapse="\\s?"),
                    paste(c("",paste(c("(",paste(word1,collapse="|"),")"),collapse=""),""),collapse="\\s*"),
                    paste(rep("(\\w+[[:graph:]]*)?",Nmid),collapse="\\s?"),
                    paste(c("",paste( c("(",paste(word2,collapse="|"),")"),collapse=""),""),collapse="\\s*"),
                    paste(rep("(\\w+[[:graph:]]*)?",Npos),collapse="\\s?"),"\\b"),collapse="")
  regstr2 <-paste(c("\\b", paste(rep("(\\w+[[:graph:]]*)?",Npre),collapse="\\s?"),
                    paste(c("",paste(c("(",paste(word2,collapse="|"),")"),collapse=""),""),collapse="\\s*"),
                    paste(rep("(\\w+[[:graph:]]*)?",Nmid),collapse="\\s?"),
                    paste(c("",paste( c("(",paste(word1,collapse="|"),")"),collapse=""),""),collapse="\\s*"),
                    paste(rep("(\\w+[[:graph:]]*)?",Npos),collapse="\\s?"),"\\b"),collapse="")
  regstr<-paste("(",regstr1,")*|(",regstr2,")*",sep="")
  print(regstr) 
  return(regstr)
})
}