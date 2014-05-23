# sciCure (http://fredhasselman.github.io/scicuRe/)
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
# 
# source_https <- function(url, ...) {
#   require(RCurl)
#   # parse and evaluate each .R script
#   sapply(c(url, ...), function(u) {
#     eval(parse(text = getURL(u, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), envir = .GlobalEnv)
#   })
# }
# # Source one of my tool collections
# source_https("https://raw2.github.com/FredHasselman/toolboxR/master/BTBTB.R")


# Packages in the list argument need will be installed if necessary and loaded
inIT <- function(need){
  ip <- .packages(all.available=T)
  if(any((need %in% ip)==F)){install.packages(need[!(need %in% ip)])}
  ok <- sapply(1:length(need),function(p) require(need[[p]],character.only=T))
}

# Packages in the list argument loose will be unloaded, if necessary (del=T) uninstalled
unIT <- function(loose,del=F){
  dp <- .packages()
  if(any(loose %in% dp)){
    for(looseLib in loose[(loose %in% dp)]){detach(paste0("package:",looseLib),unload=T,character.only=T)}
  }
  rm(dp)
  if(del==T){
    dp <- .packages(all.available=T) 
    if(any(loose %in% dp)){remove.packages(loose[(loose %in% dp)])}
  }
}


# META functions ----------------------------------------------------------
#IT <- TXTfiles

structureIT <- function(IT,chunk=10,saveOnly=FALSE){# It's a bit of a monster, but gets the job done
  require(plyr)
  # Load predefined search patterns
  structfind <- loadPatterns()
  element    <- c("header","section","table","figure")
  citeSearch <- FALSE
  pageSearch <- FALSE
  
  fnames     <- attributes(IT)$fullName
  pageSeq    <- attributes(IT)$pageSeq
  citation   <- attributes(IT)$citation
  
  if(any(is.null(citation))){citeSearch <- TRUE}
  if(any(is.null(pageSeq))){pageSearch <- TRUE}
  
  # Go chunky?
  if(length(fnames)>chunk){
    chunks <- unique(c(seq(chunk,length(fnames),by=chunk),length(fnames)))
    cat("\n\nThere are",length(fnames),"files to be analysed: Will save a temporary file after each:",chunks,"files...\n")
  } else {
    chunks <- length(fnames)
  }
  chunks <- c(0, chunks)
  
  for(chnk in 1:(length(chunks)-1)){
    ITbit     <- IT[(seq_along(IT)<=chunks[chnk+1])&(seq_along(IT)>chunks[chnk])]
    if(any(lapply(ITbit,length)>1)){
      remID <- which(lapply(ITbit,length)>1)
      CHNK <- seq((chunks[chnk]+1),chunks[chnk+1],by=1)[-remID]
    } else {
      CHNK <- seq((chunks[chnk]+1),chunks[chnk+1],by=1)
    }
    ITbit     <- IT[CHNK]
    CHNKnames <- fnames[CHNK] 
    CHNKcite  <- citation[CHNK] 
    CHNKpages <- pageSeq[CHNK]
    ITbitINFO <- list()
    
    cat("\nStarting chunk,",chnk,"...")  
    
    if(citeSearch){
      cat("\n\nTrying to find the publication on PubMed and checking whether file was dowloaded from JSTOR...\n\n")
      for(c in seq_along(CHNKcite)){
        ifelse(is.null(CHNKcite[c]),{
          ITbitINFO[[c]] <- getPUBREF(ITbit[[c]])
        },{
          ITbitINFO[[c]] <- list(Source=basename(CHNKnames[info]),pageSeq=CHNKpages[[info]],reference=CHNKcite[info],kind="User Provided")
        })
      }
    } else {
      ITbitINFO <- llply(seq_along(ITbit),function(info) return(list(Source=basename(CHNKnames[info]),pageSeq=CHNKpages[[info]],reference=CHNKcite[info],kind="User Provided")))
    }
    
    # Build a skeleton for each article
    # 1. Get location of article structure indices
    #     if(ITbitINFO$kind=="JSTOR"){
    #       firstP <- locIT("(http://(www.)?jstor.org[^/])[[:print:]]+\\n[\\n\\f]*",ITbit)[[1]]
    #       if(firstP!=-1){ITbit <- substr(ITbit,(firstP+attr(firstP,"match.length"))[[1]],nchar(ITbit))}
    #     }
    # 2. Find header info
    #     skeleton$header <- NA
    #     skeleton$header[skeleton$location!=-1] <- unlist(llply(seq_along(element), function(e) matchIT(structfind[["skeletonPatterns"]][[element[e]]],ITbit,ic=T,cln=F)))
    
    skeleton  <- getSkelet(element,ITbit)
    skeleton  <- subset(skeleton, location!=-1)
    
    # 3. Find/infer and add pagenumbers to the skeleton
    cat("\nFinding page numbers in file:\n")
    pagenumbers <- getPageSkelet(ITbit,ITbitINFO)
    pagenumbers <- subset(pagenumbers, location!=-1)
    
    cat("\nCollecting results...")
    skeleton        <- rbind(skeleton, pagenumbers)
    # 4. Add names, tidy up, sort results within each article
    #     skeleton$name   <- factor(skeleton$pdf.id,levels=seq_along(CHNKnames),labels=names(unlist(ITbit,recursive=F)))
    #     skeleton        <- subset(skeleton, location!=-1)
    #    skeleton$pdf.id   <- as.numeric(as.vector(skeleton$pdf.id))
    skeleton$location <- as.numeric(as.vector(skeleton$location))
    skeleton          <- arrange(skeleton,file.name,location) #skeleton[ do.call(order,skeleton), ] 
    
    cat("\nSaving chunk,",chnk,"\n")
    
    saveRDS(skeleton,file=paste("scicuRe_chunk",chnk,".rds",sep=""))
    saveRDS(ITbit,file=paste("scicuRe_text",chnk,".rds",sep=""))
    rm(skeleton,ITbit,pagenumbers,CHNKnames,CHNKcite,CHNKpages)
  }
  ifelse(saveOnly,{return(list("Skeletons and text for each file were saved as chunks: scicuRe_chunk##.rds"))},{
    skeleton <- ldply(seq_along(chunks[-1]),function(chnk) readRDS(paste("scicuRe_chunk",chnk,".rds",sep="")))
    ftext    <- unlist(llply(seq_along(chunks[-1]),function(chnk) readRDS(paste("scicuRe_text",chnk,".rds",sep=""))))
    return(list(skeleton=skeleton,TXT=ftext))
  })
}

getSkelet <- function(element,ITbit){
  skelet     <- list(-1)
  structfind <- loadPatterns()
  for(e in seq_along(element)){
    cat("\nFinding element\t<",element[[e]],">")
    locs <- matchIT(structfind[["skeletonPatterns"]][[element[e]]],ITbit)  
    if(any(sapply(locs,length)==0)){
      noloc <- which(sapply(locs,length)==0)
      locs[noloc] <- "Not found"
      for(loc in noloc){attr(locs[[loc]],"location") <- -1}
    }
    #(([[:punct:]]([[:blank:]]*[^[:print:]])+)+|[^[:print:]]{1,3}
    #matchIT(paste("([^[:print:]]{1,3}|\\b)+",structfind[[7]]["real_pat"],"([[:punct:]]*[[:blank:]]*[^[:print:]]){1,1}",sep=""),eLOC$value,cln=F)
    eLOC <- ldply(locs,function(l) cbind(element=element[[e]],location=attr(l,"location"),melt(l),stringsAsFactors=F))
    id   <- sort(c(which(eLOC$value==""))) #,grep(paste("([^[:print:]]{1,4}|\\b)+",structfind[[7]]["real_pat"],"(([[:punct:]]*[[:blank:]]*[^[:print:]])+|",structfind[[7]]["real_pat"],")+(\\n?\\d+(\\.)?)*",sep=""),eLOC$value,perl=TRUE)))
    if(length(id)>0){eLOC<-eLOC[-id, ]}
    id   <- sort(c(which(eLOC$value==0),grep(paste("[^[:print:]]{1,4}[\\d]+",sep=""),eLOC$value,perl=TRUE)))
    if(length(id)>0){eLOC<-eLOC[-id, ]}
    names(eLOC) <- c("file.name","element","location","header")
    skelet[[e]] <- eLOC
  }
  return(ldply(skelet))
}

getPageSkelet <- function(ITbit,ITbitINFO){
  pagePat <- loadPatterns()[[8]]["page_pat"]
  pat     <- loadPatterns()[["skeletonPatterns"]][["page"]]
  # Seperate the look ahead/behind conditions
  prepost <- strsplit(pat,pagePat,fixed=T)[[1]]
  pskelet  <- list(0)
  
  for(j in seq_along(ITbit)){
    cat(names(ITbit)[j],"\n")
    if(ITbitINFO[[j]]$kind!="Not Found"){
      loc   <- matchIT(wordXX(pre=paste(prepost[1],"(",sep=""),post=paste(")",prepost[2],sep=""),x=ITbitINFO[[j]]$pageSeq),ITbit[[j]])[[1]]
      found <- evalSeq(ITbitINFO[[j]]$pageSeq,loc)
      closeFITdiff <- 1
      tmp    <- list(0)
      tmpLoc <- loc
      while(found!="Exact"){
        if(found=="Exact+"){
          tmp   <- splitSeq(tmpLoc)
          tmp   <- rangeCheck(tmp,length(ITbitINFO[[j]]$pageSeq),closeFITdiff)
          found <- evalSeq(ITbitINFO[[j]]$pageSeq,tmp) 
        }
        if(found=="PartSeq+"){
          # If it was a partial match, longer than the target sequence try a less greedy search pattern
          tmpLoc <- gsub("(\\D)","",matchIT(wordXnYr(x="[^[:print:]]{0,1}",y="[^[:print:]]{1,2}",mid=wordXX(x=ITbitINFO[[j]]$pageSeq)),ITbit[[j]])[[1]])
          tmp   <- splitSeq(tmpLoc)
          idT <- numeric(0)
          idL <- numeric(0)
          for(seq in seq_along(tmp)){
            if(all(ITbitINFO[[j]]$pageSeq%in%tmp[[seq]])){
              idT <- cbind(idT,seq)
            }
          }
          
          if(length(tmp)>10){
            idL <- which(llply(tmp,length)>1) 
          }
          
          if(length(idL)>0|length(idT)>0){
            id <- as.vector(unique(cbind(idL,idT)))
            tmptmp <- tmp[id]
            attributes(tmptmp)$location <- attributes(tmp)$location[id]
            tmp <- tmptmp
          }
          tmp   <- rangeCheck(tmp,length(ITbitINFO[[j]]$pageSeq),closeFITdiff)
          found <- evalSeq(ITbitINFO[[j]]$pageSeq,tmp) 
        }
        if(found=="PartSeq-"){
          # If it was a partial match, shorter than the target sequence try a more greedy search pattern
          tmpLoc   <- gsub("(\\D)","",matchIT(wordXnYr(x="[^[:print:]]{0,3}",y="[^[:print:]]{1,3}",mid=wordXX(x=ITbitINFO[[j]]$pageSeq)),ITbit[[j]])[[1]])
          tmp   <- splitSeq(tmpLoc)
          idT <- numeric(0)
          idL <- numeric(0)
          for(seq in seq_along(tmp)){
            if(all(ITbitINFO[[j]]$pageSeq%in%tmp[[seq]])){
              idT <- cbind(idT,seq)
            }
          }
          
          if(length(tmp)>10){
            idL <- which(llply(tmp,length)>1) 
          }
          
          if(length(idL)>0|length(idT)>0){
            id <- as.vector(unique(cbind(idL,idT)))
            tmptmp <- tmp[id]
            attributes(tmptmp)$location <- attributes(tmp)$location[id]
            tmp <- tmptmp
          }
          tmp   <- rangeCheck(tmp,length(ITbitINFO[[j]]$pageSeq),closeFITdiff)
          found <- evalSeq(ITbitINFO[[j]]$pageSeq,tmp) 
        }
        
        if(found=="None"){
          # If it was a partial match, shorter than the target sequence try a more greedy search pattern
          tmpLoc   <- gsub("(\\D)","",matchIT(wordXnYr(x="[^[:print:]]{0,3}",y=paste0("[^[:print:]]{1,3}(?=",wordXX(x=ITbitINFO[[j]]$pageSeq),"?[^[:print:]]{0,3})?"),mid=wordXX(x=ITbitINFO[[j]]$pageSeq)),ITbit[[j]])[[1]])
          tmp   <- splitSeq(tmpLoc)
          
          idT <- numeric(0)
          idL <- numeric(0)
          for(seq in seq_along(tmp)){
            if(all(ITbitINFO[[j]]$pageSeq%in%tmp[[seq]])){
              idT <- cbind(idT,seq)
            }
          }
          
          if(length(tmp)>10){
            idL <- which(llply(tmp,length)>1) 
          }
          
          if(length(idL)>0|length(idT)>0){
            id <- as.vector(unique(cbind(idL,idT)))
            tmptmp <- tmp[id]
            attributes(tmptmp)$location <- attributes(tmp)$location[id]
            #             tmptmp <- tmp[id]
            #             for(d in seq_along(id)){attr(tmptmp[[d]],"location") <- attr(tmp,"location")[[id[d]]]}
            tmp <- tmptmp
          }
          tmp   <- rangeCheck(tmp,length(ITbitINFO[[j]]$pageSeq),closeFITdiff)
          found <- evalSeq(ITbitINFO[[j]]$pageSeq,tmp) 
        }
        
        closeFITdiff <- closeFITdiff+1
        if(closeFITdiff>10){
          # Break out of the while loop
          found <- "Exact"
          break
        }
      }
      # Check if found is really Exact, or if we broke out of the loop
      found <- evalSeq(ITbitINFO[[j]]$pageSeq,tmp) 
      if(found=="Exact"){
        pskelet[[j]] <- data.frame(file.name=names(ITbit)[j],element="page",location=attr(tmp,"location"),header=as.character(tmp)) 
      } else {
        pskelet[[j]] <- data.frame(file.name=names(ITbit)[j],element="page",location=-1,header="000") #guessSeq(ITbit[[j]],ITbitINFO[[j]]$pageSeq,loc)
      }
    } else {
      pskelet[[j]] <- data.frame(file.name=names(ITbit)[j],element="page",location=-1,header="000")
    }
  }
  return(ldply(pskelet))
}


placeIT <- function(IT,skeleton){
  require(reshape2)
  out <- vector("list")
  cnt=0
  for(p in seq_along(IT)){
    if(length(IT[[p]])>0){
      tmp1 <- melt(IT[[p]])
      names(tmp1)[1]<- "header"
      tmp1$location <- attr(IT[[p]],"location")
      tmp1$element  <- "metastat"
      skel.tmp      <- rbind(subset(skeleton, file.id==skeleton$file.name[p], select=c(element,location,header)))
      tmp2          <- rbind(skel.tmp,tmp1)
      tmp2$location <- as.numeric(as.vector(tmp2$location))
      tmp2          <- tmp2[ do.call(order,tmp2), ]
      metaID        <- which(tmp2$element=="metastat")
      
      for(i in metaID){
        idAbove   <- which(tmp2$location<tmp2$location[i])
        idHeader  <- idAbove[which(tmp2$element[idAbove]=="header")]
        # Prefer Headers with numbers
        if(length(idHeader) > 0){
          if(grepl("\\d+\\w*",tmp2$header[max(idHeader)])){
            idHeader<-max(idHeader)
          } else {
            if(any(grepl("\\d+\\w*",tmp2$header[idHeader]))){
              idHeader<-max(idHeader[grepl("\\d+\\w*",tmp2$header[idHeader])])
            } else {
              idHeader<-max(idHeader)
            }  
          }
        } else {idHeader<-nrow(tmp2)+1}
        idSection <- idAbove[which(tmp2$element[idAbove]=="section")]
        if(length(idSection) > 0){idSection<-max(idSection)} else {idSection<-nrow(tmp2)+1}
        idTable   <- idAbove[which(tmp2$element[idAbove]=="table")]
        if(length(idTable) > 0){idTable<-max(idTable)} else {idTable<-nrow(tmp2)+1}
        idBelow  <- which(tmp2$location>tmp2$location[i])
        idPage   <- idBelow[which(tmp2$element[idBelow]=="page")]
        if(length(idPage) > 0){idPage<-min(idPage)} else {idPage<-nrow(tmp2)+1}
        
        tmp <- rbind(tmp2,cbind(location=-1,element="none",header="NOT FOUND"))
        cnt=cnt+1
        out[[cnt]] <- data.frame(file.id=which(levels(skeleton$file.id)==skeleton$file.name[p]),location=tmp$location[i],stat=tmp$header[i],page=tmp$header[idPage],section=tmp$header[idSection],header=tmp$header[idHeader],table=tmp$header[idTable],name=names(IT)[p],row.names=paste(p,cnt,sep="."),stringsAsFactors = FALSE)
        rm(idAbove,idHeader,idSection,idTable,idBelow,idPage,tmp)
      }
    } else {
      cnt=cnt+1
      out[[cnt]] <- data.frame(file.id=which(levels(skelet$file.id)==skelet$file.name[p]),location=-1,stat=paste(p,sep=""),page="",section="",header="",table="",name=names(IT)[p],row.names=paste(p,cnt,sep="."),stringsAsFactors = FALSE)
    }
  }
  return(ldply(out))
}

# Extraction functions ----------------------------------------------------
getTXTtxt <- function(TXTnames,del=TRUE,enc=NULL){
  require(plyr) 
  txt <- vector("list")
  cnt <- 0
  id  <- NULL
  for(i in seq_along(TXTnames)){
    if(is.null(enc)){
      if(file.exists(Sys.which("file"))){
        report <- system(paste("file ",TXTnames[i]),intern=T)
        str(report)
        enc    <- iconvlist()[sapply(iconvlist(),function(enc) grepl(enc,report))]
        if(is.null(enc)){
          warning(paste("Unknown text encoding in:",TXTnames[i]))
        }
      } else {
        cat("The UNIX system command 'file' was not found... (need it to get the encoding of the text file)\nPlease provide a valid encoding as a character vector in the argument 'enc'\n")
      }
    }
    tmptxt <- try.CATCH(readLines(file(TXTnames[i],encoding=enc),n=-1,warn=F,encoding="UTF-8"))
    closeAllConnections()
    if(length(tmptxt$warning)!=0){
      if(del){
        cnt         <- cnt+1
        txt[cnt]  <- paste(tmptxt$value,collapse="\n")
        id          <- c(id,i)
      } else {
        txt[i]  <- paste(tmptxt$value,collapse="\n")
        id        <- c(id,i)
      }
    }
    rm(tmptxt)
  }
  txt <- lapply(seq_along(id), function(t) structure(txt[t], names = gsub(".txt","",basename(TXTnames[id[t]])), fullName = TXTnames[id[t]]))
  return(txt)
}


getPDFtxt <- function(PDFnames,del=TRUE){
  require(plyr) 
  txt <- vector("list")
  cnt <- 0
  id  <- NULL
  for(i in seq_along(PDFnames)){
    tmptxt <- getPDFto(PDFnames[i],del=del)
    if(del){
      if(!(length(tmptxt)==0)){
        cnt         <- cnt+1
        txt[cnt]    <- tmptxt
        id          <- c(id,i)
      }
    } else {
      txt[i]  <- tmptxt
      id      <- c(id,i)
    }
    rm(tmptxt)
  }
  txt <- lapply(seq_along(id), function(t) structure(txt[t], names = gsub(".pdf","",basename(PDFnames[id[t]])), fullName = PDFnames[id[t]]))
  return(txt)
}

getPDFto <- function(PDFname,del=TRUE){
  if(file.exists(Sys.which("pdftotext"))){ 
    p2t <- try.CATCH(system(paste('pdftotext -q -enc "UTF-8" "',PDFname,'"',sep=""),intern=TRUE))
    if(length(p2t$warning)>0){
      txtfile <- list(value=c("Exit1"),warning=c("SUGGESTION: Change security settings in PDF properties."))
    } else {
      fileName <- gsub("\\.pdf$","\\.txt",PDFname)
      txtfile  <- try.CATCH(readChar(fileName, file.info(fileName)$size))
      if(length(txtfile$warning)>0){
        txtfile <- list(value=c("Exit2"),warning=c("SUGGESTION: Check the file for (encoding) anomalies."))
      }
    }
  }
  txtfile <- checkPDFtoTXT(txtfile,PDFname)
  if(del&txtfile$warning){
    txtfile <- NULL
  } else {
    txtfile        <- txtfile$value
    #names(txtfile) <- gsub(".pdf","",basename(PDFname))
  }
  return(txtfile)
}

checkPDFtoTXT <- function(txtfile,PDFname){
  if(!(grepl("Exit\\d{1}",txtfile$value)[[1]])){
    if(nchar(txtfile$value)>0){
      feeds <- locIT("(\\f)+",txtfile$value)[[1]]
      if((nchar(txtfile$value)/sum(attr(feeds,"match.length")))<100){
        txtfile$value   <- "Exit3"
        txtfile$warning <- "SUGGESTION: Apply text recognition (OCR) to the PDF file"
      }
    } else {
      txtfile$value   <- "Exit4"
      txtfile$warning <- "SUGGESTION: No content in imported text file, check if file exists"
    }
  }
  switch(txtfile$value,
         Exit1 = cat("\n\nThe file [",basename(PDFname),"] could not be converted by 'pdftotext', common cause are PDF security settings\n-> Change in PDF properties\n\n"),
         Exit2 = cat("\n\nThe text file converted from [",basename(PDFname),"] could not be read by 'readChar()' \n-> Check the file for (encoding) anomalies.\nn"),
         Exit3 = cat("\n\nThe text file converted from [",basename(PDFname),"] contains mostly non-printing characters.\n-> PDF content is likely image based. Apply text recognition software (e.g., free online OCR: http://www.newocr.com)\nn"),
         Exit4 = cat("\n\nThe text file converted from [",basename(PDFname),"] is an empty character vector.\n-> PDF content could based on 1 image, or the imported text file has conversion errors or does not exist.\nn")
  )
  if(grepl("Exit\\d{1}",txtfile$value)[[1]]){
    txtfile$value   <- txtfile$warning
    txtfile$warning <- TRUE
  } else {
    txtfile$value   <- enc2utf8(txtfile$value)
    txtfile$warning <- FALSE
  }
  return(txtfile)
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
      stop("Could not find PDF file...")
    }
  } else {
    stop("Could not find 'pdfinfo'...")
  }
  return(as.numeric(sub("\\D+","",x[grep("(?<=Pages\\:)(\\s*)(\\d+)",x,perl=T)],perl=T)))
}


getREFS <- function(IT){
  #Find REFlist
  reflistLOC <- locIT(loadPatterns(c(rep(FALSE,times=5),TRUE))[["bibPatterns"]][["reflist"]],IT[[1]])[[1]]
  if((reflistLOC>0)&(min(reflistLOC/nchar(IT))>.50)){
    reflistLOC<-min(reflistLOC)
  } else {
    reflistLOC<-round(nchar(IT)*.95)
  }
  ITnoref <- substr(IT[[1]],1,(reflistLOC-1))
  ITrefs  <- substr(IT[[1]],reflistLOC,nchar(IT))
  matchIT(loadPatterns(c(rep(FALSE,times=5),TRUE))[["bibPatterns"]][["reference"]],ITrefs,cln=F,ic=FALSE)
  matchIT("(?<authors>([[:upper:]]{1,2}[[:lower:]]+[\\,])+([[:upper:]]?|[[:punct:]]|[[:blank:]])+){1,6}(?<year>([(](\\d{4}\\w{0,1})|((in\\s*press)|(n\\.d\\.)|(submitted))[)])+[[:blank:]]*[[:punct:]]*){1}",ITrefs,cln=F,ic=F)
  return(list(locREFS=reflistLOC,ITnoREFS=ITnoref,REFSnoIT=ITrefs))
}

getQUERY <- function(IT,nw=150){
  stopwords <- readRDS("stopwords_en.rds")
  stopwords <- c(stopwords,"university","institute","hospital","school","centre","center","medical","lab","laboratory","college")
  
  # Find author
  #   matchIT(paste("(((corresponding|first)*(\\s)*author)|contact|correspondence)*",loadPatterns(c(rep(FALSE,times=5),TRUE))[["bibPatterns"]][["authors"]],substr(IT,1,round(nchar(IT)*.01)))
  
  # Find doi
  # matchIT(paste("(doi|digital\\sobject\\sidentifier)+([[:punct:]]|[[:blank:]])*([[:alnum:]]|[[:punct:]])+"),substr(IT,1,round(nchar(IT)*.02)))
  
  # Find abstract text
  breaks    <- gregexpr("[^[:print:]]*",substr(IT,1,round(nchar(IT)*.02)))[[1]]
  start     <- breaks[attr(breaks,"match.length")>0]
  if(length(start)>1){
    start <- start[length(start)]   
    start <- start + attr(breaks,"match.length")[which(breaks==start)]
  } else {
    start <- start[1]
  }
  keys     <- locIT(loadPatterns(c(rep(FALSE,times=5),TRUE))[["bibPatterns"]][["keywords"]],IT)[[1]]
  if((keys!=-1)&(min(keys/nchar(IT))<.05)){
    end    <- min(keys)-1
  } else {
    end      <- round(nchar(IT)*.02)
  }
  
  abstract <- substr(IT,start,end)
  breaks   <- gregexpr("[^[:print:]]*",abstract)[[1]]
  tempend  <- breaks[attr(breaks,"match.length")>0]
  while(length(tempend)>1){
    if(tempend[2]<round(nchar(abstract)*.5)){
      start <- tempend[2]
    }
    if(tempend[length(tempend)-1]>round(nchar(abstract)*.5)){
      end <- tempend[length(tempend)-1]
    }
    abstract <- substr(abstract,start,end)
    breaks   <- gregexpr("[^[:print:]]*",abstract)[[1]]
    tempend  <- breaks[attr(breaks,"match.length")>0]
  }
  if(length(tempend)==1){
    if(tempend[1]<round(nchar(abstract)*.5)){
      start <- tempend[1] 
      start <- start + attr(breaks,"match.length")[which(breaks==start)]
    } else {
      end  <- tempend[1]-1
    }
    abstract <- substr(abstract,start,end)
  }
  abstract <- tolower(abstract)
  terms    <- matchIT("\\b\\w+\\b",abstract)[[1]]
  if(length(terms)<100){
    abstract <-  tolower(substr(IT,attr(terms,"location")[5],round(nchar(IT)*.02)))
  }
  # Remove unwanted characters  
  realPat  <- loadPatterns(rep(FALSE,times=6))[[1]]["real_pat"]
  abstract <- gsub(realPat," ",abstract)
  abstract <- gsub("[[:punct:]]"," ",abstract)
  abstract <- gsub("(\\s+)"," ",abstract)
  abstract <- gsub("^(\\s+)|(\\s+)$","",abstract)
  terms    <- matchIT("\\b\\w+\\b",abstract)[[1]]
  terms    <- unique(terms[!terms%in%stopwords])
  return(terms)
}

#   TMP <- abstract
#   while(any(grepl("(\\n|\\f|\\r)+",TMP))){
#     loc <- locIT("(\\n|\\f|\\r)+",TMP)
#     pos <- max(loc[[1]])+attr(loc[[1]],"match.length")[which.max(loc[[1]])]
#     if(pos < round(nchar(TMP)*.8)){
#       TMP <- substr(TMP,pos,nchar(TMP))
#     } else {
#       TMP <- substr(TMP,start,(nchar(TMP)-pos))
#     }
#   }
#   abstract <- TMP
#  abstract <- substr(IT,start,attr(terms,"location")[nw])
#    keys     <- attr(terms,"location")[nw]
#   if((keys!=-1)&(min(keys/nchar(IT))<.1)){
#     keys <- min(keys)
#     abstract <-  substr(IT,start,(keys-1))
#     terms    <- matchIT("\\b\\w+\\b",abstract)[[1]]
#     if(length(terms)>nw){
#       start <- attr(terms,"location")[length(terms)-nw]
#     }
# #     else {
# #       start <- 10
# #     }
#   } else {
#     #terms <- matchIT("\\b\\w+\\b",substr(IT,start,5000))[[1]]
#    
#   }
#   abstract <- substr(IT,start,(keys-1)
#   abstract <- gsub("(\\s+[[:print:]]{1,15}\\s*)$","",abstract)

checkJSTOR <- function(IT){
  ITbit    <- substr(IT[[1]],1,round(nchar(IT)*.11))
  if(grepl("JSTOR",ITbit)){
    url <- matchIT("(stable url[\\:][\\n\\f\\s\\t\\r])(http)(\\S)+",ITbit,cln=F)[[1]]
    txt <- substr(ITbit[[1]],1,(attr(url,"location")-1)[1])
    pages <- matchIT("(?<=(p[p][\\.][\\s]))[[:print:]]+(?=[\\.][\n])",txt,cln=F)[[1]]
    signPat <- sub("\\?","",loadPatterns(rep(FALSE,times=6))[[1]]["sign_pat"])
    pageMnMx<- strsplit(pages,signPat)[[1]]
    if(nchar(pageMnMx[1])!=nchar(pageMnMx[2])){
      pages <- paste(pageMnMx[1],"-",substr(pageMnMx[1],start=1,stop=nchar(pageMnMx[1])-nchar(pageMnMx[2])),pageMnMx[2],sep="")
    }
    pageSeq <- try.CATCH(seq(as.numeric(strsplit(pages,signPat)[[1]][1]),as.numeric(strsplit(pages,signPat)[[1]][2])))
    if(grepl("Error",strsplit(as.character(pageSeq$value),"\\s"))[[1]]){pageSeq <- pages} else {pageSeq <- pageSeq$value}
    JSTORref <- return(list(Source=gsub("[[:print:]]+(?=(http))","",url[1],perl=T),pageSeq=pageSeq,reference=paste(gsub("[\\n\\f]*","",txt,perl=T),url),kind="JSTOR"))
  } else {
    JSTORref <- NULL
  }
  return(JSTORref)
}

getPUBREF <- function(IT){
  require(RISmed)
  
  if(is.null(attributes(IT))){
    IT <- IT[[1]]
  }
  
  res <- readLines(paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=",paste(sample(abstract,5),collapse="+"),"refers&retmax=10&tool=scicuRe",sep=""))
  #   matchIT("(?<=\\<Id\\>)[[:print:]]*(?=\\<\\/Id\\>)",res,cln=F)
  #   content(res)
  #   matchIT("(?<=<IdList>)[[:print:]]*(?=</IdList>)",as.character(
  #     
  #     qr <- matchIT("(?<=\\t\\<CorrectedQuery\\>)[[:print:]]*(?=\\<\\/CorrectedQuery\\>)",readLines(paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/espell.fcgi?db=pubmed&term=",paste(abstract,collapse="+"),"&refers&retmax=10&tool=scicuRe",sep="")),cln=F,ic=T)
  #     qr <-unlist(qr)
  #     res <- readLines(paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=",paste(qr,collapse="+OR+"),"&refers&retmax=10&tool=scicuRe",sep=""))
  #     ids <-  unlist(matchIT("(?<=\\<Id\\>)[[:print:]]*(?=\\<\\/Id\\>)",res,cln=F))
  #  
  #     getIDs <- readLines(paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=",paste(ids,collapse=","),"&retmode=text&rettype=abstract",sep=""))
  #     
  #     qr <- getNodeSet(xmlTreeParse(res,useInternalNodes=T),"/eSearchResult/IdList/Id") # function(id) gsub("<(\\/)?Id>",aid),"")))
  #     gsub("<(\\/)?Id>","",qr[[1]])
  #     
  #     xmlElementsByTagName(xmlRoot(qr)[[1]], "IdList")
  #     ,"IdList")
  #     ))
  
  EUtilsQuery(paste(sample(abstract,5),collapse=" "))
  EUtilsGet
  trace(EUtilsSummary,edit=T)
  # First check JSTOR and return if something was found
  JSTORcite <- checkJSTOR(IT)
  if(is.null(JSTORcite)){
    JSTORfound <- FALSE
  } else {
    cat("Found citation for file [",names(IT),"]\n--",JSTORcite$reference,"\n\n")
    return(JSTORcite)
  }
  
  if(JSTORfound==FALSE){  
    cnt <- 0
    pp  <- seq(1,round(nchar(IT)*.3),by=round(nchar(IT)*.02))
    tryREFS <- FALSE
    # Try to get citation from PubMed based on text that is likely from the Abstract
    abstract <- getABSTRACT(IT[[1]])
    res      <- EUtilsSummary(paste(sample(abstract,5),collapse=" "),type="esearch",db="pubmed",encoding="UTF-8")
    while(res@count!=1){
      # Try different term lengths
      for(nwrd in seq(5,length(abstract)-1,by=5)){
        res      <- EUtilsSummary(paste(sample(abstract,nwrd),collapse=" "),type="esearch",db="pubmed",encoding="UTF-8")
      }
      # Try different pages
      cnt      <- cnt+1
      ITbit    <- substr(IT[[1]],pp[cnt],nchar(IT))
      abstract <- getABSTRACT(ITbit)
      #res      <- EUtilsSummary(paste(abstract,sep=""),type="esearch",db="pubmed",encoding="UTF-8")
      
      if(cnt==length(pp)){
        tryREFS=TRUE
        break
      }
    }
    
    if(tryREFS==TRUE){
      refs <- getREFS(IT)["REFSnoIT"]
      refs <- substr(refs,1,round(nchar(refs)/2))
      res  <- EUtilsSummary(paste(abstract,sep=""),type="esearch",db="pubmed",encoding="UTF-8")
    }
    
    if((res@count==0)&(JSTORfound==FALSE)){
      cat("\nNo references found... checked PubMed (online) and JSTOR (in PDF)\n\n")
      return(list(Source=NA,pageSeq=NA,reference=NA,kind="Not Found"))
    }
    if(res@count>1){cat("NOTE: Found more than 1 article on PubMed, using the first hit\n")}
    
    # Get the pub
    pub <- EUtilsGet(res@id[1])
    # Get the info
    authors <- Author(pub)[[1]]
    authorString <- list(character(0))
    for(a in authors[,"order"]){
      authorString[[a]] <- paste(authors[a,"LastName"], gsub("(\\w)","\\U\\1\\.",paste(authors[a,"Initials"]),perl=T),"",sep=", ")
      if(a>=6){
        authorString[[a]] <- paste(authorString[[a]],", et al. (",Year(pub),")",sep="")
        break
      } else {
        if(a==(max(authors[,"order"])-1)){authorString[[a]] <- paste(authorString[[a]],"& ",sep="")}
        if(a==max(authors[,"order"])){authorString[[a]] <- gsub("(\\s\\w\\.)(\\,\\s)$",paste("\\1 (",Year(pub),"). ",sep=""),paste(authorString[[a]],sep=""))}
      }
    }
    signPat <- sub("\\?","",loadPatterns(rep(FALSE,times=6))[[1]]["sign_pat"])
    pageMnMx<- strsplit(MedlinePgn(pub),signPat)[[1]]
    if(nchar(pageMnMx[1])!=nchar(pageMnMx[2])){
      pages <- paste(pageMnMx[1],"-",substr(pageMnMx[1],start=1,stop=nchar(pageMnMx[1])-nchar(pageMnMx[2])),pageMnMx[2],sep="")
    } else {
      pages <- MedlinePgn(pub)
    }
    pageSeq <- try.CATCH(seq(as.numeric(strsplit(pages,signPat)[[1]][1]),as.numeric(strsplit(pages,signPat)[[1]][2])))
    if(grepl("Error",strsplit(as.character(pageSeq$value),"\\s"))[[1]]){pageSeq <- pages} else {pageSeq <- pageSeq$value}
    APAcite <- paste(paste(authorString,collapse=""),ArticleTitle(pub)," ",Title(pub),", ",Volume(pub),"(",Issue(pub),"), ",pages,".",sep="")
    cat("Found citation for file [",names(IT),"]\n--",APAcite,"\n")
    return(list(Source=pub,pageSeq=pageSeq,reference=APAcite,kind="PubMed"))
  }
}

evalSeq <- function(targetSeq,foundSeq){
  # Decide to what extent the page sequence was found
  found <- ""
  if(any(targetSeq%in%unique(foundSeq))){
    if(all(targetSeq%in%unique(foundSeq))){
      if(length(foundSeq)==length(targetSeq)){found <- "Exact"} else {found <- "Exact+"}
    } else {
      if(length(foundSeq)<=length(targetSeq)){found <- "PartSeq-"} else {found <- "PartSeq+"}
    }
  } else {
    found <- "None"
  }
  return(found)
}

loadPatterns <- function(skelet=T,stats=T,evid=T,eff=T,desc=T,bib=T,punct=T,cln=F){
  # Ready to use regexpr patterns to find common sections of articles, stats, etc.
  require(Unicode)
  sup2   <- paste0("([2",intToUtf8(u_char_from_name("SUPERSCRIPT TWO")),"])")
  min    <- paste0("([",intToUtf8(u_char_from_name("HYPHEN-MINUS")),intToUtf8(u_char_from_name("HYPHEN")),intToUtf8(u_char_from_name("MINUS SIGN")),"])")
  
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
  
  # Numbers
  int_pat  <- paste("(\\d+(","(\\,)?","\\d+){0,2})",sep="[[:blank:]]?")
  real_pat <- paste("(\\d*(","(\\,|\\.)?","\\d+){0,2})",sep="[[:blank:]]?")
  # Characters encountered after 'pdftotxt' encoding for minus sign
  ifelse(cln,{sign_pat <- paste0("(",gsub("\\s","[[:blank:]]",loadCleanPatterns()[["cleanPunctuationPatterns"]]["sign"]),")?")},{sign_pat <- "([-–ϪÀ±~])?"})
  # Characters encountered after 'pdftotxt' encoding for (in)equality
  ifelse(cln,{eq_pat <- paste0("(([<>≤≥≦≧≨≩≪≫≭≮≯≰≱≲≳≴~≈≉≠≢])|(",gsub("\\s","[[:blank:]]",loadCleanPatterns()[["cleanPunctuationPatterns"]]["equal"]),"))")},{eq_pat <- "((?<unequal>[<>≤≥≦≧≨≩≪≫≭≮≯≰≱≲≳≴~≈≉≠≢])|(?<equal>([=≡≣ϭϭ])|([[:blank:]]5[[:blank:]])|(\\xf4\\x8f\\xb0\\x81)|(\\x{03fd})|(\\x{03fe})))"})
  # Characters encountered after 'pdftotxt' encoding for non significant
  ns_pat  <- "((n(\\.)?s(\\.)?)|(no(n|\\s+)*sig(\\.)?))"
  # Combinations of punctuation, nonprinting and numbers
  ifelse(cln,{test_pre <- "("},{test_pre <- "(?<test>"})
  test_pat<- paste0(test_pre,"(\\(([[:blank:]]*",int_pat,"+[[:blank:]]*(\\,)?[[:blank:]]*([nN])?[[:blank:]]*(",eq_pat,")?[[:blank:]]*){1,2}\\))([[:blank:]]*",eq_pat,"[[:blank:]]*",sign_pat,"[[:blank:]]*(",real_pat,"|",ns_pat,")+){1})")
  ifelse(cln,{nodf_pre <- "("},{nodf_pre <- "(?<nodf>"})
  nodf_pat<- paste0(nodf_pre,"[[:blank:]]*",eq_pat,"[[:blank:]]*",sign_pat,"[[:blank:]]*(",real_pat,"|",ns_pat,")+){1}")
  page_pat<- paste("[^[:print:]]*((\\d+)",")[^[:print:]]*",sep="\\w{0,2}")
  head_pat<- paste("((((?-i)(I|V|X|L){1,5}(?i))|\\d{0,2}|\\w{0,1}){1,2}(","[\\.:-]",")?)",sep="[[:blank:]]*")
  
  skeletonPatterns <- list(
    # Ordered header indicators..."Experiment 1", "Study 4b", "Appendix B"
    header=paste("(?<=[^[:print:]]{1})(?<header>",head_pat,"?((?i)(Replication|Meta[-]?|Data|Exploratory|Confirmatory|Strateg(y|i[ec])|Trial|Comput(er|ation(al)*)|Model|Network|Experiment|Stud(y|ie)|Analys(i|e)|Simulation|Appendi(x|ce))(?-i)(s)?(","(and|&)?",")){1,2}",head_pat,"?)(?=[^[:print:]]{1})+",sep="[[:blank:]]*"),
    # Common section / paragraph indicators..."Method", "Procedure", "Results"
    section=paste("(?<=[^[:print:]]{1})(?<section>",head_pat,"?((General|Final|Overall|Summar(y|i[sz]ing))?","(Abstract|Summary|Introduction|Method|Procedure|Material|Instrument|Participant|Subject|Sample|Patient|Result|Conclusion|Discussion)+(s)?(","(and|&)?",")){1,2}",head_pat,"?){1}(?=[^[:print:]]{1})",sep="[[:blank:]]*"),
    # Tables + titles
    table=paste("(?<=[^[:print:]]{1})(?<table>",head_pat,"?(Table(s)?){1}",paste0(head_pat,"+"),")",sep="[[:blank:]]*"),
    #(?<title>[[:print:]]+)?(?=[^[:print:]]{1}){1}",sep="[[:blank:]]*"),
    # Figures + titles
    figure=paste("(?<=[^[:print:]]{1})(?<figure>",head_pat,"?((Figure|Box|Diagram|Frame)(s)?){1}",paste0(head_pat,"+"),")",sep="[[:blank:]]*"),
    #(?<title>[[:print:]]+)?(?=[^[:print:]]{1}){1}",sep="[[:blank:]]*"),
    # Page numbers
    page=paste("(?<=[^[:print:]]{1})",page_pat,"(?=[^[:print:]]{1})",sep="[[:blank:]]*")
  )
  
  statsPatterns <- list(
    #     # Regexpr patterns to find statistics
    #     t = paste("((?<t>[tT])+","[s]?",test_pat,"){1}",sep="[[:blank:]]*"),
    #     F = paste("((?<F>[fF])+","[s]?",test_pat,"){1}",sep="[[:blank:]]*"),
    #     r = paste("((?<r>[rRρΡ])+","[s]?",test_pat,"){1}",sep="[[:blank:]]*"),
    #     Z = paste("((?<!(Sobel))","(?<Z>\\b[zZ]\\b)+","[s]?",nodf_pat,"){1}",sep="[[:blank:]]*"),
    #     sz= paste("((?<sz>(Sobel)+","(['`]","[s])?","([zZ])+)+","[s]?",nodf_pat,"){1}",sep="[[:blank:]]*"),
    #     W = paste("((?<W>[wW](ilcoxon)?)+",nodf_pat,"){1}",sep="[[:blank:]]*"),
    #     U = paste("((((Mann)?[-](Whitney)?)+","(?<U>[uU])+)+",nodf_pat,"){1}",sep="[[:blank:]]*"),
    #     chi2 = paste("((?<chi2>[cCχΧx](hi)?",sup2,"+)+","[s]?",test_pat,"){1}",sep="[[:blank:]]*"),
    #     B    = paste("((?<B>\\b[bB])+","[s]?",test_pat,"){1}",sep="[[:blank:]]*"),
    #     beta = paste("((?<beta>[β␤]([bB]eta)?)+","[s]?",test_pat,"){1}",sep="[[:blank:]]*"),
    #     alpha = paste("((?<alpha>[αaΑ](lpha)?)+","[s]?",nodf_pat,"){1}",sep="[[:blank:]]*")
    # Regexpr patterns to find statistics
    t = paste("(?<t>\\b[tT]+","[s]?",")",sep="[[:blank:]]*"),
    F = paste("(?<F>\\b[fF]+(1|2)?","[s]?",")",sep="[[:blank:]]*"),
    r = paste("(?<r>\\b(Pearson|Spearman)?","[s]?","[rRρΡ]+","[s]?",")",sep="[[:blank:]]*"),
    Z = paste("(?<Z>(?<!Sobel)","[s]?","\\b[zZ]\\b","[s]?",")",sep="[[:blank:]]*"),
    sz= paste("(?<sz>Sobel","['`]?","[s]?","[zZ]","[s]?",")",sep="[[:blank:]]*"),
    W = paste("(?<W>\\b[wW](ilcoxon)?",")",sep="[[:blank:]]*"),
    U = paste("(?<U>[Mann[-]Whitney]?","\\b[uU]",")",sep="[[:blank:]]*"),
    chi2 = paste("(?<chi2>\\b([wWcCχΧx](hi)?",sup2,"([\\.-]sq[uare])?","[s]?",")|(\u24392)|(chi-square))",sep="[[:blank:]]*"),
    B    = paste("(?<B>\\b[bB]+","[s]?",")",sep="[[:blank:]]*"),
    beta = paste("(?<beta>\\b([β␤]|[bB]eta)+","[s]?",")",sep="[[:blank:]]*"),
    alpha = paste("(?<alpha>\\b([α]|[aΑ]lpha)+","[s]?",")",sep="[[:blank:]]*")
  )
  
  evidencePatterns <- list(
    #    p   = paste("(?<p>\\b([pP])+","((rep)|[[:print]]{1,3})?",nodf_pat,"){1}",sep="[[:blank:]]*")
    #     ll  = paste("((?<ll>",sign_pat,"2[[:print:]]{0,2}","([lL]og)?","([lL]ik(elihood)?)?|([lL]){2}","[d]?)",nodf_pat,"){1}",sep="[[:blank:]]*"),
    #     inf = paste("((?<inf>[ABD]IC|(information","criterion)",")",nodf_pat,"){1}",sep="[[:blank:]]*"),
    #     eig = paste("((?<eig>[eE]igen","(value)?",")",nodf_pat,"){1}",sep="[[:blank:]]*")
    #     
    p  = paste("(?<!",sup2,")[[:blank:]]*(?<p>\\b[pP]+","[[:blank:]]*([\\.]?rep)?","[[:blank:]]*[s]?)",sep=""),
    ll  = paste("(?<ll>",sign_pat,"2[[:print:]]{0,2}","([lL]og)?","([lL]ik(elihood)?)?|([lL]){2}","[d]?)",sep="[[:blank:]]*"),
    inf = paste("(?<inf>[ABD]IC|(information","criterion)",")",sep="[[:blank:]]*"),
    eig = paste("(?<eig>[eE]igen","(value)?",")",sep="[[:blank:]]*")
  )
  
  effectPatterns <- list( 
    # Effect sizes
    #     MSE  = paste("((?<MSE>MSE)+",nodf_pat,"){1}",sep="[[:blank:]]*"),
    #     omega2 = paste("((?<omega2>([ωΩ])+",sup2,"+)",nodf_pat,"){1}",sep="[[:blank:]]*"),
    #     eta2 = paste("((?<=[[:punct:]])","(?<eta2>(partial)?","(eta)?","([ηΗgzZ])+","([pₚ])?",sup2,"+([pₚ])?",")+",nodf_pat,"){1}",sep="[[:blank:]]*"),
    #     eta = paste("((?<=[[:punct:]])","(?<eta>(eta)?",paste0("([ηΗgzZ])(?!",sup2,")+)+"),nodf_pat,"){1}",sep="[[:blank:]]*"),
    #     phi = paste("((?<phi>([cC]ramer","(['`]","[s])?)","([φΦ]|[pP]hi)+","[cC]?)+",nodf_pat,"){1}",sep="[[:blank:]]*"),
    #     ICC = paste("((?<ICC>[iIcC])+",nodf_pat,"){1}",sep="[[:blank:]]*"),
    #     psi = paste("((?<psi>([Ψψ]|[pP]s(i|y)))+",nodf_pat,"){1}",sep="[[:blank:]]*"),
    #     r2  = paste("((?<r2>([rRρΡ])+","(change)?",sup2,"+){1}",nodf_pat,"){1}",sep="[[:blank:]]*"),
    #     res = paste("((?<res>(effect","size)?",paste0("([rRρΡ])(?!",sup2,"))+"),nodf_pat,"){1}",sep="[[:blank:]]*"),
    #     OR = paste("((?<OR>l(n|og)[(]?","[oO]([rR]|dds)","[)]?)+",nodf_pat,"){1}",sep="[[:blank:]]*"),
    #     RR = paste("((?<RR>l(n|og)[(]?","[rR]","[)])+",nodf_pat,"){1}",sep="[[:blank:]]*"),
    #     gd = paste("((?<gd>(([gG]lass){1}","(['`]","[s])?)+","([δΔd])+)",nodf_pat,"){1}",sep="[[:blank:]]*"),
    #     hv = paste("((?<hv>([hH]edges","(['`]","[s])?)+","([gG])+)",nodf_pat,"){1}",sep="[[:blank:]]*"),
    #     cd = paste("((?<cd>([cC]ohen","(['`]","[s])?)?","([δΔd])+)",nodf_pat,"){1}",sep="[[:blank:]]*"),
    #     cf = paste("((?<cf>([cC]ohen","(['`]","[s])?)?","([fF])+","([2²])+)+",nodf_pat,"){1}",sep="[[:blank:]]*")
    omega2 = paste("(?<omega2>[ωΩ]",sup2,"){1}",sep="[[:blank:]]*"),
    eta2 = paste("((?<eta2>(eta)?","([ηΗgzZ]){1}","(p)?",sup2,")|(\u24292",sup2,"?[p]?)","){1}",sep="[[:blank:]]*"),
    eta =  paste("((?<eta>(eta)?","([ηΗgzZ]){1}",paste0(gsub("\\[","[^",sup2,perl=T)),")|(\u24292)){1}",sep="[[:blank:]]*"),
    phi = paste("(?<phi>([cC]ramer","([']","[s])?)?","([φΦ]|([pP]hi)){1}","[cC]?){1}",sep="[[:blank:]]*"),
    #    ICC = paste("(?<ICC>[iIcC]){1}",sep="[[:blank:]]*"),
    psi = paste("(?<psi>([Ψψ]|([pP]s(i|y)))){1}",sep="[[:blank:]]*"),
    r2  = paste("(?<r2>([rRρΡ])","(change)?",sup2,"){1}",sep="[[:blank:]]*"),
    #    res = paste("(?<res>(effect","size)?",paste0("([rRρΡ])(?!",sup2),")){1}",sep="[[:blank:]]*"),
    #OR = paste("(?<OR>((l(n|og))[(]?","([oO]([rR]|dds)){1}","[)]?)|([oO]dds\\s[rR]atio)){1}",sep="[[:blank:]]*"),
    OR = paste("(?<OR>([oO]dds","[rR]atio)){1}",sep="[[:blank:]]*"),
    #     RR = paste("(?<RR>(l(n|og))[(]?","([rR]){1}","[)]){1}",sep="[[:blank:]]*"),
    #     gd = paste("(?<gd>([gG]lass","([']","[s])?)?","([δΔd])){1}",sep="[[:blank:]]*"),
    #     hv = paste("(?<hv>([hH]edges","([']","[s])?)?","([gG])){1}",sep="[[:blank:]]*"),
    #     cd = paste("(?<cd>([cC]ohen","([']","[s])?)?","([δΔd])){1}",sep="[[:blank:]]*"),
    #     cf = paste("(?<cf>([cC]ohen","([']","[s])?)?","([fF]){1}",sup2,"){1}",sep="[[:blank:]]*")
    #     gd = paste("((?<gd>(([gG]lass){1}","([']","[s])?)+","([δΔd])+)",nodf_pat,"){1}",sep="[[:blank:]]*"),
    #     hv = paste("((?<hv>([hH]edges","([']","[s])?)+","([gG])+)",nodf_pat,"){1}",sep="[[:blank:]]*"),
    cd = paste("((?<cd>([cC]ohen","([']","[s])?)?","([δΔd])+)",nodf_pat,"){1}",sep="[[:blank:]]*")
    #    cf = paste("((?<cf>([cC]ohen","([']","[s])?)?","([fF])+","([2²])+)+",nodf_pat,"){1}",sep="[[:blank:]]*")
    
  )
  
  # Regexpr patterns to find information about analyses
  descriptivePatterns <- list(
    # Sample descriptives
    sde = paste("(([σΣ])|([sS](\\.)?([dD]|[eE])(\\.)?))+",nodf_pat, sep="[[:blank:]]*"),
    var = paste("(([σΣ])|([sS](\\.)?[dD](\\.)?)[2²]?)+",nodf_pat, sep="[[:blank:]]*"),
    MSE  = paste("(?<MSE>MSE){1}",nodf_pat, sep="[[:blank:]]*"),
    M   = paste("(([μΜ])|([mM](\\.)?[nN]?(\\.)?)|([mM]ean|[aA]verage|[mM]edian))+[s]?",nodf_pat,sep="[[:blank:]]*"),
    N = paste("([nN]",eq_pat,")",int_pat,sep="[[:blank:]]*"),
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
    authors  = "(?<authors>([[:upper:]]{1,2}[[:lower:]]+[\\,])+([[:upper:]]?|[[:punct:]]|[[:blank:]]|(et\\s*al[\\.]))+){1,6}",
    citation = "",
    reflist  = "\\b((reference(s)*(list)*)|bibliography)(?=\\n)",
    reference= paste("(?<=[\\d\\.\\n\\s])((?-i)[[:upper:]](?i)[[:print]]+([(]\\d{4}[)]){1})[[:print]]+(",page_pat,sign_pat,page_pat,")?[\\.]{1}([[:blank:]])*",sep=""),
    pubyear = "(?<year>([(](\\d{4}\\w{0,1})|((in\\s*press)|(n[\\.]?d[\\.]?)|(submitted))[)])+[[:blank:]]*[[:punct:]]*){1}",
    keywords = "((key[-]?|index)*(word|term|topic)(s)*[[:blank:]]*[\\:]?)",
    authorloc  = "",
    authornote = "",
    funding  = ""
  )
  
  # Punctuation patterns for cleaning purposes
  punctuationPatterns <- list(
    space_pat_rep = paste("(?<space_rep>(","{2,}))",sep="[[:blank:]]"),
    space_pat_del = paste("(?<space_del>((?<=[(])","+)|(","+(?=[),])|(^","+)|(","+$)))",sep="[[:blank:]]"),
    sign_pat_cln  = paste("(?<sign>",gsub("\\?","",sign_pat),"[[:blank:]]*)(?=",real_pat,")",sep=""),
    eq_pat_cln    = paste("(?<=[[:graph:]])",paste0("(",paste0(strsplit(x=eq_pat,"\\|",perl=T)[[1]][-1],collapse="|")),paste0("(?=[[:digit:]\\.\\,",gsub("((\\?)|(\\(\\[)|(\\]\\)))+","",sign_pat),intToUtf8(u_char_from_name("MINUS SIGN")),"nNsS])"),sep="[[:blank:]]*"),
    punct_pat_del = paste("(?<punct_del>((^","+)|(","+$)))",sep="[[:punct:]]")
  )
  
  return(list(skeletonPatterns=skeletonPatterns,
              statsPatterns=statsPatterns,
              evidencePatterns=evidencePatterns,
              effectPatterns=effectPatterns,
              descriptivePatterns=descriptivePatterns,
              bibPatterns=bibPatterns,
              punctuationPatterns=punctuationPatterns,
              c(int_pat=int_pat,real_pat=real_pat,sign_pat=sign_pat,eq_pat=eq_pat,ns_pat=ns_pat,page_pat=page_pat,test_pat=test_pat,nodf_pat=nodf_pat,head_pat=head_pat)
  )[c(skelet,stats,evid,eff,desc,bib,punct,TRUE)]
  )
}


loadCleanPatterns <- function(skelet=T,stats=T,evid=T,eff=T,desc=T,bib=T,punct=T){
  require(Unicode)
  
  cleanSkeleton <- list()
  cleanStats <- list(
    # Regexpr patterns to find statistics
    t = "t",
    F = "F",
    r = "r",
    Z = "Z",
    sz= "Sobel's z",
    W = "W",
    U = "U",
    chi2 = "χ²",
    B    = "B",
    beta = "β",
    alpha = "α"
  )
  
  cleanEvidence <- list(
    p = "p",
    ll = "Log Likelihood",
    inf = "Information Criterion",
    eig = "Eigenvalue"
  )
  
  cleanEffects <- list(
    MSE  = "MSE",
    omega2 = "ω²",
    eta2 = "η²",
    eta = "η",
    phi = "Cramer's Φ",
    ICC = "ICC",
    psi = "ψ",
    r2 = "r²",
    res = "r_es",
    OR = "OR",
    RR = "RR",
    gd = "Glass' Δ",
    hv = "Hedges' g",
    cd = "Cohen's d",
    cf = "Cohen's f²"
  )
  
  cleanDescriptives <- list()
  
  cleanBib <- list()
  
  cleanPunct <- list(
    space_rep = intToUtf8(u_char_from_name("SPACE")),
    space_del = "",
    sign      = paste0(intToUtf8(u_char_from_name("SPACE")),intToUtf8(u_char_from_name("MINUS SIGN"))),
    equal     = paste0(intToUtf8(u_char_from_name("SPACE")),intToUtf8(u_char_from_name("EQUALS SIGN")),intToUtf8(u_char_from_name("SPACE"))),
    punct_del = ""
  )
  
  return(list(cleanSkeletonPatterns=cleanSkeleton,
              cleanStatsPatterns=cleanStats,
              cleanEvidencePatterns=cleanEvidence,
              cleanEffectPatterns=cleanEffects,
              cleanDescriptivePatterns=cleanDescriptives,
              cleanBibPatterns=cleanBib,
              cleanPunctuationPatterns=cleanPunct
  )[c(skelet,stats,evid,eff,desc,bib,punct)]
  )
}


splitSeq <- function(seq){# Split a vector of numbers into sequences of lag = 1 (if any)
  seqLoc <- attr(seq,"location")
  seq    <- as.numeric(gsub("\\D","0",seq)) #unique(sort(as.numeric(gsub("\\D","0",seq))))
  seq    <- seq[seq>0]
  cutat  <- c(1,(which(diff(seq,lag=1)!=1)+1)) #c(1,(which(diff(sort(seq),lag=1)!=1)+1))
  if(length(cutat)>0){
    cutsz   <- diff(c(cutat,(length(seq)+1)))
    splitz  <- split(seq,unlist(lapply(seq_along(cutat),function(s) rep(cutat[s],times=cutsz[s]))))
    #split(sort(seq),unlist(lapply(seq_along(cutat),function(s) rep(cutat[s],times=cutsz[s]))))
    splitzLoc <- split(seqLoc,unlist(lapply(seq_along(cutat),function(s) rep(cutat[s],times=cutsz[s]))))
    seqname   <- sapply(seq_along(splitz),function(s) paste("Split sequence ",s,sep=""))
    names(splitz) <- seqname
    attr(splitz,"location") <- splitzLoc
  } else {
    splitz <- seq
    attr(splitz,"location") <- seqLoc
    if(length(cutat)==1){
      names(splitz) <- "Ordered vector is 1 sequence!"
    }
    names(splitz) <- "No sequences in this vector"
  }
  return(splitz)
}

# SEARCHERS ---------------------------------------------------------------

subIT <- function(clnIT,pat,cleanpat=NULL){
  require(plyr)
  if(is.null(cleanpat)){cleanpat <- loadCleanPatterns()}
  if(is.list(clnIT)){clnIT<-unlist(clnIT,recursive=F,use.names=F)}
  
  cleanTYPE <- gsub("\\w+(\\.){1}","",unlist(names(unlist(cleanpat))))
  logicTab  <- matrix(laply(pat,grepl,x=clnIT,perl=TRUE),nrow=length(pat),ncol=length(clnIT))
  
  if(any(rowSums(logicTab)>0)){
    for(p in which(rowSums(logicTab)>0)){
      LOC    <- laply(clnIT[which(logicTab[p,])],locIT,str=pat[[p]])
      idTYPE <- laply(seq_along(LOC),function(l) which(cleanTYPE%in%attr(LOC[[l]],"capture.names")))
      idLOC  <- laply(seq_along(LOC),function(l) which(attr(LOC[[l]],"capture.names")%in%cleanTYPE[idTYPE]))                                       
      LOCs   <- llply(seq_along(LOC),function(id) as.data.frame(cbind(start=attr(LOC[[id]],"capture.start")[,idLOC[id]],stop=attr(LOC[[id]],"capture.length")[,idLOC[id]])))
      
      ifelse(length(which(logicTab[p,]))>1,{tmp <- clnIT[which(logicTab[p,])]},{tmp <- clnIT[[which(logicTab[p,])]]})
      
      if(length(LOCs)>0){
        for(lc in seq_along(LOCs)){
          for(c in 1:nrow(LOCs[[lc]])){
            if(LOCs[[lc]][["start"]][c]!=-1){
              substr(tmp[c], start=LOCs[[lc]][["start"]][c], stop=LOCs[[lc]][["start"]][c]+(LOCs[[lc]][["stop"]][c]-1)) <- paste0(rep("#",times=LOCs[[lc]][["stop"]][c]),collapse="")}
          } 
        } 
        clnIT[which(logicTab[p,])] <- laply(seq_along(tmp),function(t) gsub("[#]+",paste0(unlist(cleanpat,use.names=F)[idTYPE[t]]),tmp[t]))
      }
      rm(tmp, LOC, idLOC, idTYPE, LOCs)  
    }
  }
  return(clnIT) 
}

cleanIT <- function(IT,searchPat="punctuationPatterns",cleanPat=paste0("clean",gsub("(\\w)(\\w+)","\\U\\1\\E\\2",searchPat,perl=T))){
  Patterns      <- loadPatterns()[[searchPat]]
  cleanPatterns <- loadCleanPatterns()[[cleanPat]] 
  namedPatterns <- gsub("\\w+(\\.){1}","",unlist(names(unlist(cleanPatterns))))
  Patterns      <- gsub("\\(\\?\\<(nodf|test)\\>([[:print:]]|[^[:print:]])+",")",Patterns,perl=T)
  
  cat("\nCleaning up",searchPat,"in text extracted from ",length(IT)," items\n")
  if(length(IT)>10){cat("This may take a while...\n")}
  for(p in seq_along(Patterns)){
    locs <- gregexpr(pattern=Patterns[p],text=IT,perl=TRUE)
    regmatches(IT,locs) <- rep(cleanPatterns[p],length(locs)) 
    #result <- gregexpr(pattern=Patterns[p],text=IT,perl=TRUE)
    #subs <- regmatches(IT,result)
  }
  return(IT)
}

# cleanIT <- function(IT){# Clean up the output from matchIT()
#   regpat   <- loadPatterns()
#   cleanpat <- loadCleanPatterns()
#   
#   for(i in seq_along(IT)){
#     # Clean the minus and equal sign
#     #     IT[i] <- gsub(paste("(?<sign>",gsub("\\?","",regpat[[7]]["sign_pat"]),")(?=",regpat[[7]]["real_pat"],")",sep="[[:blank:]]*")," -", IT[i], perl=TRUE)
#     IT[[i]] <- subIT(IT[i],pat=regpat[["punctuationPatterns"]],cleanpat[["cleanPunctuationPatterns"]])
#     IT[[i]] <- subIT(IT[i],pat=regpat[["statsPatterns"]],cleanpat[["cleanStatsPatterns"]])
#     IT[[i]] <- subIT(IT[i],pat=regpat[["effectPatterns"]],cleanpat[["cleanEffectPatterns"]])
#     #IT[[i]] <- subIT(IT[i],pat=regpat[["punctuationPatterns"]])
#     
#   }
#   return(list(IT))
# }

tidyIT <- function(IT,cln=FALSE){# Tidy up an IT)
  if(cln){cleantIT(IT)}
  for(i in seq_along(IT)){
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

matchIT <- function(pat,IT, ic=TRUE,cln=NULL,loc=TRUE) {# Return matching strings in source using regmatches
  require(reshape2)
  n <- names(IT)
  m <- gregexpr(pat,IT,perl=TRUE,ignore.case=ic) 
  #  t <- regmatches(IT,m)
  t <- list()
  for(i in seq_along(IT)){
    tmp <- regmatches(IT[[i]],m[i])
    #cat(n[i],"\n")
    ifelse(is.null(cln),{t[[i]] <- tmp[[1]]},{t[[i]] <- cleanIT(tmp[[1]],searchPat=cln)})
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

citedCheck <- function(signPat,pagePat,ppMX,ppMN,IT){
  # Try to infer page range from entries such as 'XXX-YYY = #Pages' cited in the PDF text
  require(plyr)
  tmpSeq  <- NULL 
  reflist <- locIT("\\b((reference(s)*(list)*)|bibliography)(?=\\n)",IT[[1]])[[1]]
  if(length(reflist)>0){
    reflist<-min(reflist)
  } else {
    reflist<-round(nchar(IT)*.95)
  }
  ITnoref <- substr(IT[[1]],1,reflist)
  nrs     <- gsub("[[:alpha:]]","",gsub(paste(signPat,sep=")"),"-",unlist(matchIT(paste(pagePat,signPat,pagePat,sep=""),ITnoref,cln=F))))
  if(length(nrs)!=0){
    diffP <- (abs(ldply(as.quoted(nrs),eval))+1)
    if(length(nrs)!=0){ids <- which(((diffP>=ppMN)&(diffP<=ppMX)))} else {ids <- integer(0)}
  }
  if(length(ids)!=0){
    tmpSeq   <- sapply(nrs[ids],strsplit,split=signPat)
    tmpSeq   <- sapply(tmpSeq,function(sq) seq(as.numeric(sq[1]),as.numeric(sq[2]),by=1))
  } else {
    tmpSeq <- NULL
  }
  return(tmpSeq)
}
# 
# evaluatePages <- function(tmpSeqs,IT,nPages){
#   # First try to find the major body of text (without references)
#   reflist <- locIT("\\b((reference(s)*(list)*)|bibliography)(?=\\n)",IT[[1]])[[1]]
#   if(length(reflist)>0){
#     reflist<-min(reflist)
#   } else {
#     reflist<-round(nchar(IT)*.95)
#   }
#   ITnoref <- substr(IT[[1]],1,reflist)
#   
#   # Average expected page length
#   avg_page1 <- round(nchar(IT)/nPages)
#   avg_page2 <- round(nchar(ITnoref)/nPages)
#   
#   # Now find the numbers in tmpseq surrounded by formatting characters
#   
#   for(s in seq_along(tmpSeqs)){
#     # last page in the sequence?
#     
#     lastPage <- min(abs(attr(strings[[s]],"location")-nchar(ITnoref)))
#     
#     if(
#       mean(unique(attr(strings[[2]],"location")))<avg_page1
#     )
#       
#   }
# }


rangeCheck <- function(splitsFound,nExpected,closeFITdiff=2){
  # 1. Are certain seq combinations a close fit with expected page range: ((sum(sequence lengths) - nPages) <= closeFITdiff)?
  # 2. Are certain ranges not found in an ordered sequence?
  if(length(splitsFound) > 1){
    combi    <- unlist(llply(2:(length(splitsFound)-1), function(c) combn(seq_along(splitsFound),c,simplify=F)),recursive=F)
  } else {
    combi <- list(1)
  }
  dCombi   <- ldply(combi, function(c){vals<-unlist(splitsFound[c]); return(nExpected-length(which(diff(vals)==1)))}) 
  if(any(dCombi$V1<=closeFITdiff)){
    tempSeqID<- which(dCombi$V1<=min(dCombi$V1))
  } else {
    tempSeqID<-integer(0)
  }
  if(length(tempSeqID)>0){
    tempSeqID <- tempSeqID[which.min(laply(combi[tempSeqID],length))]
    tmpSeqs <- unlist(splitsFound[combi[[tempSeqID]]])
    tmpLocs <- attributes(splitsFound)$location[combi[[tempSeqID]]] 
  } else {
    tmpSeqs <- unlist(splitsFound)
    tmpLocs <- attributes(splitsFound)$location
  }
  names(tmpSeqs) <- seq_along(tmpSeqs)
  tmpLocs        <- unlist(tmpLocs)
  names(tmpLocs) <- tmpSeqs
  attr(tmpSeqs,"location") <- tmpLocs
  return(tmpSeqs)
}

# Get location of page nrs. in units of IT size
locatePages  <- function(seqF,prepost,ITbit,norm=nchar(ITbit)){
  locs <- vector("list")
  for(i in seq_along(seqF)){locs[[i]] <- 
                              ldply(seqF,function(n) as.data.frame(cbind(page=n, 
                                                                         loc=as.vector(locIT(wordXX(pre=paste(prepost[1],"(",sep=""),post=paste(")",prepost[2],sep=""),x=seqF),ITbit)[[1]]/nchar(ITbit))
                              )))
  }
  names(locs) <- names(seqF)
  return(ldply(locs))
}

# 
# str
#PDFname
# PDFid <- 7
# ITbit <- txtfiles[[PDFid]]
# closeFITdiff=NA
# closeFITmiss=NA

getPageSeq <- function(ITbit,PDFid,ITbitINFO,closeFITdiff=NA,closeFITmiss=NA){
  # Lot of seting up to do...
  # Define some tools (requires plyr)
  require(plyr)
  f_l <- each(length)
  f_m <- each(min,max)
  
  PDFname <- attr(ITbit,"fullName")
  #ITbit <- ITbit[[1]]
  
  # Get sign regex pattern, without the ?
  signPat <- sub("\\?","",loadPatterns(rep(FALSE,times=6))[[1]]["sign_pat"])
  # Get the page nr. pattern and search string
  pagePat <- loadPatterns(rep(FALSE,times=6))[[1]]["page_pat"]
  str     <- loadPatterns(c(TRUE,rep(FALSE,times=5)))[["skeletonPatterns"]][["page"]]
  # Seperate the look ahead/behind conditions
  prepost <- strsplit(str,pagePat,fixed=T)[[1]]
  
  #Absolute min and max nr. of pages
  ppMX <- 50
  ppMN <- 1
  
  # Min and Max perc. values for first and last pagenumber locations
  minP <- .11
  maxP <- .94
  
  # Do we have a pagerange?
  if(length(ITbitINFO[[PDFid]]$pageSeq)==0){
    
    # Get info on number of pages in the pdf from 'pdfinfo'
    nPages <- try.CATCH(getPDFpageinfo(PDFname))
    if((length(nPages$warning)>0)|!(is.numeric(nPages$value))){
      cat("\n\nPDFinfo failed")
    } else {
      nPages <- nPages$value
    }
    # Try to get more info, sometimes 'pdfinfo' will return nPages that is much larger than the citation (e.g. suppl. material) 
    # If such is the case, the article citation will often be somewhere in the pdf
    tmpSeqs <- citedCheck(signPat,pagePat,ppMX,ppMN,ITbit)
    
  } else {
    newSeq <- ITbitINFO$pageSeq
    
    #Absolute min and max nr. of pages
    ppMX <- max(newSeq)
    ppMN <- min(newSeq)
    
    pageSeq <- splitSeq(newSeq)
    locSplits   <- locatePages(pageSeq,prepost,ITbit[[1]],norm=nchar(ITbit[[1]]))
    matchIT(wordXX(pre=paste(prepost[1],"(",sep=""),post=paste(")",prepost[2],sep=""),x=locSplits$page),ITbit[[1]])
    locSplits$loc*10
  }
  
  fst_page  <- round(minP*nchar(ITbit))
  lst_page  <- round(maxP*nchar(ITbit))
  avg_page  <- round(nchar(ITbit)/nPages)
  
  # Set how many page numbers can be missing (e.g., the odd pages may be missing)
  if(is.na(closeFITmiss)){if(nPages>10){closeFITmiss=5} else {closeFITmiss=4}}
  # Set how many values can be missing consequtively from a sequence
  if(is.na(closeFITdiff)){if(nPages>10){closeFITdiff=3} else {closeFITdiff=2}}
  # Set a range of sequence lengths that qualifies as close, but...
  NOcigar <- ((nPages-(closeFITmiss%/%2)):(nPages+(closeFITmiss%/%2)))
  
  # Sequence (s) found in IT using str
  sFound  <- sort(unique(as.numeric(gsub("\\D","0",matchIT(str,ITbit[[1]],cln=F)[[1]])),na.last=NA))
  sFound  <- sFound[sFound>0]
  # Lengths (n) of found sequences
  nFound  <- f_l(sFound)
  if(nFound==0){stop("No page references found in PDF using current search pattern.")}
  # Range (r) of found sequence
  rFound  <- f_m(sFound)
  
  # Check true if any missing info is found
  foundD  <- FALSE
  foundR  <- FALSE
  foundS  <- FALSE
  
  # Range within allowed range?
  if(length(rFound)%in%NOcigar){foundD=TRUE}
  SeqMiss <- diff(sFound)
  if((length(SeqMiss)>0)&(all(SeqMiss<=closeFITmiss))){foundD=TRUE}
  
  # See if multiple (sensible) sequences returned in sFound
  splitsFound <- splitSeq(sFound)
  ranges      <- ldply(splitsFound,f_m)
  d           <- (ranges$max-ranges$min)
  d[d==0]     <- 1
  splitsFound <- splitsFound[(d<=ppMX)&(d>=ppMN)]
  
  # Get locations
  locSplits   <- locatePages(splitsFound,prepost,ITbit,norm=nchar(ITbit))
  
  # Do a first check on the page ranges found 
  tmpSeqs <- rangeCheck(splitsFound,locSplits,nPages,closeFITdiff)
  if(length(tmpSeqs)>0){
    locSplits <- locatePages(tmpSeqs,prepost,ITbit,norm=nchar(ITbit))
  } else {
    tmpSeqs <- splitsFound
  }
  newSeq    <- min(unlist(tmpSeqs)):max(unlist(tmpSeqs))
  if(length(newSeq)>ppMX*2){stop("No page references found in PDF using current search pattern.")}
  
  # BEGIN CHECKS
  #ptm <- proc.time()[3]
  while(!(length(newSeq)%in%NOcigar)){
    
    # First try to find a range of numbers based on page locations in the PDF file
    # Is a likely first or last pagenr. in the sequence?
    if(length(unique(locSplits$.id))>1) {ids <- laply(unique(locSplits$.id),function(u) locSplits$.id%in%u)
    } else {
      ids <- 1
    }
    if(length(ids)>0){
      firstP <- unlist(llply(1:length(ids), function(s) locSplits$page[ids[s]][locSplits$loc[ids[s]]<=minP]))
      if(length(firstP>0)){firstP<-min(firstP)} else {firstP<-integer(0)}
      lastP  <- unlist(llply(1:length(ids), function(s) locSplits$page[ids[s]][locSplits$loc[ids[s]]>=maxP]))
      if(length(lastP>0)){lastP<-min(lastP)} else {lastP<-integer(0)}
    }
    
    # Evaluate results
    if((length(firstP)!=0)&(length(lastP)!=0)){ 
      # Range was likely found
      foundR <- TRUE
      newSeq <- (firstP:lastP)
      #tmpSeqs <- c(splitsFound[firstP[[".id"]]],splitsFound[lastP[[".id"]]])
    } else {
      if(length(firstP)==0){
        # Guess firstP based on average expected pagesize
        firstP  <- locatePages(min(unlist(tmpSeqs)),prepost,ITbit,norm=1)
        extra   <- round((firstP$loc-round((minP/2)*nchar(ITbit)))/avg_page)
        if(extra>0){
          newSeq  <- ((firstP$page-extra):max(unlist(tmpSeqs)))
          tmpSeqs <- list((min(unlist(tmpSeqs))-extra):(min(unlist(tmpSeqs))-1),as.vector(unlist(tmpSeqs)[unlist(tmpSeqs)%in%newSeq]))
        } else {
          tmpSeqs <- list(NA,unlist(tmpSeqs))
        }
      }
      if(length(lastP)==0){
        # Guess lastP
        lastP  <- locatePages(max(unlist(tmpSeqs),na.rm=T),prepost,ITbit,norm=1)
        extra  <- round((nchar(ITbit)-lastP$loc)/avg_page)
        if(extra>0){
          newSeq  <- (min(unlist(tmpSeqs),na.rm=T)):(max(unlist(tmpSeqs),na.rm=T)+extra)
          if(is.na(tmpSeqs[[1]])){
            tmpSeqs <- list(tmpSeqs[[1]],as.vector(unlist(tmpSeqs[-1])[unlist(tmpSeqs[-1])%in%newSeq]),(max(unlist(tmpSeqs[-1])+1,na.rm=T):(max(unlist(tmpSeqs[-1]))+extra))) 
          } else {
            tmpSeqs <- list((min(unlist(tmpSeqs))-extra):(min(unlist(tmpSeqs))-1),as.vector(unlist(tmpSeqs)[unlist(tmpSeqs)%in%newSeq]))
          }  
        } else {
          tmpSeqs <- list(tmpSeqs[[1]],unlist(tmpSeqs[-1]),NA)
        }
      } 
    }
    
    attr(tmpSeqs,"Close")   <- (f_l(unique(unlist(tmpSeqs)))%in%NOcigar)
    attr(tmpSeqs,"SeqMiss") <- diff(newSeq[(newSeq%in%unique(unlist(tmpSeqs)))])
    
    if(attr(tmpSeqs,"Close")){foundD=TRUE} 
    if((length(attr(tmpSeqs,"SeqMiss"))>0)&(all(attr(tmpSeqs,"SeqMiss")<=closeFITmiss))){foundR=TRUE}
    
    if(foundD&foundR){
      foundS=TRUE
      break
    } else {
      # TRY ALTERNATIVE STRATEGIES
      # Store current
      vault <- list(newSeq,tmpSeqs)
      
      # 1st repeat range check with new sequences
      ranges      <- ldply(tmpSeqs,f_m)
      d           <- (ranges$max-ranges$min)
      d[d==0]     <- 1
      tmpSeqs     <- tmpSeqs[(d<=ppMX)&(d>=ppMN)]
      locSplits <- locatePages(unlist(tmpSeqs),prepost,ITbit,norm=nchar(ITbit))
      tmpSeqs   <- rangeCheck(tmpSeqs,locSplits,nPages,closeFITdiff)
      if(length(tmpSeqs)>0){
        locSplits <- locatePages(tmpSeqs,prepost,ITbit,norm=nchar(ITbit))
      } else {
        tmpSeqs <- vault[[2]]
      }
      newSeq    <- min(unlist(tmpSeqs)):max(unlist(tmpSeqs))
    }
    # Note results in attributes of tmpSeqs for evaluation
    attr(tmpSeqs,"Close")   <- (f_l(unique(unlist(tmpSeqs)))%in%NOcigar)
    attr(tmpSeqs,"SeqMiss") <- diff(newSeq[(newSeq%in%unique(unlist(tmpSeqs)))])
    
    # Evaluate
    if(attr(tmpSeqs,"Close")){foundD=TRUE} 
    if((length(attr(tmpSeqs,"SeqMiss"))>0)&(all(attr(tmpSeqs,"SeqMiss")<=closeFITmiss))){foundR=TRUE} 
    
    if(foundD&foundR){
      foundS=TRUE
    } else {
      # Exit the while loop
      foundD <- TRUE
      foundR <- TRUE
      foundS <- TRUE
      break
    }
    
  } # While any
  
  # Fill in the blanks
  ranges   <- ldply(tmpSeqs,f_m)
  d        <- (ranges$max-ranges$min)
  d[d==0]  <- 1
  tmpSeqs  <- tmpSeqs[(d<=ppMX)&(d>=ppMN)]
  newSeq   <- min(unlist(tmpSeqs)):max(unlist(tmpSeqs))
  pagelocs <- locatePages(newSeq,prepost,ITbit,norm=1)
  
  if(pagelocs$loc[length(pagelocs$loc)]==-1){
    pagelocs$loc[length(pagelocs$loc)]<-lst_page
  }
  if(pagelocs$loc[1]==-1){
    pagelocs$loc[1]<-fst_page
  }
  
  pagelocs$loc[pagelocs$loc<=0] <- NA
  pagelocs$loc[which(is.na(pagelocs$loc))] <- round(approx(x=pagelocs$page,y=pagelocs$loc,xout=pagelocs$page[which(is.na(pagelocs$loc))])$y)
  
  pagelocs$keep <- rep(TRUE,times=length(pagelocs$page))
  ids <- splitSeq(which(diff(sort(as.numeric(gsub("\\D","0",pagelocs$page))))>closeFITmiss))
  if(length(ids)>=1){
    if(length(ids)==1){ids <- unlist(ids)
    } else {
      ids  <- c(ids[[1]],unlist(ids[-1])+1)
    }
    pagelocs$keep[ids] <- FALSE
    pagelocs <- subset(pagelocs,keep)
  }
  
  # Create a join-ready structure
  pageSkeleton <- data.frame(location=pagelocs$loc,pdf.id=PDFid,name="",element="page",header=pagelocs$page)
  pageSkeleton <- arrange(pageSkeleton,pdf.id,location)
  return(pageSkeleton)
} # getPage


# ESCI functions ----------------------------------------------------------

# Conversion formula's from Friedman(1982) and Wolf(1986)
# Also see http://www.soph.uab.edu/Statgenetics/People/MBeasley/Courses/EffectSizeConversion.pdf
f_d <- function(f,df1,df2){ifelse(df1>1,{2*sqrt(df1*f/df2)},{2*sqrt(f/df2)})}
f_r <- function(f,df1,df2){ifelse(df1>1,{sqrt((df1*f)/((df1*f) + df2))},{sqrt(f/(f + df2))})}
t_d <- function(t,df,n1=1,n2=1){((n1+n2)*t)/(sqrt(df)*sqrt(n1*n2))}
t_r <- function(t,df){sqrt(t^2/(t^2 + df))}
z_d <- function(z,N){(2*z)/sqrt(N)}
z_r <- function(z,N){sqrt(z^2/(z^2 + N))}
X_d <- function(X,N,df=1){ifelse(df>1,{2*sqrt(X/N)},{2*sqrt(X/(N-X))})}
X_r <- function(X,N,df=1){ifelse(df>1,{sqrt(X/(X+N))},{sqrt(X/N)})}
r_d <- function(r){sqrt((4*(r^2))/(1-r^2))}  
d_r <- function(d){sqrt(d^2/(4+d^2))}  


del2lam <- function(delta,N){# Change an observed t-value (delta) to a standardised mean (difference)
  if(length(N)>1){
    n1 <- N[1]
    n2 <- N[2]
  } else {
    n1 <- floor(N/2)
    n2 <- ceiling(N/2)
    if(N%%2>0){n1<-n1+1}
  }
  return(delta * sqrt((n1 * n2)/(n1 + n2)))
}

lam2del <- function(lambda,N){# Change an observed standardised mean (lambda) to a t-value (delta)
  if(length(N)>1){
    n1 <- N[1]
    n2 <- N[2]
  } else {
    n1 <- floor(N/2)
    n2 <- ceiling(N/2)
    if(N%%2>0){n1<-n1+1}
  }
  return(lambda * sqrt((n1 + n2)/(n1 * n2)))
}
# 
# p2del <- function(p,N=NULL, df, prediction="not equal"){
#   if(!is.null(N)){
#   if(length(N)>1){
#     df <- sum(N)-2
#   } else {
#     df <- N-1
#   }} 
#   if(prediction=="not equal") p<-p/2 
#   return(delta = qt(p, df))
# }

sev.t <- function(t.obs,t.df,t.crit,discrepancy){
  require(MBESS)
  
  l  <- sapply(list(t.obs,discrepancy),length)
  if(l[2]!=1&l[2]<=10){
    cat("\nExpected a discrepancy vector of length == 1, or length >= 2\nWill use a 100 step sequence representing 95%CI for: t(",t.df,") = ",t.obs,"\n")
    CI <- conf.limits.nct(ncp=t.obs, df=t.df, conf.level=.95)
    discrepancy <- seq(1.5*CI$Lower.Limit,1.5*CI$Upper.Limit,length=100)
    l[2] <- 100
  }
  if(l[1]<l[2]){
    t.obs  <- rep(t.obs,length(discrepancy))
  }
  return(pt(t.obs,Ori.df,ncp=discrepancy,lower.tail=(t.crit<0)))
}


php.t = function(P, df, prediction="not equal", alpha=.05) {# 'p-value based' post-hoc power
  # Lenth, R. (2007). Post hoc power: tables and commentary. http://www.stat.uiowa.edu/files/stat/techrep/tr378.pdf
  ifelse(prediction=="not equal",two.tailed=TRUE,two.tailed=FALSE)
  
  if (two.tailed) {
    delta = qt(1 - P/2, df)
    cv = qt(1 - alpha/2, df)
    power = 1 - pt(cv, df, delta) + pt(-cv, df, delta)
  }
  else {
    delta = qt(1 - P, df)
    cv = qt(1 - alpha, df)
    power = 1 - pt(cv, df, delta)
  }
  power
}


posthocPOWer <- function(stat.type,inference,stat.df){
  # Check input arguments
  stat.type<-tolower(stat.type)
  if(stat.type%in%tolower(c("X^2","X2","chi.sq"))){stat.type<-"chisq"}
  
  if(length(stat.type)!=1){
    stop('\nChoose 1 distribution out of: ',paste(stat.type,collapse="  "))
  } else {
    if(!any(tolower(c("z","t","f","chisq","chi.sq","X^2","X2"))%in%stat.type)){ 
      warning(paste("Unknown test statistic. Choose from: Z, t, F, chisq"),immediate.=T)
      posthocPOW<-as.data.frame(t(rep(NA,1)))
      colnames(posthocPOW) <- list("posthocPOW")
      return(posthocPOW)
    }
  }
  
  if(length(inference$prediction)>1){
    cat("\nArgument prediction > 1.\nAssuming undirected test...\n")
    inference$prediction <- "not equal"
  } else {
    inference$prediction <- tolower(inference$prediction)
    if(!any(c("not equal","greater","less")%in%inference$prediction)){
      stop("Unknown prediction, choose from: not equal, greater, less")
    }
  }
  
  if(stat.type!="f"){
    stat.df <- c(stat.df[1],NA)
  }
  
  if(stat.type=="t"){
  posthocPOW <- if(inference$prediction=="not equal"){
        (1 - pt(abs(inference$stat.crit), stat.df[1],ncp=inference$stat.ncp) + pt(-1*abs(inference$stat.crit), stat.df[1],ncp=inference$stat.ncp))
      } else {
        (1 - pt(abs(inference$stat.crit), stat.df[1]))
      }
  }
   if(stat.type=="f"){
      if(length(stat.df)==2){
        posthocPOW  <- (1-pf(abs(inference$stat.crit), stat.df[1], stat.df[2]))
      } else {
        stop("The F distribution requires 2 degrees of freedom.")
      }
   }
     if(stat.type=="chisq"){
        posthocPOW  <- (1 - pchisq(abs(inference$stat.crit), stat.df[1])) 
      }
  return(data.frame(posthocPOW=posthocPOW))
}

convertES <- function(stat.type,inference,stat.N,stat.df){
  
  # Check input arguments
  stat.type<-tolower(stat.type)
  if(stat.type%in%tolower(c("X^2","X2","chi.sq"))){stat.type<-"chisq"}
  
  if(length(stat.type)!=1){
    stop('\nChoose 1 distribution out of: ',paste(stat.type,collapse="  "))
  } else {
    if(!any(tolower(c("z","t","f","chisq","chi.sq","X^2","X2"))%in%stat.type)){ 
      warning(paste("Unknown test statistic. Choose from: Z, t, F, chisq"),immediate.=T)
      ESconvert<-as.data.frame(t(rep(NA,6)))
      colnames(ESconvert) <- list("ES.d","ES.d.ciL","ES.d.ciU","ES.r","ES.r.ciL","ES.r.ciU")
      return(ESconvert)
    }
  }
  
  if(length(inference$prediction)>1){
    cat("\nArgument prediction > 1.\nAssuming undirected test...\n")
    inference$prediction <- "not equal"
  } else {
    inference$prediction <- tolower(inference$prediction)
    if(!any(c("not equal","greater","less")%in%inference$prediction)){
      stop("Unknown prediction, choose from: not equal, greater, less")
    }
  }
  
  if(stat.type!="f"){
    stat.df <- c(stat.df[1],NA)
  }
  
  if(stat.type=="t"){
    ES.d     <- t_d(inference$stat.ncp, stat.df[1])
    ES.d.ciL <- t_d(inference$stat.ncp.ciL, stat.df[1])
    ES.d.ciU <- t_d(inference$stat.ncp.ciU, stat.df[1])
    ES.r     <- t_r(inference$stat.ncp, stat.df[1])
    ES.r.ciL <- t_r(inference$stat.ncp.ciL, stat.df[1])
    ES.r.ciU <- t_r(inference$stat.ncp.ciU, stat.df[1])
  }
  if(stat.type=="f"){
    ES.d     <- f_d(inference$stat.ncp, stat.df[1], stat.df[2])
    ES.d.ciL <- f_d(inference$stat.ncp.ciL, stat.df[1], stat.df[2])
    ES.d.ciU <- f_d(inference$stat.ncp.ciU, stat.df[1], stat.df[2])
    ES.r     <- f_r(inference$stat.ncp, stat.df[1], stat.df[2])
    ES.r.ciL <- f_r(inference$stat.ncp.ciL, stat.df[1], stat.df[2])
    ES.r.ciU <- f_r(inference$stat.ncp.ciU, stat.df[1], stat.df[2])
  }
  
  if(stat.type=="chisq"){
    ES.d     <- X_d(inference$stat.ncp, stat.N, stat.df[1])
    ES.d.ciL <- X_d(inference$stat.ncp.ciL,stat.N, stat.df[1])
    ES.d.ciU <- X_d(inference$stat.ncp.ciU,stat.N, stat.df[1])
    ES.r     <- X_r(inference$stat.ncp, stat.N, stat.df[1])
    ES.r.ciL <- X_r(inference$stat.ncp.ciL,stat.N, stat.df[1])
    ES.r.ciU <- X_r(inference$stat.ncp.ciU,stat.N, stat.df[1])
  }
  tmp<-list(ES.d,ES.d.ciL,ES.d.ciU,ES.r,ES.r.ciL,ES.r.ciU)
  tmp[sapply(tmp,length)==0]<-NA
  ESconvert<-as.data.frame(t(unlist(tmp)))
  colnames(ESconvert) <- list("ES.d","ES.d.ciL","ES.d.ciU","ES.r","ES.r.ciL","ES.r.ciU")

  return(ESconvert)
}

decideNP <- function(stat.type=c("z","t","f","chisq"), stat.ncp, stat.df, stat.N, alpha=0.05, CL=0.95, prediction="not equal"){
  
  # Check input arguments
  stat.type<-tolower(stat.type)
  if(stat.type%in%tolower(c("X^2","X2","chi.sq"))){stat.type<-"chisq"}
  
  if(length(stat.type)!=1){
    stop('\nChoose 1 distribution out of: ',paste(stat.type,collapse="  "))
  } else {
    if(!any(tolower(c("z","t","f","chisq","chi.sq","X^2","X2"))%in%stat.type)){ 
      warning(paste("Unknown test statistic. Choose from: Z, t, F, chisq"),immediate.=T)
      return(decide <- data.frame(stat.ncp=stat.ncp,
                                  stat.ncp.ciL=NA,
                                  stat.ncp.ciU=NA,
                                  ci.type="unknown",
                                  stat.crit=NA,
                                  stat.ncp.p=NA,
                                  prediction=prediction,
                                  stat.crit.p=alpha,
                                  H0=NA,
                                  H1=NA,
                                  decide="unknown"))
    }
  }
  
  if(length(prediction)>1){
    cat("\nArgument prediction > 1.\nAssuming undirected test...\n")
    prediction <- "not equal"
  } else {
    prediction <- tolower(prediction)
    if(!any(c("not equal","greater","less")%in%prediction)){
      stop("Unknown prediction, choose from: not equal, greater, less")
    }
  }
  
  if(stat.type!="f"){
    stat.df <- c(stat.df[1],NA)
  }
  
  ID    <- NULL
  alpha <- c(alpha,(1-(alpha)))
  
  if(stat.type=="t"){
    if(prediction=="not equal"){alpha=c(alpha[1]/2,(1-(alpha[1]/2)))}
    x.crit <- qt(alpha,stat.df[1])
    stat.p <- pt(stat.ncp,stat.df[1])
    CI <- conf.limits.nct(ncp=stat.ncp, df=stat.df[1], conf.level=CL)
    stat.ncp.ciL <- CI$Lower.Limit
    stat.ncp.ciU <- CI$Upper.Limit
    ci.type <- "symmetric"
  }
  if(stat.type=="f"){
    if(length(stat.df)==2){
      x.crit <- qf(alpha,df1=stat.df[1],df2=stat.df[2])
      stat.p <- pf(stat.ncp,df1=stat.df[1],df2=stat.df[2])
      if(stat.ncp>x.crit[2]){
        CI      <- conf.limits.ncf(F.value = stat.ncp, conf.level = CL, df.1 = stat.df[1], df.2 = stat.df[2])
        ci.type <- "symmetric"
      } else {
        CI     <- conf.limits.ncf(F.value = stat.ncp, conf.level = NULL, alpha.lower=0, alpha.upper=alpha[1], df.1 = stat.df[1], df.2 = stat.df[2])
        CI$Lower.Limit <- stat.ncp
        ci.type <- "asymmetric"
      }
      if(length(CI$Lower.Limit)==0){CI$Lower.Limit <- NA}
      stat.ncp.ciL <- CI$Lower.Limit
      stat.ncp.ciU <- CI$Upper.Limit
      prediction <- "greater"
    } else {
      stop("The F distribution requires 2 degrees of freedom.")
    }
  }
  if(stat.type=="chisq"){
    x.crit <- qchisq(alpha,df=stat.df[1])
    stat.p <- pchisq(stat.ncp,df=stat.df[1])
    if(stat.ncp>(x.crit[2]+1)){
      CI     <- conf.limits.nc.chisq(Chi.Square=stat.ncp, conf.level=CL, df=stat.df[1])  
      ci.type <- "symmetric"
    } else {
      CI     <- conf.limits.nc.chisq(Chi.Square=stat.ncp, conf.level=NULL, alpha.lower=0, alpha.upper=alpha[1], df=stat.df[1])
      CI$Lower.Limit <- stat.ncp
      ci.type <- "asymmetric"
    }
    if(length(CI$Lower.Limit)==0){CI$Lower.Limit <- NA}
    stat.ncp.ciL  <- CI$Lower.Limit
    stat.ncp.ciU  <- CI$Upper.Limit
    prediction <- "greater"
  }
  
  if(prediction=="not equal"){
    H0       <- (stat.ncp >=  x.crit[1] & stat.ncp <=  x.crit[2])
    H1       <- (stat.ncp <   x.crit[1] | stat.ncp >   x.crit[2])
    ifelse(stat.ncp < 0,{ID<-1},{ID<-2})
  }
  if(prediction=="greater"){
    H0     <- (stat.ncp <= x.crit[2])
    H1     <- (stat.ncp >  x.crit[2])
    ID <- 2
  }
  if(prediction=="lesser"){
    H0     <- (stat.ncp >= x.crit[1])
    H1     <- (stat.ncp <  x.crit[1])
    ID <- 1
  }
  
  ifelse(H0, say<-"Accept H0",say<-"Reject H0")
  decide <- data.frame(stat.ncp=stat.ncp,
                       stat.ncp.ciL=stat.ncp.ciL,
                       stat.ncp.ciU=stat.ncp.ciU,
                       ci.type=ci.type,
                       stat.crit=x.crit,
                       stat.ncp.p=stat.p,
                       prediction=prediction,
                       stat.crit.p=alpha,
                       H0=H0,
                       H1=H1,
                       decide=say)
  return(decide[ID, ])
}

sev.data <- function(RPPdata,study){
  SEV.ori<- with(RPPdata, sev.info(stat.type=stat.ori.type[[study]],
                                   stat.ncp=as.numeric(stat.ori.ncp[[study]]),
                                   stat.df=c(as.numeric(stat.ori.df1[[study]]), as.numeric(stat.ori.df2[[study]])),
                                   stat.N=as.numeric(stat.ori.N[[study]]), prediction=prediction.ori[[study]],
                                   mus=list(SEV.ori.rep=as.numeric(stat.rep.ncp[[study]])),
                                   compare=list(SEV.comp.df=c(as.numeric(stat.rep.df1[[study]]),as.numeric(stat.rep.df2[[study]]))))
  )
  SEV.rep <- with(RPPdata, sev.info(stat.type=stat.rep.type[[study]],
                                    stat.ncp=as.numeric(stat.rep.ncp[[study]]),
                                    stat.df=c(as.numeric(stat.rep.df1[[study]]), as.numeric(stat.rep.df2[[study]])), 
                                    stat.N=as.numeric(stat.rep.N[[study]]), 
                                    prediction=prediction.rep[[study]],
                                    mus=list(SEV.ori.rep=as.numeric(stat.ori.ncp[[study]])),
                                    compare=list(SEV.comp.df=c(as.numeric(stat.ori.df1[[study]]),as.numeric(stat.ori.df2[[study]]))))
  )
  
  return(list(SEV.ori,SEV.rep))
}

sev.info <- function(stat.type=c("Z","t","F","chisq"),stat.ncp, stat.df, stat.N ,CL=.95, alpha=.05, prediction=c("not equal","greater","less"),mus=list(SEV.0=0),compare=list(comp.df=1000)){
  require(MBESS)
  #   stat.ncp <- as.numeric(stat.ncp) 
  #   stat.df[[1]]  <- as.numeric(stat.df[1])
  #   stat.df[[2]] <- as.numeric(stat.df[2])
  #   stat.N <- as.numeric(stat.N)
  #   CL       <- as.numeric(CL)
  #   alpha    <- as.numeric(alpha)
  
  stat.type<-tolower(stat.type)
  if(stat.type%in%tolower(c("X^2","X2","chi.sq"))){stat.type<-"chisq"}
  
  if(length(stat.type)!=1){
    stop('\nChoose 1 distribution out of: ',paste(stat.type,collapse="  "))
  } else {
    prediction <- tolower(prediction)
    if(!any(tolower(c("z","t","f","chisq","chi.sq","X^2","X2"))%in%stat.type)){
      stop("Unknown test statistic, choose from: Z, t, F, chisq")
    }
    
    if(length(prediction)>1){
      cat("\nArgument prediction > 1.\nAssuming undirected test...\n")
      prediction <- "not equal"
    } else {
      if(!any(c("not equal","greater","less")%in%prediction)){
        stop("Unknown prediction, choose from: not equal, greater, less")
      }
    }
    
    if(stat.type!="f"){
      stat.df <- c(stat.df[1],NA)
    }
    
    infer      <- decideNP(stat.type=stat.type,stat.ncp=stat.ncp,stat.df=stat.df,alpha=alpha,prediction=prediction)
    infer.comp <- decideNP(stat.type=stat.type,stat.ncp=mus[[1]],stat.df=compare[[1]],alpha=alpha,prediction=prediction)
    
    if(is.na(infer$stat.ncp.ciL)){infer$stat.ncp.ciL<-0}
    if(is.na(infer$stat.ncp.ciU)){infer$stat.ncp.ciU<-0}
    x     <- seq(infer$stat.ncp.ciL,infer$stat.ncp.ciU,by=.01)
    x     <- sort(c(x,infer$stat.ncp.ciL,infer$stat.ncp.ciU))
    
    #     if(stat.type=="z"){
    #       CI    <- c(Lower.Limit=-1*infer$stat.crit,Upper.Limit=infer$stat.crit)
    #       x     <- seq(CI$Lower.Limit*1.5,CI$Upper.Limit*1.5,by=.01)
    #       x     <- sort(c(x,stat.ncp,CI$Lower.Limit,CI$Upper.Limit))
    #       y.ncl <- dt(x, dnorm=(ncp=CI$Lower.Limit)
    #       y.ncu <- dt(x, dnorm=(stat.df[1]),ncp=CI$Upper.Limit)
    #       y.ncp <- dt(x, df=(stat.df[1]),ncp=stat.ncp)
    #       POW   <- ifelse(prediction=="not equal",{(1 - pt(abs(infer$stat.crit), stat.df[1], x) + pt(-1*abs(infer$stat.crit), stat.df[1], x))},{(1 - pt(abs(infer$stat.crit), stat.df[1], x))})
    #       SEV.crit <- pt(stat.ncp,stat.df[1],ncp=infer$stat.crit,lower.tail=(infer$stat.crit<0))
    #       SEV.obs  <- pt(stat.ncp,stat.df[1],ncp=stat.ncp, lower.tail=(infer$stat.crit<0))  
    #     }
    
    if(stat.type=="t"){
      x.d   <- t_d(x,stat.df[1])
      x.r   <- t_r(x,stat.df[1])
      y.ncl <- dt(x, df=(stat.df[1]),ncp=infer$stat.ncp.ciL)
      y.ncu <- dt(x, df=(stat.df[1]),ncp=infer$stat.ncp.ciU)
      y.ncp <- dt(x, df=(stat.df[1]),ncp=stat.ncp)
      y.SEV.obs  <- pt(rep(stat.ncp,length(x)),stat.df[1],ncp=x,lower.tail=infer$H1)
      y.POW.post <- if(prediction=="not equal"){
        1 - pt(rep(abs(infer$stat.crit),length(x)), stat.df[1], ncp=x) + pt(rep(-1*abs(infer$stat.crit),length(x)), stat.df[1], ncp=x)
      } else {
        (1 - pt(rep(abs(infer$stat.crit),length(x)), stat.df[1], ncp=x))
      }
      POW.post <- if(prediction=="not equal"){
        (1 - pt(abs(infer$stat.crit), stat.df[1],ncp=stat.ncp) + pt(-1*abs(infer$stat.crit), stat.df[1],ncp=stat.ncp))
      } else {
        (1 - pt(abs(infer$stat.crit), stat.df[1]))
      }
      
      SEV.x   <- c(SEV.obs=stat.ncp,SEV.ciL=infer$stat.ncp.ciL,SEV.ciU=infer$stat.ncp.ciU,SEV.crit=infer$stat.crit,SEV.50= x[which.min(y.SEV.obs>=.5)],unlist(mus))
      SEV.d   <- t_d(SEV.x, stat.df[1])
      SEV.r   <- t_r(SEV.x, stat.df[1])
      SEV.y   <- sapply(SEV.x, function(mu) pt(stat.ncp, stat.df[1], ncp=mu, lower.tail=infer$H1))
      if(length(compare)==1){
        SEV.comp.x  <- qt(p=pt(SEV.x,stat.df[1],lower.tail=infer$H1), compare[[1]][1], lower.tail=infer.comp$H1)
        SEV.comp.y  <- sapply(SEV.comp.x,function(mu) pt(stat.ncp, compare[[1]][1], ncp=mu, lower.tail=infer.comp$H1))
      }
    }
    
    if(stat.type=="f"){
      prediction <- "greater"
      ifelse(length(stat.df)==2,{
        x.d   <- f_d(x,stat.df[1],stat.df[2])
        x.r   <- f_r(x,stat.df[1],stat.df[2])
        y.ncl  <- df(x, df1=(stat.df[1]), df2=(stat.df[2]), ncp=infer$stat.ncp.ciL)
        y.ncu  <- df(x, df1=(stat.df[1]), df2=(stat.df[2]), ncp=infer$stat.ncp.ciU)
        y.ncp  <- df(x, df1=(stat.df[1]), df2=(stat.df[2]), ncp=stat.ncp)
        y.SEV.obs  <- pf(rep(stat.ncp,length(x)),stat.df[1],stat.df[2],ncp=x,lower.tail=infer$H1)
        y.POW.post <- (1 - pf(rep(abs(infer$stat.crit),length(x)), stat.df[1], stat.df[2],ncp=x))
        POW.post   <- (1 - pf(abs(infer$stat.crit), stat.df[1], stat.df[2]))
        
        SEV.x   <- c(SEV.obs=stat.ncp,SEV.ciL=infer$stat.ncp.ciL,SEV.ciU=infer$stat.ncp.ciU,SEV.crit=infer$stat.crit,SEV.50= x[which.min(y.SEV.obs>=.5)],unlist(mus))
        SEV.d   <- f_d(SEV.x,stat.df[1],stat.df[2])
        SEV.r   <- f_r(SEV.x,stat.df[1],stat.df[2])        
        SEV.y   <- sapply(SEV.x,function(mu) pf(stat.ncp,stat.df[1],stat.df[2],ncp=mu, lower.tail=infer$H1))
        
        if(length(compare)==1){
          SEV.comp.x  <- qf(p=pf(SEV.x,stat.df[1],stat.df[2],lower.tail=infer$H1),df1=compare[[1]][1],df2=compare[[1]][2], lower.tail=infer.comp$H1)
          SEV.comp.y  <- sapply(SEV.comp.x,function(mu) pf(stat.ncp,df1=compare[[1]][1],df2=compare[[1]][2],ncp=mu, lower.tail=infer.comp$H1))
        }
      },{
        stop("The F distribution requires 2 degrees of freedom.")
      })
    }
    
    if(stat.type=="chisq"){
      prediction <- "greater" 
      x.d   <- X_d(x,stat.N,stat.df[1])
      x.r   <- X_r(x,stat.N,stat.df[1])
      y.ncl <- dchisq(x, df=(stat.df[1]),ncp=infer$stat.ncp.ciL)
      y.ncu <- dchisq(x, df=(stat.df[1]),ncp=infer$stat.ncp.ciU)
      y.ncp <- dchisq(x, df=(stat.df[1]),ncp=stat.ncp)
      y.SEV.obs  <- pchisq(rep(stat.ncp,length(x)),stat.df[1],ncp=x,lower.tail=infer$H1)
      
      y.POW.post <- (1 - pchisq(rep(abs(infer$stat.crit),length(x)), stat.df[1],ncp=x))
      POW.post   <- (1 - pchisq(abs(infer$stat.crit), stat.df[1], x))
      
      SEV.x   <- c(SEV.obs=stat.ncp,SEV.ciL=infer$stat.ncp.ciL,SEV.ciU=infer$stat.ncp.ciU,SEV.crit=infer$stat.crit,SEV.50= x[which.min(y.SEV.obs>=.5)],unlist(mus))
      SEV.d   <-  X_d(x,stat.N,stat.df[1])
      SEV.r   <-  X_r(x,stat.N,stat.df[1])
      SEV.y   <- sapply(SEV.x,function(mu) pchisq(stat.ncp,stat.df[1],ncp=mu, lower.tail=infer$H1))
      if(length(compare)==1){
        SEV.comp.x  <- qchisq(p=pchisq(SEV.x,stat.df[1],lower.tail=infer$H1),compare[[1]][1],lower.tail=infer.comp$H1)
        SEV.comp.y  <- sapply(SEV.comp.x,function(mu) pchisq(stat.ncp,compare[[1]][1],ncp=mu,  lower.tail=infer.comp$H1))
      }
    }
    
  }
  
  return(list(inference= data.frame(stat=stat.type,df1=stat.df[1],df2=stat.df[2],alpha=alpha,stat.CI=CL,infer,POW.post=POW.post), 
              severity = data.frame(cbind(SEV.x,SEV.d,SEV.r,SEV.y,SEV.comp.x,SEV.comp.y)),
              curves   = data.frame(stat.x=x,stat.d=x.d,stat.r=x.r,y.ncl=y.ncl,y.ncu=y.ncu,y.ncp=y.ncp,y.POW.post=y.POW.post,y.SEV.obs=y.SEV.obs)))
}

plotReplication <- function(SEV.ori,SEV.rep, d.axis=c("stat","d","r"), studyname="Study 1", pl.sevlabels=TRUE, pl.sevlabels.comp=FALSE, pl.power=TRUE, pl.crit=TRUE, pl.labels=TRUE, pl.connect=FALSE, pl.ci=FALSE){
  require(ggplot2)
  
  if(length(d.axis)!=1){
    cat("\nArgument discrepancy axis not 1.\nAssuming original test statistic test...\n")
    d.axis <- "stat"
  } else {
    if(!any(c("stat","d","r")%in%d.axis)){
      stop("Unknown discrepancy axis, choose from: stat, d, r")
    }
  }
  
  ori.dat <- SEV.ori[["curves"]]
  ori.t   <- SEV.ori[["inference"]]
  ori.sev <- SEV.ori[["severity"]]
  
  rep.dat <- SEV.rep[["curves"]]
  rep.t   <- SEV.rep[["inference"]]
  rep.sev <- SEV.rep[["severity"]]
  
  switch(d.axis,
         "stat" = dID <- 1,
         "d"    = dID <- 2,
         "r"    = dID <- 3
  )
  
  
  df <- data.frame(rbind(ori=ori.dat,rep=rep.dat),study=c(rep(paste("Original:",ori.t$decide),times=nrow(ori.dat)),rep(paste("Replication:",rep.t$decide),nrow(rep.dat))),power=c(rep(paste("Original: POW =",round(ori.t$POW.post,digits=2)),times=nrow(ori.dat)),rep(paste("Replication: POW =",round(rep.t$POW.post,digits=2)),nrow(rep.dat))),row.names=NULL,stringsAsFactors=F)
  
  test <- c(paste0("Original d(x0) = ",round(ori.sev[1,dID],digits=2),"\nCrit. ",ori.t$stat,"(",ori.t$df1,if(!is.na(ori.t$df2)){paste0(", ",ori.t$df2)},") = ",round(ori.t$stat.crit,digits=2)),paste0("Replication d(x0) = ",round(rep.sev[1,dID],digits=2),"\nCrit. ",rep.t$stat,"(",rep.t$df1,if(!is.na(rep.t$df2)){paste0(", ",rep.t$df2)},") = ",round(rep.t$stat.crit,digits=2)))
  
  df.sevr <- data.frame(rbind(ori=ori.sev,rep=rep.sev),study=c(rep(paste("Original:",ori.t$decide),times=nrow(ori.sev)),rep(paste("Replication:",rep.t$decide),nrow(rep.sev))),yend=c(ori.sev$SEV.y,rep.sev$SEV.y),yend.comp=c(ori.sev$SEV.comp.y,rep.sev$SEV.comp.y) )
  xmin   <- min(df[,1])
  
  if(pl.sevlabels){
    df.sev <- df.sevr[c(1,6,7,12,4,10), ]
    mlt    <- grep("Original",df.sev$study,fixed=T)
    xmin   <- min(df[,1])-((max(df[,1])-min(df[,1]))/5)
    df.sev$xend <- c(xmin,(xmin+min(df[,1]))/2,xmin,(xmin+min(df[,1]))/2,xmin,(xmin+min(df[,1]))/2)[rank(df.sev$yend)]
    col <- 1:6
    col[ mlt] <- rep("grey50",length(mlt))
    col[-mlt] <- rep("black",length(mlt))
    #     df.sev$xend        <- df.sev[,1]
    #     df.sev$xend[ mlt]  <- xmin
    #     df.sev$xend[-mlt]  <- (xmin+min(df[,1]))/2
    df.sev$dfrom       <- c(test,rev(test),test)
    pre <- character(nrow(df.sev))
    ifelse(ori.t$H1,{pre[ mlt]<-"SEV(mu>"},{pre[ mlt]<-"SEV(mu<="})
    ifelse(rep.t$H1,{pre[-mlt]<-"SEV(mu>"},{pre[-mlt]<-"SEV(mu<="})
    
    #labels <- paste0(pre,round(df.sev[ ,dID],digits=2),")==",round(df.sev$SEV.y,digits=2))
    labels <- paste0(pre,round(df.sev[ ,dID],digits=2),")")
  }
  
  if(pl.sevlabels.comp){
    df.sev.comp <- df.sevr[c(6,12), ]
    mlt    <- grep("Original",df.sev.comp$study,fixed=T)
    dst    <- abs(mean(c(ori.sev[,1]-rep.sev[,1])))*1.2
    df.sev.comp$xend        <- df.sev.comp$SEV.comp.x
    df.sev.comp$xend[mlt]   <- df.sev.comp$xend[mlt]+(-1.5*dst)
    df.sev.comp$xend[-mlt]  <- df.sev.comp$xend[-mlt]+(1.5*dst)
    df.sev.comp$dfrom       <- rev(test)
  }  
  
  if(pl.crit){
    df.crit<- df.sevr[c(4,10), ]
    df.crit$dfrom <- test 
  }
  
  # Discrepance axis will be labelled according to d.axis
  df$disc     <- df[ ,1]
  df.sev$disc <- df.sev[ ,1]
  df.sev.comp$disc <- df.sev.comp$SEV.comp.x
  df.crit$disc<- df.crit[ ,1]
  d.label <- paste0("Discrepancy in units of ",list(ori.t$stat,"Cohen's d","Effect Size r")[[dID]])
  verdict <- ifelse(rep.t$H0,paste0("ReplicationInference: mu <= mu[1]"),paste0("ReplicationInference: mu > mu[1]"))
  
  repP <-  ggplot(df,aes(group=study)) + 
    geom_line(data=df,aes(x=disc,y=y.SEV.obs,color=study),size=2) + 
    geom_point(data=df.sev,aes(x=disc, y=SEV.y,shape=dfrom,fill=dfrom),alpha=.5,size=5)
  
  if(pl.sevlabels){repP <- repP + 
                     annotate("text",  x = df.sev$xend+.05, y = df.sev$yend, label = labels, size=2, parse=T, colour=col)}
  #annotate("text",  x = c(xmin,(xmin+min(df[,1]))/2), y = c(1,1), label = c("Original","Replication"), cex=4,fontface=2)}
  #==",round(df.sev.comp$SEV.comp.y,digits=2)
  if(pl.sevlabels.comp){repP <- repP + 
                          annotate("text",  x = df.sev.comp$xend, y = df.sev.comp$SEV.comp.y, label = paste0("SEV(mu>",round(df.sev.comp$SEV.comp.x,digits=2),")"),parse=T,cex=3) +
                          geom_point(data=df.sev.comp,aes(x=disc, y=SEV.comp.y,shape=dfrom,fill=dfrom),alpha=.5,size=5) +
                          geom_segment(data=df.sev.comp[1:4, ],aes(x=SEV.x, y=SEV.y, xend=SEV.comp.x, yend=SEV.comp.y),color="grey80")}
  
  if(pl.power){repP <- repP + geom_line(data=df,aes(x=disc,y=y.POW.post,group=power,colour=study,linetype=power))}
  
  if(pl.connect){repP <- repP + geom_segment(data=df.sev[1:4,],aes(x=disc, y=SEV.y, xend=rev(disc), yend=rev(SEV.y)),color="blue",alpha=.5)}
  
  if(pl.crit){repP <- repP + geom_point(data=df.crit,aes(x=disc, y=SEV.y,shape=dfrom),color="red",fill="red",alpha=.6,size=2)}
  
  if(pl.ci){repP <- repP + geom_point(data=df.crit,aes(x=disc, y=SEV.y,shape=dfrom),color="red",fill="red",alpha=.6,size=2)}
  
  if(pl.labels){repP <- repP +
                  ggtitle(paste(studyname)) + 
                  annotate("text",x=(max(df$disc)-xmin)/3,xmin=xmin,xmax=max(df$disc), y = 1.1, label = verdict, size=4,parse=T) +
                  ylab("POWER / SEVERITY") + 
                  xlab(paste0(d.label))}
  #seq(round(min(df$disc),digits=1),round(max(df$disc)
  repP <- repP +  
    scale_x_continuous(breaks=round(c(min(df$disc),df.sev$disc,max(df$disc)),digits=2),labels= round(c(min(df$disc),df.sev[ ,dID],max(df$disc)),digits=2),limits=c(xmin,max(df$disc))) +
    scale_y_continuous(breaks=c(0,df.sev$SEV.y,1),labels= round(c(0,df.sev$SEV.y,1),digits=2),limits=c(0,1.1)) +
    scale_color_manual(values=c("grey","grey30"),guide=guide_legend("Study: N-P decision")) +
    scale_shape_manual(values=c(21:25),guide=guide_legend(expression(paste("Evaluate ",d(x[0]),":",mu[1] == (mu[0] + gamma))))) +
    scale_linetype_manual(values=c(2,2),guide=guide_legend("Post-hoc power")) +
    scale_size_manual(values=c(1,2),guide=F) +
    scale_fill_manual(values=c("grey","black"),guide=guide_legend(expression(paste("Evaluate ",d(x[0]),":",mu[1] == (mu[0] + gamma))))) +
    theme_bw(base_size = 10, base_family = "")
  
  return(repP)
}


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

wordXX   <- function(x,clps="|",pre="(",post=")") {# ( X1|X2 ) Example: wordXX(c("yes","no","maybe"))
  paste(pre, paste(x,collapse=clps), post, sep="")
}

# wordXXp  <- function(x,clps=")|(",pre="\\b((",post="))\\b") {# Same as wordXX + escape puntuation
#   x <- gsub("(?=[[:punct:]])","\\",x,perl=T)
#   paste(pre, paste(x,collapse=clps), post, sep="")
# }

wordXnY   <- function(x,y,clps="|",pre="(",post=")",mid="[[:print:]]",n=1,m=1){# ( XanyY ) Example: wordXnY("yes","no")
  paste(pre, paste(x,collapse=clps),")",mid,"{",n,",",m,"}(", paste(y,collapse=clps), post, sep="")
}

wordXnYr  <- function(x,y,clps="|",pre="(",post=")",mid="[[:print:]]",n=1,m=1){# ( XanyY | YanyX ) Example: wordXnYr("yes","no")
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

wordXorY  <- function(x,y,clps="|",pre="\\b",post="\\b",mid=""){# ( XanyY? | Y?anyX ) Example: wordXnYr("yes","no")
  paste(paste(pre,paste(x,collapse=clps),mid,paste(y,collapse=clps),post,sep=""),
        "|",paste(pre,paste(y,collapse=clps),mid,paste(x,collapse=clps),post,sep=""), sep="")
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


PwP   <- function(pre="(\\w+[[:graph:]]*)?", word=NULL, post="(\\w+[[:graph:]]*)?"){
  # ( pre Word post ) Example: 
  # txt="Love thy neigbours"  
  # matchIT(PwP("\\s","thy","\\d"),txt)
  ifelse(all(is.null(pre),is.null(word),is.null(post)),
         print("PwordP: One or more arguments missing!"),
{
  #word<-gsub("[[:punct:]]+","",word)
  regstr <-paste(c(paste(pre,collapse="[[:blank:]]?"),
                   paste(c("",paste(c("(",paste(word,collapse="|"),")"),collapse=""),""),collapse="[[:blank:]]?"),
                   paste(post,collapse="[[:blank:]]?")),collapse="")
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


# TRY … CATCH -------------------------------------------------------------------------------------------------------------------------

##================================================================##
###  In longer simulations, aka computer experiments,            ###
###  you may want to                                             ###
###  1) catch all errors and warnings (and continue)             ###
###  2) store the error or warning messages                      ###
###                                                              ###
###  Here's a solution  (see R-help mailing list, Dec 9, 2010):  ###
##================================================================##

# Catch *and* save both errors and warnings, and in the case of
# a warning, also keep the computed result.
#
# @title tryCatch both warnings (with value) and errors
# @param expr an \R expression to evaluate
# @return a list with 'value' and 'warning', where
#   'value' may be an error caught.
# @author Martin Maechler;
# Copyright (C) 2010-2012  The R Core Team
# 
try.CATCH <- function(expr){
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                   warning = w.handler),
       warning = W)
}
