# Extraction functions ----------------------------------------------------

####
# Based on the function "getPDF()" in the package "statscheck" (https://github.com/MicheleNuijten/statcheck)
####

# Inner function to read pdf:
getPDFs <- function(x)
{
  txtfiles <- character(length(x))
  for (i in 1:length(x))
  {
    system(paste('pdftotext -q -enc "UTF-8" "',x[i],'"',sep=""))
    if (file.exists(gsub("\\.pdf","\\.txt",x[i])))
    {
      fileName <- gsub("\\.pdf","\\.txt",x[i])
      txtfiles[i] <- readChar(fileName, file.info(fileName)$size)
    } else
    {
      warning(paste("Failure in file",x[i]))
      txtfiles[i] <- ""
    }
  }
  return(enc2utf8(txtfiles))
}

####
# Based on the function "checkPDFdir()" from package "statchek" (https://github.com/MicheleNuijten/statcheck)
####

# Function to check directory of PDFs:
returnPDFdir <- function(dir, ...){
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


loadPatterns <- function(pubs=TRUE,stat=TRUE,info=TRUE){
  #Load 
  
  # Regexpr patterns to find common sections of articles
  #           small   capital     alt.
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
  
  pubsPatterns <- list(
    # Ordered header indicators..."Experiment 1", "Study 4b", "Appendix B"
    hdr="\n+[^[:print:]]*(((Experiment(s*)|Stud(y|ies))\\s*\\d{1,3})|(Appendi(x|ces)\\s*\\w{1,2}))[[:print:]]{0,100}\\n{1,2}",
    # Common section / paragraph indicators..."Method", "Procedure", "Results"
    sec="\n+[^[:print:]]*((Abstract|Summary)|((General[[:blank:]]*)?Introduction)|(Method(s*))|(Procedure(s*))|(Participants|Subjects|Sample(s)*|Patients|Clients)|Material(s)*|Instrumen(s)*|(Replication(s)*[[:blank:]]*(Stud(y|ies))?)|Results|((Meta[-]?|Data|Exploratory|Confirmatory)?[[:blank:]]*Analysis)|(Simulation(s)*)|((General[[:blank:]]*)?Conclusion(s)*[[:blank:]]*(and|&)*[[:blank:]]*(Discussion)*)|((General|Result(s)*)[[:blank:]]*(and|&)?[[:blank:]]*)?(Discussion)[[:blank:]]*)[^[:print:]]*",
    # Tables + titles
    tab="\\n+[^[:print:]]*[[:blank:]]*(Table(s)*[[:blank:]]*\\w*\\d{1,3}\\w*[[:blank:]]*)[[:print:]]{0,150}[^[:print:]]+",
    # Figures + titles
    fig="\\n+[^[:print:]]*[[:blank:]]*((Figure(s)*|Box)[[:blank:]]*\\w*\\d{1,3}\\w*[[:blank:]]*)[[:print:]]{0,150}[^[:print:]]+",
    # Page numbers
    pag="(\\n+[^[:print:]]*[[:blank:]]*\\d+[[:blank:]]*[^[:print:]]*\\n+)"
  )
  # Regexpr patterns to find statistics
  statPatterns <- list(
    # Search patterns and extraction strategy for t, F, r, Z, Wald and Chi^2 
    # adapted from package "statcheck" (https://github.com/MicheleNuijten/statcheck) 
    t="t\\s?\\(\\s?\\d*\\.?\\d+\\s?\\)\\s?.?\\s?\\D{0,3}\\s?\\d*,?\\d*\\.?\\d+\\s?,\\s?(ns|p\\s?.?\\s?\\d?\\.?\\d+)",
    F="F\\s?\\(\\s?\\d*\\.?\\d+\\s?,\\s?\\d*\\.?\\d+\\s?\\)\\s?.?\\s?\\d*,?\\d*\\.?\\d+\\s?,\\s?(ns|p\\s?.?\\s?\\d?\\.\\d+)",
    r="[rρΡ]\\s?\\(\\s?\\d*\\.?\\d+\\s?\\)\\s?.?\\s?\\D{0,3}\\s?\\d*\\.?\\d+\\s?,\\s?(ns|p\\s?.?\\s?\\d?\\.\\d+)",
    Z="[^a-z]?(z|Z)\\s?.?\\s?\\D{0,3}\\s?\\d*,?\\d*\\.?\\d+\\s?,\\s?(ns|p\\s?.?\\s?\\d?\\.\\d+)",
    W="[^a-z]?([wW]ald)\\s?\\D?\\s?\\D{0,3}\\s?\\d*,?\\d*\\.?\\d+\\s?,\\s?(ns|p\\s?.?\\s?\\d?\\.\\d+)",
    x2="[χΧx][[:blank:]]*[2²][[:blank:]]*\\(\\s?\\d*\\.?\\d+\\s?\\)\\s?.?\\s?\\s?\\d*\\,?\\d*\\.?\\d+\\s?,\\s?(ns|p\\s?.?\\s?\\d?\\.\\d+)",
    # Effect sizes
    g2="[(,]?[[:blank:]]*((partial)?[[:blank:]]*[ηΗg][[:blank:]]*[p]?[2²]?)[[:blank:]]*[=]?[[:blank:]]*(\\d*\\.+\\d+)[[:blank:]]*[)]?[[:blank:]]* ",
    d="[(,]?[[:blank:]]*([cC]ohen[[:blank:]]*(['`]s)?)?((\\s|[(])+[[:blank:]]*[,]?[[:blank:]]*[δΔd][[:blank:]]*[=]?[[:blank:]]*)+(\\d?\\.\\d+)+[[:blank:]]*[)]?",
    o2="[ωΩ][2²]?",
    phi="[φΦ]",
    # Sample info
    sd="[(,]?[[:blank:]]*(([σΣ])|([sS][.]?[dD][.]?))[[:blank:]]*[=]?[[:blank:]]*(\\d*[.,]?\\d*[.]?\\d+)+[[:blank:]]*[)]?",
    M="[(,]?[[:blank:]]*(([μΜ])|([mM][.]?[nN]?[.]?)|([mM]ean|[aA]verage|[mM]edian))[[:blank:]]*[=]?[[:blank:]]*(\\d*[.,]?\\d*[.]?\\d+)+[[:blank:]]*[)]?",
    N=c("[(,]?[[:blank:]]*(participant(s)*|undergraduate(s)*|graduate(s)*|subject(s)*|client(s)*|patient(s)*|student(s)*|adult(s)*|male(s)?|female(s)?|child(ren)*|kindergartner(s)*|preschooler(s)*|adolescent(s)*|toddler(s)*|infant(s)*|neonate(s)*|m(e|a)n|wom(e|a)n|girl(s)*|boy(s)*)[[:blank:]]*[,)]?[[:blank:]]*", "[[:blank:]]*[(,]?[[:blank:]]*[nN]?[[:blank:]]*[=]?[[:blank:]]*(\\d*[,]?\\d+)[[:blank:]]*[)]?[[:blank:]]*")
  )
  # Regexpr patterns to find information about analyses
  infoPatterns <- list(
    # Adjustments for multiple comparisons
    adj="[[:blank:]]*[(]?[[:blank:]]*(Bonferroni|Holm|Hochberg|Hommel|Benjamini|Yekutieli|Sidak|Shaffer|Simes|Familywise|FWE|LSD|(False\\s*Discovery)|FDR|(Multiple\\s*Comparisons)|((Adjusted|Corrected)\\s*P[-]?\\s*value(s)*))[[:blank:]]*[)]?[[:blank:]]*",
    # Direction / tailed / sided -ness of tests
    dir="[[:blank:]]*[(]?[[:blank:]]*((Un)*Directed|((\\d+|one|two)?\\s*[-]?\\s*sided)|((\\d+|one|two)?\\s*[-]?\\s*tailed))[[:blank:]]*[)]?[[:blank:]]*",
    dep="[[:blank:]]*[(]?[[:blank:]]*((((In)*Dependent|Paired|Matched)\\s*(Sample(s)*|Group(s)*|Pairs|Observation(s)*|Measurement(s)*|(\\w*[-]?test(s)*))+)|(Repeated\\s*Measure(s|ment(s)*)*)|((within|between)\\s*(subject(s)*|participant(s)*|factor(s)*|design|predictor(s)*|covariate(s)*|level(s)*|effect(s)*)))[[:blank:]]*[)]?[[:blank:]]*",
    rnd="[[:blank:]]*[(]?[[:blank:]]*((Assign(ed|ment)*|Select(ed|ion)*)?(At)*Random(ly)*\\s*(Assign(ed|ment)*|Select(ed|ion)*)?)[[:blank:]]*[)]?[[:blank:]]*"
  )
  return(list(pubsPatterns,statPatterns,infoPatterns)[c(pubs,stat,info)])
}


# Find labels such as study header and page number for each statistic 
getLabels <- function(raw_pos,hdr_pos,page_pos,tab_pos,txt,title){
  
  Ori <- names(raw_pos)
  page_label <- rep(NA,times=length(Ori))
  hdr_label  <- rep(NA,times=length(Ori))
  tab_label  <- rep(NA,times=length(Ori))
  
  if(attr(page_pos,"match.length")!=-1){
  # Get page numbers, discard any numbers that destroy sequential order
  Numbers  <- gsub("(^\\D*)|(\\D*$)","",substring(txt,page_pos,page_pos+attr(page_pos,"match.length")-1))
  Pages    <- Numbers[which(abs(diff(as.numeric(Numbers)))<=2)]
  while(any(as.numeric(Numbers)>max(as.numeric(Pages)))){
    Numbers  <- c(Pages,Numbers[as.numeric(Numbers)>max(as.numeric(Pages))])
    id       <- which(abs(diff(as.numeric(Numbers)))<=2)
    Pages    <- Numbers[c(id,max(id)+1)]
  }
  Numbers  <- gsub("(^\\D*)|(\\D*$)","",substring(txt,page_pos,page_pos+attr(page_pos,"match.length")-1))
  page_pos <- page_pos[which(Numbers%in%Pages)]
  
  # Decide on which page a statistic is printed
  page_diff  <- as.matrix(sapply(raw_pos,function(d) d-page_pos))
  page_label <- sapply(1:ncol(page_diff), function(d) ifelse(any(page_diff[,d]<0),{
    Pages[which(page_diff[,d]==(max(page_diff[page_diff[,d]<0,d])))]},{
      Pages[which(page_diff[,d]==(min(page_diff[page_diff[,d]>=0,d])))]
    }))
  }
  
  if(attr(hdr_pos,"match.length")!=-1){
  # Decide which Experiment or Study header belongs to the statistics
  Headers <- gsub("(^[^[:print:]]*)|([^[:print:]]*$)","",substring(txt,hdr_pos,hdr_pos+attr(hdr_pos,"match.length")-1))  
  # Decide which section a statistic belongs to
  hdr_diff  <- as.matrix(sapply(raw_pos,function(d) d-hdr_pos))
  hdr_label <- sapply(1:ncol(hdr_diff), function(d)  ifelse(any(hdr_diff[,d]>=0),{
    Headers[which(hdr_diff[,d]==min(hdr_diff[hdr_diff[,d]>=0,d]))]},{
      Headers[which(hdr_diff[,d]==min(hdr_diff[hdr_diff[,d]<0,d]))]
    }))
  }
  
  if(attr(tab_pos,"match.length")!=-1){
  # Decide if it is likely stats were extracted from a table
  Tables <- gsub("(^[^[:print:]]*)|([^[:print:]]*$)","",substring(txt,tab_pos,tab_pos+attr(tab_pos,"match.length")-1))  
  # Decide which table a statistic may have come from
  tab_diff  <- as.matrix(sapply(raw_pos,function(d) d-tab_pos))
  tab_label <- sapply(1:ncol(tab_diff), function(d)  ifelse(any(tab_diff[,d]>=0),{
      Tables[which(tab_diff[,d]==min(tab_diff[tab_diff[,d]>=0,d]))]},{
        rep(NA,times=length(Tables))
        }))
  }
 
  return(as.data.frame(cbind(title,Ori,page_label,hdr_label,tab_label)))
}


# REGEXPR tools -----------------------------------------------------------

# SEARCHERS
# A wrapper for grep, searches for string in source
regIT   <- function(str,src) {
  grep(str,src,ignore.case=TRUE,value=TRUE)
}

# Return matching positions in source using regmatches
posIT <- function(str,src) {
  m<-gregexpr(str,src,perl=TRUE,ignore.case=TRUE) 
  #print(regmatches(src,m))
  return(m)
}

# Return matching strings in source using regmatches
matIT <- function(str,src) {
  m<-gregexpr(str,src,perl=TRUE,ignore.case=F) 
  #print(regmatches(src,m))
  return(regmatches(src,m))
}

# Find tags returned from grepexpr in the source and return unique matches as a regex pattern to use in subsequent searches
tagIT<-function(tags,src){
    tmp<-unique(sub("\\s$*","",unlist(regmatches(src,tags))))
    tmp<-sub("^\\s*","",tmp)
    regmatches(tmp,gregexpr("\\s",tmp))<-rep("\\s",length(tags))
    return(tmp)
  }

# SUBSTITUTE (Corpus)

# Here we actually need to change the corpus content, solution: pass as a list.
# Other complications, regmatches turns the corpus into a character list, solution: PlainTextDocument
# A call might look like this, where TMcorpus:  TMcorpus<-subIT(src<-list(tags=myTAGS,str="mine",cor=TMcorpus))

subIT  <- function(src) {
  if(any(which(src$tags==""))){src$tags<-src$tags[which(src$tags!="",arr.ind=T)]}
  ifelse((length(src$tags)==0),{
    print("No tags in src$tags !!!")
    return(src$cor)
  },{
    Ddata<-DMetaData(src$cor)
    m <-gregexpr(wordXX(src$tags),src$cor,perl=TRUE)
    subs<-regmatches(src$cor,m)
    regmatches(src$cor,m)<-rep(src$str,length(m))
    print(paste("Changed ",length(which(unlist(m)>0))," strings into: ",src$str))
    src$cor<-Corpus(VectorSource(src$cor))
    oldname<-names(Ddata)
    Ddata<-data.frame(cbind(Ddata,paste(m),paste(subs)),stringsAsFactors=options(stringsAsFactors=FALSE))
    names(Ddata)<-c(oldname,paste(src$str,"~pos",sep=""),paste(src$str,"~tag",sep=""))
    DMetaData(src$cor)<-Ddata
    return(src$cor)
  })
}

subNIT  <- function(src) {
  if(any(which(src$tags==""))){src$tags<-src$tags[which(src$tags!="",arr.ind=T)]}
  ifelse((length(src$tags)==0),{
    print("No tags in src$tags !!!")
    return(src$cor)
  },{
    Ddata<-DMetaData(src$cor)
    m <-gregexpr(wordXX(src$tags),src$cor,perl=TRUE)
    subs<-regmatches(src$cor,m,invert=TRUE)
    regmatches(src$cor,m,invert=TRUE)<-rep(src$str,length(m))
    print(paste("Changed ",length(which(unlist(m)>0))," strings into: ",src$str))
    src$cor<-Corpus(VectorSource(src$cor))
    oldname<-names(Ddata)
    Ddata<-data.frame(cbind(Ddata,paste(m),paste(subs)),stringsAsFactors=options(stringsAsFactors=FALSE))
    names(Ddata)<-c(oldname,paste(src$str,"~pos",sep=""),paste(src$str,"~tag",sep=""))
    DMetaData(src$cor)<-Ddata
    return(src$cor)
  })
}

subNIT2  <- function(src) {
  if(any(which(src$tags==""))){src$tags<-src$tags[which(src$tags!="",arr.ind=T)]}
  ifelse((length(src$tags)==0),{
    print("No tags in src$tags !!!")
    return(src$cor)
  },{
    Ddata<-DMetaData(src$cor)
    m <-gregexpr(src$tags,src$cor,perl=TRUE)
    subs<-regmatches(src$cor,m,invert=TRUE)
    regmatches(src$cor,m,invert=TRUE)<-rep(src$str,length(m))
    print(paste("Changed ",length(which(unlist(m)>0))," strings into: ",src$str))
    src$cor<-Corpus(VectorSource(src$cor))
    oldname<-names(Ddata)
    Ddata<-data.frame(cbind(Ddata,paste(m),paste(subs)),stringsAsFactors=options(stringsAsFactors=FALSE))
    names(Ddata)<-c(oldname,paste(src$str,"~pos",sep=""),paste(src$str,"~tag",sep=""))
    DMetaData(src$cor)<-Ddata
    return(src$cor)
  })
}


# PATTERN GENERATORS 

# ( X1|X2 ) Example: wordXX(c("yes","no","maybe"))
wordXX   <- function(x,clps="|",pre="\\b(",post=")\\b") {
  paste(pre, paste(x,collapse=clps), post, sep="")
}


# ( XspaceY ) Example: wordXsY("yes","no")
wordXsY   <- function(x,y,clps="|",pre="\\b(",post=")\\b") {
  paste(pre, paste(x,collapse=clps), ")[[:blank:]]*(", paste(y,collapse=clps), post, sep="")
}

# ( XspaceY? ) Example: wordXsY("yes","no")
wordXsYo   <- function(x,y,clps="|",pre="\\b(",post=")?\\b") {
  paste(pre, paste(x,collapse=clps), ")[[:blank:]]*(", paste(y,collapse=clps), post, sep="")
}


# ( XspaceY | YspaceX ) Example: wordXsYr("yes","no")
wordXsYr  <- function(x,y,clps="|",pre="\\b(",post=")\\b") {
  paste("(",paste(pre,paste(x,collapse=clps),")[[:blank:]]*(",paste(y,collapse=clps),")\\b",sep=""),
    ")|(",paste("\\b(",paste(y,collapse=clps),")[[:blank:]]*(",paste(x,collapse=clps),post,sep=""),")", sep="")
}
# 
# # ( XspaceYspaceZ | YspaceZspaceX |  ZspaceXspaceY ) Example: wordXsYsZr("yes","no","maybe")
# wordXsYsZr <- function(x,y,z) {
#   paste("(",paste("\\b(",paste(x,collapse="|"),")\\s(",paste(y,collapse="|"),")\\s(",paste(z,collapse="|"), ")\\b", sep=""),
#     ")|(",paste("\\b(",paste(y,collapse="|"),")\\s(",paste(z,collapse="|"),")\\s(",paste(x,collapse="|"), ")\\b", sep=""),
#     ")|(",paste("\\b(",paste(z,collapse="|"),")\\s(",paste(x,collapse="|"),")\\s(",paste(y,collapse="|"), ")\\b", sep=""),")")
# }

# (XspaceY)? | (YspaceY)?
wordXsYro <- function(x,y,clps="|",pre="\\b(",post=")?\\b") {
  paste("(",paste(pre, paste(x,collapse=clps),")[[:blank:]]*(",paste(y,collapse=clps),post,sep=""),
    ")|(",paste(pre, paste(y,collapse=clps),")[[:blank:]]*(",paste(x,collapse=clps),post,sep=""),")",sep="")
}

# LOOK FOR NEIGHBOURS (WORDS)

# ( Npre X Npos ) Example: txt="Love thy neigbours"  
#                           matIT(NwN(1,"thy",0),txt)

NwN   <- function(Npre=NULL, word=NULL, Npos=NULL){
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

# ( Npre X Nmid Y Npos ) Example:  txt="Didn't you know? Love thy neigbours' neigbours, too!"  
#                                  matIT(NwNwN(0,c("Love","hate","you"),2,c("neigbours","friends","family"),0),gsub("[[:punct:]]+","",txt))  
#                                  matIT(NwNwN(1,"Love",3,"too",0),gsub("[[:punct:]]+","",txt))                                   
NwNwN   <- function(Npre=NULL, word1=NULL, Nmid=NULL, word2=NULL, Npos=NULL){
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

# (Npre X Nmid Y Npos) | (Npre Y Nmid X Npos) Example:    txt="Didn't you know? Your neigbours love you!"  
#                                                         matIT(NwNwNr(1,"know",2,"neigbours",1),gsub("[[:punct:]]+","",txt))  
#                                                         matIT(NwNwNr(0,"Your",2,"love",0),gsub("[[:punct:]]+","",txt))

NwNwNr  <- function(Npre=NULL, word1=NULL, Nmid=NULL, word2=NULL, Npos=NULL){
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

# Function to get ESCIs based on the test statistic
getESCI <- function(data, CL=.95){
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
