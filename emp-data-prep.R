# Read in empirical data. 
# Requires SNP data with three columns (x,y,z) = (NS,S,intergenic) numbers (saved in folder Bacterial.Datasets)
#     also requires a data frame of three columns (organism, g, pn) values for each bacterial organism under investigation (gpn.csv) 
#       - found using "calc-pN.R"

read.emp.data <- function(gpn,foldname="Bacterial.Datasets"){
  require(dplyr)
  data <- c()
  
  files <- list.files(path=foldname, pattern="*.csv", full.names=TRUE, recursive=FALSE)
  serial <- read.csv(files[grep("serial",files)])
  files <- files[- grep("serial",files)]
  files <- files[- grep("empiricalsnps",files)]
  
  data <- lapply(files, function(x) {
    t <- read.csv(x) # load file
    t <- as.data.frame(apply(t(t), 1, function(i) as.numeric(gsub(",","",i)))) # get rid of any commas and make numeric
    # apply function
    org <- strsplit(x,split="[/,.]")[[1]][3]
    gp <- gpn[which(grepl(org,gpn$Organism)),]
    d <- data.frame(org = as.character(org),g = gp$g, pN = unlist(gp$pN), t,row.names = NULL,stringsAsFactors=FALSE)
    return(d)
    }) %>% bind_rows()

  serial <- cbind(serial,gpn[match(serial$org,gpn$Organism),2:3],row.names=NULL)
  

  orglist <- as.character(unlist(lapply(data$org,function(i) gpn$Organism[which(grepl(i,gpn$Organism))])))
  data$org <- orglist
  
  list( data = data, 
        serial = serial
  )
  
}


# Empirical data requires some restructuring for plotting, in order to get confidence intervals

# To plot empirical data, need limits in dataset of certain structure
emp.ci.prep <- function(data,citype="agresti-coull") {

  test <- lapply(unique(data$pN),fun.ci,type=citype)
  names(test) <- unique(data$org)
  dd<- melt(test)
  # in dd: L1 = org, L2 = bound value
  tmp <- ddply(dd, .(L1, L2), transform, newid = paste(L1, seq_along(L2)))
  out<- dcast(tmp,L1+newid ~ L2,value.var="value")
  out<-out[with(out,order(L1,x1)),]
  
  out<- out %>% select(-2)
  colnames(out)<- c("org","x1","x2","y1","y2")

  
  ids <- out$org
  pns <- lapply(ids,function(i) gpn[which(grepl(i,gpn$Organism)),])
  pns <- lapply(1:length(ids),function(x) pns[[x]][3])
  pns<- as.vector(unlist(pns))
  
  out <- cbind(out,pN = pns)
  
  out
}



