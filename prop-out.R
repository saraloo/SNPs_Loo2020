
# Calculate proportion of simulated points outside of bounds


pn.bounds <- function(ns,s,pN = 0.73,citype="binom"){

  n <- ns + s
  
  zal1 <- qnorm(1-0.05/2)
  zal2 <- qnorm(0.05/2)
  
  if (citype == "wald"){
    p.up <- pN + zal1 * sqrt(pN/n)
    p.down <- pN + zal2 * sqrt(pN/n)
  }
  else if (citype == "wald-cc"){
    # Wald with continuity correction
    p.up <- pN + zal1 * sqrt(pN*(1-pN)/n) + 1/(2*n)
    p.down <- pN + zal2 * sqrt(pN*(1-pN)/n) - 1/(2*n)
  }
  else if (citype == "agresti-coull"){
    # Agresti-Coull
    p.up <- (pN + zal1^2/(2*n) + zal1*sqrt(pN*(1-pN)/n + zal1^2/(4*n^2)))/(1+zal1^2/n)
    p.down <- (pN + zal2^2/(2*n) + zal2*sqrt(pN*(1-pN)/n + zal2^2/(4*n^2)))/(1+zal2^2/n)
  }
  
  phat <- ns/(ns+s)
  
  df <- data.frame(phat=phat,p1=p.down,p2=p.up)
  
}


sim.prop.outbounds <- function(data, pN = 0.73,type="binom",get="prop"){
  
  if(type=="binom"){
    data$n <- data$ns.all+data$s.all
    data$phat <- data$ns.all/(data$ns.all+data$s.all)
    data$p1 <- sapply(data$n,function(x) qbinom(p=0.025,prob=pN,size=x)/x)
    data$p2 <- sapply(data$n,function(x) qbinom(p=0.975,prob=pN,size=x)/x)
    data$c <- rep(1,length(data$ns.all))
    
    sdata <- data[is.finite(data$phat),]
    
    for (i in 1:length(sdata$phat)){
      if(sdata$phat[i] >= sdata$p1[i] && sdata$phat[i] <= sdata$p2[i]){
        sdata$c[i] <- 0
      }
    }
  }else{
    
    thispn <- pN
    a<- do.call(function(ns.all,s.all,pN,...) pn.bounds(ns.all,s.all,pN=thispn,citype=type), data)
    a1 <- a[is.finite(rowSums(a)),]
    sdata <- data.frame(data[is.finite(rowSums(a)),],phat=a1$phat,p1=a1$p1,p2=a1$p2,c=1)
    
    for (i in 1:length(sdata$phat)){
      if(sdata$phat[i] > sdata$p1[i] & sdata$phat[i] < sdata$p2[i]){
        sdata$c[i] <- 0
      }
    }
    
  }
  
  
  if(get == "prop"){
    prop.out <- sum(sdata$c)/length(sdata$c)
    list(c(prop.out = prop.out))
  }else if(get == "df"){
    return(sdata)
  }
  
  
}

power.table <- function(data,type="binom"){
  
  dtest <- sim.prop.outbounds(data,get="df",type=type)
  
  power <- aggregate(c~gens+s1+s2,data=dtest,FUN = function(x)  sum(x)/length(x))
  total.snps <- aggregate(ns.all+s.all~gens+s1+s2,data=dtest,FUN= sum)
  
  tab.out <- data.frame(total = total.snps[,4],
                        power)
  
  return(tab.out)
  
}

ptab <- power.table(sim.data)[c("gens", "s1", "s2","total","c")]
ptab$total <- ptab$total / 100
ptab <- ptab[order( ptab[,1] ),]
print(xtable(ptab),include.rownames=FALSE)
