## Function for obtaining confidence intervals given pN confidence interval type

fun.ci <- function(pN,type = "agresti-coull"){ 
  
  # x is synonymous, y is non-synonymous
  x0<- seq(from=0.2, to=10000, length.out = 30000)
  y0<- x0*pN/(1-pN)
  n <- x0+y0
  z1 <- qnorm(0.975)
  z2 <- qnorm(0.025)
  
  if (type != "binom"){
    if (type == "wald"){
      # Wald
      p.up <- pN + z1 * sqrt(pN*(1-pN)/n) 
      p.down <- pN + z2 * sqrt(pN*(1-pN)/n)
    }
    else if (type == "wald-cc"){
      # Wald with continuity correction
      p.up <- pN + z1 * sqrt(pN*(1-pN)/n) +1/(2*n)
      p.down <- pN + z2 * sqrt(pN*(1-pN)/n) -1/(2*n)
    }
    else if (type == "agresti-coull"){
      # Agresti-Coull
      p.up <- (pN + (z1^2)/(2*n) + z1*sqrt(pN*(1-pN)/n + (z1^2)/(4*n^2)))/(1+(z1^2)/n)
      p.down <- (pN + z2^2/(2*n) + z2*sqrt(pN*(1-pN)/n + z2^2/(4*n^2)))/(1+z2^2/n)
    }
    else if (type == "clopper-pearson"){
      # Clopper-Pearson
      F.1 <- qf(1-0.05/2,df1=2*(y0+1),df2=2*(n-y0))
      p.up <- (y0+1)*F.1/(n-y0+(y0+1)*F.1)
      F.2 <- qf(1-0.05/2,df1=2*(n-y0+1),df2=2*y0)
      p.down <- y0/(y0+(n-y0+1)*F.2)
    }
    
    y1 <- n*p.up
    x1 <- n-y1 
    y2 <- n*p.down
    x2 <- n-y2
    
  }else if (type == "binom"){
    n <- seq(1,10000)
    y2 <- sapply(n,function(i) qbinom(p=0.025,prob=pN,size=i))
    x2 <- n-y2
    y1 <- sapply(n,function(i) qbinom(p=0.975,prob=pN,size=i))
    x1 <- n-y1
  }

  list(x1 = x1, y1 = y1, x2 = x2, y2 = y2)
  
}
