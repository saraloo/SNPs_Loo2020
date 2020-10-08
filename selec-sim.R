## This model assumes that only nonsynonymous mutations are under selection
##  - functions for simulation model
## To run model: 
#   - sim.samples(p=setp(s1=0,s2=0), fname="Simulation.data/neutral.csv", sims=20, sample.size=100)



## Set the parameters of the model 
setp <- function(
  N = 1e4, # popsize 
  mu = 0.01,
  s1 = 0.0,  # positive sel effect 
  s2 = 0.0,  # negative sel effect
  pNS= 0.75,  # proportion of mutations that are nonsynonymous
  pa = 0.05,  # prop advantageous given nonsyn
  pd = 0.8,  # prop deleterious given nonsyn 
  gens = 200,  # num generations to run sim
  sample.size=200){
  list(
    N=N,
    mu=mu,
    s1=s1,
    s2 =s2, 
    pNS=pNS, 
    pa =pa,
    pd = pd,
    gens = gens,
    sample.size=sample.size
  )
}

## Compute vector of fitnesses times n (counts) 
fitness <- function(p, pop) {
  #   genos <- length(pop$n)
  wn <- NULL
  for (i in 1:length(pop$n)) {
    wn[i] <- pop$n[i]*(1+p$s1)^pop$g1[i] * (1-p$s2)^pop$g2[i]
  }
  wn
}


## Now do mutations
mutate <- function(p,pop) {
  ## Components of genotype, newly made by mutation (appended to list)
  # numbers of each category
  new.g0 <- c()  # number of neutral nonsyn muts  ... ns.n? NSN
  new.g1 <- c()  # num pos selected nonsyn muts       ns.a  NSA
  new.g2 <- c()  # num neg sel nonsyn muts            ns.d  NSD
  new.gs  <- c() # num synonymous                     syn   S
  new.n <- c()   # new entries to be appended to n
  vec.n <- c()
  new.nlist <- c()
  
  if (length(pop$n) == 1){
    min.mut.num <- 0
  }else{
    min.mut.num <- max(as.numeric(gsub(".*?([0-9]+).*", "\\1", unlist(pop$nlist))))
  }
  
  nn <- pop$n  # make a copy of n in pop 
  for (i in 1:length(pop$n)) {
    muts <- rbinom(1, pop$n[i], p$mu)  # how many mutants? 
    new.muts.nums <- (min.mut.num+1):(muts+min.mut.num)  # set up vector of nums of new mutants
    new.n <- c(new.n, rep(1,muts))  # add muts 1s to vector of new n elements (densities of new mutants)
    vec.n <- c(vec.n, new.muts.nums)  # add mut nums to vector of new n elements (mutation number)
    if (muts>0) {
      for (j in 1:muts) {
        ##            new.n <- c(new.n,1)
        if (runif(1)<p$pNS) { # if nonsynonymous 
          ran01 <- runif(1)
          if (ran01<p$pa) { # advantageous nonsyn mutation
            new.g0 <- c(new.g0, pop$g0[i])    # inherit parent's g0 state
            new.g1 <- c(new.g1, pop$g1[i]+1)  # append new type with 1+
            new.g2 <- c(new.g2, pop$g2[i])    # inherit parent's g2 state
            new.gs <- c(new.gs, pop$gs[i])    # inherit parent's gs state
            new.nlist <- append(new.nlist,list(c(pop$nlist[[i]],paste0("A",vec.n[j]))))
          } else if (ran01<(p$pa+p$pd)) { # deleterious nonsynon 
            new.g0 <- c(new.g0, pop$g0[i])
            new.g1 <- c(new.g1, pop$g1[i])
            new.g2 <- c(new.g2, pop$g2[i]+1)
            new.gs <- c(new.gs, pop$gs[i])   
            new.nlist <- append(new.nlist,list(c(pop$nlist[[i]],paste0("D",vec.n[j]))))
          } else {   # neutral nonsynonymous
            new.g0 <- c(new.g0, pop$g0[i]+1)
            new.g1 <- c(new.g1, pop$g1[i])
            new.g2 <- c(new.g2, pop$g2[i])
            new.gs <- c(new.gs, pop$gs[i])   
            new.nlist <- append(new.nlist,list(c(pop$nlist[[i]],paste0("N",vec.n[j]))))
          }
        } else {   # synonymous mut 
          new.g0 <- c(new.g0, pop$g0[i])
          new.g1 <- c(new.g1, pop$g1[i])
          new.g2 <- c(new.g2, pop$g2[i])
          new.gs <- c(new.gs, pop$gs[i]+1) 
          new.nlist <- append(new.nlist,list(c(pop$nlist[[i]],paste0("S",vec.n[j]))))
          
        }
      }
      min.mut.num <- max(as.numeric(gsub(".*?([0-9]+).*", "\\1", unlist(new.nlist))))
      }
    nn[i] <- nn[i] - muts   # take away muts
  }
  #   cat("sum of n is now", sum(c(nn,new.n)), "\n")
  list(n = c(nn, new.n),
       g0 = c(pop$g0, new.g0),
       g1 = c(pop$g1, new.g1),
       g2 = c(pop$g2, new.g2),
       gs = c(pop$gs, new.gs),
       nlist = c(pop$nlist, new.nlist)
  )
}

## Compress the population object to remove extinct genotypes
compress <- function(pop) {
  with(pop, {
    list(
      n = n[n>0], g0 = g0[n>0], g1 = g1[n>0], g2 = g2[n>0], gs=gs[n>0],nlist=nlist[n>0]
    )
  })
}

## Simulate population over p$gens generations 
sim <- function (p) { # p = parameters
  ## initialise pop
  pop <- list(n=c(p$N),  # start with a monomorphic pop 
              g0=c(0),   # number of neutral nonsynon muts
              g1=c(0),   # number of advantageous nonsynon mutations
              g2=c(0),   # number of deleterious nonsynon muts
              gs=c(0),    # number of synonymous muts
              nlist=list(NULL)
  ) 
  
  totns <- 0 
  tots <- 0
  
  for (i in 1:p$gens) {
    probs <- fitness(p, pop)        # get w*n vector 
    nprime <- rmultinom (1, p$N, probs)   # reproduce, with fitness effects
    pop <- list(n=nprime, g0=pop$g0, g1=pop$g1, g2=pop$g2, gs=pop$gs,nlist=pop$nlist) # update n in pop
    pop <- compress(mutate(p,pop)) # then mutate
  }
  pop
  ## Return the pop 
  pop
}


## Random multivariate hypergeometric
##  INPUT
##    vector = vector of integers (balls in an urn) 
##    size = number of balls to take out the urn 
rmvh <- function(size,vector){
  #   cat("size,vector in rmvh: ",size, vector,"\n")
  if (size>sum(vector)|size<0) {
    stop("size out of bounds.")
  }
  L <- length(vector)
  drawn <- 0
  out <- c()
  if (L>1) {
    for (i in 1:(L-1)) {
      out[i] <- rhyper(1,vector[i],sum(vector[(i+1):L]),size-drawn)
      drawn <- drawn+out[i]
    }
    out[L] <- size-drawn
  }
  else if (L==1)   # only one element 
    out[1] <- size
  else # if (L<1) 
    out[1] <- 0   
  out
}


## Take sample and compute stats
sample.stats <- function(sample.size, pop) {
  samp.n <- rmvh(sample.size, pop$n) 
  ## Get corresponding genotypes; compress in final list
  sample <- list(g0=pop$g0[samp.n>0],
                 g1=pop$g1[samp.n>0],
                 g2=pop$g2[samp.n>0],
                 gs=pop$gs[samp.n>0],
                 n =samp.n[samp.n>0],
                 nlist=pop$nlist[samp.n>0])
  gen.no <- length(sample$nlist) # no. different genomes sampled
  all.snps <- unlist(sample$nlist)
  all.snps <- table(all.snps)[which(table(all.snps)<gen.no)]
  snp.names <-as.data.frame(all.snps)$all.snps
  
  ns.d <-    snp.names[grep("D",snp.names)]
  ns.a <-  snp.names[grep("A",snp.names)]
  ns.n <-  snp.names[grep("N",snp.names)]
  syn <-  snp.names[grep("S",snp.names)]
  ns.all <- length(ns.d)+length(ns.a)+length(ns.n)  # unique snps (assuming no polymorphics)
  ns.snps <- sample$g0 + sample$g1 + sample$g2 # num of ns snps per sequence
  cod.snps <- sample$g0 + sample$g1 + sample$g2 +sample$gs
  var.ns.snps <- var(ns.snps)  # variance in num ns snps per sequence
  var.s.snps <- var(sample$gs)  # variance in num ns snps per sequence
  ns.prop <- ns.all / (ns.all + length(syn))
  s.all <- length(syn)
  g0.all <- sum(sample$g0 * sample$n)
  g1.all <- sum(sample$g1 * sample$n)
  g2.all <- sum(sample$g2 * sample$n)
  
  c(ns.all=ns.all,
    ns.prop=ns.prop,
    max.snps=max(cod.snps),
    max.ns.snps=max(ns.snps),
    var.ns.snps=var.ns.snps,
    ns.range = max(ns.snps) - min(ns.snps),
    s.all=s.all,
    gneut.all=g0.all,
    gadv.all=g1.all,
    gdel.all=g2.all,
    var.s.snps = var.s.snps
  )
}


## Now simulate pop then sample
sim.samples <- function(p, fname="simdata", sims=20, sample.size=100) {
  sink(fname)
  cat("#  N mu s1 s2 pNS pa pd gens sample.size \n") 
  cat("# ", as.vector(unlist(p)), "\n")
  cat(" ns.all  ns.prop  max.snps  max.ns.snps  var.ns.snps  ns.range  s.all gneut.all gadv.all gdel.all var.s.snps \n ") 
  cat("# ns.all s.all gneut.all gadv.all gdel.all \n ") 
  for (i in 1:sims) {
    pop <- sim(p)               # run the population simulation  (p is parameters)
    ss <- sample.stats(sample.size, pop)
    cat(ss, "\n")
  }
  sink()
}
