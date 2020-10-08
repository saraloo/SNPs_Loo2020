## This file calculates the probability of non-synonymous SNP given a bacterial codon frequency (saved as a .csv)
## Codon frequencies for organisms investigated in Loo et al 2020 are saved in folder Org.data


# Functions to obtain rate of synonymous SNP ------------------------------------------------

get.codons <- function(){
  
  # build vector of 64 codons
  nucleotides =  strsplit("TCAG","")[[1]]
  codons=c()
  for (n1 in nucleotides) {
    for (n2 in nucleotides) {
      for (n3 in nucleotides) {
        codons = c(codons,paste(n1,n2,n3, sep=""))
      }
    }
  }
  
  # animo acids
  AA = c(
    "F", "F", "L", "L", "S", "S", "S", "S", "Y", "Y", "STOP", "STOP", "C", "C",  "STOP", "W",
    "L", "L", "L", "L", "P", "P", "P", "P", "H", "H", "Q", "Q", "R", "R", "R", "R",
    "I", "I", "I", "M", "T", "T", "T", "T", "N", "N", "K", "K", "S", "S", "R", "R",
    "V", "V", "V", "V", "A", "A", "A", "A", "D", "D", "E", "E", "G", "G", "G", "G")
  
  # associate each codon to its animo acid
  names(AA) = codons
  
  subs.AA <- AA[which(AA != "STOP")]
  subs.codons <- codons[which(AA != "STOP")]
  
  # given a nucleotide return the transition
  transition = function(l) {
    list("T"="C", "C"="T", "A" = "G", "G"="A")[[l]]
  }
  
  # function to get mutants of codon, 3 first are transitions
  SNP = function (codon) {
    
    # split condon in 3 nucleotides
    cod = strsplit(codon, "")[[1]]
    
    # transitions
    mutants=c()
    for (i in 1:3) {
      mutant = cod
      mutant[i] = transition(cod[i])
      mutants = c(mutants, paste(mutant, collapse=""))
    }
    
    # easier to redo all of them at the same time
    # than exclude transitions
    for (i in 1:3) {
      for (ldest in nucleotides[nucleotides != cod[i]]) {
        mutant = cod
        mutant[i] = ldest
        mutants = c(mutants, paste(mutant, collapse=""))
      }
    }
    
    
    # identify synonimous mutations
    as.integer(subs.AA[unique(mutants)] == subs.AA[[codon]])
  }
  
  # build 61 x 9 matrix of codons x snps (without stop codons)
  # 1 indicates a synonymous mutation
  M = t(sapply(subs.codons, SNP))
  
  
  return(M)
  
}


read.bac.get.synrate <- function(fname="Org.data/Mtuberc_codonfreq.csv",
                                 prb.trans = 0.5){
  
  M <- get.codons()
  
  # get codon frequencies
  dat <- read.csv(fname)
  colnames(dat)[1]<- c("Codon")
  cds = as.character(dat$Codon)
  cds = sapply(strsplit(cds," "), paste, collapse="")
  cds = gsub("U", "T", cds)
  f <- dat$f[match(rownames(M),cds)] #freq/num in same order as rows of M
  f <- f/sum(f)
  
  # Codon to amino acid map
  
  synrate <- function(p){

    a <- c()
    b <- c()
    cod.p <- rep(0,61)
    
    for(i in 1:length(f)){
      avec <- M[i,1:3]
      bvec <- M[i,4:9]
        
      a[i] <- length(avec[which(is.na(avec) != TRUE)])
      b[i] <- length(bvec[which(is.na(bvec) != TRUE)])
        
      cod.p[i] <- f[i] * M[i,][complete.cases(M[i,])] %*% c(rep(p/a[i],a[i]), rep((1-p)/b[i],b[i]))
        
    }

    return(sum(cod.p))
    
  }
  
  
  
  # distribution of probabilities of syn mutation given bacterial codon frequency f 
  prob.synonymous <- sapply(0:100 / 100, synrate)
  # probability of syn mutation 
  synrate.at.prbtrans <- sapply(prb.trans, synrate)
  
  pN <- 1 - synrate.at.prbtrans
  
  list(prob.synonymous = prob.synonymous, pN = pN)
  
}




# Get g and pN from codon frequencies ------------------------------------------------

library(tidyverse)
library(readxl)
library(knitr)

tab <- read_excel("Org.data/org_traits.xlsx")
kable(tab,caption = "Organism traits from literature")

M <- read.bac.get.synrate(fname="Org.data/Mtuberc_codonfreq.csv",
                          prb.trans = tab$`Probability of transition`[which(tab$Organism=="M. tuberculosis")])

E <- read.bac.get.synrate(fname="Org.data/ecoli_codonfreq.csv",
                          prb.trans = tab$`Probability of transition`[which(tab$Organism=="E. coli")])

S <- read.bac.get.synrate(fname="Org.data/styphim_codonfreq.csv",
                          prb.trans = tab$`Probability of transition`[which(tab$Organism=="S. typhimurium")])

P <- read.bac.get.synrate(fname="Org.data/paerug_codonfreq.csv",
                          prb.trans = tab$`Probability of transition`[which(tab$Organism=="P. aeruginosa")])

B <- read.bac.get.synrate(fname="Org.data/bpertuss_codonfreq.csv",
                          prb.trans = tab$`Probability of transition`[which(tab$Organism=="B. pertussis")])

params <- data.frame(Organism=as.character(tab$Organism),g=tab$`Coding region (%)`)
rownames(params)<-sapply(1:length(params$Organism), function(i) strsplit(as.character(params$Organism),". ")[[i]][1])
pN <- data.frame(E=E$pN,B=B$pN,M=M$pN,S=S$pN,P=P$pN)
params$pN <- t(pN[match(rownames(params),colnames(pN))])

write.csv(params,file="Org.data/gpn.csv",row.names = FALSE)


# Plot change in non-synonymous rate as probability of transition changes for all organisms ------------------------------------------------

# Figure 1

pdf("nonsyn-vs-trans.pdf",height=7,width=7)
op <- par(mar=c(4,4,2,1))

# setup data frame of all synonymous mutation rates, with organism as variable
df <- data.frame(x=seq(0,1,length.out=101),mtb=M$prob.synonymous,ec=E$prob.synonymous,st=S$prob.synonymous,paer=P$prob.synonymous,bp=B$prob.synonymous)
mdf <- reshape2::melt(df,id.var="x")

# setup data frame to plot vertical lines to show prob of transition from literature
a <- c(tab$`Probability of transition`[which(tab$Organism=="M. tuberculosis")],tab$`Probability of transition`[which(tab$Organism=="E. coli")],tab$`Probability of transition`[which(tab$Organism=="S. typhimurium")],tab$`Probability of transition`[which(tab$Organism=="P. aeruginosa")],tab$`Probability of transition`[which(tab$Organism=="B. pertussis")])
vlines <- data.frame(xint = c(a),variable = unique(mdf$variable))

labl <- list(expression(italic('M. tuberculosis')), expression(italic('E. coli')),expression(italic('S. typhimurium')),expression(italic('P. aeruginosa')),expression(italic('B. pertussis'))) 

gg <- ggplot(mdf,aes(x = x, y = 1-value, colour = variable)) + 
  labs(x = "Probability of transition", y = "Probability of non-synonymous mutations") +
  geom_line(size=1) +
  scale_colour_discrete(name="Organism", labels = labl) +
  geom_vline(data = vlines,aes(xintercept = xint,colour = variable), linetype = "dashed") +
  theme_bw()+theme(aspect.ratio=1,text = element_text(size=16)) + coord_cartesian(xlim=c(0,1))

print(gg)

dev.off()



