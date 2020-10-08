source("prop-out.R")
data <- read.csv("Simulation.data/fullsimdat.csv")

data <- data[,-c(1)]
ptab <- power.table(data)[c("gens", "s1", "s2","total","c")]

xpos <- c()
ypos <- c()
for(i in 1:dim(ptab)[1]){
  if(ptab[i,]$s2 == 5e-2){
    xpos[i] <- max(data[which(data$gens==ptab[i,]$gens & data$s1 == ptab[i,]$s1 & data$s2 == ptab[i,]$s2),]$s.all)
    ypos[i] <- min(data[which(data$gens==ptab[i,]$gens & data$s1 == ptab[i,]$s1 & data$s2 == ptab[i,]$s2),]$ns.all)
  }else{
    xpos[i] <- min(data[which(data$gens==ptab[i,]$gens & data$s1 == ptab[i,]$s1 & data$s2 == ptab[i,]$s2),]$s.all)
    ypos[i] <- max(data[which(data$gens==ptab[i,]$gens & data$s1 == ptab[i,]$s1 & data$s2 == ptab[i,]$s2),]$ns.all)
  }
}

pt <- cbind(ptab,x=xpos,y=ypos)

# Manually adjusting location for proportions (visualisation)
pt[4,6] <- 40
pt[9,6:7] <- c(10^0.5,10^1.6)
pt[12,6:7] <- c(10^1.25,10^2.15)
pt[15,6:7] <- c(10^1.1,10^2)
pt[16,6:7] <- c(10^1.75,10^1.5)
pt[19,6:7] <- c(10^1.25,10^2.1)
pt[20,6:7] <- c(10^1.76,10^1.625)
pt[23,6:7] <- c(10^1,10^2.1)
pt[24,6:7] <- c(10^1.7,10^1.75)
pt[35,6:7] <- c(10^1.65,10^1.45)
pt[36,6:7] <- c(10^1.77,10^1.65)

pt$s1 <- factor(pt$s1,labels=c("s[a]: 10^{-5}","s[a]: 0.01","s[a]: 0.05"))
pt$s2 <- factor(pt$s2,labels=c("s[d]: 10^{-5}","s[d]: 0.01","s[d]: 0.05"))

data$s1 <- factor(data$s1,labels=c("s[a]: 10^{-5}","s[a]: 0.01","s[a]: 0.05"))
data$s2 <- factor(data$s2,labels=c("s[d]: 10^{-5}","s[d]: 0.01","s[d]: 0.05"))

gg <- ggplot(data,aes(s.all,ns.all))+
  geom_point(size=1.3,aes((s.all), (ns.all), colour = as.factor(gens),fill=as.factor(gens)),alpha=0.4) +
  geom_abline(aes(intercept=log10(pN/(1-pN)),slope=1),colour="black") +  
  geom_line(data=d,aes(x1,y1),colour="black",linetype="dashed") +
  geom_line(data=d,aes(x2,y2),colour="black",linetype="dashed") +
  labs(x = "Synonymous", y = "Non-synonymous") +
  theme_bw() +theme(legend.justification=c(0,0))+guides(fill=FALSE)+
  scale_color_discrete(name="Time \n (generations)") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))

gg1 <- gg + 
  geom_text(data = pt, aes(x = x,  y = y, label = as.character(c),color=as.factor(gens)),show.legend=F)+
  facet_grid(s2~s1,labeller = label_parsed) +
  theme(aspect.ratio=1,text = element_text(size=24),panel.spacing = unit(1.2, "lines"),
        legend.text = element_text(size = 16)) +
  coord_cartesian(xlim=c(10^0.25,10^2),ylim=c(10^0.75,10^2.5),expand=FALSE)

  
pdf("panels500sims.pdf",onefile = TRUE,height=10,width=12)

print(gg1)
dev.off()
