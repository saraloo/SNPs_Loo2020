## Visualisation for Rule of thumb guideline Figure 2

# Neutral simulations plot (Figure 2A) -----------------------------------------------------------------
library(scales)
source("conf-int.R")

sim.data <- read.csv("Simulation.data/neutral.csv")
sim.data <- sim.data[,-c(1)]
pN <- 0.73
d <- as.data.frame(lapply(pN,fun.ci,type="binom"))

gg <- ggplot(sim.data,aes(s.all,ns.all))+
  geom_point(size=1.5,aes((s.all), (ns.all),color=gens),alpha=0.8) +
  geom_abline(aes(intercept=log10(pN/(1-pN)),slope=1),colour="black") +  
  geom_line(data=d,aes(x1,y1),colour="black",linetype="dashed") +
  geom_line(data=d,aes(x2,y2),colour="black",linetype="dashed") +
  labs(x = "Synonymous", y = "Non-synonymous") +
  theme_bw() +theme(legend.justification=c(0,0))+
  scale_color_viridis_c(name="Time\n (generations)")+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))
gg1 <- gg +
  coord_cartesian(xlim=c(10^0.25,10^2.2),ylim=c(10,10^2.5),expand=FALSE)+
  theme(aspect.ratio=1,text = element_text(size=16),
        legend.position=c(0.7,0.05),legend.text = element_text(size = 10),legend.title=element_text(size=12)) 


# Rule of thumb -----------------------------------------------------------
library(reshape2)

tot.snps <- 20:250

## Binomial confidence limits
pl <- function(phat=0.73,n) {
  phat - 1.96*sqrt(phat*(1-phat)/n)
}
pu <- function(phat=0.73,n) {
  phat + 1.96*sqrt(phat*(1-phat)/n)
}


df <- data.frame(pl=pl(0.73,n=tot.snps),pu=pu(0.73,n=tot.snps),pN=0.73,n=tot.snps)
df <- rbind(df,data.frame(pl=pl(0.72,n=tot.snps),pu=pu(0.72,n=tot.snps),pN=0.72,n=tot.snps))
df <- rbind(df,data.frame(pl=pl(0.75,n=tot.snps),pu=pu(0.75,n=tot.snps),pN=0.75,n=tot.snps))
df$pN <- as.character(df$pN)
poly.df <- melt(subset(df,pN=="0.73"),id=c("n","pN"))
poly.df <- rbind(poly.df,
                 c(n=tail(tot.snps,1),pN=0.73,variable="pu",value=1),
                 c(n=head(tot.snps,1),pN=0.73,variable="pu",value=1),
                 c(n=tail(tot.snps,1),pN=0.73,variable="pl",value=0),
                 c(n=head(tot.snps,1),pN=0.73,variable="pl",value=0))
poly.df$value <- as.numeric(poly.df$value)
poly.df$n <- as.numeric(poly.df$n)
ends <- df[which(df$n==tail(tot.snps,1)),]

p1 <- ggplot(data=melt(df,id=c("pN","n")))+   
  geom_polygon(data=poly.df,mapping=aes(x=n, y=value, group=variable),fill="grey64")+
  geom_line(aes(x=n,y=value,color=variable,linetype=pN),size=1)+
  scale_colour_manual(" ",values = c("#F8766D","#F8766D","#F8766D"))+
  scale_linetype_manual(" ",values = c("solid","dashed","dotdash") )+
  annotate("text", x = 150, y = 0.87, label = "Positive selection")+
  annotate("text", x = 150, y = 0.73, label = "Consistent with neutrality")+
  annotate("text", x = 150, y = 0.57, label = "Purifying selection")+
  geom_label(data=ends, aes( x=225, y=pl-c(0.001,0.003,0.001), label=c("0.73","0.72","0.75")), 
             size=3.5 , angle=45,label.size=0,label.padding=unit(0.1, "lines"))+
  geom_label(data=ends, aes( x=225, y=pu+c(0.005,0.001,0.003), label=c("0.73","0.72","0.75")), 
             size=3.5 , angle=45,label.size=0,label.padding=unit(0.1, "lines"))+
  scale_fill_manual(values=rep("white",3))+
  theme_bw()+
  guides(color=FALSE,linetype=FALSE,fill=FALSE)+
  coord_cartesian(xlim=c(min(tot.snps),max(tot.snps)),ylim=c(min(df$pl),max(df$pu)),expand=FALSE)+
  labs(x="Total coding SNPs", y = "Observed non-synonymous proportion") +
  theme(aspect.ratio=1,legend.justification=c(0,0),
        text = element_text(size=16))




# Combine -----------------------------------------------------------------
library(cowplot)

pdf("Fig2.pdf", height=6,width=12)
plot_grid(gg1,p1, labels = "AUTO",scale = 0.9)
dev.off()
      