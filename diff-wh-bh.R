## Visualisation and quantification of differences between serial data and transmission-level data

# Plots violin plot for comparing serial against transmission level data --------

split.violin <- function(){
  library(coin)
  gpn <- read.csv("Org.data/gpn.csv")
  '%ni%' <- Negate('%in%')
  source('emp-data-prep.R')
  
  GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                             draw_group = function(self, data, ..., draw_quantiles = NULL) {
                               data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                               grp <- data[1, "group"]
                               newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                               newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                               newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                               
                               if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                                 stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                           1))
                                 quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                                 aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                                 aesthetics$alpha <- rep(1, nrow(quantiles))
                                 both <- cbind(quantiles, aesthetics)
                                 quantile_grob <- GeomPath$draw_panel(both, ...)
                                 ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                               }
                               else {
                                 ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                               }
                             })
  
  geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                                draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                                show.legend = NA, inherit.aes = TRUE) {
    layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
          position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
          params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
  }
  
  gpn <- read.csv("Org.data/gpn.csv")
  data <- read.csv("Bacterial.Datasets/empiricalsnps.csv",stringsAsFactors = FALSE) # load file
  data <- subset(data,N>=25)
  
  data$datpN <- data$x / (data$x + data$y)
  
  data$labels <-  str_replace(paste0("italic(", data$org, ")"), " ", "~")
  
  data <- lapply(1:dim(data)[1], function(x) {
    # apply function
    org <- data$org[x]
    gp <- gpn[which(grepl(org,gpn$Organism)),]
    d <- data.frame(data[x,],g = gp$g, pN = unlist(gp$pN) ,row.names = NULL,stringsAsFactors=FALSE)
    return(d)
  }) %>% bind_rows()
  
  all.data <- data
  
  pdf("splitviolin.pdf",  width=7.58, height=5.69)
  comp.data <- all.data[which(all.data$org %ni% c("B. pertussis", "E. coli")),]
  g1 <- ggplot(comp.data)+
    geom_hline(aes(yintercept=pN) ,size=0.5)+
    geom_split_violin(alpha=0.5,aes(x=as.factor(org),y=x/(x+y),fill=type)) +
    geom_boxplot(aes(x=as.factor(org),y=x/(x+y),color=type),width=0.2,fill="white",position = "dodge2",size=1,outlier.size = 2.5)+ 
    facet_wrap(~labels,scales="free_x",labeller = label_parsed)+
    theme_bw() +
    labs(x="Organism", y = "Proportion non-synonymous") +
    scale_color_manual(name="Type",labels=c("Population","Serial"),values=c("#F8766D","#00B0F6"))+
    scale_fill_manual(name="Type",labels=c("Population","Serial"),values=c("#F8766D","#00B0F6")) +
    theme(legend.position ="bottom",
          text = element_text(size=18),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.text = element_text(size=14),
          legend.title=element_text(size=14),
          legend.text.align = 0)+
    guides(fill=FALSE,linetype=FALSE)
  
  print(g1)
  
  dev.off()
}


# Quantifying differences between groups using independence test ----------------------------------------------

## Checking differences between groups (general vs. serial)
library(coin)

gpn <- read.csv("Org.data/gpn.csv")
data <- read.csv("Bacterial.Datasets/empiricalsnps.csv",stringsAsFactors = FALSE) # load file
data <- subset(data,N>=25)

data <- cbind(data,gpn[match(data$org,gpn$Organism),2:3],row.names=NULL)

comb.e.data <- data
comb.e.data$type <- as.factor(comb.e.data$type)
comb.e.data$datpN <- comb.e.data$x/(comb.e.data$x + comb.e.data$y )

independence_test(datpN ~ type, data=comb.e.data)
independence_test(datpN ~ type, data=comb.e.data[which(comb.e.data$org=="P. aeruginosa"),])
independence_test(datpN ~ type, data=comb.e.data[which(comb.e.data$org=="M. tuberculosis"),])
independence_test(datpN ~ type, data=comb.e.data[which(comb.e.data$org=="S. typhimurium"),])

independence_test(datpN ~ as.factor(org), data=comb.e.data[which(comb.e.data$type=="serial"),])
independence_test(datpN ~ as.factor(org), data=comb.e.data[which(comb.e.data$type=="popn"),])
independence_test(datpN ~ as.factor(org), data=comb.e.data)

