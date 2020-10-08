## Code for plotting the visualisation given in Figure 5

# Tee plots ---------------------------------------------------------------
library(ggplot2)
library(ggforce)
library(wesanderson)
library(DescTools)
library(PropCIs)
library(scales)
library(plyr)
library(reshape2)
library(tidyverse)
library(ggplotify)
library(grid)

# 3D for ellipses for data 
conf <- function(x,y,z,type = "agresti-coull") {

  a <- data.frame(BinomCI(x=x,n=x+y,method=type,conf.level=0.95))
  b <- data.frame(BinomCI(x=x+y,n=x+y+z,method=type,conf.level=0.95))
  
  aup <- a$upr.ci
  alow <- a$lwr.ci
  bup <- b$upr.ci
  blow <- b$lwr.ci
  
  list(pnlength = abs(aup-alow)/2,
       glength = abs(bup-blow)/2)
  
}


ci.vect <- function(x,y,z,citype="wald"){
  hor.axis <- c()
  vert.axis <- c()
  old.hor <- c()
  old.vert <- c()
  for (i in 1:length(x)){
    old.hor[i] <- conf(x[i],y[i],z[i],type=citype)$pnlength
    old.vert[i] <- conf(x[i],y[i],z[i],type=citype)$glength
    hor.axis[i] <- 1.96 * sqrt(x[i]*y[i]/(x[i]+y[i])^3)
    vert.axis[i] <- 1.96 * sqrt((x[i]+y[i])*z[i]/(x[i]+y[i]+z[i])^3)
  }
  list(hor.axis = hor.axis,
       vert.axis = vert.axis,
       old.hor = old.hor, old.vert=old.vert)
}


# Plot empirical data with confidence intervals (ellipses) on t-diagram
plot.tee.ci.ellipse <- function(data, 
                                citype = "agresti-coull"){
  
  if ("org" %in% colnames(data)){
    
    if (deparse(substitute(data)) == "serial"){
      labl <- as.expression(lapply(levels(data$org), function(x) bquote(italic(.(x)))))
      
      data <- data.frame(data,col=data$org)
      CI <- ci.vect(data$x,data$y,data$z,citype=citype)
      
      gg <- ggplot(data,aes(x/(x+y), (x+y)/(x+y+z))) + 
        coord_cartesian(xlim=c(0,1),ylim = c(0,1),expand=FALSE)
      ellipse_plot <- gg + 
        geom_point(aes(color=org,shape=col),alpha=0.5,stroke=0.9)+
        scale_shape_manual(values=c(16,15,17),breaks=unique(data$org),name="Organism",labels=labl)+
        geom_segment(aes(x=0,y=g,xend=1,yend=g,color=org))+
        geom_segment(aes(x=pN,y=0,xend=pN,yend=g,color=org)) +
        geom_ellipsis(aes(x0=x/(x+y),y0=(x+y)/(x+y+z),
                          a = CI$hor.axis, b = CI$vert.axis,        
                          fill=org,angle=0),
                      linetype="blank",alpha=0.2,n=200)+
        labs(x = "Nonsynonymous/synonymous", y = "Coding/intergenic") + 
        scale_colour_discrete(name="Organism", labels = labl) +
        theme_bw() + theme(aspect.ratio = 1,text = element_text(size=14),legend.text.align = 0) + 
        guides(fill=FALSE) 
      ellipse_plot
    }
    else{
      data$labels <- str_replace(paste0("italic(", data$org, ")"), " ", "~")
      
      data <- data.frame(data,col=rep(1,length(data$x)))
      CI <- ci.vect(data$x,data$y,data$z,citype=citype)
      
      gg <- ggplot(data,aes(x/(x+y), (x+y)/(x+y+z))) + 
        coord_cartesian(xlim=c(0,1),ylim = c(0,1),expand=FALSE)
      ellipse_plot <- gg + 
        geom_point(alpha=0.6)+
        geom_segment(aes(x=0,y=g,xend=1,yend=g))+
        geom_segment(aes(x=pN,y=0,xend=pN,yend=g)) +
        geom_ellipsis(aes(x0=x/(x+y),y0=(x+y)/(x+y+z),
                          a = CI$hor.axis, b = CI$vert.axis,        
                          fill=col,angle=0),
                      linetype="blank",alpha=0.3,n=200)+
        facet_wrap(~labels,labeller=label_parsed)+
        labs(x = "Nonsynonymous/synonymous", y = "Coding/intergenic") + 
        scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1),
                           labels=c("0","0.25","0.5","0.75","1"))+
        scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),
                           labels=c("0","0.25","0.5","0.75","1"))+
        theme_bw() + theme(aspect.ratio = 1,text = element_text(size=14),
                           panel.spacing = unit(0.75, "lines")) + guides(fill=FALSE) 
      
      ellipse_plot
    }
    
  }
  
}



plot.tee.combine <- function(citype="agresti-coull"){

  gpn <- read.csv("Bacterial.Datasets/gpn.csv")
  data <- read.csv("Bacterial.Datasets/empiricalsnps.csv",stringsAsFactors = FALSE) # load file
  data <- subset(data,N>=25)
  
  data <- cbind(data,gpn[match(data$org,gpn$Organism),2:3],row.names=NULL)
   # Add colour and box values
  data <- data.frame(data,col="1",box=data$org,stringsAsFactors = FALSE)
  data[which(data$type=="serial"),c("col","box")]<- data.frame(col=data[which(data$type=="serial"),"org"],box="serial",stringsAsFactors = FALSE)
  
  cols <- c("type", "assembly", "col","box")
  data %<>% mutate_at(cols, funs(factor(.)))
  
  data$labels <- ifelse(data$box != "serial",
                        str_replace(paste0("bold(",LETTERS[as.numeric(data$box)],"):","italic(", data$org, ")"), " ", "~"),
                        str_replace(data$box, "serial", paste0("bold(",LETTERS[nlevels(data$box)],"):Serial")))
  data$labels[which(data$labels=="bold(E):italic(S.~typhimurium)")] <- "bold(E):italic(S.)~Typhimurium"
  colors <- c("black","#F8766D","#00BA38","#619CFF")
  colors2 <- c("gray","#F8766D","#00BA38","#619CFF")
  labl <- as.expression(lapply(levels(fct_inorder(data[which(data$type=="serial"),"org"])), function(x) bquote(italic(.(x)))))
  m<- rbind.fill(data)
  m$stroke <- ifelse(m$box != "serial",
                     0.01,
                     0.9)
  
  CI <- ci.vect(data$x,data$y,data$z,citype=citype)
  data$hor.ci <- CI$hor.axis
  data$vert.ci <- CI$vert.axis
  
  labl <- as.expression(lapply(levels(fct_inorder(data[which(data$type=="serial"),"org"])), function(x) bquote(italic(.(x)))))
  labl[3] <- expression(italic(S.)~Typhimurium)
    
  gg1 <- ggplot(data,aes(x/(x+y), (x+y)/(x+y+z))) + 
    coord_cartesian(xlim=c(0,1),ylim = c(0,1),expand=FALSE)
  p1 <- gg1 +
    geom_segment(aes(x=0,y=g,xend=1,yend=g,color=col))+
    geom_segment(aes(x=pN,y=0,xend=pN,yend=1,color=col)) +
    scale_colour_manual(values=colors,breaks=unique(m$col)[-1],name="Organism",labels = labl) +
    geom_point(aes(color=col,fill=col),alpha=0.7,stroke=0.9)+
    scale_shape_manual(values=c(16,9),breaks=unique(m$assembly)[1:2],labels=c("De novo","Reference strain"),name="Assembly")+
    geom_ellipse(aes(x0=x/(x+y),y0=(x+y)/(x+y+z),
                      a = hor.ci, b = vert.ci,        
                      fill=col,angle=0),
                  linetype="blank",alpha=0.1,n=200)+
    scale_fill_manual(values=colors) +
    facet_wrap(~labels,labeller=label_parsed)+
    labs(x = "Proportion of non-synonymous mutations", y = "Proportion of SNPs in coding regions") + 
    theme_bw() + 
    scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1),
                       labels=c("0","0.25","0.5","0.75","1"))+
    scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),
                       labels=c("0","0.25","0.5","0.75","1"))+
    theme(aspect.ratio = 1,legend.position = c(.775,.125) ,
          text = element_text(size=12),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          panel.spacing = unit(0.75, "lines"),
          legend.text = element_text(size=5),
          legend.title=element_text(size=5),
          legend.text.align = 0)+
    guides(fill=FALSE,shape=FALSE)

  pdf("test-ellipse-new.pdf",  width=7, height=5)
  
  
  print(p1)
  
  dev.off()

  
}
