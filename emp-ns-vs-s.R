## Given empirical SNP data, functions for plotting non-synonymous vs. synonymous SNPs


plot.ns.s.emp <- function(data,type="agresti-coull"){
  
  cis <- emp.ci.prep(data,citype=type)

  if (deparse(substitute(data)) == "serial"){
    
    labl <- as.expression(lapply(levels(serial$org), function(x) bquote(italic(.(x)))))

    gg<- ggplot(data) + geom_point(aes(y,x,colour=org),data=data) +
      geom_abline(aes(intercept=log10(pN/(1-pN)),slope=1,colour=org),alpha=0.7) +
      geom_line(data=cis,aes(x1,y1,color=org),linetype="dashed",alpha=0.5)+
      geom_line(data=cis,aes(x2,y2,color=org),linetype="dashed",alpha=0.5)+
      labs(x = "Synonymous mutations", y="Non-synonymous mutations") + 
      theme_bw() + theme(aspect.ratio=1,text = element_text(size=14),legend.text.align = 0) +
      scale_colour_discrete(name="Organism", labels = labl) +
      scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                    labels = trans_format("log10", math_format(10^.x))) +
      scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                    labels = trans_format("log10", math_format(10^.x))) +
      coord_cartesian(xlim=c(10^0,10^4.2),ylim=c(10^0.25,10^4),expand=FALSE)
    quartz()
    print(gg)
    
  }
  else{
    
    data$labels <- str_replace(paste0("italic(", data$org, ")"), " ", "~")
    cis$labels <- str_replace(paste0("italic(", cis$org, ")"), " ", "~")
    
    # log scale
    gg<- ggplot(data) + geom_point(aes(y,x),data=data) +
      geom_abline(aes(intercept=log10(pN/(1-pN)),slope=1),colour="#F8766D") +
      geom_line(data=cis,aes(x1,y1),colour="#F8766D",linetype="dashed")+
      geom_line(data=cis,aes(x2,y2),colour="#F8766D",linetype="dashed")+
      labs(x = "Synonymous mutations", y="Non-synonymous mutations") + 
      facet_wrap(~labels,labeller = label_parsed) + theme_bw() + 
      theme(aspect.ratio=1,text = element_text(size=14),legend.text.align = 0) +
      scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                    labels = trans_format("log10", math_format(10^.x))) +
      scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                    labels = trans_format("log10", math_format(10^.x))) +
      coord_cartesian(xlim=c(10^0,10^4.2),ylim=c(10^0.25,10^4),expand=FALSE)
    quartz()
    print(gg)
    
  }
  
}

## Plots Figure 4
plot.ns.s.combine <- function(type="binom"){
  require(magrittr)
  require(dplyr)

  gpn <- read.csv("Org.data/gpn.csv")
  data <- read.csv("Bacterial.Datasets/empiricalsnps.csv",stringsAsFactors = FALSE) # load file 
  data <- subset(data,N>=25)
  
  data <- cbind(data,gpn[match(data$org,gpn$Organism),2:3],row.names=NULL)
  d.data <- emp.ci.prep(data[which(data$type!="serial"),],citype=type)
  d.serial <- emp.ci.prep(data[which(data$type=="serial"),],citype=type)
  # Add colour and box values
  data <- data.frame(data,col="1",box=data$org,stringsAsFactors = FALSE)
  data[which(data$type=="serial"),c("col","box")]<- data.frame(col=data[which(data$type=="serial"),"org"],box="serial",stringsAsFactors = FALSE)

  cols <- c("type", "assembly", "col","box")
  data %<>% mutate_at(cols, funs(factor(.)))
  
  data$labels <- ifelse(data$box != "serial",
                        str_replace(paste0("bold(",LETTERS[as.numeric(data$box)],"):","italic(", data$org, ")"), " ", "~"),
                        str_replace(data$box, "serial", paste0("bold(",LETTERS[nlevels(data$box)],"):Serial")))
  data$labels[which(data$labels=="bold(E):italic(S.~typhimurium)")] <- "bold(E):italic(S.)~Typhimurium"
  
  d.data <- data.frame(d.data,col="1",box=d.data$org)
  d.serial <- data.frame(d.serial,col=d.serial$org,box="serial")
  d.data <- rbind(d.data,d.serial)
  d.data$labels <- ifelse(d.data$box != "serial",
                          str_replace(paste0("bold(",LETTERS[as.numeric(d.data$box)],"):","italic(", d.data$org, ")"), " ", "~"),
                          str_replace(d.data$box, "serial", paste0("bold(",LETTERS[nlevels(d.data$box)],"):Serial")))
  d.data$labels[which(d.data$labels=="bold(E):italic(S.~typhimurium)")] <- "bold(E):italic(S.)~Typhimurium"

  colors <- c("black","#F8766D","#00BA38","#619CFF")
  colors2 <- c("gray","#F8766D","#00BA38","#619CFF")
  
  m<- rbind.fill(data,d.data)
  m$stroke <- ifelse(m$box != "serial",
                     0.01,
                     0.9)
  
  labl <- as.expression(lapply(levels(fct_inorder(data[which(data$type=="serial"),"org"])), function(x) bquote(italic(.(x)))))
  labl[3] <- expression(italic(S.)~Typhimurium)
  
  # Plot
  gg1 <- ggplot(data=m,aes(y,x)) 
  
  # Plot with two legends
  gg1 <- ggplot(data=m,aes(y,x)) 
  p1 <- gg1 +
    geom_point(aes(y,x,color=col,fill=col),size=2) +
    scale_shape_manual(values=c(16,9),breaks=unique(m$assembly)[1:2],labels=c("De novo","Reference strain"),name="Assembly")+
    scale_fill_manual(values=colors) +
    geom_abline(aes(intercept=log10(pN/(1-pN)),slope=1,colour=col)) +
    geom_line(data=subset(m,x1>0),aes(x1,y1,colour=col),linetype="dashed")+
    geom_line(aes(x2,y2,colour=col),linetype="dashed")+
    scale_color_manual(values=colors,name="Organism", breaks=unique(m$col)[-1],labels = labl) +
    labs(x = "Synonymous mutations", y="Non-synonymous mutations") +
    facet_wrap(~labels,labeller=label_parsed)+
    theme_bw() +
    theme(aspect.ratio = 1,legend.position = c(.91,.13) ,
          text = element_text(size=12),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          panel.spacing = unit(0.75, "lines"),
          legend.text = element_text(size=5),
          legend.title=element_text(size=5),
          legend.text.align = 0)+
          # ,plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    coord_cartesian(xlim=c(10^0,10^4.2),ylim=c(10^0.25,10^4),expand=FALSE)  +
    guides(fill="none",shape="none")
  
  pdf("nsvss-panels.pdf",  width=7, height=5)

  print(p1)
  
  dev.off()
  
}

