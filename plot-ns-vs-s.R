# plot Ns vs S with neutral line and CIs - for simulations
## data = "Simulation.data/fullsimdat.csv"

plot.ns.s <- function(data, 
                      pN = 0.73, 
                      citype = "binom",
                      fname="test.pdf"){
  
  pdf(fname, onefile = TRUE,height=7.2,width=7.2)
  
  d <- as.data.frame(lapply(pN,fun.ci,type=citype))

    if(is.null(data$s1) == FALSE){
      if(length(unique(data$s1)) != 1){
        
        data$s1 <- factor(data$s1,labels=c("s[a]: 10^{-5}","s[a]: 0.01","s[a]: 0.05"))
        data$s2 <- factor(data$s2,labels=c("s[d]: 10^{-5}","s[d]: 0.01","s[d]: 0.05"))
        
        gg <- ggplot(data %>%
                       arrange(gens),aes(s.all,ns.all))+
          geom_point(size=1.3,aes((s.all), (ns.all), colour = as.factor(gens),fill=as.factor(gens)),alpha=0.35) +
          geom_abline(aes(intercept=log10(pN/(1-pN)),slope=1),colour="black") +  
          geom_line(data=d,aes(x1,y1),colour="black",linetype="dashed") +
          geom_line(data=d,aes(x2,y2),colour="black",linetype="dashed") +
          labs(x = "Synonymous", y = "Non-synonymous") +
          theme_bw() +theme(legend.justification=c(0,0))+guides(fill=FALSE)+
          scale_color_discrete(name="Time") +
          scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                        labels = trans_format("log10", math_format(10^.x))) +
          scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                        labels = trans_format("log10", math_format(10^.x)))
        gg1 <- gg + 
          facet_grid(s2~s1,labeller = label_parsed) +
          theme(aspect.ratio=1,text = element_text(size=18),panel.spacing = unit(1.3, "lines"),legend.text = element_text(size = 10)) +
          coord_cartesian(xlim=c(10^0.25,10^2.15),ylim=c(10^0.75,10^2.6),expand=FALSE)
        print(gg1)
      }
      if(length(unique(data$s1)) == 1){

        gg <- ggplot(data,aes(s.all,ns.all))+
          geom_point(size=1.5,aes((s.all), (ns.all),color=as.factor(gens)),alpha=0.8) +
          geom_abline(aes(intercept=log10(pN/(1-pN)),slope=1),colour="black") +  
          geom_line(data=d,aes(x1,y1),colour="black",linetype="dashed") +
          geom_line(data=d,aes(x2,y2),colour="black",linetype="dashed") +
          labs(x = "Synonymous", y = "Non-synonymous") +
          theme_bw() +theme(legend.justification=c(0,0))+
          scale_color_discrete(name="Time") +
          scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                        labels = trans_format("log10", math_format(10^.x))) +
          scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                        labels = trans_format("log10", math_format(10^.x)))
        gg1 <- gg +
          coord_cartesian(xlim=c(10^0.25,10^2.2),ylim=c(10,10^2.5),expand=FALSE)+
          theme(aspect.ratio=1,text = element_text(size=18),legend.position=c(0.8,0.05),legend.text = element_text(size = 14)) 
        print(gg1)
      }
    }else {
      
      gg <- ggplot(data,aes(s.all,ns.all))+
        geom_point(size=1.3,aes((s.all), (ns.all),color=as.factor(gens))) +
        geom_abline(aes(intercept=log10(pN/(1-pN)),slope=1),colour="black") +  
        geom_line(data=d,aes(x1,y1),colour="black",linetype="dashed") +
        geom_line(data=d,aes(x2,y2),colour="black",linetype="dashed") +
        labs(x = "Synonymous", y = "Non-synonymous") +
        theme_bw() +theme(legend.justification=c(0,0))+
        scale_color_discrete(name="Time") +
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x))) +
        scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x)))
      
      gg1 <- gg +
        coord_cartesian(xlim=c(1,10^2.5),ylim=c(10,10^3),expand=FALSE)
      print(gg1)
    }
  
  dev.off()
  
  return(gg1)
  
}




