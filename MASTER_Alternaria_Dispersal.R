
#######################################################
library(bbmle) 
library(glmmTMB)
library(dplyr)
library(fitdistrplus)
library(sjPlot)
library(AICcmodavg)
library(Rfast)
library(splines)
library(mgcv)
library(glmmTMB)
library(dplyr)
library(ggplot2)
library(multcomp)
library(multcomp)
library(FSA)
library(dunn.test)

## SECTION A: summarizing count data
################
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
r <- rectGrob(gp=gpar(fill="white", lwd=0)) # 1
A <- rectGrob(gp=gpar(fill="#6666ff", lwd=0)) # 2
B <-rectGrob(gp=gpar(fill="yellow", lwd=0)) # 3
C <-rectGrob(gp=gpar(fill="grey", lwd=0)) # 4



lay2 <- rbind(c(1,2,2,2,2,2,2,2,2,2,2,5,5,5,5,5,5,5,5,5,5),
              c(1,2,2,2,2,2,2,2,2,2,2,5,5,5,5,5,5,5,5,5,5),
              c(1,2,2,2,2,2,2,2,2,2,2,5,5,5,5,5,5,5,5,5,5),
              c(1,2,2,2,2,2,2,2,2,2,2,5,5,5,5,5,5,5,5,5,5),
              c(1,2,2,2,2,2,2,2,2,2,2,5,5,5,5,5,5,5,5,5,5),
              c(1,2,2,2,2,2,2,2,2,2,2,5,5,5,5,5,5,5,5,5,5),
              c(1,2,2,2,2,2,2,2,2,2,2,5,5,5,5,5,5,5,5,5,5),
              c(1,2,2,2,2,2,2,2,2,2,2,5,5,5,5,5,5,5,5,5,5),
              c(1,2,2,2,2,2,2,2,2,2,2,5,5,5,5,5,5,5,5,5,5),
              c(1,2,2,2,2,2,2,2,2,2,2,5,5,5,5,5,5,5,5,5,5),
              c(1,2,2,2,2,2,2,2,2,2,2,5,5,5,5,5,5,5,5,5,5),
              c(1,3,3,3,3,3,4,4,4,4,4,6,6,6,6,6,7,7,7,7,7),
              c(1,3,3,3,3,3,4,4,4,4,4,6,6,6,6,6,7,7,7,7,7),
              c(1,3,3,3,3,3,4,4,4,4,4,6,6,6,6,6,7,7,7,7,7),
              c(1,3,3,3,3,3,4,4,4,4,4,6,6,6,6,6,7,7,7,7,7),
              c(1,3,3,3,3,3,4,4,4,4,4,6,6,6,6,6,7,7,7,7,7),
              c(1,3,3,3,3,3,4,4,4,4,4,6,6,6,6,6,7,7,7,7,7),
              c(11,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8),
              c(11,12,12,12,12,12,12,12,12,12,12,13,13,13,13,13,13,13,13,13,13),
              c(11,9,9,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,10,10)
)

mrgn<-c(.25,.25,.25,.25)

adoniram<-theme(legend.position = "none", 
                axis.title.x=element_blank(),
                # axis.text.x=element_blank(),
                # axis.ticks.x=element_blank(),
                axis.title.y =element_blank(),
                # axis.text.y=element_blank(),
                # axis.ticks.y=element_blank(),
                plot.margin = unit(mrgn, "cm"))



#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}



box_data<-readRDS("Fig3_Box_data.rds")
  
ALTleg<-g_legend(
  ggplot(box_data[[1]][box_data[[1]]$lambda %in% c("UVA") & box_data[[1]]$hour.rescale %in% c(4),], aes(x=temp, y=mean, fill=as.factor(RH)))+
    geom_boxplot()+scale_fill_manual(name="RH%" ,values=c("#68287e","#6666ff", "#a95fa6", "#555a9a"))+
    theme(legend.position = "bottom",
          legend.key=element_blank(),
          panel.background = element_rect(fill = "white", colour = "white"),
          panel.grid=element_blank(),
          legend.title=element_text(size=8),
          legend.box.background = element_rect(color="white"),
          legend.box.margin = margin(10, 50, 10, 50)))

SOLleg<-g_legend(
  ggplot(box_data[[2]][box_data[[2]]$lambda %in% c("UVA") & box_data[[2]]$hour.rescale %in% c(4),], aes(x=temp, y=mean, fill=as.factor(RH)))+
    geom_boxplot()+scale_fill_manual(name="RH%",values=c("#dc4048","#f6821f", "#feb913", "#ed9fa3"))+
    theme(legend.position = "bottom",
          legend.key=element_blank(),
          panel.background = element_rect(fill = "white", colour = "white"),
          panel.grid=element_blank(),
          legend.title=element_text(size=8),
          legend.box.background = element_rect(color="white"),
          legend.box.margin = margin(10, 50, 10, 50)))


boxoutcl<-grid.arrange(
  grobTree(rectGrob(gp=gpar(col=NA,fill="white", alpha=1)), textGrob("Franction of Spores Germinated (at 96 hours)", rot=90)),
  ggplot(box_data[[1]][box_data[[1]]$lambda %in% c("UVA") & box_data[[1]]$hour.rescale %in% c(4),], aes(x=temp, y=mean, fill=as.factor(RH)))+
    geom_boxplot(width=.8, outlier.alpha = 0, lwd=.3,position = position_dodge2(preserve = "single"))+theme_bw()+
    ylim(0,1)+adoniram+scale_fill_manual(values=c("#68287e","#6666ff", "#a95fa6", "#555a9a"))+scale_x_discrete(labels = paste0(c(10,15,20,25), "\u00b0")),
  ggplot(box_data[[1]][box_data[[1]]$lambda %in% c("UVB") & box_data[[1]]$hour.rescale %in% c(4),], aes(x=temp, y=mean, fill=as.factor(RH)))+
    geom_boxplot(width=.8, outlier.alpha = 0, lwd=.3,position = position_dodge2(preserve = "single"))+theme_bw()+
    ylim(0,1)+adoniram+scale_fill_manual(values=c("#68287e","#6666ff", "#a95fa6", "#555a9a"))+scale_x_discrete(labels = paste0(c(10,15,20,25), "\u00b0")),
  ggplot(box_data[[1]][box_data[[1]]$lambda %in% c("b21") & box_data[[1]]$hour.rescale %in% c(4),], aes(x=temp, y=mean, fill=as.factor(RH)))+
    geom_boxplot(width=.8, outlier.alpha = 0, lwd=.3,position = position_dodge2(preserve = "single"))+theme_bw()+
    ylim(0,1)+adoniram+scale_fill_manual(values=c("#68287e","#6666ff", "#a95fa6", "#555a9a"))+scale_x_discrete(labels = paste0(c(10,15,20,25), "\u00b0")),
  ggplot(box_data[[2]][box_data[[2]]$lambda %in% c("UVA") & box_data[[2]]$hour.rescale %in% c(4),], aes(x=temp, y=mean, fill=as.factor(RH)))+
    geom_boxplot(width=.8, outlier.alpha = 0, lwd=.3,position = position_dodge2(preserve = "single"))+theme_bw()+
    ylim(0,1)+adoniram+scale_fill_manual(values=c("#dc4048","#f6821f", "#feb913", "#ed9fa3"))+scale_x_discrete(labels = paste0(c(10,15,20,25), "\u00b0")),
  ggplot(box_data[[2]][box_data[[2]]$lambda %in% c("UVB") & box_data[[2]]$hour.rescale %in% c(4),], aes(x=temp, y=mean, fill=as.factor(RH)))+
    geom_boxplot(width=.8, outlier.alpha = 0, lwd=.3,position = position_dodge2(preserve = "single"))+theme_bw()+
    ylim(0,1)+adoniram+scale_fill_manual(values=c("#dc4048","#f6821f", "#feb913", "#ed9fa3"))+scale_x_discrete(labels = paste0(c(10,15,20,25), "\u00b0")),
  ggplot(box_data[[2]][box_data[[2]]$lambda %in% c("b21") & box_data[[2]]$hour.rescale %in% c(4),], aes(x=temp, y=mean, fill=as.factor(RH)))+
    geom_boxplot(width=.8, outlier.alpha = 0, lwd=.3,position = position_dodge2(preserve = "single"))+theme_bw()+
    ylim(0,1)+adoniram+scale_fill_manual(values=c("#dc4048","#f6821f", "#feb913", "#ed9fa3"))+scale_x_discrete(labels = paste0(c(10,15,20,25), "\u00b0")),
  grobTree(rectGrob(gp=gpar(col=NA,fill="white", alpha=1)), textGrob("Temperature (ºC)", rot=0)),
  grobTree(rectGrob(gp=gpar(col=NA,fill="white", alpha=1)), ALTleg),
  grobTree(rectGrob(gp=gpar(col=NA,fill="white", alpha=1)), SOLleg),  
  grobTree(rectGrob(gp=gpar(col=NA,fill="white", alpha=1))),
  grobTree(rectGrob(gp=gpar(col=NA,fill="white", alpha=1)), textGrob(expression(italic("A. alternata")), rot=0)),
  grobTree(rectGrob(gp=gpar(col=NA,fill="white", alpha=1)), textGrob(expression(italic("A. solani")), rot=0)),
  layout_matrix=lay2
)


# Suplemental Figure 6

sup_count_dat<-readRDS("sup_count_dat.rds")
library(ggrepel)

ALTert.lm<-lm(mean~hour.rescale, data=sup_count_dat[[1]])
SOLert.lm<-lm(mean~hour.rescale, data=sup_count_dat[[2]])

repentir<-ggplot(data=NULL)+
  geom_jitter(
    data=sup_count_dat[[1]], 
    aes(x=as.character(hour.rescale), y=mean),color="blue")+
  geom_jitter(
    data=sup_count_dat[[2]], 
    aes(x=as.character(hour.rescale), y=mean),color="red")+
  scale_x_discrete(labels = paste0(c(0,24,72,96), "hr"))+
  theme_bw()+
  adoniram+
  geom_abline(slope = coef(ALTert.lm)[["hour.rescale"]], 
              intercept = coef(ALTert.lm)[["(Intercept)"]], color="blue")+
  geom_abline(slope = coef(SOLert.lm)[["hour.rescale"]], 
              intercept = coef(SOLert.lm)[["(Intercept)"]], color="red")+
  annotate("text", label="ALT: UVA, RH=90%, T=15ºC", x=0.5, y=0.25, hjust=0, color="blue")+
  annotate("text", label="SOL: UVA, RH=90%, T=20ºC", x=0.5, y=0.2, hjust=0, color="red")+
  theme(axis.title.y = element_text(angle=90, margin=margin(r=10)))+
  labs(y="Proportion Germinated Spores")


# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# █░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█
# █░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█
# █░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
###################################################################
#                BOXPLOT Grid for Suppls                          #
###################################################################

# BOXPLOT Grid for Suppls
###################
mytheme<-theme(legend.position = "none", 
               axis.title.x=element_blank(),
               axis.title.y =element_blank(),
               plot.margin = unit(c(.1, .1, .1,.1), "cm"),
               panel.background = element_blank(),
               panel.grid.minor.y = element_line(),
               panel.grid.major = element_line(),
               plot.background = element_rect(),
               axis.text.x = element_text(size=7),
               axis.text.y = element_text(size=7))


altgrid<-vector("list")
x<-1
while(x<length(unique(box_data[[1]]$experiment))+1){
  altgrid[[x]]<-ggplot(box_data[[1]][box_data[[1]]$hour.rescale %in% c(1,4) & 
                                            box_data[[1]]$experiment %in% unique(box_data[[1]]$experiment)[[x]],], 
                       aes(x=as.character(hour.rescale), y=mean, fill=lambda))+
    stat_smooth(method = "lm", se=F, position = position_dodge2(width=2),
                linetype="dotted",lwd=0.5,aes(group=lambda, color=lambda))+
    geom_boxplot(outlier.shape = NA, lwd=.3)+
    theme_bw()+mytheme+ylim(0,1)+
    scale_color_manual(breaks=c("1","4"), labels=c("24hrs", "96hrs"), values=c("#A0A0A0","#6666ff","#33337f"))+
    scale_fill_manual(values=c("#A0A0A0","#6666ff","#33337f"))+
    scale_x_discrete(labels=c("1"="24hrs","4"="96hrs"))
  x<-x+1
}
names(altgrid)<-unique(box_data[[1]]$experiment)




solgridCL<-vector("list")
x<-1
while(x<length(unique(box_data[[2]]$experiment))+1){
  
  solgridCL[[x]]<-ggplot(box_data[[2]][box_data[[2]]$hour.rescale %in% c(1,4) & 
                                              box_data[[2]]$experiment %in% unique(box_data[[2]]$experiment)[[x]],], 
                         aes(x=as.character(hour.rescale), y=mean, fill=lambda))+
    stat_smooth(method = "lm", se=F, position = position_dodge2(width=2),
                linetype="dotted",lwd=0.5,aes(group=lambda, color=lambda))+
    geom_boxplot(outlier.shape = NA, lwd=.3)+
    theme_bw()+mytheme+ylim(0,1)+
    scale_color_manual(breaks=c("1","4"), labels=c("24hrs", "96hrs"), values=c("#A0A0A0","#ff3232","#661414"))+
    scale_fill_manual(values=c("#A0A0A0","#ff3232","#661414"))+
    scale_x_discrete(labels=c("1"="24hrs","4"="96hrs"))
  x<-x+1
}

names(solgridCL)<-unique(box_data[[2]]$experiment)





ALTblocksub<-box_data[[1]]
ALTblocksub$lambda <- factor(ALTblocksub$lambda, levels=c("UVA", "UVB", "b21"), labels=c("UVA", "UVB", "Dark"))
ALTgridleg<-g_legend(ggplot(ALTblocksub[ALTblocksub$hour.rescale %in% c(1,4) & 
                                          ALTblocksub$experiment %in% c("Exp1.1"),], 
                            aes(x=as.character(hour.rescale), y=mean, fill=lambda))+
                       geom_boxplot(outlier.shape = NA, lwd=.3)+
                       labs(fill="Lightsource")+theme(legend.key=element_blank())+
                       scale_fill_manual(breaks=c("UVA","UVB","Dark"), labels=c("UVA","UVB","Dark"),values=c("#6666ff","#33337f","#A0A0A0")))


SOLblocksub<-box_data[[2]]
SOLblocksub$lambda <- factor(SOLblocksub$lambda, levels=c("UVA", "UVB", "b21"), labels=c("UVA", "UVB", "Dark"))
SOLgridleg<-g_legend(ggplot(SOLblocksub[SOLblocksub$hour.rescale %in% c(1,4) & 
                                          SOLblocksub$experiment %in% c("Exp1.1"),], 
                            aes(x=as.character(hour.rescale), y=mean, fill=lambda))+
                       geom_boxplot(outlier.shape = NA, lwd=.3)+
                       labs(fill="Lightsource")+theme(legend.key=element_blank())+
                       scale_fill_manual(breaks=c("UVA","UVB","Dark"), labels=c("UVA","UVB","Dark"),values=c("#ff3232","#661414","#A0A0A0")))


lay3<-rbind(t(matrix(rep(c(12,rep(21,15),rep(7,5)),5), nrow=21, ncol=5)),
            t(matrix(rep(c(13,rep(1,5),rep(23,5),rep(4,5),rep(8,5)),5), nrow=21, ncol=5)),
            t(matrix(rep(c(14,rep(2,5),rep(22,5),rep(5,5),rep(9,5)),5), nrow=21, ncol=5)),
            t(matrix(rep(c(15,rep(20,5),rep(3,5),rep(6,5),rep(10,5)),5), nrow=21, ncol=5)),
            t(matrix(rep(c(11,rep(16,5),rep(17,5),rep(18,5),rep(19,5)),1), nrow=21, ncol=1))
)

sup.alt<-grid.arrange(
  altgrid$Exp1.1,
  altgrid$Exp2.1,
  altgrid$Exp3,
  altgrid$Exp4,
  altgrid$Exp5,
  altgrid$Exp6,
  altgrid$Exp7,
  altgrid$Exp8,
  altgrid$Exp9,
  altgrid$Exp10,
  grobTree(rectGrob(gp=gpar(col=NA,fill="NA", alpha=1))),
  grobTree(rectGrob(gp=gpar(col="white",fill="#a3a3ff", alpha=1,lwd=5)), textGrob("25ºC", rot=90)),
  grobTree(rectGrob(gp=gpar(col="white",fill="#a3a3ff", alpha=1,lwd=5)), textGrob("20ºC", rot=90)),
  grobTree(rectGrob(gp=gpar(col="white",fill="#a3a3ff", alpha=1,lwd=5)), textGrob("15ºC", rot=90)),
  grobTree(rectGrob(gp=gpar(col="white",fill="#a3a3ff", alpha=1,lwd=5)), textGrob("10ºC", rot=90)),
  grobTree(rectGrob(gp=gpar(col="white",fill="#a3a3ff", alpha=1,lwd=5)), textGrob("50%", rot=0)),
  grobTree(rectGrob(gp=gpar(col="white",fill="#a3a3ff", alpha=1,lwd=5)), textGrob("60%", rot=0)),
  grobTree(rectGrob(gp=gpar(col="white",fill="#a3a3ff", alpha=1,lwd=5)), textGrob("75%", rot=0)),
  grobTree(rectGrob(gp=gpar(col="white",fill="#a3a3ff", alpha=1,lwd=5)), textGrob("90%", rot=0)),
  grobTree(rectGrob(gp=gpar(col=NA,fill="white", alpha=1)), ALTgridleg),
  grobTree(rectGrob(gp=gpar(col=NA,fill="white", alpha=1)), textGrob(expression(italic("A. alternata")),gp=gpar(cex=2))),
  grobTree(rectGrob(gp=gpar(col=NA,fill="white", alpha=1))),
  grobTree(rectGrob(gp=gpar(col=NA,fill="white", alpha=1))),
  layout_matrix=lay3)


sup.solCL<-grid.arrange(
  solgridCL$Exp1.1,
  solgridCL$Exp2.1,
  solgridCL$Exp3,
  solgridCL$Exp4,
  solgridCL$Exp5,
  solgridCL$Exp6,
  solgridCL$Exp7,
  solgridCL$Exp8,
  solgridCL$Exp9,
  solgridCL$Exp10,
  rectGrob(gp=gpar(fill=NA, lwd=0)),
  grobTree(rectGrob(gp=gpar(col="white",fill="#ff8484", alpha=1,lwd=5)), textGrob("25ºC", rot=90)),
  grobTree(rectGrob(gp=gpar(col="white",fill="#ff8484", alpha=1,lwd=5)), textGrob("20ºC", rot=90)),
  grobTree(rectGrob(gp=gpar(col="white",fill="#ff8484", alpha=1,lwd=5)), textGrob("15ºC", rot=90)),
  grobTree(rectGrob(gp=gpar(col="white",fill="#ff8484", alpha=1,lwd=5)), textGrob("10ºC", rot=90)),
  grobTree(rectGrob(gp=gpar(col="white",fill="#ff8484", alpha=1,lwd=5)), textGrob("50%", rot=0)),
  grobTree(rectGrob(gp=gpar(col="white",fill="#ff8484", alpha=1,lwd=5)), textGrob("60%", rot=0)),
  grobTree(rectGrob(gp=gpar(col="white",fill="#ff8484", alpha=1,lwd=5)), textGrob("75%", rot=0)),
  grobTree(rectGrob(gp=gpar(col="white",fill="#ff8484", alpha=1,lwd=5)), textGrob("90%", rot=0)),
  grobTree(rectGrob(gp=gpar(col=NA,fill="white", alpha=1)), SOLgridleg),
  grobTree(rectGrob(gp=gpar(col=NA,fill="white", alpha=1)), textGrob(expression(italic("A. solani")),gp=gpar(cex=2))),
  grobTree(rectGrob(gp=gpar(col=NA,fill="white", alpha=1))),
  grobTree(rectGrob(gp=gpar(col=NA,fill="white", alpha=1))),
  layout_matrix=lay3)





# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# █░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█
# █░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█
# █░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# █░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█
# █░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█
# █░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# █░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█
# █░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█
# █░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████

## SECTION B: get best models and do post hoc tests
# Best and sensible models
# post hoc analyses
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# █░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█
# █░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█
# █░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ******************************************************** Part I

# Here are the models already performed
## Note there are two SOL models as described in the methods (second best selected to match ALT model)
MOMA<-readRDS("Models_PartI.rds")
# ALT...MOMA[[1]] >>> AliveSpores ~ hour.rescale + lambda + RH + temp + (1 | experiment/cart:base) +  offset(log(TotalSpores))
# SOL...MOMA[[2]] >>> AliveSpores ~ hour.rescale + lambda + RH + temp + (1 | experiment/cart:base) +  offset(log(TotalSpores))
# SOL...MOMA[[3]] >>> AliveSpores ~ hour.rescale + lambda + RH + temp + height + (1 | experiment/cart:base) + offset(log(TotalSpores))

MOMA90<-readRDS("Models_PartI_90RH.rds")
# ALT...MOMA90[[1]] >>> AliveSpores ~ hour.rescale + lambda + temp + (1 | experiment/cart:base) + offset(log(TotalSpores))
# SOL...MOMA90[[2]] >>> AliveSpores ~ hour.rescale + lambda + temp + (1 | experiment/cart:base) + offset(log(TotalSpores))
# SOL...MOMA90[[3]] >>> AliveSpores ~ hour.rescale + lambda + temp + height + (1 | experiment/cart:base) + offset(log(TotalSpores))



glht_glmmTMB <- function (model, ..., component="cond") {
  glht(model, ...,
       coef. = function(x) fixef(x)[[component]],
       vcov. = function(x) vcov(x)[[component]],
       df = NULL)
}




tem90<-list()
lam90<-list()
hig90<-list()
for (f in 1:length(MOMA90)){
  tryCatch({
    tem90[[f]]<-glht_glmmTMB(MOMA90[[f]], linfct = mcp(temp = "Tukey"))
    lam90[[f]]<-glht_glmmTMB(MOMA90[[f]], linfct = mcp(lambda = "Tukey"))
    hig90[[f]]<-glht_glmmTMB(MOMA90[[f]], linfct = mcp(height = "Tukey"))
  },
  error=function(e){})
}



hum<-list()
tem<-list()
lam<-list()
hig<-list()
for (f in 1:length(MOMA)){
  tryCatch({
    hum[[f]]<-glht_glmmTMB(MOMA[[f]], linfct = mcp(RH = "Tukey"))
    tem[[f]]<-glht_glmmTMB(MOMA[[f]], linfct = mcp(temp = "Tukey"))
    lam[[f]]<-glht_glmmTMB(MOMA[[f]], linfct = mcp(lambda = "Tukey"))
    hig[[f]]<-glht_glmmTMB(MOMA[[f]], linfct = mcp(height = "Tukey"))
  },
  error=function(e){})
}

# _ _ _ _ _ _ _ _ _ _ _ _ _
##### MODEL TABLES ##### MODEL TABLES ##### MODEL TABLES - Part 1

EEE.alt<-as.data.frame(summary(MOMA[[1]])$coef$cond)
EEE.alt$exp_trans<-exp(EEE.alt$Estimate)

EEE.solcl<-as.data.frame(summary(MOMA[[2]])$coef$cond)
EEE.solcl$exp_trans<-exp(EEE.solcl$Estimate)

# _ _ _ _ _ _ _ _ _ _ _ 

EEE90.alt<-as.data.frame(summary(MOMA90[[1]])$coef$cond)
EEE90.alt$exp_trans<-exp(EEE90.alt$Estimate)

EEE90.solcl<-as.data.frame(summary(MOMA90[[2]])$coef$cond)
EEE90.solcl$exp_trans<-exp(EEE90.solcl$Estimate)

# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# █░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█
# █░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█
# █░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
#####################################################################
################             KW-Dunn             ####################

spore_means<-readRDS("spore_count_means.rds")

ALTmean_s<-spore_means[[1]][spore_means[[1]]$hour.rescale ==4,]
SOLCLmean_s<-spore_means[[2]][spore_means[[2]]$hour.rescale ==4,]

ALTmean_s$AvD<-ALTmean_s$AliveSpores/ALTmean_s$TotalSpores
SOLCLmean_s$AvD<-SOLCLmean_s$AliveSpores/SOLCLmean_s$TotalSpores

arabia<-list(
  ALTmean_s,
  SOLCLmean_s)

kruskal.test(AvD ~ lambda, 
             data = ALTmean_s)
dunn.test(ALTmean_s$AvD,
          ALTmean_s$lambda, method="bh", list=TRUE)


kruskal.test(AvD ~ lambda, 
             data = SOLCLmean_s)
dunn.test(SOLCLmean_s$AvD,
          SOLCLmean_s$lambda, method="bh", list=TRUE)


#BY T #BY T #BY T #BY T #BY T #BY T #BY T #BY T #BY T #BY T #BY T
LLL<-c("UVA","UVB","b21")
KW.t<-list()
DUNN.t<-list()
for (v in 1:length(arabia)){
  KW.t[[v]]<-list()
  DUNN.t[[v]]<-list()
  for (u in 1:length(LLL)){
    KW.t[[v]][[u]]<-kruskal.test(AvD ~ temp, 
                                 data = arabia[[v]][arabia[[v]]$lambda == LLL[u],])
    DUNN.t[[v]][[u]]<-dunn.test(arabia[[v]][arabia[[v]]$lambda == LLL[u],]$AvD,
                                arabia[[v]][arabia[[v]]$lambda == LLL[u],]$temp, method="bh", list=TRUE)
  }
}


KW.t90<-list()
DUNN.t90<-list()
for (v in 1:length(arabia)){
  KW.t90[[v]]<-list()
  DUNN.t90[[v]]<-list()
  for (u in 1:length(LLL)){
    KW.t90[[v]][[u]]<-kruskal.test(AvD ~ temp, 
                                   data = arabia[[v]][arabia[[v]]$lambda == LLL[u] & arabia[[v]]$RH == 90,])
    DUNN.t90[[v]][[u]]<-dunn.test(arabia[[v]][arabia[[v]]$lambda == LLL[u] & arabia[[v]]$RH == 90,]$AvD,
                                  arabia[[v]][arabia[[v]]$lambda == LLL[u] & arabia[[v]]$RH == 90,]$temp, method="bh", list=TRUE)
  }
}

#BY RH #BY RH #BY RH #BY RH #BY RH #BY RH #BY RH #BY RH #BY RH #BY RH
LLL<-c("UVA","UVB","b21")
KW.rh<-list()
DUNN.rh<-list()
for (v in 1:length(arabia)){
  KW.rh[[v]]<-list()
  DUNN.rh[[v]]<-list()
  for (u in 1:length(LLL)){
    KW.rh[[v]][[u]]<-kruskal.test(AvD ~ RH, 
                                  data = arabia[[v]][arabia[[v]]$lambda == LLL[u],])
    DUNN.rh[[v]][[u]]<-dunn.test(arabia[[v]][arabia[[v]]$lambda == LLL[u],]$AvD,
                                 arabia[[v]][arabia[[v]]$lambda == LLL[u],]$RH, method="bh", list=TRUE)
  }
}


#BY EXP BY EXP BY EXP BY EXP BY EXP BY EXP BY EXP BY EXP BY EXP BY EXP
LLL<-c("UVA","UVB","b21")
KW.ex<-list()
DUNN.ex<-list()
for (v in 1:length(arabia)){
  KW.ex[[v]]<-list()
  DUNN.ex[[v]]<-list()
  for (u in 1:length(LLL)){
    KW.ex[[v]][[u]]<-kruskal.test(AvD ~ experiment, 
                                  data = arabia[[v]][arabia[[v]]$lambda == LLL[u],])
    DUNN.ex[[v]][[u]]<-dunn.test(arabia[[v]][arabia[[v]]$lambda == LLL[u],]$AvD,
                                 arabia[[v]][arabia[[v]]$lambda == LLL[u],]$experiment, method="bh", list=TRUE)
  }
}

# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# █░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█
# █░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█
# █░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
######################     PUT POSTHOC TOGETHER   ##############################
HARIF<-data.frame(
  comp=DUNN.t[[1]][[1]]$comparisons,
  p.aj.a.a=DUNN.t[[1]][[1]]$P.adjusted,
  p.aj.b.a=DUNN.t[[1]][[2]]$P.adjusted,
  p.aj.d.a=DUNN.t[[1]][[3]]$P.adjusted,
  p.aj.a.c=DUNN.t[[2]][[1]]$P.adjusted,
  p.aj.b.c=DUNN.t[[2]][[2]]$P.adjusted,
  p.aj.d.c=DUNN.t[[2]][[3]]$P.adjusted
)


HARIF90<-data.frame(
  comp=DUNN.t90[[1]][[1]]$comparisons,
  p.aj.a.a=DUNN.t90[[1]][[1]]$P.adjusted,
  p.aj.b.a=DUNN.t90[[1]][[2]]$P.adjusted,
  p.aj.d.a=DUNN.t90[[1]][[3]]$P.adjusted,
  p.aj.a.c=DUNN.t90[[2]][[1]]$P.adjusted,
  p.aj.b.c=DUNN.t90[[2]][[2]]$P.adjusted,
  p.aj.d.c=DUNN.t90[[2]][[3]]$P.adjusted
)

RATUV<-data.frame(
  comp=DUNN.rh[[1]][[1]]$comparisons,
  p.aj.a.a=DUNN.rh[[1]][[1]]$P.adjusted,
  p.aj.b.a=DUNN.rh[[1]][[2]]$P.adjusted,
  p.aj.d.a=DUNN.rh[[1]][[3]]$P.adjusted,
  p.aj.a.c=DUNN.rh[[2]][[1]]$P.adjusted,
  p.aj.b.c=DUNN.rh[[2]][[2]]$P.adjusted,
  p.aj.d.c=DUNN.rh[[2]][[3]]$P.adjusted
)

EXP<-data.frame(
  comp=DUNN.ex[[1]][[1]]$comparisons,
  p.aj.a.a=DUNN.ex[[1]][[1]]$P.adjusted,
  p.aj.b.a=DUNN.ex[[1]][[2]]$P.adjusted,
  p.aj.d.a=DUNN.ex[[1]][[3]]$P.adjusted,
  p.aj.a.c=DUNN.ex[[2]][[1]]$P.adjusted,
  p.aj.b.c=DUNN.ex[[2]][[2]]$P.adjusted,
  p.aj.d.c=DUNN.ex[[2]][[3]]$P.adjusted
)

# BIG OLE MATRIX OF PPPPP for EXP

EXPunge<-as.data.frame(cbind(str_split_fixed(EXP$comp, " - ", 2),EXP[,-1]))
colnames(EXPunge)[c(1,2)]<-c("X","Y")
EXPunge$X<-as.character(EXPunge$X)
EXPunge$Y<-as.character(EXPunge$Y)
EXPunge[,1:2] <- data.frame(lapply(EXPunge[,1:2], function(x) {
  gsub(".1","",gsub("Exp", "", x))
}))
EXPunge$indx<-paste(paste("Exp",EXPunge$X,".1",sep=""),
                    paste("Exp",EXPunge$Y,sep=""),
                    sep=":")
EXPunge$indx2<-paste(paste("Exp",EXPunge$X,".2",sep=""),
                     paste("Exp",EXPunge$Y,sep=""),
                     sep=":")
EXPunge$indx3<-paste(paste("Exp",EXPunge$X,".3",sep=""),
                     paste("Exp",EXPunge$Y,sep=""),
                     sep=":")

df<-t(data.frame(
  row.names = c(rbind(paste("Exp",seq(1.1, 10.1,by=1),sep=""),
                      paste("Exp",seq(1.2, 10.2,by=1),sep=""),
                      paste("Exp",seq(1.3, 10.3,by=1),sep=""))),
  Exp1=rep("A",30),
  Exp2=rep("A",30),
  Exp3=rep("A",30),
  Exp4=rep("A",30),
  Exp5=rep("A",30),
  Exp6=rep("A",30),
  Exp7=rep("A",30),
  Exp8=rep("A",30),
  Exp9=rep("A",30),
  Exp10=rep("A",30)
))
df[] <- paste(col(df, TRUE), row(df, TRUE), sep = ":")

df.a<-df
df.s<-df
alter<-colnames(EXPunge[,3:5])
soler<-colnames(EXPunge[,6:8])


for(x in 1:nrow(EXPunge)){
  test<-paste(paste("Exp",EXPunge$X[x],".1",sep=""),
              paste("Exp",EXPunge$Y[x],sep=""),
              sep=":")
  test2<-paste(paste("Exp",EXPunge$X[x],".2",sep=""),
               paste("Exp",EXPunge$Y[x],sep=""),
               sep=":")
  test3<-paste(paste("Exp",EXPunge$X[x],".3",sep=""),
               paste("Exp",EXPunge$Y[x],sep=""),
               sep=":")
  
  df.a[which(df.a==test, arr.ind=TRUE)[1],which(df.a==test, arr.ind=TRUE)[2] ]<-formatC( EXPunge[EXPunge$indx == test,alter[1]],format = "e", digits = 2)
  df.s[which(df.s==test, arr.ind=TRUE)[1],which(df.s==test, arr.ind=TRUE)[2] ]<-formatC( EXPunge[EXPunge$indx == test,soler[1]],format = "e", digits = 2)
  
  df.a[which(df.a==test2, arr.ind=TRUE)[1],which(df.a==test2, arr.ind=TRUE)[2] ]<-formatC( EXPunge[EXPunge$indx2 == test2,alter[2]],format = "e", digits = 2)
  df.s[which(df.s==test2, arr.ind=TRUE)[1],which(df.s==test2, arr.ind=TRUE)[2] ]<-formatC( EXPunge[EXPunge$indx2 == test2,soler[2]],format = "e", digits = 2)
  
  df.a[which(df.a==test3, arr.ind=TRUE)[1],which(df.a==test3, arr.ind=TRUE)[2] ]<-formatC( EXPunge[EXPunge$indx3 == test3,alter[3]],format = "e", digits = 2)
  df.s[which(df.s==test3, arr.ind=TRUE)[1],which(df.s==test3, arr.ind=TRUE)[2] ]<-formatC( EXPunge[EXPunge$indx3 == test3,soler[3]],format = "e", digits = 2)
}


# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# █░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█
# █░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█
# █░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ******************************************************** Part II

modII <- readRDS("Model_PartII.rds")
summary(modII)
glht_glmmTMB(modII, linfct = mcp(SPECIES = "Tukey"))
summary(glht_glmmTMB(modII, linfct = mcp(SPECIES = "Tukey")))
parte2<-as.data.frame(summary(modII)$coef$cond)
parte2$exp_trans<-exp(parte2$Estimate)
plot_model(modII,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")



# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# █░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█
# █░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█
# █░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
######################################################################################################
# Compile Tables for Pub ## Compile Tables for Pub ## Compile Tables for Pub ## Compile Tables for Pub #
######################################################################################################


A = matrix(c(2, 4, 3, 1, 5, 7), nrow=2, ncol=3,  byrow = TRUE)        # fill matrix by rows 
B = matrix(c(89, 84, 23, 11, 58, 27), nrow=2, ncol=3,  byrow = TRUE)        # fill matrix by rows 
matrix(paste(A, B, sep = "|"), nrow=2, ncol=3, byrow=T)

## Summary of best models ## Summary of best models ## Summary of best models ## Summary of best models
## Summary of best models ## Summary of best models ## Summary of best models ## Summary of best models
## Summary of best models ## Summary of best models ## Summary of best models ## Summary of best models

# Comapre some of the better odels for ALT, SOL, ALT:RH=90, and SOL:RH=90


MOODOO<-readRDS("compare_top_models.PartI.rds")
names(MOODOO)

#ALT 
cbind(
  data.frame(dAIC=AICtab(MOODOO[[1]],
                         MOODOO[[2]],
                         MOODOO[[3]])),
  
  params=c(names(MOODOO[1]),
           names(MOODOO[2]),
           names(MOODOO[3]))
)

#SOL cl
cbind(
  data.frame(dAIC=AICtab(MOODOO[[4]],
                         MOODOO[[5]],
                         MOODOO[[6]])),
  
  params=c(names(MOODOO[4]),
           names(MOODOO[5]),
           names(MOODOO[6]))
)
# $ $ $ $  90 
#ALT 
cbind(
  data.frame(dAIC=AICtab(MOODOO[[7]],
                         MOODOO[[8]],
                         MOODOO[[9]])),
  
  params=c(names(MOODOO[7]),
           names(MOODOO[8]),
           names(MOODOO[9]))
)
#SOL cl
cbind(
  data.frame(dAIC=AICtab(MOODOO[[10]],
                         MOODOO[[11]],
                         MOODOO[[12]])),
  
  params=c(names(MOODOO[10]),
           names(MOODOO[11]),
           names(MOODOO[12]))
)

# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# █░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█
# █░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█
# █░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
#### TABLE OF CHOSEN MODEL SUPPL TABLE OF CHOSEN MODEL SUPPL TABLE OF CHOSEN MODEL SUPPL TABLE OF CHOSEN MODEL
## {FIGURE 4}
ukelele<-theme(
  plot.background = element_blank(),
  panel.background = element_rect(fill = "white", colour = "grey50"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)

gore<-gridExtra::grid.arrange(
  plot_model(MOODOO[[1]],
             show.values = TRUE, 
             value.offset = .4, 
             dot.size = 1, 
             value.size = 3,
             transform="exp", 
             vline.color = "grey", 
             type="est",
             axis.labels = c("25º","20º","15º","90%","75%","60%","UVB","UVA","Day"),
             title="",
             width=.2,
             line.size=.2,
             colors = c("#000080", "#306EFF"))+
    ukelele+
    labs(y="Exponentiated Effect Size")+ylim(0.02,2),
  
  plot_model(MOODOO[[5]],
             show.values = TRUE, 
             value.offset = .4, 
             dot.size = 1, 
             value.size = 3,
             transform="exp", 
             vline.color = "grey", 
             type="est",
             axis.labels = c("25º","20º","15º","90%","75%","60%","UVB","UVA","Day"),
             title="",
             width=.2,
             line.size=.2,
             colors = c("#FF0000", "#FF7F50"))+
    ukelele+
    labs(y="Exponentiated Effect Size")+ylim(0.02,2),
  
  
  plot_model(MOODOO[[7]],
             show.values = TRUE, 
             value.offset = .4, 
             dot.size = 1, 
             value.size = 3,
             transform="exp", 
             vline.color = "grey", 
             type="est",
             axis.labels = c("25º","20º","15º","UVB","UVA","Day"),
             title="",
             width=.2,
             line.size=.2,
             colors = c("#000080", "#306EFF"))+
    ukelele+
    labs(y="Exponentiated Effect Size")+ylim(0.02,2),
  
  
  plot_model(MOODOO[[11]],
             show.values = TRUE, 
             value.offset = .4, 
             dot.size = 1, 
             value.size = 3,
             transform="exp", 
             vline.color = "grey", 
             type="est",
             axis.labels = c("25º","20º","15º","UVB","UVA","Day"),
             title="",
             width=.2,
             line.size=.2,
             colors = c("#FF0000", "#FF7F50"))+
    ukelele+
    labs(y="Exponentiated Effect Size")+ylim(0.02,2),
  ncol=2, nrow=2)


# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# █░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█
# █░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█
# █░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
#### POSHOC POSHOC POSHOC POSHOC POSHOC POSHOC POSHOC POSHOC POSHOC POSHOC POSHOC POSHOC POSHOC POSHOC 

coefp1<-data.frame(coef=round(summary(tem[[1]])[["test"]][["coefficients"]],digits=2),
                   pvaj=round(summary(tem[[1]])[["test"]][["pvalues"]],digits=2))
coefp1$fancy<-paste(coefp1[,1],paste("(",coefp1[,2],")",sep=""), sep=" ")

coefp3<-data.frame(coef=round(summary(tem[[3]])[["test"]][["coefficients"]],digits=2),
                   pvaj=round(summary(tem[[3]])[["test"]][["pvalues"]],digits=2))
coefp3$fancy<-paste(coefp3[,1],paste("(",coefp3[,2],")",sep=""), sep=" ")

maria<-c("10","15","20","25")
mother<-matrix(levels(interaction(maria,maria,sep=' - ')), nrow=4,ncol=4,byrow = T)
colnames(mother)=rownames(mother)=maria


tem.tri.pre1<-matrix(coefp1$fancy[match(mother,rownames(coefp1))],nrow=4)
tem.tri.pre3<-matrix(coefp3$fancy[match(mother,rownames(coefp3))],nrow=4)
tem.tri<-matrix(paste(tem.tri.pre1, tem.tri.pre3, sep=" | "),nrow=4)
colnames(tem.tri.pre1)=colnames(tem.tri.pre3)=row.names(tem.tri.pre1)=row.names(tem.tri.pre1)=maria

new.tem <- matrix(NA, nrow = 4, ncol = 4)
new.tem[upper.tri(new.tem)] <- tem.tri.pre1[upper.tri(tem.tri.pre1)]
new.tem[lower.tri(new.tem)] <- tem.tri.pre3[upper.tri(tem.tri.pre3)]
rownames(new.tem)=colnames(new.tem)=maria


# RH RH
cohum1<-data.frame(coef=round(summary(hum[[1]])[["test"]][["coefficients"]],digits=2),
                   pvaj=round(summary(hum[[1]])[["test"]][["pvalues"]],digits=2))
cohum1$fancy<-paste(cohum1[,1],paste("(",cohum1[,2],")",sep=""), sep=" ")

cohum3<-data.frame(coef=round(summary(hum[[3]])[["test"]][["coefficients"]],digits=2),
                   pvaj=round(summary(hum[[3]])[["test"]][["pvalues"]],digits=2))
cohum3$fancy<-paste(cohum3[,1],paste("(",cohum3[,2],")",sep=""), sep=" ")

susan<-c("50","60","75","90")
father<-matrix(levels(interaction(susan,susan,sep=' - ')), nrow=4,ncol=4,byrow = T)
colnames(father)=rownames(father)=susan

hum.tri.pre1<-matrix(cohum1$fancy[match(father,rownames(cohum1))],nrow=4)
hum.tri.pre3<-matrix(cohum3$fancy[match(father,rownames(cohum3))],nrow=4)
hum.tri<-matrix(paste(hum.tri.pre1, hum.tri.pre3, sep=" | "),nrow=4)
colnames(hum.tri.pre1)=row.names(hum.tri.pre1)=colnames(hum.tri.pre3)=row.names(hum.tri.pre3)=susan

new.hum <- matrix(NA, nrow = 4, ncol = 4)
new.hum[upper.tri(new.hum)] <- hum.tri.pre1[upper.tri(hum.tri.pre1)]
new.hum[lower.tri(new.hum)] <- hum.tri.pre3[upper.tri(hum.tri.pre3)]
rownames(new.hum)=colnames(new.hum)=susan

# TO GO TO GO

roses<-matrix(paste(new.hum, new.tem, sep=" | "), nrow=4, ncol=4, byrow = F)
#write.csv(roses, "MULCOMP.csv")

newb <- matrix(NA, nrow = 4, ncol = 4)
newb[upper.tri(newb)] <- tem.tri[upper.tri(tem.tri)]
newb[lower.tri(newb)] <- hum.tri[upper.tri(hum.tri)]
#write.csv(newb, "MULCOMP_tri.csv")


# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# █░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█
# █░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█░░▄▀▄▀▄▀▄▀▄▀░░█
# █░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█░░░░░░░░░░░░░░█
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████
# ██████████████████████████████████████████████████████████████████████████████████████████████████████████

################################################################################
#  Part II box plot in Figure 5

SUMMA.sumsum<-readRDS("PartII_data.summary.rds")

vindication<-grid.arrange(
  ggplot(SUMMA.sumsum[SUMMA.sumsum$SPECIES %in% "ALT" & SUMMA.sumsum$IRR.cat %in% "LIGHT",],
         aes(x=as.factor(HOUR), y=meansum))+geom_bar(stat="identity", fill="#6666ff")+theme_bw()+margot+
    geom_errorbar(aes(ymin=meansum-sebysum, ymax=meansum+sebysum), width=.2,
                  position=position_dodge(.9)) + scale_x_discrete(labels=c("0hrs","24hrs","72hrs","144hrs","216hrs","288hrs"))+
    ylim(0,1)+
    annotation_custom(ggplotGrob(
      ggplot(SUMMA.sumsum[SUMMA.sumsum$SPECIES %in% "ALT" & SUMMA.sumsum$IRR.cat %in% "DARK",],
             aes(x=as.factor(HOUR), y=meansum))+geom_bar(stat="identity", fill="darkgrey")+theme_bw()+margot+
        scale_x_discrete(labels=c("0hrs","24hrs","72hrs","144hrs","216hrs","288hrs"))+
        theme(panel.background = element_rect(fill = "white"),
              plot.background = element_rect(
                fill = "white",
                colour = "black",
                size = 0),
              axis.text.x = element_text(size=8, angle=45, hjust=1, vjust=1.1),
              axis.text.y = element_text(size=8, angle=0))+
        geom_errorbar(aes(ymin=meansum-sebysum, ymax=meansum+sebysum), width=.2,
                      position=position_dodge(.9)) + ylim(0,1)   ), 
      xmin = 4.5, xmax = 6.5, ymin = 0.75,  ymax = 1.04),
  
  
  
  ggplot(SUMMA.sumsum[SUMMA.sumsum$SPECIES %in% "SOL" & SUMMA.sumsum$IRR.cat %in% "LIGHT",],
         aes(x=as.factor(HOUR), y=meansum))+geom_bar(stat="identity", fill="#ff4c4c")+theme_bw()+margot+
    geom_errorbar(aes(ymin=meansum-sebysum, ymax=meansum+sebysum), width=.2,
                  position=position_dodge(.9)) + scale_x_discrete(labels=c("0hrs","24hrs","72hrs","144hrs","216hrs","288hrs"))+
    ylim(0,1)+
    annotation_custom(ggplotGrob(
      ggplot(SUMMA.sumsum[SUMMA.sumsum$SPECIES %in% "SOL" & SUMMA.sumsum$IRR.cat %in% "DARK",],
             aes(x=as.factor(HOUR), y=meansum))+geom_bar(stat="identity", fill="darkgrey")+theme_bw()+margot+
        scale_x_discrete(labels=c("0hrs","24hrs","72hrs","144hrs","216hrs","288hrs"))+
        theme(panel.background = element_rect(fill = "white"),
              plot.background = element_rect(
                fill = "white",
                colour = "black",
                size = 0),
              axis.text.x = element_text(size=8, angle=45, hjust=1, vjust=1.1),
              axis.text.y = element_text(size=8, angle=0))+
        geom_errorbar(aes(ymin=meansum-sebysum, ymax=meansum+sebysum), width=.2,
                      position=position_dodge(.9)) + ylim(0,1)   ), 
      xmin = 4.5, xmax = 6.5, ymin = 0.75,  ymax = 1.04),
  
  nrow=1,ncol=2)



