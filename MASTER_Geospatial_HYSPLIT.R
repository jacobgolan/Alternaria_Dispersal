library(data.table)
library(parallel)
library(future.apply)
library(dplyr)
library(ggmap)
library(alphahull)
library(tidyverse)
library(sp)
library(geosphere)
library(aspace)
library(ggplot2)
library(rgdal)
library(raster)
library(ggalt)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(ggforce)
library(sf)
library(ggrepel)
library(scales)
library(pbapply)
library(RColorBrewer)
library(emmeans)
library(rgeos)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)

## DECOMPRESS HYSPLIT_RAW.zip
# 1  #"ALTE_BH500000.rdsSPEED.csv" 
# 2  ## CLOSURE MODEL NOT USED::: "ALTE_KC500000.rdsSPEED.csv" 
# 3  #"SOLA_BH500000.rdsSPEED.csv" 
# 4  ## CLOSURE MODEL NOT USED::: "SOLA_KC500000.rdsSPEED.csv"

mafuiles<-list.files() #whichever of the above files you wish to analyse...

for (bbb in 1:4){
  mmx<-list()

  jj<-bbb
  mmx<-as.data.frame(fread(mafuiles[jj], verbose=F))
  colnames(mmx)[c(4,5)]<-c("FT","MX")
  mmx$VOO06<-floor(mmx$FT/360)
  mmx$VOO12<-floor(mmx$FT/720)
  mmx$VOO24<-floor(mmx$FT/1440)

  ####################################################################################################
  #################################       CENTROID         ###########################################
  ####################################################################################################

  mods<-list()
  mods[[1]]<-mmx%>%
    group_by(HOUR) %>%
    group_split()
  
  
  bellaciao<-mods[[1]]
  entregalo<-4
  
  pb <- txtProgressBar(min = 0, max = 3 , style = 3)
  centrir<-list()
  centrirkey<-list()
  for (z in 1:length(bellaciao)){
    centrir[[z]] <- bellaciao[[z]] %>%   
      group_by(VOO06) %>%  #YEAR,MONTH,DAY,VOOINT
      group_split() 
    centrirkey[[z]] <- bellaciao[[z]] %>%
      group_by(VOO06) %>% #YEAR,MONTH,DAY,HOUR,VOOINT
      group_keys()
    setTxtProgressBar(pb, z)
  }
 
  for (w in 1:length(centrir)){
    centrir[[w]]<-centrir[[w]][sapply(centrir[[w]],nrow)>3]
  }
  
  jejum<-"+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
  zaza<-"+init=epsg:2288 +proj=laea +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
  
  ramses<-list()
  pwantz<-list()
  pb <- txtProgressBar(min = 0, max = unname(cumsum(lapply(centrir, function(o) length(o)))[3]), style = 3)
  
  for (c in 1:length(centrir)){
    ramses[[c]]<-list()
    pwantz[[c]]<-list()
    for (r in 1:length(centrir[[c]])){
      tmp<-centrir[[c]][[r]]
      coordinates(tmp) <- c("LONG", "LAT")
      proj4string(tmp) <- CRS(zaza) #local projection
      calc_sde(id=paste("kat",tmp@data[1,11],sep="."),
               points=as.data.frame(tmp@coords, verbose=F)) 
      A<-sdeatt
      B<-sdeloc
      ramses[[c]][[r]]<-A
      pwantz[[c]][[r]]<-data.frame(id=B$id,
                                   x=as.numeric(as.character(B$x)),
                                   y=as.numeric(as.character(B$y)))
      rm(sdeatt)
      rm(sdeloc)
      unlink("SDE_Output.txt")
      setTxtProgressBar(pb, c*r)
    }
  }
  
  ############## Mapping
  
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  states<- ne_states(returnclass = "sf",country = c("United States of America","Canada"))
  lakes110 <- ne_download(type = 'lakes', category = 'physical',returnclass = "sf")
  rivers110 <- ne_download(type = 'rivers_lake_centerlines', category = 'physical',returnclass = "sf")
  
  
  epopeia<-list()
  for (i in 1:length(ramses)){
    epopeia[[i]]<-do.call(rbind, ramses[[i]])
  }
  
  junkie<-list()
  for (u in 1:length(epopeia)){
    doodie<-data.frame(long=epopeia[[u]][,"CENTRE.x"],
                       lat =epopeia[[u]][,"CENTRE.y"])
    coordinates(doodie) <- c("long", "lat")
    proj4string(doodie) <- CRS(jejum) # WGS 84 - just the baseline expected coord system
    CRS.new <- CRS(zaza)
    cura <- spTransform(doodie, CRS.new)
    junkie[[u]]<-data.frame(long=cura@coords[,"long"],
                            lat=cura@coords[,"lat"],
                            six=epopeia[[u]]$Sigma.x,
                            siy=epopeia[[u]]$Sigma.y,
                            theta=epopeia[[u]]$Theta,
                            ID=epopeia[[u]]$id)
  }
  # ZOOOOOOM :::: https://github.com/tidyverse/ggplot2/issues/2090
  # the bouding box polygon in long/lat projection, i.e. axis-aligned
  bb <- st_sfc(
    st_polygon(list(cbind(
      c(-100,-10,-100,-10,-100), # x-coordinates (longitudes) of points A,B,C,D,A
      c(20,20,65,65,20)     # y-coordinates (latitudes) of points A,B,C,D,A
    ))),
    crs = jejum)
  
  laeabb <- st_transform(bb, crs = zaza)
  b <- st_bbox(laeabb)
  
  ephestos<-list()
  for (u in 1:length(pwantz)){
    gherbi<-do.call(rbind, pwantz[[u]])
    doodie<-data.frame(long=(gherbi[,"x"]),
                       lat =(gherbi[,"y"]))
    doodie$long<-as.numeric(as.character(doodie$long))
    doodie$lat<- as.numeric(as.character(doodie$lat))
    coordinates(doodie) <- c("long", "lat")
    proj4string(doodie) <- CRS(jejum) # WGS 84 - just the baseline expected coord system
    CRS.new <- CRS(zaza)
    cura <- spTransform(doodie, CRS.new)
    ephestos[[u]]<-data.frame(long=cura@coords[,"long"],
                              lat=cura@coords[,"lat"],
                              ID=gherbi$id,
                              alpha=as.integer(as.factor(gherbi$id)))
  }
  
  qww<-list()
  for (e in 1:length(ephestos)){
    chingamono<-c(4,12,24,36,48)

    if(length(unique(ephestos[[e]]$ID))>=48){
      surv<-unique(ephestos[[e]]$ID)[chingamono]
    } else if(length(unique(ephestos[[e]]$ID))>=36){
      surv<-unique(ephestos[[e]]$ID)[chingamono][c(1:4)]
    } else if(length(unique(ephestos[[e]]$ID))>=24){
      surv<-unique(ephestos[[e]]$ID)[chingamono][c(1:3)]
    } else if(length(unique(ephestos[[e]]$ID))>=12){
      surv<-unique(ephestos[[e]]$ID)[chingamono][c(1:2)]
    } else if(length(unique(ephestos[[e]]$ID))>=4){
      surv<-unique(ephestos[[e]]$ID)[chingamono][c(1)]
    }
    ww<-ephestos[[e]][ephestos[[e]]$ID %in% surv,]
    g<-list()
    t<-list()
    for (i in 1:length(surv)){
      h<-sample(nrow(ww[ww$ID %in% surv[i],]),1)
      g[i]<-ww[h,][1]
      t[i]<-ww[h,][2]
    }
    inoc<-c("24hrs","72hrs","144hrs","216hrs","218hrs")
    labs<-inoc[c(1:length(surv))]#[chingamono[chingamono <= length(unique(ephestos[[e]]$ID))]/4]
    qww[[e]]<-data.frame(row.names = labs[!is.na(labs)],
                         long=unlist(g),
                         lat=unlist(t)
    )
  }

  if (jj<=2){
    # #ALT BLUE
    colfunc <- colorRampPalette(c( "#191970", "#63D1F4"))
    nas<- "#001399" #ALT
  } else if(jj>=3){
    #SOL RED
    colfunc <- colorRampPalette(c( "#B2022F", "#FF3ADB"))
    nas<- "#B2022F" #SOL
  }
  
  n=49
  colfunc(n)

  major<-"#B1B1B1"
  minor<-"#B1B1B1"
  lines<-"#FFFFFF"

  LDIY="12" #line/dot spacing
  TH=0.1    #ellipse line thickness
  X=Y=7 #pdf size

  
  zitti1<-ggplot(data = world) +
    geom_sf(color = lines, fill = minor)+
    geom_sf(data = states, fill = NA, color=lines,lwd=.25)+
    geom_sf(data = lakes110, fill = lines, color=lines, lwd=.25)+
    coord_sf(crs = zaza, xlim = c(b["xmin"], b["xmax"]), ylim = c(b["ymin"], b["ymax"]))+
    geom_polygon(data=ephestos[[1]], aes(x=long, y=lat, fill=ID, alpha=alpha),show.legend = F,lwd=.25)+
    geom_polygon(data=ephestos[[1]][ephestos[[1]]$ID %in% "kat.3",] , aes(x=long, y=lat),color=nas,fill=NA,lty=c(LDIY),show.legend = F,lwd=TH)+
    geom_polygon(data=ephestos[[1]][ephestos[[1]]$ID %in% "kat.11",], aes(x=long, y=lat),color=nas,fill=NA,lty=c(LDIY),show.legend = F,lwd=TH)+
    geom_polygon(data=ephestos[[1]][ephestos[[1]]$ID %in% "kat.23",], aes(x=long, y=lat),color=nas,fill=NA,lty=c(LDIY),show.legend = F,lwd=TH)+
    geom_polygon(data=ephestos[[1]][ephestos[[1]]$ID %in% "kat.35",], aes(x=long, y=lat),color=nas,fill=NA,lty=c(LDIY),show.legend = F,lwd=TH)+
    geom_polygon(data=ephestos[[1]][ephestos[[1]]$ID %in% "kat.47",], aes(x=long, y=lat),color=nas,fill=NA,lty=c(LDIY),show.legend = F,lwd=TH)+
    theme(panel.grid.major = element_line(color = major, linetype = "dashed", size = .25), 
          panel.background = element_rect(fill = "white"), 
          panel.border = element_rect(color="black",fill = NA),
          axis.text=element_text(size=14),
          axis.title=element_text(size=14))+
    scale_alpha(range = c(0.04, 0.04))+
    scale_fill_manual(values=colfunc(n))+
    labs(x="Longitude",y="Latitude")
  
  if (jj==1){
    ggsave("ALT_BH50_hr05.pdf", plot = zitti1, device="pdf",width = X, height = Y, units = c("in"),dpi = 300)
  } else if(jj==2){
    ggsave("ALT_KC50_hr05.pdf", plot = zitti1, device="pdf",width = X, height = Y, units = c("in"),dpi = 300)
  } else if(jj==3){
    ggsave("SOL_BH50_hr05.pdf", plot = zitti1, device="pdf",width = X, height = Y, units = c("in"),dpi = 300)
  } else if(jj==4){
    ggsave("SOL_KC50_hr05.pdf", plot = zitti1, device="pdf",width = X, height = Y, units = c("in"),dpi = 300)
  }
  
  zitti2<-ggplot(data = world) +
    geom_sf(color = lines, fill = minor)+
    geom_sf(data = states, fill = NA, color=lines,lwd=.25)+
    geom_sf(data = lakes110, fill = lines, color=lines, lwd=.25)+
    coord_sf(crs = zaza, xlim = c(b["xmin"], b["xmax"]), ylim = c(b["ymin"], b["ymax"]))+
    geom_polygon(data=ephestos[[2]], aes(x=long, y=lat, fill=ID, alpha=alpha),show.legend = F,lwd=.25)+
    geom_polygon(data=ephestos[[2]][ephestos[[2]]$ID %in% "kat.3",] , aes(x=long, y=lat),color=nas,fill=NA,lty=c(LDIY),show.legend = F,lwd=TH)+
    geom_polygon(data=ephestos[[2]][ephestos[[2]]$ID %in% "kat.11",], aes(x=long, y=lat),color=nas,fill=NA,lty=c(LDIY),show.legend = F,lwd=TH)+
    geom_polygon(data=ephestos[[2]][ephestos[[2]]$ID %in% "kat.23",], aes(x=long, y=lat),color=nas,fill=NA,lty=c(LDIY),show.legend = F,lwd=TH)+
    geom_polygon(data=ephestos[[2]][ephestos[[2]]$ID %in% "kat.35",], aes(x=long, y=lat),color=nas,fill=NA,lty=c(LDIY),show.legend = F,lwd=TH)+
    geom_polygon(data=ephestos[[2]][ephestos[[2]]$ID %in% "kat.47",], aes(x=long, y=lat),color=nas,fill=NA,lty=c(LDIY),show.legend = F,lwd=TH)+
    theme(panel.grid.major = element_line(color = major, linetype = "dashed", size = .25), 
          panel.background = element_rect(fill = "white"), 
          panel.border = element_rect(color="black",fill = NA),
          axis.text=element_text(size=14),
          axis.title=element_text(size=14))+
    scale_alpha(range = c(0.04, 0.04))+
    scale_fill_manual(values=colfunc(n))+
    labs(x="Longitude",y="Latitude")
  
  if (jj==1){
    ggsave("ALT_BH50_hr15.pdf", plot = zitti2, device="pdf",width = X, height = Y, units = c("in"),dpi = 300)
  } else if(jj==2){
    ggsave("ALT_KC50_hr15.pdf", plot = zitti2, device="pdf",width = X, height = Y, units = c("in"),dpi = 300)
  } else if(jj==3){
    ggsave("SOL_BH50_hr15.pdf", plot = zitti2, device="pdf",width = X, height = Y, units = c("in"),dpi = 300)
  } else if(jj==4){
    ggsave("SOL_KC50_hr15.pdf", plot = zitti2, device="pdf",width = X, height = Y, units = c("in"),dpi = 300)
  }
  
  zitti3<-ggplot(data = world) +
    geom_sf(color = lines, fill = minor)+
    geom_sf(data = states, fill = NA, color=lines,lwd=.25)+
    geom_sf(data = lakes110, fill = lines, color=lines, lwd=.25)+
    coord_sf(crs = zaza, xlim = c(b["xmin"], b["xmax"]), ylim = c(b["ymin"], b["ymax"]))+
    geom_polygon(data=ephestos[[3]], aes(x=long, y=lat, fill=ID, alpha=alpha),show.legend = F,lwd=.25)+
    geom_polygon(data=ephestos[[3]][ephestos[[3]]$ID %in% "kat.3",] , aes(x=long, y=lat),color=nas,fill=NA,lty=c(LDIY),show.legend = F,lwd=TH)+
    geom_polygon(data=ephestos[[3]][ephestos[[3]]$ID %in% "kat.11",], aes(x=long, y=lat),color=nas,fill=NA,lty=c(LDIY),show.legend = F,lwd=TH)+
    geom_polygon(data=ephestos[[3]][ephestos[[3]]$ID %in% "kat.23",], aes(x=long, y=lat),color=nas,fill=NA,lty=c(LDIY),show.legend = F,lwd=TH)+
    geom_polygon(data=ephestos[[3]][ephestos[[3]]$ID %in% "kat.35",], aes(x=long, y=lat),color=nas,fill=NA,lty=c(LDIY),show.legend = F,lwd=TH)+
    geom_polygon(data=ephestos[[3]][ephestos[[3]]$ID %in% "kat.47",], aes(x=long, y=lat),color=nas,fill=NA,lty=c(LDIY),show.legend = F,lwd=TH)+
    theme(panel.grid.major = element_line(color = major, linetype = "dashed", size = .25), 
          panel.background = element_rect(fill = "white"), 
          panel.border = element_rect(color="black",fill = NA),
          axis.text=element_text(size=14),
          axis.title=element_text(size=14))+
    scale_alpha(range = c(0.04, 0.04))+
    scale_fill_manual(values=colfunc(n))+
    labs(x="Longitude",y="Latitude")
  
  if (jj==1){
    ggsave("ALT_BH50_hr19.pdf", plot = zitti3, device="pdf",width = X, height = Y, units = c("in"),dpi = 300)
  } else if(jj==2){
    ggsave("ALT_KC50_hr19.pdf", plot = zitti3, device="pdf",width = X, height = Y, units = c("in"),dpi = 300)
  } else if(jj==3){
    ggsave("SOL_BH50_hr19.pdf", plot = zitti3, device="pdf",width = X, height = Y, units = c("in"),dpi = 300)
  } else if(jj==4){
    ggsave("SOL_KC50_hr19.pdf", plot = zitti3, device="pdf",width = X, height = Y, units = c("in"),dpi = 300)
  }
  
  setwd("~/Dropbox/Madison/Biotron/Alternaria_HYSPLIT/HYSPLIT_sims/for Jacob")
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
#  CURVES CURVES CURVES CURVES CURVES CURVES CURVES CURVES CURVES CURVES CURVES CURVES


ngola<-list()
for (yyy in 1:4){
  #pb <- txtProgressBar(min = 0, max = length(mafuiles), style = 3)
  mmx<-list()
  ii<-yyy
  mmx<-as.data.frame(fread(mafuiles[[ii]], verbose=F))
  colnames(mmx)[c(4,5)]<-c("FT","MX")
  mmx$VOO06<-floor(mmx$FT/360)
  mmx$VOO12<-floor(mmx$FT/720)
  mmx$VOO24<-floor(mmx$FT/1440)
  
  
  mods<-list()
  mods[[1]]<-mmx%>%
    group_by(HOUR) %>%
    group_split()
  bellaciao<-mods[[1]]
  #how many pieces am I splitting each day?
  entregalo<-4
  
  pb <- txtProgressBar(min = 0, max = 3 , style = 3)
  centrir<-list()
  centrirkey<-list()
  for (z in 1:length(bellaciao)){
    centrir[[z]] <- bellaciao[[z]] %>%   
      group_by(VOO06) %>%  #YEAR,MONTH,DAY,VOOINT
      group_split() 
    centrirkey[[z]] <- bellaciao[[z]] %>%
      group_by(VOO24) %>% #YEAR,MONTH,DAY,HOUR,VOOINT
      group_keys()
    setTxtProgressBar(pb, z)
  }

  #CHRONOS CHRONOS CHRONOS CHRONOS CHRONOS CHRONOS CHRONOS CHRONOS CHRONOS CHRONOS 

  if(ii <= 2){
    #set colors ALT
    A<-"#4F97A3"
    B<-"#000080"
    C<-"#0080FF"    
  } else {
    #set colors SOL
    A<-"#7E191B"
    B<-"#ED2929"
    C<-"#FF2400" 
  }
  
  #double
  Aa<-"#4F97A3"
  Ba<-"#000080"
  Ca<-"#0080FF"  
  As<-"#7E191B"
  Bs<-"#ED2929"
  Cs<-"#FF2400" 
  #alpha
  K<-25
  
  #line
  M<-"black"
  
  #line weight
  W<-.05
  
  #mins
  mins<-720
  
  tepido<-theme(panel.border = element_rect(linetype = 1, fill = NA),
                panel.grid.major = element_line(size=0.05,colour = "grey50"),
                panel.grid.minor = element_line(size=0.05,colour = "grey50"),
                panel.background = element_rect(fill = NA),
                axis.text.x = element_text(angle=0),
                axis.text.y = element_text( angle=-0))
  teruego<-theme(panel.border = element_rect(linetype = 1, fill = NA),
                 panel.grid.major = element_line(size=0.05,colour = "grey50"),
                 panel.grid.minor = element_line(size=0.05,colour = "grey50"),
                 panel.background = element_rect(fill = NA),
                 axis.text.x = element_text(angle=0,size=3),
                 axis.text.y = element_text( angle=-0,size=3),
                 axis.title = element_text(size=4))
  
  fancy_scientific <- function(l) {
    # turn in to character string in scientific notation
    l <- formatC(l, format = "e", digits = 2)
    # quote the part before the exponent to keep all the digits
    l <- gsub("^(.*)e", "'\\1'e", l)
    # turn the 'e+' into plotmath format
    l <- gsub("e", "%*%10^", l)
    # return this as an expression
    parse(text=l)
  }
  
  #################################
  kumura<-list()
  for (c in 1:length(centrir))
    kumura[[c]] <-pbapply::pblapply(centrir[[c]], function(x) 
      x %>%
        dplyr::group_by(VOO06) %>%
        dplyr::summarise(
          n=n(),
          mpH=mean(DIST/1000),
          sdH=sd(DIST/1000)
        )  %>%
        dplyr::mutate(seH=sdH/sqrt(n)))
  
  sodastars<-pbapply::pblapply(kumura, function(x) bind_rows(x))
  
  habibi<-list()
  for (c in 1:length(centrir))
    habibi[[c]] <-pbapply::pblapply(centrir[[c]], function(x) 
      x %>%
        dplyr::group_by(FT) %>%
        dplyr::summarise(
          n=n(),
          mpH=mean(DIST/1000),
          sdH=sd(DIST/1000),
          mpH.sub=mean(DIST),
          sdH.sub=sd(DIST)
        )  %>%
        dplyr::mutate(seH=sdH/sqrt(n)) %>%
        dplyr::mutate(seH.sub=sdH.sub/sqrt(n)))
  
  oumri<-pbapply::pblapply(habibi, function(x) bind_rows(x))
  
  
  kmvt<-ggplot(data=NULL)+
    geom_ribbon(data=sodastars[[1]],aes(x=VOO06,ymin = mpH - 5*seH, ymax = mpH + 5*seH), fill = A,alpha=1/K) +
    geom_ribbon(data=sodastars[[1]],aes(x=VOO06,ymin = mpH - 4*seH, ymax = mpH + 4*seH), fill = A,alpha=2/K) +
    geom_ribbon(data=sodastars[[1]],aes(x=VOO06,ymin = mpH - 3*seH, ymax = mpH + 3*seH), fill = A,alpha=3/K) +
    geom_ribbon(data=sodastars[[1]],aes(x=VOO06,ymin = mpH - 2*seH, ymax = mpH + 2*seH), fill = A,alpha=4/K) +
    geom_ribbon(data=sodastars[[1]],aes(x=VOO06,ymin = mpH - seH, ymax = mpH + seH), fill = A,alpha=5/K) +
    geom_ribbon(data=sodastars[[2]],aes(x=VOO06,ymin = mpH - 5*seH, ymax = mpH + 5*seH), fill = B,alpha=1/K) +
    geom_ribbon(data=sodastars[[2]],aes(x=VOO06,ymin = mpH - 4*seH, ymax = mpH + 4*seH), fill = B,alpha=2/K) +
    geom_ribbon(data=sodastars[[2]],aes(x=VOO06,ymin = mpH - 3*seH, ymax = mpH + 3*seH), fill = B,alpha=3/K) +
    geom_ribbon(data=sodastars[[2]],aes(x=VOO06,ymin = mpH - 2*seH, ymax = mpH + 2*seH), fill = B,alpha=4/K) +
    geom_ribbon(data=sodastars[[2]],aes(x=VOO06,ymin = mpH - seH, ymax = mpH + seH), fill = B,alpha=5/K) +
    geom_ribbon(data=sodastars[[3]],aes(x=VOO06,ymin = mpH - 5*seH, ymax = mpH + 5*seH), fill = C,alpha=1/K) +
    geom_ribbon(data=sodastars[[3]],aes(x=VOO06,ymin = mpH - 4*seH, ymax = mpH + 4*seH), fill = C,alpha=2/K) +
    geom_ribbon(data=sodastars[[3]],aes(x=VOO06,ymin = mpH - 3*seH, ymax = mpH + 3*seH), fill = C,alpha=3/K) +
    geom_ribbon(data=sodastars[[3]],aes(x=VOO06,ymin = mpH - 2*seH, ymax = mpH + 2*seH), fill = C,alpha=4/K) +
    geom_ribbon(data=sodastars[[3]],aes(x=VOO06,ymin = mpH - seH, ymax = mpH + seH), fill = C,alpha=5/K) +
    geom_line(data=sodastars[[1]],aes(x=VOO06, y=mpH), linetype=1, lwd=W, color=M)+
    geom_line(data=sodastars[[2]],aes(x=VOO06, y=mpH), linetype=1, lwd=W, color=M)+
    geom_line(data=sodastars[[3]],aes(x=VOO06, y=mpH), linetype=1, lwd=W, color=M)+
    tepido+
    scale_x_continuous(name="Time in Flight (hrs)",
                       breaks=c(0,4,12,24,36,48),
                       labels=c(6,24,72,144,216,288),
                       limits = c(0,48))+
    scale_y_continuous(name="Distance Traveled (km)",labels=fancy_scientific)+
    coord_cartesian(ylim = c(0, 3500)) 
  
  kmvt.10<-ggplot(data=NULL)+
    geom_ribbon(data=oumri[[1]],aes(x=FT,ymin = mpH - 5*seH, ymax = mpH + 5*seH), fill = A,alpha=1/K) +
    geom_ribbon(data=oumri[[1]],aes(x=FT,ymin = mpH - 4*seH, ymax = mpH + 4*seH), fill = A,alpha=2/K) +
    geom_ribbon(data=oumri[[1]],aes(x=FT,ymin = mpH - 3*seH, ymax = mpH + 3*seH), fill = A,alpha=3/K) +
    geom_ribbon(data=oumri[[1]],aes(x=FT,ymin = mpH - 2*seH, ymax = mpH + 2*seH), fill = A,alpha=4/K) +
    geom_ribbon(data=oumri[[1]],aes(x=FT,ymin = mpH - seH, ymax = mpH + seH), fill = A,alpha=5/K) +
    geom_ribbon(data=oumri[[2]],aes(x=FT,ymin = mpH - 5*seH, ymax = mpH + 5*seH), fill = B,alpha=1/K) +
    geom_ribbon(data=oumri[[2]],aes(x=FT,ymin = mpH - 4*seH, ymax = mpH + 4*seH), fill = B,alpha=2/K) +
    geom_ribbon(data=oumri[[2]],aes(x=FT,ymin = mpH - 3*seH, ymax = mpH + 3*seH), fill = B,alpha=3/K) +
    geom_ribbon(data=oumri[[2]],aes(x=FT,ymin = mpH - 2*seH, ymax = mpH + 2*seH), fill = B,alpha=4/K) +
    geom_ribbon(data=oumri[[2]],aes(x=FT,ymin = mpH - seH, ymax = mpH + seH), fill = B,alpha=5/K) +
    geom_ribbon(data=oumri[[3]],aes(x=FT,ymin = mpH - 5*seH, ymax = mpH + 5*seH), fill = C,alpha=1/K) +
    geom_ribbon(data=oumri[[3]],aes(x=FT,ymin = mpH - 4*seH, ymax = mpH + 4*seH), fill = C,alpha=2/K) +
    geom_ribbon(data=oumri[[3]],aes(x=FT,ymin = mpH - 3*seH, ymax = mpH + 3*seH), fill = C,alpha=3/K) +
    geom_ribbon(data=oumri[[3]],aes(x=FT,ymin = mpH - 2*seH, ymax = mpH + 2*seH), fill = C,alpha=4/K) +
    geom_ribbon(data=oumri[[3]],aes(x=FT,ymin = mpH - seH, ymax = mpH + seH), fill = C,alpha=5/K) +
    geom_line(data=oumri[[1]],aes(x=FT, y=mpH), linetype=1, lwd=W, color=M)+
    geom_line(data=oumri[[2]],aes(x=FT, y=mpH), linetype=1, lwd=W, color=M)+
    geom_line(data=oumri[[3]],aes(x=FT, y=mpH), linetype=1, lwd=W, color=M)+
    tepido+
    scale_x_continuous(name="Time in Flight (hrs)",
                       breaks=c(0,1440,4320,8640,12960,17280),
                       labels=c(6,24,72,144,216,288),
                       limits = c(0,17280))+
    scale_y_continuous(name="Distance Traveled (km)",labels=fancy_scientific)+
    coord_cartesian(ylim = c(0, 3500)) 
  
  kmvt.sub<-ggplot(data=NULL)+
    geom_ribbon(data=oumri[[1]][oumri[[1]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 5*seH.sub, ymax = mpH.sub + 5*seH.sub), fill = A,alpha=1/K) +
    geom_ribbon(data=oumri[[1]][oumri[[1]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 4*seH.sub, ymax = mpH.sub + 4*seH.sub), fill = A,alpha=2/K) +
    geom_ribbon(data=oumri[[1]][oumri[[1]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 3*seH.sub, ymax = mpH.sub + 3*seH.sub), fill = A,alpha=3/K) +
    geom_ribbon(data=oumri[[1]][oumri[[1]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 2*seH.sub, ymax = mpH.sub + 2*seH.sub), fill = A,alpha=4/K) +
    geom_ribbon(data=oumri[[1]][oumri[[1]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - seH.sub, ymax = mpH.sub + seH.sub), fill = A,alpha=5/K) +
    geom_line(data=oumri[[1]][oumri[[1]]$FT <=mins,],aes(x=FT, y=mpH.sub), linetype=1, lwd=W, color=M)+
    geom_ribbon(data=oumri[[2]][oumri[[2]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 5*seH.sub, ymax = mpH.sub + 5*seH.sub), fill = B,alpha=1/K) +
    geom_ribbon(data=oumri[[2]][oumri[[2]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 4*seH.sub, ymax = mpH.sub + 4*seH.sub), fill = B,alpha=2/K) +
    geom_ribbon(data=oumri[[2]][oumri[[2]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 3*seH.sub, ymax = mpH.sub + 3*seH.sub), fill = B,alpha=3/K) +
    geom_ribbon(data=oumri[[2]][oumri[[2]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 2*seH.sub, ymax = mpH.sub + 2*seH.sub), fill = B,alpha=4/K) +
    geom_ribbon(data=oumri[[2]][oumri[[2]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - seH.sub, ymax = mpH.sub + seH.sub), fill = B,alpha=5/K) +
    geom_ribbon(data=oumri[[3]][oumri[[3]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 5*seH.sub, ymax = mpH.sub + 5*seH.sub), fill = C,alpha=1/K) +
    geom_ribbon(data=oumri[[3]][oumri[[3]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 4*seH.sub, ymax = mpH.sub + 4*seH.sub), fill = C,alpha=2/K) +
    geom_ribbon(data=oumri[[3]][oumri[[3]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 3*seH.sub, ymax = mpH.sub + 3*seH.sub), fill = C,alpha=3/K) +
    geom_ribbon(data=oumri[[3]][oumri[[3]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 2*seH.sub, ymax = mpH.sub + 2*seH.sub), fill = C,alpha=4/K) +
    geom_ribbon(data=oumri[[3]][oumri[[3]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - seH.sub, ymax = mpH.sub + seH.sub), fill = C,alpha=5/K) +
    geom_line(data=oumri[[1]][oumri[[1]]$FT <=mins,],aes(x=FT, y=mpH.sub), linetype=1, lwd=W, color=M)+
    geom_line(data=oumri[[2]][oumri[[2]]$FT <=mins,],aes(x=FT, y=mpH.sub), linetype=1, lwd=W, color=M)+
    geom_line(data=oumri[[3]][oumri[[3]]$FT <=mins,],aes(x=FT, y=mpH.sub), linetype=1, lwd=W, color=M)+
    teruego+
    scale_x_continuous(name="Time in Flight (mins)",
                       breaks=seq(0,720,120),
                       limits=c(0,720))+
    scale_y_continuous(name="Distance Traveled (m)",labels=fancy_scientific)+
    coord_cartesian(ylim = c(0, 3.5e5)) 
  
  
  ############################
  diamonds<-list()
  for (c in 1:length(centrir))
    diamonds[[c]] <-pbapply::pblapply(centrir[[c]], function(x) 
      x %>%
        dplyr::group_by(VOO06,YEAR) %>%
        dplyr::summarise(
          n=n()) )
  
  daemons<-pbapply::pblapply(diamonds, function(x) bind_rows(x))
  
  for (u in 1:length(daemons)){
    colnames(daemons[[u]])[3]<-"spores"
  }
  
  fatass<-pblapply(daemons, function(x) 
    x %>%
      dplyr::group_by(YEAR) %>%
      mutate(cumdump = cumsum(spores)))
  
  katsura<-pblapply(fatass, function(x) 
    x %>%
      dplyr::group_by(VOO06) %>%
      dplyr::summarize(
        n=n(),
        meancum=mean(cumdump/20e5),
        sdcum=sd(cumdump/20e5)
      ) %>%
      mutate(secum=sdcum/sqrt(n)))
  
  ngola[[yyy]]<-katsura
  
  # SMALL SMALL SMALL SMALL
  emeralds<-list()
  for (c in 1:length(centrir))
    emeralds[[c]] <-pbapply::pblapply(centrir[[c]], function(x) 
      x %>%
        dplyr::group_by(FT,YEAR) %>%
        dplyr::summarise(
          n=n()) )
  
  kraken<-pbapply::pblapply(emeralds, function(x) bind_rows(x))
  
  for (u in 1:length(kraken)){
    colnames(kraken[[u]])[3]<-"spores"
  }
  
  stab<-pblapply(kraken, function(x) 
    x %>%
      dplyr::group_by(YEAR) %>%
      mutate(cumdump = cumsum(spores)))
  
  fuckthepolice<-pblapply(stab, function(x) 
    x %>%
      dplyr::group_by(FT) %>%
      dplyr::summarize(
        n=n(),
        meancum=mean(cumdump/20e5),
        sdcum=sd(cumdump/20e5)
      ) %>%
      mutate(secum=sdcum/sqrt(n)))
  
  dep<-ggplot(data=NULL)+
    geom_ribbon(data=katsura[[1]],aes(x=VOO06,ymin = meancum - 5*secum, ymax = meancum + 5*secum), fill = A,alpha=1/K) +
    geom_ribbon(data=katsura[[1]],aes(x=VOO06,ymin = meancum - 4*secum, ymax = meancum + 4*secum), fill = A,alpha=2/K) +
    geom_ribbon(data=katsura[[1]],aes(x=VOO06,ymin = meancum - 3*secum, ymax = meancum + 3*secum), fill = A,alpha=3/K) +
    geom_ribbon(data=katsura[[1]],aes(x=VOO06,ymin = meancum - 2*secum, ymax = meancum + 2*secum), fill = A,alpha=4/K) +
    geom_ribbon(data=katsura[[1]],aes(x=VOO06,ymin = meancum - secum,   ymax = meancum + secum), fill = A,alpha=5/K) +
    geom_ribbon(data=katsura[[2]],aes(x=VOO06,ymin = meancum - 5*secum, ymax = meancum + 5*secum), fill = B,alpha=1/K) +
    geom_ribbon(data=katsura[[2]],aes(x=VOO06,ymin = meancum - 4*secum, ymax = meancum + 4*secum), fill = B,alpha=2/K) +
    geom_ribbon(data=katsura[[2]],aes(x=VOO06,ymin = meancum - 3*secum, ymax = meancum + 3*secum), fill = B,alpha=3/K) +
    geom_ribbon(data=katsura[[2]],aes(x=VOO06,ymin = meancum - 2*secum, ymax = meancum + 2*secum), fill = B,alpha=4/K) +
    geom_ribbon(data=katsura[[2]],aes(x=VOO06,ymin = meancum - secum,   ymax = meancum + secum), fill = B,alpha=5/K) +
    geom_ribbon(data=katsura[[3]],aes(x=VOO06,ymin = meancum - 5*secum, ymax = meancum + 5*secum), fill = C,alpha=1/K) +
    geom_ribbon(data=katsura[[3]],aes(x=VOO06,ymin = meancum - 4*secum, ymax = meancum + 4*secum), fill = C,alpha=2/K) +
    geom_ribbon(data=katsura[[3]],aes(x=VOO06,ymin = meancum - 3*secum, ymax = meancum + 3*secum), fill = C,alpha=3/K) +
    geom_ribbon(data=katsura[[3]],aes(x=VOO06,ymin = meancum - 2*secum, ymax = meancum + 2*secum), fill = C,alpha=4/K) +
    geom_ribbon(data=katsura[[3]],aes(x=VOO06,ymin = meancum - secum,   ymax = meancum + secum), fill = C,alpha=5/K) +
    geom_line(data=katsura[[1]],aes(x=VOO06,y=meancum),linetype=1, lwd=W, color=M)+
    geom_line(data=katsura[[2]],aes(x=VOO06,y=meancum),linetype=1, lwd=W, color=M)+
    geom_line(data=katsura[[3]],aes(x=VOO06,y=meancum),linetype=1, lwd=W, color=M)+
    tepido+
    scale_x_continuous(name="Time in Flight (hrs)",
                       breaks=c(0,4,12,24,36,48),
                       labels=c(6,24,72,144,216,288),
                       limits=c(0,48))+
    scale_y_reverse(name="Spores Remaining Airborne",labels=fancy_scientific)+
    coord_cartesian(ylim = c(0, 1), xlim=c(0,48)) 
  
  dep.10<-ggplot(data=NULL)+
    geom_ribbon(data=fuckthepolice[[1]],aes(x=FT,ymin = meancum - 5*secum, ymax = meancum + 5*secum), fill = A,alpha=1/K) +
    geom_ribbon(data=fuckthepolice[[1]],aes(x=FT,ymin = meancum - 4*secum, ymax = meancum + 4*secum), fill = A,alpha=2/K) +
    geom_ribbon(data=fuckthepolice[[1]],aes(x=FT,ymin = meancum - 3*secum, ymax = meancum + 3*secum), fill = A,alpha=3/K) +
    geom_ribbon(data=fuckthepolice[[1]],aes(x=FT,ymin = meancum - 2*secum, ymax = meancum + 2*secum), fill = A,alpha=4/K) +
    geom_ribbon(data=fuckthepolice[[1]],aes(x=FT,ymin = meancum - secum,   ymax = meancum + secum), fill=    A,alpha=5/K) +
    geom_ribbon(data=fuckthepolice[[2]],aes(x=FT,ymin = meancum - 5*secum, ymax = meancum + 5*secum), fill = B,alpha=1/K) +
    geom_ribbon(data=fuckthepolice[[2]],aes(x=FT,ymin = meancum - 4*secum, ymax = meancum + 4*secum), fill = B,alpha=2/K) +
    geom_ribbon(data=fuckthepolice[[2]],aes(x=FT,ymin = meancum - 3*secum, ymax = meancum + 3*secum), fill = B,alpha=3/K) +
    geom_ribbon(data=fuckthepolice[[2]],aes(x=FT,ymin = meancum - 2*secum, ymax = meancum + 2*secum), fill = B,alpha=4/K) +
    geom_ribbon(data=fuckthepolice[[2]],aes(x=FT,ymin = meancum - secum,   ymax = meancum + secum), fill =   B,alpha=5/K) +
    geom_ribbon(data=fuckthepolice[[3]],aes(x=FT,ymin = meancum - 5*secum, ymax = meancum + 5*secum), fill = C,alpha=1/K) +
    geom_ribbon(data=fuckthepolice[[3]],aes(x=FT,ymin = meancum - 4*secum, ymax = meancum + 4*secum), fill = C,alpha=2/K) +
    geom_ribbon(data=fuckthepolice[[3]],aes(x=FT,ymin = meancum - 3*secum, ymax = meancum + 3*secum), fill = C,alpha=3/K) +
    geom_ribbon(data=fuckthepolice[[3]],aes(x=FT,ymin = meancum - 2*secum, ymax = meancum + 2*secum), fill = C,alpha=4/K) +
    geom_ribbon(data=fuckthepolice[[3]],aes(x=FT,ymin = meancum - secum,   ymax = meancum + secum), fill = C,alpha=  5/K) +
    geom_line(data=fuckthepolice[[1]],aes(x=FT,y=meancum),linetype=1, lwd=W, color=M)+
    geom_line(data=fuckthepolice[[2]],aes(x=FT,y=meancum),linetype=1, lwd=W, color=M)+
    geom_line(data=fuckthepolice[[3]],aes(x=FT,y=meancum),linetype=1, lwd=W, color=M)+
    tepido+
    scale_x_continuous(name="Time in Flight (hrs)",
                       breaks=c(0,1440,4320,8640,12960,17280),
                       labels=c(6,24,72,144,216,288),
                       limits=c(0,17280))+
    scale_y_reverse(name="Spores Remaining Airborne",labels=fancy_scientific)+
    coord_cartesian(ylim = c(0, 1)) 
  
  
  dep.sub<-ggplot(data=NULL)+
    geom_ribbon(data=fuckthepolice[[1]][fuckthepolice[[1]]$FT <= mins,],aes(x=FT,ymin = meancum - 5*secum, ymax = meancum + 5*secum), fill = A,alpha=1/K) +
    geom_ribbon(data=fuckthepolice[[1]][fuckthepolice[[1]]$FT <= mins,],aes(x=FT,ymin = meancum - 4*secum, ymax = meancum + 4*secum), fill = A,alpha=2/K) +
    geom_ribbon(data=fuckthepolice[[1]][fuckthepolice[[1]]$FT <= mins,],aes(x=FT,ymin = meancum - 3*secum, ymax = meancum + 3*secum), fill = A,alpha=3/K) +
    geom_ribbon(data=fuckthepolice[[1]][fuckthepolice[[1]]$FT <= mins,],aes(x=FT,ymin = meancum - 2*secum, ymax = meancum + 2*secum), fill = A,alpha=4/K) +
    geom_ribbon(data=fuckthepolice[[1]][fuckthepolice[[1]]$FT <= mins,],aes(x=FT,ymin = meancum - secum,   ymax = meancum + secum), fill = A,alpha=5/K) +
    geom_ribbon(data=fuckthepolice[[2]][fuckthepolice[[2]]$FT <= mins,],aes(x=FT,ymin = meancum - 5*secum, ymax = meancum + 5*secum), fill = B,alpha=1/K) +
    geom_ribbon(data=fuckthepolice[[2]][fuckthepolice[[2]]$FT <= mins,],aes(x=FT,ymin = meancum - 4*secum, ymax = meancum + 4*secum), fill = B,alpha=2/K) +
    geom_ribbon(data=fuckthepolice[[2]][fuckthepolice[[2]]$FT <= mins,],aes(x=FT,ymin = meancum - 3*secum, ymax = meancum + 3*secum), fill = B,alpha=3/K) +
    geom_ribbon(data=fuckthepolice[[2]][fuckthepolice[[2]]$FT <= mins,],aes(x=FT,ymin = meancum - 2*secum, ymax = meancum + 2*secum), fill = B,alpha=4/K) +
    geom_ribbon(data=fuckthepolice[[2]][fuckthepolice[[2]]$FT <= mins,],aes(x=FT,ymin = meancum - secum,   ymax = meancum + secum), fill = C,alpha=5/K) + 
    geom_ribbon(data=fuckthepolice[[3]][fuckthepolice[[3]]$FT <= mins,],aes(x=FT,ymin = meancum - 5*secum, ymax = meancum + 5*secum), fill = C,alpha=1/K) +
    geom_ribbon(data=fuckthepolice[[3]][fuckthepolice[[3]]$FT <= mins,],aes(x=FT,ymin = meancum - 4*secum, ymax = meancum + 4*secum), fill = C,alpha=2/K) +
    geom_ribbon(data=fuckthepolice[[3]][fuckthepolice[[3]]$FT <= mins,],aes(x=FT,ymin = meancum - 3*secum, ymax = meancum + 3*secum), fill = C,alpha=3/K) +
    geom_ribbon(data=fuckthepolice[[3]][fuckthepolice[[3]]$FT <= mins,],aes(x=FT,ymin = meancum - 2*secum, ymax = meancum + 2*secum), fill = C,alpha=4/K) +
    geom_ribbon(data=fuckthepolice[[3]][fuckthepolice[[3]]$FT <= mins,],aes(x=FT,ymin = meancum - secum,   ymax = meancum + secum), fill = C,alpha=5/K) +
    geom_line(data=fuckthepolice[[1]][fuckthepolice[[1]]$FT <= mins,],aes(x=FT,y=meancum),linetype=1, lwd=W, color=M)+
    geom_line(data=fuckthepolice[[2]][fuckthepolice[[2]]$FT <= mins,],aes(x=FT,y=meancum),linetype=1, lwd=W, color=M)+
    geom_line(data=fuckthepolice[[3]][fuckthepolice[[3]]$FT <= mins,],aes(x=FT,y=meancum),linetype=1, lwd=W, color=M)+
    teruego+
    scale_x_continuous(name="Time in Flight (mins)",
                       breaks=seq(0,720,120),
                       limits = c(0,720))+
    scale_y_reverse(name="Spores Remaining Airborne",labels=fancy_scientific)+
    coord_cartesian(ylim = c(0, 1)) 
  
  ##################################
  rooski<-list()
  for (c in 1:length(centrir))
    rooski[[c]] <-pbapply::pblapply(centrir[[c]], function(x) 
      x %>%
        dplyr::group_by(VOO06) %>%
        dplyr::summarise(
          n=n(),
          mean.mx=mean(MX/1000),
          sd.mx=sd(MX/1000)
        )  %>%
        dplyr::mutate(se.mx=sd.mx/sqrt(n)))
  
  stalingrad<-pbapply::pblapply(rooski, function(x) bind_rows(x))
  
  ##### FINE SCALE
  
  USSR<-list()
  for (c in 1:length(centrir))
    USSR[[c]] <-pbapply::pblapply(centrir[[c]], function(x) 
      x %>%
        dplyr::group_by(FT) %>%
        dplyr::summarise(
          n=n(),
          mean.mx=mean(MX),
          sd.mx=sd(MX),
          mean.mx.sub=mean(MX/1000),
          sd.mx.sub=sd(MX/1000)
        )  %>%
        dplyr::mutate(se.mx=sd.mx/sqrt(n)) %>%
        dplyr::mutate(se.mx.sub=sd.mx.sub/sqrt(n)))
  
  leningrad<-pbapply::pblapply(USSR, function(x) bind_rows(x))
  ###
  
  hvt<-ggplot(data=NULL)+
    geom_ribbon(data=stalingrad[[1]],aes(x=VOO06,ymin = mean.mx - 5*se.mx, ymax = mean.mx + 5*se.mx), fill = A,alpha=1/K) +
    geom_ribbon(data=stalingrad[[1]],aes(x=VOO06,ymin = mean.mx - 4*se.mx, ymax = mean.mx + 4*se.mx), fill = A,alpha=2/K) +
    geom_ribbon(data=stalingrad[[1]],aes(x=VOO06,ymin = mean.mx - 3*se.mx, ymax = mean.mx + 3*se.mx), fill = A,alpha=3/K) +
    geom_ribbon(data=stalingrad[[1]],aes(x=VOO06,ymin = mean.mx - 2*se.mx, ymax = mean.mx + 2*se.mx), fill = A,alpha=4/K) +
    geom_ribbon(data=stalingrad[[1]],aes(x=VOO06,ymin = mean.mx - se.mx, ymax   = mean.mx + se.mx),   fill = A,alpha=5/K) +
    geom_ribbon(data=stalingrad[[2]],aes(x=VOO06,ymin = mean.mx - 5*se.mx, ymax = mean.mx + 5*se.mx), fill = B,alpha=1/K) +
    geom_ribbon(data=stalingrad[[2]],aes(x=VOO06,ymin = mean.mx - 4*se.mx, ymax = mean.mx + 4*se.mx), fill = B,alpha=2/K) +
    geom_ribbon(data=stalingrad[[2]],aes(x=VOO06,ymin = mean.mx - 3*se.mx, ymax = mean.mx + 3*se.mx), fill = B,alpha=3/K) +
    geom_ribbon(data=stalingrad[[2]],aes(x=VOO06,ymin = mean.mx - 2*se.mx, ymax = mean.mx + 2*se.mx), fill = B,alpha=4/K) +
    geom_ribbon(data=stalingrad[[2]],aes(x=VOO06,ymin = mean.mx - se.mx, ymax   = mean.mx + se.mx),   fill = B,alpha=5/K) +
    geom_ribbon(data=stalingrad[[3]],aes(x=VOO06,ymin = mean.mx - 5*se.mx, ymax = mean.mx + 5*se.mx), fill = C,alpha=1/K) +
    geom_ribbon(data=stalingrad[[3]],aes(x=VOO06,ymin = mean.mx - 4*se.mx, ymax = mean.mx + 4*se.mx), fill = C,alpha=2/K) +
    geom_ribbon(data=stalingrad[[3]],aes(x=VOO06,ymin = mean.mx - 3*se.mx, ymax = mean.mx + 3*se.mx), fill = C,alpha=3/K) +
    geom_ribbon(data=stalingrad[[3]],aes(x=VOO06,ymin = mean.mx - 2*se.mx, ymax = mean.mx + 2*se.mx), fill = C,alpha=4/K) +
    geom_ribbon(data=stalingrad[[3]],aes(x=VOO06,ymin = mean.mx - se.mx, ymax   = mean.mx + se.mx),   fill = C,alpha=5/K) +
    geom_line(data=stalingrad[[1]], aes(x=VOO06, y=mean.mx),linetype=1, lwd=W, color=M)+
    geom_line(data=stalingrad[[2]], aes(x=VOO06, y=mean.mx),linetype=1, lwd=W, color=M)+
    geom_line(data=stalingrad[[3]], aes(x=VOO06, y=mean.mx),linetype=1, lwd=W, color=M)+
    tepido+
    scale_x_continuous(name="Time in Flight (hrs)",
                       breaks=c(0,4,12,24,36,48),
                       labels=c(0,24,72,144,216,288),
                       limits = c(0,48))+
    scale_y_continuous(name="Altitude (km)",labels=fancy_scientific)+
    coord_cartesian(ylim = c(0, 7),xlim=c(0,48)) 
  
  
  hvt.10<-ggplot(data=NULL)+
    geom_ribbon(data=leningrad[[1]],aes(x=FT,ymin = mean.mx - 5*se.mx, ymax = mean.mx + 5*se.mx), fill = A,alpha=1/K) +
    geom_ribbon(data=leningrad[[1]],aes(x=FT,ymin = mean.mx - 4*se.mx, ymax = mean.mx + 4*se.mx), fill = A,alpha=2/K) +
    geom_ribbon(data=leningrad[[1]],aes(x=FT,ymin = mean.mx - 3*se.mx, ymax = mean.mx + 3*se.mx), fill = A,alpha=3/K) +
    geom_ribbon(data=leningrad[[1]],aes(x=FT,ymin = mean.mx - 2*se.mx, ymax = mean.mx + 2*se.mx), fill = A,alpha=4/K) +
    geom_ribbon(data=leningrad[[1]],aes(x=FT,ymin = mean.mx - se.mx, ymax   = mean.mx + se.mx),   fill = A,alpha=5/K) +
    geom_ribbon(data=leningrad[[2]],aes(x=FT,ymin = mean.mx - 5*se.mx, ymax = mean.mx + 5*se.mx), fill = B,alpha=1/K) +
    geom_ribbon(data=leningrad[[2]],aes(x=FT,ymin = mean.mx - 4*se.mx, ymax = mean.mx + 4*se.mx), fill = B,alpha=2/K) +
    geom_ribbon(data=leningrad[[2]],aes(x=FT,ymin = mean.mx - 3*se.mx, ymax = mean.mx + 3*se.mx), fill = B,alpha=3/K) +
    geom_ribbon(data=leningrad[[2]],aes(x=FT,ymin = mean.mx - 2*se.mx, ymax = mean.mx + 2*se.mx), fill = B,alpha=4/K) +
    geom_ribbon(data=leningrad[[2]],aes(x=FT,ymin = mean.mx - se.mx, ymax   = mean.mx + se.mx),   fill = B,alpha=5/K) +
    geom_ribbon(data=leningrad[[3]],aes(x=FT,ymin = mean.mx - 5*se.mx, ymax = mean.mx + 5*se.mx), fill = C,alpha=1/K) +
    geom_ribbon(data=leningrad[[3]],aes(x=FT,ymin = mean.mx - 4*se.mx, ymax = mean.mx + 4*se.mx), fill = C,alpha=2/K) +
    geom_ribbon(data=leningrad[[3]],aes(x=FT,ymin = mean.mx - 3*se.mx, ymax = mean.mx + 3*se.mx), fill = C,alpha=3/K) +
    geom_ribbon(data=leningrad[[3]],aes(x=FT,ymin = mean.mx - 2*se.mx, ymax = mean.mx + 2*se.mx), fill = C,alpha=4/K) +
    geom_ribbon(data=leningrad[[3]],aes(x=FT,ymin = mean.mx - se.mx, ymax   = mean.mx + se.mx),   fill = C,alpha=5/K) +
    geom_line(data=leningrad[[1]], aes(x=FT, y=mean.mx),linetype=1, lwd=W, color=M)+
    geom_line(data=leningrad[[2]], aes(x=FT, y=mean.mx),linetype=1, lwd=W, color=M)+
    geom_line(data=leningrad[[3]], aes(x=FT, y=mean.mx),linetype=1, lwd=W, color=M)+
    tepido+
    scale_x_continuous(name="Time in Flight (hrs)",
                       breaks=c(0,1440,4320,8640,12960,17280),
                       labels=c(6,24,72,144,216,288),
                       limits=c(0,17280))+
    scale_y_continuous(name="Altitude (km)",labels=fancy_scientific)+
    coord_cartesian(ylim = c(0, 8e3), xlim=c(0,max(leningrad[[2]]$FT)))
  
  
  hvt.sub<-ggplot(data=NULL)+
    geom_ribbon(data=leningrad[[1]][leningrad[[1]]$FT <= mins,],aes(x=FT,ymin = mean.mx.sub - 5*se.mx.sub, ymax = mean.mx.sub + 5*se.mx.sub), fill = A,alpha=1/K) +
    geom_ribbon(data=leningrad[[1]][leningrad[[1]]$FT <= mins,],aes(x=FT,ymin = mean.mx.sub - 4*se.mx.sub, ymax = mean.mx.sub + 4*se.mx.sub), fill = A,alpha=2/K) +
    geom_ribbon(data=leningrad[[1]][leningrad[[1]]$FT <= mins,],aes(x=FT,ymin = mean.mx.sub - 3*se.mx.sub, ymax = mean.mx.sub + 3*se.mx.sub), fill = A,alpha=3/K) +
    geom_ribbon(data=leningrad[[1]][leningrad[[1]]$FT <= mins,],aes(x=FT,ymin = mean.mx.sub - 2*se.mx.sub, ymax = mean.mx.sub + 2*se.mx.sub), fill = A,alpha=4/K) +
    geom_ribbon(data=leningrad[[1]][leningrad[[1]]$FT <= mins,],aes(x=FT,ymin = mean.mx.sub - se.mx.sub, ymax   = mean.mx.sub + se.mx.sub),   fill = A,alpha=4/K) +
    geom_ribbon(data=leningrad[[2]][leningrad[[2]]$FT <= mins,],aes(x=FT,ymin = mean.mx.sub - 5*se.mx.sub, ymax = mean.mx.sub + 5*se.mx.sub), fill = A,alpha=1/K) +
    geom_ribbon(data=leningrad[[2]][leningrad[[2]]$FT <= mins,],aes(x=FT,ymin = mean.mx.sub - 4*se.mx.sub, ymax = mean.mx.sub + 4*se.mx.sub), fill = A,alpha=2/K) +
    geom_ribbon(data=leningrad[[2]][leningrad[[2]]$FT <= mins,],aes(x=FT,ymin = mean.mx.sub - 3*se.mx.sub, ymax = mean.mx.sub + 3*se.mx.sub), fill = A,alpha=3/K) +
    geom_ribbon(data=leningrad[[2]][leningrad[[2]]$FT <= mins,],aes(x=FT,ymin = mean.mx.sub - 2*se.mx.sub, ymax = mean.mx.sub + 2*se.mx.sub), fill = A,alpha=4/K) +
    geom_ribbon(data=leningrad[[2]][leningrad[[2]]$FT <= mins,],aes(x=FT,ymin = mean.mx.sub - se.mx.sub, ymax   = mean.mx.sub + se.mx.sub),   fill = A,alpha=4/K) +
    geom_ribbon(data=leningrad[[3]][leningrad[[3]]$FT <= mins,],aes(x=FT,ymin = mean.mx.sub - 5*se.mx.sub, ymax = mean.mx.sub + 5*se.mx.sub), fill = A,alpha=1/K) +
    geom_ribbon(data=leningrad[[3]][leningrad[[3]]$FT <= mins,],aes(x=FT,ymin = mean.mx.sub - 4*se.mx.sub, ymax = mean.mx.sub + 4*se.mx.sub), fill = A,alpha=2/K) +
    geom_ribbon(data=leningrad[[3]][leningrad[[3]]$FT <= mins,],aes(x=FT,ymin = mean.mx.sub - 3*se.mx.sub, ymax = mean.mx.sub + 3*se.mx.sub), fill = A,alpha=3/K) +
    geom_ribbon(data=leningrad[[3]][leningrad[[3]]$FT <= mins,],aes(x=FT,ymin = mean.mx.sub - 2*se.mx.sub, ymax = mean.mx.sub + 2*se.mx.sub), fill = A,alpha=4/K) +
    geom_ribbon(data=leningrad[[3]][leningrad[[3]]$FT <= mins,],aes(x=FT,ymin = mean.mx.sub - se.mx.sub, ymax   = mean.mx.sub + se.mx.sub),   fill = A,alpha=4/K) +
    geom_line(data=leningrad[[1]][leningrad[[1]]$FT <= mins,], aes(x=FT, y=mean.mx.sub),linetype=1, lwd=W, color=M)+
    geom_line(data=leningrad[[2]][leningrad[[2]]$FT <= mins,], aes(x=FT, y=mean.mx.sub),linetype=1, lwd=W, color=M)+
    geom_line(data=leningrad[[3]][leningrad[[3]]$FT <= mins,], aes(x=FT, y=mean.mx.sub),linetype=1, lwd=W, color=M)+
    teruego+ 
    scale_x_continuous(name="Time in Flight (mins)",
                       breaks=seq(0,720,120),
                       limits=c(0,720))+
    scale_y_continuous(name="Altitude (m)",labels=fancy_scientific)+
    coord_cartesian(ylim = c(0, max(leningrad[[3]][leningrad[[3]]$FT <= mins,]$mean.mx.sub))) 
  
  #### Ok, so I am gonna take kmtv, kmtv.sub, dep.10, dep.sub, hvt, hvt.sub
  X=4.5
  Y=4
  p<-0.5

  if (ii==1){
    ggsave("ALT_kmvt_BH50.pdf", plot = kmvt, device="pdf",width = X, height = Y, units = c("in"),dpi = 300)
    ggsave("ALT_dep10_BH50.pdf", plot = dep.10, device="pdf",width = X, height = Y, units = c("in"),dpi = 300)
    ggsave("ALT_hvt_BH50.pdf", plot = hvt, device="pdf",width = X, height = Y, units = c("in"),dpi = 300)
    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    ggsave("ALT_kmvt.sub_BH50.pdf", plot = kmvt.sub, device="pdf",width = p*X, height = p*Y, units = c("in"),dpi = 300)
    ggsave("ALT_dep10.sub_BH50.pdf", plot = dep.sub, device="pdf",width = p*X, height = p*Y, units = c("in"),dpi = 300)
    ggsave("ALT_hvt.sub_BH50.pdf", plot = hvt.sub, device="pdf",width = p*X, height = p*Y, units = c("in"),dpi = 300)
  } else if(ii==2){
    ggsave("ALT_kmvt_KC50.pdf", plot = kmvt, device="pdf",width = X, height = Y, units = c("in"),dpi = 300)
    ggsave("ALT_dep10_KC50.pdf", plot = dep.10, device="pdf",width = X, height = Y, units = c("in"),dpi = 300)
    ggsave("ALT_hvt_KC50.pdf", plot = hvt, device="pdf",width = X, height = Y, units = c("in"),dpi = 300)
    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    ggsave("ALT_kmvt.sub_KC50.pdf", plot = kmvt.sub, device="pdf",width = p*X, height = p*Y, units = c("in"),dpi = 300)
    ggsave("ALT_dep10.sub_KC50.pdf", plot = dep.sub, device="pdf",width = p*X, height = p*Y, units = c("in"),dpi = 300)
    ggsave("ALT_hvt.sub_KC50.pdf", plot = hvt.sub, device="pdf",width = p*X, height = p*Y, units = c("in"),dpi = 300)
  } else if(ii==3){
    ggsave("SOL_kmvt_BH50.pdf", plot = kmvt, device="pdf",width = X, height = Y, units = c("in"),dpi = 300)
    ggsave("SOL_dep10_BH50.pdf", plot = dep.10, device="pdf",width = X, height = Y, units = c("in"),dpi = 300)
    ggsave("SOL_hvt_BH50.pdf", plot = hvt, device="pdf",width = X, height = Y, units = c("in"),dpi = 300)
    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    ggsave("SOL_kmvt.sub_BH50.pdf", plot = kmvt.sub, device="pdf",width = p*X, height = p*Y, units = c("in"),dpi = 300)
    ggsave("SOL_dep10.sub_BH50.pdf", plot = dep.sub, device="pdf",width = p*X, height = p*Y, units = c("in"),dpi = 300)
    ggsave("SOL_hvt.sub_BH50.pdf", plot = hvt.sub, device="pdf",width = p*X, height = p*Y, units = c("in"),dpi = 300)
  } else if(ii==4){
    ggsave("SOL_kmvt_KC50.pdf", plot = kmvt, device="pdf",width = X, height = Y, units = c("in"),dpi = 300)
    ggsave("SOL_dep10_KC50.pdf", plot = dep.10, device="pdf",width = X, height = Y, units = c("in"),dpi = 300)
    ggsave("SOL_hvt_KC50.pdf", plot = hvt, device="pdf",width = X, height = Y, units = c("in"),dpi = 300)
    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    ggsave("SOL_kmvt.sub_KC50.pdf", plot = kmvt.sub, device="pdf",width = p*X, height = p*Y, units = c("in"),dpi = 300)
    ggsave("SOL_dep10.sub_KC50.pdf", plot = dep.sub, device="pdf",width = p*X, height = p*Y, units = c("in"),dpi = 300)
    ggsave("SOL_hvt.sub_KC50.pdf", plot = hvt.sub, device="pdf",width = p*X, height = p*Y, units = c("in"),dpi = 300)
  }
  

}



# ███████████████████████████████████████████████████████████████████████████████████████████████████
# ███████████████████████████████████████████████████████████████████████████████████████████████████
# ███████████████████████████████████████████████████████████████████████████████████████████████████
# ███████████████████████████████████████████████████████████████████████████████████████████████████



ALMOXARIFE<-list()
for (yyy in 1:4){
  #pb <- txtProgressBar(min = 0, max = length(mafuiles), style = 3)
  #mmx<-list()
  ii<-yyy
  mmx<-as.data.frame(fread(mafuiles[[ii]], verbose=F))
  colnames(mmx)[c(4,5)]<-c("FT","MX")
  mmx$VOO06<-floor(mmx$FT/360)
  mmx$VOO12<-floor(mmx$FT/720)
  mmx$VOO24<-floor(mmx$FT/1440)
  
  
  mods<-list()
  mods[[1]]<-mmx%>%
    group_by(HOUR) %>%
    group_split()
  bellaciao<-mods[[1]]
  #how many pieces am I splitting each day?
  entregalo<-4
  
  pb <- txtProgressBar(min = 0, max = 3 , style = 3)
  centrir<-list()
  centrirkey<-list()
  for (z in 1:length(bellaciao)){
    centrir[[z]] <- bellaciao[[z]] %>%   
      group_by(VOO06) %>%  #YEAR,MONTH,DAY,VOOINT
      group_split() 
    centrirkey[[z]] <- bellaciao[[z]] %>%
      group_by(VOO24) %>% #YEAR,MONTH,DAY,HOUR,VOOINT
      group_keys()
    setTxtProgressBar(pb, z)
  }
  
  
  ##################################
  
  kumura<-list()
  for (c in 1:length(centrir))
    kumura[[c]] <-pbapply::pblapply(centrir[[c]], function(x) 
      x %>%
        dplyr::group_by(VOO06) %>%
        dplyr::summarise(
          n=n(),
          mpH=mean(DIST/1000),
          sdH=sd(DIST/1000)
        )  %>%
        dplyr::mutate(seH=sdH/sqrt(n)))
  
  sodastars<-pbapply::pblapply(kumura, function(x) bind_rows(x))
  
  habibi<-list()
  for (c in 1:length(centrir))
    habibi[[c]] <-pbapply::pblapply(centrir[[c]], function(x) 
      x %>%
        dplyr::group_by(FT) %>%
        dplyr::summarise(
          n=n(),
          mpH=mean(DIST/1000),
          sdH=sd(DIST/1000),
          mpH.sub=mean(DIST),
          sdH.sub=sd(DIST)
        )  %>%
        dplyr::mutate(seH=sdH/sqrt(n)) %>%
        dplyr::mutate(seH.sub=sdH.sub/sqrt(n)))
  
  oumri<-pbapply::pblapply(habibi, function(x) bind_rows(x))
  
  #########################
  diamonds<-list()
  for (c in 1:length(centrir))
    diamonds[[c]] <-pbapply::pblapply(centrir[[c]], function(x) 
      x %>%
        dplyr::group_by(VOO06,YEAR) %>%
        dplyr::summarise(
          n=n()) )
  
  daemons<-pbapply::pblapply(diamonds, function(x) bind_rows(x))
  
  for (u in 1:length(daemons)){
    colnames(daemons[[u]])[3]<-"spores"
  }
  
  fatass<-pblapply(daemons, function(x) 
    x %>%
      dplyr::group_by(YEAR) %>%
      mutate(cumdump = cumsum(spores)))
  
  #########################
  
  katsura<-pblapply(fatass, function(x) 
    x %>%
      dplyr::group_by(VOO06) %>%
      dplyr::summarize(
        n=n(),
        meancum=mean(cumdump/20e5),
        sdcum=sd(cumdump/20e5)
      ) %>%
      mutate(secum=sdcum/sqrt(n)))
  
  # SMALL SMALL SMALL SMALL
  emeralds<-list()
  for (c in 1:length(centrir))
    emeralds[[c]] <-pbapply::pblapply(centrir[[c]], function(x) 
      x %>%
        dplyr::group_by(FT,YEAR) %>%
        dplyr::summarise(
          n=n()) )
  
  kraken<-pbapply::pblapply(emeralds, function(x) bind_rows(x))
  
  for (u in 1:length(kraken)){
    colnames(kraken[[u]])[3]<-"spores"
  }
  
  stab<-pblapply(kraken, function(x) 
    x %>%
      dplyr::group_by(YEAR) %>%
      mutate(cumdump = cumsum(spores)))
  
  fuckthepolice<-pblapply(stab, function(x) 
    x %>%
      dplyr::group_by(FT) %>%
      dplyr::summarize(
        n=n(),
        meancum=mean(cumdump/20e5),
        sdcum=sd(cumdump/20e5)
      ) %>%
      mutate(secum=sdcum/sqrt(n)))
  
 ############
  rooski<-list()
  for (c in 1:length(centrir))
    rooski[[c]] <-pbapply::pblapply(centrir[[c]], function(x) 
      x %>%
        dplyr::group_by(VOO06) %>%
        dplyr::summarise(
          n=n(),
          mean.mx=mean(MX/1000),
          sd.mx=sd(MX/1000)
        )  %>%
        dplyr::mutate(se.mx=sd.mx/sqrt(n)))
  
  stalingrad<-pbapply::pblapply(rooski, function(x) bind_rows(x))
  
  ##### FINE SCALE
  
  USSR<-list()
  for (c in 1:length(centrir))
    USSR[[c]] <-pbapply::pblapply(centrir[[c]], function(x) 
      x %>%
        dplyr::group_by(FT) %>%
        dplyr::summarise(
          n=n(),
          mean.mx=mean(MX),
          sd.mx=sd(MX),
          mean.mx.sub=mean(MX/1000),
          sd.mx.sub=sd(MX/1000)
        )  %>%
        dplyr::mutate(se.mx=sd.mx/sqrt(n)) %>%
        dplyr::mutate(se.mx.sub=sd.mx.sub/sqrt(n)))
  
  leningrad<-pbapply::pblapply(USSR, function(x) bind_rows(x))
  
  almoxarifado<-list()
  almoxarifado[[1]]<-sodastars
  almoxarifado[[2]]<-oumri
  
  almoxarifado[[3]]<-katsura
  almoxarifado[[4]]<-fuckthepolice
  
  almoxarifado[[5]]<-stalingrad
  almoxarifado[[6]]<-leningrad
  
  
  ALMOXARIFE[[yyy]]<-almoxarifado
}

# ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░
# ▄ ▄ ▄ ▄ ▄ ▄ ▄ ▄ ▄ ▄ ▄ ▄ ▄ ▄ ▄ ▄ ▄ ▄ ▄ ▄ ▄ ▄ ▄ ▄ ▄ ▄ ▄ ▄ ▄ ▄ ▄ ▄ ▄ ▄ ▄ ▄ 

for (gg in 1:2){
  #double
  Aa<-"#4F97A3"
  Ba<-"#000080"
  Ca<-"#0080FF"  
  As<-"#7E191B"
  Bs<-"#ED2929"
  Cs<-"#FF2400" 
  #alpha
  K<-25
  
  #line
  M<-"black"
  
  #line weight
  W<-.1
  n<-0.25
  
  #mins
  mins<-720
  
  tepido<-theme(panel.border = element_rect(linetype = 1, fill = NA),
                panel.grid.major = element_line(size=0.05,colour = "grey50"),
                panel.grid.minor = element_line(size=0.05,colour = "grey50"),
                panel.background = element_rect(fill = NA),
                axis.text.x = element_text(angle=0),
                axis.text.y = element_text( angle=0),
                axis.ticks.length=unit(0, "cm"))

  teruego<-theme(panel.border = element_rect(linetype = 1, fill = NA),
                 panel.grid.major = element_line(size=0.05,colour = "grey50"),
                 panel.grid.minor = element_line(size=0.05,colour = "grey50"),
                 panel.background = element_rect(fill = NA),
                 axis.text.x = element_text(angle=45, vjust=-0, hjust=0, size=7),
                 axis.text.y = element_text(angle=0, size=7),
                 axis.ticks.length=unit(0, "cm"))
  
  fancy_scientific <- function(l) {
    # turn in to character string in scientific notation
    l <- formatC(l, format = "e", digits = 2)
    # quote the part before the exponent to keep all the digits
    l <- gsub("^(.*)e", "'\\1'e", l)
    # turn the 'e+' into plotmath format
    l <- gsub("e", "%*%10^", l)
    # return this as an expression
    parse(text=l)
  }
  
  
  if(gg ==1){
    z<-1
    v<-3
  } else{
    z<-2
    v<-4
  }
  
  KMVT<-ggplot(data=NULL)+
    geom_ribbon(data=ALMOXARIFE[[v]][[1]][[1]],aes(x=VOO06,ymin = mpH - 5*seH, ymax = mpH + 5*seH), fill = As,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[1]][[1]],aes(x=VOO06,ymin = mpH - 4*seH, ymax = mpH + 4*seH), fill = As,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[1]][[1]],aes(x=VOO06,ymin = mpH - 3*seH, ymax = mpH + 3*seH), fill = As,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[1]][[1]],aes(x=VOO06,ymin = mpH - 2*seH, ymax = mpH + 2*seH), fill = As,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[1]][[1]],aes(x=VOO06,ymin = mpH - seH, ymax = mpH + seH), fill = As,alpha=5/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[1]][[2]],aes(x=VOO06,ymin = mpH - 5*seH, ymax = mpH + 5*seH), fill = Bs,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[1]][[2]],aes(x=VOO06,ymin = mpH - 4*seH, ymax = mpH + 4*seH), fill = Bs,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[1]][[2]],aes(x=VOO06,ymin = mpH - 3*seH, ymax = mpH + 3*seH), fill = Bs,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[1]][[2]],aes(x=VOO06,ymin = mpH - 2*seH, ymax = mpH + 2*seH), fill = Bs,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[1]][[2]],aes(x=VOO06,ymin = mpH - seH, ymax = mpH + seH), fill = Bs,alpha=5/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[1]][[3]],aes(x=VOO06,ymin = mpH - 5*seH, ymax = mpH + 5*seH), fill = Cs,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[1]][[3]],aes(x=VOO06,ymin = mpH - 4*seH, ymax = mpH + 4*seH), fill = Cs,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[1]][[3]],aes(x=VOO06,ymin = mpH - 3*seH, ymax = mpH + 3*seH), fill = Cs,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[1]][[3]],aes(x=VOO06,ymin = mpH - 2*seH, ymax = mpH + 2*seH), fill = Cs,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[1]][[3]],aes(x=VOO06,ymin = mpH - seH, ymax = mpH + seH), fill = Cs,alpha=5/K) +
    
    geom_ribbon(data=ALMOXARIFE[[z]][[1]][[1]],aes(x=VOO06,ymin = mpH - 5*seH, ymax = mpH + 5*seH), fill = Aa,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[1]][[1]],aes(x=VOO06,ymin = mpH - 4*seH, ymax = mpH + 4*seH), fill = Aa,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[1]][[1]],aes(x=VOO06,ymin = mpH - 3*seH, ymax = mpH + 3*seH), fill = Aa,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[1]][[1]],aes(x=VOO06,ymin = mpH - 2*seH, ymax = mpH + 2*seH), fill = Aa,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[1]][[1]],aes(x=VOO06,ymin = mpH - seH, ymax = mpH + seH), fill = Aa,alpha=5/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[1]][[2]],aes(x=VOO06,ymin = mpH - 5*seH, ymax = mpH + 5*seH), fill = Ba,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[1]][[2]],aes(x=VOO06,ymin = mpH - 4*seH, ymax = mpH + 4*seH), fill = Ba,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[1]][[2]],aes(x=VOO06,ymin = mpH - 3*seH, ymax = mpH + 3*seH), fill = Ba,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[1]][[2]],aes(x=VOO06,ymin = mpH - 2*seH, ymax = mpH + 2*seH), fill = Ba,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[1]][[2]],aes(x=VOO06,ymin = mpH - seH, ymax = mpH + seH), fill = Ba,alpha=5/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[1]][[3]],aes(x=VOO06,ymin = mpH - 5*seH, ymax = mpH + 5*seH), fill = Ca,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[1]][[3]],aes(x=VOO06,ymin = mpH - 4*seH, ymax = mpH + 4*seH), fill = Ca,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[1]][[3]],aes(x=VOO06,ymin = mpH - 3*seH, ymax = mpH + 3*seH), fill = Ca,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[1]][[3]],aes(x=VOO06,ymin = mpH - 2*seH, ymax = mpH + 2*seH), fill = Ca,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[1]][[3]],aes(x=VOO06,ymin = mpH - seH, ymax = mpH + seH), fill = Ca,alpha=5/K) +
    
    geom_line(data=ALMOXARIFE[[z]][[1]][[1]],aes(x=VOO06, y=mpH), linetype=1, lwd=W, color=M)+
    geom_line(data=ALMOXARIFE[[z]][[1]][[2]],aes(x=VOO06, y=mpH), linetype=1, lwd=W, color=M)+
    geom_line(data=ALMOXARIFE[[z]][[1]][[3]],aes(x=VOO06, y=mpH), linetype=1, lwd=W, color=M)+
    geom_line(data=ALMOXARIFE[[v]][[1]][[1]],aes(x=VOO06, y=mpH), linetype=1, lwd=W, color=M)+
    geom_line(data=ALMOXARIFE[[v]][[1]][[2]],aes(x=VOO06, y=mpH), linetype=1, lwd=W, color=M)+
    geom_line(data=ALMOXARIFE[[v]][[1]][[3]],aes(x=VOO06, y=mpH), linetype=1, lwd=W, color=M)+
    tepido+
    scale_x_continuous(name="Time in Flight [hrs]",
                       breaks=c(0,4,12,24,36,48),
                       labels=c(6,24,72,144,216,288),
                       limits = c(0,48))+
    scale_y_continuous(name="Distance Traveled [km]",labels=fancy_scientific)+
    coord_cartesian(ylim = c(0, 3500)) 
  
  KMVT.SUB<-ggplot(data=NULL)+
    geom_ribbon(data=ALMOXARIFE[[v]][[2]][[1]][ALMOXARIFE[[v]][[2]][[1]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 5*seH.sub, ymax = mpH.sub + 5*seH.sub), fill = As,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[2]][[1]][ALMOXARIFE[[v]][[2]][[1]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 4*seH.sub, ymax = mpH.sub + 4*seH.sub), fill = As,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[2]][[1]][ALMOXARIFE[[v]][[2]][[1]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 3*seH.sub, ymax = mpH.sub + 3*seH.sub), fill = As,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[2]][[1]][ALMOXARIFE[[v]][[2]][[1]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 2*seH.sub, ymax = mpH.sub + 2*seH.sub), fill = As,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[2]][[1]][ALMOXARIFE[[v]][[2]][[1]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - seH.sub, ymax = mpH.sub + seH.sub), fill = As,alpha=5/K) +
    
    geom_ribbon(data=ALMOXARIFE[[v]][[2]][[2]][ALMOXARIFE[[v]][[2]][[2]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 5*seH.sub, ymax = mpH.sub + 5*seH.sub), fill = Bs,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[2]][[2]][ALMOXARIFE[[v]][[2]][[2]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 4*seH.sub, ymax = mpH.sub + 4*seH.sub), fill = Bs,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[2]][[2]][ALMOXARIFE[[v]][[2]][[2]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 3*seH.sub, ymax = mpH.sub + 3*seH.sub), fill = Bs,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[2]][[2]][ALMOXARIFE[[v]][[2]][[2]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 2*seH.sub, ymax = mpH.sub + 2*seH.sub), fill = Bs,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[2]][[2]][ALMOXARIFE[[v]][[2]][[2]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - seH.sub, ymax = mpH.sub + seH.sub), fill = Bs,alpha=5/K) +
    
    geom_ribbon(data=ALMOXARIFE[[v]][[2]][[3]][ALMOXARIFE[[v]][[2]][[3]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 5*seH.sub, ymax = mpH.sub + 5*seH.sub), fill = Cs,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[2]][[3]][ALMOXARIFE[[v]][[2]][[3]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 4*seH.sub, ymax = mpH.sub + 4*seH.sub), fill = Cs,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[2]][[3]][ALMOXARIFE[[v]][[2]][[3]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 3*seH.sub, ymax = mpH.sub + 3*seH.sub), fill = Cs,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[2]][[3]][ALMOXARIFE[[v]][[2]][[3]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 2*seH.sub, ymax = mpH.sub + 2*seH.sub), fill = Cs,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[2]][[3]][ALMOXARIFE[[v]][[2]][[3]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - seH.sub, ymax = mpH.sub + seH.sub), fill = Cs,alpha=5/K) +
    
    geom_ribbon(data=ALMOXARIFE[[z]][[2]][[1]][ALMOXARIFE[[z]][[2]][[1]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 5*seH.sub, ymax = mpH.sub + 5*seH.sub), fill = Aa,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[2]][[1]][ALMOXARIFE[[z]][[2]][[1]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 4*seH.sub, ymax = mpH.sub + 4*seH.sub), fill = Aa,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[2]][[1]][ALMOXARIFE[[z]][[2]][[1]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 3*seH.sub, ymax = mpH.sub + 3*seH.sub), fill = Aa,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[2]][[1]][ALMOXARIFE[[z]][[2]][[1]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 2*seH.sub, ymax = mpH.sub + 2*seH.sub), fill = Aa,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[2]][[1]][ALMOXARIFE[[z]][[2]][[1]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - seH.sub, ymax = mpH.sub + seH.sub), fill = Aa,alpha=5/K) +
    
    geom_ribbon(data=ALMOXARIFE[[z]][[2]][[2]][ALMOXARIFE[[z]][[2]][[2]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 5*seH.sub, ymax = mpH.sub + 5*seH.sub), fill = Ba,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[2]][[2]][ALMOXARIFE[[z]][[2]][[2]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 4*seH.sub, ymax = mpH.sub + 4*seH.sub), fill = Ba,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[2]][[2]][ALMOXARIFE[[z]][[2]][[2]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 3*seH.sub, ymax = mpH.sub + 3*seH.sub), fill = Ba,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[2]][[2]][ALMOXARIFE[[z]][[2]][[2]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 2*seH.sub, ymax = mpH.sub + 2*seH.sub), fill = Ba,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[2]][[2]][ALMOXARIFE[[z]][[2]][[2]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - seH.sub, ymax = mpH.sub + seH.sub), fill = Ba,alpha=5/K) +
    
    geom_ribbon(data=ALMOXARIFE[[z]][[2]][[3]][ALMOXARIFE[[z]][[2]][[3]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 5*seH.sub, ymax = mpH.sub + 5*seH.sub), fill = Ca,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[2]][[3]][ALMOXARIFE[[z]][[2]][[3]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 4*seH.sub, ymax = mpH.sub + 4*seH.sub), fill = Ca,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[2]][[3]][ALMOXARIFE[[z]][[2]][[3]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 3*seH.sub, ymax = mpH.sub + 3*seH.sub), fill = Ca,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[2]][[3]][ALMOXARIFE[[z]][[2]][[3]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - 2*seH.sub, ymax = mpH.sub + 2*seH.sub), fill = Ca,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[2]][[3]][ALMOXARIFE[[z]][[2]][[3]]$FT <=mins,],aes(x=FT,ymin = mpH.sub - seH.sub, ymax = mpH.sub + seH.sub), fill = Ca,alpha=5/K) +
    
    geom_line(data=ALMOXARIFE[[v]][[2]][[1]][ALMOXARIFE[[v]][[2]][[1]]$FT <=mins,],aes(x=FT, y=mpH.sub), linetype=1, lwd=W*n, color=M)+
    geom_line(data=ALMOXARIFE[[v]][[2]][[2]][ALMOXARIFE[[v]][[2]][[2]]$FT <=mins,],aes(x=FT, y=mpH.sub), linetype=1, lwd=W*n, color=M)+
    geom_line(data=ALMOXARIFE[[v]][[2]][[3]][ALMOXARIFE[[v]][[2]][[3]]$FT <=mins,],aes(x=FT, y=mpH.sub), linetype=1, lwd=W*n, color=M)+
    
    geom_line(data=ALMOXARIFE[[z]][[2]][[1]][ALMOXARIFE[[z]][[2]][[1]]$FT <=mins,],aes(x=FT, y=mpH.sub), linetype=1, lwd=W*n, color=M)+
    geom_line(data=ALMOXARIFE[[z]][[2]][[2]][ALMOXARIFE[[z]][[2]][[2]]$FT <=mins,],aes(x=FT, y=mpH.sub), linetype=1, lwd=W*n, color=M)+
    geom_line(data=ALMOXARIFE[[z]][[2]][[3]][ALMOXARIFE[[z]][[2]][[3]]$FT <=mins,],aes(x=FT, y=mpH.sub), linetype=1, lwd=W*n, color=M)+
    
    teruego+
    scale_x_continuous(name="",
                       breaks=seq(0,720,120),
                       labels=paste0(seq(0,720,120),rep("min",720/120)),
                       limits=c(0,720),
                       position = "top")+
    scale_y_continuous(name="[m]",labels=fancy_scientific)+
    coord_cartesian(ylim = c(0, 3.5e5)) 
  
  
  # ▀▄░▄▀ ▀▄░▄▀ ▀▄░▄▀ ▀▄░▄▀ ▀▄░▄▀ ▀▄░▄▀ ▀▄░▄▀ ▀▄░▄▀ ▀▄░▄▀ ▀▄░▄▀ ▀▄░▄▀ 
  # ─░█── ─░█── ─░█── ─░█── ─░█── ─░█── ─░█── ─░█── ─░█── ─░█── ─░█── 
  # ▄▀░▀▄ ▄▀░▀▄ ▄▀░▀▄ ▄▀░▀▄ ▄▀░▀▄ ▄▀░▀▄ ▄▀░▀▄ ▄▀░▀▄ ▄▀░▀▄ ▄▀░▀▄ ▄▀░▀▄
  
  
  
  DEP.10<-ggplot(data=NULL)+
    geom_ribbon(data=ALMOXARIFE[[v]][[4]][[1]],aes(x=FT,ymin = meancum - 5*secum, ymax = meancum + 5*secum), fill = As,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[4]][[1]],aes(x=FT,ymin = meancum - 4*secum, ymax = meancum + 4*secum), fill = As,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[4]][[1]],aes(x=FT,ymin = meancum - 3*secum, ymax = meancum + 3*secum), fill = As,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[4]][[1]],aes(x=FT,ymin = meancum - 2*secum, ymax = meancum + 2*secum), fill = As,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[4]][[1]],aes(x=FT,ymin = meancum - secum,   ymax = meancum + secum), fill=    As,alpha=5/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[4]][[2]],aes(x=FT,ymin = meancum - 5*secum, ymax = meancum + 5*secum), fill = Bs,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[4]][[2]],aes(x=FT,ymin = meancum - 4*secum, ymax = meancum + 4*secum), fill = Bs,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[4]][[2]],aes(x=FT,ymin = meancum - 3*secum, ymax = meancum + 3*secum), fill = Bs,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[4]][[2]],aes(x=FT,ymin = meancum - 2*secum, ymax = meancum + 2*secum), fill = Bs,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[4]][[2]],aes(x=FT,ymin = meancum - secum,   ymax = meancum + secum), fill =   Bs,alpha=5/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[4]][[3]],aes(x=FT,ymin = meancum - 5*secum, ymax = meancum + 5*secum), fill = Cs,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[4]][[3]],aes(x=FT,ymin = meancum - 4*secum, ymax = meancum + 4*secum), fill = Cs,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[4]][[3]],aes(x=FT,ymin = meancum - 3*secum, ymax = meancum + 3*secum), fill = Cs,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[4]][[3]],aes(x=FT,ymin = meancum - 2*secum, ymax = meancum + 2*secum), fill = Cs,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[4]][[3]],aes(x=FT,ymin = meancum - secum,   ymax = meancum + secum), fill = Cs,alpha=  5/K) +
    
    geom_ribbon(data=ALMOXARIFE[[z]][[4]][[1]],aes(x=FT,ymin = meancum - 5*secum, ymax = meancum + 5*secum), fill = Aa,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[4]][[1]],aes(x=FT,ymin = meancum - 4*secum, ymax = meancum + 4*secum), fill = Aa,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[4]][[1]],aes(x=FT,ymin = meancum - 3*secum, ymax = meancum + 3*secum), fill = Aa,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[4]][[1]],aes(x=FT,ymin = meancum - 2*secum, ymax = meancum + 2*secum), fill = Aa,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[4]][[1]],aes(x=FT,ymin = meancum - secum,   ymax = meancum + secum), fill=    Aa,alpha=5/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[4]][[2]],aes(x=FT,ymin = meancum - 5*secum, ymax = meancum + 5*secum), fill = Ba,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[4]][[2]],aes(x=FT,ymin = meancum - 4*secum, ymax = meancum + 4*secum), fill = Ba,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[4]][[2]],aes(x=FT,ymin = meancum - 3*secum, ymax = meancum + 3*secum), fill = Ba,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[4]][[2]],aes(x=FT,ymin = meancum - 2*secum, ymax = meancum + 2*secum), fill = Ba,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[4]][[2]],aes(x=FT,ymin = meancum - secum,   ymax = meancum + secum), fill =   Ba,alpha=5/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[4]][[3]],aes(x=FT,ymin = meancum - 5*secum, ymax = meancum + 5*secum), fill = Ca,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[4]][[3]],aes(x=FT,ymin = meancum - 4*secum, ymax = meancum + 4*secum), fill = Ca,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[4]][[3]],aes(x=FT,ymin = meancum - 3*secum, ymax = meancum + 3*secum), fill = Ca,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[4]][[3]],aes(x=FT,ymin = meancum - 2*secum, ymax = meancum + 2*secum), fill = Ca,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[4]][[3]],aes(x=FT,ymin = meancum - secum,   ymax = meancum + secum), fill = Ca,alpha=  5/K) +
    
    geom_line(data=ALMOXARIFE[[v]][[4]][[1]],aes(x=FT,y=meancum),linetype=1, lwd=W, color=M)+
    geom_line(data=ALMOXARIFE[[v]][[4]][[2]],aes(x=FT,y=meancum),linetype=1, lwd=W, color=M)+
    geom_line(data=ALMOXARIFE[[v]][[4]][[3]],aes(x=FT,y=meancum),linetype=1, lwd=W, color=M)+
    geom_line(data=ALMOXARIFE[[z]][[4]][[1]],aes(x=FT,y=meancum),linetype=1, lwd=W, color=M)+
    geom_line(data=ALMOXARIFE[[z]][[4]][[2]],aes(x=FT,y=meancum),linetype=1, lwd=W, color=M)+
    geom_line(data=ALMOXARIFE[[z]][[4]][[3]],aes(x=FT,y=meancum),linetype=1, lwd=W, color=M)+
    tepido+
    scale_x_continuous(name="Time in Flight [hrs]",
                       breaks=c(0,1440,4320,8640,12960,17280),
                       labels=c(6,24,72,144,216,288),
                       limits=c(0,17280))+
    scale_y_reverse(name="Spores Remaining Airborne",labels=fancy_scientific)+
    coord_cartesian(ylim = c(0, 1)) 
  
  
  DEP.SUB<-ggplot(data=NULL)+
    geom_ribbon(data=ALMOXARIFE[[v]][[4]][[1]][ALMOXARIFE[[v]][[4]][[1]]$FT <=mins,],aes(x=FT,ymin = meancum - 5*secum, ymax = meancum + 5*secum), fill = As,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[4]][[1]][ALMOXARIFE[[v]][[4]][[1]]$FT <=mins,],aes(x=FT,ymin = meancum - 4*secum, ymax = meancum + 4*secum), fill = As,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[4]][[1]][ALMOXARIFE[[v]][[4]][[1]]$FT <=mins,],aes(x=FT,ymin = meancum - 3*secum, ymax = meancum + 3*secum), fill = As,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[4]][[1]][ALMOXARIFE[[v]][[4]][[1]]$FT <=mins,],aes(x=FT,ymin = meancum - 2*secum, ymax = meancum + 2*secum), fill = As,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[4]][[1]][ALMOXARIFE[[v]][[4]][[1]]$FT <=mins,],aes(x=FT,ymin = meancum - secum,   ymax = meancum + secum), fill=    As,alpha=5/K) +
    
    geom_ribbon(data=ALMOXARIFE[[v]][[4]][[2]][ALMOXARIFE[[v]][[4]][[2]]$FT <=mins,],aes(x=FT,ymin = meancum - 5*secum, ymax = meancum + 5*secum), fill = Bs,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[4]][[2]][ALMOXARIFE[[v]][[4]][[2]]$FT <=mins,],aes(x=FT,ymin = meancum - 4*secum, ymax = meancum + 4*secum), fill = Bs,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[4]][[2]][ALMOXARIFE[[v]][[4]][[2]]$FT <=mins,],aes(x=FT,ymin = meancum - 3*secum, ymax = meancum + 3*secum), fill = Bs,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[4]][[2]][ALMOXARIFE[[v]][[4]][[2]]$FT <=mins,],aes(x=FT,ymin = meancum - 2*secum, ymax = meancum + 2*secum), fill = Bs,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[4]][[2]][ALMOXARIFE[[v]][[4]][[2]]$FT <=mins,],aes(x=FT,ymin = meancum - secum,   ymax = meancum + secum), fill =   Bs,alpha=5/K) +
    
    geom_ribbon(data=ALMOXARIFE[[v]][[4]][[3]][ALMOXARIFE[[v]][[4]][[3]]$FT <=mins,],aes(x=FT,ymin = meancum - 5*secum, ymax = meancum + 5*secum), fill = Cs,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[4]][[3]][ALMOXARIFE[[v]][[4]][[3]]$FT <=mins,],aes(x=FT,ymin = meancum - 4*secum, ymax = meancum + 4*secum), fill = Cs,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[4]][[3]][ALMOXARIFE[[v]][[4]][[3]]$FT <=mins,],aes(x=FT,ymin = meancum - 3*secum, ymax = meancum + 3*secum), fill = Cs,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[4]][[3]][ALMOXARIFE[[v]][[4]][[3]]$FT <=mins,],aes(x=FT,ymin = meancum - 2*secum, ymax = meancum + 2*secum), fill = Cs,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[4]][[3]][ALMOXARIFE[[v]][[4]][[3]]$FT <=mins,],aes(x=FT,ymin = meancum - secum,   ymax = meancum + secum), fill = Cs,alpha=  5/K) +
    
    geom_ribbon(data=ALMOXARIFE[[z]][[4]][[1]][ALMOXARIFE[[z]][[4]][[1]]$FT <=mins,],aes(x=FT,ymin = meancum - 5*secum, ymax = meancum + 5*secum), fill = Aa,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[4]][[1]][ALMOXARIFE[[z]][[4]][[1]]$FT <=mins,],aes(x=FT,ymin = meancum - 4*secum, ymax = meancum + 4*secum), fill = Aa,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[4]][[1]][ALMOXARIFE[[z]][[4]][[1]]$FT <=mins,],aes(x=FT,ymin = meancum - 3*secum, ymax = meancum + 3*secum), fill = Aa,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[4]][[1]][ALMOXARIFE[[z]][[4]][[1]]$FT <=mins,],aes(x=FT,ymin = meancum - 2*secum, ymax = meancum + 2*secum), fill = Aa,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[4]][[1]][ALMOXARIFE[[z]][[4]][[1]]$FT <=mins,],aes(x=FT,ymin = meancum - secum,   ymax = meancum + secum), fill=    Aa,alpha=5/K) +
    
    geom_ribbon(data=ALMOXARIFE[[z]][[4]][[2]][ALMOXARIFE[[z]][[4]][[2]]$FT <=mins,],aes(x=FT,ymin = meancum - 5*secum, ymax = meancum + 5*secum), fill = Ba,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[4]][[2]][ALMOXARIFE[[z]][[4]][[2]]$FT <=mins,],aes(x=FT,ymin = meancum - 4*secum, ymax = meancum + 4*secum), fill = Ba,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[4]][[2]][ALMOXARIFE[[z]][[4]][[2]]$FT <=mins,],aes(x=FT,ymin = meancum - 3*secum, ymax = meancum + 3*secum), fill = Ba,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[4]][[2]][ALMOXARIFE[[z]][[4]][[2]]$FT <=mins,],aes(x=FT,ymin = meancum - 2*secum, ymax = meancum + 2*secum), fill = Ba,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[4]][[2]][ALMOXARIFE[[z]][[4]][[2]]$FT <=mins,],aes(x=FT,ymin = meancum - secum,   ymax = meancum + secum), fill =   Ba,alpha=5/K) +
    
    geom_ribbon(data=ALMOXARIFE[[z]][[4]][[3]][ALMOXARIFE[[z]][[4]][[3]]$FT <=mins,],aes(x=FT,ymin = meancum - 5*secum, ymax = meancum + 5*secum), fill = Ca,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[4]][[3]][ALMOXARIFE[[z]][[4]][[3]]$FT <=mins,],aes(x=FT,ymin = meancum - 4*secum, ymax = meancum + 4*secum), fill = Ca,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[4]][[3]][ALMOXARIFE[[z]][[4]][[3]]$FT <=mins,],aes(x=FT,ymin = meancum - 3*secum, ymax = meancum + 3*secum), fill = Ca,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[4]][[3]][ALMOXARIFE[[z]][[4]][[3]]$FT <=mins,],aes(x=FT,ymin = meancum - 2*secum, ymax = meancum + 2*secum), fill = Ca,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[4]][[3]][ALMOXARIFE[[z]][[4]][[3]]$FT <=mins,],aes(x=FT,ymin = meancum - secum,   ymax = meancum + secum), fill = Ca,alpha=  5/K) +
    
    geom_line(data=ALMOXARIFE[[v]][[4]][[1]][ALMOXARIFE[[v]][[4]][[1]]$FT <=mins,],aes(x=FT,y=meancum),linetype=1, lwd=W, color=M)+
    geom_line(data=ALMOXARIFE[[v]][[4]][[2]][ALMOXARIFE[[v]][[4]][[2]]$FT <=mins,],aes(x=FT,y=meancum),linetype=1, lwd=W, color=M)+
    geom_line(data=ALMOXARIFE[[v]][[4]][[3]][ALMOXARIFE[[v]][[4]][[3]]$FT <=mins,],aes(x=FT,y=meancum),linetype=1, lwd=W, color=M)+
    
    geom_line(data=ALMOXARIFE[[z]][[4]][[1]][ALMOXARIFE[[z]][[4]][[1]]$FT <=mins,],aes(x=FT,y=meancum),linetype=1, lwd=W*n, color=M)+
    geom_line(data=ALMOXARIFE[[z]][[4]][[2]][ALMOXARIFE[[z]][[4]][[2]]$FT <=mins,],aes(x=FT,y=meancum),linetype=1, lwd=W*n, color=M)+
    geom_line(data=ALMOXARIFE[[z]][[4]][[3]][ALMOXARIFE[[z]][[4]][[3]]$FT <=mins,],aes(x=FT,y=meancum),linetype=1, lwd=W*n, color=M)+
    teruego+
    scale_x_continuous(name="",
                       breaks=seq(0,720,120),
                       labels=paste0(seq(0,720,120),rep("min",720/120)),
                       limits = c(0,720),
                       position="bottom")+
    scale_y_reverse(name="Spores...",labels=fancy_scientific)+
    coord_cartesian(ylim = c(0, 1)) 
  
  
  # ▀▄░▄▀ ▀▄░▄▀ ▀▄░▄▀ ▀▄░▄▀ ▀▄░▄▀ ▀▄░▄▀ ▀▄░▄▀ ▀▄░▄▀ ▀▄░▄▀ ▀▄░▄▀ ▀▄░▄▀ 
  # ─░█── ─░█── ─░█── ─░█── ─░█── ─░█── ─░█── ─░█── ─░█── ─░█── ─░█── 
  # ▄▀░▀▄ ▄▀░▀▄ ▄▀░▀▄ ▄▀░▀▄ ▄▀░▀▄ ▄▀░▀▄ ▄▀░▀▄ ▄▀░▀▄ ▄▀░▀▄ ▄▀░▀▄ ▄▀░▀▄
  
  
  HVT<-ggplot(data=NULL)+
    geom_ribbon(data=ALMOXARIFE[[v]][[5]][[1]],aes(x=VOO06,ymin = mean.mx - 5*se.mx, ymax = mean.mx + 5*se.mx), fill = As,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[5]][[1]],aes(x=VOO06,ymin = mean.mx - 4*se.mx, ymax = mean.mx + 4*se.mx), fill = As,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[5]][[1]],aes(x=VOO06,ymin = mean.mx - 3*se.mx, ymax = mean.mx + 3*se.mx), fill = As,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[5]][[1]],aes(x=VOO06,ymin = mean.mx - 2*se.mx, ymax = mean.mx + 2*se.mx), fill = As,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[5]][[1]],aes(x=VOO06,ymin = mean.mx - se.mx, ymax = mean.mx + se.mx), fill = As,alpha=5/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[5]][[2]],aes(x=VOO06,ymin = mean.mx - 5*se.mx, ymax = mean.mx + 5*se.mx), fill = Bs,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[5]][[2]],aes(x=VOO06,ymin = mean.mx - 4*se.mx, ymax = mean.mx + 4*se.mx), fill = Bs,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[5]][[2]],aes(x=VOO06,ymin = mean.mx - 3*se.mx, ymax = mean.mx + 3*se.mx), fill = Bs,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[5]][[2]],aes(x=VOO06,ymin = mean.mx - 2*se.mx, ymax = mean.mx + 2*se.mx), fill = Bs,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[5]][[2]],aes(x=VOO06,ymin = mean.mx - se.mx, ymax = mean.mx + se.mx), fill = Bs,alpha=5/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[5]][[3]],aes(x=VOO06,ymin = mean.mx - 5*se.mx, ymax = mean.mx + 5*se.mx), fill = Cs,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[5]][[3]],aes(x=VOO06,ymin = mean.mx - 4*se.mx, ymax = mean.mx + 4*se.mx), fill = Cs,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[5]][[3]],aes(x=VOO06,ymin = mean.mx - 3*se.mx, ymax = mean.mx + 3*se.mx), fill = Cs,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[5]][[3]],aes(x=VOO06,ymin = mean.mx - 2*se.mx, ymax = mean.mx + 2*se.mx), fill = Cs,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[5]][[3]],aes(x=VOO06,ymin = mean.mx - se.mx, ymax = mean.mx + se.mx), fill = Cs,alpha=5/K) +
    
    geom_ribbon(data=ALMOXARIFE[[z]][[5]][[1]],aes(x=VOO06,ymin = mean.mx - 5*se.mx, ymax = mean.mx + 5*se.mx), fill = Aa,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[5]][[1]],aes(x=VOO06,ymin = mean.mx - 4*se.mx, ymax = mean.mx + 4*se.mx), fill = Aa,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[5]][[1]],aes(x=VOO06,ymin = mean.mx - 3*se.mx, ymax = mean.mx + 3*se.mx), fill = Aa,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[5]][[1]],aes(x=VOO06,ymin = mean.mx - 2*se.mx, ymax = mean.mx + 2*se.mx), fill = Aa,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[5]][[1]],aes(x=VOO06,ymin = mean.mx - se.mx, ymax = mean.mx + se.mx), fill = Aa,alpha=5/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[5]][[2]],aes(x=VOO06,ymin = mean.mx - 5*se.mx, ymax = mean.mx + 5*se.mx), fill = Ba,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[5]][[2]],aes(x=VOO06,ymin = mean.mx - 4*se.mx, ymax = mean.mx + 4*se.mx), fill = Ba,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[5]][[2]],aes(x=VOO06,ymin = mean.mx - 3*se.mx, ymax = mean.mx + 3*se.mx), fill = Ba,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[5]][[2]],aes(x=VOO06,ymin = mean.mx - 2*se.mx, ymax = mean.mx + 2*se.mx), fill = Ba,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[5]][[2]],aes(x=VOO06,ymin = mean.mx - se.mx, ymax = mean.mx + se.mx), fill = Ba,alpha=5/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[5]][[3]],aes(x=VOO06,ymin = mean.mx - 5*se.mx, ymax = mean.mx + 5*se.mx), fill = Ca,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[5]][[3]],aes(x=VOO06,ymin = mean.mx - 4*se.mx, ymax = mean.mx + 4*se.mx), fill = Ca,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[5]][[3]],aes(x=VOO06,ymin = mean.mx - 3*se.mx, ymax = mean.mx + 3*se.mx), fill = Ca,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[5]][[3]],aes(x=VOO06,ymin = mean.mx - 2*se.mx, ymax = mean.mx + 2*se.mx), fill = Ca,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[5]][[3]],aes(x=VOO06,ymin = mean.mx - se.mx, ymax = mean.mx + se.mx), fill = Ca,alpha=5/K) +
    
    geom_line(data=ALMOXARIFE[[z]][[5]][[1]],aes(x=VOO06, y=mean.mx), linetype=1, lwd=W, color=M)+
    geom_line(data=ALMOXARIFE[[z]][[5]][[2]],aes(x=VOO06, y=mean.mx), linetype=1, lwd=W, color=M)+
    geom_line(data=ALMOXARIFE[[z]][[5]][[3]],aes(x=VOO06, y=mean.mx), linetype=1, lwd=W, color=M)+
    geom_line(data=ALMOXARIFE[[v]][[5]][[1]],aes(x=VOO06, y=mean.mx), linetype=1, lwd=W, color=M)+
    geom_line(data=ALMOXARIFE[[v]][[5]][[2]],aes(x=VOO06, y=mean.mx), linetype=1, lwd=W, color=M)+
    geom_line(data=ALMOXARIFE[[v]][[5]][[3]],aes(x=VOO06, y=mean.mx), linetype=1, lwd=W, color=M)+
    tepido+
    scale_x_continuous(name="Time in Flight (hrs)",
                       breaks=c(0,4,12,24,36,48),
                       labels=c(0,24,72,144,216,288),
                       limits = c(0,48))+
    scale_y_continuous(name="Altitude [km]",labels=fancy_scientific)+
    coord_cartesian(ylim = c(0, 7),xlim=c(0,48)) 
  
  
  HVT.SUB<-ggplot(data=NULL)+
    geom_ribbon(data=ALMOXARIFE[[v]][[6]][[1]][ALMOXARIFE[[v]][[6]][[1]]$FT <=mins,],aes(x=FT,ymin = mean.mx.sub - 5*se.mx.sub, ymax = mean.mx.sub + 5*se.mx.sub), fill = As,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[6]][[1]][ALMOXARIFE[[v]][[6]][[1]]$FT <=mins,],aes(x=FT,ymin = mean.mx.sub - 4*se.mx.sub, ymax = mean.mx.sub + 4*se.mx.sub), fill = As,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[6]][[1]][ALMOXARIFE[[v]][[6]][[1]]$FT <=mins,],aes(x=FT,ymin = mean.mx.sub - 3*se.mx.sub, ymax = mean.mx.sub + 3*se.mx.sub), fill = As,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[6]][[1]][ALMOXARIFE[[v]][[6]][[1]]$FT <=mins,],aes(x=FT,ymin = mean.mx.sub - 2*se.mx.sub, ymax = mean.mx.sub + 2*se.mx.sub), fill = As,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[6]][[1]][ALMOXARIFE[[v]][[6]][[1]]$FT <=mins,],aes(x=FT,ymin = mean.mx.sub - se.mx.sub, ymax = mean.mx.sub + se.mx.sub), fill = As,alpha=5/K) +
    
    geom_ribbon(data=ALMOXARIFE[[v]][[6]][[2]][ALMOXARIFE[[v]][[6]][[2]]$FT <=mins,],aes(x=FT,ymin = mean.mx.sub - 5*se.mx.sub, ymax = mean.mx.sub + 5*se.mx.sub), fill = Bs,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[6]][[2]][ALMOXARIFE[[v]][[6]][[2]]$FT <=mins,],aes(x=FT,ymin = mean.mx.sub - 4*se.mx.sub, ymax = mean.mx.sub + 4*se.mx.sub), fill = Bs,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[6]][[2]][ALMOXARIFE[[v]][[6]][[2]]$FT <=mins,],aes(x=FT,ymin = mean.mx.sub - 3*se.mx.sub, ymax = mean.mx.sub + 3*se.mx.sub), fill = Bs,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[6]][[2]][ALMOXARIFE[[v]][[6]][[2]]$FT <=mins,],aes(x=FT,ymin = mean.mx.sub - 2*se.mx.sub, ymax = mean.mx.sub + 2*se.mx.sub), fill = Bs,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[6]][[2]][ALMOXARIFE[[v]][[6]][[2]]$FT <=mins,],aes(x=FT,ymin = mean.mx.sub - se.mx.sub, ymax = mean.mx.sub + se.mx.sub), fill = Bs,alpha=5/K) +
    
    geom_ribbon(data=ALMOXARIFE[[v]][[6]][[3]][ALMOXARIFE[[v]][[6]][[3]]$FT <=mins,],aes(x=FT,ymin = mean.mx.sub - 5*se.mx.sub, ymax = mean.mx.sub + 5*se.mx.sub), fill = Cs,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[6]][[3]][ALMOXARIFE[[v]][[6]][[3]]$FT <=mins,],aes(x=FT,ymin = mean.mx.sub - 4*se.mx.sub, ymax = mean.mx.sub + 4*se.mx.sub), fill = Cs,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[6]][[3]][ALMOXARIFE[[v]][[6]][[3]]$FT <=mins,],aes(x=FT,ymin = mean.mx.sub - 3*se.mx.sub, ymax = mean.mx.sub + 3*se.mx.sub), fill = Cs,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[6]][[3]][ALMOXARIFE[[v]][[6]][[3]]$FT <=mins,],aes(x=FT,ymin = mean.mx.sub - 2*se.mx.sub, ymax = mean.mx.sub + 2*se.mx.sub), fill = Cs,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[v]][[6]][[3]][ALMOXARIFE[[v]][[6]][[3]]$FT <=mins,],aes(x=FT,ymin = mean.mx.sub - se.mx.sub, ymax = mean.mx.sub + se.mx.sub), fill = Cs,alpha=5/K) +
    
    geom_ribbon(data=ALMOXARIFE[[z]][[6]][[1]][ALMOXARIFE[[z]][[6]][[1]]$FT <=mins,],aes(x=FT,ymin = mean.mx.sub - 5*se.mx.sub, ymax = mean.mx.sub + 5*se.mx.sub), fill = Aa,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[6]][[1]][ALMOXARIFE[[z]][[6]][[1]]$FT <=mins,],aes(x=FT,ymin = mean.mx.sub - 4*se.mx.sub, ymax = mean.mx.sub + 4*se.mx.sub), fill = Aa,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[6]][[1]][ALMOXARIFE[[z]][[6]][[1]]$FT <=mins,],aes(x=FT,ymin = mean.mx.sub - 3*se.mx.sub, ymax = mean.mx.sub + 3*se.mx.sub), fill = Aa,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[6]][[1]][ALMOXARIFE[[z]][[6]][[1]]$FT <=mins,],aes(x=FT,ymin = mean.mx.sub - 2*se.mx.sub, ymax = mean.mx.sub + 2*se.mx.sub), fill = Aa,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[6]][[1]][ALMOXARIFE[[z]][[6]][[1]]$FT <=mins,],aes(x=FT,ymin = mean.mx.sub - se.mx.sub, ymax = mean.mx.sub + se.mx.sub), fill = Aa,alpha=5/K) +
    
    geom_ribbon(data=ALMOXARIFE[[z]][[6]][[2]][ALMOXARIFE[[z]][[6]][[2]]$FT <=mins,],aes(x=FT,ymin = mean.mx.sub - 5*se.mx.sub, ymax = mean.mx.sub + 5*se.mx.sub), fill = Ba,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[6]][[2]][ALMOXARIFE[[z]][[6]][[2]]$FT <=mins,],aes(x=FT,ymin = mean.mx.sub - 4*se.mx.sub, ymax = mean.mx.sub + 4*se.mx.sub), fill = Ba,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[6]][[2]][ALMOXARIFE[[z]][[6]][[2]]$FT <=mins,],aes(x=FT,ymin = mean.mx.sub - 3*se.mx.sub, ymax = mean.mx.sub + 3*se.mx.sub), fill = Ba,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[6]][[2]][ALMOXARIFE[[z]][[6]][[2]]$FT <=mins,],aes(x=FT,ymin = mean.mx.sub - 2*se.mx.sub, ymax = mean.mx.sub + 2*se.mx.sub), fill = Ba,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[6]][[2]][ALMOXARIFE[[z]][[6]][[2]]$FT <=mins,],aes(x=FT,ymin = mean.mx.sub - se.mx.sub, ymax = mean.mx.sub + se.mx.sub), fill = Ba,alpha=5/K) +
    
    geom_ribbon(data=ALMOXARIFE[[z]][[6]][[3]][ALMOXARIFE[[z]][[6]][[3]]$FT <=mins,],aes(x=FT,ymin = mean.mx.sub - 5*se.mx.sub, ymax = mean.mx.sub + 5*se.mx.sub), fill = Ca,alpha=1/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[6]][[3]][ALMOXARIFE[[z]][[6]][[3]]$FT <=mins,],aes(x=FT,ymin = mean.mx.sub - 4*se.mx.sub, ymax = mean.mx.sub + 4*se.mx.sub), fill = Ca,alpha=2/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[6]][[3]][ALMOXARIFE[[z]][[6]][[3]]$FT <=mins,],aes(x=FT,ymin = mean.mx.sub - 3*se.mx.sub, ymax = mean.mx.sub + 3*se.mx.sub), fill = Ca,alpha=3/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[6]][[3]][ALMOXARIFE[[z]][[6]][[3]]$FT <=mins,],aes(x=FT,ymin = mean.mx.sub - 2*se.mx.sub, ymax = mean.mx.sub + 2*se.mx.sub), fill = Ca,alpha=4/K) +
    geom_ribbon(data=ALMOXARIFE[[z]][[6]][[3]][ALMOXARIFE[[z]][[6]][[3]]$FT <=mins,],aes(x=FT,ymin = mean.mx.sub - se.mx.sub, ymax = mean.mx.sub + se.mx.sub), fill = Ca,alpha=5/K) +
    
    geom_line(data=ALMOXARIFE[[v]][[6]][[1]][ALMOXARIFE[[v]][[6]][[1]]$FT <=mins,],aes(x=FT, y=mean.mx.sub), linetype=1, lwd=W*n, color=M)+
    geom_line(data=ALMOXARIFE[[v]][[6]][[2]][ALMOXARIFE[[v]][[6]][[2]]$FT <=mins,],aes(x=FT, y=mean.mx.sub), linetype=1, lwd=W*n, color=M)+
    geom_line(data=ALMOXARIFE[[v]][[6]][[3]][ALMOXARIFE[[v]][[6]][[3]]$FT <=mins,],aes(x=FT, y=mean.mx.sub), linetype=1, lwd=W*n, color=M)+
    
    geom_line(data=ALMOXARIFE[[z]][[6]][[1]][ALMOXARIFE[[z]][[6]][[1]]$FT <=mins,],aes(x=FT, y=mean.mx.sub), linetype=1, lwd=W*n, color=M)+
    geom_line(data=ALMOXARIFE[[z]][[6]][[2]][ALMOXARIFE[[z]][[6]][[2]]$FT <=mins,],aes(x=FT, y=mean.mx.sub), linetype=1, lwd=W*n, color=M)+
    geom_line(data=ALMOXARIFE[[z]][[6]][[3]][ALMOXARIFE[[z]][[6]][[3]]$FT <=mins,],aes(x=FT, y=mean.mx.sub), linetype=1, lwd=W*n, color=M)+
    
    teruego+ 
    scale_x_continuous(name="",
                       breaks=seq(0,720,120),
                       labels=paste0(seq(0,720,120),rep("min",720/120)),
                       limits=c(0,720),
                       position="top")+
    scale_y_continuous(name="[m]",labels=fancy_scientific)+
    coord_cartesian(ylim = c(0, max(leningrad[[3]][leningrad[[3]]$FT <= mins,]$mean.mx.sub))) 
  

  #### Ok, so I am gonna take kmtv, kmtv.sub, dep.10, dep.sub, hvt, hvt.sub
  X=4.5
  Y=4
  p<-3/5#0.37
  q<-1/3

  if (gg==1){
    ggsave("BO_kmvt_BH50.pdf", plot = KMVT, device="pdf",width = X, height = Y, units = c("in"),dpi = 300)
    ggsave("BO_dep10_BH50.pdf", plot = DEP.10, device="pdf",width = X, height = Y, units = c("in"),dpi = 300)
    ggsave("BO_hvt_BH50.pdf", plot = HVT, device="pdf",width = X, height = Y, units = c("in"),dpi = 300)
    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    ggsave("BO_kmvt.sub_BH50.pdf", plot = KMVT.SUB, device="pdf",width = p*X, height = q*Y, units = c("in"),dpi = 300)
    ggsave("BO_dep10.sub_BH50.pdf", plot = DEP.SUB, device="pdf",width = p*X, height = q*Y, units = c("in"),dpi = 300)
    ggsave("BO_hvt.sub_BH50.pdf", plot = HVT.SUB, device="pdf",width = p*X, height = q*Y, units = c("in"),dpi = 300)
  } else if(gg==2){
    ggsave("BO_kmvt_KC50.pdf", plot = KMVT, device="pdf",width = X, height = Y, units = c("in"),dpi = 300)
    ggsave("BO_dep10_KC50.pdf", plot = DEP.10, device="pdf",width = X, height = Y, units = c("in"),dpi = 300)
    ggsave("BO_hvt_KC50.pdf", plot = HVT, device="pdf",width = X, height = Y, units = c("in"),dpi = 300)
    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    ggsave("BO_kmvt.sub_KC50.pdf", plot = KMVT.SUB, device="pdf",width = p*X, height = q*Y, units = c("in"),dpi = 300)
    ggsave("BO_dep10.sub_KC50.pdf", plot = DEP.SUB, device="pdf",width = p*X, height = q*Y, units = c("in"),dpi = 300)
    ggsave("BO_hvt.sub_KC50.pdf", plot = HVT.SUB, device="pdf",width = p*X, height = q*Y, units = c("in"),dpi = 300)
  } 
  
}















