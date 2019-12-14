rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")


library(metaMS)
library(XML)
library(xcms)
library(ggplot2)
library(scales)
library(ggthemes)

#Function below taken from the xcms (version 2.1) wrapper on Galaxy (workflow4metabolomics), it parses the mzML file and grabs the BPC values.
getBPC <- function(file,rtcor=NULL, ...) {
  object <- xcms::xcmsRaw(file)
  sel <- xcms::profRange(object, ...)
  cbind(if (is.null(rtcor)) object@scantime[sel$scanidx] else rtcor ,xcms:::colMax(object@env$profile[sel$massidx,sel$scanidx,drop=FALSE]))
}

#cleans up the scientific notation on the y-axis..a little buggy in R-studio but may be due to my version of MacOS. 

scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))
}

#edit the files below with the locations of your .mzML files
#code below creates a data frame from the output of getBPC, adds the species labels, adds time

AminoAcid.MS<-getBPC("data/mzML_Files_forR/LEN_MAR1016_05.mzML")
AAnameList<-rep("AA~~Standard",length(AminoAcid.MS))
AminoAcid.MS<-cbind(AminoAcid.MS,AAnameList)
colnames(AminoAcid.MS)<-c("time","intensity","sample")
AminoAcid.MSdata<-as.data.frame(AminoAcid.MS)

Pseudomonas.MS<-getBPC("data/mzML_Files_forR/LEN_FEB2316_12.mzML")
PsnameList<-rep("atop(italic('Pseudomonas'),'sp. KBS0802')",length(Pseudomonas.MS))
Pseudomonas.MS<-cbind(Pseudomonas.MS,PsnameList)
colnames(Pseudomonas.MS)<-c("time","intensity","sample")
Pseudomonas.MSdata<-as.data.frame(Pseudomonas.MS)

Pedobacter.MS<-getBPC("data/mzML_Files_forR/LEN_FEB2316_13.mzML")
PenameList<-rep("atop(italic('Pedobacter'),'sp. KBS0701')",length(Pedobacter.MS))
Pedobacter.MS<-cbind(Pedobacter.MS,PenameList)
colnames(Pedobacter.MS)<-c("time","intensity","sample")
Pedobacter.MSdata<-as.data.frame(Pedobacter.MS)

Janthinobacterium.MS<-getBPC("data/mzML_Files_forR/LEN_FEB2316_14.mzML")
JanameList<-rep("atop(italic('Janthinobacterium'),'sp. KBS0711')",length(Janthinobacterium.MS))
Janthinobacterium.MS<-cbind(Janthinobacterium.MS,JanameList)
colnames(Janthinobacterium.MS)<-c("time","intensity","sample")
Janthinobacterium.MSdata<-as.data.frame(Janthinobacterium.MS)

Variovorax.MS<-getBPC("data/mzML_Files_forR/LEN_FEB2316_16.mzML")
VanameList<-rep("atop(italic('Variovorax'),'sp. KBS0712')",length(Variovorax.MS))
Variovorax.MS<-cbind(Variovorax.MS,VanameList)
colnames(Variovorax.MS)<-c("time","intensity","sample")
Variovorax.MSdata<-as.data.frame(Variovorax.MS)

Burkholdia.MS<-getBPC("data/mzML_Files_forR/LEN_FEB2316_17.mzML")
BunameList<-rep("atop(italic('Burkholdia'),'sp. KBS0801')",length(Burkholdia.MS))
Burkholdia.MS<-cbind(Burkholdia.MS,BunameList)
colnames(Burkholdia.MS)<-c("time","intensity","sample")
Burkholdia.MSdata<-as.data.frame(Burkholdia.MS)

#combine all data sets and change values from factor to numeric or character.
MassSpec<-rbind(AminoAcid.MSdata,Pseudomonas.MSdata,Pedobacter.MSdata,Janthinobacterium.MSdata,Variovorax.MS,Burkholdia.MSdata)
MassSpec$intensity<-as.numeric(as.character(MassSpec$intensity))
MassSpec$time<-as.numeric(as.character(MassSpec$time))
MassSpec$sample<-as.character(MassSpec$sample)


#if ggplot returns an error: error in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : polygon edge not found
#then remove the label=scientific from scale_y_continuous and replot. You will get the same graph but with the regular r encoded scientific notation

metab.fig <- ggplot()+
geom_line(data=MassSpec,aes(x=time,y=intensity, color=sample))+
scale_color_manual(values=c("#999999",  "#BBDF27FF", "#482576FF", "#481467FF", "#463480FF","#2F6C8EFF"))+theme_bw()+theme(legend.position="none")+
xlab("Time (seconds)")+ylab("Base Peak Chromatogram")+scale_x_continuous(limits=c(360,1600))+scale_y_continuous(label=scientific)+
facet_grid(rows=vars(sample), labeller = label_parsed)

ggsave(file="figs/metabolomics.png", metab.fig, width=12,height=9, units='in', dpi=600)


