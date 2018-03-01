## This is an R script to produce the figures that are attached to the qcML format
library("ggplot2")
library(dplyr)

options(warn=-1) #suppress warnings

#options
options(digits=10)

file<-commandArgs(TRUE)[1]
post<-commandArgs(TRUE)[2]
#file<-"~/Desktop/BDID/bladder-zh02-example/150209_DK_BD-ZH02_Bladder_W_20perc_Rep#1_50cm195min3s_msms4.seperformance-xtandem.csv"
knime.in<-read.csv(file=file,head=TRUE,sep="\t")

randfile <- file.path(dirname(post), paste(toString(sample(111111:999999, 1)), ".png"))
png(randfile) #averts nasty file names R png cannot cope with
######################################
###Binder counts
######################################
df<-knime.in
levels(df$allele)[levels(df$allele)==""] <- "nonbinder"
ggplot(data=df[df$target == 'target' & df$qvalue < 0.05, ],  aes(allele)) + 
  geom_histogram() +
  xlab("Allele") + 
  ggtitle("Identified sequences distribution over sample typing")
######################################
garbage<-dev.off()
file.rename(randfile, post) #but we have to give the requested out file
