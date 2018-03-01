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
df1 <- df[df$target == 'target' & df$binder != 'no' & df$qvalue < 0.05, ] %.% group_by(qvalue) %.% summarise(n = n()) %.% mutate(n = cumsum(n))
df1$netMHC = "binder"
df2 <- df[df$target == 'target' & df$binder == 'no' & df$qvalue < 0.05, ] %.% group_by(qvalue) %.% summarise(n = n()) %.% mutate(n = cumsum(n))
df2$netMHC = "non-binder"
dfu <- df[df$target == 'target' & df$binder != 'no' & df$qvalue < 0.05, ]
dfu <- dfu[order(dfu$qvalue),] 
dfu <- distinct(dfu, sequence)
df3 <- dfu[dfu$target == 'target' & dfu$binder != 'no' & dfu$qvalue < 0.05, ]  %.% group_by(qvalue) %.% summarise(n = n()) %.% mutate(n = cumsum(n))
df3$netMHC = "unique binder"
dfn <- rbind(df1,df2)
dfn <- rbind(dfn,df3)
ggplot(data=dfn,  aes(x=qvalue, y=n, color=netMHC)) + 
  ylab("target identifications") + xlab("FDR") + scale_x_reverse() +
  geom_line() + ggtitle("(Non-)Binder count over FDR")
######################################
garbage<-dev.off()
file.rename(randfile, post) #but we have to give the requested out file
