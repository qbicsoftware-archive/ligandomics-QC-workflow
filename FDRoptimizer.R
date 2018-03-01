## This is an R script to produce the figures that are attached to the qcML format
library("ggplot2")
library(dplyr)
library(scales)
library("numDeriv")

options(warn=-1) #suppress warnings

#options
options(digits=10)

file<-commandArgs(TRUE)[1]
post<-commandArgs(TRUE)[2]
out<-commandArgs(TRUE)[3]
#file<-"~/Desktop/BDID/bladder-zh02-example/150209_DK_BD-ZH02_Bladder_W_20perc_Rep#1_50cm195min3s_msms4.seperformance-xtandem.csv"
#out<-"/tmp/test.csv"
knime.in<-read.csv(file=file,head=TRUE,sep="\t")

randfile <- file.path(dirname(post), paste(toString(sample(111111:999999, 1)), ".png"))
png(randfile) #averts nasty file names R png cannot cope with
######################################
###Binder cut on unique identifications over FDR
######################################
df<-knime.in
dfu <- df[df$target == 'target', ]
dfu <- dfu[order(dfu$qvalue),] 
dfu <- distinct(dfu, sequence)

df1 <- dfu[dfu$binder != 'no', ] %.% group_by(qvalue) %.% summarise(b = n()) %.% mutate(b = cumsum(b))
df2 <- dfu %.% group_by(qvalue) %.% summarise(n = n()) %.% mutate(n = cumsum(n))
df3 <- df[df$target == 'target', ] %.% group_by(qvalue) %.% summarise(a = n()) %.% mutate(a = cumsum(a))
dfn <- merge(x=df2, y=df1, by.x="qvalue", by.y="qvalue")
dfn <- merge(x=dfn, y=df3, by.x="qvalue", by.y="qvalue")
dfn$bp <- 100/dfn$n * dfn$b 
dfn$ap <- 100/dfn$a * dfn$b 

#ggplot(data=dfn) + 
#  ylab("binder percentage of all target identifications") + 
#  xlab("FDR") + 
#  #scale_x_reverse() +
#  geom_line(aes(x=qvalue, y=bp, color="Percentage of unique binders from the\nunique identifications at given thresholds")) + 
#  geom_line(aes(x=qvalue, y=ap, color="Percentage of unique binders from all\ntarget identifications at given thresholds")) + 
#  ggtitle("Binder percentage of all target identifications over FDR") + 
#  geom_line(aes(x=qvalue, y=bp, color="smoothed uniques"), stat="smooth", method = "loess", formula=y~x, span=0.05, alpha=0.3) + 
#  scale_x_continuous(breaks=pretty_breaks()) +
#  scale_colour_manual(name="Line\ncolor", values=c("black","blue","red")) + 
#  theme(legend.position="bottom")


###find a good fdr part
loe <- loess(bp~qvalue,dfn, span=0.05)
#ggplot(dfn, aes(x=qvalue,y=bp), geom='smooth', span =0.5)

#When the fit was made using surface = "interpolate" (the default), predict.loess will not extrapolate – so points outside an axis-aligned hypercube enclosing the original data will have missing (NA) predictions and standard errors
loe <- loess(bp~qvalue,dfn, span=0.05, control=loess.control(surface="direct"))
#Another disadvantage of LOESS is the fact that it does not produce a regression function that is easily represented by a mathematical formula. This can make it difficult to transfer the results of an analysis to other people. In order to transfer the regression function to another person, they would need the data set and software for LOESS calculations.


func0 <- function(x){ predict(loe, x) }
dfn$grad <-grad(func0, dfn$qvalue)

#ggplot(dfn) + geom_line(aes(x=qvalue, y=grad, color="gradient")) + 
#  geom_line(aes(x=qvalue, y=bp, color="smoothed uniques"), stat="smooth", method = "loess", formula=y~x, span=0.05, alpha=0.3) +
#  geom_line(aes(x=qvalue, y=bp, color="Percentage of unique binders from the\nunique identifications at given thresholds")) +
#  xlim(c(0,.1))

#min(dfn$grad)
mx <- dfn[which.min(dfn$grad), ]$qvalue
my <- dfn[which.min(dfn$grad), ]$bp

tresh <- FALSE
if (mx > 0.05){ 
  mx <- 0.05
  my <- dfn[min(which(dfn$qvalue > 0.05)), ]$bp
  tresh <- TRUE
}

ggplot(data=dfn) + 
  ylab("binder percentage of all target identifications") + 
  xlab("FDR") + 
  #scale_x_reverse() +
  geom_line(aes(x=qvalue, y=bp, color="Percentage of unique binders from the\nunique identifications at given thresholds")) + 
  geom_line(aes(x=qvalue, y=ap, color="Percentage of unique binders from all\ntarget identifications at given thresholds")) + 
  ggtitle("Binder percentage of all target identifications over FDR") + 
  #geom_line(aes(x=qvalue, y=bp, color="smoothed uniques"), stat="smooth", method = "loess", formula=y~x, span=0.05, alpha=0.3)+
  geom_point(size = 3, aes(x=mx, y=my, colour = "lowest relative binder increase")) + 
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_colour_manual(name="Line\ncolor", values=c("black","blue","red")) + 
  { if(tresh) geom_vline(xintercept = 0.05, colour="darkred", linetype = "longdash") else geom_vline(xintercept = mx, colour="green", linetype = "longdash")} +
  annotate("text", x=mx+0.015, y=0, label=round(mx, digits = 3), color = "black", size=4) +
  theme(legend.position="bottom") 

ofo <- data.frame("optimal fdr"=numeric(), "binder percentage"=numeric(), "file"=character()) 
ofo[1,] <- c(mx,my,basename(file))
ofo$file <- basename(file)
write.csv(ofo, file = out, row.names=FALSE, quote=FALSE)

####old percent plot - still there is information in it
#reversed fdr scale, indicates the increase of binders or unique binder in percentage to all idents
#ab <- nrow(df[df$target == 'target' & df$binder != 'no', ])
#df_binder <- df[df$target == 'target' & df$binder != 'no' & df$qvalue < 0.05, ] %.% group_by(qvalue) %.% summarise(n = n()) %.% mutate(n = 100/ab * cumsum(n))
#df_binder$sequences = "binder"
#df_uniquebinder <- dfu[dfu$target == 'target' & dfu$binder != 'no' & dfu$qvalue < 0.05, ] %.% group_by(qvalue) %.% summarise(n = n()) %.% mutate(n = 100/ab * cumsum(n))
#df_uniquebinder$sequences = "uniques"
#df_sum <- rbind(df_binder,df_uniquebinder)

#ggplot(data=df_sum,  aes(x=qvalue, y=n, color=sequences)) +
#  ylab("percentage of all target binder identifications") + xlab("FDR") + scale_x_reverse() +
#  geom_line() +
#  geom_vline(xintercept = mx, colour="green", linetype = "longdash") +
#  ggtitle("Percentage of binder over FDR")

######################################
garbage<-dev.off()
file.rename(randfile, post) #but we have to give the requested out file
