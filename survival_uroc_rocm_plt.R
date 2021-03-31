#------------------------------------------------------------------------------#
# This code was used to compute the survival example in 

# Tilmann Gneiting and Eva-Maria Walz. 
# Receiver operating characteristic (ROC) movies, universal ROC (UROC) curves, and coefficient of predictive ability (CPA)
# URL https://arxiv.org/abs/1912.01956

# First UROC curve for survival exmaple is computed
# Then two corresponding ROCMs are created
#------------------------------------------------------------------------------#


library(SMPracticals)
library(dplyr)
library(uroc)
library(ggplot2)
library(ggpubr)
library(gridExtra)

data(pbc)

# filter censoring and NAs in t5
pbc_filter <- pbc %>% filter(status==1) 

response <- pbc_filter$time
alb <- pbc_filter$alb
bili <- pbc_filter$bili
bili <- max(bili)-bili

name_alb <- "Albumin"
name_bili <- "Bilirubin"

# uroc plot 
uroc_alb <- uroc(response, alb)
uroc_bili <- uroc(response, bili)

cpa_alb <- round(cpa(response, alb),2)
cpa_bili <- round(cpa(response, bili),2)

uroc_far <- c(uroc_alb$farate,uroc_bili$farate)
uroc_hit <- c(uroc_alb$hitrate,uroc_bili$hitrate)

type = c(rep(name_alb,length(uroc_alb$farate)),rep(name_bili,length(uroc_bili$farate)))

df <- data.frame(Far = uroc_far, Hitrate = uroc_hit, Marker=as.character(type))

df$Marker <- factor(df$Marker, levels = c(name_alb, name_bili))

uroc_survival = ggplot()+geom_line(data=data.frame(x=c(0,1),y=c(0,1)),aes(x=x,y=y),col=grey(0.3),lty=2) + 
  theme_minimal(base_size = 22)+
  xlab(label="1 - Specificity")+
  ylab(label = "Sensitivity")+
  scale_x_continuous(breaks = seq(0, 1, by = 1/6) , labels = c("0", "1/6", "1/3", "1/2", "2/3", "5/6", "1"))+
  scale_y_continuous(breaks = seq(0, 1, by = 1/6), labels = c("0", "1/6", "1/3", "1/2", "2/3", "5/6", "1"))+
  geom_line(data = df, aes(x=Far, y=Hitrate, group = Marker, colour=Marker), lwd=1.0)+
  scale_color_manual(values=c("#225EA8","#00BA38"))+
  annotate("text",x=0.8,y=0.205,label=paste("CPA:",format(cpa_alb,nsmall = 2)),size=7, col="#225EA8")+
  annotate("text",x=0.8,y=0.1275,label=paste("CPA:",format(cpa_bili,nsmall = 2)),size=7, col="#00BA38")+
  annotate("text", x=0.01, y=1.03, label=paste("UROC curves"), size=7, hjust = 0)

print(uroc_survival)

# animation of rocm
rocm(response, alb)
rocm(response, bili)