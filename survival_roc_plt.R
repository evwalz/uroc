#------------------------------------------------------------------------------#
# This code was used to compute the survival example in 

# Tilmann Gneiting and Eva-Maria Walz. 
# Receiver operating characteristic (ROC) movies, universal ROC (UROC) curves, and coefficient of predictive ability (CPA)
# URL https://arxiv.org/abs/1912.01956

# Histogram of two biochemical markers for survival and non survival with corresponding ROC curve
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

# thresholding survival time
response_bin = as.numeric(response>1462)
auc_alb = round(cpa(response_bin, alb),2)
auc_bili = round(cpa(response_bin, bili),2)

# define two thresholds
thres_alb1 = 3.815
thres_alb2 = 3.125

thres_bili1 = 21.35
thres_bili2 = 26.75

# Histogram
df_alb = data.frame(x = alb, survival = as.factor(response_bin==0))
df_bili = data.frame(x = bili, survival = as.factor(response_bin==0))

hist_alb = ggplot(data =df_alb, aes(x=x,color=survival,fill=survival))+geom_histogram(aes(y = ..count..),position="identity",binwidth=0.1, alpha = 0.4) +
  theme_classic(base_size = 20)+
  scale_fill_manual(values = c("#1D91C0","#081D58"), labels=c("Survival","Non-Survival"))+
  scale_color_manual(values = c("#1D91C0","#081D58"),labels=c("Survival","Non-Survival"))+
  theme(legend.title = element_blank())+
  xlab(label="Serum Albumin (mg/dl)")+
  ylab(label="Count")+
  labs(caption = "its a test")+
  geom_segment(aes(x=thres_alb1,y=0,xend=thres_alb1,yend=11), size=1, lty=5,color="#DD3497")+
  geom_segment(aes(x=thres_alb2,y=0,xend=thres_alb2,yend=11), size=1, lty=5,color="#FD8D3C")+
  theme(legend.position = c(0.18,0.8),
        axis.line = element_line(size=0.5),
        axis.ticks = element_line(size = 0.5),
        plot.caption = element_text(size = 55,color = "white", face = "italic"))


hist_bili = ggplot(data =df_bili, aes(x=x,color=survival,fill=survival))+geom_histogram(aes(y = ..count..),position="identity",binwidth=1, alpha = 0.4) +
  theme_classic(base_size = 20)+
  scale_fill_manual(values = c("#74C476","#01665E"), labels=c("Survival","Non-Survival"))+
  scale_color_manual(values = c("#74C476","#01665E"),labels=c( "Survival","Non-Survival"))+
  scale_x_continuous(breaks = seq(8,28,by=10), labels = c(20, 10, 0))+
  scale_y_continuous(breaks = seq(0, 25, by = 5))+
  theme(legend.title = element_blank())+
  xlab(label="Serum Bilirubin (mg/dl)")+
  ylab(label="Count")+
  labs(caption = "its a test")+
  theme(legend.position = c(0.18,0.8),
        axis.line = element_line(size=0.5),
        axis.ticks = element_line(size = 0.5),
        plot.caption = element_text(size = 55,color = "white", face = "italic"))+
  geom_segment(aes(x=thres_bili1,y=0,xend=thres_bili1,yend=22), size=1, lty=5,color="#DD3497")+
  geom_segment(aes(x=thres_bili2,y=0,xend=thres_bili2,yend=22), size=1, lty=5,color="#FD8D3C")

# ROC curve
roc_alb = uroc(response_bin, alb)
roc_bili = uroc(response_bin, bili)

response_order <- order(response_bin, decreasing=FALSE)
response <- response[response_order]
thresholds <- unique(response)[-1]

idx_alb1 = which(thresholds==thres_alb1)
idx_alb2 = which(thresholds==thres_alb2)

points_alb1=data.frame(x = (roc_alb$farate[idx_alb1]) ,y=roc_alb$hitrate[idx_alb1])
points_alb2=data.frame(x = (roc_alb$farate[idx_alb2]) ,y=roc_alb$hitrate[idx_alb2])

idx_bili1 = which(thresholds==thres_bili1)
idx_bili2 = which(thresholds==thres_bili2)

points_bili1=data.frame(x = (roc_bili$farate[idx_bili1]) ,y=roc_bili$hitrate[idx_bili1])
points_bili2=data.frame(x = (roc_bili$farate[idx_bili2]) ,y=roc_bili$hitrate[idx_bili2])


df_roc = data.frame(far =c(roc_alb$farate, roc_bili$farate) , hit= c(roc_alb$hitrate, roc_bili$hitrate), Marker = c(rep("Albumin",length(roc_alb$farate)), rep("Bilirubin",length(roc_bili$farate))))

roc_survival = ggplot()+geom_line(data = df_roc, aes(x=far,y=hit, group=Marker, col=Marker), lwd=1)+
  theme_minimal(base_size = 20)+
  geom_line(data=data.frame(x=c(0,1),y=c(0,1)),aes(x=x,y=y),col="grey40",lty=2)+
  xlab(label="1 - Specificity")+
  ylab(label = "Sensitivity")+
  scale_x_continuous(breaks = seq(0, 1, by = 1/6) , labels = c("0", "1/6", "1/3", "1/2", "2/3", "5/6", "1"))+
  scale_y_continuous(breaks = seq(0, 1, by = 1/6), labels = c("0", "1/6", "1/3", "1/2", "2/3", "5/6", "1"))+
  scale_color_manual(values = c("#225EA8", "#00BA38"))+
  annotate("text",x=0.79,y=0.135,label=paste("AUC:",format(auc_alb,nsmall = 2)),size=6.5, col="#225EA8")+
  annotate("text",x=0.79,y=0.20,label=paste("AUC:",format(auc_bili,nsmall = 2)),size=6.5, col = "#00BA38")+
  geom_point(data = points_alb1,aes(x=x,y=y),col="#DD3497", size=3, shape = 4, stroke = 2)+
  geom_point(data = points_alb2,aes(x=x,y=y),col="#FD8D3C",size=3, shape = 4, stroke = 2)+
  geom_point(data = points_bili1,aes(x=x,y=y),col="#DD3497", size=3, shape = 4, stroke = 2)+
  geom_point(data = points_bili2,aes(x=x,y=y),col="#FD8D3C",size=3, shape = 4, stroke = 2)+
  theme(plot.margin = unit(c(1,1,0.3,1.5),"cm"),
        legend.position = "bottom",
        legend.title=element_text(size=22), 
        legend.text=element_text(size=22))


plt_alb =  ggarrange(hist_alb, labels=c("A"),ncol=1,nrow=1, font.label = list(size=22))
plt_roc = ggarrange(roc_survival ,labels=c("B"),ncol=1,nrow=1, font.label = list(size=22))
plt_bili = ggarrange(hist_bili ,labels=c("C"),ncol=1,nrow=1, font.label = list(size=22))

pdf("survival_roc.pdf", width=21, height = 7)
grid.arrange(plt_alb, plt_roc, plt_bili,layout_matrix=rbind(c(1,2,3)))
dev.off()