#############################################
#############################################
# load data
load("./a18.RData")
library(magrittr)
library(dplyr)
library(ggplot2)
a18 <- a18 %>% 
  mutate(sido=substr(reg, 1, 2)) %>%
  filter(sido=="서울") %>%
  mutate(gu=gsub(pattern="서울 ",  replace="", reg))

##########################
# selection of data
xg <- a18
# hour의 factor화
xg <- xg %>% mutate(fhour=as.factor(hour))

# boxplot
xg %>% 
  filter(gu=="강남구") %>% 
  ggplot(mapping=aes(x=as.factor(hour), y=so2)) + 
  geom_boxplot() #+ facet_grid(~gu)
# mean plot according to hour
#filter(so2 < quantile(so2, prob=0.9)) %>%
tmp <- xg %>% 
  filter(gu=="강남구") %>% 
  group_by(hour) %>%
  mutate(mso2=mean(so2)) %>%
  ggplot(mapping=aes(x=hour, y=mso2)) + 
  geom_line(linetype='dotted', cex=3) 

xg %>% 
  filter(gu=="강남구") %>% 
  group_by(hour) %>%
  filter(o3 < quantile(o3, prob=0.9)) %>%
  mutate(mo3=mean(o3)) %>%
  ggplot(mapping=aes(x=hour, y=mo3)) + 
  geom_line(aes(linetype='dotted'))

xg %>% 
  filter(gu=="강남구") %>% 
  group_by(hour) %>%
  filter(pm10 < quantile(pm10, prob=0.99)) %>%
  mutate(mpm10=mean(pm10)) %>%
  ggplot(mapping=aes(x=hour, y=mpm10)) + 
  geom_line(linetype='dotted') 

xg %>% 
  filter(gu=="강남구") %>% 
  group_by(hour) %>%
  filter(pm25 < quantile(pm25, prob=0.99)) %>%
  mutate(mpm25=mean(pm25)) %>%
  ggplot(mapping=aes(x=hour, y=mpm25)) + 
  geom_line(linetype='dotted') 

xg %>% 
  ggplot(mapping=aes(x=gu, y=so2)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5))   

xg %>% 
  ggplot(mapping=aes(x=forcats::fct_reorder(gu, so2, .fun=quantile, 0.5, .desc=TRUE), y=so2)) + 
  geom_boxplot() + 
  theme(axis.text.x=element_text(angle=90)) + xlab("")

####
mxg <- xg %>% 
  group_by(hour) %>% 
  summarize(so2=mean(so2), 
            no2=mean(no2), co=mean(co), o3=mean(o3),
            pm10=mean(pm10), pm25=mean(pm25))

mxg <- xg %>% 
  group_by(hour) %>% 
  select(hour, so2:pm25) %>%
  summarise_all(funs(m="mean")) %>%
  rename(so2=so2_m, no2=no2_m, co=co_m, o3=o3_m, pm10=pm10_m, pm25=pm25_m)

gg <- mxg %>% ggplot() 
gg1 <- gg + theme(legend.position = "right") +
  geom_line(mapping=aes(x=hour, y=so2*1e4/2), linetype="dotted") +
  geom_line(mapping=aes(x=hour, y=no2*2e3/2), colour="green") +
  geom_line(mapping=aes(x=hour, y=co*1e2/2), colour="blue") +
  geom_line(mapping=aes(x=hour, y=pm10/2), colour="black") +
  geom_line(mapping=aes(x=hour, y=pm25), colour="red") +
  geom_vline(xintercept=c(9,12), color="black", linetype=2, cex=0.5)
gg1

# long type: absolute scale & simultaneously
long_mxg <- mxg %>% 
  tidyr::gather(so2:pm25, key="source", value="quan")  

long_mxg %>% 
  ggplot(mapping=aes(x=as.factor(hour), y=quan, group=source, color=source)) +
  geom_line()

# long type: relative scale
#long_mxg <- mxg %>% 
#  mutate(so2=so2*1e4/2, no2=no2*2e3/2, o3=o3*2e3, co=co*1e2/2, 
#         pm10=pm10/2, pm25=pm25) %>% 
#  tidyr::gather(so2:pm25, key="source", value="quan")  
long_mxg <- mxg %>% 
  mutate(so2=scale(so2), no2=scale(no2), o3=scale(o3), co=scale(co), 
         pm10=scale(pm10), pm25=scale(pm25)) %>% 
  tidyr::gather(so2:pm25, key="source", value="quan")  

#heatmap of 6 sources
heatmap(cor(mxg[,-1]), scale="none")

# along hour
gg1 <- long_mxg %>% 
  ggplot(mapping=aes(x=hour, y=quan, group=source, color=source)) +
  geom_line(cex=1.1,aes(linetype=source)) + 
  ylab("") + 
  theme(axis.text.y=element_blank()) + xlim(1,24) +
  labs(title="source", subtitle="25 regions in Seoul")
gg1

# theme 위치
## legend.position = "none", "top", "bottom", "left"
gg1 + theme(legend.position = "right", legend.text=element_text(size=15, colour=1:6))

# along hour(only pm10, pm25)
gg2 <- long_mxg %>% filter(source %in% c("pm10", "pm25")) %>%
  ggplot(mapping=aes(x=hour, y=quan, group=source, color=source)) +
  geom_line(cex=1.1,aes(linetype=source)) + 
  ylab("") + 
  theme(axis.text.y=element_blank()) 
gg2
#heatmap of 2 sources
heatmap(cor(mxg[,-1]), scale="none")
library(GGally)
par(las=1)
GGally::ggpairs(mxg,diag=list(continuous="density", alpha=0.5),axisLabels="show")

####

rm(list=ls())
library(Matrix)
library(Rcpp)
library(RcppArmadillo)
library(devtools)
library(Rglpk)
library(MethylCapSig)
install_github("glmgen/genlasso")
library(genlasso)
#setwd("~/github/ComLasso")
sourceCpp('./src/inner.cpp')
source("./src/ComLassoC.R")
# income
rX<- read.csv("./data/capital_season.csv")
ry<- read.csv("./data/income.csv")[,-1]
y = unlist(ry[5,]/ry[1,])
y = y[-length(y)] # 2003-1~2019-1
X = as.matrix(rX[,130:197]) # 2002-1~2018-4
X = t(X)
x = sweep(X,1,rowSums(X),"/")
idx = 1:nrow(x)
y
x1 = x[tail(idx, 65),]
x2 = x[tail(idx, 65)-1,]
x3 = x[tail(idx, 65)-2,]
x4 = x[tail(idx, 65)-3,]
xx = cbind(x1,x2,x3,x4)
xx = log(xx)
pk = rep(8,4)
fit = comLassoC(xx,y,pk=pk,lam_min=0,tol=1e-08,KKT_check=FALSE)
dim(fit$coefMat)
coefMat = fit$coefMat
est.var = sum((y-cbind(1,xx)%*%coefMat[nrow(coefMat),])^2)/(length(y)-1)

aic.vec = c()
i = 1
for (i in 1:nrow(coefMat))
{
  coefvec = coefMat[i,-1]
  df = 0
  for (j in 1:4)
  {
    subvec = coefvec[((j-1)*8+1):(j*8)]
    if (sum(subvec !=0)>0) df <- df + (sum(subvec !=0))-1
  }
  aic.vec[i]<-
    sum((y-cbind(1,xx)%*%coefMat[i,])^2)/est.var + log(length(y)) + 2*df
}
plot(aic.vec)
which.min(aic.vec)
library(ggplot2)
data.frame(ID=1:length(aic.vec), aic=aic.vec) %>% 
  ggplot(aes(x=ID,y=aic))+ geom_line() + 
  xlab("") + 
  geom_vline(aes(xintercept=17, color="red"), linetype=2) + 
  theme(legend.position = "none")

## 17
par(mfrow=c(2,2))
for(k in 1:4){
  #k = 1
  for(i in 1:8)
  {
    if (i == 1)
      plot(fit$coefMat[,i+8*(k-1)+1], ylim = c(-10,10), col = i, 
           type = 'l', ylab="coefficients", xlab="")
    if (i>1)
      lines(fit$coefMat[,i+8*(k-1)+1], ylim = c(-10,10),
            col = i)
    abline(v=17, lty=2, lwd=0.6)
  }
  mtext(paste0(k,"-lagged year"), side=3)
}
fit$coefMat[17,-1]
library(xtable)
mat <- matrix(fit$coefMat[17,-1],byrow=T,ncol=8)
xtable(round(mat,3),4)
mat

#for
rX<- read.csv("./data/capital_season.csv", stringsAsFactors=FALSE, header=TRUE)
ry<- read.csv("./data/income.csv")[,-1]
y = unlist(ry[5,]/ry[1,])
y = y[-length(y)] # 2003-1~2019-1
X = as.matrix(rX[,130:197]) # 2002-1~2018-4
X = t(X)
x = sweep(X,1,rowSums(X),"/")
idx = 1:nrow(x)
y
x1 = x[tail(idx, 65),]

library(ggplot2)
mx1 <- reshape2::melt(x1)
names(mx1)
library(magrittr)
mx1 %<>% 
  #  dplyr::mutate(Var1=as.character(Var1)) %>%
  dplyr::rename(time=Var1, asset=Var2, prop=value) %>%
  dplyr::mutate(asset=factor(asset, levels=1:8, labels=c("RS","NS","IF","TE","ME","IP","R&D","Others")))

g1 <- data.frame(mx1) %>% 
  ggplot(mapping=aes(x=time, y=prop, group=asset)) + 
  geom_area(aes(fill=asset)) + xlab("")
g1

dd <- c()
for(i in 2003:2018){
  if(i %% 2==1){
    #dd <- c(dd, c(paste0(i, " - Q1"), "-Q2", "-Q3", "-Q4"))
    dd <- c(dd, c(paste0(i, " - 1Q lagged"), "", "", ""))
  }
  #else{
  #  dd <- c(dd, c(paste0(i, " - 1Q lagged"), "", "", ""))
  #}
}
dd <- c(dd, "2019 - 1Q lagged")

#g2 <- g1 +  scale_x_discrete(breaks=seq(1,65, by=2), labels=dd) + 
#  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust =1, size=7)) +
#  labs(x="", y="Proportion") 
xlab_title <- levels(mx1[,1])
xlab_title[c(1,9,17,25,33,41,49,57,65)] <- 
  sapply(strsplit(xlab_title[c(1,9,17,25,33,41,49,57,65)], "X"), function(x){x[2]})
xlab_title[-c(1,9,17,25,33,41,49,57,65)] <- ""
mx1$time <- as.character(mx1$time)
mx1$time[1] <- "time"
g2 <- g1 +  scale_x_discrete(mx1$time, labels=xlab_title) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust =1, size=15)) +
  theme(axis.text.y = element_text(size=15)) +
  theme(legend.text = element_text(size=20), 
        legend.title=element_text(size=20)) +
  theme(legend.position = "top") + 
  xlab("") + 
  labs(y="Proportion")
g2

ey <- rep(y, times=8)
ma <- 7.920426
mi <- 5.191999
g3 <- g2 + geom_line(aes(x=time, y=(ey-mi)/(ma-mi), group=asset))
g3

g3 + scale_y_continuous(
  sec.axis = sec_axis(~ ./1.5 + 0.2, name = "income inequality", 
                      labels = function(b) { round((ma-mi)*(b)+mi,2) }
  )
) + theme(axis.title.y = element_text(size=15)) 




