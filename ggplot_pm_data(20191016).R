cal_aic <- function(y, x, obj){
  bm <- rbind(obj$beta0.rec, obj$beta.rec[[1]])
  #n by step  
  n <- nrow(x)
  nstep <- ncol(bm)
  dfv = colSums(bm!=0)
  pred_m <- cbind(1,x) %*% bm
  loss <- apply(pred_m, 2, function(xx){sum((xx-y)^2)})
  mse <- loss[nstep]/(n-dfv[nstep])
  rloss <- loss/2/n/mse + 2*(dfv)/n
  return(rloss)  
}

matplot.comlasso.rev <- function(lam, coef1, cnums1, cnums2, cnums3, breaks = FALSE, ...) 
{
  S <- sum(abs(coef1[nrow(coef1), ]))
  s1 <- rowSums(abs(coef1))/S
  col_names <- colnames(coef1)
  matplot(s1, coef1, ..., type = "b", pch = "*")
  abline(h = 0, lty = 3)
  axis(4, at = coef1[nrow(coef1), cnums1], 
       labels = col_names[cnums1], cex = 0.8, cex.lab = 0.3, 
       line = 0)
  axis(4, at = coef1[nrow(coef1), cnums2], 
       labels = col_names[cnums2], cex = 0.8, 
       cex.lab = 0.3, tick = FALSE, line = 1)
  axis(4, at = coef1[nrow(coef1), cnums3], 
       labels = col_names[cnums3], cex = 0.8, 
       cex.lab = 0.3, tick = FALSE, line = 0.5)
  stepid = trunc(as.numeric(dimnames(coef1)[[1]]))
  if (breaks) {
    axis(3, at = s1, labels = paste(stepid), cex = 0.8)
    abline(v = s1)
  }
}
#############################################
#############################################
#############################################
#############################################
# load data
load("./a18_comp.RData");
library(comlasso)
library(magrittr)
library(dplyr)
library(ggplot2)
a18_comp <- a18_comp %>% 
  mutate(sido=substr(reg, 1, 2)) %>%
  filter(sido=="서울") %>%
  mutate(gu=gsub(pattern="서울 ",  replace="", reg))
a18_comp[which(a18_comp==0, arr.ind=TRUE),4] <- 1e-8
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


