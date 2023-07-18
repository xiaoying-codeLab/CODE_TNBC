
rm(list = ls())  
options(stringsAsFactors = F) 

library(data.table)
library(dplyr)
library(survminer) 
library(survival)
library(loose.rock)
library(futile.logger) 
library(glmSparseNet)
library(ggrisk)  
library(pheatmap) 
library(ggplot2)
library(ggsci)
library(ggstatsplot)


# 0. load data----
load('04-1-doModel-lasso//out-data/coefs.v_lasso_model.Rdata')
xdata[1:4, 1:4]
head(ydata)
# 1. multivarible-cox----
## 1.1 tidy data for multicox----
dat <- ydata[, -3]
head(dat)  # include survdata information
## z-score
gsub('-','_',colnames(n))
colnames(xdata) = gsub('-','_',colnames(xdata))
colnames(xdata) 

cg <- sort(names(coefs.v))
cg <- names(coefs.v)
cg = gsub('-','_',cg)


n <- apply(xdata,2,scale)[,cg] %>%
  as.data.frame()
head(n)
rownames(n) <- rownames(xdata)
boxplot(n)
dat_cox <- cbind(dat,n)
head(dat_cox)

## prepare variale for Surv(), paste gene names
multivariate <- paste(sort(colnames(n)), collapse = '+') 

## 1.2 multi-cox----
attach(dat_cox)
s <-  paste0(' Surv(time, event) ~  ', multivariate )
model <- coxph(as.formula(s), data = dat_cox )
summary(model, data = dat_cox)
risk_score <- predict(model, type = 'risk', data = dat_cox)
dat_cox$Risk_score <- risk_score 


# 2.  visualization----
## 2.1 forest----
library(survminer)
library(tidyverse)
library(ggplot2)
library(ggstatsplot)
library(survival)
library(stringr)
library(viridis)
library(scales)
options(scipen=1)
dat_forest <- cbind(dat,n)
colnames(dat_forest)
Coxoutput=data.frame()
for(i in colnames(dat_forest[,3:ncol(dat_forest)])){
  cox <- coxph(Surv(time, event) ~ dat_forest[,i], data = dat_forest)
  coxSummary = summary(cox)
  Coxoutput=rbind(Coxoutput,cbind(gene=i,HR=coxSummary$coefficients[,"exp(coef)"],
                                  z=coxSummary$coefficients[,"z"],
                                  pvalue=coxSummary$coefficients[,"Pr(>|z|)"],
                                  lower=coxSummary$conf.int[,3],
                                  upper=coxSummary$conf.int[,4]))
}
for(i in c(2:6)){
  Coxoutput[,i] <- as.numeric(as.vector(Coxoutput[,i]))
}
Coxoutput <- arrange(Coxoutput,pvalue)
write.csv(Coxoutput,file = "04-1-doModel-lasso/out-data/forest.cox.csv")
plotCoxoutput <- filter(Coxoutput,HR <=0.92 | HR>= 1.15)  #
ggplot(data=plotCoxoutput,aes(x=HR,y=gene,color=pvalue))+
  geom_errorbarh(aes(xmax=upper,xmin=lower),color='black',height=0,size=1.0)+
  geom_point(aes(x=HR,y=gene),size=5,shape=18)+   #
  geom_vline(xintercept = 1,linetype='dashed',size=1.2)+
  scale_x_continuous(breaks = c(0.75,1,1.30))+
  coord_trans(x='log2')+ 
  ylab("Gene")+  #
  xlab("Hazard ratios of lncRNA in FUSCC")+ 
  labs(color="P value",title ="" )+
  scale_color_gradient2(low = muted("#DB423E"),mid ="white",high =muted("#008ECA"),midpoint = 0.15)+ #
  theme_bw(base_size = 12)+   #
  theme(panel.grid =element_blank(),  #
        axis.text.x = element_text(face="bold", color="black", size=9),    #
        axis.text.y = element_text(face="bold",  color="black", size=9),
        axis.title.x = element_text(face="bold", color="black", size=11),
        axis.title.y = element_text(face="bold",color="black", size=11),
        legend.text= element_text(face="bold", color="black", size=9),
        legend.title = element_text(face="bold", color="black", size=11),
        panel.border = element_rect(colour = 'black',size=1.4))   #
ggsave('04-1-doModel-lasso//out-plot/multicox_forest_3.pdf', height = 9, width = 8)


## 2.2 nomo----
library(rms)
library(openxlsx)
phe.2 = read.xlsx("data_input/FUSCC/mmc2.xlsx",startRow = 2)
rownames(dat_cox)
rownames(phe.2) = phe.2$Project_ID
table(rownames(dat_cox) %in% rownames(phe.2))
identical(rownames(dat_cox), rownames(phe.2))
phe.2 = phe.2[rownames(dat_cox),]
identical(rownames(dat_cox), rownames(phe.2))
colnames(phe.2)
dat_cox = cbind(dat_cox,phe.2)
colnames(dat_cox)
multivariate <- paste(c("Age_at_surgery","Grade","T","Ki67","Size_cm",
                        "sTILs","iTILs","ERBB2_IHC_score","Risk_score"), collapse = '+') 
s <-  paste0(' Surv(time, event) ~  ', multivariate )
library(rms)
dc <- datadist(dat_cox);dc
options(datadist="dc")
boxplot(dat_cox$time)
# Cox Proportional Hazards Model and Extensions , cph {rms}
f <- cph(as.formula(s),
         x=T, y=T, surv=T, 
         data=dat_cox, time.inc=15)
summary(f) 

# Fit Proportional Hazards Regression Model , coxph {survival} 
surv<- Survival(f)
nom <- nomogram(f, fun=list(function(x) surv(1, x), 
                            function(x) surv(3, x),
                            function(x) surv(5, x)), 
                lp=F, funlabel=c("1-year survival", 
                                 "3-yearsurvival",
                                 "5-year survival"))
# maxscale=10, 
plot(nom, tcl=-0.5, label.every=2)
pdf('04-1-doModel-lasso//out-plot/multicox_nomogram_2.pdf',width = 13)
plot(nom, tcl=-0.5, label.every=2)
dev.off()



## 3.3 ROC----
library(timeROC)
head(dat_cox)
new_dat <- dat_cox[, c('event', 'time','Risk_score')]
head(new_dat)
## need 3 cols，time、event and risk scores
result <- with(new_dat, timeROC(T=time,
                                delta=event,
                                marker=Risk_score,
                                cause=1,
                                times = c(1, 3, 5),
                                iid = TRUE))
#identical(c(result$TP[,1],result$TP[,2],result$TP[,3]),as.numeric(result$TP))
dat = data.frame(fpr = as.numeric(result$FP),
                 tpr = as.numeric(result$TP),
                 time = rep(as.factor(c(1,3,5)),each = nrow(result$TP)))
ggplot() + 
  geom_smooth(data = dat, 
              aes(x = fpr, y = tpr, color = time), 
              size = 1,
              method = "loess",
              se = FALSE) + 
  scale_color_manual(name = NULL,
                     values = c("#92C5DE", "#F4A582", "#66C2A5"),
                     labels = paste0("AUC of ",c(1,3,5),"-year survival: ",
                                     format(round(result$AUC,2),nsmall = 2)))+
  geom_line(aes(x = c(0, 1), y = c(0,1)), 
            color = "grey",
            linetype = 'dotdash')+
  theme_ggstatsplot()+
  theme(axis.text = element_text(size = 10, face = 'bold'),
        axis.line = element_line(linetype = 1),
        panel.grid = element_blank(),
        legend.background = element_rect(linetype = 2, size = 0.2, colour = "black"),
        legend.position = c(0.665,0.135))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  coord_fixed()
ggsave('04-1-doModel-lasso//out-plot/smooth_ROC.pdf')
save(model, dat_cox, new_dat, file = '04-1-doModel-lasso//out-data/multicox_model.Rdata')



## 3.4 km----
risk_level <- as.factor(ifelse(new_dat$Risk_score > median(new_dat$Risk_score),'High','Low'))
new_dat$Risk_level <- risk_level 
sfit <- survfit(Surv(time, event)~Risk_level, data=new_dat)
sfit
summary(sfit)
save(new_dat,file ='04-1-doModel-lasso//out-data/basel-risk.Rdata')  #
identical(rownames(phe.2),rownames(new_dat))
phe.2$Risk_level = new_dat$Risk_level
write.csv(phe.2,file = "04-1-doModel-lasso/out-data/phe_risk.csv")




## more complicate figures.
survp=ggsurvplot(
  sfit,                     # survfit object with calculated statistics.
  legend.title = 'Risk level', 
  # legend = "top",#
  legend.labs = c('High', 'Low'),
  pval = T, #
  risk.table = TRUE, 
  risk.table.y.text = F,#
  xlab = "Time in years", #
  # xlim = c(0, 10), #
  break.time.by = 1, #
  size = 1.5, #
  ggtheme = theme_ggstatsplot(),
  palette="nejm", #
)
print(survp)
pdf('04-1-doModel-lasso//out-plot/multicox_KM.pdf', onefile = F)
print(survp)
dev.off()

## 3.5 ggrisk----
## https://cran.r-project.org/web/packages/ggrisk/ggrisk.pdf
## https://cloud.tencent.com/developer/article/1765625
library(ggrisk)
# save in pdf----
pdf('04-1-doModel-lasso//out-plot/riskscore.pdf', onefile = F)
ggrisk(model,
       color.A = c(low = "#0B5E9D", high = "#EA5205"),
       color.B = c(code.0 = "#0B5E9D", code.1 = "#EA5205"),
       color.C = c(low = "#0B5E9D", median = "white", high = "#EA5205"))
dev.off() 
save(model, dat_cox, new_dat, file = '04-1-doModel-lasso//out-data/multicox_model.Rdata')
# system('mv *.pdf Plot')
# system('mv *.Rdata Rdata')



