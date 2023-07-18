# Type first
##### A model using the previous model. Modeling this data is a high-low risk
# Settings
rm(list=ls())
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
load('04-1-doModel-lasso/out-data/multicox_model.Rdata')
model
rm(dat_cox,new_dat)
#   xdata
load('04-3-TCGA_check_model//out_data/input.Rdata')
library(limma)
exp.tcga.basal = normalizeBetweenArrays(exp.tcga.basal)
exp.tcga.basal[1:4,1:4]
cg = row.names(as.data.frame( model$coefficients)) 
cg
cg = gsub('_','-',cg)
cdat <- exp.tcga.basal[cg, ]
cdat <- na.omit(cdat)
xdata <- t(cdat)
xdata[1:4,1:4]
xdata= as.data.frame(xdata)
colnames(xdata) = gsub('-','_',colnames(xdata))  
colnames(xdata)

#   ydata
ydata <- phe.basal[,c("time","status")]
head(ydata)
identical(rownames(xdata), rownames(ydata))
## time transform to year
head(ydata)
boxplot(ydata$time)
kp <- ydata$time >0
table(kp)
xdata <- xdata[kp,]
ydata <- ydata[kp,]
names(ydata)[names(ydata) == 'status'] <- 'event'

cg = colnames(xdata)
n <- apply(xdata,2,scale)[,cg] %>%
  as.data.frame()
head(n)
rownames(n) <- rownames(xdata)
boxplot(n)
dat_cox <- cbind(ydata,n)
head(dat_cox)
## prepare variale for Surv(), paste gene names
multivariate <- paste(sort(colnames(n)), collapse = '+') 
## 1.2 multi-cox----
s <-  paste0(' Surv(time, event) ~  ', multivariate )
risk_score <- predict(model, type = 'risk', newdata = dat_cox)
dat_cox$Risk_score <- risk_score
risk_level <- as.factor(ifelse(dat_cox$Risk_score > median(dat_cox$Risk_score),'High','Low'))
dat_cox$Risk_level <- risk_level 


## 2.2 nomo----
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
# fun.at=c(0.95,0.9, 0.85, 0.8, 0.75, 0.7, 0.6, 0.5)
pdf('04-3-TCGA_check_model/out_plot//multicox_nomogram.pdf',width = 13)
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
ggsave('04-3-TCGA_check_model/out_plot//smooth_ROC.pdf')


## 3.4 km----
risk_level <- as.factor(ifelse(new_dat$Risk_score > median(new_dat$Risk_score),'High','Low'))
new_dat$Risk_level <- risk_level 
sfit <- survfit(Surv(time, event)~Risk_level, data=new_dat)
sfit
summary(sfit)
save(new_dat,file ='04-3-TCGA_check_model/out_data//basel-risk.Rdata') 
phe=fread('data_input/TCGA//TCGA.BRCA.sampleMap_BRCA_clinicalMatrix',data.table = F)
phe[1:4,1:4]
rownames(phe) = phe$sampleID
phe = phe[rownames(new_dat),]
identical(rownames(phe),rownames(new_dat))
phe$Risk_level = new_dat$Risk_level
write.csv(phe,file = "04-3-TCGA_check_model/out_data//phe_risk_tcga.csv")



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
pdf('04-3-TCGA_check_model/out_plot//multicox_KM.pdf', onefile = F)
print(survp)
dev.off()


##################
########################################################
# rm(list=ls())
cg <- read.table('04-1-doModel-lasso/out-data//model_lasso_cg_genes.txt')[, 1]
cg   

#  risk
data_fit <- cbind(ydata,n)
s <-  paste0('Surv(time, event) ~  ', multivariate )
fit <- coxph(as.formula(s), data = data_fit )
summary(model, data = data_fit)
library(rms)
dc <- datadist(data_fit);dc
ggrisk(fit)
## 3.5 ggrisk----
## https://cran.r-project.org/web/packages/ggrisk/ggrisk.pdf
## https://cloud.tencent.com/developer/article/1765625
library(ggrisk)
pdf('04-3-TCGA_check_model/out_plot//riskscore.pdf', onefile = F)
ggrisk(fit,
       color.A = c(low = "#0B5E9D", high = "#EA5205"),
       color.B = c(code.0 = "#0B5E9D", code.1 = "#EA5205"),
       color.C = c(low = "#0B5E9D", median = "white", high = "#EA5205"))
dev.off() 
save(model, dat_cox, new_dat, file = '04-3-TCGA_check_model/out_data//multicox_model.Rdata')
# system('mv *.pdf Plot')
# system('mv *.Rdata Rdata')

