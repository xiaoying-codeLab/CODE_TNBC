# Settings
# Objective: To do COX on genes in batches
# Data source: Rdata
# Save data:
rm(list = ls())  
options(stringsAsFactors = F) 
library(survminer) 
library(survival)


load(file = '02-doCOX/out-data//2.input_survival_data.Rdata' )
expr = expr.BLIS
phe = phe_BLIS
expr[1:4,1:4]
head(phe) 

kp=apply(expr,1, function(x) sum(x>1) > 30)
table(kp)
expr <- expr[kp,] 
expr[1:4,1:4]
head(phe) 
dim(expr)
expr = expr[apply(expr, MARGIN=1, FUN =  function(xxx) {
  (sum(xxx==0)/length(xxx))<=0.2
}),]
dim(expr)


# 1. prepare data for coxph----
mySurv <- with(phe, Surv(time, event))
cox_results <-apply(expr , 1 , function(gene){
  # gene= as.numeric(expr[1,])
  group=ifelse(gene>median(gene),'high','low') 
  if( length(table(group))<2)
    return(NULL)
  survival_dat <- data.frame(group=group,# stage=phe$stage,
                             stringsAsFactors = F)
  m=coxph(mySurv ~ group, 
          # mySurv ~  stage+ group, 
          data =  survival_dat)
  
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)
  HRse <- HR * se
  
  #summary(m)
  tmp <- round(cbind(coef = beta, se = se, z = beta/se, p = 1 - pchisq((beta/se)^2, 1),
                     HR = HR, HRse = HRse,
                     HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                     HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
  return(tmp['grouplow',])
  
})

# 2. specify the value----
cox_results=t(cox_results)
head(cox_results)
table(cox_results[,4]<0.01)
table(cox_results[,4]<0.05)

save(cox_results, 
     file = '02-doCOX/out-data//2.batch_cox_results.Rdata')


