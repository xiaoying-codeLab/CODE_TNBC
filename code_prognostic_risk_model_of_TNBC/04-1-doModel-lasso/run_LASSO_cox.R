lasso_cox <- function(xdata, ydata){
# 1. LASSO----   
library(glmSparseNet)
## 1.1 just see---
model_lasso <- glmnet(xdata,
                      Surv(ydata$time, ydata$event),
                      family = "cox",
                      alpha = 1)
plot(model_lasso, label = T)
plot(model_lasso, xvar = "lambda", label = T)  # lasso1
 
## 1.2 creat model---- 
fit <- cv.glmHub(xdata, 
                 Surv(ydata$time, ydata$event),
                 family  = 'cox',
                 lambda = buildLambda(1),
                 network = 'correlation',
                 network.options = networkOptions(cutoff = .3,
                                                  min.degree = .2))
print(fit)
plot(fit)

## 1.3 verity----
coefs.v <- coef(fit, s = 'lambda.min')[,1] %>% { .[. != 0]} 
coefs.v %>% {
  data.frame(gene.name   = names(.),
             coefficient = .,
             stringsAsFactors = FALSE)
} %>%
  arrange(gene.name) %>%
  knitr::kable()

coefs.v

ydata$status = ydata$event

separate2GroupsCox(as.vector(coefs.v),
                   xdata[, names(coefs.v)],
                   ydata,
                   plot.title = 'Full dataset', 
                   legend.outside = FALSE)
ggsave('04-1-doModel-lasso/out-plot//model_lasso_min.pdf')

save(coefs.v, fit, ydata, xdata,
     file = '04-1-doModel-lasso/out-data//coefs.v_lasso_model.Rdata')

write.table(sort(names(coefs.v)),
            file = '04-1-doModel-lasso/out-data//model_lasso_cg_genes.txt',
            row.names = F,col.names = F,quote = F)

}




