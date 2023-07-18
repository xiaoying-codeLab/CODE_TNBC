# Settings
# Objective: Visualize the gsva results and draw a graph
# Data Source:
# Save data:

rm(list = ls())  
options(stringsAsFactors = F) 
library(pheatmap) 


# 0. 
input <- "06-Immune-ssGSEA/out-data/step1_ssgseaOut.txt"  
immune <- read.table( input ,sep="\t",header=T,row.names=1,check.names=F)
immune[1:4, 1:4]

load('04-1-doModel-lasso/out-data///basel-risk.Rdata')
risk = new_dat
dim(risk)
risk[1:4, 1:4]


# 1. overlap----
a <- intersect(colnames(immune),rownames(risk))
risk <- risk[a,]
risk <- risk[order(risk$Risk_score),]  
immune <- immune[,rownames(risk)]
immune[1:4, 1:4]
table(colnames(immune)==rownames(risk))
Type <- risk$Risk_level
names(Type) <- colnames(immune)
Type <- as.data.frame(Type)


# 2. heatmap----
n <- t(scale(t(immune))) 
n[n>2]=2 
n[n< -2]= -2
n[1:4,1:4]


p1 <-  pheatmap(n, 
         annotation = Type, 
         show_colnames =F,
         cluster_cols =F,
         fontsize = 8,
         fontsize_row=8,
         fontsize_col=5)
p1
ggsave(p1,filename = '06-Immune-ssGSEA/out-plot/heatmap.pdf',height=6,width=15)

# 3. vioplot----
outTab=data.frame()
table(Type$Type)
lowNum=70
highNum=69
immune <- data.frame(t(immune),check.names = F)

pdf("06-Immune-ssGSEA/out-plot/vioplot_last_t.test.pdf",height=7,width=15)
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(immune))
y=c(1:ncol(immune))
plot(x,y,
     xlim=c(0,72),ylim=c(min(immune),max(immune)+0.02),
     main="",xlab="", ylab="Fraction",
     pch=24,
     col="white",
     xaxt="n")

library(vioplot)
pFilter=0.05
for(i in 1:ncol(immune)){
  #i = 1
  if(sd(immune[1:lowNum,i])==0){
    immune[1,i]=0.001
  }
  if(sd(immune[(lowNum+1):(lowNum+highNum),i])==0){
    immune[(lowNum+1),i]=0.001
  }
  lowData=immune[1:lowNum,i]
  highData=immune[(lowNum+1):(lowNum+highNum),i]
  vioplot(lowData,at=3*(i-1),lty=1,add = T,col = '#00A087')
  vioplot(highData,at=3*(i-1)+1,lty=1,add = T,col = '#E64B35')
  wilcoxTest=t.test(lowData,highData)
  p=wilcoxTest$p.value
  if(p<pFilter){
    cellPvalue=cbind(Cell=colnames(immune)[i],pvalue=p)
    outTab=rbind(outTab,cellPvalue)
  }
  mx=max(c(lowData,highData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}
legend("topright", 
       c("low", "high"),
       lwd=5,bty="n",cex=1.5,
       col=c("#00A087","#E64B35"))
text(seq(1,70,3),-0.1,xpd = NA,labels=colnames(immune),cex = 1,srt = 45,pos=2)
dev.off()




