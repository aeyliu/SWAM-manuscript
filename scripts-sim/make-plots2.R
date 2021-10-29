#same as make-plots.R but now includes the extra 9000 runs

library(ggplot2)

n.tissues = c(2,4,6,8,10,20,30,40,50)
n.ss = c(50,100,150,200,250,300,350,400,450,500)
r2 = c(0.01,0.02,0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5)
r2.shared = c(0.0001,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.9999)
n.tissues.relevant = c(0,1,2,3,4,5,6,7,8,9)

swam.n.tissues = vector(mode="list",length=length(n.tissues))
swam.n.ss = vector(mode="list",length=length(n.ss))
swam.r2 = vector(mode="list",length=length(r2))
swam.r2.shared = vector(mode="list",length=length(r2.shared))
swam.n.tissues.relevant = vector(mode="list",length=length(n.tissues.relevant))

utmost.n.tissues = vector(mode="list",length=length(n.tissues))
utmost.n.ss = vector(mode="list",length=length(n.ss))
utmost.r2 = vector(mode="list",length=length(r2))
utmost.r2.shared = vector(mode="list",length=length(r2.shared))
utmost.n.tissues.relevant = vector(mode="list",length=length(n.tissues.relevant))


for(i in 1:length(n.tissues))
{
 temp9000 = read.table(paste("/net/snowwhite/home/aeyliu/pima/prediXcan/simulations/utmost-comparison/swam-results/n-tissues/n-tissues-",n.tissues[i],".txt",sep=""),header=T)
 temp1000 = read.table(paste("/net/snowwhite/home/aeyliu/pima/prediXcan/simulations/utmost-comparison/swam-results-1000/n-tissues/n-tissues-",n.tissues[i],".txt",sep=""),header=T)
 swam.n.tissues[[i]] = rbind(temp9000,temp1000)
 temp9000 = read.table(paste("/net/snowwhite/home/aeyliu/pima/prediXcan/simulations/utmost-comparison/results/n-tissues-results-",n.tissues[i],".txt",sep=""),header=F)
 temp1000 = read.table(paste("/net/snowwhite/home/aeyliu/pima/prediXcan/simulations/utmost-comparison/results-1000/n-tissues-results-",n.tissues[i],".txt",sep=""),header=F)
 utmost.n.tissues[[i]] = rbind(temp9000,temp1000)
}


for(i in 1:length(n.ss))
{
 temp9000 = read.table(paste("/net/snowwhite/home/aeyliu/pima/prediXcan/simulations/utmost-comparison/swam-results/sample-size/sample-size-",n.ss[i],".txt",sep=""),header=T)
 temp1000 = read.table(paste("/net/snowwhite/home/aeyliu/pima/prediXcan/simulations/utmost-comparison/swam-results-1000/sample-size/sample-size-",n.ss[i],".txt",sep=""),header=T)
 swam.n.ss[[i]] = rbind(temp9000,temp1000)
 temp9000 = read.table(paste("/net/snowwhite/home/aeyliu/pima/prediXcan/simulations/utmost-comparison/results/sample-size-results-",n.ss[i],".txt",sep=""),header=F)
 temp1000 = read.table(paste("/net/snowwhite/home/aeyliu/pima/prediXcan/simulations/utmost-comparison/results-1000/sample-size-results-",n.ss[i],".txt",sep=""),header=F)
 utmost.n.ss[[i]] = rbind(temp9000,temp1000)

}


for(i in 1:length(r2))
{
 temp9000 = read.table(paste("/net/snowwhite/home/aeyliu/pima/prediXcan/simulations/utmost-comparison/swam-results/r2/r2-",r2[i],".txt",sep=""),header=T)
 temp1000 = read.table(paste("/net/snowwhite/home/aeyliu/pima/prediXcan/simulations/utmost-comparison/swam-results-1000/r2/r2-",r2[i],".txt",sep=""),header=T)
 swam.r2[[i]] = rbind(temp9000,temp1000)

 temp9000 = read.table(paste("/net/snowwhite/home/aeyliu/pima/prediXcan/simulations/utmost-comparison/results/r2-results-",r2[i],".txt",sep=""),header=F)
 temp1000 = read.table(paste("/net/snowwhite/home/aeyliu/pima/prediXcan/simulations/utmost-comparison/results-1000/r2-results-",r2[i],".txt",sep=""),header=F)
 utmost.r2[[i]] = rbind(temp9000,temp1000)

}

options("scipen"=10000)

for(i in 1:length(r2.shared))
{
 temp9000 = read.table(paste("/net/snowwhite/home/aeyliu/pima/prediXcan/simulations/utmost-comparison/swam-results/r2-shared/r2-shared-",r2.shared[i],".txt",sep=""),header=T)
 temp1000 = read.table(paste("/net/snowwhite/home/aeyliu/pima/prediXcan/simulations/utmost-comparison/swam-results-1000/r2-shared/r2-shared-",r2.shared[i],".txt",sep=""),header=T)
 swam.r2.shared[[i]] = rbind(temp9000,temp1000)

 temp9000 = read.table(paste("/net/snowwhite/home/aeyliu/pima/prediXcan/simulations/utmost-comparison/results/r2-shared-results-",r2.shared[i],".txt",sep=""),header=F)
 temp1000 = read.table(paste("/net/snowwhite/home/aeyliu/pima/prediXcan/simulations/utmost-comparison/results-1000/r2-shared-results-",r2.shared[i],".txt",sep=""),header=F)
 utmost.r2.shared[[i]] = rbind(temp9000,temp1000)
}


for(i in 1:length(n.tissues.relevant))
{
 temp9000 = read.table(paste("/net/snowwhite/home/aeyliu/pima/prediXcan/simulations/utmost-comparison/swam-results/relevant-tissues/relevant-tissues-",n.tissues.relevant[i],".txt",sep=""),header=T)
 temp1000 = read.table(paste("/net/snowwhite/home/aeyliu/pima/prediXcan/simulations/utmost-comparison/swam-results-1000/relevant-tissues/relevant-tissues-",n.tissues.relevant[i],".txt",sep=""),header=T)
 swam.n.tissues.relevant[[i]] = rbind(temp9000,temp1000)

 temp9000 = read.table(paste("/net/snowwhite/home/aeyliu/pima/prediXcan/simulations/utmost-comparison/results/relevant-tissues-results-",n.tissues.relevant[i],".txt",sep=""),header=F)
 temp1000 = read.table(paste("/net/snowwhite/home/aeyliu/pima/prediXcan/simulations/utmost-comparison/results-1000/relevant-tissues-results-",n.tissues.relevant[i],".txt",sep=""),header=F)
 utmost.n.tissues.relevant[[i]] = rbind(temp9000,temp1000)
}


check.threshold=function(vec,thres)
{
 return(length(which(vec>thres)))
}

#swam = swam.n.tissues
#utmost = utmost.n.tissues
#outpar = n.tissues


pv.to.cor = function(pv,n=200)
{
 val = qt((1-pv/2),n-2)
 cor = sqrt((val^2)/(val^2+n-2))
 return(cor)
}


pv.thresholds = c(0.1,0.05,0.01,0.0001,1e-5,1e-6)
cor.thresholds = pv.to.cor(pv.thresholds,200) #correlation thresholds for the p-values

results.list.sim = vector(mode="list",length=length(cor.thresholds))

for(i in 1:length(results.list.sim))
{
 results.list.sim[[i]] = compile.results(sim,cor.threshold=cor.thresholds[i])
}


plot.results = function(swam,utmost,outpar,outpar.name="")
{
 swam.counts = vector(length=length(swam))
 naive.counts = vector(length=length(swam))
 best.tissue.counts = vector(length=length(swam))
 utmost.counts = vector(length=length(utmost))
 single.counts = vector(length=length(swam))
 cor.thres = pv.to.cor(0.05,200)

 for(i in 1:length(swam))
 {
  swam.counts[i] = apply(swam[[i]],2,check.threshold,thres=cor.thres)['meta']
  single.counts[i] = apply(swam[[i]],2,check.threshold,thres=cor.thres)['single']
  naive.counts[i] = apply(swam[[i]],2,check.threshold,thres=cor.thres)['naive']
  best.tissue.counts[i] = apply(swam[[i]],2,check.threshold,thres=cor.thres)['best.tissue']
  utmost.counts[i] = apply(utmost[[i]],2,check.threshold,thres=cor.thres)
 }
 df = data.frame(rep(outpar,5),c(rep("SWAM",length(swam)),rep("Naive.Average",length(swam)),rep("Best.Tissue",length(swam)),rep("UTMOST",length(utmost)),rep("single",length(swam))),
			c(swam.counts,naive.counts,best.tissue.counts,utmost.counts,single.counts)/dim(swam[[1]])[1])
 names(df) = c("x","method","y")
 
 p=ggplot(df,aes(x=x,y=y,group=method))+geom_line(aes(color=method))+geom_point(aes(color=method))
 p=p+xlab(outpar.name)+ylab("Proportion of Predictable Genes")+ylim(0,1)
 print(p)
 
}

pdf("/net/snowwhite/home/aeyliu/pima/prediXcan/simulations/utmost-comparison/plots/sample-size.pdf")
plot.results(swam.n.ss,utmost.n.ss,n.ss,"Sample Size of Target Tissue")
dev.off()

pdf("/net/snowwhite/home/aeyliu/pima/prediXcan/simulations/utmost-comparison/plots/n-tissues-10k.pdf")
plot.results(swam.n.tissues,utmost.n.tissues,n.tissues,"Number of Total Tissues")

dev.off()

pdf("/net/snowwhite/home/aeyliu/pima/prediXcan/simulations/utmost-comparison/plots/r2-10k.pdf")
plot.results(swam.r2,utmost.r2,r2,"Genetic Heritability")
dev.off()

pdf("/net/snowwhite/home/aeyliu/pima/prediXcan/simulations/utmost-comparison/plots/r2.shared.pdf")
plot.results(swam.r2.shared,utmost.r2.shared,r2.shared,"Cross-tissue heritability")
dev.off()

pdf("/net/snowwhite/home/aeyliu/pima/prediXcan/simulations/utmost-comparison/plots/n-tissues-relevant.pdf")
plot.results(swam.n.tissues.relevant,utmost.n.tissues.relevant,n.tissues.relevant,"Number of Relevant Tissues")

dev.off()
