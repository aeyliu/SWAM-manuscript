args = vector(length=1)
args[1] = "/net/snowwhite/home/aeyliu/pima/prediXcan/simulations/utmost-comparison/results/relevant-tissues"



library(glmnet)
library(MASS)
library(ggplot2)
library(leaps)
pkgs <- list("glmnet", "doParallel", "foreach", "pROC")
lapply(pkgs, require, character.only = T)
registerDoParallel(cores = 4)

#functions
cor.to.pv = function(cors,n=465)
{
 val = abs(cors/sqrt((1-cors^2)/(n-2)))
 return(2*(1-pt(val,n-2)))
}

check.threshold = function(v,threshold=0.1965)
{
 v = na.omit(v)
 temp = sum(v>threshold)/length(v)
 return(temp)
}


cor.to.pv = function(cors,n=465)
{
 val = abs(cors/sqrt((1-cors^2)/(n-2)))
 return(2*(1-pt(val,n-2)))
}

pv.to.cor = function(pv,n=200)
{
 val = qt((1-pv/2),n-2)
 cor = sqrt((val^2)/(val^2+n-2))
 return(cor)
}

#generate eQTL tissue structure
#the structure variable specifies how related the tissues are
gen.tissues.eqtl = function(n.tissues,n.snps,structure)
{
 eqtl.mat = matrix(nrow=n.tissues,ncol=n.snps)
 #generate first tissue (target)
 eqtl.mat[1,]=rbinom(n.snps,1,structure[1])
 for(i in 2:n.tissues)
 {
  temp = rbinom(n.snps,1,structure[i])
  for(j in 1:n.snps)
  {
   if(temp[j]==1) #if they match
   {
    eqtl.mat[i,j]=eqtl.mat[1,j] #shared eqtl
   }
   if(temp[j]==0) #if they don't match
   {
    eqtl.mat[i,j]=1-eqtl.mat[1,j] #shared eqtl
   }
  }
 }
 return(eqtl.mat)
}

#generate X matrices
gen.X = function(n,mafs)
{
 X = matrix(nrow=n,ncol=length(mafs))
 for(i in 1:length(mafs))
 {
  X[,i] = rbinom(n,2,mafs[i])
 } 
 return(X)
}
gen.Y = function(beta,X,env)
{
 Y = matrix(nrow=dim(X)[1],ncol=dim(beta)[1])
 for(i in 1:dim(beta)[1])
 {
  Y[,i] = X%*%beta[i,] + rnorm(dim(X)[1],0,env[i])
 }
 return(Y)
}

calc.models.enet = function(Y,X)
{
 models = vector(mode="list",length=dim(Y)[2])
 for(i in 1:length(models))
 {
  y = Y[,i]
  a <- seq(0.1, 0.9, 0.05)
  search <- foreach(j = a, .combine = rbind) %dopar% {
    cv <- cv.glmnet(X, y, family = "gaussian", nfold = 10, type.measure = "deviance", paralle = TRUE, alpha = j)
    data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = j)
  }
  cv3 <- search[search$cvm == min(search$cvm), ]
  model.temp = glmnet(X,y,family="gaussian",alpha=cv3$alpha,lambda=cv3$lambda.1se)
  models[[i]] = as.numeric(coef(model.temp))[-1]
 }
 return(models)
}


calc.models.leaps = function(Y,X,sample.sizes,max.val,iter,outpar)
{
 models = vector(mode="list",length=dim(Y)[2])
 n = dim(X)[1]
 for(i in 1:length(models))
 {
  ind = sample(1:n,floor(sample.sizes[i]*n)) #sample based on sample size
  y = Y[ind,i]
  X.temp = X[ind,]
  #write out data for utmost to analyze
  y.utmost = data.frame(ind,y)
  write.table(y.utmost,paste(args[1],"/relevant-tissues-",outpar,"/run",iter,"/Y/Y",i,".txt",sep=""),row.names=F,quote=F,sep=" ",col.names=F)


  df.test = data.frame(y,X.temp)
  fit = regsubsets(y~.,data=df.test,nbest=5,nvmax=max.val)
  vals = apply(summary(fit,matrix.logical=TRUE)$outmat,2,sum)
  reg.index = which(vals>0)
  X.names = names(reg.index)
  model.temp = lm(y~X.temp[,reg.index],data=df.test)
  names(model.temp$coefficients)[2:length(model.temp$coefficients)]=X.names
  models[[i]] = summary(model.temp)
  #fit = lm(y~.,data=df.test)
  #step=stepAIC(fit,direction="both")
 }
 return(models)
}

calc.models.ideal = function(Y,X,sample.sizes,eqtl.index)
{
 models = vector(mode="list",length=dim(Y)[2])
 n = dim(X)[1]
 for(i in 1:length(models))
 {
  ind = sample(1:n,floor(sample.sizes[i]*n)) #sample based on sample size
  y = Y[ind,i]
  X.temp = X[ind,]
  X.ind = which(as.vector(eqtl.index[,i])==1)
  df.temp = data.frame(y,X.temp)
  model.temp = lm(y~X.temp[,X.ind],data=df.temp)
  names(model.temp$coefficients)[2:length(model.temp$coefficients)] = names(data.frame(X.temp))[X.ind]
  models[[i]] = summary(model.temp)
 }
 return(models)
}


calc.models = function(Y,X,sample.sizes)
{
 models = vector(mode="list",length=dim(Y)[2])
 n = dim(X)[1]
 for(i in 1:length(models))
 {
  ind = sample(1:n,floor(sample.sizes[i]*n)) #sample based on sample size
  y = Y[ind,i]
  X.temp = X[ind,] 
  models[[i]] = summary(lm(y~X.temp))
 }

 return(models)
}

calc.models.target = function(target.index,Y)
{
 ind = c(1:dim(Y)[2])
 remaining.index = ind[-target.index]
 tau = vector(length=length(ind)) #cross tissue bias
 for(i in remaining.index)
 {
  temp=aov(lm(Y[,target.index]~Y[,i]))
  tau[i] = summary(temp)[[1]][2,3]
 }
 return(tau)
}

pred.expression = function(models,X)#models is in list form
{
 Y.pred = matrix(nrow=dim(X)[1],ncol=length(models))
 for(i in 1:length(models))
 {
  n.temp = dim(models[[i]]$coefficients)[1]
  n.temp1 = models[[i]]$df[3] #should be safer

  sig.index = which(models[[i]]$coefficients[2:n.temp,4]<0.3)
  if(length(sig.index)==0)
  {
   temp = models[[i]]$coefficients[-1,]
   temp = temp[order(temp[,4]),]
   temp = na.omit(temp)
   sig.index = which(models[[i]]$coefficients[2:n.temp,4]%in%temp[c(1,2,3),4]) #take the top three signals, regardless of p-value
  }
  beta.hat = rep(0,n.temp1-1)
  beta.hat[sig.index] = models[[i]]$coefficients[sig.index+1,1]

  Y.pred[,i] = X%*%beta.hat
 }
 return(Y.pred)
}

pred.expression.stepwise = function(models,X) #models is in list form, this needs to be the models obtained from the variable selection
{
 Y.pred = matrix(nrow=dim(X)[1],ncol=length(models))
 for(i in 1:length(models))
 {
  n.temp = dim(models[[i]]$coefficients)[1]
  n.temp1 = models[[i]]$df[3] #should be safer

  x.index = row.names(models[[i]]$coefficients)[-1]
  x.index = as.numeric(as.character(substr(x.index,2,nchar(x.index))))

  beta.hat = rep(0,dim(X)[2])
  beta.hat[x.index] = models[[i]]$coefficients[2:n.temp,1]

  Y.pred[,i] = X%*%beta.hat
 }
 return(Y.pred)
}


calc.aggregate = function(Y,weights)
{
 Y.agg = Y%*%weights
 return(Y.agg)
}

get.sigmas = function(mods)
{
 sigs = vector(length=length(mods))
 for(i in 1:length(mods))
 {
  sigs[i] = mods[[i]]$sigma
 }
 return(sigs)
}

sim.pred = function(n.tissue=10,n.tissue.irrelevant=0,n.samples=100,n.snps=100,n.snps.shared=1,n.snps.specific=1,r2=0.1,r2.cross=0.1,mean.effect.size=0.1,sample.sizes=rep(1,50),diag.val = 1,iter,outpar)
{
 redo.all=TRUE
 while(redo.all)
 {
  #generate training set
  #minor allele frequencies
  mafs = rbeta(n.snps,1,3)
  #effect sizes
  #try using fixed effect size first
  beta = rep(1,n.snps)
 
  #otherwise use random effect size
  #beta = rexp(n.snps,1/mean.effect.size)*(2*(rbinom(n.snps,1,0.5)-0.5))
  #beta = rexp(n.snps,1/mean.effect.size)
 
  redo = TRUE
  #genotypes
  while(redo)
  {
   X=gen.X(n.samples,mafs)
   if(!0%in%apply(X,2,sum))
   {
    redo=FALSE
   }
  }
  #shared SNPs between tissues
  shared.snp.index = sample(1:n.snps,n.snps.shared) #sample list of shared SNPs
  remaining.index = (1:n.snps)[-shared.snp.index] #candidate for remaining tissue-specific eQTLs
 
  #set SNP values to 1 for eqtl.index matrix
  eqtl.index = matrix(0,nrow=n.tissue,ncol=n.snps)
  eqtl.index.shared = matrix(0,nrow=n.tissue,ncol=n.snps)
  eqtl.index.specific = matrix(0,nrow=n.tissue,ncol=n.snps)
  eqtl.index.shared[,shared.snp.index] = rep(1,n.tissue)
  if(n.tissue.irrelevant>0) #if we have irrelevant tissues
  {
   eqtl.index.shared[(n.tissue-n.tissue.irrelevant+1):n.tissue,] = 0
  }

  eqtl.shared = t(eqtl.index.shared)*beta
  Y.genetic.shared = X%*%eqtl.shared
  
  sigma.specific = vector(length=n.tissue)
  sigma.irrelevant = vector(length=n.tissue)
  
  #now simulate the tissue-specific eQTLs
  for(i in 1:n.tissue)
  {
   if(n.tissue.irrelevant==0)
   {
    j = sample(remaining.index,n.snps.specific)
    eqtl.index.specific[i,j] = 1
    sigma.specific[i] = sum(2*(1-mafs[j])*mafs[j]) #sum up the individual variances
   }
   if(n.tissue.irrelevant>0)
   {
    if(i<(n.tissue-n.tissue.irrelevant+1)) #relevant tissues
    {
     j = sample(remaining.index,n.snps.specific)
     eqtl.index.specific[i,j] = 1
     sigma.specific[i] = sum(2*(1-mafs[j])*mafs[j]) #sum up the individual variances
    }
    if(i>=(n.tissue-n.tissue.irrelevant+1)) #irrelevant tissues
    {
     j = sample(remaining.index,n.snps.specific+n.snps.shared)
     eqtl.index.specific[i,j] = 1
     sigma.specific[i] = sum(2*(1-mafs[j])*mafs[j]) #sum up the individual variances
    }
   }
  }
 
  #now, calculate expected variance based on mafs
  sigma.shared = sum(2*(1-mafs[shared.snp.index])*mafs[shared.snp.index])
 
  #scaling factor to make expected heritability match specified parameters
  k = ((r2.cross*sigma.specific)/((1-r2.cross)*sigma.shared))^(-1)
  var.gen.expected = sigma.shared+(sigma.specific*k)[1]

  if(n.tissue.irrelevant>0)
  {
   k[(n.tissue-n.tissue.irrelevant+1):n.tissue]=(var.gen.expected/(sigma.specific))[(n.tissue-n.tissue.irrelevant+1):n.tissue]
  }

  #calculate tissue specific values
  eqtl.specific=t(eqtl.index.specific)*beta
  Y.genetic.specific = X%*%eqtl.specific
  eqtl.specific.scaled = eqtl.specific
  for(i in 1:n.tissue)
  {
   Y.genetic.specific[,i] = Y.genetic.specific[,i]*sqrt(k[i]) #scale by sqrt(k[i]) as per the equation
   eqtl.specific.scaled[,i] = eqtl.specific[,i]*sqrt(k[i])
  }

  Y.genetic = Y.genetic.specific+Y.genetic.shared

  #check that the heritability is what we want it to be
  #shared.var.gen = apply(Y.genetic.shared,2,var)
  #specific.var.gen = apply(Y.genetic.specific,2,var)
  var.gen = apply(Y.genetic,2,var)


  #shared.var.gen/var.gen #this should be approximately r2.cross

  #now calculate environmental effect
  var.env = (1-r2)*var.gen.expected/r2
  Y.env = matrix(ncol=n.tissue,rnorm(n.samples*n.tissue,0,sqrt(var.env)))

  Y = Y.genetic+Y.env

  #calculate eqtl.index for later Y simulation
  eqtl.effects = eqtl.shared + eqtl.specific.scaled #this should produce the same value of Y.genetic when doing Y.genetic = X%*%eqtl.effects

  redo.all=FALSE

  #generate "prediction" models 
  models = calc.models.leaps(Y=Y,X=X,sample.sizes=sample.sizes,max.val=n.snps.shared+n.snps.specific,iter,outpar)
  for(i in 1:length(models))
  {
   if(dim(models[[i]]$coefficients)[1]<2) #if the regression didn't pick up any SNPs 
   {
    redo.all=TRUE
   }
  }
 
  #generate "true" prediction models - knowing where the eQTLs are
  eqtl.index = eqtl.effects
  eqtl.index[eqtl.index>0] = 1
  models.ideal = calc.models.ideal(Y=Y,X=X,sample.sizes=sample.sizes,eqtl.index)
  for(i in 1:length(models.ideal))
  {
   if(dim(models.ideal[[i]]$coefficients)[1]<2) #if the regression didn't pick up any SNPs 
   {
    redo.all=TRUE
   }
  }

 X.utmost = data.frame(1:dim(X)[1],X)
 write.table(X.utmost,paste(args[1],"/relevant-tissues-",outpar,"/run",iter,"/X",".txt",sep=""),row.names=F,quote=F,sep=" ",col.names=F)

 }
 

 #calculate weights based on known models

 #predict expression for all tissues
 Y.pred = pred.expression.stepwise(models,X)

 Y.pred.ideal = pred.expression.stepwise(models.ideal,X)

 #correlation based weights
 weights.cor = vector(length=n.tissue)
 weights.cor.ideal = vector(length=n.tissue)
 weights.cor.best = vector(length=n.tissue)

 
 X.cv = gen.X(n.samples,mafs)
 Y.env = matrix(ncol=n.tissue,rnorm(n.tissue*n.samples,0,sqrt(var.env)))
 Y.cv = X.cv%*%eqtl.effects+Y.env
 models1=vector(mode="list",length=2)
 models1[[1]]=models[[1]]
 models1[[2]]=models[[2]]
 Y.pred.cv = pred.expression.stepwise(models1,X.cv)
 models2=vector(mode="list",length=2)
 models2[[1]]=models.ideal[[1]]
 models2[[2]]=models.ideal[[2]]
 Y.pred.cv.ideal = pred.expression.stepwise(models2,X.cv) #ideal model
 #calculated "Cross-validated" R-squared for target tissue
 
 weights.cor[1] = sqrt(models[[1]]$r.squared)
 weights.cor.ideal[1] = sqrt(models.ideal[[1]]$r.squared)
 weights.cor.best[1] = cor(Y.cv[,1],Y.pred.cv[,1])

 for(i in 2:n.tissue)
 {
  weights.cor[i] = cor(Y.pred[,i],Y[,1])
  weights.cor.best[i] = cor(Y.pred[,i],Y[,1])
  weights.cor.ideal[i] = cor(Y.pred.ideal[,i],Y[,1])
 }

 #weights based on correlation and tissue structure
 C.y = cov(Y.pred)
 weights.cor.tissue = as.numeric(ginv(C.y+diag(diag.val,n.tissue))%*%weights.cor)
 weights.cor.tissue[weights.cor.tissue<0]=0

 C.y.ideal = cov(Y.pred.ideal)
 weights.cor.tissue.ideal = as.numeric(ginv(C.y.ideal+diag(diag.val,n.tissue))%*%weights.cor.ideal)
 weights.cor.tissue.ideal[weights.cor.tissue.ideal<0]=0

 weights.cor.best.tissue = as.numeric(ginv(C.y+diag(diag.val,n.tissue))%*%weights.cor.best)
 weights.naive = rep(1/n.tissue,n.tissue)
 best.index =  which(weights.cor.best.tissue==max(weights.cor.best.tissue))
 weights.best.tissue = rep(0,n.tissue)
 weights.best.tissue[best.index] = 1
 weights.target = c(1,rep(0,n.tissue-1))

 #generate the new data (test set)
 X.test = gen.X(n.samples,mafs)
 Y.pred.test = pred.expression.stepwise(models,X.test)
 Y.env = matrix(ncol=n.tissue,rnorm(n.tissue*n.samples,0,sqrt(var.env)))
 Y.test = X.test%*%eqtl.effects + Y.env

 #write out utmost data
 X.utmost.test = data.frame(1:dim(X.test)[1],X.test)
 write.table(X.utmost.test,paste(args[1],"/relevant-tissues-",outpar,"/run",iter,"/X.test",".txt",sep=""),row.names=F,quote=F,sep=" ",col.names=F)

 for(index in 1:dim(Y.test)[2])
 {
  Y.utmost.test = data.frame(1:dim(Y.test)[1],Y.test[,index])
  write.table(Y.utmost.test,paste(args[1],"/relevant-tissues-",outpar,"/run",iter,"/Y.test/Y",index,".txt",sep=""),row.names=F,quote=F,sep=" ",col.names=F)

 }
 #ideal setting
 Y.pred.test.ideal = pred.expression.stepwise(models.ideal,X.test)

 n.results=6
 Y.naive = calc.aggregate(Y.pred.test,weights.naive)
 Y.best.tissue = calc.aggregate(Y.pred.test,weights.best.tissue)
 Y.meta = calc.aggregate(Y.pred.test,weights.cor.tissue)
 Y.single = calc.aggregate(Y.pred.test,weights.target)
 Y.naive.ideal = calc.aggregate(Y.pred.test.ideal,weights.naive)
 Y.meta.ideal = calc.aggregate(Y.pred.test.ideal,weights.cor.tissue.ideal)

 names.results = c("naive","best.tissue","meta","single","naive.ideal","meta.ideal")
 results = vector(length=(n.results))
 
 results[1] = cor(Y.test[,1],Y.naive)
 results[2] = cor(Y.test[,1],Y.best.tissue)
 results[3] = cor(Y.test[,1],Y.meta)
 results[4] = cor(Y.test[,1],Y.single)
 results[5] = cor(Y.test[,1],Y.naive.ideal)
 results[6] = cor(Y.test[,1],Y.meta.ideal)
 
 target.proportion.weight = weights.cor.tissue[1]/sum(weights.cor.tissue)

 names(results) = names.results


 #return results and information
 info = list(mafs,X,Y,models,results,target.proportion.weight)

 return(info)
}



#####################
#run the simulation
run.sim=function(n.times=1000,n.tissue=10,n.tissue.irrelevant=0,n.samples=100,n.snps=100,n.snps.shared=1,n.snps.specific=1,r2=0.1,r2.cross=0.1,mean.effect.size=0.1,sample.sizes=rep(1,50),diag.val=1,outpar)
{
 results = matrix(nrow=n.times,ncol=6)
 target.proportion.weight = vector(length=n.times)
 for(i in 1:n.times)
 {
  if(i%%100==0)
  {
   print(paste("Processing ", i, " of ", n.times, sep=""))
  }
  test = sim.pred(n.tissue=n.tissue,n.tissue.irrelevant=n.tissue.irrelevant,n.samples=n.samples,n.snps=n.snps,n.snps.shared=n.snps.shared,n.snps.specific=n.snps.specific,r2=r2,r2.cross=r2.cross,mean.effect.size=mean.effect.size,sample.sizes=sample.sizes,diag.val = diag.val,iter=i,outpar=outpar)
  results[i,]=test[[5]]
  target.proportion.weight[i] = test[[6]]
 }
 output = data.frame(results,target.proportion.weight)
 names(output) = c("naive","best.tissue","meta","single","naive.ideal","meta.ideal","target.weight")
 return(output)
}



compile.results=function(results.list,cor.threshold=0.1)
{
 threshold.results = matrix(nrow=length(results.list),ncol=6)
 mean.results = matrix(nrow=length(results.list),ncol=6)
 for(i in 1:length(results.list))
 {
  threshold.results[i,] = as.numeric(apply(na.omit(results.list[[i]][,1:6]),2,check.threshold,threshold=cor.threshold))
  mean.results[i,] = as.numeric(apply(na.omit(results.list[[i]][,1:6]),2,mean))
 }
 return(list(threshold.results,mean.results))
}

n.tiss=10
n.tiss.irrelevant=c(9,8,7,6,5,4,3,2,1,0)
n.snps=35 #total number of potential eQTL SNPs
n.snps.shared=5 #number of shared eQTL across all tissues
n.snps.specific=5
r2=0.1 #total genetic r2
r2.c=0.5 #cross-tissue r2
mean.effect.size=0.1
dv = 3
N=length(n.tiss.irrelevant)
sim.x = 9-n.tiss.irrelevant
ss.in = c(0.5,rep(1,n.tiss-1))

sim = vector(mode="list",length=N)
for(i in 1:N)
{
 sim[[i]] = run.sim(n.times=9000,n.tissue=n.tiss,n.tissue.irrelevant=n.tiss.irrelevant[i],n.samples=200,n.snps=n.snps,n.snps.shared=1,n.snps.specific=1,r2=r2,r2.cross=r2.c,mean.effect.size=0.1,sample.sizes=ss.in,diag.val=dv,outpar=(9-n.tiss.irrelevant[i]))
}

pv.thresholds = c(0.1,0.05,0.01,0.0001,1e-5,1e-6)
cor.thresholds = pv.to.cor(pv.thresholds,200) #correlation thresholds for the p-values

results.list.sim = vector(mode="list",length=length(cor.thresholds))

for(i in 1:length(results.list.sim))
{
 results.list.sim[[i]] = compile.results(sim,cor.threshold=cor.thresholds[i])
}


#write out results for SWAM
for(i in 1:N)
{
 write.table(sim[[i]],paste("/net/snowwhite/home/aeyliu/pima/prediXcan/simulations/utmost-comparison/swam-results/relevant-tissues/relevant-tissues-",(9-n.tiss.irrelevant[i]),".txt",sep=""),quote=F,row.names=F,sep="\t")
}
