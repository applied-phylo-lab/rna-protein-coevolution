library(PCMBase,lib.loc = "~/R_dev")
library(PCMBaseCpp)
library(LaplacesDemon,lib.loc = "~/R_dev/")
require(coda)
library(tidyverse)
library(phytools)
library(extraDistr)
library(parallel)
source("local_functions.R")

prior_ln <- function(params)
{
  params.present <- names(params)
  x0.index <- which(str_detect(params.present,"(X0)|(Theta)"))
  h.index <- which(str_detect(params.present,"(H)|(Sigma)"))
  pr <- 0
  if (length(x0.index) > 0)
  {
    pr <- pr + sum(dnorm(params[x0.index],mean=0,sd=1,log=T)) # assuming same prior
  }
  if (length(h.index) > 0)
  {
    pr  <- pr + sum(dlnorm(params[h.index],meanlog = 0.25,sdlog = 1.5,log=T))
  }
  return(pr)
}

PGF <- function(Data)
{
  H_A <- rnorm(1,mean = 0.25,sd = 1.5)
  H_BA <- rnorm(1,mean = 0.25,sd = 1.5)
  H_AB <- rnorm(1,mean = 0.25,sd = 1.5)
  H_B <- rnorm(1,mean = 0.25,sd = 1.5)
  Sigma_A <- rnorm(1,mean = 0.25,sd = 1.5)
  Sigma_AB <- rnorm(1,mean = 0.25,sd = 1.5)
  Sigma_B <- rnorm(1,mean = 0.25,sd = 1.5)
  Theta <- rnorm(n=Data$num.genes*2,mean=0,sd=1)
  return(c(H_A,H_BA,H_AB,H_B,Sigma_A,Sigma_AB,Sigma_B,Theta))
}

Model_attempt <- function(parm, Data)
{
  ### Parameters
  h <- exp(parm[Data[["pos.h"]]])
  sigma <- exp(parm[Data[["pos.sigma"]]])
  #prior <- Data[["prior"]]
  ll_fun <- Data[["ll_fun"]]
  prop <- Data[["prop"]]
  gene.ll <- Data[["Gene.LL"]]
  if (prop[1] >= 8 && !is.null(gene.ll))
  {
    gene.index <- floor((prop[1] - 7 - 1)/2) + 1
    other.ll <- sum(gene.ll[-gene.index])
    current.ll <- ll_fun[[gene.index]](c(parm[prop[1]],parm[prop[2]],h,parm[prop[1]],parm[prop[2]],sigma))
    LL <- current.ll + other.ll
    gene.ll[gene.index] <- current.ll
    
  } else {
    num.genes <- (length(parm)-7)/2
    gene.ll <- unlist(lapply(1:num.genes, function(current.loci)
    {
      par.index <- 7 + 2*(current.loci-1) + 1
      ll <- ll_fun[[current.loci]](c(parm[par.index],parm[par.index+1],h,parm[par.index],parm[par.index+1],sigma))
      ll
    })
    )
    LL <- sum(gene.ll)
  }
  hyper.pr <-  sum(dlnorm(c(h,sigma),meanlog = 0.25,sdlog = 1.5,log=T))
  pr <- sum(dnorm(parm[8:length(parm)],mean = 0,sd = 1,log=T))
  LP <- LL + pr + hyper.pr# + 
  if (Data$correct.transform)
  {
    LP <- LP + sum(parm[Data[["pos.h"]]]) + sum(parm[Data[["pos.sigma"]]])
  }
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=1, parm=parm,Parameter.LL=gene.ll)
  return(Modelout)
}

Model <- function(parm, Data)
{
  ### Parameters
  h <- exp(parm[Data[["pos.h"]]])
  sigma <- exp(parm[Data[["pos.sigma"]]])
  #prior <- Data[["prior"]]
  ll_fun <- Data[["ll_fun"]]
  prop <- Data[["prop"]]
  #gene.ll <- Data[["Gene.LL"]]
  
  num.genes <- (length(parm)-7)/2
  gene.ll <- unlist(lapply(1:num.genes, function(current.loci)
    {
      par.index <- 7 + 2*(current.loci-1) + 1
      ll <- ll_fun[[current.loci]](c(parm[par.index],parm[par.index+1],h,parm[par.index],parm[par.index+1],sigma))
      ll
    })
  )
  
  LL <- sum(gene.ll)
  
  hyper.pr <-  sum(dlnorm(c(h,sigma),meanlog = 0.25,sdlog = 1.5,log=T))
  pr <- sum(dnorm(parm[8:length(parm)],mean = 0,sd = 1,log=T))
  LP <- LL + pr + hyper.pr
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=1, parm=parm,Parameter.LL=LL)
  return(Modelout)
}

set.seed(20)
ba.tree <- read.nexus("../Data/tree_11sp_noGpig.nex")
ba.tree <- force.ultrametric(ba.tree,method = "extend")

n <- 50 # number of loci
# Simulate the root state of the 2 traits across 10 loci
x0.a <- rnorm(n=n,mean=0,sd=1)
x0.b <- rnorm(n=n,mean=0,sd=1)

# Assume at stationarity
theta.a <- x0.a
theta.b <- x0.b

# Use same alpha and sigma as above

sigma.a <- 0.5
sigma.ab <- 0.1
sigma.b <- 0.25
alpha.a <- 0.75
alpha.b <- 0.5
alpha.ab <- 0.2
alpha.ba <- 0.4


data.pcm <- vector(mode="list",length=n)
likFun.list <- vector(mode="list",length=n)
for (i in 1:n)
{
  true.model <- PCM(model="OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",k=2)
  true.param <- c(x0.a[i],x0.b[i],alpha.a,alpha.ba,alpha.ab,alpha.b,theta.a[i],theta.b[i],sigma.a,sigma.ab,sigma.b)
  PCMParamLoadOrStore(true.model, true.param, offset=0, load=T)
  
  data.pcm[[i]] <- PCMSim(tree = ba.tree,model = true.model,X0 = true.model$X0)
  
  modelOU <- PCM(model="OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Global_Theta__Global_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",k=2)
  
  likFun.list[[i]] <- PCMCreateLikelihood(data.pcm[[i]], ba.tree, modelOU,metaI = PCMInfoCpp)
}



true.param <- c(alpha.a,
                alpha.ba,
                alpha.ab,
                alpha.b,
                sigma.a,
                sigma.ab,
                sigma.b)
names.true.param <- c("H_A","H_BA","H_AB","H_B","Sigma_A","Sigma_AB","Sigma_B")

for (i in 1:n)
{
  true.param <- c(true.param,theta.a[i],theta.b[i])
  names.true.param <- c(names.true.param,paste0(c("Theta_A","Theta_B"),"_",i))
}
names(true.param) <- names.true.param

pos.h <- grep("H",names(true.param)) 
pos.sigma <- grep("Sigma",names(true.param))
pos.theta <- grep("Theta",names(true.param))

MyData_hier <- list(ll_fun = likFun.list,
                    mon.names = "LP",
                    parm.names=names(true.param), 
                    pos.h = pos.h, 
                    pos.sigma=pos.sigma, 
                    pos.theta=pos.theta,
                    prop = 1:length(true.param),
                    prior=prior_ln,
                    num.genes = n,
                    correct.transform = F,
                    N=length(ba.tree$tip.label)*n*2)


blockwise.sample.list <- list()
blockwise.sample.list[[1]] <- c(1,5)
blockwise.sample.list[[2]] <- c(4,7)
blockwise.sample.list[[3]] <- c(2,3,6)
#blockwise.sample.list[[4]] <- pos.theta
num.other <- length(blockwise.sample.list)
for (i in 1:n)
{
  pos.theta.i <- grep(paste0("Theta_[AB]_",i,"$"),names(true.param))
  blockwise.sample.list[[i+num.other]] <- pos.theta.i
}




log.true.param <- c(log(true.param[1:7]),true.param[8:length(true.param)])

fit.RAM.hier.local <- LaplacesDemon_local(Model_attempt,
                                          MyData_hier,
                                          Initial.Values = log.true.param,
                                          Iterations = 20000,
                                          Algorithm = "RAM",
                                          Thinning = 1,
                                          Debug = list(DB.Model = F,
                                                       DB.chol = T),
                                          Specs = list(alpha.star = 0.44,
                                                       B = blockwise.sample.list,
                                                       Dist="t",
                                                       gamma = 2/3,
                                                       n=0
                                          ))

save(fit.RAM.hier.local,file="mcmc_11_species_50_genes_0.44_acceptance_20000_blockwise_lognormal_prior_do_not_apply_jacobian.Rda")
pdf("/data2/cope/rna-protein-coevolution/MCMC/example_mcmc_laplace_demon_hier_blockwise_11_species_50_genes_0.44_acceptance_20000_lognormal_prior_do_not_apply_jacobian.pdf")
plot(fit.RAM.hier.local$Monitor[,"LP"],type="l",ylab="Log(Posterior)")
plot(fit.RAM.hier.local$Posterior1[,"H_A"],type="l",ylab=expression(alpha["A"]))
abline(h=log.true.param[1])
plot(fit.RAM.hier.local$Posterior1[,"H_BA"],type="l",ylab=expression(alpha["BA"]))
abline(h=log.true.param[2])
plot(fit.RAM.hier.local$Posterior1[,"H_AB"],type="l",ylab=expression(alpha["AB"]))
abline(h=log.true.param[3])
plot(fit.RAM.hier.local$Posterior1[,"H_B"],type="l",ylab=expression(alpha["B"]))
abline(h=log.true.param[4])
plot(fit.RAM.hier.local$Posterior1[,"Sigma_A"],type="l",ylab=expression(sigma["A"]))
abline(h=log.true.param[5])
plot(fit.RAM.hier.local$Posterior1[,"Sigma_AB"],type="l",ylab=expression(sigma["AB"]))
abline(h=log.true.param[6])
plot(fit.RAM.hier.local$Posterior1[,"Sigma_B"],type="l",ylab=expression(sigma["B"]))
abline(h=log.true.param[7])

for (i in 8:ncol(fit.RAM.hier.local$Posterior1))
{
  plot(fit.RAM.hier.local$Posterior1[,i],type="l",ylab=expression(theta))
  abline(h=log.true.param[i])
}
dev.off()
