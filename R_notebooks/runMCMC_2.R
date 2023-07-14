library(PCMBase,lib.loc = "~/R_dev")
library(PCMBaseCpp)
library(LaplacesDemon,lib.loc = "~/R_dev/")
require(coda)
library(tidyverse)
library(phytools)
source("local_functions.R")



Model_OUOU <- function(parm, Data)
{
  ### Parameters
  pos.psi <- Data[["pos.psi"]]
  h <- exp(parm[Data[["pos.h"]]])
  b <- parm[Data[["pos.b"]]]
  sigma <- exp(parm[Data[["pos.sigma"]]])
  ll_fun <- Data[["ll_fun"]]
  prop <- Data[["prop"]]
  gene.ll <- Data[["Gene.LL"]]
  if (prop[1] >= pos.psi[1] && !is.null(gene.ll))
  {
    gene.index <- floor((prop[1] - 5 - 1)/2) + 1
    other.ll <- sum(gene.ll[-gene.index])
    current.ll <- ll_fun[[gene.index]](c(parm[prop[1]] + b*parm[prop[2]],
                                         parm[prop[2]],
                                         h[1],
                                         0,
                                         -h[1] * b,
                                         h[2],
                                         parm[prop[1]] + b*parm[prop[2]],
                                         parm[prop[2]],
                                         sigma))
    LL <- current.ll + other.ll
    gene.ll[gene.index] <- current.ll
    
  } else {
    num.genes <- (length(parm)-5)/2
    gene.ll <- unlist(lapply(1:num.genes, function(current.loci)
    {
      par.index <- 5 + 2*(current.loci-1) + 1
      ll <- ll_fun[[current.loci]](c(parm[par.index] + b*parm[par.index+1],
                                     parm[par.index+1],
                                     h[1],
                                     0,
                                     -h[1] * b,
                                     h[2],
                                     parm[par.index] + b*parm[par.index+1],
                                     parm[par.index+1],
                                     sigma))
      ll
    })
    )
    LL <- sum(gene.ll)
  }
  hyper.pr <-  sum(dlnorm(c(h,sigma),meanlog = 0.25,sdlog = 1.5,log=T)) + dnorm(b,mean = 0,sd = 1,log=T)
  pr <- sum(dnorm(parm[pos.psi[1]:length(parm)],mean = 0,sd = 1,log=T))
  LP <- LL + pr + hyper.pr
  Modelout <- list(LP=LP,Dev=-2*LL, Monitor=c(LP,LL), yhat=1, parm=parm,Parameter.LL=gene.ll)
  return(Modelout)
}

set.seed(20)
ba.tree <- read.nexus("../Data/tree_11sp_noGpig.nex")
ba.tree <- force.ultrametric(ba.tree,method = "extend")
# tip.heights <- ips::tipHeights(ba.tree)
# ba.tree <- pbtree(n = 50,scale = unname(tip.heights[1]))


n <- 200 # number of loci

psi.a <- rnorm(n=n,mean=0,sd=1)
psi.b <- rnorm(n=n,mean=0,sd=1)

# Use same alpha and sigma as above

sigma.a <- 0.5
sigma.b <- 0.25
alpha.a <- 0.75
alpha.b <- 0.5
b <- 0.5

## Get root state
x0.b <- psi.b
x0.a <- psi.a + b*x0.b

model.type <- "OU__Global_X0__Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x"
#model.type <- "OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x"


data.pcm <- vector(mode="list",length=n)
likFun.list <- vector(mode="list",length=n)
for (i in 1:n)
{
  true.model <- PCM(model=model.type,k=2)
  true.param <- c(x0.a[i],
                  x0.b[i],
                  alpha.a,
                  0,
                  -alpha.a * b,
                  alpha.b,
                  x0.a[i],
                  x0.b[i],
                  sigma.a,
                  sigma.b)
  PCMParamLoadOrStore(true.model, true.param, offset=0, load=T)
  
  data.pcm[[i]] <- PCMSim(tree = ba.tree,model = true.model,X0 = true.model$X0)
  
  modelOU <- PCM(model=model.type,k=2)
  
  likFun.list[[i]] <- PCMCreateLikelihood(data.pcm[[i]], ba.tree, modelOU,metaI = PCMInfoCpp)
}

true.param <- c(alpha.a,
                alpha.b,
                b,
                sigma.a,
                sigma.b)
names.true.param <- c("H_A",
                      "H_B",
                      "Q",
                      "Sigma_A",
                      "Sigma_B")


for (i in 1:n)
{
  true.param <- c(true.param,psi.a[i],psi.b[i])
  names.true.param <- c(names.true.param,paste0(c("Psi_A","Psi_B"),"_",i))
}
names(true.param) <- names.true.param

pos.h <- grep("H",names(true.param)) 
pos.b <- grep("Q",names(true.param))
pos.sigma <- grep("Sigma",names(true.param))
pos.psi <- grep("Psi",names(true.param))

blockwise.sample.list <- list()
blockwise.sample.list[[1]] <- c(2,5)
blockwise.sample.list[[2]] <- c(3)
blockwise.sample.list[[3]] <- c(1,4)
num.block <- length(blockwise.sample.list)
for (i in 1:n)
{
  pos.theta.i <- grep(paste0("Psi_[AB]_",i,"$"),names(true.param))
  blockwise.sample.list[[i+num.block]] <- pos.theta.i
}

MyData_OUOU <- list(ll_fun = likFun.list,
                    mon.names = c("LP","LL"),
                    parm.names=names(true.param), 
                    pos.h = pos.h, 
                    pos.b = pos.b,
                    pos.sigma=pos.sigma, 
                    pos.psi=pos.psi,
                    prop = 1:length(true.param),
                    N=length(ba.tree$tip.label)*n*2)


log.true.param <- c(log(true.param[1:2]),true.param[3],log(true.param[4:5]),true.param[6:length(true.param)])
fit.RAM.model <- LaplacesDemon_local(Model_OUOU,
                                     MyData_OUOU,
                                     Initial.Values = log.true.param,
                                     Iterations = 50000,
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

save(fit.RAM.model,file="mcmc_11_species_200_genes_0.44_acceptance_50000_blockwise_lognormal_prior_ouou_predict_b_0.5_no_decompose.Rda")
pdf("/data2/cope/rna-protein-coevolution/MCMC/example_mcmc_laplace_demon_hier_blockwise_11_species_200_genes_0.44_acceptance_50000_lognormal_prior_ouou_predict_b_0.5_no_decompose.pdf")
plot(fit.RAM.model$Monitor[,"LP"],type="l",ylab="Log(Posterior)")
plot(fit.RAM.model$Monitor[,"LL"],type="l",ylab="Log(Likelihood)")
plot(fit.RAM.model$Posterior1[,"H_A"],type="l",ylab=expression(alpha["A"]))
abline(h=log.true.param[1])
plot(fit.RAM.model$Posterior1[,"H_B"],type="l",ylab=expression(alpha["B"]))
abline(h=log.true.param[2])
plot(fit.RAM.model$Posterior1[,"Q"],type="l",ylab="Q")
abline(h=log.true.param[3])
plot(fit.RAM.model$Posterior1[,"Sigma_A"],type="l",ylab=expression(sigma["A"]))
abline(h=log.true.param[4])
plot(fit.RAM.model$Posterior1[,"Sigma_B"],type="l",ylab=expression(sigma["B"]))
abline(h=log.true.param[5])

for (i in 6:ncol(fit.RAM.model$Posterior1))
{
  plot(fit.RAM.model$Posterior1[,i],type="l",ylab=expression(theta))
  abline(h=log.true.param[i])
}
dev.off()