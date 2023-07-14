library(PCMBase,lib.loc = "~/R_dev")
library(PCMBaseCpp)
library(LaplacesDemon,lib.loc = "~/R_dev/")
require(coda)
library(tidyverse)
library(phytools)
library(argparse)
library(bayestestR)
source("/data2/cope/rna-protein-coevolution/R_notebooks/local_functions.R")

options(PCMBase.Threshold.EV = 1e-07)
options(PCMBase.Threshold.SV = 1e-08)
parser <- ArgumentParser()
parser$add_argument("-i","--input",help="tsv file with data input",type="character",default="./")
parser$add_argument("-o","--output",help="Directory of where to put results. Will automatically generate lowest-level directory if not already generated.",type="character",default="./")
parser$add_argument("-s","--samples",help="Number of samples for parameter estimation (RWM)",type="integer",default=50000)
parser$add_argument("--adapt_samples",help="Number of samples to perform adapting (RAM)",type="integer",default=10000)
parser$add_argument("-t","--thin",help="Thinning value. Total number of iterations will be samples * thinning",type="integer",default=5)
parser$add_argument("--alpha_star",type="double",default=0.44)
parser$add_argument("--tree",type="character",default="../Data/tree_11sp_noGpig.nex")
parser$add_argument("--prev_adapt_run",type="character",default=NULL)
parser$add_argument("--prev_rwm_run",type="character",default=NULL)
parser$add_argument("--rwm_run_number",type="integer",default=1)
parser$add_argument("--newick",action="store_true")
parser$add_argument("--independent_evo",action="store_true")
parser$add_argument("--reverse_trait_order",action="store_true")
parser$add_argument("--omit_root",action="store_true")
parser$add_argument("--no_std_err",action="store_true")


args <- parser$parse_args()
print(args)
input <- args$input
directory <- args$output
thinning <- args$thin
samples <- args$samples
alpha.star <- args$alpha_star
tree.file <- args$tree
adapt.samples <- args$adapt_samples
newick <- args$newick
independent.evo <- args$independent_evo
reverse.trait.order <- args$reverse_trait_order
omit.root <- args$omit_root
no.std.err <- args$no_std_err
prev.adapt.run <- args$prev_adapt_run
prev.rwm.run <- args$prev_rwm_run
rwm.run.number <- args$rwm_run_number


Model_OUOU_independent_evolution <- function(parm, Data)
{
  ### Parameters
  pos.psi <- Data[["pos.psi"]]
  h <- exp(parm[Data[["pos.h"]]])
  sigma <- exp(parm[Data[["pos.sigma"]]])
  ll_fun <- Data[["ll_fun"]]
  prop <- Data[["prop"]]
  gene.ll <- Data[["Gene.LL"]]
  if (prop[1] >= pos.psi[1] && !is.null(gene.ll))
  {
    gene.index <- floor((prop[1] - 4 - 1)/2) + 1
    other.ll <- sum(gene.ll[-gene.index])
    current.ll <- ll_fun[[gene.index]](c(parm[prop[1]],
                                         parm[prop[2]],
                                         h[1],
                                         h[2],
                                         parm[prop[1]],
                                         parm[prop[2]],
                                         sigma))
    LL <- current.ll + other.ll
    gene.ll[gene.index] <- current.ll
    
  } else {
    num.genes <- (length(parm)-4)/2
    gene.ll <- unlist(lapply(1:num.genes, function(current.loci)
    {
      par.index <- 4 + 2*(current.loci-1) + 1
      ll <- ll_fun[[current.loci]](c(parm[par.index] + b*parm[par.index+1],
                                     parm[par.index+1],
                                     h[1],
                                     h[2],
                                     parm[par.index] + b*parm[par.index+1],
                                     parm[par.index+1],
                                     sigma))
      ll
    })
    )
    LL <- sum(gene.ll)
  }
  hyper.pr <-  sum(dlnorm(c(h,sigma),meanlog = 0.25,sdlog = 1.5,log=T))
  pr <- sum(dnorm(parm[pos.psi[1]:length(parm)],mean = 0,sd = 1,log=T))
  LP <- LL + pr + hyper.pr
  Modelout <- list(LP=LP,Dev=-2*LL, Monitor=c(LP,LL), yhat=1, parm=parm,Parameter.LL=gene.ll)
  return(Modelout)
}
  
Model_OUOU_correlated_evolution <- function(parm, Data)
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



dir.create(directory)
cmd <- paste(commandArgs(T),collapse = " ")
cmd <- paste("Rscript --vanilla R_scripts/runMCMC_OUOU.R",cmd,sep=" ")
readme <- paste("Model was run with bash command\n```\n",cmd,"\n```",sep="")
write(readme,file.path(directory,"README.md"))

if (!independent.evo)
{
  model.type <- "OU__Global_X0__Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x"
  Model_OUOU <- Model_OUOU_correlated_evolution
} else {
  model.type <- "OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x"
  Model_OUOU <- Model_OUOU_independent_evolution
}

if (newick)
{
  tree <- read.tree(tree.file)
} else {
  tree <- read.nexus(tree.file)
}
if (!is.ultrametric(tree))
{
  tree <- force.ultrametric(tree,method = "extend")
}

data <- read_tsv(input)

if (reverse.trait.order)
{
  print("Reversing order of traits...")
  data <- data %>% 
    relocate(V2,.before=V1)
}


data.gene <- data %>%
  group_by(Locus) %>%
  group_split()


## Format data to work with PCMBase
data.pcmbase.format <- purrr::map(data.gene,function(x){
                                x %>% 
                                  dplyr::select(-Locus) %>%
                                  dplyr::slice(match(tree$tip.label,Species)) %>%
                                  column_to_rownames("Species") %>%
                                  t()
      }
)
n <- length(data.pcmbase.format)


# Create likelihood function to speed up computation of likelihood                  
likFun.list <- vector(mode="list",length=n)
for (i in 1:n)
{
  modelOU <- PCM(model= model.type,k=2)
  likFun.list[[i]] <- PCMCreateLikelihood(data.pcmbase.format[[i]], tree, modelOU,metaI = PCMInfoCpp)
}




if (!independent.evo)
{
  # Get starting values. For now, use the mean of the prior for H and Sigma, 0 for Q (no affect of X on Y), and mean of each trait for optimum
  # Root state will be determined by optimum values for Y and X, as well as Q. 
  start.values <- c(log(0.25),log(0.25),0,log(0.25),log(0.25)) # will propose H, Sigma on log-scale
  names(start.values) <- c("H_Y","H_X","Q","Sigma_Y","Sigma_X")
  optimum.values <- unlist(
                       lapply(1:length(data.pcmbase.format), function(i)
                       {
                         tmp <- data.pcmbase.format[[i]]
                         trait.mean <- rowMeans(tmp)
                         names(trait.mean) <- paste0(c("Psi_Y_","Psi_X_"),i)
                         return(trait.mean)
                       }
  
        )
  )
  
  start.values <- c(start.values,optimum.values)
  
  pos.h <- grep("H",names(start.values)) 
  pos.b <- grep("Q",names(start.values))
  pos.sigma <- grep("Sigma",names(start.values))
  pos.psi <- grep("Psi",names(start.values))
  
  MyData_OUOU <- list(ll_fun = likFun.list,
                      mon.names = c("LP","LL"),
                      parm.names=names(start.values), 
                      pos.h = pos.h, 
                      pos.b = pos.b,
                      pos.sigma=pos.sigma, 
                      pos.psi=pos.psi,
                      prop = 1:length(start.values),
                      N=length(tree$tip.label)*n*2)
  
  
  
  blockwise.sample.list <- list()
  blockwise.sample.list[[1]] <- c(2,5)
  blockwise.sample.list[[2]] <- c(3)
  blockwise.sample.list[[3]] <- c(1,4)
  num.block <- length(blockwise.sample.list)
  for (i in 1:n)
  {
    pos.psi.i <- grep(paste0("Psi_[YX]_",i,"$"),names(start.values))
    blockwise.sample.list[[i+num.block]] <- pos.psi.i
  }
} else {
  # Get starting values. For now, use the mean of the prior for H and Sigma. Assuming independent evolution
  start.values <- c(log(0.25),log(0.25),log(0.25),log(0.25)) # will propose H, Sigma on log-scale
  names(start.values) <- c("H_Y","H_X","Sigma_Y","Sigma_X")
  optimum.values <- unlist(
    lapply(1:length(data.pcmbase.format), function(i)
    {
      tmp <- data.pcmbase.format[[i]]
      trait.mean <- rowMeans(tmp)
      names(trait.mean) <- paste0(c("Psi_Y_","Psi_X_"),i)
      return(trait.mean)
    }
    
    )
  )
  
  start.values <- c(start.values,optimum.values)
  
  pos.h <- grep("H",names(start.values)) 
  pos.sigma <- grep("Sigma",names(start.values))
  pos.psi <- grep("Psi",names(start.values))
  
  MyData_OUOU <- list(ll_fun = likFun.list,
                      mon.names = c("LP","LL"),
                      parm.names=names(start.values), 
                      pos.h = pos.h, 
                      pos.sigma=pos.sigma, 
                      pos.psi=pos.psi,
                      prop = 1:length(start.values),
                      N=length(tree$tip.label)*n*2)
  
  
  
  blockwise.sample.list <- list()
  blockwise.sample.list[[1]] <- c(2,5)
  blockwise.sample.list[[2]] <- c(1,4)
  num.block <- length(blockwise.sample.list)
  for (i in 1:n)
  {
    pos.psi.i <- grep(paste0("Psi_[YX]_",i,"$"),names(start.values))
    blockwise.sample.list[[i+num.block]] <- pos.psi.i
  }
}


if (adapt.samples > 0 && is.null(prev.adapt.run) && is.null(prev.rwm.run))
{
  print("Beginning adaptive MCMC...")
  adapt.fit <- LaplacesDemon_local(Model_OUOU,
                                   MyData_OUOU,
                                   Initial.Values = start.values,
                                   Iterations = adapt.samples,
                                   Algorithm = "RAM",
                                   Thinning = thinning,
                                   Debug = list(DB.Model = F,
                                                DB.chol = T),
                                   Specs = list(alpha.star = alpha.star,
                                                B = blockwise.sample.list,
                                                Dist="t",
                                                gamma = 2/3,
                                                n=0
                                   ))
  
  save(adapt.fit,file=file.path(directory,"adaptive_mcmc.Rda"))
  
  burnin<-1
  nsample <- nrow(adapt.fit$Monitor)
  range.to.plot <- burnin:nsample
  pdf(file.path(directory,"adaptive_mcmc_traces.pdf"))
  plot(adapt.fit$Monitor[range.to.plot,"LP"],type="l",ylab="Log(Posterior)")
  plot(adapt.fit$Monitor[range.to.plot,"LL"],type="l",ylab="Log(Likelihood)")
  plot(adapt.fit$Posterior1[range.to.plot,"H_Y"],type="l",ylab=expression(alpha["Y"]))
  plot(adapt.fit$Posterior1[range.to.plot,"H_X"],type="l",ylab=expression(alpha["X"]))
  if (!independent.evo)
  {
    plot(adapt.fit$Posterior1[range.to.plot,"Q"],type="l",ylab="Q")
  }
  plot(adapt.fit$Posterior1[range.to.plot,"Sigma_Y"],type="l",ylab=expression(sigma["Y"]))
  plot(adapt.fit$Posterior1[range.to.plot,"Sigma_X"],type="l",ylab=expression(sigma["X"]))
  
  dev.off()
} else if(!is.null(prev.adapt.run)) {
  print("Loading previous adaptive run")
  load(prev.adapt.run)
} else if(!is.null(prev.rwm.run)) {
  print("Loading previous RWM run")
  load(prev.rwm.run)
  adapt.fit <- fit
}
num.prev.sample <- nrow(adapt.fit$Posterior1)


# num.prev.sample <- nrow(adapt.fit$Posterior1)
# fit <- LaplacesDemon_local(Model_OUOU,
#                                  MyData_OUOU,
#                                  Initial.Values = adapt.fit$Posterior1[num.prev.sample,],
#                                  Covar=adapt.fit$Covar,
#                                  Iterations = samples,
#                                  Algorithm = "RWM",
#                                  Thinning = thinning,
#                                  Debug = list(DB.Model = F,
#                                               DB.chol = T),
#                                  Specs = list(B = blockwise.sample.list)
#                            )
# 
# save(fit,file=file.path(directory,"rwm_mcmc.Rda"))
# 
# 
# burnin<-1
# nsample <- nrow(fit$Monitor)
# range.to.plot <- burnin:nsample
# pdf(file.path(directory,"mcmc_traces.pdf"))
# plot(fit$Monitor[range.to.plot,"LP"],type="l",ylab="Log(Posterior)")
# plot(fit$Monitor[range.to.plot,"LL"],type="l",ylab="Log(Likelihood)")
# plot(fit$Posterior1[range.to.plot,"H_Y"],type="l",ylab=expression(alpha["Y"]))
# plot(fit$Posterior1[range.to.plot,"H_X"],type="l",ylab=expression(alpha["X"]))
# if (!independent.evo)
# {
#   plot(fit$Posterior1[range.to.plot,"Q"],type="l",ylab="Q")
# }
# plot(fit$Posterior1[range.to.plot,"Sigma_Y"],type="l",ylab=expression(sigma["Y"]))
# plot(fit$Posterior1[range.to.plot,"Sigma_X"],type="l",ylab=expression(sigma["X"]))
# dev.off()
# 
# 
