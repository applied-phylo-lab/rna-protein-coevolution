library(PCMBase,lib.loc = "~/R_dev")
library(PCMBaseCpp)
library(LaplacesDemon,lib.loc = "~/R_dev/")
require(coda)
library(tidyverse)
library(phytools)
library(argparse)
library(extraDistr)
source("/data2/cope/rna-protein-coevolution/R_notebooks/local_functions.R")

options(PCMBase.Threshold.EV = 1e-07)
options(PCMBase.Threshold.SV = 1e-08)
# options(PCMBase.ParamValue.LowerLimit=-20)
# options(PCMBase.ParamValue.UpperLimit=20)
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
  
  hy.mean.prior <- parm[Data[["pos.hy.mean"]]]
  hy.std.prior <- exp(parm[Data[["pos.hy.std"]]])
  hx.mean.prior <- parm[Data[["pos.hx.mean"]]]
  hx.std.prior <- exp(parm[Data[["pos.hx.std"]]])

  sigmay.mean.prior <- parm[Data[["pos.sigmay.mean"]]]
  sigmay.std.prior <- exp(parm[Data[["pos.sigmay.std"]]])
  sigmax.mean.prior <- parm[Data[["pos.sigmax.mean"]]]
  sigmax.std.prior <- exp(parm[Data[["pos.sigmax.std"]]])
  
  hy <- exp(parm[Data[["pos.hy"]]])
  hx <- exp(parm[Data[["pos.hx"]]])

  sigmay <- exp(parm[Data[["pos.sigmay"]]])
  sigmax <- exp(parm[Data[["pos.sigmax"]]])
  pos.psi <- Data[["pos.psi"]]
  
  ll_fun <- Data[["ll_fun"]]
  prop <- Data[["prop"]]
  gene.ll <- Data[["Gene.LL"]]
  ## Can use pos.psi[1] as check because parm should first have hyperparameters and then gene-specific parameters
  if (prop[1] >= pos.psi[1] && !is.null(gene.ll))
  {
    gene.index <- floor((prop[1] - 8 - 1)/6) + 1
    other.ll <- sum(gene.ll[-gene.index])
    current.ll <- ll_fun[[gene.index]](c(parm[prop[1]],
                                         parm[prop[2]],
                                         parm[prop[3]],
                                         parm[prop[4]],
                                         parm[prop[1]],
                                         parm[prop[2]],
                                         parm[prop[5]],
                                         parm[prop[6]]))
    
    LL <- current.ll + other.ll
    gene.ll[gene.index] <- current.ll
    
  } else {
    num.genes <- (length(parm)-8)/6
    gene.ll <- unlist(lapply(1:num.genes, function(current.loci)
    {
      par.index <- 8 + 6*(current.loci-1) + 1
      ll <- ll_fun[[current.loci]](c(parm[par.index],
                                     parm[par.index+1],
                                     parm[par.index+2],
                                     parm[par.index+3],
                                     parm[par.index],
                                     parm[par.index+1],
                                     parm[par.index+4],
                                     parm[par.index+5]))
      ll
    })
    )
    LL <- sum(gene.ll)
  }
  
  PR <- sum(dlnorm(hy,hy.mean.prior,hy.std.prior,log=T)) +
    sum(dlnorm(hx,hx.mean.prior,hx.std.prior,log=T)) +
    sum(dlnorm(sigmay,sigmay.mean.prior,sigmay.std.prior,log=T)) +
    sum(dlnorm(sigmax,sigmax.mean.prior,sigmax.std.prior,log=T)) +
    sum(dnorm(parm[pos.psi],mean = 0,sd = 1,log=T))
  
  HYP <- sum(dhcauchy(c(sigmay.std.prior,
                        sigmax.std.prior,
                        hx.std.prior,
                        hy.std.prior),
                      sigma=1,log=T)) +
    sum(dunif(c(hy.mean.prior,
                hx.mean.prior,
                sigmay.mean.prior,
                sigmax.mean.prior),min=log(0.008),max=0,log=T))
  LP.unc <- LL + PR + HYP
  ## Apply correction for log-transformation of hyperparameters. Using lognormal for other parameters already corrects for this, but currently using uniform for hyperparameters
  LP <- LP.unc + sum(log(c(hy.std.prior,
                           hx.std.prior,
                           b.std.prior,
                           sigmay.std.prior,
                           sigmax.std.prior)))# + sum(c(hy,hx,sigmay,sigmax))
  ## Parameter.LL is a bit of misnomer, is actually log posterior for each gene (no hyperpriors)
  Modelout <- list(LP=LP,Dev=-2*LL, Monitor=c(LP,LP.unc,LL,PR,HYP), yhat=1, parm=parm,Parameter.LL=gene.ll)
  return(Modelout)
}

Model_OUOU_correlated_evolution <- function(parm, Data)
{
  
  ### Parameters
  ## Means are for a lognormal distribution, represent mean on the log-scale. Can be negative.
  hy.mean.prior <- parm[Data[["pos.hy.mean"]]]
  hy.std.prior <- exp(parm[Data[["pos.hy.std"]]])
  
  hx.mean.prior <- parm[Data[["pos.hx.mean"]]]
  hx.std.prior <- exp(parm[Data[["pos.hx.std"]]])
  
  b.mean.prior <- parm[Data[["pos.b.mean"]]]
  b.std.prior <- exp(parm[Data[["pos.b.std"]]])
  
  sigmay.mean.prior <- parm[Data[["pos.sigmay.mean"]]]
  sigmay.std.prior <- exp(parm[Data[["pos.sigmay.std"]]])
  
  sigmax.mean.prior <- parm[Data[["pos.sigmax.mean"]]]
  sigmax.std.prior <- exp(parm[Data[["pos.sigmax.std"]]])
  
  hy <- exp(parm[Data[["pos.hy"]]])
  hx <- exp(parm[Data[["pos.hx"]]])
  b <- parm[Data[["pos.b"]]]
  sigmay <- exp(parm[Data[["pos.sigmay"]]])
  sigmax <- exp(parm[Data[["pos.sigmax"]]])
  pos.psi <- Data[["pos.psi"]]
  
  ll_fun <- Data[["ll_fun"]]
  prop <- Data[["prop"]]
  gene.ll <- Data[["Gene.LL"]]
  offset <- 10
  ## Can use pos.psi[1] as check because parm should first have hyperparameters and then gene-specific parameters
  if (prop[1] >= pos.psi[1] && !is.null(gene.ll))
  {
    gene.index <- floor((prop[1] - offset - 1)/7) + 1
    other.ll <- sum(gene.ll[-gene.index])
    current.ll <- ll_fun[[gene.index]](c(parm[prop[1]] + parm[prop[5]]*parm[prop[2]],
                                         parm[prop[2]],
                                         exp(parm[prop[3]]),
                                         0,
                                         -exp(parm[prop[3]]) * parm[prop[5]],
                                         exp(parm[prop[4]]),
                                         parm[prop[1]] + parm[prop[5]]*parm[prop[2]],
                                         parm[prop[2]],
                                         exp(parm[prop[6]]),
                                         exp(parm[prop[7]])))
  
    LL <- current.ll + other.ll
    gene.ll[gene.index] <- current.ll
    
  } else {
    num.genes <- (length(parm)-offset)/7
    gene.ll <- unlist(lapply(1:num.genes, function(current.loci)
    {
      par.index <- offset + 7*(current.loci-1) + 1
      ll <- ll_fun[[current.loci]](c(parm[par.index] + parm[par.index+4]*parm[par.index+1],
                                     parm[par.index+1],
                                     exp(parm[par.index+2]),
                                     0,
                                     -exp(parm[par.index+2]) * parm[par.index+4],
                                     exp(parm[par.index+3]),
                                     parm[par.index] + parm[par.index+4]*parm[par.index+1],
                                     parm[par.index+1],
                                     exp(parm[par.index+5]),
                                     exp(parm[par.index+6])))

      ll
    })
    )
    LL <- sum(gene.ll)
  }
  PR <- sum(dlnorm(hy,hy.mean.prior,hy.std.prior,log=T)) +
    sum(dlnorm(hx,hx.mean.prior,hx.std.prior,log=T)) +
    sum(dlnorm(sigmay,sigmay.mean.prior,sigmay.std.prior,log=T)) +
    sum(dlnorm(sigmax,sigmax.mean.prior,sigmax.std.prior,log=T)) +
    sum(dnorm(b,mean = b.mean.prior,sd = b.std.prior,log=T)) +
    sum(dnorm(parm[pos.psi],mean = 0,sd = 1,log=T))
  
  HYP <- sum(dhcauchy(c(sigmay.std.prior,
                        sigmax.std.prior,
                        b.std.prior,
                        hx.std.prior,
                        hy.std.prior),
                      sigma=1,log=T)) +
    sum(dunif(c(hy.mean.prior,
                hx.mean.prior,
                sigmay.mean.prior,
                sigmax.mean.prior),min=log(0.008),max=0,log=T)) +
   dnorm(b.mean.prior,mean=0,sd=1,log=T)
  LP.unc <- LL + PR + HYP
  ## Apply correction for log-transformation of hyperparameters. Using lognormal for other parameters already corrects for this, but currently using uniform for hyperparameters
  LP <- LP.unc + sum(log(c(hy.std.prior,
                           hx.std.prior,
                           b.std.prior,
                           sigmay.std.prior,
                           sigmax.std.prior)))
  
  ## Parameter.LL is a bit of misnomer, is actually log posterior for each gene (no hyperpriors)
  Modelout <- list(LP=LP,Dev=-2*LL, Monitor=c(LP,LP.unc,LL,PR,HYP), yhat=1, parm=parm,Parameter.LL=gene.ll)
  return(Modelout)
}


dir.create(directory)
cmd <- paste(commandArgs(T),collapse = " ")
cmd <- paste("Rscript --vanilla R_scripts/runMCMC_OUOU_full_model.R",cmd,sep=" ")
readme <- paste("Model was run with bash command\n```\n",cmd,"\n```",sep="")
write(readme,file.path(directory,"README.md"))

if (!independent.evo)
{
  model.type <- "OU__Global_X0__Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x"
  #model.type <- "OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x"
  if (!omit.root)
  {
    Model_OUOU <- Model_OUOU_correlated_evolution
  } else {
    Model_OUOU <- Model_OUOU_correlated_evolution_omit_root
  }
} else {
  model.type <- "OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x"
  #model.type <- "OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x"
  if (!omit.root)
  {
    Model_OUOU <- Model_OUOU_independent_evolution 
  } else {
    Model_OUOU <- Model_OUOU_independent_evolution_omit_root
  }
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
colnames(data)[1] <- "Gene_ID"

if (no.std.err)
{
  data <- data[,c(1,2,3,4)] 
}

# the order of the traits will determine direction of "causality". First trait evolves towards optimum determined by the second trait
if (reverse.trait.order)
{
  print("Reversing order of traits...")
  # data <- data %>% 
  #   relocate(Mean_Protein,.before=Mean_RNA) %>% 
  #   relocate(SE_Protein,.before=SE_RNA)
  if (no.std.err)
  {
    data <- data[,c(1,2,4,3)]
  } else {
    data <- data[,c(1,2,4,3,6,5)]
  }
}

if (!no.std.err)
{
  
  data.se <- data[,c(1,2,5,6)] %>%
    group_by(Gene_ID) %>%
    group_split()
  
  data.pcmbase.format.se <- purrr::map(data.se,function(x){
    x %>% 
      dplyr::select(-Gene_ID) %>%
      dplyr::slice(match(tree$tip.label,Species)) %>%
      column_to_rownames("Species") %>%
      t()
  }
  )
} else{
  data.se <- NULL
  data.pcmbase.format.se <- NULL
}

data.gene <- data[,c(1,2,3,4)] %>%
  group_by(Gene_ID) %>%
  group_split()


## Format data to work with PCMBase
data.pcmbase.format.trait <- purrr::map(data.gene,function(x){
  x %>% 
    dplyr::select(-Gene_ID) %>%
    dplyr::slice(match(tree$tip.label,Species)) %>%
    column_to_rownames("Species") %>%
    t()
}
)


n <- length(data.pcmbase.format.trait)

# Create likelihood function to speed up computation of likelihood                  
likFun.list <- vector(mode="list",length=n)
for (i in 1:n)
{
  modelOU <- PCM(model= model.type,k=2)
  if (!no.std.err)
  {
    likFun.list[[i]] <- PCMCreateLikelihood(data.pcmbase.format.trait[[i]], tree,SE=data.pcmbase.format.se[[i]],modelOU,metaI = PCMInfoCpp)
  } else {
    likFun.list[[i]] <- PCMCreateLikelihood(data.pcmbase.format.trait[[i]], tree,modelOU,metaI = PCMInfoCpp)
    
  }
}

tip.heights <- ips::tipHeights(tree)[1]
start.alpha <- log(2)/(tip.heights/4)

if (!independent.evo)
{
  # Get starting values. For now, use the mean of the prior for H and Sigma, 0 for Q (no affect of X on Y), and mean of each trait for optimum
  # Root state will be determined by optimum values for Y and X, as well as Q. 
  #start.values <- c(log(0.25),log(0.25),0,log(0.25),log(0.25)) # will propose H, Sigma on log-scale
  start.values <- c(log(start.alpha),log(1.5),log(start.alpha),log(1.5),1,log(0.25),log(start.alpha),log(1.5),log(start.alpha),log(1.5)) # will propose H, Sigma on log-scale
  names(start.values) <- c("H_Y_Mean","H_Y_Std","H_X_Mean","H_X_Std","Q_Mean","Q_Std","Sigma_Y_Mean","Sigma_Y_Std","Sigma_X_Mean","Sigma_X_Std")
  optimum.values <- unlist(
    lapply(1:length(data.pcmbase.format.trait), function(i)
    {
      tmp <- data.pcmbase.format.trait[[i]]
      trait.mean <- c(rowMeans(tmp),log(0.25),log(0.25),0,log(0.25),log(0.25))
      names(trait.mean) <- paste0(c("Psi_Y_","Psi_X_","H_Y_","H_X_","Q_","Sigma_Y_","Sigma_X_"),i)
      return(trait.mean)
    }
    
    )
  )
  
  start.values <- c(start.values,optimum.values)
  
  pos.hy.mean <- grep("H_Y_Mean",names(start.values)) 
  pos.hy.std <- grep("H_Y_Std",names(start.values)) 
  pos.hx.mean <- grep("H_X_Mean",names(start.values)) 
  pos.hx.std <- grep("H_X_Std",names(start.values)) 
  
  pos.b.mean <- grep("Q_Mean",names(start.values))
  pos.b.std <- grep("Q_Std",names(start.values))
  
  pos.sigmay.mean <- grep("Sigma_Y_Mean",names(start.values)) 
  pos.sigmay.std <- grep("Sigma_Y_Std",names(start.values)) 
  pos.sigmax.mean <- grep("Sigma_X_Mean",names(start.values)) 
  pos.sigmax.std <- grep("Sigma_X_Std",names(start.values)) 
  
  
  pos.hx <- grep("H_X_[0-9]+",names(start.values)) 
  pos.hy <- grep("H_Y_[0-9]+",names(start.values)) 
  pos.b <- grep("Q_[0-9]+",names(start.values))
  pos.sigmax <- grep("Sigma_X_[0-9]+",names(start.values))
  pos.sigmay <- grep("Sigma_Y_[0-9]+",names(start.values))
  
  pos.psi <- grep("Psi",names(start.values))
  
  MyData_OUOU <- list(ll_fun = likFun.list,
                      mon.names = c("LP","LP.unc","LL","PR","HYP"),
                      parm.names=names(start.values), 
                      pos.hy.mean = pos.hy.mean,
                      pos.hy.std = pos.hy.std,
                      pos.hx.mean = pos.hx.mean,
                      pos.hx.std = pos.hx.std,
                      pos.b.mean = pos.b.mean,
                      pos.b.std = pos.b.std,
                      pos.sigmay.mean = pos.sigmay.mean,
                      pos.sigmay.std = pos.sigmay.std, 
                      pos.sigmax.mean = pos.sigmax.mean,
                      pos.sigmax.std = pos.sigmax.std,
                      pos.hy = pos.hy, 
                      pos.hx = pos.hx, 
                      pos.b = pos.b,
                      pos.sigmax=pos.sigmax,
                      pos.sigmay=pos.sigmay,
                      pos.psi=pos.psi,
                      prop = 1:length(start.values),
                      N=length(tree$tip.label)*n*2)
  
  
  
  blockwise.sample.list <- list()
  blockwise.sample.list[[1]] <- c(3,4,9,10)
  blockwise.sample.list[[2]] <- c(5,6)
  blockwise.sample.list[[3]] <- c(1,2,7,8)
  num.block <- length(blockwise.sample.list)
  for (i in 1:n)
  {
    pos.psi.i <- grep(paste0("Psi_[YX]_",i,"$"),names(start.values))
    pos.h.i <- grep(paste0("H_[YX]_",i,"$"),names(start.values))
    pos.q.i <- grep(paste0("Q_",i,"$"),names(start.values))
    pos.sigma.i <- grep(paste0("Sigma_[YX]_",i,"$"),names(start.values))
    blockwise.sample.list[[i+num.block]] <- c(pos.psi.i,pos.h.i,pos.q.i,pos.sigma.i)
  }
  
} else {
  # Get starting values. For now, use the mean of the prior for H and Sigma, 0 for Q (no affect of X on Y), and mean of each trait for optimum
  # Root state will be determined by optimum values for Y and X, as well as Q. 
  #start.values <- c(log(0.25),log(0.25),0,log(0.25),log(0.25)) # will propose H, Sigma on log-scale
  start.values <- c(log(start.alpha),log(1.5),log(start.alpha),log(1.5),log(start.alpha),log(1.5),log(start.alpha),log(1.5)) # will propose H, Sigma on log-scale
  names(start.values) <- c("H_Y_Mean","H_Y_Std","H_X_Mean","H_X_Std","Sigma_Y_Mean","Sigma_Y_Std","Sigma_X_Mean","Sigma_X_Std")
  optimum.values <- unlist(
    lapply(1:length(data.pcmbase.format.trait), function(i)
    {
      tmp <- data.pcmbase.format.trait[[i]]
      trait.mean <- c(rowMeans(tmp),log(0.25),log(0.25),log(0.25),log(0.25))
      names(trait.mean) <- paste0(c("Psi_Y_","Psi_X_","H_Y_","H_X_","Sigma_Y_","Sigma_X_"),i)
      return(trait.mean)
    }
    
    )
  )
  
  start.values <- c(start.values,optimum.values)
  pos.hy.mean <- grep("H_Y_Mean",names(start.values)) 
  pos.hy.std <- grep("H_Y_Std",names(start.values)) 
  pos.hx.mean <- grep("H_X_Mean",names(start.values)) 
  pos.hx.std <- grep("H_X_Std",names(start.values)) 
 
  
  pos.sigmay.mean <- grep("Sigma_Y_Mean",names(start.values)) 
  pos.sigmay.std <- grep("Sigma_Y_Std",names(start.values)) 
  pos.sigmax.mean <- grep("Sigma_X_Mean",names(start.values)) 
  pos.sigmax.std <- grep("Sigma_X_Std",names(start.values)) 
  
  
  pos.hx <- grep("H_X_[0-9]+",names(start.values)) 
  pos.hy <- grep("H_Y_[0-9]+",names(start.values)) 
  pos.sigmax <- grep("Sigma_X_[0-9]+",names(start.values))
  pos.sigmay <- grep("Sigma_Y_[0-9]+",names(start.values))
  
  pos.psi <- grep("Psi",names(start.values))
  
  MyData_OUOU <- list(ll_fun = likFun.list,
                      mon.names = c("LP","LP.unc","LL","PR","HYP"),
                      parm.names=names(start.values), 
                      pos.hy.mean = pos.hy.mean,
                      pos.hy.std = pos.hy.std,
                      pos.hx.mean = pos.hx.mean,
                      pos.hx.std = pos.hx.std,
                      pos.sigmay.mean = pos.sigmay.mean,
                      pos.sigmay.std = pos.sigmay.std, 
                      pos.sigmax.mean = pos.sigmax.mean,
                      pos.sigmax.std = pos.sigmax.std,
                      pos.hy = pos.hy, 
                      pos.hx = pos.hx, 
                      pos.sigma=pos.sigma,
                      pos.sigmay=pos.sigmay,
                      pos.psi=pos.psi,
                      prop = 1:length(start.values),
                      N=length(tree$tip.label)*n*2)
  
  
  
  blockwise.sample.list <- list()
  blockwise.sample.list[[1]] <- c(3,4,7,8)
  blockwise.sample.list[[2]] <- c(1,2,5,6)
  num.block <- length(blockwise.sample.list)
  for (i in 1:n)
  {
    pos.psi.i <- grep(paste0("Psi_[YX]_",i,"$"),names(start.values))
    pos.h.i <- grep(paste0("H_[YX]_",i,"$"),names(start.values))
    pos.sigma.i <- grep(paste0("Sigma_[YX]_",i,"$"),names(start.values))
    blockwise.sample.list[[i+num.block]] <- c(pos.psi.i,pos.h.i,pos.sigma.i)
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
  plot(adapt.fit$Posterior1[range.to.plot,"H_Y_Mean"],type="l",ylab="Estimate",main="H_Y_Mean")
  plot(adapt.fit$Posterior1[range.to.plot,"H_Y_Std"],type="l",ylab="Estimate",main="H_Y_Std")
  plot(adapt.fit$Posterior1[range.to.plot,"H_X_Mean"],type="l",ylab="Estimate",main="H_X_Mean")
  plot(adapt.fit$Posterior1[range.to.plot,"H_X_Std"],type="l",ylab="Estimate",main="H_X_Std")
  if (!independent.evo)
  {
    plot(adapt.fit$Posterior1[range.to.plot,"Q_Mean"],type="l",ylab="Estimate",main="Q_Mean")
    plot(adapt.fit$Posterior1[range.to.plot,"Q_Std"],type="l",ylab="Estimate",main="Q_Std")
  }
  plot(adapt.fit$Posterior1[range.to.plot,"Sigma_Y_Mean"],type="l",ylab="Estimate",main="Sigma_Y_Mean")
  plot(adapt.fit$Posterior1[range.to.plot,"Sigma_Y_Std"],type="l",ylab="Estimate",main="Sigma_Y_Std")
  plot(adapt.fit$Posterior1[range.to.plot,"Sigma_X_Mean"],type="l",ylab="Estimate",main="Sigma_X_Mean")
  plot(adapt.fit$Posterior1[range.to.plot,"Sigma_X_Std"],type="l",ylab="Estimate",main="Sigma_X_Std")
  
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


# print("Beginning Random-Walk Metropolis MCMC for parameter estimation...")
# fit <- LaplacesDemon_local(Model_OUOU,
#                            MyData_OUOU,
#                            Initial.Values = adapt.fit$Posterior1[num.prev.sample,],
#                            Covar=adapt.fit$Covar,
#                            Iterations = samples,
#                            Algorithm = "RWM",
#                            Thinning = thinning,
#                            Debug = list(DB.Model = F,
#                                         DB.chol = T),
#                            Specs = list(B = blockwise.sample.list)
# )
# 
# save(fit,file=file.path(directory,paste0("rwm_mcmc_",rwm.run.number,".Rda")))
# 
# if (!is.null(prev.rwm.run))
# {
#   fit <- Combine(list(adapt.fit,fit),Data=MyData_OUOU)
#   save(fit,file=file.path(directory,"rwm_mcmc_total.Rda"))
# }
# 
# 
# burnin<-1
# nsample <- nrow(fit$Monitor)
# range.to.plot <- burnin:nsample
# pdf(file.path(directory,"mcmc_traces.pdf"))
# plot(fit$Monitor[range.to.plot,"LP"],type="l",ylab="Log(Posterior)")
# plot(fit$Monitor[range.to.plot,"LL"],type="l",ylab="Log(Likelihood)")
# plot(fit$Posterior1[range.to.plot,"H_Y_Mean"],type="l",ylab="Estimate",main="H_Y_Mean")
# plot(fit$Posterior1[range.to.plot,"H_Y_Std"],type="l",ylab="Estimate",main="H_Y_Std")
# plot(fit$Posterior1[range.to.plot,"H_X_Mean"],type="l",ylab="Estimate",main="H_X_Mean")
# plot(fit$Posterior1[range.to.plot,"H_X_Std"],type="l",ylab="Estimate",main="H_X_Std")
# if (!independent.evo)
# {
#   plot(fit$Posterior1[range.to.plot,"Q_Mean"],type="l",ylab="Estimate",main="Q_Mean")
#   plot(fit$Posterior1[range.to.plot,"Q_Std"],type="l",ylab="Estimate",main="Q_Std")
# }
# plot(fit$Posterior1[range.to.plot,"Sigma_Y_Mean"],type="l",ylab="Estimate",main="Sigma_Y_Mean")
# plot(fit$Posterior1[range.to.plot,"Sigma_Y_Std"],type="l",ylab="Estimate",main="Sigma_Y_Std")
# plot(fit$Posterior1[range.to.plot,"Sigma_X_Mean"],type="l",ylab="Estimate",main="Sigma_X_Mean")
# plot(fit$Posterior1[range.to.plot,"Sigma_X_Std"],type="l",ylab="Estimate",main="Sigma_X_Std")
# 
# dev.off()
# 
# 
