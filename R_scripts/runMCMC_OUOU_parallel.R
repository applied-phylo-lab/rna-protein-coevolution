library(PCMBase,lib.loc = "~/R_dev")
library(PCMBaseCpp)
library(LaplacesDemon,lib.loc = "~/R_dev/")
library(tidyverse)
library(phytools)
library(argparse)
library(parallel)
source("/data2/cope/Active_projects/Cope_Schraiber_Pennell/rna-protein-coevolution/R_notebooks/local_functions.R")


parser <- ArgumentParser()
parser$add_argument("-i","--input",help="tsv file with data input",type="character",default="./")
parser$add_argument("-o","--output",help="Directory of where to put results. Will automatically generate lowest-level directory if not already generated.",type="character",default="./")
parser$add_argument("-s","--samples",help="Number of samples for parameter estimation (RWM)",type="integer",default=50000)
parser$add_argument("--burnin_samples",help="Number of samples to perform RWM before starting adapting. May not be necessary for RAM.",type="integer",default=1000)
parser$add_argument("--adapt_samples",help="Number of samples to perform adapting (RAM)",type="integer",default=10000)
parser$add_argument("-t","--thin",help="Thinning value. Total number of iterations will be samples * thinning",type="integer",default=5)
parser$add_argument("--alpha_star",type="double",default=0.44)
parser$add_argument("--tree",type="character",default="../Data/tree_11sp_noGpig.nex")
parser$add_argument("--prev_adapt_run",type="character",default=NULL)
parser$add_argument("--prev_rwm_run",type="character",default=NULL)
parser$add_argument("--rwm_run_number",type="integer",default=1)
parser$add_argument("--start_values_file",type="character",default=NULL,help="A two column file specifying starting values for model parameters. ")
parser$add_argument("--randomize_psi_start_sd",type="double",default=0,help="Psi values, by default, start at the average values of the traits across species. This adds random variation N(0,sd) around the traits, where sd is the value specified here. Default = 0 (no random variation)")
parser$add_argument("--newick",action="store_true")
parser$add_argument("--independent_evo",action="store_true",help="Traits are assumed to be evolving independently, i.e. all off-diagonals of matrices are 0")
parser$add_argument("--reverse_trait_order",action="store_true",help="Reverse the order of the traits in the data matrices. This determines the response and predictor trait.")
parser$add_argument("--sigma_off_diagonal",action="store_true")
parser$add_argument("--omit_root",action="store_true")
parser$add_argument("--no_std_err",action="store_true",help="Ignore standard error estimates.")
parser$add_argument("--num_cores",type="integer",default=12)
parser$add_argument("--PCM_options_file",type="character",default=NULL,help="R file that will be sourced. Should contain commands for setting PCM options. Please provide the absolute path or path relative to this file. Default=NULL.")
parser$add_argument("--drop_species",type="character",default=NULL,help="A semi-colon delimited list with species to drop from data and tree")


args <- parser$parse_args()
print(args)
input <- args$input
directory <- args$output
thinning <- args$thin
samples <- args$samples
alpha.star <- args$alpha_star
tree.file <- args$tree
burnin.samples <- args$burnin_samples
adapt.samples <- args$adapt_samples
newick <- args$newick
independent.evo <- args$independent_evo
reverse.trait.order <- args$reverse_trait_order
sigma.off.diagonal <- args$sigma_off_diagonal
omit.root <- args$omit_root
prev.adapt.run <- args$prev_adapt_run
prev.rwm.run <- args$prev_rwm_run
rwm.run.number <- args$rwm_run_number
no.std.err <- args$no_std_err
num.cores <- args$num_cores
start.values.file <- args$start_values_file
PCM.options.file <- args$PCM_options_file
psi.random.sd <- args$randomize_psi_start_sd
drop.species <- args$drop_species

if(!is.null(PCM.options.file))
{
  source(PCM.options.file)
}

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
    gene.ll <- unlist(mclapply(1:num.genes, function(current.loci)
    {
      par.index <- 4 + 2*(current.loci-1) + 1
      ll <- ll_fun[[current.loci]](c(parm[par.index],
                                     parm[par.index+1],
                                     h[1],
                                     h[2],
                                     parm[par.index],
                                     parm[par.index+1],
                                     sigma))
      ll
    },mc.cores=Data$num.cores)
    )
    LL <- sum(gene.ll)
  }
  HYP <-  sum(dlnorm(h,meanlog = log(log(2)/(173/4)),sdlog = 1.5,log=T)) + 
    sum(dlnorm(sigma,meanlog = log(0.25),sdlog = 1.5,log=T))
  PR <- sum(dnorm(parm[pos.psi[1]:length(parm)],mean = 0,sd = 1.0,log=T))
  LP.unc <- LL + PR + HYP
  LP <- LP.unc + sum(log(c(h,sigma)))
  Modelout <- list(LP=LP,Dev=-2*LL, Monitor=c(LP,LP.unc,LL,PR,HYP), yhat=1, parm=parm,Parameter.LL=gene.ll)
  return(Modelout)
}

Model_OUOU_independent_evolution_omit_root <- function(parm, Data)
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
    current.ll <- ll_fun[[gene.index]](c(NA,
                                         NA,
                                         h[1],
                                         h[2],
                                         parm[prop[1]],
                                         parm[prop[2]],
                                         sigma))
    LL <- current.ll + other.ll
    gene.ll[gene.index] <- current.ll
    
  } else {
    num.genes <- (length(parm)-4)/2
    gene.ll <- unlist(mclapply(1:num.genes, function(current.loci)
    {
      par.index <- 4 + 2*(current.loci-1) + 1
      ll <- ll_fun[[current.loci]](c(NA,
                                     NA,
                                     h[1],
                                     h[2],
                                     parm[par.index],
                                     parm[par.index+1],
                                     sigma))
      ll
    },mc.cores=Data$num.cores)
    )
    LL <- sum(gene.ll)
  }
  HYP <-  sum(dlnorm(h,meanlog = log(log(2)/(173/4)),sdlog = 1.5,log=T)) + 
    sum(dlnorm(sigma,meanlog = log(0.25),sdlog = 1.5,log=T))
  PR <- sum(dnorm(parm[pos.psi[1]:length(parm)],mean = 0,sd = 1,log=T))
  LP.unc <- LL + PR + HYP
  LP <- LP.unc + sum(log(c(h,sigma)))
  Modelout <- list(LP=LP,Dev=-2*LL, Monitor=c(LP,LP.unc,LL,PR,HYP), yhat=1, parm=parm,Parameter.LL=gene.ll)
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
    gene.ll <- unlist(mclapply(1:num.genes, function(current.loci)
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
    },mc.cores=Data$num.cores)
    )
    LL <- sum(gene.ll)
  }
  HYP <-  sum(dlnorm(h,meanlog = log(log(2)/(173/4)),sdlog = 1.5,log=T)) +
    sum(dlnorm(sigma,meanlog = log(0.25),sdlog = 1.5,log=T)) +
    dnorm(b,mean = 0,sd = 1,log=T)
  PR <- sum(dnorm(parm[pos.psi[1]:length(parm)],mean = 0,sd = 1,log=T))
  LP.unc <- LL + PR + HYP
  LP <- LP.unc + sum(log(c(h,sigma)))
  Modelout <- list(LP=LP,Dev=-2*LL, Monitor=c(LP,LP.unc,LL,PR,HYP), yhat=1, parm=parm,Parameter.LL=gene.ll)
  return(Modelout)
}

Model_OUOU_correlated_evolution_omit_root <- function(parm, Data)
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
    current.ll <- ll_fun[[gene.index]](c(NA,
                                         NA,
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
    gene.ll <- unlist(mclapply(1:num.genes, function(current.loci)
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
    },mc.cores=Data$num.cores)
    )
    LL <- sum(gene.ll)
  }
  HYP <-  sum(dlnorm(h,meanlog = log(log(2)/(173/4)),sdlog = 1.5,log=T)) + 
    sum(dlnorm(sigma,meanlog = log(0.25),sdlog = 1.5,log=T)) + 
    dnorm(b,mean = 0,sd = 1,log=T)
  PR <- sum(dnorm(parm[pos.psi[1]:length(parm)],mean = 0,sd = 1,log=T))
  LP.unc <- LL + PR + HYP
  LP <- LP.unc + sum(log(c(h,sigma)))
  Modelout <- list(LP=LP,Dev=-2*LL, Monitor=c(LP,LP.unc,LL,PR,HYP), yhat=1, parm=parm,Parameter.LL=gene.ll)
  return(Modelout)
}




Model_OUOU_correlated_evolution_sigma_off <- function(parm, Data)
{
  ### Parameters
  pos.psi <- Data[["pos.psi"]]
  h <- exp(parm[Data[["pos.h"]]])
  b <- parm[Data[["pos.b"]]]
  c.val <- parm[Data[["pos.c"]]]
  sigma <- exp(parm[Data[["pos.sigma"]]])
  ll_fun <- Data[["ll_fun"]]
  prop <- Data[["prop"]]
  gene.ll <- Data[["Gene.LL"]]
  if (prop[1] >= pos.psi[1] && !is.null(gene.ll))
  {
    gene.index <- floor((prop[1] - 6 - 1)/2) + 1
    other.ll <- sum(gene.ll[-gene.index])
    current.ll <- ll_fun[[gene.index]](c(parm[prop[1]] + b*parm[prop[2]],
                                         parm[prop[2]],
                                         h[1],
                                         0,
                                         h[2] * c.val - h[1] * b,
                                         h[2],
                                         parm[prop[1]] + b*parm[prop[2]],
                                         parm[prop[2]],
                                         sigma[1],
                                         c.val*sigma[2],
                                         sigma[2]))
    LL <- current.ll + other.ll
    gene.ll[gene.index] <- current.ll
    
  } else {
    num.genes <- (length(parm)-6)/2
    gene.ll <- unlist(mclapply(1:num.genes, function(current.loci)
    {
      par.index <- 6 + 2*(current.loci-1) + 1
      ll <- ll_fun[[current.loci]](c(parm[par.index] +  b * parm[par.index+1],
                                     parm[par.index+1],
                                     h[1],
                                     0,
                                     h[2] * c.val - h[1] * b, # c should be multiplying alpha of predictor, b should be multiplying alpha of response
                                     h[2],
                                     parm[par.index] + b*parm[par.index+1],
                                     parm[par.index+1],
                                     sigma[1],
                                     c.val*sigma[2], # c should be multiplying sigma of predictor
                                     sigma[2]))
      ll
    },mc.cores=Data$num.cores)
    )
    LL <- sum(gene.ll)
  }
  HYP <-  sum(dlnorm(h,meanlog = log(log(2)/(173/4)),sdlog = 1.5,log=T)) +
    sum(dlnorm(sigma,meanlog = log(0.25),sdlog = 1.5,log=T)) +
    sum(dnorm(c(c.val,b),mean = 0,sd = 1,log=T))
  PR <- sum(dnorm(parm[pos.psi],mean = 0,sd = 1,log=T))
 
  LP.unc <- LL + PR + HYP
  LP <- LP.unc + sum(log(c(h,sigma)))
  Modelout <- list(LP=LP,Dev=-2*LL, Monitor=c(LP,LP.unc,LL,PR,HYP), yhat=1, parm=parm,Parameter.LL=gene.ll)
  return(Modelout)
}

Model_OUOU_correlated_evolution_sigma_off_omit_root <- function(parm, Data)
{
  ### Parameters
  pos.psi <- Data[["pos.psi"]]
  h <- exp(parm[Data[["pos.h"]]])
  b <- parm[Data[["pos.b"]]]
  c.val <- parm[Data[["pos.c"]]]
  sigma <- exp(parm[Data[["pos.sigma"]]])
  ll_fun <- Data[["ll_fun"]]
  prop <- Data[["prop"]]
  gene.ll <- Data[["Gene.LL"]]
  if (prop[1] >= pos.psi[1] && !is.null(gene.ll))
  {
    gene.index <- floor((prop[1] - 6 - 1)/2) + 1
    other.ll <- sum(gene.ll[-gene.index])
    current.ll <- ll_fun[[gene.index]](c(#parm[prop[1]] + (b - (h[2]*c.val)/h[1])*parm[prop[2]],
      NA,
      NA,
      h[1],
      0,
      h[2] * c.val - h[1] * b,
      h[2],
      #parm[prop[1]] + (b - (h[2]*c.val)/h[1])*parm[prop[2]],
      parm[prop[1]] + b*parm[prop[2]],
      parm[prop[2]],
      sigma[1],
      c.val*sigma[2],
      sigma[2]))
    LL <- current.ll + other.ll
    gene.ll[gene.index] <- current.ll
    
  } else {
    num.genes <- (length(parm)-6)/2
    gene.ll <- unlist(mclapply(1:num.genes, function(current.loci)
    {
      par.index <- 6 + 2*(current.loci-1) + 1
      ll <- ll_fun[[current.loci]](c(NA,
        NA,
        h[1],
        0,
        h[2] * c.val - h[1] * b, # c should be multiplying alpha of predictor, b should be multiplying alpha of response
        h[2],
        parm[par.index] + b*parm[par.index+1],
        parm[par.index+1],
        sigma[1],
        c.val*sigma[2], # c should be multiplying sigma of predictor
        sigma[2]))
      ll
    },mc.cores=Data$num.cores)
    )
    LL <- sum(gene.ll)
  }
  HYP <-  sum(dlnorm(h,meanlog = log(log(2)/(173/4)),sdlog = 1.5,log=T)) + 
    sum(dlnorm(sigma,meanlog = log(0.25),sdlog = 1.5,log=T)) + 
    sum(dnorm(c(c.val,b),mean = 0,sd = 1,log=T))
  PR <- sum(dnorm(parm[pos.psi[1]:length(parm)],mean = 0,sd = 1,log=T))
  LP.unc <- LL + PR + HYP
  LP <- LP.unc + sum(log(c(h,sigma)))
  Modelout <- list(LP=LP,Dev=-2*LL, Monitor=c(LP,LP.unc,LL,PR,HYP), yhat=1, parm=parm,Parameter.LL=gene.ll)
  return(Modelout)
}



dir.create(directory)
cmd <- paste(commandArgs(T),collapse = " ")
cmd <- paste("Rscript --vanilla R_scripts/runMCMC_OUOU_parallel.R",cmd,sep=" ")
readme <- paste("Model was run with bash command\n```\n",cmd,"\n```",sep="")
write(readme,file.path(directory,"README.md"))


if (independent.evo)
{
  model.type <- "OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x"
  if (!omit.root)
  {
    Model_OUOU <- Model_OUOU_independent_evolution 
  } else {
    Model_OUOU <- Model_OUOU_independent_evolution_omit_root
  }
  
} else if (sigma.off.diagonal) {
  model.type <- "OU__Global_X0__Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x"
  if (!omit.root)
  {
    Model_OUOU <- Model_OUOU_correlated_evolution_sigma_off
  } else {
    Model_OUOU <- Model_OUOU_correlated_evolution_sigma_off_omit_root
  }
  
} else {
  model.type <- "OU__Global_X0__Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x"
  if (!omit.root)
  {
    Model_OUOU <- Model_OUOU_correlated_evolution
  } else {
    Model_OUOU <- Model_OUOU_correlated_evolution_omit_root
  }
}

if (newick)
{
  tree <- read.tree(tree.file)
} else {
  tree <- read.nexus(tree.file)
}
# The tree we are working with appears to have a rounding error, making it not ultrametric
if (!is.ultrametric(tree))
{
  tree <- force.ultrametric(tree,method = "extend")
}

data <- read_tsv(input)
colnames(data)[1] <- "Gene_ID"

if (!is.null(drop.species))
{
  drop.species <- unlist(str_split(drop.species,pattern = ";"))
  data <- data %>% 
    filter(!Species %in% drop.species)
  tree <- drop.tip(tree,tip = drop.species)
  print(paste0("Species after drop: ",tree$tip.label))
  print(unique(data$Species))
}

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

print(colnames(data.pcmbase.format.trait[[1]]) == tree$tip.label)

n <- length(data.pcmbase.format.trait)

# Create likelihood function to speed up computation of likelihood                  
likFun.list <- vector(mode="list",length=n)
for (i in 1:n)
{
  if (omit.root)
  {
    positiveValueGuard <- 0
  } else {
    positiveValueGuard <- Inf
  }
  modelOU <- PCM(model= model.type,k=2)
  if (!no.std.err)
  {
    metaICPP <- PCMInfoCpp(X=data.pcmbase.format.trait[[i]],tree=tree,SE=data.pcmbase.format.se[[i]],model=modelOU)
    likFun.list[[i]] <- PCMCreateLikelihood(data.pcmbase.format.trait[[i]], tree,SE=data.pcmbase.format.se[[i]],modelOU,metaI = metaICPP,positiveValueGuard = positiveValueGuard)
  } else {
    metaICPP <- PCMInfoCpp(X=data.pcmbase.format.trait[[i]],tree=tree,model=modelOU)
    likFun.list[[i]] <- PCMCreateLikelihood(data.pcmbase.format.trait[[i]], tree,modelOU,metaI = metaICPP,positiveValueGuard = positiveValueGuard)
  }
}

min.edge <- tree$edge.length[which.min(tree$edge.length)]
alpha.start.range <- c(log(2)/(ips::tipHeights(tree)[1]/4),
                       log(2)/min.edge)




if (!independent.evo & !sigma.off.diagonal)
{
  # Get starting values. For now, use the mean of the prior for H and Sigma, 0 forQ (no affect of X on Y), and mean of each trait for optimum
  # Root state will be determined by optimum values for Y and X, as well as Q. 
  #start.values <- c(log(alpha.start),log(alpha.start),0,log(0.25),log(0.25)) # will propose H, Sigma on log-scale
  start.values <- c(log(runif(n=1,min=alpha.start.range[1],max=alpha.start.range[2])),
                    log(runif(n=1,min=alpha.start.range[1],max=alpha.start.range[2])),
                    rnorm(n=1,mean=0,sd=1),
                    log(runif(n=1,min=0,max=5)),
                    log(runif(n=1,min=0,max=5))) # will propose H, Sigma on log-scale
  
  
  
  names(start.values) <- c("H_Y","H_X","Q","Sigma_Y","Sigma_X")
  optimum.values <- unlist(
    lapply(1:length(data.pcmbase.format.trait), function(i)
    {
      tmp <- data.pcmbase.format.trait[[i]]
      #trait.mean <- rowMeans(tmp)
      trait.mean <- c(runif(1,min=-2,max=2),
                      runif(1,min=-2,max=2)
      )
      names(trait.mean) <- paste0(c("Psi_Y_","Psi_X_"),i)
      return(trait.mean)
    }
    
    )
  )
  
  start.values <- c(start.values,optimum.values + rnorm(n=length(optimum.values),0,psi.random.sd))
  
  pos.h <- grep("H",names(start.values)) 
  pos.b <- grep("Q",names(start.values))
  pos.sigma <- grep("Sigma",names(start.values))
  pos.psi <- grep("Psi",names(start.values))
  
  MyData_OUOU <- list(ll_fun = likFun.list,
                      mon.names = c("LP","LP.unc","LL","PR","HYP"),
                      parm.names=names(start.values), 
                      pos.h = pos.h, 
                      pos.b = pos.b,
                      pos.sigma=pos.sigma, 
                      pos.psi=pos.psi,
                      prop = 1:length(start.values),
                      num.cores = num.cores,
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
  
} else if (sigma.off.diagonal)
{
  #start.values <- c(log(alpha.start),log(alpha.start),0,0,log(0.25),log(0.25)) # will propose H, Sigma on log-scale
  start.values <- c(log(runif(n=1,min=alpha.start.range[1],max=alpha.start.range[2])),
                    log(runif(n=1,min=alpha.start.range[1],max=alpha.start.range[2])),
                    rnorm(n=1,mean=0,sd=1),
                    rnorm(n=1,mean=0,sd=1),
                    log(runif(n=1,min=0,max=5)),
                    log(runif(n=1,min=0,max=5))) # will propose H, Sigma on log-scale
  
  
                    
  
  names(start.values) <- c("H_Y","H_X","Q","C","Sigma_Y","Sigma_X")
  optimum.values <- unlist(
    lapply(1:length(data.pcmbase.format.trait), function(i)
    {
      tmp <- data.pcmbase.format.trait[[i]]
      #trait.mean <- rowMeans(tmp)
      trait.mean <- c(runif(1,min=-2,max=2),
                      runif(1,min=-2,max=2)
      )
      names(trait.mean) <- paste0(c("Psi_Y_","Psi_X_"),i)
      return(trait.mean)
    }
    
    )
  )
  start.values <- c(start.values,optimum.values + rnorm(n=length(optimum.values),0,psi.random.sd))

  
  
  pos.h <- grep("H",names(start.values)) 
  pos.b <- grep("Q",names(start.values))
  pos.c <- grep("C",names(start.values))
  pos.sigma <- grep("Sigma",names(start.values))
  pos.psi <- grep("Psi",names(start.values))
  MyData_OUOU <- list(ll_fun = likFun.list,
                      mon.names = c("LP","LP.unc","LL","PR","HYP"),
                      parm.names=names(start.values), 
                      pos.h = pos.h,
                      pos.b = pos.b,
                      pos.c = pos.c,
                      pos.sigma=pos.sigma,
                      pos.psi=pos.psi,
                      prop = 1:length(start.values),
                      num.cores = num.cores,
                      N=length(tree$tip.label)*n*2)
  
  
  blockwise.sample.list <- list()
  blockwise.sample.list[[1]] <- c(2,6)
  blockwise.sample.list[[2]] <- c(3,4)
  blockwise.sample.list[[3]] <- c(1,5)
  # blockwise.sample.list[[1]] <- c(2)
  # blockwise.sample.list[[2]] <- c(6)
  # blockwise.sample.list[[3]] <- c(3)
  # blockwise.sample.list[[4]] <- c(4)
  # blockwise.sample.list[[5]] <- c(1)
  # blockwise.sample.list[[6]] <- c(5)
  num.block <- length(blockwise.sample.list)
  for (i in 1:n)
  {
    pos.psi.i <- grep(paste0("Psi_[YX]_",i,"$"),names(start.values))
    blockwise.sample.list[[i+num.block]] <- pos.psi.i
  }

  
} else {
  # Get starting values. For now, use the mean of the prior for H and Sigma. Assuming independent evolution
  #start.values <- c(log(alpha.start),log(alpha.start),log(0.25),log(0.25)) # will propose H, Sigma on log-scale
  start.values <- c(log(runif(n=1,min=alpha.start.range[1],max=alpha.start.range[2])),
                    log(runif(n=1,min=alpha.start.range[1],max=alpha.start.range[2])),
                    log(runif(n=1,min=0,max=5)),
                    log(runif(n=1,min=0,max=5))) # will propose H, Sigma on log-scale
  
  names(start.values) <- c("H_Y","H_X","Sigma_Y","Sigma_X")
  optimum.values <- unlist(
    lapply(1:length(data.pcmbase.format.trait), function(i)
    {
      tmp <- data.pcmbase.format.trait[[i]]
      #trait.mean <- rowMeans(tmp)
      trait.mean <- c(runif(1,min=-2,max=2),
                      runif(1,min=-2,max=2)
      )
      names(trait.mean) <- paste0(c("Psi_Y_","Psi_X_"),i)
      return(trait.mean)
    }
    
    )
  )
  
  start.values <- c(start.values,optimum.values + rnorm(n=length(optimum.values),0,psi.random.sd))
  
  
  pos.h <- grep("H",names(start.values)) 
  pos.sigma <- grep("Sigma",names(start.values))
  pos.psi <- grep("Psi",names(start.values))
  
  MyData_OUOU <- list(ll_fun = likFun.list,
                      mon.names = c("LP","LP.unc","LL","PR","HYP"),
                      parm.names=names(start.values), 
                      pos.h = pos.h, 
                      pos.sigma=pos.sigma, 
                      pos.psi=pos.psi,
                      prop = 1:length(start.values),
                      num.cores = num.cores,
                      N=length(tree$tip.label)*n*2)
  
  
  
  blockwise.sample.list <- list()
  blockwise.sample.list[[1]] <- c(2,4)
  blockwise.sample.list[[2]] <- c(1,3)
  num.block <- length(blockwise.sample.list)
  for (i in 1:n)
  {
    pos.psi.i <- grep(paste0("Psi_[YX]_",i,"$"),names(start.values))
    blockwise.sample.list[[i+num.block]] <- pos.psi.i
  }
}

if (!is.null(start.values.file))
{
  start.values.from.file <- read_tsv(start.values.file) %>%
    column_to_rownames("Parameter") %>%
    as.matrix() 
  colnames(start.values.from.file) <- NULL
  start.values.from.file <- t(start.values.from.file)[1,]
  start.values.from.file <- start.values.from.file[which(names(start.values.from.file) %in% names(start.values))]
  start.values[names(start.values.from.file)] <- start.values.from.file[names(start.values.from.file)]
  print(start.values)
}

burnin.covar.matrix.list <- vector(mode="list",length = length(blockwise.sample.list))
for (i in 1:length(burnin.covar.matrix.list))
{
  burnin.covar.matrix.list[[i]] <- diag(1,nrow=length(blockwise.sample.list[[i]]))
}



if (is.null(prev.adapt.run) && is.null(prev.rwm.run) & burnin.samples > 0)
{
  print("Beginning burn-in...")
  burnin <- LaplacesDemon_local(Model_OUOU,
                                MyData_OUOU,
                                Initial.Values = start.values,
                                Covar = burnin.covar.matrix.list,
                                Iterations = burnin.samples,
                                Algorithm = "RWM",
                                Thinning = 1,
                                Debug = list(DB.Model = F,
                                             DB.chol = T),
                                Specs = list(B = blockwise.sample.list)
  )
  start.values <- burnin$Posterior1[burnin.samples,]
  rm(burnin)
}

if(!is.null(prev.adapt.run)) {
  print("Loading previous adaptive run")
  load(prev.adapt.run)
  
} else if(!is.null(prev.rwm.run)) {
  print("Loading previous RWM run")
  load(prev.rwm.run)
  adapt.fit <- fit
}

if (adapt.samples > 0 && is.null(prev.adapt.run) && is.null(prev.rwm.run))
{
  print("Beginning adaptive MCMC...")
  adapt.fit <- LaplacesDemon_local(Model_OUOU,
                                   MyData_OUOU,
                                   Initial.Values = start.values,
                                   Iterations = adapt.samples,
                                   Algorithm = "RAM",
                                   Status = 100,
                                   Thinning = thinning,
                                   Debug = list(DB.Model = T,
                                                DB.eigen=T,
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
    if (sigma.off.diagonal)
    {
      plot(adapt.fit$Posterior1[range.to.plot,"C"],type="l",ylab="C")
    }
  }
  plot(adapt.fit$Posterior1[range.to.plot,"Sigma_Y"],type="l",ylab=expression(sigma["Y"]))
  plot(adapt.fit$Posterior1[range.to.plot,"Sigma_X"],type="l",ylab=expression(sigma["X"]))
  dev.off()
  
} else if (adapt.samples > 0 & !is.null(prev.adapt.run) && is.null(prev.rwm.run)) {
  print("Resuming adaptive mcmc...")
  prev.samp <- nrow(adapt.fit$Posterior1)
  prev.iter <- adapt.fit$Iterations
  prev.covar <- adapt.fit$Covar
  start.values <- adapt.fit$Posterior1[prev.samp,]
  covar.matrix <- adapt.fit$Covar
  adapt.fit.new <- LaplacesDemon_local(Model_OUOU,
                                   MyData_OUOU,
                                   Initial.Values = start.values,
                                   Covar=covar.matrix,
                                   Iterations = adapt.samples,
                                   Algorithm = "RAM",
                                   Thinning = thinning,
                                   Debug = list(DB.Model = T,
                                                DB.eigen=T,
                                                DB.chol = T),
                                   Specs = list(alpha.star = alpha.star,
                                                B = blockwise.sample.list,
                                                Dist="t",
                                                gamma = 2/3,
                                                n=prev.iter
                                   ))
  
  adapt.fit <- Combine(list(adapt.fit,adapt.fit.new),Data=MyData_OUOU)
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
    if (sigma.off.diagonal)
    {
      plot(adapt.fit$Posterior1[range.to.plot,"C"],type="l",ylab="C")
    }
  }
  plot(adapt.fit$Posterior1[range.to.plot,"Sigma_Y"],type="l",ylab=expression(sigma["Y"]))
  plot(adapt.fit$Posterior1[range.to.plot,"Sigma_X"],type="l",ylab=expression(sigma["X"]))
  dev.off()
}
num.prev.sample <- nrow(adapt.fit$Posterior1)
start.values <- adapt.fit$Posterior1[num.prev.sample,]
covar.matrix <- adapt.fit$Covar

if (samples > 0)
{
  print("Beginning Random-Walk Metropolis MCMC for parameter estimation...")
  fit <- LaplacesDemon_local(Model_OUOU,
                             MyData_OUOU,
                             Initial.Values = start.values,
                             Covar=covar.matrix,
                             Iterations = samples,
                             Algorithm = "RWM",
                             Thinning = thinning,
                             Debug = list(DB.Model = T,
                                          DB.eigen=T,
                                          DB.chol = T),
                             Specs = list(B = blockwise.sample.list)
  )
  
  save(fit,file=file.path(directory,paste0("rwm_mcmc_",rwm.run.number,".Rda")))
  
  if (!is.null(prev.rwm.run))
  {
    fit <- Combine(list(adapt.fit,fit),Data=MyData_OUOU)
    save(fit,file=file.path(directory,"rwm_mcmc_total.Rda"))
  }
  
  
  burnin<-1
  nsample <- nrow(fit$Monitor)
  range.to.plot <- burnin:nsample
  pdf(file.path(directory,"mcmc_traces.pdf"))
  plot(fit$Monitor[range.to.plot,"LP"],type="l",ylab="Log(Posterior)")
  plot(fit$Monitor[range.to.plot,"LL"],type="l",ylab="Log(Likelihood)")
  plot(fit$Posterior1[range.to.plot,"H_Y"],type="l",ylab=expression(alpha["Y"]))
  plot(fit$Posterior1[range.to.plot,"H_X"],type="l",ylab=expression(alpha["X"]))
  if (!independent.evo)
  {
    plot(fit$Posterior1[range.to.plot,"Q"],type="l",ylab="Q")
    if (sigma.off.diagonal)
    {
      plot(fit$Posterior1[range.to.plot,"C"],type="l",ylab="C")
    }
  }
  plot(fit$Posterior1[range.to.plot,"Sigma_Y"],type="l",ylab=expression(sigma["Y"]))
  plot(fit$Posterior1[range.to.plot,"Sigma_X"],type="l",ylab=expression(sigma["X"]))
  dev.off()
}


