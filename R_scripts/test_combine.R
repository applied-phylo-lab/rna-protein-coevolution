library(PCMBase,lib.loc = "~/R_dev")
library(PCMBaseCpp)
library(LaplacesDemon,lib.loc = "~/R_dev/")
require(coda)
library(tidyverse)
library(phytools)
library(argparse)
library(bayestestR)
source("/data2/cope/rna-protein-coevolution/R_notebooks/local_functions.R")


parser <- ArgumentParser()
parser$add_argument("-i","--input",help="tsv file with data input",type="character",default="./")
parser$add_argument("-o","--output",help="Directory of where to put results. Will automatically generate lowest-level directory if not already generated.",type="character",default="./")
parser$add_argument("--tree",type="character",default="../Data/tree_11sp_noGpig.nex")
parser$add_argument("--newick",action="store_true")
parser$add_argument("--independent_evo",action="store_true")
parser$add_argument("--reverse_trait_order",action="store_true")


args <- parser$parse_args()
print(args)
input <- args$input
directory <- args$output
tree.file <- args$tree
newick <- args$newick
independent.evo <- args$independent_evo
reverse.trait.order <- args$reverse_trait_order


if (!independent.evo)
{
  model.type <- "OU__Global_X0__Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x"
  #model.type <- "OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x"
} else {
  model.type <- "OU__Global_X0__Diagonal_WithNonNegativeDiagonal_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x"
  #model.type <- "OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x"

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

# the order of the traits will determine direction of "causality". First trait evolves towards optimum determined by the second trait
if (reverse.trait.order)
{
  print("Reversing order of traits...")
  data <- data %>% 
    relocate(Mean_Protein,.before=Mean_RNA) %>% 
    relocate(SE_Protein,.before=SE_RNA)
}

data.se <- data %>%
  dplyr::select(Gene_ID,Species,contains("SE")) %>%
  group_by(Gene_ID) %>%
  group_split()
data.gene <- data %>%
  dplyr::select(Gene_ID,Species,contains("Mean")) %>%
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

data.pcmbase.format.se <- purrr::map(data.se,function(x){
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
  likFun.list[[i]] <- PCMCreateLikelihood(data.pcmbase.format.trait[[i]], tree,SE=data.pcmbase.format.se[[i]],modelOU,metaI = PCMInfoCpp)
}

if (!independent.evo)
{
  # Get starting values. For now, use the mean of the prior for H and Sigma, 0 for Q (no affect of X on Y), and mean of each trait for optimum
  # Root state will be determined by optimum values for Y and X, as well as Q. 
  start.values <- c(log(0.25),log(0.25),0,log(0.25),log(0.25)) # will propose H, Sigma on log-scale
  names(start.values) <- c("H_Y","H_X","Q","Sigma_Y","Sigma_X")
  optimum.values <- unlist(
    lapply(1:length(data.pcmbase.format.trait), function(i)
    {
      tmp <- data.pcmbase.format.trait[[i]]
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
  
} else {
  # Get starting values. For now, use the mean of the prior for H and Sigma. Assuming independent evolution
  start.values <- c(log(0.25),log(0.25),log(0.25),log(0.25)) # will propose H, Sigma on log-scale
  names(start.values) <- c("H_Y","H_X","Sigma_Y","Sigma_X")
  optimum.values <- unlist(
    lapply(1:length(data.pcmbase.format.trait), function(i)
    {
      tmp <- data.pcmbase.format.trait[[i]]
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
  
}

# load(file.path(directory,"rwm_mcmc.Rda"))
# fit.1 <- fit
load(file.path(directory,"rwm_mcmc_2.Rda"))
fit.2 <- fit
load(file.path(directory,"rwm_mcmc_3.Rda"))
fit.3 <- fit
fit.total <- Combine(list(fit.2,fit.3),Data=MyData_OUOU)
save(fit.total,file=file.path(directory,"rwm_mcmc_2_3.Rda"))