knitr::opts_chunk$set(echo = TRUE)
library(PCMBaseCpp)
library(adaptMCMC)
require(coda)
library(tidyverse)
library(cowplot)
library(ggpubr)
library(tidyverse)
library(argparse)

prior_BM <- function(x)
{
  x0.A <- dnorm(x[1],mean = 0, sd = 1, TRUE) # prior for the X0.A
  x0.B <- dnorm(x[2],mean = 0, sd = 1, TRUE) # prior for X0.B
  sigma.A <- dlnorm(x[3],meanlog = 0.25,sdlog = 1.5, log= TRUE) # prior for the standard deviations of trait A
  sigma.AB <- dlnorm(x[4],meanlog = 0.25,sdlog = 1.5, log= TRUE) # prior for the off-diagonal of symmetric sigma matrix
  sigma.B <- dlnorm(x[5],meanlog = 0.25,sdlog = 1.5, log= TRUE) # prior for the standard deviations of trait B
  return(x0.A + x0.B + sigma.A + sigma.AB + sigma.B)
}


# define the log-likelihood distribution
log_post_BM<- function(par,ll_fun)
{
  ll <- ll_fun(par)
  pr <- prior_BM(par)
  return(ll+pr)
}

prior_OU <- function(x)
{
  x0.A <- dnorm(x[1],mean = 0, sd = 1, TRUE) # prior for the X0.A
  x0.B <- dnorm(x[2],mean = 0, sd = 1, TRUE) # prior for X0.B
  alpha.A <- dlnorm(x[3],meanlog = 0.25,sdlog = 1.5, log= TRUE) # prior for the alphaof trait A
  alpha.BA <- dlnorm(x[4],meanlog = 0.25,sdlog = 1.5, log= TRUE) # prior for the off-diagonal of non-symmetric H matrix
  alpha.AB <- dlnorm(x[5],meanlog = 0.25,sdlog = 1.5, log= TRUE) # prior for the off-diagonal of non-symmetric H matrix
  alpha.B <- dlnorm(x[6],meanlog = 0.25,sdlog = 1.5, log= TRUE) # prior for the alpha of trait B
  theta.A <- dnorm(x[7],mean = 0, sd = 1, TRUE) # prior for the optimum of trait A
  theta.B <- dnorm(x[8],mean = 0, sd = 1, TRUE) # prior for the optimum of trait B
  sigma.A <- dlnorm(x[9],meanlog = 0.25,sdlog = 1.5, log= TRUE) # prior for the standard deviations of trait A
  sigma.AB <- dlnorm(x[10],meanlog = 0.25,sdlog = 1.5, log= TRUE) # prior for the off-diagonal of symmetric sigma matrix
  sigma.B <- dlnorm(x[11],meanlog = 0.25,sdlog = 1.5, log= TRUE) # prior for the standard deviations of trait B
  return(x0.A + x0.B + alpha.A + alpha.BA + alpha.AB + alpha.B + theta.A + theta.B + sigma.A + sigma.AB + sigma.B)
}


# define the log-likelihood distribution
log_post_OU<- function(par,ll_fun){
  ll <- ll_fun(par)
  pr <- prior_OU(par)
  return(ll+pr)
}


hyper_prior_bayou <- function(x)
{
  alpha.A <- dlnorm(x[1],meanlog = 0.25,sdlog = 1.5, log= TRUE) # prior for the alphaof trait A
  alpha.BA <- dlnorm(x[2],meanlog = 0.25,sdlog = 1.5, log= TRUE) # prior for the off-diagonal of non-symmetric H matrix
  alpha.AB <- dlnorm(x[3],meanlog = 0.25,sdlog = 1.5, log= TRUE) # prior for the off-diagonal of non-symmetric H matrix
  alpha.B <- dlnorm(x[4],meanlog = 0.25,sdlog = 1.5, log= TRUE) # prior for the alpha of trait B
  sigma.A <- dlnorm(x[5],meanlog = 0.25,sdlog = 1.5, log= TRUE) # prior for the standard deviations of trait A
  sigma.AB <- dlnorm(x[6],meanlog = 0.25,sdlog = 1.5, log= TRUE) # prior for the off-diagonal of symmetric sigma matrix
  sigma.B <- dlnorm(x[7],meanlog = 0.25,sdlog = 1.5, log= TRUE) # prior for the standard deviations of trait B
  return(alpha.A + alpha.BA + alpha.AB + alpha.B + sigma.A + sigma.AB + sigma.B)
}

prior_OU_loci <- function(x)
{
  x0.A <- dnorm(x[1],mean = 0, sd = 1, TRUE) # prior for the X0.A
  x0.B <- dnorm(x[2],mean = 0, sd = 1, TRUE) # prior for X0.B
  theta.A <- dnorm(x[3],mean = 0, sd = 1, TRUE) # prior for the optimum of trait A
  theta.B <- dnorm(x[4],mean = 0, sd = 1, TRUE) # prior for the optimum of trait B
  return(x0.A + x0.B + theta.A + theta.B)
}

# Treat the first 7 parameters as the H and Sigma matrix. Remaining are X0_A, X0_B, Theta_A, Theta_B for each loci
# | Alpha A | Alpha BA | Alpha AB | Alpha B | Sigma A | Sigma AB | Sigma B |
log_post_hier<- function(par,ll_fun)
{
  lp <- 0
  hypr <- hyper_prior_bayou(par[1:7])
  num.genes <- (length(par)-7)/4 
  lp <- sum(unlist(lapply(1:num.genes, function(current.loci)
  {
    par.index <- 7 + 4*(current.loci-1) + 1
    tmp <- c(par[par.index],par[par.index+1],par[1],par[2],par[3],par[4],par[par.index+2],par[par.index+3],par[5],par[6],par[7])
    ll <- ll_fun[[current.loci]](tmp)
    pr <- prior_OU_loci(c(par[par.index],par[par.index+1],par[par.index+2],par[par.index+3]))
    ll + pr
    
  })))
  lp <- lp + hypr
  return(lp)
}


# Treat the first 7 parameters as the H and Sigma matrix. Remaining are X0_A, X0_B, Theta_A, Theta_B for each loci
# | Alpha A | Alpha BA | Alpha AB | Alpha B | Sigma A | Sigma AB | Sigma B |
log_post_hier_parallel<- function(par,ll_fun,n.cores)
{
  lp <- 0
  hypr <- hyper_prior_bayou(par[1:7])
  num.genes <- (length(par)-7)/4 
  lp <- sum(unlist(mclapply(1:num.genes, function(current.loci)
  {
    par.index <- 7 + 4*(current.loci-1) + 1
    tmp <- c(par[par.index],par[par.index+1],par[1],par[2],par[3],par[4],par[par.index+2],par[par.index+3],par[5],par[6],par[7])
    ll <- ll_fun[[current.loci]](tmp)
    pr <- prior_OU_loci(c(par[par.index],par[par.index+1],par[par.index+2],par[par.index+3]))
    ll + pr
    
  },mc.cores = n.cores)))
  lp <- lp + hypr
  return(lp)
}


tracePlot <- function(model.fit,truth,type="X0_A")
{
  num.genes <- (ncol(model.fit$samples)-7)/4 
  if (type == "X0_A")
  {
    index <- seq(8,ncol(model.fit$samples),by=4)
  } else if (type == "X0_B"){
    index <- seq(9,ncol(model.fit$samples),by=4)
  } else if (type == "Theta_A"){
    index <- seq(10,ncol(model.fit$samples),by=4)
  } else if (type == "Theta_B"){
    index <- seq(11,ncol(model.fit$samples),by=4)
  }
  param.matrix <- model.fit$samples[,index]
  
  post.means <- colMeans(param.matrix)
  param.df <- data.frame(True.Values = truth[index],
                         Posterior.Mean = post.means)
  p <- ggplot(param.df,aes(x=True.Values,y=Posterior.Mean)) +
    geom_point() +
    xlab("Truth") +
    ylab("Posterior Means") +
    ggtitle(type) +
    theme_cowplot() +
    theme(aspect.ratio=1) +
    stat_cor() +
    geom_abline(intercept=0,slope = 1,linetype="dashed")
  return(p)
  
}


createStartValues <- function(data,model)
{
  
}



parser <- ArgumentParser()
parser$add_argument("--data",help="Tab-separated file (.tsv) storing trait data. First column should represent the higher-level trait (e.g. a gene), second column should be the species, and the remaining columns should represent traits.",type="character",default="./")
parser$add_argument("--tree",help="File storing phylogenetic tree. Specify format using --tree_format")
parser$add_argument("--tree_format",help="Format of tree. Valid options are newick and nexus",default="newick")
parser$add_argument("--output",help="Directory of where to put results. Will automatically generate lowest-level directory if not already generated.",type="character",default="./MCMC")
parser$add_argument("--model",help="a string describing the parameterization of a PCM model in the PCMBase format",default=NULL)
parser$add_argument("--model_file",help="a file path that leads to a file containing a PCM model in PCMBase format as a single line. Supersedes --model.",default = NULL)
parser$add_argument("--adapt",help="Will turn on adapting via the robust adaptive MCMC algorithm",action="store_true")
parser$add_argument("--accept",help="Target acceptance rate",type="double")
parser$add_argument("-n","--threads",help="Number of threads to use for MCMC",type="integer",default=1)
parser$add_argument("--iterations",help="Number of iteration for MCMC",type="integer",default=10000)




args <- parser$parse_args()
directory <- args$output
data.file <- args$data
tree.file <- args$tree
tree.format <- args$tree_format
directory <- args$output
model <- args$model
model.file <- args$model_file
adapt <- args$adapt
accept <- args$accept
iterations <- args$iterations
n.cores <- args$threads


dir.create(directory)

cmd <- paste(commandArgs(T),collapse = " ")
cmd <- paste("Rscript --vanilla R_scripts/runMCMC.R",cmd,sep=" ")
readme <- paste("Model was run with bash command\n```\n",cmd,"\n```",sep="")
write(readme,file.path(directory,"README.md"))

data <- read_tsv(data.file)
data.split <- data %>%
  group_by(Higher.Trait) %>%
  group_split() %>%
  purrr::map(.x %>% 
               dplyr::select(-Higher.Trait) %>% 
               column_to_rownames(Species) %>% 
               as.matrix() %>%
               t())

if (is.null(tree))
{
  stop("ERROR: file containing phylogenetic tree not specified")
} else {
  if (tree.format == "newick")
  {
    tree <- read.tree(tree.file)
  } else if (tree.format == "nexus"){
    tree <- read.nexus(tree.file)
  } else {
    stop("ERROR: please specify if tree is in Newick (newick) or Nexus (nexus) format")
  }
}
if (!is.null(model.file))
{
  model <- readLines(model.file)
} else if (is.null(model))
{
  stop("ERROR: model not specified")
}

num.traits <- length(data.split)
llFun.list <- vector(mode = "list", length=num.traits)
model.obj <- PCM(model=model,k=nrow(data.split[[1]]))

for (i in 1:num.traits)
{
  likFun.list[[i]] <- PCMCreateLikelihood(data.split[[i]], tree, model.obj)
}

start.val <- createStartValues()

mcmc.fit <- MCMC(log_post_hier_parallel, iterations, init = start.val, scale = rep(1,length(start.val)),
                               adapt = adapt, acc.rate = acc.rate, gamma = 2/3, ll_fun=likFun.list)

save(mcmc.fit,file.path(directory,"mcmc.Rda"))

