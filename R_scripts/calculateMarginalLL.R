library(PCMBase,lib.loc = "~/R_dev")
library(PCMBaseCpp)
library(LaplacesDemon,lib.loc = "~/R_dev/")
require(coda)
library(tidyverse)
library(phytools)
library(argparse)
source("/data2/cope/rna-protein-coevolution/R_notebooks/local_functions.R")

LML <- function(Model=NULL, Data=NULL, Modes=NULL, theta=NULL, LL=NULL,
                Covar=NULL, method="NSIS")
{
  LML.out <- list(LML=NA, VarCov=NA)
  if(!is.null(Modes) & !is.null(theta) & !is.null(LL) & method == "GD") {
      Interval <- 1.0e-6
      parm.len <- length(Modes)
      if(is.null(Covar)) {
        eps <- Interval * Modes
        Approx.Hessian <- Hessian(Model, Modes, Data)
        Covar <- try(-as.inverse(Approx.Hessian), silent=TRUE)
        if(!inherits(VarCov, "try-error"))
          diag(Covar)[which(diag(Covar) <= 0)] <- .Machine$double.eps
        else {
          cat("\nWARNING: Failure to solve matrix inversion of ",
              "Approx. Hessian in LML.\n", sep="")
          cat("NOTE: Identity matrix is supplied instead.\n")
          Covar <- diag(parm.len)
        }
    }
    else {
      VarCov <- Covar
      rm(Covar)
    }
    log.g.theta <- dmvn(theta, Modes, Covar/1.5, log=TRUE)
    #LML <- log(1 / mean(exp(log.g.theta - LL)))
    N <- length(log.g.theta)
    x <- log.g.theta - LL
    x.max <- x[which.max(x)]
    x.rescale <- x - x.max
    LML <- log(1) - ((log(1) - log(N)) + (x.max + log(sum(exp(x.rescale)))))
    ### Output
    LML.out <- list(LML=LML, VarCov=NA)
  }
  if(!is.null(LL) & method == "HME") {
    med <- median(LL)
    LML <- med - log(mean(exp(-LL + med)))
    ### Output
    LML.out <- list(LML=LML, VarCov=NA)
  }
  if(!is.null(Model) & !is.null(Data) & !is.null(Modes) &
     (method == "LME")) {
    Interval <- 1.0e-6
    parm.len <- length(Modes)
    if(is.null(Covar)) {
      eps <- Interval * Modes
      Approx.Hessian <- Hessian(Model, Modes, Data)
      VarCov <- try(-as.inverse(Approx.Hessian), silent=TRUE)
      if(!inherits(VarCov, "try-error"))
        diag(VarCov)[which(diag(VarCov) <= 0)] <- .Machine$double.eps
      else {
        cat("\nWARNING: Failure to solve matrix inversion of ",
            "Approx. Hessian in LML.\n", sep="")
        cat("NOTE: Identity matrix is supplied instead.\n")
        VarCov <- diag(parm.len)}
    }
    else {
      VarCov <- Covar
      rm(Covar)}
    ### Logarithm of the Marginal Likelihood
    LML <- NA
    options(warn=-1)
    LML.test <- try(parm.len/2 * log(2*pi) + 0.5*logdet(VarCov) +
                      as.vector(Model(Modes, Data)[["LP"]]), silent=TRUE)
    if(!inherits(LML.test, "try-error")) LML <- LML.test[1]
    options(warn=0)
    ### Output
    LML.out <- list(LML=LML, VarCov=VarCov)
  }
  if(!is.null(theta) & !is.null(LL) & (method =="NSIS")) {
    if(!is.matrix(theta)) stop("theta must be a matrix.")
    thetacol <- ncol(theta)
    LL <- as.vector(LL)
    LLlen <- length(LL)
    if(nrow(theta) != LLlen)
      stop("The number of rows in theta differs from the ",
           "length of LL.")
    if(LLlen < 301) {
      cat("\nWARNING: At least 301 samples are required for NSIS in LML.\n")
      return(list(LML=NA, VarCov=NA))}
    if(thetacol > round(LLlen / 2)) {
      cat("\nWARNING: The number of parameters,", thetacol,
          ",\n exceeds half the number of stationary samples,",
          round(LLlen / 2), ",\n required for NSIS.\n")
      return(list(LML=NA, VarCov=NA))}
    cov.prob <- 0.5
    bounds <- matrix(c(-Inf, Inf), 2, ncol(theta))
    .GetID <- function(point, center, width)
    {return(paste(floor((point - center) / width),
                  collapse=" "))}
    .PopulateHist <- function(points, heights, width)
    {
      max.point <- points[heights == max(heights), ,
                          drop=FALSE][1, ]
      center <- max.point - width / 2
      hist.env <- new.env(hash=TRUE, parent=emptyenv())
      for (i in seq(along=heights)) {
        id <- .GetID(points[i, ], center, width)
        if(exists(id, envir=hist.env))
          cur.tuple <- get(id, envir=hist.env)
        else cur.tuple <- c(Inf, 0)
        assign(id, c(min(heights[i], cur.tuple[1]),
                     cur.tuple[2] + 1), envir=hist.env)
      }
      return(list(env=hist.env, center=center, width=width,
                  dim=length(center)))
    }
    .Profile <- function(hist)
    {
      ids <- ls(hist$env)
      counts <- sapply(ids, function(id) get(id,
                                             envir=hist$env)[2])
      names(counts) <- ids
      return(sort(counts, decreasing=TRUE))
    }
    .Lookup <- function(point, hist)
    {
      id <- .GetID(point, hist$center, hist$width)
      return(ifelse(exists(id, envir=hist$env),
                    get(id, envir=hist$env)[1], -Inf))
    }
    .SetHistNorm <- function(hist)
    {
      ids <- ls(hist$env)
      heights <- sapply(ids, function(id) get(id,
                                              envir=hist$env)[1])
      max.height <- max(heights)
      hist$norm <- (log(hist$width ^ hist$dim * sum(exp(heights - max.height))) + max.height)
      return(hist)
    }
    .Coverage <- function(points, hist)
    {
      heights <- apply(points, 1, function(point) .Lookup(point,
                                                          hist))
      return(length(heights[heights > -Inf]) / length(heights))
    }
    .Dist <- function(p1, p2)
    {
      return(max(abs(p1 - p2)))
    }
    .GetWidth <- function(theta.hist, theta.width, LL.hist,
                          opt.prob=0.5)
    {
      low.point <- theta.hist[which(LL.hist == min(LL.hist))[1], ,
                              drop=FALSE]
      high.point <- theta.hist[which(LL.hist == max(LL.hist))[1], ,
                               drop=FALSE]
      minLL.dists <- apply(theta.hist, 1, function(theta) {
        .Dist(theta, low.point)})
      maxLL.dists <- apply(theta.hist, 1, function(theta) {
        .Dist(theta, high.point)})
      small.dist <- 0.5 * min(maxLL.dists[maxLL.dists > 0])
      big.dist <- 2 * max(minLL.dists)
      F <- function(width) {
        hist <- .PopulateHist(theta.hist, LL.hist, width)
        return(.Coverage(theta.width, hist) - opt.prob)
      }
      options(warn=-1)
      solve.test <- try(uniroot(F, lower=small.dist,
                                upper=big.dist, tol=min(.05, 2/nrow(theta.width))),
                        silent=TRUE)
      options(warn=0)
      if(!inherits(solve.test, "try-error")) return(solve.test$root)
      else return(1)
    }
    .GetLimits <- function(hist)
    {
      split.ids <- strsplit(ls(hist$env), " ")
      coords <- matrix(NA, ncol=hist$dim, nrow=length(split.ids))
      for (i in seq(along=split.ids))
        coords[i, ] <- as.integer(split.ids[[i]])
      max.coords <- apply(coords, 2, max)
      min.coords <- apply(coords, 2, min)
      result <- rbind(apply(coords, 2, min) * hist$width +
                        hist$center, (apply(coords, 2, max) + 1) * hist$width +
                        hist$center)
      rownames(result) <- c("min", "max")
      return(result)
    }
    .MargLL <- function(hist, theta.imp, LL.imp)
    {
      n <- length(LL.imp)
      samples <- rep(0, n)
      for (i in seq(length=n))
      {
        samples[i] <- .Lookup(theta.imp[i, ],hist) - LL.imp[i]
      }
      max.sample <- samples[which.max(samples)]
      samples.rescaled <- samples - max.sample
      return(list(samples=samples, 
                  mll=log(1) - ((log(1) - log(n)) + (max.sample + sum(exp(samples.rescaled))))
            ))
    }
    .SplitComponents <- function(theta, LL)
    {
      tot.count <- length(LL)
      hist.count <- min(0.2 * tot.count,
                        round(sqrt(tot.count) * 2))
      width.count <- 40
      imp.count <- tot.count - hist.count - width.count
      return(list(theta.hist=theta[1:hist.count, , drop=FALSE],
                  theta.bw=theta[(hist.count + 1):(hist.count +
                                                     width.count), , drop=FALSE],
                  theta.imp=theta[(tot.count - imp.count + 1):tot.count, ,
                                  drop=FALSE],
                  LL.hist=LL[1:hist.count],
                  LL.imp=LL[(tot.count - imp.count + 1):tot.count]))
    }
    .CheckBounds <- function(hist, bounds)
    {
      limits <- hist$limits <- .GetLimits(hist)
      if(!all((bounds[1, ] == -Inf) |
              (bounds[1, ] < limits[1, ])))
        stop(paste("Bounds error: histogram lower limits (",
                   paste(limits[1, ], collapse=" "),
                   ") lower than specified bounds (",
                   paste(bounds[1, ], collapse=" "), ").", sep=""))
      if(!all((bounds[2, ] == Inf) | (bounds[2, ] > limits[2, ])))
        stop(paste("Bounds error: histogram upper limits (",
                   paste(limits[2, ], collapse=" "),
                   ") greater than specified bounds (",
                   paste(bounds[2, ], collapse=" "), ").", sep=""))
      return(hist)
    }
    options(warn=-1)
    comps <- .SplitComponents(theta, LL)
    width <- .GetWidth(comps$theta.hist, comps$theta.bw,
                       comps$LL.hist, cov.prob)
    hist <- .PopulateHist(comps$theta.hist, comps$LL.hist, width)
    hist <- .SetHistNorm(hist)
    hist <- .CheckBounds(hist, bounds)
    ml.out <- .MargLL(hist, comps$theta.imp, comps$LL.imp)
    ml.out$samples <- ml.out$samples[is.finite(ml.out$samples)]
    sd.samples <- sd(ml.out$samples) / sqrt(length(ml.out$samples))
    conf.interval <- qnorm(c(0.975, 0.025)) * sd.samples +
      mean(ml.out$samples)
    conf.interval[which(conf.interval < .Machine$double.eps)] <- .Machine$double.eps
    conf.interval <- -log(conf.interval)
    options(warn=0)
    LML <- ml.out$mll + hist$norm
    LML[which(!is.finite(LML))] <- NA
    ### Output
    LML.out <- list(LML=LML, VarCov=NA)
  }
  return(LML.out)
}

#End

parser <- ArgumentParser()
parser$add_argument("-i","--input",help="tsv file with data input",type="character",default="./")
parser$add_argument("--tree",type="character",default="../Data/tree_11sp_noGpig.nex")
parser$add_argument("--prev_rwm_run",type="character",default=NULL)
parser$add_argument("--newick",action="store_true")
parser$add_argument("--independent_evo",action="store_true")
parser$add_argument("--reverse_trait_order",action="store_true")
parser$add_argument("--omit_root",action="store_true")


args <- parser$parse_args()
print(args)
input <- args$input
tree.file <- args$tree
newick <- args$newick
independent.evo <- args$independent_evo
reverse.trait.order <- args$reverse_trait_order
omit.root <- args$omit_root
prev.rwm.run <- args$prev_rwm_run


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
      ll <- ll_fun[[current.loci]](c(parm[par.index],
                                     parm[par.index+1],
                                     h[1],
                                     h[2],
                                     parm[par.index],
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
    gene.ll <- unlist(lapply(1:num.genes, function(current.loci)
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
    gene.ll <- unlist(lapply(1:num.genes, function(current.loci)
    {
      par.index <- 5 + 2*(current.loci-1) + 1
      ll <- ll_fun[[current.loci]](c(NA,
                                     NA,
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



load(prev.rwm.run)
posterior.means <- colMeans(fit$Posterior2)
#covar <- cov(fit$Posterior2)


lml <- LML(Model=Model_OUOU,
    Data=MyData_OUOU,
    Modes = posterior.means,
    theta=fit$Posterior2,
    LL = fit$Monitor[,"LP"], #documentation says loglikelihood, but really want logposterior, see https://doi.org/10.1016/j.jmp.2017.09.005
    Covar = NULL,
    method="GD")

# lml <- LML(theta=fit$Posterior2,
#            LL = fit$Monitor[,"LL"],
#            method="NSIS")

# fit$LML <- lml$LML
print(fit$LML)

