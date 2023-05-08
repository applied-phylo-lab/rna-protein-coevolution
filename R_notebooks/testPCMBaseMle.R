---
title: "Test PCMBase Likelihood"
author: "Alex Cope"
date: '2023-05-06'
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(ape)
library(PCMBase)
library(PCMBaseCpp)
library(tidyverse)
library(cowplot)
```

```{r}
negL <- function(par,
                 likFun,
                 fixed=c(),
                 param.order = c("X0_A","X0_B","H_A","H_BA","H_AB","H_B","Theta_A","Theta_B","Sigma_A","Sigma_AB","Sigma_B"))
{
  all.par <- c(par,fixed)
  all.par <- all.par[par.order]
  -likFun(all.par,log = T)
}

generateRandomParams <- function()
{
  x0 <- rnorm(n=2,mean=0,sd=1)
  theta <- x0
  alpha <- rlnorm(n=4,meanlog = 0.25,sdlog = 1.5)
  sigma <- rlnorm(n=3,meanlog = 0.25,sdlog = 1.5)
  return(data.frame(X0_A = x0[1],
                    X0_B = x0[2],
                    H_A = alpha[1],
                    H_BA = alpha[2],
                    H_AB = alpha[3],
                    H_B = alpha[4],
                    Theta_A = theta[1],
                    Theta_B = theta[2],
                    Sigma_A = sigma[1],
                    Sigma_AB = sigma[2],
                    Sigma_B = sigma[3]))
}

generateSimData <- function(model.type,tree,params,root.state,k=2,offset=0)
{
  model <- PCM(model=model.type,k=k)
  PCMParamLoadOrStore(model, params, offset=offset, load=T)
  data.pcm <- PCMSim(tree = tree,model = model,X0 = root.state)
  return(data.pcm)
}
ba.tree <- read.nexus("Data/tree_11sp_noGpig.nex")
plot(ba.tree)
ba.tree <- force.ultrametric(ba.tree,method = "extend")
plot(ba.tree)
```


```{r}
df <- lapply(1:50,function(x) 
  {
  generateRandomParams()
}) %>% bind_rows()

model.type <- "OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x"

sim.data <- apply(df,
                  MARGIN = 1,
                  FUN = generateSimData,
                  model.type = model.type,
                  )




```
