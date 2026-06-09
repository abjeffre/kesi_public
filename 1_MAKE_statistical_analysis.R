###################################
######## BUILD HELPER FUNCTIONS ###

library(rethinking)
library(cmdstanr)
library(posterior)
library(matrixStats)
library(shape)
library(stringr)

###################################
######## BUILD HELPER FUNCTIONS ###

create_empty_like <- function(obj) {
  if (is.vector(obj)) {
    return(vector(mode = typeof(obj), length = length(obj)))
  } else if (is.matrix(obj)) {
    return(matrix(NA, nrow = nrow(obj), ncol = ncol(obj)))
  } else if (is.data.frame(obj)) {
    return(obj[0, ])
  } else if (is.array(obj)) {
    return(array(NA, dim = dim(obj)))
  } else if (is.list(obj)) {
    empty_list <- vector("list", length(obj))
    names(empty_list) <- names(obj)
    return(empty_list)
  } else {
    stop("Unsupported object type")
  }
}

get_init_list<-function(pf){
  stanfit = posterior::as_draws_rvars(pf)
  init_list <- list()
  cnt <- 1
  for(i in names(stanfit)){
    print(i)
    means <- summarise_draws(as_draws_array(stanfit[[i]]))
    temp<-create_empty_like(stanfit[[i]])
    for(j in 1:length(means$mean)){
      temp[j] = c(means$mean[j])
    }
    init_list[[cnt]]<-as.array(temp)
    cnt <- cnt+1
  }
  
  names(init_list) <- names(stanfit)
  for(i in names(init_list)){
    if(dims(init_list[[i]])==1 & length(init_list[[i]]) ==1) init_list[[i]] <- init_list[[i]][1]
  }
  return(init_list)
}

extract.samples2 <- function(x){
  output <- list()
  stanfit = posterior::as_draws_rvars(x)
  for(i in names(stanfit)){
    output[[i]] = posterior::draws_of(stanfit[[i]])
  }
  return(output)
}  



##################################
########## LOAD DATA #############


data <- readRDS("~/kesi/data/data_kesi2025-09-22.RDS")

##########################
######### FULL MODEL #####

full_24_model <- cmdstan_model("kesi/code/stan_models/main_model.stan")

full_hmc <- full_24_model$sample(
  data = data,
  iter_sampling = 500,
  iter_warmup = 500,
  chains = 4,
  init = 0,
  parallel_chains = 4,
  refresh = 2
)
post <- extract.samples(full_hmc)
saveRDS(post, "data/full_hmc.RDS")



