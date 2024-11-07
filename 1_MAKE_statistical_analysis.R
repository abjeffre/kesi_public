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


data <- readRDS("data/data.RDS")


##################################
########## RUN MODELS ############

# Set up pathfinder for efficent sampling

m1 <- cmdstan_model("code/stan_models/model_sectors.stan")
pf<-m1$pathfinder(data = data, init = 0, num_paths = 4)

# Extract Initalization values
init_list<-get_init_list(pf)
# Fit the main model 
fit_sectors <- m1$sample(
  data = data,
  iter_sampling = 250,
  iter_warmup = 250,
  chains = 8,
  init = list(init_list, init_list, init_list, init_list, init_list, init_list, init_list, init_list),  # one list for each chain
  parallel_chains = 8,
  refresh = 10
)
# Save
saveRDS(fit_sectors, "data/fit_sectors.RDS")


m2 <- cmdstan_model("code/stan_models/model_aggregate.stan")
pf<-m2$pathfinder(data = data, init = 0, num_paths = 4)

# Extract Initalization values
init_list2<-get_init_list(pf)
# Fit the main model 
fit_aggregate <- m2$sample(
  data = data,
  iter_sampling = 250,
  iter_warmup = 250,
  chains = 8,
  init = list(init_list2, init_list2, init_list2, init_list2, init_list2, init_list2, init_list2, init_list2),  # one list for each chain
  parallel_chains = 8,
  refresh = 10
)
# Save
saveRDS(fit_aggregate, "data/fit_aggregate.RDS")



###################################
####### MAKE PLOTS ################
# fit_sectors <-readRDS("data/fit_sectors.RDS")

post <- extract.samples2(fit_sectors)
post2 <- extract.samples2(fit_aggregate)

make_layout <- function(){
  
  layout_matrix <- rbind(
    c(1, 1, 2, 2),    # First row: two plots, each taking half the width
    c(3, 4, 5, 6),    # Second row: four plots
    c(7, 8, 9, 10),   # Third row: four plots
    c(11, 12, 13, 14) # Fourth row: four plots
  )
  
  # Adjust the heights of the rows
  layout_heights <- c(2, 1, 1, 1)  # First row shorter, next rows taller
  
  # Set up the layout
  layout(mat = layout_matrix, heights = layout_heights)
  par(oma = c(2,2,0, 0))
}

font_size = 1.8
lwidth = 3.5
sectors<- c("housing", "retail",  "fishing", "manufacturing", "agriculture",   "cloves",  "ntfp",    "government",    "livestock",     "agroforestry",  "seaweed", "services"  )

pdf("figures/main.pdf", width = 11, height = 14)
make_layout()
par(oma = c(2,2,0, 0))
par(mar=c(5.2,4.2,3,1))
source("code/plotting/base_plot.R")
par(mar=c(3,3,2,1))
source("code/plotting/earnings_on_kesi.R")
dev.off()




