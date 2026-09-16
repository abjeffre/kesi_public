###################################################################
############ SHARED HELPERS FOR THE KESI REPLICATION ##############
# Sourced by every MAKE script after its library() calls.

# Test mode: KESI_TEST=1 runs every model with a handful of iterations so the
# whole pipeline can be exercised in minutes. Unset it for the paper's settings.
KESI_TEST <- nzchar(Sys.getenv("KESI_TEST"))
ITER_FULL <- 500
ITER_TEST <- 20
SAMPLER_ITER <- if (KESI_TEST) ITER_TEST else ITER_FULL
# One gradient of the main model costs about 0.5 s, so unadapted warmup at the
# default tree depth of 10 runs for hours. Test mode caps the tree depth.
MAX_TREEDEPTH_FULL <- 10
MAX_TREEDEPTH_TEST <- 5
SAMPLER_MAX_TREEDEPTH <- if (KESI_TEST) MAX_TREEDEPTH_TEST else MAX_TREEDEPTH_FULL
SAMPLER_CHAINS <- 4
SAMPLER_REFRESH <- 100
TEST_SWEEP_CELLS <- 4

# figures/ is not in git; every R step writes a PDF there.
dir.create("figures", showWarnings = FALSE)

`%||%` <- function(a, b) if (!is.null(a)) a else b

###################################################################
############ POSTERIOR HELPERS ####################################

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

# Posterior means of a pathfinder (or any draws object), shaped as a Stan init list.
get_init_list <- function(pf) {
  stanfit <- posterior::as_draws_rvars(pf)
  init_list <- list()
  cnt <- 1
  for (i in names(stanfit)) {
    means <- posterior::summarise_draws(posterior::as_draws_array(stanfit[[i]]))
    temp <- create_empty_like(stanfit[[i]])
    for (j in 1:length(means$mean)) {
      temp[j] <- c(means$mean[j])
    }
    init_list[[cnt]] <- as.array(temp)
    cnt <- cnt + 1
  }
  names(init_list) <- names(stanfit)
  for (i in names(init_list)) {
    if (length(dim(init_list[[i]])) == 1 & length(init_list[[i]]) == 1) init_list[[i]] <- init_list[[i]][1]
  }
  return(init_list)
}

extract.samples2 <- function(x) {
  output <- list()
  stanfit <- posterior::as_draws_rvars(x)
  for (i in names(stanfit)) {
    output[[i]] <- posterior::draws_of(stanfit[[i]])
  }
  return(output)
}

# Pathfinder-based inits with a fallback to init = 0 when pathfinder fails.
sample_with_pathfinder_inits <- function(model, data, refresh) {
  init <- tryCatch({
    init_list <- get_init_list(model$pathfinder(data = data))
    rep(list(init_list), SAMPLER_CHAINS)
  }, error = function(e) {
    message("Pathfinder failed (", conditionMessage(e), "); sampling from init = 0.")
    0
  })
  model$sample(parallel_chains = SAMPLER_CHAINS,
               chains = SAMPLER_CHAINS,
               data = data,
               init = init,
               iter_warmup = SAMPLER_ITER, iter_sampling = SAMPLER_ITER,
               max_treedepth = SAMPLER_MAX_TREEDEPTH, refresh = refresh)
}

###################################################################
############ GAUSSIAN PROCESS DISCRETISATION ######################

# Bin a continuous predictor into L categories and return the bin index,
# bin means and break points; the means feed makeGPDmat.
makeGPCat <- function(x, y = NULL, L, min_buffer = 0, max_buffer = 0) {
  breaks <- seq(min(min(y), min(x)) - min_buffer, max(max(x), max(y)) + min_buffer, length.out = L + 1)
  means <- rep(NA, L)
  for (i in 1:(L)) {
    means[i] <- mean(breaks[i], breaks[i + 1])
  }
  cat <- cut(x, breaks = breaks, right = FALSE)
  return(list(as.integer(cat), means, breaks))
}

makeGPDmat <- function(means) {
  L <- length(means)
  Dmat <- matrix(NA, L, L)
  for (i in 1:L) {
    for (j in 1:L) {
      Dmat[i, j] <- abs(means[i] - means[j])
    }
  }
  return(Dmat)
}

###################################################################
############ PARALLEL SWEEP RUNNER ################################

# Runs one_cell(row_index) over the sweep rows on a PSOCK cluster, which works
# on Windows as well as Unix. KESI_CORES caps the workers; KESI_TEST keeps the
# first few cells only. A cell whose fit fails is reported and skipped, so the
# heatmaps are drawn from the cells that finished.
# Compiles each Stan file once on the master before workers start; parallel
# workers compiling the same file concurrently corrupt each other's build.
compile_stan_models <- function(stan_files) {
  invisible(lapply(stan_files, cmdstanr::cmdstan_model))
}

report_failed_cells <- function(results) {
  failed <- vapply(results, inherits, logical(1), what = "sweep_cell_failure")
  for (result in results[failed]) message(result$message)
  results
}

run_sweep_cells <- function(sweep_list, one_cell) {
  rows <- seq_len(nrow(sweep_list))
  if (KESI_TEST) rows <- head(rows, TEST_SWEEP_CELLS)
  guarded_cell <- function(row) {
    tryCatch(one_cell(row), error = function(e) {
      structure(list(row = row, message = paste0("Sweep cell ", row, " failed: ", conditionMessage(e))),
                class = "sweep_cell_failure")
    })
  }
  cores_requested <- as.integer(Sys.getenv("KESI_CORES", unset = max(1L, parallel::detectCores() - 1L)))
  workers <- max(1L, min(cores_requested, length(rows)))
  if (workers == 1L) return(report_failed_cells(lapply(rows, guarded_cell)))
  cluster <- parallel::makeCluster(workers, rscript_args = "--vanilla")
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  parallel::clusterEvalQ(cluster, {
    library(readr)
    library(posterior)
    library(cmdstanr)
  })
  parallel::clusterExport(cluster, c("sweep_list", "one_cell", "observed_years"), envir = parent.frame())
  parallel::clusterCall(cluster, function() source("code/functions/utility.R"))
  report_failed_cells(parallel::parLapply(cluster, rows, guarded_cell))
}
