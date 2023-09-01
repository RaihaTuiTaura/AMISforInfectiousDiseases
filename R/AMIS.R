#' Run the AMIS algorithm to fit a transmission model to a map
#' 
#' For details of the algorithm, see 
#' Integrating geostatistical maps and infectious disease transmission models 
#' using adaptive multiple importance sampling.
#' Renata Retkute, Panayiota Touloupou, Maria-Gloria Basanez,
#' T. Deirdre Hollingsworth, Simon E.F. Spencer
#' Ann. Appl. Stat. 15 (4) 1980 - 1998, December 2021.
#' DOI: https://doi.org/10.1214/21-AOAS1486
#'
#' @param prevalence_map An L x M matrix containing samples from the fitted prevalence map, where L is the number of locations and M the number of samples.
#' The location names are inherited from \code{rownames(prevalence_map)} if possible. Alternatively, a list with one entry for each timepoint.
#' Each entry must be a list containing objects \code{data} (an L x M matrix of data as above);
#' and \code{likelihood} a function taking arguments \code{param} (model parameters used in the simulation), \code{data} (a matrix of data as above),
#' \code{prevalence} (a matrix of output from the transmission model) and optional logical \code{log}, which returns the vector of (log)-likelihoods.
#' If a likelihood is not specified then it is assumed that
#' the data consist of samples from a geo-statistical model and empirical methods are used.  
#' @param transmission_model A function taking a vector of n seeds and an n x d matrix of parameter vectors as inputs
#'  and producing a n x timepoints MATRIX of prevalences as output (it must be a matrix even when timepoints==1).
#' @param prior A list containing the functions \code{dprior} and \code{rprior} (density and RNG).
#' Argument to \code{dprior} must be a d-length vector of model parameters and binary indicator \code{log} to indicate whether to calculate log density or not. 
#' Argument to \code{rprior} must be a single integer \code{n} that determines the number of samples to draw. 
#' \code{rprior} must produce an n by d MATRIX of parameters, even when d=1.
#' parameter names are inherited from the \code{colnames} from the output of \code{rprior} if possible.
#' @param amis_params A list containing the control parameters for the AMIS algorithm
#' \describe{
#' \item{\code{delta} the smoothing parameter in the empirical RN derivative (usually 0.01).}
#' \item{\code{nsamples} the number of new samples drawn within each AMIS iteration.}
#' \item{\code{mixture_samples} the number of samples used to represent the weighted parameters in the mixture fitting.}
#' \item{\code{df} the degrees of freedom in the t-distributions, used to yield a heavy tailed proposal.}
#' \item{\code{target_ess} the target effective sample size.}
#' \item{\code{log} logical indicating if calculations are to be performed on log scale.} 
#' \item{\code{max_iters} maximum number of AMIS iterations.}
#' \item{\code{breaks} optional vector specifying the breaks for the histogram. Last entry must be strictly larger than the largest possible prevalence. }
#' \item{\code{bayesian} optional logical for whether to perform Bayesian updated step (does not divide by the induced prior 
#' over the simulated prevalence when calculating weights).}
#' }
#' @param seed Optional seed for the random number generator
#' @param initial_amis_vals optional list containing the object created by setting intermittent_output=TRUE. 
#' Necessary to initialise the algorithm from intermittent output from a previous run (where at least one iteration was successful!)
#' @return A dataframe of the sampled parameters, simulation seed, and weight in each location.
#' @export
amis <- function(prevalence_map, transmission_model, prior, amis_params, seed = NULL, initial_amis_vals = NULL) {

  if(amis_params[["intermittent_output"]]){
    message("Saving output after each iteration (this will increase memory usage). \n")
  }

  save_output <- function(){
    
    if (is.null(initial_amis_vals)){
      allseeds <- 1:(niter * nsamples) 
    } else {
      allseeds <- c(1:(initial_amis_vals$amis_bits$last_simulation_seed + (niter * nsamples))) 
    }
    
    ret <- data.frame(allseeds, param, simulated_prevalences, weight_matrix)
    if (is.null(rownames(prevalence_map[[1]]$data))) {
      iunames<-sapply(1:dim(weight_matrix)[2], function(idx) sprintf("iu%g", idx))
    } else {
      iunames<-rownames(prevalence_map[[1]]$data)
    }
    if (is.null(colnames(param))) {
      paramnames<-sapply(1:dim(param)[2], function(idx) sprintf("param%g", idx))
    } else {
      paramnames<-colnames(param)
    }
    if (is.null(colnames(simulated_prevalences))) {
      prevnames<-sapply(1:dim(simulated_prevalences)[2], function(idx) sprintf("prev%g", idx))
    } else {
      prevnames<-paste0("prev",colnames(simulated_prevalences))
    }
    colnames(ret) <- c("seeds",paramnames,prevnames,iunames)
    
    return(list(results=ret, amis_bits = list(simulated_prevalences=simulated_prevalences, 
                                              likelihoods=likelihoods, 
                                              first_weight=first_weight,
                                              weight_matrix=weight_matrix, 
                                              ess=ess, 
                                              components=components, 
                                              param=param, 
                                              prior_density=prior_density,
                                              last_simulation_seed=max(allseeds))))
  }
  
  # Checks
  if (!is.null(amis_params[["breaks"]])){print("Using empirical weights from user-defined histogram. Ensure last entry is strictly larger than the largest possible prevalence.\n")}
  if (tryCatch(!(is.numeric(amis_params$delta) || is.numeric(amis_params$nsamples || is.numeric(amis_params$mixture_samples) || is.numeric(amis_params$df) || is.numeric(amis_params$target_ess) || is.logical(amis_params$log || is.numeric(amis_params$max_iters)))),error=function(e) return(TRUE))) {stop("Error in arguments for amis_params. Refer to function documentation for correct specification. For defaults set amis_params = default_amis_params(). \n")}
  if (tryCatch(amis_params[["nsamples"]] <= 1,error=function(e) return(TRUE))) {stop("Number of samples must be greater than 1. \n")}  
  if (tryCatch(amis_params[["df"]] < 1,error=function(e) return(TRUE))) {stop("Degrees of freedom must be greater than 0. \n")}
  if (tryCatch(amis_params[["target_ess"]] < 1,error=function(e) return(TRUE))) {stop("Target ESS must be greater than 0. \n")}
  if (tryCatch(amis_params[["max_iters"]] <= 1,error=function(e) return(TRUE))) {stop("Maximum number of iterations must be greater than 1. \n")}
  if (!(is.matrix(prevalence_map) || is.data.frame(prevalence_map) || is.list(prevalence_map))) {stop("prevalence_map must be a matrix or data frame of size #locations by #samples (for one timepoint) or a list of matrices (for >1 timepoints). \n")}
  if (tryCatch(!(is.function(prior$rprior) || is.function(prior$dprior)),error=function(e) return(TRUE))) {stop("prior must be a list containing functions named rprior and dprior. \n")}
  if (!is.function(transmission_model)) {stop("transmission_model must be a function. \n")}
  if (is.list(prevalence_map)) {
    spatial_units = sapply(1:length(prevalence_map), function(t) nrow(prevalence_map[[t]]))
    if (!(length(unique(spatial_units)) == 1)){
      stop("Number of spatial units at each time point in prevalence_map must be equal. If data for some locations are missing at a timepoint set to NA.\n")
    }
  }
  
  # Formatting
  if (is.matrix(prevalence_map) || is.data.frame(prevalence_map)) {prevalence_map=list(list(data=prevalence_map))}
  if (is.null(amis_params[["bayesian"]])) {amis_params[["bayesian"]]<-FALSE}
  
  # Initialise
  if(!is.null(seed)) set.seed(seed)
  nsamples <- amis_params[["nsamples"]]
  
  # Define first_weight object 
  first_weight = rep(1-amis_params[["log"]], nsamples)

  t=1
  if(is.null(initial_amis_vals)){
    cat("AMIS iteration 1\n")
    # Sample first set of parameters from the prior
    param <- prior$rprior(nsamples)
    if(!is.matrix(param)) {stop("rprior function must produce a MATRIX of size #simulations by #parameters, even when #parameters is equal to 1. \n")}
    if(length(c(which(is.na(param)),which(is.nan(param))))>0) {warning("Greater than 1 sample from the prior was an NA or NaN value. \n")}
    if(ncol(param)==1 & amis_params[["bayesian"]]==F) {warning("Currently running with amis_params[['bayesian']]==FALSE. 
                                                                              For models with only one parameter it is recommended to set 
                                                                              amis_params[['bayesian']] = TRUE for prior to influence the weights calculation.\n")}
    # to avoid duplication, evaluate prior density now.
    prior_density<-sapply(1:nsamples,function(b) {prior$dprior(param[b,],log=amis_params[["log"]])})
    if(!(length(prior_density)==nsamples)) {stop("Output from dprior function must have length 1. \n")}
    if(length(c(which(is.na(prior_density)),which(is.nan(prior_density))))>0) {warning("Greater than 1 prior density evaluation was an NA or NaN value. \n")}
    # Simulate from transmission model
    simulated_prevalences <- transmission_model(seeds = 1:nsamples, param)
    if(!is.matrix(simulated_prevalences)) {warning("Unless specifying a bespoke likelihood function, transmission_model function should produce a MATRIX of size #simulations by #timepoints, even when #timepoints is equal to 1. \n")}
    if(!(nrow(param) == nrow(simulated_prevalences))) {warning("Unless specifying a bespoke likelihood function, number of rows in matrices from transmission_model and rprior functions must be equal (#simulations). \n")}
    if(!(length(prevalence_map) == ncol(simulated_prevalences))) {warning("Unless specifying a bespoke likelihood function, number of timepoints in prevalence_map and the number of columns in output from transmission_model function must be equal to #timepoints. \n")}
    if(length(c(which(is.na(simulated_prevalences)),which(is.nan(simulated_prevalences))))>0) {warning("Output from transmission model produced greater than 1 NA or NaN value. \n")}
    
    # to avoid duplication, evaluate likelihood now.
    likelihoods<- compute_likelihood(param,prevalence_map,simulated_prevalences,amis_params)
    if(length(c(which(is.nan(likelihoods))))>0) {warning("Likelihood evaluation produced greater than 1 NaN value. \n")}
  
    # Determine first time each location appears in the data
    n_tims = dim(likelihoods)[1]
    n_locs = dim(likelihoods)[2]
    if (n_locs == 1) {
      lik_matrix <- function(l) as.matrix(l)
    } else {
      lik_matrix <- function(l) t(l)
    }
    locations_first_t = rep(NA,n_locs)
    for (j in 1:n_tims) {
      lik_mat = lik_matrix(likelihoods[j,,])
      for (l in 1:n_locs){
      locations_first_t[which(!is.na(lik_mat[1,]) & is.na(locations_first_t))] = j
      }
    }
     
    weight_matrix <- compute_weight_matrix(
      likelihoods,
      simulated_prevalences,
      amis_params,
      first_weight = first_weight,
      locations_first_t = locations_first_t
    )
    if(length(c(which(is.na(weight_matrix)), which(is.nan(weight_matrix))))>0) {warning("Weight matrix contains greater than 1 NA or NaN value. \n")}
    
    ess <- calculate_ess(weight_matrix,amis_params[["log"]])
    cat("  min ESS:",round(min(ess))," mean ESS:",round(mean(ess))," max ESS:",round(max(ess)),"\n")
    cat(" ",length(which(ess<amis_params[["target_ess"]])),"locations are below the target ESS.\n")
    # Make object to store the components of the AMIS mixture.
    components <- list(
      G = c(0), # number of mixture components from proposal for each iteration (zero is for prior)
      Sigma = list(), # list of covariance matrices for each component
      Mean = list(), # list of means for each component
      probs = list() # probability of each component (unnormalised)
    )
    seeds <- function(t) ((t - 1) * nsamples + 1):(t * nsamples)  #function to calculate the seeds for iteration t.
    niter <- 1 # number of completed iterations 
    if (amis_params[["intermittent_output"]]){
      res = save_output()
      assign("intermittent_output",res,.GlobalEnv)
    }
  } else {     
    cat("AMIS iteration 1\n")
    message("Initialising algorithm from a previous run provided by the user. \n")
    check_initial_vals = function(d){
      if(!is.null(initial_amis_vals$amis_bits[[d]])){
        initial_amis_vals$amis_bits[[d]]
      } else {
        stop(paste0("Cannot find object 'initial_amis_vals$amis_bits[['",d,"']]'. To initialise from a previous run, use object saved by setting amis_params[['intermittent output']]=TRUE.\n"))
      }
    }
    
    simulated_prevalences = check_initial_vals("simulated_prevalences")
    likelihoods= check_initial_vals("likelihoods")
    weight_matrix= check_initial_vals("weight_matrix")
    ess = check_initial_vals("ess")
    components = check_initial_vals("components")
    param = check_initial_vals("param")
    prior_density = check_initial_vals("prior_density")
    first_weight = check_initial_vals("first_weight")
    n_tims = dim(likelihoods)[1]
    n_locs = dim(likelihoods)[2]
    if (n_locs == 1) {
      lik_matrix <- function(l) as.matrix(l)
    } else {
      lik_matrix <- function(l) t(l)
    }
    
    locations_first_t = rep(NA,n_locs)
    for (j in 1:n_tims) {
      lik_mat = lik_matrix(likelihoods[j,,])
      for (l in 1:n_locs){
        locations_first_t[which(!is.na(lik_mat[1,]) & is.na(locations_first_t))] = j
      }
    }
    
    seeds <- function(t) ((t - 2) * nsamples + 1):((t - 1) * nsamples) + initial_amis_vals$amis_bits$last_simulation_seed #function to calculate the seeds for iteration t.
    niter <- 0 # number of completed iterations 
  }
  
  # More checks
  if (length(which(is.na(locations_first_t)))>0){
    warning(paste0(length(which(is.na(locations_first_t))), " location(s) provided with no data. Using prior information to determine weights for these locations.\n"))
  }
  
  # Continue if target_ess not yet reached
  if (min(ess) < amis_params[["target_ess"]]){
    for (t in 2:amis_params[["max_iters"]]) {
      cat("AMIS iteration",t,"\n")
      mean_weights <- update_according_to_ess_value(weight_matrix, ess, amis_params[["target_ess"]],amis_params[["log"]])
      if ((amis_params[["log"]] && max(mean_weights)==-Inf) || (!amis_params[["log"]] && max(mean_weights)==0)) {stop("No weight on any particles for locations in the active set.\n")}
      mixture <- weighted_mixture(param, amis_params[["mixture_samples"]], mean_weights, amis_params[["log"]])
      cat("  A",mixture$G,"component mixture has been fitted.\n")
      components <- update_mixture_components(mixture, components, t)
      new_params <- sample_new_parameters(mixture, nsamples, amis_params[["df"]], prior, amis_params[["log"]])
      if(length(c(which(is.na(new_params$params)), which(is.nan(new_params$params))))>0) {warning("Greater than 1 sample from the proposal after the first iteration of AMIS was an NA or NaN value. \n")}
      param <- rbind(param, new_params$params)
      prior_density <- c(prior_density,new_params$prior_density)
      new_prevalences <- transmission_model(seeds(t), new_params$params)
      if(length(c(which(is.na(new_prevalences)), which(is.nan(simulated_prevalences))))>0) {warning("Output from transmission model produced greater than 1 NA or NaN value. \n")}
      simulated_prevalences <- rbind(simulated_prevalences,new_prevalences)
      likelihoods <- compute_likelihood(new_params$params,prevalence_map,new_prevalences,amis_params,likelihoods)
      if(length(c(which(is.nan(likelihoods))))>0) {warning("Likelihood evaluation produced greater than 1 NaN value. \n")}
      first_weight <- compute_prior_proposal_ratio(components, param, prior_density, amis_params[["df"]], amis_params[["log"]]) # Prior/proposal
      weight_matrix <- compute_weight_matrix(likelihoods, simulated_prevalences, amis_params, first_weight,locations_first_t) # RN derivative (shd take all amis_params)
      if(length(c(which(is.na(weight_matrix)),which(is.nan(weight_matrix))))>0) {warning("Weight matrix contains greater than 1 NA or NaN value. \n")}
      ess <- calculate_ess(weight_matrix,amis_params[["log"]])
      cat("  min ESS:",round(min(ess))," mean ESS:",round(mean(ess))," max ESS:",round(max(ess)),"\n")
      cat(" ",length(which(ess<amis_params[["target_ess"]])),"locations are below the target ESS.\n")
      niter <- niter + 1
      if (amis_params[["intermittent_output"]]){
        res = save_output()
        assign("intermittent_output",res,.GlobalEnv)
      }
      if (min(ess) >= amis_params[["target_ess"]]) break
    }
  }

  if(niter == amis_params[["max_iters"]] && min(ess) < amis_params[["target_ess"]]) {
    msg <- sprintf(
      "Some locations did not reach target ESS (%g) after %g iterations",
      amis_params[["target_ess"]], niter
      )
    warning(msg)
  }
  
  # Save output
  res = save_output()
  
  # Calculate model evidence if using Bayesian updating
  if (amis_params[["bayesian"]]){
    model_evidence = compute_model_evidence(likelihoods, simulated_prevalences, amis_params, first_weight)
    return(list(sample=res$results, evidence=model_evidence))
  } else {
    return(list(sample=res$results,evidence=NULL))
  }
}
