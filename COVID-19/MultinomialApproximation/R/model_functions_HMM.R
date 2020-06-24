# one predictive step
predictive_step_one_particle <- function( t_start, t_end, dt, theta, filtering_prev, simzetaA, travelF ){
  
  # y has to be in the form of a matrix 10x10
  # filtering_prev has to be in the form of a matrix 10x10
  
  pitt = matrix(colSums(filtering_prev), 10, 1)
  
  # pop size
  N = theta[["pop_travel"]]
  
  # rates
  inf_rate <- (simzetaA*sum(pitt[c(4, 5)]))
  inc_rate <- theta[["incubation"]]*2
  rec_rate <- theta[["recover"]]*2
  
  travel_frac <- travelF  
  
  # scale transitions
  P_t = 1 - exp(-dt*inf_rate)
  P_C = 1 - exp(-dt*inc_rate)
  P_R = 1 - exp(-dt*rec_rate)

  # set the transition kernel
  Kappa = matrix(0, 10, 10)
  Kappa[1,1] = 1 - P_t
  Kappa[1,2] = P_t*(1 - travel_frac)
  Kappa[1,6] = P_t*(travel_frac)
  
  Kappa[2,3] = P_C
  Kappa[2,2] = 1 - P_C
  Kappa[3,4] = P_C
  Kappa[3,3] = 1 - P_C
  Kappa[6,7] = P_C
  Kappa[6,6] = 1 - P_C
  Kappa[7,8] = P_C
  Kappa[7,7] = 1 - P_C
  
  Kappa[4,5]  = P_R
  Kappa[4,4]  = 1 - P_R
  Kappa[5,10] = P_R
  Kappa[5,5]  = 1 - P_R
  Kappa[8,9]  = P_R
  Kappa[8,8]  = 1 - P_R
  Kappa[9,10] = P_R
  Kappa[9,9]  = 1 - P_R
  
  Kappa[10,10] =  1
  
  # intermediate steps
  pitt_next = pitt
  for(ii in seq((t_start+dt),t_end,dt) ){
    # rates
    inf_rate <- (simzetaA*sum(pitt_next[c(4, 5)]))
    
    # scale transitions
    P_t = 1 - exp(-dt*inf_rate)
    
    Kappa[1,1] = 1 - P_t
    Kappa[1,2] = P_t*(1 - travel_frac)
    Kappa[1,6] = P_t*(travel_frac)
    
    pitt_expanded_prev   = Kappa*matrix(rep(pitt_next, 10), 10, 10)
    pitt_next = matrix(colSums(pitt_expanded_prev), 10, 1)
        
  }
  
  pitt_expanded_prev 

}

# one filtering step
correction_step_one_particle <- function( theta, pitt_expanded_prev, y_current, Q ){
  
  # pop size
  N = theta[["pop_travel"]]
    
  rho_vec              = pitt_expanded_prev*(1 - Q)
  rho_vec              = rho_vec/sum(rho_vec)
  
  filtering_new = y_current/N + ( 1 - sum(y_current)/N )*rho_vec
  
  filtering_new
}

# one step smoothing
smoothing_step_one_particle <- function( dt, smoothing_t, filtering_tminus1 ){
  
  # y has to be in the form of a matrix 10x10
  # filtering_prev has to be in the form of a matrix 10x10
  smoothing_t_next = smoothing_t
  
  for(ii in seq((0+dt),1,dt) ){
    smoothing_short_tminus1 = rowSums(smoothing_t_next)
  
    L = t(filtering_tminus1)/(colSums(filtering_tminus1) %o% matrix(1, nrow=10))[,,1] 
    L[is.na(L)] = 0

    smoothing_t_next   = (matrix(1, nrow=10) %o% smoothing_short_tminus1 )[,1,]*t(L)
  }
  
  list( smoothing_short = smoothing_short_tminus1, smoothing = smoothing_t_next )
}

# sample from smoothing
smoothing_sample <- function(dt, sample_short_t, filtering_t ){

  L = t(filtering_t)/(colSums(filtering_t) %o% matrix(1, nrow=10))[,,1] 
  L[is.na(L)] = 0
  
  Flag = TRUE
  for(ii in seq((0+dt),1,dt) ){
    
    if(Flag){Flag = FALSE}
    
    else{
      L = L %*% L 
      L[is.na(L)] = 0
    }
    
  } 
  
  multinom_param = (L)
  
  sampling_function = function(i){if(sample_short_t[i]!=0){rmultinom(1, size= sample_short_t[i], prob = multinom_param[i,])}else{rep(0, 10)}}

  sapply(seq(1,10), sampling_function)
}



# SMC function --------------------------------------------

smc_model <- function(theta, nn, dt, q, q_local){
  
  # nn = 100;   dt <- 0.25
  
  # Assumptions - using daily growth rate
  ttotal <- t_period
  t_length <- ttotal
  

  
  #simzeta <- matrix(rlnorm(nn*t_length, mean = -theta[["betavol"]]^2/2, sd = theta[["betavol"]]),ncol=ttotal)
  simzeta_sim <- matrix(rnorm(nn*t_length, mean = 0, sd = theta[["betavol"]]),nrow=ttotal)
  simzeta <- matrix(rnorm(nn*t_length, mean = 0, sd = 0),nrow=ttotal)
  simzeta[1,] <- exp(simzeta_sim[1,])*theta[["beta"]] # define IC
  
  # Latent variables
  beta_traj = matrix(NA,ncol=1,nrow=ttotal)
  filt_traj = array(NA, dim = c(ttotal, 10, 10))
  smoo_traj = array(NA, dim = c(ttotal, 10, 10))
  
  sample_from_smoothing = array(NA, dim = c(ttotal, 10))
  sample_from_full_smoothing = array(NA, dim = c(ttotal-1, 10, 10))
  
  #names = c("S", "E_1W", "E_2W", "I_1W", "I_2W", "E_1T", "E_2T", "I_1T", "I_2T", "R")
  prediction = array(0,dim=c(t_length, 10, 10, nn ))#, dimnames = list(NULL, names, names, NULL ))
  filtering = array(0,dim=c(t_length, 10, 10, nn ))#, dimnames = list(NULL, names, names, NULL ))
  w <- matrix(NA,nrow=nn,ncol=ttotal); w[,1] <- 1  # weights
  W <- matrix(NA,nrow=nn,ncol=ttotal)
  A <- matrix(NA,nrow=nn,ncol=ttotal) # particle parent matrix
  l_sample <- rep(NA,ttotal)
  lik_values <- rep(NA,ttotal)
  
  #initial population and infected
  N = theta[["pop_travel"]] 
  s_init    = theta[["pop_travel"]] - theta[["init_cases"]]
  inf1_init = theta[["init_cases"]]/2
  inf2_init = theta[["init_cases"]]/2

  filtering[1,1,1,] = s_init/N
  filtering[1,4,4,] = inf1_init/N
  filtering[1,5,5,] = inf2_init/N
  
  prediction[1,,,] = filtering[1,,,]
  
  # Iterate through steps
  
  for(tt in 2:ttotal){
    
    # Add geom random walk on transmission 
    simzeta[tt,] <- simzeta[tt-1,]*exp(simzeta_sim[tt,])
    
    # travel restrictions in place?
    if(tt<wuhan_travel_time){travelF <- theta[["travel_frac"]]}else{travelF <- 0}
    
    # run filtering
    y_current = matrix(0, 10, 10)
    
    Q = matrix(0, 10, 10)
    
    if(!is.na(data_list$local_case_data_onset[tt])){
      y_current[3,4] = data_list$local_case_data_onset[tt]
      
      Q[3,4] = q_local
    }
    
    if (!is.na(data_list$int_case_onset[tt])){
      y_current[7,8] = data_list$int_case_onset[tt]
      
      Q[7,8] = q
    }
    
    predictive_per_particle = function(particle_num){ predictive_step_one_particle( tt-1, tt, dt, theta, filtering[tt-1,,,particle_num], simzeta[tt, particle_num], travelF ) }
    prediction[tt,,,] <- array(sapply(seq(1,nn,1), predictive_per_particle), dim = c(10, 10, nn))
    
    correction_per_particle = function(particle_num){ correction_step_one_particle( theta, prediction[tt,,,particle_num], y_current, Q ) }
    filtering[tt,,,] <- array(sapply(seq(1,nn,1), correction_per_particle), dim = c(10, 10, nn))

    AssignWeights = function(particle_num){AssignWeights_per_particle( N, y_current, Q, prediction[tt,,,particle_num])}
    
    # calculate weights
    w[,tt] <- sapply( seq(1,nn,1), AssignWeights )
    
    
    # check likelihood isn't NA
    if(is.na(max(w[1:nn,tt])) | max(w[1:nn,tt]) == 0){
      likelihood0 = -Inf
      return(list(lik=likelihood0 ))
    }
    
    # normalise particle weights
    sum_weights <- sum(w[1:nn,tt])
    W[1:nn,tt] <- w[1:nn,tt]/sum_weights
    
    # resample particles by sampling parent particles according to weights:
    A[, tt] <- sample(1:nn,prob = W[1:nn,tt],replace = T)
    
    # DEPRECATED
    # for (j in 1:nn){
    #   locs <- pickA[cumsum_W >= rand_vals[j]]; A[j, tt] <- locs[1]
    # }
    
    # Resample particles for corresponding variables
    simzeta[tt,] <- simzeta[tt, A[, tt]] #- needed for random walk on beta
    prediction[tt,,,] = prediction[tt,,,A[, tt]]
    filtering[tt,,,]  = filtering[tt,,,A[, tt]]
    
    
  } # END PARTICLE LOOP
  
  
  # Estimate likelihood:
  for(tt in 1:ttotal){
    
    lik_values[tt] = log(sum(w[1:nn,tt])) # log-likelihoods
  }
  
  likelihood0 = -ttotal*log(nn)+ sum(lik_values) # add full averaged log-likelihoods
  
  # Sample latent variables:
  locs <-  sample(1:nn,prob = W[1:nn,tt],replace = T)
  l_sample[ttotal] <- locs[1]
  filt_traj[ttotal,,] <- filtering[ttotal,,,l_sample[ttotal]] 
  beta_traj[ttotal,] <- simzeta[ttotal,l_sample[ttotal]]
  
  smoo_traj[ttotal,,] <- filtering[ttotal,,,l_sample[ttotal]]

  smoothing_short_t = matrix(colSums(smoo_traj[ttotal,,]), 10, 1)
  
  sample_from_smoothing[ttotal,] = rmultinom(1, N, smoothing_short_t)

  # this recover the beta trajectory of a given particle. This is like sampling the traj from
  # the full posterior 
  # Note given beta the filtering algorithm is deterministic, hence we can recover it 
  # from the ancestors as well
  for(ii in seq(ttotal,2,-1)){
    #print(ii)
    l_sample[ii-1] <- A[l_sample[ii],ii] # have updated indexing
    filt_traj[ii-1,,] <- filtering[ii-1,,,l_sample[ii-1]]
    beta_traj[ii-1,] <- simzeta[ii-1,l_sample[ii-1]]
    
    if(ii<wuhan_travel_time){travelF <- theta[["travel_frac"]]} else{travelF <- 0}
    
    smoothing_step = smoothing_step_one_particle(dt, smoo_traj[ii,,], filt_traj[ii-1,,] )
    smoothing_short_t = smoothing_step$smoothing_short
    smoo_traj[ii-1,,] = smoothing_step$smoothing
    
    sample_from_full_smoothing[ii-1,,] = smoothing_sample(dt, sample_from_smoothing[ii,], filt_traj[ii-1,,] )
    sample_from_smoothing[ii-1,] = rowSums(sample_from_full_smoothing[ii-1,,])
 
  }

  ESS = (colSums(w[1:nn,]))^2/(colSums(w[1:nn,]^2))
  
  # just smoothing here with the final beta

  return(list(ESS = ESS, filtering_trace=filt_traj, smoothing_trace=smoo_traj, beta_trace=beta_traj, sample_full_smoothing = sample_from_full_smoothing, sample_smoothing = sample_from_smoothing, lik=likelihood0 ))
  
}


# Likelihood calc for SMC --------------------------------------------

AssignWeights_per_particle <- function( N, y_current, Q, prediction_current){
  
  if(sum(prediction_current<0)==0){
    
    multinomial_prop = c( array(prediction_current*Q, dim = 100), 1-sum(prediction_current*Q))
    
    sample = c( array(y_current, dim = 100), N-sum(y_current))
    
    dmultinom(sample, size = NULL, multinomial_prop, log =FALSE)    
  }
  else{0}

}



# MLE grid search -  2D ---------------------------------------------------------
# grid search on q and q_local
MLE_check_2D <- function(q_local_grid, q_grid, nn, dt,
                         filename = 1){
  
  store_lik <- NULL
  
  for(ii in 1:length(q_local_grid)){
    print(ii)
    for(jj in 1:length(q_grid)){
      
      # Run SMC and output likelihooda
      q       = q_grid[jj]
      q_local = q_local_grid[ii]
      output_smc <- smc_model(theta, nn, dt, q, q_local)
      
      store_lik <- rbind(store_lik,c(q, q_local,output_smc$lik))
      
    } 
  }
  
  colnames(store_lik) <- c("q", "q_local", "lik")
  store_lik <- as_tibble(store_lik)
  
  write_csv(store_lik, paste0("outputs/param_search_",filename,".csv"))
  store_lik
  
}

