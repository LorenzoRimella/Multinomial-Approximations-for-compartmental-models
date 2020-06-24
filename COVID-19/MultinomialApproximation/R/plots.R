# Run SMC to get bootstrap estimates --------------------------------------------

run_fits <- function(rep_plot,nn,cut_off, dt, q, q_local, filename="1"){
  
  out_rep <- foreach(kk = 1:rep_plot) %dopar% {
    print(kk)
    output_smc <- smc_model(theta,nn,dt, q, q_local)
    output_smc
  }
  
  R0_plot = matrix(NA,ncol=rep_plot,nrow=t_period)
  smoothing_plot = array(NA, dim = c(t_period, 10, 10, rep_plot))
  ESS_plot = matrix(NA,ncol=rep_plot,nrow=t_period)
  lik_plot = matrix(NA,ncol=rep_plot,nrow=t_period)
  sample_smoothing_plot = array(NA, dim = c(t_period-1, 10, 10, rep_plot))
  
  for(kk in 1:rep_plot){
    output_smc <- out_rep[[kk]]
    if(output_smc$lik != - Inf){
      
      smoothing_plot[,,,kk] = output_smc$smoothing_trace
      R0_plot[,kk] <- output_smc$beta_trace/(theta[["recover"]])
      ESS_plot[,kk] <- output_smc$ESS
      lik_plot[,kk] <- output_smc$lik
      sample_smoothing_plot[,,,kk] <- output_smc$sample_full_smoothing
    }
  }
  
  save(
    smoothing_plot,
    R0_plot,
    ESS_plot,
    lik_plot,
    sample_smoothing_plot,
    file=paste0("outputs/bootstrap_fit_",filename,".RData")) 
  
}

# Plot outputs from SMC --------------------------------------------


plot_outputs <- function(rep_plot, cut_off, filename="1"){
  
  load(paste0("outputs/bootstrap_fit_",filename,".RData"))
  
  par(mfrow=c(1,1))
  plot(ESS_plot[,1], type ="l", ylim =  c(min(ESS_plot),3000), xlab = "Time step", ylab = "ESS")
  for(i in 2:rep_plot){lines(ESS_plot[,i])}

  # Remove NA fits
  if(sum(is.na(smoothing_plot[t_period,,,]))!=0){print("There are NA's in the filtering")}

  # Compute the evolving in time reproduction number and its quantiles
  R0_plot = R0_plot[,!is.na(R0_plot[t_period,])]
  
  # Calculate quantiles R_0
  R0_quantile <- apply(R0_plot,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975))}) 
  
  # Remove final few points (as estimation less reliable)
  R0_quantileA <- R0_quantile[,1:(ncol(R0_quantile)-cut_off)]
  date_rangeA <- date_range[1:(length(date_range)-cut_off)]
  
  # forecast window
  date_rangeF <- date_range[date_range>max(cases_Wuhan$date)]
  yyF <- rep(0,length(date_rangeF))
  
  # - - - - - - - 
  # - - - - - - - 
  # Plot the reproduction number, the new cases onset locally and internationally
  # the reported cases locally and internationally
  par(mar=c(2,3,1,1),mgp=c(2,0.55,0)) #mfrow=c(4,2),
  layout(matrix(c(1,1,1,1,2,2,3,3,4,4,5,5), 3, 4, byrow = TRUE))
  
  # Plot outputs
  a_col <- 0.4 # alpha
  xMin1 <- as.Date("2019-12-15") #min(as.Date("2019-12-01")) 
  xMin <- xMin1 #min(as.Date("2020-01-01"))
  xMax <- end_date-1 #max(date_range)
  
  
  # Plot reproduction number
  xMaxR <- as.Date("2020-02-05")
  date_rangeB <- date_rangeA#[date_rangeA>as.Date("2019-12-15")]
  R0_quantileB <- R0_quantileA#[,date_rangeA>as.Date("2019-12-15")]
  xMax1 <- xMax #as.Date("2020-02-01") #xMax #xMax #
  
  plot(date_rangeB,R0_quantileB[1,],col="white",ylim=c(0,8),xlim=c(xMin1,xMaxR),xlab="",ylab=expression(paste(R[t])))
  polygon(c(date_rangeF,rev(date_rangeF)),c(yyF,rev(yyF+1e5)),lty=0,col=rgb(0.9,0.9,0.9))
  
  polygon(c(date_rangeB,rev(date_rangeB)),c(R0_quantileB[2,],rev(R0_quantileB[4,])),lty=0,col=rgb(0,0.3,1,0.35))
  polygon(c(date_rangeB,rev(date_rangeB)),c(R0_quantileB[1,],rev(R0_quantileB[5,])),lty=0,col=rgb(0,0.3,1,0.2))
  lines(date_rangeB,R0_quantileB[3,],type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  lines(date_rangeB,1+0*R0_quantileB[3,],lty=2)
  
  lines(c(wuhan_travel_restrictions,wuhan_travel_restrictions),c(0,10),col="red")
  
  # - - - - - - - 
  # Compute new local infectious onset and its quantiles
  sample_population_from_smoothing = function(p_i){rbinom(1, size = theta[["pop_travel"]], prob = p_i)}
  smoothing_population = matrix(sapply(smoothing_plot[,3,4,], sample_population_from_smoothing), nrow=82)
  
  samplebinom = function(N_i){rbinom(1, size = N_i, prob = q_local)}
  Case_local_onsetA = matrix(sapply(smoothing_population, samplebinom), nrow=82)
  
  Case_local_quantile_onsetA = apply(Case_local_onsetA,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975))})
  
  # Plot local case onsets
  ym1 <- 25
  plot(date_rangeA,Case_local_quantile_onsetA[1,],col="white",ylim=c(0,ym1),xlim=c(xMin1,xMax),xlab="",ylab="New onsets in Wuhan")
  polygon(c(date_rangeF,rev(date_rangeF)),c(yyF,rev(yyF+1e5)),lty=0,col=rgb(0.9,0.9,0.9))
  
  polygon(c(date_rangeA,rev(date_rangeA)),c(Case_local_quantile_onsetA[2,],rev(Case_local_quantile_onsetA[4,])),lty=0,col=rgb(0,0.3,1,0.35))
  polygon(c(date_rangeA,rev(date_rangeA)),c(Case_local_quantile_onsetA[1,],rev(Case_local_quantile_onsetA[5,])),lty=0,col=rgb(0,0.3,1,0.2))
  lines(date_rangeA,Case_local_quantile_onsetA[3,],type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  
  lines(c(wuhan_travel_restrictions,wuhan_travel_restrictions),c(0,1e6),col="red")

  legend(date_rangeA[23], 25, legend=c("Data", "Posterior mean", "50% credible interval", "95% credible interval"), 
         box.lwd = 1, box.col = "grey",bg = "white",
         col=c("black", "blue", rgb(0,0.3,1,0.65), rgb(0,0.3,1,0.2)), lty= c(NA, 1, 1, 1), lwd = c(NA, 2, 8, 10), pch=c(19, NA, NA, NA),
         seg.len=0.25, y.intersp=0.65, x.intersp=0.25, cex=0.8,
         text.width = 15)
  
  points(case_data_wuhan$date,case_data_wuhan$number,pch=19,cex=0.8)
  points(case_data_china$date,case_data_china$number,pch=19,cex=0.8)
  
  
  # - - -
  # Compute new cases onset internationally
  sample_population_from_smoothing = function(p_i){rbinom(1, size = theta[["pop_travel"]], prob = p_i)}
  smoothing_population = matrix(sapply(smoothing_plot[,7,8,], sample_population_from_smoothing), nrow=82)

  samplebinom = function(N_i){rbinom(1, size = N_i, prob = q)}
  Case_A = matrix(sapply(smoothing_population, samplebinom), nrow=82)
  
  Case_quantile = apply(Case_A,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975))})
  
  # Plot international cases onsets
  ym1 <- 10
  plot(date_range,case_time,pch=19,ylim=c(0,ym1),xlim=c(xMin1,xMax),ylab="New international onsets",col="white")
  polygon(c(date_rangeF,rev(date_rangeF)),c(yyF,rev(yyF+1e5)),lty=0,col=rgb(0.9,0.9,0.9))
  
  polygon(c(date_range,rev(date_range)),c(Case_quantile[2,],rev(Case_quantile[4,])),lty=0,col=rgb(0,0.3,1,0.35))
  polygon(c(date_range,rev(date_range)),c(Case_quantile[1,],rev(Case_quantile[5,])),lty=0,col=rgb(0,0.3,1,0.2))
  lines(date_range,Case_quantile[3,],type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  points(date_range,case_data_onset_time,pch=19)
  
  legend(date_rangeA[23], 10, legend=c("Data", "Posterior mean", "50% credible interval", "95% credible interval"), 
         box.lwd = 1, box.col = "grey",bg = "white",
         col=c("black", "blue", rgb(0,0.3,1,0.65), rgb(0,0.3,1,0.2)), lty= c(NA, 1, 1, 1), lwd = c(NA, 2, 8, 10), pch=c(19, NA, NA, NA),
         seg.len=0.25, y.intersp=0.65, x.intersp=0.25, cex=0.8,
         text.width = 15)
 
  lines(c(wuhan_travel_restrictions,wuhan_travel_restrictions),c(0,1e6),col="red")

  # - - - - - - - 
  # Compute the new reported cases locally
  # use the estimates from kucharski
  sample_new_I_plot = array(0, dim =c(t_period, rep_plot))
  sample_new_I_plot[2:t_period,] = sample_smoothing_plot[,3,4,]
  
  report_rec__rate = exp(-theta[["report_local"]]*theta[["recover"]])
  report_rec_prob = 1 - exp(- report_rec__rate)
  
  sample_new_Q = array(rbinom(t_period*rep_plot, size = sample_new_I_plot, prob = array(report_rec_prob, dim(sample_new_I_plot))), dim =c(t_period, rep_plot))
  
  report_rate = theta[["report_local"]]
  report_prob = 1 - exp(- report_rate)
  Q = array(0, dim = c(t_period, rep_plot))
  C = array(0, dim = c(t_period, rep_plot))
  for(tt in 2:t_period){
    new_C = rbinom(rep_plot, Q[tt-1,], rep(report_prob, rep_plot))
    C[tt,] = new_C
    Q[tt,] = Q[tt-1,] + sample_new_Q[tt, ] - new_C
  }
   
  inf_local_quantile = apply(C,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975))})

  # Plot predicted reported cases locally
  ym1 <- 10000
  plot(date_range,inf_local_quantile[1,],col="white",ylim=c(0,ym1),xlim=c(xMin1,xMax),xlab="",ylab="New cases in Wuhan")
  
  polygon(c(date_rangeF,rev(date_rangeF)),c(yyF,rev(yyF+1e5)),lty=0,col=rgb(0.9,0.9,0.9))
  
  polygon(c(date_range,rev(date_range)),c(inf_local_quantile[2,],rev(inf_local_quantile[4,])),lty=0,col=rgb(0,0.3,1,0.35))
  polygon(c(date_range,rev(date_range)),c(inf_local_quantile[1,],rev(inf_local_quantile[5,])),lty=0,col=rgb(0,0.3,1,0.2))
  lines(date_range,inf_local_quantile[3,],type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")
 
  # Plot data on a different scale
  par(new=TRUE)
  ym1a <- 5000
  plot(cases_Wuhan$date,cases_Wuhan$new_case,pch=19,xaxt="n",xlim=c(xMin1,xMax),bty="l",yaxt="n",xaxt="n",bty="l",yaxt="n",xlab="",ylab="",ylim=c(0,ym1a))
  axis(4)
  
  mtext("confirmed", side=4, cex=0.7,line=-1)
  
  lines(c(wuhan_travel_restrictions,wuhan_travel_restrictions),c(0,1e6),col="red")
  
  legend(date_rangeA[23], 5000, legend=c("Data (out-of-sample)", "Posterior mean", "50% credible interval", "95% credible interval"), 
         box.lwd = 1, box.col = "grey",bg = "white",
         col=c("black", "blue", rgb(0,0.3,1,0.65), rgb(0,0.3,1,0.2)), lty= c(NA, 1, 1, 1), lwd = c(NA, 2, 8, 10), pch=c(19, NA, NA, NA),
         seg.len=0.25, y.intersp=0.65, x.intersp=0.25, cex=0.8,
         text.width = 15)
  

  # - - - - - - - 
  # Compute the new reported cases internationally
  sample_new_I_plot = array(0, dim =c(t_period, rep_plot))
  sample_new_I_plot[2:t_period,] = sample_smoothing_plot[,7,8,]
  
  report_rec__rate = exp(-theta[["report"]]*theta[["recover"]])
  report_rec_prob = 1 - exp(- report_rec__rate)
  
  sample_new_Q = array(rbinom(t_period*rep_plot, size = sample_new_I_plot, prob = array(report_rec_prob, dim(sample_new_I_plot))), dim =c(t_period, rep_plot))
  
  report_rate = theta[["report"]]
  report_prob = q
  Q = array(0, dim = c(t_period, rep_plot))
  C = array(0, dim = c(t_period, rep_plot))
  for(tt in 2:t_period){
    new_C = rbinom(rep_plot, Q[tt-1,], rep(report_prob, rep_plot))
    C[tt,] = new_C
    Q[tt,] = Q[tt-1,] + sample_new_Q[tt, ] - new_C
  }
  
  inf_quantile = apply(C,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975))})

  # Plot international cases confirmed
  ym1 <- 12
  plot(date_range,case_time,pch=19,ylim=c(0,ym1),xlim=c(xMin1,xMax),ylab="",col="white")
  polygon(c(date_rangeF,rev(date_rangeF)),c(yyF,rev(yyF+1e5)),lty=0,col=rgb(0.9,0.9,0.9))
  
  polygon(c(date_range,rev(date_range)),c(inf_quantile[2,],rev(inf_quantile[4,])),lty=0,col=rgb(0,0.3,1,0.35))
  polygon(c(date_range,rev(date_range)),c(inf_quantile[1,],rev(inf_quantile[5,])),lty=0,col=rgb(0,0.3,1,0.2))
  lines(date_range,inf_quantile[3,],type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  points(date_range,case_time,pch=19)

  title(ylab="New international exports confirmed", line=1.5, cex.lab=1)
  
  legend(date_rangeA[23], 12, legend=c("Data (out-of-sample)", "Posterior mean", "50% credible interval", "95% credible interval"), 
         box.lwd = 1, box.col = "grey",bg = "white",
         col=c("black", "blue", rgb(0,0.3,1,0.65), rgb(0,0.3,1,0.2)), lty= c(NA, 1, 1, 1), lwd = c(NA, 2, 8, 10), pch=c(19, NA, NA, NA),
         seg.len=0.25, y.intersp=0.65, x.intersp=0.25, cex=0.8,
         text.width = 15) 

  
  lines(c(wuhan_travel_restrictions,wuhan_travel_restrictions),c(0,1e6),col="red")

  # - - - - - - - 
  # save the plot
  # dev.copy(pdf,paste("../stoch_model_V2_paper_HMM/Figure_2.pdf",sep=""),width=8,height=8)
  # dev.off()
  
  
}
