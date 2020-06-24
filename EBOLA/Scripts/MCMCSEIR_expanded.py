import numpy as np
from scipy.stats import gamma as RVgamma
from scipy.stats import beta as RVbeta
from scipy.stats import dirichlet as RVdiri
from scipy.stats import norm as RVnorm
from scipy.stats import foldnorm as RVfoldnorm
from scipy.stats import truncnorm as RVtruncnorm

# the gamma distribution consider a varying shape parameter and a scale parameter equal to 1

from Scripts.HMMSEIR_expanded import *

class MCMC_SEIR_expanded:

    def __init__( self, N, t_star, beta_prior, rho_prior, gamma_prior, q_prior,
                        q_r_prior, eta_zero_prior,
                        beta_proposal_variance, rho_proposal_variance, gamma_proposal_variance, q_proposal_variance,
                        q_r_proposal_variance):

        self.N      = N
        self.t_star = t_star

        self.beta_prior = beta_prior
        self.rho_prior = rho_prior
        self.gamma_prior = gamma_prior
        self.q_prior = q_prior
        self.eta_zero_prior = eta_zero_prior
        self.q_r_prior = q_r_prior

        self.beta_proposal_variance  = beta_proposal_variance
        self.rho_proposal_variance   = rho_proposal_variance
        self.gamma_proposal_variance = gamma_proposal_variance

        self.q_proposal_variance   = q_proposal_variance
        self.q_r_proposal_variance = q_r_proposal_variance


    # we choose a gamma prior on all the parameter:
    # Gamma prior for beta,gamma,rho
    # Beta prior for the non zero element of q
    # Dirichlet prior for eta_0

    # priors
    def prior_ratio_beta(self, beta_proposed): 

        beta_ratio  = np.log(( RVgamma.pdf( beta_proposed,   self.beta_prior[0] , 0, 1/self.beta_prior[1] ) ))

        return beta_ratio


    def prior_ratio_rho(self, rho_proposed): 

        rho_ratio = np.log( RVgamma.pdf( rho_proposed,    self.rho_prior[0], 0, 1/self.rho_prior[1]  ) ) 

        return rho_ratio


    def prior_ratio_gamma(self, gamma_proposed): 

        gamma_ratio = np.log( RVgamma.pdf( gamma_proposed,  self.gamma_prior[0], 0, 1/self.gamma_prior[1]) ) 

        return gamma_ratio


    def prior_ratio_q_r(self, q_r_proposed): 

        q_r_ratio = np.log( RVgamma.pdf( q_r_proposed,  self.q_r_prior[0], 0, 1/self.q_r_prior[1]) ) 

        return q_r_ratio


    def prior_ratio_q(self, q_proposed): 

        q_ratio_12 = np.log( RVbeta.pdf( q_proposed[1,2], self.q_prior[0], self.q_prior[1] ) )
        q_ratio_23 = np.log( RVbeta.pdf( q_proposed[2,3], self.q_prior[0], self.q_prior[1] ) )

        return q_ratio_12+q_ratio_23

    # likelihood
    def likelihood_ratio(self, y,  beta_proposed, rho_proposed, gamma_proposed, q_proposed, eta_proposed, q_r_proposed ):

        T  = np.size( y, 2 )
        C1 = np.size( y, 0 )
        C2 = np.size( y, 1 )

        HMM1 = HMM_approxSEIR_expanded( self.N, beta_proposed, rho_proposed, gamma_proposed, q_proposed, eta_proposed, q_r_proposed, self.t_star )

        _, _, _, pitt_prev1 = HMM1.filtering(y)

        # use log to avoid over/under flow problems
        # posterior predictive on the observed compartments
        y_from1               = y[:, : , 1:]
        y_from1_index_nozeros = (y_from1!=0)

        pitt_prev_q = q_proposed.reshape(4,4,1)*pitt_prev1
        pitt_prev_q_index_nozeros = (pitt_prev_q[:,:,1:]!=0)

        non_zero_index = y_from1_index_nozeros*pitt_prev_q_index_nozeros

        likelihood_observed_compartment = np.sum( (y_from1[non_zero_index])*np.log( pitt_prev_q[:, :, 1:][non_zero_index] ) )

        # posterior predictive on the unobserved compartments
        pitt_prev_q_complement = ( 1 - np.sum( np.sum( pitt_prev_q[:, :, 1:], 0 ), 0 ) )
        likelihood_unobserved_compartment = np.sum( np.log( pitt_prev_q_complement**( self.N - np.sum( np.sum( y_from1, 0), 0 ) ) ) ) 

        likelihood = likelihood_observed_compartment + likelihood_unobserved_compartment
        return likelihood

    # Proposal distributions
    def propose_parameters_beta(self, beta_old ):

        beta_proposed  = np.random.normal( beta_old,  self.beta_proposal_variance  )

        return beta_proposed


    def propose_parameters_rho(self, rho_old ):

        rho_proposed   = np.random.normal( rho_old,   self.rho_proposal_variance   )

        return rho_proposed


    def propose_parameters_gamma(self, gamma_old ):

        gamma_proposed = np.random.normal( gamma_old, self.gamma_proposal_variance )

        return gamma_proposed


    def propose_parameters_q_r(self, q_r_old ):

        q_r_proposed = np.random.normal( q_r_old, self.q_r_proposal_variance )

        return q_r_proposed


    def propose_parameters_q(self, q_old ):

        q_proposed     = np.zeros((4,4))

        q_proposed[1, 2]  = np.random.normal( q_old[1, 2], self.q_proposal_variance )
        q_proposed[2, 3]  = np.random.normal( q_old[2, 3], self.q_proposal_variance )

        return q_proposed


    # mcmc sampling
    def posterior_sample(self, sample_size, y, beta_0, rho_0, gamma_0, eta_zero_0, x_0, q_0, q_r_0 ):

#         ############################################################################
#         title = "Check-Chain-expected-prop-LearnQ-"
#         title = title+"Prior"+str(self.beta_prior[0])+"-"+str(self.beta_prior[1])
#         title = title+"-Proposal"+"-"+str(self.beta_proposal_variance)+"-"+str(self.rho_proposal_variance)+"-"+str(self.gamma_proposal_variance)+"-"+str(self.q_r_proposal_variance)+"-"+str(self.q_proposal_variance)+".txt"       
#         f= open(title,"a")
#         f.close()
#         ############################################################################

        # initialization of the mcmc
        beta_chain  = np.zeros((sample_size))
        rho_chain   = np.zeros((sample_size))
        gamma_chain = np.zeros((sample_size))
        q_r_chain   = np.zeros((sample_size))

        eta_zero_chain = np.zeros((4, sample_size))
        x_chain        = np.zeros((4, sample_size))

        q_chain     = np.zeros(( 4, 4, sample_size))

        # initial conditions
        beta_chain[0]   = beta_0
        rho_chain[0]    = rho_0
        gamma_chain[0]  = gamma_0
        q_r_chain[0]    = q_r_0

        eta_zero_chain[:, 0] = eta_zero_0
        x_chain[:,0]         = x_0

        q_chain[:, :, 0]        = q_0

        # initialization of the acceptance ratio
        beta_acceptance_ratio  = 0
        rho_acceptance_ratio   = 0
        gamma_acceptance_ratio = 0
        q_acceptance_ratio     = 0
        eta_acceptance_ratio   = 0
        q_r_acceptance_ratio   = 0

        ####################################################################
        # FIRST ITERATION
        ####################################################################
        iter = 1
        # BETA
        prior_ratio_beta_old    = self.prior_ratio_beta( beta_chain[iter-1] )
        likelihood_ratio_old = self.likelihood_ratio( y,  beta_chain[iter-1], rho_chain[iter-1], gamma_chain[iter-1], q_chain[:, :, iter-1], eta_zero_chain[:, iter-1], q_r_chain[iter-1] )

        beta_proposed = self.propose_parameters_beta( beta_chain[iter-1] )
        if beta_proposed>0:

            prior_ratio_beta_new    = self.prior_ratio_beta( beta_proposed )
            
            likelihood_ratio_beta_new = self.likelihood_ratio( y,  beta_proposed,      rho_chain[iter-1], gamma_chain[iter-1], q_chain[:, :, iter-1], eta_zero_chain[:, iter-1], q_r_chain[iter-1] )      

            ratio_beta = ( likelihood_ratio_beta_new + prior_ratio_beta_new ) - ( likelihood_ratio_old + prior_ratio_beta_old )

            if (ratio_beta) >= np.log(np.random.uniform(0,1)):
                beta_chain[iter]      = beta_proposed
                beta_acceptance_ratio = beta_acceptance_ratio + 1

                # proposal_ratio_beta_old   = proposal_ratio_beta_new
                prior_ratio_beta_old      = prior_ratio_beta_new

                likelihood_ratio_old = likelihood_ratio_beta_new

            else:
                beta_chain[iter] = beta_chain[iter-1]
        else:
            beta_chain[iter] = beta_chain[iter-1]
            

        #RHO
        prior_ratio_rho_old = self.prior_ratio_rho( rho_chain[iter-1] )

        rho_proposed       = self.propose_parameters_rho( rho_chain[iter-1] )
        if (rho_proposed)>0:

            prior_ratio_rho_new = self.prior_ratio_rho( rho_proposed )
            
            likelihood_ratio_rho_new   = self.likelihood_ratio( y,  beta_chain[iter], rho_proposed,      gamma_chain[iter-1], q_chain[:, :, iter-1], eta_zero_chain[:, iter-1], q_r_chain[iter-1] )
            
            ratio_rho = ( prior_ratio_rho_new + likelihood_ratio_rho_new) - ( prior_ratio_rho_old + likelihood_ratio_old)

            if (ratio_rho) >= np.log(np.random.uniform(0,1)):
                rho_chain[iter]         = rho_proposed
                rho_acceptance_ratio = rho_acceptance_ratio + 1

                prior_ratio_rho_old      = prior_ratio_rho_new

                likelihood_ratio_old = likelihood_ratio_rho_new

            else:
                rho_chain[iter] = rho_chain[iter-1]
        else:
            rho_chain[iter] = rho_chain[iter-1]

        # GAMMA
        prior_ratio_gamma_old    = self.prior_ratio_gamma( gamma_chain[iter-1] )
        
        gamma_proposed = self.propose_parameters_gamma( gamma_chain[iter-1] )
        if (gamma_proposed)>0:

            prior_ratio_gamma_new    = self.prior_ratio_gamma( gamma_proposed )     

            likelihood_ratio_gamma_new = self.likelihood_ratio( y, beta_chain[iter], rho_chain[iter], gamma_proposed,      q_chain[:, :, iter-1], eta_zero_chain[:, iter-1], q_r_chain[iter-1] )
            
            ratio_gamma = ( likelihood_ratio_gamma_new + prior_ratio_gamma_new ) - ( likelihood_ratio_old + prior_ratio_gamma_old )

            if ratio_gamma >= np.log(np.random.uniform(0,1)):
                gamma_chain[iter] = gamma_proposed
                gamma_acceptance_ratio = gamma_acceptance_ratio + 1

                prior_ratio_gamma_old      = prior_ratio_gamma_new
                likelihood_ratio_old = likelihood_ratio_gamma_new

            else:
                gamma_chain[iter] = gamma_chain[iter-1]
        else:
            gamma_chain[iter]     = gamma_chain[iter-1]

        # Q_2 IS FIXED FOR EBOLA
        q_proposed       = self.propose_parameters_q( q_chain[:,:, iter-1] )
        prior_ratio_q_new = self.prior_ratio_q( q_proposed )
        prior_ratio_q_old = self.prior_ratio_q( q_chain[:, :, iter-1] )

        if (q_proposed< 1).all() and (q_proposed>=0).all():

            likelihood_ratio_q_new = self.likelihood_ratio( y,  beta_chain[iter], rho_chain[iter], gamma_chain[iter], q_proposed,            eta_zero_chain[:, iter-1], q_r_chain[iter-1] )

            ratio_q = ( likelihood_ratio_q_new + prior_ratio_q_new )- ( likelihood_ratio_old + prior_ratio_q_old )

            if ratio_q >= np.log(np.random.uniform(0,1)):
                q_chain[:, :, iter] = q_proposed
                q_acceptance_ratio = q_acceptance_ratio + 1

                prior_ratio_q_old      = prior_ratio_q_new
                likelihood_ratio_old = likelihood_ratio_q_new

            else:
                q_chain[:, :, iter] = q_chain[:, :, iter-1]
        else:
            q_chain[:, :, iter] = q_chain[:, :, iter-1]

        # ETA_0 is fixed for ebola
        eta_zero_chain[:, iter] = eta_zero_chain[:, iter-1]

        # Q_r
        prior_ratio_q_r_old = self.prior_ratio_q_r( q_r_chain[iter-1] )

        q_r_proposed  = self.propose_parameters_q_r( q_r_chain[iter-1] )
        if q_r_proposed>0:

            prior_ratio_q_r_new = self.prior_ratio_q_r( q_r_proposed )

            likelihood_ratio_q_r_new = self.likelihood_ratio( y,  beta_chain[iter], rho_chain[iter], gamma_chain[iter], q_chain[:, :, iter], eta_zero_chain[:, iter-1], q_r_proposed )

            ratio_q_r = ( likelihood_ratio_q_r_new+ prior_ratio_q_r_new ) - ( likelihood_ratio_old + prior_ratio_q_r_old )

            if (ratio_q_r) >= np.log(np.random.uniform(0,1)):
                q_r_chain[iter] = q_r_proposed
                q_r_acceptance_ratio = q_r_acceptance_ratio + 1

                prior_ratio_q_r_old      = prior_ratio_q_r_new
                likelihood_ratio_old = likelihood_ratio_q_r_new

            else:
                q_r_chain[iter] = q_r_chain[iter-1]
        else:
            q_r_chain[iter] = q_r_chain[iter-1]


        for iter in range(2, sample_size):

            if iter%100==0:

#                 string1 = ["Iteration: "+ str(iter), "\n"]
#                 string2 = ["Acceptance ratio beta: "+ str(beta_acceptance_ratio/iter), "\n"]
#                 string3 = ["Acceptance ratio rho: "+ str(rho_acceptance_ratio/iter), "\n"]
#                 string4 = ["Acceptance ratio gamma: "+ str(gamma_acceptance_ratio/iter), "\n"]
#                 string5 = ["Acceptance ratio q_r: "+ str(q_r_acceptance_ratio/iter), "\n"]
#                 string6 = ["Acceptance ratio q: "+ str(q_acceptance_ratio/iter), "\n"]
#                 string7 = ["Acceptance ratio eta: "+ str(eta_acceptance_ratio/iter), "\n"]

#                 f= open(title,"a")
#                 f.writelines(string1)
#                 f.writelines(string2)
#                 f.writelines(string3)
#                 f.writelines(string4)
#                 f.writelines(string5)
#                 f.writelines(string6)
#                 f.writelines(string7)
#                 f.close()
                print("Iteration: "+ str(iter))
                print("Acceptance ratio beta: "+ str(beta_acceptance_ratio/iter),)
                print("Acceptance ratio rho: "+ str(rho_acceptance_ratio/iter), )
                print("Acceptance ratio gamma: "+ str(gamma_acceptance_ratio/iter),)
                print("Acceptance ratio q_r: "+ str(q_r_acceptance_ratio/iter), )
                print("Acceptance ratio q: "+ str(q_acceptance_ratio/iter), )
                print("Acceptance ratio eta: "+ str(eta_acceptance_ratio/iter),)

            # Beta
            beta_proposed = self.propose_parameters_beta( beta_chain[iter-1] )
            if beta_proposed>0:

                prior_ratio_beta_new    = self.prior_ratio_beta( beta_proposed )

                likelihood_ratio_beta_new = self.likelihood_ratio( y,  beta_proposed, rho_chain[iter-1], gamma_chain[iter-1], q_chain[:, :, iter-1], eta_zero_chain[:, iter-1], q_r_chain[iter-1] )

                ratio_beta = ( likelihood_ratio_beta_new + prior_ratio_beta_new ) - ( likelihood_ratio_old + prior_ratio_beta_old )

                if (ratio_beta) >= np.log(np.random.uniform(0,1)):
                    beta_chain[iter]      = beta_proposed
                    beta_acceptance_ratio = beta_acceptance_ratio + 1

                    prior_ratio_beta_old      = prior_ratio_beta_new

                    likelihood_ratio_old = likelihood_ratio_beta_new

                else:
                    beta_chain[iter] = beta_chain[iter-1]
            else:
                beta_chain[iter]     = beta_chain[iter-1]

            #RHO
            rho_proposed       = self.propose_parameters_rho( rho_chain[iter-1] )
            if (rho_proposed)>0:

                prior_ratio_rho_new = self.prior_ratio_rho( rho_proposed )

                likelihood_ratio_rho_new   = self.likelihood_ratio( y,  beta_chain[iter], rho_proposed, gamma_chain[iter-1], q_chain[:, :, iter-1], eta_zero_chain[:, iter-1], q_r_chain[iter-1] )

                ratio_rho = ( prior_ratio_rho_new + likelihood_ratio_rho_new) - ( prior_ratio_rho_old + likelihood_ratio_old)

                if (ratio_rho) >= np.log(np.random.uniform(0,1)):
                    rho_chain[iter]         = rho_proposed
                    rho_acceptance_ratio = rho_acceptance_ratio + 1

                    prior_ratio_rho_old      = prior_ratio_rho_new

                    likelihood_ratio_old = likelihood_ratio_rho_new

                else:
                    rho_chain[iter] = rho_chain[iter-1]
            else:
                rho_chain[iter] = rho_chain[iter-1]

            # GAMMA
            gamma_proposed       = self.propose_parameters_gamma( gamma_chain[iter-1] )
            if (gamma_proposed)>0:

                prior_ratio_gamma_new    = self.prior_ratio_gamma( gamma_proposed )

                likelihood_ratio_gamma_new = self.likelihood_ratio( y, beta_chain[iter], rho_chain[iter], gamma_proposed, q_chain[:, :, iter-1], eta_zero_chain[:, iter-1], q_r_chain[iter-1] )

                ratio_gamma = ( likelihood_ratio_gamma_new + prior_ratio_gamma_new ) - ( likelihood_ratio_old + prior_ratio_gamma_old )

                if ratio_gamma >= np.log(np.random.uniform(0,1)):
                    gamma_chain[iter] = gamma_proposed
                    gamma_acceptance_ratio = gamma_acceptance_ratio + 1

                    prior_ratio_gamma_old      = prior_ratio_gamma_new
                    likelihood_ratio_old = likelihood_ratio_gamma_new

                else:
                    gamma_chain[iter] = gamma_chain[iter-1]
            else:
                gamma_chain[iter]     = gamma_chain[iter-1]

            # Q_2 IS FIXED FOR EBOLA
            q_proposed       = self.propose_parameters_q( q_chain[:,:, iter-1] )
            if (q_proposed< 1).all() and (q_proposed>=0).all():

                prior_ratio_q_new = self.prior_ratio_q( q_proposed )

                likelihood_ratio_q_new = self.likelihood_ratio( y,  beta_chain[iter], rho_chain[iter], gamma_chain[iter], q_proposed, eta_zero_chain[:, iter-1], q_r_chain[iter-1] )

                ratio_q = ( likelihood_ratio_q_new + prior_ratio_q_new ) - ( likelihood_ratio_old + prior_ratio_q_old )

                if ratio_q >= np.log(np.random.uniform(0,1)):
                    q_chain[:, :, iter] = q_proposed
                    q_acceptance_ratio = q_acceptance_ratio + 1

                    prior_ratio_q_old      = prior_ratio_q_new
                    likelihood_ratio_old = likelihood_ratio_q_new

                else:
                    q_chain[:, :, iter] = q_chain[:, :, iter-1]
            else:
                q_chain[:, :, iter] = q_chain[:, :, iter-1]

            # ETA_0 IS FIXED FOR EBOLA
            eta_zero_chain[:, iter] = eta_zero_chain[:, iter-1]


            # Q_r
            q_r_proposed  = self.propose_parameters_q_r( q_r_chain[iter-1] )
            if q_r_proposed>0:

                prior_ratio_q_r_new = self.prior_ratio_q_r( q_r_proposed )

                likelihood_ratio_q_r_new = self.likelihood_ratio( y,  beta_chain[iter], rho_chain[iter], gamma_chain[iter], q_chain[:, :, iter], eta_zero_chain[:, iter-1], q_r_proposed )

                ratio_q_r = ( likelihood_ratio_q_r_new + prior_ratio_q_r_new ) - ( likelihood_ratio_old + prior_ratio_q_r_old )

                if (ratio_q_r) >= np.log(np.random.uniform(0,1)):
                    q_r_chain[iter] = q_r_proposed
                    q_r_acceptance_ratio = q_r_acceptance_ratio + 1

                    prior_ratio_q_r_old      = prior_ratio_q_r_new
                    likelihood_ratio_old = likelihood_ratio_q_r_new

                else:
                    q_r_chain[iter] = q_r_chain[iter-1]
            else:
                q_r_chain[iter] = q_r_chain[iter-1]


        return beta_chain, rho_chain, gamma_chain, eta_zero_chain, q_chain, x_chain, q_r_chain




