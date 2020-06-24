import numpy as np
from Scripts.HMMSEIR_expanded import HMM_approxSEIR_expanded

class EM_approxSEIR_expanded:
    
    def __init__( self, N, beta, rho, gamma, q, init_pop, q_r, t_star, y, h_fac, beta_prior, rho_prior, gamma_prior):
    
        self.N            = N
        self.beta_init    = beta
        self.rho_init     = rho
        self.gamma_init   = gamma
        self.q_init       = q
        self.init_pop     = init_pop
        self.q_r     = q_r
        self.t_star  = t_star
        self.y            = y
        self.h_fac        = h_fac
        self.beta_prior = beta_prior
        self.rho_prior = rho_prior
        self.gamma_prior = gamma_prior

    def iterate( self, n_iter ):
        
        beta_est = np.zeros(n_iter+1)
        beta_est[0] = self.beta_init
        rho_est = np.zeros(n_iter+1)
        rho_est[0] = self.rho_init
        gamma_est = np.zeros(n_iter+1)
        gamma_est[0] = self.gamma_init
        q_est = np.zeros([4, 4, n_iter+1])
        q_est[:,:,0] = self.q_init
        init_pop_est = np.zeros([4,n_iter+1])
        init_pop_est[:,0] = self.init_pop

       
        HMM = HMM_approxSEIR_expanded( self.N, self.beta_init, self.rho_init, self.gamma_init, self.q_init, \
                                      self.init_pop,  self.q_r, self.t_star)
        
        T = np.size(self.y[0, 0, :])
        h = 1/T * 1/self.N  * self.h_fac#step size for beta gradient ascent
        #beta_grid = np.linspace(0.01,3,num=20)
        #beta_obj = np.zeros(20)
        
        for i in range(n_iter):
    
            #forward-backward
            pitt, Kappa, pitt_expanded, _ = HMM.filtering(self.y)
            #pist, Kappa_hat, pist_expanded = HMM.smoothing(pitt, Kappa)
            pist, pist_expanded = HMM.smoothing(pitt_expanded, pitt)
            
            pist = pist / np.outer(np.ones(4),pist.sum(axis=0))
            
            #print(pist_expanded)
        
            #parameter updates
            
            #HMM.rho = np.log( 1 + np.sum( pist[2,1:] * Kappa_hat[2,1,:]) 
            #                        / (self.rho_prior[1]/self.N + np.sum( pist[1,1:] * Kappa_hat[1,1,:])))
            HMM.rho = np.log( 1 + np.sum( pist_expanded[1,2,:]) 
                                    / (self.rho_prior[1]/self.N + np.sum( pist_expanded[1,1,:])))
            
            #HMM.gamma = np.log( 1 + np.sum( pist[3,1:] * Kappa_hat[3,2,:]) 
            #                        / (self.gamma_prior[1]/self.N + np.sum( pist[2,1:] * Kappa_hat[2,2,:]))) 
            HMM.gamma = np.log( 1 + np.sum( pist_expanded[2,3,:]) 
                                    / (self.gamma_prior[1]/self.N + np.sum( pist_expanded[2,2,:]))) 
            
            
            #need to double check the following and the calculation of pist_expanded
            HMM.q[1,2] = np.minimum(1/self.N * np.sum(self.y[1,2,1:]) / np.sum(pist_expanded[1,2,:]),1)
            HMM.q[2,3] = np.minimum(1/self.N * np.sum(self.y[2,3,1:]) / np.sum(pist_expanded[2,3,:]),1)
            
            #HMM.init_pop = pist[:,0] #NB does not necessarily give integers
            
                    
                            
            #GRADIENT CALCULATION NEEDS TO BE CHECKED FOR EBOLA MODEL
            if self.h_fac > 0 :   #stochastic gradient update for beta                                          
                beta_grad = (self.beta_prior[0]-1)/HMM.beta - self.beta_prior[1] #assuming that the alpha parameter in the prior for beta is 1
                for t in range(1,T-1):
                    Z_prob_mat =  np.outer(np.ones(4),pist[:,t]) * Kappa_hat[:,:,t-1].transpose()
                    Z_draw = np.random.multinomial(self.N,Z_prob_mat.reshape(16),size=1).reshape([4,4])
                    X_draw = np.sum(Z_draw,axis=1)
                    decay_factor = 1*(t< self.t_star) + (np.exp(-self.q_r*(t-self.t_star)))*(t>= self.t_star)
                    if X_draw[2] > 0:
                        beta_grad = beta_grad - (self.N - 1) * Z_prob_mat[0,0] * pist[2,t-1] * decay_factor + \
                        1/self.N * Z_draw[0,1] * X_draw[2] * decay_factor / (np.exp(HMM.beta * X_draw[2]/self.N * decay_factor)-1) 
                    else:
                        beta_grad = beta_grad - (self.N - 1) * Z_prob_mat[0,0] * pist[2,t-1] * decay_factor
                print('beta_grad =',  beta_grad)
                HMM.beta = np.maximum(HMM.beta + h * beta_grad, 0.01) #project on to allowed values for beta
                                          
                                          
            #store the new params in param_est
            beta_est[i+1] = HMM.beta
            rho_est[i+1] = HMM.rho
            gamma_est[i+1] = HMM.gamma
            q_est[:,:,i+1] = HMM.q
            init_pop_est[:,i+1] = HMM.eta_zero
            
        
        #and evaluate the marginal likelihood for the final parameter values
        _, _, _ , pitt_prev = HMM.filtering(self.y)
        T  = np.size( self.y, 2 )
        C1 = np.size( self.y, 0 )
        C2 = np.size( self.y, 1 )
        
        #from Lorenzo's MCMC code
        y_from1 = self.y[:, : , 1:]
        y_from1_index_nozeros = (y_from1!=0)
        y_pi_from1_index_nozeros = (pitt_prev[:, :, 1:][y_from1_index_nozeros]!=0)
        m_groups_prob1   = np.sum( np.sum( np.sum( y_from1[y_from1_index_nozeros][y_pi_from1_index_nozeros]*np.log( ( np.repeat( HMM.q, T-1).reshape( C1, C2, T-1)*pitt_prev[:, :, 1:])[y_from1_index_nozeros][y_pi_from1_index_nozeros] ), 0 ), 0 ), 0)
        mp1_groups_prob1 = np.sum( ( self.N - np.sum(np.sum(self.y[:, :, 1:], 0), 0) )*np.log( (1 - np.sum( ( np.repeat(HMM.q, T-1).reshape( C1, C2, T-1)*pitt_prev[:, :, 1:] ), 0 ) ) ) )
        log_like = m_groups_prob1+mp1_groups_prob1
        
        return beta_est, rho_est, gamma_est, q_est, init_pop_est, log_like