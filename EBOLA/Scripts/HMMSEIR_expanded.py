import numpy as np
from scipy.stats import gamma as RVgamma
# the gamma distribution consider a varying shape parameter and a scale parameter equal to 1

class HMM_approxSEIR_expanded:

    def __init__( self, N, beta, rho, gamma, q, eta_zero, q_r, t_star ):

        self.N        = N
        self.beta     = beta
        self.rho      = rho
        self.gamma    = gamma
        self.q        = q
        self.eta_zero = eta_zero
        self.q_r     = q_r
        self.t_star  = t_star



    def eta_computation(self, T):

        eta       = np.zeros((4, T))
        eta[:, 0]      = self.eta_zero

        pC = 1 - np.exp(-self.rho)
        pR = 1 - np.exp(-self.gamma)

        for t in range(1, T):
            Kappa_eta_prev = np.array([[ np.exp(-self.beta*eta[2,t-1]), 1 - np.exp(-self.beta*eta[2,t-1]), 0, 0 ], [ 0, 1 - pC, pC, 0 ], [ 0, 0, 1 - pR, pR ], [ 0, 0, 0, 1 ]])
            eta[:, t] = eta[:, t-1] @ Kappa_eta_prev

        return eta


    def filtering(self, y):

        T = np.size(y[0, 0, :])

        pC = 1 - np.exp(-self.rho)
        pR = 1 - np.exp(-self.gamma)

        pitt = np.zeros([4,T])
        pitt[:,0]= self.eta_zero

        pitt_expanded  = np.zeros((4, 4, T))
        pitt_expanded[0, :, 0] = pitt[:,0] 

        pitt_prev_expanded  = np.zeros((4, 4, T))
        pitt_prev_expanded[0, :, 0] = pitt[:,0]
        
        Kappa = np.zeros([4,4,T-1])
        
        pitt_expanded_q = np.zeros([4,4,T])

        for t in range(1, T):

            beta_restr = self.beta*(t< self.t_star) + self.beta*(np.exp(-self.q_r*(t-self.t_star)))*(t>= self.t_star)

            Kappa_eta_prev = np.array([[ np.exp(-beta_restr*pitt[2,t-1]), 1 - np.exp(-beta_restr*pitt[2,t-1]), 0, 0 ], [ 0, 1 - pC, pC, 0 ], [ 0, 0, 1 - pR, pR ], [ 0, 0, 0, 1 ]])

            Kappa[:,:,t-1] = Kappa_eta_prev
            
            pitt_prev_expanded[:,:, t]   = Kappa_eta_prev*( np.sum(pitt_expanded[:, :, t-1], 0) ).reshape(4,1)
            #rho_vec              = pitt_prev_expanded[:,:, t]*(1-self.q)
            #rho_vec              = rho_vec/np.sum(rho_vec)
            
            pitt_expanded_q[:,:,t] = pitt_prev_expanded[:,:, t]*(1-self.q)
            pitt_expanded_q[:,:,t] = pitt_expanded_q[:,:,t]/np.sum(pitt_expanded_q[:,:,t])

            pitt_expanded[:,:, t]       = y[:,:, t]/self.N + ( 1 - (np.sum( y[:,:, t] ))/(self.N) )*pitt_expanded_q[:,:,t]
            pitt[:,t] = np.sum( pitt_expanded[:,:, t], 0 )

        return pitt, Kappa, pitt_expanded, pitt_prev_expanded
    
    
    def smoothing(self,pitt_expanded, pitt):
    
        T = np.size(pitt_expanded[1,1,:])
        
        pist = np.zeros((4, T))
        pist[:,T-1] = np.sum(pitt_expanded[:,:,T-1],0)
        L = np.zeros((4,4))
        pist_expanded  = np.zeros((4, 4, T))
        pist_expanded[:,:,T-1] = pitt_expanded[:,:,T-1]

        for t in range(T-1,1,-1):

            pist[:,t-1] = np.sum(pist_expanded[:,:,t],1)
            L[np.outer(pitt[:,t-1],np.ones(4))!=0] = np.transpose(pitt_expanded[:,:,t-1])[np.outer(pitt[:,t-1],np.ones(4))!=0] / np.outer(pitt[:,t-1],np.ones(4))[np.outer(pitt[:,t-1],np.ones(4))!=0] 

            pist_expanded[:,:,t-1] = np.outer(np.ones(4),pist[:,t-1]) * np.transpose(L) 
         
        pist[:,0] = np.sum(pist_expanded[:,:,1],1)
        pist_expanded[0, :, 0] = pist[:,0]
            
        return pist, pist_expanded