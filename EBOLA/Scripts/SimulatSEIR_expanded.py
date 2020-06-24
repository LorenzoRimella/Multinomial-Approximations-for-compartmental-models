import numpy as np

class SEIR_expanded:

    def __init__( self, N, beta, rho, gamma, q_r, t_star ):

        self.N       = N
        self.beta    = beta
        self.rho     = rho
        self.gamma   = gamma
        self.q_r     = q_r
        self.t_star  = t_star

    def Evolution( self, eta_zero, T, q = (np.zeros((4,4))+1) ):

        if isinstance(eta_zero, np.ndarray)==False:
            raise ValueError("The initial condition should be a dictionary wit keys: S,E,I,R.")

        if sum(eta_zero) != 1:
            raise ValueError("The total population is not consistent with the total population number.")

        SEIR = np.zeros((4, T))
        SEIR[:, 0] = eta_zero*self.N

        y    = np.zeros((4, 4, T))
        y[:, :, 0]    = np.zeros((4, 4))

        for t in range(1, T):

            beta_restr = self.beta*(t< self.t_star) + self.beta*(np.exp(-self.q_r*(t-self.t_star)*(t>= self.t_star)))*(t>= self.t_star)

            Pt = 1 - np.exp(-(beta_restr/self.N)*SEIR[2,t-1] )
            pC = 1 - np.exp(-self.rho)
            pR = 1 - np.exp(-self.gamma)

            Bt = np.random.binomial( SEIR[0, t-1], Pt)
            Ct = np.random.binomial( SEIR[1, t-1], pC)
            Dt = np.random.binomial( SEIR[2, t-1], pR)

            SEIR[0, t] = (SEIR[0, t-1] - Bt)
            SEIR[1, t] = (SEIR[1, t-1] + Bt - Ct)
            SEIR[2, t] = (SEIR[2, t-1] + Ct - Dt)
            SEIR[3, t] = (SEIR[3, t-1] + Dt)

            y[0, 1, t] = np.random.binomial( Bt, q[0,1] )
            y[1, 2, t] = np.random.binomial( Ct, q[1,2] )
            y[2, 3, t] = np.random.binomial( Dt, q[2,3] )

        return SEIR, y