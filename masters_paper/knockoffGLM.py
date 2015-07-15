""" This module is meant to implement and test the knockoff procedure on L1 regularized GLM models. The method was originally designed for linear least squares regression, and there is no theoretical guarantee it will work on other GLMS. Hopefully, with simulation we can discover if it breaks, if so where it breaks, and maybe how to fix it in those situations.

More can be read about the original procedure on least squares linear models here: http://web.stanford.edu/~candes/Knockoffs/index.html

A lot of the code is stolen/modified from the knockoff R package

The method for binary models is my development, based on the paper "On parametric families for sampling binary data with specified mean and correlation", 2012, Christian Schafer
"""

# imports 
from __future__ import division
import numpy as np
import numpy.linalg as nplin
import numpy.random as npran
import scipy as sp
import scipy.linalg as splin
import scipy.stats as spstat
import matplotlib.pyplot as plt
import cvxpy as cvx
import rpy2 as rp 
import rpy2.robjects.numpy2ri
import statsmodels.api as sm 
from rpy2.robjects.packages import importr
from glmnet import LogisticNet, ElasticNet
from sklearn.preprocessing import normalize
from scipy.special import expit as invlogit
from scipy.special import logit
from sys import float_info

def dinvlogit(x):
    '''Derivative of logit function at each point of vextor x'''
    return invlogit(x)*(1-invlogit(x))

def get_index_z(coef_vector):
    """This returns the first index for a variable from the glmnet coef matrix which isn't zero"""
    n     = coef_vector.size
    index = 0 
    while index<n: 
        if abs(coef_vector[index])>1e-6: return index
        index+=1
    return n+1

def solve_sdp(G):
    """ Solves the SDP problem:
    
    maximize    1' * s
    subject to  0 <= s <= 1
                [G, G - diag(s); G - diag(s), G] >= 0
    stolen directly from knockoff R package - solve_sdp.py
    """
    assert G.ndim == 2 and G.shape[0] == G.shape[1]
    p = G.shape[0]
    
    # Define the problem.
    s = cvx.Variable(p)
    objective = cvx.Maximize(cvx.sum_entries(s))
    constraints = [
        0 <= s, s <= 1,
        2*G - cvx.diag(s) == cvx.Semidef(p),
    ]
    prob = cvx.Problem(objective, constraints)
    
    # Solve the problem.
    prob.solve()
    assert prob.status == cvx.OPTIMAL
    
    # Return as array, not as matrix.
    return np.asarray(s.value).flatten()

class knockoff_net(object):
    """ Parent class for knockoff lasso and knockoff logistic regression. Defines everything besides fitting the particular GlmNet object """
    def __init__(self,y,X,q,knockoff='binary',cov_method='equicor',randomize=False,MCsize=10000,tol=1E-5,maxiter=40,full_A=False):
        self.y      = y
        self.X      = normalize(X.astype(float),norm='l2',axis=0)
        self.X_orig = X
        self.q      = q                # Level we want to control FDR at
        self.n,self.p = self.X.shape
        self.knockoff_type = knockoff  # what type of knockoff are we generating - binary or the original deterministic ones?
        self.cov_type = cov_method     # how are we generating the desired covariance matrix - equicor or SDP?
        self.randomize=False
        self.tol = tol
        self.maxiter = maxiter
        self.MCsize = MCsize   # how many values for Monte Carlo simulation
        self.zerotol = 1E-5
        self.full_A = full_A

    def _create_SDP(self):
        """ Creates the SDP knockoff of X"""
 
        # Check for rank deficiency (will add later).
 
        # SVD and come up with perpendicular matrix
        U, d, V = nplin.svd(self.X,full_matrices=True) 
        d[d<0] = 0
        U_perp = U[:,self.p:(2*self.p)]
        if self.randomize:
            U_perp = np.dot(U_perp,splin.orth(npran.randn(self.p,self.p)))
 
        # Compute the Gram matrix and its (pseudo)inverse.
        G     = np.dot(V.T * d**2 ,V)
        G_inv = np.dot(V.T * d**-2,V)
 
        # Optimize the parameter s of Equation 1.3 using SDP.
        s = solve_sdp(G)
        s[s <= self.zerotol] = 0
 
        # Construct the knockoff according to Equation 1.4:
        C_U,C_d,C_V = nplin.svd(2*np.diag(s) - (s * G_inv.T).T * s)
        C_d[C_d < 0] = 0
        self.X_ko = self.X - np.dot(self.X,G_inv*s) + np.dot(U_perp*np.sqrt(C_d),C_V)

    def _create_equicor(self):
        """ Creates the equal correlation knockoff of X"""
        # Check for rank deficiency (will add later).
 
        # SVD and come up with perpendicular matrix
        U, d, V = nplin.svd(self.X,full_matrices=True) 
        d[d<0] = 0
        U_perp = U[:,self.p:(2*self.p)]
        U = U[:,:self.p]
        if self.randomize:
            U_perp = np.dot(U_perp,splin.orth(npran.randn(self.p,self.p)))
 
        # Set s = min(2 * smallest eigenvalue of X'X, 1), so that all the correlations
        # have the same value 1-s.
        lambda_min = min(d)**2
        s = min(2*lambda_min,1)
 
        # Construct the knockoff according to Equation 1.4.
        s_diff = 2*s - (s/d)**2
        s_diff[s_diff<0]=0 # can be negative due to numerical error
        self.X_ko = np.dot(U*(d-s/d) + U_perp*(np.sqrt(s_diff)) , V)

    def _original_knockoff(self):
        """ Creates the original style knockoffs"""
        if self.cov_type == 'equicor':
            self._create_equicor()
        else: 
            self._create_SDP()

    def _binary_knockoff(self):
        ''' This creates the new binary knockoffs, which are random multivariate bernoulli which should have, in expectation,
        the same first two moments as X. Only will work if X is all binaray '''

        ###############################################
        # Figure out desired cross moment matrix
        ###############################################

        # First, calculate s
        # Largely replicates begining of self._create_equicor()/self._create_SDP()

        # SVD and come up with perpendicular matrix
        U, d, V = nplin.svd(self.X,full_matrices=True) 
        d[d<0] = 0
        if   self.cov_type == 'SDP': 
            # Compute the Gram matrix and its (pseudo)inverse.
            G     = np.dot(V.T * d**2 ,V)
            G_inv = np.dot(V.T * d**-2,V)
         
            # Optimize the parameter s of Equation 1.3 using SDP.
            s = solve_sdp(G)
            s[s <= self.zerotol] = 0

        # Same as first part self._create_equicor 
        elif self.cov_type == 'equicor':
            # Set s = min(2 * smallest eigenvalue of X'X, 1), so that all the correlations
            # have the same value 1-s.
            lambda_min = min(d)**2
            s = np.ones(self.p)*min(2*lambda_min,1)

        # Return to the original X
        self.X = self.X_orig              # Now that we've adjusted the variance, revert back to the original X 

        # Calculate first moments 
        mu       = np.mean(self.X,axis=0)
        mu_lrg   = np.concatenate((mu,mu))

        # Scale s up to original variance

        # Generate the covariance matrix Sigma
        M    = np.dot(self.X.T,self.X)/self.n
        Sigma = M - np.outer(mu,mu)
        s     = s*np.diag(Sigma)        # scale s back up to original variance


        # Generate the large Sigma matrix, and keep shrinking s until its a valid crossmoment matrix
        test = False 
        shrink = 1 
        while not test:
            # s will be shrunk a little more each time the test is failed
            s = s*shrink

            # Generate large Sigma and crossmoment matrix
            Sigma_lrg = np.hstack(((np.vstack((Sigma,Sigma-np.diag(s)))),np.vstack((Sigma-np.diag(s),Sigma))))
            self.M   = Sigma_lrg + np.outer(mu_lrg,mu_lrg)

            # Check condition 2.1 from Schafer - make sure its a valid crossmoment matrix
            cond1 = self.M <= np.minimum(np.outer(mu_lrg,np.ones(2*self.p)),np.outer(np.ones(2*self.p),mu_lrg))
            cond2 = np.outer(mu_lrg,np.ones(2*self.p))+np.outer(np.ones(2*self.p),mu_lrg) -1 <= self.M
            testmat = np.logical_or(np.logical_and(cond1,cond2),np.diag(np.ones(2*self.p)))  # ignore diagonal
            test = np.min(testmat)

            shrink_s = False
            if shrink_s:
                if shrink<.25:
                    test = True
                
                if test and shrink<1:
                    print "s had to be shrunk by a factor of %.2f" % shrink

                shrink = shrink*.90
            elif not test:
                test = True 
                print "Cross moment test failed, didn't do anything about it"

        ####################################################
        # Fit the A matrix for the original variables first. 
        ####################################################

        A = np.zeros((2*self.p,2*self.p))

        # we don't actually need to fit A on the original data; we use the original points to simulate
        if self.full_A == True:
            # We can fit on the actual data, making this easier
            A[1,1] = logit(mu[0])
            for i in np.arange(0,self.p):
            # inject the paramters from the logit X_i ~ X_1 + ... + X_(i-1) + 1 into the ith row of A
                A[i,0:(i+1)] = sm.GLM(self.X[:,i],np.hstack((self.X[:,0:i],np.ones((self.n,1)))),family=sm.families.Binomial()).fit().params

        ###################################################
        # Derive remaing A from Newton-Raphson
        ###################################################

        # Largely from 5.1 in Schafer, including notation

        # the current value and the derivatves are derived by simulation

        # Rather than simulate entirely new X at each stage, I will used a fixed set of X_1 ... X_(i-1)
        # This definitely makes sense for the orignal X vars (why simulate when we already have it), but possibly less sense for the knockoffs
        # To get the desired size of montecarlo simulation, I will replicate X until it has at least self.MCsize rows
        repl = np.min((self.MCsize//self.n,1))
        X_fix = np.repeat(self.X,repl,0)
        nMC   = X_fix.shape[0]
        ones  = np.ones((nMC,1))

        # sequence of portions between 0 and 1 to deal with case of ill conditioned hessian
        por_seq = np.append(0,invlogit(np.arange(-5,5,.5)))
        por_seq = np.append(por_seq,1)
        self.por = np.empty((1,0))

        for i in np.arange(self.p,2*self.p):
            # m are the cross moments we are trying to fit
            m = self.M[i,0:(i+1)]
            # X_tmp is the list of (x_k,1) vectors
            X_tmp     = np.hstack((X_fix,ones))
            # X_out_tmp is the array of (x_k,1)'(x_k,1) matricies (emperical cross moments for the vector)
            # I shouldn't have to do the full multiplication each time (only the X_fix*ones part is new). improve later
            X_out_tmp = X_tmp[:,:,np.newaxis]*X_tmp[:,np.newaxis,:]

            # Now, the Newton-Raphson steps
            # If the hessian becomes singular, we will relax the cross moment requirements, as described in Schafer 5.1 point 2
            # the idea is that the problem is relaxed until X_i is independent of all prior vars
            # as por increases, covariance drops
            for por in por_seq:
                # m are the cross moments we are trying to fit
                m = (1-por)*self.M[i,0:(i+1)] + por*self.M[i,i]*np.append(np.diag(self.M)[0:i],1)

                # a is the row we are adding to A. Initialize with values as if independent all other vars
                a = np.append(np.zeros(i),logit(mu_lrg[i]))
                error, counter = np.Inf,0

                # If we fit the given m without running into ill-conditioned hessians, we can stop
                illcond = False
                while error > self.tol and self.maxiter>counter and not illcond:
                    # mu is the probability of X_i=1|X_1 ... X_(i-1) times the normalizing constant
                    mu = invlogit(np.dot(X_tmp,a))[:,np.newaxis]
                    # fa is the current estimate of the cross moments
                    fa = np.mean(mu*X_tmp      ,axis=0)
                    # mup is the derivative of mu with respect to np.dot(X_tmp,a)
                    mup = dinvlogit(np.dot(X_tmp,a))[:,np.newaxis,np.newaxis]
                    # fa is the Hessian of fa with respect to a
                    fap = np.mean(mup*X_out_tmp,axis=0)
                    if nplin.cond(fap) < 1/float_info.epsilon:
                        a = a - np.dot(nplin.inv(fap),fa-m)
                        counter += 1
                        error = nplin.norm(fa-m)  # how close am I?
                    else:
                        illcond = True

                # Stop once have gone through Newton-Raphson without running into ill conditioning 
                if not illcond:
                    self.por = np.append(self.por,por)
                    break

            # put a into A matrix, draw X_i for 'fixed' matrix, update X_out_fix
            A[i,0:(i+1)] = a
            X_fix = np.hstack((X_fix,npran.binomial(1,invlogit(np.dot(X_tmp,a)))[:,np.newaxis]))

        ##############################################
        # Wrapup and get X_ko
        ##############################################

        # since we've been drawing the X along the way, can subset X_fix to get X_ko
        self.X_ko = X_fix[0::repl,self.p:]

        # hang onto A
        self.A = A

    def _get_z(self): 
        """ Given the coefficient matrix from a glmnet object, returns the Z, the first non-zero entries"""
        self.z_rank  = np.array(map(get_index_z,self.coef_matrix))
        self.z_value = np.array(map(lambda x: self.lambdas[x] if x<self.lambdas.size else 0,self.z_rank))

    def _get_w(self): 
        """Produces the w values using the z values"""
        self.w_filter = np.sign(self.z_rank[self.p:(2*self.p)]-self.z_rank[0:self.p])
        self.w_rank   = self.w_filter * np.nanmin((self.z_rank[self.p:(2*self.p)],self.z_rank[:self.p]),axis=0)
        self.w_value  = self.w_filter * np.max((self.z_value[self.p:(2*self.p)],self.z_value[:self.p]),axis=0)

    def _get_T(self):
        """ Calculates the data-dependent threshold for the statistics """
        Ws = set(np.abs(self.w_value))
        self.T = np.amin(np.append([t for t in Ws if np.sum(self.w_value<=-t)/max(1,np.sum(self.w_value>=t)) <= self.q],np.inf))

    def _get_S(self):
        """ Calculate the vector selecting features S """
        self.S = self.w_value>=self.T

class knockoff_logit(knockoff_net):
    """ Preforms the knockoff technique with logistic regression """

    def fit(self):
        """ Generates the knockoffs, fits the regression, and performs the FDR calculations """
        # Generate knockoff as inherited from knockoff_net
        if   self.knockoff_type == 'original': self._original_knockoff()
        elif self.knockoff_type == 'binary':   self._binary_knockoff()

        X_lrg = np.concatenate((self.X,self.X_ko), axis=1)

        # initialize and fit the glmnet object
        self.lognet = LogisticNet(alpha=1) 
        self.lognet.fit(X_lrg,self.y,normalize=False,include_intercept=False)

        # pull out some values from the glmnet object and clean
        self.lambdas = self.lognet.out_lambdas
        self.var_index_ent = np.sort(self.lognet._indices)
        self.coef_matrix = np.zeros((2*self.p,self.lognet.n_lambdas))
        self.coef_matrix[self.var_index_ent] = self.lognet._comp_coef.squeeze()[self.lognet._indices]

        # figure out when different variables entered the model
        self.var_entered = np.zeros(2*self.p).astype(bool)
        self.var_entered[self.var_index_ent] = True

        # Preform all the FDR calculations as inherited from knockoff_net
        self._get_z()
        self._get_w()
        self._get_T()
        self._get_S()


class knockoff_lasso(knockoff_net):
    """ Preforms the knockoff technique with lasso """
    
    def fit(self):
        """ Generates the knockoffs, fits the regression, and performs the FDR calculations """
        # Generate knockoff as inherited from knockoff_net
        if   self.knockoff_type == 'original': self._original_knockoff()
        elif self.knockoff_type == 'binary':   self._binary_knockoff()

        X_lrg = np.concatenate((self.X,self.X_ko), axis=1)

        # initialize the glmnet object
        self.elasticnet = ElasticNet(alpha=1) 
        self.elasticnet.fit(X_lrg,self.y,normalize=False,include_intercept=False)

        # pull out some values from the glmnet object and clean
        self.lambdas = self.elasticnet.out_lambdas
        self.var_index_ent = np.sort(self.elasticnet._indices)
        self.coef_matrix = np.zeros((2*self.p,self.elasticnet.n_lambdas))
        self.coef_matrix[self.var_index_ent] = self.elasticnet._comp_coef.squeeze()[self.elasticnet._indices]

        # figure out when different variables entered the model
        self.var_entered = np.zeros(2*self.p).astype(bool)
        self.var_entered[self.var_index_ent] = True

        # Preform all the FDR calculations as inherited from knockoff_net
        self._get_z()
        self._get_w()
        self._get_T()
        self._get_S()

def main():
    print 'nobody lives here right now'

if __name__ == '___main__':
    status = main()
    sys.exit(status)

