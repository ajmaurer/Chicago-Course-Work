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
#import cvxpy as cvx 
import statsmodels.api as sm
from glmnet import LogisticNet, ElasticNet
from sklearn.preprocessing import normalize
from scipy.special import expit as invlogit
from scipy.special import logit
from sys import float_info
from scipy.optimize import root 

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

def draw_random_binary(n,A):
    ''' If p is the size of the square, lower triangular A, creates a n by p matrix of random binary vectors using Schafer's method '''
    m,p = A.shape
    ones   = np.ones((n,1))
    output = np.empty((n,p))
    for i in np.arange(0,p):
        output[:,i] = npran.binomial(1,invlogit(np.dot(np.hstack((output[:,0:i],ones)),A[i,0:(i+1)])))
    return output

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
    def __init__(self,y,X,q,knockoff='binary',cov_method='equicor',randomize=False,MCsize=10000,fresh_sim=True,pseudocount=0,intercept=False):
        self.y      = y
        self.X      = normalize(X.astype(float),norm='l2',axis=0)
        self.X_orig = X
        self.q      = q                # Level we want to control FDR at
        self.n,self.p = self.X.shape
        self.knockoff_type = knockoff  # what type of knockoff are we generating - binary or the original deterministic ones?
        self.cov_type = cov_method     # how are we generating the desired covariance matrix - equicor or SDP?
        self.randomize=False
        self.MCsize = MCsize   # how many values for Monte Carlo simulation
        self.zerotol = 1E-5
        self.fresh_sim = fresh_sim 
        self.pseudocount=pseudocount
        self.intercept = intercept 

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
        self.s = solve_sdp(G)
        self.s[s <= self.zerotol] = 0
 
        # Construct the knockoff according to Equation 1.4:
        C_U,C_d,C_V = nplin.svd(2*np.diag(s) - (self.s * G_inv.T).T * self.s)
        C_d[C_d < 0] = 0
        X_ko = self.X - np.dot(self.X,G_inv*s) + np.dot(U_perp*np.sqrt(C_d),C_V)
        self.X_lrg = np.concatenate((self.X_orig,X_ko), axis=1)

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
        self.s = min(2*lambda_min,1)
 
        # Construct the knockoff according to Equation 1.4.
        s_diff = 2*self.s - (self.s/d)**2
        s_diff[s_diff<0]=0 # can be negative due to numerical error
        X_ko = np.dot(U*(d-self.s/d) + U_perp*(np.sqrt(s_diff)) , V)
        self.X_lrg = np.concatenate((self.X_orig,X_ko), axis=1)

    def _original_knockoff(self):
        """ Creates the original style knockoffs"""
        if self.cov_type == 'equicor':
            self._create_equicor()
        else: 
            self._create_SDP()
        self.emp_ko_corr= np.dot(self.X_lrg.T,self.X_lrg)[:self.p,self.p:2*self.p][np.identity(self.p)==1]/self.n

    def _derive_crossmoments(self):
        """Figures out desired cross momemnt matrix matrix for binary knockoffs"""

        ###############################################
        # Figure out desired cross moment matrix
        ###############################################

        # First, calculate s
        # Largely replicates begining of self._create_equicor()/self._create_SDP()
        # However, we now want to perform calculations on the covariance matrix
        # The goal is minimize the correlation between X_i and ~X_i
        # Generate the second moment matrix Sigma
        mu   = np.mean(self.X_orig,axis=0)
        relaxer = np.where(np.identity(self.p),np.diag(mu),np.outer(mu,mu))
        M    = (np.dot(self.X_orig.T,self.X_orig) + self.pseudocount*relaxer)/(self.n + self.pseudocount)
        Cov  = M - np.outer(mu,mu)
        #Cov  = np.cov(self.X_orig,rowvar=0,bias=1)      # I'm doing redundant calculations here
        Corr = Cov/np.sqrt(np.outer(mu*(1-mu),mu*(1-mu)))
        #Corr = np.corrcoef(self.X_orig,rowvar=0,bias=1)

        # SVD and come up with perpendicular matrix
        d, V = nplin.eig(Corr)
        #U, d, V = nplin.svd(self.X,full_matrices=True) 
        d[d<0] = 0
        if   self.cov_type == 'SDP': 
         
            # Optimize the parameter s of Equation 1.3 using SDP.
            self.s = solve_sdp(Corr)
            self.s[s <= self.zerotol] = 0

        # Same as first part self._create_equicor 
        elif self.cov_type == 'equicor':
            # Set s = min(2 * smallest eigenvalue of X'X, 1), so that all the correlations
            # have the same value 1-s.
            lambda_min = min(d)**2
            self.s = np.ones(self.p)*min(2*lambda_min,1)

        # Scale s up to original variance
        self.s     = self.s*np.diag(Cov)        # scale s back up to original variance

        # Calculate first moments 
        self.mu_lrg   = np.concatenate((mu,mu))

        # Generate the large Sigma matrix, and keep shrinking s until its a valid crossmoment matrix
        test1, test2 = True, True
        shrink = 1
        while test1 or test2:
            # s will be shrunk a little more each time the test is failed
            self.s = self.s*shrink

            # Generate large Sigma and crossmoment matrix
            Cov_lrg = np.hstack(((np.vstack((Cov,Cov-np.diag(self.s)))),np.vstack((Cov-np.diag(self.s),Cov))))
            self.M   = Cov_lrg + np.outer(self.mu_lrg,self.mu_lrg)

            # Check condition 2.1 from Schafer - make sure its a valid crossmoment matrix
            test1, test2 = False, False
            cond1 = self.M >= np.minimum(np.outer(self.mu_lrg,np.ones(2*self.p)),np.outer(np.ones(2*self.p),self.mu_lrg))
            cond1 = np.logical_and(np.diag(1-np.ones(2*self.p)),cond1)
            if (np.sum(cond1)) > 0:
                print "%d cross moments which were too high" % (np.sum(cond1))
                test1 = True
            cond2 = np.outer(self.mu_lrg,np.ones(2*self.p))+np.outer(np.ones(2*self.p),self.mu_lrg) -1 >= self.M
            cond2 = np.logical_and(np.diag(1-np.ones(2*self.p)),cond2)
            if (np.sum(cond2)) > 0:
                print "%d cross moments which were too low" % (np.sum(cond2))
                test2 = True

            if test1 or test2:
                if shrink>.1:
                    shrink = .9*shrink
                else:
                    test1, test2 = False, False
                    "Shrinking didn't solve the issue"

            if not (test1 or test2) and shrink<1:
                print "s had to be shrunk by a factor of %.2f" % shrink

    def _binary_knockoff(self):
        ''' This creates the new binary knockoffs, which are random multivariate bernoulli which should have, in expectation,
        the same first two moments as X. Only will work if X is all binaray '''

        self._derive_crossmoments()

        ####################################################
        # Get the data corresponding to the original x for the simulations 
        ####################################################

        A = np.zeros((2*self.p,2*self.p))

        # Simulate fresh x based on the original data
        if self.fresh_sim: 
            # Fit the upperhalf of A on the actual data, making this easier
            A[1,1] = logit(self.mu_lrg[0])
            for i in np.arange(0,self.p):
            # inject the paramters from the logit X_i ~ X_1 + ... + X_(i-1) + 1 into the ith row of A
                A[i,0:(i+1)] = sm.GLM(self.X_orig[:,i],np.hstack((self.X_orig[:,0:i],np.ones((self.n,1)))),family=sm.families.Binomial()).fit().params

            # Then draw the X
            X_fix = draw_random_binary(self.MCsize,A[:self.p,:self.p])
            nMC = self.MCsize
        
        # just repeat X a bunch of times
        else:
            # Rather than simulate entirely new X at each stage, I will used a fixed set of X_1 ... X_(i-1)
            # This definitely makes sense for the orignal X vars (why simulate when we already have it), but possibly less sense for the knockoffs
            # To get the desired size of montecarlo simulation, I will replicate X until it has at least self.MCsize rows
            repl = np.min((self.MCsize//self.n,1))
            X_fix = np.repeat(self.X_orig,repl,0)
            nMC   = X_fix.shape[0]

        ###################################################
        # Derive remaing A from Newton-Raphson
        ###################################################

        # Largely from 5.1 in Schafer, including notation

        # the current value and the derivatves are derived by simulation

        # sequence of portions between 0 and 1 to deal with case of ill conditioned hessian
        por_seq = np.arange(0,1,.25)
        self.por = np.empty((1,0))

        for i in np.arange(self.p,2*self.p):
            # m are the cross moments we are trying to fit
            m = self.M[i,0:(i+1)]
            # X_tmp is the list of (x_k,1) vectors
            X_fix     = np.hstack((X_fix,np.ones((nMC,1))))

            #X_tmp_out = X_fix[:,:,np.newaxis]*X_fix[:,np.newaxis,:]

            # Now, the Newton-Raphson steps
            # If the hessian becomes singular, we will relax the cross moment requirements, as described in Schafer 5.1 point 2
            # the idea is that the problem is relaxed until X_i is independent of all prior vars
            # as por increases, covariance drops
            # a is the row we are adding to A. Initialize with values as if independent all other vars
            a = np.append(np.zeros(i),logit(self.mu_lrg[i]))

            for por in por_seq:
                # m are the cross moments we are trying to fit
                m = (1-por)*self.M[i,0:(i+1)] + por*self.M[i,i]*np.append(np.diag(self.M)[0:i],1)

                # Minimize the actual difference vector
                opt = root(self._vector_objective,
                        a,
                        args=(X_fix,m),
                        method='anderson',
                        options = {'maxiter':(i*2+150),'fatol':1E-5,'jac_options':{'M':20}}
                        )

                # update a to most recent estimate, even without convergence
                a = opt.x

                # Stop once optimal has been reached
                if opt.success:
                    self.por = np.append(self.por,por)
                    if por>0:
                        print "Variable %d relaxed by tau=%.2f" % (i-self.p+1,por)
                    break

            if not opt.success:
                a = np.append(np.zeros(i),logit(self.mu_lrg[i]))
                self.por = np.append(self.por,1)
                print "Variable %d fully relaxed" % (i-self.p+1)

            # put a into A matrix, draw X_i for 'fixed' matrix, update X_out_fix
            A[i,0:(i+1)] = a
            X_fix[:,-1] = npran.binomial(1,invlogit(np.dot(X_fix,a)))

        ##############################################
        # Wrapup and get X_ko
        ##############################################


        # If we freshly simulated x, we need to draw ~x based on x
        if self.fresh_sim: 
            self.X_lrg = np.hstack((self.X_orig,np.empty((self.n,self.p))))
            for i in np.arange(self.p,2*self.p):
                # need to make sure the knockoff isn't uniformly 0 or 1
                count = 0
                j=0
                while count==0 or count==self.n:
                    # first five times we try regenerating
                    if j<5: 
                        self.X_lrg[:,i] = npran.binomial(1,invlogit(np.dot(np.hstack((self.X_lrg[:,0:i],np.ones((self.n,1)))),A[i,0:(i+1)])))
                        if j>0:
                            print "Knockoff regenerated to avoid constant value"
                    # otherwise, just randomly flip a few bits
                    else:
                        print "Random noise added to knockoff to avoid constant value"
                        self.X_lrg[:,i] = np.where(npran.binomial(1,.01*np.ones(self.n)),1-self.X_lrg[:,i],self.X_lrg[:,i])
                    count = np.sum(self.X_lrg[:,i])
                    j += 1

        else:
            # since we've been drawing the X along the way, can subset X_fix to get X_ko
            self.X_lrg = np.concatenate((self.X_orig,X_fix[0::repl,self.p:]), axis=1)

        # hang onto A
        self.A = A

        # Evaluate how close we are emperically to M
        self.M_distortion = nplin.norm(self.M-np.dot(self.X_lrg.T,self.X_lrg)/self.n)/nplin.norm(self.M)

        self.emp_ko_corr= np.corrcoef(self.X_lrg,rowvar=0,bias=1)[:self.p,self.p:2*self.p][np.identity(self.p)==1]
        if np.sum(np.isnan(self.emp_ko_corr))>0:
            print "There were %d out of %d variables who had missing correlation" % (np.sum(np.isnan(self.emp_ko_corr)),self.p)

    def _vector_objective(self,a,X_fix,m):
        return np.mean(invlogit(np.dot(X_fix,a))[:,np.newaxis]*X_fix,axis=0) - m

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

    def fit(self,X_lrg=None):
        """ Generates the knockoffs, fits the regression, and performs the FDR calculations """
        # Generate knockoff as inherited from knockoff_net
        if X_lrg is None:
            if   self.knockoff_type == 'original': self._original_knockoff()
            elif self.knockoff_type == 'binary':   self._binary_knockoff()
        else:
            self.X_lrg = X_lrg

        # initialize and fit the glmnet object
        self.lognet = LogisticNet(alpha=1) 
        self.lognet.fit(self.X_lrg,self.y,normalize=False,include_intercept=self.intercept)

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
    
    def fit(self,X_lrg=None):
        """ Generates the knockoffs, fits the regression, and performs the FDR calculations """
        # Generate knockoff as inherited from knockoff_net
        if X_lrg is None:
            if   self.knockoff_type == 'original': self._original_knockoff()
            elif self.knockoff_type == 'binary':   self._binary_knockoff()
        else:
            self.X_lrg = X_lrg

        # initialize the glmnet object
        self.elasticnet = ElasticNet(alpha=1) 
        self.elasticnet.fit(self.X_lrg,self.y,normalize=False,include_intercept=self.intercept)

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

