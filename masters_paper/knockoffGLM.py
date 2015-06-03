""" This module is meant to implement and test the knockoff procedure on L1 regularized GLM models. The method was originally designed for linear least squares regression, and there is no theoretical guarantee it will work on other GLMS. Hopefully, with simulation we can discover if it breaks, if so where it breaks, and maybe how to fix it in those situations.

More can be read about the original procedure on least squares linear models here: http://web.stanford.edu/~candes/Knockoffs/index.html

Most of the code is stolen/modified from the knockoff R package

"""

# imports 
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
from rpy2.robjects.packages import importr
from glmnet import LogisticNet
from sklearn.preprocessing import normalize

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

class knockoff_logit(object):
    def __init__(self,y,X,knockoff='equicor',randomize=False):
        self.y      = y
        self.X      = normalize(X.astype(float),norm='l2',axis=0)
        self.X_orig = X
        self.n,self.p = self.X.shape
        self.knockoff_type = knockoff
        self.randomize=False
        self.tol = 1E-5

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
        s[s <= self.tol] = 0
 
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

    def _binary_knockoff(self):
        ''' This creates knockoffs which are random multivariate bernoulli which should have, in expectation,
        the same first two moments as X. Only will work if X is all binaray '''
        # First, calculate s
        # Largely replicates begining of self._create_equicor()/self._create_SDP()

        # SVD and come up with perpendicular matrix
        U, d, V = nplin.svd(self.X,full_matrices=True) 
        d[d<0] = 0
        if   self.knockoff_type == 'SDP': 
            # Compute the Gram matrix and its (pseudo)inverse.
            G     = np.dot(V.T * d**2 ,V)
            G_inv = np.dot(V.T * d**-2,V)
         
            # Optimize the parameter s of Equation 1.3 using SDP.
            s = solve_sdp(G)
            s[s <= self.tol] = 0

        # Same as first part self._create_equicor 
        elif self.knockoff_type == 'equicor':
            # Set s = min(2 * smallest eigenvalue of X'X, 1), so that all the correlations
            # have the same value 1-s.
            lambda_min = min(d)**2
            s = np.ones(self.p)*min(2*lambda_min,1)

        # Generate the large covariance matrix we want
        self.X = self.X_orig              # Now that we've adjusted the variance, revert back to the original X 
        Sigma     = np.cov(self.X,rowvar=0)
        s         = s*np.diag(Sigma)        # scale s back up to original variance
        Sigma_lrg = np.hstack(((np.vstack((Sigma,Sigma-np.diag(s)))),np.vstack((Sigma-np.diag(s),Sigma))))

        # Generate the mean vector
        mu       = np.mean(self.X,axis=0)
        mu_lrg   = np.concatenate((mu,mu))

        # Get a covariance matrix so we can generate random normals that when thresholded provide random binaries with mean and covariance matchin mu_lrg, Sigma_lrg 
        # Use the r package bindata
        rpy2.robjects.numpy2ri.activate()
        bindata    = importr('bindata')
        commonprob = bindata.bincorr2commonprob(mu_lrg,Sigma_lrg)
        Sigma_nor  = np.asarray(bindata.commonprob2sigma(commonprob))
        mu_nor     = spstat.norm.ppf(mu_lrg)

        lam,v      = nplin.eig(Sigma_nor)
        print     lam

        # Now, draw X_ko such that (X,X_ko) has first and second moments mu_lrg, Sigma_lrg
        # This is accomplished by performing gibbs sampling based on Z~N(mu_nor,mu_lrg), with the first p terms conditioned to reflect z[i]<mu[i] if x[i]=0, or z[i]>mu[i] if x[i]=1. 
        Z = mu_lrg + np.concatenate((self.X,self.X), axis=1) - .5     # inital value, set Z to .5 above/below upper/lower bound
        # Can precompute the matricies for taking the conditional normal
        S_yxS_xx = []       # Sigma_YX * Sigma_XX^-1
        sig_c    = []       # Sigma_Y|X
        for k in range(2*self.p):
            nk   = np.arange(2*self.p)!=k
            S_yx = Sigma_nor[k,nk]
            S_xx = nplin.inv(Sigma_nor[nk,:][:,nk])    # not sure why subsetting rows and columns doesn't work
            S_yxS_xx.append(np.dot(S_yx,S_xx))
            sig_c.append(Sigma_nor[k,k] - np.dot(S_yx,S_yxS_xx[k]))

        # Loop over iterations
        for i in range(self.gibbs_iter):
            # loop through variables
            for k in range(self.p):
                nk   = np.arange(2*self.p)!=k
                mu_cond  = mu_nor[k] + np.dot(S_yxS_xx[k],(Z[:,nk]-mu_nor[nk]).T)
                # If they are originals, we need to condition truncate above/below mu_nor[k]
                if k//self.p==0:
                    # The npran.random is uniform on [0,1), so need to always map 0 to mu_lrg (the mean of the original bernoulli)
                    qcut   = spstat.norm.cdf((mu_nor[k]-mu_cond)/sig_c[k])
                    quan   = qcut + npran.random(size=self.n)*(self.X[:,k] - qcut)
                    quan   = np.clip(quan,0+self.tol,1-self.tol)      # Make sure to avoid numerical issues
                    Z[:,k] = spstat.norm.ppf(quan)*sig_c[k] + mu_cond
                    Z[:,k] = np.where(self.X[:,k],
                                        np.clip(Z[:,k],mu_nor[k]+self.tol,np.inf            ),
                                        np.clip(Z[:,k],-np.inf           ,mu_nor[k]-self.tol)
                                     )        # more numerical safeguards
                # If they aren't original, no truncation
                else:
                    Z[:,k] = npran.normal(size=self.n)*sig_c[k] + mu_cond
        self.X_ko = Z[:,self.p:]>0



    def _fit_lognet(self):
        X_lrg = np.concatenate((self.X,self.X_ko), axis=1)

        # initialize the glmnet object
        self.lognet = LogisticNet(alpha=1) 
        self.lognet.fit(X_lrg,self.y,normalize=False,include_intercept=False)

        self.lambdas = self.lognet.out_lambdas
        self.var_index_ent = np.sort(self.lognet._indices)
        self.coef_matrix = np.zeros((2*self.p,self.lognet.n_lambdas))
        self.coef_matrix[self.var_index_ent] = self.lognet._comp_coef.squeeze()[self.lognet._indices]

        self.var_entered = np.zeros(2*self.p).astype(bool)
        self.var_entered[self.var_index_ent] = True

    def _get_z(self): 
        """ Given the coefficient matrix from a glmnet object, returns the Z, the first non-zero entries"""
        self.z_rank  = np.array(map(get_index_z,self.coef_matrix))
        self.z_value = np.array(map(lambda x: self.lambdas[x] if x<self.lambdas.size else 0,self.z_rank))

    def _get_w(self): 
        """Produces the w values using the z values"""
        self.w_filter = np.sign(self.z_rank[self.p:(2*self.p)]-self.z_rank[0:self.p])
        self.w_rank   = self.w_filter * np.nanmin((self.z_rank[self.p:(2*self.p)],self.z_rank[:self.p]),axis=0)
        self.w_value  = self.w_filter * np.max((self.z_value[self.p:(2*self.p)],self.z_value[:self.p]),axis=0)

    def fit(self):
        if   self.knockoff_type == 'equicor': self._create_equicor()
        elif self.knockoff_type == 'SDP'    : self._create_SDP()
        self._fit_lognet()
        self._get_z()
        self._get_w()

def main(n,p,q):
    print 'nobody lives here right now'

if __name__ == '___main__':
    status = main()
    sys.exit(status)

