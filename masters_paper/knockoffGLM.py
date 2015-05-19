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
from glmnet import LogisticNet

def get_index_z(coef_vector):
    """This returns the first index for a variable from the glmnet coef matrix which isn't zero"""
    n     = coef_vector.size
    index = 0 
    while index<n: 
        if coef_vector[index]!=0: return index
        index+=1
    return np.nan 

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
        self.y = y
        self.X = X
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

    def _fit_lognet(self):
        X_lrg = np.concatenate((self.X,self.X_ko), axis=1)

        # initialize the glmnet object
        self.lognet = LogisticNet(alpha=1) 
        self.lognet.fit(X_lrg,self.y)

        self.lambdas = self.lognet.out_lambdas
        self.var_index_ent = np.sort(self.lognet._indices)
        self.coef_matrix = np.zeros((2*self.p,self.lognet.n_lambdas))
        self.coef_matrix[self.var_index_ent] = self.lognet._comp_coef.squeeze()[self.lognet._indices]

        self.var_entered = np.zeros(2*self.p).astype(bool)
        self.var_entered[self.var_index_ent] = True

    def _get_z(self): 
        """ Given the coefficient matrix from a glmnet object, returns the Z, the first non-zero entries"""
        self.z_rank  = np.array(map(get_index_z,self.coef_matrix))
        self.z_value = np.array(map(lambda x: self.lambdas[x] if not np.isnan(x) else 0,self.z_rank))

    def _get_w(self): 
        """Produces the w values using the z values"""
        self.w_filter = np.sign(self.z_value[0:self.p]-self.z_value[self.p:(2*self.p)])
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

