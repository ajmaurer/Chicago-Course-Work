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
from sklearn import linear_model as lm
from sklearn.svm import l1_min_c 

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

def create_SDP(X,randomize):
    """ Creates the SDP knockoff of X"""

    # Check for rank deficiency (will add later).
    tol = 1e-5

    # SVD and come up with perpendicular matrix
    n,p = X.shape
    U, d, V = nplin.svd(X,full_matrices=True) 
    d[d<0] = 0
    U_perp = U[:,p:(2*p)]
    if randomize:
        U_perp = np.dot(U_perp,splin.orth(npran.randn(p,p)))

    # Compute the Gram matrix and its (pseudo)inverse.
    G     = np.dot(V.T * d**2 ,V)
    G_inv = np.dot(V.T * d**-2,V)

    # Optimize the parameter s of Equation 1.3 using SDP.
    s = solve_sdp(G)
    s[s <= tol] = 0

    # Construct the knockoff according to Equation 1.4:
    C_U,C_d,C_V = nplin.svd(2*np.diag(s) - (s * G_inv.T).T * s)
    C_d[C_d < 0] = 0
    X_ko = X - np.dot(X,G_inv*s) + np.dot(U_perp*np.sqrt(C_d),C_V)
    return(X_ko)

def create_equicor(X, randomize):
    """ Creates the equal correlation knockoff of X"""
    # Check for rank deficiency (will add later).
    tol = 1e-5

    # SVD and come up with perpendicular matrix
    n,p = X.shape
    U, d, V = nplin.svd(X,full_matrices=True) 
    d[d<0] = 0
    U_perp = U[:,p:(2*p)]
    U = U[:,:p]
    if randomize:
        U_perp = np.dot(U_perp,splin.orth(npran.randn(p,p)))

    # Set s = min(2 * smallest eigenvalue of X'X, 1), so that all the correlations
    # have the same value 1-s.
    lambda_min = min(d)**2
    s = min(2*lambda_min,1)

    # Construct the knockoff according to Equation 1.4.
    s_diff = 2*s - (s/d)**2
    s_diff[s_diff<0]=0 # can be negative due to numerical error
    X_ko = np.dot(U*(d-s/d) + U_perp*(np.sqrt(s_diff)) , V)
    return(X_ko)

def find_w(y,X,X_ko):
    """ Generate the W statistics given the X matrix, X knockoff, and y """
    n,p = X.shape

    X_lrg = np.concatenate((X,X_ko), axis=1)

    clen = p*10
    cs = l1_min_c(X_lrg, y, loss='log') * 10 * np.logspace(0, 3,clen)
    clf = lm.LogisticRegression(C=1.0, penalty='l1', tol=1e-6)
    coefs_ = []
    for c in cs:
        clf.set_params(C=c)
        clf.fit(X_lrg, y)
        coefs_.append(clf.coef_.ravel().copy())

    coefsMat = np.vstack(coefs_)
    coefsMat[abs(coefsMat)<1e-6] = 0
    Z   = np.zeros(2*p) 
    for i in range(clen):
        Z[np.all(np.vstack((Z==0,coefsMat[i,:]!=0)),axis=0)] = cs[i]
    w = np.min((Z[:p],Z[p:2*p]),axis=0) * np.sign(Z[p:2*p]-Z[:p])

    return w 

def generate_logit_w(y,X,knockoff='equicor',randomize=False):
    """ Generates the w statistic for a logistic model predicting y~X"""

    # normalize X
    X = X/np.linalg.norm(X,ord=2,axis=0)

    # creat knockoff
    if   knockoff == 'equicor': X_ko = create_equicor(X,randomize)
    elif knockoff == 'SDP'    : X_ko =     create_SDP(X,randomize)

    return find_w(y,X,X_ko)

def iid_normal_null_test(n,p,q):
    """ Test run with n observations. p covariates, each iid Normal. Base rate q"""
    y = npran.binomial(1,.2,n)
    X = npran.randn(n,p)
    w = generate_logit_w(y,X)
    print '%.2f of null knockoffs came in first; %.2f of variables never entered' % (np.sum(w>0)/float(np.sum(w!=0)),np.sum(w==0)/float(p))

def main(n,p,q):
    print 'nobody lives here right now'

if __name__ == '___main__':
    status = main()
    sys.exit(status)

