#imports
from __future__ import division
import numpy as np
import numpy.linalg as nplin
import numpy.random as npran
import scipy as sp
import scipy.linalg as splin
import scipy.stats as spstat1gt
import matplotlib.pyplot as plt
import knockoffGLM as ko
import simulate as sim
import statsmodels.api as sm
from statsmodels.tools.sm_exceptions import PerfectSeparationError
from multiprocessing import Pool
from scipy.special import expit as invlogit
from scipy.special import logit

def cutoff(array,min=-50,max=50):
    shape = array.shape
    cut_array = np.min((max*np.ones(shape),array),axis=0)
    cut_array = np.max((min*np.ones(shape),cut_array),axis=0)
    return cut_array

def make_X_ind(X,q=.05):
    """ Add some random noise to make X linearly independent """
    n,p = X.shape
    for i in range(p):
        while True:
            try:
                model = sm.GLM(self.X_orig[:,i],np.hstack((np.ones((n,1)),X[:,:i])),family=sm.families.Binomial()).fit()
            except:
                # Keep flipping a few bits until we have linear independence
                X[:,i] = np.where(npran.binomial(1,q*np.ones(n)),1-X[:,i],X[:,i])
                v = np.min(nplin.svd(np.hstack((np.ones((n,1)),X[:,:(i+1)])))[1])
                break
    return X

def ising_X(p,n,A_base_diag=-1,A_sd=.2):
    """ generate X from ising model """
    A = npran.normal(0,A_sd,(p,p))+np.diag(A_base_diag*np.ones(p))
    m,p = A.shape
    ones   = np.ones((n,1))
    X = np.empty((n,p))
    for i in np.arange(0,p):
        X[:,i] = npran.binomial(1,invlogit(np.dot(np.hstack((X[:,0:i],ones)),A[i,0:(i+1)])))
    return X

def given_X(p,n,data):
    h,w = data.shape
    X      = data[npran.choice(h,n),:][:,npran.choice(w,p)]
    return X

def norm_y(X,p1,sd=1,beta_sd=1):
    """ Generate a normal Y from the X """
    n,p = X.shape
    X_1    = X[:,:p1]
    beta   = npran.randn(p1)*beta_sd
    if p1>0:
        eta    = np.dot(X_1,beta)
        y      = eta + npran.normal(0,sd,n)
    else:
        y      = npran.normal(0,sd,n)
    return y

def bern_y(X,p1,base_prob=.25,beta_sd=1):
    n,p = X.shape
    X_1    = X[:,:p1]
    v = 0 
    while v<1E-5:
        beta   = npran.randn(p1)*beta_sd
        if p1>0:
            eta    = cutoff(np.dot(X_1,beta)+logit(base_prob))
            y      = npran.binomial(1,invlogit(eta),n)
        else:
            y      = npran.binomial(1,base_prob,n)
        v = np.min(nplin.svd(np.hstack((X,y[:,np.newaxis])))[1])
    return y

def generate(input):
    """ Testing function, all parameters are in a single tuple """
    seed,gen,p0,p1,n,method,MCsize= input
    npran.seed(seed)
    if gen.lower() == "ising":
        X    = ising_X(p1+p0,n)
        ynor = norm_y(X,p1)
    elif gen.lower() == "genetic":
        genes = np.genfromtxt('data/SNPdata.txt', delimiter=',')
        np.place(genes,genes!=0,1)
        X = given_X(p1+p0,n,genes)
    ybin = bern_y(X,p1)
    ynor = norm_y(X,p1)

    # Logit
    bin_logit = ko.knockoff_logit(ybin,X,.2,
                              knockoff='binary',
                              method=method,
                              MCsize=MCsize,
                              intercept=True
                              )
    bin_logit.fit()
    ori_logit = ko.knockoff_logit(ybin,X,.2,
                              knockoff='original',
                              intercept=False
                              )
    ori_logit.fit()
    trueS = (np.arange(p0+p1)<p1).astype(int)
    bin_FDR   = np.dot(bin_logit.S,1-trueS)/max(np.sum(bin_logit.S),1)
    bin_power = np.dot(bin_logit.S,trueS)  /max(p1,1)
    ori_FDR   = np.dot(ori_logit.S,1-trueS)/max(np.sum(ori_logit.S),1)
    ori_power = np.dot(ori_logit.S,trueS)  /max(p1,1)
    corr      = np.corrcoef(ori_logit.S,bin_logit.S)[0,1]
    ko_corr = [cor for cor in bin_logit.emp_ko_corr if not np.isnan(cor)]

    with open('data/logit_test_'+str(p0+p1)+'_w_n.txt','a') as f:
        f.write("%d\t%s\t%d\t%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n" % (seed, gen, p1, n, bin_logit.M_distortion, np.mean(ko_corr), bin_FDR, bin_power, np.mean(ori_logit.emp_ko_corr), ori_FDR, ori_power, corr))

    # LASSO
    bin_lasso = ko.knockoff_lasso(ynor,X,.2,
                              knockoff='binary',
                              method=method,
                              MCsize=MCsize,
                              intercept=True
                              )
    bin_lasso.fit(bin_logit.X_lrg)
    ori_lasso = ko.knockoff_lasso(ynor,X,.2,
                              knockoff='original',
                              intercept=False
                              )
    ori_lasso.fit()
    trueS = (np.arange(p0+p1)<p1).astype(int)
    bin_FDR   = np.dot(bin_lasso.S,1-trueS)/max(np.sum(bin_lasso.S),1)
    bin_power = np.dot(bin_lasso.S,trueS)  /max(p1,1)
    ori_FDR   = np.dot(ori_lasso.S,1-trueS)/max(np.sum(ori_lasso.S),1)
    ori_power = np.dot(ori_lasso.S,trueS)  /max(p1,1)
    corr      = np.corrcoef(ori_lasso.S,bin_lasso.S)[0,1]

    with open('data/lasso_test_'+str(p0+p1)+'_w_n.txt','a') as f:
        f.write("%d\t%s\t%d\t%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n" % (seed, gen, p1, n, bin_logit.M_distortion, np.mean(ko_corr), bin_FDR, bin_power, np.mean(ori_lasso.emp_ko_corr), ori_FDR, ori_power, corr))

def batch_compare(b,p,p1s,gens,start_seed,method='fresh_sim',MCsize=100000,procs=4):
    with open('data/backup_seeds.txt','r') as f:
        seeds = [int(seed) for seed in f.read().split()]
    inputs = []
    i = 0
    for b in range(b):
        for p1 in p1s:
            for gen in gens:
                inputs.append((seeds[i+start_seed],gen,p-p1,p1,method,MCsize))
                i += 1

    pool = Pool(processes=procs)
    pool.map(generate,inputs)
    pool.close()

def distortion(input):
    """ Testing function, all parameters are in a single tuple """
    seed,p,n= input
    npran.seed(seed)
    X    = ising_X(p,n)
    ybin = bern_y(X,1)

    if   p<=10:
        MCsize =10000
    elif p<=20:
        MCsize =30000
    elif p<=30:     
        MCsize =70000
    elif p<=40:    
        MCsize =100000
    else:    
        MCsize =150000

    # Logit
    bin_logit = ko.knockoff_logit(ybin,X,.2,
                              knockoff='binary',
                              method='fresh_sim',
                              MCsize=MCsize,
                              intercept=True
                              )
    bin_logit._binary_knockoff()

    with open('data/distortion.txt','a') as f:
        f.write("%d\t%d\t%d\t%.5f\n" % (seed, p, n, bin_logit.M_distortion))

def batch_distortion(b,ps,ns,start_seed,procs=4):
    with open('data/backup_seeds.txt','r') as f:
        seeds = [int(seed) for seed in f.read().split()]
    inputs = []
    i = 0
    for b in range(b):
        for p in ps:
            for n in ns:
                inputs.append((seeds[i+start_seed],p,n))
                i += 1

    pool = Pool(processes=procs)
    pool.map(distortion,inputs)
    pool.close()



