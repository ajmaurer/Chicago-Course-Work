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
    A = npran.normal(0,.2,(p,p))-np.diag(A_base_diag*np.ones(p))
    m,p = A.shape
    ones   = np.ones((n,1))
    X = np.empty((n,p))
    for i in np.arange(0,p):
        X[:,i] = npran.binomial(1,invlogit(np.dot(np.hstack((X[:,0:i],ones)),A[i,0:(i+1)])))
    return make_X_ind(X)

def given_X(p,n,data):
    h,w = data.shape
    X      = data[npran.choice(h,n),:][:,npran.choice(w,p)]
    return make_X_ind(X)

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

def gen_logit(input):
    """ Testing function, all parameters are in a single tuple """
    seed,gen,p0,p1 = input
    npran.seed(seed)
    if gen.lower() == "ising":
        X = ising_X(p1+p0,1000)
        y = bern_y(X,p1)
    elif gen.lower() == "genetic":
        genes = np.genfromtxt('data/SNPdata.txt', delimiter=',')
        np.place(genes,genes!=0,1)
        X = given_X(p1+p0,1000,genes)
        y = bern_y(X,p1)

    model = ko.knockoff_logit(y,X,.2,
                              knockoff='binary',
                              MCsize=25000,
                              tol=1E-5,
                              fresh_sim=True,
                              opt_method='anderson',
                              intercept=True
                              )
    model.fit()
    trueS = (np.arange(p0+p1)<p1).astype(int)
    FDR   = np.dot(model.S,1-trueS)/max(np.sum(model.S),1)
    power = np.dot(model.S,trueS)  /max(p1,1)

    with open('data/logit_test_'+str(p0+p1)+'.txt','a') as f:
        f.write("%d\t%s\t%d\t%.5f\t%.5f\t%.5f\t%.5f\n" % (seed,gen,p1,model.M_distortion,np.mean(model.emp_ko_corr),FDR,power))

def gen_lasso(input):
    """ Testing function, all parameters are in a single tuple """
    seed,gen,p0,p1 = input
    npran.seed(seed)
    if gen.lower() == "ising":
        X = ising_X(p1+p0,1000)
        y = norm_y(X,p1)
    elif gen.lower() == "genetic":
        genes = np.genfromtxt('data/SNPdata.txt', delimiter=',')
        np.place(genes,genes!=0,1)
        X = given_X(p1+p0,1000,genes)
        y = norm_y(X,p1)

    bin_model = ko.knockoff_lasso(y,X,.2,
                              knockoff='binary',
                              MCsize=25000,
                              tol=1E-5,
                              fresh_sim=True,
                              opt_method='anderson',
                              intercept=True
                              )
    bin_model.fit()
    ori_model = ko.knockoff_lasso(y,X,.2,
                              knockoff='original',
                              tol=1E-5,
                              intercept=False
                              )
    ori_model.fit()
    trueS = (np.arange(p0+p1)<p1).astype(int)
    bin_FDR   = np.dot(bin_model.S,1-trueS)/max(np.sum(bin_model.S),1)
    bin_power = np.dot(bin_model.S,trueS)  /max(p1,1)
    ori_FDR   = np.dot(ori_model.S,1-trueS)/max(np.sum(ori_model.S),1)
    ori_power = np.dot(ori_model.S,trueS)  /max(p1,1)
    corr      = np.corrcoef(ori_model.S,bin_model.S)[0,1]

    with open('data/lasso_test_'+str(p0+p1)+'.txt','a') as f:
        f.write("%d\t%s\t%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n" % (seed, gen, p1, bin_model.M_distortion, np.mean(bin_model.emp_ko_corr), bin_FDR, bin_power, np.mean(ori_model.emp_ko_corr), ori_FDR, ori_power, corr))

def batch_lasso(b,p,p1s,gens,start_seed,procs=4):
    with open('data/backup_seeds.txt','r') as f:
        seeds = [int(seed) for seed in f.read().split()]
    inputs = []
    i = 0
    for b in range(b):
        for p1 in p1s:
            for gen in gens:
                inputs.append((seeds[i+start_seed],gen,p-p1,p1))
                i += 1

    pool = Pool(processes=procs)
    pool.map(gen_lasso,inputs)
    pool.close()

def batch_logit(b,p,p1s,gens,start_seed,procs=4):
    with open('data/backup_seeds.txt','r') as f:
        seeds = [int(seed) for seed in f.read().split()]
    inputs = []
    i = 0
    for b in range(b):
        for p1 in p1s:
            for gen in gens:
                inputs.append((seeds[i+start_seed],gen,p-p1,p1))
                i += 1

    pool = Pool(processes=procs)
    pool.map(gen_logit,inputs)
    pool.close()



