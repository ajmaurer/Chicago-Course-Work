# Imports
from __future__ import division
import numpy as np
import numpy.linalg as nplin
import numpy.random as npran
import scipy as sp
import scipy.linalg as splin
import scipy.stats as spstat
import matplotlib.pyplot as plt
import urllib3, certifi
from scipy.special import expit as invlogit
from scipy.special import logit

from multiprocessing import Pool
from bs4 import BeautifulSoup

import knockoffGLM as ko

def draw_random_binary(n,A):
    ''' If p is the size of the square, lower triangular A, creates a n by p matrix of random binary vectors using Schafer's method '''
    m,p = A.shape
    ones   = np.ones((n,1))
    output = np.empty((n,p))
    for i in np.arange(0,p):
        output[:,i] = npran.binomial(1,invlogit(np.dot(np.hstack((output[:,0:i],ones)),A[i,0:(i+1)])))
    return output

def getRandomIntegers(num=1000,min=1,max=1000000000,base=10):
    """
    This will query random.org to get a bunch of truly random integers. I want this so I can get seeds for parallel pseudo-random number generation, since you don't want to use a pseudo-random numbers or a sequence as seeds.
    
    Since I'm not using random.org's API (do I look like I'm made out of money?!?), a maximum of 10,000 random numbers can be generated, and they min and max must be within +/- 1,000,000,000
    """

    http = urllib3.PoolManager(
        cert_reqs='CERT_REQUIRED', # Force certificate check.
        ca_certs=certifi.where(),  # Path to the Certifi bundle.
    )
    address = 'https://www.random.org/integers/?num=%d&min=%d&max=%d&col=1&base=%d&format=html&rnd=new' % (num,min,max,base)
    r       = http.request('get',address)
    uniList = BeautifulSoup(r.text).find('pre',class_='data').text.split()
    numList = map(int,uniList)
    return numList

def cutoff(array,min=-50,max=50):
    shape = array.shape
    cut_array = np.min((max*np.ones(shape),array),axis=0)
    cut_array = np.max((min*np.ones(shape),cut_array),axis=0)
    return cut_array

def genXy_given_X_norm_y(seed,data,n,p1,pnull,sd=1,beta_sd=1):
    ''' X is binary from the isling model, with the coefficients drawn from a normal. Y is binary, with beta's coefficients also from a normal '''
    if not seed == None:
        npran.seed(seed)
    p = p1 + pnull
    h,w = data.shape
    rows = npran.choice(h,n)
    X      = data[rows,:][:,npran.choice(w,p)]
    X_1    = X[:,:p1]
    X_null = X[:,p1:]
    beta   = npran.randn(p1)*beta_sd
    if p1>0:
        eta    = np.dot(X_1,beta)
        y      = eta + npran.normal(0,sd,n)
    else:
        y      = npran.normal(0,sd,n)
    return X,y

def genXy_binary_X_norm_y(seed,n,p1,pnull,sd=1,beta_sd=1,A_base_diag=-1,A_sd=.2):
    ''' X is binary from the isling model, with the coefficients drawn from a normal. Y is binary, with beta's coefficients also from a normal '''
    if not seed == None:
        npran.seed(seed)
    p = p1 + pnull
    A = npran.normal(0,.2,(p,p))-np.diag(A_base_diag*np.ones(p))
    X = draw_random_binary(n,A)
    X_1    = X[:,:p1]
    X_null = X[:,p1:]
    beta   = npran.randn(p1)*beta_sd
    if p1>0:
        eta    = np.dot(X_1,beta)
        y      = eta + npran.normal(0,sd,n)
    else:
        y      = npran.normal(0,sd,n)
    return X,y




def genXy_given_X_norm_beta(seed,data,n,p1,pnull,base_prob=.25,beta_sd=1):
    ''' X is binary from the isling model, with the coefficients drawn from a normal. Y is binary, with beta's coefficients also from a normal '''
    if not seed == None:
        npran.seed(seed)
    p = p1 + pnull
    h,w = data.shape
    rows = npran.choice(h,n)
    X      = data[rows,:][:,npran.choice(w,p)]
    X_1    = X[:,:p1]
    X_null = X[:,p1:]
    beta   = npran.randn(p1)*beta_sd
    if p1>0:
        eta    = cutoff(np.dot(X_1,beta)+logit(base_prob))
        y      = npran.binomial(1,invlogit(eta),n)
    else:
        y      = npran.binomial(1,base_prob,n)
    return X,y

def genXy_binary_X_norm_beta(seed,n,p1,pnull,base_prob=.25,beta_sd=1,A_base_diag=-1,A_sd=.2):
    ''' X is binary from the isling model, with the coefficients drawn from a normal. Y is binary, with beta's coefficients also from a normal '''
    if not seed == None:
        npran.seed(seed)
    p = p1 + pnull
    A = npran.normal(0,.2,(p,p))-np.diag(A_base_diag*np.ones(p))
    X = draw_random_binary(n,A)
    X_1    = X[:,:p1]
    X_null = X[:,p1:]
    beta   = npran.randn(p1)*beta_sd
    if p1>0:
        eta    = cutoff(np.dot(X_1,beta)+logit(base_prob))
        y      = npran.binomial(1,invlogit(eta),n)
    else:
        y      = npran.binomial(1,base_prob,n)
    return X,y

def genXy_binary_X_norm_beta(seed,n,p1,pnull,base_prob=.25,beta_sd=1,A_base_diag=-1,A_sd=.2):
    ''' X is binary from the isling model, with the coefficients drawn from a normal. Y is binary, with beta's coefficients also from a normal '''
    if not seed == None:
        npran.seed(seed)
    p = p1 + pnull
    A = npran.normal(0,.2,(p,p))-np.diag(A_base_diag*np.ones(p))
    X = draw_random_binary(n,A)
    X_1    = X[:,:p1]
    X_null = X[:,p1:]
    beta   = npran.randn(p1)*beta_sd
    if p1>0:
        eta    = cutoff(np.dot(X_1,beta)+logit(base_prob))
        y      = npran.binomial(1,invlogit(eta),n)
    else:
        y      = npran.binomial(1,base_prob,n)
    return X,y

def genXy_normal_X_beta(seed,n,p1,pnull,base_prob=.25,beta_sd=1):
    """ The X are normal. p1 predictive vars, pnull null vars. beta on the p1 vars is ~normal(0,beta_sd) and the intercept is logit(base_prob)"""
    if not seed == None:
        npran.seed(seed)
    X_1    = npran.randn(n,p1)
    X_null = npran.randn(n,pnull)
    X      = np.concatenate((X_1,X_null),axis=1)
    beta   = npran.randn(p1)*beta_sd
    if p1>0:
        eta    = cutoff(np.dot(X_1,beta)+logit(base_prob))
        y      = npran.binomial(1,invlogit(eta),n)
    else:
        y      = npran.binomial(1,base_prob,n)
    return X,y

def genXy_bern_X_norm_beta(seed,n,p1,pnull,x_prob=.25,base_prob=.25,beta_sd=1):
    """ The X are normal. p1 predictive vars, pnull null vars. beta on the p1 vars is ~normal(0,beta_sd) and the intercept is logit(base_prob)"""
    if not seed == None:
        npran.seed(seed)
    X_1    = npran.binomial(1,x_prob,(n,p1))
    X_null = npran.binomial(1,x_prob,(n,pnull))
    X      = np.concatenate((X_1,X_null),axis=1)
    beta   = npran.randn(p1)*beta_sd
    if p1>0:
        eta    = cutoff(np.dot(X_1,beta)+logit(base_prob))
        y      = npran.binomial(1,invlogit(eta),n)
    else:
        y      = npran.binomial(1,base_prob,n)
    return X,y

# these helpers are here because pool.map() is stupid and won't take additional arguments
def merge_two_dicts(x,y):
    '''Given two dicts, merge them into a new dict as a shallow copy.'''
    z = x.copy()
    z.update(y)
    return z

def removekeys(d, keys):
    r = dict(d)
    for key in keys:
        del r[key]
    return r

# generates a whole lot of w statistics from different simulations
def gen_w_s(kwargs):
    func = kwargs['func']
    X,y = func(**removekeys(kwargs,['func','q','knockoff','model']))
    if kwargs['knockoff']=='both':
        if kwargs['model']=='logit':
            bin_knockoff_model = ko.knockoff_logit(y,X,kwargs['q'],knockoff='binary')
            ori_knockoff_model = ko.knockoff_logit(y,X,kwargs['q'],knockoff='original')
        elif kwargs['model']=='lasso':
            bin_knockoff_model = ko.knockoff_lasso(y,X,kwargs['q'],knockoff='binary')
            ori_knockoff_model = ko.knockoff_lasso(y,X,kwargs['q'],knockoff='original')
        bin_knockoff_model.fit()
        ori_knockoff_model.fit()
        return ((bin_knockoff_model.w_value,ori_knockoff_model.w_value),(bin_knockoff_model.S,ori_knockoff_model.S))
    else:
        if kwargs['model']=='logit':
            knockoff_model = ko.knockoff_logit(y,X,kwargs['q'],knockoff=kwargs['knockoff'])
        elif kwargs['model']=='lasso':
            knockoff_model = ko.knockoff_lasso(y,X,kwargs['q'],knockoff=kwargs['knockoff'])
        knockoff_model.fit()
        return (knockoff_model.w_value, knockoff_model.S.astype(int))

def batch_gen_w_s(procs,seeds,**kwargs):
    ''' Will calculate a w statistic using the key words arguments passed to gen_w for each seed using procs processors'''
    pool = Pool(processes=procs)
    input = [merge_two_dicts({'seed':seed},kwargs) for seed in seeds]
    outlist   = pool.map(gen_w_s,input)
    pool.close()
    Ws    = [w for w,s in outlist]
    Ss    = [s for w,s in outlist]
    return (Ws,Ss)

class ko_test(object):
    def __init__(self,procs,seeds,**kwargs):
        self.k      = len(seeds)
        self.kwargs = kwargs
        # This setup is akward, but keeping most of the arguments in the dictionary saves time with passing them
        # to batch_gen_w
        self.n      = kwargs['n']
        self.p1     = kwargs['p1']
        self.pnull  = kwargs['pnull']
        self.p      = self.p1+self.pnull
        # The Ws - matrix, with the first row the value, second the rank by absolute value,
        #    and the third the sign (did knockoff or original come in first?), and the forth whether it is a true predictor
        Ws,Ss   = batch_gen_w_s(procs,seeds,**kwargs)
        if kwargs['knockoff']!='both':
            self.Ws =   [
                            np.vstack(
                                (
                                    W,
                                    np.abs(W).argsort()[::-1].argsort(),
                                    np.sign(W),
                                    np.arange(self.p)<self.p1
                                )
                            )
                        for W in Ws]
        self.Ss = Ss

    def knockoff_rank_rate(self):
        if self.pnull>0:
            # Null knockoffs
            null_ko_w_ranks = np.concatenate([W[1,np.all((W[2,:]<0,W[3,:]==0),axis=0)] for W in self.Ws]).astype(int)
            self.null_ko_rank_ct = np.cumsum(np.bincount(null_ko_w_ranks,minlength=self.p))
            # Null Originals
            null_or_w_ranks = np.concatenate([W[1,np.all((W[2,:]>0,W[3,:]==0),axis=0)] for W in self.Ws]).astype(int)
            self.null_or_rank_ct = np.cumsum(np.bincount(null_or_w_ranks,minlength=self.p))

            # The cumulitive rate of originals coming in first up to a given rank
            self.null_orig_rank_rate = self.null_or_rank_ct/(self.null_or_rank_ct+self.null_ko_rank_ct)

        if self.p1>0:
            # True knockoffs
            true_ko_w_ranks = np.concatenate([W[1,np.all((W[2,:]<0,W[3,:]==1),axis=0)] for W in self.Ws]).astype(int)
            self.true_ko_rank_ct = np.cumsum(np.bincount(true_ko_w_ranks,minlength=self.p))
            # True Originals
            true_or_w_ranks = np.concatenate([W[1,np.all((W[2,:]>0,W[3,:]==1),axis=0)] for W in self.Ws]).astype(int)
            self.true_or_rank_ct = np.cumsum(np.bincount(true_or_w_ranks,minlength=self.p))

            # The cumulitive rate of originals coming in first up to a given rank
            self.true_orig_rank_rate = self.true_or_rank_ct/(self.true_or_rank_ct+self.true_ko_rank_ct)

    def plot_ko_rank_rate(self,title='Percentage Original Variables',plot=None,xlabel=True,ylabel=True):
        self.knockoff_rank_rate()
        if plot==None:
            if self.p1>0:
                plt.plot(np.arange(self.p)+1,self.true_orig_rank_rate,'b-', label = 'True Predictors')
            if self.pnull>0:
                plt.plot(np.arange(self.p)+1,self.null_orig_rank_rate,'g-', label = 'Null Predictors')
            end_val = self.null_orig_rank_rate[-1]
            if xlabel:
                plt.xlabel('Absolute Rank W')
            if ylabel:
                plt.ylabel('Cumulative % Positive W')
            plt.title(title+'\n %.3f Total Positive W' % end_val)
            plt.axis([1,self.p,0,1])
            plt.legend(loc='lower right')
            plt.show()

        else:
            if self.p1>0:
                plot.plot(np.arange(self.p)+1,self.true_orig_rank_rate,'b-', label = 'True Predictors')
            if self.pnull>0:
                plot.plot(np.arange(self.p)+1,self.null_orig_rank_rate,'g-', label = 'Null Predictors')
            end_val = self.null_orig_rank_rate[-1]
            if xlabel:
                plot.set_xlabel('Absolute Rank W')
            if ylabel:
                plot.set_ylabel('Cumulative % Positive W')
            plot.set_title(title+'\n %.3f Total Positive W' % end_val)
            plot.set_xlim(1,self.p)
            plot.set_ylim(0,1)
            plot.legend(loc='lower right')

    def calc_fdr_power(self):
        trueS = (np.arange(self.p)<self.p1).astype(int)
        if self.kwargs['knockoff']=='both':
            self.binFDRs  = map(lambda (Sbin,Sori): np.dot(Sbin,1-trueS)/max(np.sum(Sbin),1),self.Ss)
            self.binpowers= map(lambda (Sbin,Sori): np.dot(Sbin,trueS)  /max(self.p1,1),self.Ss)
            self.oriFDRs  = map(lambda (Sbin,Sori): np.dot(Sori,1-trueS)/max(np.sum(Sori),1),self.Ss)
            self.oripowers= map(lambda (Sbin,Sori): np.dot(Sori,trueS)  /max(self.p1,1),self.Ss)
            self.corr     = map(lambda (Sbin,Sori): np.corrcoef(Sbin,Sori)[0,1],self.Ss)
            self.corr     = [ c for c in self.corr if not np.isnan(c)]
            return (np.mean(self.binFDRs),np.mean(self.oriFDRs),np.mean(self.binpowers),np.mean(self.oripowers),np.mean(self.corr))


        else:
            self.FDRs  = map(lambda S: np.dot(S,1-trueS)/max(np.sum(S),1),self.Ss)
            self.powers  = map(lambda S: np.dot(S,trueS)/max(self.p1,1),self.Ss)
            return (np.mean(self.FDRs),np.mean(self.powers))

    






 
 
        


        
        
    

