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
from multiprocessing import Pool

def gen_logit(input):
    """ Testing function, all parameters are in a single tuple """
    seed,gen,p0,p1 = input
    if gen.lower() == "ising":
        X,y = sim.genXy_binary_X_norm_beta(seed,1000,p1,p0,base_prob=.25,beta_sd=1,A_base_diag=-1,A_sd=.2)
    elif gen.lower() == "genetic":
        genes = np.genfromtxt('data/SNPdata.txt', delimiter=',')
        np.place(genes,genes!=0,1)
        X,y = sim.genXy_given_X_norm_beta(seed,genes,1000,p1,p0,base_prob=.25,beta_sd=1)

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
    if gen.lower() == "ising":
        X,y = sim.genXy_binary_X_norm_y(seed,1000,p1,p0,beta_sd=1,A_base_diag=-1,A_sd=.2)
    elif gen.lower() == "genetic":
        genes = np.genfromtxt('data/SNPdata.txt', delimiter=',')
        np.place(genes,genes!=0,1)
        X,y = sim.genXy_given_X_norm_y(seed,genes,1000,p1,p0,beta_sd=1)

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

def batch_lasso(b,p,p1s,procs=4):
    with open('data/backup_seeds.txt','r') as f:
        seeds = [int(seed) for seed in f.read().split()]
    inputs = []
    i = 0
    for b in range(b):
        for p1 in p1s:
            for gen in ['ising','genetic']:
                inputs.append((seeds[i],gen,p-p1,p1))
                i += 1

    pool = Pool(processes=procs)
    pool.map(gen_lasso,inputs)
    pool.close()

def batch_logit(b,p,p1s,procs=4):
    with open('data/backup_seeds.txt','r') as f:
        seeds = [int(seed) for seed in f.read().split()]
    inputs = []
    i = 0
    for b in range(b):
        for p1 in p1s:
            for gen in ['ising','genetic']:
                inputs.append((seeds[i],gen,p-p1,p1))
                i += 1

    pool = Pool(processes=procs)
    pool.map(gen_logit,inputs)
    pool.close()



