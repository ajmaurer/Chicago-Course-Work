#imports
from __future__ import division
import numpy as np
import numpy.linalg as nplin
import numpy.random as npran
import scipy as sp
import scipy.linalg as splin
import scipy.stats as spstat1gt
import matplotlib.pyplot as plt
import glmnet
import knockoffGLM as ko
import simulate as sim

def logit_compare(id,n,b,p,q,p1s,new_seeds=False,procs=4,plot=False):
    # since I'm parrallel computing, need to make sure I'm using different seeeds to generate random numbers
    if new_seeds:
        # these come from random.org, and should be truely random
        # the function is a tad spotty due parsing a web page isntead of using the API, so may need to be run a few times
        # also, it errors out when the request is rejected due to server load
        seeds = sim.getRandomIntegers(num=10000,min=1,max=1000000000)
    else:
        # failing that, backup seeds gotten by good ol copy and paste
        f = open('data/backup_seeds.txt','r')
        seeds = [int(seed) for seed in f.read().split()]

    # Load the genetic data
    genes = np.genfromtxt('data/SNPdata.txt', delimiter=',')
    # replace all the 2s with 1s
    np.place(genes,genes!=0,1)

    isl_data_logit = np.empty((2,len(p1s)))
    gen_data_logit = np.empty((2,len(p1s)))

    for i in range(len(p1s)):
        p1 = p1s[i]
        isl_data_logit[:,i] = sim.ko_test(procs=procs,
                                          seeds=seeds[(i*b):(i*b+b)],
                                          q=q,
                                          knockoff='binary',
                                          model='logit',
                                          func=sim.genXy_binary_X_norm_beta,
                                          n=n,
                                          p1=p1,
                                          pnull=p-p1,
                                          beta_sd=1,
                                          A_base_diag=-1,
                                          A_sd=.2
                                ).calc_fdr_power()
        # save data
        fi = open('data/logit_ising_' +id + '.npy','w')
        np.save(fi,isl_data_logit)
        fi.close()
        print "Done Ising %d" % (i+1)

        gen_data_logit[:,i] = sim.ko_test(procs=procs,
                                          seeds=seeds[(i*b):(i*b+b)],
                                          q=q,
                                          knockoff='binary',
                                          model='logit',
                                          func=sim.genXy_given_X_norm_beta,
                                          n=n,
                                          p1=p1,
                                          pnull=p-p1,
                                          data=genes,
                                          beta_sd=1
                               ).calc_fdr_power()
        print "Done Genetic %d" % (i+1)
        # save data 
        fg = open('data/logit_ising_' +id + '.npy','w')
        np.save(fg,gen_data_logit)
        fg.close()

    # FDR and Power
    f,subplts = plt.subplots(1,2)
    f.set_size_inches(8,4)

    # FDR Plot
    subplts[0].set_xlim(0,max(p1s))
    subplts[0].set_ylim(0,q*2)
    subplts[0].set_xlabel('Sparsity')
    subplts[0].set_ylabel('FDR')
    subplts[0].plot(p1s,isl_data_logit[0,:],label="Ising, Binary", linestyle='--', marker='v')
    subplts[0].plot(p1s,gen_data_logit[0,:],label="Genetic, Binary", linestyle='-', marker='o')
    subplts[0].plot((0,max(p1s)),(q,q),'k--')

    # power plot
    subplts[1].set_xlim(0,max(p1s))
    subplts[1].set_ylim(0,1)
    subplts[1].set_xlabel('Sparsity')
    subplts[1].set_ylabel('Power')
    subplts[1].plot(p1s,isl_data_logit[1,:],label="Ising, Binary", linestyle='--', marker='v')
    subplts[1].plot(p1s,gen_data_logit[1,:],label="Genetic, Binary", linestyle='-', marker='o')

    # Legend
    subplts[1].legend(bbox_to_anchor=(1.05, 0), loc='lower left', borderaxespad=0.)

    plt.savefig('images/logit_FDR_power_' +id +'.pdf',bbox_inches='tight') 


def main():
    logit_compare(n=1000,b=100,p=100,q=.2,p1s=[5,10,15,20,25])

if __name__ == '___main__':
    status = main()
    sys.exit(status)


