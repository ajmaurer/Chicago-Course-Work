#imports
from __future__ import division
import numpy as np
import numpy.linalg as nplin
import numpy.random as npran
import scipy as sp
import scipy.linalg as splin
import scipy.stats as spstat1gt
import matplotlib.pyplot as plt
import json
import inspect
import glmnet
import knockoffGLM as ko
import simulate as sim

def lasso_compare(n,b,p,q,p1s):
    # since I'm parrallel computing, need to make sure I'm using different seeeds to generate random numbers
    # these come from random.org, and should be truely random
    # the function is a tad spotty due parsing a web page isntead of using the API, so may need to be run a few times
    # also, it errors out when the request is rejected due to server load
    seeds = sim.getRandomIntegers(num=10000,min=1,max=1000000000)

    # Load the genetic data
    genes = np.genfromtxt('data/SNPdata.txt', delimiter=',')
    # replace all the 2s with 1s
    np.place(genes,genes!=0,1)

    isl_data_lasso = np.empty((5,len(p1s)))
    gen_data_lasso = np.empty((5,len(p1s)))

    for i in range(len(p1s)):
        p1 = p1s[i]
        isl_data_lasso[:,i] = sim.ko_test(procs=4,
                                          seeds=seeds[(i*b):(i*b+b)],
                                          q=q,
                                          knockoff='both',
                                          model='lasso',
                                          func=sim.genXy_binary_X_norm_y,
                                          n=n,
                                          p1=p1,
                                          pnull=p-p1,
                                          beta_sd=1,
                                          A_base_diag=-1,
                                          A_sd=.2
                                ).calc_fdr_power()
        print "Done Ising %d" % (i+1)
        gen_data_lasso[:,i] = sim.ko_test(procs=4,
                                          seeds=seeds[(i*b):(i*b+b)],
                                          q=q,
                                          knockoff='both',
                                          model='lasso',
                                          func=sim.genXy_given_X_norm_y,
                                          n=n,
                                          p1=p1,
                                          pnull=p-p1,
                                          data=genes,
                                          beta_sd=1
                               ).calc_fdr_power()
        print "Done Genetic %d" % (i+1)

    # save data 
    fg = open('data/lasso_genetic.json','w')
    json.dump(gen_data_lasso,fg)
    fg.close()

    fi = open('data/lasso_ising.json','w')
    json.dump(isl_data_lasso,fi)
    fi.close()

    # FDR and Power
    f,subplts = plt.subplots(1,2)
    f.set_size_inches(8,4)

    # FDR Plot
    subplts[0].set_xlim(0,max(p1s))
    subplts[0].set_ylim(0,q*2)
    subplts[0].set_xlabel('Sparsity')
    subplts[0].set_ylabel('FDR')
    subplts[0].plot(p1s,isl_data_lasso[0,:],label="Ising, Binary", linestyle='--', marker='v', color='g')
    subplts[0].plot(p1s,isl_data_lasso[1,:],label="Ising, Original", linestyle='--', marker='v', color='r')
    subplts[0].plot(p1s,gen_data_lasso[0,:],label="Genetic, Binary", linestyle='-', marker='o', color='g')
    subplts[0].plot(p1s,gen_data_lasso[1,:],label="Genetic, Original", linestyle='-', marker='o', color='r')
    subplts[0].plot((0,max(p1s)),(q,q),'k--')

    # power plot
    subplts[1].set_xlim(0,max(p1s))
    subplts[1].set_ylim(0,1)
    subplts[1].set_xlabel('Sparsity')
    subplts[1].set_ylabel('Power')
    subplts[1].plot(p1s,isl_data_lasso[2,:],label="Ising, Binary", linestyle='--', marker='v', color='g')
    subplts[1].plot(p1s,isl_data_lasso[3,:],label="Ising, Original", linestyle='--', marker='v', color='r')
    subplts[1].plot(p1s,gen_data_lasso[2,:],label="Genetic, Binary", linestyle='-', marker='o', color='g')
    subplts[1].plot(p1s,gen_data_lasso[3,:],label="Genetic, Original", linestyle='-', marker='o', color='r')

    # Legend
    subplts[1].legend(bbox_to_anchor=(1.05, 0), loc='lower left', borderaxespad=0.)

    plt.savefig('images/lasso_FDR_power.pdf',bbox_inches='tight')      
     
    # Correlation plot  
    f,ax = plt.subplots(1)
    f.set_size_inches(4,4)

    ax.set_xlim(0,max(p1s))
    ax.set_ylim(0,1)
    ax.set_xlabel('Sparsity')
    ax.set_ylabel('Correlation')
    ax.plot(p1s,isl_data_lasso[4,:],linestyle='--', marker='v',label="Ising")
    ax.plot(p1s,gen_data_lasso[4,:],linestyle='-', marker='o',label="Genetic")
    ax.legend(bbox_to_anchor=(1.05, 0), loc='lower left', borderaxespad=0.)

    subplts[1].legend(bbox_to_anchor=(1.05, 0), loc='lower left', borderaxespad=0.)

    plt.savefig('images/lasso_corr.pdf',bbox_inches='tight')  


def main():
    lasso_compare(n=1000,b=100,p=100,q=.2,p1s=[5,10,15,20,25])

if __name__ == '___main__':
    status = main()
    sys.exit(status)


