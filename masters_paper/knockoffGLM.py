""" This module is meant to implement and test the knockoff procedure on L1 regularized GLM models. The method was originally designed for linear least squares regression, and there is no theoretical guarantee it will work on other GLMS. Hopefully, with simulation we can discover if it breaks, if so where it breaks, and maybe how to fix it in those situations.

More can be read about the original procedure on least squares linear models here: http://web.stanford.edu/~candes/Knockoffs/index.html

"""

# imports 
import numpy as np
import matplotlib.pyplot as plt
import cvxpy as cvx
from sklearn import linear_model as lm
from sklearn.svm import l1_min_c 


if __name__ == '___main__':
    status = main()
    sys.exit(status)

