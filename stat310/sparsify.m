% STAT 310 Winter 2015
% Programming Assignment, Problem 5
% Version 1.2

rand('state',1);
randn('state',1);
n = 500;
k = 20;
F = 1.5*rand(n,k)-0.5;
Sigma = 0.5*F*diag(rand(k,1))*F' + diag(0.1*rand(n,1));
mu0 = 0.2;
SR = 0.4;
mu = mu0 + SR*sqrt(diag(Sigma));
c = 5 + exp(2*randn(n,1));
