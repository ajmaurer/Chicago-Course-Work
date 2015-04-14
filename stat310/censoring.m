% STAT 310 Winter 2015
% Programming Assignment, Problem 4
% Version 1.1

randn('state',0);
d = 20;
k = 25;
n = 100;
x_true = randn(d,1);
A = randn(d,n);
b = A'*x_true + 0.1*(sqrt(d))*randn(n,1);

[b, sort_ind] = sort(b);
A = A(:,sort_ind);
beta = (b(k)+b(k+1))/2;
b = b(1:k);