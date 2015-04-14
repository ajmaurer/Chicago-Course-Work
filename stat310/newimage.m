% STAT 310 Winter 2015
% Programming Assignment, Problem 3
% Version 1.1

X = double(imread('newimage.png'));
[m, n] = size(X);
rand('state', 1029);
S = rand(m,n) > 0.5;

%imagesc(X);
%colormap gray;