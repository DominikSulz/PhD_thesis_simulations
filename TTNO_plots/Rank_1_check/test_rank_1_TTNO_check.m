% Check if the construction is right

clear all; clc; close all;

addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\TTNO')

% physical dimension
n = 10;
tmp = ones(n, 1);
A = diag(2*tmp,0) - diag(tmp(1:n-1),-1) - diag(tmp(1:n-1),1);
I = eye(n,n);

% problem
d = 8;
v = 0:(1/(d-1)):1;

% tree structure
[X,tau] = init_spin_all_dim_same_rank(4,1,d); % binary tree
% [X,tau] = init_spin_all_dim_diff_rank_TT(d); % tensor train
% [X,tau] = init_spin_all_dim_diff_rank_tree4(d); % tree with 4 subtrees (only for d=16)

% linearisation
B = linearisation_vv(A,d,v);

% exact operator
tic
TTNO_exact = make_operator(X,B,tau,ones(d,1)*n);
TTNO_exact2 = rounding(TTNO_exact,tau);
time_full_construction = toc
% TTNO_exact = truncate(TTNO_exact,10^-12,100,2);

% rank_1_operator
tic
TTNO = rank1_TTNO_vv(X,v,A,I,d);
time_rank1_construction = toc

% check
tmp = TTNO;
tmp{end} = -tmp{end};
E = Add_TTN(TTNO_exact,tmp,tau);
err = sqrt(abs(Mat0Mat0(E,E)))