clear all; clc; close all;

addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\TTNO')
addpath('C:\Users\Dominik\Documents\MATLAB\Matlab toolboxes\hm-toolbox-master')
% addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\rank_adaptive_integrator_for_TTN')
% addpath('C:\Users\Dominik\Documents\MATLAB\Matlab toolboxes\tensor_toolbox-master')

%% initializations
d = 2^4;           % number of particles
l = log(d)/log(2); % number of layers
n = 2;             % physical dimension

% tree structures
% [X,tau] = init_spin_all_dim_same_rank(4,1,d); % binary tree - needed for construct exact operator
% [X,tau] = init_spin_all_dim_diff_rank_tree4(d); % tree with 4 subtrees (d=16)
% [X,tau] = init_spin_all_dim_diff_rank_tree4_d32(d); % tree with 4 subtrees (d=32)
[X,tau] = init_spin_all_dim_diff_rank_TT(d); % TT 


A = cell(1,d);
% tmp = rand(n,n);
for jj=1:d
    A{jj} = rand(n,n);     % matrix acting on each site
%     A{jj} = tmp;
end

% interaction matrix
tmp = ones(d,1);
A_int = diag(2*tmp,0) - diag(tmp(1:d-1),-1) - diag(tmp(1:d-1),1);
V = inv(A_int);

% % interaction matrix 2
% k = 5;
% V = zeros(d,d);
% for ii=1:k
%     v = rand(d,1);
%     w = rand(d,1);
%     V = V + v*w';
% end
% % V = 0.5*(V + V.'); % symmetrischer Fall

%% HSS 
% cl = 1:1:d;
% H = hss(V,'cluster',cl);

%% construction of TTNO with help of HSS
tic
TTNO = TTNO_no_structure_arbitrary_tree(A,V,X,l,l,n*ones(d,1),1:d);
toc 
% TTNO = rounding(TTNO,tau);

%% exact TTNO
B = cell(1,d);
count = 1;
for ii=1:d
    for jj=1:d
        if ii < jj
            B{count,ii} = V(ii,jj)*A{ii};
            B{count,jj} = A{jj};
            count = count + 1;
%             V(ii,jj)
        end
    end
end
tic
TTNO_exact = make_operator(X,B,tau,n*ones(d,1));
TTNO_exact = rounding(TTNO_exact,tau);
toc

% error check
tmp = TTNO;
tmp{end} = -tmp{end};
E = Add_TTN(TTNO_exact,tmp,tau);
err = sqrt(abs(Mat0Mat0(E,E)))