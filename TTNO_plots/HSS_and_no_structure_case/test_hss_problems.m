clear all; clc; close all;

addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\TTNO')
addpath('C:\Users\Dominik\Documents\MATLAB\Matlab toolboxes\hm-toolbox-master\hm-toolbox-master')
addpath('C:\Users\Dominik\Documents\MATLAB\Matlab toolboxes\tensor_toolbox-master')

%% initializations
d = 2^3;           % number of particles
l = log(d)/log(2); % number of layers
n = 2;             % physical dimension

[X,tau] = init_spin_all_dim_same_rank(4,1,d); % initial data - needed for construct exact operator

A = cell(1,d);
for jj=1:d
    A{jj} = rand(n,n);     % matrix acting on each site
end

% interaction matrix 
tmp = ones(d,1);
A_int = diag(2*tmp,0) - diag(tmp(1:d-1),-1) - diag(tmp(1:d-1),1);
% A_int2 = diag(2*tmp,0) - diag(4*tmp(1:d-1),-1) - diag(3*tmp(1:d-1),1);
% % V = inv(A_int) + inv(A_int2);


% V = inv(A_int);
V = A_int;

% V = 0.5*(V + V.'); % symmetrischer Fall


%% HSS 
% first idea
H = hss(V,'cluster',1:d);
% H = proper(H);
% H = compress(H, 10^-15);
k = hssrank(H);


%% construction of TTNO with help of HSS
TTNO = TTNO_HSS_construct(H,A,l,l,n*ones(1,d),1:d,k);

% %% exact TTNO without using HSS structure (reference solution)
% TTNO_exact = TTNO_no_structure(A,V,l,l,n*ones(d,1),1:d);
% 
% %% error check
% tmp = TTNO;
% tmp{end} = -tmp{end};
% E = Add_TTN(TTNO_exact,tmp,tau);
% err = sqrt(abs(Mat0Mat0(E,E)))
% 
% err_scaled = sqrt(abs(Mat0Mat0(E,E)))/sqrt(abs(Mat0Mat0(TTNO_exact,TTNO_exact)))

