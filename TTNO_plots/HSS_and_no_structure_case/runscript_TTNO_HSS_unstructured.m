clear all; clc; close all;

addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\TTNO')
addpath('C:\Users\Dominik\Documents\MATLAB\Matlab toolboxes\hm-toolbox-master\hm-toolbox-master')
% addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\rank_adaptive_integrator_for_TTN')
addpath('C:\Users\Dominik\Documents\MATLAB\Matlab toolboxes\tensor_toolbox-master')

%% initializations
d = 2^8;           % number of particles
l = log(d)/log(2); % number of layers
n = 2;             % physical dimension

[X,tau] = init_spin_all_dim_same_rank(4,1,d); % initial data - needed for construct exact operator

A = cell(1,d);
% tmp = rand(n,n);
for jj=1:d
    A{jj} = rand(n,n);     % matrix acting on each site
%     A{jj} = tmp;
%     A{jj} = ones(n,n);
end

% interaction matrix (has hss_rank 1)
tmp = ones(d,1);
A_int = diag(2*tmp,0) - diag(tmp(1:d-1),-1) - diag(tmp(1:d-1),1);
V = inv(A_int);

% % interaction matrix 2
% k = 4;
% V = zeros(d,d);
% for ii=1:k
%     v = rand(d,1);
%     w = rand(d,1);
%     V = V + v*w';
% end
% % V = 0.5*(V + V.'); % symmetrischer Fall

%% HSS 
cl = 1:1:d;
H = hss(V,'cluster',cl);

%% construction of TTNO with help of HSS
tic
TTNO = TTNO_HSS_construct_unstructured(H,A,l,l,n*ones(1,d),1:d);
toc

%% exact TTNO without using HSS structure
tic
TTNO_exact = TTNO_no_structure(A,V,l,l,n*ones(d,1),1:d);
toc

% error check
tmp = TTNO;
tmp{end} = -tmp{end};
E = Add_TTN(TTNO_exact,tmp,tau);
err = sqrt(abs(Mat0Mat0(E,E)))

err_scaled = sqrt(abs(Mat0Mat0(E,E)))/sqrt(abs(Mat0Mat0(TTNO_exact,TTNO_exact)))


% % exact TTNO by old construction
% B = cell(1,d);
% count = 1;
% for ii=1:d
%     for jj=1:d
%         if ii < jj
%             B{count,ii} = V(ii,jj)*A{ii};
%             B{count,jj} = A{jj};
%             count = count + 1;
% %             V(ii,jj)
%         end
%     end
% end
% tic
% TTNO_exact = make_operator(X,B,tau,n*ones(d,1));
% TTNO_exact = rounding(TTNO_exact,tau);
% toc