clc; clear all;
% This script solves a ODE of a low rank approximation for tensors with
% a constant F(Y(t))= A2-A1, Y(t_0)= A1
% where A1,A2 are TTN's, via the new integrator
clear; clc;
addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\unconventional_TTN_integrator')
addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\rank_adaptive_integrator_for_TTN')
addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\parallel_TTN_integrator')


%global d a b K lamda r q r_tau tau A
% global d

% profile on

n_l = 10*ones(1,6);
r_l = 5*ones(1,6);
r_tau = [5,5,r_l(end)];    % if r_l is a leaf, correspondinig r_tau must be equal
t0 = 0;
T_end = 1;
h = 1;
d = 6;
r_max = 10;
r_min = 2;

% creates the initial data
[A1,tau] = rand_TTN_complex(n_l,r_l,r_tau);
[A2,~] = rand_TTN_complex(n_l,r_l,r_tau);

% d = 3;
% [A1,tau] = rand_Tucker_complex();
% [A2,~] = rand_Tucker_complex();

% test 2
% d = 9; % DOF, i.e. number of particles
% K = 31*ones(1,d); % number of basis functions in the Fourier base
% r = 4*ones(1,d); % rank of the system
% r_tau = 4;
% q = 2*ones(1,d); % positions of each 1-d Gaussain
% [A1,tau] = init_Gaussian_binary_tree_all_dim();
% [A2,~] = init_Gaussian_binary_tree_all_dim();


% Computes Delta_A = A2-A1
dum = A1;
dum{end} = -dum{end};
Delta_A = Add_TTN(A2,dum,tau);
clear dum;
F_Delta_A = @(t,Y) Delta_A;

%% TTN integration
t_steps = ceil((T_end - t0)/h);
tol = 10^-4;

for i=1:t_steps
%     tic
%     Y_new = TTN_integrator_complex(tau,A1,F_Delta_A,(i-1)*h,i*h);
%     toc
    
    Y_new = TTN_integrator_complex_rank_adapt(tau,A1,F_Delta_A,(i-1)*h,i*h);
    tmp = Y_new;
    Y_new = truncate(Y_new,tol,r_max,r_min);
    
    % rank-adaptive TTN integrator
    Y_new_ad = TTN_integrator_complex_rank_adapt_nonglobal(tau,Y0,F_tau,t0,t1,A,d,r_min);
    Y_new_ad = truncate_old(Y_new_ad,tol,r_max);
    

    
%     %%% Difference between augmented and truncated TTN
%     E = Y_new;
%     E{end} = -E{end};
%     F = Add_TTN(tmp,E,tau);
%     err = sqrt(Mat0Mat0(F,F))

%     [Y_new,~] = TTN_integrator_complex_rank_adapt_trunc2(tau,A1,F_Delta_A,(i-1)*h,i*h,tol);
    Y_start = Y_new;
end


%% Error plots
dum = A2;
dum{end} = -dum{end};
Diff = Add_TTN(dum,Y_new,tau);
Mat0 = Mat0Mat0(Diff,Diff);
Err = abs(sqrt(Mat0Mat0(Diff,Diff)))


% profile off
% profile viewer
profile clear
