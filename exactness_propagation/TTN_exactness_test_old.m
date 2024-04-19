% This script solves a ODE of a low rank approximation for tensors with
% a constant F(Y(t))= A2-A1, Y(t_0)= A1
% where A1,A2 are TTN's, via the new integrator
clear; clc;

% profile on

% Initialisations
% n_l = [10,11,12,13,14,15];
% r_l = [2,3,4,5,6,7];
% r_tau = [8,9,7];    % if r_l is a leaf, r_tau must be equal
n_l = [10,10,10,10,10,10];
r_l = [5,5,5,5,5,5];
r_tau = [5,5,r_l(end)];    % if r_l is a leaf, correspondinig r_tau must be equal
t0 = 0;
T_end = 1;
h = 1;

% creates the initial data
[A1,tau] = rand_TTN(n_l,r_l,r_tau);
[A2,~] = rand_TTN(n_l,r_l,r_tau);

% Computes Delta_A = A2-A1
dum = A1;
dum{end} = -dum{end};
Delta_A = Add_TTN(A2,dum,tau);
clear dum;
F_Delta_A = @(t,Y) Delta_A;

%% TTN integration
t_steps = ceil((T_end - t0)/h);

for i=1:t_steps
    Y_new = TTN_integrator(tau,A1,F_Delta_A,(i-1)*h,i*h);
    Y_start = Y_new;
end


%% Error plots
dum = A2;
dum{end} = -dum{end};
Diff = Add_TTN(dum,Y_new,tau);
Mat0 = Mat0Mat0(Diff,Diff)
Err = sqrt(Mat0Mat0(Diff,Diff))


% profile off
% profile viewer
profile clear
