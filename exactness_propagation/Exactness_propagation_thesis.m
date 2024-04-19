% This script describes the propagation of a time dependent TTN, given as
% in the TTN paper of CLW20.

clear all; clc; close all;
global W_l W_tau 

% Initialisations
n_l = [1000,1000,1000,1000,1000,1000]; % make it 1000 for plots
r_l = [10,10,10,10,10,10];
% n_l = [10,10,10,10,10,10];
% r_l = [5 5 5 5 5 5 ];
r_tau = [10,10,r_l(end)];    % if r_l is a leaf, correspondinig r_tau must be equal
T_start = 0;
T_end = 1;
h = [0.1, 0.01, 0.001]; % h = [0.1, 0.01, 0.001];

r_min = r_l(1);
r_max = 20;
d = 6;
tol = 10^-8;

% creates the initial data
[A0,tau] = rand_TTN_complex(n_l,r_l,r_tau);

%% skew symmetric matrices
W_l = cell(1,length(n_l));

for i=1:length(n_l)
    dum = rand(r_l(i),r_l(i));
    dum = dum - dum.';
    N = norm(dum(:));
    W_l{i} = (1/N)*dum;
end

W_tau = cell(1,length(r_tau)-1);
for i=1:(length(r_tau)-1)
    dum = rand(r_tau(i),r_tau(i));
    dum = dum - dum.';
    N = norm(dum(:));
    W_tau{i} = (1/N)*dum;
end

%% TTN integration
Err_ad(:,1) = zeros(1,3);
Err_fix(:,1) = zeros(1,3);
Err_par(:,1) = zeros(1,3);
t(:,1) = zeros(1,3);

for j=1:length(h)
    t_steps = ceil((T_end - T_start)/h(j));
    A1 = A0;
    Y_new_ad = A0;
    Y_new_fix = A0;
    Y_new_par = A0;
    for i=1:t_steps
        
        %% calculate A(t_i) at time t_i
        t0 = (i-1)*h(j);
        t1 = i*h(j);
        A2 = FF(t1,A0);
        
        % Add A(t1) - A(t0) = A2 - A1
        dum = A1;
        dum{end} = -dum{end};
        Delta_A = Add_TTN(A2,dum,tau);
        
        F = @(t,Y,A,d) Delta_A;

        %% rank-adaptive TTN integrator
        Y_new_ad = TTN_integrator_complex_rank_adapt_nonglobal(tau,Y_new_ad,F,0,1,Delta_A,d,r_min);
        Y_new_ad = truncate(Y_new_ad,tol,r_max,r_min);
        A1 = A2;
        % calulate error
        dum = A2;
        dum{end} = -A2{end};
        Diff = Add_TTN(Y_new_ad,dum,tau);
        Diff = Mat0Mat0(Diff,Diff);
        Err_ad(i+1,j) = sqrt(Diff);
        t(i+1,j) = t1;
        
        %% fixed-rank BUG
        Y_new_fix = TTN_integrator_complex_nonglobal(tau,Y_new_fix,F,0,1,Delta_A,d);
        % calulate error
        dum = A2;
        dum{end} = -A2{end};
        Diff = Add_TTN(Y_new_fix,dum,tau);
        Diff = Mat0Mat0(Diff,Diff);
        Err_fix(i+1,j) = sqrt(Diff);
        
        %% parallel TTN integrator
        Y_new_par = TTN_integrator_complex_parallel_nonglobal(tau,Y_new_par,F,0,1,Delta_A,d,r_min);
        Y_new_par = truncate(Y_new_par,tol,r_max,r_min);
        % calulate error
        dum = A2;
        dum{end} = -A2{end};
        Diff = Add_TTN(Y_new_par,dum,tau);
        Diff = Mat0Mat0(Diff,Diff);
        Err_par(i+1,j) = sqrt(Diff);
        
        i
    end    
end
%% plot the error
% rank adaptive BUG
figure(1)
subplot(1,2,1)
marker = {'o-','*-','x-'};
for i=1:length(h)
    steps = 1/h(i) + 1;
    Err2 = Err_ad(1:steps,i);
    tt = t(1:steps,i);
   
    semilogy(tt,abs(Err2),marker{i},'Linewidth',2)
    hold on
end
xlabel('Time','FontSize',20);
ylabel('Error','FontSize',20);
title('Rank-adaptive BUG','FontSize',20)
legend('h=0.1','h=0.01','h=0.001','FontSize',20)
axis([0 1.2 0 7*10^-12])

% fixed-rank BUG
subplot(1,2,2)
marker = {'o-','*-','x-'};
for i=1:length(h)
    steps = 1/h(i) + 1;
    Err2 = Err_fix(1:steps,i);
    tt = t(1:steps,i);
   
    semilogy(tt,abs(Err2),marker{i},'Linewidth',2)
    hold on
end
xlabel('Time','FontSize',20);
ylabel('Error','FontSize',20);
title('Fixed rank BUG','FontSize',20)
legend('h=0.1','h=0.01','h=0.001','FontSize',20)
axis([0 1.2 0 7*10^-12])

% parallel BUG
figure(3)
marker = {'o-','*-','x-'};
for i=1:length(h)
    steps = 1/h(i) + 1;
    Err2 = Err_par(1:steps,i);
    tt = t(1:steps,i);
   
    semilogy(tt,abs(Err2),marker{i},'Linewidth',2)
    hold on
end
xlabel('Time','FontSize',20);
ylabel('Error','FontSize',20);
legend('h=0.1','h=0.01','h=0.001','FontSize',20)
% axis([0 1.2 0 1.2*10^-12])


% % Test
% T1 = double(full_tensor(A2));
% T2 = double(full_tensor(Y_new));
% E = T2 - T1;
% norm(E(:))
% %%%

function [X] = FF(t,Y)
global W_l W_tau
X = Y;

% tau1
dum = Y{1};
for j=1:(length(dum)-2)
    dum{j} = dum{j}*expm(t*W_l{j});
end
dum{end} = ttm(dum{end},expm(t*W_tau{1}),length(dum)-1);
X{1} = dum;

%tau2
dum = Y{2};
for j=1:(length(dum)-2)
    dum{j} = dum{j}*expm(t*W_l{j+3});
end
dum{end} = ttm(dum{end},expm(t*W_tau{2}),length(dum)-1);
X{2} = dum;

%tau3 
X{3} = X{3}*expm(t*W_l{6});

end