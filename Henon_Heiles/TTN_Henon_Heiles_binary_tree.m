% This script computes an approximation to the solution of the Schrödinger
% equation, with H = \Delta + Henon-Heiles potential.

clear; clc; close all;

addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\unconventional_TTN_integrator')
addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\rank_adaptive_integrator_for_TTN')
addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\parallel_TTN_integrator')
addpath(genpath('tensor_toolbox-master'))
addpath(genpath('tensorlab_2016-03-28'))
addpath(genpath('hm-toolbox-master'))

% global a b d K lamda r q r_tau tau A D M1 M2 M3 W
% profile on

% Initialisations
a = -10;
b = 10;
lamda = 0.111803;

d = 4;            % DOF, i.e. number of particles
K = 31*ones(1,d); % number of basis functions in the Fourier base
q = 2*ones(1,d);  % positions of each 1-d Gaussain

tol = 10^-8;

t0 = 0;
T_end = 1; 
dt = 0.005;
r_max = 10;
r_min = 4;

r = 10*ones(1,d);  % rank of the system at the leaves
r_tau = r_max;        % ranks at the core tensors

% initial data
[X0,tau] = init_Gaussian_binary_tree_all_dim(K,r,d,a,b,q,r_tau);

%% TTN integration
t_steps = ceil((T_end - t0)/dt);
% pre_calculations_exp_composed(a,b,K,d);
[D,M1,M2,M3,W] = pre_calculations_exp(a,b,K,d);
% pre_calculations_exp_diff_CAP(a,b,K,d);

% represent the Hamiltonian as a linear operator
B = linearisation(D,M1,M2,M3,W,d,lamda);
A = make_operator(X0,B,tau,K);
A{end} = -1i*A{end};

% test to make the operator with my funciton
% B = linearisation_test();
% A = make_operator(Y0,B,tau,K);

ac = zeros(t_steps+1,1);
Norm_con = zeros(t_steps+1,1);
Energy_con = zeros(t_steps+1,1);
ac_ad = zeros(t_steps+1,1);
Norm_con_ad = zeros(t_steps+1,1);
Energy_con_ad = zeros(t_steps+1,1);
ac_BUG = zeros(t_steps+1,1);
Norm_con_BUG = zeros(t_steps+1,1);
Energy_con_BUG = zeros(t_steps+1,1);

ac(1) = Mat0Mat0_prod(X0,X0,a,b,d);
Norm_con(1) = sqrt(ac(1));
F = F_Hamiltonian(t0,X0,A,d);
Energy_con(1) = 1i*Mat0Mat0_prod(X0,F,a,b,d); % 1i* as in F_Ham I have 1/i *H

ac_ad(1) = ac(1);
ac_BUG(1) = ac(1);
Norm_con_BUG(1) = Norm_con(1);
Norm_con_ad(1) = Norm_con(1);
Energy_con_BUG(1) = Energy_con(1);
Energy_con_ad(1) = Energy_con(1);

X_start = X0;
X_start_ad = X0;
X_start_BUG = X0;

for it=1:t_steps
    t0 = (it-1)*dt;
    t1 = it*dt;
    time(it) = t1;
    
    % rank-adaptive BUG
    X_new_ad = TTN_integrator_complex_rank_adapt_nonglobal(tau,X_start_ad,@F_Hamiltonian,t0,t1,A,d,r_min);
    X_new_ad = truncate(X_new_ad,tol,r_max,r_min);
    
    % fixed-rank BUG
    X_new_BUG = TTN_integrator_complex_nonglobal_spin(tau,X_start_BUG,@F_Hamiltonian,t0,t1,A,d);
    
    % parallel BUG
    [X_new,~,~] = TTN_integrator_complex_parallel_nonglobal(tau,X_start,@F_Hamiltonian,t0,t1,A,d,r_min);
    X_new = truncate(X_new,tol,r_max,r_min);
        
    % norm and energy
    Norm_con_ad(it+1) = sqrt(Mat0Mat0_prod(X_new_ad,X_new_ad,a,b,d));
    F = F_Hamiltonian(t1,X_new_ad,A,d);
    Energy_con_ad(it+1) = 1i*Mat0Mat0_prod(X_new_ad,F,a,b,d); % 1i* as in F_Ham I have 1/i *H
    
    Norm_con_BUG(it+1) = sqrt(Mat0Mat0_prod(X_new_BUG,X_new_BUG,a,b,d));
    F = F_Hamiltonian(t1,X_new_BUG,A,d);
    Energy_con_BUG(it+1) = 1i*Mat0Mat0_prod(X_new_BUG,F,a,b,d); % 1i* as in F_Ham I have 1/i *H
    
    Norm_con(it+1) = sqrt(Mat0Mat0_prod(X_new,X_new,a,b,d));
    F = F_Hamiltonian(t1,X_new,A,d);
    Energy_con(it+1) = 1i*Mat0Mat0_prod(X_new,F,a,b,d); % 1i* as in F_Ham I have 1/i *H
    
    % compute overlap with initial data
    ac_ad(it+1) = Mat0Mat0_prod(X_new_ad,X0,a,b,d);
    ac_BUG(it+1) = Mat0Mat0_prod(X_new_BUG,X0,a,b,d);
    ac(it+1) = Mat0Mat0_prod(X_new,X0,a,b,d);
    
    % set for next time step
    X_start_ad = X_new_ad;
    X_start_BUG = X_new_BUG;
    X_start = X_new;

    fprintf('Rank-adaptive t=%f, E=%f, norm=%f \n', it*dt,  abs(Energy_con_ad(it+1,1)), abs(Norm_con_ad(it+1,1)))
    fprintf('Fixed-rank t=%f, E=%f, norm=%f \n', it*dt,  abs(Energy_con_BUG(it+1,1)), abs(Norm_con_BUG(it+1,1)))
    fprintf('Parallel t=%f, E=%f, norm=%f \n', it*dt,  abs(Energy_con(it+1,1)), abs(Norm_con(it+1,1)))

end


%% compute vibrational spectrum
omega = 0:0.001:30; 

% compute the integral of the FT via the trapezoidal rule rank-adaptive BUG
FT3_ra = zeros(length(omega),1);
time = 0:dt:T_end;
for i=1:length(omega)
    for k=1:length(ac)
        if k==1 || k==length(ac)
            FT3_ra(i) = FT3_ra(i) + 0.5*dt*ac(k)*exp(-1i*omega(i)*time(k));
        else
            FT3_ra(i) = FT3_ra(i) + dt*ac(k)*exp(-1i*omega(i)*time(k));
        end
    end
end
FT3_ra = (1/pi)*real(FT3_ra);

% compute the integral of the FT via the trapezoidal rule fixed-rank BUG
FT3_BUG = zeros(length(omega),1);
time = 0:dt:T_end;
for i=1:length(omega)
    for k=1:length(ac)
        if k==1 || k==length(ac)
            FT3_BUG(i) = FT3_BUG(i) + 0.5*dt*ac(k)*exp(-1i*omega(i)*time(k));
        else
            FT3_BUG(i) = FT3_BUG(i) + dt*ac(k)*exp(-1i*omega(i)*time(k));
        end
    end
end
FT3_BUG = (1/pi)*real(FT3_BUG);

% compute the integral of the FT via the trapezoidal rule parallel BUG
FT3 = zeros(length(omega),1);
time = 0:dt:T_end;
for i=1:length(omega)
    for k=1:length(ac)
        if k==1 || k==length(ac)
            FT3(i) = FT3(i) + 0.5*dt*ac(k)*exp(-1i*omega(i)*time(k));
        else
            FT3(i) = FT3(i) + dt*ac(k)*exp(-1i*omega(i)*time(k));
        end
    end
end
FT3 = (1/pi)*real(FT3);


% Plot the results
figure(1)
subplot(1,3,1); plot(omega,abs(FT3_ra))
legend('$S(\omega )$','Interpreter','latex','FontSize',20,'Linewidth',2);
title('Rank-adaptive BUG')
subplot(1,3,2); plot(omega,abs(FT3_BUG))
legend('$S(\omega )$','Interpreter','latex','FontSize',20,'Linewidth',2)
title('Fixed-rank BUG')
subplot(1,3,3); plot(omega,abs(FT3))
legend('$S(\omega )$','Interpreter','latex','FontSize',20,'Linewidth',2)
title('Parallel BUG')

figure(2)
subplot(2,3,1); plot(0:dt:T_end,abs(Norm_con_ad))
title('Rank-adaptive BUG')
legend('$ \vert \vert X(t) \vert \vert$','Interpreter','latex','FontSize',14)
subplot(2,3,4); plot(0:dt:T_end,abs(Energy_con_ad))
legend('$ < X(t) \vert H \vert X(t) > $','Interpreter','latex','FontSize',14)

subplot(2,3,2); plot(0:dt:T_end,abs(Norm_con_BUG))
title('Fixed-rank BUG')
legend('$ \vert \vert X(t) \vert \vert$','Interpreter','latex','FontSize',14)
subplot(2,3,5); plot(0:dt:T_end,abs(Energy_con_BUG))
legend('$ < X(t) \vert H \vert X(t) > $','Interpreter','latex','FontSize',14)

subplot(2,3,3); plot(0:dt:T_end,abs(Norm_con))
title('Parallel BUG')
legend('$ \vert \vert X(t) \vert \vert$','Interpreter','latex','FontSize',14)
subplot(2,3,6); plot(0:dt:T_end,abs(Energy_con))
legend('$ < X(t) \vert H \vert X(t) > $','Interpreter','latex','FontSize',14)

% profile off
% profile viewer