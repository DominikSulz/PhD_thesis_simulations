% This script computes an approximation to the solution of the Schrödinger
% equation, with H = \Delta + Henon-Heiles potential using a Strang
% splitting. The Laplacian part is then treated exactly, while the
% potential is evolved with the BUG integrators

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
lambda = 0.111803;

d = 2;            % DOF, i.e. number of particles
K = 31*ones(1,d); % number of basis functions in the Fourier base
q = 2*ones(1,d);  % positions of each 1-d Gaussain

tol = 10^-8;

t0 = 0;
T_end = 60; 
dt = 0.01;
r_max = 4;
r_min = 4;

r = 4*ones(1,d);  % rank of the system at the leaves
r_tau = r_max;        % ranks at the core tensors

% initial data
[X0,tau] = init_Gaussian_binary_tree_all_dim(K,r,d,a,b,q,r_tau);

%% TTN integration
t_steps = ceil((T_end - t0)/dt);
% pre_calculations_exp_composed(a,b,K,d);
[D,M1,M2,M3,W] = pre_calculations_exp(a,b,K,d);
% pre_calculations_exp_diff_CAP(a,b,K,d);

% represent the full Hamiltonian as a linear operator
% B = linearisation(D,M1,M2,M3,W,d,lamda);
% A = make_operator(Y0,B,tau,K);
% A{end} = (1/1i)*A{end};

% % represent the potential as a linear operator, without the Laplacian
% DD = D;
% for ii=1:d
%     DD{ii} = 0*D{ii};
% end
% B = linearisation(DD,M1,M2,M3,W,d,lambda);
% A = make_operator(X0,B,tau,K);
% A = rounding(A,tau);
% A = truncate(A,10^-14,50,2);
% A{end} = -1i*A{end};
% A1 = A;
% 
% B = linearisation(D,M1,M2,M3,W,d,lambda);
% A_full = make_operator(X0,B,tau,K);
% A_full = rounding(A_full,tau);
% A_full = truncate(A_full,10^-14,50,2);
% A_full{end} = -1i*A_full{end};
% A_full1 = A_full;

%% make operators with Kressner method
DD = D;
for ii=1:d
    DD{ii} = 0*D{ii};
end
B = linearisation_kressner(DD,M1,M2,M3,W,d,lambda);
A = make_operator_kressner(X0,B,d);
A{end} = -1i*A{end};

B = linearisation_kressner(D,M1,M2,M3,W,d,lambda);
A_full = make_operator_kressner(X0,B,d);
A_full{end} = -1i*A_full{end};

%%

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
F = F_Hamiltonian(t0,X0,A_full,d);
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

rank_max = zeros(t_steps+1,1);
rank_max(1) = max_rank(X0);
rank_max_ad = zeros(t_steps+1,1);
rank_max_ad(1) = max_rank(X0);
rank_max_BUG = zeros(t_steps+1,1);
rank_max_BUG(1) = max_rank(X0);

for it=1:t_steps
    t0 = (it-1)*dt;
    t1 = it*dt;
    t_half = (t1 + t0)/2;
    time(it) = t1;
    
    % rank-adaptive BUG
    X_new_ad = TTN_integrator_complex_rank_adapt_nonglobal(tau,X_start_ad,@F_Hamiltonian,t0,t_half,A,d,r_min);
    X_new_ad = truncate(X_new_ad,tol,r_max,r_min); % half-step with potential
    X_new_ad = exact_Laplace(X_new_ad,dt,D,d);     % full step with the Laplacian
    X_new_ad = TTN_integrator_complex_rank_adapt_nonglobal(tau,X_new_ad,@F_Hamiltonian,t_half,t1,A,d,r_min);
    X_new_ad = truncate(X_new_ad,tol,r_max,r_min); % half-step with potential
    
    % fixed-rank BUG
    X_new_BUG = TTN_integrator_complex_nonglobal(tau,X_start_BUG,@F_Hamiltonian,t0,t_half,A,d); % half-step with potential
    X_new_BUG = exact_Laplace(X_new_BUG,dt,D,d);     % full step with the Laplacian
    X_new_BUG = TTN_integrator_complex_nonglobal(tau,X_new_BUG,@F_Hamiltonian,t_half,t1,A,d); % half-step with potential
    
    % parallel BUG
    [X_new,~,~] = TTN_integrator_complex_parallel_nonglobal(tau,X_start,@F_Hamiltonian,t0,t_half,A,d,r_min);
    X_new = truncate(X_new,tol,r_max,r_min); % half-step with potential
    X_new = exact_Laplace(X_new,dt,D,d);     % full step with the Laplacian
    [X_new,~,~] = TTN_integrator_complex_parallel_nonglobal(tau,X_new,@F_Hamiltonian,t_half,t1,A,d,r_min);
    X_new = truncate(X_new,tol,r_max,r_min); % half-step with potential
        
    % norm and energy
    Norm_con_ad(it+1) = sqrt(Mat0Mat0_prod(X_new_ad,X_new_ad,a,b,d));
    F = F_Hamiltonian(t1,X_new_ad,A_full,d);
    Energy_con_ad(it+1) = 1i*Mat0Mat0_prod(X_new_ad,F,a,b,d); % 1i* as in F_Ham I have 1/i *H
    
    Norm_con_BUG(it+1) = sqrt(Mat0Mat0_prod(X_new_BUG,X_new_BUG,a,b,d));
    F = F_Hamiltonian(t1,X_new_BUG,A_full,d);
    Energy_con_BUG(it+1) = 1i*Mat0Mat0_prod(X_new_BUG,F,a,b,d); % 1i* as in F_Ham I have 1/i *H
    
    Norm_con(it+1) = sqrt(Mat0Mat0_prod(X_new,X_new,a,b,d));
    F = F_Hamiltonian(t1,X_new,A_full,d);
    Energy_con(it+1) = 1i*Mat0Mat0_prod(X_new,F,a,b,d); % 1i* as in F_Ham I have 1/i *H
    
    % compute overlap with initial data
    ac_ad(it+1) = Mat0Mat0_prod(X_new_ad,X0,a,b,d);
    ac_BUG(it+1) = Mat0Mat0_prod(X_new_BUG,X0,a,b,d);
    ac(it+1) = Mat0Mat0_prod(X_new,X0,a,b,d);
    
    % compute max. ranks 
    rank_max_ad(it+1) = max_rank(X_new_ad);
    rank_max_BUG(it+1) = max_rank(X_new_BUG);
    rank_max(it+1) = max_rank(X_new);
    
    % set for next time step
    X_start_ad = X_new_ad;
    X_start_BUG = X_new_BUG;
    X_start = X_new;

    fprintf('Rank-adaptive t=%f, E=%f, norm=%f \n', it*dt,  abs(Energy_con_ad(it+1,1)), abs(Norm_con_ad(it+1,1)))
    fprintf('Fixed-rank t=%f, E=%f, norm=%f \n', it*dt,  abs(Energy_con_BUG(it+1,1)), abs(Norm_con_BUG(it+1,1)))
    fprintf('Parallel t=%f, E=%f, norm=%f \n', it*dt,  abs(Energy_con(it+1,1)), abs(Norm_con(it+1,1)))

end


%% compute vibrational spectrum
omega = 0:0.001:20; 

% compute the integral of the FT via the trapezoidal rule rank-adaptive BUG
FT3_ra = zeros(length(omega),1);
time = 0:dt:T_end;
for i=1:length(omega)
    for k=1:length(ac_ad)
        if k==1 || k==length(ac_ad)
            FT3_ra(i) = FT3_ra(i) + 0.5*dt*ac_ad(k)*exp(-1i*omega(i)*time(k));
        else
            FT3_ra(i) = FT3_ra(i) + dt*ac_ad(k)*exp(-1i*omega(i)*time(k));
        end
    end
end
FT3_ra = (1/pi)*real(FT3_ra);

% compute the integral of the FT via the trapezoidal rule fixed-rank BUG
FT3_BUG = zeros(length(omega),1);
time = 0:dt:T_end;
for i=1:length(omega)
    for k=1:length(ac_BUG)
        if k==1 || k==length(ac_BUG)
            FT3_BUG(i) = FT3_BUG(i) + 0.5*dt*ac_BUG(k)*exp(-1i*omega(i)*time(k));
        else
            FT3_BUG(i) = FT3_BUG(i) + dt*ac_BUG(k)*exp(-1i*omega(i)*time(k));
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
subplot(1,3,1); plot(omega,abs(FT3_ra),'Linewidth',1)
legend('$S(\omega )$','Interpreter','latex','FontSize',20,'Linewidth',2);
title('Rank-adaptive BUG')
subplot(1,3,2); plot(omega,abs(FT3_BUG),'Linewidth',1)
legend('$S(\omega )$','Interpreter','latex','FontSize',20,'Linewidth',2)
title('Fixed-rank BUG')
subplot(1,3,3); plot(omega,abs(FT3),'Linewidth',1)
legend('$S(\omega )$','Interpreter','latex','FontSize',20,'Linewidth',2)
title('Parallel BUG')

figure(2)
subplot(2,3,1); plot(0:dt:T_end,abs(Norm_con_ad),'Linewidth',2)
title('Rank-adaptive BUG')
legend('$ \vert \vert X(t) \vert \vert$','Interpreter','latex','FontSize',14)
subplot(2,3,4); plot(0:dt:T_end,abs(Energy_con_ad),'Linewidth',2)
legend('$ < X(t) \vert H \vert X(t) > $','Interpreter','latex','FontSize',14)

subplot(2,3,2); plot(0:dt:T_end,abs(Norm_con_BUG),'Linewidth',2)
title('Fixed-rank BUG')
legend('$ \vert \vert X(t) \vert \vert$','Interpreter','latex','FontSize',14)
subplot(2,3,5); plot(0:dt:T_end,abs(Energy_con_BUG),'Linewidth',2)
legend('$ < X(t) \vert H \vert X(t) > $','Interpreter','latex','FontSize',14)

subplot(2,3,3); plot(0:dt:T_end,abs(Norm_con),'Linewidth',2)
title('Parallel BUG')
legend('$ \vert \vert X(t) \vert \vert$','Interpreter','latex','FontSize',14)
subplot(2,3,6); plot(0:dt:T_end,abs(Energy_con),'Linewidth',2)
legend('$ < X(t) \vert H \vert X(t) > $','Interpreter','latex','FontSize',14)

figure(3)
plot(0:dt:T_end,rank_max_ad,'Linewidth',2)
hold on
plot(0:dt:T_end,rank_max_BUG,'Linewidth',2)
plot(0:dt:T_end,rank_max,'Linewidth',2)
legend('Rank-adaptive BUG','Fixed-rank BUG','Parallel BUG','FontSize', 20)

% profile off
% profile viewer