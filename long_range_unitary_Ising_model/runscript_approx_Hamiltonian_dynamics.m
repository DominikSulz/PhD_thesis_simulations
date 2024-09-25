clear all; close all; clc
addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\unconventional_TTN_integrator')
addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\parallel_TTN_integrator')
addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\rank_adaptive_integrator_for_TTN')
addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\TTNO\HSS_and_no_structure_case')
addpath('C:\Users\Dominik\Documents\MATLAB\Matlab toolboxes\hm-toolbox-master')

% number particles
d = 4;
l = log(d)/log(2); % number of layers
n = 2;             % physical dimension

hss_tols = [10^-2 10^-4 10^-6 10^-8 10^-10 10^-12];

% parameters of the model
nu = 2;
Delta = 1;
Omega = 1;
alpha = 1;
% time step size and final time
T_end = 1;
dt = 0.01;
% rank of initial data at bottom layer
r = 2;
% bond dimension operator
r_op_max = 40;
r_op_min = 2;

% for rank-adaptive integrator
tol = 10^-8;
r_min = 2;
r_max = 30;

% constant for step rejection
const = 10;

%%% time-step
Iter=T_end/dt;

% initialisations
sx=[0,1;1,0];  %% Pauli Matrix x
sy=[0,-1i;1i,0]; %% Pauli Matrix y
sz=[1,0;0,-1]; %% Pauli Matrix z
n_Pu=[1,0;0,0];    %% Projector onto the excited state Pu=(sz+id)/2;

% % % create initial data binary tree
% [X,tau] = init_spin_all_dim_diff_rank(r,2,d);
% X = rounding(X,tau);
% X = truncate(X,10^-10,30,2);
% [X_test,tau_test] = init_spin_all_dim_diff_rank1(1,d);
% [X,tau] = init_spin_all_dim_diff_rank2(d);
% X = rounding(X,tau);
% [X,tau] = init_spin_all_dim_rank2(r,2,d);
r_vec = [r_max 16 4 2];  % r_vec = [r_max 16 4 2];          % adjust to each TTN!
[X,tau] = init_unconventional(r_vec,2,d);


%% construction of the reference TTNO
hssoption('threshold',10^-16);
% interaction matrix - single particle
A_single = cell(1,d);
for ii=1:d
    A_single{ii} = Omega*sx + Delta*n_Pu;
end
V_single = eye(d,d);

% interaction matrix - long-range interactions
V_int = zeros(d,d);
A_int = cell(1,d);
for ii=1:d
    for jj=1:d
        if ii ~= jj
            V_int(ii,jj) = nu*(1/abs(ii-jj))^alpha;
        end
    end
    A_int{ii} = n_Pu;
end

% HSS
H_single = hss(V_single,'cluster',1:d);
H_single = adjust_H(H_single);

H_int = hss(V_int,'cluster',1:d);
H_int = adjust_H(H_int);

k_full = hssrank(H_int) + 1;
k_save = [];

% construction of TTNO with help of HSS
TTNO_single = TTNO_HSS_construct_single(H_single,A_single,l,l,n*ones(1,d),1:d);
TTNO_int = TTNO_HSS_construct(H_int,A_int,l,l,n*ones(1,d),1:d);

A = Add_TTN(TTNO_single,TTNO_int,tau);
A = rounding(A,tau);
A{end} = -1i*A{end};

%% make operator magnetisation
B2 = cell(d,d);
for ii=1:d
    B2{ii,ii} = sz;
end
Mag = make_operator(X,B2,tau,2*ones(d,1));

%% time evolution
time = [];
obs_sz = [];
obs_sz_ad = [];
obs_sz_BUG = [];
obs_sz_eps = zeros(length(hss_tols),Iter);
obs_sz_ad_eps = zeros(length(hss_tols),Iter);
obs_sz_BUG_eps = zeros(length(hss_tols),Iter);

err_par = zeros(length(hss_tols),Iter);
err_ad = zeros(length(hss_tols),Iter);
err_BUG = zeros(length(hss_tols),Iter);

%% \eps-approximated computations
for kk=1:length(hss_tols)
    % construct the approximated TTNO
    hssoption('threshold',hss_tols(kk));
    % HSS
    H_single_eps = hss(V_single,'cluster',1:d);
    H_single_eps = adjust_H(H_single_eps);
    
    H_int_eps = hss(V_int,'cluster',1:d);
    H_int_eps = adjust_H(H_int_eps);
    
    k_save(kk) = hssrank(H_int_eps) + 1;
    
    % construction of TTNO with help of HSS
    TTNO_single_eps = TTNO_HSS_construct_single(H_single_eps,A_single,l,l,n*ones(1,d),1:d);
    TTNO_int_eps = TTNO_HSS_construct(H_int_eps,A_int,l,l,n*ones(1,d),1:d);
    
    A_eps = Add_TTN(TTNO_single_eps,TTNO_int_eps,tau);
    A_eps = rounding(A_eps,tau);
    A_eps{end} = -1i*A_eps{end};
    %%%%%
    
    % set inital values for dynamics
    X_start = X;
    X_start_ad = X;
    X_start_BUG = X;
    X_start_eps = X;
    X_start_ad_eps = X;
    X_start_BUG_eps = X;
    
    for it=1:Iter
        t0 = (it-1)*dt;
        t1 = it*dt;
        time(it) = t1;
        
        %% reference solution
        % parallel without step rejection
        [X_new,~,~] = TTN_integrator_complex_parallel_nonglobal(tau,X_start,@F_Ising,t0,t1,A,d,r_min);
        X_new = truncate(X_new,tol,r_max,r_min);
        % rank-adaptive
        X_new_ad = TTN_integrator_complex_rank_adapt_nonglobal_spin(tau,X_start_ad,@F_Ising,t0,t1,A,d,r_min);
        X_new_ad = truncate(X_new_ad,tol,r_max,r_min);
        % fixed-rank BUG
        X_new_BUG = TTN_integrator_complex_nonglobal_spin(tau,X_start_BUG,@F_Ising,t0,t1,A,d);
        
        %% approximated dynamics
        % parallel without step rejection
        [X_new_eps,~,~] = TTN_integrator_complex_parallel_nonglobal(tau,X_start_eps,@F_Ising,t0,t1,A_eps,d,r_min);
        X_new_eps = truncate(X_new_eps,tol,r_max,r_min);
        % rank-adaptive
        X_new_ad_eps = TTN_integrator_complex_rank_adapt_nonglobal_spin(tau,X_start_ad_eps,@F_Ising,t0,t1,A_eps,d,r_min);
        X_new_ad_eps = truncate(X_new_ad_eps,tol,r_max,r_min);
        % fixed-rank BUG
        X_new_BUG_eps = TTN_integrator_complex_nonglobal_spin(tau,X_start_BUG_eps,@F_Ising,t0,t1,A_eps,d);
        
        % setting for next time step
        X_start = X_new;
        X_start_ad = X_new_ad;
        X_start_BUG = X_new_BUG;
        X_start_eps = X_new_eps;
        X_start_ad_eps = X_new_ad_eps;
        X_start_BUG_eps = X_new_BUG_eps;
        
        
        % compute difference between solutions
        tmp = X_new_eps;
        tmp{end} = -tmp{end};
        E = Add_TTN(X_new,tmp,tau);
        err_par(kk,it) = abs(sqrt(Mat0Mat0(E,E)));
        
        tmp = X_new_ad_eps;
        tmp{end} = -tmp{end};
        E = Add_TTN(X_new_ad,tmp,tau);
        err_ad(kk,it) = abs(sqrt(Mat0Mat0(E,E)));
        
        tmp = X_new_BUG_eps;
        tmp{end} = -tmp{end};
        E = Add_TTN(X_new_BUG,tmp,tau);
        err_BUG(kk,it) = abs(sqrt(Mat0Mat0(E,E)));
        
        % compute magnetisation in z-direction
        tmp = apply_operator_nonglobal(X_new,Mag,d);
        obs_sz(it) = (1/d)*Mat0Mat0(X_new,tmp);
        
        tmp = apply_operator_nonglobal(X_new_ad,Mag,d);
        obs_sz_ad(it) = (1/d)*Mat0Mat0(X_new_ad,tmp);
        
        tmp = apply_operator_nonglobal(X_new_BUG,Mag,d);
        obs_sz_BUG(it) = (1/d)*Mat0Mat0(X_new_BUG,tmp);
        
        % approximated observables
        tmp = apply_operator_nonglobal(X_new_eps,Mag,d);
        obs_sz_eps(kk,it) = (1/d)*Mat0Mat0(X_new_eps,tmp);
        
        tmp = apply_operator_nonglobal(X_new_ad_eps,Mag,d);
        obs_sz_ad_eps(kk,it) = (1/d)*Mat0Mat0(X_new_ad_eps,tmp);
        
        tmp = apply_operator_nonglobal(X_new_BUG_eps,Mag,d);
        obs_sz_BUG_eps(kk,it) = (1/d)*Mat0Mat0(X_new_BUG_eps,tmp);
        
        
        it
        
    end

end

hssoption('threshold',10^-12) % set all hss option values back to default

if d<=8
    [obs_ref,sol] = ref_sol(d,Omega,nu,Delta,alpha,dt,T_end);
end

%% observable

markers = {'o','*','x','+','.','s','d','^','v','>','<','p','h'};

figure(1)
subplot(1,3,1)
plot(time,obs_sz_ad,'Linewidth',2)
hold on
for jj=1:length(hss_tols)
    plot(time,obs_sz_ad_eps(jj,:),'Marker',markers{jj},'Linewidth',2,'MarkerSize',6)
end
xlabel('Time')
title('Rank-adaptive BUG')
legend('Exact Hamiltonian','10^{-2}','10^{-4}','10^{-6}','10^{-8}','10^{-10}','10^{-12}')

subplot(1,3,2)
plot(time,obs_sz_BUG,'Linewidth',2)
hold on
for jj=1:length(hss_tols)
    plot(time,obs_sz_BUG_eps(jj,:),'Marker',markers{jj},'Linewidth',2,'MarkerSize',6)
end
xlabel('Time')
title('Fixed-rank BUG')
legend('Exact Hamiltonian','10^{-2}','10^{-4}','10^{-6}','10^{-8}','10^{-10}','10^{-12}')

subplot(1,3,3)
plot(time,obs_sz,'Linewidth',2)
hold on
for jj=1:length(hss_tols)
    plot(time,obs_sz_eps(jj,:),'Marker',markers{jj},'Linewidth',2,'MarkerSize',6)
end
xlabel('Time')
title('Parallel BUG')
legend('Exact Hamiltonian','10^{-2}','10^{-4}','10^{-6}','10^{-8}','10^{-10}','10^{-12}')

%% error observable (seems like for tol=10^-12 the error is 0 all the time)
figure(2)
subplot(1,3,1)
for jj=1:length(hss_tols)
    semilogy(time,abs(obs_sz_ad_eps(jj,:) - obs_sz_ad),'Marker',markers{jj},'Linewidth',2,'MarkerSize',6)
    hold on
end
xlabel('Time')
title('Rank-adaptive BUG')
legend('10^{-2}','10^{-4}','10^{-6}','10^{-8}','10^{-10}','10^{-12}')

subplot(1,3,2)
for jj=1:length(hss_tols)
    semilogy(time,abs(obs_sz_BUG_eps(jj,:) - obs_sz_BUG),'Marker',markers{jj},'Linewidth',2,'MarkerSize',6)
    hold on
end
xlabel('Time')
title('Fixed-rank BUG')
legend('10^{-2}','10^{-4}','10^{-6}','10^{-8}','10^{-10}','10^{-12}')

subplot(1,3,3)
for jj=1:length(hss_tols)
    semilogy(time,abs(obs_sz_eps(jj,:) - obs_sz),'Marker',markers{jj},'Linewidth',2,'MarkerSize',6)
    hold on
end
xlabel('Time')
title('Parallel BUG')
legend('10^{-2}','10^{-4}','10^{-6}','10^{-8}','10^{-10}','10^{-12}')

%% error over time to ref solution
figure(3)
subplot(1,3,1)
for jj=1:length(hss_tols)
    semilogy(time,err_ad(jj,:),'Marker',markers{jj},'Linewidth',2,'MarkerSize',6)
    hold on
end
xlabel('Time')
title('Rank-adaptive BUG')
legend('10^{-2}','10^{-4}','10^{-6}','10^{-8}','10^{-10}','10^{-12}')

subplot(1,3,2)
for jj=1:length(hss_tols)
    semilogy(time,err_BUG(jj,:),'Marker',markers{jj},'Linewidth',2,'MarkerSize',6)
    hold on
end
xlabel('Time')
title('Fixed-rank BUG')
legend('10^{-2}','10^{-4}','10^{-6}','10^{-8}','10^{-10}','10^{-12}')

subplot(1,3,3)
for jj=1:length(hss_tols)
    semilogy(time,err_par(jj,:),'Marker',markers{jj},'Linewidth',2,'MarkerSize',6)
    hold on
end
xlabel('Time')
title('Parallel BUG')
legend('10^{-2}','10^{-4}','10^{-6}','10^{-8}','10^{-10}','10^{-12}')



function [obs_exact,psi] = ref_sol(d,Omega,nu,Delta,alpha,dt,T_end)

%%%% dimension
L=d;
h=Omega;

FinalTime=T_end;
Iter=FinalTime/dt;


sx=[0,1;1,0];
sy=[0,-1i;1i,0];
sz=[1,0;0,-1];
n_Pu=[1,0;0,0];

B = linearisation_long_range_unitary_full(L,sx,n_Pu,nu,Delta,Omega,alpha);

H_sys = superkron(L,B);

psi=zeros(2^L,1);
psi(1,1)=1;

tmp2=diag(sparse(ones(2^(L-1),1)));
X{1,1}=kron(sx,tmp2);
X{L,1}=kron(tmp2,sx);

Z{1,1}=kron(sz,tmp2);
Z{L,1}=kron(tmp2,sz);


for i1=2:L-1
    tmp1=diag(sparse(ones(2^(i1-1),1)));
    tmp2=diag(sparse(ones(2^(L-i1),1)));
    X{i1,1}=kron(kron(tmp1,sx),tmp2);
    Z{i1,1}=kron(kron(tmp1,sz),tmp2);
end

Mag=0*X{1,1};
H_X=0*X{1,1};

for i1=1:L
    Mag=Mag+Z{i1,1};
    H_X=H_X+X{i1,1};
end

Vt=expm(-1i*dt*H_sys); % Vt=expm(-1i*dt*H_sys);

for tst=1:Iter
    time(tst)=tst*dt;
    psi=Vt*psi;
    psi=psi/norm(psi);
    
    p1(tst)=dot(psi,Mag*psi)/L;
    energy(tst)=dot(psi,H_sys*psi)/L;
    
    mag(tst)=dot(psi,Mag*psi)/L;
end


%%%% diagonalization

Ground_State_Energy=min(real(eig(H_sys)))/L;
%%%%% diagonalization of H_eff
[v2,d2]=eig(full(H_sys));
[x2,y2]=sort(diag(d2),'ascend');


%%% positive eigenvalue
ind_p=y2(end);


norm(H_sys*v2(:,ind_p)-d2(ind_p,ind_p)*v2(:,ind_p));
Ground_state=v2(:,ind_p);
Ground_state=Ground_state/norm(Ground_state);


Mag_Ground_state=dot(Ground_state,Mag*Ground_state)/L;


obs_exact = mag;

end


function [M] = superkron(L,B)
% B is a cell array
[m,n] = size(B);

M = sparse(2^L,2^L);
for jj=1:m
    tmp = 1;
    for ii=n:-1:1
        if isempty(B{jj,ii}) == 1
            tmp = kron(tmp,speye(2,2));
        else
            tmp = kron(tmp,B{jj,ii});
        end
    end
    M = M + tmp;
end

end