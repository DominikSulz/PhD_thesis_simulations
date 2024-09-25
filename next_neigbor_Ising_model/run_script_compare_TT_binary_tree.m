clear; clc; close all;

addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\unconventional_TTN_integrator')
addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\rank_adaptive_integrator_for_TTN')
addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\parallel_TTN_integrator')
addpath(genpath('tensor_toolbox-master'))
addpath(genpath('tensorlab_2016-03-28'))
addpath(genpath('hm-toolbox-master'))

% number particles
d = 8;
% parameters of the model
Omega = 1;
% time step size and final time
T_end = 1;
dt = 0.01;
% rank of initial data at bottom layer
r = 2;
% bond dimension operator
r_op_max = 20;
r_op_min = 3;

% for rank-adaptive integrator
tol = 10^-8;
r_min = 2;
r_max = 150;

%%% time-step
Iter=T_end/dt;

% initialisations
sx=[0,1;1,0];  %% Pauli Matrix x
sy=[0,-1i;1i,0]; %% Pauli Matrix y
sz=[1,0;0,-1]; %% Pauli Matrix z

% create initial data binary tree
[X_bin,tau_bin] = init_spin_all_dim_same_rank(r,2,d);

% create initial data tensor train
[X_TT,tau_TT] = init_spin_all_dim_same_rank_TT(d);

% create cell array for Hamiltonian
B = linearisation_Ising(Omega*sx,sz,d);

% make operator of Ising model in binary tree representation
A = make_operator(X_bin,B,tau_bin,2*ones(d,1));
A{end} = -1i*A{end};
A = rounding(A,tau_bin);
A = truncate(A,10^-14,r_op_max,r_op_min);

% make operator of Ising model in tensor train representation
A_TT = make_operator(X_TT,B,tau_TT,2*ones(d,1));
A_TT{end} = -1i*A_TT{end};
A_TT = rounding(A_TT,tau_TT);
A_TT = truncate(A_TT,10^-14,r_op_max,r_op_min);

% make operator magnetisation 
B2 = cell(d,d);
for ii=1:d
    B2{ii,ii} = sz;
end
O_Mag_bin = make_operator(X_bin,B2,tau_bin,2*ones(d,1)); % binary tree
O_Mag_bin = rounding(O_Mag_bin,tau_bin);
O_Mag_TT = make_operator(X_TT,B2,tau_TT,2*ones(d,1)); % tensor train
O_Mag_TT = rounding(O_Mag_TT,tau_TT);

tmp = apply_operator_nonglobal(X_bin,O_Mag_bin,d);
obs_sz_bin = (1/d)*Mat0Mat0(X_bin,tmp);
tmp = apply_operator_nonglobal(X_TT,O_Mag_TT,d);
obs_sz_TT = (1/d)*Mat0Mat0(X_TT,tmp);
Diff_mag = abs(obs_sz_bin - obs_sz_TT);

tmp = apply_operator_nonglobal(X_bin,O_Mag_bin,d);
obs_sz_bin_par = (1/d)*Mat0Mat0(X_bin,tmp);
tmp = apply_operator_nonglobal(X_TT,O_Mag_TT,d);
obs_sz_TT_par = (1/d)*Mat0Mat0(X_TT,tmp);
Diff_mag_par = abs(obs_sz_bin_par - obs_sz_TT_par);

% time evolution
time = [];
obs_sz = [];
time_bin = [];
time_TT = [];

% max rank
r_max_bin = max_rank(X_bin);
r_max_TT = max_rank(X_TT);

r_max_bin_par = max_rank(X_bin);
r_max_TT_par = max_rank(X_TT);

X_start_bin = X_bin;
X_start_TT = X_TT;

X_start_bin_par = X_bin;
X_start_TT_par = X_TT;

for it=1:Iter
    t0 = (it-1)*dt;
    t1 = it*dt;
    time(it)=t1;

    % rank-adaptive 
    % time integration with balanced binary tree
    X_new_bin = TTN_integrator_complex_rank_adapt_nonglobal_spin(tau_bin,X_start_bin,@F_Ising,t0,t1,A,d,r_min);
    X_new_bin = truncate(X_new_bin,tol,r_max,r_min);
  
    % time integration with tensor train
    X_new_TT = TTN_integrator_complex_rank_adapt_nonglobal_spin(tau_TT,X_start_TT,@F_Ising,t0,t1,A_TT,d,r_min);
    X_new_TT = truncate(X_new_TT,tol,r_max,r_min);
    
    % parallel integrator
    % time integration balanced binary tree
    X_new_bin_par = TTN_integrator_complex_parallel_nonglobal(tau_bin,X_start_bin_par,@F_Ising,t0,t1,A,d,r_min);
    X_new_bin_par = truncate(X_new_bin_par,tol,r_max,r_min);
  
    % time integration with tensor train
    X_new_TT_par = TTN_integrator_complex_parallel_nonglobal(tau_TT,X_start_TT_par,@F_Ising,t0,t1,A_TT,d,r_min);
    X_new_TT_par = truncate(X_new_TT_par,tol,r_max,r_min);
    
    % memory
    S = whos('X_new_bin');
    mem(it) = S.bytes;
    S = whos('X_new_TT');
    mem_TT(it) = S.bytes;
    
    S = whos('X_new_bin_par');
    mem_par(it) = S.bytes;
    S = whos('X_new_TT_par');
    mem_TT_par(it) = S.bytes;
  
    % setting for next time step
    X_start_bin = X_new_bin;
    X_start_TT = X_new_TT;
    
    X_start_bin_par = X_new_bin_par;
    X_start_TT_par = X_new_TT_par;
    
    % compute magnetisation
    tmp = apply_operator_nonglobal(X_new_bin,O_Mag_bin,d);
    obs_sz_bin(it+1) = (1/d)*Mat0Mat0(X_new_bin,tmp);
    tmp = apply_operator_nonglobal(X_new_TT,O_Mag_TT,d);
    obs_sz_TT(it+1) = (1/d)*Mat0Mat0(X_new_TT,tmp);
    Diff_mag(it+1) = abs(obs_sz_bin(it+1) - obs_sz_TT(it+1));
    
    tmp = apply_operator_nonglobal(X_new_bin_par,O_Mag_bin,d);
    obs_sz_bin_par(it+1) = (1/d)*Mat0Mat0(X_new_bin_par,tmp);
    tmp = apply_operator_nonglobal(X_new_TT_par,O_Mag_TT,d);
    obs_sz_TT_par(it+1) = (1/d)*Mat0Mat0(X_new_TT_par,tmp);
    Diff_mag_par(it+1) = abs(obs_sz_bin_par(it+1) - obs_sz_TT_par(it+1));
    
    % max. rank
    r_max_bin(it+1) = max_rank(X_new_bin);
    r_max_TT(it+1) = max_rank(X_new_TT);
    
    r_max_bin_par(it+1) = max_rank(X_new_bin_par);
    r_max_TT_par(it+1) = max_rank(X_new_TT_par);
    
    fprintf('t=%f, Diff_mag=%f \n', it*dt,  Diff_mag(it+1))
    
end
figure(1)
subplot(1,2,1)
plot(0:dt:T_end,r_max_bin,'o-','Linewidth',2)
hold on
plot(0:dt:T_end,r_max_TT,'*-','Linewidth',2)
legend('Balaced binary tree','Tensor Train','FontSize',20)
title('Rank-adaptive BUG')
xlabel('Time')

subplot(1,2,2)
plot(0:dt:T_end,r_max_bin_par,'o-','Linewidth',2)
hold on
plot(0:dt:T_end,r_max_TT_par,'*-','Linewidth',2)
legend('Balaced binary tree','Tensor Train','FontSize',20)
title('Parallel BUG')
xlabel('Time')

figure(2)
subplot(1,2,1)
plot(dt:dt:T_end,mem,'o-','Linewidth',2)
hold on
plot(dt:dt:T_end,mem_TT,'*-','Linewidth',2)
legend('Balaced binary tree','Tensor Train')
title('Rank-adaptive BUG')
xlabel('Time')

subplot(1,2,2)
plot(dt:dt:T_end,mem_par,'o-','Linewidth',2)
hold on
plot(dt:dt:T_end,mem_TT_par,'*-','Linewidth',2)
legend('Balaced binary tree','Tensor Train')
title('Parallel BUG')
xlabel('Time')

if d<=8
    max_Diff_mag = max(Diff_mag);
    [obs_exact,psi] = ref_sol(d,Omega,dt,T_end);
    figure(3)
    plot(dt:dt:T_end,obs_exact,'Linewidth',2)
    hold on
    plot(0:dt:T_end,obs_sz_bin,'Linewidth',2)
    plot(0:dt:T_end,obs_sz_TT,'Linewidth',2)
    legend('Exact','Balaced binary tree','Tensor Train','FontSize',20)
    figure(4)
    plot(dt:dt:T_end,abs(obs_exact - obs_sz_bin(2:end)),'Linewidth',2)
    hold on
    plot(dt:dt:T_end,abs(obs_exact - obs_sz_TT(2:end)),'Linewidth',2)
end



function [obs_exact,psi] = ref_sol(d,Omega,dt,T_end)

%%%% dimension
L=d;
h=Omega;

FinalTime=T_end;
Iter=FinalTime/dt;


sx=[0,1;1,0];
sy=[0,-i;i,0];
sz=[1,0;0,-1];

iden=[1 0;0 1];




Proj1=[1 0;0 0];
Proj2=[0 0;0 1];

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

H_NN=0*X{1,1};
for i1=1:L-1
    H_NN=H_NN+Z{i1,1}*Z{i1+1,1};
end

H_sys=-h*H_X-H_NN;

psi=zeros(2^L,1);
psi(1,1)=1;

tic
Vt=expm(-1i*dt*H_sys); % Vt=expm(-1i*dt*H_sys);
toc

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

% function [r] = max_rank(X)
% 
% m = length(X) - 2;
% s = size(X{end});
% r = max(s);
% 
% for ii=1:m
%     if 1==iscell(X{ii})
%         r_i = max_rank(X{ii});
%         if r_i > r
%             r = r_i;
%         end
%     end
% end
% 
% end
