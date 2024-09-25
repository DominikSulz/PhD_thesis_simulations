clear; clc; close all;

addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\unconventional_TTN_integrator')
addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\rank_adaptive_integrator_for_TTN')
addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\parallel_TTN_integrator')
addpath(genpath('tensor_toolbox-master'))
addpath(genpath('tensorlab_2016-03-28'))
addpath(genpath('hm-toolbox-master'))
addpath('C:\Users\Dominik\Documents\MATLAB\Matlab toolboxes\hm-toolbox-master')
addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\TTNO\HSS_and_no_structure_case')


% number particles
d = 8;
l = log(d)/log(2); % number of layers
n = 2;             % physical dimension

% parameters of the model
nu = 2;
Delta = 1;
Omega = 1;
alpha = 1;
% time step size and final time
T_end = 1;
dt = 0.005;
% rank of initial data at bottom layer
r = 2;
% bond dimension operator
r_op_max = 30;
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
n_Pu=[1,0;0,0];    %% Projector onto the excited state Pu=(sz+id)/2;

% create initial data binary tree
[X_bin,tau_bin] = init_spin_all_dim_same_rank(r,2,d);

% create initial data tensor train
[X_TT,tau_TT] = init_spin_all_dim_same_rank_TT(d);

%% construction of the TTNO
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

k = hssrank(H_int); 

% construction of TTNO with help of HSS
TTNO_single = TTNO_HSS_construct_single(H_single,A_single,l,l,n*ones(1,d),1:d);
TTNO_int = TTNO_HSS_construct(H_int,A_int,l,l,n*ones(1,d),1:d);

A = Add_TTN(TTNO_single,TTNO_int,tau_bin);
A = rounding(A,tau_bin);
A{end} = -1i*A{end};

% construction of TT represntation with help of HSS
cluster = hss();
cluster = cluster_rec_HSS(cluster,d);
cluster.topnode = 1;

H_single_TT = hss(V_single,'cluster',cluster);
H_int_TT = hss(V_int,'cluster',cluster);

TTNO_TT_single = TTNO_HSS_construct_single_TT(H_single_TT,A_single,d-1,d-1,n*ones(1,d),1:d);
TTNO_TT_int = TTNO_HSS_construct_TT(H_int_TT,A_int,d-1,d-1,n*ones(1,d),1:d);

A_TT = Add_TTN(TTNO_TT_single,TTNO_TT_int,tau_TT);
A_TT = rounding(A_TT,tau_TT);
A_TT{end} = -1i*A_TT{end};

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
    
    % rank-adaptive integrator
    % time integration balanced binary tree
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
  
    % setting for next time step
    X_start_bin = X_new_bin;
    X_start_TT = X_new_TT;
    
    X_start_bin_par = X_new_bin_par;
    X_start_TT_par = X_new_TT_par;
    
    % memory
    S = whos('X_new_bin');
    mem(it) = S.bytes;
    S = whos('X_new_TT');
    mem_TT(it) = S.bytes;
    
    S = whos('X_new_bin_par');
    mem_par(it) = S.bytes;
    S = whos('X_new_TT_par');
    mem_TT_par(it) = S.bytes;
    
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
    [obs_exact,psi] = ref_sol(d,Omega,nu,Delta,alpha,dt,T_end);
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

function[cluster] = cluster_rec_HSS(cluster,d)

cluster.leafnode = 0;
cluster.topnode = 0;
cluster.ml = 1; cluster.nl = 1;
cluster.mr = d-1; cluster.nr = d-1;

cluster.A11 = hss();
cluster.A11.D = 0;
cluster.A11.leafnode = 1;
cluster.A11.topnode = 0;
cluster.A22.topnode = 0;

if d==2
    cluster.A22 = hss();
    cluster.A22.D = 0;
    cluster.A22.leafnode = 1;
else
    d = d-1;
    cluster.A22 = hss(); cluster.A22 = cluster_rec_HSS(cluster.A22,d);
    cluster.A22.D = zeros(d);
end

end