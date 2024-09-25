clear all; close all; clc
addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\unconventional_TTN_integrator')
addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\parallel_TTN_integrator')
addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\rank_adaptive_integrator_for_TTN')
addpath('C:\Users\Dominik\Documents\MATLAB\Matlab toolboxes\hm-toolbox-master')
addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\TTNO\HSS_and_no_structure_case')

% number particles
d = 8;
l = log(d)/log(2); % number of layers
% parameters of the model
Omega = 1;

% time step size and final time
T_end = 1;
dt_save = [0.1 0.05 0.025 0.01 0.005 0.0025 0.001];

% rank of initial data at bottom layer
r = 2;
% bond dimension operator
r_op_max = 10;
r_op_min = 2;

% for rank-adaptive integrator
tol = 10^-8;
r_min = 2;
r_max = 20;

% constant for step rejection
const = 10;

% initialisations
sx=[0,1;1,0];  %% Pauli Matrix x
sy=[0,-1i;1i,0]; %% Pauli Matrix y
sz=[1,0;0,-1]; %% Pauli Matrix z

% % % create initial data binary tree
% [X,tau] = init_spin_all_dim_diff_rank(r,2,d);
% [X_test,tau_test] = init_spin_all_dim_diff_rank1(1,d);
% [X,tau] = init_spin_all_dim_diff_rank2(d);
% X = rounding(X,tau);
r_vec = [16 4 2];  % r_vec = [r_max 16 4 2];          % adjust to each TTN! 
[X,tau] = init_unconventional(r_vec,2,d);

% initial data Tucker tensor
% [X,tau] = init_Tucker(d);

% create cell array for Hamiltonian
B = linearisation_Ising(Omega*sx,sz,d);

% make operator of Ising model in TTN representation -> change that to HSS construction
A = make_operator(X,B,tau,2*ones(d,1));
A{end} = -A{end};
A{end} = -1i*A{end};
A = rounding(A,tau);
A = truncate(A,10^-14,r_op_max,r_op_min);

% time evolution
time = [];
% obs_sz = [];
% obs_sz_ad = [];
% obs_sz_BUG = [];

for kk=1:length(dt_save)
    dt = dt_save(kk);
    %%% time-step
    Iter=T_end/dt;
    X_start = X;
    X_start_ad = X;
    X_start_BUG = X;

    tic
    for it=1:Iter
        t0 = (it-1)*dt;
        t1 = it*dt;
        time(it) = t1;
        
        % time integration with parallel integrator
        %     X_new = TTN_integrator_complex_parallel_nonglobal(tau,X_start,@F_Ising,t0,t1,A,d,r_min);
        X_new = TTN_integrator_complex_parallel_nonglobal(tau,X_start,@F_Ising,t0,t1,A,d,r_min);
        X_new = truncate(X_new,tol,r_max,r_min);

%         rej = 1;
%         while rej == 1
%             [X_new,U_tilde,~] = TTN_integrator_complex_parallel_nonglobal(tau,X_start,@F_Ising,t0,t1,A,d,r_min);
%             X_new = truncate(X_new,tol,r_max,r_min);
%             rej_rk = rejection_check(X_start,X_new);
%             U_tilde_SR = U_tilde_TTN(X_start,X_new);
%             rej_eta = eta_check(X_start,X_new,U_tilde_SR,F_Ising(t0,X_new,A,d),const,dt,tol);
%             if (rej_rk == 1) || (rej_eta == 1)
%                 X_start = augment_zero(X_start,X_new); % Kann evtl. paralleler BUG nicht mit 0-Zeilen umghen?
%                 X_start = orthogonalize(X_start);
%             else
%                 rej = 0;
%             end
%         end

        % rank-adaptive
        X_new_ad = TTN_integrator_complex_rank_adapt_nonglobal_spin(tau,X_start_ad,@F_Ising,t0,t1,A,d,r_min);
        X_new_ad = truncate(X_new_ad,tol,r_max,r_min);
        
        % BUG integrator
        X_new_BUG = TTN_integrator_complex_nonglobal_spin(tau,X_start_BUG,@F_Ising,t0,t1,A,d);
        
        % setting for next time step
        X_start = X_new;
        X_start_ad = X_new_ad;
        X_start_BUG = X_new_BUG;
        
    end

    kk
    if d<=8
        [obs_ref,sol] = ref_sol(d,Omega,dt,T_end);
        
        sol_ra_BUG = double(full_tensor(X_new_ad));
        sol_par = double(full_tensor(X_new));
        sol_BUG = double(full_tensor(X_new_BUG));
        err_ad(kk) = norm(sol - sol_ra_BUG(:));
        err_par(kk) = norm(sol - sol_par(:));
        err_BUG(kk) = norm(sol - sol_BUG(:));
    end
    
end

figure(1)
loglog(dt_save,err_ad,'*-','Linewidth',2)
hold on 
loglog(dt_save,err_BUG,'x-','Linewidth',2)
loglog(dt_save,err_par,'o-','Linewidth',2)
loglog(dt_save,0.001*dt_save,'--','Linewidth',2)
% loglog(dt_save,0.001*dt_save.^2,'Linewidth',2)
% loglog(dt_save,0.001*dt_save.^3,'Linewidth',2)
grid on 
set(gca,'YMinorGrid','off')
% title('Error','FontSize', 14)
xlabel('Time step size','FontSize', 20)
ylabel('Error','FontSize', 20)
% legend('Error Parallel Integrator','Error Rank-adaptive BUG','Error fixed-rank BUG','h','h^2','h^3','FontSize', 20)
legend('Rank-adaptive BUG','Fixed-rank BUG','Parallel BUG','h','FontSize', 20)

%save('workspace_errorplot_d4_Tend1')

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



% % TTNO via HSS-ansatz
% %% interaction matrix - single particle
% A_single = cell(1,d);
% for ii=1:d
%     A_single{ii} = -Omega*sx;
% end
% V_single = eye(d,d);
% 
% %% interaction matrix - next-neighbor interaction
% V_int = zeros(d,d);
% A_int = cell(1,d);
% for ii=1:d-1
%     V_int(ii,ii+1) = -1;
%     A_int{ii} = sz;
% end
% A_int{end} = sz;
% 
% %% HSS decomposition of the interaction matrices
% H_single = hss(V_single,'cluster',1:d);
% H_single = adjust_H(H_single);
% 
% H_int = hss(V_int,'cluster',1:d);
% H_int = adjust_H(H_int);
% 
% % construction of TTNs via HSS
% TTNO_single = TTNO_HSS_construct_single(H_single,A_single,l,l,2*ones(1,d),1:d);
% TTNO_int = TTNO_HSS_construct(H_int,A_int,l,l,2*ones(1,d),1:d);
% 
% TTNO = Add_TTN(TTNO_single,TTNO_int,tau);
% TTNO = rounding(TTNO,tau);
